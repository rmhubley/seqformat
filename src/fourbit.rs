//! BWA/BAM-style 4-bit sequence codec.
//!
//! The *nibble encoding* is the de-facto "4-bit" packing of the BWA/samtools
//! ecosystem (htslib `seq_nt16`): each base is one nibble indexing the alphabet
//! `=ACMGRSVTWYHKDBN` (0..=15), high nibble first, so two bases share a byte and
//! every IUB degenerate code (including `N`) is representable directly — no side
//! tables needed.
//!
//! Raw `.pac`/htslib blobs carry no names or framing, so for a self-contained,
//! round-trippable test file we wrap the packed nibbles in a tiny container:
//!
//! ```text
//! header
//!   magic    u32  0x54494234  ("4BIT")
//!   version  u32  1
//!   seqCount u32
//!   reserved u32  0
//! per sequence
//!   nameSize u8
//!   name     nameSize bytes
//!   length   u64        (number of bases)
//!   packed   ceil(length/2) bytes (4 bits/base, first base in the high nibble)
//! ```
//!
//! Note: 4-bit packing does not carry soft-mask (case) information, so a FASTA
//! round-trip through this format is case-insensitive (output is upper-case).

use crate::error::{fmt_err, Result};
use crate::io::{peek_u32, put_u32, put_u64, put_u8, Reader};
use crate::seq::{base_to_nibble, nibble_to_base, Sequence};
use std::fs;
use std::path::Path;

const MAGIC: u32 = 0x5449_4234; // "4BIT" little-endian
const VERSION: u32 = 1;

pub fn pack(bases: &[u8]) -> Vec<u8> {
    let mut packed = vec![0u8; bases.len().div_ceil(2)];
    for (i, &b) in bases.iter().enumerate() {
        let nib = base_to_nibble(b);
        if i % 2 == 0 {
            packed[i / 2] |= nib << 4; // high nibble = first base
        } else {
            packed[i / 2] |= nib;
        }
    }
    packed
}

pub fn unpack(packed: &[u8], len: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(len);
    for i in 0..len {
        let byte = packed[i / 2];
        let nib = if i % 2 == 0 { byte >> 4 } else { byte & 0x0F };
        out.push(nibble_to_base(nib));
    }
    out
}

pub fn to_bytes(seqs: &[Sequence]) -> Result<Vec<u8>> {
    for s in seqs {
        if s.name.len() > 255 {
            return fmt_err(format!(
                "sequence name {:?} is too long ({} > 255 bytes)",
                s.name,
                s.name.len()
            ));
        }
    }
    let mut out = Vec::new();
    put_u32(&mut out, MAGIC);
    put_u32(&mut out, VERSION);
    put_u32(&mut out, seqs.len() as u32);
    put_u32(&mut out, 0);
    for s in seqs {
        put_u8(&mut out, s.name.len() as u8);
        out.extend_from_slice(s.name.as_bytes());
        put_u64(&mut out, s.bases.len() as u64);
        out.extend_from_slice(&pack(&s.bases));
    }
    Ok(out)
}

pub fn from_bytes(data: &[u8]) -> Result<Vec<Sequence>> {
    if data.len() < 16 || peek_u32(data, 0, true) != Some(MAGIC) {
        return fmt_err("not a 4-bit file (bad magic)");
    }
    let mut r = Reader::new(data, true);
    let _magic = r.u32()?;
    let version = r.u32()?;
    if version != VERSION {
        return fmt_err(format!("unsupported 4-bit version {version}"));
    }
    let count = r.u32()? as usize;
    let _reserved = r.u32()?;

    let mut seqs = Vec::with_capacity(count);
    for _ in 0..count {
        let name_len = r.u8()? as usize;
        let name = String::from_utf8_lossy(r.take(name_len)?).into_owned();
        let len = r.u64()? as usize;
        let packed = r.take(len.div_ceil(2))?;
        seqs.push(Sequence::new(name, unpack(packed, len)));
    }
    Ok(seqs)
}

/// Random-access reader for the 4-bit container. The container has no offset
/// table, so `open` scans the records once to build an in-memory index; after
/// that `extract` decodes only the requested window.
pub struct FourBitReader {
    data: Vec<u8>,
    /// name -> (packed_start_byte, length_in_bases)
    by_name: std::collections::HashMap<String, (usize, usize)>,
    order: Vec<String>,
}

impl FourBitReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Self::from_vec(fs::read(path)?)
    }

    pub fn from_vec(data: Vec<u8>) -> Result<Self> {
        if data.len() < 16 || peek_u32(&data, 0, true) != Some(MAGIC) {
            return fmt_err("not a 4-bit file (bad magic)");
        }
        let mut r = Reader::new(&data, true);
        let _magic = r.u32()?;
        let version = r.u32()?;
        if version != VERSION {
            return fmt_err(format!("unsupported 4-bit version {version}"));
        }
        let count = r.u32()? as usize;
        let _reserved = r.u32()?;

        let mut by_name = std::collections::HashMap::with_capacity(count);
        let mut order = Vec::with_capacity(count);
        for _ in 0..count {
            let name_len = r.u8()? as usize;
            let name = String::from_utf8_lossy(r.take(name_len)?).into_owned();
            let len = r.u64()? as usize;
            let packed_start = r.pos();
            r.seek(packed_start + len.div_ceil(2)); // skip packed bytes
            by_name.insert(name.clone(), (packed_start, len));
            order.push(name);
        }
        let _ = r;
        Ok(FourBitReader { data, by_name, order })
    }

    pub fn names(&self) -> &[String] {
        &self.order
    }

    pub fn sequence_infos(&self) -> Vec<(String, usize)> {
        self.order
            .iter()
            .map(|n| (n.clone(), self.by_name[n].1))
            .collect()
    }

    pub fn extract(&self, name: &str, start: usize, end: Option<usize>) -> Result<Vec<u8>> {
        let &(packed_start, len) = self
            .by_name
            .get(name)
            .ok_or_else(|| crate::error::Error::Format(format!("no sequence named {name:?}")))?;
        let end = end.unwrap_or(len).min(len);
        let start = start.min(end);
        let mut out = Vec::with_capacity(end - start);
        for i in start..end {
            let byte = self.data[packed_start + i / 2];
            let nib = if i % 2 == 0 { byte >> 4 } else { byte & 0x0F };
            out.push(nibble_to_base(nib));
        }
        Ok(out)
    }
}

pub fn is_fourbit(data: &[u8]) -> bool {
    data.len() >= 4 && peek_u32(data, 0, true) == Some(MAGIC)
}

pub fn write_file(path: impl AsRef<Path>, seqs: &[Sequence]) -> Result<()> {
    fs::write(path, to_bytes(seqs)?)?;
    Ok(())
}

pub fn read_file(path: impl AsRef<Path>) -> Result<Vec<Sequence>> {
    from_bytes(&fs::read(path)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fourbit_roundtrip_all_codes() {
        // Upper-cased because 4-bit drops case.
        let s = vec![
            Sequence::new("c", b"ACGTRYSWKMBDHVN".to_vec()),
            Sequence::new("d", b"ACGTA".to_vec()), // odd length
        ];
        let bytes = to_bytes(&s).unwrap();
        let out = from_bytes(&bytes).unwrap();
        assert_eq!(out, s);
    }

    #[test]
    fn lowercase_folds_to_upper() {
        let s = vec![Sequence::new("c", b"acgtryswn".to_vec())];
        let out = from_bytes(&to_bytes(&s).unwrap()).unwrap();
        assert_eq!(out[0].bases, b"ACGTRYSWN");
    }
}
