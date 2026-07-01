//! **2be** ("twoBit-edit") — a from-scratch redesign of twoBit (no backward
//! compatibility), addressing the two scaling problems of the original:
//!
//! 1. **High sequence count.** The flat TOC is replaced by an on-disk
//!    [`crate::bptree`] B+ tree (name → record offset): O(log N) lookup with no
//!    full index load.
//! 2. **IUB codes as runs *and* points.** The single N-block table is replaced by
//!    one per-sequence **merged, position-sorted tagged-edit stream** whose
//!    entries are typed: `N_RUN` / `MASK_RUN` (run-length, for gaps & masking),
//!    `IUB_POINT` (an isolated degenerate code), and `IUB_RUN` (a degenerate
//!    run). N gaps stay cheap; isolated IUB cost ~one entry each and never
//!    pollute a separate gap table.
//!
//! Emitting a sequence is just sweeping that one stream alongside the 2-bit
//! array. A tiny per-record index of the run entries makes random sub-range
//! extraction handle straddling runs without scanning past the (possibly many)
//! point entries.
//!
//! ```text
//! HEADER (32 bytes, little-endian)
//!   magic     u32   "2BE1"
//!   version   u16 ; flags u16
//!   seqCount  u64
//!   tocOffset u64   -> B+ tree blob
//!   reserved  u64
//! SEQUENCE RECORD (self-contained; coordinates are sequence-local u32)
//!   dnaSize   u32
//!   editCount u32
//!   edits     editCount × EDIT(10B): pos:u32, len:u32, type:u8, code:u8
//!   runCount  u32
//!   runIndex  runCount × u32      (indices, into edits, of non-point entries)
//!   packedDna ceil(dnaSize/4)     (2-bit, T=00 C=01 A=10 G=11)
//! B+ TREE TOC  (name -> absolute record offset)
//! ```

use crate::bptree;
use crate::error::{fmt_err, Error, Result};
use crate::io::{peek_u32, put_u16, put_u32, put_u64};
use crate::seq::{base_to_twobit, is_iub_degenerate, twobit_to_base, Sequence};
use crate::source::Source;
use std::fs;
use std::path::Path;

const MAGIC: u32 = 0x3145_4232; // "2BE1" little-endian
const HEADER: usize = 32;
const EDIT_SIZE: usize = 10;

// Edit types.
const N_RUN: u8 = 0;
const IUB_POINT: u8 = 1;
const IUB_RUN: u8 = 2;
const MASK_RUN: u8 = 3;

struct Edit {
    pos: u32,
    len: u32,
    etype: u8,
    code: u8,
}

#[derive(PartialEq, Clone, Copy)]
enum Class {
    Plain,
    N,
    Iub(u8),
}

fn base_class(b: u8) -> Class {
    let u = b.to_ascii_uppercase();
    if matches!(u, b'A' | b'C' | b'G' | b'T') {
        Class::Plain
    } else if is_iub_degenerate(u) {
        Class::Iub(u)
    } else {
        Class::N // bare N and any other junk
    }
}

/// Build the merged, position-sorted edit stream for one sequence.
fn build_edits(bases: &[u8]) -> Vec<Edit> {
    let n = bases.len();
    let mut edits = Vec::new();

    // Base layer: runs of equal class. N (and junk) → N_RUN; identical IUB code
    // → IUB_RUN (len>1) or IUB_POINT (len==1); plain ACGT contributes nothing.
    let mut i = 0;
    while i < n {
        let c = base_class(bases[i]);
        let mut j = i + 1;
        while j < n && base_class(bases[j]) == c {
            j += 1;
        }
        match c {
            Class::Plain => {}
            Class::N => edits.push(Edit { pos: i as u32, len: (j - i) as u32, etype: N_RUN, code: 0 }),
            Class::Iub(code) => {
                let (etype, _) = if j - i == 1 { (IUB_POINT, ()) } else { (IUB_RUN, ()) };
                edits.push(Edit { pos: i as u32, len: (j - i) as u32, etype, code });
            }
        }
        i = j;
    }

    // Case layer: runs of lower-case (soft-mask), independent of base identity.
    let mut i = 0;
    while i < n {
        if bases[i].is_ascii_lowercase() {
            let mut j = i + 1;
            while j < n && bases[j].is_ascii_lowercase() {
                j += 1;
            }
            edits.push(Edit { pos: i as u32, len: (j - i) as u32, etype: MASK_RUN, code: 0 });
            i = j;
        } else {
            i += 1;
        }
    }

    edits.sort_by_key(|e| (e.pos, e.etype));
    edits
}

fn pack_dna(bases: &[u8]) -> Vec<u8> {
    let mut packed = vec![0u8; bases.len().div_ceil(4)];
    for (i, &b) in bases.iter().enumerate() {
        packed[i / 4] |= base_to_twobit(b) << (6 - 2 * (i % 4));
    }
    packed
}

fn encode_record(bases: &[u8]) -> Vec<u8> {
    let edits = build_edits(bases);
    // Run index: every non-point edit (the few entries that can straddle a
    // query start), so extraction can find straddlers without walking points.
    let run_idx: Vec<u32> = edits
        .iter()
        .enumerate()
        .filter(|(_, e)| e.etype != IUB_POINT)
        .map(|(i, _)| i as u32)
        .collect();

    let mut out = Vec::new();
    put_u32(&mut out, bases.len() as u32);
    put_u32(&mut out, edits.len() as u32);
    for e in &edits {
        put_u32(&mut out, e.pos);
        put_u32(&mut out, e.len);
        out.push(e.etype);
        out.push(e.code);
    }
    put_u32(&mut out, run_idx.len() as u32);
    for &ri in &run_idx {
        put_u32(&mut out, ri);
    }
    out.extend_from_slice(&pack_dna(bases));
    out
}

/// Serialize sequences to a 2be byte image.
pub fn to_bytes(seqs: &[Sequence]) -> Result<Vec<u8>> {
    let records: Vec<Vec<u8>> = seqs.iter().map(|s| encode_record(&s.bases)).collect();

    let mut offset = HEADER;
    let mut toc_items = Vec::with_capacity(seqs.len());
    for (s, r) in seqs.iter().zip(&records) {
        toc_items.push((s.name.clone(), offset as u64));
        offset += r.len();
    }
    let toc_offset = offset;
    let blob = bptree::build(toc_items);

    let mut out = Vec::with_capacity(toc_offset + blob.len());
    put_u32(&mut out, MAGIC);
    put_u16(&mut out, 1); // version
    put_u16(&mut out, 0); // flags
    put_u64(&mut out, seqs.len() as u64);
    put_u64(&mut out, toc_offset as u64);
    put_u64(&mut out, 0); // reserved
    for r in &records {
        out.extend_from_slice(r);
    }
    out.extend_from_slice(&blob);
    Ok(out)
}

pub fn is_twobyte(data: &[u8]) -> bool {
    data.len() >= HEADER && peek_u32(data, 0, true) == Some(MAGIC)
}

/// Random-access reader over a 2be file (local file, memory, or remote HTTP).
pub struct TwoByteReader {
    src: Source,
    toc_offset: usize,
}

impl TwoByteReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Self::from_source(Source::from_path(path)?)
    }

    pub fn from_vec(data: Vec<u8>) -> Result<Self> {
        Self::from_source(Source::Mem(data))
    }

    /// Open a remote `http(s)://` 2be over HTTP range reads (UDC-style). The B+
    /// tree resolves a name in ~3 range reads; see [`TwoByteReader::http_stats`].
    pub fn open_url(url: &str) -> Result<Self> {
        Self::from_source(Source::from_url(url)?)
    }

    pub(crate) fn from_source(src: Source) -> Result<Self> {
        if src.len() < HEADER || src.u32_at(0, true)? != MAGIC {
            return fmt_err("not a 2be file (bad magic)");
        }
        let toc_offset = src.u64_at(16, true)? as usize;
        if toc_offset > src.len() {
            return fmt_err("2be: TOC offset past end of file");
        }
        Ok(TwoByteReader { src, toc_offset })
    }

    /// `(requests, bytes)` issued so far over HTTP; `None` for local/in-memory.
    pub fn http_stats(&self) -> Option<(u64, u64)> {
        self.src.http_stats()
    }

    pub fn names(&self) -> Vec<String> {
        match bptree::iter_all_src(&self.src, self.toc_offset) {
            Ok(v) => v.into_iter().map(|(n, _)| n).collect(),
            Err(_) => Vec::new(),
        }
    }

    /// `(name, dnaSize, editCount)` per sequence, in name order.
    pub fn sequence_stats(&self) -> Result<Vec<(String, usize, usize)>> {
        let mut v = Vec::new();
        for (name, off) in bptree::iter_all_src(&self.src, self.toc_offset)? {
            let off = off as usize;
            let dna = self.src.u32_at(off, true)? as usize;
            let edits = self.src.u32_at(off + 4, true)? as usize;
            v.push((name, dna, edits));
        }
        Ok(v)
    }

    /// Extract `[start, end)` of `name` (0-based, half-open).
    pub fn extract(&self, name: &str, start: usize, end: Option<usize>) -> Result<Vec<u8>> {
        let src = &self.src;
        let rec = bptree::find_src(src, self.toc_offset, name)?
            .ok_or_else(|| Error::Format(format!("no sequence named {name:?}")))?
            as usize;

        let dna_size = src.u32_at(rec, true)? as usize;
        let edit_count = src.u32_at(rec + 4, true)? as usize;
        let edits_off = rec + 8;
        let run_count_off = edits_off + edit_count * EDIT_SIZE;
        let run_count = src.u32_at(run_count_off, true)? as usize;
        let run_idx_off = run_count_off + 4;
        let packed_off = run_idx_off + run_count * 4;

        let end = end.unwrap_or(dna_size).min(dna_size);
        let start = start.min(end);
        if start >= end {
            return Ok(Vec::new());
        }

        // Unpack the window from the 2-bit array (read only the spanning bytes).
        let first_byte = start / 4;
        let last_byte = (end - 1) / 4;
        let packed = src.bytes_at(packed_off + first_byte, last_byte - first_byte + 1)?;
        let mut seq = Vec::with_capacity(end - start);
        for i in start..end {
            let byte = packed[i / 4 - first_byte];
            seq.push(twobit_to_base((byte >> (6 - 2 * (i % 4))) & 0x03));
        }

        // Read one EDIT (10B): (pos, len, type, code).
        let edit_at = |k: usize| -> Result<(usize, usize, u8, u8)> {
            let b = src.bytes_at(edits_off + k * EDIT_SIZE, EDIT_SIZE)?;
            Ok((
                u32::from_le_bytes(b[0..4].try_into().unwrap()) as usize,
                u32::from_le_bytes(b[4..8].try_into().unwrap()) as usize,
                b[8],
                b[9],
            ))
        };

        let mut put = |pos: usize, len: usize, val: u8| {
            let lo = pos.max(start);
            let hi = (pos + len).min(end);
            if lo < hi {
                seq[lo - start..hi - start].fill(val);
            }
        };

        // --- phase 1: base edits ---
        // IUB points that start within the window (binary search the stream by pos).
        let pos_at = |k: usize| -> Result<usize> {
            Ok(src.u32_at(edits_off + k * EDIT_SIZE, true)? as usize)
        };
        let lo = lower_bound_res(edit_count, |k| Ok(pos_at(k)? < start))?;
        let hi = lower_bound_res(edit_count, |k| Ok(pos_at(k)? < end))?;
        for k in lo..hi {
            let (p, _l, t, c) = edit_at(k)?;
            if t == IUB_POINT {
                put(p, 1, c);
            }
        }
        // Base runs overlapping the window (from the small run index).
        for r in 0..run_count {
            let ri = src.u32_at(run_idx_off + r * 4, true)? as usize;
            let (p, l, t, c) = edit_at(ri)?;
            match t {
                N_RUN => put(p, l, b'N'),
                IUB_RUN => put(p, l, c),
                _ => {}
            }
        }

        // --- phase 2: mask runs (lower-case, applied last) ---
        for r in 0..run_count {
            let ri = src.u32_at(run_idx_off + r * 4, true)? as usize;
            let (p, l, t, _c) = edit_at(ri)?;
            if t == MASK_RUN {
                let lo2 = p.max(start);
                let hi2 = (p + l).min(end);
                if lo2 < hi2 {
                    for b in &mut seq[lo2 - start..hi2 - start] {
                        b.make_ascii_lowercase();
                    }
                }
            }
        }

        Ok(seq)
    }

    pub fn read_all(&self) -> Result<Vec<Sequence>> {
        bptree::iter_all_src(&self.src, self.toc_offset)?
            .into_iter()
            .map(|(name, _)| self.extract(&name, 0, None).map(|b| Sequence::new(name, b)))
            .collect()
    }
}

impl crate::SeqReader for TwoByteReader {
    fn names(&self) -> Vec<String> {
        TwoByteReader::names(self)
    }
    fn extract(&self, name: &str, start: usize, end: Option<usize>) -> Result<Vec<u8>> {
        TwoByteReader::extract(self, name, start, end)
    }
    fn http_stats(&self) -> Option<(u64, u64)> {
        self.src.http_stats()
    }
}

/// Fallible binary-search lower bound: first index where `pred` is false.
fn lower_bound_res(n: usize, pred: impl Fn(usize) -> Result<bool>) -> Result<usize> {
    let (mut lo, mut hi) = (0, n);
    while lo < hi {
        let mid = (lo + hi) / 2;
        if pred(mid)? {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    Ok(lo)
}

pub fn write_file(path: impl AsRef<Path>, seqs: &[Sequence]) -> Result<()> {
    fs::write(path, to_bytes(seqs)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_mixed_content() {
        // points, an IUB run (WWW / RRRRR), N run, soft-mask, plain.
        let s = vec![Sequence::new(
            "mix",
            b"ACGTNNNNNNNNacgtRYSWWWKMacgtRRRRRACGT".to_vec(),
        )];
        let rd = TwoByteReader::from_vec(to_bytes(&s).unwrap()).unwrap();
        let whole = rd.extract("mix", 0, None).unwrap();
        assert_eq!(whole, s[0].bases);
        // every sub-window matches a plain slice of the full decode
        let n = s[0].bases.len();
        for a in 0..=n {
            for b in a..=n {
                assert_eq!(rd.extract("mix", a, Some(b)).unwrap(), &whole[a..b], "{a}-{b}");
            }
        }
    }

    #[test]
    fn bptree_many_sequences() {
        // 600 sequences forces a multi-level B+ tree (fan-out 256).
        let s: Vec<Sequence> = (0..600)
            .map(|i| {
                Sequence::new(
                    format!("seq{i:04}"),
                    format!("ACGT{}NNNacgtRY", "ACG".repeat(i % 7 + 1)).into_bytes(),
                )
            })
            .collect();
        let rd = TwoByteReader::from_vec(to_bytes(&s).unwrap()).unwrap();
        assert_eq!(rd.names().len(), 600);
        // random-ish lookups resolve to the right sequence
        for orig in s.iter().step_by(37) {
            assert_eq!(rd.extract(&orig.name, 0, None).unwrap(), orig.bases, "{}", orig.name);
        }
        assert!(rd.extract("nonexistent", 0, None).is_err());
    }
}
