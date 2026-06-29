//! UCSC twoBit codec — standard, long, and the IUB extension.
//!
//! # On-disk layout
//!
//! ```text
//! header (16 bytes)
//!   signature   u32  0x1A412743  (native byte order; doubles as a BOM)
//!   version     u32  0 = standard, 1 = long
//!   seqCount    u32
//!   reserved    u32  0
//! index (one entry per sequence)
//!   nameSize    u8
//!   name        nameSize bytes (ASCII)
//!   offset      u32 (standard) | u64 (long)   -> absolute offset of the record
//! record (at `offset`)
//!   dnaSize      u32
//!   nBlockCount  u32
//!   nBlockStarts u32[nBlockCount]
//!   nBlockSizes  u32[nBlockCount]
//!   maskCount    u32
//!   maskStarts   u32[maskCount]
//!   maskSizes    u32[maskCount]
//!   reserved     u32  0
//!   packedDna    ceil(dnaSize/4) bytes  (T=00,C=01,A=10,G=11; first base hi bits)
//!   ── IUB extension (only when the file carries the trailer below) ──
//!   iubCount     u32
//!   iubStarts    u32[iubCount]
//!   iubSizes     u32[iubCount]
//!   iubCodes     u8[iubCount]   (ASCII IUB letter for each run)
//! ```
//!
//! ## Why the IUB extension is backward compatible
//!
//! A legacy reader locates each record purely through the index `offset` and
//! reads exactly `ceil(dnaSize/4)` packed bytes — it never reads past the DNA of
//! a record. We tuck the per-record IUB table into the gap *after* the packed
//! DNA and *before* the next record (the offsets we write already account for
//! it). The presence of the extension is advertised only by an 8-byte magic
//! **trailer at the very end of the file**, which legacy readers never look at.
//! So an old tool reads the file as an ordinary 2bit, reporting every degenerate
//! position as `N`; our reader additionally recovers the exact IUB code.
//!
//! The IUB table does not duplicate `N`-blocks: a degenerate run is *also*
//! present in the standard N-block table (so old readers mask it as `N`), and
//! the IUB table merely overrides those positions with the precise code. Bare
//! `N` runs are left out of the IUB table entirely.
//!
//! ## The backward-compatible name index (`--index`)
//!
//! The stock 2bit index is a run of *variable-length* `nameSize/name/offset`
//! entries in file order, so a reader that wants a sequence by name must load
//! every entry (UCSC builds a hash; we build a `HashMap`) — O(N) work on every
//! open, painful for files with millions of short sequences.
//!
//! To get O(log N) lookup *without* disturbing old readers we append, after all
//! records, a fixed-width **sorted pointer array**: one `u64` per sequence,
//! sorted by name, each holding the absolute file offset of that sequence's
//! existing index entry (the `nameSize` byte). The names themselves are **not**
//! duplicated — a probe follows a pointer into the original index, reads the
//! name there to compare, and on a hit reads the record offset that already sits
//! after the name. So the only added bytes are 8 per sequence plus a footer.
//!
//! ```text
//! nameIndex (only when the footer below is present)
//!   ptr          u64[seqCount]   sorted by name; ptr[i] = absolute offset of
//!                                the index entry for the i-th name in order
//! footer (24 bytes, at EOF)
//!   indexOffset  u64   absolute offset of the ptr array
//!   seqCount     u64   number of pointers
//!   flags        u32   bit 0 = records carry IUB tables
//!   magic        u32   INDEX_FOOTER_MAGIC ("INDX")
//! ```
//!
//! The footer replaces the 8-byte IUB trailer (its `flags` advertise whether IUB
//! tables are present), and like that trailer it lives past the last record, so
//! a legacy reader — which reaches records purely through index offsets — never
//! sees it.

use crate::error::{fmt_err, Result};
use crate::io::{peek_u32, peek_u64, put_u32, put_u64, put_u8, Reader};
use crate::seq::{
    base_to_twobit, is_acgt, is_iub_degenerate, twobit_to_base, Sequence,
};
use std::cell::RefCell;
use std::fs;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

const SIGNATURE: u32 = 0x1A41_2743;
/// Magic for the EOF trailer that flags an IUB-extended file.
/// As little-endian bytes this reads `2 B 2 U` ("2bit + 2-bit IUB").
const IUB_TRAILER_MAGIC: u32 = 0x5532_4232;
const IUB_EXT_VERSION: u32 = 1;

/// Magic for the EOF footer of a name-indexed file. As little-endian bytes this
/// reads `I N D X`. It must differ from `IUB_EXT_VERSION` (the last 4 bytes of
/// the plain IUB trailer) so the two trailers never alias.
const INDEX_FOOTER_MAGIC: u32 = 0x5844_4E49;
/// Footer flag: the records carry per-record IUB tables.
const INDEX_FLAG_IUB: u32 = 1;
/// Size of the name-index footer in bytes: indexOffset u64, count u64,
/// flags u32, magic u32.
const INDEX_FOOTER_SIZE: usize = 24;

/// A parsed twoBit file plus the flavour flags we detected.
#[derive(Debug, Clone)]
pub struct TwoBitFile {
    pub long: bool,
    pub iub: bool,
    /// True when the file carries the backward-compatible sorted-name index.
    pub indexed: bool,
    pub sequences: Vec<Sequence>,
}

// ---------------------------------------------------------------------------
// Block helpers
// ---------------------------------------------------------------------------

/// Collect maximal runs over `seq` where `key(b)` is `Some(k)`; a run breaks
/// when the key changes or becomes `None`. Returns `(start, size, key)`.
fn collect_runs<K, F>(seq: &[u8], key: F) -> Vec<(usize, usize, K)>
where
    K: PartialEq,
    F: Fn(u8) -> Option<K>,
{
    let mut runs = Vec::new();
    let mut cur: Option<(usize, K)> = None;
    for (i, &b) in seq.iter().enumerate() {
        match key(b) {
            Some(k) => match &cur {
                Some((_, ck)) if *ck == k => {}
                _ => {
                    if let Some((start, ck)) = cur.take() {
                        runs.push((start, i - start, ck));
                    }
                    cur = Some((i, k));
                }
            },
            None => {
                if let Some((start, ck)) = cur.take() {
                    runs.push((start, i - start, ck));
                }
            }
        }
    }
    if let Some((start, ck)) = cur.take() {
        runs.push((start, seq.len() - start, ck));
    }
    runs
}

/// N-blocks: every position whose base is not A/C/G/T (this includes bare N and
/// all degenerate codes — exactly what stock faToTwoBit records).
fn n_blocks(seq: &[u8]) -> Vec<(usize, usize)> {
    collect_runs(seq, |b| if is_acgt(b) { None } else { Some(()) })
        .into_iter()
        .map(|(s, l, _)| (s, l))
        .collect()
}

/// Soft-mask blocks: runs of lower-case letters.
fn mask_blocks(seq: &[u8]) -> Vec<(usize, usize)> {
    collect_runs(seq, |b| if b.is_ascii_lowercase() { Some(()) } else { None })
        .into_iter()
        .map(|(s, l, _)| (s, l))
        .collect()
}

/// IUB blocks: runs of identical specific degenerate codes (case-folded to
/// upper). Bare N and plain bases are excluded.
fn iub_blocks(seq: &[u8]) -> Vec<(usize, usize, u8)> {
    collect_runs(seq, |b| {
        if is_iub_degenerate(b) {
            Some(b.to_ascii_uppercase())
        } else {
            None
        }
    })
}

// ---------------------------------------------------------------------------
// Encoding
// ---------------------------------------------------------------------------

fn pack_dna(seq: &[u8]) -> Vec<u8> {
    let mut packed = vec![0u8; seq.len().div_ceil(4)];
    for (i, &b) in seq.iter().enumerate() {
        let code = base_to_twobit(b);
        let shift = 6 - 2 * (i % 4);
        packed[i / 4] |= code << shift;
    }
    packed
}

/// Serialise one sequence record (optionally with the IUB table appended).
fn encode_record(bases: &[u8], iub: bool) -> Vec<u8> {
    let mut out = Vec::new();
    let nb = n_blocks(bases);
    let mb = mask_blocks(bases);

    put_u32(&mut out, bases.len() as u32);

    put_u32(&mut out, nb.len() as u32);
    for &(s, _) in &nb {
        put_u32(&mut out, s as u32);
    }
    for &(_, l) in &nb {
        put_u32(&mut out, l as u32);
    }

    put_u32(&mut out, mb.len() as u32);
    for &(s, _) in &mb {
        put_u32(&mut out, s as u32);
    }
    for &(_, l) in &mb {
        put_u32(&mut out, l as u32);
    }

    put_u32(&mut out, 0); // reserved
    out.extend_from_slice(&pack_dna(bases));

    if iub {
        let ib = iub_blocks(bases);
        put_u32(&mut out, ib.len() as u32);
        for &(s, _, _) in &ib {
            put_u32(&mut out, s as u32);
        }
        for &(_, l, _) in &ib {
            put_u32(&mut out, l as u32);
        }
        for &(_, _, code) in &ib {
            put_u8(&mut out, code);
        }
    }
    out
}

/// Serialise a full twoBit file to bytes (always little-endian).
///
/// * `long`: write 64-bit index offsets (version 1).
/// * `iub`:  append per-record IUB tables and the discovery trailer.
pub fn to_bytes(seqs: &[Sequence], long: bool, iub: bool) -> Result<Vec<u8>> {
    build(seqs, long, iub, false)
}

/// Like [`to_bytes`] but also append the backward-compatible sorted-name index
/// (and an index footer in place of the plain IUB trailer). See the module docs.
pub fn to_bytes_indexed(seqs: &[Sequence], long: bool, iub: bool) -> Result<Vec<u8>> {
    build(seqs, long, iub, true)
}

fn build(seqs: &[Sequence], long: bool, iub: bool, index: bool) -> Result<Vec<u8>> {
    for s in seqs {
        if s.name.len() > 255 {
            return fmt_err(format!(
                "sequence name {:?} is {} bytes; twoBit names must be <= 255",
                s.name,
                s.name.len()
            ));
        }
        if s.bases.len() > u32::MAX as usize {
            return fmt_err(format!(
                "sequence {:?} has {} bases; a single twoBit sequence is limited to 2^32-1",
                s.name,
                s.bases.len()
            ));
        }
    }

    let records: Vec<Vec<u8>> = seqs.iter().map(|s| encode_record(&s.bases, iub)).collect();

    let off_size = if long { 8 } else { 4 };
    let header_size = 16usize;
    let index_size: usize = seqs.iter().map(|s| 1 + s.name.len() + off_size).sum();

    let trailer = if index {
        seqs.len() * 8 + INDEX_FOOTER_SIZE
    } else if iub {
        8
    } else {
        0
    };
    let mut out = Vec::with_capacity(
        header_size + index_size + records.iter().map(|r| r.len()).sum::<usize>() + trailer,
    );

    // Header.
    put_u32(&mut out, SIGNATURE);
    put_u32(&mut out, if long { 1 } else { 0 });
    put_u32(&mut out, seqs.len() as u32);
    put_u32(&mut out, 0);

    // Index — compute absolute record offsets, remembering where each entry
    // begins so the optional name index can point back at it.
    let mut entry_offsets = Vec::with_capacity(seqs.len());
    let mut offset = header_size + index_size;
    for (s, rec) in seqs.iter().zip(&records) {
        entry_offsets.push(out.len());
        put_u8(&mut out, s.name.len() as u8);
        out.extend_from_slice(s.name.as_bytes());
        if long {
            put_u64(&mut out, offset as u64);
        } else {
            put_u32(&mut out, offset as u32);
        }
        offset += rec.len();
    }

    // Records.
    for rec in &records {
        out.extend_from_slice(rec);
    }

    if index {
        // Sorted pointer array: one u64 per sequence, ordered by name, pointing
        // at the index entry above. Names live only in that entry.
        let index_offset = out.len();
        let mut order: Vec<usize> = (0..seqs.len()).collect();
        order.sort_by(|&a, &b| seqs[a].name.as_bytes().cmp(seqs[b].name.as_bytes()));
        for &i in &order {
            put_u64(&mut out, entry_offsets[i] as u64);
        }
        // Footer (discovery + IUB advertisement).
        put_u64(&mut out, index_offset as u64);
        put_u64(&mut out, seqs.len() as u64);
        put_u32(&mut out, if iub { INDEX_FLAG_IUB } else { 0 });
        put_u32(&mut out, INDEX_FOOTER_MAGIC);
    } else if iub {
        // Discovery trailer for the IUB extension.
        put_u32(&mut out, IUB_TRAILER_MAGIC);
        put_u32(&mut out, IUB_EXT_VERSION);
    }

    Ok(out)
}

/// Inspect the EOF trailer/footer. Returns `(iub, name_index)` where
/// `name_index` is `Some((ptr_array_offset, count))` for an indexed file.
fn detect_trailer(data: &[u8], little: bool) -> (bool, Option<(usize, usize)>) {
    let n = data.len();
    if n >= INDEX_FOOTER_SIZE
        && peek_u32(data, n - 4, little) == Some(INDEX_FOOTER_MAGIC)
    {
        let index_offset = peek_u64(data, n - INDEX_FOOTER_SIZE, little).unwrap_or(0) as usize;
        let count = peek_u64(data, n - 16, little).unwrap_or(0) as usize;
        let flags = peek_u32(data, n - 8, little).unwrap_or(0);
        return (flags & INDEX_FLAG_IUB != 0, Some((index_offset, count)));
    }
    let iub = n >= 8 && peek_u32(data, n - 8, little) == Some(IUB_TRAILER_MAGIC);
    (iub, None)
}

// ---------------------------------------------------------------------------
// Decoding
// ---------------------------------------------------------------------------

fn unpack_dna(packed: &[u8], dna_size: usize) -> Vec<u8> {
    let mut s = Vec::with_capacity(dna_size);
    for i in 0..dna_size {
        let byte = packed[i / 4];
        let shift = 6 - 2 * (i % 4);
        s.push(twobit_to_base((byte >> shift) & 0x03));
    }
    s
}

fn read_u32_array(r: &mut Reader, n: usize) -> Result<Vec<u32>> {
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        v.push(r.u32()?);
    }
    Ok(v)
}

fn decode_record(r: &mut Reader, iub: bool) -> Result<Vec<u8>> {
    let dna_size = r.u32()? as usize;

    let n_count = r.u32()? as usize;
    let n_starts = read_u32_array(r, n_count)?;
    let n_sizes = read_u32_array(r, n_count)?;

    let m_count = r.u32()? as usize;
    let m_starts = read_u32_array(r, m_count)?;
    let m_sizes = read_u32_array(r, m_count)?;

    let _reserved = r.u32()?;

    let packed = r.take(dna_size.div_ceil(4))?;
    let mut seq = unpack_dna(packed, dna_size);

    // Apply N-blocks.
    for (&s, &l) in n_starts.iter().zip(&n_sizes) {
        fill_range(&mut seq, s as usize, l as usize, b'N')?;
    }

    // Apply IUB overrides (extension only).
    if iub {
        let i_count = r.u32()? as usize;
        let i_starts = read_u32_array(r, i_count)?;
        let i_sizes = read_u32_array(r, i_count)?;
        let i_codes = r.take(i_count)?.to_vec();
        for ((&s, &l), &code) in i_starts.iter().zip(&i_sizes).zip(&i_codes) {
            fill_range(&mut seq, s as usize, l as usize, code)?;
        }
    }

    // Apply soft-mask last so it lower-cases whatever residue is present.
    for (&s, &l) in m_starts.iter().zip(&m_sizes) {
        let (start, end) = clamp_range(s as usize, l as usize, seq.len())?;
        for b in &mut seq[start..end] {
            b.make_ascii_lowercase();
        }
    }

    Ok(seq)
}

fn fill_range(seq: &mut [u8], start: usize, len: usize, val: u8) -> Result<()> {
    let (s, e) = clamp_range(start, len, seq.len())?;
    seq[s..e].fill(val);
    Ok(())
}

fn clamp_range(start: usize, len: usize, total: usize) -> Result<(usize, usize)> {
    let end = start
        .checked_add(len)
        .filter(|&e| e <= total)
        .ok_or_else(|| crate::error::Error::Format(format!(
            "block [{start}, {start}+{len}) is out of bounds for sequence of length {total}"
        )))?;
    Ok((start, end))
}

/// Detect endianness from the signature. Returns `little` or a format error.
fn detect_endianness(data: &[u8]) -> Result<bool> {
    match peek_u32(data, 0, true) {
        Some(SIGNATURE) => Ok(true),
        _ => match peek_u32(data, 0, false) {
            Some(SIGNATURE) => Ok(false),
            _ => fmt_err("not a twoBit file (bad signature)"),
        },
    }
}

/// Cheap signature check (either endianness) for format dispatch — only the
/// first 4 bytes are inspected, so a small file prefix suffices.
pub fn is_twobit(data: &[u8]) -> bool {
    peek_u32(data, 0, true) == Some(SIGNATURE) || peek_u32(data, 0, false) == Some(SIGNATURE)
}

/// Parse a full twoBit file from bytes (auto-detects endianness, long, and the
/// IUB extension).
pub fn from_bytes(data: &[u8]) -> Result<TwoBitFile> {
    if data.len() < 16 {
        return fmt_err("file too small to be a twoBit");
    }
    let little = detect_endianness(data)?;

    // Detect the IUB trailer / name-index footer at EOF.
    let (iub, name_index) = detect_trailer(data, little);
    let indexed = name_index.is_some();

    let mut r = Reader::new(data, little);
    let _sig = r.u32()?;
    let version = r.u32()?;
    if version > 1 {
        return fmt_err(format!(
            "unsupported twoBit version {version} (this tool handles 0 and 1)"
        ));
    }
    let long = version == 1;
    let count = r.u32()? as usize;
    let _reserved = r.u32()?;

    // Index.
    let mut index: Vec<(String, usize)> = Vec::with_capacity(count);
    for _ in 0..count {
        let name_len = r.u8()? as usize;
        let name = String::from_utf8_lossy(r.take(name_len)?).into_owned();
        let offset = if long { r.u64()? as usize } else { r.u32()? as usize };
        index.push((name, offset));
    }

    // Records.
    let mut sequences = Vec::with_capacity(count);
    for (name, offset) in index {
        r.seek(offset);
        let bases = decode_record(&mut r, iub)?;
        sequences.push(Sequence::new(name, bases));
    }

    Ok(TwoBitFile { long, iub, indexed, sequences })
}

// ---------------------------------------------------------------------------
// Random-access reader (decode only the requested region)
// ---------------------------------------------------------------------------

/// Backing store for the random-access reader: either an in-memory buffer or a
/// seekable file we `seek`+`read`. The file variant — which `open` uses — never
/// loads the whole file: every access is a small positioned read of just the
/// header, an index-probe path, or the requested window, matching how UCSC
/// `twoBitToFa` (`lseek`+`read` via its UDC layer) touches only what it needs.
enum Source {
    Mem(Vec<u8>),
    File { file: RefCell<fs::File>, len: usize },
}

impl Source {
    fn len(&self) -> usize {
        match self {
            Source::Mem(d) => d.len(),
            Source::File { len, .. } => *len,
        }
    }

    /// Read exactly `buf.len()` bytes starting at `off`.
    fn read_at(&self, off: usize, buf: &mut [u8]) -> Result<()> {
        match self {
            Source::Mem(d) => {
                let end = off
                    .checked_add(buf.len())
                    .filter(|&e| e <= d.len())
                    .ok_or_else(|| crate::error::Error::Format(format!(
                        "read of {} bytes at offset {off} is past the {}-byte buffer",
                        buf.len(),
                        d.len()
                    )))?;
                buf.copy_from_slice(&d[off..end]);
                Ok(())
            }
            Source::File { file, .. } => {
                let mut f = file.borrow_mut();
                f.seek(SeekFrom::Start(off as u64))?;
                f.read_exact(buf)?;
                Ok(())
            }
        }
    }

    fn bytes_at(&self, off: usize, n: usize) -> Result<Vec<u8>> {
        let mut v = vec![0u8; n];
        self.read_at(off, &mut v)?;
        Ok(v)
    }

    fn u8_at(&self, off: usize) -> Result<u8> {
        let mut b = [0u8; 1];
        self.read_at(off, &mut b)?;
        Ok(b[0])
    }

    fn u32_at(&self, off: usize, little: bool) -> Result<u32> {
        let mut b = [0u8; 4];
        self.read_at(off, &mut b)?;
        Ok(if little { u32::from_le_bytes(b) } else { u32::from_be_bytes(b) })
    }

    fn u64_at(&self, off: usize, little: bool) -> Result<u64> {
        let mut b = [0u8; 8];
        self.read_at(off, &mut b)?;
        Ok(if little { u64::from_le_bytes(b) } else { u64::from_be_bytes(b) })
    }
}

/// The flat TOC parsed into a name→offset map. Built lazily, because a
/// name-indexed file can answer single lookups by binary search without ever
/// paying this O(N) scan — the whole point of the index.
struct FlatToc {
    by_name: std::collections::HashMap<String, usize>,
    order: Vec<String>,
}

/// Parses only the header up front and then `seek`s to decode an arbitrary
/// sub-range without touching the rest of the file — the access pattern UCSC
/// `twoBitToFa` uses. Created from a path via [`TwoBitReader::open`] (seek-based,
/// never slurps the file) or from an in-memory buffer via
/// [`TwoBitReader::from_vec`]. When the file carries a name index, lookups go
/// through it (O(log N)); otherwise the flat TOC is parsed on first use.
pub struct TwoBitReader {
    src: Source,
    little: bool,
    pub long: bool,
    pub iub: bool,
    count: usize,
    /// `(ptr_array_offset, count)` when the file carries the sorted-name index.
    name_index: Option<(usize, usize)>,
    flat: std::cell::OnceCell<FlatToc>,
}

impl TwoBitReader {
    /// Open by path using `seek`+`read`: only the header, an index probe, and
    /// the requested window are read — never the whole file.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let file = fs::File::open(path)?;
        let len = file.metadata()?.len() as usize;
        Self::from_source(Source::File { file: RefCell::new(file), len })
    }

    pub fn from_vec(data: Vec<u8>) -> Result<Self> {
        Self::from_source(Source::Mem(data))
    }

    fn from_source(src: Source) -> Result<Self> {
        if src.len() < 16 {
            return fmt_err("file too small to be a twoBit");
        }
        // Endianness: read the 4 signature bytes once and try both orders.
        let mut sig = [0u8; 4];
        src.read_at(0, &mut sig)?;
        let little = if u32::from_le_bytes(sig) == SIGNATURE {
            true
        } else if u32::from_be_bytes(sig) == SIGNATURE {
            false
        } else {
            return fmt_err("not a twoBit file (bad signature)");
        };

        let (iub, name_index) = detect_trailer_src(&src, little)?;

        let version = src.u32_at(4, little)?;
        if version > 1 {
            return fmt_err(format!("unsupported twoBit version {version}"));
        }
        let long = version == 1;
        let count = src.u32_at(8, little)? as usize;

        Ok(TwoBitReader {
            src,
            little,
            long,
            iub,
            count,
            name_index,
            flat: std::cell::OnceCell::new(),
        })
    }

    /// Whether this file carries the backward-compatible name index.
    pub fn has_name_index(&self) -> bool {
        self.name_index.is_some()
    }

    /// Parse the flat TOC into a name→offset map on first use.
    fn flat(&self) -> Result<&FlatToc> {
        if let Some(f) = self.flat.get() {
            return Ok(f);
        }
        let mut by_name = std::collections::HashMap::with_capacity(self.count);
        let mut order = Vec::with_capacity(self.count);
        if self.count > 0 {
            // The index runs from offset 16 to the first record. Records are
            // written in index order with ascending offsets, so entry 0's offset
            // marks the index end — letting us pull the whole index in one read
            // instead of one syscall per entry.
            let first_name_len = self.src.u8_at(16)? as usize;
            let first_off_pos = 16 + 1 + first_name_len;
            let index_end = if self.long {
                self.src.u64_at(first_off_pos, self.little)? as usize
            } else {
                self.src.u32_at(first_off_pos, self.little)? as usize
            };
            let buf = self.src.bytes_at(16, index_end.saturating_sub(16))?;
            let mut r = Reader::new(&buf, self.little);
            for _ in 0..self.count {
                let name_len = r.u8()? as usize;
                let name = String::from_utf8_lossy(r.take(name_len)?).into_owned();
                let offset = if self.long { r.u64()? as usize } else { r.u32()? as usize };
                by_name.insert(name.clone(), offset);
                order.push(name);
            }
        }
        // OnceCell::set fails only if already set, which the early return rules
        // out; ignore the returned Err in that impossible case.
        let _ = self.flat.set(FlatToc { by_name, order });
        Ok(self.flat.get().unwrap())
    }

    /// Resolve a sequence name to its record offset. Uses the sorted pointer
    /// array (O(log N)) when present, else the flat TOC map.
    fn record_offset(&self, name: &str) -> Result<usize> {
        if let Some((idx_off, n)) = self.name_index {
            return self
                .lookup_indexed(idx_off, n, name)?
                .ok_or_else(|| crate::error::Error::Format(format!("no sequence named {name:?}")));
        }
        self.flat()?
            .by_name
            .get(name)
            .copied()
            .ok_or_else(|| crate::error::Error::Format(format!("no sequence named {name:?}")))
    }

    /// Binary-search the sorted pointer array. Each probe dereferences a pointer
    /// into the original index entry to read the name (for comparison) and, on a
    /// hit, the record offset that follows it — a handful of small reads total.
    fn lookup_indexed(&self, idx_off: usize, n: usize, name: &str) -> Result<Option<usize>> {
        let little = self.little;
        let key = name.as_bytes();
        let (mut lo, mut hi) = (0usize, n);
        while lo < hi {
            let mid = (lo + hi) / 2;
            let entry = self.src.u64_at(idx_off + mid * 8, little)? as usize;
            let nlen = self.src.u8_at(entry)? as usize;
            let entry_name = self.src.bytes_at(entry + 1, nlen)?;
            match entry_name.as_slice().cmp(key) {
                std::cmp::Ordering::Less => lo = mid + 1,
                std::cmp::Ordering::Greater => hi = mid,
                std::cmp::Ordering::Equal => {
                    let off_pos = entry + 1 + nlen;
                    let rec = if self.long {
                        self.src.u64_at(off_pos, little)? as usize
                    } else {
                        self.src.u32_at(off_pos, little)? as usize
                    };
                    return Ok(Some(rec));
                }
            }
        }
        Ok(None)
    }

    pub fn names(&self) -> &[String] {
        // If parsing the flat TOC ever fails the file is corrupt; treat it as no
        // names rather than panicking from this infallible accessor.
        match self.flat() {
            Ok(f) => &f.order,
            Err(_) => &[],
        }
    }

    /// `(name, dnaSize)` for every sequence, in file order.
    pub fn sequence_infos(&self) -> Result<Vec<(String, usize)>> {
        let flat = self.flat()?;
        let mut v = Vec::with_capacity(flat.order.len());
        for name in &flat.order {
            let off = flat.by_name[name];
            v.push((name.clone(), self.src.u32_at(off, self.little)? as usize));
        }
        Ok(v)
    }

    /// Decode `[start, end)` of `name` (end clamped to the sequence length;
    /// `end == None` means to the end). Only the requested window is decoded.
    pub fn extract(&self, name: &str, start: usize, end: Option<usize>) -> Result<Vec<u8>> {
        let off = self.record_offset(name)?;
        let little = self.little;

        let dna_size = self.src.u32_at(off, little)? as usize;
        let end = end.unwrap_or(dna_size).min(dna_size);
        let start = start.min(end);

        // Walk the record header by *offset arithmetic* — never materialise the
        // (potentially millions of) block entries. Each table is a pair of
        // parallel, ascending, disjoint u32 arrays we binary-search in place.
        let n_count = self.src.u32_at(off + 4, little)? as usize;
        let n_starts_off = off + 8;
        let n_sizes_off = n_starts_off + n_count * 4;
        let m_count_off = n_sizes_off + n_count * 4;

        let m_count = self.src.u32_at(m_count_off, little)? as usize;
        let m_starts_off = m_count_off + 4;
        let m_sizes_off = m_starts_off + m_count * 4;
        let packed_start = m_sizes_off + m_count * 4 + 4; // +4 = reserved

        // Decode just the window: read only the packed bytes that span it.
        let mut seq = Vec::with_capacity(end - start);
        if start < end {
            let first_byte = start / 4;
            let last_byte = (end - 1) / 4;
            let packed = self
                .src
                .bytes_at(packed_start + first_byte, last_byte - first_byte + 1)?;
            for i in start..end {
                let byte = packed[i / 4 - first_byte];
                let shift = 6 - 2 * (i % 4);
                seq.push(twobit_to_base((byte >> shift) & 0x03));
            }
        }

        // N-blocks → 'N'.
        for_overlapping_blocks(
            &self.src, little, n_starts_off, n_sizes_off, n_count, start, end,
            |_, lo, hi| seq[lo..hi].fill(b'N'),
        );

        // IUB overrides → the exact code.
        if self.iub {
            let packed_len = dna_size.div_ceil(4);
            let i_count_off = packed_start + packed_len;
            let i_count = self.src.u32_at(i_count_off, little)? as usize;
            let i_starts_off = i_count_off + 4;
            let i_sizes_off = i_starts_off + i_count * 4;
            let i_codes_off = i_sizes_off + i_count * 4;
            for_overlapping_blocks(
                &self.src, little, i_starts_off, i_sizes_off, i_count, start, end,
                |k, lo, hi| {
                    let code = self.src.u8_at(i_codes_off + k).unwrap_or(b'N');
                    seq[lo..hi].fill(code);
                },
            );
        }

        // Soft-mask → lower-case whatever is present.
        for_overlapping_blocks(
            &self.src, little, m_starts_off, m_sizes_off, m_count, start, end,
            |_, lo, hi| seq[lo..hi].make_ascii_lowercase(),
        );
        Ok(seq)
    }
}

/// Visit the blocks of a sorted, disjoint block table (parallel u32 `starts`/
/// `sizes` arrays in `src`) that overlap the window `[start, end)`, calling
/// `f(block_index, rel_lo, rel_hi)` with window-relative coordinates. Uses
/// binary search + a short forward scan, so cost is O(log n + overlaps) — for a
/// file source that's a handful of 4-byte reads, not a full table scan.
fn for_overlapping_blocks(
    src: &Source,
    little: bool,
    starts_off: usize,
    sizes_off: usize,
    count: usize,
    start: usize,
    end: usize,
    mut f: impl FnMut(usize, usize, usize),
) {
    if count == 0 || start >= end {
        return;
    }
    let bstart = |i: usize| src.u32_at(starts_off + i * 4, little).unwrap_or(0) as usize;
    let bsize = |i: usize| src.u32_at(sizes_off + i * 4, little).unwrap_or(0) as usize;

    // First index whose block start is >= `start`.
    let mut lo = 0usize;
    let mut hi = count;
    while lo < hi {
        let mid = (lo + hi) / 2;
        if bstart(mid) < start {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    // Step back one so a block that begins before `start` but extends into the
    // window is not missed (blocks are disjoint, so only one can).
    let mut i = lo.saturating_sub(1);
    while i < count {
        let bs = bstart(i);
        if bs >= end {
            break;
        }
        let block_lo = bs.max(start);
        let block_hi = (bs + bsize(i)).min(end);
        if block_lo < block_hi {
            f(i, block_lo - start, block_hi - start);
        }
        i += 1;
    }
}

/// Source-backed twin of [`detect_trailer`] for the seek reader: reads only the
/// last few bytes rather than requiring the whole file in memory.
fn detect_trailer_src(src: &Source, little: bool) -> Result<(bool, Option<(usize, usize)>)> {
    let len = src.len();
    if len >= INDEX_FOOTER_SIZE
        && src.u32_at(len - 4, little)? == INDEX_FOOTER_MAGIC
    {
        let index_offset = src.u64_at(len - INDEX_FOOTER_SIZE, little)? as usize;
        let count = src.u64_at(len - 16, little)? as usize;
        let flags = src.u32_at(len - 8, little)?;
        return Ok((flags & INDEX_FLAG_IUB != 0, Some((index_offset, count))));
    }
    if len >= 8 && src.u32_at(len - 8, little)? == IUB_TRAILER_MAGIC {
        return Ok((true, None));
    }
    Ok((false, None))
}

// ---------------------------------------------------------------------------
// File-path convenience wrappers
// ---------------------------------------------------------------------------

pub fn write_file(path: impl AsRef<Path>, seqs: &[Sequence], long: bool, iub: bool) -> Result<()> {
    fs::write(path, to_bytes(seqs, long, iub)?)?;
    Ok(())
}

/// Like [`write_file`] but append the backward-compatible sorted-name index.
pub fn write_file_indexed(
    path: impl AsRef<Path>,
    seqs: &[Sequence],
    long: bool,
    iub: bool,
) -> Result<()> {
    fs::write(path, to_bytes_indexed(seqs, long, iub)?)?;
    Ok(())
}

pub fn read_file(path: impl AsRef<Path>) -> Result<TwoBitFile> {
    from_bytes(&fs::read(path)?)
}

/// Per-sequence statistics for `info`.
pub struct SeqStats {
    pub name: String,
    pub len: usize,
    pub n_blocks: usize,
    pub mask_blocks: usize,
    pub iub_blocks: usize,
}

pub fn stats(file: &TwoBitFile) -> Vec<SeqStats> {
    file.sequences
        .iter()
        .map(|s| SeqStats {
            name: s.name.clone(),
            len: s.bases.len(),
            n_blocks: n_blocks(&s.bases).len(),
            mask_blocks: mask_blocks(&s.bases).len(),
            iub_blocks: iub_blocks(&s.bases).len(),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seq(name: &str, bases: &str) -> Sequence {
        Sequence::new(name, bases.as_bytes().to_vec())
    }

    fn roundtrip(seqs: &[Sequence], long: bool, iub: bool) -> Vec<Sequence> {
        let bytes = to_bytes(seqs, long, iub).unwrap();
        from_bytes(&bytes).unwrap().sequences
    }

    #[test]
    fn standard_roundtrip_basic() {
        let s = vec![seq("chr1", "ACGTACGTACGT"), seq("chr2", "TTTTGGGGCCCCAAAA")];
        assert_eq!(roundtrip(&s, false, false), s);
    }

    #[test]
    fn n_and_mask_blocks() {
        // lower-case = soft mask, N runs collapse the degenerate-free way.
        let s = vec![seq("c", "ACGTnnnNNNacgtACGT")];
        let out = roundtrip(&s, false, false);
        // Without the extension, degenerate codes would already be N; here input
        // only has plain n/N so it round-trips exactly (case preserved).
        assert_eq!(out, s);
    }

    #[test]
    fn long_format_roundtrip() {
        let s = vec![seq("chrLong", "ACGTACGTNNNNacgtRYSW")];
        // In long mode without iub, the R/Y/S/W degrade to N (standard behaviour).
        let out = roundtrip(&s, true, false);
        assert!(from_bytes(&to_bytes(&s, true, false).unwrap()).unwrap().long);
        assert_eq!(out[0].bases, b"ACGTACGTNNNNacgtNNNN");
    }

    #[test]
    fn degenerate_lost_without_extension() {
        let s = vec![seq("c", "ACGTRYSWKMBDHVNacgt")];
        let out = roundtrip(&s, false, false);
        // Every degenerate code becomes N; soft-mask + plain bases survive.
        assert_eq!(out[0].bases, b"ACGTNNNNNNNNNNNacgt");
    }

    #[test]
    fn degenerate_preserved_with_extension() {
        let s = vec![
            seq("c", "ACGTRYSWKMBDHVNacgtryswn"),
            seq("d", "NNNRRRYYYacgt"),
        ];
        let out = roundtrip(&s, false, true);
        assert_eq!(out, s, "IUB extension must preserve degenerate codes & case");
    }

    #[test]
    fn extension_is_backward_compatible() {
        // A reader that doesn't know about the extension (iub=false on decode)
        // must still parse the file as a plain 2bit, seeing N for degenerates.
        let s = vec![seq("c", "ACGTRYSWacgtN"), seq("d", "RRRRGGGG")];
        let bytes = to_bytes(&s, false, true).unwrap();

        // Simulate a legacy reader: same bytes, but force iub=false decoding by
        // stripping the trailer so detection fails.
        let legacy = &bytes[..bytes.len() - 8];
        let parsed = from_bytes(legacy).unwrap();
        assert!(!parsed.iub);
        assert_eq!(parsed.sequences[0].bases, b"ACGTNNNNacgtN");
        assert_eq!(parsed.sequences[1].bases, b"NNNNGGGG");

        // And our own reader recovers everything.
        let full = from_bytes(&bytes).unwrap();
        assert!(full.iub);
        assert_eq!(full.sequences, s);
    }

    #[test]
    fn random_access_matches_full_decode() {
        // A sequence with mixed plain bases, soft-mask, scattered N and IUB runs.
        let s = vec![
            seq("c", "ACGTNNNryacgtRRRYYYacgtNacgtKMBacgtACGTnnnNNN"),
            seq("d", "RRRRRRRRRRacgtACGTNNNN"),
        ];
        for &iub in &[false, true] {
            let bytes = to_bytes(&s, false, iub).unwrap();
            let full = from_bytes(&bytes).unwrap().sequences;
            let rd = TwoBitReader::from_vec(bytes).unwrap();
            for (whole, src) in full.iter().zip(&s) {
                let n = src.bases.len();
                // Probe every sub-window [a, b) and compare to the full decode.
                for a in 0..=n {
                    for b in a..=n {
                        let got = rd.extract(&src.name, a, Some(b)).unwrap();
                        assert_eq!(
                            got,
                            &whole.bases[a..b],
                            "iub={iub} {}:{a}-{b}",
                            src.name
                        );
                    }
                }
                // end=None means to the end.
                assert_eq!(rd.extract(&src.name, 3, None).unwrap(), &whole.bases[3..]);
            }
        }
    }

    #[test]
    fn empty_and_edge_sizes() {
        let s = vec![seq("empty", ""), seq("one", "A"), seq("five", "ACGTN")];
        assert_eq!(roundtrip(&s, false, true), s);
        assert_eq!(roundtrip(&s, true, true), s);
    }

    #[test]
    fn indexed_roundtrip_full_decode() {
        // Names deliberately out of sorted order, mixed lengths.
        let s = vec![
            seq("chrM", "ACGTRYSWacgtNNN"),
            seq("chr1", "RRRRGGGGacgtKMB"),
            seq("chr10", "ACGTNNNryacgtRR"),
            seq("chr2", "TTTT"),
        ];
        for &long in &[false, true] {
            let bytes = to_bytes_indexed(&s, long, true).unwrap();
            let f = from_bytes(&bytes).unwrap();
            assert!(f.indexed && f.iub && f.long == long);
            assert_eq!(f.sequences, s);
        }
    }

    #[test]
    fn indexed_is_backward_compatible() {
        // Strip the index footer + pointer array: the prefix must still parse as
        // an ordinary 2bit (degenerates show as N), proving old readers cope.
        let s = vec![seq("b", "ACGTRYSW"), seq("a", "RRRRGGGG")];
        let bytes = to_bytes_indexed(&s, false, true).unwrap();
        let footer_off = {
            let little = detect_endianness(&bytes).unwrap();
            peek_u64(&bytes, bytes.len() - INDEX_FOOTER_SIZE, little).unwrap() as usize
        };
        let legacy = &bytes[..footer_off]; // drop ptr array + footer
        let parsed = from_bytes(legacy).unwrap();
        assert!(!parsed.iub && !parsed.indexed);
        assert_eq!(parsed.sequences[0].bases, b"ACGTNNNN");
        assert_eq!(parsed.sequences[1].bases, b"NNNNGGGG");
    }

    #[test]
    fn indexed_reader_uses_binary_search() {
        let s = vec![
            seq("seqZ", "ACGTNNNryacgtRRRYYYacgtNacgtKMB"),
            seq("seqA", "RRRRRRRRRRacgtACGTNNNN"),
            seq("seqM", "ACGTacgtACGT"),
        ];
        let bytes = to_bytes_indexed(&s, false, true).unwrap();
        let full = from_bytes(&bytes).unwrap().sequences;
        let rd = TwoBitReader::from_vec(bytes).unwrap();
        assert!(rd.has_name_index());

        for (whole, src) in full.iter().zip(&s) {
            let n = src.bases.len();
            for a in 0..=n {
                for b in a..=n {
                    assert_eq!(
                        rd.extract(&src.name, a, Some(b)).unwrap(),
                        &whole.bases[a..b],
                        "{}:{a}-{b}",
                        src.name
                    );
                }
            }
        }
        // Missing name is an error, not a panic or wrong hit.
        assert!(rd.extract("nope", 0, None).is_err());
        // Lazy flat TOC still answers names()/sequence_infos() on demand.
        let mut names = rd.names().to_vec();
        names.sort();
        assert_eq!(names, vec!["seqA", "seqM", "seqZ"]);
        assert_eq!(rd.sequence_infos().unwrap().len(), 3);
    }

    #[test]
    fn file_backed_reader_matches_mem() {
        // The seek+read (File) source must decode byte-for-byte like the in-memory
        // source across every window, for plain / long / iub / indexed flavours.
        let s = vec![
            seq("seqZ", "ACGTNNNryacgtRRRYYYacgtNacgtKMBacgtACGTnnnNNN"),
            seq("seqA", "RRRRRRRRRRacgtACGTNNNN"),
            seq("seqM", "ACGTacgtACGTnnnNNNacgtryswKMB"),
        ];
        let combos = [
            (false, false, false),
            (false, true, false),
            (true, true, true),
            (false, true, true),
            (false, false, true),
        ];
        for (ci, &(long, iub, index)) in combos.iter().enumerate() {
            let bytes = build(&s, long, iub, index).unwrap();
            let mem = TwoBitReader::from_vec(bytes.clone()).unwrap();

            let path = std::env::temp_dir().join(format!(
                "seqformat_filetest_{}_{ci}.2bit",
                std::process::id()
            ));
            std::fs::write(&path, &bytes).unwrap();
            let file = TwoBitReader::open(&path).unwrap();

            assert_eq!(file.iub, mem.iub);
            assert_eq!(file.long, mem.long);
            assert_eq!(file.has_name_index(), index);

            for src in &s {
                let n = src.bases.len();
                for a in 0..=n {
                    for b in a..=n {
                        assert_eq!(
                            file.extract(&src.name, a, Some(b)).unwrap(),
                            mem.extract(&src.name, a, Some(b)).unwrap(),
                            "combo {ci} {}:{a}-{b}",
                            src.name
                        );
                    }
                }
            }
            assert!(file.extract("missing", 0, None).is_err());
            let mut names = file.names().to_vec();
            names.sort();
            assert_eq!(names, vec!["seqA", "seqM", "seqZ"]);
            assert_eq!(file.sequence_infos().unwrap().len(), 3);

            std::fs::remove_file(&path).ok();
        }
    }
}
