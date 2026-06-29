//! The "samtools format": FASTA indexed by `.fai`, optionally BGZF-compressed
//! (`.fa.gz`) with a `.gzi` block index. This is an *index-over-FASTA* approach,
//! not a self-contained container — the bases stay as FASTA text; the sidecar
//! files just enable random access (and bgzip adds compression).
//!
//! What this module implements, all interoperable with `samtools faidx`:
//!   * write plain FASTA + `.fai`
//!   * write BGZF FASTA + `.fai` + `.gzi`
//!   * random-access extraction from either (reads `.fai`; for BGZF it scans the
//!     self-describing block headers, so a `.gzi` is not required to read)
//!
//! `.fai` line: `name<TAB>length<TAB>offset<TAB>linebases<TAB>linewidth` where
//! `offset` is the byte offset of the first base in the *uncompressed* FASTA,
//! `linebases` is bases per line and `linewidth` is bytes per line (incl. `\n`).
//!
//! BGZF is gzip with a `BC` extra field carrying each block's size; see
//! <https://samtools.github.io/hts-specs/SAMv1.pdf> §4.

use crate::error::{fmt_err, Error, Result};
use crate::seq::Sequence;
use libdeflater::{CompressionLvl, Compressor, Crc, Decompressor};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

/// Max uncompressed bytes per BGZF block (htslib default, guarantees the
/// compressed block stays under the 64 KiB BSIZE limit).
const BGZF_BLOCK: usize = 0xff00;

/// The standard 28-byte empty BGZF block that marks end-of-file.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

// ---------------------------------------------------------------------------
// .fai records
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct FaiRec {
    pub name: String,
    pub length: usize,
    pub offset: usize,
    pub linebases: usize,
    pub linewidth: usize,
}

/// Build FASTA bytes together with the matching `.fai` records, so the recorded
/// offsets are exact by construction. `width` 0 means one line per sequence.
fn build_fasta_with_fai(seqs: &[Sequence], width: usize) -> (Vec<u8>, Vec<FaiRec>) {
    let mut fasta = Vec::new();
    let mut recs = Vec::with_capacity(seqs.len());
    for s in seqs {
        fasta.push(b'>');
        fasta.extend_from_slice(s.name.as_bytes());
        fasta.push(b'\n');
        let offset = fasta.len();

        let linebases = if width == 0 || s.bases.len() <= width {
            s.bases.len()
        } else {
            width
        };
        if width == 0 || s.bases.is_empty() {
            fasta.extend_from_slice(&s.bases);
            fasta.push(b'\n');
        } else {
            for chunk in s.bases.chunks(width) {
                fasta.extend_from_slice(chunk);
                fasta.push(b'\n');
            }
        }
        recs.push(FaiRec {
            name: s.name.clone(),
            length: s.bases.len(),
            offset,
            linebases,
            linewidth: linebases + 1,
        });
    }
    (fasta, recs)
}

fn fai_bytes(recs: &[FaiRec]) -> Vec<u8> {
    let mut s = String::new();
    for r in recs {
        s.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\n",
            r.name, r.length, r.offset, r.linebases, r.linewidth
        ));
    }
    s.into_bytes()
}

fn parse_fai(text: &str) -> Result<Vec<FaiRec>> {
    let mut recs = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 {
            return fmt_err(format!(".fai line {} has {} fields, need 5", i + 1, f.len()));
        }
        let p = |idx: usize| -> Result<usize> {
            f[idx]
                .parse()
                .map_err(|_| Error::Format(format!("bad integer in .fai: {:?}", f[idx])))
        };
        recs.push(FaiRec {
            name: f[0].to_string(),
            length: p(1)?,
            offset: p(2)?,
            linebases: p(3)?,
            linewidth: p(4)?,
        });
    }
    Ok(recs)
}

// ---------------------------------------------------------------------------
// BGZF write
// ---------------------------------------------------------------------------

fn deflate_raw(data: &[u8], level: u32) -> Result<Vec<u8>> {
    let lvl = CompressionLvl::new(level as i32)
        .map_err(|e| Error::Format(format!("invalid deflate level {level}: {e:?}")))?;
    let mut c = Compressor::new(lvl);
    let mut out = vec![0u8; c.deflate_compress_bound(data.len())];
    let n = c
        .deflate_compress(data, &mut out)
        .map_err(|e| Error::Format(format!("deflate failed: {e:?}")))?;
    out.truncate(n);
    Ok(out)
}

fn inflate_raw(data: &[u8], expected: usize) -> Result<Vec<u8>> {
    if expected == 0 {
        return Ok(Vec::new());
    }
    let mut d = Decompressor::new();
    let mut out = vec![0u8; expected];
    let n = d
        .deflate_decompress(data, &mut out)
        .map_err(|e| Error::Format(format!("inflate failed: {e:?}")))?;
    out.truncate(n);
    Ok(out)
}

fn bgzf_block(uncomp: &[u8], level: u32) -> Result<Vec<u8>> {
    let cdata = deflate_raw(uncomp, level)?;
    let bsize = 18 + cdata.len() + 8 - 1; // total block length minus 1
    if bsize > 0xffff {
        return fmt_err("BGZF block overflow (compressed chunk too large)");
    }
    let mut b = Vec::with_capacity(bsize + 1);
    b.extend_from_slice(&[
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    ]);
    b.extend_from_slice(&(bsize as u16).to_le_bytes());
    b.extend_from_slice(&cdata);
    let mut crc = Crc::new();
    crc.update(uncomp);
    b.extend_from_slice(&crc.sum().to_le_bytes());
    b.extend_from_slice(&(uncomp.len() as u32).to_le_bytes());
    Ok(b)
}

/// BGZF-compress `data`, returning the compressed bytes and the `(compressed,
/// uncompressed)` start offset of each *data* block (one per block, the first
/// being `(0,0)`) — exactly what the `.gzi` is derived from. The EOF marker is
/// appended but not recorded as a boundary (matching htslib).
fn bgzf_compress(data: &[u8], level: u32) -> Result<(Vec<u8>, Vec<(u64, u64)>)> {
    let mut out = Vec::new();
    let mut boundaries = Vec::new();
    let mut uoff: u64 = 0;
    for chunk in data.chunks(BGZF_BLOCK) {
        boundaries.push((out.len() as u64, uoff));
        out.extend_from_slice(&bgzf_block(chunk, level)?);
        uoff += chunk.len() as u64;
    }
    out.extend_from_slice(&BGZF_EOF);
    Ok((out, boundaries))
}

/// `.gzi`: u64 count, then for every data-block boundary *after* the first a
/// `(compressed_offset, uncompressed_offset)` u64 pair (htslib layout). A
/// single-block file therefore has count 0.
fn gzi_bytes(boundaries: &[(u64, u64)]) -> Vec<u8> {
    let entries = boundaries.get(1..).unwrap_or(&[]);
    let mut v = Vec::with_capacity(8 + entries.len() * 16);
    v.extend_from_slice(&(entries.len() as u64).to_le_bytes());
    for &(c, u) in entries {
        v.extend_from_slice(&c.to_le_bytes());
        v.extend_from_slice(&u.to_le_bytes());
    }
    v
}

fn sidecar(path: &Path, ext: &str) -> PathBuf {
    let mut s = path.as_os_str().to_os_string();
    s.push(ext);
    PathBuf::from(s)
}

/// Write `seqs` as a samtools-style indexed FASTA. With `bgzip`, the output is
/// BGZF-compressed and a `.gzi` is written alongside the `.fai`.
pub fn write_file(
    path: impl AsRef<Path>,
    seqs: &[Sequence],
    width: usize,
    bgzip: bool,
    level: u32,
) -> Result<()> {
    let path = path.as_ref();
    let (fasta, recs) = build_fasta_with_fai(seqs, width);
    fs::write(sidecar(path, ".fai"), fai_bytes(&recs))?;
    if bgzip {
        let (comp, boundaries) = bgzf_compress(&fasta, level)?;
        fs::write(path, comp)?;
        fs::write(sidecar(path, ".gzi"), gzi_bytes(&boundaries))?;
    } else {
        fs::write(path, fasta)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// BGZF read + random access
// ---------------------------------------------------------------------------

struct Block {
    c_start: usize,
    c_len: usize,
    u_start: usize,
    u_len: usize,
}

pub fn is_bgzf(d: &[u8]) -> bool {
    d.len() >= 18 && d[0] == 0x1f && d[1] == 0x8b && d[2] == 0x08 && d[3] == 0x04 && d[12] == b'B'
        && d[13] == b'C'
}

pub fn is_fasta(d: &[u8]) -> bool {
    d.first() == Some(&b'>')
}

/// Random-access reader over a samtools-style indexed FASTA (plain or BGZF).
pub struct FaidxReader {
    data: Vec<u8>,
    bgzf_blocks: Option<Vec<Block>>, // None => plain FASTA
    by_name: HashMap<String, FaiRec>,
    order: Vec<String>,
}

impl FaidxReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let data = fs::read(path)?;
        let bgzf_blocks = if is_bgzf(&data) {
            Some(scan_bgzf(&data)?)
        } else if is_fasta(&data) {
            None
        } else {
            return fmt_err("not a FASTA or BGZF file");
        };

        // Prefer the sidecar .fai; otherwise build one by scanning.
        let fai_path = sidecar(path, ".fai");
        let recs = if fai_path.exists() {
            parse_fai(&fs::read_to_string(fai_path)?)?
        } else {
            build_fai_by_scan(&data, &bgzf_blocks)?
        };

        let mut by_name = HashMap::with_capacity(recs.len());
        let mut order = Vec::with_capacity(recs.len());
        for r in recs {
            order.push(r.name.clone());
            by_name.insert(r.name.clone(), r);
        }
        Ok(FaidxReader { data, bgzf_blocks, by_name, order })
    }

    pub fn names(&self) -> &[String] {
        &self.order
    }

    pub fn sequence_infos(&self) -> Vec<(String, usize)> {
        self.order
            .iter()
            .map(|n| (n.clone(), self.by_name[n].length))
            .collect()
    }

    /// Read uncompressed FASTA bytes `[a, b)` (transparently across BGZF blocks).
    fn read_uncompressed(&self, a: usize, b: usize) -> Result<Vec<u8>> {
        if a >= b {
            return Ok(Vec::new());
        }
        match &self.bgzf_blocks {
            None => Ok(self.data[a..b.min(self.data.len())].to_vec()),
            Some(blocks) => {
                let mut out = Vec::with_capacity(b - a);
                // First block whose end is past `a`.
                let mut i = blocks.partition_point(|blk| blk.u_start + blk.u_len <= a);
                while i < blocks.len() && blocks[i].u_start < b {
                    let blk = &blocks[i];
                    let cdata = &self.data[blk.c_start + 18..blk.c_start + blk.c_len - 8];
                    let u = inflate_raw(cdata, blk.u_len)?;
                    let lo = a.max(blk.u_start) - blk.u_start;
                    let hi = b.min(blk.u_start + blk.u_len) - blk.u_start;
                    out.extend_from_slice(&u[lo..hi]);
                    i += 1;
                }
                Ok(out)
            }
        }
    }

    /// Extract `[start, end)` of `name` (0-based, half-open).
    pub fn extract(&self, name: &str, start: usize, end: Option<usize>) -> Result<Vec<u8>> {
        let rec = self
            .by_name
            .get(name)
            .ok_or_else(|| Error::Format(format!("no sequence named {name:?}")))?;
        let end = end.unwrap_or(rec.length).min(rec.length);
        let start = start.min(end);
        if start == end {
            return Ok(Vec::new());
        }
        let lb = rec.linebases.max(1);
        let byte_off = |p: usize| rec.offset + (p / lb) * rec.linewidth + (p % lb);
        let raw = self.read_uncompressed(byte_off(start), byte_off(end))?;
        // Strip line breaks to recover the bases.
        Ok(raw
            .into_iter()
            .filter(|&c| c != b'\n' && c != b'\r')
            .collect())
    }
}

fn scan_bgzf(data: &[u8]) -> Result<Vec<Block>> {
    let mut blocks = Vec::new();
    let mut c = 0usize;
    let mut u = 0usize;
    while c + 18 <= data.len() {
        if !(data[c] == 0x1f && data[c + 1] == 0x8b) {
            break;
        }
        let bsize = u16::from_le_bytes([data[c + 16], data[c + 17]]) as usize;
        let blen = bsize + 1;
        if c + blen > data.len() {
            return fmt_err("truncated BGZF block");
        }
        let isize = u32::from_le_bytes([
            data[c + blen - 4],
            data[c + blen - 3],
            data[c + blen - 2],
            data[c + blen - 1],
        ]) as usize;
        if isize == 0 {
            break; // EOF marker / empty block
        }
        blocks.push(Block { c_start: c, c_len: blen, u_start: u, u_len: isize });
        c += blen;
        u += isize;
    }
    Ok(blocks)
}

/// Fallback when no `.fai` sidecar exists: decompress (if needed) and scan the
/// FASTA to recover the index.
fn build_fai_by_scan(data: &[u8], bgzf_blocks: &Option<Vec<Block>>) -> Result<Vec<FaiRec>> {
    // Materialise the uncompressed FASTA.
    let text: Vec<u8> = match bgzf_blocks {
        None => data.to_vec(),
        Some(blocks) => {
            let mut out = Vec::new();
            for blk in blocks {
                let cdata = &data[blk.c_start + 18..blk.c_start + blk.c_len - 8];
                out.extend_from_slice(&inflate_raw(cdata, blk.u_len)?);
            }
            out
        }
    };
    // Single pass over lines, tracking byte offsets.
    let mut recs: Vec<FaiRec> = Vec::new();
    let mut pos = 0usize;
    let mut cur: Option<FaiRec> = None;
    for line in text.split_inclusive(|&b| b == b'\n') {
        let start = pos;
        pos += line.len();
        let content = line.strip_suffix(b"\n").unwrap_or(line);
        if content.first() == Some(&b'>') {
            if let Some(r) = cur.take() {
                recs.push(r);
            }
            let name = String::from_utf8_lossy(&content[1..])
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            cur = Some(FaiRec {
                name,
                length: 0,
                offset: start + line.len(),
                linebases: 0,
                linewidth: 0,
            });
        } else if let Some(r) = cur.as_mut() {
            if r.linebases == 0 {
                r.linebases = content.len();
                r.linewidth = line.len();
            }
            r.length += content.len();
        }
    }
    if let Some(r) = cur.take() {
        recs.push(r);
    }
    Ok(recs)
}

pub fn read_all(path: impl AsRef<Path>) -> Result<Vec<Sequence>> {
    let rd = FaidxReader::open(path)?;
    rd.names()
        .iter()
        .map(|n| rd.extract(n, 0, None).map(|b| Sequence::new(n.clone(), b)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seqs() -> Vec<Sequence> {
        vec![
            Sequence::new("chr1", b"ACGTNNNNRYSWacgtKMBDHVacgtACGTACGTAC".to_vec()),
            Sequence::new("chr2", b"AAAACCCCGGGGTTTT".to_vec()),
            Sequence::new("short", b"ACG".to_vec()),
        ]
    }

    fn roundtrip(bgzip: bool, width: usize) {
        let dir = std::env::temp_dir().join(format!("seqfmt_fai_{}_{}", bgzip, width));
        let _ = fs::create_dir_all(&dir);
        let path = dir.join(if bgzip { "t.fa.gz" } else { "t.fa" });
        let s = seqs();
        write_file(&path, &s, width, bgzip, 6).unwrap();
        let rd = FaidxReader::open(&path).unwrap();
        for orig in &s {
            // whole sequence
            assert_eq!(rd.extract(&orig.name, 0, None).unwrap(), orig.bases, "whole {}", orig.name);
            // every sub-window
            let n = orig.bases.len();
            for a in 0..=n {
                for b in a..=n {
                    assert_eq!(
                        rd.extract(&orig.name, a, Some(b)).unwrap(),
                        &orig.bases[a..b],
                        "{} {a}-{b} bgzip={bgzip} width={width}",
                        orig.name
                    );
                }
            }
        }
    }

    #[test]
    fn plain_roundtrip() {
        roundtrip(false, 10);
        roundtrip(false, 0);
    }

    #[test]
    fn bgzf_roundtrip() {
        roundtrip(true, 10);
        roundtrip(true, 60);
    }

    #[test]
    fn bgzf_self_describing_without_fai() {
        // Remove the .fai and ensure the scan fallback reconstructs it.
        let dir = std::env::temp_dir().join("seqfmt_fai_nofai");
        let _ = fs::create_dir_all(&dir);
        let path = dir.join("t.fa.gz");
        let s = seqs();
        write_file(&path, &s, 10, true, 6).unwrap();
        fs::remove_file(sidecar(&path, ".fai")).unwrap();
        let rd = FaidxReader::open(&path).unwrap();
        assert_eq!(rd.extract("chr1", 4, Some(12)).unwrap(), b"NNNNRYSW");
    }
}
