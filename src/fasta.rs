//! Minimal multi-record FASTA reader/writer (used for testing & conversion).

use crate::error::{fmt_err, Result};
use crate::seq::Sequence;
use std::fs;
use std::path::Path;

/// Parse FASTA from a byte buffer. The record name is the header text up to the
/// first run of whitespace; the remainder of the header line is ignored. Case
/// and degenerate codes are preserved verbatim.
pub fn parse(data: &[u8]) -> Result<Vec<Sequence>> {
    let mut seqs: Vec<Sequence> = Vec::new();
    let mut cur: Option<Sequence> = None;

    for raw in data.split(|&b| b == b'\n') {
        // Trim a trailing '\r' (CRLF files) and surrounding ASCII whitespace.
        let line = trim(raw);
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            if let Some(s) = cur.take() {
                seqs.push(s);
            }
            let header = &line[1..];
            let name_end = header
                .iter()
                .position(|b| b.is_ascii_whitespace())
                .unwrap_or(header.len());
            let name = String::from_utf8_lossy(&header[..name_end]).into_owned();
            cur = Some(Sequence::new(name, Vec::new()));
        } else {
            match cur.as_mut() {
                Some(s) => {
                    // Drop any internal whitespace just in case.
                    s.bases
                        .extend(line.iter().copied().filter(|b| !b.is_ascii_whitespace()));
                }
                None => return fmt_err("FASTA data does not start with a '>' header"),
            }
        }
    }
    if let Some(s) = cur.take() {
        seqs.push(s);
    }
    Ok(seqs)
}

fn trim(mut s: &[u8]) -> &[u8] {
    while let [first, rest @ ..] = s {
        if first.is_ascii_whitespace() {
            s = rest;
        } else {
            break;
        }
    }
    while let [rest @ .., last] = s {
        if last.is_ascii_whitespace() {
            s = rest;
        } else {
            break;
        }
    }
    s
}

/// Serialise sequences to FASTA with the given line width (0 = no wrapping).
pub fn write_bytes(seqs: &[Sequence], width: usize) -> Vec<u8> {
    let mut out = Vec::new();
    for s in seqs {
        out.push(b'>');
        out.extend_from_slice(s.name.as_bytes());
        out.push(b'\n');
        if width == 0 || s.bases.is_empty() {
            out.extend_from_slice(&s.bases);
            out.push(b'\n');
        } else {
            for chunk in s.bases.chunks(width) {
                out.extend_from_slice(chunk);
                out.push(b'\n');
            }
        }
    }
    out
}

pub fn read_file(path: impl AsRef<Path>) -> Result<Vec<Sequence>> {
    parse(&fs::read(path)?)
}

pub fn write_file(path: impl AsRef<Path>, seqs: &[Sequence], width: usize) -> Result<()> {
    fs::write(path, write_bytes(seqs, width))?;
    Ok(())
}
