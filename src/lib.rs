//! `seqformat` — a small, dependency-free library for reading and writing
//! genomic sequence container formats:
//!
//! * **twoBit** (UCSC `.2bit`) — both the *standard* (version 0, 32-bit offsets)
//!   and *long* (version 1, 64-bit offsets) variants.
//! * **twoBit + IUB extension** — a backward-compatible superset of twoBit that
//!   additionally records the exact IUB/IUPAC degenerate codes (R, Y, S, W, K,
//!   M, B, D, H, V) that the stock format collapses into plain `N`. Old readers
//!   are completely unaware of the extra data (see [`twobit`] for the on-disk
//!   layout and why it is invisible to legacy tools).
//! * **4-bit** — a BWA/BAM-style 4-bits-per-base packing (`=ACMGRSVTWYHKDBN`
//!   nibble code) that natively represents every IUB code including `N`.
//!
//! FASTA reading/writing is provided for testing and round-tripping.

pub mod error;
pub mod io;
pub mod seq;
pub mod fasta;
pub(crate) mod source;
pub mod twobit;
pub mod fourbit;
pub mod samtools;
pub mod bptree;
pub mod twobyte;
pub mod generate;

pub use error::{Error, Result};
pub use seq::Sequence;

/// A uniform random-access view over any supported container, used by the URL
/// path so callers don't dispatch on format themselves.
pub trait SeqReader {
    /// Sequence names in file order.
    fn names(&self) -> Vec<String>;
    /// Decode `[start, end)` of `name` (0-based, half-open; `end == None` = to end).
    fn extract(&self, name: &str, start: usize, end: Option<usize>) -> Result<Vec<u8>>;
    /// `(requests, bytes)` if reading over HTTP; `None` otherwise.
    fn http_stats(&self) -> Option<(u64, u64)>;
}

/// Open any supported container served over `http(s)://` for UDC-style range
/// reads, auto-detecting the format from a small prefix. One [`source::Source`]
/// is built, peeked, then handed to the matching reader (no re-probe).
pub fn open_url(url: &str) -> Result<Box<dyn SeqReader>> {
    let src = source::Source::from_url(url)?;
    let prefix = src.bytes_at(0, src.len().min(64))?;
    if twobit::is_twobit(&prefix) {
        Ok(Box::new(twobit::TwoBitReader::from_source(src)?))
    } else if twobyte::is_twobyte(&prefix) {
        Ok(Box::new(twobyte::TwoByteReader::from_source(src)?))
    } else if fourbit::is_fourbit(&prefix) {
        Ok(Box::new(fourbit::FourBitReader::from_source(src)?))
    } else if samtools::is_bgzf(&prefix) || samtools::is_fasta(&prefix) {
        Ok(Box::new(samtools::FaidxReader::from_source(src)?))
    } else {
        Err(Error::Format(format!("{url}: unrecognized sequence format")))
    }
}
