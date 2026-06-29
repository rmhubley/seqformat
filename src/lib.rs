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
pub mod twobit;
pub mod fourbit;
pub mod samtools;
pub mod bptree;
pub mod twobyte;
pub mod generate;

pub use error::{Error, Result};
pub use seq::Sequence;
