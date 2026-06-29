//! The in-memory sequence model plus the IUB/IUPAC code tables shared by the
//! twoBit and 4-bit codecs.

/// A single named sequence. `bases` holds the raw residues exactly as they
/// appear in FASTA, i.e. case is preserved (lower-case = soft-masked) and
/// degenerate IUB codes are kept verbatim.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    pub name: String,
    pub bases: Vec<u8>,
}

impl Sequence {
    pub fn new(name: impl Into<String>, bases: impl Into<Vec<u8>>) -> Self {
        Sequence {
            name: name.into(),
            bases: bases.into(),
        }
    }

    pub fn len(&self) -> usize {
        self.bases.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }
}

/// Returns true for the four unambiguous bases (case-insensitive).
pub fn is_acgt(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

/// The ten "specific" IUB degenerate codes — everything that is neither a plain
/// base nor a bare `N`. These are what the standard twoBit format throws away
/// (it stores them all as `N`) and what our extension table preserves.
pub fn is_iub_degenerate(b: u8) -> bool {
    matches!(
        b.to_ascii_uppercase(),
        b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V'
    )
}

/// The BWA/BAM 4-bit nibble alphabet. Index 0..=15 maps to these symbols; this
/// is the same table htslib/seqtk use for `seq_nt16`. It covers every IUB code,
/// so a 4-bit file needs no side tables for degenerate bases.
pub const NIBBLE_ALPHABET: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// Map a base character to its 4-bit nibble. Unknown characters become `N` (15).
pub fn base_to_nibble(b: u8) -> u8 {
    let up = b.to_ascii_uppercase();
    match NIBBLE_ALPHABET.iter().position(|&c| c == up) {
        Some(i) => i as u8,
        None => 15, // N
    }
}

/// Map a 4-bit nibble back to its (upper-case) base character.
pub fn nibble_to_base(n: u8) -> u8 {
    NIBBLE_ALPHABET[(n & 0x0F) as usize]
}

/// 2-bit code for a base, per the twoBit spec: T=0, C=1, A=2, G=3.
/// Anything that isn't A/C/G/T packs as 0 (its real identity is recorded
/// out-of-band in the N-block / IUB tables).
pub fn base_to_twobit(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'T' => 0,
        b'C' => 1,
        b'A' => 2,
        b'G' => 3,
        _ => 0,
    }
}

/// Inverse of [`base_to_twobit`].
pub fn twobit_to_base(code: u8) -> u8 {
    match code & 0x03 {
        0 => b'T',
        1 => b'C',
        2 => b'A',
        _ => b'G',
    }
}
