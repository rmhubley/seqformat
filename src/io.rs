//! Tiny binary read/write helpers.
//!
//! twoBit files are written in the machine's native byte order, with the
//! file signature acting as a byte-order mark: a reader that sees the signature
//! byte-swapped must swap every multi-byte field. We always *write* little-endian
//! (by far the common case), but we can *read* either endianness.

use crate::error::{fmt_err, Result};

/// A forward/seekable cursor over an in-memory buffer that decodes integers
/// with a fixed endianness.
pub struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
    little: bool,
}

impl<'a> Reader<'a> {
    pub fn new(data: &'a [u8], little: bool) -> Self {
        Reader { data, pos: 0, little }
    }

    pub fn pos(&self) -> usize {
        self.pos
    }

    pub fn seek(&mut self, pos: usize) {
        self.pos = pos;
    }

    pub fn remaining(&self) -> usize {
        self.data.len().saturating_sub(self.pos)
    }

    /// Borrow the next `n` bytes and advance.
    pub fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        if self.pos + n > self.data.len() {
            return fmt_err(format!(
                "unexpected end of data: wanted {n} bytes at offset {}, only {} remain",
                self.pos,
                self.remaining()
            ));
        }
        let s = &self.data[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }

    pub fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }

    pub fn u32(&mut self) -> Result<u32> {
        let b = self.take(4)?;
        let a = [b[0], b[1], b[2], b[3]];
        Ok(if self.little {
            u32::from_le_bytes(a)
        } else {
            u32::from_be_bytes(a)
        })
    }

    pub fn u64(&mut self) -> Result<u64> {
        let b = self.take(8)?;
        let a = [b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]];
        Ok(if self.little {
            u64::from_le_bytes(a)
        } else {
            u64::from_be_bytes(a)
        })
    }
}

/// Read a u32 at an absolute offset (used to peek at magics/trailers).
pub fn peek_u32(data: &[u8], off: usize, little: bool) -> Option<u32> {
    let b = data.get(off..off + 4)?;
    let a = [b[0], b[1], b[2], b[3]];
    Some(if little {
        u32::from_le_bytes(a)
    } else {
        u32::from_be_bytes(a)
    })
}

pub fn peek_u16(data: &[u8], off: usize, little: bool) -> Option<u16> {
    let b = data.get(off..off + 2)?;
    let a = [b[0], b[1]];
    Some(if little {
        u16::from_le_bytes(a)
    } else {
        u16::from_be_bytes(a)
    })
}

pub fn peek_u64(data: &[u8], off: usize, little: bool) -> Option<u64> {
    let b = data.get(off..off + 8)?;
    let mut a = [0u8; 8];
    a.copy_from_slice(b);
    Some(if little {
        u64::from_le_bytes(a)
    } else {
        u64::from_be_bytes(a)
    })
}

// ---- writer helpers (always little-endian) ----

pub fn put_u8(v: &mut Vec<u8>, x: u8) {
    v.push(x);
}

pub fn put_u16(v: &mut Vec<u8>, x: u16) {
    v.extend_from_slice(&x.to_le_bytes());
}

pub fn put_u32(v: &mut Vec<u8>, x: u32) {
    v.extend_from_slice(&x.to_le_bytes());
}

pub fn put_u64(v: &mut Vec<u8>, x: u64) {
    v.extend_from_slice(&x.to_le_bytes());
}
