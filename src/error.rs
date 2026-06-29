//! Error type shared across the crate.

use std::fmt;

/// Crate result alias.
pub type Result<T> = std::result::Result<T, Error>;

/// Anything that can go wrong while parsing/serialising a sequence file.
#[derive(Debug)]
pub enum Error {
    /// Wraps an underlying I/O failure.
    Io(std::io::Error),
    /// The bytes don't match the expected format (bad magic, truncation, ...).
    Format(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(e) => write!(f, "I/O error: {e}"),
            Error::Format(m) => write!(f, "format error: {m}"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Io(e) => Some(e),
            Error::Format(_) => None,
        }
    }
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e)
    }
}

/// Convenience for building a [`Error::Format`].
pub(crate) fn fmt_err<T>(msg: impl Into<String>) -> Result<T> {
    Err(Error::Format(msg.into()))
}
