//! Backing store shared by every random-access reader: an in-memory buffer, a
//! seekable local file, or a remote file fetched with HTTP range requests
//! (seqformat's analogue of UCSC's UDC layer). Each access is a small positioned
//! read of just the header, an index probe, or the requested window — never the
//! whole file — exactly how UCSC `twoBitToFa` (`lseek`+`read`, or UDC over HTTP)
//! touches only what it needs.

use crate::error::{fmt_err, Error, Result};
use std::cell::{Cell, RefCell};
use std::collections::HashMap;
use std::fs;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

pub(crate) enum Source {
    Mem(Vec<u8>),
    File { file: RefCell<fs::File>, len: usize },
    /// Remote file fetched with HTTP range requests — the UDC analogue.
    Http(HttpSource),
}

impl Source {
    /// Open a local file for seek+read (never slurps).
    pub(crate) fn from_path(path: impl AsRef<Path>) -> Result<Source> {
        let file = fs::File::open(path)?;
        let len = file.metadata()?.len() as usize;
        Ok(Source::File { file: RefCell::new(file), len })
    }

    /// Open a remote `http(s)://` file for range reads (UDC-style).
    pub(crate) fn from_url(url: &str) -> Result<Source> {
        Ok(Source::Http(HttpSource::open(url)?))
    }

    pub(crate) fn len(&self) -> usize {
        match self {
            Source::Mem(d) => d.len(),
            Source::File { len, .. } => *len,
            Source::Http(h) => h.len,
        }
    }

    /// The resolved URL for an HTTP source (e.g. to derive a sidecar like `.fai`).
    pub(crate) fn url(&self) -> Option<&str> {
        match self {
            Source::Http(h) => Some(&h.url),
            _ => None,
        }
    }

    /// `(requests, bytes)` issued so far over HTTP; `None` for local/in-memory.
    pub(crate) fn http_stats(&self) -> Option<(u64, u64)> {
        match self {
            Source::Http(h) => Some((h.requests.get(), h.bytes.get())),
            _ => None,
        }
    }

    /// Read exactly `buf.len()` bytes starting at `off`.
    pub(crate) fn read_at(&self, off: usize, buf: &mut [u8]) -> Result<()> {
        match self {
            Source::Mem(d) => {
                let end = off
                    .checked_add(buf.len())
                    .filter(|&e| e <= d.len())
                    .ok_or_else(|| Error::Format(format!(
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
            Source::Http(h) => h.read_at(off, buf),
        }
    }

    pub(crate) fn bytes_at(&self, off: usize, n: usize) -> Result<Vec<u8>> {
        let mut v = vec![0u8; n];
        self.read_at(off, &mut v)?;
        Ok(v)
    }

    pub(crate) fn u8_at(&self, off: usize) -> Result<u8> {
        let mut b = [0u8; 1];
        self.read_at(off, &mut b)?;
        Ok(b[0])
    }

    pub(crate) fn u16_at(&self, off: usize, little: bool) -> Result<u16> {
        let mut b = [0u8; 2];
        self.read_at(off, &mut b)?;
        Ok(if little { u16::from_le_bytes(b) } else { u16::from_be_bytes(b) })
    }

    pub(crate) fn u32_at(&self, off: usize, little: bool) -> Result<u32> {
        let mut b = [0u8; 4];
        self.read_at(off, &mut b)?;
        Ok(if little { u32::from_le_bytes(b) } else { u32::from_be_bytes(b) })
    }

    pub(crate) fn u64_at(&self, off: usize, little: bool) -> Result<u64> {
        let mut b = [0u8; 8];
        self.read_at(off, &mut b)?;
        Ok(if little { u64::from_le_bytes(b) } else { u64::from_be_bytes(b) })
    }
}

/// Default cache/fetch block size for the HTTP source. Matches UCSC UDC's
/// default `udcBlockSize` of 8 KiB so round-trip counts are comparable.
/// Overridable with the `SEQFORMAT_HTTP_BLOCK` environment variable.
const HTTP_BLOCK_DEFAULT: usize = 8192;

/// HTTP range-read transport with a UDC-style block cache. A single pooled
/// `ureq::Agent` keeps **one** connection alive across all range requests — the
/// way UCSC's UDC reuses a connection — so per-fetch cost reflects round-trips +
/// bytes, not repeated TLS handshakes. Each cache miss issues one ranged GET for
/// the covering block(s); hits are free. We tally requests and bytes so a
/// benchmark can see how the access pattern maps to network cost.
pub(crate) struct HttpSource {
    agent: ureq::Agent,
    url: String, // resolved (post-redirect) URL
    len: usize,
    block: usize,
    cache: RefCell<HashMap<usize, Vec<u8>>>,
    requests: Cell<u64>,
    bytes: Cell<u64>,
}

impl HttpSource {
    /// Probe the URL once: follow redirects, learn the total length from the
    /// `Content-Range` of a 1-byte ranged GET, and remember the final URL so
    /// subsequent block fetches skip the redirect and reuse the connection.
    fn open(url: &str) -> Result<Self> {
        let block = std::env::var("SEQFORMAT_HTTP_BLOCK")
            .ok()
            .and_then(|s| s.parse().ok())
            .filter(|&b| b > 0)
            .unwrap_or(HTTP_BLOCK_DEFAULT);

        let agent = ureq::AgentBuilder::new().build();
        let resp = agent
            .get(url)
            .set("Range", "bytes=0-0")
            .call()
            .map_err(|e| Error::Format(format!("{url}: {e}")))?;
        let eff = resp.get_url().to_string();
        let len = resp
            .header("Content-Range")
            .and_then(|v| v.rsplit('/').next())
            .and_then(|n| n.trim().parse::<usize>().ok())
            .ok_or_else(|| Error::Format(format!(
                "{url}: server did not report a Content-Range; range requests unsupported?"
            )))?;
        let _ = std::io::copy(&mut resp.into_reader(), &mut std::io::sink());
        Ok(HttpSource {
            agent,
            url: eff,
            len,
            block,
            cache: RefCell::new(HashMap::new()),
            requests: Cell::new(0),
            bytes: Cell::new(0),
        })
    }

    /// Fetch blocks `[a, b)` in a *single* range request and split the result
    /// into per-block cache entries. Coalescing a contiguous run into one
    /// request mirrors UDC streaming a run of blocks over one kept-alive
    /// connection — so a big sequential read costs one request of many bytes,
    /// not one request per block.
    fn fetch_run(&self, a: usize, b: usize) -> Result<()> {
        let start = a * self.block;
        let end = (b * self.block).min(self.len); // exclusive
        if start >= end {
            for idx in a..b {
                self.cache.borrow_mut().insert(idx, Vec::new());
            }
            return Ok(());
        }
        let range = format!("bytes={start}-{}", end - 1);
        let resp = self
            .agent
            .get(&self.url)
            .set("Range", &range)
            .call()
            .map_err(|e| Error::Format(format!("{}: {e}", self.url)))?;
        let mut data = Vec::with_capacity(end - start);
        resp.into_reader().read_to_end(&mut data).map_err(Error::Io)?;
        if data.len() != end - start {
            return fmt_err(format!(
                "{}: range {range} returned {} bytes, expected {}",
                self.url,
                data.len(),
                end - start
            ));
        }
        self.requests.set(self.requests.get() + 1);
        self.bytes.set(self.bytes.get() + data.len() as u64);
        let mut cache = self.cache.borrow_mut();
        for idx in a..b {
            let lo = idx * self.block - start;
            let hi = ((idx + 1) * self.block).min(end) - start;
            cache.insert(idx, data[lo..hi].to_vec());
        }
        Ok(())
    }

    fn read_at(&self, off: usize, buf: &mut [u8]) -> Result<()> {
        let end = off.checked_add(buf.len()).filter(|&e| e <= self.len).ok_or_else(|| {
            Error::Format(format!(
                "read of {} bytes at offset {off} is past the {}-byte remote file",
                buf.len(),
                self.len
            ))
        })?;
        if buf.is_empty() {
            return Ok(());
        }
        let b0 = off / self.block;
        let b1 = (end - 1) / self.block; // inclusive
        // Fetch each maximal run of missing blocks in one request.
        let mut idx = b0;
        while idx <= b1 {
            if self.cache.borrow().contains_key(&idx) {
                idx += 1;
                continue;
            }
            let run_start = idx;
            while idx <= b1 && !self.cache.borrow().contains_key(&idx) {
                idx += 1;
            }
            self.fetch_run(run_start, idx)?;
        }
        // Assemble from the (now-cached) blocks.
        let cache = self.cache.borrow();
        let mut pos = off;
        while pos < end {
            let bidx = pos / self.block;
            let blk = &cache[&bidx];
            let in_blk = pos - bidx * self.block;
            let take = (blk.len() - in_blk).min(end - pos);
            if take == 0 {
                return fmt_err(format!("short remote block at offset {pos}"));
            }
            buf[pos - off..pos - off + take].copy_from_slice(&blk[in_blk..in_blk + take]);
            pos += take;
        }
        Ok(())
    }
}
