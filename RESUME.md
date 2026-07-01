# Resume notes — UDC (HTTP range-read) for all formats

Snapshot for picking this back up in a fresh (e.g. `screen`) session. See
`STATUS.md` for the full project state; this file is the short "where we are".

## What was done this session
- **Shared `Source`** (`src/source.rs`): `Mem` / `File` / `Http` backing store
  with typed positioned reads. `Http` = pooled `ureq` agent (one kept-alive
  connection) + UDC-style 8 KiB block cache + request/byte counters.
- **All readers are now `Source`-backed** with `open` / `from_vec` / `open_url`
  / `from_source`: `twobit`, `twobyte` (2be), `fourbit`, `samtools`.
- **Seek-based B+ tree**: `bptree::find_src` / `iter_all_src` walk nodes over a
  `Source` (one node per level → ~3 requests remote).
- **Unified URL entry point**: `seqformat::open_url(url) -> Box<dyn SeqReader>`
  (lib.rs) auto-detects format from a 64-byte prefix. `extract <url>` in
  `main.rs` routes here; `--http-stats` prints requests/bytes.
- **Dependency added**: `ureq = "2"` (rustls) in Cargo.toml.
- **Benchmark**: `bench/webseq.sh` (configurable `FORMATS="label=url ..."`).
- Docs updated: README "Serving over the web" section + STATUS.md.

## Verified
- `cargo test --release` → 19 passed.
- Remote == local (byte-identical) for 2bit-idx, 2be, 4bit, faidx (plain+BGZF).
- Web per-fetch (500k×300, real host): std 217 ms / 6.7 MiB, idx 51 ms / 170 KiB,
  **2be 16 ms / 49 KiB** (7 requests — the B+ tree's ~3-node lookup wins).
- 4bit = O(N) open (~10k req / 84 MiB, no index); faidx loads whole `.fai`.

## Test files served (already uploaded to /home/rhubley/public_html)
- `seqformat-std.2bit`, `seqformat-idx.2bit`, `seqformat.2be`, `seqformat.4bit`
  (all 500k×300, seed 1)
- `seqformat.fa.gz` + `.fai`, `seqformat.fa` + `.fai` (50k×300, seed 7)
- Base URL: `https://repeatmasker.org/~rhubley/`

## Re-validate quickly
```sh
cargo build --release && cargo test --release
BIN=./target/release/seqformat B=https://repeatmasker.org/~rhubley
$BIN extract $B/seqformat.2be seq250000 --http-stats >/dev/null   # ~7 req / ~50 KiB
FETCHES=15 bash bench/webseq.sh
```

## Pending / candidate next steps
- **Interrupted**: the final all-formats validation sweep did not finish running
  (was cut off to save state) — re-run the block above to confirm.
- Gate `ureq` behind an optional `http` cargo feature so the default build stays
  std-only + libdeflater (cfg the URL paths + `open_url`).
- Re-run `manyseq.sh` to refresh **local** 2be per-fetch (should now be ~1 ms
  since 2be seeks instead of slurping; old tables say ~37 ms).
- `.fai` fetch for remote faidx uses a *separate* `Source`, so its bytes aren't
  in `http_stats` — consider folding it in for honest accounting.
- Optional: `.gzi`-based BGZF seek (avoid the O(blocks) header scan).

## Git
Work committed on branch `udc-all-formats` (off `main`).
