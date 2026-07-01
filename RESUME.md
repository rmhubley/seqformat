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

## Done since first checkpoint
- **New format `--bpt`** (`fa2twobit --iub --bpt`): 2bit + IUB + a full B+ tree
  appended as a *complete TOC duplicate* (`BPT_FOOTER_MAGIC`, names inline).
  Fully twoBit backward compatible (verified `twoBitToFa` reads it). Reader uses
  `bptree::find_src` (seek/HTTP). Web: **18 ms / 6.6 req / 52 KiB — matches 2be**
  while staying compatible, vs idx 83/22/170. Local ties idx/2be (~1.3 ms).
  Size 75.8 MiB (+4.3 over idx = duplicated names). Caveat baked into docs: bench
  names are `seq+index` (4–9 B) → best case; realistic names ~3× the size gap and
  ~10–13 req. 20 tests pass. In manyseq.sh + webseq.sh. User declined a
  realistic-name rerun for now.
- Final all-formats validation sweep: **passed** (remote==local for all).
- `cmd_extract` now routes single-region **2be** through the seek path (was
  slurping) → 2be local per-fetch **37 ms → 1.4 ms**; BGZF few-big **95 → 8 ms**.
- Remote **`.fai` fetch folded into `http_stats`** (FaidxReader.fai_http) — faidx
  web cost now honestly includes the index load (~1.37 MiB plain 50k).
- Benchmarks re-run (manyseq, benchmark, webseq); README + STATUS tables refreshed.

## Pending / candidate next steps
- Gate `ureq` behind an optional `http` cargo feature so the default build stays
  std-only + libdeflater (cfg the URL paths + `open_url`).
- Minor: 2be *bulk* (slurp path) ~0.5 s vs ~0.1 s for 2bit family —
  `Source::Mem` allocates per positioned read. Seek path unaffected. A Mem
  fast-path returning slices would fix it.
- Optional: `.gzi`-based BGZF seek (avoid the O(blocks) header scan).

## Git
Work committed on branch `udc-all-formats` (off `main`).
