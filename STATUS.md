# seqformat — project status

Checkpoint: 2026-06-29. Build clean, **19 tests pass** (`cargo test`), release
binary at `target/release/seqformat`. Not a git repo.

## What this is
A dependency-light Rust tool + library for genomic sequence container formats,
with FASTA conversion and a benchmark suite. Only external dep is `libdeflater`
(C `libdeflate`, used solely by the samtools/BGZF module — needs a C compiler).

## Formats implemented (all round-trip + cross-validated against UCSC/samtools)
- **twoBit standard** (v0) and **long** (v1, 64-bit offsets). Byte-identical to
  `faToTwoBit`. `src/twobit.rs`.
- **twoBit + IUB extension** — backward-compatible: degenerate codes (R Y S W K M
  B D H V) preserved in a per-record table after packedDNA + an EOF magic
  trailer; old readers see `N`. Double-coded (also in N-block table). `--iub`.
- **twoBit + IUB + name index** (`--iub --index`, or `--index` alone) — adds a
  backward-compatible sorted **pointer array** (one u64/seq, sorted by name,
  pointing at the existing flat-TOC entry) + a 24-byte EOF footer that subsumes
  the IUB trailer. Names are NOT duplicated. Gives O(log N) binary-search lookup
  (reader skips the O(N) flat-TOC load). Cross-validated: UCSC `twoBitToFa`/
  `twoBitInfo` read it as plain 2bit. `to_bytes_indexed`/`write_file_indexed`,
  `TwoBitReader::lookup_indexed`. manyseq per-fetch: 3.6 ms vs flat 12.5 ms,
  ~matches 2be's B+ tree (3.3 ms).
- **4-bit** — BWA/BAM `seq_nt16` nibble packing. `src/fourbit.rs`.
- **samtools format** — FASTA + `.fai`, optionally BGZF `.fa.gz` + `.fai`/`.gzi`.
  Our `.fai`/`.gzi` are **byte-identical to samtools**; BGZF output byte-identical
  to `bgzip`. `src/samtools.rs`.
- **2be** (experimental, no backward compat) — B+ tree TOC (`src/bptree.rs`) for
  O(log N) name lookup + per-sequence merged tagged-edit stream (N_RUN /
  IUB_POINT / IUB_RUN / MASK_RUN) + run-index. `src/twobyte.rs`.

## CLI (`src/main.rs`)
`fa2twobit`/`twobit2fa` (`--long`,`--iub`,`--index`), `fa2fourbit`/`fourbit2fa`,
`fa2faidx` (`--bgzip`), `fa2be`/`be2fa`, `extract` (region or `--seq-list`,
auto-detects format), `random` (test data; `--n-runs`/`--iub-runs` cluster
ambiguity), `info`.

## Benchmarks (`bench/`, gzip removed)
- `benchmark.sh` — few big seqs (default 3×100 Mbp, N clustered `N_RUNS=3`, IUB
  scattered `IUB_RUNS=0`). Reports size/bits-base, encode, bulk-extract,
  per-fetch latency.
- `manyseq.sh` — many short seqs (default 500k×300 bp). Same metrics.
- `webseq.sh` — **remote** per-fetch over HTTP range requests (the UDC case).
  Fetches single seqs from files served by URL and reports ms/fetch + HTTP
  requests + bytes/fetch. Configurable `FORMATS="label=url ..."`; defaults to
  the 500k `2bit std` / `2bit idx` / `2be` comparison.
- Both report the same columns/labels (`2bit (UCSC kentsrc)`, `BGZF (samtools)`,
  etc.); latest numbers + analysis are in README.

## Key findings
- Clustered-N genome data: standard 2bit is the most compact (2.315 b/base) and
  beats BGZF; with scattered IUB, `twoBitToFa` bulk extract is slow (~17 s/20k)
  because it reads the per-fragment N-block list — our binary-searched readers
  stay ~0.13 s.
- The **2bit-family reader is now seek+read** (`TwoBitReader::open`, a `Source`
  enum: `RefCell<File>` seek/`read_exact` vs in-mem `Vec`; `from_vec` keeps the
  mem path). A single fetch reads only header + index probe + window — never the
  whole file. Few-big per-fetch: ~45–53 ms → **~1.1 ms**, now *beating*
  `twoBitToFa` (8.2 ms). `cmd_extract` peeks a 64-byte prefix (`twobit::is_twobit`)
  to pick the seek path, and slurps via `from_vec` only for >1024 regions /
  extract-all (so bulk-20k stays ~0.12 s). Only 4bit still slurps (no index);
  2be/BGZF/samtools now seek too — see the Remote/UDC bullet below.
- **2bit+IUB+index and 2be tie the many-short per-fetch** at **1.4 ms** vs
  flat-TOC 173 ms, `twoBitToFa` 89 ms, samtools 245 ms. (2be was 37 ms when it
  slurped; now that it seeks, its B+ tree lands at 1.4 ms — the pointer array and
  the tree are a wash locally, the tree wins remotely. See the Remote/UDC bullet.)
- **Remote/UDC (new)**: a shared `Source` (`src/source.rs`; Mem/File/Http) gives
  **every** reader an HTTP range-read path — `ureq` agent with a pooled
  connection + UDC-style 8 KiB block cache. `seqformat::open_url()` auto-detects
  format and returns a `Box<dyn SeqReader>`; `extract <url> --http-stats` reports
  requests/bytes. Web per-fetch on 500k×300 (served from repeatmasker.org):

  | format | ms/fetch | req/fetch | KiB/fetch |
  |---|--:|--:|--:|
  | 2bit std (flat TOC)          | 242 | 3.9  | 6743 |
  | 2bit idx (sorted array)      |  83 | 21.7 |  170 |
  | **2bit +bptree (`--bpt`)**   | **18** | **6.6** | **52** |
  | 2be (B+ tree)                |  20 | 7.0  |   49 |
  | 4bit (no index) ¹            |  —  | 10241 | 84 MiB |
  | faidx plain (`.fai`) ¹       |  —  | 3    | 14.6 MiB |
  | faidx bgzf ¹                 |  —  | 2431 | 34.5 MiB |

  ¹ single fetch (O(N) on open; a timed loop would be dominated by these). The
  story from the theory holds: flat TOC pulls O(N) bytes (whole TOC) every open;
  the sorted index does ~log₂N *scattered* probes (~15–25 requests); the **B+
  tree (fan-out 256, `bptree::find_src`) hits the ideal ~3-node lookup**. The new
  **`--bpt`** format (2bit + IUB + a full B+ tree appended as a TOC duplicate,
  `BPT_FOOTER_MAGIC`) **matches 2be remotely while staying twoBit backward
  compatible** — at the cost of duplicated names (+4.3 MiB over `--index` on 500k
  short names; more on realistic names). 4bit scans every record header (~whole
  file); faidx pulls its whole `.fai` (14.6 MiB at 500k — folded into
  `--http-stats`) then a window read / BGZF block scan. `bench/webseq.sh`.
  NOTE: names are `seq0…seq499999` (4–9 B) → best case for `--bpt`'s size + node
  packing; realistic 15–25 B names widen the gap ~3× and push tree nodes toward
  2 blocks/level (~7 → ~10–13 req).

  **Local per-fetch** now matches: `cmd_extract` routes single-region 2be through
  the seek path too (was slurping), so 2be local dropped **37 ms → 1.4 ms**,
  tying the sorted index (1.4 ms); BGZF (seqformat) few-big dropped **95 → 8 ms**.
  4bit stays slurp-based (no index → one sequential read beats O(N) seeks).

## Known limitations / candidate next steps
- **All readers are now Source-backed (seek + HTTP range reads)** — `twobit`,
  `twobyte` (2be), `fourbit`, and `samtools`/BGZF share `src/source.rs`. The old
  whole-file slurp on open is gone; local single-fetch is seek-based and remote
  is range-based.
- HTTP `Source` adds a `ureq` (rustls) dependency; could be gated behind an
  optional `http` cargo feature to keep the default build std-only + libdeflater.
  Remote is `http(s)` only (no `ftp`/`s3`). faidx requires a `.fai` sidecar
  remotely (the scan fallback would defeat range access); remote BGZF scans block
  headers (O(blocks)) rather than reading a `.gzi`.
- Minor: 2be *bulk* (>1024 regions, slurp path) is ~0.5 s vs ~0.1 s for the 2bit
  family — `Source::Mem` allocates a Vec per positioned read, so the many-edit 2be
  decode does more small allocs. Per-fetch (seek) path is unaffected.
- Flat-TOC (non-indexed) 2bit still builds a `String` HashMap of all names when a
  name is looked up without the index → 173 ms on 500k seqs. The `--index` format
  removes it; plain 2bit could also lazily binary-search a seek-read TOC.
- 2be is a prototype: no `--long`-style 64-bit sequence-local coords (u32 caps a
  single sequence at 4 Gbp), no mask/IUB-only selective indexes.
- `cmd_extract` slurp/seek cutoff is a fixed region count (>1024); a size-aware
  threshold would avoid slurping a multi-GB file for a mid-size batch.
- Possible: clustered-IUB single-code mode in the generator.
