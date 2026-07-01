# seqformat

A small Rust tool + library for reading and writing several genomic sequence
container formats, with FASTA conversion for testing. It supports:

- **twoBit standard** (version 0, 32-bit index offsets)
- **twoBit long** (version 1, 64-bit index offsets — the UCSC `faToTwoBit -long`
  format)
- **twoBit + IUB extension** — a *backward-compatible* superset that preserves
  the exact IUB/IUPAC degenerate codes (`R Y S W K M B D H V`) that stock twoBit
  throws away as `N`
- **twoBit + IUB + name index** — the IUB superset *plus* a backward-compatible
  **sorted pointer array** at EOF that gives O(log N) name lookup (binary search,
  no full index load) without duplicating any names. Still reads as plain 2bit in
  UCSC tools. `--index`; see below.
- **4-bit** — BWA/BAM `seq_nt16` nibble packing (`=ACMGRSVTWYHKDBN`), which
  represents every IUB code (incl. `N`) directly
- **samtools format** — FASTA indexed by `.fai`, optionally **BGZF**-compressed
  (`.fa.gz`) with a `.gzi`. The writer produces `.fai`/`.gzi` **byte-identical to
  samtools**, and reads/writes interoperate with `samtools faidx` both ways.
- **2be** (experimental) — a from-scratch twoBit redesign (no backward compat):
  a **B+ tree TOC** for O(log N) name lookup at high sequence counts, and a
  per-sequence **merged tagged-edit stream** (N runs, IUB points, IUB runs, mask
  runs) so isolated degenerate codes cost ~one entry each instead of polluting a
  gap table. See `src/twobyte.rs`.

The twoBit/4bit/FASTA codecs are pure-std. The samtools/BGZF module needs a
DEFLATE codec, so the crate depends on `libdeflater` — Rust bindings to the C
`libdeflate` library that htslib itself uses — for that one module. This needs a
C compiler at build time, and makes our BGZF output **byte-identical to
`bgzip`**.

## Build & test

```sh
cargo build --release
cargo test
./target/release/seqformat help
```

The only dependency is `libdeflater` (used solely by the samtools/BGZF module,
and requiring a C compiler); the twoBit, 4bit and FASTA codecs are
standard-library only.

## Usage

```
seqformat fa2twobit  <in.fa> <out.2bit> [--long] [--iub] [--index | --bpt]
seqformat twobit2fa  <in.2bit> <out.fa> [--width N]
seqformat fa2fourbit <in.fa> <out.4bit>
seqformat fourbit2fa <in.4bit> <out.fa> [--width N]
seqformat fa2faidx   <in.fa> <out.fa> [--bgzip] [--width N] [--level N]
seqformat extract    <file> [region | --seq-list <f>] [--out <fa>]  # 2bit/4bit/fasta/bgzf
seqformat info       <file>            # auto-detects 2bit / 4bit / fasta / bgzf
```

`fa2faidx` writes `out.fa` + `out.fa.fai`; with `--bgzip` it writes a
BGZF-compressed `out.fa.gz` + `.fai` + `.gzi`. `extract` and `info` auto-detect
2bit, 4bit, plain FASTA (`.fai`) and BGZF FASTA (`.fai`/`.gzi`).

Example round-trips:

```sh
seqformat fa2twobit genome.fa genome.2bit            # standard
seqformat fa2twobit genome.fa genome.2bit --long     # 64-bit offsets
seqformat fa2twobit genome.fa genome.2bit --iub      # preserve degenerate codes
seqformat fa2twobit genome.fa genome.2bit --iub --index  # + O(log N) name index
seqformat twobit2fa genome.2bit out.fa
seqformat info genome.2bit
```

## Benchmarking

`bench/benchmark.sh` compares the formats on randomly generated data (default:
3 × 100 Mbp, ~1% `N` **clustered into 3 assembly-gap-style runs** per sequence,
plus ~0.5% IUB codes **scattered** as isolated bases — controlled by `N_RUNS` /
`IUB_RUNS`) against the UCSC tools and `samtools faidx`. It reports
(1) storage size, (2) encoding time, and (3) random-access extraction time of
many regions, with a warmup plus several timed iterations (min/median/mean/sd).
It also cross-checks correctness:

- our standard 2bit ≡ what UCSC `twoBitToFa` reads (degenerates as `N`)
- our **2bit+IUB+index** still reads as a plain 2bit through `twoBitToFa`
  (degenerates as `N`) — the index footer is invisible to it
- our IUB 2bit and 4bit ≡ what `samtools faidx` reads from the source (codes
  preserved)

```sh
bash bench/benchmark.sh                      # 3 x 100 Mbp (several minutes)
SEQS=2 LEN=1000000 bash bench/benchmark.sh   # quick run
bash bench/manyseq.sh                        # 500k x 300 bp (default) — sequence-count axis
```

Requires `faToTwoBit`, `twoBitToFa`, `gzip`, `bgzip`, `samtools` on `PATH`
(missing tools are skipped).

### Representative results (3 × 100 Mbp; ~1% N in 3 runs/seq + ~0.5% IUB scattered)

N is clustered into 3 assembly-gap runs; IUB codes are scattered as isolated
bases (~490 k/seq). Because standard twoBit stores **every** non-ACGT position as
an N-block, the scattered IUB dominate the N-block table (~490 k blocks/seq).

| format | size | bits/base | encode (median) | bulk 20k (1 proc) | per-fetch (separate procs) |
|---|---:|---:|---:|---:|---:|
| FASTA (raw) | 290.9 MiB | 8.133 | — | — | — |
| 2bit (UCSC kentsrc) | 82.8 MiB | 2.315 | 2.54 s | 19.7 s ¹ | 8.8 ms/fetch |
| BGZF (samtools) | 90.8 MiB | 2.538 | 9.96 s | 5.93 s | 5.5 ms/fetch |
| BGZF (seqformat) | 90.8 MiB | 2.538 | 8.72 s ² | 2.33 s | 7.9 ms/fetch ⁵ |
| **2bit (seqformat)** | **82.8 MiB** | **2.315** | 2.39 s | **0.117 s** | **1.3 ms/fetch** ⁵ |
| 2bit+IUB (seqformat) | 95.6 MiB | 2.672 | 3.52 s | 0.136 s | **1.3 ms/fetch** ⁵ |
| 2bit+IUB+index (seqformat) | 95.6 MiB ⁴ | 2.672 | 3.53 s | 0.141 s | **1.3 ms/fetch** ⁵ |
| 2bit+IUB+bptree (seqformat) | 95.6 MiB ⁴ | 2.672 | 3.74 s | 0.141 s | **1.4 ms/fetch** ⁵ |
| 4bit (seqformat) | 143.1 MiB | 4.000 | 1.57 s | 0.134 s | 81.1 ms/fetch ³ |
| 2be (seqformat) | 85.7 MiB | 2.396 | 2.60 s | 0.512 s | **1.9 ms/fetch** ⁵ |

Notes:
- Even with ~490 k IUB N-blocks per sequence, **standard 2bit is still the
  smallest format** (2.315 bits/base) — beating BGZF (2.538) — because the `N`
  gaps are clustered and the bases pack at 2 bits.
- ¹ `2bit (UCSC kentsrc)` = the `faToTwoBit`/`twoBitToFa` tools. In **bulk**,
  `twoBitToFa` reads the full block list per fragment, so the ~490 k scattered IUB
  blocks slow it to ~1 ms/region (19.7 s for 20 k); our binary-searched reader
  stays at **0.117 s (~170× faster)**. Clustering IUB too (`IUB_RUNS=3`) or a
  realistic sparse IUB rate collapses this to ~6 blocks/seq and `twoBitToFa`
  becomes sub-second.
- ³ **4bit is the only slurper left, by design.** It carries no offset table, so
  a lookup must touch every record header — one sequential whole-file read beats
  O(N) scattered seeks, so the reader slurps 143 MB (81 ms single-fetch). Every
  other seqformat reader now seeks (note ⁵).
- ⁴ Both name indexes are negligible at N=3: the sorted array adds 8 bytes/seq +
  a 24-byte footer (40 bytes total), and the `--bpt` B+ tree adds only ~40 bytes
  more (the tree over 3 names) — so both are `95.6 MiB` / `2.672`, unchanged to
  two decimals, and identical to the plain `2bit+IUB`. With so few sequences there
  is no O(N) index load to skip, so neither helps nor hurts per-fetch here; the
  index *structure* only separates them in the many-short and web tables, where
  the duplicated-name cost of `--bpt` (and its remote payoff) actually appear.
- ⁵ **All readers are now `Source`-backed `seek`+`read`** (like `twoBitToFa`'s
  `lseek`/`read`; `src/source.rs`): a single fetch reads only the header, the
  index probe path, and the requested window — never the whole file. This is the
  change that dropped `BGZF (seqformat)` from **95 → 7.9 ms** and `2be` from
  **46 → 1.9 ms** (both used to slurp); the 2bit family sits at ~1.3 ms, faster
  than `twoBitToFa` itself (8.8 ms — our process startup is lighter). Bulk
  extraction (>1024 regions) amortizes a single whole-file read instead. The same
  `Source` also does HTTP range reads — see "Serving over the web" below.
- Everything cross-validates: our standard 2bit is **byte-identical** to
  `faToTwoBit`'s, our BGZF + `.fai`/`.gzi` are **byte-identical** to `bgzip`/
  `samtools`', and `samtools faidx` agrees with our 2bit+IUB / 4bit / BGZF / 2be
  decode.
- ² `BGZF (seqformat)` uses `libdeflate` (the same codec as htslib `bgzip`), so
  output is **byte-identical to `bgzip`**; encode matches (~9 s) and bulk
  extraction beats `samtools faidx` (2.33 s vs 5.93 s).

### Many short sequences — the sequence-count axis (`bench/manyseq.sh`)

`bench/manyseq.sh` stresses sequence *count* instead of size — **500,000 × 300 bp**
here — and measures build, size, bulk extraction, and per-fetch latency. This is
where a name index earns its keep: O(log N) lookup instead of loading a
500k-entry table on every open.

Note `bits/base` is *higher* than the few-big scenario for the same formats —
short sequences carry more per-record overhead (index entry, record header, block
tables) relative to their 300 bp of data; that overhead is exactly what differs
between the two tables.

| format | size | bits/base | build (median) | bulk 20k (1 proc) | **per-fetch (separate procs)** |
|---|---:|---:|---:|---:|---:|
| FASTA (raw) | 150.6 MiB | 8.421 | — | — | — |
| 2bit (UCSC kentsrc) | 59.4 MiB | 3.320 | 1.69 s | 0.151 s | 89.2 ms/fetch |
| BGZF (samtools) | 49.2 MiB | 2.754 | 5.43 s | 0.975 s | 244.8 ms/fetch ³ |
| 2bit (seqformat, flat TOC) | 59.4 MiB | 3.320 | 1.60 s | 0.276 s | 173.1 ms/fetch ⁶ |
| **2bit+IUB+index (seqformat, ptr array)** | 71.5 MiB | 3.997 | 2.42 s | **0.107 s** | **1.4 ms/fetch** ⁵ |
| **2bit+IUB+bptree (seqformat, B+ tree)** | 75.8 MiB | 4.239 | 2.62 s | **0.114 s** | **1.3 ms/fetch** ⁵ |
| 4bit (seqformat, no index) | 80.0 MiB | 4.474 | 0.94 s | 0.314 s | 244.6 ms/fetch ³ |
| **2be (seqformat, B+ tree)** | 63.4 MiB | 3.545 | 1.72 s | **0.104 s** | **1.4 ms/fetch** ⁵ |

The per-fetch column (open + fetch one sequence in a *separate* process — the
real "grab a contig from a huge multi-FASTA" pattern) is the headline, and all
three indexed formats **tie at ~1.4 ms** — the sorted pointer array, the appended
B+ tree, and 2be's B+ tree all make the name lookup O(log N) (no 500k-entry TOC
to load), and the `Source` seek+read reader (note ⁵) means the lookup probes and
the one record are the only bytes touched. That's ~120× faster than our own
flat-TOC 2bit reader, ~64× faster than `twoBitToFa`, and ~175× faster than
`samtools faidx`. Every non-indexed reader pays an O(N) cost per call — the
flat-TOC 2bit rebuilds the whole 500k-entry TOC (173 ms ⁶), 4bit scans all 500k
interleaved record headers (245 ms), `twoBitToFa` loads its index into a hash
(89 ms), and `samtools faidx` loads its `.fai` (245 ms).

Three notes on the comparison:

- **Locally, the index *structure* doesn't matter — all three tie at ~1.4 ms.**
  The measured per-fetch is dominated by process spawn + open, not the lookup
  (which is microseconds either way). The pointer array's ~19 probes vs the B+
  tree's ~3 node reads is invisible against that fixed overhead. The structure
  only separates them **remotely**, where each access is a network round-trip —
  see "Serving over the web" above.
- **The `--bpt` variant duplicates the names to get 2be's remote speed *with*
  backward compatibility.** It appends a full B+ tree (names inline) after the
  standard TOC; legacy `twoBitToFa` ignores it. That costs **+4.3 MiB over
  `--index`** here (the duplicated names) — and, being short `seq######` names,
  this is a *best case*: realistic 15–25 byte names widen that penalty ~3× and
  raise its remote request count (see the web section's caveat).
- **The pointer array's size here is the *IUB* tables, not the index.**
  2bit+IUB+index is 71.5 MiB, but that bulk is the scattered-degenerate N-blocks
  (which 2be's merged stream avoids); the pointer array itself is only ~3.8 MiB.
  Dropping `--iub` brings the file back near plain 2bit.

The takeaway: for **local** exact-name lookup, all three index structures are a
wash — the flat 8-byte pointer array (backward compatible, smallest) is the
obvious pick. The structure only earns its keep **over the web**, where the B+
tree's fan-out cuts round-trips — and `--bpt` gets that remotely without giving
up twoBit compatibility, at the cost of the duplicated name bytes.

> ⁶ The flat-TOC 2bit reader is seek-based too (note ⁵), so it no longer reads the
> whole 59 MB file — but with no name index it must still load the entire 500k-entry
> TOC to resolve one name, and that O(N) load is the 173 ms (provable here: the
> *indexed* reader is the same reader minus this step and runs in 1.4 ms).
> `twoBitToFa` pays the **same** O(N) load — kent's `twoBitOpen` also reads the whole
> index into a name hash on open — so it too is in this slow regime at 89 ms, ~64×
> the indexed reader. It's ~1.9× faster than *our* flat reader only by constant
> factor: our `flat()` heap-allocates each name twice (the map key plus the `order`
> vec) and uses Rust's SipHash, where kent keeps each name once in a `localmem`
> arena with a lighter hash. The name index removes the load entirely; it does not
> merely out-implement it.

### A note on block-table overhead

twoBit stores `N` **and every degenerate base** (all "non-ACGT") as run-length
**blocks**. That is extremely compact for real genomes, where `N`s come in long
runs (assembly gaps). But scattered ambiguity yields one length-1 block each, at
8 bytes — millions of them at high rates. The generator models both via
`--n-runs` / `--iub-runs` (benchmark `N_RUNS` / `IUB_RUNS`):

| config | N-blocks / 100 Mbp seq | 2bit standard | `twoBitToFa` 20k extract |
|---|---|---|---|
| N + IUB scattered (`0 / 0`) | ~1.5 M | 2.945 b/base (105 MiB) | 60–100 s |
| N clustered, IUB scattered (`3 / 0`, default) | ~490 k | 2.315 b/base (83 MiB) | ~20 s |
| N + IUB clustered (`3 / 3`) | 6 | 2.000 b/base (75 MiB) | sub-second |

Key point: because standard twoBit lumps `N` and every IUB code into one block
table, the **IUB rate/layout drives the block count** as much as `N` does. Our
reader stays fast in all three thanks to binary-searched lookup; `twoBitToFa`
(linear per fragment) tracks the block count. For a realistic genome, isolated
IUB codes are sparse (hundreds, not ~1.5 M), so the default's ~490 k blocks is a
deliberately heavy case.

### Serving over the web (HTTP range reads / UDC)

Every reader can open an `http(s)://` URL and fetch only the bytes it needs via
HTTP range requests — seqformat's analogue of UCSC's UDC layer (a shared
`Source` in `src/source.rs`, a pooled `ureq` connection, an 8 KiB block cache).
Just pass a URL where you'd pass a path:

```sh
seqformat extract https://host/genome.2be chr1:0-1000       # only header+index+window cross the wire
seqformat extract https://host/genome.idx.2bit chr7 --http-stats   # print requests + bytes
```

This is where the **name-lookup index earns its keep remotely**, because over a
network the cost is *round-trips and bytes*, not comparisons. Per single-sequence
fetch on 500k × 300 bp (`bench/webseq.sh`, served from a real host; `ms/fetch`
varies with the network, but `requests` and `bytes` are architecture-determined
and stable):

| format | index structure | ms/fetch | HTTP requests | bytes/fetch |
|---|---|--:|--:|--:|
| 2bit standard | flat, unsorted TOC | 242 | 3.9 | 6.7 MiB |
| 2bit + index | sorted pointer array | 83 | 21.7 | 170 KiB |
| **2bit + B+ tree** (`--bpt`) | **B+ tree, TOC duplicate** | **18** | **6.6** | **52 KiB** |
| 2be | on-disk B+ tree (fan-out 256) | 20 | 7.0 | 49 KiB |
| 4bit | none (interleaved) | — ¹ | 10241 | 84 MiB |
| faidx (plain) | `.fai` sidecar | — ¹ | 3 | 14.6 MiB |
| faidx (BGZF) | `.fai` + block scan | — ¹ | 2431 | 34.5 MiB |

- **Standard 2bit** has no ordered in-file index, so resolving one name pulls the
  **entire TOC** (O(N) bytes) — 6.7 MiB every open, bandwidth-bound.
- **2bit + index** binary-searches its sorted array — O(log₂N), but the probes
  scatter across the TOC (a pointer then a distant name), so ~22 small reads.
- **`2bit + B+ tree`** appends a full B+ tree (names inline, a TOC duplicate);
  a lookup touches ~3 nodes → **it matches 2be remotely while staying fully
  twoBit backward compatible** (legacy `twoBitToFa` ignores the appended blob).
  The price is storage: the duplicated names (+4.3 MiB over `--index` here).
- **2be**'s B+ tree gives the same ~7-request shape, but is not twoBit compatible
  (different record encoding). `--bpt` is how you get 2be's remote lookup on a
  standard 2bit.
- ¹ **4bit and faidx are O(N) on open, by design** (shown as a single fetch — a
  timed loop would be dominated by these). 4bit has no index, so a lookup scans
  every interleaved record header (~84 MiB, ~the whole file). faidx first pulls
  its entire `.fai` sidecar — **14.6 MiB** at 500k seqs — which dominates; BGZF
  adds an O(blocks) header scan. `--http-stats` folds the `.fai` fetch into the
  total, so faidx isn't flattered by hiding its index load. A remote faidx needs
  its `.fai` served alongside the file.

**Caveat — these numbers assume short names.** The benchmark uses `seq0…seq499999`
(4–9 bytes). The B+ tree pads keys to the longest name and inlines them, so with
short names its duplication penalty (+4.3 MiB) and its node packing (one 8 KiB
block per level → ~7 requests) are a *best case*. Realistic 15–25 byte names
(`scaffold_000123`, accessions) roughly **triple** the `--bpt`-vs-`--index`
storage gap and can push tree nodes across two blocks (~7 → ~10–13 requests). The
ranking holds, but the exact tree numbers are name-length-sensitive.

```sh
bash bench/webseq.sh    # FORMATS="label=url ..." FETCHES=15 to customize
```

## Format reference

### twoBit (standard / long)

```
header (16 bytes)
  signature   u32  0x1A412743   (native byte order — also a byte-order mark)
  version     u32  0 = standard, 1 = long
  seqCount    u32
  reserved    u32  0
index (per sequence)
  nameSize    u8
  name        nameSize bytes
  offset      u32 (standard) | u64 (long)
record (at offset)
  dnaSize      u32
  nBlockCount  u32 ;  nBlockStarts u32[] ;  nBlockSizes u32[]
  maskCount    u32 ;  maskStarts   u32[] ;  maskSizes   u32[]
  reserved     u32  0
  packedDna    ceil(dnaSize/4) bytes   (T=00, C=01, A=10, G=11; first base in
                                        the most-significant bits)
```

`-long` widens **only the index offsets** to 64 bits (matching current UCSC
kent source); everything else is identical. Endianness is auto-detected from the
signature on read; files are always written little-endian.

### IUB extension (backward compatible)

The standard format stores every non-ACGT residue — `N` *and* all degenerate
codes — as an `N`-block, so `R/Y/S/...` are indistinguishable from `N` on
readback. The extension recovers them **without breaking legacy readers**.

Two changes, both invisible to old tools:

1. **Per-record IUB table**, appended *after* `packedDna`:

   ```
   iubCount   u32
   iubStarts  u32[iubCount]
   iubSizes   u32[iubCount]
   iubCodes   u8[iubCount]    (ASCII IUB letter per run)
   ```

   A legacy reader finds each record via the index `offset` and reads exactly
   `ceil(dnaSize/4)` packed bytes — it never reads past the DNA. The IUB table
   lives in the gap *after* the packed DNA and *before* the next record (the
   offsets we write already skip over it), so legacy readers never see it.

2. **An 8-byte EOF trailer** that advertises the extension:

   ```
   magic       u32  0x55324232   ("2 B 2 U")
   extVersion  u32  1
   ```

   Legacy readers never look at the end of the file, so this is invisible too.
   Our reader checks the last 8 bytes to decide whether to parse IUB tables.

The IUB table does **not** duplicate `N`-blocks: a degenerate run is *also*
recorded as a normal `N`-block (so old readers mask it as `N`), and the IUB
table simply overrides those positions with the precise code. Bare `N` runs are
omitted from the IUB table.

Net effect: an old `twoBit` reader opens an `--iub` file as an ordinary 2bit and
reports degenerate positions as `N`; this tool additionally recovers the exact
code. `--iub` composes with `--long`.

> **Why the trailer exists at all.** The magic and `extVersion` (and likewise the
> name-index footer's magic, below) are there only to make the extension
> *discoverable out-of-band*. The natural place to advertise a format extension
> is the header `version` field — but bumping it would break every existing
> reader, since the kent code (and ours) rejects unknown versions. So instead the
> signal is tucked at EOF, where legacy readers never look, which is exactly what
> keeps the file 100% readable as a plain 2bit. **If UCSC agreed to adopt the IUB
> extension and the name index under a new official `version` number, that header
> field would signal their presence directly and the EOF magic/`extVersion` would
> be unnecessary.** The trade-off is the usual one: a version bump makes the file
> opaque to pre-existing tools (they'd refuse it), whereas the EOF trailer is the
> price of backward compatibility *without* a coordinated spec change.

### Sorted-name index (backward compatible)

The stock twoBit index is a run of **variable-length** `nameSize/name/offset`
entries in file order, so finding a sequence by name means loading every entry
(UCSC builds a hash; we build a `HashMap`) — O(N) work on every open, which
dominates per-fetch latency on files with hundreds of thousands of short
sequences (see the many-short table above). Variable-length, unsorted entries
can't be binary-searched in place, so O(log N) lookup needs an auxiliary
fixed-width structure.

`--index` appends one, **after all records**, that adds no name duplication:

```
nameIndex
  ptr          u64[seqCount]   sorted by name; ptr[i] = absolute file offset of
                               the index entry for the i-th name in sorted order
footer (24 bytes, at EOF)
  indexOffset  u64   absolute offset of the ptr array
  seqCount     u64   number of pointers
  flags        u32   bit 0 = records carry IUB tables
  magic        u32   0x58444E49  ("INDX")
```

A lookup binary-searches the pointer array: each probe follows a pointer into the
*original* index entry, reads the name there to compare, and on a hit reads the
record offset that already sits after the name. The names live exactly once, in
the existing index; the only added bytes are **8 per sequence** plus the footer.

The 24-byte footer **replaces** the 8-byte IUB trailer (its `flags` advertise
whether IUB tables are present), and like that trailer it sits past the last
record — a legacy reader, which reaches records purely through index offsets,
never sees it. So `--iub --index` is the full "transitional" format: exact
degenerate codes *and* O(log N) lookup, while UCSC `twoBitToFa`/`twoBitInfo`
still read it as an ordinary 2bit. `--index` composes with `--long`, and may be
used without `--iub` (a plain indexed 2bit).

This is the backward-compatible counterpart to 2be's B+ tree: a flat pointer
array gives the same O(log N) exact-name lookup at lower redundancy (no padded
key copies), trading away the B+ tree's range scans and disk-paging locality —
neither of which exact-name fetch needs.

**`--bpt` — the B+ tree variant.** For the *remote* case (see "Serving over the
web"), the pointer array's disk locality *does* matter: each of its ~log₂N probes
is a separate scattered HTTP range read (a pointer, then a distant name). `--bpt`
instead appends a **full [`bptree`] blob** — the same on-disk B+ tree 2be uses,
with names stored inline — as a *complete duplicate* of the TOC, under a distinct
footer magic (`0x58545042`, "BPTX"). Because keys are co-located per node, a
lookup touches ~log₂₅₆N ≈ 3 nodes instead of ~19 scattered probes, cutting remote
requests ~22 → ~7 (it matches 2be). It stays fully twoBit compatible — the blob
sits past the last record, invisible to legacy readers — at the cost of the
duplicated name bytes (**~name-length per sequence**; +4.3 MiB on 500k short
names, proportionally more on realistic ones). Use `--bpt` instead of `--index`
when serving over HTTP; use `--index` (or nothing) for local-only files, where
the structure makes no measurable difference.

### 4-bit (BWA/BAM style)

The nibble alphabet `=ACMGRSVTWYHKDBN` (index 0..15, high nibble first) is the
htslib `seq_nt16` encoding. Raw `.pac`/htslib blobs have no framing, so for a
self-contained, round-trippable test file the nibbles are wrapped in a minimal
container:

```
header:  magic u32 0x54494234 ("4BIT") | version u32 1 | seqCount u32 | reserved u32 0
per seq: nameSize u8 | name | length u64 | packed ceil(length/2) bytes
```

4-bit packing carries no case information, so FASTA → 4-bit → FASTA is
case-insensitive (output is upper-case).

### samtools format (indexed / BGZF FASTA)

This is an *index over FASTA*, not a self-contained container — the bases stay as
FASTA text and small sidecar files add random access (and bgzip adds
compression):

- **`.fai`** — plain-text TSV, one line per sequence:
  `name<TAB>length<TAB>offset<TAB>linebases<TAB>linewidth`, where `offset` is the
  byte offset of the first base in the *uncompressed* FASTA, `linebases` is bases
  per line and `linewidth` is bytes per line incl. the newline. The byte offset
  of base `p` is `offset + (p/linebases)*linewidth + p%linebases`.
- **BGZF** (`.fa.gz`) — gzip split into ≤64 KiB blocks, each carrying a `BC`
  extra-field with its size, so any block can be found and inflated alone.
- **`.gzi`** — `u64 count`, then a `(compressed_offset, uncompressed_offset)`
  u64 pair for each block boundary after the first; maps uncompressed positions
  to compressed blocks.

`seqformat`'s `.fai` and `.gzi` are **byte-identical** to those `samtools faidx`
produces, and BGZF files round-trip through `samtools`/`bgzip` both ways. (The
reader doesn't actually need the `.gzi` — BGZF block headers are self-describing,
so it can scan them — but the writer emits a correct one for interop.)

### 2be (experimental twoBit redesign)

`2be` drops backward compatibility to fix twoBit's two scaling limits. Layout:

```text
HEADER (32 B)   magic "2BE1" | version u16 | flags u16 | seqCount u64
                tocOffset u64 | reserved u64
SEQUENCE RECORD (self-contained; sequence-local u32 coordinates)
  dnaSize    u32
  editCount  u32
  edits      editCount × 10 B:  pos u32, len u32, type u8, code u8
                                 type ∈ {N_RUN, IUB_POINT, IUB_RUN, MASK_RUN}
  runCount   u32
  runIndex   runCount × u32      indices (into edits) of the non-point entries
  packedDna  ceil(dnaSize/4)     2-bit, T=00 C=01 A=10 G=11
B+ TREE TOC     name → absolute record offset
```

Two ideas:

- **B+ tree TOC** ([`src/bptree.rs`], bigBed-style: fixed-width keys, fan-out
  256, blob-relative node offsets) → name lookup is O(log N) with no full index
  load, so it scales to millions of sequences.
- **One merged, position-sorted tagged-edit stream per sequence.** `N` and mask
  stay run-length (cheap for gaps); isolated IUB codes are `IUB_POINT` entries
  (~one each) and degenerate runs are `IUB_RUN`. Emitting a sequence is a single
  sweep of this stream over the 2-bit array. The small `runIndex` lets a
  sub-range fetch find straddling runs without scanning past the (possibly many)
  point entries; points in range are found by binary search.

It round-trips exactly (IUB + soft-mask), and region extraction matches
`samtools faidx` of the source. Its size sits just above standard 2bit (the edit
stream + tiny B+ tree); its real advantage — scalable name lookup — shows at high
sequence counts. Convert with `fa2be` / `be2fa`; `extract` and `info` auto-detect
it.

## Library

The same functionality is exposed as a library (`seqformat::twobit`,
`seqformat::fourbit`, `seqformat::fasta`) with in-memory `to_bytes` / `from_bytes`
entry points and on-disk `read_file` / `write_file` wrappers. The name-indexed
variant adds `to_bytes_indexed` / `write_file_indexed`; `TwoBitReader` reads it
transparently and uses the pointer array (`has_name_index()` reports its
presence).

`TwoBitReader::open(path)` is **seek-based** — it parses only the header and then
`seek`s to read each requested window (and, for an indexed file, a handful of
index-probe bytes), never loading the whole file, which is what gives the ~1 ms
single-fetch latency. `TwoBitReader::from_vec(bytes)` keeps an in-memory buffer
for callers that already hold the file (the `extract` CLI uses this for large
batches, where one whole-file read beats thousands of seeks). See the unit tests
in `src/twobit.rs` for round-trip examples, including the backward-compatibility,
binary-search, and file-vs-memory parity checks.
