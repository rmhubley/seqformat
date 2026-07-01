#!/usr/bin/env bash
#
# Web (HTTP-range) benchmark — the remote analogue of manyseq.sh's per-fetch
# column. Fetches single sequences from *remote* files over HTTP range requests
# (seqformat's UDC-style transport) and reports, per format:
#   - wall-clock latency per fetch (separate process each, like manyseq.sh)
#   - HTTP requests and bytes transferred per fetch (from --http-stats)
#
# It shows how the name-lookup index maps to network cost:
#   2bit std   flat unsorted TOC   -> pull whole TOC (O(N) bytes) per open
#   2bit idx   sorted name index   -> binary search  (~log2 N scattered reads)
#   2be        on-disk B+ tree     -> ~log_256 N node reads (~3), the ideal
#   4bit       no index            -> O(N): scan every record header
#   faidx      .fai sidecar        -> load whole .fai (O(N)); plain window read
#
# Config via env:
#   FORMATS="label=url label=url ..."   (default: the four 500k core formats)
#   FETCHES=15                          single-seq fetches per format
#   COUNT=500000                        sequence count (for picking names)
#   BIN=./target/release/seqformat
#   SEQFORMAT_HTTP_BLOCK=8192           cache/fetch block size (UDC default)

set -euo pipefail

BASE=${BASE:-https://repeatmasker.org/~rhubley}
FORMATS=${FORMATS:-"2bit-std=$BASE/seqformat-std.2bit 2bit-idx=$BASE/seqformat-idx.2bit 2be=$BASE/seqformat.2be"}
FETCHES=${FETCHES:-15}
COUNT=${COUNT:-500000}
BIN=${BIN:-./target/release/seqformat}

now() { date +%s.%N; }
[ -x "$BIN" ] || { echo "building..."; cargo build --release >/dev/null 2>&1; }

names=(); step=$((COUNT/FETCHES)); [ "$step" -lt 1 ] && step=1
for ((i=0;i<FETCHES;i++)); do names+=("seq$((i*step))"); done

echo "web per-fetch benchmark ($FETCHES fetches, separate processes, block=${SEQFORMAT_HTTP_BLOCK:-8192})"
echo "----------------------------------------------------------------------"
printf '  %-12s %13s %14s %14s\n' "format" "ms/fetch" "req/fetch" "KiB/fetch"

run() {  # run LABEL URL
  local label="$1" url="$2" t0 t1 reqs=0 bytes=0 line b
  t0=$(now)
  for nm in "${names[@]}"; do
    line=$("$BIN" extract "$url" "$nm" --http-stats 2>&1 >/dev/null | sed -n 's/^http: //p')
    reqs=$((reqs + ${line%% requests*}))
    b=${line#*requests, }; bytes=$((bytes + ${b% bytes}))
  done
  t1=$(now)
  awk -v l="$label" -v a="$t0" -v b="$t1" -v f="$FETCHES" -v r="$reqs" -v by="$bytes" 'BEGIN{
    printf "  %-12s %10.1f    %10.1f     %10.1f\n", l, 1000*(b-a)/f, r/f, by/f/1024 }'
}

for entry in $FORMATS; do
  run "${entry%%=*}" "${entry#*=}"
done
echo "----------------------------------------------------------------------"
echo "done."
