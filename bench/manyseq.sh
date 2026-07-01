#!/usr/bin/env bash
#
# Many-short-sequences benchmark — stresses the *sequence-count* axis (the one
# the 2be B+ tree addresses) rather than per-sequence size. Generates a large
# number of short sequences and measures, across formats:
#   1) build time  (FASTA -> format)
#   2) storage size
#   3) bulk extraction  (fetch a big sample of whole sequences in ONE process)
#   4) per-fetch latency (fetch single sequences in SEPARATE processes — each
#      pays the index-open cost, so flat TOCs pay O(N) every call, the 2be B+
#      tree pays O(log N))
#
# Config via env (defaults shown):
#   COUNT=500000 LEN=300 NFRAC=0.01 IUBFRAC=0.005 SEED=1
#   BULK=20000   FETCHES=300   ENC_ITERS=2
#   WORKDIR=./bench/manywork  BIN=./target/release/seqformat

set -euo pipefail

COUNT=${COUNT:-500000}
LEN=${LEN:-300}
NFRAC=${NFRAC:-0.01}
IUBFRAC=${IUBFRAC:-0.005}
SEED=${SEED:-1}
BULK=${BULK:-20000}
FETCHES=${FETCHES:-300}
ENC_ITERS=${ENC_ITERS:-2}
WORKDIR=${WORKDIR:-./bench/manywork}
BIN=${BIN:-./target/release/seqformat}

FA="$WORKDIR/many.fa"
have() { command -v "$1" >/dev/null 2>&1; }
hr() { printf '%s\n' "----------------------------------------------------------------------"; }
now() { date +%s.%N; }
mib() { awk -v b="$1" 'BEGIN{printf "%.2f", b/1048576}'; }
bpb() { awk -v s="$1" -v n="$((COUNT*LEN))" 'BEGIN{printf "%.3f", s*8/n}'; }
sizeof() { stat -c %s "$1" 2>/dev/null || echo 0; }

# bench LABEL ITERS -- cmd...
bench() {
  local label="$1" iters="$2"; shift 2
  if ! "$@" >/dev/null 2>&1; then printf '  %-26s FAILED\n' "$label"; return 0; fi
  local times=() t0 t1
  for ((i=0;i<iters;i++)); do t0=$(now); "$@" >/dev/null 2>&1; t1=$(now)
    times+=("$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.6f",b-a}')"); done
  printf '%s\n' "${times[@]}" | sort -n | awk -v l="$label" '
    {v[NR]=$1} END{n=NR; med=(n%2)?v[(n+1)/2]:(v[n/2]+v[n/2+1])/2;
      printf "  %-26s min %8.3fs  median %8.3fs\n", l, v[1], med}'
}

echo "many-short-sequences benchmark"
echo "  $COUNT sequences x $LEN bp = $((COUNT*LEN)) bp;  N $NFRAC, IUB $IUBFRAC scattered"
echo "  bulk sample $BULK seqs (1 process);  latency $FETCHES fetches (separate processes)"
for t in faToTwoBit twoBitToFa bgzip samtools; do
  have "$t" && echo "  tool : $t" || echo "  tool : $t MISSING (skipped)"; done
[ -x "$BIN" ] || { echo "building..."; cargo build --release >/dev/null 2>&1; }
mkdir -p "$WORKDIR"

hr; echo "Generating $COUNT short sequences..."
"$BIN" random "$FA" --seqs "$COUNT" --length "$LEN" --n-frac "$NFRAC" \
  --iub-frac "$IUBFRAC" --n-runs 1 --seed "$SEED" >/dev/null

# ---- format builders -------------------------------------------------------
enc_ucsc()  { faToTwoBit "$FA" "$WORKDIR/m.2bit"; }
enc_std()   { "$BIN" fa2twobit "$FA" "$WORKDIR/m.std.2bit"; }
# 2bit + IUB + sorted-name index: backward-compatible, O(log N) lookup.
enc_idx()   { "$BIN" fa2twobit "$FA" "$WORKDIR/m.idx.2bit" --iub --index; }
# 2bit + IUB + full B+ tree name index (TOC duplicate); backward compatible.
enc_bpt()   { "$BIN" fa2twobit "$FA" "$WORKDIR/m.bpt.2bit" --iub --bpt; }
enc_4bit()  { "$BIN" fa2fourbit "$FA" "$WORKDIR/m.4bit"; }
enc_2be()   { "$BIN" fa2be "$FA" "$WORKDIR/m.2be"; }
enc_sam()   { rm -f "$WORKDIR/m.fa.gz" "$WORKDIR/m.fa.gz".{fai,gzi};
              bgzip -c "$FA" > "$WORKDIR/m.fa.gz"; samtools faidx "$WORKDIR/m.fa.gz"; }

echo "Building all formats once..."
have faToTwoBit && enc_ucsc; enc_std; enc_idx; enc_bpt; enc_4bit; enc_2be
have bgzip && have samtools && enc_sam || true

# ---- region lists (whole short sequences) ----------------------------------
# reg0 = 0-based half-open (twoBitToFa / seqformat); reg1 = 1-based (samtools).
awk -v c="$COUNT" -v L="$LEN" -v b="$BULK" -v seed="$SEED" \
    -v f0="$WORKDIR/reg0.txt" -v f1="$WORKDIR/reg1.txt" '
  BEGIN{ srand(seed); step=(c>b)?int(c/b):1;
    for(i=0;i<c && n<b;i+=step){ printf "seq%d:%d-%d\n", i, 0, L > f0;
                                 printf "seq%d:%d-%d\n", i, 1, L > f1; n++ } }'

hr; echo "(1) BUILD TIME   (FASTA -> format)"
have faToTwoBit && bench "2bit (UCSC kentsrc)" "$ENC_ITERS" enc_ucsc
bench "2bit (seqformat)" "$ENC_ITERS" enc_std
bench "2bit+IUB+index (seqformat)" "$ENC_ITERS" enc_idx
bench "2bit+IUB+bptree (seqformat)" "$ENC_ITERS" enc_bpt
bench "4bit (seqformat)" "$ENC_ITERS" enc_4bit
bench "2be (seqformat)"  "$ENC_ITERS" enc_2be
have bgzip && have samtools && bench "BGZF (samtools)" "$ENC_ITERS" enc_sam

hr; echo "(2) STORAGE SIZE"
printf '  %-30s %12s %10s %12s\n' "format" "bytes" "MiB" "bits/base"
rs() { printf '  %-30s %12s %10s %12s\n' "$1" "$2" "$(mib "$2")" "$(bpb "$2")"; }
rs "FASTA (raw)" "$(sizeof "$FA")"
have faToTwoBit && rs "2bit (UCSC kentsrc)" "$(sizeof "$WORKDIR/m.2bit")"
[ -f "$WORKDIR/m.fa.gz" ] && rs "BGZF (samtools)" "$(sizeof "$WORKDIR/m.fa.gz")"
rs "2bit (seqformat)" "$(sizeof "$WORKDIR/m.std.2bit")"
rs "2bit+IUB+index (seqformat)" "$(sizeof "$WORKDIR/m.idx.2bit")"
rs "2bit+IUB+bptree (seqformat)" "$(sizeof "$WORKDIR/m.bpt.2bit")"
rs "4bit (seqformat)" "$(sizeof "$WORKDIR/m.4bit")"
rs "2be (seqformat)" "$(sizeof "$WORKDIR/m.2be")"

hr; echo "(3) BULK EXTRACT   ($BULK whole seqs, ONE process — open cost amortized)"
have twoBitToFa && bench "2bit (UCSC kentsrc)" 2 \
  bash -c "twoBitToFa -seqList='$WORKDIR/reg0.txt' '$WORKDIR/m.2bit' '$WORKDIR/o.ucsc.fa'"
bench "2bit (seqformat)" 2 "$BIN" extract "$WORKDIR/m.std.2bit" --seq-list "$WORKDIR/reg0.txt" --out "$WORKDIR/o.std.fa"
bench "2bit+IUB+index (seqformat)" 2 "$BIN" extract "$WORKDIR/m.idx.2bit" --seq-list "$WORKDIR/reg0.txt" --out "$WORKDIR/o.idx.fa"
bench "2bit+IUB+bptree (seqformat)" 2 "$BIN" extract "$WORKDIR/m.bpt.2bit" --seq-list "$WORKDIR/reg0.txt" --out "$WORKDIR/o.bpt.fa"
bench "4bit (seqformat)" 2 "$BIN" extract "$WORKDIR/m.4bit" --seq-list "$WORKDIR/reg0.txt" --out "$WORKDIR/o.4bit.fa"
bench "2be (seqformat)"  2 "$BIN" extract "$WORKDIR/m.2be" --seq-list "$WORKDIR/reg0.txt" --out "$WORKDIR/o.2be.fa"
have samtools && bench "BGZF (samtools)" 2 \
  bash -c "samtools faidx -r '$WORKDIR/reg1.txt' '$WORKDIR/m.fa.gz' -o '$WORKDIR/o.sam.fa'"

# interop on the bulk output (sequence-only md5; coords differ but bases match)
so() { grep -v '^>' "$1" | tr -d '\n' | md5sum | cut -d' ' -f1; }
if have samtools; then
  [ "$(so "$WORKDIR/o.2be.fa")" = "$(so "$WORKDIR/o.sam.fa")" ] && echo "  interop: 2be == samtools  OK" || echo "  interop: 2be != samtools"
fi

hr; echo "(4) PER-FETCH LATENCY   ($FETCHES single-seq fetches, SEPARATE processes)"
# Sample names spread across the file.
names=(); step=$((COUNT/FETCHES)); [ "$step" -lt 1 ] && step=1
for ((i=0;i<FETCHES;i++)); do names+=("seq$((i*step))"); done

lat() {  # lat LABEL cmd_template (use {N} for name, {S}/{E} for 0-based, {S1} for 1-based start)
  local label="$1"; shift
  local t0 t1; t0=$(now)
  for nm in "${names[@]}"; do "$@" "$nm" >/dev/null 2>&1 || true; done
  t1=$(now)
  awk -v l="$label" -v a="$t0" -v b="$t1" -v f="$FETCHES" \
    'BEGIN{printf "  %-26s total %7.3fs  (%.3f ms/fetch)\n", l, b-a, 1000*(b-a)/f}'
}
fetch_std()  { "$BIN" extract "$WORKDIR/m.std.2bit" "$1:0-$LEN" >/dev/null; }
fetch_idx()  { "$BIN" extract "$WORKDIR/m.idx.2bit" "$1:0-$LEN" >/dev/null; }
fetch_bpt()  { "$BIN" extract "$WORKDIR/m.bpt.2bit" "$1:0-$LEN" >/dev/null; }
fetch_4bit() { "$BIN" extract "$WORKDIR/m.4bit" "$1:0-$LEN" >/dev/null; }
fetch_2be()  { "$BIN" extract "$WORKDIR/m.2be" "$1:0-$LEN" >/dev/null; }
fetch_ucsc() { twoBitToFa -seq="$1" "$WORKDIR/m.2bit" /dev/stdout >/dev/null; }
fetch_sam()  { samtools faidx "$WORKDIR/m.fa.gz" "$1:1-$LEN" >/dev/null; }

have twoBitToFa && lat "2bit (UCSC kentsrc)" fetch_ucsc
lat "2bit (seqformat, flat TOC)" fetch_std
lat "2bit+IUB+idx (seqformat, ptr array)" fetch_idx
lat "2bit+IUB+bpt (seqformat, B+ tree)" fetch_bpt
lat "4bit (seqformat, scan TOC)" fetch_4bit
lat "2be (seqformat, B+ tree)"  fetch_2be
have samtools && lat "BGZF faidx (samtools)" fetch_sam

hr; echo "done. artifacts in $WORKDIR"
