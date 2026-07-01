#!/usr/bin/env bash
#
# Benchmark the sequence container formats produced by `seqformat` against the
# UCSC tools (faToTwoBit / twoBitToFa) and samtools faidx (bgzip). seqformat also
# implements the samtools BGZF format and the experimental 2be format natively.
#
# Reports, for randomly generated data with ~1% N and ~0.5% IUB codes:
#   1) storage size of each format
#   2) encoding time (warmup + N timed iterations -> min/median/mean/sd)
#   3) random-access extraction time of many regions
#
# Config via env (defaults shown):
#   SEQS=3 LEN=100000000 NFRAC=0.01 IUBFRAC=0.005 SEED=1
#   ENC_ITERS=3 EXT_ITERS=3 NREGIONS=20000 REGION_LEN=500
#   WORKDIR=./bench/work  BIN=./target/release/seqformat

set -euo pipefail

SEQS=${SEQS:-3}
LEN=${LEN:-100000000}
NFRAC=${NFRAC:-0.01}
IUBFRAC=${IUBFRAC:-0.005}
# Cluster ambiguity into this many random-sized runs per sequence (assembly-gap
# style). 0 = scatter per-base (millions of length-1 blocks — the old default).
N_RUNS=${N_RUNS:-3}
IUB_RUNS=${IUB_RUNS:-0}
SEED=${SEED:-1}
ENC_ITERS=${ENC_ITERS:-3}
EXT_ITERS=${EXT_ITERS:-3}
# twoBitToFa is ~5 ms/region on this pathological (millions of scattered
# N-blocks) input, so its full extraction run is minutes; give it fewer iters.
SLOW_EXT_ITERS=${SLOW_EXT_ITERS:-1}
NREGIONS=${NREGIONS:-20000}
REGION_LEN=${REGION_LEN:-500}
FETCHES=${FETCHES:-300}   # single-region fetches in separate processes (per-fetch latency)
WORKDIR=${WORKDIR:-./bench/work}
BIN=${BIN:-./target/release/seqformat}

TOTAL_BP=$((SEQS * LEN))
FA="$WORKDIR/seq.fa"

have() { command -v "$1" >/dev/null 2>&1; }
hr() { printf '%s\n' "----------------------------------------------------------------------"; }
now() { date +%s.%N; }
mib() { awk -v b="$1" 'BEGIN{printf "%.2f", b/1048576}'; }
sizeof() { stat -c %s "$1" 2>/dev/null || echo 0; }
bpb() { awk -v s="$1" -v n="$TOTAL_BP" 'BEGIN{printf "%.3f", s*8/n}'; }

bench() {
  local label="$1" iters="$2"; shift 2
  if ! "$@" >/dev/null 2>&1; then printf '  %-28s FAILED\n' "$label"; return 0; fi
  local times=() t0 t1
  for ((i=0; i<iters; i++)); do
    t0=$(now); "$@" >/dev/null 2>&1; t1=$(now)
    times+=("$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.6f", b-a}')")
  done
  printf '%s\n' "${times[@]}" | sort -n | awk -v l="$label" '
    {v[NR]=$1; s+=$1}
    END{ n=NR; mean=s/n; min=v[1]; med=(n%2)?v[(n+1)/2]:(v[n/2]+v[n/2+1])/2;
         for(i=1;i<=n;i++){d=v[i]-mean; ss+=d*d} sd=(n>1)?sqrt(ss/(n-1)):0;
         printf "  %-28s min %7.3fs  median %7.3fs  mean %7.3fs  sd %6.3fs\n", l, min, med, mean, sd }'
}
seqonly_md5() { grep -v '^>' "$1" | tr -d '\n' | md5sum | cut -d' ' -f1; }

echo "seqformat benchmark"
echo "  data : $SEQS x $LEN bp = $TOTAL_BP bp;  N $NFRAC, IUB $IUBFRAC, seed $SEED"
echo "  ambig: N runs/seq=$N_RUNS, IUB runs/seq=$IUB_RUNS  (0 = scattered per-base)"
echo "  iters: encode $ENC_ITERS, extract $EXT_ITERS (slow extract $SLOW_EXT_ITERS)"
echo "  extr : $NREGIONS regions x $REGION_LEN bp;  bin $BIN"
for t in faToTwoBit twoBitToFa bgzip samtools; do
  have "$t" && echo "  tool : $t" || echo "  tool : $t MISSING (skipped)"
done

[ -x "$BIN" ] || { echo "building..."; cargo build --release >/dev/null 2>&1; }
mkdir -p "$WORKDIR"
hr; echo "Generating test data..."
"$BIN" random "$FA" --seqs "$SEQS" --length "$LEN" --n-frac "$NFRAC" \
  --iub-frac "$IUBFRAC" --n-runs "$N_RUNS" --iub-runs "$IUB_RUNS" --seed "$SEED" >/dev/null

# ---- encoders / extractors -------------------------------------------------
enc_ucsc_2bit()  { faToTwoBit "$FA" "$WORKDIR/ucsc.2bit"; }
enc_our_std()    { "$BIN" fa2twobit "$FA" "$WORKDIR/our.std.2bit"; }
enc_our_iub()    { "$BIN" fa2twobit "$FA" "$WORKDIR/our.iub.2bit" --iub; }
enc_our_idx()    { "$BIN" fa2twobit "$FA" "$WORKDIR/our.idx.2bit" --iub --index; }
enc_our_bpt()    { "$BIN" fa2twobit "$FA" "$WORKDIR/our.bpt.2bit" --iub --bpt; }
enc_our_4bit()   { "$BIN" fa2fourbit "$FA" "$WORKDIR/our.4bit"; }
enc_our_bgzf()   { "$BIN" fa2faidx "$FA" "$WORKDIR/our.fa.gz" --bgzip; }
enc_our_2be()    { "$BIN" fa2be "$FA" "$WORKDIR/our.2be"; }
enc_sam()        { rm -f "$WORKDIR/sam.fa.gz" "$WORKDIR/sam.fa.gz".{fai,gzi};
                   bgzip -c "$FA" > "$WORKDIR/sam.fa.gz"; samtools faidx "$WORKDIR/sam.fa.gz"; }

ext_ucsc()       { twoBitToFa -seqList="$WORKDIR/reg0.txt" "$WORKDIR/ucsc.2bit" "$WORKDIR/o.ucsc.fa"; }
ext_our_std()    { "$BIN" extract "$WORKDIR/our.std.2bit" --seq-list "$WORKDIR/reg0.txt" --out "$WORKDIR/o.std.fa"; }
ext_our_iub()    { "$BIN" extract "$WORKDIR/our.iub.2bit" --seq-list "$WORKDIR/reg0.txt" --out "$WORKDIR/o.iub.fa"; }
ext_our_idx()    { "$BIN" extract "$WORKDIR/our.idx.2bit" --seq-list "$WORKDIR/reg0.txt" --out "$WORKDIR/o.idx.fa"; }
ext_our_bpt()    { "$BIN" extract "$WORKDIR/our.bpt.2bit" --seq-list "$WORKDIR/reg0.txt" --out "$WORKDIR/o.bpt.fa"; }
ext_our_4bit()   { "$BIN" extract "$WORKDIR/our.4bit"     --seq-list "$WORKDIR/reg0.txt" --out "$WORKDIR/o.4bit.fa"; }
ext_our_bgzf()   { "$BIN" extract "$WORKDIR/our.fa.gz"    --seq-list "$WORKDIR/reg0.txt" --out "$WORKDIR/o.ourbgzf.fa"; }
ext_our_2be()    { "$BIN" extract "$WORKDIR/our.2be"      --seq-list "$WORKDIR/reg0.txt" --out "$WORKDIR/o.2be.fa"; }
ext_sam()        { samtools faidx -r "$WORKDIR/reg1.txt" "$WORKDIR/sam.fa.gz" -o "$WORKDIR/o.sam.fa"; }

echo "Building all formats once..."
have faToTwoBit && enc_ucsc_2bit
enc_our_std; enc_our_iub; enc_our_idx; enc_our_bpt; enc_our_4bit; enc_our_bgzf; enc_our_2be
have bgzip && have samtools && enc_sam || true

# ---- block counts (shows how clustering affects the 2bit block tables) ------
hr; echo "BLOCK COUNTS  (per-sequence N-blocks / mask / IUB, from the IUB 2bit)"
"$BIN" info "$WORKDIR/our.iub.2bit" | sed -n '3,$p'

# ---- (1) storage -----------------------------------------------------------
hr; echo "(1) STORAGE SIZE   (total sequence = $TOTAL_BP bp)"
printf '  %-30s %12s %10s %12s\n' "format" "bytes" "MiB" "bits/base"
rs() { printf '  %-30s %12s %10s %12s\n' "$1" "$2" "$(mib "$2")" "$(bpb "$2")"; }
rs "FASTA (raw, uppercase)" "$(sizeof "$FA")"
have faToTwoBit && rs "2bit (UCSC kentsrc)" "$(sizeof "$WORKDIR/ucsc.2bit")"
if [ -f "$WORKDIR/sam.fa.gz" ]; then
  b=$(sizeof "$WORKDIR/sam.fa.gz")
  i=$(( $(sizeof "$WORKDIR/sam.fa.gz.fai") + $(sizeof "$WORKDIR/sam.fa.gz.gzi") ))
  rs "BGZF (samtools)" "$b"
  printf '  %-30s %12s %10s\n' "  + .fai/.gzi" "$i" "$(mib "$i")"
fi
if [ -f "$WORKDIR/our.fa.gz" ]; then
  b=$(sizeof "$WORKDIR/our.fa.gz")
  i=$(( $(sizeof "$WORKDIR/our.fa.gz.fai") + $(sizeof "$WORKDIR/our.fa.gz.gzi") ))
  rs "BGZF (seqformat)" "$b"
  printf '  %-30s %12s %10s\n' "  + .fai/.gzi" "$i" "$(mib "$i")"
fi
rs "2bit (seqformat)" "$(sizeof "$WORKDIR/our.std.2bit")"
rs "2bit+IUB (seqformat)" "$(sizeof "$WORKDIR/our.iub.2bit")"
rs "2bit+IUB+index (seqformat)" "$(sizeof "$WORKDIR/our.idx.2bit")"
rs "2bit+IUB+bptree (seqformat)" "$(sizeof "$WORKDIR/our.bpt.2bit")"
rs "4bit (seqformat)" "$(sizeof "$WORKDIR/our.4bit")"
rs "2be (seqformat)" "$(sizeof "$WORKDIR/our.2be")"

# ---- regions + interop -----------------------------------------------------
hr; echo "INTEROP / CORRECTNESS CHECK"
awk -v k="$SEQS" -v L="$LEN" -v w="$REGION_LEN" -v n="$NREGIONS" -v seed="$SEED" \
    -v f0="$WORKDIR/reg0.txt" -v f1="$WORKDIR/reg1.txt" '
  BEGIN{ srand(seed);
    for(i=0;i<n;i++){ s=int(rand()*k); st=int(rand()*(L-w)); en=st+w;
      printf "seq%d:%d-%d\n", s, st,   en > f0;     # 0-based half-open
      printf "seq%d:%d-%d\n", s, st+1, en > f1; } }' # 1-based inclusive

ext_our_std; ext_our_iub; ext_our_idx; ext_our_4bit; ext_our_bgzf; ext_our_2be
have twoBitToFa && ext_ucsc || true
have samtools && ext_sam || true
check() { [ "$2" = "$3" ] && echo "  PASS  $1" || echo "  FAIL  $1"; }
ourstd=$(seqonly_md5 "$WORKDIR/o.std.fa")
if have twoBitToFa; then
  check "seqformat 2bit == UCSC kentsrc (degenerate=N)" "$ourstd" "$(seqonly_md5 "$WORKDIR/o.ucsc.fa")"
  # The indexed file must still read as a plain 2bit through the UCSC tool.
  twoBitToFa -seqList="$WORKDIR/reg0.txt" "$WORKDIR/our.idx.2bit" "$WORKDIR/o.idx.ucsc.fa" 2>/dev/null \
    && check "UCSC twoBitToFa reads 2bit+IUB+index (degenerate=N)" "$ourstd" "$(seqonly_md5 "$WORKDIR/o.idx.ucsc.fa")"
  twoBitToFa -seqList="$WORKDIR/reg0.txt" "$WORKDIR/our.bpt.2bit" "$WORKDIR/o.bpt.ucsc.fa" 2>/dev/null \
    && check "UCSC twoBitToFa reads 2bit+IUB+bptree (degenerate=N)" "$ourstd" "$(seqonly_md5 "$WORKDIR/o.bpt.ucsc.fa")"
fi
if have samtools; then
  sam=$(seqonly_md5 "$WORKDIR/o.sam.fa")
  check "seqformat 2bit+IUB == samtools faidx " "$(seqonly_md5 "$WORKDIR/o.iub.fa")" "$sam"
  check "seqformat 2bit+IUB+index== samtools faidx" "$(seqonly_md5 "$WORKDIR/o.idx.fa")" "$sam"
  check "seqformat 4bit      == samtools faidx" "$(seqonly_md5 "$WORKDIR/o.4bit.fa")" "$sam"
  check "seqformat BGZF      == samtools faidx" "$(seqonly_md5 "$WORKDIR/o.ourbgzf.fa")" "$sam"
  check "seqformat 2be       == samtools faidx" "$(seqonly_md5 "$WORKDIR/o.2be.fa")" "$sam"
  # byte-identical index check: let samtools re-index OUR bgzf file
  cp "$WORKDIR/our.fa.gz" "$WORKDIR/reidx.fa.gz"; samtools faidx "$WORKDIR/reidx.fa.gz"
  cmp -s "$WORKDIR/our.fa.gz.fai" "$WORKDIR/reidx.fa.gz.fai" && echo "  PASS  our .fai byte-identical to samtools" || echo "  FAIL  .fai"
  cmp -s "$WORKDIR/our.fa.gz.gzi" "$WORKDIR/reidx.fa.gz.gzi" && echo "  PASS  our .gzi byte-identical to samtools" || echo "  FAIL  .gzi"
fi

# ---- (2) encoding ----------------------------------------------------------
hr; echo "(2) ENCODING TIME   (FASTA -> format)"
have faToTwoBit && bench "2bit (UCSC kentsrc)" "$ENC_ITERS" enc_ucsc_2bit
have bgzip && have samtools && bench "BGZF (samtools)" "$ENC_ITERS" enc_sam
bench "BGZF (seqformat)"       "$ENC_ITERS" enc_our_bgzf
bench "2bit (seqformat)"       "$ENC_ITERS" enc_our_std
bench "2bit+IUB (seqformat)"   "$ENC_ITERS" enc_our_iub
bench "2bit+IUB+index (seqformat)" "$ENC_ITERS" enc_our_idx
bench "2bit+IUB+bptree (seqformat)" "$ENC_ITERS" enc_our_bpt
bench "4bit (seqformat)"       "$ENC_ITERS" enc_our_4bit
bench "2be (seqformat)"        "$ENC_ITERS" enc_our_2be

# ---- (3) bulk extraction (one process, open cost amortized) ----------------
hr; echo "(3) BULK EXTRACT   ($NREGIONS regions x $REGION_LEN bp, ONE process)"
have twoBitToFa && bench "2bit (UCSC kentsrc)" "$SLOW_EXT_ITERS" ext_ucsc
have samtools && bench "BGZF (samtools)" "$EXT_ITERS" ext_sam
bench "BGZF (seqformat)"      "$EXT_ITERS" ext_our_bgzf
bench "2bit (seqformat)"      "$EXT_ITERS" ext_our_std
bench "2bit+IUB (seqformat)"  "$EXT_ITERS" ext_our_iub
bench "2bit+IUB+index (seqformat)" "$EXT_ITERS" ext_our_idx
bench "2bit+IUB+bptree (seqformat)" "$EXT_ITERS" ext_our_bpt
bench "4bit (seqformat)"      "$EXT_ITERS" ext_our_4bit
bench "2be (seqformat)"       "$EXT_ITERS" ext_our_2be

# ---- (4) per-fetch latency (one region per separate process) ---------------
hr; echo "(4) PER-FETCH LATENCY   ($FETCHES single-region fetches, SEPARATE processes)"
mapfile -t R0 < <(head -n "$FETCHES" "$WORKDIR/reg0.txt")   # 0-based half-open
mapfile -t R1 < <(head -n "$FETCHES" "$WORKDIR/reg1.txt")   # 1-based inclusive
lat() {  # lat LABEL fetchfn   (fetchfn takes a region index)
  local label="$1" fn="$2" t0 t1; t0=$(now)
  for ((i=0; i<${#R0[@]}; i++)); do "$fn" "$i" >/dev/null 2>&1 || true; done
  t1=$(now)
  awk -v l="$label" -v a="$t0" -v b="$t1" -v f="${#R0[@]}" \
    'BEGIN{printf "  %-28s total %7.3fs  (%.3f ms/fetch)\n", l, b-a, 1000*(b-a)/f}'
}
pf_ucsc()  { local s="${R0[$1]}" nm rng a b; nm="${s%%:*}"; rng="${s#*:}"; a="${rng%-*}"; b="${rng#*-}";
             twoBitToFa -seq="$nm" -start="$a" -end="$b" "$WORKDIR/ucsc.2bit" /dev/stdout >/dev/null; }
pf_sam()   { samtools faidx "$WORKDIR/sam.fa.gz" "${R1[$1]}" >/dev/null; }
pf_obgzf() { "$BIN" extract "$WORKDIR/our.fa.gz"    "${R0[$1]}" >/dev/null; }
pf_std()   { "$BIN" extract "$WORKDIR/our.std.2bit" "${R0[$1]}" >/dev/null; }
pf_iub()   { "$BIN" extract "$WORKDIR/our.iub.2bit" "${R0[$1]}" >/dev/null; }
pf_idx()   { "$BIN" extract "$WORKDIR/our.idx.2bit" "${R0[$1]}" >/dev/null; }
pf_bpt()   { "$BIN" extract "$WORKDIR/our.bpt.2bit" "${R0[$1]}" >/dev/null; }
pf_4bit()  { "$BIN" extract "$WORKDIR/our.4bit"     "${R0[$1]}" >/dev/null; }
pf_2be()   { "$BIN" extract "$WORKDIR/our.2be"      "${R0[$1]}" >/dev/null; }
have twoBitToFa && lat "2bit (UCSC kentsrc)" pf_ucsc
have samtools && lat "BGZF (samtools)" pf_sam
[ -f "$WORKDIR/our.fa.gz" ] && lat "BGZF (seqformat)" pf_obgzf
lat "2bit (seqformat)" pf_std
lat "2bit+IUB (seqformat)" pf_iub
lat "2bit+IUB+index (seqformat)" pf_idx
lat "2bit+IUB+bptree (seqformat)" pf_bpt
lat "4bit (seqformat)" pf_4bit
lat "2be (seqformat)" pf_2be

hr; echo "done. artifacts in $WORKDIR"
