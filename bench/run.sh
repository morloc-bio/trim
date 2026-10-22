#!/usr/bin/env bash
# Time, peak memory and storage of trimming simulated read pairs: trim on
# morloc stream files and on FASTQ text, the C++ trimmer bioproject, and
# Trim Galore 2.3.0 (the version nf-core/rnaseq 3.27.0 pins), each with Trim
# Galore's defaults, on every core.
#
#   GENOME_SIZE   genome size in bases (default 100000000)
#   COVERAGES     read coverages (default "3 9"; 3 of 100 Mb is 1M pairs)
#
# Columns: tool, step, pairs, seconds, peak_rss_mib, peak_disk_mib,
# written_mib, kept_pairs.
set -euo pipefail
cd "$(dirname "$0")"
L=../..
P=$L/bioprojects
W=work
RES=results.tsv
HIST=history.tsv
MEASURE="python3 $P/tools/measure.py"
export MEASURE_DIRS="$PWD/$W/out"
SIZE=${GENOME_SIZE:-100000000}
for x in "$P/synth/synth" "$P/trimmer/trimmer"; do
  [ -x "$x" ] || { echo "bench: $x missing" >&2; exit 1; }
done
[ -x ./trim-bench ] || morloc make -o trim-bench bench.loc
command -v trim_galore >/dev/null || { echo "bench: trim_galore is not on the PATH (bench/package.yaml)" >&2; exit 1; }
mkdir -p "$W/out"
$P/synth/synth -f packet -z 3 -o "$W/g.mlc" genome -n 200 -s "$SIZE" --seed 3
row() { local IFS=$'\t'; echo "$*" >> "$RES"; }
pairs_in() { echo $(( $(wc -l < "$1") / 4 )); }
printf 'tool\tstep\tpairs\tseconds\tpeak_rss_mib\tpeak_disk_mib\twritten_mib\tkept_pairs\n' > "$RES"
for c in ${COVERAGES:-3 9}; do
  $P/trimmer/trimmer simulate "$W/g.mlc" -c "$c" --r1="$W/r1.fq" --r2="$W/r2.fq" --no-stdout
  n=$(pairs_in "$W/r1.fq")
  rm -rf "$W/out"/*
  row trim ingest "$n" "$($MEASURE ./trim-bench ingestFastqPairs --compression-level 3 "$W/r1.fq" "$W/r2.fq" "$W/out/r.mlc")" NA
  mv "$W/out/r.mlc" "$W/r.mlc"
  row trim trim-morloc "$n" "$($MEASURE ./trim-bench trimMlc "$W/r.mlc" "$W/out/t.mlc")" NA
  rm -rf "$W/out"/*
  row trim trim-fastq "$n" "$($MEASURE ./trim-bench trimText "$W/r1.fq" "$W/r2.fq" "$W/out/t1.fq" "$W/out/t2.fq")" "$(pairs_in "$W/out/t1.fq")"
  rm -rf "$W/out"/*
  row trimmer-cpp trim-fastq "$n" "$($MEASURE $P/trimmer/trimmer trimText "$W/r1.fq" "$W/r2.fq" -O 1 --r1="$W/out/t1.fq" --r2="$W/out/t2.fq" --no-stdout)" "$(pairs_in "$W/out/t1.fq")"
  rm -rf "$W/out"/*
  row trim_galore-2.3.0 trim-fastq "$n" "$($MEASURE trim_galore --paired --cores "$(getconf _NPROCESSORS_ONLN)" -o "$W/out" "$W/r1.fq" "$W/r2.fq")" "$(pairs_in "$W/out/r1_val_1.fq")"
  rm -rf "$W/out"/* "$W/r.mlc" "$W"/*.fq
done
commit=$(git -C /work/morloc-workspace/compiler/morloc rev-parse --short HEAD 2>/dev/null || echo NA)
[ -f "$HIST" ] || printf 'date\tcompiler\tlabel\t%s\n' "$(head -n 1 "$RES")" > "$HIST"
tail -n +2 "$RES" | while IFS= read -r line; do
  printf '%s\t%s\t%s\t%s\n' "$(date +%F)" "$commit" "${BENCH_LABEL:-run}" "$line" >> "$HIST"
done
rm -rf "$W"
cat "$RES"
