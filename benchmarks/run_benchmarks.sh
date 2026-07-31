#!/bin/bash
# Phase 5 benchmark runner (plan.md) - the single reproducible command
# for both the primitive-level and end-to-end old-vs-new comparisons.
#
# Usage:
#   BUILD_DIR=/path/to/build ./run_benchmarks.sh
# BUILD_DIR defaults to ../build relative to this script (the repo's
# conventional build directory).
#
# Writes benchmarks/results_primitives.csv and
# benchmarks/results_end_to_end.csv, and prints a summary.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${BUILD_DIR:-$SCRIPT_DIR/../build}"
FASTPROT="$BUILD_DIR/src/c++/fastprot"
BENCH_PRIMITIVES="$BUILD_DIR/src/c++/bench_primitives"
DATA_DIR="$SCRIPT_DIR/data"
REPS="${REPS:-7}" # odd, so the median is an actual sample

if [ ! -x "$FASTPROT" ] || [ ! -x "$BENCH_PRIMITIVES" ]; then
    echo "Building fastprot and bench_primitives in $BUILD_DIR ..."
    cmake --build "$BUILD_DIR" --target fastprot --target bench_primitives -j"$(sysctl -n hw.ncpu 2>/dev/null || nproc)"
fi

mkdir -p "$DATA_DIR"

echo "=== Primitive-level micro-benchmark (count_id_dist, count_replacements) ==="
"$BENCH_PRIMITIVES" | tee "$SCRIPT_DIR/results_primitives.csv"

echo
echo "=== End-to-end fastprot benchmark (JCK/Kimura, old vs new, $REPS reps/case) ==="

median_of() {
    # reads whitespace-separated numbers from stdin, prints the median
    python3 -c "
import sys
xs = sorted(float(x) for x in sys.stdin.read().split())
n = len(xs)
print(xs[n // 2])
"
}

echo "n_seqs,length,old_median_s,new_median_s,speedup" > "$SCRIPT_DIR/results_end_to_end.csv"

# (n_seqs, length) - short/typical/long lengths x small/large dataset
# sizes, matching plan.md's Phase 5 requirement to vary both.
for case in "150 100" "150 500" "600 300" "1500 300"; do
    read -r n_seqs length <<< "$case"
    dataset="$DATA_DIR/bench_${n_seqs}x${length}.fasta"
    if [ ! -f "$dataset" ]; then
        python3 "$SCRIPT_DIR/gen_dataset.py" "$n_seqs" "$length" "$dataset"
    fi

    old_times=""
    new_times=""
    for i in $(seq 1 "$REPS"); do
        t0=$(python3 -c "import time;print(time.time())")
        "$FASTPROT" -I fasta "$dataset" -D JCK -O phylip -o /dev/null
        t1=$(python3 -c "import time;print(time.time())")
        old_times="$old_times $(python3 -c "print($t1-$t0)")"

        t0=$(python3 -c "import time;print(time.time())")
        FASTPROT_NEW_PROT_CODE=1 "$FASTPROT" -I fasta "$dataset" -D JCK -O phylip -o /dev/null
        t1=$(python3 -c "import time;print(time.time())")
        new_times="$new_times $(python3 -c "print($t1-$t0)")"
    done

    old_median=$(echo "$old_times" | median_of)
    new_median=$(echo "$new_times" | median_of)
    speedup=$(python3 -c "print($old_median/$new_median)")

    echo "n_seqs=$n_seqs length=$length: old=${old_median}s new=${new_median}s speedup=${speedup}x"
    echo "$n_seqs,$length,$old_median,$new_median,$speedup" >> "$SCRIPT_DIR/results_end_to_end.csv"
done

echo
echo "Results written to $SCRIPT_DIR/results_primitives.csv and $SCRIPT_DIR/results_end_to_end.csv"
