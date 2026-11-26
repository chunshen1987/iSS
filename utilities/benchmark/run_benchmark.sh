#!/usr/bin/env bash
# Time iSS executables on the benchmark surface, interleaving them so that
# drifting CPU clocks (laptops) affect all of them alike.
#
#   utilities/benchmark/run_benchmark.sh WORK_DIR OUT_DIR EXE [EXE ...]
#
# Run from the iSS repo root (the table paths are relative). WORK_DIR holds
# surface.dat and music_input (make_benchmark_surface.py). For each
# executable, bulk_deltaf_kind and repeat this writes one line to
# OUT_DIR/timings.txt: wall and CPU time, peak memory and the md5 of the
# samples, and keeps the summary OUT_DIR/<label>_df<kind>_<repeat>.json
# (sample_stats.py) instead of the samples.
#
# Environment: KINDS (default "20 21"), REPEATS (3), EVENTS (1000),
# SEED (1), OMP_NUM_THREADS (passed through).
#
# The executables write particle_samples.bin into the current directory.
# Use executables built with -Dlink_with_lib=OFF when comparing versions:
# the default build loads libiSS.so from the repo's lib/, so a copied
# iSS.e runs whatever library was installed last.
set -u
if [ $# -lt 3 ]; then
    sed -n '2,20p' "$0"
    exit 1
fi
work_dir=$1
out_dir=$2
shift 2
kinds=${KINDS:-"20 21"}
repeats=${REPEATS:-3}
events=${EVENTS:-1000}
seed=${SEED:-1}
here=$(cd "$(dirname "$0")" && pwd)

if [ ! -f iSS_parameters.dat ] || [ ! -d iSS_tables ]; then
    echo "Run from the iSS repo root." >&2
    exit 1
fi
if [ -e particle_samples.bin ]; then
    echo "particle_samples.bin exists in $(pwd); move it away first." >&2
    exit 1
fi
mkdir -p "$out_dir"
timings=$out_dir/timings.txt

for rep in $(seq 1 "$repeats"); do
    for kind in $kinds; do
        for exe in "$@"; do
            label=$(basename "$(dirname "$(dirname "$(dirname "$exe")")")")
            [ "$label" = "." ] && label=$(basename "$exe")
            log=$out_dir/${label}_df${kind}_${rep}.log
            /usr/bin/time -f "%e %U %M" -o "$out_dir/time.tmp" \
                "$exe" iSS_parameters.dat "$work_dir" surface.dat \
                bulk_deltaf_kind="$kind" number_of_repeated_sampling="$events" \
                randomSeed="$seed" use_binary_format=1 > "$log" 2>&1 || {
                    echo "$exe failed, see $log" >&2
                    exit 1
                }
            md5=$(md5sum particle_samples.bin | cut -c1-12)
            python3 "$here/sample_stats.py" summarize particle_samples.bin \
                "$out_dir/${label}_df${kind}_${rep}.json" > /dev/null
            rm -f particle_samples.bin
            read -r wall cpu rss < "$out_dir/time.tmp"
            printf "%-12s df=%s rep=%s threads=%s wall=%6.2fs cpu=%6.2fs rss=%4dMB md5=%s\n" \
                "$label" "$kind" "$rep" "${OMP_NUM_THREADS:-default}" \
                "$wall" "$cpu" $((rss/1024)) "$md5" | tee -a "$timings"
        done
    done
done
rm -f "$out_dir/time.tmp"
