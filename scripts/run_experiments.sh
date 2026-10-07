#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p build/tmp
export TMPDIR="$PWD/build/tmp"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
ctest --test-dir build --output-on-failure
build/valving --test > results/verification.txt
python3 scripts/download_data.py
python3 scripts/prepare_archive.py
python3 scripts/initial_profiles.py
for n in 8192 16384 32768; do
  build/valving --n "$n" --years 12 --T 1e8 \
    --profile inputs/T1e8_restart.csv --output "results/restart_T1e8_n$n" \
    > "results/restart_T1e8_n$n.log" 2>&1
done
for exponent in 7 9 10; do
  duration=25
  if [ "$exponent" = 7 ]; then duration=12; fi
  build/valving --n 16384 --years "$duration" --T "1e$exponent" \
    --profile "inputs/T1e${exponent}_restart.csv" \
    --output "results/restart_T1e${exponent}_n16384" \
    > "results/restart_T1e${exponent}_n16384.log" 2>&1
done
python3 scripts/plot_figures.py
python3 scripts/analyze_cpp.py
(cd report && latexmk -pdf -interaction=nonstopmode -halt-on-error reproduction.tex)
