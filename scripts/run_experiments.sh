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
python3 scripts/check_archive_units.py
python3 scripts/initial_profiles.py
phase_end() {
  python3 - "$1" "$2" <<'PY'
import json,sys
from pathlib import Path
p=Path('results/run_provenance.json')
r=json.loads(p.read_text()).get('optimized_continuations',{}).get('runs',{}) if p.exists() else {}
print(r[sys.argv[1]]['start_year'] if sys.argv[1] in r else sys.argv[2])
PY
}
thread_end() {
  python3 - "$1" "$2" <<'PY'
import json,sys
from pathlib import Path
p=Path('results/run_provenance.json')
r=json.loads(p.read_text()).get('threaded_continuation',{}) if p.exists() else {}
print(r['start_year'] if r.get('run')==sys.argv[1] else sys.argv[2])
PY
}
for n in 8192 16384 32768; do
  build/valving --n "$n" --years 12 --T 1e8 --rk-max-dt 1 \
    --profile inputs/T1e8_restart.csv --output "results/restart_T1e8_n$n" \
    > "results/restart_T1e8_n$n.log" 2>&1
done
for n in 8192 16384 32768; do
  cp "results/restart_T1e8_n$n/final_state.bin" "results/restart_T1e8_n$n/phase1_state.bin"
  end=35
  if [ "$n" = 32768 ]; then end=16; fi
  split=$(phase_end "restart_T1e8_n$n" "$end")
  build/valving --n "$n" --years "$split" --T 1e8 --rk-max-dt 1 --profile inputs/T1e8_restart.csv \
    --state "results/restart_T1e8_n$n/phase1_state.bin" --start-year 12 --append \
    --output "results/restart_T1e8_n$n" >> "results/restart_T1e8_n$n.log" 2>&1
  if [ "$split" != "$end" ]; then
    cp "results/restart_T1e8_n$n/checkpoint_state.bin" "results/restart_T1e8_n$n/phase2_state.bin"
    build/valving --n "$n" --years "$end" --T 1e8 --rk-max-dt 3000 --profile inputs/T1e8_restart.csv \
      --state "results/restart_T1e8_n$n/phase2_state.bin" --start-year "$split" --append \
      --output "results/restart_T1e8_n$n" >> "results/restart_T1e8_n$n.log" 2>&1
  fi
done
build/valving --n 16384 --years 1 --T 1e8 --rtol 2.5e-5 --rk-max-dt 1 \
  --profile inputs/T1e8_restart.csv --output results/time_refinement \
  > results/time_refinement.log 2>&1
for tolerance in 1e-4 2.5e-5; do
  name=seismic_time_baseline
  if [ "$tolerance" = 2.5e-5 ]; then name=seismic_time_refinement; fi
  build/valving --n 16384 --years 14.7 --T 1e8 --rtol "$tolerance" --rk-max-dt 1 \
    --profile inputs/T1e8_restart.csv --state results/restart_T1e8_n16384/phase1_state.bin \
    --start-year 12 --output "results/$name" > "results/$name.log" 2>&1
done
build/valving --n 16384 --years 14.7 --T 1e8 --rtol 1e-4 --rk-max-dt 3000 \
  --profile inputs/T1e8_restart.csv --state results/restart_T1e8_n16384/phase1_state.bin \
  --start-year 12 --output results/time_stepping_candidate > results/time_stepping_candidate.log 2>&1
for exponent in 7 9 10; do
  duration=25
  if [ "$exponent" = 7 ]; then duration=12; fi
  build/valving --n 16384 --years "$duration" --T "1e$exponent" --rk-max-dt 1 \
    --profile "inputs/T1e${exponent}_restart.csv" \
    --output "results/restart_T1e${exponent}_n16384" \
    > "results/restart_T1e${exponent}_n16384.log" 2>&1
  cp "results/restart_T1e${exponent}_n16384/final_state.bin" "results/restart_T1e${exponent}_n16384/phase1_state.bin"
  case "$exponent" in
    7) end=48;;
    9) end=80;;
    10) end=180;;
  esac
  split=$(phase_end "restart_T1e${exponent}_n16384" "$end")
  build/valving --n 16384 --years "$split" --T "1e$exponent" --rk-max-dt 1 \
    --profile "inputs/T1e${exponent}_restart.csv" \
    --state "results/restart_T1e${exponent}_n16384/phase1_state.bin" --start-year "$duration" --append \
    --output "results/restart_T1e${exponent}_n16384" \
    >> "results/restart_T1e${exponent}_n16384.log" 2>&1
  if [ "$split" != "$end" ]; then
    cp "results/restart_T1e${exponent}_n16384/checkpoint_state.bin" "results/restart_T1e${exponent}_n16384/phase2_state.bin"
    thread_split=$(thread_end "restart_T1e${exponent}_n16384" "$end")
    build/valving --n 16384 --years "$thread_split" --T "1e$exponent" --rk-max-dt 3000 \
      --profile "inputs/T1e${exponent}_restart.csv" \
      --state "results/restart_T1e${exponent}_n16384/phase2_state.bin" --start-year "$split" --append \
      --output "results/restart_T1e${exponent}_n16384" >> "results/restart_T1e${exponent}_n16384.log" 2>&1
    if [ "$thread_split" != "$end" ]; then
      cp "results/restart_T1e${exponent}_n16384/checkpoint_state.bin" "results/restart_T1e${exponent}_n16384/phase3_state.bin"
      build/valving --n 16384 --years "$end" --T "1e$exponent" --rk-max-dt 3000 \
        --profile "inputs/T1e${exponent}_restart.csv" \
        --state "results/restart_T1e${exponent}_n16384/phase3_state.bin" --start-year "$thread_split" --append \
        --output "results/restart_T1e${exponent}_n16384" >> "results/restart_T1e${exponent}_n16384.log" 2>&1
    fi
  fi
done
python3 scripts/plot_figures.py
python3 scripts/analyze_cpp.py
(cd report && latexmk -pdf -interaction=nonstopmode -halt-on-error reproduction.tex)
python3 scripts/check_outputs.py
