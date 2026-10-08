#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
root="$PWD"
mkdir -p .tmp .build-deps .deps
export TMPDIR="$root/.tmp"
if [[ ! -f .deps/lib/libfftw3.so ]]; then
  archive="$root/.build-deps/fftw-3.3.10.tar.gz"
  curl -L --fail --retry 2 https://www.fftw.org/fftw-3.3.10.tar.gz -o "$archive"
  printf '%s  %s\n' 56c932549852cddcfafdab3820b0200c7742675be92179e59e6215b340e26467 "$archive" | sha256sum -c -
  mkdir -p .build-deps/fftw
  tar xzf "$archive" -C .build-deps/fftw --strip-components=1
  (
    cd .build-deps/fftw
    ./configure --prefix="$root/.deps" --enable-shared --enable-sse2 --enable-avx
    make -j8
    make install
  ) > .build-deps/fftw-build.log 2>&1
fi
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
ctest --test-dir build --output-on-failure
python3 -c 'import numpy, matplotlib; print("Plotting dependencies:", numpy.__version__, matplotlib.__version__)'
