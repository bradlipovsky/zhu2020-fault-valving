#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$root/build/fftw-source"
cd "$root/build/fftw-source"
curl -fL https://www.fftw.org/fftw-3.3.10.tar.gz -o fftw.tar.gz
echo '56c932549852cddcfafdab3820b0200c7742675be92179e59e6215b340e26467  fftw.tar.gz' | sha256sum -c -
tar -xzf fftw.tar.gz
cd fftw-3.3.10
./configure --prefix="$root/build/fftw" --enable-shared --disable-static --disable-fortran
make -j8
make install
