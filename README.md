# Earthquake cycles, permeability evolution, and fault valving

Reproducibility study of Zhu, Allison, Dunham, and Yang (2020),
[Fault valving and pore pressure evolution in simulations of earthquake sequences and aseismic slip](https://doi.org/10.1038/s41467-020-18598-z).

The project includes a new C++ implementation of the coupled model, a reconstruction
of the six main and six supplementary figures from the published equations and
the authors' archived fields, and a compiled [report](report/reproduction.pdf)
with its [LaTeX source](report/reproduction.tex). The report distinguishes
independent calculations from replotting archived simulations. Redrawing a
published earthquake sequence does not validate the new solver's prediction of
that sequence.

## Model

`src/valving.cpp` implements quasi-dynamic antiplane elasticity, regularized
rate-and-state friction with the aging law, vertical Darcy flow, stress-dependent
permeability, and slip enhancement and time-dependent healing of permeability.
All solver quantities use SI units. FFTW evaluates the exact continuum elastic
Dirichlet-to-Neumann map on a rectangular domain with traction-free top and
bottom boundaries. Fluid transport uses conservative second-order finite volumes.
Adaptive midpoint integration uses step doubling to control all four evolving
fields. At long time steps, a Newton/CG solve treats mechanics implicitly;
short steps use an explicit mechanical midpoint. Pressure is implicit in both.

The actual Scycle inputs are preserved in `reference/original_inputs/` under
their original license. They contain the friction profiles, a 1 MPa effective
stress floor, and a taper in permeability bounds between 30 and 60 km. These
details are part of the original implementation. The archived depth grid ends
at 578.952681 km; the input file's nominal depth is 500 km. The new solver
defaults to the archived extent and accepts `--height 500000` for the nominal
domain. Its discretization and initialization differ from Scycle; see the report.

## Build and verify

Dependencies: a C++17 compiler, CMake, FFTW3, Python packages listed in
`requirements.txt`, and a LaTeX installation with `latexmk`.

```sh
# If FFTW3 development files are unavailable:
bash scripts/install_fftw.sh
mkdir -p build/tmp
export TMPDIR="$PWD/build/tmp"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
ctest --test-dir build --output-on-failure
build/valving --test
```

The verification checks elastic Fourier modes, scalar friction inversion,
steady flow preservation, a manufactured diffusion solution, discrete fluid
mass conservation, and second-order time convergence in a spring-slider limit.

## Rebuild the archived figures

```sh
python3 -m pip install -r requirements.txt
python3 scripts/download_data.py
python3 scripts/prepare_archive.py
python3 scripts/plot_figures.py
```

The download is about 5 GB. Every original file is checked against the SHA-256
value supplied by [OSF project 9YGRP](https://doi.org/10.17605/OSF.IO/9YGRP).
Raw files and plotting caches are excluded from Git; the manifest, conversion
code, generated figures, and measured diagnostics are included. In particular,
OSF permeability values must be multiplied by 10^6 to convert km² to m².
Pressure and effective stress in the archive are in MPa, depth is in km, slip
is in m, and time is in s.

Figure 1 is evaluated independently from the constitutive equations. Figures
2–5 and S1–S6 are reconstructed from archived outputs with revised layouts and
explicitly labelled provenance. Figure 6 uses an explicitly defined velocity
threshold and time windows because these plotting choices were not deposited.
The plot cache retains every saved time sample and every fifth depth node
(50 m spacing above 20 km); the underlying archive resolves 10 m there.
The horizontal sample-index panels count archived samples, which were saved
every ten integration steps.

## Run the C++ model

```sh
# Fresh initial steady sliding, with a specified small state perturbation.
build/valving --n 16384 --years 120 --T 1e8 --output results/fresh_T1e8
build/valving --n 16384 --years 120 --T 1e8 --fixed --output results/fixed_T1e8

# Approximate restart from the original archive, with reconstructed aging state.
python3 scripts/initial_profiles.py
build/valving --n 16384 --years 12 --T 1e8 \
  --profile inputs/T1e8_restart.csv --output results/restart_T1e8_n16384
```

`--n` specifies the number of depth intervals; 16384 and 32768 give 35.34 and
17.67 m spacing across the archived extent. Other supported options include
`--rtol`, `--max-dt`, `--height`, `--perturb`, and `--stride`.
Healing times are 1e7, 1e8, 1e9, or 1e10 s. The 1e7 s case uses the original
input flux 3.3764e-10 m/s; the paper rounds this to 3.3e-10 m/s. Other cases
use 3e-9 m/s. Earthquake cycles with a 2 mm state evolution distance are
expensive: grid refinement and event timing must be assessed before treating
a new sequence as a quantitative reproduction.

The archive omits frictional state and shear traction. Restart profiles
therefore reconstruct the aging state from slip history and infer shear
traction from frictional balance. These are approximate, archive-assisted
initial conditions, not the authors' original checkpoints. Their assumptions
and state-memory sensitivity are recorded beside each CSV.

Each run writes metadata, a progress history, and binary field records. Each
record contains a float64 time followed by seven float64 depth vectors: slip,
velocity, effective stress, permeability, upward face flux, frictional state,
and reference permeability. The metadata gives vector lengths and units.
Output is restricted to the upper 30 km; the final state retains the whole grid.

## Report and provenance

```sh
cd report
latexmk -pdf -interaction=nonstopmode -halt-on-error reproduction.tex
```

`AGENTS.md` and `brad-lipovsky-academic-style-guide.md` were imported from
[bradlipovsky/vdv-damage](https://github.com/bradlipovsky/vdv-damage).
The original solver branch is
[Scycle / Zhu_et_al_2020](https://bitbucket.org/kallison/scycle/src/Zhu_et_al_2020/),
commit `d65e248900e193add0cc901ee7ee4d8ee7a2d50a`.
The article and supplement are CC BY 4.0; original Scycle inputs retain their
MIT license in `reference/SCYCLE_LICENSE`. This project's source is MIT licensed.
FFTW is an external dependency with its own license and is not vendored.
