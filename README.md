# Independent fault-valving calculation

This revision implements the equations in Zhu et al. (2020),
*Fault valving and pore pressure evolution in simulations of earthquake sequences
and aseismic slip*, doi:10.1038/s41467-020-18598-z, from the published description.
It replaces the earlier reconstruction, which used archived authors' results.
The earlier attempt remains in the repository history and v1.0.0 release and is
not an input to this calculation.

The C++ model couples finite-domain antiplane elasticity, regularized rate-and-state
friction, one-dimensional fault-zone Darcy flow, pore pressure, and permeability
evolution. FFTW supplies cosine transforms. Python supplies output analysis and
plotting. All physical constants, assumptions, and numerical settings are written
in `configs/*.cfg`; see `documentation/independence.md`.

## Regeneration

Requirements: Linux, a C++17 compiler, CMake, OpenMP, Python 3 with NumPy and Matplotlib,
and LaTeX with `latexmk`. `scripts/bootstrap.sh` installs a pinned FFTW release
under `.deps` if required. It needs internet access only for that numerical library.
The model and figure pipeline do not download scientific data.

```
python3 -m pip install -r requirements.txt
python3 scripts/reproduce.py
```

This compiles and verifies the solver, runs every main case and the declared
validation cases, generates panel data and figures, and compiles
`report/reproduction.pdf`. Complete runs are reused only when their source and
configuration hashes match. `--fresh` moves previous generated cases to `.tmp`
before starting a new calculation. Do not use it while another run is active.
The command also creates `report/latex-source.zip` and verifies that its extracted
LaTeX sources and generated figures compile in an isolated directory.
`report/reproduction-complete.tex` is the expanded single-file LaTeX source.

```
# A short smoke calculation
bash scripts/bootstrap.sh
build/valving configs/smoke.cfg data/smoke

# A single physical experiment
build/valving configs/baseline.cfg data/baseline
python3 scripts/analyze.py baseline

# Rebuild figures and report after the simulations have finished
python3 scripts/reproduce.py --postprocess-only
```

`data/run_status.json` and per-case `run.log` expose live progress. The numerical
data are generated entirely by this executable; the paper is a comparison target.
Results and reproduction status will be tabulated panel by panel in the report.
No published panel is substituted for a generated numerical result.

The complete software pipeline can also be checked with
`python3 scripts/smoke_pipeline.py`. It creates its own directory under `.tmp`,
uses deliberately short and coarse simulations, and verifies all panel products
and both report builds. Those test outputs are never used as scientific results.
