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
The fault mesh resolves the rupture region at 3.81 m and coarsens at depth.
`documentation/graded_mesh.md` describes its conservative fluid operator,
energy-preserving elastic projection, and refinement checks.

The [interim scientific report](report/interim.pdf), dated 9 October 2026 (UTC),
contains all 49 generated panel packages from eleven completed cases and direct
constitutive calculations. Two constitutive-law panels are independently
reproduced; 47 are partial comparisons, including the newly drawn schematic.
The main and tighter-tolerance 200-year cases have finished. The fine-mesh run
continues; this is not the final reproduction report.
Its [LaTeX source](report/interim.tex) and
[self-contained source archive](report/interim-latex-source.zip) are included.
Run `python3 scripts/interim_report.py` to rebuild this report from completed,
audited outputs. The [build record](data/interim_report_verification.json)
includes input hashes and an isolated source-archive compilation.
The main case's final complete interval is 31.44 years versus approximately
32 years in the paper, but its event sequence, stress variation, and migration
rates differ. Its [supplementary diagnostics](data/baseline/supplementary_verification.json)
also show that a surface-reaching 9.48 km rupture is excluded by the declared
10 km large-event cutoff; recurrence statistics depend on that event definition.
The [full temporal-tolerance comparison](data/tolerance_200yr_verification.json)
finds close agreement through 100 years, but 60 versus 58 resolved ruptures over
200 years. Large-event counts are also sensitive to the span cutoff. These checks
do not demonstrate convergence of the individual late events used in the figures.

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
The final audit checks every accepted history row and every saved field profile,
including times outside the figure windows. It records per-case results and raw
file hashes in `data/*/output_verification.json`. These consistency checks do not
establish numerical convergence or agreement with the paper. Run them separately
after simulations finish with `python3 scripts/verify_outputs.py` (or append case
names to check selected completed runs).

The committed `data/panels/*.npz` files contain the plotted numerical arrays.
Full field files, accepted-step histories, and checkpoints are generated locally
by the numerical runs. The supplied panel files support plotting directly with
`python3 scripts/plot.py all`.

When a fresh clone contains completed-run metadata but no raw field file, the
runner archives that metadata and starts a new simulation. An isolated
[regeneration check](data/cache_regeneration_verification.json) verified this
path and obtained bit-identical field and history files for the short software
fixture; those test outputs are excluded from scientific results.

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
The fixed-pressure reference runs for 400 years: an initial 200-year trial
contained only one large earthquake and no recurrence interval. Its duration
check is recorded in `documentation/reference_duration.json`; the shorter trial
is excluded from the final figure inputs. Concurrent case batches can use
`--status-file` to keep their progress records separate.
Results and reproduction status will be tabulated panel by panel in the report.
No published panel is substituted for a generated numerical result.

The complete software pipeline can also be checked with
`python3 scripts/smoke_pipeline.py`. It creates its own directory under `.tmp`,
uses deliberately short and coarse simulations, and verifies all panel products
and both report builds. Those test outputs are never used as scientific results.
`python3 tests/output_audit.py <fixture-directory>` then checks that the raw-output
audit rejects deliberately damaged copies of that fixture. The test records its
results in `data/output_audit_verification.json` and leaves the fixture intact.
