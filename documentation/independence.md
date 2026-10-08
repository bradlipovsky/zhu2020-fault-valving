# Scope and provenance

This implementation was started in an empty worktree on 8 October 2026 UTC.
The previous attempt's source, restart states, simulation results, plotting code,
and figure products are excluded. The only files retained from that tree are
the requested AGENTS.md and academic style guide. Git history and the v1.0.0
release preserve the earlier attempt, which must not be used to regenerate
this revision.

Scientific inputs consulted:

* Zhu, Allison, Dunham, and Yang (2020), Nature Communications 11, 4833,
  doi:10.1038/s41467-020-18598-z: published PDF and article XML from Europe PMC.
* The published Supplementary Information PDF (MOESM1): six figure captions
  and figures, inspected only as reference targets.
* Allison and Dunham (2018), Tectonophysics 733, 232–256,
  doi:10.1016/j.tecto.2017.10.021: written discussion of the friction profile.

No authors' code, scripts, simulation archives, source-data files, or digitized
output curves are inputs to the new solver or plotting pipeline. Reference
documents are held outside the tracked build inputs in `.reference/`.
FFTW is an ordinary numerical transform library, independently downloaded from
fftw.org. NumPy and Matplotlib are used only for handling and plotting our output.

The attachment defining the revised task ends with “LaTeX source, and”. The
compiled PDF remains part of the deliverable under the original request.

## Decisions required by incomplete documentation

All units in configurations and data are SI, except names explicitly ending
in `_years`. A year is 365.25 days.

1. The normal stress gradient is not tabulated. We assume 22 MPa/km, consistent
   with the straight normal-stress line displayed in Fig. 1d. This is a model
   input assumption, not a numerical reproduction of that line.
2. The exact a(z), b(z) arrays are not provided. We specify a piecewise-linear
   approximation to the parameter diagram in Fig. 1a: a=0.01 at the surface,
   0.03 at 15 km, and slope 0.003/km below; a-b=-0.01 to 13.5 km and increases
   linearly with slope 0.01/3.5 km below, crossing zero at 17 km. Both profiles
   continue to the 500 km domain boundary. These are declared assumptions.
3. Initial k* is its analytic steady-sliding value. Initial p is the discrete
   steady Darcy solution with prescribed bottom influx and zero surface p.
4. Initially the fault is nearly locked above 17 km and creeps at Vp below,
   with a tanh transition of width 1.5 km and a shallow rate of 0.0001 Vp.
   The initial state variable is the steady aging-law value at Vp. The initial
   elastic prestress balances the specified initial velocity and radiation damping.
   A small, explicitly configured Gaussian decrease in the state variable near
   12 km starts the instability; no earthquake times or rupture footprints are
   prescribed. The initial displacement giving this prestress exists through
   inversion of the positive elastic stiffness operator.
   An exploratory all-depth steady-sliding start was also checked: shallow
   elastic stiffness severely restricted the explicit mechanical time step.
   Its outputs are not used in the manuscript panels. The declared locked-start
   profile is a physical initial-condition assumption, not an imposed constraint
   on subsequent slip. Steady sliding remains a solver verification case.
5. We solve the documented rectangular elastic boundary-value problem using
   its cosine eigenfunctions, rather than reproduce the fourth-order SBP
   discretization. Cell centers avoid the degenerate zero-effective-stress
   surface node. Pressure uses a conservative second-order finite-volume
   scheme and harmonic face permeability. These numerical changes require
   verification and refinement, not a claim of identical trajectories.
6. Time integration uses the independently implemented second-order ARS(2,2,2)
   IMEX Runge–Kutta method, with nonlinear implicit pressure stages and explicit
   slip, state, and k*. Step doubling estimates local error in all four fields.
7. Signed velocity is allowed by regularized friction. State evolution and
   permeability enhancement use slip speed |V|, extending the paper's positive
   sliding convention to possible reverse slip. No velocity, effective stress,
   or permeability floor is imposed beyond the published k_min law. States
   with nonpositive effective stress are reported as failures; opening physics
   is outside the documented model.
8. Event definitions, analysis windows, and front tracking are analysis choices,
   explicitly recorded with the generated metrics. First cycles are treated as
   transients. Differences in event phase cannot be eliminated using an authors'
   restart, data-derived forcing, or manual trajectory adjustment.
