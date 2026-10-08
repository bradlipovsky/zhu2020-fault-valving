# Equations and numerical choices

The implemented physical equations are Zhu et al. (2020), Eqs. (1)–(9), with
parameters in their Table 1 and the explicit assumptions in `independence.md`.
There is no thermal pressurization, porosity evolution, poroelastic traction, or
inelastic country-rock deformation in the published model or this implementation.

## Elasticity

Let d(z,t) be slip accumulated since initialization minus Vp t. Expand d in
cos(j pi z/Lz). The harmonic continuation into 0<y<Ly for each nonconstant mode
is proportional to sinh[kappa(Ly-y)]/sinh(kappa Ly), with kappa=j pi/Lz.
It satisfies the prescribed side displacements and traction-free top and base.
The fault stiffness eigenvalues are

    K_0 = mu/(2 Ly)
    K_j = mu kappa coth(kappa Ly)/2, j>0.

Thus tau_qs=tau_initial-Kd. DCT-II and DCT-III implement this operator on cell
centers z_i=(i+1/2)Lz/N. Their composition is normalized by 2N. The factor 1/2
is the conversion from fault slip to one-sided displacement. This is the
finite-domain operator, not an infinite-half-space approximation. The initial
prestress is equivalent to an initial slip field -K^{-1} tau_initial with zero
remote displacement; only subsequent accumulated slip is plotted.

The production fault mesh is graded. Linear interpolation P maps its slip to
the auxiliary cosine grid. Tractions return by the weighted adjoint
W^(-1) P^T h, where W contains physical control-volume widths. This preserves
constant traction and weighted elastic reciprocity. The complete construction,
mesh-refinement family, and verification are given in `graded_mesh.md`.

## Friction and state

At each explicit stage solve the scalar, monotone equation

    tau_qs = N a asinh[V exp(psi/a)/(2 V0)] + eta_rad V.

The velocity root is bracketed and solved in log|V|. Stable asinh(exp(x)) and
log(sinh(x)) expressions avoid exponential overflow. State evolves according to
the aging law, expressed in the paper's psi variable, with slip speed used in
the reverse-slip extension documented separately. The k* state is represented
by u=(k*-kmin)/(kmax-kmin), for which

    du/dt = |V|(1-u)/L - u/T.

The velocity and effective stress are not clipped. Invalid trial states reject
the step; an invalid accepted state or a time step below 1e-10 s terminates
the calculation and cannot produce a `completed.json` record.

## Fluid flow

We evolve excess pressure e=p-rho g z, giving q=(k/eta) de/dz. Integrating
fluid conservation over a cell gives n beta width_i de_i/dt=q_{i+1/2}-q_{i-1/2}.
Interior face permeability is the harmonic mean of adjacent cell values.
Pressure gradients use the actual distance between adjacent physical nodes.
At the surface the pressure is zero at the boundary face, at its actual distance
from the first unknown; the face permeability is the harmonic mean of the first
cell and its extrapolated zero-effective-stress boundary value k*_0. This is
the half-cell distance for a uniform mesh. At the bottom the
flux is the prescribed q0. There are no internal fluid sources.

Each implicit pressure stage uses fixed-point iteration: form k(e), assemble
and solve the tridiagonal backward-Euler operator, and repeat until the largest
pressure update is below 0.001 Pa. Iterations are damped after eight updates.
A failed nonlinear solve rejects the time step. Because fluxes telescope, the
discrete scheme conserves storage to the accuracy of the implicit solve.

The initial Darcy profile is solved sequentially from the surface by bisection
of each nonlinear face-flux relation, imposing q=q0 at every face. It is a
discrete steady solution of the same operator used in time integration.

An independent analytical check avoids testing that discrete solution only
against its own flux formula. Let G=normal_gradient-rho*g, C=kstar-kmin,
keq=viscosity*influx/G, and D=keq-kmin. With constant kstar, the continuous
steady-flow equations give dN/dz=G-viscosity*influx/[kmin+C exp(-N/stress_scale)].
Separating variables and imposing N(0)=0 yields

    z = N/G - stress_scale*keq/(G*D)
                  * log[(C-D exp(N/stress_scale))/(C-D)].

The tested parameters satisfy C>D>0. In the zero-floor limit this inverts to

    N(z) = -stress_scale * log[keq/C + (1-keq/C) exp(-G*z/stress_scale)].

`tests/steady.cpp` compares the production initial-pressure solver with both
solutions on four meshes. Maximum pressure errors, normalized by stress_scale,
converge at second order. This checks the nonlinear hydraulic resistance and
surface treatment against a continuum solution, independently of time stepping.

## Time integration

The production method is ARK4(3)6L[2]SA, a six-stage fourth-order additive
Runge–Kutta scheme with a third-order embedded estimate. Its exact rational
coefficients are written in `struct ARK4` in `src/model.cpp`. They are from
Kennedy and Carpenter (2001), NASA/TM-2001-211038, Appendix C, checked against
the SUNDIALS 7.7 written Butcher tables:
https://sundials.readthedocs.io/en/v7.7.0/arkode/Butcher_link.html.
We implement the table directly; no earthquake simulator is a source.

Let F contain slip, state, and k* rates, and G contain the pressure rate.
For explicit coefficients Ae and implicit coefficients Ai, each stage satisfies

    Y_s = Y_0 + dt sum_{j<s} [Ae_sj F(Y_j) + Ai_sj G(Y_j)]
                + dt Ai_ss G(Y_s).

F has zero pressure component and G has zero mechanical/state components.
Thus the implicit stage is only a tridiagonal nonlinear pressure solve; all
other stage fields are already known. The pressure derivative is evaluated
from conservative face fluxes at the converged stage. The accepted solution
uses the fourth-order weights b4; the estimate is dt sum (b4-b3)(F+G).
Error is bounded in the maximum norm after scaling slip by dc, psi by a,
normalized k* by 0.01+|u|, and pressure by sigma*. The controller uses
0.9 error^(-1/4), bounded between 0.15 and 2. Configurations state the tolerance
and maximum step. Invalid final values reject the trial just as invalid stages do.

ARS(2,2,2) with step doubling is retained as an independent, lower-order
comparison, and for simple limiting tests. It is not used for final figure data.
The coupled verification test measures fourth-order convergence of the production
method and second-order convergence of that comparison. The exploratory
rupture benchmark and the reason for changing the method are recorded in
`integrator_benchmark.json`; all final production runs start anew.

This integrator is an explicit deviation from the paper's stage-split adaptive
Runge–Kutta/backward-Euler method. Convergence is checked for constant-coefficient
diffusion, the constitutive ODEs, and a fully coupled nonlinear calculation.
Earthquake sequences additionally require the separate mesh and time-tolerance
experiments; equation verification alone is not trajectory convergence.

## Stored data and analysis

`fields.bin` starts with ZHUIND01, an unsigned 64-bit depth count, and little-endian
64-bit depth coordinates. Each record holds time (float64), accepted macro-step
(int64), and six float32 depth profiles in order: cumulative slip, signed velocity,
effective normal stress, permeability, k*, upward flux. The numerical integration
itself uses float64 throughout. These snapshots cover the upper 30 km with about
30 m output spacing; spatial output subsampling does not change the solution grid.
Every accepted step also writes a text history with maxima and four depth traces.
The trace-column names denote nominal depths of 5, 10, 15, and 20 km. Each trace
uses the nearest physical model node, choosing the deeper node in an exact tie.
The S1 panel files retain both these nominal depths and the actual node centers;
their sidecars record the same coordinates. `common.py` reconstructs the mesh
from the run configuration for this metadata. Diagnostics computed from binary
profiles instead use the nearest stored profile coordinate, which can differ
because profiles are spatially subsampled.
`data/depth_coordinate_verification.json` records exact agreement of every
reconstructed physical node with the C++ mesh for the 12 scientific cases and
the two software test configurations.
It additionally records the minimum numbers of physical cell widths across
Lb and h* in the velocity-weakening region, using the current effective stress.
The text histories use 17 significant decimal digits, preserving float64 times.
Binary snapshot times likewise retain full float64 precision.

Nominal snapshot intervals are 0.025 yr during slow slip and 0.5 s during an
earthquake; each snapshot is taken at the first accepted step reaching its
target interval, so the actual spacing can exceed it by one step.
Additional snapshots capture crossings of 1e-3 m/s and changes
exceeding 0.2 in log10 maximum speed. The exact snapshot times and macro-step
indices are retained. Coseismic profile contours use nearest snapshots to a
one-second sequence, not a claimed exact one-second sampling of the integrator.

`analyze.py` explicitly defines seismic events, large ruptures, selected cycles,
and migration fits. Events use the maximum velocity recorded at every accepted
step. Their rupture footprints use at least 1 cm of accumulated slip between
the saved first-above-threshold and first-below-threshold profiles; a local velocity maximum can be missed by
0.5 s snapshots and is therefore not an additional footprint requirement.
Each crossing is saved explicitly by the executable. The catalog also retains
the onset-time bracket between the preceding accepted step and the first step
above threshold. Slip during that initial bracket is not included in the footprint;
the exact measured slip interval is recorded separately. Using an earlier profile
would incorrectly include nucleation creep in a brief threshold excursion.
Large events have a connected footprint reaching above 2 km and spanning more
than 10 km. The last complete large-event cycle is selected by a fixed
rule. The published time origins and restarts are never imported or fitted.
Threshold events truncated at either output boundary are retained with
complete=false and cannot close a cycle or enter complete-event comparisons.
The report compares the first two complete partial events in chronological
order, using each event's longest connected footprint. Maximum local slip does
not select the events because a later swarm event can exceed an earlier
partial rupture. The full catalog remains available in each analysis file.
Small-rupture counts require a resolved 1 cm slip footprint. Complete threshold
excursions without that footprint are counted separately, not discarded or merged.
This distinguishes brief speed-threshold recrossings from additional resolved
ruptures in the sequence comparison; no detection threshold is adjusted.
The time of maximum effective stress near 10 km provides a separately defined
drainage diagnostic; it is not equated automatically to a published phase boundary.
Refinement and input-sensitivity comparisons use the common elapsed duration
available in every validation run, including startup. Complete seismic-event
counts and large-event recurrence intervals exclude events ending beyond that
duration. Large-event timing differences pair events in chronological order,
without fitting a time shift or choosing the closest event. All paired differences
are retained, and event counts expose unequal catalog lengths.
The same chronological large-event pairs also compare duration, peak speed, and
maximum slip. The report gives their largest absolute relative differences and
the numerical summary retains every signed difference. These measure sensitivity
of the sequence, including later ruptures; event-number pairing does not establish
physical correspondence when event counts or rupture patterns differ.
Figure 6 uses the independently computed |V|=Vp contour; its definition is an
analysis assumption because the paper does not specify its exact extraction rule.
The quantitative tracker retains two choices at every saved time: the shallowest
and deepest crossings where speed increases through Vp with depth, between
2 and 23 km, interpolated in log speed. The shallowest (leading) choice describes
the outer boundary; the deepest can follow a secondary pulse behind it.
Both are reported for every case without choosing the result closest to a
published speed. Figure 6 itself retains all contours. The two tracked extrema
can change branches; fit segments split at depth increments outside -250 to
+30 m, gaps longer than 0.25 year, or seismic intervals. Thus these diagnostics
do not assign persistent identities to individual propagating pulses.
Deep migration fits use 13--20 km; separate shallow fits use 2--10 km for
comparison with the printed 4.56 km/year shallow-front annotation. Both require
at least 1 km of upward motion over 0.1 year and report the median of segments
with R-squared at least 0.8 for each contour choice. All segment fits are retained
in the analysis file, and `front.csv` stores both depth series. A missing fit does
not establish that no propagating feature exists; it means the chosen contour
and fit criteria do not provide a qualifying measurement.
Slow-slip diagnostics at 15, 18, and 20 km count complete local intervals above
1.1 Vp, longer than 0.01 yr, during which the maximum slip speed anywhere remains
below 1e-3 m/s. Truncated intervals at a selected window edge are excluded.
The seismicity exclusion uses the whole-fault maximum at every accepted step,
including depths beyond the saved profile range and times between snapshots.
The threshold is near the plate rate because the article describes centimetres
of slip over about a year; a 10 Vp threshold would exclude that amplitude scale.
This is a declared measurement convention, not a prescribed condition in the model.
