# Graded fault discretization

Allison and Dunham (2018), Section 4, discuss variable spatial spacing and the
frictional length scales h*=mu*dc/[N*(b-a)] and Lb=mu*dc/(N*b). Their Eq. (21)
is the more restrictive criterion. The initial uniform 15.26 m spacing left
approximately one cell across Lb for some long-healing cases. Those exploratory
runs are excluded from the final panels. Finer spacing is needed to assess the
rupture behavior rather than merely draw a smoother figure.

The physical fault nodes form a subset of an auxiliary uniform cosine grid.
The finest spacing is h=Lz/n. Starting from depth (j+1/2)h, advance by the nearest
integer to max[1, surface_ratio*exp(-z/surface_scale)]. Below fine_depth, instead
advance by the nearest integer to min[maximum_spacing/h,
exp((z-fine_depth)/stretch_scale)]. Thus the shallow transition and the deep
stretching prescribe numerical spacing, not physical constitutive changes.
The first index is floor[(surface_ratio-1)/2]. The final interval is combined
with the preceding interval if otherwise shorter than half its target width.
The construction is explicit in `mesh_indices` in `src/model.cpp`.

The production auxiliary grid has 131,072 points (h=3.8147 m). The fault has
finest spacing through the seismogenic depths and down to 30 km. Surface spacing
decreases from about 8h to h over the upper 2 km; it coarsens below 30 km with
stretch_scale=10 km toward the target maximum_spacing=1 km. The final interval
can be up to about 1.5 times that target to avoid a tiny interval at the bottom.
The refinement cases halve/double both h and this spacing target, preserving
the mesh family. Surface spacing
therefore also refines; no fixed stress or velocity floor replaces resolution.

Let P linearly interpolate physical slip to the auxiliary grid, extending its
first value to the surface. Let W be the diagonal matrix of control-volume
widths, with boundaries at adjacent-node midpoints and at the domain ends.
The physical stiffness operator is

    K_physical = W^(-1) P^T h K_cosine P.

Here K_cosine is the finite-domain operator already derived in `numerics.md`.
For this subset mesh, W_ii=sum_j h P_ji. Consequently a constant slip has exactly
the original constant traction. Also W K_physical is symmetric positive definite:
its elastic work equals the work of the interpolated slip on the cosine grid.
This is a Galerkin projection with lumped mass, not point sampling of tractions.
Only the numerical approximation changes; the same elastic boundary conditions
and rate-and-state laws apply at all physical nodes.

For fluid flow, replace uniform h in the storage balance by each control-volume
width, and replace h in the pressure gradient by the distance between adjacent
nodes. The first pressure gradient uses the actual surface-to-node distance.
Harmonic face permeability and the prescribed bottom flux retain the conservative
telescoping balance. The initial nonlinear steady flux is solved on this same mesh.

`tests/graded.cpp` checks constant traction, weighted elastic reciprocity,
positive elastic energy, agreement with uniform-grid tractions from localized
slip, nonlinear fluid storage conservation, and spatial convergence to the
analytical steady-pressure and transient diffusion solutions. These checks do
not establish convergence of nonlinear earthquake sequences. The separate
coarse, fine, and tighter-time-tolerance experiments address that question.

At every accepted step, `history.csv` records the smallest Lb/width and
h*/width in the velocity-weakening region. `completed.json` retains their
minima over the entire run. These are resolution diagnostics rather than
constraints imposed on the physical solution.

The optional `mesh_compare` target compares two independently generated final
checkpoints at exactly the same physical time. It reports maximum differences
over a declared depth interval, linearly interpolating the second spatial mesh
only where coordinates differ. Interseismic benchmarks cannot replace rupture
convergence. Benchmark timings exclude initial setup and depend on hardware load.

## Controlled comparisons

The recorded experiments are in `graded_benchmark.json`. Each begins from the
independently specified initial conditions. No checkpoint is imported to start
an experiment. Checkpoints are read only to compare its completed final state.

At 15.26 m finest spacing and tolerance 0.001, the uniform and deep-graded
meshes completed five years in 522.38 and 36.66 s, respectively. Between 3 and
25 km depth, their maximum differences were 2.63 micrometres of slip, 179.3 Pa
of pressure, and 9.49e-12 m/s in velocity. This is an interseismic comparison;
it does not establish earthquake convergence.

At 1.91 m finest spacing and tolerance 1e-5, adding shallow grading reduced
the one-year cost from 213.46 to 34.39 s. Over the same 3--25 km interval the
maximum differences were 0.156 micrometres of slip, 2.21 Pa of pressure, and
2.86e-13 m/s in velocity. The two tight calculations use the final production
source, as their recorded SHA-256 hashes confirm.

Reducing tolerance from 0.001 to 1e-5 on the surface-graded mesh changed slip
by 31.7 micrometres and velocity by 3.04e-11 m/s over one year. The runtimes
were almost equal (34.00 and 34.39 s). This motivates the production tolerance
of 1e-5 and the additional 2.5e-6 sensitivity run. The larger loose-tolerance
difference must not be attributed entirely to mesh geometry.

To repeat the six controlled experiments, run

    python3 scripts/mesh_benchmark.py --directory .tmp/repeated-mesh-benchmark

The directory must be new. This optional calculation is separate from the
main figure pipeline and writes its own configurations, logs, and comparisons.
Its timings will vary with hardware and concurrent work.
