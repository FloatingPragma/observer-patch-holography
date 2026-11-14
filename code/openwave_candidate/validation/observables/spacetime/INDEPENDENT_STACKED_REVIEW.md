# Independent stacked-ring Bloch and observable audit

The extension is a supplied periodic stack of canonical ring cells. It does
not derive the existence of a third spatial direction or a gravitational
interaction. Uniform modular Gauss-fiber reduction applies to a connected
stacked graph with the same specified local-coarse laws; the demonstration
selects zero commit rate and the positive acoustic band. It does not substitute
this preparation for definite records.

A cell contains 24 ports and 63 horizontal seams. Adding one positive-z seam
for each port adds 24 transit components, giving 111 internal positions.
With inverse Fourier waves exp(+i q.n), a vertical edge from R to R+z gives
`B_z=(1-exp(-i q3)) I_24`. Therefore

```text
L_3(q) = L_ring(q1,q2) + 2(1-cos(q3)) I_24.
```

For horizontal glue A(R) to B(R-t_i), the same cell gauge gives the B-endpoint
incidence entry `-exp(+i q.t_i)`, where t_i=(0,0),(1,0),(0,1). The physical
B-carrier offset is D=(sqrt(3)/2,1/2,0). A different phase convention can be
unitarily equivalent, but its Fourier sign, seam anchoring and position
offsets must transform consistently.

Let u be a normalized lowest port-band eigenvector, lambda its nonnegative
3D Laplacian eigenvalue, E=sqrt(m²+lambda), with m>0. The normalized positive
H_m spinor is

```text
P=sqrt((E+m)/(2E)) u;
Q=B^dagger u/sqrt(2E(E+m)).
```

Independent tests verify the full 111-component eigen-equation and norm.
A packet needs a specified phase gauge across q; arbitrary independent
signs of eigenvectors preserve these pointwise tests but can destroy spatial
localization. The preparation phase, envelope and gauge must be retained.

Write the physical cell-coordinate matrix A with columns
`a1=(sqrt(3),0,0)`, `a2=(sqrt(3)/2,3/2,0)`, `a3=(0,0,1)`.
Then q=A^T k and r=A n+internal_offset. Consequently

```text
v_physical = A gradient_q E;
Hessian_k(E²) = A Hessian_q(E²) A^T.
```

The horizontal leading form is
`(15/352)*(4/9)*(q1²-q1*q2+q2²)`, while the vertical leading term is q3².
Thus physical coefficients are diag(15/352,15/352,1), with anisotropy retained.
Independent finite differences verify rank3; deleting vertical links gives
rank2 and zero z coefficient. This does not establish a shared isotropic
relativistic light cone.

A Gaussian momentum envelope evolved by exp(-iEt) gives a probability packet.
Before periodic wrap, its centroid follows the envelope-weighted group
velocity, subject to the declared position readout and gauge. The center-q
ray is a finite-width approximation. A complete position readout includes
ports and every horizontal and vertical transit component. Seam positions
are their declared midpoints, including cross-cell offsets. Norm-preserving
FFT normalization is required. A centroid is not a sampled hidden-particle
trajectory.

Finite periodic grids require either a circular/unwrap convention or a
quantified edge/wrap diagnostic. A naive arithmetic centroid changes when
probability crosses the torus boundary. Grid refinement, packet widths and
ray deviation are numerical checks, not proof of a continuum geometry. The
absent-z control should retain a stationary z marginal and zero vertical
group velocity under its separately prepared positive band.

`independent_stacked_reference.py` reconstructs the Bloch incidence without
importing the demonstration. `test_independent_stacked.py` checks the vertical
identity, spinor and physical metric. Direct demo tests also check its Bloch matrices, phase-normalized positive-band
projectors, Hellmann-Feynman velocity, orthonormal FFT, missing-z marginal and
saved full position apparatus. Ten tests pass. The supplied seven-offset-bin
distributions reconstruct every reported centroid and variance, and sum to the
cell marginal. The original cell-only output omission is resolved. No reserved acceptance configuration is used.

The default n=32 output has probability norm error below 4e-15 and a reported
maximum edge-shell probability 0.0004100120. Its maximum arithmetic-centroid
distance from the weighted group-velocity ray is 0.0011758164. These finite
periodic effects are measured, without a certified infinite-volume bound.
The central-momentum velocity differs from the weighted velocity, especially
in z (.3748876092 versus .2467860627), so the weighted comparison is necessary
for this packet width. The root's reference flat metric and ray are explicitly
supplied interpretations, without derived gravity, a hidden particle trajectory
or a modular-to-proper-time conversion.

Final source guards reject invalid masses, momenta, grids, widths and time samples.
A nonzero first sample verifies the relative-time ray anchor. The producer states
a finite-grid phase convention rather than a global continuous-gauge theorem.
The retained source and final default artifacts are frozen for this review.
