# Wavepacket on an explicitly stacked ring graph

This experiment propagates a probability packet on a declared three-dimensional
periodic graph. The original canonical ring screen is two-dimensional. The
third direction is supplied by stacking copies of that screen and adding one
vertical seam per port. Its existence, spacing and unit hopping weight are
construction choices. The calculation does not derive three dimensions,
spacetime geometry, gravity, or a hidden particle trajectory.

Each periodic cell contains two twelve-port carriers, 63 horizontal seams and
24 vertical seams. Ports and seam-transit states give 111 internal components.
Carrier A is at the cell origin and carrier B at
D=(sqrt(3)/2,1/2,0). The lattice vectors are

```text
a1=(sqrt(3),0,0), a2=(sqrt(3)/2,3/2,0), a3=(0,0,1).
```

The Fourier variables q are cell covectors. With the inverse transform
exp(+i q.n), horizontal glue from A(R) to B(R-t) has the appropriate endpoint
phase exp(+i q.t). Vertical seams from R to R+a3 give incidence block
`B_z=(1-exp(-i q3)) I_24`. Thus the port Laplacian is exactly

```text
L_3(q)=L_ring(q1,q2)+2(1-cos(q3)) I_24.
```

Deleting vertical seams removes all q3 dependence. Independent tests verify
the missing-z control, including its stationary z probability marginal.

The coherent generator is the supplied port/transit form
`H_m=[[m,B],[B^dagger,-m]]`, with m=0.23. If u is the normalized lowest
port-band eigenvector, lambda its Laplacian eigenvalue and
E=sqrt(m²+lambda), the positive-band spinor is

```text
P=sqrt((E+m)/(2E)) u,
Q=B^dagger u/sqrt(2E(E+m)).
```

The code checks the full eigen-equation and normalization. Only the positive
acoustic branch is populated; other positive, negative and flat branches exist
in the full theory. Zero commit rate is selected. A local native commit can
mix bands, so this preparation is not an invariant sector of a nonzero-rate
experiment by assertion alone.

Registers have modulus 6 and the preparation uses uniform modular Gauss
fibers. The complete local-coarse reduction theorem applies to the connected
augmented graph. The huge record fibers are represented implicitly, without
an executable preparation procedure for their global coherent superposition.
This is distinct from definite-zero-record preparation. No record energy,
magnetic field, interaction, or empirical mass/length/time anchor is selected.

A Gaussian envelope is specified in cell momentum, centered at
q0=(0.75,0.45,0.10), with sigma=0.20. A finite-grid component phase convention
fixes the Bloch eigenvector gauge: port component 0 is made real positive,
with a declared fallback pivot if it is too small. This is a preparation
choice, without a global smooth-gauge theorem. The packet is centered in a
32-by-32-by-32 periodic cell array. Evolution multiplies each prepared
momentum spinor by exp(-i E t), followed by an orthonormal inverse FFT.

The position apparatus places internal ports and carrier-internal seams at
their carrier center and glue/vertical transits at whole-seam midpoints.
There are seven distinct offsets per cell. The saved `spatial_probability`
array retains every probability for this named coarse geometric apparatus;
colocated internal labels are explicitly aggregated. The 111-to-seven map
and offsets are retained. A separate cell marginal is used for the cloud
plot. Independent tests reconstruct all reported centroids and variances from
these saved bins.

Each displayed time uses a separately and identically prepared ensemble.
The figure does not depict sequential position measurements on one load.
The centroid is the first moment of the full declared position distribution,
not an observed point-particle path.

Let A have the lattice vectors as columns. Physical coordinates and wavevectors
satisfy r=A n+offset and q=A^T k. Therefore the physical group velocity is
`A gradient_q E`, and the quadratic form transforms as
`Hessian_k(E²)=A Hessian_q(E²) A^T`. Independent reconstruction gives

```text
E²=m²+(15/352)(kx²+ky²)+kz²+high-order lattice terms.
```

The horizontal acoustic speed squared is 15/352 and the vertical value is 1.
This anisotropy is retained. No common isotropic relativistic light cone is
established. The width of the momentum packet matters: its centroid is
compared with the envelope-weighted group velocity. The central-q ray is
reported separately and differs appreciably in z. The reference ray is
anchored at the first sampled time and position.

The reference metric is flat and is supplied as a long-wavelength
interpretation of this symbol. Its straight constant-velocity ray is therefore
a geodesic of that reference metric. Agreement of a finite wavepacket centroid
with it does not derive curved spacetime, a gravitational force, or a modular
to proper-time conversion. The native evolution parameter has inverse supplied
hop-amplitude units, without an empirical physical-clock calibration.

Periodic wrap is monitored with an edge-shell probability. Arithmetic
centroids can change when probability crosses a torus boundary. The reported
centroid-ray error and edge weight describe this finite run; no certified
infinite-volume error bound is supplied. The code rejects invalid positive-band
masses, grid sizes, momentum vectors, packet widths and time samples.

The two-dimensional force-kernel results do not transfer automatically to this
three-dimensional graph. Changing spectral dimension changes Green-function
asymptotics, and a proposed force additionally needs a source coupling,
interaction operator, state and readout. In particular a two-dimensional
inverse-square-root kernel cannot be called the same three-dimensional static
force after adding z links without a separate derivation. For a low-momentum
quadratic Laplacian, Fourier scaling gives an inverse-square-root kernel
proportional to r^(1-d): 1/r in two dimensions and 1/r² in three. An ordinary
inverse Laplacian has a 1/r kernel in three dimensions. These conditional
kernel scalings alone do not define an interaction or a measured force.

Files: `stacked_wavepacket.py`, `STACKED_WAVEPACKET.json`,
`position_distributions.npz`, PNG/SVG plots, and independent reference/tests/
review. The receipt binds source identity and declares the preparation,
geometry, readouts, parameters and limitations. No OpenWave cell or full
physical criterion is awarded by this demonstration.
