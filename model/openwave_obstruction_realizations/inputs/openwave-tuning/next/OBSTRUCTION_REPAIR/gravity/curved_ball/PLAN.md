# Fixed-volume curved-ball development test

Use the declared quartic source from the frozen parent packet, with the local
Hermitian coframe rule `H=sum alpha_j{v,D_j}/2+beta R^3(v L v)^2`.
At v=1 the geometry tangent is
`V_sigma(k+q,k)=-(A(k+q)+A(k))/2-beta R^3(lambda(k+q)+lambda(k))^2`.

The real periodic profile is `sigma(x)=sum_j[mean_A cos(q x_j)-cos(q x_j)]`.
Its mean on A is exactly zero, so the discrete proper volume
`sum_(x in A) exp(3 epsilon sigma(x))` has zero derivative at epsilon=0.
For small q it is the radial profile `q^2(r^2-mean_A r^2)/2`.
No subsystem is added or removed; in the radial continuum metric its boundary
is a geodesic sphere. The finite-q profile is not exactly radial, and the voxel
ball is not a smooth sphere; these errors remain declared.

Compute actual Gaussian regional entropy variation from the filled projector's
Frechet derivative. Compare to the independently measured planar area
coefficient and the continuum fixed-volume sphere integral. Measure the
matter modular coefficient on the same ball and source. Initial configurations:
R=1, L=24 with ball radius 3 and 4; R=1, L=40 with radius4. These are development
calculations, not prefiled acceptance. Retain nonmatches. A nonzero response
alone does not count as the Einstein area identity or stationarity.

Controls: uniform full-H scaling must give zero projector response;
the coframe tangent is checked against a direct finite matrix derivative on
a small lattice; the fixed-coordinate comparison is preserved. Numerical
entropy-eigenvalue floors must be reported. No native law, Newton strength,
second metric potential, or all-observer conclusion is inferred.
