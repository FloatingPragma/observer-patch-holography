# Same-source curved-ball and bulk calculation

This extends the frozen parent candidate beyond a uniform area dilation.
The graph and regional subsystem remain fixed. The local coframe tangent has
zero regional mean, so the proper volume is fixed to first order. The actual
filled Gaussian state's regional entropy is differentiated; no area weight
is applied to the entropy afterwards.

The complete local Hermitian source is

    H[v] = sum_j alpha_j {v,D_j}/2 + beta R^3 (v L v)^2,
    v=exp(-epsilon sigma),

where L is the positive lattice Laplacian. Its tangent is

    delta H = -sum_j alpha_j {sigma,D_j}/2
              -beta R^3 (L^2 sigma+2 L sigma L+sigma L^2).

`dense_control.py` constructs this full 256-dimensional Hamiltonian,
diagonalizes its actual perturbed ground projectors, and compares the regional
entropy difference to `run.py`'s Bloch response. The difference is
`3.14e-9` at step `1e-5`; the independent operator tangent also passes.
Uniform scaling of the complete Hamiltonian is retained as a null control.

For the initial periodic perturbation

    sigma(x)=sum_j [mean_A cos(q x_j)-cos(q x_j)],

`sum_A sigma=0` exactly algebraically. At small q this is the radial profile
`q^2(r^2-mean_A r^2)/2`. A separate post-diagnostic calculation uses the
exact signed-coordinate quadratic profile, periodically continued beyond the
region, to remove the cosine's local angular terms. Its Fourier components
are summed without discarding the high modes. The small residual difference
between lattice sums and smooth proper volumes remains explicit.

For a smooth sphere of volume `4 pi ell^3/3`, the cosine profile's first-order
area change at fixed volume is

    delta A = 24 pi ell^2 [f_volume(q ell)-sinc(q ell)],
    f_volume(z)=3(sin z-z cos z)/z^3.

The discrete-volume comparator substitutes the actual regional cosine mean.
For the radial quadratic, the smooth result is `16 pi ell^4/5`; its
discrete-volume counterpart is `8 pi ell^2(ell^2-mean_A r^2)`.
The planar eta is measured independently from the same source, never fitted
to the ball result.

At L=24, R=1:

| Ball | Cells | Actual curved entropy response | Smooth geometric ratio | Discrete-volume geometric ratio | Same-ball bulk/continuum ratio |
|---|---:|---:|---:|---:|---:|
| radius 3 | 123 | 2.1053083 | 1.04270 | 1.05871 | 0.76878 |
| radius 4 | 257 | 5.4295478 | 1.03734 | 1.03150 | 0.87747 |

The same-ball bulk entry is the actual modular energy of the lowest complete
particle/hole pair, divided by `2 pi C_ell`, with `C_ell=4 pi ell^4/15`.
It has no separately fitted temperature. Complete projector and four spinor
components are retained.

These results supply a correct-sign curvature response on a declared fixed
graph source with a bulk response approaching its required value. They do
not supply a convergence bound, generalized-entropy stationarity, a dynamical
geometry, ten observers, or the physical value of G. In particular, no
Poisson equation is said to follow merely from the positive response. If
stationarity were supplied, the measured ratios would determine its finite
error; imposing stationarity remains an additional assumption.

The original F087 source and its failed handles are unchanged. This is a
different candidate, with a proper-regulator rule which itself needs a native
OPH derivation. `RESULT_*.json` retain all numerical outcomes, including
failures of any proposed area or bulk equality. `SUMMARY.json` consolidates
the completed runs and normalization sensitivity when available.
