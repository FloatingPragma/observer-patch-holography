# Contact interaction, particle pole and complete source readout

This CP-06 component uses the same supplied two-component Weyl/contact
Hamiltonian as [the complete Fock field](../full_fock_field/README.md):

```
H = dGamma(h) + U sum_x n_(x,up) n_(x,down),
h(k) = -c sum_a sigma_a sin(k_a),  c,U >= 0.
```

It connects the existing isolated pair to the complete two-fermion
interaction response, constructs a controlled fixed-broadening volume
limit, and executes a finite-time interaction read with an untouched
reference. The parameter t is Hamiltonian model time; z=E+i eta is its
complex frequency and eta>0 is a declared broadening. The action, lattice
spacing, input preparation and physical clock dictionary are supplied.
The subsequent [on-shell theorem](ONSHELL.md) constructs complete
two-fermion lattice wave operators and a unitary scattering matrix. An
exact rational certificate bounds its eight-node contact threshold; the
interacting density reaches zero energy. No particle species, physical mass,
chiral spatial continuum or physical cross section is assigned.

## Complete contact response

At total momentum K, the four-spin relative-momentum operator is

```
H0,K(k) = h(k) tensor I + I tensor h(K-k),
H_K = H0,K + U |d><d|,
d(k) = singlet / sqrt(L^3).
```

The rank-one method is the standard resolvent construction described in
[rank-one perturbation theory](https://arxiv.org/abs/math/9601206); the explicit
Weyl/contact and source-operation maps are the attachment established here.

The finite physical subspace satisfies `psi(K-k)=-SWAP psi(k)`. Both H and d
preserve that subspace. All relative momenta and four spin entries are
retained before this exact fermion restriction. The finite untwisted torus
requires reciprocal-grid K. General K in the quadrature routine denotes a
Brillouin integral approximation, and is labeled accordingly.

For z in the upper half plane define

```
R0(z) = (z-H0)^(-1),   g0(z) = <d,R0(z)d>,
tau(z) = U / (1-U g0(z)),   T(z) = tau(z) |d><d|.
```

The rank-one resolvent identity gives, on every state/reference column,

```
R(z) = R0(z) + R0(z) T(z) R0(z),
gU(z) = <d,R(z)d> = g0(z)/(1-U g0(z)).
```

The spectral theorem gives `-Im g0=eta ||R0 d||^2>0`, so the denominator
cannot vanish for eta>0 and U>0. The same identity applies to gU. Direct
substitution then gives the optical identity

```
-Im tau = |tau|^2 (-Im g0),
T-T* = T* (R0-R0*) T.
```

The density `-Im gU/pi` is the complete contact spectral measure convolved
with a Lorentz kernel of width eta. T and R are response operators, not
trace-preserving instruments or normalized probabilities. The code
compares them with an independent full-matrix inverse, including the
complete physical exchange projector and arbitrary reference columns.
Zero interaction gives T=0. Contact-orthogonal vectors remain in the
carrier even when T annihilates them; that kernel is not asserted to be a
protected subspace of the time evolution.

The existing [bound-pair pole](../bound_pair/README.md) obeys exactly the
same denominator equation `1-U g0(E_B)=0` outside the free band. Its
scattering residue is `Res_(E_B) tau=U^2 Z_B`, where Z_B is the previously
constructed contact weight. An independent complete physical eigensystem
checks that relation. The isolated pair and this interaction response
therefore use one Hamiltonian, contact vector and pole. The routines also
admit weak U outside the separate strong-coupling bound-pair hypothesis;
nonzero interaction alone does not assert an isolated particle pole.

## Exponential volume control at fixed broadening

The infinite relative-coordinate fiber is the multiplication operator
H0,K on `L2(T^3,C^4)` plus the rank-one contact. Use normalized Brillouin
measure and constant normalized singlet d. This bounded self-adjoint
operator defines g0 and gU at every eta>0. The construction keeps lattice
spacing fixed.

For a complex momentum shift with `|Im k_a|<=a`,

```
||H0,K(k+i y)-H0,K(k)|| <= 6c(exp(a)-1).
```

Choose `a=log(1+eta/(12c))` for c>0. The change is at most eta/2, so the
resolvent exists throughout this closed strip and has norm at most 2/eta.
The Fourier argument is the periodic analytic-strip trapezoidal method
of [Trefethen and Weideman](https://people.maths.ox.ac.uk/trefethen/publication/PDF/2014_149.pdf),
with the three-dimensional alias sum made explicit here. Its scalar Fourier coefficients obey
`|g_hat(n)| <= (2/eta) exp(-a |n|_1)`. The L-grid quadrature aliases only
Fourier indices in `L Z^3`; absolute convergence gives

```
q = (12c/(12c+eta))^L,
epsilon_g = (2/eta) ( ((1+q)/(1-q))^3 - 1 ),
|g0,L-g0,infinity| <= epsilon_g.
```

For c=0, the integrand is constant and the error is zero. All displayed
bounds are exact rational arithmetic when c,U,eta are rational. The
bounded calculator permits L through 4096; the argument holds for every
positive finite L. Streaming numerical quadrature permits L through 256
without allocating a global many-body state.

The identity `tau=U+U^2 gU` and `|gU|<=1/eta` give
`|tau|<=M=U+U^2/eta`. Therefore

```
|tau_L-tau_infinity| <= M^2 epsilon_g,
|gU,L-gU,infinity| <= (1+U/eta)^2 epsilon_g.
```

These bounds are uniform in real E and total K. At c=1/8, U=1, eta=1/4,
L=128, the exact rational bounds give
`epsilon_g < 1.295e-7` and `epsilon_tau < 3.236e-6`.
For the illustrative z=1/8+i/4, the evaluated amplitude is approximately
`0.2169001462 - 0.4165668493 i`. Floating evaluation error is separate from
the analytic volume bound. No interchange with eta tending to zero,
on-shell S matrix, asymptotic completeness, physical detector bandwidth or
spatial continuum limit is claimed.

## A complete source-operated interaction experiment

The local read uses two sites joined along a declared coordinate axis,
with all four field modes and one initially blank pointer. A supplied
normalized preparation has two separated fermions, entangled with a
reference qubit:

```
(|modes 0,3> |0>_R + i |modes 0,2> |1>_R) / sqrt(2).
```

Both terms initially occupy mode zero. Execute the actual source words for
`exp(-it H_U)`, then the controlled free reverse `exp(+it H_0)`, and the
occupation pointer for mode zero. The executor uses every complete
six-dimensional source primitive and its sender-reset export instrument.
Every first-failure Kraus group and both pointer outcomes are retained;
there is no success normalization or conditioning. This is a finite-time
interaction echo, not a source preparation of asymptotic incoming particles.
Initialization, reference storage and the reverse-free control are supplied.
This local echo does not prepare asymptotic scattering states or implement
a detector for the infinite-volume T matrix.

On the z-bond use the normalized sums
`S=(|0,3>+|1,2>)/sqrt(2)` and `D=(|0,1>+|2,3>)/sqrt(2)`. Their block is
`[[0,ic],[-ic,U]]`; the orthogonal separated combination has energy zero,
and the other prepared spin channel `|0,2>` is blocked. Let
`Omega=sqrt(c^2+U^2/4)`. Its exact empty-pointer probability is

```
p_empty = (1-Re alpha)/4,
alpha = exp(-iUt/2) [ cos(ct) cos(Omega t)
        + (c/Omega) sin(ct) sin(Omega t)
        + i (U/(2 Omega)) cos(ct) sin(Omega t) ].
```

At c=1/8, U=1, t=4, the rational certificate brackets sqrt(17), evaluates
sin/cos Taylor polynomials through degree 39 with remainder bounded by
`|x|^40/40!`, and bounds the dependence on the root by a Lipschitz constant
less than one for root>=4. Specifically, writing r for the root gives
`|d cos(r/2)/dr|<=1/2` and
`|d(sin(r/2)/r)/dr|<=1/(2r)+1/r^2<=3/16`, so
`|dp_empty/dr|<=23/64`. Each remaining product has at most three
trigonometric factors, whose polynomial approximations have magnitude
less than 1.01. Product errors are therefore bounded conservatively by
100 times the largest Taylor remainder. Outward rational rounding gives

```
0.046942 <= p_empty <= 0.046943.
```

The complete 32-step source schedule has ideal half-diamond error at most
`11/1024`, using the existing overlapping-support product-formula bound;
the free reverse is exact on this disjoint-hopping bond in ideal arithmetic.
Thus the implemented ideal-word empty-pointer probability exceeds
`579197/16000000 = 0.0361998125`. The U=0 control has exact empty probability
zero, with the same initialization and requested source structure.

Literal numerical source execution gives approximately 0.04691534 and
0.95308466 for the two pointer outcomes, with retained abort probability
about `5.7e-29` and reference Gram error about `1.2e-12`. It requests
**22,531 pulses and 45,062 load/export service events**. Numerical admission
has its declared `2e-8` tolerance; the rational ideal-word bound does not
include floating or hardware error. Independent Jordan--Wigner/exponential
controls check all three spatial axes. Additional input controls populate
every field number sector and retain their coherences and references.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.contact_scattering
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/contact_scattering
```
