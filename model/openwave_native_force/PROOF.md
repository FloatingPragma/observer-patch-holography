# Native charges in a finite elastic geometry

## Local source and declared spatial extension

The actual three-oscillator, cutoff-two reciprocal scalar action is
`H(r)=sqrt(14) M(r)`, where

```
M(r)=a(r)D-C/(7r),
a(r)=r/4+1/(7r)+3r²/28,
D_mask,mask=2 popcount(mask)+3,    C=XXI+IXX.             (1)
```

The producer rebuilds this matrix from the original
`native_equilibration.vacuum_repair` coefficient expressions and compressed
oscillator matrices. No gravitational equation enters (1). The common
factor sqrt(14) is used as the energy unit below. At reference geometry
`r=1`,

```
M0=D/2-C/7,     J=dM/dr=9D/28+C/7.                      (2)
```

The dimensionless geometry charge of a specified prepared-state contrast
delta rho is `j=Tr(delta rho J)`. The contrast between `|001>` and `|000>`
has excess energy one and charge `9/14`. Both numbers follow from the
same source matrices. The prepared sources are externally pinned states;
they are not identified freely evolving particles or stationary vacua.

Declare a cubic spatial extension with geometry `r_x=1+q_x`, zero
Dirichlet boundary values, and positive nearest-neighbour elastic energy

```
E_lin(q)= (kappa/2) sum_edges(q_x-q_y)² + sum_x j_x q_x.
```

The spatial gluing, stiffness kappa, absence of an onsite geometry mass,
boundary condition and minimum-energy preparation are architecture inputs.
They are not consequences of the two-valued reciprocal thermal update.
The local matter term is the expectation of (1) relative to its displayed
reference preparation, with its r=1 constant removed. A zero-reference
subtraction and fixed source preparations are explicit choices of this
static-response experiment. No inference about the source's autonomous
vacuum or a full Einstein limit follows.

Let `L` be the positive Dirichlet graph Laplacian. Completing the square
derives, rather than independently postulating,

```
q_star=-L^-1 j/kappa,
E_lin,min=-j^T L^-1 j/(2kappa).                         (3)
```

The matrix L has diagonal six and nearest-neighbour entries minus one.
Its inverse is positive, by the discrete maximum principle. Subtracting
the two one-source minimum energies from the pair minimum therefore gives

```
V_int(x,y)=-j_x j_y G_D(x,y)/kappa.                     (4)
```

Equal-sign charges have negative interaction energy. The executed central
differences at separations 4, 6 and 8 on 31, 63 and 95 point boxes give
attraction and the long-range profile `j²/(4pi kappa r²)` with stated
finite-box and lattice errors. This is the mutual interaction after
subtracting boundary-dependent self energies; it is not the total force
including each source's own boundary response. The sine transform in the
code inverts the specified quadratic energy, without inserting a `1/r`
array. Doubling kappa halves the response. Adding an allowed positive
onsite geometry term screens it. Neither strength nor infinite range is
selected by the local matter Hamiltonian alone.

## A completely finite geometry carrier

At every site q is restricted to the dyadic grid of spacing `2^-48` within
`[-1/50000,1/50000]`. This gives a finite Cartesian carrier. The spatial
source uses mixture strength `1/10000`, hence `j=9/140000` at each occupied
source site. The finite grid, source mixture and box are declared controls,
not fitted physical scales.

The unconstrained minimizer in (3) is used to find a dyadic trial q.
Every trial is checked to lie strictly inside the finite geometry window.
Its energy and residual `b=Lq+j` are evaluated with exact integers and
rational arithmetic. Since

```
lambda_min(L)=12 sin²(pi/(2(n+1))) >=12/(n+1)²,
E_lin(q)-||b||²/(2 lambda_min_lower) <= E_lin,min <= E_lin(q), (5)
```

the same interval encloses the finite-grid minimum: the unrestricted
minimum is a lower bound and the retained dyadic q supplies an upper bound.
No claim of having found the exact discrete minimizer is needed. The
residual norm is bounded by `n³ ||b||_infinity²`, using an integer overflow
guard. The proof of (5) is the exact positive quadratic identity
`E(q)-E(q_star)=b^T L^-1 b/2`.

The nonlinear native matter term can also be retained without minimizing
an enlarged nonlinear carrier. For this occupation contrast the exact
remainder is

```
2a(1+q)-1-(9/14)q = q²(7+3q)/(14(1+q)).               (6)
```

It is positive throughout the declared window. Its coefficient is at most
`B=(7-3h)/(14(1-h))`, where `h=1/50000`. Each source therefore adds a
remainder between zero and `mixture*B*h²` on the ENTIRE finite carrier.
The full nonlinear minimum for a configuration with k sources lies between
the linear lower endpoint and the linear upper endpoint plus
`k*mixture*B*h²`. For a pair interaction, the two single-source subtractions
give an absolute cross-energy error at most `2*mixture*B*h²`.

At separation six on the 31-point box the full nonlinear mutual-force
central difference lies in

```
9.6427225513e-12 < attractive magnitude < 9.7227234669e-12. (7)
```

The underlying endpoints are exact rationals in the receipt. This bounds
the full finite-carrier native matter energy, rather than applying a
linear Taylor approximation at one chosen minimizer. The minima refer to
the declared static energy experiment; a native controller preparing them,
a physical length unit and a dynamical spacetime interpretation are not
constructed.

## A universal positive-mass interpretation fails

The same matrices (2) give exact contrasts relative to `|000>`:

| Prepared change | Excess energy | Geometry charge | Charge/energy |
|---|---:|---:|---:|
| occupation 001 | 1 | 9/14 | 9/14 |
| normalized sum of occupations 000 and 011 | 6/7 | 11/14 | 11/12 |
| normalized difference of occupations 000 and 011 | 8/7 | 1/2 | 7/16 |
| normalized 000 minus one tenth of 011 | 34/707 | -11/707 | -11/34 |

All four excess energies are positive, while charge/energy varies and
one charge has the opposite sign. This prevents one universal
strength from identifying this force with positive mass on that declared
preparation family. The reference `|000>` is not the source ground; that
limitation is not hidden.

A separate certified counterexample uses the ACTUAL ground. Exact rational
inertia isolates all eight eigenvalues of M0, including its simple lowest
level. A binary rational trial ground v has Rayleigh value E_v and
infidelity bounded by

```
1-|<g,v>|² <= (E_v-E0_lower)/(E1_lower-E0_upper).
```

The norm bound `||J||<=89/28` bounds its charge-expectation error by
`2*(89/28)*sqrt(infidelity)`. Choose the two binary rational vectors near
`v +/- (1/100) (J-<J>)v/||(J-<J>)v||`. Their exact rational Rayleigh
energies lie strictly above the upper bound on the true ground energy,
while their charges relative to the true ground have strictly opposite
signs, even after the full reference-error bound. The receipt retains the
vectors and every interval. Floating eigensolvers propose the vectors;
rational spectral and Rayleigh bounds certify the claim.

Thus the sign obstruction persists with a true ground reference. It rules
out universal attractive positive-energy gravity for THIS scalar coupling
and its admitted state family. It does not exclude other geometry sources,
restricted physically selected particles, or a spin-two Einstein
construction. The concrete requirements exposed here are a justified
spatial geometry law, a physical stiffness/scale, range protection, and
universal source coupling on the actual particle sector.
