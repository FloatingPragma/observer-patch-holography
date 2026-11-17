# Four-cell interacting wall source and bounded charged response

## Carrier, references and resources

The source has four spatial cells, two regulator slices, one physical
spin and sixteen complete CAR modes. It retains every number sector
N=0 through16. All spatial momenta, both energy signs and every wall
bulk coordinate belong to the operator. Its Wilson coefficient is1,
its cell clock is h=1/8 and its electric coefficient is1/4. The bare
mass is either zero or the actual432 electron coefficient. Each mass
case has its own source state. There is no multi-family vacuum claim.

The compact one-dimensional gauge field has one holonomy after a tree
gauge choice. At vertex x the reference is −2, with an additional
charge −(N−8) at vertex zero. For every matter occupation w and integer
ell, define

```
q_x(w)=sum_(v<=x)n_v(w)+2(3−x)−N,
E_x=ell+q_x(w), H_electric=sum_x E_x^2/4.
```

The last q is zero. These fields solve every local Gauss constraint.
Interior hops update their occupation-dependent fluxes; the closing
hop also increments ell. Local charged probes at vertex zero change
the reference, so their complete action preserves Gauss law. The
reference has zero declared energy. This U1 loop control has no proved
reducing-subspace transfer from the full nonabelian U3 gauge vacuum.

The full source is the direct sum of finite CAR spaces tensored with
the full rotor. Electric energy has compact resolvent, and the complete
matter multiplication operator is bounded. The self-adjoint Gauss
Hamiltonian therefore has an attained ground energy. At flux cutoffK4
the full direct sum has589824 states; its largest N8 sector has115830.
The matrix-free action retains all of them. The declared numerical
ceilings are16 CAR modes, K<=4 and250000 states per individual sector.
No dense-spectrum requirement is silently replaced by a reduced band.

## Infinite-rotor lower bounds for every sector

At a fixed holonomy theta the matter operator has exactly four spatial
momenta, `k_j=(theta+2pi j)/4`, j=0,1,2,3. Put `W=1−cos(k_j)`. For
two wall slices the complete positive bands in cell units satisfy

```
e_±^2=2W+1/2 ± sqrt(1+4W^2)/2.
```

Both negative bands are retained. Multiplication by8 converts to the
declared source-energy clock. The fixed-N matter lower edge V_N(theta)
is the sum of the N lowest entries in this full sixteen-state spectrum.
Because electric energy is positive, `H_N>=min_theta V_N(theta)`.
This inequality applies to the actual fluctuating infinite rotor; it
does not assume a classical holonomy or a source-selected free vacuum.

The closing-hop operator has fixed-N norm at most
`c_N=8 sqrt(floor((n+1)^2/4))`, n=min(2,N,16−N). Its two rank-one
channels form two raising operators, with other modes retained as
spectators. Therefore V_N is Lipschitz with constant2c_N. The executable
uses1024 holonomy cells and interval arithmetic at40 decimal digits.
Every momentum and band is evaluated at each grid point, and the
global correction2c_N pi/1024 is subtracted. Endpoint conversion to
binary floating point is outward. All seventeen sector bounds are
stored, including sectors excluded by their lower edge.

There is also a lower bound that retains electric zero-point energy.
Write `qbar=sum q_x/4`. Exactly,

```
sum_x(E+q_x)^2/4=(E+qbar)^2+sum_x(q_x−qbar)^2/4.
```

For a vector-valued wavefunction Psi(theta),
`|d ||Psi||/dtheta|<=||(d/dtheta+i qbar)Psi||`. The remaining variance
term is positive. Thus the complete source is bounded below by the
scalar periodic operator `−d^2/dtheta^2+V_N(theta)`, acting on the
wavefunction norm. This vector diamagnetic inequality uses the actual
Gauss connection, including all occupation charges.

For any positive periodic f and real periodic r, integration by parts
gives

```
int (|r'|^2+V r^2)
 =int [V−f''/f]r^2 + int f^2 |(r/f)'|^2.
```

The retained verification functions are `f=exp(a cos theta)`, with
a=0, ±1/4, ±1/2, ±1, ±2. They alter no Hamiltonian coefficient.
Every resulting Barta bound is retained. Its sampled function is
`V_N+a cos(theta)−a^2 sin^2(theta)`, with global Lipschitz constant
at most `2c_N+|a|+a^2`. The interval grid therefore supplies a complete
lower bound for each test. Taking the largest verified lower bound
is a verification operation, not a physical parameter fit.

At zero bare mass the bounds include

```
H_7,H_9 >= −89.996420264674,
H_8 >= −93.209400822772.
```

All other number sectors have higher lower bounds. A boundary input
mass changes each local CAR contact by norm at most |m|, hence the
full four-cell shift is at most4|m|. That bound transfers the sector
inequalities to the declared electron-input case. Empty and full
sectors have zero contact shift.

## What the variational state certifies

Every preparation is a normalized vector on the complete finite N8
occupation/flux space, embedded by zero extension into the infinite
rotor. Its full complex amplitudes, all probabilities, occupation labels
and Gauss-prefix charges are retained in the preparation NPZ files.
The actual infinite-source Rayleigh expectation equals its finite
expectation because exterior hops are orthogonal to that preparation.
The energy upper bound includes the stated numerical margin.

For m=0, K4, the variational energy is approximately−92.7441895061.
It lies below every competing number-sector lower bound. Therefore
the actual full-rotor vacuum number is N8, and every N7/N9 state lies
at least2.7477 source-energy units above that variational reference.
Since the true ground energy is no larger than the variational one,
this also gives a positive charged-sector separation above the true
finite-volume vacuum. It is not an observed lepton mass.

The neutral energy excess of the preparation is bounded by its upper
Rayleigh value minus the N8 Barta lower bound, about0.4653 at K4. Sparse
neutral eigenvectors do not establish neutral spectral order or exact
ground-state fidelity. Those acceptance fields are false. A warm-seed
control returned inconsistent second Ritz levels across nearby masses;
its complete results and vectors are preserved under
`warm_seed_control/`. The comparison run adds an explicit generic
component1/16 to each warm seed. That improves the numerical search
without changing the source or turning its neutral order into a theorem.

The K4 preparation has `<cos(theta)>` approximately−0.657. Unlike the
one-cell source, every holonomy has a lower free momentum satisfying
`W_low<=1−cos(pi/4)≈0.2929`. The measured event `W_low<=1/4` has
probability about0.784. It is computed by analytic Fourier integration
of the complete rotor density matrix. This is a localization-hypothesis
diagnostic, not a decomposition of the interacting state into free
particles. The source-selected neutral wavefunction and its charged
responses remain interacting objects.

## A controlled complete charged resolvent

Use both local boundary Dirac probes at vertex zero. Their coefficient
isometry is `E=|wall0>P+ +|wall1>P−`, and `E*E=I2`. The particle and
hole currents act on all occupations and preserve Gauss law through
the reference translation. Their full CAR sum is checked, including
both energy signs. Exact particle-hole conjugation, with Dirac X and
ell reversal, maps the complete N9 operator to N7. Its canonical CAR
signs transport the particle current; no equality of particle/hole
residues is assumed.

Let Eref be the declared variational reference, b one of these complete
currents, and v=b psi or b* psi. The sector bounds prove
`H_charged−Eref>=g I`, with g>0. If a finite vector x solves its truncated
linear system approximately, its full infinite-source residual consists
of the stored interior residual and exactly two exterior flux blocks:

```
r_in=(H_K−Eref)x−v,
r_out,+=C_N x_K, r_out,−=C_N* x_−K,
||(H_infinite−Eref)^−1 v−x|| <= sqrt(||r_in||^2+||r_out||^2)/g.
```

This is an operator-norm resolvent statement on the entire infinite
rotor. It requires neither a spectral-pole truncation nor a guess about
unrecorded high-energy modes. Every interior coefficient and both
exterior residuals are retained in the resolvent NPZ files. The declared
floating-point residual margin is included in the stored upper bound.

The two particle and two hole solutions give the complete local
prepared-state Green matrix. Cauchy-Schwarz transfers their solution
errors to a matrix-norm error. At m=0, K4, the full infinite-rotor Green
error is bounded by approximately4.12e−5. The K2 bound is about0.052,
which fails the attempted inverse-response admission. That failed
control is retained. At K4 the inverse perturbation bound is finite:
if `||G_K^−1|| epsilon<1`, then

```
||G_infinite^−1−G_K^−1||
 <= ||G_K^−1||^2 epsilon/(1−||G_K^−1|| epsilon).
```

Its K4 bound is about0.074. The beta coefficient of the local inverse
response is a specified prepared-state observable; it is not identified
with a propagating mass. Spatially integrated local Green functions,
finite charge references and residue normalization prevent that
identification without an additional pole analysis.

The reference energy is variational and the state is an explicit
preparation. These bounded responses are distinguished from the exact
neutral-vacuum propagator. The neutral-state fidelity theorem and its
pole/current transfer are absent. The two bare-mass cases do not form
a common interacting three-family source, and no muon-ratio repair is
claimed.

## Numerical and physical scope

The interval band and Barta bounds use `mpmath.iv` with outward endpoint
conversion and analytic global Lipschitz estimates. Sparse vectors and
matrix products use double precision, with source-energy upper margin
1e−6, full-residual margin1e−8 and Green-product margin1e−8. The finite
source has bounded row complexity, source-energy row bounds and the
declared dimension ceiling. These are conservative standard-IEEE
admissions, not an interval sparse-linear-algebra claim. Independent
review reconstructs the operator from its literal coefficients and
checks the complete retained vectors, currents and residuals without
importing the producer source module.

This construction supplies a full finite-volume Gauss source, a proved
vacuum number sector and a bounded complete charged response of explicit
preparations. It does not establish a three-dimensional continuum,
physical charged-family poles, neutral-vacuum fidelity, a dynamical
coframe join or a full native execution. The free spatial-refinement
result in `wall_ring_control` preserves its own scope. These interacting
finite-volume bounds neither replace that result nor exclude its
infrared limit.
