# Position and joint crystal momentum on the same finite field

The periodic gauge-coordinate construction has exact conserved spatial
translations and complete position records. The translation moves every
matter and link register. Its physical measurement on the referenced mixed
state also moves the compensating reference. The resulting conserved
crystal momentum belongs to this whole system; it is not the momentum of
an identified particle. All fermion numbers, including vacuum, are retained.

## 1. Reused results and the specified field

Use the [mixed gauge preparation](../particle_local_gauge_state/PROOF.md)
on a cubic torus of side s and V=s³ cells, with sixty native Dirac/internal
CAR modes per cell. At fixed spacing a, polynomial order m, plaquette
coefficient kappa and word cutoff K, the actual field Hamiltonian is

```
H_K=P_K[H_E+dGamma(h_m)-kappa sum_p Re Tr R(U_p)/15]P_K.
```

The full fermionic Fock space is retained. P_K acts on the coordinate
word space; it does not impose a number cutoff. All construction
parameters are the ones declared by that source. The same clock is
`exp(-it H_K)`; no laboratory duration or empirical mass is inserted.

The existing [periodic momentum theorem](../../Lean/QFT/PeriodicMomentum.lean)
transports commutation with a unitary to its Fourier polynomials. Its
Weyl/contact source is distinct; only that finite algebra is reused after
the present source attachment. The existing
[finite state Robertson adapter](../../Lean/EventAlgebra/Robertson.lean)
applies to actual density matrices and Hermitian matrices on a finite
carrier. The [complete instrument contract](../particle_observer_outputs/READOUT_CONTRACT.md)
supplies primitive or compiled realizations with all outcomes and
references. None of these results requires rebuilding its proof or
identifying this field with the separate ordered-spin position limit.

## 2. Translation of every field register

For b in `(Z/sZ)^3`, define T_F(b) by its exterior-power action
`T_F(b)c_(x,A)^dagger T_F(b)^dagger=c_(x+b,A)^dagger`, fixing vacuum.
On an ordered occupied list its sign is the parity of the permutation
needed to sort the translated list. This is an exact group representation,
not an unsigned permutation of Fock basis vectors. The coordinate unitary
relabels every link `(x,mu)` to `(x+b,mu)`. Their product is T_S(b).

All elementary translations commute and `T_S(e_i)^s=I`. Translation
preserves the coordinate-word filtration and its Gram inner product, hence
commutes with P_K. It permutes electric Casimirs and plaquette terms.
The covariant Wilson operator transforms by the same site permutation;
functional calculus gives the same covariance for its polynomial h_m.
Second quantization preserves that identity on the entire Fock space.
Therefore `[T_S(b),H_K]=0` on every retained input. The product
preparation and site-gauge twirl are translation equivariant, so its
source density rho also commutes with T_S(b).

The executable checks all 480 one-body columns on a side-two cubic
carrier with the actual four Clifford and fifteen native representation
rows, at m=0. Nontrivial rational weak links are translated with matter.
Holding those links fixed while translating only matter fails the same
check. Higher polynomial orders, electric terms, all occupations and
word-cutoff covariance follow from the displayed operator construction;
their exponential-size matrices are not materialized.

## 3. Gauge constraints and the translated reference

Translation normalizes the site-gauge action:
`T_S(b)U_S(g)T_S(b)^dagger=U_S(g_b)`, where `g_b(x)=g(x-b)`.
It preserves the source Gauss-zero subspace, but need not commute with
each gauge generator on a charged ambient carrier. The mixed preparation
contains nonsinglet field sectors. Measuring T_S alone while leaving its
compensating reference fixed is therefore not asserted to be an admitted
Gauss-preserving instrument.

Use the source's actual canonical purification `|sqrt(rho)>` with a
conjugate reference. Translation preserves its support. Define

```
T_J(b)=T_S(b) tensor conjugate(T_S(b)),
G_joint=G_source+G_reference,       H_joint=H_K tensor I.
```

The reference Hamiltonian is zero. T_J normalizes the combined gauge
action and preserves its invariant subspace. Both the Hamiltonian and
T_J preserve that subspace. Canonical vectorization gives
`T_J(b)|sqrt(rho)>=|sqrt(rho)>`, because rho commutes with T_S(b).
The same remains true after any duration of the stated Hamiltonian.

The full finite Fourier projections are

```
Pi_k=(1/V) sum_b exp(-2pi i k.b/s) T_J(b),
k in (Z/sZ)^3.                                             (1)
```

Character orthogonality proves `Pi_k Pi_l=delta_(k,l)Pi_k`,
`Pi_k^dagger=Pi_k` and `sum_k Pi_k=I`. Zero-rank labels, when present,
are retained. They commute with H_joint. They preserve the joint
Gauss-zero subspace, and their complements are not discarded. They do
not in general commute with gauge transformations on every charged
ambient input. This is a complete instrument on the full carrier whose
physical preparation domain is stated explicitly. The optional joint
Gauss-zero/outside flag commutes with the projections and can be retained
without accepting only one flag.

Choose real representatives kappa_i(k) of the residues and define
`P_i=(2pi hbar/(as))sum_k kappa_i(k)Pi_k`. These commuting Hermitian
operators are conserved, with `exp(i a P_i/hbar)=T_J(e_i)`.
hbar here is a declared unit conversion. The logarithm branch and lattice
spacing must be recorded; no physical continuum momentum is inferred.
The canonical input has exactly momentum k=0, at every uninterrupted
time. This is total field/reference crystal momentum, not proof that its
field excitations have zero momentum.

## 4. A complete position instrument on all occupations

Let `n_x=sum_A c_(x,A)^dagger c_(x,A)` and `N=sum_x n_x`.
The same Hamiltonian conserves N. Define its inverse to be zero on N=0:

```
F_x=n_x N^(-1) 1_(N>0),       F_vac=1_(N=0).                 (2)
```

Every F_x is positive and gauge invariant. The complete sum is identity:
`sum_x F_x+F_vac=I`. It describes selecting a uniformly random occupied
fermion mode and reporting its site; it does not select one species or
assert a particle position. The square-root Kraus maps give a specified
complete instrument, including the vacuum record. No number sector is
filtered out. Tensoring any untouched reference preserves completeness.

Translation transports these same effects:
`T_J(b)F_x T_J(b)^dagger=F_(x+b)`, while F_vac is fixed.
The initial maximally mixed full Fock law gives

```
p_vac=2^(-60V),      p_x=(1-p_vac)/V.                        (3)
```

The same probabilities hold throughout uninterrupted evolution: N=0
probability is conserved, and translation symmetry makes the V site
probabilities equal. This does not make the full multi-time occupation
law static. The [local-history calculation](../particle_local_history_family/PROOF.md#6-a-nonzero-source-local-dynamical-diagnostic)
derives its nonzero initial squared-increment curvature on the same
order-zero member.

## 5. Complete position-then-momentum records and backaction

The canonical state is translation invariant as a vector, but a retained
position intervention need not leave it in that momentum eigenspace.
Set `A=(1-2^(-60V))/V` and, for two distinct sites, define

```
C=2^(-60V) sum_(i,j=1..60) binom(60,i)binom(60,j)sqrt(ij)
                * sum_(l=0..60(V-2)) binom(60(V-2),l)/(i+j+l).
```

This is the exact full-Fock expectation of `sqrt(n_x n_y)/N`.
It is independent of the link state. The complete initial two-menu law is

```
Prob(vac,k=0)=2^(-60V),        Prob(vac,k!=0)=0,
Prob(x,k=0)=[A+(V-1)C]/V,
Prob(x,k!=0)=(A-C)/V.                                      (4)
```

Each nonzero k has the displayed weight. To derive (4), insert (1)
between the two square-root maps and move T_J through the invariant
input vector. Its identity term gives A; every nonidentity translation
gives C. The sum of a nontrivial character over all translations is zero.
Thus each site's full momentum sum is A, and all records sum to one.

Positivity is explicit:
`A-C=(1/2)E[(sqrt(n_x)-sqrt(n_y))²/N;N>0]`.
The configurations with precisely one of the sixty modes occupied at x
or y give `A-C>=60*2^(-60V)>0`. Every nonzero-momentum record after a
site outcome therefore has positive original weight. There is no
normalization by a selected detector success probability. The expression
for C is a finite computable algebraic sum; the native Fock matrix and
all terms of a large example are not allocated.

At a later time the same character expression holds with
`C_b(t)=Tr[rho(t)sqrt(F_x)sqrt(F_(x+b))]`. The distinct displacements
need not have one common value C. Equation (4) is the exact specified
initial experiment, not a solution of those interacting time correlations.

## 6. Cyclic position, Weyl covariance and finite uncertainty

Assign site coordinates `x_i in {0,...,s-1}` and the vacuum coordinate
zero. The first-moment observable is `X_i=a sum_x x_i F_x`.
With `Q=I-F_vac`, its exact boundary-wrap law is

```
T_J(e_i)^dagger X_i T_J(e_i)
 =X_i+a Q-a s sum_(x_i=s-1)F_x.                             (5)
```

The boundary term cannot be removed on the finite carrier. A distinct
unitary coordinate phase uses the total occupation dipole
`D_i=sum_x x_i n_x`, `Z_i=exp(2pi i D_i/s)`. Integer occupations give

```
T_J(e_i)^dagger Z_i T_J(e_i)=exp(2pi i N/s) Z_i.              (6)
```

Equation (6) is an exact finite Weyl relation. Its phase depends on N;
it is the usual root-of-unity relation on N=1. Z_i is not the exponential
of the random-occupied-site POVM's first moment on every number sector.
Both statements concern the actual translated field/reference carrier.

For any density on this finite carrier, the existing Robertson adapter
applies to the actual Hermitian X_i and P_j:

```
Var(X_i) Var(P_j) >= |Tr rho[X_i,P_j]|²/4.                   (7)
```

Every matrix has the whole finite space as its domain. The variance of
the recorded site coordinate is at least Var(X_i), since its additional
noise is `sum_x (a x_i)² F_x-X_i²>=0`; this follows directly from the
complete position isometry and its positive orthogonal complement.
Consequently (7) also lower-bounds the product using the actual position
record variance and momentum variance. Vacuum records remain included.

Finite trace gives `Tr[X_i,P_j]=0`, excluding `[X_i,P_j]=i hbar I`
for nonzero hbar. The initial canonical state has zero momentum variance
and zero commutator expectation. No universal `hbar/2` lower bound follows.
A physical position/momentum uncertainty relation requires the same
operators' controlled limiting domain, in addition to this finite result.

## 7. Complete finite controls and remaining interpretation

A distinct two-site, one-mode Z2 gauge control checks the reference
requirement on all eight field states and all sixty-four joint states.
Its canonical field density is I/8. Translating matter alone violates a
combined Gauss constraint; simultaneous conjugate-reference translation
preserves the invariant domain. The complete position/momentum rows are
`(site0,0)=1/4`, `(site0,1)=1/8`, the same two values for site1,
`(vac,0)=1/4`, `(vac,1)=0`. Its nontrivial number sector and all zeros
are retained. It is an exact structural control, not an SU(3) source
substitution.

An admitted diagonal coordinate-phase preparation in that control gives
a nonzero commutator readout and verifies (7) with the actual complete
POVM noise. For unit lattice spacing and unit hbar the exact values are
`Var(X)=11/64`, recorded coordinate variance `15/64`,
`Var(P)=3pi^2/16` and `Tr rho[X,P]=i pi/4`. The operator variance
product is `33pi^2/1024`, exceeding the Robertson lower bound
`pi^2/64`; the extra position-record noise is `1/16`.
The full native cubic test separately verifies the source
representation, translation of links and all one-body columns.

The established outputs are a complete position instrument, a conserved
joint crystal-momentum instrument, their exact initial joint records,
cyclic/Weyl covariance and state-specific finite uncertainty. These are
global finite apparatus choices, with compilation errors added through
the owning instrument contract. They do not identify any protected
defect with a particle, remove a charged reference, select a Lorentz
vacuum, derive scattering or assign a measured mass.
