# Stationary low-energy states of the spatial gauge-coordinate source

## 1. The source, cutoff and common clock

Use the native fifteen-dimensional representation of
`G=S(U(3)×U(2))`, four Dirac components and all sixty canonical
anticommutation-relation modes per spatial cell. A directed link has
coordinate Hilbert space `L²(G,Haar)`. Its finite subspace `V_ell` is the
span of products of at most `ell` native matrix coefficients or their
complex conjugates, including the constant function. Endpoint gauge
transformations preserve this subspace. Electric Casimirs preserve it.

The regulator here is the **product of the per-link word spaces**,

```
P_ell = tensor_links P_(V_ell),
H_(Lambda,ell) = P_ell [H_E + dGamma(h) + H_B] P_ell,
H_E = e0/(2a) sum_links Casimir_trace,
H_B = -kappa0/a sum_plaquettes Re Tr_M15(U_p)/15.        (1)
```

The full fermionic Fock factor is unchanged. Products in (1) are formed
in the parent coordinate algebra before compression. Substituting products
of compressed link operators defines a different expression in general.
The product cutoff differs from the total-word cutoff in
[the field route](../particle_observer_outputs/FIELD_ROUTE.md).
It is an explicitly chosen regulator. It provides compatible finite
local factors when spatial volume grows at fixed `ell`.

For the finite-range Wilson member,

```
h = (beta_D + H_W)/a = beta_D D_W,
H_W = beta_D(aD_W-I),    aD_W=3I-sum_(j=1)^3 W_j.      (2)
```

Here `beta_D` is the Dirac matrix, not inverse temperature. Each `W_j`
is a unitary covariant shift resolved by orthogonal spin projections,
so `||h||<=6/a`. The same bound holds with open boundaries by compression
of the one-body shifts. Every cell owns its three outgoing link factors.
Open-box Hamiltonians retain precisely those interaction terms whose
support lies inside the box; dangling outgoing links retain their
electric factor. This makes the sequence restrictions of one fixed local
interaction. Periodic boxes have the same conservative finite bounds.

The smooth alternative is the declared action
`h_delta=(beta_D+sqrt(1+delta²)H_W(H_W²+delta²I)^(-1/2))/a`
of [the smooth-field proof](../particle_smooth_sign_field/PROOF.md).
Its bound is `(2+delta)/a`. Every finite-path polynomial truncation in
that proof obeys this bound as well. The finite-volume cooling estimates
below apply to each of these specified actions. The spatial ground-state
subsequence established in section6 uses the finite-range Wilson action.

Time throughout is `tau_t(A)=exp(itH) A exp(-itH)` for the selected
finite-cutoff Hamiltonian in (1). Subtracting a scalar inside a preparation
polynomial changes neither this clock nor any transition energy.

## 2. Source dimension and spectral bounds

The actual native trace metric gives conventional SU(3) and SU(2) indices
2 and integer hypercharge `q=6Y` index120. The source packet and
[`native_casimirs`](../particle_field_preparation/preparation.py) verify
these values. The generator-norm proof in
[the preparation proof, section5](../particle_field_preparation/PROOF.md)
gives `sum_A ||t_A||²<=29/5`. On words of length at most `ell`, the
Casimir is bounded above by `(29/5)ell²`.

For `V` cells, `E<=3V` link factors and `P<=3V` plaquettes,

```
D <= 2^(60V) 451^(ell E) <= 2^b,    b=60V+9ell E,
log D <= b,
C_E = (29/10)(e0/a) E ell²,
B = 360V/a + (kappa0/a) P,
-B I <= H <= (C_E+B) I,    W=C_E+2B.                  (3)
```

Each of `ell` word slots has451 choices: identity,225 matrix entries,
or225 conjugates. Linear dependence only decreases the dimension.
For the smooth action replace `360V/a` by `60V(2+delta)/a`.
These are bounds on the complete carrier, including every occupation
and every admitted link word. There is no projection onto a chosen
particle number or a Gauss-singlet source sector.

The source metric is retained when selecting coefficients. For example,
`e0=105, kappa0=1/7` is one declared member. The cooling proof does not
infer any continuum speed or measured coupling from these choices.

## 3. Gap-free cooling and deterministic canonical preparation

For any nonempty finite-dimensional Hermitian `H`, write `E0=min spec H`
and `rho_beta=exp(-beta H)/Tr exp(-beta H)`, where `beta>0`. Spectral
diagonalization reduces its entropy to the finite Gibbs entropy. The
identities `gibbsEntropy_normalization` and
`gibbsEntropy_le_log_card` in
[`GibbsEntropyLimit.lean`](../../Lean/Thermodynamics/GibbsEntropyLimit.lean)
give

```
S(rho_beta)=beta Tr(rho_beta H)+log Z <= log D.
```

Since the ground term belongs to `Z`, `log Z>=-beta E0`. Consequently

```
0 <= Tr(rho_beta H)-E0 <= log D/beta <= b/beta.        (4)
```

No spectral gap, unique ground, fitted mass or oracle for `E0` enters
(4). A family whose excitation gaps tend to zero obeys the same bound.
[`GapFreeCooling.lean`](../../Lean/Thermodynamics/GapFreeCooling.lean)
formalizes the attained-minimum energy bound directly from the cited
finite entropy identities.

Take a certified interval `[L,U]`, `W=U-L`, and a dyadic integer
`s>=max(1,beta W)`. Define the finite polynomial

```
X=beta(U I-H)/(2s),    0<=X<=yI,    y=beta W/(2s)<=1/2,
T_m(X)=sum_(j=0)^m X^j/j!,    Q=T_m(X)^s,
rho_Q=Q²/Tr Q²,    Psi_Q=vec(Q)/||Q||_F.              (5)
```

The vectorization uses the full source carrier and its conjugate reference.
`Q` is positive definite and commutes exactly with `H`, all exact source
symmetries and parity. Thus `rho_Q` is exactly stationary. The reference
has zero Hamiltonian. The joint purification need not be stationary under
`H tensor I`; all stationarity claims concern the source marginal.

There is a deterministic finite unitary preparation. Let `b0` be the
empty-Fock/constant-link source vector and its conjugate reference. Both
`b0` and `Psi_Q` are invariant under joint gauge transformations and even
parity, and their overlap is real positive. If they differ, put
`v=(b0-Psi_Q)/||b0-Psi_Q||` and `R=I-2|v><v|`; then `R b0=Psi_Q`.
If they coincide, use `R=I`. Extend this operator by the identity on the
orthogonal symmetry sectors. This retains the full reference and every
source state, with no successful-branch conditioning. It establishes
joint gauge invariance, not a claim that the source density is supported
in a no-external-charge Gauss sector. A singlet-only source ground would
be a different minimum and preparation problem.

For the Wilson source, finite Haar moments, representation matrices and
orthogonalization give algebraic finite matrix entries. Equation(5) and
the Householder preparation are finite algebraic constructions. The
finite instrument compiler consumes their matrix description to any
specified output accuracy. Exact smooth-coordinate compression can have
computable-real entries; its finite-path polynomial action supplies an
algebraic alternative with its own stationary state. A finite-precision
approximation of a target state is not called exactly stationary here.
Its additional output error must be included separately.

## 4. Certified error and finite resource plan

For `0<=z<=y`, the positive Taylor tail is bounded by

```
r_m = y^(m+1)/(m+1)! / (1-y/(m+2)).
```

Since `exp(z)>=1`, `0<=1-T_m(z)/exp(z)<=r_m`. For `s` powers,
Bernoulli's inequality gives relative error at most `s r_m`.
All factors are functions of the same `H`, so in the eigenbasis

```
||Q-exp(beta(U I-H)/2)||_F
 <= s r_m ||exp(beta(U I-H)/2)||_F.
```

Normalizing two nonzero vectors changes this bound by at most a factor2.
The half trace distance of pure canonical states is at most their vector
norm distance; partial trace contracts it. Thus

```
D_half(Psi_Q,Psi_beta) <= eta=min(1,2s r_m),
D_half(rho_Q,rho_beta) <= eta,
Tr(rho_Q H)-E0 <= b/beta + W eta.                     (6)
```

The last term follows by centering `H` at `(L+U)/2`; the difference of
normalized densities has trace zero. It is `W eta`, with no dependence
on an arbitrary energy zero.

For any rational requested `epsilon>0`, choose
`beta=2b/epsilon`, the least dyadic `s>=max(1,beta W)`, and increase `m`
until `eta<=min(1/4,epsilon/(2W))` (the scalar case `W=0` is immediate).
Factorial decay makes this loop terminate. Equation(6) then gives global
mean ground-energy excess at most `epsilon` on the actual complete source.
The planner requires no diagonalization or ground gap.

The polynomial degree is `ms`. Its factored evaluation uses at most
`m+log2(s)` matrix multiplications before forming the density. This is
not an algorithm efficient in spatial volume: storing a generic `Q`
requires `D²` entries, dense multiplication costs order `D³` arithmetic
operations, and rational coefficient bit lengths grow with the powers.
The canonical vector has `D²` amplitudes; its generic preparation has
exponential resources. `spatial_plan` reports bounds and a finite
preparation specification without allocating these matrices. This
existence-and-budget construction supplies no practical large-volume
vacuum simulation or state convergence rate.

## 5. A local ground inequality for the same stationary state

Let `[rho,H]=0`, `Tr rho=1`, and `Tr rho(H-E0)<=epsilon`. For arbitrary
`A`, including charged odd probes with compensating references,

```
omega(A* [H,A])
 = Tr rho A* (H-E0) A - Tr rho A*A (H-E0)
 >= -||A||² epsilon.                                 (7)
```

The first term is nonnegative. Commutation makes `rho(H-E0)` positive,
so cyclicity bounds the second term by `||A||² Tr rho(H-E0)`.
Both terms are real. Dropping stationarity would invalidate this step.
The polynomial state (5) supplies it exactly. The Gibbs approximation
error is used to bound energy excess, not to substitute approximate
commutation into (7).

## 6. A fixed-cutoff spatial algebraic ground subsequence

Fix `a>0`, finite `ell` and finite coefficients `e0,kappa0`. Give each
cell its full sixty-mode CAR factor and three outgoing link word factors.
Keep odd CAR observables and even Hamiltonian interactions. For the Wilson
source, every interaction has at most four cell factors and diameter at
most two lattice edges. A conservative incident interaction bound is

```
J = 900/a + (87/10)e0 ell²/a + 12kappa0/a,
sum_(Z contains x) ||Phi_Z|| <= J.                    (8)
```

The onsite term has norm at most `180/a`. Each of at most six incident
fermion hopping bonds has norm at most `120/a`, from its full120-mode
one-body support. The three electric factors give the second term.
Twelve plaquette terms per cell is a conservative upper count.
Compression does not increase these norms. All terms are even, so
disjoint graded observables commute with them.

This supplies the quasi-local dynamics without transferring numerical
bounds from a different contact model. For `A` on `s0` cells, iterated
commutators satisfy

```
||ad_H^n(A)|| <= ||A|| (2J)^n product_(j=0)^(n-1)(s0+3j). (9)
```

At each commutator, the interaction must meet the accumulated support,
its incident norm sum is bounded by `J` per cell, and it adds at most
three cells. The exponential series therefore converges uniformly in
volume for `|t|<1/(6J)` when `J>0`; `J=0` is trivial. Each fixed-order
coefficient stabilizes once its finite support lies inside the box.
Dominated convergence gives the local infinite-volume evolution on this
time interval. Isometry, density and composition extend it to a strongly
continuous automorphism group for all finite times. This is the same
local-algebra method used for the distinct two-mode contact action in
[thermodynamic_ground](../thermodynamic_ground/README.md), with the actual
quantum-link interaction bounds (8).

Take growing open cubes and the explicit cooling plans with
`epsilon_R=1/R`. Extend each exactly stationary finite density outside
its cube by the fixed even, gauge-invariant tracial product state. Local
matrix algebras have a countable dense collection of entries. Bounded
state values admit a diagonal subsequence converging on those entries;
positivity, normalization and the norm bound extend its limit to a
state on the completed quasi-local algebra.

For any fixed local `A`, all terms meeting its support lie inside every
sufficiently large cube. Equation(7) therefore reads
`omega_R(A*[H,A])>=-||A||²/R` for a fixed local expression. Passing to the
subsequence gives `omega(A*[H,A])>=0`. Uniform local dynamics convergence
and exact finite-volume stationarity pass stationarity to the limit.
The state is a fixed-spacing, fixed-link-cutoff algebraic ground state
of this specified Wilson quantum-link action.

Spatial averaging can additionally give a translation-invariant ground
state: average the extended states over the cube before extracting the
subsequence. Interior translates obey (7); the fraction of boundary
translates tends to zero, and their contribution has norm at most
`2s0 J ||A||²`. The same interior/boundary argument applies to local
time evolution. This is the averaging construction in the cited ground
proof with bound (8), including the quantum links and all CAR sectors.

Compactness selects no computable subsequence or unique ground. This
result has no `a->0` or `ell->infinity` convergence assertion, Lorentz
covariance, chiral mirror removal, measured mass identification or
Einstein dynamics. Those statements cannot be inferred from small mean
energy excess or from a three-dimensional lattice alone.

## 7. An exact cooled native control

`prepared_control(beta=4,order=8)` uses the unchanged208-state native
singlet/quantum-flux Hamiltonian: four Dirac modes, all sixteen occupation
states and all thirteen flux labels. Its complete 82 reducing blocks have
maximum dimension 12. The source interval is `[-4,22]`, so the construction
uses `s=128`, `y=13/32`, polynomial degree1024 and

```
eta <= 10604499373/47847722383835136 < 2.217e-7.
```

Every amplitude block is evaluated exactly over rational matrices in the
retained real CAR phase chart; the source phases are returned with it.
`native_control_summary` combines exact mean-energy arithmetic with the
independent rational positive-definiteness ground brackets in
[the spectral record package](../particle_ground_spectral_records/PROOF.md).
It certifies mean ground-energy excess below `0.157` source units;
the mean energy is approximately `-2.00394106127` and the full ground is
approximately `-2.16034436`. Decimal values are diagnostics. The receipt
contains outward rational intervals and the rational upper excess.
This is a complete one-cell control for cooling and detector records,
not a three-dimensional continuum or a rest-mass determination.
