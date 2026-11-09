# An interacting field state with complete thermal records

The unchanged native charge-six quantum-flux control has a complete
finite Gibbs target, a deterministic reference-preserving preparation
with an exact error bound, and complete two-time flux records. Its
thermodynamic identities concern the same Hamiltonian that drives those
records. A separate exact variational comparison shows that the spatial
degree-zero mixed preparation is not a finite ground state. Neither
result identifies a physical vacuum or a relativistic continuum.

## 1. Same finite source, occupation space and clock

Use `field_route.exact_rotor(6)` from the
[field source](../particle_observer_outputs/FIELD_ROUTE.md). The native
singlet has charge `q=6`. There are four source Dirac CAR modes, all
sixteen Fock states, and every electric flux `e=-6,...,6`. Thus the
carrier has dimension208. Let S be the compression of multiplication
by the charge-six link character; it shifts flux by six, with zero
outside this declared finite carrier. Write B for the second
quantization of the source beta matrix and K for its source hopping
matrix. The actual action is

```
H=I tensor E²/2+B tensor I+K tensor S+K† tensor S†,
G=6(N-2I) tensor I,       [H,G]=0,
J=i[H,E]=6i(K† tensor S†-K tensor S).
```

This is exactly the existing source, including its cutoff and all
occupation sectors. No state is admitted or removed by an eigenvalue
test. Its clock is `exp(-itH)` in the declared source units. The
reference Hamiltonian is zero. The selected cutoff-six control is not
called an approximation to an unbounded-flux thermal state.

For comparison with coefficient response define the exact Hermitian
components

```
H_E=I tensor E²/2,
H_W=B tensor I+(K+K†)/2 tensor (S+S†),
H_D=(K-K†)/2 tensor (S-S†),       H=H_E+H_W+H_D.
```

The suffixes identify electric, Wilson and Dirac terms. Geometric
interpretation of a coefficient variation requires its declared metric
dictionary; the decomposition alone is an operator identity.

The uncompressed one-body coordinate symbol obeys
`h(theta)²=4 sin²(theta/2)I`, with positive and negative multiplicity
two. Its full second quantization has norm at most4. Compression does
not increase that norm, and `0<=H_E<=18I`. Hence
`-4I<=H<=22I`. The code uses exact connected components only to evaluate
the full matrix economically: dimensions `1:52,4:20,6:4,8:5,12:1`
sum to208. They are reducing blocks, without a selected-sector state.

## 2. A declared Gibbs target and a complete preparation

Choose the dimensionless inverse temperature `b=1/64` and target

```
rho_b=exp(-bH)/Tr exp(-bH),
psi_b=vec(exp(-bH/2))/sqrt(Tr exp(-bH)).
```

Temperature is a declared preparation parameter, without a laboratory
calibration or autonomous cooling claim. This target differs from the
spatial gauge-twirled product with a tracial Fock factor. A declaration
of a preparation primitive is allowed by the construction contract;
the explicit approximation below gives finite rational instructions
and an error budget for that primitive.

Put `A=b(22I-H)/2`; then `0<=A<=xI` with `x=13b`. For integer M with
`x<M+2`, define

```
Q_M=sum_(j=0)^M A^j/j!,
epsilon_M=x^(M+1)/(M+1)! / (1-x/(M+2)),
rho_M=Q_M²/Tr Q_M²,       psi_M=vec(Q_M)/sqrt(Tr Q_M²).
```

All matrices Q_M and rho_M have exact Gaussian rational entries on this source.
Every term is positive, so `Q_M>=I`; `rho_M` is faithful and normalized.
The exponential tail is bounded by the first omitted term times the
displayed geometric majorant. Thus `||Q_M-exp(A)||<=epsilon_M`.
The harmless factor `exp(11b)` cancels in the normalized target.

For nonzero vectors u,v,
`||u/||u||-v/||v||||<=2||u-v||/||v||`.
Here the Frobenius numerator is at most `sqrt(208)epsilon_M` and
`||exp(A)||_F>=sqrt(208)`. Consequently the half-trace distance of
the canonical pure states, and hence their source marginals, is at
most `min(1,2epsilon_M)`. This is independent of a reference dimension
and is an analytic bound, not a floating diagonalization comparison.
The default M=12 makes this error much smaller than the subsequent
time-evolution allowance.

For efficient exact arithmetic the code uses the fixed diagonal CAR
basis change `P|mask,e>=i^(number of occupied negative-beta modes)|mask,e>`.
The matrix `P^dagger H P` is real rational. Every connected block is
computed in that chart and the returned density is restored as
`P rho_real P^dagger`. This is a basis change, not a changed action.
Flux projections commute with P, so complete flux probabilities are
unchanged. The canonical reference carries the conjugate basis change.

The preparation has a complete deterministic unitary implementation.
Take a blank matched source/reference basis vector e0 in the joint
Gauss-zero and even-parity subspace. If `v=e0-psi_M` is nonzero, the
Householder operator `W=I-2vv†/(v†v)` maps e0 to psi_M; otherwise use I.
The chosen vector component at e0 is real, so this formula applies
without a phase adjustment. W is defined on every input, with no
success flag or accepted outcome. Since Q_M commutes with G and parity,
psi_M is invariant under `U_G tensor conjugate(U_G)` and has even total
parity. Both e0 and v lie in that invariant subspace, so W preserves the
joint gauge constraint and parity on its full carrier. Dense matrices
of dimension208² are specified by this formula, not materialized.

The same statements hold for psi_b. The source densities are stationary
under H because they are functions of H. Their purifications generally
are not stationary under `H tensor I`. The reference is retained as a
preparation resource, not silently made into a thermal bath.

## 3. Complete rational flux histories

Let P_f be all thirteen original flux projections, including every
occupation. The two-time instrument is `P_g U P_f`, with reference
identity, where U is the same source-time evolution. It is complete:

```
sum_(f,g) (P_g U P_f)†(P_g U P_f)=I.
```

For an exactly unitary rational evolution define, at J positive integer
substeps and total duration t,

```
U_(t,J)=[(I-itH/(2J))(I+itH/(2J))^(-1)]^J.
```

Hermiticity makes the two factors adjoints and commuting. The inverse
exists for every real t, so U_(t,J) is unitary on all208 states and
commutes with the same Gauss operator. It is a compiled approximation
to the original clock, not a redefinition of the exact Hamiltonian.
The scalar inequality `|atan z-z|<=|z|³/3` gives

```
||U_(t,J)-exp(-itH)|| <=min(2, |t|³22³/(12J²)).
```

For the actual prepared source/reference state the original joint
weights are computed exactly as

```
p_M(f,g)=||P_g U_(t,J) P_f Q_M||_F²/Tr Q_M².
```

The certificate retains all169 `(f,g)` labels, including zeros; their
exact rational weights sum to one. It uses every reducing block and
every occupation, with no renormalization after discarding a record.
Contraction under the full instrument and the unitary comparison bound
give half-trace error, including final conditioned states and the
untouched reference, at most

```
min(1, 2epsilon_M+|t|³22³/(12J²))
```

relative to the Gibbs target under the exact clock. The default
`t=1/64,J=4` has evolution allowance `1331/6291456`. More substeps
reduce this explicit error. Measuring the first flux disturbs the
stationary source Gibbs state; its subsequent history is calculated
with that backaction retained.

## 4. Thermodynamics on this source algebra

The target is Gibbs for the actual H, rather than Gibbs for an
unrelated modular operator chosen after the fact. In its spectral
basis the existing
[Gibbs energy identification theorem](../../Lean/Thermodynamics/GibbsReferenceEnergyIdentification.lean)
applies to all eigenvalues with multiplicity. Functional calculus gives
its matrix form

```
-log rho_b=bH+log Z I,      S(rho_b)=b<E>+log Z.
```

The finite KMS identity follows by cyclicity of trace, for arbitrary
source A,B and `alpha_z(B)=exp(izH)B exp(-izH)`:
`Tr rho_b A alpha_(t+ib)(B)=Tr rho_b alpha_t(B) A`.
At fixed H, `d<E>/db=-Var(H)` and heat capacity in the declared units
is `b² Var(H)>=0`. The stationary current mean vanishes, whereas its
second moment can be nonzero. That does not suppress the two-time
record dynamics.

For a cyclic source unitary U the Gibbs state is passive. In the energy
basis its populations decrease with energy; the doubly stochastic
matrix `|U_ij|²` cannot lower their energy-weighted sum. The same
rearrangement proves source passivity for rho_M, because the positive
polynomial Q_M is monotonically increasing in `22-H`. Both states
are stationary. Rho_M is generally not exact Gibbs or KMS at b, and
no complete-passivity claim is made for rho_M. Source passivity does
not apply to arbitrary joint unitaries using the purifying reference;
that reference contains accessible correlations and purity.

For a differentiable same-source coefficient family H(s), put
`K=H'(s)` and `rho_b(s)=exp(-bH(s))/Z(s)`. With
`Cov_KM(A,B)=integral_0^1 Tr[rho^u(A-<A>)rho^(1-u)(B-<B>)]du`,
the finite Duhamel derivative of the exponential gives

```
F'=<K>,       F=-(log Z)/b,
E'=<K>-b Cov_KM(H,K),
S'=-b² Cov_KM(H,K),
d<A>/ds=<A'>-b Cov_KM(A,K).
```

These are target Gibbs identities. In particular, `S'=b<E'>` is false
for varying H unless its work term `<K>` vanishes. At a fixed target
state, any trace-zero state tangent obeys `dS=b Tr(H d rho)`. This
supplies the finite entropy first law with the actual source energy.
The [existing finite work identity](../../Lean/Thermodynamics/FirstLawIdentity.lean)
keeps its bilinear cross term for finite strokes; at first order it is
`dE=Tr(H d rho)+Tr(rho dH)`. Bounded record expectations transfer to
rho_M with error at most twice their operator norm times the certified
half-trace distance. Derivatives and logarithms do not transfer through
that bound without their own estimates.

The [finite low-temperature theorem](../../Lean/Thermodynamics/LowTemperatureControl.lean)
also applies in this exact source energy basis. A ground gap and
dimension bound give excited weight at most `208 exp(-b gap)`, and
the target concentrates on its full minimum space as b increases.
The present default b is not a ground-state preparation. No gap
uniform in volume or spacing, autonomous cooling dynamics, physical
temperature scale, or interacting spatial vacuum follows here.

## 5. A concrete boundary of the spatial mixed preparation

Consider the separate degree-zero Legendre field on a periodic cubic
carrier of side at least two. Its full gauge-state factor is the one
specified in the [local state construction](../particle_local_gauge_state/PROOF.md),
and its fermion factor is tracial on all60V modes. The actual h0 has
onsite coefficient `(25/(17a)) beta tensor I_15`; other terms connect
distinct cells. The mixed fermion covariance is I/2, so the full kinetic
mean is zero.

Replace only that fermion factor by the product Slater state filling
the two negative-beta spin rows for every native internal row at every
site. This occupies thirty modes per site and is a local gauge singlet.
The semisimple representation determinants are one and the sum of all
fifteen native integer hypercharges is zero. Its occupied determinant
character is therefore `(det R)^2=1`. The two states have the same
gauge density, electric mean and plaquette mean. Gauge twirling leaves
the alternative fermion factor invariant. The coordinate-word
compression acts only on gauge registers and changes none of these
comparisons.
The gauge factor can itself carry charge. Both complete densities have
canonical purifications on the same conjugate reference carrier, with
zero reference Hamiltonian and the same joint Gauss constraint. The
energy comparison is therefore admissible on that common referenced
domain; it does not assume that the gauge factor is a source-only singlet.

Both states have zero intersite one-body coherence, so every hopping
mean vanishes. In the filled state each cell instead has onsite kinetic
energy `-30*25/(17a)`. Thus its total source mean energy is lower by

```
750 V/(17a)>0.
```

This exact same-H comparison proves that the mixed spatial preparation
is not a finite ground state. It uses side at least two to exclude
self-link terms. Its coefficient belongs to h0 and is not assigned to
the positive-smoothing action. It is not a theorem excluding every
physical continuum interpretation of the mixed family, and it does
not prove that the Slater alternative is the interacting ground.

## 6. Executable checks and scientific boundary

The code reconstructs the actual source operators, rational positive
preparation and complete Cayley records. Separate floating diagnostics
check target Gibbs energy, variance and stationary current on all208
spectral rows; they are not interval certificates. Tests include
exact normalization, operator identities, reference/pointer scope,
complete labels, and independent matrix comparisons.

This result supplies finite same-source thermodynamic premises and
an explicit state-preparation path. It supplies no thermodynamic
volume limit of these Gibbs states, geometric diamond modular kernel,
Lorentz vacuum, particle species, GeV mass or Einstein dynamics.
