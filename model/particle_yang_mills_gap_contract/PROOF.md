# Euclidean decay, the physical clock and a continuum spectral gap

## 1. What this criterion supplies

An exponential bound on a dense family of centered Euclidean correlations
excludes low-energy spectrum of their reconstructed Hamiltonian. The
argument needs neither equality of prelimit and limit spectra nor convergence
of every finite eigenvalue. It uses the same limiting correlation functions
and time translations that define the quantum theory.

This strengthens the usable interface of
[`RepairGapChain.mass_gap`](../../Lean/ObserverPatchHolography/RepairGapChain.lean).
That logical theorem consumes an assumed repair gap and an assumed equality
of nonzero spectra. The criterion below gives a direct route to spectral
exclusion from local correlation estimates. It does not establish those
estimates or reconstruct a four-dimensional Yang--Mills theory.

The [Jaffe--Witten problem](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf),
sections 3--5, requires a nontrivial theory for each compact simple group on
four-dimensional flat spacetime, axiomatic quantum-field properties, the
Yang--Mills short-distance behavior, and a positive finite vacuum threshold.
A stable isolated one-particle state and confinement are further questions.
The theorem below addresses the spectral threshold after the required
quantum theory and source identification have been established.

## 2. Dense local correlations determine spectral exclusion

Let a Hilbert space carry a strongly continuous unitary time group with
self-adjoint generator K>=0 and normalized vector Omega satisfying
K Omega=0. Let A be a unital bounded local star algebra whose vectors
A Omega are dense. No Lorentz covariance or Yang--Mills interpretation is
inferred from these assumptions. They name one reconstructed source.

For a in A define

```
v_a=(a-<Omega,a Omega>I)Omega,
C_a(t)=<v_a,exp(-tK)v_a>.
```

It suffices that a collection of such vectors has dense linear span in
Omega-perp. Fix h>0 and 0<q<1. Suppose each vector in this collection has
a finite constant B_a such that

```
0 <= C_a(kh) <= B_a q^k          for every integer k>=0.       (1)
```

Then, with delta=-log(q)/h,

```
K >= delta (I-|Omega><Omega|).                               (2)
```

In particular the vacuum is unique in this cyclic Hilbert space and
Spec(K) has no points in (0,delta). A shared q is required; B_a may depend
on the observable.

Proof. Let mu_a be the positive spectral measure of K in v_a. For
0<=E<delta, spectral calculus gives

```
exp(-khE) mu_a([0,E]) <= C_a(kh) <= B_a q^k,
mu_a([0,E]) <= B_a exp[-kh(delta-E)] -> 0.                    (3)
```

Thus the spectral projection 1_[0,E](K) annihilates every v_a. Its
boundedness makes it zero on their closed linear span Omega-perp.
It acts as the identity on Omega. Increasing E to delta proves (2),
including exclusion of an additional zero-energy vector. No assumption
about isolated eigenvalues enters this proof.

Conversely, (2) gives (1) with B_a=||v_a||² by spectral calculus. Thus
the dense decay criterion is equivalent to the gap inequality in this
positive-energy cyclic representation. A bound for one selected probe
certifies only its spectral measure, unless its generated vectors are
proved dense.

## 3. Fixed-time limits and complete errors

At regulator n, declare one Hamiltonian H_n, exact ground Omega_n and
ground energy E_n. Declare s_n>0 as the number of source time units per
one physical time unit, before comparing spectra. The Euclidean generator
and transfer are

```
K_n=s_n(H_n-E_n I),       T_n=exp(-h K_n).                    (4)
```

The ground-energy subtraction is explicit. It leaves source Heisenberg
observations unchanged; reference gates, absolute phase conventions and
the retained source energy zero are separate pieces of the experiment.
The source preparation in (4) cannot be replaced by a different warm
state, purification or reference clock without a proved comparison.

Suppose local representatives a_n of each fixed continuum test satisfy

```
C_{a,n}(kh) -> C_a(kh)                       for each fixed k,
0 <= C_{a,n}(kh) <= B_a q^k + epsilon_{a,n}(k),
epsilon_{a,n}(k) -> 0                        for each fixed k. (5)
```

The limiting C_a must be the diagonal correlations of the very same
positive-energy source and time translations in section 2. This is the
attachment premise. If computed or apparatus estimates replace C_{a,n},
their complete absolute error must be included in epsilon. Failed
outcomes and outside-code weights cannot be discarded. Exponential
imaginary-time operators are not automatically admitted real-time
detectors; a measurement or analytic reconstruction needs its own error
bound.

First take n to infinity at a fixed k in (5). This proves (1). Then
take k to infinity in (3). Neither uniform convergence for unbounded
time nor equality of finite and infinite spectra is required. This
order of limits is essential: epsilon may be amplified by exp(khE).

A sufficient special case is the uniform finite-source inequality

```
K_n >= delta(I-|Omega_n><Omega_n|),
sup_n ||v_{a,n}||² <= B_a.                                  (6)
```

It gives (5) with no error for exact matrices. The more general (5) can
hold while global finite gaps close: a mode of energy tending to zero
can have original weight tending to zero in every limiting local test.
For example the exact three-level controls use transfer eigenvalues
1, 1-1/(n+1), 1/2 and centered weights 0, 1/n, 1-1/n. Their correlation
is at most 2^(-k)+1/n and tends to 2^(-k) at each fixed k. The limiting
cyclic sector is gapped, although the full prelimit gap tends to zero.
This example prevents replacing (5) by an asserted necessary uniform
gap of every full prelimit carrier.

## 4. Noncollapse gives a finite upper threshold

Suppose one nonzero centered limiting vector v lies in the form domain
of K and satisfies

```
||v||² >= v0 > 0,       ||K^(1/2)v||² <= M < infinity.       (7)
```

Then the non-vacuum spectrum is nonempty and its infimum m obeys

```
delta <= m <= M/v0 < infinity.                             (8)
```

Indeed, the probability measure mu_v/||v||² is supported on the
non-vacuum spectrum, so its mean energy is at least m. Equation(7)
then bounds that mean by M/v0. This is the finite-energy noncollapse
condition in item(Y6) of the
[Yang--Mills continuum certificate](../../extra/yang_mills_gap_clay_problem.tex).
There are two sufficient ways to derive (7) from finite stages. Norm
convergence of centered vectors preserves their variance floor, while
an appropriate form lower semicontinuity transfers the uniform energy
upper bound. Alternatively, no Hilbert-space identifications are needed
if the same centered correlations also converge at time zero and at
a fixed sequence t_j decreasing to zero. Suppose

```
C_{a,n}(0)>=v0,       <v_{a,n},K_n v_{a,n}><=M.
```

For every j, the scalar inequality 1-exp(-t_j E)<=t_j E gives

```
0 <= [C_{a,n}(0)-C_{a,n}(t_j)]/t_j <= M.
```

Pass n to infinity for this fixed j. Then let j tend to infinity.
The spectral functions (1-exp(-t_j E))/t_j increase to E for E>=0.
Monotone convergence yields ||K^(1/2)v_a||²<=M, including membership
in the form domain. Time-zero convergence gives ||v_a||²>=v0.
This proves (7) from correlation and energy controls of the same
source, without Mosco convergence or a preassigned limiting vector
identification. Weak state compactness alone supplies neither route.

A sequence K_n=n diag(0,1) has gaps diverging to infinity and fixed
nonzero centered variance. Its positive-time correlations tend to zero,
while their time-zero value stays positive. This limiting function is
discontinuous at zero, so it cannot be a correlation of a strongly
continuous nontrivial quantum time group. The finite vector energies
also diverge. A lower gap bound alone does not prevent this collapse.

## 5. Exact clock accounting for the selected Wilson source

The [primary contract](../particle_common_contract/PROOF.md) fixes
a_n=2^(-n), L_n=2^(2n+1), ell_n=n+1, e0=105 and kappa0=1/7.
Every term of its compressed Hamiltonian has one factor 1/a_n:

```
H_n = a_n^(-1) Hhat(L_n,ell_n;105,1/7).                      (9)
```

The word projections and Gauss projector are independent of a at fixed
L and ell, so restriction preserves this exact identity. Hhat changes
with volume and cutoff. For its dimensionless physical-sector gap
g_n and the declared clock exp(-it H_n), the physical gap is g_n/a_n.
There is no additional fitted time conversion in (9).

If the uniform finite-gap route (6) and a uniform finite-stage witness
with variance at least v0 and physical energy at most M both apply,
then necessarily

```
a_n delta <= g_n <= a_n M/v0.                              (10)
```

In this route g_n has the scale of a_n. A constant positive lower
bound on g_n instead drives g_n/a_n to infinity and is incompatible
with that noncollapse witness. Changing source time after seeing this
behavior defines a different clock and requires a fresh spacetime
identification. Equation(10) is not asserted as necessary for the
vanishing-residue alternative of section 3.

The selected coefficients are fixed along n; a trajectory of running
bare couplings and asymptotically free short-distance correlations is
not supplied. The source also retains sixty matter modes per cell and
a gauge representation of a product group, rather than a declared pure
compact-simple Yang--Mills theory. A pure gauge comparison must name
its own action, carrier and vacuum, or prove a reducing sector of this
one. Choosing one factor of the Lie algebra does not prove that its
Hilbert space reduces the coupled M15 Wilson action. Fixing other link
variables to the identity is not normalizable in Haar L². Their
constant-word subspaces need not be invariant under plaquette
multiplication either. A separate pure-gauge source or a proved
decoupling map is required. No positive global gap of a matter-containing
Standard Model source follows from a pure Yang--Mills gap.

The selected polynomial preparation has a ground-energy excess bound.
An excess upper bound is neither a spectral-gap lower bound nor a
bound on distance to a unique ground without further hypotheses.
Its exact stationarity does not replace (4) or establish (5). For
example H=diag(0,1), rho=diag(1-epsilon,epsilon) and A=|0><1| give
the stationary ordered correlation

```
Tr[rho A* exp(-tH) A exp(tH)] = epsilon exp(t).              (11)
```

Its ground-energy excess is epsilon, but it contains a negative energy
transfer. No uniform decaying vacuum bound follows. If epsilon tends
to zero it gives a vanishing error for each fixed t; that limit order
is permitted in (5). An energy-excess upper bound alone does not
control exponential moments of arbitrary negative transfers when
source spectral widths grow. The required time-dependent error must
therefore be proved for the actual prepared correlation, with its
retained reference, rather than inferred from stationarity or cooling.

For a distinct native screen source, the retained diagnostic eigenvalues
lambda_l of its port Laplacian give the per-attempt relaxation scale
lambda_l/(2M_l), where M_l=630*4^l. The declared diffusive clock has
1260*16^l attempts per comparison unit. Hence its scale is exactly
4^l lambda_l in that clock. The three retained decimal diagnostics
become approximately 0.04414, 0.04126, 0.04087, while their per-attempt
values fall by factors near 16. These are floating eigenvalue diagnostics
from the two-dimensional comparison source, not interval certificates
or mass gaps of the selected Wilson source.

## 6. The remaining Yang--Mills premises and the rate comparison

The [existing Yang--Mills paper](../../extra/yang_mills_gap_clay_problem.tex)
proves a compact cylinder/GNS subsequence extraction and a conditional
spectral-gap transfer. Its continuum certificate separately requires
renormalized local convergence, reflection positivity, four-dimensional
Euclidean covariance and regularity, the common transfer clock and
vacuum identification, and finite-energy noncollapse. The criterion
here can replace its spectral-gap transfer step if (5) is proved for
the reconstructed local algebra. It supplies none of the preceding
quantum-field existence or Yang--Mills identification premises.

One useful source comparison is established in
[`dirichlet_domination_family`](../../Lean/Screen/KogutSusskindFiberRateComparison.lean):
for pi-symmetric idempotent conditional expectations E_l and positive
fiber-constant rates c_l>=c*, the actual variable-rate form dominates
c* sum_l ||(I-E_l)f||². A variance inequality with one constant A on
the same non-product gauge quotient would therefore give gap>=c*/A.
In the finite Z2 Kogut--Susskind control the exact rate is
lambda(r+1/r)>=2lambda. Scalar cross-fiber equality is unnecessary.
The independent-product Walsh theorem supplies no such variance
inequality for the interacting quotient law. Neither that finite
discrete-group control nor a local electric Casimir gap proves (5)
uniformly along a four-dimensional compact-simple continuum family.

## 7. Executed checks

`finite_transfer` uses exact algebraic spectral projectors of actual
positive finite matrices. It retains all original centered spectral
weights, including zero weights, and checks whether the declared probe
vectors span the full vacuum complement. The tests calculate powers
of non-diagonal transfer matrices, distinguish fast probes from slower
unobserved modes, reject a changed clock and a collapsing gap floor,
and retain the vanishing-residue and divergent-energy controls.

The receipt authenticates two selected primary members and their
unchanged source clocks. It computes no primary-field spectral gap,
continuum Euclidean decay rate, or particle mass.
