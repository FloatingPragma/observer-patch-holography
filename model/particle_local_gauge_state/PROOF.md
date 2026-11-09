# A referenced gauge state with uniform local energy densities

The construction gives a family of finite gauge-coordinate experiments
with one local word degree and one plaquette coefficient for every volume.
A mixed state obtained by gauge twirling has explicit electric and magnetic
energy-density bounds. Its retained conjugate reference makes the joint
preparation gauge invariant. Conservation of the same Hamiltonian bounds
the probability that any specified plaquette is rough, for every model
time. This establishes the density hypotheses of the
[local roughness theorem](../particle_field_preparation/REFINEMENT.md#conditional-local-roughness-from-energy-density).
It gives neither a global Wilson-gap probability nor the missing local
low-mode/sign-response bound, and assigns no particle mass.

## 1. Independently normalized finite link polynomials

Use the actual fifteen-dimensional representation R of
`G=S(U(3)×U(2))`, its normalized Haar measure and source trace metric.
The [owning finite preparation](../particle_field_preparation/PROOF.md)
proves the Haar-ball and Casimir bounds used below. Define

```
s(U)=(1+Re Tr R(U)/15)/2=1-||R(U)-I||_F²/60,
L=n²,    Z1=integral s(U)^(2L) dU,    f_L=s^L/sqrt(Z1).
```

Here `0<=s<=1`, and f_L belongs to the matrix-word space of degree
at most L. The Frobenius norm is unnormalized. The source image lies in
the radius-sqrt15 sphere about zero in real ambient dimension 450.
The existing packing proof gives Haar mass at least `(1+8n)^(-450)`
to its identity ball of Frobenius radius 1/n. On that ball,
`s^L>=(1-1/(60n²))^(n²)>=59/60`. Integrating its square yields

```
Z1 >= (59/60)²(1+8n)^(-450) >0.                         (1)
```

This is a one-link normalization. No power of the number of links
appears in a local marginal or its probability bound.

For `0<delta<=1`, a link with represented operator distance
`||R(U)-I||>delta/4` has Frobenius distance at least that large. Thus
`s<=1-delta²/960`, and its normalized single-link probability obeys

```
q_L(delta) <= (60/59)²(1+8n)^450 exp[-n² delta²/480].     (2)
```

For `delta=2^(-r)`, the code uses only integer arithmetic:

```
b_link=2+450 ceil(log2(1+8n))-floor(n²/(480*4^r)),
q_L(delta)<=min(1,2^b_link).                            (3)
```

Indeed `(60/59)²<4`, the covering factor is bounded by its binary
ceiling, and `exp(-x)<=2^(-floor x)` for x>=0. The ceiling is exactly
`bit_length(8n)`. Doubling n eventually satisfies every fixed requested
tail because n² grows faster than log n. No floating exponential supplies
a proof premise.

## 2. Site-gauge twirling of the full density

On a periodic cubic lattice with V cells and E=3V positive-direction
links, form the normalized pure link product `F_L=product_e f_L(U_e)`.
The fermion input is the maximally mixed state on the entire
`2^(60V)` Dirac/M15 Fock carrier. Declare the source density

```
rho_V = integral_site gauges dg
        U(g) [|F_L><F_L| tensor I_Fock/2^(60V)] U(g)†.   (4)
```

This is a trace-preserving mixture of unitary conjugations. It is not
projection of the amplitude onto its singlet component. Its trace is one
without conditioning or a volume-dependent normalization. The independent
link product, uniform Fock marginal and gauge-averaging operation are
translation equivariant, so rho_V is translation invariant. It also
commutes with every site gauge transformation.

The physical gauge/reference distinction matters here. Commutation with
the gauge action makes rho_V invariant; it does not imply support only
on field singlets. Its nontrivial gauge-charge sectors are retained.
A complete joint preparation uses the finite canonical purification

```
|sqrt(rho_V)> = sum_ij sqrt(rho_V)_ij |i> |j*>.           (5)
```

The reference transforms in the conjugate of the source representation.
Since sqrt(rho_V) commutes with U(g), (5) is invariant under
`U(g) tensor conjugate(U(g))`. The source and reference can carry opposite
charges. The reference is an indispensable component of this declared
Gauss-invariant experiment; it is not silently projected out or called
the vacuum. The reference Hamiltonian is zero during observation, so
`H_source tensor I_reference` preserves the joint gauge constraint on
the same source clock. A field-only Gauss-singlet preparation is a
separate problem.

All source and reference spaces in (5) are finite. Gauge transformations
preserve the per-link word spaces, so (4) has total degree at most E L.
Its matrix entries are finite Haar integrals of native representation
coefficients. These can be computed by compact-group invariant-tensor
decomposition and finite Gram quotients. Their positive square root
gives a finite computable preparation target. A complete unitary taking
an invariant blank to (5), extended by identity on the orthogonal
complement of their at-most-two-dimensional span, is an admitted
apparatus primitive. The existing
[coherent-instrument compiler](../particle_observer_outputs/READOUT_CONTRACT.md)
provides the alternative implementation/error contract. This certificate
does not materialize the large density or execute a compiled pulse tape.
It claims no autonomous repair or cooling preparation.

Site gauges correlate adjacent links. The periodic density in (4) is
not replaced by a tensor product of independently endpoint-twirled link
densities. The latter appears only as the distinct single-link control
in section 6.

## 3. Uniform initial densities from actual invariant observables

For a plaquette with four unitary link factors, telescoping gives

```
||R(U_p)-I|| <= sum_(four factors) ||R(U_e)-I||.
```

Adjoints have the same distance. If the left side exceeds delta, at
least one factor exceeds delta/4. A four-event union bound and (2) give
`Prob(||R(U_p)-I||>delta)<=4 q_L(delta)` in the normalized product state.
This implication does not require four independent edges and applies
to repeated edges on a one-cell torus.

The plaquette cost is
`s_p=1-Re Tr R(U_p)/15=||R(U_p)-I||_F²/30`, between zero and two.
On its smooth event, `||R(U_p)-I||_F²<=15 delta²`; thus

```
<s_p>_initial <= delta²/2+8q_L(delta),
s_0=<S>_initial/V <=3delta²/2+24q_L(delta),
S=sum_p s_p.                                          (6)
```

Each plaquette cost and roughness event is gauge invariant, so the
density twirl in (4) preserves its expectation exactly. Electric
Casimirs are also gauge invariant. The inherited source-metric bound
on each degree-L link gives

```
<H_E>_initial/V <=3*(29/10)L²=(87/10)L²=:e_E.            (7)
```

Equations(6)–(7) are independent of V. They use finite word degree for
electric energy, rather than dividing a derivative estimate by a tiny
global normalization. The maximally mixed full Fock input is a permitted
state choice; no ground or thermal interpretation is assigned to it.

## 4. One action for all volumes and an all-time local bound

Fix a, polynomial order m, and the local coefficient kappa>0. On each
periodic parent use exactly

```
H_m=H_E+K_m-kappa sum_p Re Tr_M15(U_p)/15,
K_m=dGamma(h_m),    ||K_m||<=bV,
b=60(2+epsilon_m)/a,    epsilon_m=16(15/17)^(m+1).
```

This full-Fock norm bound retains all occupations. Constant-shifting
the action gives `H'_m=H_E+K_m+kappa S`. The source marginal is
translation invariant and that symmetry is preserved by the Hamiltonian.
Energy conservation and H_E>=0 yield

```
<S>_t/V <= s_0+(e_E+2b)/kappa.                         (8)
```

For each fixed orientation, translation symmetry makes its V plaquette
means equal. Its orientation sum is bounded by the total positive S,
so each particular plaquette mean is bounded by (8). No division by
three is used. From `R_p=1_{||R(U_p)-I||>1/40}<=48000 s_p`,

```
Prob_t(rough plaquette p)
 <=min(1,48000[s_0+(e_E+2b)/kappa])                    (9)
```

at every model time. This is the same source/state/clock energy argument
as the [conditional density lemma](../particle_field_preparation/REFINEMENT.md#conditional-local-roughness-from-energy-density),
with its initial density premises realized by (4).

The exact planner takes eta=2^(-b_target). It chooses delta=2^(-r) and
a dyadic q such that `3delta²/2<=eta/(4*48000)` and
`24q<=eta/(4*48000)`, then doubles n until (3)<=q. It fixes

```
L=n²,    kappa=96000(e_E+2b)/eta.                       (10)
```

Then (6)<=eta/96000, and (9)<=eta. Both L and kappa depend on target
accuracy, a and m, and are independent of V. They are not changed during
observation. Their large values are architectural inputs with no
experimental particle fitting. The construction supplies no physical
normalization or claim that these coefficients reproduce observed
couplings.

## 5. Finite compression and complete effects

For any K>=E L, let P_K be the total matrix-word projection and retain
the entire Fock carrier. The actual finite Hamiltonian is
`H_(m,K)=P_K H_m P_K`. P_K commutes with gauge transformations and
translations, contains the initial support and preserves positivity of
electric compression. Kinetic compression has the same bV norm bound;
the initial energy and magnetic expectations are unchanged. The finite
energy-conservation argument proves (9) directly for

```
E_rough=P_K R_p P_K,    E_smooth=I_code-E_rough.          (11)
```

These positive effects sum to identity, including on inputs outside the
chosen state's support. They may be realized by the standard complete
square-root instrument with a retained two-label pointer. Neither effect
is replaced by the roughness of a separately truncated link or Wilson
matrix. The canonical reference is untouched during observation, and
the full joint output is retained. No smooth outcome is accepted and
renormalized. Initial zero-support directions remain valid later input
and detector directions.

The finite matrices in (11) are defined by integrals of bounded
semialgebraic events and finite representation coefficients. This
analytic certificate bounds them without materializing their entries.
An approximate detector implementation must supply its operator or
complete-channel error; no numerical detector accuracy is assumed here.
Preparation errors must similarly be propagated into actual energy
allowances, or into finite-channel comparisons for a specified finite
observation, rather than discarded by calling a prepared state exact.

## 6. An exact native one-link normalization and charge control

For L=1, the native representation has five distinct irreducible blocks
of dimensions 6,3,3,2,1, with ten distinct blocks when their conjugates
are included. In normalized Haar character coordinates,

```
s(U)=1/2+(1/60)sum_(ten blocks) chi_r(U),
Z1=1/4+10/3600=91/360.
```

For a single open link with two independent endpoint gauges, the
complete twirl replaces the state in each irrep's `d_r²` coefficient
space by its scalar density. The trivial block has original probability
90/91; each of the ten nontrivial blocks has probability 1/910, spread
over its entire d_r² space. All 119 positive-support states are retained,
as is the zero-weight outside-support outcome. Original probabilities
sum to one. The exact native Casimir table gives mean electric energy
`33/9100` on this link. Its nontrivial charge probability is 1/91, making
explicit why an invariant mixed field state requires the conjugate
reference for a joint singlet interpretation.

This control verifies native normalization, representation and reference
conventions. It does not stand in for the correlated periodic site twirl
and does not certify the planner's roughness target at L=1.

## 7. Exact boundary and the next operator theorem

The result supplies a nonshrinking-in-volume family of local action
parameters, uniform initial densities and all-time local roughness
bounds. The family differs from the coherently projected filled-origin
preparation. No state equality or local-record comparison between those
two preparations is asserted. The source marginal in (4) is not a
field-only Gauss-singlet ground.

A rough plaquette is not a local Wilson low-mode. The
[missing operator estimate](../particle_field_preparation/REFINEMENT.md#the-missing-local-operator-estimate)
requires local weight/decay of `1_{|H_W|<1/2}` and a volume-uniform
complete-record Duhamel comparison between exact-sign and polynomial
evolution on this quantum state. Its right side must go to zero with
the proved local roughness and approximation error, at fixed nonzero
observation duration. The local density bound cannot replace that
operator theorem. Summing its roughness allowance over the entire
lattice restores a volume factor.

An interacting relativistic limit additionally requires estimates
uniform under spacing, polynomial order, state and coupling refinement,
with surviving connected interactions and the correct chiral content.
The existing [native spin density theorem](../native_repair_spin_limit/PROOF.md#exact-conservation-and-the-energy-estimate)
and [ground/thermal construction](../spin_reference_limit/PROOF.md#charge-density-from-source-energy)
offer a local proof pattern on their own two-flavor action; their
spectral conclusions are not imported into this coordinate gauge source.
