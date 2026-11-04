# Finite modular flow across a spatial boundary

The complete interacting plaquette has a faithful reference and a modular
clock, but that clock moves local observables outside every proper site
region. Exact certificates cover all fourteen regions and all 256 occupation
states. The exclusion also holds for the literal prepared source reference,
using its certified logarithm error. A separate
[Lean theorem](../../Lean/QFT/FiniteModularInclusion.lean) proves that any
finite-dimensional half-sided invariant algebra is invariant in both time
directions. A nontrivial null-horizon construction therefore requires a
controlled limiting inclusion, rather than exact finite compression.

## Dimension theorem

Let `S` be a finite-dimensional vector subspace and `T` an invertible linear
map on its ambient space. If `T(S) subset S`, the restriction `T|S:S->S` is
injective. Equal finite dimensions make it surjective, so `T(S)=S` and
`T^-1(S)=S`. No compactness, continuity or numerical approximation is used.
The ambient space itself need not have finite dimension.

Apply this to the vector space underlying a finite matrix subalgebra `N`
and the modular automorphism

```
alpha_t(A) = rho^(it) A rho^(-it),   rho > 0.
```

If `alpha_t(N) subset N` for every `t>=0`, each such inclusion is equality;
inverting it supplies invariance for every negative parameter as well.
Strict half-sided compression is impossible on this finite algebra.
`inverse_preserves`, `image_eq` and `half_sided_is_two_sided` formalize
this argument. The theorem does not assume the desired inclusion, then
declare its physical interpretation: the following source calculation
tests that inclusion and rejects it for the declared spatial cuts.

This dimension argument does not exclude a nontrivial half-sided inclusion
after a limit changes the algebra. It also does not exclude proving a
finite stress Ward identity by another route. No identification of a
modular-generator difference with a positive null-translation generator is
made here.

## The same complete interacting source

The source is the existing open `2 x 2 x 1` Weyl/contact action with
`c=1/8`, `U=1`, eight CAR modes, and its Gibbs reference at `beta=1/32`.
The action, temperature, spatial dictionary, native controls, observer
scheduler and record memory are specified services. Their specification is
available to the observer. No physical calibration follows from reading it.

For a nonempty proper site subset `R`, let `N_R` be the full local CAR
algebra, embedded by the signed Fock factorization. It contains every
regional number sector and odd operator. All `A in N_R` commute with the
number projection `n_q` of every exterior mode `q`. In particular,
`<p|A|q>=0` if the one-particle mode `p` is inside and `q` is outside.
This separation tests the entire algebra, without sampling an operator basis.

Every proper region has a crossing kinetic link. Choose linked modes `p`
inside and `q` outside, so `|H_pq|=1/16`, and use the norm-one local probe
`n_p`. On the complete Fock carrier,

```
[H,n_p]_(p,q) = -H_(p,q),
[n_q,[H,n_p]] != 0.
```

The exact Gaussian-integer checker constructs both commutators on all
256 states. Thus the modular derivative `-i beta [H,n_p]` lies outside
`N_R`. The one-particle matrix element witnesses an identity of complete
operators; it is not a restriction of the reference or dynamics to a
one-particle sector.

The signed source row sum gives `||H-I||<=3`. Unitary conjugation and
twice integrating its second derivative give

```
||alpha_t(n_p)-n_p+i beta t [H,n_p]|| <= 18 beta^2 t^2.
```

No exponential majorant is needed: the real-parameter unitary flow
preserves the operator norm of each derivative. Taking the forbidden
matrix element yields, for every `B in N_R` and `t>=0`,

```
||alpha_t(n_p)-B|| >= beta t/16 - 18 beta^2 t^2.
```

At `beta=t=1/32` the right side is exactly `23/524288`, greater than
`4.38e-5`. Every one of the fourteen proper site regions therefore fails
positive-time inclusion with a uniform quantitative margin. All statements
also survive tensoring with an untouched reference and replacing the
regional algebra by `N_R tensor B(H_reference)`, since the same source
matrix element between equal reference vectors separates the algebras.

## Transfer to the literal prepared reference

The [certified source preparation](../certified_reference/README.md)
retains every classical column record and export abort. Its successful
field block `sigma` is positive and subnormalized. The retained certificate
proves, on the complete field carrier,

```
||-log(sigma) - (-log(rho_Gibbs))|| <= d,
d < 11/10^10.
```

The checker verifies the inherited certificate's source hashes, independently
recomputes this logarithm bound from its rational complete-output error and
the polynomial design theorem, and retains the exact rational value of `d`.
It does not repeat the 65,535-node native control proof; that proof and its
replay command remain with the bound source certificate.

Duhamel's identity for two self-adjoint generators gives

```
||sigma^(it) n_p sigma^(-it)-alpha_t(n_p)|| <= 2 |t| d.
```

Consequently the same fourteen cuts have source-reference distance at
least `23/524288-d/16 > 1/25000` at `t=1/32`. Faithfulness of `sigma` is
certified separately by its positive full-field minimum weight. Multiplying
`sigma` by any positive scalar leaves this modular action unchanged; no
success conditioning or missing-outcome interpretation is required. The
flow considered here is the field success block's modular automorphism,
not a claim about a faithful modular algebra for the entire retained
record-and-abort history. Neither this flow nor the Gibbs comparison is
identified with the irreversible repair channel.

## Controls and scope

Removing `U` leaves the crossing-current witness unchanged: the exclusion
is caused by hopping across the chosen boundary and does not require the
contact interaction. Removing every crossing hopping bond makes the Gibbs
Hamiltonian exactly `H_R tensor I + I tensor H_outside`. The signed integer
tensor check proves this identity for every cut, so its Gibbs modular flow
preserves both regional algebras at every parameter. This deletion control
concerns the Gibbs reference, not a claim that a truncated polynomial in a
disconnected Hamiltonian factorizes. Zero speed with the onsite interaction
retained also passes the tensor-sum check. At `beta=0` the Gibbs modular flow
is identity. The empty and full algebras are invariant.

This closes the finite feasibility test for an exact half-sided modular
attachment on the specified source and all its proper site cuts. It does
not fill the common Einstein register's null-translation supply, identify
geometric diamonds, derive area/edge entropy or Newton normalization, or
prove a common interacting continuum. The continuum route needs a sequence
of source algebras, references and inclusion errors controlled in the
limiting representation. A finite boundary-current failure alone supplies
no lower bound uniform in every allowed continuum scaling.

## Reproduction

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.finite_modular_inclusion --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/finite_modular_inclusion
cd Lean
lake env lean QFT/FiniteModularInclusion.lean
```

`--write` regenerates the source-bound rational receipt. Its fresh-process
replay computes the full source commutators, all cut controls and both
positive-time bounds. Independent tests assemble the field through the
separate full-Fock evaluator, use the complete signed CAR restrictions,
and compare with full-spectrum modular evolution. Floating evaluations
test the rational bounds; they are not used to establish them.
