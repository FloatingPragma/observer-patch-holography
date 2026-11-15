# The interacting field on complete Fock space

This construction lifts the supplied [Weyl/contact action](../interacting_field/README.md)
from its two-particle sector to the complete finite CAR algebra. The same
completed source hopping, contact and fermionic-swap words act on every
occupation number, including coherence between number sectors and arbitrary
untouched reference columns. An exact locality bound also constructs its
fixed-spacing thermodynamic local dynamics. This is an interacting lattice
field in three spatial coordinates and Hamiltonian model time. The lattice,
action, coherent control services and clock convention are supplied; no
Lorentzian continuum, taste removal or native physical-action selection is
established.

## One action and every particle number

For an open rectangular box `Lambda` of `V` sites, use two CAR modes per site
and the full `2^(2V)`-dimensional Fock space:

```
H_Lambda = sum_<x,y> (a_x* h_xy a_y + a_y* h_xy* a_x)
           + U sum_x n_x,up n_x,down,
h_(x,x+e_j) = i c sigma_j/2,   c>=0, U>=0.
```

There is no particle-number cutoff, Gaussian/Wick closure or momentum-pocket
projection in this definition. The complete space is the direct sum of all
exterior powers, and `[H,N]=0` on that entire space. A spatial bond contains
two hopping terms of norm `c/2`; hence its norm is at most `c`. If `E` is the
number of spatial edges, positivity of the contact terms gives

```
-c E I <= H_Lambda <= (c E + U V) I,     E <= 3V.
```

This establishes extensive stability without using a fixed particle number.
Every finite-temperature Gibbs density on this finite space is faithful.
It does not prepare that reference or establish the regional product
conditions needed for a Gibbs-preserving retraction.

The empty Fock vector has declared energy zero, but is generally **not** the
interacting ground state. The chemical potential, physical vacuum and
normal-ordering convention have not been selected by the native process.
Including this vector therefore does not by itself identify a particle mass
or remove the physical vacuum-comparison premise.

## Source execution and a complete error budget

A completed two-mode hopping or contact word acts on all four local Fock
inputs. Move a nonadjacent mode next to its partner by FSWAPs, apply the word,
then undo the route. The minus sign on double occupation is essential.
The full-Fock executor performs these routed local contractions literally,
without dropping vacuum phases, number sectors or small numerical off-sector
entries. An optional bounded `source_word` emits the actual concatenated
pulse tape. Its direct source execution agrees with local contraction, and
its pulse/load/export totals agree with the existing resource counter.
Failure branches of the underlying valid-code services remain zero in the
ideal compiler contract; numerical source leakage is a separate diagnostic.

For the symmetric split, let `B=cE+UV`, the sum of local generator norm
bounds. The third derivative of either the exact exponential or the
symmetric product has norm at most `B^3`: expand the product derivative and
use unitarity and the multinomial identity. Their first two derivatives
agree at zero. Taylor's integral remainder therefore bounds one step by
`B^3 |t|^3/3`; telescoping `n` unitary steps gives

```
||U_source,ideal(t,n)-exp(-itH)|| <= B^3 |t|^3/(3n^2).
```

The half-diamond distance of the unitary channels is no larger than this
operator bound. Summing segment errors and using instrument contractivity
gives the complete finite-history bound, with all outcomes and arbitrary
references retained. There is no success-probability division or record-count
multiplier. The **full-Fock bound depends on volume**; the volume-uniform
bound proved on the `N=2` carrier must not be reused here.
Floating source contraction and hardware errors are outside the ideal bound.
The existing local-residual sum is retained as a diagnostic, not an outward
interval certificate.

`FieldGeometry` and resource counts need no global Fock allocation. The
explicit evaluator is separately bounded to at most 16 modes and `2^20`
state/reference entries. Its three-dimensional `2x2x2` control executes a reference-entangled state
on the complete 65,536-dimensional Fock carrier and compares it with an
independent sparse CAR exponential, and retains both occupation-pointer
outputs. Larger geometries are declared schedules and mathematical families,
not claims to have materialized exponential-size state arrays.

## Local thermodynamic dynamics at fixed spacing

Finite particle count is unnecessary for the locality argument. Group the
two modes of each site, use spatial Manhattan distance, and keep their CAR
parity. Disjoint local operators commute when at least one is even; two odd
operators require the graded anticommutator formulation. The calculation
below uses ordinary commutators only under the stated commuting hypothesis.

The finite-support commutator-chain iteration is the standard fermionic
Lieb--Robinson argument of Nachtergaele, Sims and Young,
[Theorem 3.1 and equations 3.28--3.35](https://arxiv.org/html/1705.08553v3#S3).
Its use here is a concrete model mapping: the onsite contact terms are removed
by a product of local interaction-picture unitaries. Each resulting bond is
time dependent but retains its two-site support, even parity and norm bound
`c`. Thus the estimate is uniform in `U`. It does not replace fermionic
locality by tensor commutativity of odd fields.

For initially commuting disjoint `A_X,B_Y`, the first chain link has at most
`6|X|` choices and each subsequent bond at most 12 choices. Each commutator
contributes a factor two. A chain spanning distance `d>=1` has at least `d`
bonds. The resulting norm bound is

```
||[tau_t(A_X),B_Y]||
 <= ||A_X|| ||B_Y|| |X| sum_(n>=d) (24c|t|)^n/n!
 <= ||A_X|| ||B_Y|| |X| 2^(-d) [exp(48c|t|)-1].
```

The same chain count bounds the iteration remainder by a factorial tail,
which tends to zero. The estimate is uniform in finite volume and particle
number. It is also uniform in an untouched reference dimension: the proof
uses operator norms and the same field-support chains after tensoring by an
identity on the reference.

Take `||A||<=1`, supported on `|X|` sites in the cube of radius `r`, and a
strictly padded cube of radius `R>r`. In the Duhamel difference between its
dynamics and any larger cube, exterior terms commute exactly with the
smaller-box evolution. Only the `6(2R+1)^2` crossing bonds remain. Their
distance from `X` is at least `R-r`; each has norm at most `c` and is even,
so this comparison applies to arbitrary local `A`, including odd `A`.
Consequently

```
||tau_t^larger(A)-tau_t^R(A)||
 <= min(2, 6c |X| (2R+1)^2 |t| 2^(-(R-r)) [exp(48c|t|)-1]).
```

For every bounded model-time interval this tends to zero as `R` grows.
Local dynamics is therefore norm Cauchy; its limit preserves products and
adjoints. Applying the same limit at negative time gives the inverse, and
density of the local algebras gives a strongly continuous automorphism group
on the full quasi-local CAR algebra. This proves existence of the same
interacting **fixed-spacing** local dynamics, without selecting an
infinite-volume state or taking spatial spacing to zero.

`local_boundary_bound` evaluates an outward exact rational upper bound for
this expression. It bounds the exponential by a Taylor polynomial plus a
geometric tail after scaling to an argument at most one half. Its bounded
arithmetic domain (radius at most 4,096 and enclosed exponent at most 64)
is distinct from the mathematical arbitrary-radius limit. Larger exponents
return the trivial bound two.
This is a local-observable comparison, not a global-state distance between
different-sized systems. The finite-volume complete-history split bound above
is a separate result; no unrestricted thermodynamic detector-instrument
error bound is inferred from a single observable estimate.

## The limiting dynamics remains interacting

For neighboring sites, let `A=a_0* a_2`, an even hopping observable, and let
`n_1,n_3` be the opposite-spin occupations at those sites. Contact deletion
changes its local Hamiltonian derivation by exactly

```
[H_U-H_0,A] = U A (n_1-n_3).
```

This nonzero quartic operator persists inside every containing box and in
the infinite-volume local derivation. It is not a one-body generator
renamed as an interaction. The complete source-executed non-Gaussian
four-point witness on the `N=2` subspace remains a valid subspace witness in
the full theory; the full-Fock implementation separately checks its exact
`N=2` restriction against the existing executor. No taste or chirality has
been discarded. Source-selected field coefficients, physical reference and
clock preparation, the chiral continuum limit and particle assignments
are additional obligations.

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest -q -W error model/full_fock_field
```
