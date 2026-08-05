# B3 finite public/private dynamics boundary

## Formal receipts

`Dynamics/PublicMarkov.lean` uses the active-label function algebra supplied
by `EventAlgebra/PublicRecordAlgebra.lean`. Its positivity convention is
explicit: a complex record function is nonnegative when every coordinate is
real with nonnegative real part. A complex-linear map is positive when it
preserves that cone and unital when it fixes the constant function one.

Under those conventions,
`activeRecord_positive_unital_iff_stochastic` proves an equivalence. Every
positive unital map has a unique real kernel in the record basis, each kernel
entry is nonnegative, every row sums to one, and the map acts by

```text
(T f)(x) = sum_y K(x,y) f(y).
```

Conversely, every real row-stochastic kernel defines such a positive unital
map. This is the Heisenberg observable convention. It does not assert that
the kernel is doubly stochastic or that a Schrödinger-picture state map uses
the same row convention.

`ContinuousPermutationFlow.toPerm_eq_refl` proves that every coordinatewise
continuous real one-parameter permutation flow on a discrete record set is
the identity. The proof uses connectedness of the real line and discreteness
of the label set. The group law is part of the interface, although the
continuity obstruction does not need it.

`Dynamics/PrivateInner.lean` proves three private-block results:

1. `finitePrivateStarAutomorphism_inner` states that every star-algebra
   automorphism of the continuous endomorphisms of a finite-dimensional
   complex Hilbert space is conjugation by a unitary linear isometry.
2. `hasDerivAt_innerExponentialFlow` gives the exact commutator derivative of
   exponential conjugation. `hasDerivAt_vonNeumannFlow` specializes it to
   `d rho / dt = -i[H,rho]` in complex time, while
   `hasDerivAt_realVonNeumannFlow` proves the real-parameter equation without
   identifying that parameter with a physical clock.
3. `hamiltonianPropagator_mem_unitary` proves that `exp(-itH)` is unitary for
   every real `t` when `H` is self-adjoint.

The innerness proof specializes Mathlib's
`StarAlgEquiv.eq_linearIsometryEquivConjStarAlgEquiv`. Finite dimensionality
discharges its continuity premise. The exponential and derivative results use
Mathlib's Banach-algebra exponential, unitary, and derivative interfaces.

## Exact boundary

The public stochastic theorem acts on the active-label function algebra. The
star-algebra equivalence from B1 identifies that function algebra with the
projective-partition public subalgebra. A bundled ordered or positive map on
the matrix subalgebra is not constructed here. Complete positivity is not
part of the theorem interface.

The continuity obstruction is proved for label permutations. A theorem that
classifies every star-algebra automorphism of a finite function algebra as a
label permutation is not included. Consequently the receipt supports
continuous reversible dynamics once the reversible public action is supplied
as a label-permutation flow; it does not silently replace an arbitrary
automorphism flow by such data.

The private innerness theorem covers one full endomorphism algebra, hence one
simple matrix block. It does not classify an arbitrary finite-dimensional
complex C-star algebra as a direct sum of full matrix blocks. The pinned
Mathlib tree exposes no applicable block-decomposition theorem. Proving
the decomposition, transporting the topology and star structure, and
handling central-block permutations form the remaining route for the
per-central-block statement in issue 679.

For a continuous one-parameter group of private automorphisms, the innerness
theorem supplies a unitary implementer for each group element. It does not
supply a coherent continuous choice of implementers, remove the phase
ambiguity, or derive one time-independent self-adjoint generator. The pinned
Mathlib tree exposes no packaged finite-dimensional Stone theorem for this
lift. The explicit fixed-Hamiltonian flow has the unitary and von
Neumann receipts, but the converse from an arbitrary continuous automorphism
group remains open.

No source rule selects a stochastic kernel, Hamiltonian, time scale, or
physical clock. These results establish the finite algebraic dichotomy under
declared dynamical data and do not constitute a physical prediction.

## Closure classification

The public stochastic classification, public continuity obstruction, simple
private-block innerness theorem, fixed-Hamiltonian unitary flow, and von
Neumann differential equation are attained without admissions. Issue 679 is
open for classification of arbitrary public star automorphisms as label
permutations, the arbitrary finite central-block decomposition, and the
converse continuous-group-to-Hamiltonian theorem.
