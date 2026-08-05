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

`Dynamics/PublicAutomorphism.lean` strengthens the original supplied-flow
obstruction in two steps. `publicStarAutomorphism_is_labelPermutation` proves
that every star-algebra automorphism of the finite public function algebra is
uniquely pullback by a permutation of the active labels. The proof classifies
the images of the orthogonal record-basis idempotents pointwise and does not
assume positivity as extra data. `ContinuousPublicStarFlow.toAut_eq_refl`
then proves that every pointwise-continuous real-parameter group of arbitrary
public star automorphisms is the identity. Operator continuity makes the
classified label locally constant, and connectedness of the real line fixes
it at its value at zero.

`ContinuousPermutationFlow.toPerm_eq_refl` is the smaller supplied-label
interface and follows the same connectedness/discreteness principle.

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

The public classification is complete for the exact finite function
algebra, and the continuous obstruction starts from arbitrary pointwise-
continuous public star automorphisms rather than a supplied permutation flow.
Transport to the ordered matrix subalgebra uses the B1 star-algebra
equivalence; no source law or physical clock is introduced.

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
group is an open obligation.

No source rule selects a stochastic kernel, Hamiltonian, time scale, or
physical clock. These results establish the finite algebraic dichotomy under
declared dynamical data and do not constitute a physical prediction.

## Closure classification

The public stochastic classification, arbitrary-public-automorphism
classification and continuity obstruction, simple private-block innerness
theorem, fixed-Hamiltonian unitary flow, and von Neumann differential equation
are attained without admissions. Issue 679 is open only for the
arbitrary finite central-block decomposition and the converse continuous-
group-to-Hamiltonian theorem.
