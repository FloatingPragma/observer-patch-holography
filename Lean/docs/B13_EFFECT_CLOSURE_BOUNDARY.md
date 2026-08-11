# B13 finite-effect and phase-lift boundary

`EventAlgebra/FiniteEffectClosureBoundary.lean` records the exact boundary
between the available celestial binary-projector data and a Born-rule
representation.

## Exact negative result

On the C1 unit celestial sphere, Lean defines

\[
F(n)=\frac{1+n_z^3}{2}.
\]

It proves that `F` is continuous, takes values in `[0,1]`, and satisfies
`F(-n)=1-F(n)`.  It also proves that no vector `q` represents this function as
`(1+q·n)/2` on every unit direction.  The contradiction is exact and uses the
three coordinate axes together with `(3/5,0,4/5)`.

Therefore continuity, probability range, and normalization on every
antipodal binary context do **not** force affine/Born form, even when the
entire celestial sphere is available.

## Exact positive result after affinity

If the centered response is already supplied in affine form, bounds on any
dense subset of unit directions extend by continuity to the whole sphere.
Those bounds force the coefficient vector into the closed unit ball.  The
module exposes both the centered and probability-interval versions.

## Exact phase-lift target and real-closure no-go

`QFT/SourcePhaseLiftBridge.lean` starts from two noncommuting real algebraic
projector candidates obtained by applying the declared two-dimensional S3
representation to source-realized gauge labels: `recordProjector` and
`conjProjector 3`. Subtracting their normalized complex
commutator from `I/2` gives the phase lift

\[
\frac I2-\frac{2\sqrt 3}{3}i(QP-PQ)=\rho_{Y+},
\]

exactly. Lean proves that this phase lift is a projection, that its binary
complement completes a context, and that it distinguishes the two pure
Pauli-Y states which the algebraic real web cannot distinguish. The record
projector, the algebraically rotated projector, and the phase lift separate all
two-by-two matrices on every fixed-trace slice; this statement does not need
positivity or Hermiticity.

The positive algebraic closure is paired with a stronger negative result.
Even a generous phase-free grammar containing every native real outcome,
complements, real sums and scalings, and pullback by arbitrary real Kraus
matrices remains real symmetric. Every effect in that closure gives equal
Born weight to the opposite Pauli-Y states, and the complex phase lift cannot
belong to it. `sourcePhaseLift_boundary_summary` packages the complex-algebra
membership, exact projection, fixed-trace tomography, real-closure exclusion,
and the receipt fact that the committed run has outcomes only for its native
diagonal context.

## Conjugation gauge and post-hoc count-oriented completion

`QFT/ConjugationGauge.lean` makes the completion ambiguity structural.
Entrywise conjugation preserves Hermiticity, positivity, states, effects,
and events, equals transposition on Hermitian matrices, and conjugates the
trace. Simultaneously conjugating state and effect therefore preserves the
real part of their Born weight; conjugating the state alone does so only for
a conjugation-fixed effect. The
complexified real web is fixed pointwise, and the two candidate phase
completions \(\rho_{Y+}\) and \(\rho_{Y-}\) form one conjugation orbit.
Hence no conjugation-fixed frame separates any state from its conjugate,
and the choice between the two candidates is the choice of an orientation
of a two-element torsor. These statements bound the web's ambiguity from
below; the full Pauli-Y coordinate, magnitude included, is separately
hidden from the web.

`Thermodynamics/RepairCurrentOrientation.lean` mirrors the committed
payload `docs/REPAIR_CURRENT_PAYLOAD.json` of `oph-physics-sim` (schema
v3, sha256 `7f8ea7ef9c92a50e23207c2fe85d09ed2bce1c1aa539ae9914a9b9edd0df26d6`,
recounted post hoc from the retained bundle of a run executed under an
unrelated, locally hash-pinned B12 contract
`runs/b12_prereg_16k_20260806` with report-total, alphabet,
strict-integer, and npz support/recount cross-checks). Kernel decides
certify the totals, the designated pair `3 -> 4` with `1343` forward and
`0` backward counts, the designated cycle `(3,4,5)` with products
`1239691068` and `0`, the maximality and unique-orbit attainment of both
designations and lexicographically least tie breaks, the exact normalized
selected-cycle products `9391599/57188378` and `0`, the reversal-oddness of
the bit under transposition, and the synthetic symmetric-control degeneration.
Neither this statistic nor its designation rule was preregistered, and the
payload records that it is ineligible as validation.

`QFT/SourceOrientedCompletion.lean` composes the two under a declared
typed convention with an explicit applicability condition (the strict
designated-cycle inequality, kernel-decided for the committed table): the
orientation bit selects the torsor element. Under the declared pairing the
selected completion coincides with the phase lift; the coincidence is made
by the pairing, and the opposite pairing is equally admissible. The frame
completed by either torsor element identifies every certified state, so
the convention is outcome-robust in both branches, and the composition
with `finite_busch_gleason` pins the state representing any additive
effect valuation by its three oriented Born weights.

## Boundary and continuation

This is substantial bounded progress for B13, not closure. The phase lift
is an exact operator-algebra target, and the orientation is a post-hoc
raw-count diagnostic consumed under a declared convention; neither
is a source operation, public effect, instrument, or outcome receipt. The
package does not derive affinity, a source effect algebra, Gleason/Busch
hypotheses, noncontextuality across rich coexistent-effect contexts,
operational coarse graining, a physical quantum instrument, or any
arrow-of-time claim. Issue #702 remains open for a phase-sensitive
source-produced instrument with common-preparation outcomes (including the
y-magnitude readout of an unknown state), the operational additivity
bridge, and a freshly preregistered validation of the oriented readout; adding only
more phase-free real contexts or continuity to binary-projector
normalization cannot work.
