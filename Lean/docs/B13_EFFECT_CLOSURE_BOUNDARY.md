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

`QFT/SourcePhaseLiftBridge.lean` starts from the two noncommuting real
projection effects earned by the committed source payload,
`recordProjector` and `conjProjector 3`. Subtracting their normalized complex
commutator from `I/2` gives the phase lift

\[
\frac I2-\frac{2\sqrt 3}{3}i(QP-PQ)=\rho_{Y+},
\]

exactly. Lean proves that this phase lift is a projection, that its binary
complement completes a context, and that it distinguishes the two pure
Pauli-Y states which the native real web cannot distinguish. The record
projector, the realized rotated projector, and the phase lift separate all
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

## Boundary and continuation

This is substantial bounded progress for B13, not closure. The phase lift is
an exact operator-algebra target, not a source operation, public effect,
instrument, or outcome receipt. The package does not derive affinity, a
source effect algebra, Gleason/Busch hypotheses, noncontextuality across rich
coexistent-effect contexts, operational coarse graining, or a physical
quantum instrument. Issue #702 remains open for a phase-sensitive
source-produced instrument with common-preparation outcomes and the
operational additivity bridge; adding only more phase-free real contexts or
continuity to binary-projector normalization cannot work.
