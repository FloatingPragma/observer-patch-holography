# B13 finite-effect closure boundary

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

## Boundary and continuation

This is a substantial partial result for B13, not closure.  It does not
derive affinity, an effect algebra, Gleason/Busch hypotheses, source-produced
effects, noncontextuality across richer contexts, a public readback, or a
physical quantum instrument.  Issue #702 must obtain one of those stronger
bridges; adding only continuity to binary-projector normalization cannot work.
