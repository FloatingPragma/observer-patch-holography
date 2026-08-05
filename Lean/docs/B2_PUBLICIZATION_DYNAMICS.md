# B2 finite publicization dynamics boundary

Issue `#678` asks for a finite publicization map, an exact relaxation law,
and the stable algebra of the corresponding generator without treating a
displayed Kraus sum as a proved physical channel.

## Proved in Lean

`EventAlgebra/PartitionAverageCP.lean` constructs an explicit finite
Parseval-frame Kraus family for partition averaging. The family satisfies

\[
  \sum_{i,a,b} K_{iab}^{\dagger}K_{iab}=I,
  \qquad
  E_{\rm av}(X)=\sum_{i,a,b}K_{iab}XK_{iab}^{\dagger}.
\]

Together with the existing trace theorem, this is a normalized Kraus and
trace-preservation certificate. The module does not define a
complete-positivity predicate or a bundled CP/CPTP channel.

`EventAlgebra/TwoScalePublicRepair.lean` proves the exact linear relaxation

\[
  R_a(x)=E x+a(x-E x)
\]

for an idempotent map `E`. It fixes the public component, scales the residual,
and satisfies `R_a R_b = R_{ab}`. The specialization
`a = exp(-gamma t)` is an additive-parameter semigroup.

`Thermodynamics/PoissonizedRepair.lean` bundles the same closed form as a
linear map, proves its initial value, semigroup, fixed-space, and generator
flow identities, and records `L = gamma(E-I)`. The algebraic identities hold
for real parameters. Calling `gamma` a Poisson rate and `t` forward time
requires `gamma >= 0` and `t >= 0`; no physical clock calibration follows.

`Thermodynamics/PoissonizedRepairOperatorExp.lean` closes the literal
operator-exponential interface. For a bounded idempotent endomorphism of a
real Banach space it proves, using Mathlib's Banach-algebra exponential,

\[
  \operatorname{exp}\!\left(t\gamma(E-I)\right)
  =E+e^{-\gamma t}(I-E).
\]

The module also proves the generic identity
`exp(a P) = 1 + (exp(a)-1) P` for a Banach-algebra idempotent and gives exact
bridges back to the existing algebraic generator and repair map.

`Dynamics/ConditionalExpectationGenerator.lean` specializes `E` to partition
pinching. It proves the displayed projector-rate dissipator identity, that a
nonzero-rate single-collar kernel is exactly the partition commutant, and that
a multi-collar kernel is the intersection of commutants under an explicit
no-cancellation premise.

The exact-rational companion certificate in
`code/publicization/b2_exact_rational_certificate.py` checks the finite
three-dimensional matrix-unit instance and its negative controls.

## Boundary

The package does not define a CP/CPTP predicate, derive multi-collar
no-cancellation from positivity, select a source generator, calibrate a rate
or clock, or attach a physical readout. Those are nonclaims; the CP/CPTP
interface becomes mandatory only if either label is asserted. The literal
operator-exponential theorem assumes a bounded endomorphism of a complete real
normed space. The package emits no prediction.
