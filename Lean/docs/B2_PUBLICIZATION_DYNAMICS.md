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
trace-preservation certificate. `Dynamics/ChoiCPTP.lean` supplies the separate
formal complete-positivity predicate and proves this map CPTP from those
Kraus and trace identities.

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
nonzero-rate single-collar kernel is exactly the partition commutant, and the
trace residual identity

\[
  \operatorname{Tr}\bigl(X^\dagger L_c(X)\bigr)
  =-\gamma_c\,\lVert X-E_c(X)\rVert_{\mathrm{HS}}^2.
\]

`multiCollarNoCancellation_of_nonneg` uses this identity to prove that a
finite sum of these generators cannot cancel when every rate is nonnegative.
Zero-rate terms are included. With strictly positive rates,
`multiCollarGenerator_eq_zero_iff_stableIntersection_of_pos` identifies the
kernel of the sum with the intersection of the partition commutants. The
pinchings need not commute. The signed-rate theorem retains its
explicit no-cancellation premise: two identical nontrivial pinchings with
rates `+1` and `-1` cancel even when each has a nonzero residual.

`Dynamics/ChoiCPTP.lean` defines complete positivity through positivity at
every finite matrix amplification and proves the Kraus criterion. It proves
partition averaging, partition pinching, and each relaxation channel with
`gamma >= 0` and `t >= 0` CPTP. The relaxation map is identified with the
existing `publicRelaxTime`. An exact transpose-map Choi counterexample
distinguishes ordinary positivity and trace preservation from complete
positivity.

The exact-rational companion certificate in
`code/publicization/b2_exact_rational_certificate.py` checks the finite
three-dimensional matrix-unit instance and its negative controls.

## Boundary

The no-cancellation proof concerns the stated trace-orthogonal pinching
generators with nonnegative declared rates; it does not extend to arbitrary
positive maps or signed rates. The package supplies no source-selected
generator, physical rate or clock calibration, or physical readout. The
literal operator-exponential theorem assumes a bounded endomorphism of a
complete real normed space. The finite CPTP theorems supply no physical
prediction by themselves.
