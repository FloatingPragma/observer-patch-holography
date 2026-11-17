# Declared preparation works on trained regions and fails common-region consistency

The source-defined entropy response gives a unique global coframe for the
27 trained x-oriented region protocols: linear residual `4.80e-14`.
This is an entropy-based preparation rule, with no Poisson or Einstein data.
Actual finite-epsilon entropies of the exact neutral-pair mixture were then
computed in regional Fock spaces. The trained residual decreases by factors
`3.968,3.968` when epsilon halves, consistent with its expected second-order
remainder.

The held-out linear residual is `3.92014`. An overdetermined solve using all
81 orientations still leaves `2.57466`; rank is 27. Reversed-centre regions
and an unused smaller-radius protocol also change entropy. One global scalar
coframe does not make these complete readouts stationary on this small model.
At the three executed epsilons, selector Shannon curvature is large enough
that the held-out finite residual is not yet dominated by its first-order
term; the nonzero first derivative is given by the independently retained
matrix, not inferred from a supposed factor-of-two finite scaling.

This is a failure of the specific 3 cubed, R=1, one-/two-cell expected-volume
preparation architecture. There is no separation `R << ell`, and the chosen
pair is not a rotationally averaged matter ensemble. It is not a continuum
regulator no-go, nor a general obstruction to tensor geometry.

## Complete records and observer algebra

Each preparation branch is either the full filled Slater vacuum or that
vacuum with exactly one particle and one hole. Both branches have identical
total charge. Their regional density matrices are mixed before taking the
entropy. The Gaussian of the mean covariance is not used.

The external preparation flag is retained outside the field-region algebra.
The region selector is included in that algebra, so its full Shannon entropy
is part of the computed quantity. Including the preparation flag as well
instead gives a different classical-quantum entropy with an explicit
`H(epsilon)` term; this alternative is calculated in `INSPECTION.json` and
is not substituted for the selected regional marginal. No flag is erased.
Controller work/native implementation costs are not derived.

Expected proper volume is exactly 1.5 for the trained protocol. Individual
selector branches have different volumes. The construction does not realize
literal branchwise fixed proper volume.

## Response decomposition and geometry diagnostics

The Frobenius norms of the full response components are:

- actual coframe/covariance response: `1.60182`;
- selector-probability response: `1.42085`;
- selector Shannon first derivative: zero at weight 1/2.

The uniform eigenvalue is `0.217916 = 0.429724 - 0.211808`.
The zero mode was retained. Constraining its mean to zero leaves trained
residual `0.700573`. At the largest epsilon and central region the selector
Shannon contribution alone is `-0.00286652`; this explains most of the
second-order trained residual.

A geometric comparison performed only after the solve gives a 0.785 relative
remainder against a proportional lattice-Poisson response. The fitted
coefficient is 30.299, versus conditional tower value 46.318; the fit was
never used to prepare sigma. No lapse law or physical G match is supplied.

## Next constructive attempt

Use the same regulator with an isotropic exact-charge excitation ensemble,
larger digital balls, and an independently estimated proper-volume shape
tangent. Compare a common low-q coframe across unused radii, retaining both
shape-discretization errors and the zero-mode residual. A general tensor
coframe is a separate candidate if isotropic matter does not remove the
anisotropy. This failed selector architecture remains frozen.
