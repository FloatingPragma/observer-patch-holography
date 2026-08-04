# B7 conditional-history boundary

## Attained theorem packet

`InformationProjection/PathGibbs.lean` proves conditional finite-history
mathematics. Given a positive normalized projected weight, a positive
reference weight, a displayed exponential log relation, and moment matching,
it supplies the information-projection Pythagorean identity and minimizer
inequality. A supplied modal path minimizes the effective action. For the
normalized inverse-noise family it proves both a fixed positive-gap tail bound
and the actual limit in which the total mass of every strict nonminimizer
tends to zero.

`Variational/DiscreteEulerLagrange.lean` and
`Variational/DiscreteNoether.lean` separately prove scalar real-path results:
a local-action minimum under every real replacement at a supplied interior
junction satisfies the discrete Euler--Lagrange equation, and differentiable
invariance transports segment momentum. If those premises hold at every
interior record, `noether_current_constant_on_finite_chain` constructs one
scalar shared by every segment momentum. The free affine path supplies a
nonzero current witness.

`Variational/FiniteHistoryBridge.lean` proves the key interface obstruction.
For a fixed real path and site, replacement by `x : Real` is injective in
`x`; therefore no finite real-path family contains every single-site
variation. The finite-state Gibbs theorem and the universal-real-variation
Euler--Lagrange theorem cannot be silently identified.

## Open B7 obligations

The package does not:

- construct the exponential history law, constraints, or multipliers from an
  OPH Axiom-3 source object;
- prove multiplier or minimizer uniqueness;
- provide an enrichment, density, or transfer theorem from the finite state
  history type to the real differential path space;
- recover stationary saddle histories from a global Gibbs minimum;
- construct physical action, time, fields, complex amplitudes, interference,
  a continuum path space, or a laboratory current.

Issue #683 owns those open obligations. The attained results are
standard conditional helpers plus a machine-checked interface no-go. They
emit no prediction-ladder row and do not enter the zero-sum flagship.
