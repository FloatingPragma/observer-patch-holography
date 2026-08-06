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

`Variational/FiniteRealTransfer.lean` proves the composition that survives
the obstruction. Under an undercut receipt, the statement that every real
single-site variation of an embedded path is undercut in action by some
embedded finite path, finite minimality transfers to real local minimality
and to the real Euler--Lagrange packet. Membership receipts are impossible;
the committed quadratic witness carries an undercut receipt, and a one-path
control family provably carries none, with the transferred conclusion
failing there.

`InformationProjection/LogTransitionAction.lean` derives the action. For a
strictly positive row-stochastic kernel and a normalized initial law, the
Markov path law equals the exponential tilt of the step-uniform reference by
the log-transition action at multiplier one, and an action-multiplier pair
reproduces the law exactly when the multiplier-weighted action equals the
log-transition action plus a constant, so the action is unique up to the
additive-constant and multiplier-rescale gauge and the normalization
convention pins the multiplier. The committed source chain instantiates both
statements with kernel-decided receipts, and the committed repair-count
action reproduces the chain law at no multiplier.

`Variational/LegendreBridge.lean` joins the two faces. The finite Legendre
transform carries involutivity, Fenchel--Young, and degenerate controls;
discrete Euler--Lagrange transport at a junction is equivalent to one
discrete Hamilton step for the quadratic and strictly convex classes;
single-site minimizers of the log-transition local action coincide with most
probable paths of the realized chain through an exact corner identity; and
the constant chain current equals the quadratic Legendre momentum, conserved
together with the energy along the free Hamilton orbit that reproduces the
committed witness path.

## Open B7 obligations

The package does not:

- derive the reference measure; the initial law times the uniform-step
  counting weight is a declared object;
- produce the realized chain from a Hamiltonian flow; the naive
  quadratic-Gibbs identification is refuted on the chain literals;
- recover stationary saddle histories from a global Gibbs minimum;
- construct physical action, time, fields, complex amplitudes, interference,
  a continuum path space, or a laboratory current.

Issue #683 owns those open obligations. The derived-action and bridge
results enter the flagship variational passage; they emit no
prediction-ladder row.
