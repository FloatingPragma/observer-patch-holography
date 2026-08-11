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

`Variational/RealizedHistoryLegendreNoGo.lean` proves a complementary inverse
boundary for that same committed two-state source chain. The canonical
bilinear real extension of the exact four-corner log-transition table is
affine in velocity and has no global momentum solver. For every real `a`,

\[
L_a(x,y)=L_0(x,y)+\frac a2 y(y-1)
\]

agrees with the source action on every realized binary history at every
length. For `a > 0`, however, it is strictly convex in velocity and has an
explicit velocity solver and Legendre transform. The cases `a = 1` and
`a = 2` are distinct Lagrangians with distinct Hamiltonians, while their
first off-alphabet midpoint values differ by exactly `1/8`. The summary
theorem `realizedHistory_legendre_nonidentifiability_receipt` packages the
corner agreement, all-history indistinguishability, strict-convexity,
solver, midpoint, and distinct-Hamiltonian controls. The finite source law
therefore does not select a real variation curvature or Hamiltonian
enrichment.

## Reference normal form

`InformationProjection/ReferenceNormalForm.lean` characterizes the declared
reference exactly within a declared independent-target-scrambling normal
form. Among row-stochastic kernels, invariance under every independent
relabeling of transition targets at fixed source, row-constant transition weight, and
constant log-transition step action are each equivalent to the uniform
kernel (`relabel_invariant_iff_uniform`, `row_const_iff_uniform`,
`constant_step_action_iff_uniform`), and any independently target-relabeling-invariant
row-stochastic Markov reference has path law equal to `stepUniformRef`
(`unique_invariant_reference`).  A biased two-state control is
row-stochastic, strictly positive, and neither invariant nor uniform, so the
invariance premise is load-bearing. A second stay-biased control is positive,
row-stochastic, nonuniform, and invariant under ordinary simultaneous
source-target relabeling, proving that the independent scrambling condition is
strictly stronger. This fixes the reference gauge by an
invariance property, exactly as a coordinate convention is fixed; it is a
representation-level normal-form theorem, not a source selection.

## Open B7 obligations

The package does not:

- source-select the reference measure; the step-uniform reference is now
  characterized as the unique independently target-scrambling-invariant Markov reference,
  but the invariance principle itself is a stated normal-form convention,
  not a source product;
- select one real continuation away from the binary history alphabet; the
  exact source law admits the distinct strictly convex `a = 1` and `a = 2`
  enrichments above;
- produce the realized chain from a Hamiltonian flow; the naive
  quadratic-Gibbs identification is refuted on the chain literals;
- recover stationary saddle histories from a global Gibbs minimum;
- construct physical action, time, fields, complex amplitudes, interference,
  a continuum path space, or a laboratory current.

Issue #683 owns those open obligations. The derived-action and bridge
results enter the flagship variational passage; they emit no
prediction-ladder row.

## Source selection of the reference and the multiplier boundary

`InformationProjection/SourceReferenceSelection.lean` upgrades the
target-relabeling normal form to a source product.  The Axiom-3
conditional-resampling repair kernel of the trivial visible datum under
the counting reference equals the uniform kernel exactly, so the
representation theorem's step-uniform reference is the repair law's own
output.  Two controls prove both inputs load-bearing: the identity
visible datum turns the repair kernel into the identity kernel, and a
biased two-point reference weights the rows, each breaking
target-relabel invariance.  The counting reference is the committed
source counting measure of the repair-word packet, taken as a declared
identification.

The multiplier is proved not packet-determined: the exact rational
tilts of the committed chain law at weights one and one half are both
strictly positive and normalized, and their mean actions differ, so the
packet's internal receipts admit every weight and only the declared
constraint level, the committed empirical mean action, separates them.
The multiplier selection gap therefore closes onto the declared level:
the level is a source literal and the intermediate-value receipt matches
a multiplier to it, while no packet-internal rule can produce the
multiplier without that level.
