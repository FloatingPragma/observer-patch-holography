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
strictly positive row-stochastic kernel and a strictly positive normalized
initial law, the
Markov path law equals the exponential tilt of the step-uniform reference by
the log-transition action at multiplier one, and an action-multiplier pair
reproduces the law exactly when the multiplier-weighted action equals the
log-transition action plus a constant, so the action is unique up to the
additive-constant and multiplier-rescale gauge. The bare-action convention
chooses multiplier one and fixes it only when the path action is nonconstant;
a constant action leaves the multiplier invisible after normalization. The committed source chain instantiates both
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

## Conditional realization of the reference kernel and the multiplier boundary

`InformationProjection/SourceReferenceSelection.lean` conditionally realizes
the target-relabeling normal form.  If the trivial visible datum and counting
reference are supplied, their conditional-resampling repair kernel equals the
uniform transition kernel exactly.  This does not source-select those inputs:
Axiom 3 optimizes relative to an exact reference and does not manufacture that
reference.  The identity visible datum and a biased two-point reference are
separating controls showing that alternative supplied inputs can change the
kernel.

Two exact non-uniqueness controls prevent promotion.  Every positive constant
rescaling of the counting reference gives the same uniform kernel, so the
kernel does not determine the mass normalization.  Two distinct positive
normalized initial laws combined with that same uniform transition kernel
give distinct `stepUniformRef` path laws, so selecting the transition kernel
does not select the whole history reference.  The counting reference and
trivial datum remain declared identifications pending an A1-generated source
receipt.

The multiplier receipt is narrower.  The exact rational tilts of the committed
chain law at weights one and one half are both strictly positive and normalized,
and their mean actions differ.  The intermediate-value theorem matches some
multiplier to the declared empirical mean `197/1754`.  The matching quadratic
is strictly increasing on the positive exponential-parameter ray, so that
parameter is unique at the supplied target.  These results prove existence,
uniqueness at that target, and two-point sensitivity; they do not source-select
the constraint observable or level.

## Positive-Gibbs mode boundary

`Variational/StationarySaddleCoverage.lean` certifies a narrow negative for the
claim that positive-Gibbs modes recover every stationary history.  The
three-record zero history under the concave two-point Lagrangian satisfies the
discrete Euler--Lagrange stationarity condition with explicit derivative
witnesses, fails minimality against the unit variation, and carries strictly
smaller unnormalized Gibbs weight than that variation at every positive
multiplier.  This witness is a stationary maximum, not a saddle, and the
theorem does not assign it zero normalized mass or remove it from the support.

The counterexample therefore closes only the universal positive-Gibbs
mode/minimizer route.  It does not rule out genuine saddle recovery, complex or
signed stationary-phase weights, constrained ensembles, a source-selected real
enrichment, or a refinement limit.  Together with the reference-scale and
initial-law controls, it leaves source selection of the whole history
reference, the real enrichment and transfer receipt, and the physical action,
time, current, amplitude, field, and continuum attachments explicit.  The
finite theorems in this packet do not by themselves justify closing those
obligations.
