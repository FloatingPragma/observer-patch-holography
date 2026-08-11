# B4 finite-locality boundary

## Attained theorem packet

`ObserverPatchHolography/Locality/DependencyCone.lean` proves against the
repository's concrete single-site `localRepair` that:

- one move reads only the firing site's closed edge neighborhood;
- for every fixed exogenous repair word `w`, agreement on
  `ball S w.length` forces agreement of the final records on `S`; and
- changing one initial record outside that grown set cannot change the final
  readout on `S`.

`ObserverPatchHolography/Locality/NoSignalling.lean` separately proves two
generic finite algebraic identities on a supplied bipartite split:

- a row-normalized real first-factor kernel preserves the second marginal of
  any signed joint array; and
- a Kraus-complete local matrix family preserves the partial trace over the
  acted factor.

The non-normalized one-point control changes remote total mass. It shows that
row normalization cannot be dropped from the general algebraic identity; it
is not an operational signalling experiment.

All displayed declarations are admission-free and their axiom reports contain
only Lean's standard extensionality, choice, and quotient principles.

## Exact boundary

The fixed-word cone is an `n`-move, `n`-fold-neighborhood upper bound. It is
not proved minimal, influence need not reach its boundary, and no graph-radius
`nR`, distance, clock, or propagation speed is present. Both compared runs use
the same externally fixed word.

`ObserverPatchHolography/Locality/AdaptiveScheduler.lean` is a separate E2
conditional helper, not an enlargement of the B4 result. Given a supplied
adaptive scheduler `σ` and a supplied consultation region `R` satisfying
`ConsultsOnly σ R`, it proves agreement on `ball (S ∪ R) n`, one-site
no-influence outside `ball S n ∪ ball R n`, and a two-cell countermodel showing
that the `R` term cannot be dropped. Given a declared `ConeRefinement`, it also
proves cone-image inclusion and run/readback naturality. These results do not
produce `σ`, `R`, refinement data, fairness, a positive state, or a channel
from the source; nor do they have distance, clock, CPTP, spacelike, or
laboratory-channel semantics.

The no-signalling results do not construct an OPH graph-region product or
tensor factor, positivity and normalization for a source state or channel, a
bundled stochastic or CPTP laboratory operation, or spacelike separation.
Issue #692 (E1) is the B4 claim's sole live gate and owns only the finite
coverage/region-factor attachment. Source scheduler and state/channel
production remain open under #693 (E2), the operational clock to
#703 (E5), and continuum causal/time-slice and spacelike attachment to #700
(E3); those continuations are outside B4's own live gate.

This packet therefore closes only B4's finite helper deliverables. It emits no
prediction-ladder row and is not promoted as physical causality.
