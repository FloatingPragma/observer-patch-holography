import ObserverPatchHolography.Locality.DependencyCone
import ObserverPatchHolography.Locality.NoSignalling
import ObserverPatchHolography.Locality.AdaptiveScheduler

/-!
# Locality umbrella root

Reusable finite locality interfaces: a schedule-explicit dependency-cone
upper bound for the concrete local repair and generic finite bipartite
no-signalling identities.  Their attachment to physical regions, clocks,
stochastic channels, or tensor factors is deliberately not part of this
umbrella.
`AdaptiveScheduler` bounds the dependency cone of adaptively
scheduled repair by the accumulated balls of the probe and the
scheduler's consultation region, proves the consultation term
indispensable by a two-cell countermodel where adaptivity influences a
probe no fixed word can reach, and carries refinement naturality with an
identity witness.
-/
