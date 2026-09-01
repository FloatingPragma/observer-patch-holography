# A3 Consensus-Tower Interface

`Tower/ConsensusTower.lean` defines the timeless root object used by the V2
construction in the `OPH.Tower` namespace, re-exported by the `Tower`
OPHConstruction umbrella. A tower is a nonempty directed covariant refinement system
with:

- finite observer and record types at each regulator;
- observer-indexed strict record orders from the A1 ledger;
- finite private matrix algebras and commutative public star-subalgebras;
- public representatives of record labels;
- selected density states and linear generators; and
- explicit identity, composition, order, public-algebra, record-element,
  contravariant state-restriction, and generator naturality receipts.

The structure does not assert that any physical source supplies this data.
It contains no repair fixed point, worldline, causal region, Lorentzian
metric, clock, proper time, conditional expectation, semigroup, continuum
limit, or global time. In particular, a record order is not a clock and a
linear generator is not physical evolution without separate realization
receipts.

## Existing finite fiber

`ConsensusTower.constantConsensusTower` is an adaptor, not a replacement
for the event-algebra stack. Given an existing `ProjectivePartition` and an
existing certified `StateMatrix`, it constructs a constant one-regulator
tower whose public algebra is definitionally the partition's
`publicSubalgebra` and whose records are definitionally its active projector
labels. The adaptor uses the discrete record order and zero generator, so its
inhabitation proves packaging only.

## E1 wiring

E1 can define finite causal regions and local algebra assignments as a
decoration of a fixed `ConsensusTower`. Its isotony, locality or
no-signalling, overlap descent, conditional repair, and compatible regional
states remain additional fields or theorems. The A3 `public_natural` and
`generator_natural` laws provide refinement targets; they do not prove any E1
causal or channel property.

## E2 wiring

E2 should quantify over one enriched `ConsensusTower` and attach the
source-native informational-order/placement interface and the finite causal
algebra net to that same value. Source-selected spatial placement and a
separate physical-continuum certificate are still required before this data
has a geometric reading. The record-germ atlas is an optional algebraic
soldering interface; it is not an input to the source-native route. The A3
functoriality and naturality fields express the
joint finite naturality and custody obligation without constructing a second
tower. E2 requires an inhabited nonconstant example and all order, placement,
and net-attachment receipts. A locally covariant functor, causal embeddings,
time-slice property, continuum existence, and sector or field reconstruction
are endpoint gates.
