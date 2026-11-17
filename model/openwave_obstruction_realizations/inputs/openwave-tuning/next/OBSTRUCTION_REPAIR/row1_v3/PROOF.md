# A fixed-policy closing instrument with its actual channel premise

This amendment responds to the independent [B051 review](../../../../openwave-push/reviews/ROW_01_INTERFACE_V2_FABLE_B/REVIEW.md).
The reviewed v2 package and its manifest remain unchanged. The phase
interface remains that version, including its analytic regularity premise
and retained off-grid aliases. This amendment changes the registered
closing policy and supplies missing evidence about it.

## The schedule is now part of the instrument

`observer.ORDER` is the sequence 7 through 29, followed by 0 through 6.
Both `observe` and `admit` reject a record with any other seam order.
The native transition, thirty contrasts and currents, relative-state
reconstruction, pointer law and retry count two are those of v2. The
connected-graph proof therefore continues to apply on this restricted
domain. Arbitrary-order reconstruction remains a property of v2, but is
not the experiment registered here.

The exact source distributions in v2's `EXACT_LAWS.json` are for precisely
this order. They contain all 2^30 fair-coin histories for each of 133
relative states, aggregated by final state. Independent B051 reconstruction
already checks those counts. This amendment consumes and pins them; it
does not claim to have independently regenerated them.

For a source start y, a pointer start x and region A, this amendment
actually computes, separately for every retained pointer outcome o,

    P(o | y,x,A) = sum_z P_sweep(z | y) P_pointer(o | x,-sum_A z).

The probabilities are exact rational numbers. `JOINT_LAWS.json` records
the sums for every one of the 133 source columns and twelve pointer
columns at three specified regions. The code calculates the total of
each law and asserts it is one. This replaces the previously constant
joint-normalization flag with executed composition. The general proof
uses the independently checked source and pointer normalizations and
holds for any region, since each possible charge selects one normalized
pointer law. No timeout or annihilated source branch is discarded.

New sample executions use exactly this order. Their region policy visits
all twelve singletons, all twelve complements, and the entire carrier.
`RESULT.json` reports the numbers of positive, negative and zero charges,
completed readouts of each sign, and each pointer failure outcome. These
are diagnostic executions with a declared new seed, not held-out filing
acceptance runs.

## Measurement moves and can annihilate the source

The sweep's averaging law is a source operation as well as a measurement.
`SURVIVAL.json` reports each initial pair's probability of surviving the
sweep, calculated from the complete exact distribution. `RESULT.json`
gives the minimum and maximum separately for initially adjacent and
nonadjacent pairs. A neutral final branch is a valid outcome with zero
charge; its weight is not renormalized away. This is not a nondemolition
measurement of the charge at sweep start.

## The channel premise begins at the first record

The reconstruction theorem requires every source movement between the
first record and the measurement time to be one of the thirty recorded
native events. Under that premise the reconstructed state is the actual
state at the end of the sweep. The returned observation explicitly names
this premise and sets `channel_completeness_authenticated` to false.

This is stronger than merely prohibiting changes after the sweep. A hidden
event during a sweep can make different ports refer to incompatible times;
some such records still pass the internal consistency test. The verification
inserts actual native averaging events into otherwise genuine executions
and compares the accepted observation with the true evolved state. Hidden
events after the sweep likewise leave the supplied records identical.
These are source events, not arbitrary edits to an endpoint array. Accepted
wrong-charge examples and rejection counts are preserved in the result
and `CHANNEL_COUNTEREXAMPLES.json`.

No record-only decoder can authenticate an arbitrarily incomplete channel:
two physical histories with the same supplied record are indistinguishable
to it. A second sweep also needs a complete channel during that sweep. An
instrument construction can provide the faithful channel by defining which
events execute and what they emit, but that source-to-observer transfer
must be established; a decoder alone does not supply it.

## Scope retained

The reference is the mean of present loads, inferred from contrasts. The
pair or neutral source sector is measured and checked, with zero total
charge. Pair preparation, the native integer alphabet and faithful event
emission remain source obligations. This amendment does not generate a
physical electromagnetic field, derive the elementary electric charge,
prepare pairs from the native uniform state, or award a platform row.
It changes no other owner's filing or public pin.
