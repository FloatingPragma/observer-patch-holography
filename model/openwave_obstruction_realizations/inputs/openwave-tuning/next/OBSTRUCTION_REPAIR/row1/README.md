# Charge readout repairs

Two additive constructions repair specific defects. Neither changes the
official row-1 status or establishes the physical identity of the electron.
The original frozen sources and filing drivers are unchanged.

## Adaptive phase sampling

The filed diagnostic refuses any tape containing more than three turns on
one seam. A registered-kind native history reaches four turns. The new
reader chooses `n_e = 8 abs(F_e) + 8` intervals per spatial edge. For the
declared deformed field its phase derivative has magnitude at most
`2 pi abs(F_e) + 0.62 + 2 pi 0.23`. Each sampled increment is therefore less
than pi/4, so principal-phase increments recover the same winding. This
analytic bound is on the mathematical interpolant; floating-point checks
are diagnostic and retain a numerical-failure outcome.

The honest four-turn history now returns its exact degree. All 168 signed
and deformed diagnostic reads through maximum winding 128 agree with the
integer reference within 5.69e-14. Exceeding the declared sample budget gives
an explicit resource stop with no readout. A huge-integer control stops
before float conversion; fractional currents are rejected.

## Neutral preparation with a native pointer

The existing neutral-start arm removes the per-port reference but reads
terminal loads and is outside the old total-one-or-two apparatus domain.
`neutral_record.py` provides the missing finite instrument on a new named
preparation domain. It retains the old endpoint update law and native
pointer. Inputs are one true scalar neutral filling q, a recorded pair
preparation, native source and closing tapes, a region, and the pointer tape.
Source and pointer terminal snapshots are ignored by admission.

The finite invariant is explicit. Start with n_x=q, then transport one unit
across a prepared seam. The deviations consist of one +1, one -1, or none
after annihilation. A native seam event can move either unit, annihilate an
adjacent pair, or do nothing; it cannot split either unit. This exhausts
133 configurations. Translation covariance follows directly from the
integer endpoint sum: adding q to both endpoints adds q to both outputs and
leaves the current unchanged. All 23,940 state/seam/coin branches at three
fillings agree with an independently written integer update.

For every region, the accumulated current satisfies
`degree(R) = sum_R(q-n_x)`, which lies in {-1,0,1}. Six-state counters are
injective on this interval even though the total load is 12q. The observer
reconstructs solely from the preparation intervention and current/contrast
records, checks every native event and the complete closing sweep, and
commands the unchanged pointer with this bounded degree. Every possible
unknown pointer initial column is checked against the received pointer
tape. The fine-record effect is `2^(-N)` on compatible columns and zero
otherwise, where N counts all source and pointer coin attempts. A failure
returns no charge; it is never conditioned away. The preparation is an
explicit deterministic intervention, not a native spontaneous birth.

The fresh run executes 288 complete source/pointer tapes, including
routing and measurement timeouts, and obtains signed charges -1, 0, +1.
Malformed preparation, fractional currents, false closing contrasts and
wrong successful pointer commands are rejected. Source and pointer hidden
snapshot mutations leave admission unchanged. The exact complete pointer
success/failure law transfers from the existing instrument because its
law, command domain, retry bound and unknown initial columns are unchanged;
these 288 tapes are executions, not a fresh exhaustive probability sum.

This removes the per-port preparation array and terminal-load access from
the neutral-start instrument. The true scalar filling, prepared pair,
integer alphabet, faithful emitted E channel, scheduler and charge
dictionary remain declared. Whole-channel forgery is outside that emitted
record contract. The new preparation/domain must be named in any filing;
its evidence is not silently substituted into the old candidate.

## Reproduction and independent evidence

Run `verify_adaptive.py` and `verify_neutral_record.py` from this directory
or by their full paths. `ADAPTIVE_RESULT.json` and
`NEUTRAL_RECORD_RESULT.json` contain the outcomes and source pins;
`NEUTRAL_RECORD_EXAMPLES.json` retains all new complete tapes.

Separately, the frozen Fable-B held-out driver was executed only in its
recorded mode: all 171 recorded examples reproduced. Its neutral-start arm
was rerun on the already-used development-1 seed: 2,461,696 cut readings,
zero differences, all five mutations rejected as specified. Those results
are `RECORDED_REPLAY.json` and `NEUTRAL_REPLAY.json`. No held-out or reserved
acceptance seed was executed. These are second-party reruns of those two
scripts, not a full source-code audit or an official validation.
