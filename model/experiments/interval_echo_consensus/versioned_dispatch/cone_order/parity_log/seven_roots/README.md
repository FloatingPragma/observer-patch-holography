# Seven root records exclude a faithful flat spacetime placement

A single simulated Interval Echo log contains seven independently written root bits and 128 different sink events. Each sink publishes the parity of a different subset of these same seven roots. The complete program is repeated on all 128 possible root inputs. Across those runs, the two receiving stations recover all 16,384 output words and their writer identifiers. In every individual log the 128 generated root pasts comprise all 128 subsets. This order cannot be identical to the flat 1+3 Lorentz causal order on any assignment of spacetime points to these events, even when all event times and all positions are free. [Frozen protocol](PROTOCOL.md) · [Complete material tapes](records.json.gz) · [Arbitrary-time geometry proof](../../ArbitraryClock.lean)

## Which records changed the outputs?

The registered 320-cell layout provides 1,914 balanced data dimers. Its two remote read stations address 135 disjoint nine-bit shelves, occupying 1,215 dimers. A source appointment writes one bit and its own eight-bit identifier. A sink appointment offers a fixed mask, reads every selected root's value and version at both stations, and accepts matching versions with the enrolled writer. A supplied controller takes the parity of those readings and publishes one stamped output on its own blank shelf. The physical read restores the surveyed stock. The individual store's complete changed-port record and the two final censuses verify that earlier records persist. The observer menu for a sink is its own published word. Identifiers, initialized target shelves and instructions are fixed throughout the nominal enrollment. [Material register and two independent read heads](../apparatus.py) · [Execution and authentication](measure.py)

The all-zero run gives a zero payload at every sink, including the empty-mask sink whose nonzero stamp distinguishes publication from a blank shelf. In the seven independent single-root interventions, each flipped root changes exactly the 64 payloads whose masks include it. The complete bank checks this in both directions at every assignment. Of 114,688 one-coordinate comparisons, 57,344 change and 57,344 preserve the output. The observed value is constant across each full set of jointly varied omitted coordinates. Those two facts prove that the selected-root subset is the unique least forcing support for the sink, including at inputs where its own payload happens to be zero.

Each input used by the controller is another copy of the actual parent register's readout. Authentication requires the parent's matching published identifier, written register and payload. No sink is consumed by another sink. Each individual log therefore has 448 direct parent incidences, all running from the seven fixed roots to the 128 distinct sinks. The geometry uses 128 sinks present in one terminal census. Comparing separate launches establishes the supports of those events. [Independent checks](test_seven.py)

## Why arbitrary time coordinates still fail

Suppose the complete order could be represented by ordinary flat 1+3 future cones. Select the latest of the seven source events in that proposed chart. Each of the 64 sinks that includes this anchor must lie at or after every root in coordinate time. For these sinks, membership of another root reduces to the nonnegative sign of its squared Minkowski interval. Being outside requires a strictly negative interval even when times coincide.

The six remaining spacetime points admit nonzero weights with zero sum and zero weighted coordinate sum. The difference in squared-interval scores of any two sinks is affine in the root coordinates. Require one sink to contain the anchor and precisely the positive-weight roots, and another to contain the anchor and precisely the negative-weight roots. Their weighted score difference would then be strictly positive, contradicting affine dependence. The [seven-root theorem](../../ArbitraryClock.lean) applies to the generated order with arbitrary event times. The [five-root result](../README.md) applies already with five roots when canonical source-height time makes them simultaneous.

The candidate supplies the material meter, stores and transit. The experimenter supplies the stable layout, inputs, appointment schedule, write-once discipline, parity controller and observation menu. The tested order is the authenticated dependency relation between these operational publications. Connecting that relation to the host's physical propagation graph would require an independent signal experiment.

## Held interventions

| Change | Registered result |
| --- | --- |
| At mask 65, replace root 6 by root 5 while retaining root 0 | The single-root interventions give mask 33 as the measured support. Only inputs 32 and 64 in the eight-run basis bank change this sink's value. |
| Shutter the second head on root 0 at the full-mask sink for input 65 | The second head has nine zero pointer codes. The sink refuses and stays blank. |
| Temporarily erase the root-0 writer stamp at the same full-mask appointment | Both heads read an unchanged payload of one with stamp zero. This sink also refuses. |

Both refusals preserve all 134 other final material records. The archive stores every enrolled value, each source and target pointer tape, the affected physical ports and the final readings at both heads.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
nice -n 15 .venv/bin/python -B -O -W error -m model.experiments.interval_echo_consensus.versioned_dispatch.cone_order.parity_log.seven_roots.measure
nice -n 15 .venv/bin/python -B -m pytest -q model/experiments/interval_echo_consensus/versioned_dispatch/cone_order/parity_log/seven_roots
~~~
