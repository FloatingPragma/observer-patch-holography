# Receipt for Codex: strict tape validator (C040) and closing-sweep freeze (C041)

Auditor AU4/AU5 (Fable side), 2026-11-14. Directory `openwave-push/audits_row1/AU5_validators/`. Additive: no
frozen byte of the route changes; nothing was written under `openwave-tuning/` or in the theory repository.
Pins checked before and after each run: `multiload_readout.py` `442d6977…5d94eef`, `MULTILOAD_READOUT.json`
`b70871db…b73f8e1f`, `native_cycle_meter.py` `4db8d707…c1a05329`. Hashes of the two modules: `FREEZE.md`.

## A. `strict_tape_validator.validate(result) -> (ok, reasons)`

`result` is the dictionary returned by `multiload_readout.execute` with the inputs of the call attached
(`with_inputs(result, G, schedule, source_coins, region, pointer_initial, retries)`). Rules: R1 source tape
replays under the native rule to the stated final loads; R2 class rebuilt from the emitted currents only; R3
pointer tape replays as discovery, routing, requested traversal, each attempt a native event, cursor advancing
only after a moving current, tape ending where the policy ends; R4 winding = signed cycle steps / 6 of the
replayed tape = requested traversal, final pointer loads replayed; R5 fine weights 2^(−attempts), attempts
within (30 + 3 + 12)·R; R6 no winding after a timeout; R7 stated conditional law (R0: complete, typed, admitted).

| Input | Count | Result | Rules that reject |
| --- | --- | --- | --- |
| six complete executions of the frozen receipt | 6 | 6 accepted | |
| sixty pointer replays of the receipt, regenerated with the pinned meter | 60 | 60 valid, summaries equal | |
| fresh readouts of the pinned executable (degrees ±1, ±2; 60,324 pointer attempts) | 1,680 | 1,680 accepted | |
| fresh timeouts (280 per stage) and class-0 readouts | 840 + 720 | all accepted | |
| echo of the command, zero pointer attempts (the route returns a degree in all) | 1,680 | 1,680 rejected | R0, R3, R4 |
| the same with every pointer field fabricated | 1,680 | 1,680 rejected | R3, R4 |
| echo of command 0, zero attempts; fabricated variant | 720; 720 | all rejected | R0, R3, R4; R3, R4 |
| pointer tape with one move deleted (weights made consistent) | 16,320 | 16,320 rejected | R3, R4 |
| one non-native current in the pointer tape | 12,108 | 12,108 rejected | R3, R4 |
| one non-native current in the source tape | 5,280 | 5,280 rejected | R1, R2 (and R3, R4, R7 in 3,840) |
| timeout with the requested winding; decoded as 0; relabelled as success | 840 each | all rejected | R6; R6; R4 |
| regional class edited by one | 3,360 | 3,360 rejected | R2 |
| class edited by one with a genuine pointer tape for the edited class | 1,680 | 1,680 rejected | R2, R3, R4, R7 |
| one moving source record dropped, true final loads stated | 2,160 | 2,160 rejected | R1 (and R2, R3, R4, R7) |
| the same, whole claim recomputed from the shortened tape | 1,680 | **1,680 accepted** | none |

Dropped record: the validator alone rejects it only against the true final loads. For 1,680 of the 2,160 the
shortened tape is itself a native history from G and the recomputed claim is a valid readout of that shorter
history; for the other 480 the shortened tape does not replay. The 1,680 are rejected only with part B: the
closing sweep on the actual source detects all of them in all twelve orders (20,160 sweeps, nine mismatching
seams each).

## B. `closing_sweep.closing_sweep(G, retained_tape, source_state, order, coins)`

Returns the thirty records (seam, coin, contrast, current), `source_before` and `source_after`,
`reconstruction_before` (G − B F from the retained currents: the object compared against) and
`reconstruction_after`, the comparison per event, the mismatch positions and the verdict. The charge is read at
latch L0 (before the sweep); the sweep changes the source (one load: unchanged with probability 1/32 to 25/128);
counters advanced by the thirty sweep currents describe L1. Proof and hypotheses: `PROOF.md` (connected swept
graph, sweep records received faithfully, known G, total conserved). Exhaustive, no sampling, no restriction on
coin words (all 2^30 per sweep by a transfer recursion, checked against literal enumeration):

| Domain | Histories | Damaged tapes | Sweeps (pair × 12 orders) | Coin words | Without mismatch |
| --- | --- | --- | --- | --- | --- |
| one load, 1 to 3 events, 12 ports | 2,635,920 | 655,260 | 7,920 | 8,504,035,246,080 | 0 |
| two loads, 4 preparations, 1 to 2 events | 14,640 | 4,578 | 4,668 | 5,012,226,834,432 | 0 |

Nine mismatching seams in every sweep, independent of order and coins. No false alarm: 144 and 732 sweeps of
complete tapes, all coin words. Not detected, executed: two lost records with cancelling currents (10,560 tapes,
loads correct); sweep records replaced by the reconstruction's prediction (1,320 of 1,320 pass); a 21-seam sweep
avoiding the lost seam's ports (1,320 of 1,320 inconclusive). The sweep is an appended instrument, not a
completeness property of the current-only tape.

## C. Two distinctions asked for in C040

**Geometric degree-two fission against relative charge under two loads.** The route reads the relative charge
q_A = n(A) − G(A) against the retained preparation. From a doubly loaded port the first addressed incident event
sends (2,0) to (1,1) for both coins (60 seam events, 120 branches). In these variables it creates the pentagon
degrees +1 at the start port and −1 at the neighbour from the zero class; no class of degree two divides. Over
the 5,160 natively reachable pairs (G, n) the single-pentagon degrees are −1, 0, +1, and +2 only at the start
port of a doubly loaded preparation after both units have left (660 port instances); a pentagon of degree −2
never occurs. The start port reaches +2 by two departures with probability 1/2 under the addressed schedule
(240 schedules, 960 leaves, 480 with degree 2) and 1/90 under two uniform draws (40 of 3,600 histories). A
geometric fission, one pentagon of degree two dividing into two of degree one, needs a doubly charged arrival;
it occurs only in the spatial route's four-load control, outside the modulus-6 counter domain (2L + 1 = 9), and
was not executed here.

**One-event continuation table against a second measurement.** The continued route stores 27,360 rows = 456
prefixes × 30 seams × 2 coins: 720 distinct applications of the native event to the twelve one-load
configurations, each with conditional probability 1/60 under a supplied uniform draw; 2,280 carry a current. The
meter runs on one source history (eight tapes). The continuation appends one source event, does not advance the
counters and runs no meter: the table holds no second reading. A second measurement needs the counters advanced
by the new current, a second complete pointer run with its own timeout budget, and the joint law of both
readings; the first reading conditions the source (total variation up to 2.29 × 10⁻⁵ among conclusive readouts at
R = 16, up to 0.875 after a measurement timeout). The AU4 replay continued 2,736 joined executions 15,216 times
with zero differences: that confirms the one-event table after every meter leaf.
