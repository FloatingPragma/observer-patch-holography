# All 32 root pasts in one executed log

One 80-cell Interval Echo program writes five independent root bits, then leaves 32 distinct parity results on the same board. Both material receiving heads recover each result and its six-bit writer ID. In each of the 32 enrolled input runs the shelf for mask 5, for example, has precisely roots 0 and 2 in its generated past. Every one of the other 31 subsets also occurs in that same run. The records are deterministic simulated pointer readings. [Frozen protocol, commit 385506d0b](PROTOCOL.md) · [Complete receipt](records.json.gz)

## The events and their measured support

The [mounted native apparatus](../../../seam_word_meter/finite_program/layout/README.md) provides 474 balanced data dimers. The program uses 259 on 37 disjoint seven-bit shelves. Five unconditional appointments prepare one independent input bit each. Every appointment writes its unique positive identifier into the other six bit positions of its shelf. Each following mask command reads only its selected root words. Heads at ports 63 and 85 independently compare each extremal data rail against their own zero and double references. They require equal valid words with the writer enrolled at that root. A supplied parity controller combines the observed bits and commands the target's native balanced-dimer stores. The controller supplies the appointments, parity and addressing. Each store is monitored for changes outside its target, and each input survey restores the full physical state. Both heads read each sealed word again and repeat the full census after the last sink. [Instrument](apparatus.py)

For a fixed sink at an enrolled input, we hold a candidate subset of root values fixed and vary the remaining values over all admissible completions. The promise that all 37 shelves start blank and all stamps were correctly enrolled is part of this input family. The sink's public result is its own stamped word. The earlier sinks and the diagnostic tape are available to auditors but are outside that event's output menu. For all 1,024 input/sink pairs, precisely the supersets of the advertised mask force the answer. Flipping any one selected input alone reverses the parity. The unchanged target stamp leaves this distinction visible on the published low bit. [All-completion decoder](measure.py)

The producer of each parent must have written the certified register. The child must consume that parent's writer ID and exact deposited value on that register. Applying these conditions to the two physical readbacks gives 80 parent incidences per launch and 2,560 across the nominal bank. No sink reads another sink. The five immutable roots have height zero. The constant, empty-mask sink is a sixth height-zero event, and the other 31 sinks have height one. All these events coexist before the terminal read. Physical routing, controller chronology, and diagnostic rereads are available in the implementation; only authenticated semantic reads generate the stated order. [Authenticated parent definition](../../../../../../Lean/ObserverPatchHolography/Provenance/SemanticEventProvenance.lean)

## Why this one order has no faithful rank-spatial placement

A faithful placement at canonical source-height time would simultaneously place its five height-zero roots at five points in the three-dimensional spatial quotient. Each sink's accepted roots would be the five points inside a closed ball, with its own center and squared radius. A root exactly on the boundary is inside. Five points in three dimensions have a nonzero affine dependence. If closed balls realized the positive and negative coefficient groups, the difference of their two ball scores would be an affine function with strictly positive weighted sum. Affine dependence sets that sum to zero. This excludes any such placement of a single nominal log, even if all spatial positions may be chosen after the full history or permitted to vary between launches.

The checked theorem `noFaithfulSingleLog` takes the existing `SemanticEventLog`, distinct root and sink events, root height zero and the generated-past incidence. It excludes the repository's `FaithfulRankSpatialCausalPlacement`. The proof allows each sink to have its own nonnegative query time. Under source-height time the empty-mask sink has radius zero. The two independent receivers and full root intervention table establish the incidence for this program. The conclusion concerns its registered generated semantic order, metric, and source-height clock. Assignments with unequal root times or other metrics, and the order of passive host signal propagation, are separate questions. [Closed-ball and log theorems](../FiveRoots.lean)

## Controls and replay

| Intervention | Physical result |
| --- | --- |
| Sink 5 reads roots 0 and 1 while declaring 0 and 2 | The measured least support is mask 3 on all 32 inputs. Its bit differs from nominal on 16. The nominal writer-support check rejects the altered command. |
| Park the second head for root 0 at sink 31, input 21 | The first head reads `9966666`, the second reads `0000000`, and the full-mask shelf stays blank. |
| Temporarily erase the low writer bit of root 0 at that same appointment | Both heads read `9666666`. The value is still one, but the stamp is now zero. The full-mask shelf stays blank. |

Each hex digit encodes both two-reference rail responses. The least-significant dimer is first. Zero and one produce `6` and `9` respectively. Blank material reads as seven `6` digits. The shutter and stamp controls preserve the other 36 final words. The tests reconstruct every changed data port, all pairwise pivotal flips, sources and final versions from the two complete tapes. [Independent tests](test_parity.py)

From the repository root:

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
nice -n 15 .venv/bin/python -B -O -W error -m model.experiments.interval_echo_consensus.versioned_dispatch.cone_order.parity_log.measure
nice -n 15 .venv/bin/python -B -m pytest -q model/experiments/interval_echo_consensus/versioned_dispatch/cone_order/parity_log
(cd Lean && nice -n 15 lake env lean -j1 -DwarningAsError=true ../model/experiments/interval_echo_consensus/versioned_dispatch/cone_order/FiveRoots.lean)
~~~
