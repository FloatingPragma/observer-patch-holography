# How many messages survive an unknown rotation?

Two observers receive a twelve-bit message after an unreported proper rotation. With no orientation reference they can distinguish 96 messages. An ordered pair of adjacent port markers restores all 4,096. The comparison uses the sharp pointer outcomes of the registered Interval Echo law on twenty owners. Each receiving station compares two material rails per digit against its own reference. The [frozen protocol](PROTOCOL.md) specifies all preparations, rotations and readouts. The [record archive](records.json.gz) contains the two stations’ separate pointer tapes for all 245,760 launch and rotation combinations.

| Reference shared with the payload | Rotations preserving the reference | Distinguishable messages | With exactly six bright ports |
| --- | ---: | ---: | ---: |
| None | 60 | 96 | 24 |
| One named host port | 5 | 832 | 188 |
| One unordered edge | 2 | 2,080 | 472 |
| Two ordered adjacent ports | 1 | 4,096 | 924 |

Each row counts a largest zero-error code. Its logarithm gives the corresponding capacity in nats. The counts apply when the reference is a common apparatus setting, independent of the sender’s message. The sixty orientations of the ordered edge supply a frame. Counting them again as sixty independent messages would change the sender's alphabet.

## Reconstructing the count

The three material planes carry the payload, the host and the second endpoint. All mounts are drawn from the proper rotations of the twelve-port cell. We enroll the ports as named register pairs on the federation. To exchange two endpoint rails the controller routes adjacent native exchanges forward and back, restoring the intermediate stocks. Both auditors then survey all three planes. The semantic record retains the markers specified by its observation menu. No decoder receives the mount index. The simulator packs 4,096 independent sharp launches into binary columns for replay. The [independent checks](test_frames.py) read the archived pairs at every address, reconstruct the expected permutation and replay held launches with the scalar native echo and inverse. Parking crossing exchanges changes a held tape. Parking the comparison echo yields double-dark refusals. [Apparatus](apparatus.py) · [Acquisition and decoder](measure.py)

For an input word, collect the set of outputs possible under all sixty mounts. In each of the four menus every pair of sets is identical or disjoint. One representative per distinct set is therefore achievable and optimal. A separate calculation counts fixed patterns by Hamming weight. For an unmarked instrument the identity fixes all 4,096 words, fifteen half-turns fix 64 each, twenty third-turns fix 16 each and twenty-four fifth-turns fix 16 each. Their average is 96. A single host leaves five compatible mounts. An unordered edge leaves two, and ordering its endpoints fixes the mount uniquely. When the marker planes can move independently of the payload, all sixty marker readings accompany every admissible payload reading. The 3,600 independently chosen payload and marker settings follow by combining the measured marginal supports. Even the ordered reference then recovers only 96 messages. A second public copy of an identical tape preserves every ambiguity. [Lean proofs of the frame, independent-mount and copy statements](Frame.lean)

## Can a receiver reconstruct the veto?

In the [second frozen protocol](VETO_PROTOCOL.md) a controller reads the host marker and queries the data at that material port. It exchanges one of two output dimers to publish a veto or an allowance. The two auditors independently read the flags. Exactly 122,880 launches produce each decision. The table asks what an observer with the indicated reference can infer about this decision from the original message.

| Public reference | Always vetoed | Always allowed | Both decisions possible |
| --- | ---: | ---: | ---: |
| None | 1 | 1 | 94 |
| Host | 416 | 416 | 0 |
| Unordered edge | 528 | 528 | 1,024 |
| Ordered edge | 2,048 | 2,048 | 0 |

The distinction is operational. A coherent interpretation of the veto from a public record exists exactly when all runs with that record give the same decision. With an unordered edge the receiver cannot tell which end is the host when the two endpoint bits differ. Fixing the controller at the original port after the mount changes 2,048 decisions out of 4,096. An unreadable host survey leaves both flags quiet. [Controller and decoder](veto.py) · [Veto readings](veto.json.gz)

For the six branch choices in PR-128, this experiment realizes marked and unmarked readouts, a rotation quotient with measured stabilizers, a redundant public copy and a finite logarithmic code count. The veto measures a one-slot presence rule on a uniform binary preparation. Identifying that slot with a protected reserve of trace \(P/24\), choosing cell versus edge pricing and attaching the records to a physical horizon require additional evidence. The global capacity certificate therefore retains its two unselected alternatives. [Capacity branch certificate](../../../../../code/capacity_readback/manifests/n_closure_branch_certificate.json)

Replay from the repository root:

```sh
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -O -W error -m model.experiments.interval_echo_consensus.publication_tree.mounted_capacity.measure
python -O -W error -m model.experiments.interval_echo_consensus.publication_tree.mounted_capacity.veto
python -m pytest -q model/experiments/interval_echo_consensus/publication_tree/mounted_capacity
(cd Lean && lake env lean -j1 -DwarningAsError=true ../model/experiments/interval_echo_consensus/publication_tree/mounted_capacity/Frame.lean)
```
