# PR 723 follow-up: implemented experiments and review

**The follow-up experiments and alternatives are implemented and verified.
The physical requirements are not all solved; the official rows remain 🚧.**
Four Codex workers used all available concurrency slots for production and
cross-review. Everything here is private repository work, with no new upstream
submission, tag, source adoption or cell award.

The [complete maintainer review](UPSTREAM_REVIEW.md) and
[criticism-by-criticism checklist](REVIEW_CHECKLIST.md) cover the questions as
well as the three row decisions. The row packets own the executable work:

| Row | Implemented and observed | Remaining scientific requirement |
| --- | --- | --- |
| [1: charge](../../work/2026-11-15_codex_pr723_unit_selection/README.md) | Expanded multi-load source; 1,820 states and 109,200 branches; all census states have a path into the invariant unit-occupancy class. All 48 native histories settle, while 48 swap/hold controls retain their doubles. Current-only reconstruction and a general signed-integer convergence proof are checked. | The integer alphabet remains supplied. The general law selects adjacent levels; absolute 0/1 needs the density promise. Electric coupling and a non-bookkeeping topological identification remain open. |
| [2: rest energy](../../../row_02_electron_rest_energy/work/2026-11-15_codex_pr723_rest_selection/README.md) | 49,445 native calls across record-discovered localization, persistent recovery, work and pair-composition assays, and a second gate policy. Original compact state leaks by 1/315 after one unread event. The new modular gap changes from log2 to log4 with the gate policy. | Localization uses declared address feedback. Stationary transfer energy is zero; the modular gap is tunable and is not an independently predicted electron energy. A physical energy/work and clock attachment is still missing. |
| [7: double cover](../../../row_07_double_cover/work/2026-11-15_codex_pr723_rotation/README.md) | 100,908 native calls; symmetry-preserving controls, full-record antipodal tests and coherent transfer-loop alternatives. A polar program gives visibility -0.997122766 after one turn and +0.991866005 after two. Exact signs and complete rejection outcomes are checked. | Preparation, coherent eraser/analyzer and polar coupling are supplied. Parameter rotation has not been realized as physical spatial rotation. Success probabilities are approximately 10^-85 and 10^-170; tiny constant-mode contamination destroys the negative fringe. |

Failed alternatives stay in their original packets: row 2's leakage, lower
modes and pair excess; row 7's ordinary permutations, sign/axis relaxation,
failed fringe thresholds and controls. These were used to guide the new tests,
not removed after better results appeared.

## Direct answers to the questions

1. **What would count against the claims?** Each row has an explicit
   `RESPONSE_TO_REVIEW.md`. Examples are a closed unbalanced charge class or
   incorrect record reconstruction; a wrong localization/work distribution;
   and absent interference, failed two-turn return or incomplete normalization.
   Several earlier claims fail these tests. Data-corruption checks are separate
   from scientific controls, and a finite timeout is not a proof against
   almost-sure convergence.
2. **Which results require coherent memory?** The new charge selection and
   classical energy-controller/work results do not. The new interference does,
   as did the older eigenray certificate. Committed classical records cannot
   silently be reused as coherent flags; the necessary apparatus remains
   declared and uncompiled.
3. **Is the 2.45 keV excess intended?** It is a reproducible consequence of the
   original logarithmic energy under its exclusive event clock: the separated
   pair is 1024.45084870234 keV at the 511 keV anchor. No interaction supports a
   binding-energy interpretation. The new copy-clock tests distinguish
   composition conventions without changing that original result or claiming
   a new independent mass prediction.

## Review findings and fixes

The [charge cross-review](../2026-11-15_codex_pr723_cross_review/REVIEW.md)
challenged signed inputs, large background shifts, backwards observation and
timeouts. The [energy adversarial review](../../../row_02_electron_rest_energy/reviews/2026-11-15_codex_pr723_adversarial/REVIEW.md)
required actual persistent trajectories and caught an unchecked first paired
boundary; the amended reader rejects that mutation. The
[rotation adversarial review](../../../row_07_double_cover/reviews/2026-11-15_codex_pr723_adversarial/REVIEW.md)
caught missing full rejection accounting for the exact two-turn program;
the [separate supplement](../../../row_07_double_cover/work/2026-11-15_codex_pr723_two_turn_rejections/README.md)
now retains all 192 × 384 outcomes without modifying the original packet.

The [common instrument v2](COMMON_SOURCE_AMENDMENT.md) implements the bounded
positive source primitive and all conditional failure states on the original
cell graph. Its [independent complex-state review](../2026-11-15_codex_pr723_common_instrument/REVIEW.md)
checks every Kraus branch and inverse erasure. It resolves the v1 graph and
failure-state gaps; v1 remains available. This is primitive compatibility,
not a complete common native observer controller or a combined physical claim.

All these reviews are within the Codex group. They do not pass external G6 or
owner/platform G9. Original Lean builds, stopped-response proofs and producer
P1–P6 have not been retrospectively certified.

## Reproduce and inspect

For a payload that can be read in full first, each row also has one standalone
standard-library-only script. Each embeds its small inputs, reads no files,
uses no network or subprocesses and imports no author package:

| Compact script | Lines including comments | Checks |
| --- | ---: | --- |
| [Row 1](../../work/2026-11-15_codex_pr723_compact_review/compact_review.py) | 114 | Complete finite census, closure/reachability, double control and uniform-two boundary. |
| [Row 2](../../../row_02_electron_rest_energy/work/2026-11-15_codex_pr723_compact_review/compact_review.py) | 123 | Discovery, all 120 return contexts, both gate/work laws, pair composition and the old excess formula. |
| [Row 7](../../../row_07_double_cover/work/2026-11-15_codex_pr723_compact_review/compact_review.py) | 131 | Native triangle matrices, frozen 96-step input, exact one/two-turn signs, visibility, success cost and fixed control. |

Read the chosen script and run `python3 -I -S -B compact_review.py` from its
directory. Integration inspected all three complete sources and replayed them
in isolation; [row 1](COMPACT_RUN_01.json), [row 2](COMPACT_RUN_02.json) and
[row 7](COMPACT_RUN_07.json) receipts bind their exact sources and manifests.
These are verification-only arithmetic routes, not new prospective experiments,
full-record audits or automatic approval under the platform's execution policy.
The complete data archives remain substantial and separate.

For the full follow-up records and separate readers, use the command below.
From this directory, with Python 3 and already-installed NumPy:

```sh
python3 -B replay_followups.py
```

This runs nine offline commands, checks nine bound manifests before and after,
and prints a receipt without overwriting stored evidence. There is no download,
installation or execution of the archived author package. The saved
[final receipt](FINAL_RECEIPT.json) records all commands passing in about 15.6
seconds total, with no stderr. The nested Berry manifest overlaps the main
rotation manifest; their leaf counts must not be added as unique files.
Original one-shot producer scripts and `audit_row01.py` write receipts; use a
disposable copy for those. The documented replay entry points are read-only.

Local pre-run freezes, known exploratory origins, later verification amendments
and initial failures are retained. They do not imply prospective platform
registration. [Source custody](SOURCE_CUSTODY.json) preserves the exact original
329-file package at c330c6e7. The requested durable public tag remains a
[publication step](PUBLICATION.md); no remote action was taken.

The mutable row statuses, registry, tracker and generated summaries distinguish
PR 723's merged/adjudicated state from these new unfiled candidates. The final
PR head is 8f551e4b; its merge is d3029c5e. No official status is upgraded by a
successful local replay.
