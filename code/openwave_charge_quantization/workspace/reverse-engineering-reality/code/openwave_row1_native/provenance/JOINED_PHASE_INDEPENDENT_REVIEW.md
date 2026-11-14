# Independent review of the joined native readout

The frozen `../auditor_bridge/joined_readout.py` passes this independent
source, probability and controller review. The implementation advances the
actual source, updates its counter only from native current records, and
then executes the literal native pointer. It receives no precomputed charge
or spatial degree. The primary source idles while the pointer is scheduled,
so the reported charge belongs to the explicitly retained source prefix.

`review_joined_readout.py` reconstructs all 456 source histories with a
separate floor/ceiling endpoint implementation. All 564 no-move events and
324 moving events remain in their original histories. Every history has
its original exact fair-coin probability, and each of the 132 prescribed
source schedules has total probability one.

The reviewer independently splits pointer probability into conclusive,
discovery-timeout, routing-timeout and measurement-timeout outcomes. It
enumerates the actual one-load moving/holding possibilities, uses a separate
breadth-first graph-distance calculation, and forms exact rational stage
probabilities. The resulting conclusive probabilities match the producer's
reset and readout factors. All 19,008 configured source/region/pointer-input
probability budgets sum to one. No success normalization is used.

These are laws conditional on configured inputs. The contract now states
that there is no implicit random prior over preparation, schedule, region
or unknown pointer position. If a future protocol randomizes those settings,
its own declared selection weights multiply the existing conditional laws.
All native coin and meter failure outcomes are already unconditional within
each configured protocol.

As an adversarial controller check, the reviewer removes the terminal-load
and initial-background audit fields from `source_prefix` before the actual
joined policy receives its result. Positive, negative, held-zero, empty-region
and full-region cases still execute correctly. The source's retained
current-derived counter supplies the meter's input. The independent review
also checks the fine-history factor `2^(-executed native attempts)` for every
complete execution example, including discovery, routing and measurement
failure. A timeout has no decoded charge; it is not a zero-charge result.

Two reporting fixes were applied by the owner before freezing: the
configured-input probability scope is explicit, and every declared input
dependency is checked before and after execution. This review pins the
corrected code and receipt and independently verifies that they remain
unchanged during the opposing review.

The method closes the full finite readout path for the stated native source
and known-reference charge class. Its spatial interpretation is supplied by
the independently checked integral dual-cell theorem, not by the pointer's
ability to display an integer. No reserved acceptance configuration was run,
and no official platform status is awarded here.

The executable and receipt are `review_joined_readout.py` and
`JOINED_READOUT_INDEPENDENT_REVIEW.json`. The earlier coherent phase-meter
artifacts remain frozen and are not needed for this literal-native path.
