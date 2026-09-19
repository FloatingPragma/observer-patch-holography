# Objective, Deliverables and Exit audit

The acceptance authority is issue #776, including its September 10 evidence
check. The pre-execution specification is commit `b013a019`. The mapping of
every requirement to concrete evidence is in README.md. This audit does not
change the issue's law or accept a target dimension in place of execution.

## Requested second audit and corrections

The second audit re-read the live issue against the implementation, rather
than treating the first audit or green CI as proof of completeness. It found
three acceptance gaps and one missing scientific comparison:

* The q13 replay converted stored payloads with `int(...)`. A forged final
  payload increased by 0.5, stored as floats with a recomputed tape hash, passed.
  Replay now requires the declared integer dtypes and shapes before arithmetic.
* An appended read and offset, with the reported count and hash updated, also
  passed because the appended record was never visited. Replay now checks the
  complete array census, exact lengths and every offset against the metric
  graph before checking any per-event read list.
* Replacing either historical W3 or integer comparison with an empty object
  passed the receipt's comparison check: it only checked keys still present.
  The verifier now requires every declared law, budget and result field.
* The paired receipt included W3/isolated controls but omitted the record-metric
  family's existing one- and two-dimensional scientific controls. The explicit
  [post-run addendum](../../code/support_wiring/CONTROL_ADDENDUM.md) declares their
  execution without changing the frozen original experiment or claiming blind
  preregistration. Both complete tapes and all eight interval rows are now in
  the receipt, mirror, independent replay and fresh-reproduction comparison.

The first three findings were demonstrated against the old verifier, not
inferred from hypothetical failures. Regression tests reject these corruptions
after resealing their hashes, as well as missing/duplicated control families,
altered counts, false historical agreement and suppressed clipping flags.
No committed production payload or original interval count required correction.
All original W12 trace chunks, geometry, kernels, provenance rows and main q13
artifacts remain byte-identical to the first PR head.

Both added controls reproduce their archived exact counts. Their lag-four
dimensions are approximately 2.00658 and 3.04457, but **both are clipped**.
All smaller intervals remain visible; no agreement with a reference dimension
is required. The full verification below includes these controls. The targeted
suite now contains **60 passing tests** (43 support-wiring and 17 source-routing).

Closure assessment: the Objective is the executed canonical L3-L5 architecture;
Deliverables 1-2 have wiring census, confluence and kernel controls; Deliverable
3 has authenticated one-level/tower intervals with both placements; Deliverable
4 now includes the q13 readout and its scientific controls; Deliverable 5 has
independent replay, mutation tests and the byte-exact mirror. The Exit is met
by those finite measurements. W12 interior means an order interval with executed
past/future on the uncut, closed S2 support; q13 additionally reports the exact
continuum-cube containment test. These are explicitly different interior tests,
not evidence for a shared physical clock or a continuum dimension. The small
one-level samples still do not establish a stable 2+1 reading.

## Corrections and controls carried forward

* **The issue's own prior audit:** the abandoned partial provenance run used
  integer nearest agreement while canonical means were a separate calculation.
  Here all 3,870,360 provenance actions are exact canonical means. The verifier
  checks both output values against the actual old port versions. The historical
  integer result remains a separately labelled comparison only.
* **Custody and source provenance:** the previously referenced
  `oph_exact/support_wiring.py` is unavailable in the public revision. This
  package imports no such module. Capture checks four actual public Git blobs,
  reconstructs the existing L3 fixture exactly, and retains L4/L5 geometry and
  the committed refinement maps. No uncommitted local implementation is called
  publicly reproducible by naming a nearby commit.
* **Complete provenance:** all intra-carrier operations and every refinement
  copy are logged. A copy is an executed read of the parent's current version,
  not an invented causal edge. The two writes of a mean action are simultaneous;
  event-ID order is not silently promoted to causality.
* **Finite versus asymptotic confluence:** a cyclic-projection argument gives
  convergence to the component mean for the declared repeated sweeps. Actual
  four-sweep endpoints remain unequal across schedules, and are reported as
  such. The historical L6 floating W3 run is explicitly budgeted, not terminated.
* **Readback versus support geometry:** the two placements annotate the same
  authenticated order and use predetermined anchor selectors. A fixed interval's
  MM fraction cannot change under a coordinate relabelling. The receipt does
  not manufacture a metric cone to obtain a desired dimension.
* **A discovered historical presentation error:** q13 lag four already carries
  `continuum_diamond_inside_cube=false`. Its approximately 4.15 estimate was
  described by the package table as an interior-diamond dimension. That label
  is corrected to central diamond, and the paired result exposes the clipping
  flag alongside all three smaller interior intervals. Original receipt bytes
  and their hashes are unchanged.
* **Source selection and M1:** the supplied assignments, loads, clocks, event
  granularities and componentwise join remain declared. Their different measured
  dimensions do not identify wiring as the sole cause or derive M1. No dependency
  on the unmerged #777 branch is introduced.
* **Executed gates:** the independent verifier is connected to its own CI
  workflow, with a fresh full trace reproduction. The existing frozen mandatory
  suite runner is unchanged. Tests corrupt semantic evidence and then recompute
  its hashes, so passing a hash check alone cannot satisfy the Exit.

## Verification performed

The production run started with no output artifacts and executed all three
levels, all nine confluence controls, all declared response kernels, and the
q13 read log. The independent trace check authenticated **8,063,280 writes**,
verified exact conservation and quadratic descent at each matching phase,
recomputed every rank-three readback and all eight provenance interval counts.
Its q13 check reproduced every neighbour decision and the complete readback
metric identity, every versioned read and all four exact interval counts.

The independent kernel check reproduced every retained matrix using `T^n`
full-field Grams, independently of the producer's `T^(2n)` local-return
calculation. The historical L6 verifier also passed, including its actual
integer schedule replay, floating descent ledger and isolated response check.
The causal-poset verifier passed with six inventoried attachments; its five
original artifacts retain their original hashes. Adding the mirror required
updating the archive's file count, total bytes and aggregate inventory digest.

No dimension or slow-share threshold was used to choose a later horizon,
anchor, load seed, population or sample. The retained one-level intervals are
small (9 and 75 events), which limits any dimension interpretation. All count
rows, including those inconsistent with an expected dimension, remain visible.

Executed local gates:

| Command or comparison | Result |
| --- | --- |
| `python code/support_wiring/verify.py` | `SUPPORT_WIRING_VERIFIED all`: all 99 trace chunks, 135 kernel matrices, interval counts, paired receipt and mirror passed |
| `python -m pytest -q code/support_wiring/test_support_wiring.py code/source_routing/test_routing.py` | 60 passed after the requested second audit |
| Fresh `experiment.py trace`, `readouts.py provenance`, `readouts.py q13`, followed by `verify.py --part trace` | Passed; fresh trace manifest and every chunk hash also byte-identical locally |
| Historical `exact_federation_L6_canonical_20260909/verify_archive.py` | Passed, including actual integer schedule replay and independent canonical-control checks |
| `source_net_causal_poset/verify_causal_poset_archive.py` | `CAUSAL_POSET_ARCHIVE_VERIFIED 6 files` |

The dedicated CI workflow additionally regenerates all controls and kernels
and compares the complete fresh production run using `check_reproduction.py`.
Integer payloads, writer versions and membership arrays must match exactly;
only derived floating quantities allow the documented rounding tolerance.
