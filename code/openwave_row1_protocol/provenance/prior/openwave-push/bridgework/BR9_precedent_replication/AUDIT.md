# AUDIT BR9: the native row-1 chain under the maintainers' lens (questions 3 and 5)

Auditor BR9, 2026-11-14 (+07). Object: `openwave-tuning/next/ROW1_NATIVE/ROW1_CLAIM_GATE.md` and `.json` (all nine
clauses), `platform/native_charge_winding.py`, `platform/OPENWAVE_LOCAL_EXECUTION.json`. Scope of this audit:
question 3 of the audit brief ("can it fail?") and question 5 ("reading"). Questions 1, 2 and 4 are touched only
where a number was needed here.

**No OpenWave row is validated by anything in this file. The chain's own official status, untested, is correct.**

Verdict: **holds with corrections.** Every identity, count and hash that this audit recomputed agrees with the
frozen package. The corrections concern what the gates test. Of the gates that the claim and the platform replay
print, the ones that carry the word "topological" are identities of the construction; the properties of the law
that row 1 needs (currents are whole units; one unit stays one unit) are tested by no gate of the platform
replay and fail under mutated laws only in gates that this audit added.

Code: `recompute.py` (own code from the definitions; geometry from `kernel/ophk.py`), `replay.sh` with
`replay_check.py` (the frozen bundle rebuilt and run inside `./replay`), `run.py` with the script
`packet/scripts/m10_2_record_flux.py` (the neutral-start experiment of section 3, item 5). Logs: `logs/recompute.log`,
`replay/replay.log`, `run.log`. Numbers: `recompute_results.json`, `replay/replay_check.json`, `results.json`.
Nothing was written under `openwave-tuning/` or in the theory repository; every interpreter ran with `-B` and
`PYTHONDONTWRITEBYTECODE=1`; no reserved or acceptance configuration was run.

## 1. Gate table: identity or able to fail

Laws run on the same 132 shortest prescribed schedules with every coin string (456 histories each; the random
laws repeated with 20 and 3 seeds). *native*: the literal rule. *half*: j = d/2. *random fraction*: j = u·d,
u drawn from {0, ¼, ½, ¾, 1}. *random integer*: j drawn from {−1, 0, 1} whatever d is. *omitted*: native law,
the first moving record is not written to the tape the reader holds (the 228 histories that contain a move).
*no law*: 400 random integer tapes with entries in [−4, 4].

| Gate | native | half | random fraction | random integer | omitted | no law | Kind |
| --- | --- | --- | --- | --- | --- | --- | --- |
| g1 continuity, BF = G − n | 456/456 | 456/456 | 9120/9120 | 1368/1368 | 0/228 | not defined | identity of the bookkeeping for every recorded law; fails only when the tape is not the record |
| g2 phase curves glue (every F_e an integer) | 456/456 | 0/456 | 2857/9120 | 1368/1368 | 228/228 | 400/400 | **can fail**: detects fractional currents |
| g3 degree from sampled phases = signed period | 456/456 | not evaluable | 2857/2857 | 1368/1368 | 228/228 | 400/400 | identity of the curve construction for every integer tape |
| homotopy and vertex-phase invariance | pass in the bundle replay | not evaluable | not run | not run | not run | 400/400 | identity (degree of a closed circle map) |
| g4 the occupation is one whole unit | 456/456 | 0/456 | 2857/9120 | 539/1368 | 228/228 | not defined | **can fail**: detects fractional and integer fragmentation |
| g5 regional charge in {−1, 0, 1}, counter classes in {0, 1, 5} | 456/456 | 0/456 | 2857/9120 | 1081/1368 | 228/228 | not defined | **can fail**: detects loads outside the one-unit domain |
| lattice: face boundaries generate every integer cycle | 20 minors, all ±1 | | | | | | property of the icosahedron; independent of the law |

Further numbers. All schedules of length one and two (43,920 histories, native law): g1, g2, g4, g5 without
exception, the literal `native()` and the in-place rule agree on all of them. Half law, first event: loads ½ and
½, tape entry ½, phase advance +0.5 and −0.5 turns around the two ports, junction defect |e^{iπ} − 1| = 2.0;
second event: defect 1.414. Half law over the 456 histories: nonzero periods ½ (576), ¼ (576), ⅛ (192).
Random integer law: largest regional charge 3; periods of magnitude 2 in 268 port readings. Omitted record:
detected from the tape alone in 0 of 228 histories, detected against the loads in 228 of 228.

Clause by clause (the nine clauses of `ROW1_CLAIM_GATE`):

| Clause | Can it fail? | Under what |
| --- | --- | --- |
| `source` | yes, as a transcription check | a kernel that differs from the balanced endpoint form (their 50 inputs); here the literal `native()` against an in-place rule on 43,920 histories. A mutated kernel is caught by the replay of the stored histories as a data mismatch, not by a topological gate |
| `spatial_class` | no | a theorem about the icosahedron, true for every law. It changes with the complex: eleven further classes on twenty carriers |
| `charge_dictionary` | only against the loads | BF = G − n holds for every law that writes its current (half, random fraction, random integer: all pass). It fails for a tape that is not the record (omitted move: 228/228), and the failure is visible only to a checker who has the loads |
| `topological_gate` | in one part | "degrees are integral" and "1 is not ½ + ½" are theorems about circle maps. The gluing requirement can fail: half law, 456/456 |
| `native_observer` | yes | the regional bound fails under the random integer law (287/1368); under fractional laws the counter update is undefined |
| `complete_readout` | not tested by execution here | the omitted-record control is decided against the loads (F4) |
| `probability_clock` | only through an enumeration error | counts and normalizations of a complete enumeration; 456 / 888 / 564 / 324 recomputed, the pointer budgets not |
| `independent_reviews` | not a gate | |
| `platform_development` | replayed | hashes and checks equal; a development run, labelled as one |

## 2. Findings

Severity: **wrong** (a statement is false), **incomplete** (something required is absent), **wording** (true in a
narrower sense than written), **confirmed**.

**F1. confirmed. Counts, identities and hashes.** 132 schedules, 456 histories, 888 events, 564 holds, 324
moves; BF = G − n on all of them and on the 43,920 histories of length one and two; sampled degrees equal the
signed periods (456 histories, 12 ports each; largest deviation below 10⁻⁹); the cycle 0 → 1 → 8 → 6 → 2 → 10 → 0
is six seams of the carrier in the kernel's labelling; every one of the 20 face minors of the chord matrix has
determinant ±1, so the quotient Z³⁰ / im(face boundaries) is the zero-sum lattice with no torsion. The bundle
rebuilt in `./replay/bundle` has the four SHA-256 values of the frozen manifest; its run returns the six checks
true, 456 histories, 888 events, 564 holds, largest phase error 2.2 × 10⁻¹⁶; all twelve pins of
`ROW1_CLAIM_GATE.json` match their files.

**F2. incomplete. The two checks named after the topology are identities.**
`native_record_spatial_degree_identity` and `nonlinear_homotopy_and_shared_vertex_gauge_invariance` cannot fail
for any integer tape. Evidence: 400 of 400 random integer tapes that no law produced pass both, with random
vertex phases and edge deformations; a charge pattern (+5, −5) that one load cannot reach is drawn along three
seams and read back as degrees (5, …, −5). The reason is in `curve()`: the curve on a dual edge is given
`signs[edge] * tape[edge]` full turns, and `winding()` returns the sum of those integers; the largest error of
2.2 × 10⁻¹⁶ is the rounding of that sum. The comparison that carries content is the integer one, Σ ± F_e = G − n,
and that is g1, the bookkeeping identity of record law C, which holds under the half law (456/456), the random
fraction law (9120/9120) and the random integer law (1368/1368). The platform's standard for such a line
(`CONTRIBUTING.md`, "Self-checks a script prints"; `dev_docs/PR_REVIEW_STANDARDS.md` D10) is the label
*asserted*.
Replacement for the two check names in `native_charge_winding.py` and for the receipt:
`asserted_identity_degree_equals_signed_period` and `asserted_identity_homotopy_invariance`, with the sentence
"These two lines hold for every integer tape, including tapes that no native history produces; they are
identities of the phase construction and are not tests of the law."

**F3. incomplete. No gate of the platform replay runs a mutated law; the rejecting controls act on hand-made
inputs.** The four controls are `curve(0, 1, .5, False)`, `winding([1., 0., 1.])`, `winding([1., 1j, -1.])` and a
comparison of the degrees of the tape (1 on seam 0) with those of the zero tape. None replaces the law and none
is evaluated on a history. Under the platform's rows D12 ("Was the mutation arm shown to fire on the arena its
gate actually runs on?") and D13 ("Is the check's parent running the pinned law, on a case where its predicate
can discriminate?") the arms are instrument-input validation. Executed here on the arena of the gate:
(a) the half law makes g2 fail in 456 of 456 histories (defect 2.0 at the first event); passed through the
bundle's own `curve`, the first event of the half law is "rejected: terminal gluing"; (b) the random integer
law passes g1, g2, g3 in 1368 of 1368 histories while the occupation is not one unit in 829 of them; passed
through the bundle's own `dual_boundaries`, `curve` and `winding`, 300 of 300 such histories are accepted with
integral degrees equal to G − n, and in 283 of them the occupation is not one unit.
The same holds in the routes behind the claim, read at the lines named (not executed here):
`routes/auditor_flux_topology/audit_spatial_lift.py` lines 252 and 253 are `half_endpoint=3` and
`assert half_endpoint!=0`, reported in the receipt as `half_loop_does_not_close`: an assertion on a literal,
which cannot fail; `routes/spatial/native_dual_winding.py` line 113 evaluates
`abs(np.exp(2j*np.pi*.5)-1)`, a constant equal to 2, with the coefficient ½ typed in;
`routes/auditor_bridge/native_cycle_meter.py` lines 311 to 318 pass `F(1, 2)` to the meter and record the type
refusal. In none of the three is a half unit produced by a law.
Consequence: the gate that detects a unit divided into fractions is the gluing gate; the gate that detects a
unit turned into several signed integers is the occupation gate g4 (or the regional bound g5); the degree
identity detects neither.
Replacement, to be added to the package before a task is filed: two gates, "at every step the load vector is
one whole unit" and "every regional charge lies in {−1, 0, 1}", and three law-level arms on the 132 schedules,
"half law: gluing red", "random integer law: occupation red, degree gates green", "omitted record: continuity
red against the loads".

**F4. wording. "Omitting an actual record update breaks the current/degree identity and is detected"**
(`routes/spatial/METHOD.md`) **and "A deliberately omitted actual current gives a counter mismatch"**
(`routes/auditor_bridge/JOINED_METHOD.md`). The mismatch is with the source loads. An observer who holds the
shortened tape reads an integer tape with closed curves, integral degrees and zero total: 0 of 228 omissions are
detected from the tape alone, 228 of 228 against the loads. The same document states that the controller does
not read the loads ("Simulator final loads are retained only for independent checking. The controller does not
read them."). The two statements are compatible only if "detected" means detected by the checker.
Replacement: "An omitted record changes the decoded class. The omission is detected by comparison with the
source loads, which the controller does not read; the readout has no internal check for a lost record."

**F5. wording. Clause `topological_gate`:** "Closed circle-map degrees are integral; half-turn has incompatible
endpoint gluing. Charge two can split into one plus one. A degree-one class cannot split into two closed
degree-one-half classes." As a statement about circle maps this is a theorem and holds for every integer
cochain. As a statement about the process it says that the recorded currents are integers. Under the half law
the unit does split into ½ + ½ and continuity holds (456/456); what fails is the reader's requirement that
tape entries be integers.
Replacement: "Closed circle-map degrees are integral. Under a law that transfers half units the continuity
identity still holds and the phase curves do not close (junction defect 2 at the first event); the gate
therefore tests that the recorded currents are integers. That they are integers, and that one unit is never
divided, is the integer event domain of the declared law."

**F6. wording. Clause `spatial_class` and the phrase "carrier-derived spatial cohomology observable"** (last
paragraph of the claim gate). The lattice statement is confirmed (F1). It is a property of the icosahedron as a
sphere and does not involve the law. The class is computed from the carrier's graph and the tape; the phase
field is constructed by the reader from the integers of the tape (the method's own words: "The primitive2π phase
interpolation is an explicit observer map"). In the platform's vocabulary that is an encoding
(`issues/482`, maintainer answer: "the exponents are measured rather than encoded"). On more than one carrier
the statement changes: on the level-0 federation (240 ports, 630 seams) the cycle rank is 391 and the face
boundaries have rank 380, so eleven integer classes of a tape are not fixed by the charge.
Replacement for the phrase: "the class of the native current record on the carrier's dual graph, read through a
phase encoding chosen by the observer".

**F7. confirmed, with a consequence for the reading. Clause `charge_dictionary`:** "Present arrival charge
equals q only where G(A)=0." Correct and necessary. In the 228 prescribed histories without a move every degree
is zero while one load is present, and after a move the reading is +1 at the origin and −1 at the arrival port:
a unit and its own origin. The observable is the charge transferred since the tape was zeroed (the push's
referee finding H1, item 2; ledger line J2). An absolute reading needs a neutral start. The smallest such
experiment was executed here with separate code on the unchanged law (section 3, item 5).

**F8. confirmed, and under-used. Clause `native_observer`:** "bounded one-load relative charge lies in
{-1,0,1}, injectively decoded." True for the pinned law on all 456 and all 43,920 histories. This bound is the
one statement of the chain that a fragmenting integer law violates (287 of 1368 histories fail g5; largest
regional charge 3). It is a property of the law and should be a gate with its mutation arm (F3).

**F9. incomplete, format. `platform/OPENWAVE_LOCAL_EXECUTION.json` and the bundle against the platform's
conventions** (`openwave-push/platform/PRECEDENTS.md`, section 5): the module path
`openwave.xperiments.oph_row1_native_development` is outside the column's tree `openwave/xperiments/m10_oph/`
and carries no task id (gates A7 and A10); the field `isolated_checkout` holds a machine-local path (the
no-local-path rule); `source_packet.json` is 1,371,080 bytes of stored development histories that 132 schedules
regenerate (gate A1 asks for summary-scale tracked data, and the M7 record applies a 1 MB limit to raw data;
the script should enumerate the histories from the law instead of replaying stored ones, which also turns the
stored-history check into a run of the pinned law); a failed check ends the script through an exception, and
the output JSON has a map of checks and no `verdict` string; no gates were filed before the run, so the
execution is a development run and is labelled as one, correctly, by the file itself.

**F10. confirmed. The claim's own boundary statements stand:** official status untested; integrality inherited
from the integer event domain; the phase interpolation an observer map; counters, controller, pointer and
scheduler declared apparatus; the source idle during the terminal measurement; G known and retained; one load;
"native C does not read F". None of them is weakened or strengthened here.

## 3. Reading (question 5)

1. What the chain shows that a maintainer can check: on the unchanged native rule and one carrier, the complete
   record of currents has an integer class per closed cut that equals the charge transferred across the cut, in
   every history of 132 prescribed schedules; the class is the only invariant of a tape modulo closed currents
   on one carrier (the lattice statement); a readout by the same native rule exists with every failure
   retained.
2. What the dual class adds to the occupation: cut-independence (a statement about all 4,096 regions at once
   through twelve periods), the absence of hidden torsion on the sphere, and a contour form of the reading. It
   adds no second determination of the integer: every degree gate passes on tapes that no law produced and on
   histories of a law that breaks the unit into several signed integers.
3. What it does not show: that quantization is other than the integer alphabet of the declared law; an absolute
   charge (the reading is relative to the retained G); a charged object in space (one carrier; the load has no
   localized stationary state); any response to the charge.
4. Status under the platform rule, on the precedents (`MAINTAINER_LENS.md`): filed as frozen, the review would
   stop at a blocking finding of the tier "a gate that cannot fail" and no cell would move; filed with the two
   added gates, the law-level arms, a neutral-start arm and terms fixed before the run, it is a case for
   partial with four caveats in the cell; validated is not supported by the rule that binds the column.
5. Neutral-start experiment (own code, unchanged law, development configuration only; `results.json`, block
   `R2_packet_development`): uniform filling, zero tape, one recorded pair, then native events. The flux of the
   tape through a closed cut equals the present charge inside, read from the loads, on 1,228,800 of 1,228,800
   readings on the carrier (all 4,096 cuts) and on 94,440 of 94,440 on twenty carriers; the values are +1
   around the excess alone, −1 around the deficit alone and 0 otherwise. Stacks of 2, 3 and 5 units separate
   into units in every history. Five mutations of law and record turn their gates red; the half law and a law
   with coin-driven unit currents leave the continuity gate green. Table in `TASK_DRAFT.md`, section 4.
6. The M5 observable itself, transposed (degree of the direction field of the mean record on a cube shell of
   the declared cubic routing, which is not at the platform pin): degree and flux are (+1, +1) and (−1, −1)
   around the two charges and (0, 0) around both and around neither; for sources of strength 2 and 0.37 the
   degree stays ±1 while the flux is ±2 and ±0.37; around each half of a split unit the degree is +1 with flux
   ½. That degree does not read the charge; the chain's own remark on a "sign-only degree" is confirmed. In one
   history 80 per cent of the shell carries no record and no direction (10 histories with the pair present,
   none with a direction on the whole shell).

## 4. Items for the other auditors (not decided here)

- Question 2: the omitted-record detection uses the loads (F4); whether any other input reaches the controller
  was not examined here.
- Question 4: eleven classes not fixed by the charge on twenty carriers (F6); the addressed scheduler against
  the law's uniform draw is listed in the canonical ledger of the platform as declared architecture
  ("Carrier, regulator, representation, scheduler, preparation and apparatus for each source").
