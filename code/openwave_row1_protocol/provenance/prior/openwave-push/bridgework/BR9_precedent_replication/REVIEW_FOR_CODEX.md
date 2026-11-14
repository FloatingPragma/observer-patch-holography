# Review receipt for the native row-1 chain (auditor BR9, maintainer lens)

2026-11-14 (+07). Scope: questions 3 ("can it fail?") and 5 ("reading") of the audit brief, for all nine
clauses of `ROW1_CLAIM_GATE`, for `platform/native_charge_winding.py` and for
`platform/OPENWAVE_LOCAL_EXECUTION.json`. Verdict: **holds with corrections.** No official row changes; the
frozen statement "Official status: untested" is correct. Full findings: `AUDIT.md` beside this file.

## 1. Recomputed independently

Code: `recompute.py` in the auditor's directory. It imports nothing from `openwave-tuning`. Geometry from
`openwave-push/kernel/ophk.py` (`carrier().seams`, `carrier().faces`); the dual graph, the pentagon
orientations, the phase curves, the sampled winding, the law and its mutations are written from the
definitions. Arithmetic: `int` and `fractions.Fraction` for every identity; floating point only for sampled
phases.

| Quantity | Frozen value | Recomputed | Difference |
| --- | --- | --- | --- |
| prescribed schedules (ordered pairs of ports, one shortest path each) | 132 | 132 | none |
| complete coin histories | 456 | 456 | none |
| native events / holds / moves | 888 / 564 / 324 | 888 / 564 / 324 | none |
| BF = G − n | all histories | 456 of 456, and 43,920 of 43,920 histories of all schedules of length one and two | none |
| degree from sampled phases = signed period | error 2.2 × 10⁻¹⁶ | 456 of 456 histories (12 ports each), below 10⁻⁹ | none |
| cycle 0 → 1 → 8 → 6 → 2 → 10 → 0 | six seams | six seams in the kernel's labelling | none |
| face boundaries generate the integer cycle lattice | nineteen Smith factors 1 | all 20 face minors of the 20 × 19 chord matrix have determinant ±1 | none (different route) |
| regional charge on the one-load domain | in {−1, 0, 1} | largest 1 on 456 and on 43,920 histories | none |

## 2. Replayed from the auditor's directory

`replay.sh`: `platform/build_bundle.py` was run with a destination inside the auditor's directory, then the
copy `replay/bundle/native_charge_winding.py --output replay/portable_review.json`. Interpreter flags `-B` and
`PYTHONDONTWRITEBYTECODE=1`; nothing was created under `openwave-tuning/`. Log: `replay/replay.log`.

| Item | Result |
| --- | --- |
| bundle files against the manifest in `OPENWAVE_LOCAL_EXECUTION.json` | `native_kernel.py` b84e5ecd…, `source_packet.json` 07ecff52…, `native_charge_winding.py` 0aa34967…, `LICENSE` 7c5428d9…: all four equal |
| six checks | all true, equal to the frozen `scientific_checks` |
| histories / events / holds / largest phase error | 456 / 888 / 564 / 2.220446049250313e-16, equal |
| twelve pins of `ROW1_CLAIM_GATE.json` | all twelve match their files |

## 3. What the gates test (the corrections)

Mutated laws on the same 132 schedules, own code; the last two rows also through the bundle's own
`dual_boundaries`, `curve` and `winding` (copy in `replay/bundle`).

| Mutation | Continuity BF = G − n | Curves glue | Degree = period | Occupation is one unit | Regional bound |
| --- | --- | --- | --- | --- | --- |
| none (native) | 456/456 | 456/456 | 456/456 | 456/456 | 456/456 |
| half law, j = d/2 | 456/456 | 0/456 (defect 2.0 at the first event) | not evaluable | 0/456 | 0/456 |
| random fraction, j = u·d | 9120/9120 | 2857/9120 | 2857/2857 | 2857/9120 | 2857/9120 |
| random integer, j from {−1, 0, 1} | 1368/1368 | 1368/1368 | 1368/1368 | 539/1368 | 1081/1368 |
| one moving record omitted | 0/228 against the loads; nothing detectable from the tape alone (0/228) | 228/228 | 228/228 | 228/228 | 228/228 |
| 400 random integer tapes from no law | not defined | 400/400 | 400/400, also with vertex phases and deformations | not defined | not defined |
| half law, first event, through the bundle's `curve` | | "rejected: terminal gluing" | | | |
| random integer law, 300 histories, through the bundle's contour functions | | 300 accepted | 300 integral and equal to G − n | 283 of the 300 are not one unit | |

Requested changes, with replacement text (details in `AUDIT.md`, findings F2 to F6 and F9):

1. Label `native_record_spatial_degree_identity` and `nonlinear_homotopy_and_shared_vertex_gauge_invariance` as
   asserted identities. They hold for every integer tape, including tapes that no native history produces and a
   drawn pattern (+5, −5).
2. Add two gates that can fail, "at every step the load vector is one whole unit" and "every regional charge
   lies in {−1, 0, 1}", and three law-level arms on the 132 schedules: half law (gluing red), random integer law
   (occupation red while every degree gate stays green), omitted record (continuity red against the loads).
3. Replace "Omitting an actual record update breaks the current/degree identity and is detected" and "A
   deliberately omitted actual current gives a counter mismatch" by: "An omitted record changes the decoded
   class. The omission is detected by comparison with the source loads, which the controller does not read;
   the readout has no internal check for a lost record."
4. Replace the `topological_gate` clause by: "Closed circle-map degrees are integral. Under a law that
   transfers half units the continuity identity still holds and the phase curves do not close (junction defect
   2 at the first event); the gate therefore tests that the recorded currents are integers. That they are
   integers, and that one unit is never divided, is the integer event domain of the declared law."
5. Replace "carrier-derived spatial cohomology observable" by "the class of the native current record on the
   carrier's dual graph, read through a phase encoding chosen by the observer"; add that on twenty carriers
   eleven integer classes of a tape are not fixed by the charge (cycle rank 391, face rank 380).
6. For filing: move the script under `openwave/xperiments/m10_oph/research/scripts/` with a task id, remove the
   machine-local path from the receipt, let the script enumerate the 456 histories from the law instead of
   reading 1,371,080 bytes of stored histories, write a `verdict` string, and file the gates before the run.

## 4. Reading under the platform rule

The chain computes an observable on the carrier's own process, which is what the admission reply names as the
condition for a row to move ("until its observable is computed on the carrier's own process in-platform"). Three
pieces stay supplied and belong in a cell text: the integer event alphabet of the declared law (the integer is
not produced by the process from non-integer data), the phase encoding and the retained reference G (the
reading is the transfer since the tape was zeroed), and the one-carrier, one-load domain. Nothing responds to
the charge, which the claim states. On the row-1 precedents (M5 validated on a quadrature of a seeded field
with no dynamics; M7 partial with "its value is source-set: quantization not yet emergent"; M8 untested with
"not an integrated winding") a strict maintainer would treat the frozen package as a development run with a
blocking finding (gates that cannot fail), and a re-filed task with items 1 to 6 as a case for partial.
Validated is not supported: the rule asks for "measured rather than encoded".

## 5. Files read, SHA-256

| File | SHA-256 |
| --- | --- |
| `openwave-tuning/next/ROW1_NATIVE/ROW1_CLAIM_GATE.md` | b71a7bd5330fb8a6ece8a40649786bbbade81fae99eab2e602b1b7e9142c8b09 |
| `openwave-tuning/next/ROW1_NATIVE/ROW1_CLAIM_GATE.json` | 7ebeb61c10fecca35d1b60448aeecd9ce20d4611c8f4d294fdaf95a8a5c5f498 |
| `openwave-tuning/next/ROW1_NATIVE/platform/native_charge_winding.py` | 0aa349674a64f42e07b2cc0b38586319e767b5b4a662150ad5d9d24dfe903fc0 |
| `openwave-tuning/next/ROW1_NATIVE/platform/OPENWAVE_LOCAL_EXECUTION.json` | f67cc2a26f03294d42544752c27abdfd89775d0ff22db3cf97bb8014e467556a |
| `openwave-tuning/next/ROW1_NATIVE/platform/build_bundle.py` | f9302ead56948233a93446add704242f321b08fdddbb478fb24a5946ac39a45a |
| `openwave-tuning/next/ROW1_NATIVE/routes/spatial/METHOD.md` | 327f0195cd7552b4306bf0122f8ad53555e6f7f540e648b27ce31133b79fad75 |
| `openwave-tuning/next/ROW1_NATIVE/routes/auditor_bridge/NATIVE_CYCLE_METHOD.md` | 4dd7ad4f00e37637b16652601d4bc5cc33390e2836bd4eaef7db1aeb90dc7521 |
| `openwave-tuning/next/ROW1_NATIVE/routes/auditor_bridge/JOINED_METHOD.md` | c545ee324cc1681bd4dabdfc8b7c3c7acc4635886c0677b4c3b92fd7d5873351 |
| `openwave-tuning/next/ROW1_NATIVE/routes/auditor_flux_topology/SPATIAL_LIFT_PROOF.md` | 717d23fa661be830d8a8b61c149c8a2996b634d7f5e8b72afa07903d12858e91 |
| `openwave-tuning/next/ROW1_NATIVE/NEXT_TEN.md` | c2d3c3d2b357fb387c5418224fe0e9f97aa7d0a06f13a172ae60d951c48bcf12 |
| `openwave-tuning/auditor/upstream/m5_1_winding.py` (read, not executed) | 222d5d7a18db1eabdb5ba251213648e84e986e6354970082091f08b2b90421cb |
| `openwave-tuning/auditor/upstream/issue_482_comments.json` (read as data) | c1e9506cefbd2e116a75172c5608bab98653260c18fc5c05de30c08d8087d7c3 |
| `openwave-tuning/auditor/upstream/discussion_712.json` (read as data) | 8d5ee0bad9fc0bb5bf93194b1e051b44984d2467e5d256c15aa4af6449075dee |
| `reverse-engineering-reality/model/observer/native.py` (byte-identical at the platform pin and at head) | dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d |
| `openwave-push/kernel/ophk.py` | 5b7e10aa60af47fc67a90057469c2ed6fbf2706859ca41d8140c34fbf39f2be0 |

Not read: `routes/auditor_phase/`, `routes/multiload/`, `routes/continued/`, `routes/gauge_loop/`,
`routes/compiled_texture/`, the joined executor's code, the reviews under `reviews/`, the ledger entries C024 to
C029 and A016 to A024. The pointer probabilities (p_x, the timeout 0.000137320487, the 19,008 budgets, the cap
of 627 attempts) were not recomputed here.
