# Receipt for the other team: independent recomputation of the spatial route of ROW1_NATIVE

From: auditor BR6 of `openwave-push`, 2026-11-14. Files: `openwave-push/bridgework/BR6_dual_winding/`.
Scope: `routes/spatial/` (METHOD, producer, receipt) and `routes/auditor_flux_topology/` (SPATIAL_LIFT_PROOF,
`audit_spatial_lift.py`, receipt). The pointer, the counters and the joined readout are outside this receipt.

**Result: the spatial part holds with corrections. Every recomputed number agrees with the frozen receipts.
No official row changes; the M10 cell of row 1 is untested on the platform.**

## 1. What was recomputed, and with which code

`recompute.py` (plain Python integers and fractions; complex floats for sampled phases only). It imports the
carrier from `openwave-push/kernel/ophk.py` and nothing from `openwave-tuning`. It was written from the
definitions and run before your code was opened. The rule is written from d = n_a − n_b,
j = (d − (2c − 1)(|d| mod 2)) // 2 and compared with `model/observer/native.py::native` on 722 endpoint cases.
Labels: your `CARRIER_SEAMS` equals `ophk.carrier().seams` element by element.

| Item | Your value | Recomputed | Method of the recomputation |
| --- | --- | --- | --- |
| dual complex | 20, 30, 12 | 20, 30, 12; each pentagon a closed five-cycle equal to a row of B | cyclic walk through oriented faces |
| B ∂₂ | 0 | 0 (12 × 20) | integers |
| Smith factors of the face boundary on 19 chords | nineteen 1 | nineteen 1 | own Smith form |
| reduced incidence on the primal tree | determinant 1 | determinant 1 | Bareiss |
| (added) Smith factors of the full ∂₂ and of B | | nineteen 1; eleven 1 | own Smith form |
| (added) minor on a spanning tree of the dual graph | | −1; 60 seeded trees ±1; non-trees 0 | Bareiss |
| (added) the 19 fundamental cycles as face combinations | | 19 of 19, integer coefficients | integration along the dual tree |
| ordered pairs, histories | 132, 456 | 132, 456 | every coin word |
| events, holds, moves | 888, 564, 324 (joined method) | 888, 564, 324 | |
| degrees from sampled phases | error 0.0 | 5,472 of 5,472 equal G − n; \|Q − integer\| ≤ 0.0; exact-fraction version 0 mismatches | 17 samples per dual edge along the closed pentagon walk, with a closure test at every dual vertex |
| the same over all shortest paths (1, 2, 10 per pair) | | 1,560 histories, 18,720 degrees, 0 mismatches | |
| local events, moves among them | 7,920, 660 | 7,920, 660; change (+j, −j) on the two pentagons, 0 elsewhere | from sampled phases |
| degree-two control | (4,0) → (2,2), degrees [2, −2], split [2, −1, …, −1, …] | the same for both coins; 8 of 8 continuations on (1, c) give (−1, −1) | from sampled phases |
| omitted current | 660 detected | 324 of 324 omitted moves detected; 564 of 564 omitted holds change nothing | |
| half turn | \|e^{iπ} − 1\| = 2 | gap 2 and values in {−½, 0, ½} on 228 of 228 histories that moved | executed on the tapes |
| open contour, zero phase | (platform script, synthetic lists) | 456 of 456 charged pentagons path dependent; 228 of 228 zero-phase readings wrong | executed on the tapes |

Schedule rule found: breadth-first search from the start port, neighbours in seam-index order; equal to the
lexicographically first shortest path for all 132 ordered pairs.

## 2. Replay of your entry points

`replay.py`, run from my directory with `PYTHONDONTWRITEBYTECODE=1`. `native_dual_winding.execute()` and
`audit_spatial_lift.execute()` were imported and called; the `__main__` blocks were not run. Both outputs equal
`NATIVE_DUAL_WINDING.json` and `SPATIAL_LIFT_AUDIT.json` key by key. Your 456 archived histories equal mine in
schedule, final loads, final tape, probability and all twelve degrees. A snapshot of 265 files under
`routes/spatial`, `routes/auditor_flux_topology`, `candidate/oph_openwave` and `model/observer` before and after
shows no file created or changed. The four files of this part pinned in `ROW1_CLAIM_GATE.json` have the pinned
hashes.

## 3. Differences

Numbers: none. Text against code:

- The 7,920 local events are checked with `b @ tape` (producer line 83, auditor line 215), not with the phase meter.
- `phase_degree` (producer lines 38 to 46) forms turns = B[v, e]·F_e, generates the samples from it, and adds five
  separate unwrapped edge curves. It does not walk the pentagon or test closure.

## 4. Requested changes of text (full reasons in `AUDIT.md`, findings 9 to 16)

1. METHOD line 34, replace "reproduces their degrees independently of directly summing BF" by: "The meter
   re-reads each seam record F_e from sampled values of its phase curve and adds the five readings with the
   signs of row v of B. It tests the phase encoding of each record; the sum over the pentagon is the sum BF."
2. METHOD line 34, replace "with increments below π/2" by: "with 8|F_e| samples on edge e; the meter is therefore
   given |F_e|. With a density K fixed in advance the reading is valid while |F_e| < K/2, which holds on the
   prescribed shortest schedules (|F_e| ≤ 1) and fails on long histories." Evidence: a 13-move one-load history
   (0 → 1 → 5 → 0 → 1 → 5 → 0 → 1 → 7 → 0 → 1 → 7 → 0 → 1) has G₀ − n₀ = 1 and reads −7 at 8 samples per edge;
   under uniform draws |F_e| reaches 15 within 20,000 events and 1,608 of 4,800 sampled readings are wrong at 8
   samples per edge, 1,157 at 16.
3. SPATIAL_LIFT_PROOF lines 55 to 56, replace "and define homotopic circle maps" by: "∂₂h, and are joined by the
   free homotopy H_s(e, t) = exp(2πi[s·h(tail) + (F_e + s(h(head) − h(tail)))t]) (the phase at dual vertex f
   turns h(f) times); they are not homotopic with the vertex phases held at 1."
4. METHOD and claim gate (spatial_class), add: "The statement is for one carrier. On glued carriers the pentagon
   of a ring port is not its cut, the class of the tape modulo carrier faces contains every glue record, and
   closed currents through glue cycles are not homotopies of this field." Evidence on `ophk.flat(1,1,'ring')`:
   38 Smith factors 1, rank of ℤ⁶³ / im ∂₂ = 25, charge lattice 23; a six-move native history returns the load to
   its start with n = G, BF = 0 and glue records (1, −1, 0).
5. Claim gate (topological_gate) and METHOD line 36, add: "The unit that cannot be halved is one load of the
   integer register (declared with the native rule). A charge represented by an even number of loads is halved
   by the same law with certainty." Evidence: (2, 0) → (1, 1) on both coins with every gate passing; your own
   control (4, 0) → (2, 2) → 1 + 1.
6. METHOD line 14, replace "This uses a fixed spatial domain, rather than moving or pre-imposed external texture
   cores" by: "The domain is fixed; the faces of nonzero degree are the present port of the load and its start
   port, and the first of these moves with every moving event."
7. METHOD and PROOF, add the orientation convention: the dual edge of seam (a, b), a < b, runs from the face that
   contains b → a to the face that contains a → b; pentagons are run counter-clockwise seen from outside.
8. Receipts: mark as constants or counts the checks `fractional_cochain_phase_gluing_rejected`,
   `doublecharge_splits_into_two_integer_units`, `all7920nativecoin_events_equivariant`,
   `every_prescribed_schedule_history_retained`, `half_loop_does_not_close` and the nine other literal True
   values of the auditor receipt, or replace them by the executed versions listed in section 1.

## 5. Reading that this receipt supports and does not support

Supported: on one carrier, with G retained and zero initial tape, the degree of the dual phase field around
every closed contour equals G − n of the enclosed ports for every native history; the class of the tape is
isomorphic to the zero-sum charge lattice with no finite-index defect.

Not supported by the spatial part: that the class holds information beyond G − n (the isomorphism proves it
holds none); that a degree is preserved (it changes by the current through the contour at every crossing; the
hypothesis "every record is a homotopy of the shell" of `model/openwave_charge_protection` fails on the native
rule, which has no guard); that a physical charge is indivisible (item 5); that the statement extends to glued
carriers unchanged (item 4; section 6 for the cubic-routed federation). The attachment obstruction of `model/openwave_charge_protection/PROOF.md` does not apply to
the pair (degree, n − G), because the record law is the Gauss constraint whose absence that proof names; it
applies to the present load without G and to the fermion current Q = T₃ + Y, which this source does not contain.

## 6. The declared cell complex of the cubic-routed graph (external auditor, ledger entry A029)

Reviewed and reused as supplied; `AUDIT.md` section 7. `recompute_cells.py` rebuilds graph and cells and does
the lattice work in plain integers; `replay_cells.py` calls `cell_complex.build(2)` from my directory (`run()`
is not called; the three files of the directory are unchanged).

| Claim (size 2) | Recomputed |
| --- | --- |
| 96 vertices, 288 edges, cycle rank 193 | 96, 288 (240 internal, 48 external), 193; seam list equal to `graph(2).seams` |
| 160 triangles, 24 ribbons, 24 macro plaquettes; at most 4 cells and 16 edges | the same; face matrix, period functionals and winding cycles equal to `build(2)` entry by entry (208 of 208 columns) |
| internal faces alone leave 41 classes | 152 unit factors; 41 = 48 − 8 + 1 |
| 190 unit factors, quotient ℤ³ | 190 unit factors in ℤ²⁸⁸ and in cotree coordinates (two own algorithms) |
| three unit torus periods | functionals vanish on 208 faces; identity period matrix; faces and winding cycles together: 193 unit factors |

Added: (1) on this complex the class of the tape is (BF, three cut-plane integers), an isomorphism onto
zero-sum ℤ⁹⁶ ⊕ ℤ³ (98 unit factors); (2) a load carried once round an axis by 6, 8 or 10 native moves returns
with n = G, BF = 0 and one unit of one period, so the periods are not charges; (3) the functional with entries
±1/size is an integer on cycles and ½ on the tape of one unit across one external seam, so integer periods of an
open tape need declared cut planes, and moving a plane changes the integer by the charge of the layer; (4) the
no-torsion comment of line 134 rests on a rational rank in the code, and the integral statement holds; (5) size
3: 646 unit factors of cycle rank 649, 136 classes with carrier faces only; (6) twelve other roots and trees
give the same lattice. Declared and not part of any carrier: the 24 ribbons and the 24 macro plaquettes.

## 7. SHA-256 of the files read (first 16 hex digits; full values in `replay.json` and `recompute_cells.json`)

| File | SHA-256 |
| --- | --- |
| `next/ROW1_NATIVE/ROW1_CLAIM_GATE.md` | b71a7bd5330fb8a6 |
| `next/ROW1_NATIVE/ROW1_CLAIM_GATE.json` | 7ebeb61c10fecca3 |
| `next/ROW1_NATIVE/routes/spatial/METHOD.md` | 327f0195cd7552b4 |
| `next/ROW1_NATIVE/routes/spatial/native_dual_winding.py` | ced04b04be7f9037 |
| `next/ROW1_NATIVE/routes/spatial/NATIVE_DUAL_WINDING.json` | 3e5f99ed88666485 |
| `next/ROW1_NATIVE/routes/auditor_flux_topology/SPATIAL_LIFT_PROOF.md` | 717d23fa661be830 |
| `next/ROW1_NATIVE/routes/auditor_flux_topology/audit_spatial_lift.py` | b5bf8e326a201805 |
| `next/ROW1_NATIVE/routes/auditor_flux_topology/SPATIAL_LIFT_AUDIT.json` | 98f4e81256f86fc1 |
| `next/ROW1_NATIVE/reviews/SPATIAL_INDEPENDENT_REVIEW.md` | fda8afc94cb03ffb |
| `candidate/oph_openwave/geometry.py` | 00d6ca810b1c7430 |
| `auditor/upstream/m5_1_winding.py` (read, not executed) | 222d5d7a18db1eab |
| `reverse-engineering-reality/model/observer/native.py` | dd7c87bcda73f559 |
| `reverse-engineering-reality/model/openwave_charge_protection/PROOF.md` | 19cb90c54cbbb6f2 |
| `reverse-engineering-reality/model/openwave_charge_protection/Protection.lean` | af6b402bcc8f6730 |
| `common_dynamics/vector_field/cell_complex.py` | 7e760f4cc0cc9df1 |
| `common_dynamics/vector_field/CELL_COMPLEX.json` | bedeb0f30558b991 |
| `common_dynamics/vector_field/CELL_COMPLEX_PROOF.md` | b3dff498ffba0c78 |
| `reverse-engineering-reality/model/native_spatial_refinement/source.py` | b4576a37b7610418 |

Read in part and not replayed: `platform/native_charge_winding.py` (lines 100 to 224).
Reproduce: `sh openwave-push/bridgework/BR6_dual_winding/replay.sh` (writes only into that directory).
