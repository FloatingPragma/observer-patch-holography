# AUDIT BR6: the spatial map of the native row-1 chain and its topology

Auditor BR6, 2026-11-14. Object: `openwave-tuning/next/ROW1_NATIVE/` (read only), parts
`routes/spatial/METHOD.md`, `native_dual_winding.py`, `NATIVE_DUAL_WINDING.json` and
`routes/auditor_flux_topology/SPATIAL_LIFT_PROOF.md`, `audit_spatial_lift.py`, `SPATIAL_LIFT_AUDIT.json`.
Binding brief: `openwave-push/ROW1_AUDIT_BRIEF.md`; coordinator update of 10:50 (items 1 to 6).

**No OpenWave row is validated by anything in this file. The M10 cell of row 1 is untested on the platform.**

## 0. Verdict and what was run

**The spatial part holds with corrections.** Every identity, count and lattice statement was recomputed with code
written from the definitions before the audited code was read, and agrees. No number is wrong. Seven statements
need a change of wording or an added scope sentence (findings 9 to 15); finding 16 lists the gates that are
constants, literal values or loop counts. The reading of question 5 is in section 3: on one carrier the audited
isomorphism itself proves that the dual class is the function G − n and holds nothing else.

| Step | Code | Log | Result |
| --- | --- | --- | --- |
| Independent implementation (original brief; written first) | `dual.py`, `run.py` | `run.log`, `results.json` | 30 checks pass; three runs of 10⁶ uniform native events |
| Independent recomputation of the audited items | `recompute.py` | `recompute.log`, `recompute.json` | 16 checks pass; the whole wrapper `replay.sh` takes 44 s of CPU (`logs/replay_sh.log`) |
| Replay of the two frozen entry points, from this directory | `replay.py`, `replay.sh` | `replay.log`, `replay.json` | 5 checks pass; both `execute()` outputs equal the frozen receipts key by key; 265 watched files unchanged |
| Existing theorem of the repository, replayed read only | one call of `protection.texture_certificate()` | `logs/rer_texture_certificate.log` | degrees (−1, +1) on 3 states, 344 guard proposals |
| Declared cell complex of the cubic-routed graph (section 7) | `recompute_cells.py`, `replay_cells.py` | `recompute_cells.log`, `replay_cells.log` | 10 checks pass; the audited face matrix equals the rebuilt one in 208 of 208 columns |

Order of work: `dual.py` and `run.py` ran at 10:40; sections 0 to 6 of `recompute.py` (R0 to R9) were written
from the brief and ran at 10:48; the audited files were opened after that; sections 7 to 10 (R10 to R13) were
added afterwards. States: theory repository head `19288e2b5`, `openwave-tuning` head
`866d151933`; the four files of this part that the claim gate pins have the pinned SHA-256 values (`replay.log`,
P0). Labels: `CARRIER_SEAMS` of `oph_openwave/geometry.py` equals `ophk.carrier().seams` element by element;
the reference face (0, 11, 5) is face 0 of the kernel; the cycle 0 → 1 → 8 → 6 → 2 → 10 → 0 consists of six seams.
No reserved or acceptance configuration of any party was run (`ophk.sphere(0,'ring')` was used for a rank
computation only; no process was run on it). Nothing was written outside this directory. Not audited here: the
pointer, the joined readout, the six-state counters, `routes/multiload`, `routes/continued`,
`platform/native_charge_winding.py` (read at lines 100 to 224 to locate its controls; not replayed).

Conventions (theirs and mine agree): seam e = (a, b) with a < b; B has +1 at a and −1 at b; a native event with
current j gives n′ = n − jB_e, F′ = F + j e; faces are oriented outward; the dual edge of e runs from the face that
contains b → a to the face that contains a → b (the seam direction turned by +90° seen from outside); the pentagon
of port v is run counter-clockwise seen from outside; degree_v = (BF)_v = G_v − n_v; relative charge q = n − G.

## 1. Findings

Severity: **wrong** / **incomplete** / **wording** / **confirmed**.

### Confirmed

1. **confirmed. The dual complex.** 20 dual vertices of degree 3, 30 dual edges, 12 pentagons; every pentagon is
   a closed walk of five distinct dual edges through five distinct dual vertices; its signed crossings are row v
   of B; Euler characteristic 2; B ∂₂ = 0 (12 × 20 zero matrix, exact). `recompute.log` R1, R2; independently
   `run.log` A1 (the explicit loops of 116 port sets are closed walks; 94 are single closed curves).
2. **confirmed. The lattice statement is integral, not a statement of ranks.** Smith factors of the face boundary
   ∂₂ (30 × 20): nineteen 1; of B (12 × 30): eleven 1 (own implementation, exact integers). The two objects of
   the audited certificate: reduced incidence on the breadth-first spanning tree of the port graph from port 0
   (rows 1 to 11; seams 0, 1, 2, 3, 4, 7, 8, 18, 23, 12, 16): determinant 1 (theirs: 1); face boundary restricted to the 19
   chords: nineteen Smith factors 1 (theirs: the same). Three further certificates of mine: the 19 × 19 minor on a
   spanning tree of the dual graph has determinant −1; 60 seeded dual spanning trees all give ±1 and seeded
   19-subsets that are not trees give 0; each of the 19 fundamental cycles is an explicit integer combination of
   face boundaries (integration along the dual tree, residual 0). Hence ker B = im ∂₂ as lattices. R3a to R3c, R10.
3. **confirmed. The isomorphism ℤ³⁰ / im ∂₂ → {ρ ∈ ℤ¹² : Σρ = 0}, [F] ↦ BF.** Onto: the eleven tree paths map to
   e₀ − e_x. One to one: on 300 seeded tapes, F minus the tree representative of BF is an integer face
   combination. R4.
4. **confirmed. 456 one-load histories.** 132 ordered (start, arrival) pairs; 60, 60 and 12 pairs at distance 1, 2
   and 3; every coin word of the prescribed shortest schedule: 120 + 240 + 96 = 456 histories, 888 events, 564
   holds, 324 moves, 132 arrivals. The degree at each of the 12 pentagons, computed from 17 sampled phase values
   per dual edge along the closed cyclic walk, equals G_v − n_v in all 5,472 cases; largest |Q − integer| = 0.0;
   largest gap between consecutive curve ends at a dual vertex 2.4 × 10⁻¹⁶; the same with phases kept as exact
   fractions of a turn: 0 non-integers, 0 mismatches. R5. Port by port against the frozen receipt: schedules,
   final loads, final tapes, probabilities 2^(−length) and all degrees equal in 456 of 456 histories
   (`replay.log` P3).
5. **confirmed. The schedule rule.** `path()` of `native_dual_winding.py` (lines 49 to 61) is a breadth-first
   search from the start port with neighbours in seam-index order, stopped when the arrival port is dequeued.
   For all 132 ordered pairs its schedule equals the lexicographically first shortest path (P3). The schedule is a
   function of the declared start and of the chosen arrival; no hidden state enters. There are 1, 2 and 10
   shortest paths at distance 1, 2 and 3; over all of them (1,560 histories, 18,720 pentagon degrees) the result
   and the counts of holds and moves per coin word are the same (R5).
6. **confirmed. 7,920 local events.** 132 ordered (start, present) pairs with the tape of the all-moving word,
   30 seams, 2 coins; 660 moves (theirs: 660). From sampled phases: the degree after the event equals G − n′ at all
   12 pentagons; the change is (+j, −j) at the two pentagons of the seam and 0 at the other ten; 0 exceptions. The
   720 events from a load at rest on its start port: 0 exceptions. R6. The audited scripts check these events with
   the matrix sum b @ tape (producer line 83, auditor line 215), not with the phase meter.
7. **confirmed. The degree-two control and the split.** (4, 0) → (2, 2) on seam (0, 1) for both coins, tape 2,
   sampled degrees +2 at port 0 and −2 at port 1, equal to G − n. Each of the eight continuations on a seam (1, c)
   with c empty gives (1, 1): the degree −2 at port 1 becomes −1 and −1 with certainty; all 60 continuations obey
   degree = G − n. R7. Theirs: degrees [2, −2], split_degree [2, −1, 0, 0, 0, −1, …]: equal.
8. **confirmed. Replay and completeness of inputs.** `native_dual_winding.execute()` and
   `audit_spatial_lift.execute()`, called from this directory without their `__main__` blocks, reproduce
   `NATIVE_DUAL_WINDING.json` and `SPATIAL_LIFT_AUDIT.json` with no differing key (`replay.log` P1, P2). The
   meter `phase_degree(tape, vertex, incidence)` receives the tape and the incidence only; the terminal load
   enters the comparison, not the meter. A tape that forgets one moving event is detected in 324 of 324 cases;
   a forgotten hold changes nothing in 564 of 564 cases (a hold writes no current) (R8).

### Corrections

9. **wording. "Independently of directly summing BF" (METHOD line 34; `native_dual_winding.py` lines 38 to 46).**
   For each of the five seams of port v the meter forms turns = B[v, e]·F_e, generates the samples
   exp(2πi·turns·k/count) itself, unwraps them and adds the five results. It does not walk the pentagon, does not
   use the dual vertices and does not test closure. The value is Σ_e B[v, e]·decode(encode(F_e)).
   Replacement: "The meter re-reads each seam record F_e from sampled values of its phase curve and adds the five
   readings with the signs of row v of B. It tests the phase encoding of each record; the sum over the pentagon is
   the sum BF." (The cyclic walk with a closure test is done in `recompute.py` and, by its description, in
   `reviews/review_spatial_independent.py`.)
10. **incomplete. The sample density is an input taken from the tape (line 43: `count = max(8, 8*abs(turns))`;
    their review: "eleven subdivisions per absolute turn plus eleven").** A meter with a density fixed in advance
    reads the curve of F_e correctly only while |F_e| is below half that density. Executed (R12): the 13-move
    one-load history 0 → 1 → 5 → 0 → 1 → 5 → 0 → 1 → 7 → 0 → 1 → 7 → 0 → 1 has records (5, −2, −2) on the seams
    (0,1), (0,5), (0,7); G₀ − n₀ = 1; with 8 samples per edge the reading is −7; with the density taken from the
    tape it is 1. Under the law's own uniform schedule (one load, zero tape, 20,000 events, 1,682 moves) the
    identity BF = G − n holds at every event, the largest |F_e| reaches 15, and of 4,800 sampled pentagon readings
    a fixed density of 8 misreads 1,608 (first at event 1,650) and a fixed density of 16 misreads 1,157; the
    tape-dependent density misreads none. On the frozen domain the largest |F_e| is 2, so the executed cases are
    not affected. Replacement for "with increments below π/2": "with 8|F_e| samples on edge e; the meter is
    therefore given |F_e|. With a density K fixed in advance the reading is valid while |F_e| < K/2, which holds on
    the prescribed shortest schedules (|F_e| ≤ 1) and fails on long histories."
11. **wording. "Define homotopic circle maps" (SPATIAL_LIFT_PROOF lines 55 to 58; the brief's "closed currents are
    homotopies").** The step from "differ by a dual coboundary" to "homotopic" is not written out. It holds as a
    free homotopy: for F′ = F + ∂₂h the family H_s(e, t) = exp(2πi[s·h(tail) + (F_e + s(h(head) − h(tail)))t])
    is continuous on the dual graph for every s, starts at the curve of F, ends at the curve of F′ and keeps every
    pentagon degree (R11: 40 seeded pairs, nine values of s, gap at dual vertices below 3 × 10⁻¹⁵, 0 degree
    changes). During the deformation the phase at dual vertex f is exp(2πi·s·h(f)); it leaves the value 1 (280 of
    360 sampled stages). Relative to the dual vertices the curves of F and F′ are not homotopic: F_e is the degree
    of edge e relative to its end points. Replacement: "…differ only by a dual coboundary ∂₂h and are joined by the
    free homotopy H_s (the phase at dual vertex f turns h(f) times); they are not homotopic with the vertex phases
    held at 1." Two further facts belong beside it: a closed current is a property of a tape, not a native event
    (each event with j ≠ 0 changes the class); and the statement is B ∂₂ = 0 read in the dual complex.
12. **incomplete. Several carriers.** No carrier face contains a glue seam. On `ophk.flat(1,1,'ring')` (two
    carriers, three glue seams, 24 ports, 63 seams) the face boundary has 38 Smith factors 1, so ℤ⁶³ / im ∂₂ is
    free of rank 25, while the charge lattice has rank 23 (B has 23 Smith factors 1). The two extra classes are
    the cycle rank of the graph of carriers (3 − 2 + 1). In general, for c carriers and g glue seams:
    rank = 11c + g against 12c − 1, difference g − c + 1. Each glue record is a class invariant. Executed (R9): one
    load carried by six native moves round the cycle through glue seams 0 and 1 returns to its start port with
    n = G and BF = 0, and leaves the glue records (1, −1, 0); its tape is closed and is not a face combination.
    After the first glue crossing the pentagon degree of the ring port it left is 0 while G − n = 1: the cut of a
    ring port has six seams and the pentagon has five. With one cell per plaquette added (`run.log` C7) the extra
    rank is 2 on the torus `flat(3,3)` and 0 on `sphere(0)` and `sphere(1)`; on the torus the two remaining
    classes are not charges (at 8 of 9 neutral times of a 10⁶-event run every closed-cut winding was 0 and they
    were not, C9). Sentence to add to METHOD and to the claim gate, clause spatial_class: "The statement is for one
    carrier. On glued carriers the pentagon of a ring port is not its cut, the class of the tape modulo carrier
    faces contains every glue record, and closed currents through glue cycles are not homotopies of this field."
    For the cubic-routed federation with the external auditor's declared cells see section 7.
13. **incomplete. The unit of the non-halving statement (claim gate, clause topological_gate; METHOD line 36).**
    The statement is about one load of the integer register. Their control shows the complement: degree two
    becomes 1 + 1 with certainty (finding 7). Executed here (R13): a charge prepared as two loads on one port is
    halved by the unchanged law on both coins, (2, 0) → (1, 1), with a closed curve (gap 2 × 10⁻¹⁶) and degrees
    (+1, −1): every gate of the route passes. A law that moves half units (j = d/2) is rejected: the curve of
    F_e = ½ ends at −1, gap 2.0. Sentence to add: "The unit that cannot be halved is one load of the integer
    register (declared with the native rule). A charge represented by an even number of loads is halved by the
    same law with certainty."
14. **wording. "A fixed spatial domain, rather than moving or pre-imposed external texture cores" (METHOD line
    14).** The domain is fixed because all twelve faces are punctured whether charged or not. The faces that
    cannot be filled are {v : G_v ≠ n_v}; this set changes at every moving event (324 of the 888 events). From
    the moment the load leaves, its start port carries degree +1: the defect pair is the load and its start port
    (referee H1 of this push, findings 1 and 2). Replacement: "The domain is fixed; the faces of nonzero degree
    are the present port of the load and its start port, and the first of these moves with every moving event."
15. **wording. "Prevents exactness of the cocycle" (METHOD line 24) and the orientation of the dual edge.** On the
    complex with the face present the cochain is not closed there ((δF)(v*) = (BF)_v); on the punctured domain it
    is closed and is exact only if every period vanishes. The extension claims are correct: a punctured closed
    pentagon retracts radially onto its boundary, so the punctured sphere retracts onto the dual graph (first
    Betti number 30 − 20 + 1 = 11); a degree-zero boundary map extends over the disk through its real lift; a
    boundary map of nonzero degree has no continuous extension, because a circle map that extends over the disk
    is null-homotopic. Neither text states the direction of the dual edge. The producer needs none (its meter
    uses B only); the auditor uses tail = the face with coefficient −1, head = the face with +1. Reversing the
    orientation of the sphere changes ∂₂ to −∂₂ and no degree. Sentence to add: the convention of section 0.
16. **wording. Gates that are constants or loop counts.** Producer (line 115):
    `fractional_cochain_phase_gluing_rejected` is |e^{iπ} − 1| > 1.9; `doublecharge_splits_into_two_integer_units`
    is the literal True (the test is the `raise` of line 110, computed with b @ split);
    `all7920nativecoin_events_equivariant` is the loop count 7,920 (the test is the `raise` of line 83);
    `every_prescribed_schedule_history_retained` is a count above 132; `forgotten_record_mutation_detected` is a
    count above 0 of events with nonzero current, in which the old tape is compared with the new loads (660 of
    660). Auditor (lines 256 to 266): ten of fourteen checks are the literal True behind earlier `assert`
    statements, and `half_loop_does_not_close` is `3 != 0`. In `platform/native_charge_winding.py` (lines 200 to
    209) the controls half turn, zero crossing and open loop are three synthetic inputs (`curve(0, 1, .5, False)`,
    `[1, 0, 1]`, `[1, 1j, -1]`), not mutations of the law or of a tape of the 456 histories. Executed versions on
    the histories (R8, R13): omitted move 324 of 324; half turn on the 228 histories that moved: gap 2 at a dual
    vertex and values in {−½, 0, ½} in 228 of 228; the five open four-edge contours of a charged pentagon give
    different values in 456 of 456; zero phase gives 0 in 228 of 228; half-unit law rejected.

No finding of severity **wrong**.

## 2. The five questions of the brief, for the spatial part

1. **Correct?** Yes: findings 1 to 8. Differences between the recomputation and the frozen receipts: none.
   Differences between text and code: the local events are checked by matrix sums (finding 6); the phase meter
   is an encode and decode of each record (finding 9).
2. **Complete?** Inputs used and not named: the sample density, which is |F_e| (finding 10). Named and used:
   the tape; the incidence B; the faces of the carrier; the zero initial tape; G for every charge reading; the
   prescribed schedule. The terminal load is not an input of the meter. A dropped current is detected (324 of 324);
   a dropped hold is not detectable and changes nothing.
3. **Can it fail?**

   | Gate | Status | Mutation under which it fails (executed) |
   | --- | --- | --- |
   | degree_v = G_v − n_v | identity of the update rule pair (n′ = n − jB_e, F′ = F + j e); ledger M2, C3 of this push | a tape that forgets a moving event: 324 of 324; record of the coin only: 10,309 of 10,310 samples; record per port: 10,205 of 10,310 (`run.log` D1 to D3) |
   | change (+j, −j) at the two pentagons, 0 elsewhere | the same identity, per event | a dropped seam breaks exactly the loops that cross it: 0 errors on loops that do not (`run.log` D2) |
   | degree is an integer | identity for integer tapes | the half-unit law: gap 2.0, value ½ (R13) |
   | the curve glues at dual vertices | identity for integer tapes | half turn per unit: 228 of 228 (R8) |
   | lattice isomorphism | property of the icosahedron | two glued carriers: rank 25 against 23 (R9) |
   | (4, 0) → (2, 2), then 1 + 1 | consequence of the rule on even differences | none run; a rule that rounds even differences would change it |
   | constants and counts of finding 16 | cannot fail | none exists |
4. **Hidden assumptions.** (a) G is known and retained: the degree is the transfer since the tape was zeroed; a
   load at rest has degree 0 (720 at-rest events; `run.log` D6: in the two-load class the winding differs from
   the present load at 2,000 of 2,000 sampled times). (b) Zero initial tape. (c) One load: total one is not a
   configuration of the frozen local carrier class of the repository (`model/completion/program_class_v2.json`:
   "all78 configurations; no arbitrary total"), and the degree-two control uses total four; the identity itself
   holds for every integer configuration. (d) The schedule is addressed, not drawn; the identity does not depend
   on it (uniform draws: 20,000 one-load events and 3 × 10⁶ events with pairs, 0 violations). (e) The 2π turn per
   load is the observer's normalization: finding 13. (f) The sample density: finding 10. (g) The faces of the
   carrier enter as the dual complex; the degree of a port needs only its seams. (h) One carrier: finding 12.
   (i) Sign: q = n − G = −degree; a region that excludes the start port has degree −n(R).
5. **Reading.** Section 3.

## 3. Question 5: what the dual class gives a maintainer beyond the occupation n − G

**As an invariant: nothing.** Finding 3 is the proof. The class of the tape in ℤ³⁰ / im ∂₂ and the vector
BF = G − n determine each other through an isomorphism of lattices with no finite-index defect. On one carrier
the homotopy class of the phase field is the occupation difference, written in another group. The integer is the
integer current of the native rule; the sampled curves are generated from it (finding 9).

**As a procedure a maintainer can run, five things that the pair (n, G) alone does not provide:**

1. A reading of the relative charge of a region from records on its boundary seams only, without access to the
   loads inside (126 dual cycles on the carrier at every moving event; 396 and 1,323 loops on the screens;
   `run.log` C1, C8, C10). This is the Gauss identity of the records.
2. A continuous S¹-valued map on a fixed domain whose degree around a contour is computed from sampled values,
   the category of the validated precedent. In the upstream script (`auditor/upstream/m5_1_winding.py`, read, not
   executed) the field is the model's own state on a 65³ grid, seeded with hedgehogs of sign ±1, `N_RELAX = 0`,
   the degree is taken on a sphere of radius 5 voxels at the fixed grid resolution, pass at |Q − integer| < 0.10,
   and a comment in its configuration block (lines 55 to 61) records that the engine's relaxation loses the
   winding after 50 to 100 steps. There the field could have lacked the degree, and the tracker measures it at a
   resolution that does not depend on the answer.
   Here the map is built from the integers it returns, and its sampling density is taken from them (finding 10).
3. Invariance under closed currents, as an executable free homotopy (finding 11).
4. Mutation controls that fail: omitted current, half-unit law, open contour, zero phase (findings 8, 13, 16).
5. The integral statement itself: no hidden finite-index class on one carrier, two extra classes on two glued
   carriers (findings 2, 12).

**What it does not give:** an integer that was not put in by the integer event domain; a present charge without
G; a field that the law reads (native C does not read F: no energy, force or back-reaction); a degree that
survives the next event (section 4); a unit that stays indivisible when the same charge is prepared as two loads
(finding 13); validity on glued carriers without further cells (finding 12).

**Under the platform rule.** The source is the pinned native rule and the observable is computed on its own
records, so the chain is a case for partial, with the supplied pieces in the cell: known G, zero tape, one load,
addressed schedule, the 2π normalization, the sample density. Against the precedents of row 1
(`platform/PRECEDENTS.md` lines 100 to 107 and 366): it is an integrated winding, which the M8 cell is not; its
sign is set by the prepared start port, the fact recorded for M4; its value is set by the source, the caveat of
the M7 partial ("quantization not yet emergent"); it is not the M5 case, because the map is an observer
construction and not the model's field. The present M10 wording "quantization is assigned" still describes it.
The claim gate says the same of itself (integrality inherited, observer map, G retained, official status
untested); with findings 9 to 16 applied its text is accurate.

## 4. The protection theorem of the repository on the dual phase field

Existing theorem (`model/openwave_charge_protection/PROOF.md` lines 5 to 24; source law
`model/compact_defect_transport/PROOF.md` lines 15 to 23 and 44 to 60, `source.py` lines 58 to 86): on the compact
texture the two shell degrees (−1, +1) are preserved by every record of the four-proposal law, because an accepted
one-vertex replacement obeys guards that make the normalized linear interpolation a homotopy of the whole shell,
and a rejected record is the identity. The Lean file (`Protection.lean` lines 10 to 20) holds five integer
lemmas "conditional on an integer-valued shell degree"; it formalizes neither the degree theorem nor the homotopy
(PROOF lines 81 to 84). Replayed read only: degrees (−1, +1) on the 3 continuing states, 344 guard proposals.

| Hypothesis of the theorem | On the dual phase field under native events | Status |
| --- | --- | --- |
| Texture: rational directions on a compact support, north outside | the map of the dual graph into S¹ with F_e turns on dual edge e | holds as a definition; target S¹ in place of S², so the degree is a winding of a contour on a surface, not of a shell in space |
| Shell: closed oriented triangulated sphere with nonzero, pairwise non-obtuse directions on each triangle | any closed dual contour (a closed cut of ports) | holds; no admissibility condition is needed, the curve of an integer record closes at every dual vertex |
| Degree: integer, values (−1, +1) | integer; values set by the preparation | holds; integrality is the integrality of the records |
| Every record is a homotopy of the whole shell (guards) or the identity (rejection) | a hold is the identity; a move on a seam the contour does not cross leaves the field on the contour unchanged; a move with current j on a seam the contour crosses changes the degree by ±j | **fails**: the native rule has no guard. One event takes the winding of the contour round an excess unit from 1 to 0 (`run.log` F1; the 324 moves of the 456 histories; in the 10⁶-event carrier run 39,867 moves were on a seam of the contour round the port on which the latest pair had placed its excess) |
| Cores and exterior fixed; no core is created | the start-port defect is fixed by G; the load moves; no unit is created (excess and deficit counts never increase: 578 endpoint cases at 17 thresholds; reachable sets of 133, 3,103 and 3,214 configurations; `run.log` B1, B2) | half holds |
| Conclusion: the degrees are the same on every finite history | false for every contour that separates two ports | **not instantiated** |

With "preserved" read as "changed only by the current through the shell" every hypothesis holds and the
conclusion is: for a closed dual contour with integer port weights r, the degree after an event on seam (a, b)
with current j is the degree before plus j(r_a − r_b). Proved (linearity of r·BF), executed at 188,965, 12,307
and 56,828 moving events on a carrier, a torus screen and a sphere screen and on the 7,920 local events, 0
exceptions. This is the continuity equation of the records (the identity N_initial(S) − N_final(S) = J_out(S) of
`model/openwave_native_defect_records/PROOF.md`; "Adding a closed integer link chain changes no boundary charge"
of `model/native_criterion_packet/PROOF.md`), stated for contours. It is a balance law, not a protection of a
degree: the texture theorem excludes the change, the native statement accounts for it.

The discrete counterpart of the homotopy argument has three parts, each proved above: a contour moved across a
port of degree 0 keeps its degree; a record off the contour leaves the field on the contour unchanged; a face
boundary added to the tape is a free homotopy. What transfers from the Lean file without change: the two
no-equal-halves lemmas (with q the degree of a daughter contour), and `integer_curl_step_preserves_sum`
((a + j) + (b − j) = a + b), which is the event step on the two pentagons of a seam.
`integer_winding_permits_unequal_signed_split` (2 + (−1) = 1) is arithmetic that the native rule excludes by a
lemma the file does not contain: the two loads after an event lie between the two loads before it, so the sum of
any convex function of the loads does not increase and no unit of either sign is created. That lemma and the
contour form of the event step are new statements; both are integer linear arithmetic.

**The attachment obstruction** (PROOF lines 112 to 141: on the unrestricted texture carrier the same winding
carries charge 0 or 1, so no function gives Q = f(W); "The compact construction supplies no local Gauss
constraint tying its classical degree to those fermion occupations"). On the native source with zero initial
tape the record law is that constraint: n + BF = G on every reachable state, so the relative load n − G equals
minus the degree. For this pair the obstruction does not apply, and it does not apply for the reason the PROOF
names. It applies in two places. First, to the present load: the tapes of a load at rest on port 0 and of a load
at rest on port 3 are both zero, equal degrees and unequal n; the dictionary needs G. Second, to every charge that
is not the load number: the obstruction of the PROOF concerns the current Q = T₃ + Y of the fermion module, which
the native source does not contain; the equality "relative load = −degree" neither uses nor supplies the
hypercharge embedding, the normalization of the unit or the identification
`port_bracket_physical_current`.

## 5. The independent implementation (original brief), in numbers

Preparation in which the winding is sourced by present loads: uniform filling 3, zero records, excess-deficit
pairs created across seeded seams with their record (declared; the native rule creates no pair), then uniform
seam draws and fair coins. Winding w(R) = −(BF)(R) = present excess of R.

| Run | Events | Moves | Loops | Identity and locality violations |
| --- | --- | --- | --- | --- |
| carrier, 20,485 pairs in 11,706 preparations | 10⁶ | 188,965 | 126 | 0, 0; port windings −3 to 2 |
| torus screen `flat(3,3,'ring')`, 216 ports, 120 pairs | 10⁶ | 12,307 | 396 | 0, 0 |
| sphere screen `sphere(1,'ring')`, 960 ports, 2,080 pairs | 10⁶ | 56,828 | 1,323 | 0, 0; region windings −10 to 9 |

At 11,705 neutral times of the carrier run the tape was a nonzero closed current at 11,704 and every winding was
0. Registers modulo 7 decode every region; modulo 6 misreads 157 regional values (three units of one sign).
An angle field on the 20 dual vertices in Coulomb gauge, read with the principal-value rule, has vortex number
equal to the excess at 9,965 of 9,965 samples where every seam flow of the Coulomb field is below ½ (one pair:
at most 11/30) and differs at 338 of the 345 samples where it is not. Unchanged rule from one pair: 133 reachable
configurations, never two ports of winding +1. Rule without rounding: windings ½ and ½. Doubled register: ½ and ½
with certainty.

## 6. Statements used above, with their proofs

- **Update identity.** One event adds j e to F and −jB_e to n, so n + BF is the same vector after every history;
  with a zero tape it is G.
- **Contours.** A closed dual contour is Bᵀr for integer port weights r (on the sphere every dual cycle is of
  this form: ker ∂₂ᵀ has rank 11 = rank Bᵀ, and integration along a spanning tree of the port graph gives integer
  r). Its degree is ⟨F, Bᵀr⟩ = r·BF. An event on seam (a, b) with current j changes it by j(r_a − r_b). For the
  0/1 weights of a region the seams inside cancel and the degree is the sum over the enclosed pentagons.
- **Integral class.** B ∂₂ = 0; the reduced incidence on a spanning tree is unimodular, so B maps onto the
  zero-sum lattice; the face boundary restricted to the chords is onto ℤ¹⁹ and a cycle is fixed by its chord
  coefficients, so ker B = im ∂₂. Hence [F] ↦ BF is an isomorphism (the audited proof; finding 2).
- **No creation.** After an event the two loads are ⌊s/2⌋ and ⌈s/2⌉ with s = n_a + n_b (coin 0 gives the floor to
  the first port). Both lie between the two loads before the event and their sum is unchanged, so for every
  convex φ the value φ(n_a) + φ(n_b) does not increase. With φ(x) = max(x − q, 0) and φ(x) = max(q − x, 0): the
  numbers of excess and of deficit units relative to a filling q never increase. From one prepared pair at most
  one port has winding +1 at any time, and two disjoint regions never both hold it.
- **Split.** Every record is a sum of integer currents (the rule ends in a floor division), so every degree is an
  integer; two equal integers do not sum to 1 (`Protection.lean`, line 10). The unequal split 2 + (−1) is
  excluded by the preceding item, not by integrality. Each history holds the whole unit on one side; only the
  mean over the fair coin is ½ (`run.log` E1). The rule without its rounding term gives ½ and ½ (E2); the
  unchanged rule on a doubled register gives ½ and ½ of the prepared charge with certainty (E3, R13).
- **An angle field on the dual vertices.** Write F = a + ∂₂g with a = BᵀL⁺BF (the least-norm flow of the charge
  BF, L = BBᵀ) and g rational on the faces. Adding a face combination to F shifts g by integers, so θ = 2πg
  modulo 2π depends on BF only. Its vortex number around port v by the principal-value rule is −(BF)_v whenever
  |a_e| < ½ on every seam. For one pair the largest |a_e| is 11/30, 7/30 and 1/5 at distance 1, 2 and 3 (`run.log` C0). With
  two units on one port the bound can fail (first failure at event 1,746: excess 2, −1, −1 on ports 3, 4, 11,
  flow 3/5), and the vortex number then differs from the charge.

## 7. Several carriers: the declared cell complex of the cubic-routed graph (ledger entry A029)

Object (read only): `openwave-tuning/common_dynamics/vector_field/cell_complex.py`, `CELL_COMPLEX.json`,
`CELL_COMPLEX_PROOF.md` (first 16 hex digits of SHA-256: 7e760f4cc0cc9df1, bedeb0f30558b991, b3dff498ffba0c78).
Code here: `recompute_cells.py` rebuilds the graph from `links()` of `model/native_spatial_refinement/source.py`
and the cells from their definition, and does every lattice computation in plain integers (a reduction on unit
pivots and a dense Smith form, both my own); `replay_cells.py` calls the audited `build(2)` from this directory
(`run()`, which writes the receipt, is not called) and compares. The audited receipt is **confirmed**; three
statements are added (items c to e below). This section replaces nothing above: finding 12 concerns the
ring-glued screens, and its formula g − c + 1 gives the 41 of this graph (48 − 8 + 1). The audited directory is
being edited by its owner: during the replay other files were created and changed there (listed in
`replay_cells.json`); the three audited files kept their hashes, and the bytecode in that directory carries the
tag cpython-313 while the interpreter used here is cpython-312 with bytecode writing off.

| Claim of the receipt (size 2) | Recomputed |
| --- | --- |
| 96 vertices, 288 edges | 96 ports, 288 seams: 240 internal, 48 external; every port has degree 6 and exactly one external seam; the rebuilt seam list equals `graph(2).seams` |
| cycle-lattice rank 193 | connected; 288 − 96 + 1 = 193 |
| 160 internal triangles, 24 ribbons, 24 macro plaquettes | 160, 24, 24; all 208 boundaries are closed; equal to the audited matrix entry by entry in 208 of 208 columns, labels included |
| at most 4 cells and 16 edges per face | 4 and 16 (ribbon: 2 cells, at most 8 edges; plaquette: 4 cells, at most 16 edges) |
| internal faces alone leave 41 cycle classes | 152 Smith factors 1; 193 − 152 = 41 = 48 − 8 + 1, the cycle rank of the graph of cells |
| 190 unit factors with all faces | 190 factors 1, in ℤ²⁸⁸ and in cotree coordinates, by the sparse reduction and by the dense Smith form |
| ℤ³ with three unit torus periods | the period functionals vanish on all 208 faces; the three winding cycles are closed and have the identity period matrix; faces and winding cycles together have 193 Smith factors 1, so they generate the cycle lattice |

**Which cells are declared.** The 24 ribbons and the 24 macro plaquettes are not faces of any carrier. Each is a
disc glued along a closed path through two or four carriers, lifted along a breadth-first tree from port 0 inside
each carrier (`cell_complex.py` line 3: "Faces are declared geometric data"; PROOF lines 20 to 21: "they are not
deduced as a unique consequence of the bare graph"). The 160 triangles are the faces of the icosahedra. Also
declared: the cubic routing itself (`source.py` line 1: "A declared periodic cubic routing of the original
twelve-port cells"), the root and the tree of the lifts, and, for integer periods, three cut planes. The native
rule reads no face: the declared cells change no event; they fix which closed currents count as trivial.

**The dual-winding statement of row 1 on this complex.**

a. *Closed cuts: unchanged.* For every native history from a zero tape and every port set R, the record flux
   through the cut of R is G(R) − n(R). This is the update identity and needs no cell. The cut of a port has six
   seams. Every port of this graph has an external seam, so no pentagon is a cut: the single-carrier reading
   "relative charge = −(pentagon degree)" fails at any port that has exchanged a unit through its external seam.
   An internal seam lies in 2 to 14 of the 208 faces and an external seam in 1 or 5, so the complex is not a
   surface and has no dual graph; what carries over from one carrier is the lattice statement, not the circle
   field on a dual surface.
b. *The class of the tape is the charge plus three integers.* With one integer cocycle γ_k per axis (+1 on the
   positive external seams of axis k that cross one chosen plane), the map F ↦ (BF, γ₁·F, γ₂·F, γ₃·F) sends ℤ²⁸⁸
   onto {ρ ∈ ℤ⁹⁶ : Σρ = 0} ⊕ ℤ³ and its kernel is the lattice of the 208 face boundaries: 98 Smith factors 1 for
   the map, 95 for B alone, 288 − 98 = 190. Two tapes differ by face boundaries exactly when they have the same
   G − n and the same three integers. In invariant form: 0 → ℤ³ → ℤ²⁸⁸ / im ∂₂ → zero-sum ℤ⁹⁶ → 0, the homology
   sequence of the complex relative to its ports; on one carrier the left term is 0, with carrier faces only it
   is ℤ⁴¹.
c. *The three integers are not charges.* One load carried once round an axis by 6, 8 or 10 native moves returns
   with n = G and BF = 0 and leaves the integers (1,0,0), (0,1,0) or (0,0,1); its tape is the winding cycle. The
   native rule has no guard on them: they change by the current whenever a unit crosses a plane (on the torus
   screen of section 5 the corresponding records were nonzero at 8 of 9 neutral times).
d. *For a tape with nonzero charge the integers depend on the planes.* Moving the plane of axis 0 from between
   layers 1 and 0 to between layers 0 and 1 changes the integer by the charge of layer 1 (200 seeded tapes, 0
   exceptions); only the sequence in b is canonical. The audited period functional (±1/size on every external
   seam of an axis) needs no plane and equals the cut-plane integer on every cycle (193 fundamental cycles,
   values −1, 0, 1), but on an open tape it is not an integer: one unit across one external seam has period ½ at
   size 2. "Three integer periods" is exact for closed tapes; for the tape of a displaced load it needs the
   declared planes.
e. *Beyond the receipt.* The receipt's seven checks are the literal True behind `assert` statements. Line 134
   supports "no torsion" in the splitting by a rational rank; the integral statement holds (193 unit factors,
   above). Twelve seeded choices of root and tree give the same lattice of rank 190 and the same period matrix.
   The receipt executes size 2, where the positive and the negative neighbour of a cell coincide; at size 3
   (324 ports, 972 seams, 702 faces, cycle rank 649) the carrier faces alone give 513 unit factors and leave 136
   classes (162 − 27 + 1), all faces give 646 unit factors and leave ℤ³, and the class map has 326 unit factors.

Sentence for the claim gate, clause spatial_class, when the federation is used: "On the cubic-routed federation
with the declared ribbons and macro plaquettes the class of the tape is the relative charge G − n together with
three integers, the net numbers of units that crossed three declared cut planes of the torus. These are winding
records, not charges, and the native rule changes them at every crossing. Without the declared cells the class
holds 41 such integers at size 2 in place of three."
