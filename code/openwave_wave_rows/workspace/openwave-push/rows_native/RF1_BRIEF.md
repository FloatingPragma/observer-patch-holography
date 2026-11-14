# RF1: the dual repair of the native tape (rows 20, 2, 8, 5, 22)

"Native reciprocal field attachment", item 4 of the joint plan
`/Users/muellerberndt/Projects/oph-meta/openwave-tuning/next/ROW1_NATIVE/NEXT_TEN.md`, which this side leads.

Directory: `/Users/muellerberndt/Projects/oph-meta/openwave-push/rows_native/RF1_dual_repair_field/` (create it).
Write only there. Never write to the theory repository
(`/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/`, "RER", read only), to the simulator, or to
`openwave-tuning/`. Interpreter `/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/.venv/bin/python`
with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1`; each script
under five minutes of CPU (the machine is heavily loaded); never `pkill` or `killall`. Writing: statements of
fact, no narration, never the word "honest". Budget 70 to 90 minutes; report the decisive finding as soon as you
have it.

## Read first, in this order

1. `/Users/muellerberndt/Projects/oph-meta/openwave-push/REUSE_BRIEF.md` (existing-result lookup before any
   derivation: binding).
2. `openwave-tuning/next/ROW1_NATIVE/ROW1_CLAIM_GATE.md` and `routes/spatial/METHOD.md` there: the other team's
   finished row-1 chain. On the literal native law the current tape F obeys BF = G − n; read on the dual graph of
   the carrier its winding around the pentagon of port v is G_v − n_v; relative charge q = n − G.
3. `/Users/muellerberndt/Projects/oph-meta/openwave-push/issues/N20_native_green_field/SOLUTION.md`: this side's
   finished result. Under the unchanged native law the mean of the tape given the present positions is the
   lattice Coulomb field BᵀL⁺(δ_x − δ_y) only under a uniform prior over the creation; one history carries a
   divergence-free remainder of order one; on the repository's cubic routing the potential falls as 1/r with
   coefficient 1/(6π). A maintainer would refuse "an ensemble mean instead of a field in one run" and "a
   potential that nothing feels".
4. `/Users/muellerberndt/Projects/oph-meta/openwave-push/survey/REGISTER_AND_TRANSPORT.md` sections 0 and 2: where
   the simulator is on disk; its transported-contrast law (integer loads that read integer seam offsets,
   s = x_a + x_b + b_e; holonomy charges of the offsets; release ½|d̄|² = ½ qᵀD⁺q per sweep; phase 2: unit twists
   of the offsets with a declared cost); code entry points.

## The construction to test (a hypothesis: prove, correct or refute every statement)

One carrier has 12 ports, 30 seams, 20 oriented faces; ∂₂ is the face boundary (30 × 20), B the incidence
(12 × 30), B∂₂ = 0 and ker B = im ∂₂ as lattices. Apply the simulator's transported-contrast law, unchanged, to
the dual graph of the carrier, with the native tape as its offsets: integer face loads x_f on the 20 faces, one
dual seam per seam, offset b_e := F_e. Then:

1. The holonomy of the offsets around the dual plaquette of port v is (BF)_v = G_v − n_v: the topological
   charges of the dual law are exactly the native units of the row-1 chain. No new register: the offsets are the
   record the native law already writes.
2. The real-valued minimiser of |F − ∂₂x|² over x leaves the residual BᵀL⁺(G − n): the Coulomb field of the
   present charges, in one history, with no prior and no average. The divergence-free remainder that N20 found
   is exactly what the face loads absorb.
3. The minimal value is ½ qᵀL⁺q (for a unit and its origin: half the effective resistance; on the carrier
   11/60, 7/30, 1/4 at distances 1, 2, 3), which on the cubic routing behaves as a constant minus c/r.
4. A native hop of a unit across seam e changes F_e by one: it is a unit twist of the dual offsets, the
   simulator's phase-2 move. The unchanged native law is the twist rule that does not read the field (fair coin).
5. A reciprocal law needs one further rule: the hop reads the dual field. Candidates with no new constant: the
   coin of an odd-contrast event is replaced by the sign of the change of the dual release (hop accepted when it
   lowers ½|d̄|²), or by the simulator's own phase-2 acceptance rule. This changes the fair coin (registered
   premises PR-05, PR-20) and is therefore a new candidate; say so wherever it is used.

## Tasks

A. **Lookup** (reuse brief). In RER search `model/` and `claims/` for face loads, face registers, dual graph,
   2-cells, plaquette, "transverse", `maxwell_local_selection`, `native_circulation`, `native_electric_gauss`,
   `flux_response`, `openwave_current_coulomb`, `coherent_defect_wave/native_routing_join`. In the simulator read
   the transported-contrast engine and its phase-2 rule (quote the rule and the cost table as coded). Table
   "Existing results" with the action per item. In this push read `issues/G9_offset_charges/` and
   `bridgework/BR15_repair_on_defects/` if they have delivered files (two workers on the simulator's law with
   offsets as an independent register: do not repeat them; your offsets are the native tape).

B. **Exact statements on one carrier** (plain integers and Fractions;
   `sys.path.insert(0, "/Users/muellerberndt/Projects/oph-meta/openwave-push/kernel"); import ophk;
   c = ophk.carrier()` gives seams and oriented faces). Prove (1), (2), (3) as lemmas and verify them on every
   one-load history of length up to 4 from one start port and on seeded long histories: the residual after exact
   real minimisation equals BᵀL⁺(G − n) in every history; its norm depends only on the present port.

C. **The integer law.** Implement the dual transported-contrast law with integer face loads yourself (the
   simulator's rule, on the dual graph, offsets = the tape; the port law unchanged and running at the same time;
   events drawn uniformly over the 30 seams for ports and the 30 dual seams for faces; state the interleaving).
   With the unit held at separation 1, 2, 3 from its origin (port law paused) measure in the stationary state:
   the time-averaged transported contrast on every seam against BᵀL⁺(G − n) (seam by seam, mean ± standard
   error); the release per sweep against ½ qᵀL⁺q; the circulation of the face loads around the charged pentagons
   (a persistent current? its sign against the sign of the charge; its magnitude). With no charge (unit back at
   its origin after a closed walk that left loop currents on the tape) the face loads must absorb the tape
   completely: residual zero: check.

D. **Three dimensions.** The repository's cubic routing of carriers (`model/native_spatial_refinement/source.py`
   in RER: read how carriers are joined; N20's `n20lib.py` in this push builds it: reuse that builder, read
   only). State what the 2-cells are on a routed federation: carrier faces, and which closed loops through glue
   seams have no face. The tape's classes on those loops are not absorbed by face loads: count them on the
   smallest cubic block and say what they are physically. With real minimisation over the available face loads,
   compute the residual field energy of a prepared pair against separation along an axis and a diagonal on the
   largest block that fits the budget (sparse least squares), fit E(r) = E_∞ − c/r, and compare c with N20's
   coefficient 1/(6π) and with the exact lattice Green function of that block. If the missing 2-cells spoil the
   1/r law, say what fills them (which plaquettes one must add, and whether the repository has them).

E. **The reciprocal rule.** Implement candidate (5) in its two forms on one carrier and on a ring screen patch
   from `ophk.flat` (a unit and a real compensating deficit on a balanced background, so that both are native
   objects). Measure: the drift of the pair against separation (mean approach per sweep ± standard error, seeds
   stated); the sign structure (unlike attract; two excess units repel); the annihilation time against the
   unchanged law; and two controls that must show no drift: the face loads frozen at zero (held field) and the
   fair coin (zero coupling). Is the force proportional to the gradient of ½ qᵀL⁺q? State the mobility that
   relates them and whether any free constant entered (none should; if one did, name it).

F. **Rows**, each as a block "statement / executed evidence / supplied pieces by registered premise / what a
   strict maintainer would grant" (precedents in
   `/Users/muellerberndt/Projects/oph-meta/openwave-push/platform/PRECEDENTS.md`, rows 20, 2, 8, 5, 22):
   - row 20: a 1/r potential in one history on the native units; a force with the reciprocal rule; the strength
     as the pure number that comes out (compare with the values in N20's SOLUTION.md; α is not derived unless the
     number is 1/137.036 in the right unit: say which unit);
   - row 2: the self term ½ L⁺_vv and the maintenance release of one unit as a rest energy in native units:
     finite or divergent with size in two and three dimensions; what anchors it;
   - row 8: what the field energy does at annihilation;
   - row 5: the circulation of task C: does a stationary charged unit carry a persistent circulating dual
     current, with which moment;
   - row 22: the two-dimensional screen gives a logarithm, not a linear law: confirm.

G. **What is declared.** Every piece of the construction that is not the literal native law (the face loads; the
   dual application of the simulator's rule; the interleaving rate; the reciprocal rule), with the registered
   premise each falls under and whether the theory repository already admits it anywhere.

## Deliver

`SOLUTION.md` (lemmas with proofs, tables), `ROW_SHEET.md` (task F), `run.py`, `verify.py` with logs and
`results.json`, `FOR_CODEX.md` (one page for the other team: the construction, the three lemmas, the numbers,
what is declared, the file to run).

Final message to the coordinator, under 450 words: which of (1) to (5) hold, the decisive numbers (seam-by-seam
agreement, c against 1/(6π), drift against separation with controls), the circulation result, what is declared.

## Addendum (coordinator, 10:50): what the other agents have on the same item

- The external auditor has frozen a reciprocal feedback package on the literal native law:
  `/Users/muellerberndt/Projects/oph-meta/openwave-tuning/auditor/native_field_bridge/` (FINAL_MANIFEST.json,
  ROOT_UNCAPPED_REVIEW.json, ROOT_CAPPED_REVIEW.json; ledger entries A024 to A027 in
  `openwave-tuning/sync/AUDITOR.md`). It is a current-driven feedback with a declared energy-dependent scheduler
  and a local Poisson relaxation; its κ, β and controller energy are installed inputs. Read its method text
  (read only) in task A and state, in `FOR_CODEX.md`, the exact relation: which of its installed pieces the dual
  repair replaces by a quantity of the balancing rule itself (the energy ½ qᵀL⁺q; the relaxation by face loads),
  and which it does not (the rule by which a hop reads the field; any rate ratio between port and face events).
- The coordinator's check `bridgework/COORD_checks/dual_lattice_check.log` already has: ker B = im ∂₂ as
  lattices on the carrier; effective resistances 11/30, 7/15, 1/2 at distances 1, 2, 3.
- On several glued carriers the cycles through glue seams are not face boundaries. Task D has to say what
  plays the role of the face loads for those loops.
