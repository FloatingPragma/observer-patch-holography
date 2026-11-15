# Row 4. Particle stability (Derrick escape): theorems of the repository (worker TA, 2026-11-14)

Read-only collection. Paths relative to `reverse-engineering-reality/` unless they start with `openwave-`. EOC =
`paper/executable_observer_construction.tex`, FL = `flagship/from_observer_consensus_to_standard_physics.tex`, HYP =
`docs/OPH_EXECUTABLE_HYPOTHESIS.md`, DP = `openwave-push/onecandidate/DP1_retention_schedule/CANDIDATE_DPRIME.md`.
Row test: "the particle state persists in free evolution".

## 1. Chain

| # | Statement | Where | Kind |
| ---: | --- | --- | --- |
| 1 | Native event moves whole units; regional and global continuity: a unit is never lost or divided. | `Lean/Screen/RegionalContinuity.lean` `regional_continuity`, `global_continuity`; EOC 4095-4098 (prop:fields (i)) | lean + proof |
| 2 | Executed: 0 units lost or divided in 3,464 histories and 32,424,192 events. | row README; `openwave-push/packets/ROW_04` (cited) | exec |
| 3 | The one-unit position law is the reversible kernel I − L/60 (uniform seams) or I − L/12 (LOCAL12): real spectrum, every nonconstant mode decays; the unit spreads. | EOC 4095-4112 (prop:fields (iii), (iv)); `model/openwave_native_rest_energy/PROOF.md` first section | proof |
| 4 | The counting reference has constant energy: the native reference holds no energy minimum for a unit. | `Lean/Thermodynamics/NativeBinaryThermodynamics.lean` `counting_calibration_energy_constant` | lean |
| 5 | Frustration: on a fixed graph with positive mismatch weights and fixed edge offsets, nonzero cycle holonomy forces a strictly positive minimum mismatch that no vertex repair removes. | FL 845-848 (cites the consensus paper) | conditional (offsets fixed) |
| 6 | Load differences telescope round all twenty faces for all 78 states: the native gradient map has no holonomy. | EOC 5385-5387 | proof |
| 7 | Reflection walk W = S·R: a compact span with W u = v, W v = 2κv − u stays in four cells for all time; with every register read the same start leaves (1/36 at two events, review R7). | row README "Mathematics"; `openwave-push/mf/reviews/R7_row04_stability/REVIEW.md` (cited) | exec (installed reflection) |
| 8 | Installed trapped source with the full photon Fock space: an isolated spin-half band below 16 + 1/16, next level ≥ 17, Temple bound on E_∞, photon-tail bound. | `model/openwave_infinite_photon/PROOF.md` "Exact isolated spin-half band", "Exact moment shortcut" | conditional |
| 9 | Installed kinetic/analyzer source: a compact charged packet stays outside its cube with probability ≤ (√(η/(1−η)) + \|t\|b)² for finite time (no soliton premise); ≤ 1.103e-7 for \|t\| ≤ 1e9. | `model/openwave_localized_rest_packet/PROOF.md` "Persistence, volume and disturbances" | conditional |
| 10 | Installed neutral band of the link-penalty-two source: closed-source stationarity, flat at zero hopping (object of row 19). | `model/openwave_neutral_stationary/PROOF.md` "Complete records, localization and calibration" | conditional |
| 11 | HYP target: "stable massive matter as localized, topologically protected native defects"; the sector and boundary invariant under admitted repairs must be exhibited. | HYP 103-112 | premise (program statement) |

## 2. Premises

| Id | Text | State |
| --- | --- | --- |
| PR-05, PR-20 | uniform seams, fair coin | open (declared); PR-20 forced by theorem given open data |
| unregistered | fixed integer edge offsets o_e with nonzero face holonomy (theorem 5; candidate T) | open; produced only in a quench of another worker's run (BR15, row design) |
| PR-21, PR-43, DP clauses A, B | unread choice register and reflection (theorem 7) | open for that candidate; the continuous-star-flow theorem on a finite commutative algebra does not rule out an enlarged coherent memory or discrete active correction |
| unregistered | spatial trap, field strength, absolute-kernel law (theorem 8) | realized only on an installed source (`openwave_infinite_photon` PROOF.md lines 111-114) |
| unregistered | kinetic/analyzer source and positive-band preparation (theorem 9) | realized only on an installed source |
| unregistered | preparation of the compact state | realized only on an installed source (native preparation selection gate fails, row README) |
| unregistered | co-moving bookkeeping of the observer | open (design choice of the observer) |

## 3. The gap

Integrity of the unit is proved (1, 2). Under the independent uniform schedule,
the unconditioned one-unit position law spreads (3); the counting reference in
(4) is one specified reference, not all possible record energies. The spatial
persistence results (7 to 10) use installed elements. The active-memory results
below establish a different, positive kind of persistence on explicit enlarged
sources. Three previously proposed spatial constructions are:

1. **Offsets written by the process (candidate T).** Theorem to prove on a finite object: on the twelve-port
   carrier with integer offsets o of one nonzero face holonomy, the native law applied to d = n_a − n_b − o_e has a
   unique stationary law whose one-unit marginal assigns the frustrated face's three ports a probability bounded
   below by a constant c > 1/4 independent of the start, and the co-moving width stays bounded at all times.
   Then an **execution**: a record experiment in which events write and keep the offsets (no offset supplied),
   with the unfrustrated control spreading.
2. **The unread choice register (D′).** Research construction is authorized; selection of a column pin remains separate. Shared with rows 3 and 30. Even then the
   preparation of the compact span is supplied; a theorem would have to show that a native preparation from the
   draw state has nonzero overlap with the compact span that does not decay.
3. **An installed trap or kinetic source.** May be constructed, but its installed physical source retains that classification.

Charged and massive together: the localized object must carry row 1's integer and row 2's energy on one process.

## 4. Shared premises

PR-05, PR-20 with all rows. PR-21, PR-43, A and B with rows 2, 3, 29, 30. Offsets o_e (candidate T) with rows 19
and 29. Localization with rows 2 and 19. Kinetic/analyzer source with rows 2 and 30.

## Positive finite memory stack

`model/defect_protection/README.md` and `protection.py` construct the binary
66-state C/E source with one coherent ordering bit. Actual token transport
determines each branch's bit flip. Every complete branch commutes with
`Q=I66 tensor X`, so its two sharp sectors persist for every finite adaptive
seam policy. The all-seam algebra is exactly `M66 direct_sum M66`, with
commutant `span{I,Q}`; the classical 132-state cover is connected. Thus the
classical load-only invariant theorem does not exclude this quantum charge.
The uniform source with either sharp memory character is stationary. Its
source and classical tape laws agree in the two sectors, so the memory readout
and its preservation are essential. Dephasing or resetting the ordering bit
destroys this distinction. This is an active coupled memory, not yet a
spatially localized massive electron.

`model/common_source/ENCODED_MEMORY.md` supplies a second constructive result:
record-controlled recovery gives `R_h L_h V=sqrt(p_h)V`, preserving an unknown
logical state and arbitrary reference through continuing native events. Its
six-dimensional matching code attains the proved all-seam capacity. The
encoder, recovery pulses, finite numerical diagnostics and service assumptions
remain explicit; passive native averaging alone is not claimed to implement
them. `Lean/Geometry/SourceRecordProtection.lean` characterizes linear invariants
of real scalar pair averaging only and explicitly does not exclude these
active or nonlinear constructions.

`model/defect_energy/README.md` derives the state-dependent action of the same
protected ordering memory from its actual transport parity. Histories with
the same raw-event multiset can have different memory actions, so a scalar
cost per raw event loses information that the operator action retains. The
remaining join is a source-localized core, energy and clock for this same
protected sector, rather than a proof that all native memory must dissipate.
