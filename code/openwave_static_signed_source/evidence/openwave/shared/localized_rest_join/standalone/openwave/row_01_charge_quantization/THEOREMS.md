# Row 1. Charge quantization: theorems of the repository (worker TA, 2026-11-14)

Read-only collection. Paths relative to `reverse-engineering-reality/` unless they start with `openwave-`. EOC =
`paper/executable_observer_construction.tex`, FL = `flagship/from_observer_consensus_to_standard_physics.tex`,
FR/x = `paper/tex_fragments/x.tex`. Row test: "integer topological charge; a split e → e/2 + e/2 impossible".

## 1. Chain

| # | Statement | Where | Kind |
| ---: | --- | --- | --- |
| 1 | The native event on seam (a,b) with fair coin c moves j = (d − (2c−1)(\|d\| mod 2))/2 whole units, d = x_a − x_b; total load is conserved. | `model/observer/README.md` "The construction"; EOC 4095-4100 (prop:fields (i)) | proof |
| 2 | Uniform seam choice (1/60 per directed seam) and the fair coin are the reference of that event. Given A2/A3 projection data the weight 1/60 is forced. | `Lean/Screen/SeamCurrentHomogeneousAction.lean` `a2a3_directed_seam_weight_eq_one_sixtieth` (data: structure `A2A3DirectedSeamProjection`) | conditional (PR-20, PR-05) |
| 3 | Regional continuity: change of stored integer load in any port region = regional source − outward seam flux; internal seams cancel. | `Lean/Screen/RegionalContinuity.lean` `regional_continuity`, `seamBoundary_region_eq_neg_flux`, `global_continuity` | lean |
| 4 | With known start G and zero tape F, every native outcome gives F_e ↦ F_e + j, n ↦ n − jB_e, hence BF = G − n and Q(A) = Σ_{v∈A}(n_v − G_v) = −deg(z\|∂A) for the circle path z_e(t) = exp(2πiF_e t). | EOC 1003-1022 | proof (conditional on declared G and circle map) |
| 5 | ker B = im ∂₂: the homotopy class of the circle field equals the neutral integer port charge; a degree-zero boundary extends, a nonzero degree needs a puncture. | EOC 1013-1020 | proof |
| 6 | Equal daughter degrees of a primitive unit need 2d = ±1 with d integer: impossible. | `model/openwave_charge_protection/Protection.lean` `unit_degree_has_no_equal_halves`, `negative_unit_degree_has_no_equal_halves`; EOC 1021-1023; FL 2010-2013 | lean (arithmetic only; the sphere-degree theorem is not formalized, PROOF.md lines 81-84) |
| 7 | Integer degree does not forbid 1 = 2 + (−1). | `Protection.lean` `integer_winding_permits_unequal_signed_split` | lean |
| 8 | A closing sweep of thirty native contrast reads reconstructs the final loads: δ = n − n̂ is constant during the sweep and connectedness gives δ = 0, provided the true G and faithful closing contrasts. | EOC 1026-1035 | conditional |
| 9 | Executed: 5,400 native branches, 6,228 preparation/state pairs, 768 complete source/pointer laws, every timeout retained. | EOC 1057-1066; `code/openwave_charge_quantization/PROOF.md` (cited, not opened) | exec |
| 10 | Review R1: items 3 to 6 hold for every integer law (1,000 arbitrary integer tapes pass): they are identities of the alphabet, the graph and the decoder. | `openwave-push/mf/reviews/R1_charge/charge_row_packet/REVIEW.md` (cited in the row README) | exec |
| 11 | On the declared fifteen-weight native module, Q = T₃ + Y: spec dΓ(Q) ⊂ (1/3)ℤ and exp(2πi dΓ(Q))Γ(Z_colour) = I; on colour singlets the spectrum is integral. | `model/openwave_charge_protection/PROOF.md` "Current quantization on the complete fermion carrier" | proof (conditional on the declared generator) |
| 12 | Winding alone fixes no charge: on a plus-one texture W = I, Q = diag(0,1), so Q ≠ f(W) for every f. | same PROOF.md "Exact obstruction"; `Protection.lean` `same_degree_different_charge_forbids_dictionary` | proof + lean |
| 13 | The law-C integer face lift is no particle label: equal class mod 6 need not mean equal value. | `Lean/Dynamics/IntegerFaceLift.lean` `equal_class_need_not_equal_value`; charge_protection PROOF.md lines 86-110 | lean |
| 14 | On the declared Peter-Weyl source the cut flux equals enclosed central charge, F_S = ρ_S^Y + 1_{0∈S}Y_ref, and an arrival fixes ρ_y^Y = ρ_y^Q = −1. | `model/openwave_native_charge_flux/PROOF.md` eqs. (1)-(5) | proof (installed source) |
| 15 | A modulo-six accumulated current record is not a protected source charge or an identified Wilson flux. | FR/NATIVE_REPAIR_RESPONSE_AND_RECORDS 66-77 | proof |

## 2. Premises

| Id | Text | State |
| --- | --- | --- |
| PR-05 | counting path reference: uniform step weight is an input | open (declared input) |
| PR-20 | equal source-counting measure on the sixty directed seams | discharged by theorem `a2a3_directed_seam_weight_eq_one_sixtieth` given projection data; the data are open |
| unregistered | known initial occupation G (preparation) | realized only on an installed source (declared in EOC 1005); theorem 8 replaces it by faithful closing records |
| unregistered | faithful closing records | open (EOC 1032-1033: "does not authenticate replacement of the whole channel") |
| unregistered | the local circle map z_e(t) = exp(2πiF_e t) as the decoder | realized only on an installed source (declared; faithfulness proved, uniqueness not claimed, EOC 1019-1021) |
| PR-12 | matter candidate grammar; the reading unit = e (dictionary k = 6Y, cited by the row design) | open |
| unregistered | identification of the native integer load with electric charge | open (theorem 12 is an obstruction for winding; theorem 15 for the face lift) |
| unregistered | declared generator Q = T₃ + Y and colour-singlet restriction (theorem 11) | realized only on an installed source (`openwave_charge_protection`) |
| PR-54 | source gauge-field, current and action attachment (theorem 14) | open (register: obstructed on fixed native channels) |

## 3. The gap

The row's two clauses on the native process (integer, no equal halves) are theorems 3 and 6 and need no further
premise beyond PR-05/PR-20. Review R1 shows they cannot fail for any integer law. The open premises that stand
between the native records and "charge" are:

1. **Identification of the conserved integer with electric charge e.** Owner decision (which register of the
   candidate is called charge, and the unit dictionary of PR-12), followed by a theorem that the identified
   integer is the source of row 20's field on the same process. Theorem to prove on a finite object: on the
   twelve-port carrier with the native law, the mean seam current of a held unit equals −(1/6)∇(L⁺(δ_x − δ_ref))
   on every seam (EOC prop:fields (ii) composed with the mean law), so the integer of item 4 is the Gauss source of
   the row-20 record. This ties rows 1 and 20 to one integer.
2. **A present-state reading without supplied G.** Execution: the closing sweep of item 8 run on native records
   with the pinned function, all thirty contrasts read, and a forged-closing control; theorem 8 then removes G
   up to the faithful-closing premise.
3. **A clause that can fail.** Theorem to prove: there is a seam-local, load-conserving integer law that keeps
   every tape integral and gives, on the one-unit object, a different joint law of (tape, closing sweep) than the
   native law. If no such law exists among the laws that move at most \|d\| units, the clause is an identity and the
   owner decides whether to file it as one (decision, R1 finding 2).

## 4. Shared premises

PR-05 and PR-20 with rows 2, 3, 4, 19, 20, 29, 30 (every native-law statement). PR-12 and the charge
identification with row 20 ("at charge e"). PR-54 with rows 20 and 29.

## Current-only instrument and its existing source theorem

On one continuing native trajectory, repeat a seam until the independent coin changes.
Its initial difference is exactly `2 j_first + j_first_opposite`. Correct each tree
measurement by the previously recorded incidence, integrate the tree, and transport
forward by the tape to recover the terminal state modulo a common offset. A declared
zero-minimum initial promise fixes that offset; no numeric total or repeated preparation
is needed. [Written proof and replay](work/2026-11-14_codex_current_probe/PROOF.md).
The proof reuses finite regional continuity and the source endpoint analysis; it is
checked by another Codex agent, with other-party review pending. This implements a
present-state reader; it does not identify the integer with electric charge.

The two-current identity and coin-change stopping rule already occur in EOC
Theorem `thm:recovery` (vii), including an almost-sure identifying forest on the
fixed-total 66/78-state domains. The sealed replay is an implementation and a
general-integer spanning-tree extension modulo uniform shifts, not the discovery
of active tomography. Also retain `model/native_clock/configuration.py`: its
current-only controller reads the full 78-state sector, and restores or resets
binary inputs with complete failure flags. `model/native_thermal_return/PROOF.md`
composes those censuses into finite return experiments with explicit error and
native-attempt budgets. The replay's comparison with the older repeated-preparation
packet does not describe the limits of this stronger existing stack.

## Measured native relative-source topology

The [relative-source proof](../shared/row01_source_charge/PROOF.md) joins two
independent current-only restoring censuses to a common-choice native response.
The measured discrepancy `q=delta_p-delta_r` closes, in distribution, on the
signed native dipole kernel, independently of the common occupied background.
At record-certified synchronization its actual integrated current satisfies
`BF=q`. Independent branch enumeration and an exact rational reward identity
give mean response `B^T L+ q`; no field law is evolved to obtain it.

The integral face-boundary matrix has nineteen unit Smith factors and the
tree incidence has a unimodular minor. Hence `ker B=im D` over the integers,
and `C1/im D` is the rank-eleven relative group `H1(S²,V;Z)`, not the ordinary
first homology of the sphere. Regional evaluation makes a unit source class
primitive and excludes two equal sharp half classes. Native three-move loops
around a face realize its comparison relation while returning the relative
source; the common background need not return. The same-population swap/hold
control never synchronizes, and the record-read half-averaging control detects
fractional daughters. The literal static clauses hold on this native object;
the charge dictionary and in-platform adjudication remain explicit.
