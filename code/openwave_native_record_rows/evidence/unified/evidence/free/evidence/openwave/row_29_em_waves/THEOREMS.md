# Row 29. EM waves (Maxwell): theorems of the repository (worker TA, 2026-11-14)

**Current source scope (2026-11-14):** Earlier phase-access and physical-clock
obligations below belong to the named older candidates. The sealed common-v1
continuum theorem, same-history periodic first-hit join and actual common
cubic execution now supply the free-field reconstruction route. They do not
require primitive unitary updates or unread event registers. The standalone
submission bundle in row30/work/2026-11-14_codex_common_submission binds the
exact common source, full mathematical clauses, execution, controls and pin
change. Owner/platform/other-group acceptance remains open.


Read-only collection. Paths relative to `reverse-engineering-reality/` unless they start with `openwave-`. EOC =
`paper/executable_observer_construction.tex`, FL = `flagship/from_observer_consensus_to_standard_physics.tex`,
DP = `openwave-push/onecandidate/DP1_retention_schedule/CANDIDATE_DPRIME.md`. Row test: "Maxwell recovered;
transverse waves at c".

## 1. Chain

| # | Statement | Where | Kind |
| ---: | --- | --- | --- |
| 1 | Continuity of native currents; E[j] = d/2 given past and seam; mean law I − L/60 with multipliers in (0,1); "no positive rescaling of the attempt count turns this closed mean law into a wave law". | EOC 4095-4112 (prop:fields (i) to (iv)); `Lean/Screen/RegionalContinuity.lean` `regional_continuity` | proof + lean |
| 2 | The additive contrast j₀ + j₁ = n_a − n_b has zero sum round all nineteen carrier loops and twenty faces: every U(1) character of it is pure gauge. | `model/observer/README.md` "native loop-readout gate"; `model/observer/native_loop_readout.py` | exec |
| 3 | No continuous automorphism flow exists on the committed record algebra other than the identity. | `Lean/Dynamics/PublicAutomorphism.lean` `ContinuousPublicStarFlow.toAut_eq_refl` | lean |
| 4 | "Passive repair supplies no wave sector"; spectral dimension two; a bulk reading needs a radial observation map. | FL 1213-1219 | proof |
| 5 | Local Maxwell selection: a gauge-invariant, face-local, rotation-covariant quadratic seam action is const + (κ/2)\|CA\|²; for positive ε, β, a real linear autonomous evolution with Ḃ = −CE, conserved (ε\|E\|² + β\|B\|²)/2 and reversal (E even, B odd) is exactly Ė = (β/ε)CᵀB on closed flux. | FL 1124-1143 (thm:maxwell-selection); `model/maxwell_local_selection/PROOF.md` "Deriving Ampere evolution" | conditional |
| 6 | "The premises, field identification, coefficients and clock are inputs, and native repair does not select them." | FL 1150-1151 | premise statement |
| 7 | Periodic cubic stencil: C*C = \|d\|²I − dd*, rank P_T = 2 (two transverse modes), ω_a² = (4c²/a²)Σ sin²(ak_j/2), 0 ≤ c²\|k\|² − ω_a² ≤ c²a²Σk_j⁴/12; finite-cutoff quantum fields follow the cochain evolution exactly. | `model/openwave_field_observables/PROOF.md` "Finite quantum fields", "Two polarizations" | conditional (c = 1 declared, periodic gluing declared) |
| 8 | Given `LightSignalPremiseData` (oscillator lift `modeDynamics_is_lift`, `physicalFrequency_sq_eq_symbol`), the transverse mode obeys a second-order wave equation and the readout frequency is the committed nonnegative branch, one-Lipschitz, with a zero mode. | `Lean/Screen/LightSignalAdequacySurface.lean` `transverse_wave_second_order`, `physicalFrequency_eq_mode_branch`, `physicalFrequency_global_one_lipschitz` | lean (conditional on the structure) |
| 9 | The auxiliary lift q′ = p, p′ = −Lq gives q″ + Lq = 0 algebraically; it "does not provide a trajectory or identify its mnemonic prime notation with physical time". | `Lean/Screen/SeamCurrentAuxiliaryOscillatorLift.lean` `auxiliary_q_second_order` | lean |
| 10 | Temporal Maxwell packet: stationarity of the declared discrete action (kinetic term (h/2)\|E_n\|²) is equivalent to the step-scaled Ampère update and Gauss at every step. | `tracking/premise_register.json` PR-66 statement; `model/maxwell_local_selection/` | conditional |
| 11 | DP: row 29 and the Dirac field lie outside D′. | DP §4 | exec statement |

## 2. Premises

| Id | Text | State |
| --- | --- | --- |
| PR-05, PR-20 | uniform seams, fair coin (item 1) | open (declared); PR-20 forced given open data |
| PR-21 | second-order (oscillator) evolution of a transverse mode | open (field `modeDynamics_is_lift`) |
| PR-22 | readout frequency identified with the dilated carrier symbol at a chart scale | open (field `physicalFrequency_sq_eq_symbol`) |
| PR-43 | supplied continuous automorphism flow | open; impossible on committed records (item 3) |
| PR-66 | declared discrete kinetic term of the temporal Maxwell action | open (declared input) |
| unregistered | positive energy weights ε, β and Faraday kinematics Ḃ = −CE as a record law (item 5) | open: no native event writes face data (row design: "No event writes the offsets") |
| unregistered | radial observation map / cubic periodic gluing | open (shared with row 20) |
| PR-54 | source gauge field, current and action | open (obstructed on fixed native channels) |
| PR-15 | clock calibration for "at c" in physical units | open (import) |

## 3. The gap

On the committed native records the mean law is first order and real (1, 4); the loop sums are pure gauge (2);
no flow exists on committed records (3). Every Maxwell theorem on file (5, 7, 8, 10) is conditional on a field
carrier with its own kinematics. The smallest open set:

1. **A face/link register written by events and kept (Faraday kinematics, M2 of the row design).** Execution:
   on candidate T (or another candidate that writes transition data), a record experiment in which the change of a
   face circulation equals minus the circulation of the recorded currents on every event, with the unwritten-offset
   control failing. Theorem to prove first: for the native law applied to d = n_a − n_b − o_e with an event rule that
   writes o_e, the face sums of o obey ΔB_f = −(C j)_f exactly (finite incidence identity on 30 seams, 20 faces).
2. **Reversible second-order evolution (PR-21, PR-66, ε, β).** Theorem 5 reduces this to three clauses: positive
   conserved quadratic record energy, reversal, autonomy. On a dissipative committed law reversal fails (item 1), so
   either a theorem that some kept register of the candidate has a conserved positive quadratic form under the
   native events, or an **owner decision** to include a reversible field register (the D′-type choice of rows 3, 30).
3. **One speed.** Theorem to prove on a finite object: on the common candidate, the transverse speed of item 7 and
   the matter speed of row 30 (1/√48 on LOCAL12) are equal without a chosen coefficient (today a·b·ρ = 1 is declared).
4. Radial map (decision, shared with row 20) and PR-15/PR-22 for "c" in laboratory units (import).

## 4. Shared premises

PR-05, PR-20 with all rows. PR-21, PR-43 with rows 2, 3, 4, 30. PR-22 with rows 3, 30. PR-54 and the radial map with
row 20. PR-15 with rows 2, 3, 19, 20, 30. Offsets o_e (candidate T) with rows 4, 19.

## 5. Finite-window extension and operational-time scope (Codex MFR)

[New proof](work/2026-11-14_codex_record_poles/PROOF.md): for every finite graph and
nonnegative fixed total, uniform or fixed-IID seam repair has only eigenvalues
in [0,1]. Ordering configurations by Q=sum n² makes the matrix block triangular;
each same-Q seam block is an orthogonal projection. Any passive finite-window
function of committed letters has gap correlation μT_f K^nT_g1 and no new poles.
Positive one-state Feynman-Kac weights retain nonnegative-real poles; arbitrary
history-dependent weights and adaptive schedules are outside that corollary.

[Executed instance](work/2026-11-14_codex_record_poles/RECEIPT.json): 71,160 pinned
native calls, totals1,2,3, all27,240 branches and43,200 nonlinear windows retained.
Native poles are1,9/10,(55±sqrt(5))/60; the changed-law decoder detects four nonreal
poles and the independent planted control detects two. This extends E-03 and
E-08; it does not discharge PR-21, PR-43, PR-22 or the common-source bridge.

Event count is not assumed to be physical time. A relational time or an Euclidean
reconstruction must establish its time/algebra/state/readout transfer. The existing
`model/openwave_native_rest_energy/PROOF.md` already proves adaptive current-only
phase blindness from partial-injection Kraus histories, while allowing such an
alternate reconstruction. No theorem here excludes all observer-first physics.
The written new proof awaits another group's review; no new Lean build is claimed.

## 6. Positive reconstructed-observable route (finite sources inspected)

A committed diagonal record algebra and the algebra reconstructed from record
statistics are different objects. `RER/model/native_quantum_reconstruction/PROOF.md`
derives all4,356 matrix units of M66(R) from nonzero native path words (length<=4),
then its canonical complexification and the positive functional
rho_L=P^L/Tr(P^L) from closed-return probabilities. The current-only census and
signed tape attach these ideal counts with an explicit finite apparatus error;
the stated effect-statistic systematic bound is3.76596e-6. This is a source-derived
positive matrix functional, not a supplied density matrix.

`RER/model/native_modular_dynamics/PROOF.md` reconstructs ordered modular
correlations as finite signed/complex combinations of native event probabilities,
with exact binomial truncation bounds, and fixes the comparison t=Ls between
loop-length references. `RER/model/native_effective_fields/README.md` derives
G=sum_e(I-S_e), P=I-G/60 and H=-logP, including a nonquadratic interaction witness.
The spatial metric of its finite12-point support comes from current responses.
Thus linear mean-load superposition does not imply absence of every reconstructed
interaction, and a commutative-algebra no-flow theorem does not apply to M66(C).

The packages keep their remaining boundary explicit: signed reconstructed effect
weights need not be a single native measurement, the native returned working
state is diag(rho_L), and its evolution is not the modular unitary. Native attempts,
Euclidean parameter, modular parameter and physical duration are distinct. The
same-history operational algebra/state/clock/readout intertwiner remains owed.
This positive reconstruction chain should be used before any new field equation
is chosen; the fixed-IID window-pole theorem in section5 does not contradict it.

## 7. Native occupation Maxwell reconstruction

The [new written proof and finite experiment](work/2026-11-14_codex_record_maxwell/PROOF.md)
derive a Gaussian field-strength law from four native occupation-record
channels on the declared four-axis routing. The scalar Gaussian and higher
moments are source-derived by the shared occupation theorem. The local curl
decoder supplies tensor interpretation but cannot insert the nonlocal pole.
An explicit two-form reflection calculation gives a rank2 positive quotient
and exact pole E_a=2asinh(aQ/2)/a, equal to row30's massless native band.

The finite native branch-derived population covariance agrees with the target
within3.070e-5 and its fitted energy within4.694e-6 at three momenta; the
same decoder on IID-location resampling has no separated covariance. The
65,536 actual current events are separately retained and replayed. The long
population windows are analytical compression, not actual long tapes.

This discharges the assumption of independently supplied Gaussian higher
functions and oscillator evolution for this declared free-field
reconstruction. Routing, four channels, vector interpretation, local readout,
limiting regime and the physical time/intervention correspondence remain
explicit. No charged coupling, physical clock, owner pin or full row follows.
