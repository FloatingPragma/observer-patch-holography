# Row 30. Quantum wave equation (Klein-Gordon): theorems of the repository (worker TA, 2026-11-14)

**Current source scope (2026-11-14):** Earlier phase-access and physical-clock
obligations below belong to the named older candidates. The sealed common-v1
continuum theorem, same-history periodic first-hit join and actual common
cubic execution now supply the free-field reconstruction route. They do not
require primitive unitary updates or unread event registers. The standalone
submission bundle in row30/work/2026-11-14_codex_common_submission binds the
exact common source, full mathematical clauses, execution, controls and pin
change. Owner/platform/other-group acceptance remains open.


Read-only collection. Paths relative to `reverse-engineering-reality/` unless they start with `openwave-`. EOC =
`paper/executable_observer_construction.tex`, CC = `model/observer/COHERENT_CONTINUATION.md`, DP =
`openwave-push/onecandidate/DP1_retention_schedule/CANDIDATE_DPRIME.md`. Row test: "emergent ω² = k² + m² for
matter waves".

## 1. Chain

| # | Statement | Where | Kind |
| ---: | --- | --- | --- |
| 1 | Native event and continuity; mean law I − L/60 with multipliers in (0,1); no rescaling of the attempt count gives a wave law. | EOC 4095-4112 (prop:fields); `Lean/Screen/RegionalContinuity.lean` `regional_continuity` | proof + lean |
| 2 | Unread one-load channel: real strictly positive spectrum; Φ − I is self-adjoint in Hilbert-Schmidt space while −i[H,·] is skew, so Φ is no unitary step. Record laws see only occupation diagonals. | `model/openwave_native_rest_energy/PROOF.md` first two sections | proof |
| 3 | The only continuous automorphism flow on committed records is the identity. | `Lean/Dynamics/PublicAutomorphism.lean` `ContinuousPublicStarFlow.toAut_eq_refl` | lean |
| 4 | Given a supplied flow on a represented matrix algebra, it is conjugation by a self-adjoint generator unique up to a real scalar. | `Lean/Dynamics/StoneConverse.lean` `starDerivation_exists_selfAdjoint_generator`, `selfAdjoint_generator_unique_up_to_real_scalar` | lean (conditional on PR-43) |
| 5 | Reflection walk on the native branches: A†SA = K, A†WⁿA = T_n(K), U_{n+1} − 2U_n + U_{n−1} = −U_n L/30, cos ω = 1 − λ/60 (massless at p = 0). Squared amplitudes and mean loads do not obey this recurrence. | CC lines 33-74; imported Higuchi-Segawa-Suzuki theorem | exec + imported theorem |
| 6 | Nonselection: W_z share the first coherent event; a record-phase twist keeps every committed branch and changes Tr(D³) by −1/900: "even a complete committed quantum instrument does not determine its coherent continuation". | CC lines 82-114 | exec |
| 7 | With reuse (A) and blind maximal redraw (B) the redraw is R up to phase; the mass redraw Π + z(1 − Π), z = −e^{2iμ}, drops the second clause of B; μ supplied. | `openwave-push/steps/S1_reflection_origin/SOLUTION.md` Theorem 4.1; DP §2 | proof (outside RER) |
| 8 | Scalar mass-coin source W_hm = S(cos(hm)R − i sin(hm)I): RP₀ = 2ee* − P₀, RH_a + H_aR = 0, g_a*g_b = δ_ab/48, ⟨e, e^{tK_m(p)}e⟩ = cos(Et) − i(m/E) sin(Et) with E² = m² + \|p\|²/48. | EOC 1712-1723 | proof (integer algebra) |
| 9 | (∂_t² − Δ/48 + m²)φ = 0 in declared units; X = √48 x gives unit speed "without proving a shared physical clock". | EOC 1735-1742 | proof |
| 10 | Mass-channel operator H_a(p) = (v/a)Σ sin(ap_j)γ_j + m₀β₀ + Σ m_jβ_j obeys H_a² = [(v/a)²Σ sin²(ap_j) + (m₀ ± \|m\|)²]I in each sector, with 0 ≤ v²\|p\|² − (v/a)²Σ sin² ≤ v²a²Σp_j⁴/3. | `model/openwave_field_observables/PROOF.md` "Massive wave records" | conditional (supplied coefficients) |
| 11 | Kinetic/analyzer source H(p) = M + ΣΓ_j sin p_j with Clifford Γ and {Γ_j, M} = 0, ‖H‖ < 2; fourth-order integrator error ≤ t⁵2⁵/(120n⁴); the electron scale is empirical. | `model/openwave_platform_kernel/PROOF.md` | conditional |
| 12 | Lab frequency of an internal rate is not forced without a declared tick. | `Lean/Thermodynamics/PhysicalCalibrationImport.lean` `labFrequency_not_forced` | lean |

## 2. Premises

| Id | Text | State |
| --- | --- | --- |
| PR-05, PR-20 | uniform seams, fair coin; LOCAL12 scheduler for the wave rows (DP §2 table) | open (declared) |
| PR-43 | supplied continuous automorphism flow | open; impossible on committed records (item 3); realized on D′ by W |
| PR-21 | second-order evolution shape | open; realized on D′ by the reflection (item 5) |
| unregistered (DP §2 A, B) | reuse of the unread choice register; blind maximal redraw | open; given A and B, R is fixed (item 7) |
| unregistered | mass angle μ (rule μ = π/3 from ε = 1/2 is a rule, row design) | open |
| unregistered | maintenance-to-energy map (HYP 113-120) | open |
| PR-22 | readout frequency identified with a physical frequency | open |
| PR-15 | clock and energy calibration | open (import) |
| PR-64, PR-65 | coherent phase instrument (Ramsey quadratures of item 8), source-produced phase preparation | realized only on an installed source; PR-65 obstructed for committed repair |
| unregistered | kinetic/analyzer source and mass-channel coefficients (items 10, 11) | realized only on an installed source |

## 3. The gap

On committed records the law is the real kernel I − L/60 (or I − L/12): no phase per step (1 to 3). The
dispersion ω² = m² + p²/48 exists on the unread schedule with a supplied angle (5, 7, 8). The smallest open set:

1. **Unread choice register as part of the memory law (PR-43, PR-21, A, B).** By items 3 and 6 no theorem derives
   it from committed records. **Owner decision** (F0), shared with rows 2, 3, 4.
2. **The mass angle μ.** With clause B only, the law is massless. Theorem or decision as stated in row 2: one
   native quantity must fix μ (for example 1 − cos μ = ε with ε the all-read movement excess of the same object),
   proved on the 144-coordinate Bloch cell; otherwise μ stays a supplied number and the gap clause is not met.
3. **One operational speed.** Use the same record-derived rulers and clock to compare the long-wavelength
   matter dispersion slope with the light-signal slope on one candidate, with finite-cutoff errors bounded.
   A discrete graph-cone maximum is a separate causal bound; equality of that maximum to the infrared slope
   1/√48 is not required by the criterion and does not follow from a continuum wave limit. Execution:
   ω(p) from terminal phase-sensitive records at separated cells, a common light comparison, and a declared
   10 % change of μ as a sensitivity control. Rescaling a clock cannot repair a ratio of unequal speeds.
4. PR-22 and PR-15 for physical units (import).

## 4. Shared premises

PR-05, PR-20 with all rows. PR-43, PR-21, A and B with rows 2, 3, 4, 29. PR-22 with rows 3, 29. PR-15 with rows 2,
3, 19, 20, 29. PR-64, PR-65 with rows 2, 3. Mass angle / maintenance-to-energy map with rows 2, 3, 19.

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

## 7. Positive controlled slice refinement

[Conditional proof](work/2026-11-14_codex_slice_refinement/PROOF.md), checked by the
root in a [same-group review](reviews/20261114_codex_root_slice_refinement/REVIEW.md).
On E1's declared four-axis port-diagonal routing with external seam weight a^-2,
the native record slice rate is E_a=2asinh(a sqrt(lambda+q_a)/2)/a. For all bands,
0<=lambda+|p|²-E_a²<=a²[sum p_j^4+(lambda+|p|²)²]/12, with native
lambda=0,5-sqrt(5),6,5+sqrt(5). All infrared speed-squared deficits are bounded by
a²lambda/6; the normalized covariance has a uniform O(a²) scalar limit away from
the zero mode. This extends the old weight scan with explicit error estimates.

No new process was executed. The routing, refinement scheduler, coordinate ruler,
Euclidean slice and decoder remain supplied. This is a positive conditional
source-record theorem, not an operational physical-time or full matter-wave claim.
