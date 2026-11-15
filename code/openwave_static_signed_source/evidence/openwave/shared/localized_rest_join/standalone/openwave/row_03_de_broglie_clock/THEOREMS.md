# Row 3. de Broglie clock (Zitterbewegung): theorems of the repository (worker TA, 2026-11-14)

Read-only collection. Paths relative to `reverse-engineering-reality/` unless they start with `openwave-`. EOC =
`paper/executable_observer_construction.tex`, CC = `model/observer/COHERENT_CONTINUATION.md`, DP =
`openwave-push/onecandidate/DP1_retention_schedule/CANDIDATE_DPRIME.md`. Row test: "self-starting
frequency-rigid internal oscillation at ω = mc²/ℏ".

## 1. Chain

| # | Statement | Where | Kind |
| ---: | --- | --- | --- |
| 1 | Native event, continuity, counting reference in detailed balance (rows 1, 2). | `Lean/Screen/RegionalContinuity.lean` `regional_continuity`; `Lean/Thermodynamics/NativeBinaryThermodynamics.lean` `fairSwap_detailedBalance` | lean |
| 2 | Every star automorphism of a finite commutative record algebra is a label permutation; every pointwise-continuous real-parameter group of them is the identity. | `Lean/Dynamics/PublicAutomorphism.lean` `publicStarAutomorphism_is_labelPermutation`, `ContinuousPublicStarFlow.toAut_eq_refl` | lean |
| 3 | The unread one-load channel has real strictly positive eigenvalues (phase [0,0] per event); record laws are blind to offdiagonal phases. | `model/openwave_native_rest_energy/PROOF.md` first two sections | proof |
| 4 | Mean law: multipliers in (0,1); the second-difference coefficient has the sign opposite to an undamped oscillator. | EOC 4095-4112 (prop:fields (iv)) | proof |
| 5 | For the continuous-time finite Markov memory generator used in the cited construction, damping is at least frequency·tan(π/N); this is not a statement about every discrete-time finite observer. | `openwave-push/blockers/X3b_tape_momentum/VERDICT.md` line 23 | proof (outside RER) |
| 6 | Energy-invariant source/battery collisions: if ρ, U commute with H_s + H_b and every battery effect commutes with H_b, each conditional source state commutes with H_s, so every source expectation is constant in each branch. | `model/openwave_clock_phase_resource/PROOF.md` "Conditional phase obstruction"; algebra in `Lean/ObserverPatchHolography/ClockPhaseCovariance.lean` `invariant_branch`, `invariant_sum` | conditional + lean (algebraic step) |
| 7 | Declared collective collision channel: the k-th coherence band carries exactly the phase exp(−ikmδ) with a certified nonzero linewidth; from the phase-invariant vacuum ⟨a⟩ = 0 at every step. | `model/openwave_autonomous_clock_channel/PROOF.md` "Exact frequency and finite linewidth", "Onset" | conditional (exec certificate) |
| 8 | Reflection walk on the binary shell: S² = I, A†SA = K (native kernel), W = S(2Π − I), A†WⁿA = T_n(K); Laplacian modes have cos ω = 1 − λ/60. | CC lines 17-74; imported Higuchi-Segawa-Suzuki spectral mapping | exec + imported theorem |
| 9 | Nonselection: W_z = S[Π + z(1 − Π)] share the first coherent event; A†W_z²A = K² + z(I − K²). A record-phase twist keeps every committed C/E branch and changes Tr(D³) by −1/900. | CC lines 82-114 | exec (exact fractions) |
| 10 | Under slot-permutation covariance the real redraws are ±1 and ±R; with reuse (A) and blind maximal redraw (B) the redraw is R up to a phase. Under the carrier symmetry alone the commutant has dimension 60. | `openwave-push/steps/S1_reflection_origin/SOLUTION.md` Theorem 4.1; DP §2 | proof (outside RER) |
| 11 | On D′ the all-read schedule equals the literal law; the choice-unread schedule equals W = S·R; they differ from the second event (643/900 against 101/120). | DP §1 table; `openwave-push/onecandidate/DP1_retention_schedule/run_dprime.py` (cited) | exec |
| 12 | Installed kinetic/analyzer packet: phase within d/(1 − \|t\|d) of m; no passive onset. | `model/openwave_localized_rest_packet/PROOF.md` "Autonomous phase records" | conditional |
| 13 | Clock frequency in laboratory units is not forced without a declared tick. | `Lean/Thermodynamics/PhysicalCalibrationImport.lean` `labFrequency_not_forced` | lean |

## 2. Premises

| Id | Text | State |
| --- | --- | --- |
| PR-05, PR-20 | uniform seams, fair coin | open (declared); PR-20 forced by theorem given open data |
| PR-43 | pointwise-continuous public star-automorphism flow, supplied | open; theorem 2 proves it cannot be emitted by the committed record algebra; it can hold only on an algebra with unread registers |
| PR-21 | second-order (oscillator) evolution shape | open; on D′ realized by the reflection R (theorem 8) |
| unregistered (DP §2 A, B) | reuse of the choice register without reading; blind maximal redraw | open; given A and B the redraw is fixed (theorem 10); A and B themselves are design choices |
| PR-22 | identification of a readout frequency with a physical frequency | open (`Lean/Screen/LightSignalAdequacySurface.lean` field `physicalFrequency_sq_eq_symbol` of `LightSignalPremiseData`) |
| PR-15 | clock and energy calibration anchors | open (import) |
| PR-64, PR-65 | coherent phase instrument; source-produced phase preparation | realized only on an installed source; PR-65 obstructed for committed repair |
| unregistered | maintenance-to-energy map / mass angle μ (HYP 113-120) | open |
| unregistered | nonpassive resource for onset (theorem 6) | realized only on an installed source (fresh nonpassive batteries, `openwave_autonomous_clock_channel`) |

## 3. The gap

1. **Unread choice register in the memory law (PR-43, PR-21, A and B).** Theorems 2, 3 and 4 constrain continuous automorphism flow and the specified uniform repair channel; they do not exclude every oscillating statistic of a discrete or record-fed observer. Theorem 9 shows that the committed instrument does not select a coherent continuation.
   The committed output law alone does not select that continuation. Building
   an explicit memory extension is authorized; selecting and publishing a new
   column pin remains separate from this local research.
2. **Self-start.** Theorem 6 excludes onset from energy-invariant resources; theorem 7 gives zero ensemble mean
   from the vacuum. Theorem to prove on D′: from the blank draw state with no preparation pulse, the terminal
   return-record law at n events has a nonzero Fourier component at ω₀ = arccos of an eigenvalue of K for one fixed
   native preparation (exact on the 3,960-arc space). If the line needs a superposed preparation, onset is an owner
   decision on the preparation.
3. **ω = mc²/ℏ on the object that carries mass.** On D′ with clause B the rest phase is 0 (massless, DP §4). The
   mass angle μ is supplied. Shared gap with rows 2 and 30: the maintenance-to-energy map (theorem or decision as
   stated in row 2), then PR-22/PR-15 for units (import).

## 4. Shared premises

PR-43, PR-21 and the register clauses A, B with rows 2, 4, 29, 30. PR-22 with rows 29, 30. PR-15 with rows 2, 19,
20, 29, 30. PR-64, PR-65 with rows 2, 30. Maintenance-to-energy map with rows 2, 19, 30.

## Finite-stack scope correction and reusable recurrence clock

`Lean/Dynamics/PublicAutomorphism.lean` assumes a pointwise-continuous real-parameter
family of star automorphisms of a finite function algebra. Its identity theorem
does not exclude discrete cyclic permutations, stochastic damped oscillations,
or noncommuting instruments with classical output letters. The row-10
[directed C3 control](../shared/rows08_10_observer_boundaries/PROOF.md) is an explicit
committed-record damped oscillation. A deterministic discrete three-cycle has
return period three; it is also outside the continuous-time damping bound.
Neither control provides the row's self-starting mass-frequency identity.

`Lean/QFT/SourceRecurrenceClock.lean` supplies `meanReturn_values`,
`meanReturn_ratio_tick_free` and `returnCount_block_additive` for its specified
committed two-state source chain: means `61511/7155` and `61511/54356`. It also
proves `no_single_period` on that chain. Its `acceptedStep_no_cycle` applies to
strictly accepted mismatch-lowering repair, not every continuing observer event.
These statements were read at their actual source scope. They provide a genuine
clock precursor with a scale-free ratio and do not identify a rigid de Broglie
frequency or discharge a laboratory tick. The source chain is a distinct declared
construction; its theorem must not be silently transferred to A-cont.

The actual `model/defect_protection` and `model/defect_energy` sources already
give a protected coherent ordering charge and an exact history-dependent
memory action generated by native token transport. They must be included in
any clock construction: passive decay of a reduced memory average is not a
no-go for a record-fed controller on that larger source. Its action is an
operator depending on the ordered history, not just a count of repair costs.
`ProperTimeInternalAction.slope_selection` provides the conditional
proper-time/inertia join described in row 2; it does not by itself implement
the internal process or bind its proper-time reading to the source.

## Implemented modular flow already exists on an explicit source architecture

`model/native_modular_reference/full_flow.py` and its README implement the full
Palm modular action on all 66 native configurations, not just an active corner.
The native degree-count reference is `sigma=diag(degree/600)` and
`-log(sigma)=log(60)I+log(5/4)P_degree8`. Existing proper-code source operations
compile its nontrivial phase word in 2042 pulses, with complete padding retained.
Native current-only resets and an auxiliary pointer give complete Ramsey
readouts, including preparation timeouts. The stated complete-history comparison
bound below `3.4e-8` composes ideal-service reset and numerical word/pointer
bounds; it is not a hardware-error or duration guarantee.

This is a positive implemented flow on its declared coherent architecture.
Neither the passive channel's real spectrum nor nonuniqueness of a continuation
negates it. Its reference is stationary; a noncommuting probe and coherent
services are supplied. The clock parameter, the native attempt counter and
physical proper time are distinct until their relation is measured. The later
completed/history clock services retain their own continuing-reference and
preparation contracts. The construction supplies no electron rest mass by
itself, but it is an existing component to use in the energy/clock join.

The [row2 transfer-energy stack](../row_02_electron_rest_energy/THEOREMS.md#actual-transfer-energy-measurement-and-local-energy-clock) also contains an actual complete66-label energy detector and the exact prepared-reference identity `K=H/60+log(56)I`, not merely a formal proposal for a clock. It does not establish self-starting oscillation: the reference is stationary for its modular flow, and an admitted noncommuting preparation/readout is needed to observe phase contrast. The internal action, preparation and physical proper-time join must be checked on one source.
