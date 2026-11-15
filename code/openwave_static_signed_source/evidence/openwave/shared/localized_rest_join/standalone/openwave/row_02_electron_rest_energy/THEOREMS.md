# Row 2. Electron rest energy (mass): theorems of the repository (worker TA, 2026-11-14)

Read-only collection. Paths relative to `reverse-engineering-reality/` unless they start with `openwave-`. EOC =
`paper/executable_observer_construction.tex`, HYP = `docs/OPH_EXECUTABLE_HYPOTHESIS.md`, FR/x =
`paper/tex_fragments/x.tex`. Row test: "localized state anchored to 511 keV".

## 1. Chain

| # | Statement | Where | Kind |
| ---: | --- | --- | --- |
| 1 | Native event, continuity and conserved total load (row 1, items 1 to 3). | `Lean/Screen/RegionalContinuity.lean` `regional_continuity` | lean |
| 2 | The fair hold/exchange kernel is conditional resampling from the counting reference; it is stationary and in detailed balance with it. | `Lean/Thermodynamics/NativeBinaryThermodynamics.lean` `fairSwap_eq_heatBath`, `fairSwap_stationary`, `fairSwap_detailedBalance` | lean |
| 3 | The counting reference has constant modular energy; every calibration compatible with it assigns equal energy to all configurations. | same file, `countingReference_modularEnergy`, `counting_calibration_energy_constant` | lean |
| 4 | Mean load law μ_{k+1} − μ_k = −Lμ_k/60; every nonconstant multiplier lies in (0,1); no rescaling of the attempt count turns it into a wave law. | EOC 4095-4112 (prop:fields (iii), (iv)) | proof |
| 5 | The unread one-load channel Φ has real, strictly positive eigenvalues: event phase interval [0,0]. Φ(\|i⟩⟨i\|) differs from every unitary output by trace distance ≥ d_i/(2M). | `model/openwave_native_rest_energy/PROOF.md` "Source phase and Hamiltonian boundary" | proof |
| 6 | Every repair-only history map is a scalar times a partial injection: record laws depend on the occupation diagonal only; phase twins have identical record laws. | same PROOF.md "Native records do not observe an offdiagonal phase"; FR/NATIVE_REPAIR_RESPONSE_AND_RECORDS 106-108 | proof + exec (15,625 histories) |
| 7 | Calling a one-load coherence decay rate a rest mass gives rate ratios ≤ 30/23; with the electron anchor every such rate is below 0.000667 GeV. | same PROOF.md "Retrospective common-scale check"; FR/NATIVE_REPAIR_RESPONSE_AND_RECORDS 118-126 | proof |
| 8 | Measured on the native law: the unit moves in exactly 1/2 of its incident attempts; no charged object keeps its position. | row README; `openwave-push/mf/E_row02_energy/ROW_02/STATUS.md` (cited) | exec |
| 9 | No internal statement fixes a laboratory frequency of an internal rate without a declared tick. | `Lean/Thermodynamics/PhysicalCalibrationImport.lean` `labFrequency_not_forced` | lean |
| 10 | Installed kinetic/analyzer source, h = −τ₃⊗σ·D + mβ, m = e₁: a compact packet ψ = P₊f/√w has 0 ≤ ⟨E⟩_ψ − m ≤ a₂/(2m(1−η)) and outside-cube probability ≤ η/(1−η); with L = 2^30 the excess is ≤ 3.014e-9 relative. | `model/openwave_localized_rest_packet/PROOF.md` "Exact finite moments and rest-energy bounds" | conditional |
| 11 | Same packet: \|C(t) − e^{−imt}\| ≤ \|t\|d; phase frequency within d/(1 − \|t\|d) of m; no passive onset (dephasing removes contrast). | same PROOF.md "Autonomous phase records" | conditional |
| 12 | With Y_up = Y_down = Y_lepton = diag(e_k) declared: H − E_ground ≥ m\|Q\| on the whole Fock space, attained by a charged-lepton rest mode. | `model/openwave_full_higgs_floor/PROOF.md` eq. (3) | conditional |
| 13 | The 511 keV value enters as one anchor m = 0.00051099895069 GeV; "it is not a prediction of the electron's absolute mass". | localized_rest_packet PROOF.md lines 109-114; EOC 1119-1137 (held-source member, g² = 12πα_input) | premise (PR-15) |

## 2. Premises

| Id | Text | State |
| --- | --- | --- |
| PR-05, PR-20 | uniform seams and fair coin | open (declared); PR-20 forced by `a2a3_directed_seam_weight_eq_one_sixtieth` given open projection data |
| PR-15 | clock and energy calibration anchors: laboratory clock and energy enter as declared standards | open (empirical import by disposition; the row admits one anchor) |
| unregistered: `independent_physical_scale_receipt` (`claims/assumption_dictionary.md` line 113) | two independently generated source readouts fix the positive rescaling of the internal equation | open; this is the source-side energy normalization the register separates from PR-15 |
| unregistered: maintenance-to-energy map (HYP 113-120) | "The maintenance-to-energy map must be derived, including kinetic normalization and any required regulator/reference subtraction" | open |
| PR-21, PR-43 (as used by candidate D′) | reuse and blind redraw of the unread choice register; a supplied automorphism flow | open; PR-43 obstructed on the committed record algebra (`Lean/Dynamics/PublicAutomorphism.lean` `ContinuousPublicStarFlow.toAut_eq_refl`) |
| PR-64, PR-65 | coherent phase instrument and source-produced phase preparation (item 11) | realized only on an installed source; PR-65 obstructed for committed C/E repair |
| unregistered | kinetic/analyzer source, positive-band filter, filled reference (items 10, 11) | realized only on an installed source (`openwave_localized_rest_packet`) |
| unregistered | equal family Yukawa matrices (item 12) | realized only on an installed source (declared extension, full_higgs_floor PROOF.md lines 13-16) |

## 3. The gap

Theorems 2 to 7 show that the counting reference on the stated fixed-load configuration domain has constant modular energy (3), the specified repair-only instruments are blind to source off-diagonal phase (5, 6), and the tested coherence-decay rates fail the charged-lepton mass comparison (7). This does not exclude nonconstant record-cost observables on other native record algebras. The smallest open set:

1. **Maintenance-to-energy map** (unregistered). Theorem to prove, or decision. On candidate D′ (12-port
   carrier, mass redraw Π + z(1 − Π), z = −e^{2iμ}): the zero-momentum eigenphase ω₀(μ) of the one-load walk and the
   all-read movement excess ε (defect minus reference movement bits per own attempt, exactly 1/2 at the fair coin)
   satisfy 1 − cos ω₀(μ) = ε for one μ fixed by the native law alone. If μ is not fixed by any native quantity,
   the rule "energy = record rate" is an owner decision of the candidate (F0).
2. **A localized charged object** without holding apparatus. Shared with row 4; on A the unit diffuses (theorem 4).
3. **Energy units and the common clock.** The static row explicitly admits a
   511 keV anchor. An independent prediction of that absolute scale is not an
   additional row-2 requirement. The source energy observable must still be
   specified and measured before anchoring it; the anchor does not derive that
   observable. Its clock must be shared when carrying the result into row 3.
4. **Cross-row propagation join.** Measuring the same object's dispersion and
   inertial response is required for claims importing rows 3 or 30. It is not
   an extra clause in row 2's literal static test.

## 4. Shared premises

PR-05, PR-20 with all rows. PR-15 with rows 3, 19, 20, 29, 30. `independent_physical_scale_receipt` and the
maintenance-to-energy map with rows 3, 19, 30. PR-21, PR-43 with rows 3, 4, 29, 30. PR-64, PR-65 with rows 3, 30.
Localization with rows 4 and 19.

## Finite-stack scope correction and reusable record energy

The actual Lean file `Lean/Thermodynamics/NativeRecordFourLaws.lean` proves
`record_thermal`, `record_four_laws`, `record_energy_gap` and `record_thermometer`
for the stopped prefix records `0,10,11`, with counted weights `(1/2,1/4,1/4)` and
record costs `(log 2,2 log 2,2 log 2)`. This is a nonconstant native record-energy
construction in record-cost units. It does not contradict the counting-reference
result on a different algebra and must be reused rather than ruled out by it.
Physical energy normalization, common clock and its attachment to a localized
charged object remain the stated row obligations. No laboratory mass follows
from the prefix code alone. The new [load-only invariant theorem](../shared/rows08_10_observer_boundaries/PROOF.md) also explicitly excludes active memory from its scope.

A useful corollary of the load-only invariant theorem is that an additive conserved
energy depending solely on total signed load must have `E(Q)=a Q` after fixing the
vacuum zero. If it is also invariant under charge conjugation, `E(Q)=E(−Q)`, then
`a=0`. The premises matter: this rules out identifying a nonzero equal particle/
antiparticle rest energy with an additive load-only invariant. It directs the
construction toward active observer memory, interaction records and the already
proved record-energy stack; it is not a no-go for those larger architectures.

## Same-source action and the conditional inertia theorem

`model/defect_energy/README.md` derives the exact action of the protected
orientation memory on the binary C/E source: the source-dependent ordered-token
parity gives `X^p=exp(-i pi p P_minus)`. This is a complete-branch operator
identity, with its original fair-coin amplitude retained. Two seven-move
histories with identical multisets of raw event labels produce different
memory actions; therefore a scalar additive cost of those labels alone is
insufficient. This positive action construction is not a rest-energy or
physical-duration measurement.

The actual `Lean/Geometry/InternalEnergyInertia.lean` binds modular energy to
the cap first law and proves inertia `m+E` for its declared composite action.
`Lean/Geometry/ProperTimeInternalAction.lean` goes further: midpoint refinement
invariance eliminates the per-step additive term; agreement with the declared
proper-time principle selects slope one. Its `slope_selection` and
`lengthStationary_iff_eom` then give the corresponding equation of motion with
coefficient `m+E`. These are usable finite conditional results. Their
proper-time action, Lorentz attachment and ledger-to-physical-energy premises
must be bound to the same executable observer; mathematical nonuniqueness is
not a reason to forbid constructing that observer.

## Actual transfer-energy measurement and local energy clock

The positive stack extends beyond a formal logarithm. The actual
`model/native_energy_measurement/PROOF.md` constructs a complete destructive
measurement of `H=-60 log P` on the same fair native binary66 carrier, where
`P=I-G/60` and `G=sum_e(I-S_e)` are obtained from literal endpoint branches.
Its66 sequential source filters retain every energy label, processor failure,
padding output and final residual. Exact rational residuals bound the proposed
energy frame; the stored-map complete half-diamond allowance is less than
`1/12000`. The3,853,872 source pulses are a bounded composition of declared
proper-code services, not an autonomous repair operation or a hardware result.
Its energy observable is fixed independently of detector runtime.

The existing source preparation `rho=P/56` has exactly
`K=-log rho=H/60+log(56)I`. Thus its modular action equals the controlled
transfer action under the stated parameter conversion. This is a positive
same-reference energy/clock identity, with actual source services behind its
readout. It does not identify repair attempts, service pulse durations and
physical proper time. It is distinct from the diagonal Palm reference with
degrees8 and10; their respective gaps cannot be merged without a map.

`model/native_local_energy_clock/PROOF.md` then establishes a local comparison:
on a fixed-density binary federation, the counted transfer action approaches
a common rescaled exchange action in its specified reference-weighted norm,
with an explicit `O(M^-1/2)` bound for fixed local observables and complete
fixed finite detector menus. Its seam and site energy currents obey exact
continuity under the controlled exchange. This supplies a legitimate
source-derived energy/current construction. The growing reference preparation,
physical time and attachment to the same localized charged object remain
separate, explicitly stated obligations. A uniqueness theorem for all possible
clocks is not a prerequisite to testing this declared construction.

These source proofs were inspected directly in the current pass. No new
replay, Lean build or laboratory energy identification is asserted here.
