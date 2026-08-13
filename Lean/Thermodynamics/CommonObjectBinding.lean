import MixingChainRealization
import LowTemperatureControl
import InformationProjection.GlobalObjective

namespace OPH.Thermodynamics

open OPH.InformationProjection

/-!
# Degenerate common-object compatibility witness for B12 (issue #688)

This module packages four finite representation-level ingredients relevant
to the historical B12 receipts on one structure over one realized data set.
It does not discharge the replacement source/common-object receipts of issue
#725; #688 is their bounded-negative predecessor.  The data set is the exact
restricted two-state chain of `MixingChainRealization`, extracted from the
locally pre-specified and hash-pinned bounded source run
`runs/b12_prereg_16k_20260806` of the simulator
(seed `20260806`, engine commit
`b39b78faf894894ebe573571e0902ccfaaeac32a`).  The payload with the full
extraction is `oph-physics-sim/docs/B12_MIXING_CHAIN_PAYLOAD.json`
(`scripts/extract_b12_mixing_chain.py`); the pinned transition-clock
report is `finite_repair_transition_matrix_report.json`
(sha256 `e38bb28475fb7b34127b864abd20414efc920d8ef9ddbdf5c78612874ae37f2d`)
over `observer_views.jsonl`
(sha256 `3da6f04d770b81b49a082ad1e6ecb93c611517dd85ecaccd7cc89077fcd6c0ca`).

Module chain consumed here, with the finite ingredient each carries:

* `MixingChainRealization` (receipt-3 input): the exact restricted recurrent
  chain `mixingChain` with stationary law `mixingChainStationary` and
  nonconstant protected labelling `protectedRecordLabel`;
* `InformationProjection.GlobalObjective` (receipt-1 mathematical input): the global
  objective `objective π K ρ P` whose two marginal minimizations are
  the state-side and transition-side receipts, with joint minimum
  uniquely at the pinned pair;
* `LowTemperatureControl` (receipt-4 finite control): `offMinMass`, `energyGap`,
  `refGibbs`, and `UniformGapRefinement`;
* `FiniteConditionalRepair` (through the imports): the divergence stack
  `kl`/`klTerm`, the conditional-resampling kernel `heatBath`, and the
  Gibbs family `gibbs`.

## The bound data

One state space (the restricted two-state recurrent class), one pinned
reference `bindingReference` (the real cast of the chain's exact
stationary law), one transition object `bindingChain` (the real cast of
the realized chain), one protected labelling (the committed
record-checkpoint class), and one energy `bindingEnergy` derived from
the reference by the declared convention

`E x = -log (bindingReference x)`,

with Boltzmann constant one and additive constant zero, so the
normalized Gibbs law of `E` at inverse temperature `beta = 1` equals
the pinned reference exactly (`binding_gibbs_one`), and `refGibbs` with
the uniform base law reproduces the same reference at `beta = 1`
(`binding_refGibbs_uniform_one`).

## Degenerate facts, stated as theorems

The protected labelling separates the two states, so every observation
fibre is a singleton (`binding_fiber_singleton`) and the
conditional-resampling kernel of the pinned reference over the record
fibres is the identity kernel (`bindingKernel_eq_ident`).  The kernel
component of the global objective on this data is therefore degenerate:
the joint minimum sits at the pinned pair
(`binding_objective_eq_zero_iff`), and the pinned kernel of that pair
is the identity.  The realized transition chain fails domination by
this degenerate reference (`bindingChain_not_dominated`), so the
divergence of the chain from the identity reference is undefined at
this representation level.  The well-defined orientation is the
divergence of the resampling reference from the realized chain, whose
rows are strictly positive: the objective with the realized transition
chain bound in its kernel slot computes it exactly
(`binding_kernel_divergence_value`), with strictly positive value

`(7155/61511) log (1431/1324) + (54356/61511) log (508/503)`,

the stationary-expected holding surprisal of the realized chain.  This
is the divergence of the resampling reference from the realized chain
(the chain-dominated orientation; the reverse orientation is undefined
by the domination-failure theorem), a measured separation between the
mixing dynamics and the degenerate resampling projection; the same
objective
vanishes at the realized pair itself (`binding_realized_pair_zero`).
No synthetic transition counts enter: all numerical literals trace to
retained data, while the energy, singleton-fibre kernel, and refinement
tower are declared adapter constructions.

## Claim boundary

Compatibility at the representation level over committed extracted literals.  The
stationary reference of the restricted transition chain is not identified
with the separately pre-specified state-side reference.  The realized chain is not the
identity conditional-resampling kernel forced by the singleton fibres, and
the `UniformGapRefinement` inhabitant below is a constant family rather than
a nontrivially varying regulator tower.  A common source-reference theorem,
an actual source-collar/resampling identification, and nondegenerate
refinement control therefore remain open under #732.  The physical
inverse-temperature reading and energy-clock calibration remain under #732.
-/

/-! ## The bound data -/

/-- The pinned transition-side reference: the real cast of the exact
stationary law of the realized restricted chain.  No equality with the
separately pre-specified state-side source reference is asserted. -/
noncomputable def bindingReference : Fin 2 → ℝ :=
  fun i => ((mixingChainStationary i : ℚ) : ℝ)

/-- The realized transition object: the real cast of the exact
restricted recurrent chain. -/
noncomputable def bindingChain : Fin 2 → Fin 2 → ℝ :=
  fun i j => ((mixingChain i j : ℚ) : ℝ)

/-- The conditional-resampling kernel of the pinned reference over the
fibres of the protected record labelling. -/
noncomputable def bindingKernel : Fin 2 → Fin 2 → ℝ :=
  heatBath bindingReference protectedRecordLabel

/-- The derived energy, by the declared convention `E = -log π` with
Boltzmann constant one and additive constant zero. -/
noncomputable def bindingEnergy : Fin 2 → ℝ :=
  fun i => -Real.log (bindingReference i)

/-! ## Exact values of the bound data -/

theorem bindingReference_apply_0 : bindingReference 0 = 7155 / 61511 := by
  norm_num [bindingReference, mixingChainStationary]

theorem bindingReference_apply_1 : bindingReference 1 = 54356 / 61511 := by
  norm_num [bindingReference, mixingChainStationary]

theorem bindingChain_apply_0_0 : bindingChain 0 0 = 1324 / 1431 := by
  norm_num [bindingChain, mixingChain]

theorem bindingChain_apply_0_1 : bindingChain 0 1 = 107 / 1431 := by
  norm_num [bindingChain, mixingChain]

theorem bindingChain_apply_1_0 : bindingChain 1 0 = 5 / 508 := by
  norm_num [bindingChain, mixingChain]

theorem bindingChain_apply_1_1 : bindingChain 1 1 = 503 / 508 := by
  norm_num [bindingChain, mixingChain]

/-! ## Finite facts carried over from the realized chain -/

theorem bindingReference_pos : ∀ i, 0 < bindingReference i := by
  intro i
  fin_cases i
  · show (0 : ℝ) < bindingReference 0
    rw [bindingReference_apply_0]
    norm_num
  · show (0 : ℝ) < bindingReference 1
    rw [bindingReference_apply_1]
    norm_num

theorem bindingReference_sum : (∑ i, bindingReference i) = 1 := by
  rw [Fin.sum_univ_two, bindingReference_apply_0, bindingReference_apply_1]
  norm_num

theorem bindingChain_pos : ∀ i j, 0 < bindingChain i j := by
  intro i j
  fin_cases i <;> fin_cases j
  · show (0 : ℝ) < bindingChain 0 0
    rw [bindingChain_apply_0_0]; norm_num
  · show (0 : ℝ) < bindingChain 0 1
    rw [bindingChain_apply_0_1]; norm_num
  · show (0 : ℝ) < bindingChain 1 0
    rw [bindingChain_apply_1_0]; norm_num
  · show (0 : ℝ) < bindingChain 1 1
    rw [bindingChain_apply_1_1]; norm_num

theorem bindingChain_rows : ∀ i, (∑ j, bindingChain i j) = 1 := by
  intro i
  fin_cases i
  · show (∑ j, bindingChain 0 j) = 1
    rw [Fin.sum_univ_two, bindingChain_apply_0_0, bindingChain_apply_0_1]
    norm_num
  · show (∑ j, bindingChain 1 j) = 1
    rw [Fin.sum_univ_two, bindingChain_apply_1_0, bindingChain_apply_1_1]
    norm_num

/-- The pinned transition-side reference is the exact stationary law of the
realized transition chain.  This theorem relates two fields extracted from
that chain; it does not identify the reference with the independently fixed
state-side source reference required by B12. -/
theorem bindingReference_stationary :
    ∀ j, (∑ i, bindingReference i * bindingChain i j) = bindingReference j := by
  intro j
  fin_cases j
  · show (∑ i, bindingReference i * bindingChain i 0) = bindingReference 0
    rw [Fin.sum_univ_two, bindingReference_apply_0, bindingReference_apply_1,
      bindingChain_apply_0_0, bindingChain_apply_1_0]
    norm_num
  · show (∑ i, bindingReference i * bindingChain i 1) = bindingReference 1
    rw [Fin.sum_univ_two, bindingReference_apply_0, bindingReference_apply_1,
      bindingChain_apply_0_1, bindingChain_apply_1_1]
    norm_num

/-! ## The finite compatibility structure -/

/-- **Finite compatibility record.**  One state space, one transition-side
stationary reference, one transition object, one protected labelling, and
one energy defined from that reference.  The fields record positivity,
normalization, stationarity, nonconstant labelling, and the convention
`E = -log ref`.  The type does not contain the separately pre-specified
state-side reference, a source-collar identification, or a nonconstant
refinement family. -/
structure CommonObjectBinding (Ω : Type*) [Fintype Ω] [DecidableEq Ω] where
  /-- The pinned transition-side stationary reference law. -/
  ref : Ω → ℝ
  /-- The realized transition object. -/
  chain : Ω → Ω → ℝ
  /-- The protected record labelling. -/
  label : Ω → ℕ
  /-- The derived energy. -/
  energy : Ω → ℝ
  ref_pos : ∀ x, 0 < ref x
  ref_sum : (∑ x, ref x) = 1
  chain_pos : ∀ x y, 0 < chain x y
  chain_rows : ∀ x, (∑ y, chain x y) = 1
  ref_stationary : ∀ y, (∑ x, ref x * chain x y) = ref y
  label_nonconstant : ∃ x y, label x ≠ label y
  energy_convention : ∀ x, energy x = -Real.log (ref x)

/-- The finite compatibility value carrying the run's stationary law,
transition chain, protected labelling, and derived energy simultaneously.
Its inhabitation is not a B12 common-source binding theorem. -/
noncomputable def binding : CommonObjectBinding (Fin 2) where
  ref := bindingReference
  chain := bindingChain
  label := protectedRecordLabel
  energy := bindingEnergy
  ref_pos := bindingReference_pos
  ref_sum := bindingReference_sum
  chain_pos := bindingChain_pos
  chain_rows := bindingChain_rows
  ref_stationary := bindingReference_stationary
  label_nonconstant := ⟨0, 1, protectedRecordLabel_nonconstant⟩
  energy_convention := fun _ => rfl

/-- The record fields are the named realized data, definitionally. -/
theorem binding_fields :
    binding.ref = bindingReference ∧ binding.chain = bindingChain ∧
      binding.label = protectedRecordLabel ∧
      binding.energy = bindingEnergy :=
  ⟨rfl, rfl, rfl, rfl⟩

/-! ## The degenerate resampling kernel, stated exactly -/

/-- The protected labelling separates the two states. -/
theorem protectedRecordLabel_separating :
    ∀ x y : Fin 2, protectedRecordLabel y = protectedRecordLabel x → y = x := by
  decide

/-- Every observation fibre of the protected labelling is a
singleton. -/
theorem binding_fiber_singleton :
    ∀ x : Fin 2, Finset.univ.filter
      (fun y => protectedRecordLabel y = protectedRecordLabel x) = {x} := by
  decide

/-- The fibre mass of the pinned reference at each state is the
reference mass of that state alone. -/
theorem binding_fiberMass (x : Fin 2) :
    fiberMass bindingReference protectedRecordLabel x = bindingReference x := by
  unfold fiberMass
  rw [binding_fiber_singleton x, Finset.sum_singleton]

/-- **Degenerate fact.**  On singleton fibres the conditional-resampling
kernel of the pinned reference is the identity kernel. -/
theorem bindingKernel_eq_ident :
    ∀ x y : Fin 2, bindingKernel x y = if y = x then 1 else 0 := by
  intro x y
  unfold bindingKernel heatBath
  rw [binding_fiberMass x]
  by_cases hxy : y = x
  · subst hxy
    rw [if_pos rfl, if_pos rfl,
      div_self (ne_of_gt (bindingReference_pos y))]
  · rw [if_neg (fun h => hxy (protectedRecordLabel_separating x y h)),
      if_neg hxy]

theorem bindingKernel_apply_0_0 : bindingKernel 0 0 = 1 := by
  rw [bindingKernel_eq_ident, if_pos rfl]

theorem bindingKernel_apply_0_1 : bindingKernel 0 1 = 0 := by
  rw [bindingKernel_eq_ident, if_neg (by decide : ¬ (1 : Fin 2) = 0)]

theorem bindingKernel_apply_1_0 : bindingKernel 1 0 = 0 := by
  rw [bindingKernel_eq_ident, if_neg (by decide : ¬ (0 : Fin 2) = 1)]

theorem bindingKernel_apply_1_1 : bindingKernel 1 1 = 1 := by
  rw [bindingKernel_eq_ident, if_pos rfl]

/-- **Domination failure.**  The realized transition chain is
unabsorbed by the degenerate resampling kernel: the kernel vanishes off
the diagonal where the chain is strictly positive, so the divergence of
the chain from the identity reference is undefined at this
representation level. -/
theorem bindingChain_not_dominated :
    ¬ ∀ x y, bindingKernel x y = 0 → bindingChain x y = 0 := by
  intro h
  have h01 := h 0 1 bindingKernel_apply_0_1
  rw [bindingChain_apply_0_1] at h01
  norm_num at h01

/-! ## Conditional global-objective identity on the bundled data -/

/-- The global objective over the pinned pair is nonnegative on
normalized data supported inside the record fibres. -/
theorem binding_objective_nonneg
    (ρ : Fin 2 → ℝ) (P : Fin 2 → Fin 2 → ℝ)
    (hρ0 : ∀ x, 0 ≤ ρ x) (hρ1 : (∑ x, ρ x) = 1)
    (hP0 : ∀ x y, 0 ≤ P x y) (hP1 : ∀ x, (∑ y, P x y) = 1)
    (hsupp : ∀ x y, protectedRecordLabel y ≠ protectedRecordLabel x →
      P x y = 0) :
    0 ≤ objective bindingReference bindingKernel ρ P := by
  unfold bindingKernel
  exact objective_nonneg bindingReference
    (heatBath bindingReference protectedRecordLabel) ρ P hρ0 hρ1
    bindingReference_pos bindingReference_sum hP0 hP1
    (heatBath_nonneg bindingReference_pos)
    (heatBath_row_sum bindingReference_pos)
    (fun x y h => hsupp x y
      ((heatBath_eq_zero_iff bindingReference_pos x y).mp h))

/-- **Joint-minimum identity on the bundled data.**  The global objective
over the realized reference and its record-fibre resampling kernel
vanishes exactly at the pinned pair.  On this bound object the kernel
half of the equivalence is hypothesis-forced (the record labels separate
the two states, so the support and row-sum hypotheses alone pin the
kernel); the state half carries the optimization content.  Both slots use
the same declared function inside this conditional objective, but no theorem
identifies that function with the independent state-side source reference. -/
theorem binding_objective_eq_zero_iff
    (ρ : Fin 2 → ℝ) (P : Fin 2 → Fin 2 → ℝ)
    (hρ0 : ∀ x, 0 ≤ ρ x) (hρ1 : (∑ x, ρ x) = 1)
    (hP0 : ∀ x y, 0 ≤ P x y) (hP1 : ∀ x, (∑ y, P x y) = 1)
    (hsupp : ∀ x y, protectedRecordLabel y ≠ protectedRecordLabel x →
      P x y = 0) :
    objective bindingReference bindingKernel ρ P = 0
      ↔ ρ = bindingReference ∧ P = bindingKernel := by
  unfold bindingKernel
  exact conditionalResampling_objective_eq_zero_iff bindingReference
    protectedRecordLabel ρ P bindingReference_pos bindingReference_sum
    hρ0 hρ1 hP0 hP1 hsupp

/-! ## Realized-chain separation through the objective's kernel slot -/

/-- Row `0` divergence of the resampling reference from the realized
chain. -/
theorem binding_kl_kernel_row_0 :
    kl (bindingKernel 0) (bindingChain 0) = Real.log (1431 / 1324) := by
  unfold kl
  rw [Fin.sum_univ_two, bindingKernel_apply_0_0, bindingKernel_apply_0_1,
    bindingChain_apply_0_0, klTerm_zero, add_zero,
    klTerm_of_pos one_ne_zero, one_mul, one_div_div]

/-- Row `1` divergence of the resampling reference from the realized
chain. -/
theorem binding_kl_kernel_row_1 :
    kl (bindingKernel 1) (bindingChain 1) = Real.log (508 / 503) := by
  unfold kl
  rw [Fin.sum_univ_two, bindingKernel_apply_1_0, bindingKernel_apply_1_1,
    bindingChain_apply_1_1, klTerm_zero, zero_add,
    klTerm_of_pos one_ne_zero, one_mul, one_div_div]

/-- The objective with the realized transition chain bound in its
kernel slot vanishes at the realized pair. -/
theorem binding_realized_pair_zero :
    objective bindingReference bindingChain bindingReference
      bindingChain = 0 :=
  objective_reference_pair bindingReference bindingChain

/-- **Exact divergence of the resampling reference from the realized
chain.**  With the realized transition chain bound in the objective's
kernel slot, the value at the pinned reference and the degenerate
resampling kernel is the stationary-expected holding surprisal of the
chain, computed exactly on the run literals. -/
theorem binding_kernel_divergence_value :
    objective bindingReference bindingChain bindingReference bindingKernel
      = 7155 / 61511 * Real.log (1431 / 1324)
        + 54356 / 61511 * Real.log (508 / 503) := by
  unfold objective
  rw [kl_self, zero_add, Fin.sum_univ_two, binding_kl_kernel_row_0,
    binding_kl_kernel_row_1, bindingReference_apply_0,
    bindingReference_apply_1]

/-- The exact divergence is strictly positive: the run's mixing
dynamics is measurably separated from the degenerate resampling
projection. -/
theorem binding_kernel_divergence_pos :
    0 < objective bindingReference bindingChain bindingReference
      bindingKernel := by
  rw [binding_kernel_divergence_value]
  have h1 : 0 < Real.log (1431 / 1324) := Real.log_pos (by norm_num)
  have h2 : 0 < Real.log (508 / 503) := Real.log_pos (by norm_num)
  have c1 : (0 : ℝ) < 7155 / 61511 := by norm_num
  have c2 : (0 : ℝ) < 54356 / 61511 := by norm_num
  have hm1 := mul_pos c1 h1
  have hm2 := mul_pos c2 h2
  linarith

/-! ## Single-object low-temperature facts for the derived energy -/

/-- The derived energy orders the two states strictly: the heavier
reference state carries the lower energy. -/
theorem bindingEnergy_lt : bindingEnergy 1 < bindingEnergy 0 := by
  unfold bindingEnergy
  have h : Real.log (bindingReference 0) < Real.log (bindingReference 1) := by
    apply Real.log_lt_log (bindingReference_pos 0)
    rw [bindingReference_apply_0, bindingReference_apply_1]
    norm_num
  linarith

/-- The derived energy is nonconstant, because the pinned reference is
nonuniform. -/
theorem bindingEnergy_nonconstant : bindingEnergy 1 ≠ bindingEnergy 0 :=
  ne_of_lt bindingEnergy_lt

theorem binding_gibbsWeight_one :
    ∀ x, gibbsWeight bindingEnergy 1 x = bindingReference x := by
  intro x
  unfold gibbsWeight bindingEnergy
  rw [neg_one_mul, neg_neg, Real.exp_log (bindingReference_pos x)]

theorem binding_partitionZ_one : partitionZ bindingEnergy 1 = 1 := by
  unfold partitionZ
  calc (∑ x, gibbsWeight bindingEnergy 1 x)
      = ∑ x, bindingReference x :=
        Finset.sum_congr rfl fun x _ => binding_gibbsWeight_one x
    _ = 1 := bindingReference_sum

/-- **Convention identity.**  The normalized Gibbs law of the derived
energy at inverse temperature `beta = 1` is the pinned reference,
exactly. -/
theorem binding_gibbs_one :
    ∀ x, gibbs bindingEnergy 1 x = bindingReference x := by
  intro x
  unfold gibbs
  rw [binding_gibbsWeight_one, binding_partitionZ_one, div_one]

/-- **`refGibbs` instantiation on the two-state class.**  Against the
uniform base law, the reference-relative Gibbs state of the derived
energy at `beta = 1` is the pinned reference, exactly. -/
theorem binding_refGibbs_uniform_one (x : Fin 2) :
    refGibbs (fun _ : Fin 2 => (1 : ℝ) / 2) bindingEnergy 1 x
      = bindingReference x := by
  have hw : ∀ y : Fin 2,
      refGibbsWeight (fun _ : Fin 2 => (1 : ℝ) / 2) bindingEnergy 1 y
        = 1 / 2 * bindingReference y := by
    intro y
    unfold refGibbsWeight
    have h := binding_gibbsWeight_one y
    unfold gibbsWeight at h
    rw [h]
  have hZ : refPartitionZ (fun _ : Fin 2 => (1 : ℝ) / 2) bindingEnergy 1
      = 1 / 2 := by
    unfold refPartitionZ
    rw [Finset.sum_congr rfl fun y _ => hw y, ← Finset.mul_sum,
      bindingReference_sum, mul_one]
  unfold refGibbs
  rw [hw x, hZ, mul_comm, mul_div_assoc,
    div_self (by norm_num : (1 : ℝ) / 2 ≠ 0), mul_one]

/-- The exact minimum of the derived energy sits at the heavier
reference state. -/
theorem binding_minEnergy : minEnergy bindingEnergy = bindingEnergy 1 := by
  apply le_antisymm (minEnergy_le bindingEnergy 1)
  unfold minEnergy
  apply Finset.le_inf'
  intro x _
  fin_cases x
  · show bindingEnergy 1 ≤ bindingEnergy 0
    exact le_of_lt bindingEnergy_lt
  · show bindingEnergy 1 ≤ bindingEnergy 1
    exact le_rfl

/-- The off-minimum set of the derived energy is exactly the lighter
reference state. -/
theorem binding_offMin : offMin bindingEnergy = {0} := by
  ext x
  rw [mem_offMin, binding_minEnergy, Finset.mem_singleton]
  fin_cases x
  · show bindingEnergy 0 ≠ bindingEnergy 1 ↔ (0 : Fin 2) = 0
    exact iff_of_true (ne_of_gt bindingEnergy_lt) rfl
  · show bindingEnergy 1 ≠ bindingEnergy 1 ↔ (1 : Fin 2) = 0
    exact iff_of_false (fun h => h rfl) (by decide)

theorem binding_offMin_nonempty : (offMin bindingEnergy).Nonempty := by
  rw [binding_offMin]
  exact ⟨0, Finset.mem_singleton_self 0⟩

/-- The exact energy difference between the two states is the log ratio
of the stationary masses. -/
theorem bindingEnergy_gap_value :
    bindingEnergy 0 - bindingEnergy 1 = Real.log (54356 / 7155) := by
  unfold bindingEnergy
  rw [neg_sub_neg, ← Real.log_div (ne_of_gt (bindingReference_pos 1))
    (ne_of_gt (bindingReference_pos 0)), bindingReference_apply_1,
    bindingReference_apply_0]
  congr 1
  norm_num

theorem binding_gap_pos : 0 < Real.log (54356 / 7155) :=
  Real.log_pos (by norm_num)

/-- **Exact finite gap.**  The energy gap of the derived energy is
the exact log ratio of the two stationary masses (distinct from the
chain's spectral gap recorded in the mixing-chain realization). -/
theorem binding_energyGap (h : (offMin bindingEnergy).Nonempty) :
    energyGap bindingEnergy h = Real.log (54356 / 7155) := by
  apply le_antisymm
  · have h0 : (0 : Fin 2) ∈ offMin bindingEnergy := by
      rw [binding_offMin]
      exact Finset.mem_singleton_self 0
    have hle := energyGap_le bindingEnergy h h0
    rw [binding_minEnergy] at hle
    calc energyGap bindingEnergy h
        ≤ bindingEnergy 0 - bindingEnergy 1 := hle
      _ = Real.log (54356 / 7155) := bindingEnergy_gap_value
  · unfold energyGap
    apply Finset.le_inf'
    intro x hx
    rw [binding_offMin, Finset.mem_singleton] at hx
    subst hx
    rw [binding_minEnergy]
    exact le_of_eq bindingEnergy_gap_value.symm

/-- The gap dichotomy of the derived energy, with the exact gap. -/
theorem bindingEnergy_dichotomy (x : Fin 2) :
    bindingEnergy x = minEnergy bindingEnergy ∨
      minEnergy bindingEnergy + Real.log (54356 / 7155)
        ≤ bindingEnergy x := by
  fin_cases x
  · right
    show minEnergy bindingEnergy + Real.log (54356 / 7155)
      ≤ bindingEnergy 0
    rw [binding_minEnergy]
    linarith [bindingEnergy_gap_value]
  · left
    show bindingEnergy 1 = minEnergy bindingEnergy
    exact binding_minEnergy.symm

/-- **Exact finite bound.**  The off-minimum Gibbs mass of the derived
energy carries the explicit bound with the exact cardinality and the
exact gap at every real inverse temperature. -/
theorem binding_offMinMass_le (beta : ℝ) :
    offMinMass bindingEnergy beta
      ≤ 2 * Real.exp (-beta * Real.log (54356 / 7155)) := by
  have h := offMinMass_le bindingEnergy binding_offMin_nonempty beta
  rw [binding_energyGap binding_offMin_nonempty] at h
  simpa using h

/-- The off-minimum mass of the derived energy tends to zero at low
temperature. -/
theorem binding_offMinMass_tendsto_zero :
    Filter.Tendsto (fun beta => offMinMass bindingEnergy beta)
      Filter.atTop (nhds 0) :=
  offMinMass_tendsto_zero bindingEnergy binding_offMin_nonempty

/-- **Degenerate control.**  The constant family over the bound energy
inhabits the formal refinement-uniform control structure with the exact gap
and exact cardinality.  Because every regulator carries the same `Fin 2`
object and every refinement map is the identity, this does not establish a
nontrivially varying source-derived regulator refinement. -/
noncomputable def bindingTower :
    UniformGapRefinement (fun _ : ℕ => Fin 2)
      (fun _ => bindingEnergy) where
  refineMap _ x := x
  refine_min _ _ hx := hx
  gapBound := Real.log (54356 / 7155)
  gapBound_pos := binding_gap_pos
  gap_dichotomy _ x := bindingEnergy_dichotomy x
  cardBound := 2
  card_le _ := by simp

/-- One explicit bound controls every member of the bound family at
every real inverse temperature. -/
theorem bindingTower_uniform_bound (beta : ℝ) (r : ℕ) :
    offMinMass bindingEnergy beta
      ≤ 2 * Real.exp (-beta * Real.log (54356 / 7155)) := by
  have h := bindingTower.uniform_bound beta r
  simpa using h

/-- Uniform concentration of the constant bound family.  This is a formal
control on one repeated finite object, not a continuum third-law receipt. -/
theorem bindingTower_uniform_concentration :
    ∀ eps : ℝ, 0 < eps → ∃ beta0 : ℝ, ∀ beta : ℝ, beta0 ≤ beta →
      ∀ _r : ℕ, offMinMass bindingEnergy beta < eps := by
  intro eps heps
  obtain ⟨beta0, hb⟩ := bindingTower.uniform_concentration eps heps
  exact ⟨beta0, fun beta hbeta r => hb beta hbeta r⟩

/-! ## The summary finite conjunction -/

/-- **Common-object compatibility witness (issue #688).**  Four finite
representation-level facts on one realized object, with exact values.  The
conjunction is intentionally not a closure theorem for the B12 source
receipts:

1. the global objective over the pinned reference and its record-fibre
   resampling kernel is nonnegative on admissible data and vanishes
   exactly at the pinned pair;
2. the pinned reference is the exact stationary law of the realized
   transition chain, and the protected record labelling is nonconstant; no
   equality with the separately pre-specified state-side reference is
   proved;
3. the realized restricted chain is strictly positive and row
   stochastic, its record-fibre resampling kernel is the identity (the
   degenerate fact, stated), and the objective with the realized chain
   bound in its kernel slot computes the exact strictly positive
   divergence of the resampling reference from the chain;
4. the derived energy `E = -log π` is nonconstant, reproduces the
   pinned reference as the Gibbs law at `beta = 1`, has nonempty
   off-minimum set with the exact gap `log (54356/7155)`, and carries
   the explicit off-minimum bound at every real inverse temperature on
   this single finite object. -/
theorem commonObjectBinding_receipt :
    (∀ ρ P, (∀ x, 0 ≤ ρ x) → (∑ x, ρ x) = 1 →
        (∀ x y, 0 ≤ P x y) → (∀ x, (∑ y, P x y) = 1) →
        (∀ x y, protectedRecordLabel y ≠ protectedRecordLabel x →
          P x y = 0) →
        0 ≤ objective bindingReference bindingKernel ρ P ∧
          (objective bindingReference bindingKernel ρ P = 0 ↔
            ρ = bindingReference ∧ P = bindingKernel)) ∧
    ((∀ y, (∑ i, bindingReference i * bindingChain i y)
        = bindingReference y) ∧
      protectedRecordLabel 0 ≠ protectedRecordLabel 1) ∧
    ((∀ i j, 0 < bindingChain i j) ∧
      (∀ i, (∑ j, bindingChain i j) = 1) ∧
      (∀ x y, bindingKernel x y = if y = x then 1 else 0) ∧
      objective bindingReference bindingChain bindingReference
          bindingKernel
        = 7155 / 61511 * Real.log (1431 / 1324)
          + 54356 / 61511 * Real.log (508 / 503) ∧
      0 < objective bindingReference bindingChain bindingReference
        bindingKernel) ∧
    (bindingEnergy 1 ≠ bindingEnergy 0 ∧
      (∀ x, gibbs bindingEnergy 1 x = bindingReference x) ∧
      (offMin bindingEnergy).Nonempty ∧
      (∀ h : (offMin bindingEnergy).Nonempty,
        energyGap bindingEnergy h = Real.log (54356 / 7155)) ∧
      (∀ beta : ℝ, offMinMass bindingEnergy beta
        ≤ 2 * Real.exp (-beta * Real.log (54356 / 7155)))) :=
  ⟨fun ρ P hρ0 hρ1 hP0 hP1 hsupp =>
      ⟨binding_objective_nonneg ρ P hρ0 hρ1 hP0 hP1 hsupp,
        binding_objective_eq_zero_iff ρ P hρ0 hρ1 hP0 hP1 hsupp⟩,
    ⟨bindingReference_stationary, protectedRecordLabel_nonconstant⟩,
    ⟨bindingChain_pos, bindingChain_rows, bindingKernel_eq_ident,
      binding_kernel_divergence_value, binding_kernel_divergence_pos⟩,
    ⟨bindingEnergy_nonconstant, binding_gibbs_one,
      binding_offMin_nonempty, binding_energyGap,
      binding_offMinMass_le⟩⟩

end OPH.Thermodynamics

#print axioms OPH.Thermodynamics.bindingKernel_eq_ident
#print axioms OPH.Thermodynamics.bindingChain_not_dominated
#print axioms OPH.Thermodynamics.binding_objective_eq_zero_iff
#print axioms OPH.Thermodynamics.binding_kernel_divergence_value
#print axioms OPH.Thermodynamics.binding_kernel_divergence_pos
#print axioms OPH.Thermodynamics.binding_gibbs_one
#print axioms OPH.Thermodynamics.binding_refGibbs_uniform_one
#print axioms OPH.Thermodynamics.binding_energyGap
#print axioms OPH.Thermodynamics.binding_offMinMass_le
#print axioms OPH.Thermodynamics.bindingTower_uniform_bound
#print axioms OPH.Thermodynamics.bindingTower_uniform_concentration
#print axioms OPH.Thermodynamics.commonObjectBinding_receipt
