import FiniteConditionalRepair
import CapFirstLaw

namespace OPH.InformationProjection

open OPH.Thermodynamics

/-!
# Global objective for the state and transition receipts

B12 receipt 1 asks for one global objective whose two marginal
optimizations are the committed state-side and transition-side
receipts.  For a finite state space with a pinned strictly positive
reference law `π` and a pinned reference kernel `K`, the objective is

`J(ρ, P) = D(ρ ‖ π) + ∑ x, ρ x * D(P x ‖ K x)`,

equal by the chain-rule theorem to the relative entropy of the joint
one-step law `ρ ⊗ P` to the reference one-step law `π ⊗ K`, so the two
committed receipts are marginals of one divergence to one reference
object.  With the transition slot pinned to `K` the objective reduces
to the state divergence and the state marginal is minimized uniquely at
`ρ = π`, which is the information-projection receipt.  At fixed
strictly positive `ρ` the transition marginal is minimized uniquely at
`P = K`, which is the conditional-resampling receipt.  The joint
minimum is `(π, K)` with value zero, and both marginal optimizers
consume the same pinned pair.  The conditional-resampling
instantiation takes `K` to be the committed `heatBath` repair kernel
built from the same `π`, so the state target and the transition target
share one reference by construction.  A three-state rational witness
realizes the objective exactly, and a control reference with different
masses moves the unique minimizer and the kernel entries, so the
objective pins its reference.

Scope boundary: this module is the finite mathematical representation
demanded by receipt 1.  A locally hash-pinned bounded run probed the
state-side and transition-side objects but did not identify the common
reference or realize this objective. Source/objective and uniform-tail
construction remains premise work under #739; physical calibration is the
separate empirical import PR-15. No source construction, physical energy,
clock, or continuum statement is made here.
-/

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω]

section StrictDivergence

/-- Strict pointwise Gibbs bound: away from `q = p` the divergence term
strictly dominates the mass difference. -/
theorem sub_lt_klTerm_of_ne {q p : ℝ} (hq : 0 ≤ q) (hp : 0 < p)
    (hne : q ≠ p) : q - p < klTerm q p := by
  rcases eq_or_lt_of_le hq with hq0 | hqpos
  · rw [← hq0, klTerm_zero]
    linarith
  have hne' : p / q ≠ 1 := by
    intro hcontra
    apply hne
    have hqne : q ≠ 0 := ne_of_gt hqpos
    field_simp at hcontra
    linarith
  have hlog : Real.log (p / q) < p / q - 1 :=
    Real.log_lt_sub_one_of_pos (div_pos hp hqpos) hne'
  have hflip : Real.log (q / p) = -Real.log (p / q) := by
    rw [← Real.log_inv]
    congr 1
    field_simp
  have hgt : 1 - p / q < Real.log (q / p) := by
    rw [hflip]
    linarith
  have hmul := mul_lt_mul_of_pos_left hgt hqpos
  rw [klTerm_of_pos (ne_of_gt hqpos)]
  have hexp : q * (1 - p / q) = q - p := by
    field_simp
  linarith [hmul, hexp]

omit [DecidableEq Ω] in
/-- Strict positivity of relative entropy away from equality, under
shared total mass and support domination. -/
theorem kl_pos_of_ne (q p : Ω → ℝ) (hq : ∀ x, 0 ≤ q x)
    (hp : ∀ x, 0 ≤ p x) (habs : ∀ x, p x = 0 → q x = 0)
    (hsum : ∑ x, q x = ∑ x, p x) (hne : q ≠ p) : 0 < kl q p := by
  obtain ⟨x₀, hx₀⟩ := Function.ne_iff.mp hne
  have hpx₀ : 0 < p x₀ := by
    rcases eq_or_lt_of_le (hp x₀) with h0 | h
    · exfalso
      have hq0 : q x₀ = 0 := habs x₀ h0.symm
      exact hx₀ (by rw [hq0, ← h0])
    · exact h
  have hstrict : ∑ x, (q x - p x) < ∑ x, klTerm (q x) (p x) := by
    apply Finset.sum_lt_sum
    · exact fun x _ => klTerm_ge (hq x) (hp x) (habs x)
    · exact ⟨x₀, Finset.mem_univ x₀,
        sub_lt_klTerm_of_ne (hq x₀) hpx₀ hx₀⟩
  have hzero : ∑ x, (q x - p x) = 0 := by
    rw [Finset.sum_sub_distrib, hsum, sub_self]
  unfold kl
  linarith

omit [DecidableEq Ω] in
/-- Vanishing relative entropy forces equality, under shared total mass
and support domination. -/
theorem eq_of_kl_eq_zero (q p : Ω → ℝ) (hq : ∀ x, 0 ≤ q x)
    (hp : ∀ x, 0 ≤ p x) (habs : ∀ x, p x = 0 → q x = 0)
    (hsum : ∑ x, q x = ∑ x, p x) (hkl : kl q p = 0) : q = p := by
  by_contra hne
  exact absurd hkl
    (ne_of_gt (kl_pos_of_ne q p hq hp habs hsum hne))

end StrictDivergence

section GlobalObjective

/-- The joint one-step law of a state and a transition kernel. -/
noncomputable def jointLaw (ρ : Ω → ℝ) (P : Ω → Ω → ℝ) :
    Ω × Ω → ℝ :=
  fun z => ρ z.1 * P z.1 z.2

/-- **Global objective.**  For the pinned reference pair `(π, K)`, the
value on the variable pair `(ρ, P)` is the state divergence to `π` plus
the `ρ`-expected row divergence of `P` to `K`.  Both summands consume
the same pinned reference pair. -/
noncomputable def objective (π : Ω → ℝ) (K : Ω → Ω → ℝ)
    (ρ : Ω → ℝ) (P : Ω → Ω → ℝ) : ℝ :=
  kl ρ π + ∑ x, ρ x * kl (P x) (K x)

omit [DecidableEq Ω] in
/-- The decomposition of the global objective into the state-side
divergence and the `ρ`-weighted transition-side divergences, with one
shared reference pair `(π, K)` on the right. -/
theorem objective_decomposition (π : Ω → ℝ) (K : Ω → Ω → ℝ)
    (ρ : Ω → ℝ) (P : Ω → Ω → ℝ) :
    objective π K ρ P = kl ρ π + ∑ x, ρ x * kl (P x) (K x) := rfl

omit [DecidableEq Ω] in
/-- **Chain rule.**  The global objective is one relative entropy on the
joint one-step space: `J(ρ, P) = D(ρ ⊗ P ‖ π ⊗ K)`.  This exhibits the
two committed receipts as marginals of a single divergence to a single
reference object. -/
theorem jointLaw_kl_eq_objective (π : Ω → ℝ) (K : Ω → Ω → ℝ)
    (ρ : Ω → ℝ) (P : Ω → Ω → ℝ)
    (hπ : ∀ x, 0 < π x)
    (hP1 : ∀ x, ∑ y, P x y = 1)
    (habs : ∀ x y, K x y = 0 → P x y = 0) :
    kl (jointLaw ρ P) (jointLaw π K) = objective π K ρ P := by
  have hterm : ∀ x y, klTerm (ρ x * P x y) (π x * K x y)
      = P x y * klTerm (ρ x) (π x) + ρ x * klTerm (P x y) (K x y) := by
    intro x y
    by_cases hρx : ρ x = 0
    · simp [hρx]
    by_cases hPy : P x y = 0
    · simp [hPy]
    have hKy : K x y ≠ 0 := fun h => hPy (habs x y h)
    have hπx : π x ≠ 0 := ne_of_gt (hπ x)
    rw [klTerm_of_pos (mul_ne_zero hρx hPy), klTerm_of_pos hρx,
      klTerm_of_pos hPy]
    have hsplit : ρ x * P x y / (π x * K x y)
        = ρ x / π x * (P x y / K x y) :=
      (div_mul_div_comm (ρ x) (π x) (P x y) (K x y)).symm
    rw [hsplit,
      Real.log_mul (div_ne_zero hρx hπx) (div_ne_zero hPy hKy)]
    ring
  have hrow : ∀ x, ∑ y, klTerm (ρ x * P x y) (π x * K x y)
      = klTerm (ρ x) (π x) + ρ x * kl (P x) (K x) := by
    intro x
    rw [Finset.sum_congr rfl fun y _ => hterm x y,
      Finset.sum_add_distrib, ← Finset.sum_mul, hP1 x, one_mul,
      ← Finset.mul_sum]
    rfl
  calc kl (jointLaw ρ P) (jointLaw π K)
      = ∑ x, ∑ y, klTerm (ρ x * P x y) (π x * K x y) := by
        unfold kl jointLaw
        rw [Fintype.sum_prod_type]
    _ = ∑ x, (klTerm (ρ x) (π x) + ρ x * kl (P x) (K x)) :=
        Finset.sum_congr rfl fun x _ => hrow x
    _ = objective π K ρ P := by
        rw [Finset.sum_add_distrib]
        rfl

/-- Nonnegativity of the global objective on normalized data. -/
theorem objective_nonneg (π : Ω → ℝ) (K : Ω → Ω → ℝ)
    (ρ : Ω → ℝ) (P : Ω → Ω → ℝ)
    (hρ0 : ∀ x, 0 ≤ ρ x) (hρ1 : ∑ x, ρ x = 1)
    (hπ : ∀ x, 0 < π x) (hπ1 : ∑ x, π x = 1)
    (hP0 : ∀ x y, 0 ≤ P x y) (hP1 : ∀ x, ∑ y, P x y = 1)
    (hK0 : ∀ x y, 0 ≤ K x y) (hK1 : ∀ x, ∑ y, K x y = 1)
    (habs : ∀ x y, K x y = 0 → P x y = 0) :
    0 ≤ objective π K ρ P := by
  unfold objective
  have hstate : 0 ≤ kl ρ π :=
    kl_nonneg ρ π hρ0 (fun x => le_of_lt (hπ x))
      (fun x h => absurd h (ne_of_gt (hπ x))) (by rw [hρ1, hπ1])
  have htrans : 0 ≤ ∑ x, ρ x * kl (P x) (K x) :=
    Finset.sum_nonneg fun x _ => mul_nonneg (hρ0 x)
      (kl_nonneg (P x) (K x) (hP0 x) (hK0 x)
        (fun y hy => habs x y hy) (by rw [hP1 x, hK1 x]))
  linarith

/-- With the transition slot pinned to the reference kernel, the global
objective reduces to the state divergence. -/
theorem objective_fixed_kernel (π : Ω → ℝ) (K : Ω → Ω → ℝ)
    (ρ : Ω → ℝ) : objective π K ρ K = kl ρ π := by
  unfold objective
  rw [Finset.sum_eq_zero fun x _ => by
    rw [kl_self, mul_zero], add_zero]

/-- The pinned reference pair sits at objective value zero. -/
theorem objective_reference_pair (π : Ω → ℝ) (K : Ω → Ω → ℝ) :
    objective π K π K = 0 := by
  rw [objective_fixed_kernel, kl_self]

/-- **State-side marginal receipt.**  At `P = K`, the reference law
minimizes the global objective over normalized states. -/
theorem state_marginal_minimizer (π : Ω → ℝ) (K : Ω → Ω → ℝ)
    (ρ : Ω → ℝ) (hρ0 : ∀ x, 0 ≤ ρ x) (hρ1 : ∑ x, ρ x = 1)
    (hπ : ∀ x, 0 < π x) (hπ1 : ∑ x, π x = 1) :
    objective π K π K ≤ objective π K ρ K := by
  rw [objective_reference_pair, objective_fixed_kernel]
  exact kl_nonneg ρ π hρ0 (fun x => le_of_lt (hπ x))
    (fun x h => absurd h (ne_of_gt (hπ x))) (by rw [hρ1, hπ1])

/-- **State-side uniqueness.**  At `P = K`, a normalized state that
attains the marginal minimum equals the reference law. -/
theorem state_marginal_unique (π : Ω → ℝ) (K : Ω → Ω → ℝ)
    (ρ : Ω → ℝ) (hρ0 : ∀ x, 0 ≤ ρ x) (hρ1 : ∑ x, ρ x = 1)
    (hπ : ∀ x, 0 < π x) (hπ1 : ∑ x, π x = 1)
    (h : objective π K ρ K = objective π K π K) : ρ = π := by
  rw [objective_reference_pair, objective_fixed_kernel] at h
  exact eq_of_kl_eq_zero ρ π hρ0 (fun x => le_of_lt (hπ x))
    (fun x hx => absurd hx (ne_of_gt (hπ x))) (by rw [hρ1, hπ1]) h

/-- **Transition-side marginal receipt.**  At fixed nonnegative `ρ`, the
reference kernel minimizes the global objective over row-stochastic
kernels dominated by the reference rows. -/
theorem transition_marginal_minimizer (π : Ω → ℝ) (K : Ω → Ω → ℝ)
    (ρ : Ω → ℝ) (P : Ω → Ω → ℝ) (hρ0 : ∀ x, 0 ≤ ρ x)
    (hP0 : ∀ x y, 0 ≤ P x y) (hP1 : ∀ x, ∑ y, P x y = 1)
    (hK0 : ∀ x y, 0 ≤ K x y) (hK1 : ∀ x, ∑ y, K x y = 1)
    (habs : ∀ x y, K x y = 0 → P x y = 0) :
    objective π K ρ K ≤ objective π K ρ P := by
  rw [objective_fixed_kernel]
  unfold objective
  have htrans : 0 ≤ ∑ x, ρ x * kl (P x) (K x) :=
    Finset.sum_nonneg fun x _ => mul_nonneg (hρ0 x)
      (kl_nonneg (P x) (K x) (hP0 x) (hK0 x)
        (fun y hy => habs x y hy) (by rw [hP1 x, hK1 x]))
  linarith

/-- **Transition-side uniqueness.**  At fixed strictly positive `ρ`, a
row-stochastic kernel that attains the marginal minimum equals the
reference kernel. -/
theorem transition_marginal_unique (π : Ω → ℝ) (K : Ω → Ω → ℝ)
    (ρ : Ω → ℝ) (P : Ω → Ω → ℝ) (hρ : ∀ x, 0 < ρ x)
    (hP0 : ∀ x y, 0 ≤ P x y) (hP1 : ∀ x, ∑ y, P x y = 1)
    (hK0 : ∀ x y, 0 ≤ K x y) (hK1 : ∀ x, ∑ y, K x y = 1)
    (habs : ∀ x y, K x y = 0 → P x y = 0)
    (h : objective π K ρ P = objective π K ρ K) : P = K := by
  rw [objective_fixed_kernel] at h
  unfold objective at h
  have hsum0 : ∑ x, ρ x * kl (P x) (K x) = 0 := by linarith
  have hterm : ∀ x ∈ Finset.univ, 0 ≤ ρ x * kl (P x) (K x) :=
    fun x _ => mul_nonneg (le_of_lt (hρ x))
      (kl_nonneg (P x) (K x) (hP0 x) (hK0 x)
        (fun y hy => habs x y hy) (by rw [hP1 x, hK1 x]))
  have hzero := (Finset.sum_eq_zero_iff_of_nonneg hterm).mp hsum0
  funext x
  have hx := hzero x (Finset.mem_univ x)
  have hklx : kl (P x) (K x) = 0 := by
    rcases mul_eq_zero.mp hx with h' | h'
    · exact absurd h' (ne_of_gt (hρ x))
    · exact h'
  exact eq_of_kl_eq_zero (P x) (K x) (hP0 x) (hK0 x)
    (fun y hy => habs x y hy) (by rw [hP1 x, hK1 x]) hklx

/-- The pinned reference pair minimizes the global objective jointly. -/
theorem joint_minimizer (π : Ω → ℝ) (K : Ω → Ω → ℝ)
    (ρ : Ω → ℝ) (P : Ω → Ω → ℝ)
    (hρ0 : ∀ x, 0 ≤ ρ x) (hρ1 : ∑ x, ρ x = 1)
    (hπ : ∀ x, 0 < π x) (hπ1 : ∑ x, π x = 1)
    (hP0 : ∀ x y, 0 ≤ P x y) (hP1 : ∀ x, ∑ y, P x y = 1)
    (hK0 : ∀ x y, 0 ≤ K x y) (hK1 : ∀ x, ∑ y, K x y = 1)
    (habs : ∀ x y, K x y = 0 → P x y = 0) :
    objective π K π K ≤ objective π K ρ P := by
  rw [objective_reference_pair]
  exact objective_nonneg π K ρ P hρ0 hρ1 hπ hπ1 hP0 hP1 hK0 hK1 habs

/-- **Joint uniqueness.**  Under strict positivity of the pinned
reference law, the global objective vanishes exactly at the pinned
reference pair: the state optimizer and the transition optimizer are
forced simultaneously and consume the same `(π, K)`. -/
theorem objective_eq_zero_iff (π : Ω → ℝ) (K : Ω → Ω → ℝ)
    (ρ : Ω → ℝ) (P : Ω → Ω → ℝ)
    (hρ0 : ∀ x, 0 ≤ ρ x) (hρ1 : ∑ x, ρ x = 1)
    (hπ : ∀ x, 0 < π x) (hπ1 : ∑ x, π x = 1)
    (hP0 : ∀ x y, 0 ≤ P x y) (hP1 : ∀ x, ∑ y, P x y = 1)
    (hK0 : ∀ x y, 0 ≤ K x y) (hK1 : ∀ x, ∑ y, K x y = 1)
    (habs : ∀ x y, K x y = 0 → P x y = 0) :
    objective π K ρ P = 0 ↔ ρ = π ∧ P = K := by
  constructor
  · intro h
    have h' := h
    unfold objective at h'
    have hstate : 0 ≤ kl ρ π :=
      kl_nonneg ρ π hρ0 (fun x => le_of_lt (hπ x))
        (fun x hx => absurd hx (ne_of_gt (hπ x))) (by rw [hρ1, hπ1])
    have htrans : 0 ≤ ∑ x, ρ x * kl (P x) (K x) :=
      Finset.sum_nonneg fun x _ => mul_nonneg (hρ0 x)
        (kl_nonneg (P x) (K x) (hP0 x) (hK0 x)
          (fun y hy => habs x y hy) (by rw [hP1 x, hK1 x]))
    have hkl0 : kl ρ π = 0 := by linarith
    have hρπ : ρ = π :=
      eq_of_kl_eq_zero ρ π hρ0 (fun x => le_of_lt (hπ x))
        (fun x hx => absurd hx (ne_of_gt (hπ x)))
        (by rw [hρ1, hπ1]) hkl0
    have hρpos : ∀ x, 0 < ρ x := fun x => by
      rw [hρπ]
      exact hπ x
    refine ⟨hρπ, transition_marginal_unique π K ρ P hρpos
      hP0 hP1 hK0 hK1 habs ?_⟩
    rw [h, objective_fixed_kernel, hkl0]
  · rintro ⟨rfl, rfl⟩
    exact objective_reference_pair ρ P

end GlobalObjective

section ConditionalResamplingReference

variable {B : Type*} [DecidableEq B]

omit [DecidableEq Ω] in
/-- The repair kernel vanishes off the fibre of the visible datum. -/
theorem heatBath_off_fiber {π : Ω → ℝ} {b : Ω → B} {x y : Ω}
    (h : b y ≠ b x) : heatBath π b x y = 0 := by
  unfold heatBath
  rw [if_neg h]

/-- For a strictly positive reference, the repair kernel vanishes
exactly off the fibre. -/
theorem heatBath_eq_zero_iff {π : Ω → ℝ} {b : Ω → B}
    (hπ : ∀ x, 0 < π x) (x y : Ω) :
    heatBath π b x y = 0 ↔ b y ≠ b x := by
  constructor
  · intro h hb
    unfold heatBath at h
    rw [if_pos hb] at h
    exact absurd h
      (ne_of_gt (div_pos (hπ y) (fiberMass_pos hπ x)))
  · exact fun h => heatBath_off_fiber h

/-- **Common-reference instantiation.**  With the pinned kernel taken to
be the conditional-resampling repair kernel built from the same pinned
law `π`, the global objective vanishes exactly at `(π, heatBath π b)`.
The state optimizer and the transition optimizer consume one reference
`π` by construction, which is the receipt's common-reference demand at
the representation level. -/
theorem conditionalResampling_objective_eq_zero_iff
    (π : Ω → ℝ) (b : Ω → B) (ρ : Ω → ℝ) (P : Ω → Ω → ℝ)
    (hπ : ∀ x, 0 < π x) (hπ1 : ∑ x, π x = 1)
    (hρ0 : ∀ x, 0 ≤ ρ x) (hρ1 : ∑ x, ρ x = 1)
    (hP0 : ∀ x y, 0 ≤ P x y) (hP1 : ∀ x, ∑ y, P x y = 1)
    (hsupp : ∀ x y, b y ≠ b x → P x y = 0) :
    objective π (heatBath π b) ρ P = 0
      ↔ ρ = π ∧ P = heatBath π b :=
  objective_eq_zero_iff π (heatBath π b) ρ P hρ0 hρ1 hπ hπ1 hP0 hP1
    (heatBath_nonneg hπ) (heatBath_row_sum hπ)
    (fun x y h => hsupp x y ((heatBath_eq_zero_iff hπ x y).mp h))

end ConditionalResamplingReference

section Witness

/-- Witness reference law on three states: masses `1/2, 1/4, 1/4`. -/
noncomputable def witnessRef : Fin 3 → ℝ :=
  fun x => if x = 0 then 1 / 2 else 1 / 4

/-- Witness visible field: state `0` is alone in its fibre; states `1`
and `2` share one fibre. -/
def witnessField : Fin 3 → Bool := fun x => decide (x = 0)

/-- Displaced witness state: masses `1/4, 3/8, 3/8`. -/
noncomputable def witnessDisplaced : Fin 3 → ℝ :=
  fun x => if x = 0 then 1 / 4 else 3 / 8

/-- Control reference law: masses `1/4, 1/2, 1/4`. -/
noncomputable def witnessControlRef : Fin 3 → ℝ :=
  fun x => if x = 0 then 1 / 4 else if x = 1 then 1 / 2 else 1 / 4

theorem witnessRef_apply_0 : witnessRef 0 = 1 / 2 := by
  unfold witnessRef
  rw [if_pos rfl]

theorem witnessRef_apply_1 : witnessRef 1 = 1 / 4 := by
  unfold witnessRef
  rw [if_neg (by decide : ¬ (1 : Fin 3) = 0)]

theorem witnessRef_apply_2 : witnessRef 2 = 1 / 4 := by
  unfold witnessRef
  rw [if_neg (by decide : ¬ (2 : Fin 3) = 0)]

theorem witnessRef_pos : ∀ x, 0 < witnessRef x := by
  intro x
  unfold witnessRef
  split <;> norm_num

theorem witnessRef_sum : ∑ x, witnessRef x = 1 := by
  rw [Fin.sum_univ_three, witnessRef_apply_0, witnessRef_apply_1,
    witnessRef_apply_2]
  norm_num

theorem witnessDisplaced_apply_0 : witnessDisplaced 0 = 1 / 4 := by
  unfold witnessDisplaced
  rw [if_pos rfl]

theorem witnessDisplaced_nonneg : ∀ x, 0 ≤ witnessDisplaced x := by
  intro x
  unfold witnessDisplaced
  split <;> norm_num

theorem witnessDisplaced_sum : ∑ x, witnessDisplaced x = 1 := by
  rw [Fin.sum_univ_three]
  unfold witnessDisplaced
  rw [if_pos rfl, if_neg (by decide : ¬ (1 : Fin 3) = 0),
    if_neg (by decide : ¬ (2 : Fin 3) = 0)]
  norm_num

theorem witnessControlRef_apply_0 : witnessControlRef 0 = 1 / 4 := by
  unfold witnessControlRef
  rw [if_pos rfl]

theorem witnessControlRef_apply_1 : witnessControlRef 1 = 1 / 2 := by
  unfold witnessControlRef
  rw [if_neg (by decide : ¬ (1 : Fin 3) = 0), if_pos rfl]

theorem witnessControlRef_apply_2 : witnessControlRef 2 = 1 / 4 := by
  unfold witnessControlRef
  rw [if_neg (by decide : ¬ (2 : Fin 3) = 0),
    if_neg (by decide : ¬ (2 : Fin 3) = 1)]

theorem witnessControlRef_pos : ∀ x, 0 < witnessControlRef x := by
  intro x
  unfold witnessControlRef
  split <;> (try split) <;> norm_num

theorem witnessControlRef_sum : ∑ x, witnessControlRef x = 1 := by
  rw [Fin.sum_univ_three, witnessControlRef_apply_0,
    witnessControlRef_apply_1, witnessControlRef_apply_2]
  norm_num

/-- Exact fibre mass of the shared fibre under the witness reference. -/
theorem witness_fiberMass_1 :
    fiberMass witnessRef witnessField 1 = 1 / 2 := by
  unfold fiberMass
  rw [Finset.sum_filter, Fin.sum_univ_three]
  rw [if_neg (by decide : ¬ witnessField 0 = witnessField 1),
    if_pos (rfl : witnessField 1 = witnessField 1),
    if_pos (by decide : witnessField 2 = witnessField 1)]
  rw [witnessRef_apply_1, witnessRef_apply_2]
  norm_num

/-- Exact rational entry of the pinned resampling kernel. -/
theorem witness_kernel_1_2 :
    heatBath witnessRef witnessField 1 2 = 1 / 2 := by
  unfold heatBath
  rw [if_pos (by decide : witnessField 2 = witnessField 1),
    witness_fiberMass_1, witnessRef_apply_2]
  norm_num

/-- Exact fibre mass of the shared fibre under the control reference. -/
theorem control_fiberMass_1 :
    fiberMass witnessControlRef witnessField 1 = 3 / 4 := by
  unfold fiberMass
  rw [Finset.sum_filter, Fin.sum_univ_three]
  rw [if_neg (by decide : ¬ witnessField 0 = witnessField 1),
    if_pos (rfl : witnessField 1 = witnessField 1),
    if_pos (by decide : witnessField 2 = witnessField 1)]
  rw [witnessControlRef_apply_1, witnessControlRef_apply_2]
  norm_num

/-- Exact rational entry of the control resampling kernel. -/
theorem control_kernel_1_2 :
    heatBath witnessControlRef witnessField 1 2 = 1 / 3 := by
  unfold heatBath
  rw [if_pos (by decide : witnessField 2 = witnessField 1),
    control_fiberMass_1, witnessControlRef_apply_2]
  norm_num

/-- Changing the pinned reference law changes the pinned resampling
kernel entrywise: `1/3` against `1/2` on the shared fibre. -/
theorem control_kernel_shift :
    heatBath witnessControlRef witnessField 1 2
      ≠ heatBath witnessRef witnessField 1 2 := by
  rw [control_kernel_1_2, witness_kernel_1_2]
  norm_num

/-- The witness objective vanishes at its own reference pair. -/
theorem witness_objective_zero :
    objective witnessRef (heatBath witnessRef witnessField)
      witnessRef (heatBath witnessRef witnessField) = 0 :=
  objective_reference_pair _ _

/-- The witness objective is strictly positive at the displaced state,
an exact strict-positivity certificate with no logarithm evaluation. -/
theorem witness_displaced_objective_pos :
    0 < objective witnessRef (heatBath witnessRef witnessField)
      witnessDisplaced (heatBath witnessRef witnessField) := by
  have hπ := witnessRef_pos
  have h0 : 0 ≤ objective witnessRef
      (heatBath witnessRef witnessField)
      witnessDisplaced (heatBath witnessRef witnessField) :=
    objective_nonneg _ _ _ _ witnessDisplaced_nonneg
      witnessDisplaced_sum hπ witnessRef_sum
      (heatBath_nonneg hπ) (heatBath_row_sum hπ)
      (heatBath_nonneg hπ) (heatBath_row_sum hπ) (fun _ _ h => h)
  rcases eq_or_lt_of_le h0 with hz | hpos
  · exfalso
    have hiff := conditionalResampling_objective_eq_zero_iff
      witnessRef witnessField witnessDisplaced
      (heatBath witnessRef witnessField) hπ witnessRef_sum
      witnessDisplaced_nonneg witnessDisplaced_sum
      (heatBath_nonneg hπ) (heatBath_row_sum hπ)
      (fun _ _ h => heatBath_off_fiber h)
    have heq := (hiff.mp hz.symm).1
    have hval := congrFun heq 0
    rw [witnessDisplaced_apply_0, witnessRef_apply_0] at hval
    norm_num at hval
  · exact hpos

/-- The control objective vanishes at the control reference pair. -/
theorem control_objective_zero :
    objective witnessControlRef
      (heatBath witnessControlRef witnessField)
      witnessControlRef
      (heatBath witnessControlRef witnessField) = 0 :=
  objective_reference_pair _ _

/-- **Negative control.**  Against the control reference, the pair that
minimizes the witness objective scores strictly positive, so the
objective genuinely pins its reference: replacing `π` by a different
law moves the unique joint minimizer. -/
theorem control_objective_pos :
    0 < objective witnessControlRef
      (heatBath witnessControlRef witnessField)
      witnessRef (heatBath witnessRef witnessField) := by
  have hπ := witnessRef_pos
  have hπ' := witnessControlRef_pos
  have h0 : 0 ≤ objective witnessControlRef
      (heatBath witnessControlRef witnessField)
      witnessRef (heatBath witnessRef witnessField) :=
    objective_nonneg _ _ _ _ (fun x => le_of_lt (hπ x))
      witnessRef_sum hπ' witnessControlRef_sum
      (heatBath_nonneg hπ) (heatBath_row_sum hπ)
      (heatBath_nonneg hπ') (heatBath_row_sum hπ')
      (fun x y h =>
        heatBath_off_fiber ((heatBath_eq_zero_iff hπ' x y).mp h))
  rcases eq_or_lt_of_le h0 with hz | hpos
  · exfalso
    have hiff := conditionalResampling_objective_eq_zero_iff
      witnessControlRef witnessField witnessRef
      (heatBath witnessRef witnessField) hπ' witnessControlRef_sum
      (fun x => le_of_lt (hπ x)) witnessRef_sum
      (heatBath_nonneg hπ) (heatBath_row_sum hπ)
      (fun _ _ h => heatBath_off_fiber h)
    have heq := (hiff.mp hz.symm).1
    have hval := congrFun heq 0
    rw [witnessRef_apply_0, witnessControlRef_apply_0] at hval
    norm_num at hval
  · exact hpos

end Witness

end OPH.InformationProjection

#print axioms OPH.InformationProjection.sub_lt_klTerm_of_ne
#print axioms OPH.InformationProjection.kl_pos_of_ne
#print axioms OPH.InformationProjection.eq_of_kl_eq_zero
#print axioms OPH.InformationProjection.objective_decomposition
#print axioms OPH.InformationProjection.jointLaw_kl_eq_objective
#print axioms OPH.InformationProjection.objective_nonneg
#print axioms OPH.InformationProjection.objective_fixed_kernel
#print axioms OPH.InformationProjection.objective_reference_pair
#print axioms OPH.InformationProjection.state_marginal_minimizer
#print axioms OPH.InformationProjection.state_marginal_unique
#print axioms OPH.InformationProjection.transition_marginal_minimizer
#print axioms OPH.InformationProjection.transition_marginal_unique
#print axioms OPH.InformationProjection.joint_minimizer
#print axioms OPH.InformationProjection.objective_eq_zero_iff
#print axioms OPH.InformationProjection.conditionalResampling_objective_eq_zero_iff
#print axioms OPH.InformationProjection.witness_kernel_1_2
#print axioms OPH.InformationProjection.control_kernel_1_2
#print axioms OPH.InformationProjection.control_kernel_shift
#print axioms OPH.InformationProjection.witness_objective_zero
#print axioms OPH.InformationProjection.witness_displaced_objective_pos
#print axioms OPH.InformationProjection.control_objective_zero
#print axioms OPH.InformationProjection.control_objective_pos
