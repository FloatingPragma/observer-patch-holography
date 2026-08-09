import MixingChainRealization
import GreenKubo

namespace OPH.Thermodynamics

/-!
# Common-reference obstruction for the current B12 artifact

The preregistered B12 run currently supplies two different finite objects.
Its state-side conditional-resampling kernel is idempotent, while the
transition-side recurrent chain has the exact nonstationary eigenvalue
`665437 / 726948`, strictly between zero and one.  This module records the
resulting no-go: an intertwiner from the transition observable space to the
state observable space must kill that nonstationary mode.

There is also a source-level arithmetic obstruction to identifying the
transition stationary law with a deterministic pushforward of the empirical
state table.  Such a pushforward has mass `k / 16384`, whereas the first
stationary mass is `7155 / 61511`, strictly between the adjacent empirical
masses `1905 / 16384` and `1906 / 16384`.

Boundary.  These statements reject a nondegenerate dynamic intertwiner and a
deterministic empirical pushforward for the *current* objects.  They do not
exclude a newly source-produced random-scan kernel, a different common
reference, or an independently justified stochastic coupling.  Merely
inventing such a coupling would not be source evidence.
-/

noncomputable section

/-! ## Linear observable actions -/

/-- The observable action of a finite real kernel as a linear map. -/
def kernelActLinear {Omega : Type*} [Fintype Omega]
    (K : Omega → Omega → ℝ) : (Omega → ℝ) →ₗ[ℝ] (Omega → ℝ) where
  toFun := kernelAct K
  map_add' := by
    intro f g
    funext x
    simp only [kernelAct, Pi.add_apply, mul_add, Finset.sum_add_distrib]
  map_smul' := by
    intro c f
    funext x
    simp only [kernelAct, Pi.smul_apply, smul_eq_mul, RingHom.id_apply]
    calc
      (∑ y, K x y * (c * f y)) = ∑ y, c * (K x y * f y) := by
        apply Finset.sum_congr rfl
        intro y _
        ring
      _ = c * ∑ y, K x y * f y := by rw [Finset.mul_sum]

/-- The exact recurrent-chain kernel, cast from the earned rational
literals to real scalars. -/
def mixingChainReal (i j : Fin 2) : ℝ := (mixingChain i j : ℚ)

/-- Backward observable action of the exact two-state mixing chain. -/
def mixingActLinear : (Fin 2 → ℝ) →ₗ[ℝ] (Fin 2 → ℝ) :=
  kernelActLinear mixingChainReal

/-- The unique nonconstant right mode of the two-state recurrent chain,
scaled to the integer vector `(54356, -7155)`. -/
def mixingMode : Fin 2 → ℝ := ![54356, -7155]

/-- Exact nonstationary eigenpair of the current recurrent chain. -/
theorem mixingMode_eigenpair :
    mixingActLinear mixingMode = (665437 / 726948 : ℝ) • mixingMode := by
  funext i
  fin_cases i <;>
    norm_num [mixingActLinear, kernelActLinear, kernelAct, mixingChainReal,
      mixingChain, mixingMode, Fin.sum_univ_two]

theorem mixingMode_eigenvalue_pos : (0 : ℝ) < 665437 / 726948 := by
  norm_num

theorem mixingMode_eigenvalue_lt_one : (665437 / 726948 : ℝ) < 1 := by
  norm_num

/-! ## The idempotent-intertwiner obstruction -/

/-- An idempotent target action cannot retain an eigenmode whose source
eigenvalue is neither zero nor one. -/
theorem idempotent_intertwiner_kills_nonprojector_mode
    {V W : Type*} [AddCommGroup V] [Module ℝ V]
    [AddCommGroup W] [Module ℝ W]
    (H : V →ₗ[ℝ] V) (P : W →ₗ[ℝ] W) (T : W →ₗ[ℝ] V)
    (hH : ∀ x, H (H x) = H x)
    (hT : ∀ x, H (T x) = T (P x))
    {v : W} {lam : ℝ} (hv : P v = lam • v)
    (hlam0 : lam ≠ 0) (hlam1 : lam ≠ 1) :
    T v = 0 := by
  have hEig : H (T v) = lam • T v := by
    rw [hT, hv, map_smul]
  have hquad : lam ^ 2 • T v = lam • T v := by
    calc
      lam ^ 2 • T v = lam • (lam • T v) := by rw [pow_two, mul_smul]
      _ = lam • H (T v) := by rw [hEig]
      _ = H (lam • T v) := by rw [map_smul]
      _ = H (H (T v)) := by rw [hEig]
      _ = H (T v) := hH _
      _ = lam • T v := hEig
  have hcoef : lam ^ 2 - lam ≠ 0 := by
    intro h
    have : lam * (lam - 1) = 0 := by
      nlinarith
    rcases mul_eq_zero.mp this with hzero | hone
    · exact hlam0 hzero
    · exact hlam1 (sub_eq_zero.mp hone)
  have hz : (lam ^ 2 - lam) • T v = 0 := by
    rw [sub_smul, hquad, sub_self]
  exact (smul_eq_zero.mp hz).resolve_left hcoef

/-- Observable action of the state-side conditional-resampling kernel. -/
def heatBathActLinear {Omega B : Type*} [Fintype Omega]
    [DecidableEq Omega] [DecidableEq B] (pi : Omega → ℝ) (b : Omega → B) :
    (Omega → ℝ) →ₗ[ℝ] (Omega → ℝ) :=
  kernelActLinear (heatBath pi b)

/-- The current heat-bath action is idempotent whenever the supplied
reference weights are strictly positive. -/
theorem heatBathActLinear_idempotent
    {Omega B : Type*} [Fintype Omega] [DecidableEq Omega] [DecidableEq B]
    (pi : Omega → ℝ) (b : Omega → B) (hpi : ∀ x, 0 < pi x) :
    ∀ f, heatBathActLinear pi b (heatBathActLinear pi b f) =
      heatBathActLinear pi b f := by
  intro f
  exact heatBath_kernelAct_idempotent hpi f

/-- Every dynamic intertwiner from the current recurrent chain to a
positive-reference heat bath annihilates the exact nonstationary mode. -/
theorem current_heatBath_intertwiner_kills_mixingMode
    {Omega B : Type*} [Fintype Omega] [DecidableEq Omega] [DecidableEq B]
    (pi : Omega → ℝ) (b : Omega → B) (hpi : ∀ x, 0 < pi x)
    (T : (Fin 2 → ℝ) →ₗ[ℝ] (Omega → ℝ))
    (hT : ∀ f, heatBathActLinear pi b (T f) = T (mixingActLinear f)) :
    T mixingMode = 0 := by
  apply idempotent_intertwiner_kills_nonprojector_mode
    (heatBathActLinear pi b) mixingActLinear T
    (heatBathActLinear_idempotent pi b hpi) hT mixingMode_eigenpair
  · norm_num
  · norm_num

/-- No nondegenerate intertwiner can bind the current idempotent state
object to the current genuinely mixing transition object. -/
theorem no_nondegenerate_current_common_object_intertwiner
    {Omega B : Type*} [Fintype Omega] [DecidableEq Omega] [DecidableEq B]
    (pi : Omega → ℝ) (b : Omega → B) (hpi : ∀ x, 0 < pi x) :
    ¬ ∃ T : (Fin 2 → ℝ) →ₗ[ℝ] (Omega → ℝ),
      (∀ f, heatBathActLinear pi b (T f) = T (mixingActLinear f)) ∧
      T mixingMode ≠ 0 := by
  rintro ⟨T, hT, hnonzero⟩
  exact hnonzero (current_heatBath_intertwiner_kills_mixingMode pi b hpi T hT)

/-! ## Empirical pushforward obstruction -/

/-- The stationary mass lies strictly between two adjacent masses available
to a deterministic pushforward of `16384` equally weighted observations. -/
theorem stationary_mass_between_empirical_atoms :
    (1905 : ℚ) / 16384 < 7155 / 61511 ∧
      7155 / 61511 < (1906 : ℚ) / 16384 := by
  norm_num

/-- No deterministic partition of the empirical `16384`-sample table has
first mass equal to the recurrent chain's exact stationary first mass. -/
theorem no_empirical_deterministic_stationary_pushforward (k : ℕ) :
    (k : ℚ) / 16384 ≠ 7155 / 61511 := by
  intro hk
  have hlo : (1905 : ℚ) < k := by
    have := stationary_mass_between_empirical_atoms.1
    rw [← hk] at this
    exact (div_lt_div_iff_of_pos_right (by norm_num : (0 : ℚ) < 16384)).mp this
  have hhi : (k : ℚ) < 1906 := by
    have := stationary_mass_between_empirical_atoms.2
    rw [← hk] at this
    exact (div_lt_div_iff_of_pos_right (by norm_num : (0 : ℚ) < 16384)).mp this
  have hloNat : 1905 < k := by exact_mod_cast hlo
  have hhiNat : k < 1906 := by exact_mod_cast hhi
  omega

/-- Summary of the exact obstruction packet: the earned mixing eigenpair,
the empirical denominator mismatch, and the universal dynamic no-go. -/
theorem current_common_reference_obstruction_summary :
    mixingActLinear mixingMode = (665437 / 726948 : ℝ) • mixingMode ∧
    (∀ k : ℕ, (k : ℚ) / 16384 ≠ 7155 / 61511) ∧
    ∀ {Omega B : Type*} [Fintype Omega] [DecidableEq Omega] [DecidableEq B]
      (pi : Omega → ℝ) (b : Omega → B), (∀ x, 0 < pi x) →
      ¬ ∃ T : (Fin 2 → ℝ) →ₗ[ℝ] (Omega → ℝ),
        (∀ f, heatBathActLinear pi b (T f) = T (mixingActLinear f)) ∧
        T mixingMode ≠ 0 := by
  exact ⟨mixingMode_eigenpair,
    no_empirical_deterministic_stationary_pushforward,
    fun pi b hpi => no_nondegenerate_current_common_object_intertwiner pi b hpi⟩

end

#print axioms mixingMode_eigenpair
#print axioms idempotent_intertwiner_kills_nonprojector_mode
#print axioms current_heatBath_intertwiner_kills_mixingMode
#print axioms no_nondegenerate_current_common_object_intertwiner
#print axioms no_empirical_deterministic_stationary_pushforward
#print axioms current_common_reference_obstruction_summary

end OPH.Thermodynamics
