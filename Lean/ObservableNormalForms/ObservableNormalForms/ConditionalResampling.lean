import ObservableNormalForms.Stochastic

/-!
# Finite conditional resampling

This module gives an elementary finite weighted model of conditional
resampling along the fibers of an observation map.  The operator is developed
entirely as finite-sum algebra: it fixes exactly the functions that are
constant on observation fibers, is idempotent and self-adjoint for the
weighted inner product, and satisfies the corresponding Pythagorean energy
identity and squared-`L2` contraction.
-/

namespace ObservableNormalForms

universe u v

/-- A finite state space, a strictly positive weight at every state, and an
observation map.  The observation type itself need not be finite: every
occupied observation fiber is finite because the state space is finite. -/
structure FiniteWeightedObservation
    (S : Type u) (O : Type v) [Fintype S] [DecidableEq O] where
  weight : S → ℝ
  weight_pos : ∀ x : S, 0 < weight x
  observe : S → O

namespace FiniteWeightedObservation

variable {S : Type u} {O : Type v} [Fintype S] [DecidableEq O]

/-- The finite fiber through `x`. -/
def fiber (M : FiniteWeightedObservation S O) (x : S) : Finset S :=
  Finset.univ.filter (fun y => M.observe y = M.observe x)

@[simp]
theorem mem_fiber_iff (M : FiniteWeightedObservation S O) (x y : S) :
    y ∈ M.fiber x ↔ M.observe y = M.observe x := by
  simp [fiber]

@[simp]
theorem self_mem_fiber (M : FiniteWeightedObservation S O) (x : S) :
    x ∈ M.fiber x := by
  simp

theorem fiber_eq_of_observe_eq (M : FiniteWeightedObservation S O)
    {x y : S} (hxy : M.observe x = M.observe y) :
    M.fiber x = M.fiber y := by
  ext z
  simp only [mem_fiber_iff]
  constructor <;> intro hz
  · exact hz.trans hxy
  · exact hz.trans hxy.symm

/-- Total weight of the occupied observation fiber through `x`. -/
def fiberMass (M : FiniteWeightedObservation S O) (x : S) : ℝ :=
  ∑ y ∈ M.fiber x, M.weight y

theorem fiberMass_pos (M : FiniteWeightedObservation S O) (x : S) :
    0 < M.fiberMass x := by
  unfold fiberMass
  exact Finset.sum_pos' (fun y _ => (M.weight_pos y).le)
    ⟨x, M.self_mem_fiber x, M.weight_pos x⟩

theorem fiberMass_ne_zero (M : FiniteWeightedObservation S O) (x : S) :
    M.fiberMass x ≠ 0 := (M.fiberMass_pos x).ne'

theorem fiberMass_eq_of_observe_eq (M : FiniteWeightedObservation S O)
    {x y : S} (hxy : M.observe x = M.observe y) :
    M.fiberMass x = M.fiberMass y := by
  rw [fiberMass, fiberMass, M.fiber_eq_of_observe_eq hxy]

/-- The conditional-resampling kernel: from `x`, choose a state in the same
observation fiber with probability proportional to its weight. -/
noncomputable def transition (M : FiniteWeightedObservation S O) (x y : S) : ℝ :=
  if M.observe y = M.observe x then M.weight y / M.fiberMass x else 0

theorem transition_nonneg (M : FiniteWeightedObservation S O) (x y : S) :
    0 ≤ M.transition x y := by
  by_cases hxy : M.observe y = M.observe x
  · simp only [transition, hxy, if_true]
    exact div_nonneg (M.weight_pos y).le (M.fiberMass_pos x).le
  · simp [transition, hxy]

theorem transition_sum_one (M : FiniteWeightedObservation S O) (x : S) :
    ∑ y : S, M.transition x y = 1 := by
  rw [show (∑ y : S, M.transition x y) =
      (∑ y ∈ M.fiber x, M.weight y) / M.fiberMass x by
    simp only [transition, fiber, Finset.sum_filter]
    rw [Finset.sum_div]
    apply Finset.sum_congr rfl
    intro y _
    split_ifs <;> simp]
  rw [show (∑ y ∈ M.fiber x, M.weight y) = M.fiberMass x by rfl]
  exact div_self (M.fiberMass_ne_zero x)

/-- The conditional-resampling transition matrix as a finite Markov kernel. -/
noncomputable def conditionalResamplingKernel
    (M : FiniteWeightedObservation S O) : FiniteMarkovKernel S where
  probability := M.transition
  probability_nonneg := M.transition_nonneg
  probability_sum_one := M.transition_sum_one

/-- R1: transitions do not leave an observation fiber. -/
def FiberSupported (M : FiniteWeightedObservation S O)
    (K : FiniteMarkovKernel S) : Prop :=
  ∀ ⦃x y : S⦄, M.observe x ≠ M.observe y → K.probability x y = 0

/-- R2: starting states in the same observation fiber have identical rows. -/
def FiberRowConstant (M : FiniteWeightedObservation S O)
    (K : FiniteMarkovKernel S) : Prop :=
  ∀ ⦃x x' : S⦄, M.observe x = M.observe x' →
    K.probability x = K.probability x'

/-- R3: detailed balance with respect to the prescribed positive weights. -/
def WeightedDetailedBalance (M : FiniteWeightedObservation S O)
    (K : FiniteMarkovKernel S) : Prop :=
  ∀ x y : S,
    M.weight x * K.probability x y = M.weight y * K.probability y x

theorem conditionalResamplingKernel_fiberSupported
    (M : FiniteWeightedObservation S O) :
    M.FiberSupported M.conditionalResamplingKernel := by
  intro x y hxy
  have hyx : M.observe y ≠ M.observe x := fun hyx => hxy hyx.symm
  simp [conditionalResamplingKernel, transition, hyx]

theorem conditionalResamplingKernel_fiberRowConstant
    (M : FiniteWeightedObservation S O) :
    M.FiberRowConstant M.conditionalResamplingKernel := by
  intro x x' hxx'
  funext y
  show M.transition x y = M.transition x' y
  by_cases hyx : M.observe y = M.observe x
  · have hyx' : M.observe y = M.observe x' := hyx.trans hxx'
    rw [transition, transition, if_pos hyx, if_pos hyx',
      M.fiberMass_eq_of_observe_eq hxx']
  · have hyx' : M.observe y ≠ M.observe x' :=
      fun hyx' => hyx (hyx'.trans hxx'.symm)
    rw [transition, transition, if_neg hyx, if_neg hyx']

theorem conditionalResamplingKernel_weightedDetailedBalance
    (M : FiniteWeightedObservation S O) :
    M.WeightedDetailedBalance M.conditionalResamplingKernel := by
  intro x y
  show M.weight x * M.transition x y = M.weight y * M.transition y x
  by_cases hxy : M.observe x = M.observe y
  · rw [transition, transition, if_pos hxy.symm, if_pos hxy,
      M.fiberMass_eq_of_observe_eq hxy]
    ring
  · have hyx : M.observe y ≠ M.observe x := fun hyx => hxy hyx.symm
    rw [transition, transition, if_neg hyx, if_neg hxy]
    ring

/-- Conditional resampling of a real-valued function. -/
noncomputable def resample (M : FiniteWeightedObservation S O)
    (f : S → ℝ) (x : S) : ℝ :=
  ∑ y : S, M.transition x y * f y

/-- Explicit fiber-average formula for conditional resampling. -/
theorem resample_eq_fiber_average (M : FiniteWeightedObservation S O)
    (f : S → ℝ) (x : S) :
    M.resample f x =
      (∑ y ∈ M.fiber x, M.weight y * f y) / M.fiberMass x := by
  simp only [resample, transition, fiber, Finset.sum_filter]
  rw [Finset.sum_div]
  apply Finset.sum_congr rfl
  intro y _
  split_ifs <;> ring

theorem resample_nonneg (M : FiniteWeightedObservation S O)
    {f : S → ℝ} (hf : ∀ x : S, 0 ≤ f x) (x : S) :
    0 ≤ M.resample f x := by
  unfold resample
  apply Finset.sum_nonneg
  intro y _
  exact mul_nonneg (M.transition_nonneg x y) (hf y)

theorem resample_const (M : FiniteWeightedObservation S O)
    (c : ℝ) : M.resample (fun _ => c) = fun _ => c := by
  funext x
  unfold resample
  rw [← Finset.sum_mul, M.transition_sum_one, one_mul]

theorem resample_add (M : FiniteWeightedObservation S O)
    (f g : S → ℝ) :
    M.resample (fun x => f x + g x) =
      fun x => M.resample f x + M.resample g x := by
  funext x
  simp only [resample]
  rw [← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro y _
  ring

theorem resample_smul (M : FiniteWeightedObservation S O)
    (a : ℝ) (f : S → ℝ) :
    M.resample (fun x => a * f x) = fun x => a * M.resample f x := by
  funext x
  simp only [resample]
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro y _
  ring

/-- A function is observation-measurable when it is constant on each occupied
observation fiber.  This is the finite, pointwise form of measurability with
respect to the observation map. -/
def ObservationMeasurable (M : FiniteWeightedObservation S O)
    (f : S → ℝ) : Prop :=
  ∀ ⦃x y : S⦄, M.observe x = M.observe y → f x = f y

/-- Fiberwise constancy is equivalent to factoring through the observation
map, even when the observation type contains values not attained by a state. -/
theorem observationMeasurable_iff_exists_factor
    (M : FiniteWeightedObservation S O) (f : S → ℝ) :
    M.ObservationMeasurable f ↔
      ∃ g : O → ℝ, f = g ∘ M.observe := by
  classical
  constructor
  · intro hf
    let g : O → ℝ := fun o =>
      if h : ∃ x : S, M.observe x = o then f (Classical.choose h) else 0
    refine ⟨g, ?_⟩
    funext x
    let hx : ∃ y : S, M.observe y = M.observe x := ⟨x, rfl⟩
    have hchoice : M.observe (Classical.choose hx) = M.observe x :=
      Classical.choose_spec hx
    change f x = g (M.observe x)
    rw [show g (M.observe x) = f (Classical.choose hx) by
      simp only [g]
      rw [dif_pos hx]]
    exact (hf hchoice).symm
  · rintro ⟨g, rfl⟩ x y hxy
    simp only [Function.comp_apply]
    rw [hxy]

theorem resample_observationMeasurable
    (M : FiniteWeightedObservation S O) (f : S → ℝ) :
    M.ObservationMeasurable (M.resample f) := by
  intro x y hxy
  rw [resample_eq_fiber_average, resample_eq_fiber_average,
    M.fiber_eq_of_observe_eq hxy, M.fiberMass_eq_of_observe_eq hxy]

theorem resample_eq_self_of_observationMeasurable
    (M : FiniteWeightedObservation S O) {f : S → ℝ}
    (hf : M.ObservationMeasurable f) :
    M.resample f = f := by
  funext x
  rw [M.resample_eq_fiber_average]
  have hsum :
      (∑ y ∈ M.fiber x, M.weight y * f y) = M.fiberMass x * f x := by
    rw [fiberMass, Finset.sum_mul]
    apply Finset.sum_congr rfl
    intro y hy
    have hyx : M.observe y = M.observe x := (M.mem_fiber_iff x y).mp hy
    rw [hf hyx]
  rw [hsum]
  exact mul_div_cancel_left₀ (f x) (M.fiberMass_ne_zero x)

/-- Conditional resampling fixes exactly the observation-measurable
functions. -/
theorem resample_eq_self_iff_observationMeasurable
    (M : FiniteWeightedObservation S O) (f : S → ℝ) :
    M.resample f = f ↔ M.ObservationMeasurable f := by
  constructor
  · intro hfix
    rw [← hfix]
    exact M.resample_observationMeasurable f
  · exact M.resample_eq_self_of_observationMeasurable

/-- Conditional resampling is idempotent. -/
theorem resample_idempotent
    (M : FiniteWeightedObservation S O) (f : S → ℝ) :
    M.resample (M.resample f) = M.resample f :=
  M.resample_eq_self_of_observationMeasurable
    (M.resample_observationMeasurable f)

/-- The finite weighted inner product used by the resampling operator. -/
def weightedInner (M : FiniteWeightedObservation S O)
    (f g : S → ℝ) : ℝ :=
  ∑ x : S, M.weight x * f x * g x

/-- Symmetric pair weight underlying self-adjointness. -/
noncomputable def pairWeight
    (M : FiniteWeightedObservation S O) (x y : S) : ℝ :=
  M.weight x * M.transition x y

theorem pairWeight_comm (M : FiniteWeightedObservation S O) (x y : S) :
    M.pairWeight x y = M.pairWeight y x := by
  by_cases hxy : M.observe x = M.observe y
  · rw [pairWeight, pairWeight, transition, transition,
      if_pos hxy.symm, if_pos hxy]
    rw [M.fiberMass_eq_of_observe_eq hxy]
    ring
  · have hyx : M.observe y ≠ M.observe x := fun hyx => hxy hyx.symm
    rw [pairWeight, pairWeight, transition, transition,
      if_neg hxy, if_neg hyx]
    ring

theorem finiteMarkovKernel_ext
    {K L : FiniteMarkovKernel S}
    (hProbability : K.probability = L.probability) : K = L := by
  cases K
  cases L
  cases hProbability
  rfl

/-- R1--R3 recognize the conditional-resampling kernel.  Row stochasticity
and nonnegativity are carried by `FiniteMarkovKernel`; strict positivity of
the prescribed weights makes normalization on every occupied fiber unique. -/
theorem kernel_eq_conditionalResamplingKernel_of_recognition
    (M : FiniteWeightedObservation S O) (K : FiniteMarkovKernel S)
    (hSupport : M.FiberSupported K)
    (hRows : M.FiberRowConstant K)
    (hBalance : M.WeightedDetailedBalance K) :
    K = M.conditionalResamplingKernel := by
  apply finiteMarkovKernel_ext
  funext x y
  show K.probability x y = M.transition x y
  by_cases hyx : M.observe y = M.observe x
  · have hentryOf : ∀ ⦃z : S⦄, M.observe z = M.observe x →
        K.probability x z =
          M.weight z * (K.probability x x / M.weight x) := by
      intro z hzx
      have hrow_zx : K.probability z x = K.probability x x := by
        have hrow : K.probability z = K.probability x := hRows hzx
        exact congrFun hrow x
      have hbal := hBalance x z
      rw [hrow_zx] at hbal
      calc
        K.probability x z =
            (M.weight x * K.probability x z) / M.weight x :=
          (mul_div_cancel_left₀ (K.probability x z)
            (M.weight_pos x).ne').symm
        _ = (M.weight z * K.probability x x) / M.weight x := by rw [hbal]
        _ = M.weight z * (K.probability x x / M.weight x) := by ring
    have hentry := hentryOf hyx
    have hsumFiber :
        (∑ z ∈ M.fiber x, K.probability x z) = 1 := by
      calc
        (∑ z ∈ M.fiber x, K.probability x z) =
            ∑ z : S,
              if M.observe z = M.observe x then K.probability x z else 0 := by
          simp only [fiber, Finset.sum_filter]
        _ = ∑ z : S, K.probability x z := by
          apply Finset.sum_congr rfl
          intro z _
          by_cases hzx : M.observe z = M.observe x
          · rw [if_pos hzx]
          · rw [if_neg hzx]
            exact (hSupport (fun hxz => hzx hxz.symm)).symm
        _ = 1 := K.probability_sum_one x
    have hnormal :
        M.fiberMass x * (K.probability x x / M.weight x) = 1 := by
      calc
        M.fiberMass x * (K.probability x x / M.weight x) =
            ∑ z ∈ M.fiber x,
              M.weight z * (K.probability x x / M.weight x) := by
          rw [fiberMass, Finset.sum_mul]
        _ = ∑ z ∈ M.fiber x, K.probability x z := by
          apply Finset.sum_congr rfl
          intro z hz
          have hzx : M.observe z = M.observe x := (M.mem_fiber_iff x z).mp hz
          exact (hentryOf hzx).symm
        _ = 1 := hsumFiber
    have hconstant :
        K.probability x x / M.weight x = 1 / M.fiberMass x := by
      apply (eq_div_iff (M.fiberMass_ne_zero x)).2
      simpa [mul_comm] using hnormal
    rw [transition, if_pos hyx, hentry, hconstant]
    ring
  · have hxy : M.observe x ≠ M.observe y := fun hxy => hyx hxy.symm
    rw [hSupport hxy, transition, if_neg hyx]

/-- Exact matrix recognition: a finite Markov kernel is the conditional
fiber-resampling kernel if and only if it satisfies R1 fiber support, R2
identical rows inside fibers, and R3 weighted detailed balance. -/
theorem kernel_eq_conditionalResamplingKernel_iff_recognition
    (M : FiniteWeightedObservation S O) (K : FiniteMarkovKernel S) :
    K = M.conditionalResamplingKernel ↔
      M.FiberSupported K ∧
      M.FiberRowConstant K ∧
      M.WeightedDetailedBalance K := by
  constructor
  · intro hK
    subst K
    exact ⟨M.conditionalResamplingKernel_fiberSupported,
      M.conditionalResamplingKernel_fiberRowConstant,
      M.conditionalResamplingKernel_weightedDetailedBalance⟩
  · rintro ⟨hSupport, hRows, hBalance⟩
    exact M.kernel_eq_conditionalResamplingKernel_of_recognition K
      hSupport hRows hBalance

theorem weightedInner_resample_expand
    (M : FiniteWeightedObservation S O) (f g : S → ℝ) :
    M.weightedInner (M.resample f) g =
      ∑ x : S, ∑ y : S, M.pairWeight x y * f y * g x := by
  simp only [weightedInner, resample, pairWeight]
  apply Finset.sum_congr rfl
  intro x _
  rw [Finset.mul_sum, Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro y _
  ring

/-- Conditional resampling is self-adjoint for the weighted inner product. -/
theorem resample_weighted_self_adjoint
    (M : FiniteWeightedObservation S O) (f g : S → ℝ) :
    M.weightedInner (M.resample f) g =
      M.weightedInner f (M.resample g) := by
  rw [M.weightedInner_resample_expand]
  calc
    (∑ x : S, ∑ y : S, M.pairWeight x y * f y * g x)
        = ∑ y : S, ∑ x : S, M.pairWeight x y * f y * g x :=
          Finset.sum_comm
    _ = ∑ y : S, ∑ x : S, M.pairWeight y x * f y * g x := by
      apply Finset.sum_congr rfl
      intro y _
      apply Finset.sum_congr rfl
      intro x _
      rw [M.pairWeight_comm]
    _ = M.weightedInner f (M.resample g) := by
      simp only [weightedInner, resample, pairWeight]
      apply Finset.sum_congr rfl
      intro y _
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro x _
      ring

/-- Symmetry of the finite weighted inner product. -/
theorem weightedInner_comm
    (M : FiniteWeightedObservation S O) (f g : S → ℝ) :
    M.weightedInner f g = M.weightedInner g f := by
  unfold weightedInner
  apply Finset.sum_congr rfl
  intro x _
  ring

/-- Linearity in the first argument, recorded explicitly for the energy
calculation below. -/
theorem weightedInner_sub_left
    (M : FiniteWeightedObservation S O) (f g h : S → ℝ) :
    M.weightedInner (fun x => f x - g x) h =
      M.weightedInner f h - M.weightedInner g h := by
  unfold weightedInner
  rw [← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro x _
  ring

/-- Linearity in the second argument. -/
theorem weightedInner_sub_right
    (M : FiniteWeightedObservation S O) (f g h : S → ℝ) :
    M.weightedInner f (fun x => g x - h x) =
      M.weightedInner f g - M.weightedInner f h := by
  unfold weightedInner
  rw [← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro x _
  ring

/-- The resampling residual is orthogonal to every
observation-measurable function. -/
theorem residual_orthogonal_observationMeasurable
    (M : FiniteWeightedObservation S O) (f h : S → ℝ)
    (hh : M.ObservationMeasurable h) :
    M.weightedInner (fun x => f x - M.resample f x) h = 0 := by
  rw [M.weightedInner_sub_left]
  rw [← M.resample_eq_self_of_observationMeasurable hh]
  rw [M.resample_weighted_self_adjoint]
  rw [M.resample_idempotent]
  exact sub_self _

/-- Squared finite weighted `L2` energy. -/
def weightedEnergy (M : FiniteWeightedObservation S O) (f : S → ℝ) : ℝ :=
  M.weightedInner f f

theorem weightedEnergy_nonneg
    (M : FiniteWeightedObservation S O) (f : S → ℝ) :
    0 ≤ M.weightedEnergy f := by
  unfold weightedEnergy weightedInner
  apply Finset.sum_nonneg
  intro x _
  have hw := (M.weight_pos x).le
  nlinarith [sq_nonneg (f x)]

/-- Exact Pythagorean decomposition into conditional average and residual. -/
theorem resample_weighted_energy_identity
    (M : FiniteWeightedObservation S O) (f : S → ℝ) :
    M.weightedEnergy f =
      M.weightedEnergy (M.resample f) +
        M.weightedEnergy (fun x => f x - M.resample f x) := by
  let p : S → ℝ := M.resample f
  have hpidem : M.resample p = p := by
    dsimp [p]
    exact M.resample_idempotent f
  have hcross : M.weightedInner f p = M.weightedInner p p := by
    have hsa := M.resample_weighted_self_adjoint f p
    rw [hpidem] at hsa
    exact hsa.symm
  change M.weightedInner f f =
    M.weightedInner p p +
      M.weightedInner (fun x => f x - p x) (fun x => f x - p x)
  rw [M.weightedInner_sub_left, M.weightedInner_sub_right,
    M.weightedInner_sub_right, M.weightedInner_comm p f, hcross]
  ring

/-- Conditional resampling contracts squared finite weighted `L2` energy. -/
theorem resample_weighted_energy_le
    (M : FiniteWeightedObservation S O) (f : S → ℝ) :
    M.weightedEnergy (M.resample f) ≤ M.weightedEnergy f := by
  calc
    M.weightedEnergy (M.resample f) ≤
        M.weightedEnergy (M.resample f) +
          M.weightedEnergy (fun x => f x - M.resample f x) :=
      le_add_of_nonneg_right (M.weightedEnergy_nonneg _)
    _ = M.weightedEnergy f := (M.resample_weighted_energy_identity f).symm

end FiniteWeightedObservation

end ObservableNormalForms
