import ObservableNormalForms.MechanismVariants
import ObservableNormalForms.Stochastic


/-!
# Protected behavior obstructions and exact stochastic transport

Canonical two-file generic source. It constructs the finite-horizon first-hit
law from `pathWeight`, its bounded limit and derived transport interfaces,
then exposes the T0--T8, common-carrier, and quantitative public surface.
-/

namespace ObservableNormalForms.ProtectedObstructions

open scoped BigOperators

universe u

open ObservableNormalForms
open ObservableNormalForms.FiniteMarkovKernel

variable {S : Type u} [Fintype S] [DecidableEq S]

/-- The path first enters `target` exactly after `n` transitions. -/
def FirstHitAt (target : S → Prop) [DecidablePred target] :
    (n : ℕ) → S → (Fin n → S) → Prop
  | 0, x, _ => target x
  | n + 1, x, p =>
      ¬ target x ∧ FirstHitAt target n (p 0) (Fin.tail p)

private def finPiConsEquiv (S : Type u) (n : ℕ) :
    S × (Fin n → S) ≃ (Fin (n + 1) → S) where
  toFun zq := Fin.cons zq.1 zq.2
  invFun p := (p 0, Fin.tail p)
  left_inv := fun zq => by simp
  right_inv := fun p => Fin.cons_self_tail p

private theorem sum_pi_fin_succ (n : ℕ) (f : (Fin (n + 1) → S) → ℝ) :
    ∑ p : Fin (n + 1) → S, f p =
      ∑ z : S, ∑ q : Fin n → S, f (Fin.cons z q) := by
  calc
    ∑ p : Fin (n + 1) → S, f p =
        ∑ zq : S × (Fin n → S), f (Fin.cons zq.1 zq.2) :=
      ((finPiConsEquiv S n).sum_comp f).symm
    _ = ∑ z : S, ∑ q : Fin n → S, f (Fin.cons z q) := by
      rw [Fintype.sum_prod_type]

/-- Endpoint-resolved mass of paths whose first target hit is exactly at `n`. -/
noncomputable def firstHitAtMass (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (n : ℕ) (x c : S) : ℝ := by
  classical
  exact ∑ p : Fin n → S,
    if FirstHitAt target n x p ∧ pathEndpoint n x p = c
      then K.pathWeight n x p else 0

theorem firstHitAtMass_nonneg (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (n : ℕ) (x c : S) :
    0 ≤ firstHitAtMass K target n x c := by
  classical
  exact Finset.sum_nonneg fun p _ => by
    split_ifs
    · exact K.pathWeight_nonneg n x p
    · exact le_rfl

theorem firstHitAtMass_zero (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (x c : S) :
    firstHitAtMass K target 0 x c = if target x ∧ x = c then 1 else 0 := by
  classical
  simp [firstHitAtMass, FirstHitAt, pathEndpoint, FiniteMarkovKernel.pathWeight]

theorem firstHitAtMass_succ (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (n : ℕ) (x c : S) :
    firstHitAtMass K target (n + 1) x c =
      if target x then 0
      else ∑ z : S, K.probability x z * firstHitAtMass K target n z c := by
  classical
  rw [firstHitAtMass, sum_pi_fin_succ]
  by_cases hx : target x
  · simp [FirstHitAt, hx]
  · simp [FirstHitAt, hx, FiniteMarkovKernel.pathWeight_cons,
      pathEndpoint_cons, firstHitAtMass, Finset.mul_sum]

/-- Fixed-horizon cumulative endpoint mass through horizon `N`. -/
noncomputable def cumulativeEndpointMass (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (N : ℕ) (x c : S) : ℝ :=
  ∑ n ∈ Finset.range (N + 1), firstHitAtMass K target n x c

theorem cumulativeEndpointMass_nonneg (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (N : ℕ) (x c : S) :
    0 ≤ cumulativeEndpointMass K target N x c := by
  classical
  exact Finset.sum_nonneg fun n _ => firstHitAtMass_nonneg K target n x c

theorem cumulativeEndpointMass_mono (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (N : ℕ) (x c : S) :
    cumulativeEndpointMass K target N x c ≤
      cumulativeEndpointMass K target (N + 1) x c := by
  classical
  unfold cumulativeEndpointMass
  have h := le_add_of_nonneg_right
    (a := ∑ n ∈ Finset.range (N + 1), firstHitAtMass K target n x c)
    (firstHitAtMass_nonneg K target (N + 1) x c)
  simpa [Finset.sum_range_succ, Nat.add_assoc] using h

/-- Total probability of a first hit by the fixed horizon. -/
noncomputable def hitBy (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (N : ℕ) (x : S) : ℝ :=
  ∑ c : S, cumulativeEndpointMass K target N x c

theorem hitBy_nonneg (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (N : ℕ) (x : S) : 0 ≤ hitBy K target N x := by
  classical
  exact Finset.sum_nonneg fun c _ => cumulativeEndpointMass_nonneg K target N x c

theorem hitBy_mono (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (N : ℕ) (x : S) : hitBy K target N x ≤ hitBy K target (N + 1) x := by
  classical
  exact Finset.sum_le_sum fun c _ => cumulativeEndpointMass_mono K target N x c

/-- First-hit mass at an exact time, with the endpoint summed out. -/
noncomputable def firstHitTotal (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (n : ℕ) (x : S) : ℝ := by
  classical
  exact ∑ p : Fin n → S,
    if FirstHitAt target n x p then K.pathWeight n x p else 0

theorem sum_firstHitAtMass (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (n : ℕ) (x : S) :
    (∑ c : S, firstHitAtMass K target n x c) = firstHitTotal K target n x := by
  classical
  unfold firstHitAtMass firstHitTotal
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro p _
  by_cases hp : FirstHitAt target n x p
  · simp [hp]
  · simp [hp]

theorem firstHitTotal_zero (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] (x : S) :
    firstHitTotal K target 0 x = if target x then 1 else 0 := by
  classical
  simp [firstHitTotal, FirstHitAt, FiniteMarkovKernel.pathWeight]

theorem firstHitTotal_succ (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (n : ℕ) (x : S) :
    firstHitTotal K target (n + 1) x =
      if target x then 0
      else ∑ z : S, K.probability x z * firstHitTotal K target n z := by
  classical
  rw [firstHitTotal, sum_pi_fin_succ]
  by_cases hx : target x
  · simp [FirstHitAt, hx]
  · simp [FirstHitAt, hx, FiniteMarkovKernel.pathWeight_cons,
      firstHitTotal, Finset.mul_sum]

theorem hitBy_eq_sum_firstHitTotal (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (N : ℕ) (x : S) :
    hitBy K target N x =
      ∑ n ∈ Finset.range (N + 1), firstHitTotal K target n x := by
  classical
  unfold hitBy cumulativeEndpointMass
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro n hn
  exact sum_firstHitAtMass K target n x

theorem hitBy_succ_recurrence (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (N : ℕ) (x : S) :
    hitBy K target (N + 1) x =
      if target x then 1 else ∑ z : S, K.probability x z * hitBy K target N z := by
  classical
  rw [hitBy_eq_sum_firstHitTotal]
  rw [show N + 1 + 1 = (N + 1) + 1 by omega]
  rw [Finset.sum_range_succ']
  rw [firstHitTotal_zero]
  by_cases hx : target x
  · simp [hx, firstHitTotal_succ]
  · rw [if_neg hx, if_neg hx]
    simp only [firstHitTotal_succ, hx, if_false]
    rw [Finset.sum_comm]
    simp only [add_zero]
    apply Finset.sum_congr rfl
    intro z _
    rw [← Finset.mul_sum]
    congr 1
    symm
    exact hitBy_eq_sum_firstHitTotal K target N z

theorem hitBy_le_one (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] :
    ∀ (N : ℕ) (x : S), hitBy K target N x ≤ 1 := by
  intro N
  induction N with
  | zero =>
      intro x
      rw [hitBy_eq_sum_firstHitTotal]
      by_cases hx : target x <;> simp [firstHitTotal_zero, hx]
  | succ N ih =>
      intro x
      rw [hitBy_succ_recurrence]
      by_cases hx : target x
      · simp [hx]
      · rw [if_neg hx]
        calc
          ∑ z : S, K.probability x z * hitBy K target N z ≤
              ∑ z : S, K.probability x z * 1 := by
            apply Finset.sum_le_sum
            intro z _
            exact mul_le_mul_of_nonneg_left (ih z) (K.probability_nonneg x z)
          _ = 1 := by simp [K.probability_sum_one x]

theorem hitBy_monotone (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] (x : S) :
    Monotone (fun N => hitBy K target N x) :=
  monotone_nat_of_le_succ fun N => hitBy_mono K target N x

theorem hitBy_bddAbove (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] (x : S) :
    BddAbove (Set.range fun N => hitBy K target N x) := by
  refine ⟨1, ?_⟩
  rintro _ ⟨N, rfl⟩
  exact hitBy_le_one K target N x

/-- Canonical bounded monotone limit of the path-weight approximants. -/
noncomputable def hitProbability (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] (x : S) : ℝ :=
  ⨆ N, hitBy K target N x

theorem hitBy_tendsto_hitProbability (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] (x : S) :
    Filter.Tendsto (fun N => hitBy K target N x) Filter.atTop
      (nhds (hitProbability K target x)) := by
  exact tendsto_atTop_ciSup (hitBy_monotone K target x) (hitBy_bddAbove K target x)

theorem hitBy_le_hitProbability (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (N : ℕ) (x : S) : hitBy K target N x ≤ hitProbability K target x := by
  exact le_ciSup (hitBy_bddAbove K target x) N

theorem hitProbability_nonneg (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] (x : S) :
    0 ≤ hitProbability K target x :=
  (hitBy_nonneg K target 0 x).trans (hitBy_le_hitProbability K target 0 x)

theorem hitProbability_le_one (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] (x : S) :
    hitProbability K target x ≤ 1 := by
  exact ciSup_le fun N => hitBy_le_one K target N x

/-- The path-weight limit satisfies the target-boundary Bellman equation. -/
theorem hitProbability_bellman (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] (x : S) :
    hitProbability K target x =
      if target x then 1
      else ∑ z : S, K.probability x z * hitProbability K target z := by
  refine tendsto_nhds_unique (hitBy_tendsto_hitProbability K target x)
    ((Filter.tendsto_add_atTop_iff_nat 1).1 ?_)
  simp only [hitBy_succ_recurrence]
  by_cases hx : target x
  · simpa [hx] using (tendsto_const_nhds :
      Filter.Tendsto (fun _ : ℕ => (1 : ℝ)) Filter.atTop (nhds 1))
  · simp only [hx, if_false]
    exact tendsto_finset_sum Finset.univ fun z _ =>
      Filter.Tendsto.const_mul (K.probability x z)
        (hitBy_tendsto_hitProbability K target z)

theorem cumulativeEndpointMass_le_hitBy (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (N : ℕ) (x c : S) :
    cumulativeEndpointMass K target N x c ≤ hitBy K target N x := by
  classical
  unfold hitBy
  exact Finset.single_le_sum
    (fun y _ => cumulativeEndpointMass_nonneg K target N x y)
    (Finset.mem_univ c)

theorem cumulativeEndpointMass_bddAbove (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (x c : S) :
    BddAbove (Set.range fun N => cumulativeEndpointMass K target N x c) := by
  refine ⟨1, ?_⟩
  rintro _ ⟨N, rfl⟩
  exact (cumulativeEndpointMass_le_hitBy K target N x c).trans
    (hitBy_le_one K target N x)

/-- Endpoint-resolved bounded monotone first-hit law. -/
noncomputable def endpointHitProbability (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (x c : S) : ℝ :=
  ⨆ N, cumulativeEndpointMass K target N x c

theorem cumulativeEndpointMass_tendsto (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (x c : S) :
    Filter.Tendsto (fun N => cumulativeEndpointMass K target N x c) Filter.atTop
      (nhds (endpointHitProbability K target x c)) := by
  exact tendsto_atTop_ciSup
    (monotone_nat_of_le_succ fun N => cumulativeEndpointMass_mono K target N x c)
    (cumulativeEndpointMass_bddAbove K target x c)

def SupportedFinPath (K : FiniteMarkovKernel S) :
    (n : ℕ) → S → (Fin n → S) → Prop
  | 0, _, _ => True
  | n + 1, x, p =>
      0 < K.probability x (p 0) ∧ SupportedFinPath K n (p 0) (Fin.tail p)

theorem pathWeight_pos_iff_supportedFinPath (K : FiniteMarkovKernel S) :
    ∀ (n : ℕ) (x : S) (p : Fin n → S),
      0 < K.pathWeight n x p ↔ SupportedFinPath K n x p := by
  intro n
  induction n with
  | zero =>
      intro x p
      simp [FiniteMarkovKernel.pathWeight, SupportedFinPath]
  | succ n ih =>
      intro x p
      simp only [FiniteMarkovKernel.pathWeight, SupportedFinPath]
      constructor
      · intro h
        rcases (mul_pos_iff.mp h) with hpos | hneg
        · exact ⟨hpos.1, (ih (p 0) (Fin.tail p)).mp hpos.2⟩
        · exact False.elim ((not_lt_of_ge (K.probability_nonneg x (p 0))) hneg.1)
      · rintro ⟨hstep, htail⟩
        exact mul_pos hstep ((ih (p 0) (Fin.tail p)).mpr htail)

theorem firstHitAtMass_pos_iff (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (n : ℕ) (x c : S) :
    0 < firstHitAtMass K target n x c ↔
      ∃ p : Fin n → S,
        FirstHitAt target n x p ∧ pathEndpoint n x p = c ∧
          SupportedFinPath K n x p := by
  classical
  unfold firstHitAtMass
  let f : (Fin n → S) → ℝ := fun p =>
    if FirstHitAt target n x p ∧ pathEndpoint n x p = c
      then K.pathWeight n x p else 0
  have hf : ∀ p ∈ (Finset.univ : Finset (Fin n → S)), 0 ≤ f p := by
    intro p _
    dsimp [f]
    split_ifs
    · exact K.pathWeight_nonneg n x p
    · exact le_rfl
  constructor
  · intro hsum
    have hne : (∑ p : Fin n → S, f p) ≠ 0 := ne_of_gt hsum
    have hex : ∃ p ∈ (Finset.univ : Finset (Fin n → S)), f p ≠ 0 := by
      by_contra hnone
      push Not at hnone
      exact hne (Finset.sum_eq_zero fun p hp => hnone p hp)
    rcases hex with ⟨p, _, hpne⟩
    have hppos : 0 < f p := lt_of_le_of_ne (hf p (Finset.mem_univ p)) (Ne.symm hpne)
    dsimp [f] at hppos
    split_ifs at hppos with hev
    · exact ⟨p, hev.1, hev.2,
        (pathWeight_pos_iff_supportedFinPath K n x p).mp hppos⟩
    · exact False.elim (lt_irrefl 0 hppos)
  · rintro ⟨p, hfirst, hend, hsupp⟩
    apply Finset.sum_pos' hf
    refine ⟨p, Finset.mem_univ p, ?_⟩
    dsimp [f]
    rw [if_pos ⟨hfirst, hend⟩]
    exact (pathWeight_pos_iff_supportedFinPath K n x p).mpr hsupp

theorem cumulativeEndpointMass_pos_iff (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (N : ℕ) (x c : S) :
    0 < cumulativeEndpointMass K target N x c ↔
      ∃ n ∈ Finset.range (N + 1), 0 < firstHitAtMass K target n x c := by
  classical
  unfold cumulativeEndpointMass
  constructor
  · intro hsum
    have hne : (∑ n ∈ Finset.range (N + 1), firstHitAtMass K target n x c) ≠ 0 :=
      ne_of_gt hsum
    by_contra hnone
    push Not at hnone
    apply hne
    exact Finset.sum_eq_zero fun n hn =>
      le_antisymm (hnone n hn) (firstHitAtMass_nonneg K target n x c)
  · rintro ⟨n, hn, hpos⟩
    exact Finset.sum_pos'
      (fun m _ => firstHitAtMass_nonneg K target m x c) ⟨n, hn, hpos⟩

/-- T2: endpoint law is positive exactly when a positive supported first-hit path exists. -/
theorem endpointHitProbability_pos_iff (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (x c : S) :
    0 < endpointHitProbability K target x c ↔
      ∃ (n : ℕ) (p : Fin n → S),
        FirstHitAt target n x p ∧ pathEndpoint n x p = c ∧
          SupportedFinPath K n x p := by
  classical
  rw [endpointHitProbability, lt_ciSup_iff (cumulativeEndpointMass_bddAbove K target x c)]
  constructor
  · rintro ⟨N, hN⟩
    rcases (cumulativeEndpointMass_pos_iff K target N x c).mp hN with ⟨n, _, hn⟩
    rcases (firstHitAtMass_pos_iff K target n x c).mp hn with ⟨p, hp⟩
    exact ⟨n, p, hp⟩
  · rintro ⟨n, p, hp⟩
    refine ⟨n, (cumulativeEndpointMass_pos_iff K target n x c).mpr ?_⟩
    exact ⟨n, by simp, (firstHitAtMass_pos_iff K target n x c).mpr ⟨p, hp⟩⟩

abbrev SupportStep (K : FiniteMarkovKernel S) (x y : S) : Prop :=
  0 < K.probability x y

abbrev PreHitStep (K : FiniteMarkovKernel S)
    (target : S → Prop) (x y : S) : Prop :=
  ¬ target x ∧ SupportStep K x y

abbrev PreHitReach (K : FiniteMarkovKernel S)
    (target : S → Prop) (x y : S) : Prop :=
  Relation.ReflTransGen (PreHitStep K target) x y

/-- A nonempty target-free positive-support closed set reachable before hit. -/
structure ReachableClosedTrap (K : FiniteMarkovKernel S)
    (target : S → Prop) (root : S) where
  carrier : Set S
  nonempty : carrier.Nonempty
  target_free : ∀ u ∈ carrier, ¬ target u
  support_closed : ∀ u ∈ carrier, ∀ v, SupportStep K u v → v ∈ carrier
  reachable : ∃ u ∈ carrier, PreHitReach K target root u

noncomputable def deficit (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] (x : S) : ℝ :=
  1 - hitProbability K target x

theorem deficit_nonneg (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] (x : S) :
    0 ≤ deficit K target x :=
  sub_nonneg.mpr (hitProbability_le_one K target x)

theorem deficit_le_one (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] (x : S) :
    deficit K target x ≤ 1 := by
  unfold deficit
  linarith [hitProbability_nonneg K target x]

theorem deficit_bellman (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] (x : S) :
    deficit K target x =
      if target x then 0
      else ∑ z : S, K.probability x z * deficit K target z := by
  rw [deficit, hitProbability_bellman]
  by_cases hx : target x
  · simp [hx]
  · rw [if_neg hx, if_neg hx]
    unfold deficit
    calc
      1 - ∑ z : S, K.probability x z * hitProbability K target z =
          (∑ z : S, K.probability x z) -
            ∑ z : S, K.probability x z * hitProbability K target z := by
        rw [K.probability_sum_one]
      _ = ∑ z : S,
          (K.probability x z - K.probability x z * hitProbability K target z) := by
        rw [Finset.sum_sub_distrib]
      _ = ∑ z : S, K.probability x z *
          (1 - hitProbability K target z) := by
        apply Finset.sum_congr rfl
        intro z _
        ring

theorem hitProbability_eq_one_of_supportStep
    (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    {x y : S} (hx : hitProbability K target x = 1)
    (hstep : PreHitStep K target x y) :
    hitProbability K target y = 1 := by
  have hdx : deficit K target x = 0 := by simp [deficit, hx]
  have hsum : (∑ z : S, K.probability x z * deficit K target z) = 0 := by
    calc
      _ = deficit K target x := by
        rw [deficit_bellman K target x, if_neg hstep.1]
      _ = 0 := hdx
  have hterm : K.probability x y * deficit K target y = 0 :=
    (Finset.sum_eq_zero_iff_of_nonneg
      (fun z _ => mul_nonneg (K.probability_nonneg x z)
        (deficit_nonneg K target z))).mp hsum y (Finset.mem_univ y)
  have hdy : deficit K target y = 0 := (mul_eq_zero.mp hterm).resolve_left
    (ne_of_gt hstep.2)
  unfold deficit at hdy
  linarith

theorem hitProbability_eq_one_of_preHitReach
    (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    {x y : S} (hx : hitProbability K target x = 1)
    (hreach : PreHitReach K target x y) :
    hitProbability K target y = 1 := by
  induction hreach with
  | refl => exact hx
  | tail hxy hyz ih =>
      exact hitProbability_eq_one_of_supportStep K target ih hyz

theorem hitBy_eq_zero_on_closed
    (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (U : Set S)
    (htarget : ∀ u ∈ U, ¬ target u)
    (hclosed : ∀ u ∈ U, ∀ v, SupportStep K u v → v ∈ U) :
    ∀ (N : ℕ) (u : S), u ∈ U → hitBy K target N u = 0 := by
  intro N
  induction N with
  | zero =>
      intro u hu
      rw [hitBy_eq_sum_firstHitTotal]
      simp [firstHitTotal_zero, htarget u hu]
  | succ N ih =>
      intro u hu
      rw [hitBy_succ_recurrence, if_neg (htarget u hu)]
      apply Finset.sum_eq_zero
      intro v _
      by_cases hp : 0 < K.probability u v
      · rw [ih v (hclosed u hu v hp), mul_zero]
      · have hz : K.probability u v = 0 :=
          le_antisymm (le_of_not_gt hp) (K.probability_nonneg u v)
        rw [hz, zero_mul]

theorem hitProbability_eq_zero_on_closed
    (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (U : Set S)
    (htarget : ∀ u ∈ U, ¬ target u)
    (hclosed : ∀ u ∈ U, ∀ v, SupportStep K u v → v ∈ U)
    (u : S) (hu : u ∈ U) :
    hitProbability K target u = 0 := by
  refine tendsto_nhds_unique (hitBy_tendsto_hitProbability K target u) ?_
  have heq : (fun N => hitBy K target N u) = (fun _ : ℕ => (0 : ℝ)) := by
    funext N
    exact hitBy_eq_zero_on_closed K target U htarget hclosed N u hu
  rw [heq]
  exact tendsto_const_nhds

theorem reachableClosedTrap_imp_hitProbability_ne_one
    (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (root : S) (trap : ReachableClosedTrap K target root) :
    hitProbability K target root ≠ 1 := by
  intro hroot
  rcases trap.reachable with ⟨u, hu, hreach⟩
  have huone := hitProbability_eq_one_of_preHitReach K target hroot hreach
  have huzero := hitProbability_eq_zero_on_closed K target trap.carrier
    trap.target_free trap.support_closed u hu
  linarith

noncomputable def preHitReachableFinset (K : FiniteMarkovKernel S)
    (target : S → Prop) (root : S) : Finset S := by
  classical
  exact Finset.univ.filter fun u => PreHitReach K target root u

theorem mem_preHitReachableFinset (K : FiniteMarkovKernel S)
    (target : S → Prop) (root u : S) :
    u ∈ preHitReachableFinset K target root ↔ PreHitReach K target root u := by
  classical
  simp [preHitReachableFinset]

theorem reachableClosedTrap_of_hitProbability_lt_one
    (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (root : S) (hroot : hitProbability K target root < 1) :
    Nonempty (ReachableClosedTrap K target root) := by
  classical
  let R := preHitReachableFinset K target root
  have hrootR : root ∈ R := by
    rw [show R = preHitReachableFinset K target root from rfl,
      mem_preHitReachableFinset]
  obtain ⟨u, huR, hmax⟩ := R.exists_max_image (deficit K target) ⟨root, hrootR⟩
  let M := deficit K target u
  have hMpos : 0 < M := by
    have hrootpos : 0 < deficit K target root := by
      unfold deficit
      linarith
    exact hrootpos.trans_le (hmax root hrootR)
  let U : Set S := {v | v ∈ R ∧ deficit K target v = M}
  have huU : u ∈ U := ⟨huR, rfl⟩
  have target_free : ∀ v ∈ U, ¬ target v := by
    intro v hv htarget
    have hvone : hitProbability K target v = 1 := by
      simpa [htarget] using hitProbability_bellman K target v
    have hvzero : deficit K target v = 0 := by simp [deficit, hvone]
    rw [hv.2] at hvzero
    linarith
  have support_closed :
      ∀ v ∈ U, ∀ w, SupportStep K v w → w ∈ U := by
    intro v hv w hvw
    have hvreach : PreHitReach K target root v :=
      (mem_preHitReachableFinset K target root v).mp hv.1
    have hwreach : PreHitReach K target root w :=
      Relation.ReflTransGen.tail hvreach ⟨target_free v hv, hvw⟩
    have hwR : w ∈ R :=
      (mem_preHitReachableFinset K target root w).mpr hwreach
    refine ⟨hwR, ?_⟩
    have hwle : deficit K target w ≤ M := hmax w hwR
    apply le_antisymm hwle
    by_contra hnot
    have hwlt : deficit K target w < M := lt_of_not_ge hnot
    have hall : ∀ z ∈ (Finset.univ : Finset S),
        K.probability v z * deficit K target z ≤ K.probability v z * M := by
      intro z _
      by_cases hz : 0 < K.probability v z
      · have hzreach : PreHitReach K target root z :=
          Relation.ReflTransGen.tail hvreach ⟨target_free v hv, hz⟩
        have hzR : z ∈ R :=
          (mem_preHitReachableFinset K target root z).mpr hzreach
        exact mul_le_mul_of_nonneg_left (hmax z hzR)
          (K.probability_nonneg v z)
      · have hz0 : K.probability v z = 0 :=
          le_antisymm (le_of_not_gt hz) (K.probability_nonneg v z)
        simp [hz0]
    have hstrict : ∃ z ∈ (Finset.univ : Finset S),
        K.probability v z * deficit K target z < K.probability v z * M :=
      ⟨w, Finset.mem_univ w, mul_lt_mul_of_pos_left hwlt hvw⟩
    have hsumlt :
        (∑ z : S, K.probability v z * deficit K target z) <
          ∑ z : S, K.probability v z * M :=
      Finset.sum_lt_sum hall hstrict
    have hsumM : (∑ z : S, K.probability v z * M) = M := by
      rw [← Finset.sum_mul, K.probability_sum_one]
      simp
    have hbell : deficit K target v =
        ∑ z : S, K.probability v z * deficit K target z := by
      simpa [target_free v hv] using deficit_bellman K target v
    rw [hv.2] at hbell
    linarith
  refine ⟨{
    carrier := U
    nonempty := ⟨u, huU⟩
    target_free := target_free
    support_closed := support_closed
    reachable := ⟨u, huU,
      (mem_preHitReachableFinset K target root u).mp huR⟩ }⟩

/-- T3: almost-sure first hit iff no reachable target-avoiding closed trap exists. -/
theorem hitProbability_eq_one_iff_no_reachableClosedTrap
    (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (root : S) :
    hitProbability K target root = 1 ↔
      ¬ Nonempty (ReachableClosedTrap K target root) := by
  constructor
  · intro hone htrap
    exact reachableClosedTrap_imp_hitProbability_ne_one K target root htrap.some hone
  · intro hnotrap
    apply le_antisymm (hitProbability_le_one K target root)
    by_contra hnot
    have hlt : hitProbability K target root < 1 := lt_of_not_ge hnot
    exact hnotrap (reachableClosedTrap_of_hitProbability_lt_one K target root hlt)

/-- Target-boundary nonnegative Bellman fixed-point interface. -/
def IsNonnegativeBellmanFixedPoint (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] (g : S → ℝ) : Prop :=
  (∀ x, 0 ≤ g x) ∧
  ∀ x, g x = if target x then 1 else ∑ z : S, K.probability x z * g z

theorem hitBy_zero (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] (x : S) :
    hitBy K target 0 x = if target x then 1 else 0 := by
  rw [hitBy_eq_sum_firstHitTotal]
  simp [firstHitTotal_zero]

theorem hitBy_le_nonnegativeBellmanFixedPoint
    (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (g : S → ℝ) (hg : IsNonnegativeBellmanFixedPoint K target g) :
    ∀ (N : ℕ) (x : S), hitBy K target N x ≤ g x := by
  intro N
  induction N with
  | zero =>
      intro x
      rw [hitBy_zero]
      by_cases hx : target x
      · simpa [hx] using le_of_eq (hg.2 x).symm
      · simpa [hx] using hg.1 x
  | succ N ih =>
      intro x
      rw [hitBy_succ_recurrence, hg.2 x]
      by_cases hx : target x
      · simp [hx]
      · simp only [hx, if_false]
        apply Finset.sum_le_sum
        intro z _
        exact mul_le_mul_of_nonneg_left (ih z) (K.probability_nonneg x z)

theorem hitProbability_isNonnegativeBellmanFixedPoint
    (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] :
    IsNonnegativeBellmanFixedPoint K target (hitProbability K target) :=
  ⟨hitProbability_nonneg K target, hitProbability_bellman K target⟩

/-- The path-weight limit is the least nonnegative Bellman fixed point. -/
theorem hitProbability_least_nonnegative_fixedPoint
    (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (g : S → ℝ) (hg : IsNonnegativeBellmanFixedPoint K target g) :
    ∀ x, hitProbability K target x ≤ g x := by
  intro x
  exact ciSup_le fun N => hitBy_le_nonnegativeBellmanFixedPoint K target g hg N x

theorem constantOne_isNonnegativeBellmanFixedPoint
    (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target] :
    IsNonnegativeBellmanFixedPoint K target (fun _ => 1) := by
  constructor
  · intro x
    norm_num
  · intro x
    by_cases hx : target x
    · simp [hx]
    · simp [hx, K.probability_sum_one x]

/-- Closed-trap cases explicitly witness failure of general fixed-point uniqueness. -/
theorem nonunique_nonnegativeBellmanFixedPoint_of_reachableClosedTrap
    (K : FiniteMarkovKernel S)
    (target : S → Prop) [DecidablePred target]
    (root : S) (trap : ReachableClosedTrap K target root) :
    hitProbability K target ≠ (fun _ => 1) ∧
      IsNonnegativeBellmanFixedPoint K target (fun _ => 1) := by
  constructor
  · intro heq
    apply reachableClosedTrap_imp_hitProbability_ne_one K target root trap
    exact congrFun heq root
  · exact constantOne_isNonnegativeBellmanFixedPoint K target

end ObservableNormalForms.ProtectedObstructions


/-!
# Nonidentity structural exactness

The structural record contains no direct fiber, positive, almost-sure,
selection, or endpoint-law pushforward fields. Exact-time, fixed-horizon, and
limit pushforward are derived from kernel lumping and target reflection;
target coverage derives fiber nonemptiness. Observation commutation and the
`RealizedBehavior`/`BehaviorCut` consumer remain separate interfaces.
-/

namespace ObservableNormalForms.ProtectedObstructions.NonidentityExactness

open scoped BigOperators
open ObservableNormalForms
open ObservableNormalForms.FiniteMarkovKernel

universe u v w x

structure FixedBehaviorModel (O : Type u) (Q : Type v)
    [Fintype O] [DecidableEq O] [Fintype Q] [DecidableEq Q] where
  protectedSet : Finset O
  initial : O → Finset Q
  initial_nonempty : ∀ b ∈ protectedSet, (initial b).Nonempty
  target : O → Q → Bool
  kernel : FiniteMarkovKernel Q
  silent : Setoid Q

namespace FixedBehaviorModel

variable {O : Type u} {Q : Type v}
  [Fintype O] [DecidableEq O] [Fintype Q] [DecidableEq Q]

def IsTarget (M : FixedBehaviorModel O Q) (b : O) (q : Q) : Prop :=
  M.target b q = true

instance instDecidablePredIsTarget (M : FixedBehaviorModel O Q) (b : O) :
    DecidablePred (M.IsTarget b) := fun q =>
  inferInstanceAs (Decidable (M.target b q = true))

def FiberNonempty (M : FixedBehaviorModel O Q) (b : O) : Prop :=
  ∃ q, M.IsTarget b q

noncomputable def endpointMass (M : FixedBehaviorModel O Q)
    (b : O) (x c : Q) : ℝ :=
  endpointHitProbability M.kernel (M.IsTarget b) x c

noncomputable def hit (M : FixedBehaviorModel O Q) (b : O) (x : Q) : ℝ :=
  hitProbability M.kernel (M.IsTarget b) x

def Positive (M : FixedBehaviorModel O Q) (b : O) : Prop :=
  ∃ x ∈ M.initial b, 0 < M.hit b x

def AlmostSure (M : FixedBehaviorModel O Q) (b : O) : Prop :=
  ∀ x ∈ M.initial b, M.hit b x = 1

/-- Pairwise form of a nonempty singleton quotient endpoint support. -/
def EndpointUnique (M : FixedBehaviorModel O Q) (b : O) : Prop :=
  (∃ x ∈ M.initial b, ∃ c, 0 < M.endpointMass b x c) ∧
  ∀ x ∈ M.initial b, ∀ c, 0 < M.endpointMass b x c →
    ∀ y ∈ M.initial b, ∀ d, 0 < M.endpointMass b y d → M.silent.r c d

def InL0 (M : FixedBehaviorModel O Q) (b : O) : Prop :=
  b ∈ M.protectedSet ∧ M.FiberNonempty b

def InL1 (M : FixedBehaviorModel O Q) (b : O) : Prop :=
  M.InL0 b ∧ M.Positive b

def InL2 (M : FixedBehaviorModel O Q) (b : O) : Prop :=
  M.InL1 b ∧ M.AlmostSure b

def InL3 (M : FixedBehaviorModel O Q) (b : O) : Prop :=
  M.InL2 b ∧ M.EndpointUnique b

theorem endpointMass_nonneg (M : FixedBehaviorModel O Q)
    (b : O) (x c : Q) : 0 ≤ M.endpointMass b x c := by
  unfold endpointMass endpointHitProbability
  exact (cumulativeEndpointMass_nonneg M.kernel (M.IsTarget b) 0 x c).trans
    (le_ciSup (cumulativeEndpointMass_bddAbove M.kernel (M.IsTarget b) x c) 0)

end FixedBehaviorModel

variable {O₁ : Type u} {Q₁ : Type v} {O₂ : Type w} {Q₂ : Type x}
  [Fintype O₁] [DecidableEq O₁] [Fintype Q₁] [DecidableEq Q₁]
  [Fintype O₂] [DecidableEq O₂] [Fintype Q₂] [DecidableEq Q₂]

/-- Minimal fixed-behavior stochastic structure.  No semantic conclusion is a field. -/
structure StructuralMorphism
    (A : FixedBehaviorModel O₁ Q₁) (B : FixedBehaviorModel O₂ Q₂) where
  stateMap : Q₁ → Q₂
  behaviorMap : O₁ → O₂
  protected_reflect : ∀ b, b ∈ A.protectedSet ↔ behaviorMap b ∈ B.protectedSet
  protected_cover : ∀ d ∈ B.protectedSet, ∃ b ∈ A.protectedSet, behaviorMap b = d
  initial_reflect : ∀ b x, b ∈ A.protectedSet →
    (x ∈ A.initial b ↔ stateMap x ∈ B.initial (behaviorMap b))
  initial_cover : ∀ b y, b ∈ A.protectedSet → y ∈ B.initial (behaviorMap b) →
    ∃ x ∈ A.initial b, stateMap x = y
  kernel_lumping : ∀ x (φ : Q₂ → ℝ),
    B.kernel.apply φ (stateMap x) = A.kernel.apply (φ ∘ stateMap) x
  target_reflect : ∀ b x, b ∈ A.protectedSet →
    (A.IsTarget b x ↔ B.IsTarget (behaviorMap b) (stateMap x))
  target_cover : ∀ b y, b ∈ A.protectedSet → B.IsTarget (behaviorMap b) y →
    ∃ x, A.IsTarget b x ∧ stateMap x = y
  quotient_exact : ∀ c d,
    A.silent.r c d ↔ B.silent.r (stateMap c) (stateMap d)

/-- Observation ownership is deliberately outside fixed-behavior transport. -/
structure ObservationCompatibility
    {A : FixedBehaviorModel O₁ Q₁} {B : FixedBehaviorModel O₂ Q₂}
    (F : StructuralMorphism A B)
    (observeA : Q₁ → O₁) (observeB : Q₂ → O₂) : Prop where
  observation_commute : ∀ x,
    observeB (F.stateMap x) = F.behaviorMap (observeA x)

namespace StructuralMorphism

variable {A : FixedBehaviorModel O₁ Q₁} {B : FixedBehaviorModel O₂ Q₂}
  (F : StructuralMorphism A B)

noncomputable def fiberSum (d : Q₂) (μ : Q₁ → ℝ) : ℝ :=
  ∑ c : Q₁, if F.stateMap c = d then μ c else 0

private theorem fiberSum_nonneg (d : Q₂) (μ : Q₁ → ℝ)
    (hμ : ∀ c, 0 ≤ μ c) : 0 ≤ F.fiberSum d μ := by
  classical
  exact Finset.sum_nonneg fun c _ => by
    by_cases hc : F.stateMap c = d <;> simp [fiberSum, hc, hμ c]

private theorem fiberSum_pos_of (d : Q₂) (μ : Q₁ → ℝ)
    (hμ : ∀ c, 0 ≤ μ c) {c : Q₁}
    (hc : F.stateMap c = d) (hpos : 0 < μ c) :
    0 < F.fiberSum d μ := by
  classical
  unfold fiberSum
  apply Finset.sum_pos'
  · intro z _
    by_cases hz : F.stateMap z = d <;> simp [hz, hμ z]
  · exact ⟨c, Finset.mem_univ c, by simp [hc, hpos]⟩

private theorem exists_pos_of_fiberSum_pos (d : Q₂) (μ : Q₁ → ℝ)
    (hμ : ∀ c, 0 ≤ μ c) (hpos : 0 < F.fiberSum d μ) :
    ∃ c, F.stateMap c = d ∧ 0 < μ c := by
  classical
  by_contra hnone
  have hle : ∀ c, F.stateMap c = d → μ c ≤ 0 := by
    intro c hc
    exact le_of_not_gt fun hp => hnone ⟨c, hc, hp⟩
  have hzero : F.fiberSum d μ = 0 := by
    unfold fiberSum
    apply Finset.sum_eq_zero
    intro c _
    by_cases hc : F.stateMap c = d
    · have hm : μ c = 0 := le_antisymm (hle c hc) (hμ c)
      simp [hc, hm]
    · simp [hc]
  linarith

/-- Target reflection plus coverage derives the removed fiber exactness field. -/
theorem fiber_nonempty_exact_derived {b : O₁} (hb : b ∈ A.protectedSet) :
    A.FiberNonempty b ↔ B.FiberNonempty (F.behaviorMap b) := by
  constructor
  · rintro ⟨c, hc⟩
    exact ⟨F.stateMap c, (F.target_reflect b c hb).mp hc⟩
  · rintro ⟨d, hd⟩
    rcases F.target_cover b d hb hd with ⟨c, hc, _⟩
    exact ⟨c, hc⟩

/-- Exact-time endpoint pushforward, proved from target reflection and kernel lumping. -/
theorem firstHitAtMass_pushforward {b : O₁} (hb : b ∈ A.protectedSet) :
    ∀ (n : ℕ) (x : Q₁) (d : Q₂),
      firstHitAtMass B.kernel (B.IsTarget (F.behaviorMap b)) n (F.stateMap x) d =
        F.fiberSum d (fun c => firstHitAtMass A.kernel (A.IsTarget b) n x c) := by
  intro n
  induction n with
  | zero =>
      intro x d
      simp_rw [firstHitAtMass_zero]
      have ht := F.target_reflect b x hb
      classical
      by_cases hA : A.IsTarget b x
      · have hB : B.IsTarget (F.behaviorMap b) (F.stateMap x) := ht.mp hA
        unfold fiberSum
        symm
        rw [Finset.sum_eq_single x]
        · simp [hA, hB]
        · intro c _ hcx
          simp [Ne.symm hcx]
        · simp
      · have hB : ¬ B.IsTarget (F.behaviorMap b) (F.stateMap x) :=
          fun h => hA (ht.mpr h)
        simp [hA, hB, fiberSum]
  | succ n ih =>
      intro x d
      simp_rw [firstHitAtMass_succ]
      have ht := F.target_reflect b x hb
      by_cases hA : A.IsTarget b x
      · have hB : B.IsTarget (F.behaviorMap b) (F.stateMap x) := ht.mp hA
        simp [hA, hB, fiberSum]
      · have hB : ¬ B.IsTarget (F.behaviorMap b) (F.stateMap x) :=
          fun h => hA (ht.mpr h)
        simp_rw [hA, hB, if_false]
        change B.kernel.apply
            (fun z => firstHitAtMass B.kernel (B.IsTarget (F.behaviorMap b)) n z d)
            (F.stateMap x) = _
        rw [F.kernel_lumping]
        simp only [FiniteMarkovKernel.apply, Function.comp_apply, ih]
        classical
        unfold fiberSum
        calc
          (∑ y : Q₁, A.kernel.probability x y *
              ∑ c : Q₁, if F.stateMap c = d then
                firstHitAtMass A.kernel (A.IsTarget b) n y c else 0) =
              ∑ y : Q₁, ∑ c : Q₁, if F.stateMap c = d then
                A.kernel.probability x y *
                  firstHitAtMass A.kernel (A.IsTarget b) n y c else 0 := by
            apply Finset.sum_congr rfl
            intro y _
            rw [Finset.mul_sum]
            apply Finset.sum_congr rfl
            intro c _
            by_cases hc : F.stateMap c = d <;> simp [hc]
          _ = ∑ c : Q₁, ∑ y : Q₁, if F.stateMap c = d then
                A.kernel.probability x y *
                  firstHitAtMass A.kernel (A.IsTarget b) n y c else 0 := by
            rw [Finset.sum_comm]
          _ = ∑ c : Q₁, if F.stateMap c = d then
                ∑ y : Q₁, A.kernel.probability x y *
                  firstHitAtMass A.kernel (A.IsTarget b) n y c else 0 := by
            apply Finset.sum_congr rfl
            intro c _
            by_cases hc : F.stateMap c = d <;> simp [hc]

/-- Fixed-horizon cumulative endpoint pushforward. -/
theorem cumulativeEndpointMass_pushforward {b : O₁} (hb : b ∈ A.protectedSet)
    (N : ℕ) (x : Q₁) (d : Q₂) :
    cumulativeEndpointMass B.kernel (B.IsTarget (F.behaviorMap b)) N
        (F.stateMap x) d =
      F.fiberSum d (fun c => cumulativeEndpointMass A.kernel (A.IsTarget b) N x c) := by
  classical
  unfold cumulativeEndpointMass fiberSum
  calc
    (∑ n ∈ Finset.range (N + 1),
        firstHitAtMass B.kernel (B.IsTarget (F.behaviorMap b)) n (F.stateMap x) d) =
        ∑ n ∈ Finset.range (N + 1), ∑ c : Q₁,
          if F.stateMap c = d then
            firstHitAtMass A.kernel (A.IsTarget b) n x c else 0 := by
      apply Finset.sum_congr rfl
      intro n hn
      exact F.firstHitAtMass_pushforward hb n x d
    _ = ∑ c : Q₁, ∑ n ∈ Finset.range (N + 1),
          if F.stateMap c = d then
            firstHitAtMass A.kernel (A.IsTarget b) n x c else 0 := by
      rw [Finset.sum_comm]
    _ = ∑ c : Q₁, if F.stateMap c = d then
          ∑ n ∈ Finset.range (N + 1),
            firstHitAtMass A.kernel (A.IsTarget b) n x c else 0 := by
      apply Finset.sum_congr rfl
      intro c _
      by_cases hc : F.stateMap c = d <;> simp [hc]

/-- Fixed-horizon total hit transport, derived independently of endpoint coverage. -/
theorem hitBy_exact {b : O₁} (hb : b ∈ A.protectedSet) :
    ∀ (N : ℕ) (x : Q₁),
      hitBy B.kernel (B.IsTarget (F.behaviorMap b)) N (F.stateMap x) =
        hitBy A.kernel (A.IsTarget b) N x := by
  intro N
  induction N with
  | zero =>
      intro x
      rw [hitBy_zero, hitBy_zero]
      have ht := F.target_reflect b x hb
      by_cases hA : A.IsTarget b x
      · have hB := ht.mp hA
        simp [hA, hB]
      · have hB : ¬ B.IsTarget (F.behaviorMap b) (F.stateMap x) :=
          fun h => hA (ht.mpr h)
        simp [hA, hB]
  | succ N ih =>
      intro x
      rw [hitBy_succ_recurrence, hitBy_succ_recurrence]
      have ht := F.target_reflect b x hb
      by_cases hA : A.IsTarget b x
      · have hB := ht.mp hA
        simp [hA, hB]
      · have hB : ¬ B.IsTarget (F.behaviorMap b) (F.stateMap x) :=
          fun h => hA (ht.mpr h)
        rw [if_neg hA, if_neg hB]
        change B.kernel.apply
            (fun z => hitBy B.kernel (B.IsTarget (F.behaviorMap b)) N z)
            (F.stateMap x) = _
        rw [F.kernel_lumping]
        simp only [FiniteMarkovKernel.apply, Function.comp_apply, ih]

/-- Bounded-limit total hit equality. -/
theorem hitProbability_exact_derived {b : O₁} (hb : b ∈ A.protectedSet)
    (x : Q₁) :
    B.hit (F.behaviorMap b) (F.stateMap x) = A.hit b x := by
  unfold FixedBehaviorModel.hit
  apply tendsto_nhds_unique
    (hitBy_tendsto_hitProbability B.kernel (B.IsTarget (F.behaviorMap b)) (F.stateMap x))
  convert hitBy_tendsto_hitProbability A.kernel (A.IsTarget b) x using 1
  funext N
  exact F.hitBy_exact hb N x

/-- Bounded-limit endpoint pushforward; it is a theorem, never a morphism field. -/
theorem endpointHitProbability_pushforward_derived {b : O₁}
    (hb : b ∈ A.protectedSet) (x : Q₁) (d : Q₂) :
    B.endpointMass (F.behaviorMap b) (F.stateMap x) d =
      F.fiberSum d (fun c => A.endpointMass b x c) := by
  unfold FixedBehaviorModel.endpointMass
  apply tendsto_nhds_unique
    (cumulativeEndpointMass_tendsto B.kernel
      (B.IsTarget (F.behaviorMap b)) (F.stateMap x) d)
  have hsum : Filter.Tendsto
      (fun N => F.fiberSum d
        (fun c => cumulativeEndpointMass A.kernel (A.IsTarget b) N x c))
      Filter.atTop
      (nhds (F.fiberSum d (fun c => endpointHitProbability A.kernel (A.IsTarget b) x c))) := by
    classical
    unfold fiberSum
    exact tendsto_finset_sum Finset.univ fun c _ => by
      by_cases hc : F.stateMap c = d
      · simpa [hc] using cumulativeEndpointMass_tendsto A.kernel (A.IsTarget b) x c
      · simpa [hc] using (tendsto_const_nhds :
          Filter.Tendsto (fun _ : ℕ => (0 : ℝ)) Filter.atTop (nhds 0))
  convert hsum using 1
  funext N
  exact F.cumulativeEndpointMass_pushforward hb N x d

theorem positive_exact_derived {b : O₁} (hb : b ∈ A.protectedSet) :
    A.Positive b ↔ B.Positive (F.behaviorMap b) := by
  constructor
  · rintro ⟨x, hx, hpos⟩
    refine ⟨F.stateMap x, (F.initial_reflect b x hb).mp hx, ?_⟩
    rwa [F.hitProbability_exact_derived hb x]
  · rintro ⟨y, hy, hpos⟩
    rcases F.initial_cover b y hb hy with ⟨x, hx, rfl⟩
    refine ⟨x, hx, ?_⟩
    rwa [F.hitProbability_exact_derived hb x] at hpos

theorem almostSure_exact_derived {b : O₁} (hb : b ∈ A.protectedSet) :
    A.AlmostSure b ↔ B.AlmostSure (F.behaviorMap b) := by
  constructor
  · intro hA y hy
    rcases F.initial_cover b y hb hy with ⟨x, hx, rfl⟩
    rw [F.hitProbability_exact_derived hb x]
    exact hA x hx
  · intro hB x hx
    rw [← F.hitProbability_exact_derived hb x]
    exact hB (F.stateMap x) ((F.initial_reflect b x hb).mp hx)

theorem endpointUnique_exact_derived {b : O₁} (hb : b ∈ A.protectedSet) :
    A.EndpointUnique b ↔ B.EndpointUnique (F.behaviorMap b) := by
  constructor
  · rintro ⟨⟨x, hx, c, hc⟩, hunique⟩
    refine ⟨?_, ?_⟩
    · refine ⟨F.stateMap x, (F.initial_reflect b x hb).mp hx,
          F.stateMap c, ?_⟩
      rw [F.endpointHitProbability_pushforward_derived hb x (F.stateMap c)]
      exact fiberSum_pos_of F _ _ (A.endpointMass_nonneg b x) rfl hc
    · intro y hy d hd y' hy' e he
      rcases F.initial_cover b y hb hy with ⟨x, hx, rfl⟩
      rcases F.initial_cover b y' hb hy' with ⟨x', hx', rfl⟩
      rw [F.endpointHitProbability_pushforward_derived hb x d] at hd
      rw [F.endpointHitProbability_pushforward_derived hb x' e] at he
      rcases exists_pos_of_fiberSum_pos F d _ (A.endpointMass_nonneg b x) hd with
        ⟨c, hcd, hcpos⟩
      rcases exists_pos_of_fiberSum_pos F e _ (A.endpointMass_nonneg b x') he with
        ⟨c', hc'e, hc'pos⟩
      rw [← hcd, ← hc'e]
      exact (F.quotient_exact c c').mp
        (hunique x hx c hcpos x' hx' c' hc'pos)
  · rintro ⟨⟨y, hy, d, hd⟩, hunique⟩
    rcases F.initial_cover b y hb hy with ⟨x, hx, rfl⟩
    rw [F.endpointHitProbability_pushforward_derived hb x d] at hd
    rcases exists_pos_of_fiberSum_pos F d _ (A.endpointMass_nonneg b x) hd with
      ⟨c, hcd, hcpos⟩
    refine ⟨⟨x, hx, c, hcpos⟩, ?_⟩
    intro x₁ hx₁ c₁ hc₁ x₂ hx₂ c₂ hc₂
    apply (F.quotient_exact c₁ c₂).mpr
    apply hunique (F.stateMap x₁) ((F.initial_reflect b x₁ hb).mp hx₁)
        (F.stateMap c₁)
    · rw [F.endpointHitProbability_pushforward_derived hb x₁ (F.stateMap c₁)]
      exact fiberSum_pos_of F _ _ (A.endpointMass_nonneg b x₁) rfl hc₁
    · exact (F.initial_reflect b x₂ hb).mp hx₂
    · rw [F.endpointHitProbability_pushforward_derived hb x₂ (F.stateMap c₂)]
      exact fiberSum_pos_of F _ _ (A.endpointMass_nonneg b x₂) rfl hc₂

theorem layer0_exact {b : O₁} : A.InL0 b ↔ B.InL0 (F.behaviorMap b) := by
  constructor
  · rintro ⟨hb, hf⟩
    exact ⟨(F.protected_reflect b).mp hb, (F.fiber_nonempty_exact_derived hb).mp hf⟩
  · rintro ⟨hgb, hf⟩
    have hb := (F.protected_reflect b).mpr hgb
    exact ⟨hb, (F.fiber_nonempty_exact_derived hb).mpr hf⟩

theorem layer1_exact {b : O₁} : A.InL1 b ↔ B.InL1 (F.behaviorMap b) := by
  constructor
  · rintro ⟨h0, h1⟩
    exact ⟨F.layer0_exact.mp h0, (F.positive_exact_derived h0.1).mp h1⟩
  · rintro ⟨h0, h1⟩
    have hb := (F.layer0_exact.mpr h0).1
    exact ⟨F.layer0_exact.mpr h0, (F.positive_exact_derived hb).mpr h1⟩

theorem layer2_exact {b : O₁} : A.InL2 b ↔ B.InL2 (F.behaviorMap b) := by
  constructor
  · rintro ⟨h1, h2⟩
    exact ⟨F.layer1_exact.mp h1, (F.almostSure_exact_derived h1.1.1).mp h2⟩
  · rintro ⟨h1, h2⟩
    have h1A := F.layer1_exact.mpr h1
    exact ⟨h1A, (F.almostSure_exact_derived h1A.1.1).mpr h2⟩

theorem layer3_exact {b : O₁} : A.InL3 b ↔ B.InL3 (F.behaviorMap b) := by
  constructor
  · rintro ⟨h2, h3⟩
    exact ⟨F.layer2_exact.mp h2, (F.endpointUnique_exact_derived h2.1.1.1).mp h3⟩
  · rintro ⟨h2, h3⟩
    have h2A := F.layer2_exact.mpr h2
    exact ⟨h2A, (F.endpointUnique_exact_derived h2A.1.1.1).mpr h3⟩

end StructuralMorphism

end ObservableNormalForms.ProtectedObstructions.NonidentityExactness
/- A model of the two T6 path transformations. This section uses only Lean
   core constructs, so the logical argument rests on no Mathlib lemma even
   though the surrounding module imports Mathlib. -/

universe u

namespace ObservableNormalForms.ProtectedObstructions.Scheduler

variable {State : Type u}

/-- Front-recursive reflexive-transitive closure. -/
inductive PathClosure (step : State → State → Prop) : State → State → Prop where
  | refl (x) : PathClosure step x x
  | head {x y z} : step x y → PathClosure step y z → PathClosure step x z

namespace PathClosure

theorem tail {step : State → State → Prop} {x y z : State}
    (path : PathClosure step x y) (last : step y z) :
    PathClosure step x z := by
  induction path with
  | refl => exact .head last (.refl _)
  | head first _ ih => exact .head first (ih last)

end PathClosure

/-- A finite positive path whose last state is the first target. -/
inductive FirstHitPath (support : State → State → Prop) (target : State → Prop) :
    State → State → Prop where
  | done {x} : target x → FirstHitPath support target x x
  | step {x y c} : ¬ target x → support x y →
      FirstHitPath support target y c → FirstHitPath support target x c

abbrev PreHitStep (support : State → State → Prop) (target : State → Prop)
    (x y : State) : Prop :=
  ¬ target x ∧ support x y

abbrev PreHitReach (support : State → State → Prop) (target : State → Prop)
    (x y : State) : Prop :=
  PathClosure (PreHitStep support target) x y

/-- Exactly the two non-circular one-step scheduler obligations. -/
structure SchedulerInterface
    (support rewrite : State → State → Prop) (target : State → Prop)
    (initial : State → Prop) : Prop where
  supportSound : ∀ {x y},
    (∃ a, initial a ∧ PreHitReach support target a x) →
    support x y → rewrite x y ∨ x = y
  rewriteComplete : ∀ {x y},
    (∃ a, initial a ∧ PreHitReach support target a x) →
    rewrite x y → support x y

private theorem support_to_rewrite_aux
    {support rewrite : State → State → Prop} {target initial : State → Prop}
    (sched : SchedulerInterface support rewrite target initial)
    {root x c : State} (hroot : initial root)
    (reached : PreHitReach support target root x)
    (path : FirstHitPath support target x c) :
    PathClosure rewrite x c := by
  induction path generalizing root with
  | done _ => exact .refl _
  | @step x y _ hnot hxy _ ih =>
      have reachedY : PreHitReach support target root y :=
        PathClosure.tail reached ⟨hnot, hxy⟩
      have tail := ih hroot reachedY
      cases sched.supportSound ⟨root, hroot, reached⟩ hxy with
      | inl hrewrite => exact .head hrewrite tail
      | inr heq => cases heq; exact tail

theorem support_to_rewrite
    {support rewrite : State → State → Prop} {target initial : State → Prop}
    (sched : SchedulerInterface support rewrite target initial)
    {root c : State} (hroot : initial root)
    (path : FirstHitPath support target root c) :
    PathClosure rewrite root c :=
  support_to_rewrite_aux sched hroot (.refl root) path

private theorem rewrite_to_support_aux
    {support rewrite : State → State → Prop} {target initial : State → Prop}
    (sched : SchedulerInterface support rewrite target initial)
    (target_is_normal : ∀ {x}, target x → ∀ y, ¬ rewrite x y)
    {root x c : State} (hroot : initial root)
    (reached : PreHitReach support target root x) (hc : target c)
    (path : PathClosure rewrite x c) :
    FirstHitPath support target x c := by
  induction path generalizing root with
  | refl => exact .done hc
  | @head x y _ hrewrite _ ih =>
      have hnot : ¬ target x := by
        intro htarget
        exact target_is_normal htarget y hrewrite
      have hsupport : support x y :=
        sched.rewriteComplete ⟨root, hroot, reached⟩ hrewrite
      have reachedY : PreHitReach support target root y :=
        PathClosure.tail reached ⟨hnot, hsupport⟩
      exact .step hnot hsupport (ih hroot reachedY hc)

theorem rewrite_to_support
    {support rewrite : State → State → Prop} {target initial : State → Prop}
    (sched : SchedulerInterface support rewrite target initial)
    (target_is_normal : ∀ {x}, target x → ∀ y, ¬ rewrite x y)
    {root c : State} (hroot : initial root) (hc : target c)
    (path : PathClosure rewrite root c) :
    FirstHitPath support target root c :=
  rewrite_to_support_aux sched target_is_normal hroot (.refl root) hc path

end ObservableNormalForms.ProtectedObstructions.Scheduler

/-! # Public T0--T8 adapters and structural composition -/

namespace ObservableNormalForms.ProtectedObstructions.PublicAdapter

open scoped BigOperators
open Relation
open ObservableNormalForms
open ObservableNormalForms.FiniteMarkovKernel
open ObservableNormalForms.MechanismVariants
open ObservableNormalForms.ProtectedObstructions.NonidentityExactness

universe uPA vPA wPA xPA yPA zPA

noncomputable section

abbrev CoreModel := FixedBehaviorModel

private def equalitySetoid (alpha : Type*) : Setoid alpha where
  r := Eq
  iseqv := ⟨fun _ => rfl, fun h => h.symm, fun h g => h.trans g⟩

/-! ## Frozen public profile and T0/T1 -/

/-- The static/rewrite data omitted from the fixed-behavior stochastic core.
`target_iff` only identifies the canonical Boolean target with the
approved consistent observation fiber. -/
structure Profile (O : Type uPA) (Q : Type vPA)
    [Fintype O] [DecidableEq O] [Fintype Q] [DecidableEq Q] where
  core : CoreModel O Q
  protected_nonempty : core.protectedSet.Nonempty
  observe : Q → O
  consistent : Q → Prop
  target_iff : ∀ b q, core.IsTarget b q ↔ consistent q ∧ observe q = b
  initial_observe : ∀ {b q}, b ∈ core.protectedSet →
    q ∈ core.initial b → observe q = b
  rewrite : Q → Q → Prop
  rewrite_observe : ObservableNormalForms.ObservationPreserving rewrite observe
  normal_iff_consistent : ∀ q,
    ObservableNormalForms.IsNormalForm rewrite q ↔ consistent q

variable {O : Type uPA} {Q : Type vPA}
  [Fintype O] [DecidableEq O] [Fintype Q] [DecidableEq Q]

namespace Profile

def fiber (M : Profile O Q) (b : O) : Finset Q := by
  classical
  exact Finset.univ.filter (M.core.IsTarget b)

def layer0 (M : Profile O Q) : Set O := {b | M.core.InL0 b}
def layer1 (M : Profile O Q) : Set O := {b | M.core.InL1 b}
def layer2 (M : Profile O Q) : Set O := {b | M.core.InL2 b}
def layer3 (M : Profile O Q) : Set O := {b | M.core.InL3 b}

def cut0 (M : Profile O Q) : Set O := (M.core.protectedSet : Set O) \ M.layer0
def cut1 (M : Profile O Q) : Set O := (M.core.protectedSet : Set O) \ M.layer1
def cut2 (M : Profile O Q) : Set O := (M.core.protectedSet : Set O) \ M.layer2
def cut3 (M : Profile O Q) : Set O := (M.core.protectedSet : Set O) \ M.layer3

def delta0 (M : Profile O Q) : Set O := (M.core.protectedSet : Set O) \ M.layer0
def delta1 (M : Profile O Q) : Set O := M.layer0 \ M.layer1
def delta2 (M : Profile O Q) : Set O := M.layer1 \ M.layer2
def delta3 (M : Profile O Q) : Set O := M.layer2 \ M.layer3

theorem mem_layer0_iff (M : Profile O Q) {b : O} :
    b ∈ M.layer0 ↔ M.core.InL0 b := by
  rfl

theorem mem_layer1_iff (M : Profile O Q) {b : O} :
    b ∈ M.layer1 ↔ M.core.InL1 b := by
  rfl

theorem mem_layer2_iff (M : Profile O Q) {b : O} :
    b ∈ M.layer2 ↔ M.core.InL2 b := by
  rfl

theorem mem_layer3_iff (M : Profile O Q) {b : O} :
    b ∈ M.layer3 ↔ M.core.InL3 b := by
  rfl

theorem layer1_subset_layer0 (M : Profile O Q) : M.layer1 ⊆ M.layer0 := by
  intro b hb
  exact hb.1

theorem layer2_subset_layer1 (M : Profile O Q) : M.layer2 ⊆ M.layer1 := by
  intro b hb
  exact hb.1

theorem layer3_subset_layer2 (M : Profile O Q) : M.layer3 ⊆ M.layer2 := by
  intro b hb
  exact hb.1

/-- T0: the four obstruction strata and the final survivor exhaust the
protected carrier. -/
theorem t0_partition (M : Profile O Q) :
    (M.core.protectedSet : Set O) =
      M.delta0 ∪ M.delta1 ∪ M.delta2 ∪ M.delta3 ∪ M.layer3 := by
  classical
  ext b
  constructor
  · intro hp
    by_cases h0 : b ∈ M.layer0
    · by_cases h1 : b ∈ M.layer1
      · by_cases h2 : b ∈ M.layer2
        · by_cases h3 : b ∈ M.layer3 <;>
            simp [delta0, delta1, delta2, delta3, hp, h0, h1, h2, h3]
        · simp [delta0, delta1, delta2, delta3, hp, h0, h1, h2]
      · simp [delta0, delta1, delta2, delta3, hp, h0, h1]
    · simp [delta0, delta1, delta2, delta3, hp, h0]
  · intro h
    rcases h with h | h
    · rcases h with h | h
      · rcases h with h | h
        · rcases h with h | h
          · exact h.1
          · exact h.1.1
        · exact h.1.1.1
      · exact h.1.1.1.1
    · exact h.1.1.1.1

/-- T0: no behavior belongs to two distinct obstruction/survivor strata. -/
theorem t0_pairwise_disjoint (M : Profile O Q) :
    Disjoint M.delta0 M.delta1 ∧ Disjoint M.delta0 M.delta2 ∧
    Disjoint M.delta0 M.delta3 ∧ Disjoint M.delta0 M.layer3 ∧
    Disjoint M.delta1 M.delta2 ∧ Disjoint M.delta1 M.delta3 ∧
    Disjoint M.delta1 M.layer3 ∧ Disjoint M.delta2 M.delta3 ∧
    Disjoint M.delta2 M.layer3 ∧ Disjoint M.delta3 M.layer3 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
    apply Set.disjoint_left.mpr <;> intro b hleft hright
  · exact hleft.2 hright.1
  · exact hleft.2 (M.layer1_subset_layer0 hright.1)
  · exact hleft.2
      (M.layer1_subset_layer0 (M.layer2_subset_layer1
        hright.1))
  · exact hleft.2
      (M.layer1_subset_layer0 (M.layer2_subset_layer1
        (M.layer3_subset_layer2 hright)))
  · exact hleft.2 hright.1
  · exact hleft.2 (M.layer2_subset_layer1 hright.1)
  · exact hleft.2
      (M.layer2_subset_layer1 (M.layer3_subset_layer2 hright))
  · exact hleft.2 hright.1
  · exact hleft.2 (M.layer3_subset_layer2 hright)
  · exact hleft.2 hright

/-- T1, static half. -/
theorem t1_fiber_empty_iff (M : Profile O Q) (b : O) :
    M.fiber b = ∅ ↔ ¬ ∃ c, M.consistent c ∧ M.observe c = b := by
  classical
  constructor
  · intro hempty ⟨c, hc, hobs⟩
    have : c ∈ M.fiber b := by
      simp [fiber, (M.target_iff b c), hc, hobs]
    simpa [hempty] using this
  · intro hnone
    apply Finset.eq_empty_iff_forall_notMem.mpr
    intro c hc
    apply hnone
    have ht : M.core.IsTarget b c := by simpa [fiber] using hc
    exact ⟨c, (M.target_iff b c).mp ht⟩

/-- T1, rewrite half: an observation-preserving path cannot end at a normal
endpoint in an empty protected fiber. -/
theorem t1_no_normal_endpoint_of_empty_fiber (M : Profile O Q)
    {b : O} (hempty : M.fiber b = ∅) {q c : Q}
    (hqb : M.observe q = b) (hpath : ReflTransGen M.rewrite q c)
    (hnormal : ObservableNormalForms.IsNormalForm M.rewrite c) : False := by
  have hobs : M.observe q = M.observe c :=
    ObservableNormalForms.observation_eq_of_reflTransGen
      M.rewrite_observe hpath
  have hc : M.consistent c := (M.normal_iff_consistent c).mp hnormal
  exact (M.t1_fiber_empty_iff b).mp hempty ⟨c, hc, hobs ▸ hqb⟩

/-! ## Canonical-law bridges and public T2/T3/T4/T5 -/

theorem firstHitAt_target_endpoint (target : Q → Prop) [DecidablePred target] :
    ∀ {n : ℕ} {q : Q} {p : Fin n → Q},
      FirstHitAt target n q p → target (pathEndpoint n q p) := by
  intro n
  induction n with
  | zero =>
      intro q p h
      simpa [FirstHitAt, pathEndpoint] using h
  | succ n ih =>
      intro q p h
      have htail := ih h.2
      rw [← pathEndpoint_cons n q (p 0) (Fin.tail p)] at htail
      simpa only [Fin.cons_self_tail] using htail

/-- The total bounded limit is the finite sum of endpoint bounded limits. -/
theorem hit_eq_sum_endpoint (M : Profile O Q) (b : O) (q : Q) :
    M.core.hit b q = ∑ c : Q, M.core.endpointMass b q c := by
  unfold FixedBehaviorModel.hit FixedBehaviorModel.endpointMass
  apply tendsto_nhds_unique
    (hitBy_tendsto_hitProbability M.core.kernel (M.core.IsTarget b) q)
  have hsum : Filter.Tendsto
      (fun N => ∑ c : Q,
        cumulativeEndpointMass M.core.kernel (M.core.IsTarget b) N q c)
      Filter.atTop
      (nhds (∑ c : Q,
        endpointHitProbability M.core.kernel (M.core.IsTarget b) q c)) :=
    tendsto_finset_sum Finset.univ fun c _ =>
      cumulativeEndpointMass_tendsto M.core.kernel (M.core.IsTarget b) q c
  simpa [hitBy] using hsum

theorem hit_pos_iff_endpoint (M : Profile O Q) (b : O) (q : Q) :
    0 < M.core.hit b q ↔ ∃ c, 0 < M.core.endpointMass b q c := by
  classical
  rw [M.hit_eq_sum_endpoint b q,
    Finset.sum_pos_iff_of_nonneg
      (fun c _ => M.core.endpointMass_nonneg b q c)]
  simp

def PositiveFirstHitPath (M : Profile O Q) (q : Q) (b : O) (c : Q) : Prop :=
  ∃ (n : ℕ) (p : Fin n → Q),
    FirstHitAt (M.core.IsTarget b) n q p ∧
    pathEndpoint n q p = c ∧
    SupportedFinPath M.core.kernel n q p

theorem endpoint_pos_iff_path (M : Profile O Q) (q : Q) (b : O) (c : Q) :
    0 < M.core.endpointMass b q c ↔ M.PositiveFirstHitPath q b c := by
  exact endpointHitProbability_pos_iff M.core.kernel (M.core.IsTarget b) q c

theorem endpoint_pos_target (M : Profile O Q) {q : Q} {b : O} {c : Q}
    (hpos : 0 < M.core.endpointMass b q c) : M.core.IsTarget b c := by
  rcases (M.endpoint_pos_iff_path q b c).mp hpos with ⟨n, p, hfirst, hend, _⟩
  rw [← hend]
  exact firstHitAt_target_endpoint (M.core.IsTarget b) hfirst

/-- T2, source-quantified profile wrapper, derived from the canonical
endpoint law. -/
theorem t2_positive_reach_iff (M : Profile O Q) {b : O}
    (hb : b ∈ M.core.protectedSet) :
    b ∈ M.layer1 ↔
      ∃ q ∈ M.core.initial b, ∃ c, M.PositiveFirstHitPath q b c := by
  constructor
  · intro h1
    rcases h1.2 with ⟨q, hq, hhit⟩
    rcases (M.hit_pos_iff_endpoint b q).mp hhit with ⟨c, hc⟩
    exact ⟨q, hq, c, (M.endpoint_pos_iff_path q b c).mp hc⟩
  · rintro ⟨q, hq, c, hpath⟩
    have hc : 0 < M.core.endpointMass b q c :=
      (M.endpoint_pos_iff_path q b c).mpr hpath
    have htarget := M.endpoint_pos_target hc
    exact ⟨⟨hb, ⟨c, htarget⟩⟩,
      ⟨q, hq, (M.hit_pos_iff_endpoint b q).mpr ⟨c, hc⟩⟩⟩

abbrev NoClosedTrap (M : Profile O Q) (q : Q) (b : O) : Prop :=
  ¬ Nonempty (ReachableClosedTrap M.core.kernel (M.core.IsTarget b) q)

/-- T3, source-quantified profile wrapper, derived from the canonical closed
trap theorem. -/
theorem t3_almostSure_iff_noClosedTrap (M : Profile O Q) {b : O} :
    b ∈ M.layer2 ↔
      b ∈ M.layer1 ∧ ∀ q ∈ M.core.initial b, M.NoClosedTrap q b := by
  constructor
  · rintro ⟨h1, hall⟩
    exact ⟨h1, fun q hq =>
      (hitProbability_eq_one_iff_no_reachableClosedTrap
        M.core.kernel (M.core.IsTarget b) q).mp (hall q hq)⟩
  · rintro ⟨h1, hall⟩
    exact ⟨h1, fun q hq =>
      (hitProbability_eq_one_iff_no_reachableClosedTrap
        M.core.kernel (M.core.IsTarget b) q).mpr (hall q hq)⟩

/-- Structural endpoint-union certificate stated only with positive
first-hit paths.  It contains no endpoint-law, uniqueness, or layer
conclusion premise. -/
def EndpointUnionCert (M : Profile O Q) (b : O) : Prop :=
  ∃ c0,
    (∃ q ∈ M.core.initial b, M.PositiveFirstHitPath q b c0) ∧
    ∀ q ∈ M.core.initial b, ∀ c,
      M.PositiveFirstHitPath q b c → M.core.silent.r c c0

theorem endpointUnique_iff_endpointUnionCert (M : Profile O Q) (b : O) :
    M.core.EndpointUnique b ↔ M.EndpointUnionCert b := by
  constructor
  · rintro ⟨⟨q0, hq0, c0, hc0⟩, hall⟩
    exact ⟨c0, ⟨q0, hq0, (M.endpoint_pos_iff_path q0 b c0).mp hc0⟩,
      fun q hq c hc => hall q hq c
        ((M.endpoint_pos_iff_path q b c).mpr hc) q0 hq0 c0 hc0⟩
  · rintro ⟨c0, ⟨q0, hq0, hc0⟩, hall⟩
    refine ⟨⟨q0, hq0, c0, (M.endpoint_pos_iff_path q0 b c0).mpr hc0⟩, ?_⟩
    intro q hq c hc q' hq' d hd
    exact M.core.silent.iseqv.trans
      (hall q hq c ((M.endpoint_pos_iff_path q b c).mp hc))
      (M.core.silent.iseqv.symm
        (hall q' hq' d ((M.endpoint_pos_iff_path q' b d).mp hd)))

theorem endpointUnionCert_iff_path (M : Profile O Q) (b : O) :
    M.EndpointUnionCert b ↔
      ∃ c0,
        (∃ q ∈ M.core.initial b, M.PositiveFirstHitPath q b c0) ∧
        ∀ q ∈ M.core.initial b, ∀ c,
          M.PositiveFirstHitPath q b c → M.core.silent.r c c0 := by
  rfl

/-- T4. -/
theorem t4_selection_iff_endpointUnionCert (M : Profile O Q) {b : O} :
    b ∈ M.layer3 ↔ b ∈ M.layer2 ∧ M.EndpointUnionCert b := by
  change M.core.InL3 b ↔ M.core.InL2 b ∧ M.EndpointUnionCert b
  rw [FixedBehaviorModel.InL3, M.endpointUnique_iff_endpointUnionCert]

/-- T5: observable determination collapses L2 to L3, but makes no claim
about L0/L1. -/
theorem t5_observableDetermination_collapse (M : Profile O Q)
    (hdet : ∀ {b c d}, M.core.IsTarget b c → M.core.IsTarget b d →
      M.core.silent.r c d) :
    M.layer2 = M.layer3 := by
  apply Set.Subset.antisymm
  · intro b hb2
    have h1 : b ∈ M.layer1 := hb2.1
    rcases h1.2 with ⟨q0, hq0, hhit⟩
    rcases (M.hit_pos_iff_endpoint b q0).mp hhit with ⟨c0, hc0⟩
    refine ⟨hb2, ⟨⟨q0, hq0, c0, hc0⟩, ?_⟩⟩
    intro q hq c hc q' hq' d hd
    exact hdet (M.endpoint_pos_target hc) (M.endpoint_pos_target hd)
  · exact M.layer3_subset_layer2

/-! ## T6: public scheduler adapter with conclusion fields absent -/

abbrev SupportStep (M : Profile O Q) (q r : Q) : Prop :=
  0 < M.core.kernel.probability q r

/-- Only the two approved one-step obligations are fields. -/
structure SchedulerLaws (M : Profile O Q) (b : O) : Prop where
  supportSound : ∀ {q r},
    (∃ a, a ∈ M.core.initial b ∧
      Scheduler.PreHitReach M.SupportStep (M.core.IsTarget b) a q) →
    M.SupportStep q r → M.rewrite q r ∨ q = r
  rewriteComplete : ∀ {q r},
    (∃ a, a ∈ M.core.initial b ∧
      Scheduler.PreHitReach M.SupportStep (M.core.IsTarget b) a q) →
    M.rewrite q r → M.SupportStep q r

def SchedulerLaws.toInterface (M : Profile O Q) {b : O}
    (laws : M.SchedulerLaws b) :
    Scheduler.SchedulerInterface M.SupportStep M.rewrite
      (M.core.IsTarget b) (fun q => q ∈ M.core.initial b) where
  supportSound := laws.supportSound
  rewriteComplete := laws.rewriteComplete

theorem target_is_normal (M : Profile O Q) {b : O} {q : Q}
    (hq : M.core.IsTarget b q) : ∀ r, ¬ M.rewrite q r := by
  have hc : M.consistent q := (M.target_iff b q).mp hq |>.1
  exact (M.normal_iff_consistent q).mpr hc

/-- T6, support-to-rewrite direction, derived from `supportSound`. -/
theorem t6_support_to_rewrite (M : Profile O Q) {b q c}
    (laws : M.SchedulerLaws b) (hq : q ∈ M.core.initial b)
    (path : Scheduler.FirstHitPath M.SupportStep (M.core.IsTarget b) q c) :
    Scheduler.PathClosure M.rewrite q c :=
  Scheduler.support_to_rewrite (laws.toInterface M) hq path

/-- T6, rewrite-to-support direction, derived from `rewriteComplete` and
target normality. -/
theorem t6_rewrite_to_support (M : Profile O Q) {b q c}
    (laws : M.SchedulerLaws b) (hq : q ∈ M.core.initial b)
    (hc : M.core.IsTarget b c)
    (path : Scheduler.PathClosure M.rewrite q c) :
    Scheduler.FirstHitPath M.SupportStep (M.core.IsTarget b) q c :=
  Scheduler.rewrite_to_support (laws.toInterface M)
    (fun ht => M.target_is_normal ht) hq hc path

/-! ## T7: exact `RealizedBehavior`/`BehaviorCut` recovery -/

structure SupportedTrace (M : Profile O Q) where
  n : ℕ
  start : Q
  path : Fin n → Q
  start_active : start ∈ M.core.initial (M.observe start)
  supported : SupportedFinPath M.core.kernel n start path

theorem supportedFinPath_reflTransGen (M : Profile O Q) :
    ∀ {n : ℕ} {q : Q} {p : Fin n → Q},
      SupportedFinPath M.core.kernel n q p →
      ReflTransGen M.SupportStep q (pathEndpoint n q p) := by
  intro n
  induction n with
  | zero =>
      intro q p hs
      exact .refl
  | succ n ih =>
      intro q p hs
      have htail := ih hs.2
      have hchain : ReflTransGen M.SupportStep q
          (pathEndpoint n (p 0) (Fin.tail p)) :=
        Relation.ReflTransGen.head hs.1 htail
      rw [← pathEndpoint_cons n q (p 0) (Fin.tail p)] at hchain
      simpa only [Fin.cons_self_tail] using hchain

def supportedDynamics (M : Profile O Q) : Dynamics Unit where
  State := fun _ => Q
  Step := fun _ => M.SupportStep
  Initial := fun _ q => q ∈ M.core.initial (M.observe q)
  Trace := fun _ => M.SupportedTrace
  start := fun _ t => t.start
  finish := fun _ t => pathEndpoint t.n t.start t.path
  startsInitial := fun _ t => t.start_active
  traceReachable := fun _ t => M.supportedFinPath_reflTransGen t.supported

def traceObserve (M : Profile O Q) : Observation M.supportedDynamics O :=
  fun _ t => M.observe t.start

def Admissible0 (M : Profile O Q) (_ : Unit) (t : M.SupportedTrace) : Prop :=
  M.observe t.start ∈ M.core.protectedSet ∧
    M.core.FiberNonempty (M.observe t.start)

def Admissible1 (M : Profile O Q) (_ : Unit) (t : M.SupportedTrace) : Prop :=
  M.Admissible0 () t ∧
    FirstHitAt (M.core.IsTarget (M.observe t.start)) t.n t.start t.path

def Admissible2 (M : Profile O Q) (_ : Unit) (t : M.SupportedTrace) : Prop :=
  M.Admissible1 () t ∧
    ∀ q ∈ M.core.initial (M.observe t.start),
      M.NoClosedTrap q (M.observe t.start)

def Admissible3 (M : Profile O Q) (_ : Unit) (t : M.SupportedTrace) : Prop :=
  M.Admissible2 () t ∧ M.EndpointUnionCert (M.observe t.start)

private def traceOfPath (M : Profile O Q) {b q : _}
    (hq : q ∈ M.core.initial b) (hobs : M.observe q = b)
    (n : ℕ) (p : Fin n → Q) (hs : SupportedFinPath M.core.kernel n q p) :
    M.SupportedTrace where
  n := n
  start := q
  path := p
  start_active := by simpa [hobs] using hq
  supported := hs

theorem t7_realizedBehavior0 (M : Profile O Q) :
    RealizedBehavior M.supportedDynamics M.Admissible0 M.traceObserve () =
      M.layer0 := by
  ext b
  constructor
  · rintro ⟨t, ht, rfl⟩
    exact ht
  · intro hb
    have hb0 : b ∈ M.layer0 := hb
    have hbP : b ∈ M.core.protectedSet := hb0.1
    rcases M.core.initial_nonempty b hbP with ⟨q, hq⟩
    have hobs := M.initial_observe hbP hq
    let p : Fin 0 → Q := fun i => Fin.elim0 i
    refine ⟨M.traceOfPath hq hobs 0 p (by trivial), ?_, ?_⟩
    · simpa [Admissible0, traceOfPath, hobs] using hb0
    · exact hobs

theorem t7_realizedBehavior1 (M : Profile O Q) :
    RealizedBehavior M.supportedDynamics M.Admissible1 M.traceObserve () =
      M.layer1 := by
  ext b
  constructor
  · rintro ⟨t, ht, hobs⟩
    change M.observe t.start = b at hobs
    subst b
    have hbP : M.observe t.start ∈ M.core.protectedSet := ht.1.1
    apply (M.t2_positive_reach_iff hbP).mpr
    refine ⟨t.start, t.start_active, pathEndpoint t.n t.start t.path,
      t.n, t.path, ht.2, rfl, t.supported⟩
  · intro hb
    have hb1 : b ∈ M.layer1 := hb
    have hbP : b ∈ M.core.protectedSet := hb1.1.1
    rcases (M.t2_positive_reach_iff hbP).mp hb1 with
      ⟨q, hq, c, n, p, hfirst, hend, hs⟩
    have hobs := M.initial_observe hbP hq
    refine ⟨M.traceOfPath hq hobs n p hs, ?_, ?_⟩
    · refine ⟨by simpa [Admissible0, traceOfPath, hobs] using hb1.1, ?_⟩
      simpa only [traceOfPath, hobs] using hfirst
    · exact hobs

theorem t7_realizedBehavior2 (M : Profile O Q) :
    RealizedBehavior M.supportedDynamics M.Admissible2 M.traceObserve () =
      M.layer2 := by
  ext b
  constructor
  · rintro ⟨t, ht, hobs⟩
    change M.observe t.start = b at hobs
    subst b
    apply M.t3_almostSure_iff_noClosedTrap.mpr
    constructor
    · rw [← M.t7_realizedBehavior1]
      exact ⟨t, ht.1, rfl⟩
    · exact ht.2
  · intro hb
    have hb2 : b ∈ M.layer2 := hb
    have h1 : b ∈ M.layer1 := hb2.1
    have htraps := (M.t3_almostSure_iff_noClosedTrap.mp hb2).2
    have hr1 : b ∈ RealizedBehavior M.supportedDynamics M.Admissible1
        M.traceObserve () := by
      rw [M.t7_realizedBehavior1]
      exact h1
    rcases hr1 with ⟨t, ht, hobs⟩
    change M.observe t.start = b at hobs
    refine ⟨t, ⟨ht, ?_⟩, hobs⟩
    intro q hq
    rw [hobs] at hq ⊢
    exact htraps q hq

theorem t7_realizedBehavior3 (M : Profile O Q) :
    RealizedBehavior M.supportedDynamics M.Admissible3 M.traceObserve () =
      M.layer3 := by
  ext b
  constructor
  · rintro ⟨t, ht, hobs⟩
    change M.observe t.start = b at hobs
    subst b
    apply M.t4_selection_iff_endpointUnionCert.mpr
    constructor
    · rw [← M.t7_realizedBehavior2]
      exact ⟨t, ht.1, rfl⟩
    · exact ht.2
  · intro hb
    have hb3 : b ∈ M.layer3 := hb
    have hb2 := (M.t4_selection_iff_endpointUnionCert.mp hb3).1
    have hcert := (M.t4_selection_iff_endpointUnionCert.mp hb3).2
    have hr2 : b ∈ RealizedBehavior M.supportedDynamics M.Admissible2
        M.traceObserve () := by
      rw [M.t7_realizedBehavior2]
      exact hb2
    rcases hr2 with ⟨t, ht, hobs⟩
    change M.observe t.start = b at hobs
    exact ⟨t, ⟨ht, by simpa [hobs] using hcert⟩, hobs⟩

theorem t7_behaviorCut0 (M : Profile O Q) :
    BehaviorCut M.supportedDynamics (M.core.protectedSet : Set O)
      M.Admissible0 M.traceObserve () = M.cut0 := by
  simp [BehaviorCut, cut0, M.t7_realizedBehavior0]

theorem t7_behaviorCut1 (M : Profile O Q) :
    BehaviorCut M.supportedDynamics (M.core.protectedSet : Set O)
      M.Admissible1 M.traceObserve () = M.cut1 := by
  simp [BehaviorCut, cut1, M.t7_realizedBehavior1]

theorem t7_behaviorCut2 (M : Profile O Q) :
    BehaviorCut M.supportedDynamics (M.core.protectedSet : Set O)
      M.Admissible2 M.traceObserve () = M.cut2 := by
  simp [BehaviorCut, cut2, M.t7_realizedBehavior2]

theorem t7_behaviorCut3 (M : Profile O Q) :
    BehaviorCut M.supportedDynamics (M.core.protectedSet : Set O)
      M.Admissible3 M.traceObserve () = M.cut3 := by
  simp [BehaviorCut, cut3, M.t7_realizedBehavior3]

end Profile

/-! ## T8: four-coordinate product preorder -/

def ProfileLE {Variant : Type wPA} {Obs : Type uPA}
    (cutAt : Fin 4 → Variant → Set Obs) (i j : Variant) : Prop :=
  ∀ k, BehaviorLE (cutAt k) i j

def ProfileLT {Variant : Type wPA} {Obs : Type uPA}
    (cutAt : Fin 4 → Variant → Set Obs) (i j : Variant) : Prop :=
  ProfileLE cutAt i j ∧ ∃ k, BehaviorLT (cutAt k) i j

theorem t8_coordinate_strict_iff {Variant : Type wPA} {Obs : Type uPA}
    (cutAt : Fin 4 → Variant → Set Obs) (k : Fin 4) (i j : Variant) :
    BehaviorLT (cutAt k) i j ↔
      cutAt k i ⊆ cutAt k j ∧ ¬ cutAt k j ⊆ cutAt k i := by
  rfl

theorem profileLE_refl {Variant : Type wPA} {Obs : Type uPA}
    (cutAt : Fin 4 → Variant → Set Obs) (i : Variant) :
    ProfileLE cutAt i i := by
  intro k q hq
  exact hq

theorem profileLE_trans {Variant : Type wPA} {Obs : Type uPA}
    (cutAt : Fin 4 → Variant → Set Obs) {i j k : Variant} :
    ProfileLE cutAt i j → ProfileLE cutAt j k → ProfileLE cutAt i k := by
  intro hij hjk n q hq
  exact hjk n (hij n hq)

/-! ## Identity and composition for the current structural morphism -/

namespace CoreStructuralMorphism

variable {O1 : Type uPA} {Q1 : Type vPA} {O2 : Type wPA} {Q2 : Type xPA}
  {O3 : Type yPA} {Q3 : Type zPA}
  [Fintype O1] [DecidableEq O1] [Fintype Q1] [DecidableEq Q1]
  [Fintype O2] [DecidableEq O2] [Fintype Q2] [DecidableEq Q2]
  [Fintype O3] [DecidableEq O3] [Fintype Q3] [DecidableEq Q3]

def id (A : CoreModel O1 Q1) : StructuralMorphism A A where
  stateMap := fun q => q
  behaviorMap := fun b => b
  protected_reflect := fun _ => Iff.rfl
  protected_cover := fun d hd => ⟨d, hd, rfl⟩
  initial_reflect := fun _ _ _ => Iff.rfl
  initial_cover := fun _ y _ hy => ⟨y, hy, rfl⟩
  kernel_lumping := fun _ _ => rfl
  target_reflect := fun _ _ _ => Iff.rfl
  target_cover := fun _ y _ hy => ⟨y, hy, rfl⟩
  quotient_exact := fun _ _ => Iff.rfl

def comp {A : CoreModel O1 Q1} {B : CoreModel O2 Q2}
    {C : CoreModel O3 Q3}
    (G : StructuralMorphism B C) (F : StructuralMorphism A B) :
    StructuralMorphism A C where
  stateMap := G.stateMap ∘ F.stateMap
  behaviorMap := G.behaviorMap ∘ F.behaviorMap
  protected_reflect := fun b => (F.protected_reflect b).trans
    (G.protected_reflect (F.behaviorMap b))
  protected_cover := by
    intro d hd
    rcases G.protected_cover d hd with ⟨c, hc, rfl⟩
    rcases F.protected_cover c hc with ⟨b, hb, rfl⟩
    exact ⟨b, hb, rfl⟩
  initial_reflect := by
    intro b q hb
    exact (F.initial_reflect b q hb).trans
      (G.initial_reflect (F.behaviorMap b) (F.stateMap q)
        ((F.protected_reflect b).mp hb))
  initial_cover := by
    intro b r hb hr
    have hFb := (F.protected_reflect b).mp hb
    rcases G.initial_cover (F.behaviorMap b) r hFb hr with ⟨q2, hq2, rfl⟩
    rcases F.initial_cover b q2 hb hq2 with ⟨q1, hq1, rfl⟩
    exact ⟨q1, hq1, rfl⟩
  kernel_lumping := by
    intro q phi
    change C.kernel.apply phi (G.stateMap (F.stateMap q)) =
      A.kernel.apply (phi ∘ G.stateMap ∘ F.stateMap) q
    rw [G.kernel_lumping, F.kernel_lumping]
    rfl
  target_reflect := by
    intro b q hb
    exact (F.target_reflect b q hb).trans
      (G.target_reflect (F.behaviorMap b) (F.stateMap q)
        ((F.protected_reflect b).mp hb))
  target_cover := by
    intro b r hb hr
    have hFb := (F.protected_reflect b).mp hb
    rcases G.target_cover (F.behaviorMap b) r hFb hr with ⟨q2, hq2, rfl⟩
    rcases F.target_cover b q2 hb hq2 with ⟨q1, hq1, rfl⟩
    exact ⟨q1, hq1, rfl⟩
  quotient_exact := fun c d => (F.quotient_exact c d).trans
    (G.quotient_exact (F.stateMap c) (F.stateMap d))

end CoreStructuralMorphism

end

/-! ## Common-carrier generic completion -/

namespace GenericCompletion

open scoped BigOperators
open ObservableNormalForms
open ObservableNormalForms.MechanismVariants
open ObservableNormalForms.ProtectedObstructions.NonidentityExactness

universe uGC vGC wGC xGC yGC zGC

noncomputable section

variable {O : Type uGC} {Q : Type vGC}
  [Fintype O] [DecidableEq O] [Fintype Q] [DecidableEq Q]

namespace Profile

theorem cut0_subset_cut1 (M : Profile O Q) : M.cut0 ⊆ M.cut1 := by
  intro b hb
  exact ⟨hb.1, fun h1 => hb.2 (M.layer1_subset_layer0 h1)⟩

theorem cut1_subset_cut2 (M : Profile O Q) : M.cut1 ⊆ M.cut2 := by
  intro b hb
  exact ⟨hb.1, fun h2 => hb.2 (M.layer2_subset_layer1 h2)⟩

theorem cut2_subset_cut3 (M : Profile O Q) : M.cut2 ⊆ M.cut3 := by
  intro b hb
  exact ⟨hb.1, fun h3 => hb.2 (M.layer3_subset_layer2 h3)⟩

theorem delta0_eq_cut0 (M : Profile O Q) : M.delta0 = M.cut0 := rfl

theorem delta1_eq_cut1_diff_cut0 (M : Profile O Q) :
    M.delta1 = M.cut1 \ M.cut0 := by
  ext b
  constructor
  · rintro ⟨h0, hn1⟩
    exact ⟨⟨h0.1, hn1⟩, fun hk0 => hk0.2 h0⟩
  · rintro ⟨⟨hp, hn1⟩, hnk0⟩
    exact ⟨by
      by_contra hn0
      exact hnk0 ⟨hp, hn0⟩, hn1⟩

theorem delta2_eq_cut2_diff_cut1 (M : Profile O Q) :
    M.delta2 = M.cut2 \ M.cut1 := by
  ext b
  constructor
  · rintro ⟨h1, hn2⟩
    exact ⟨⟨h1.1.1, hn2⟩, fun hk1 => hk1.2 h1⟩
  · rintro ⟨⟨hp, hn2⟩, hnk1⟩
    exact ⟨by
      by_contra hn1
      exact hnk1 ⟨hp, hn1⟩, hn2⟩

theorem delta3_eq_cut3_diff_cut2 (M : Profile O Q) :
    M.delta3 = M.cut3 \ M.cut2 := by
  ext b
  constructor
  · rintro ⟨h2, hn3⟩
    exact ⟨⟨h2.1.1.1, hn3⟩, fun hk2 => hk2.2 h2⟩
  · rintro ⟨⟨hp, hn3⟩, hnk2⟩
    exact ⟨by
      by_contra hn2
      exact hnk2 ⟨hp, hn2⟩, hn3⟩

theorem cut1_eq_cut0_union_delta1 (M : Profile O Q) :
    M.cut1 = M.cut0 ∪ M.delta1 := by
  rw [delta1_eq_cut1_diff_cut0 M]
  exact (Set.union_diff_cancel (cut0_subset_cut1 M)).symm

theorem cut2_eq_cut1_union_delta2 (M : Profile O Q) :
    M.cut2 = M.cut1 ∪ M.delta2 := by
  rw [delta2_eq_cut2_diff_cut1 M]
  exact (Set.union_diff_cancel (cut1_subset_cut2 M)).symm

theorem cut3_eq_cut2_union_delta3 (M : Profile O Q) :
    M.cut3 = M.cut2 ∪ M.delta3 := by
  rw [delta3_eq_cut3_diff_cut2 M]
  exact (Set.union_diff_cancel (cut2_subset_cut3 M)).symm

end Profile

/-! ## F1 common-carrier family and its actual T7 cut rows -/

structure CommonProfileData (O : Type uGC) (Q : Type vGC) (P : Finset O)
    [Fintype O] [DecidableEq O] [Fintype Q] [DecidableEq Q] where
  initial : O → Finset Q
  initial_nonempty : ∀ b ∈ P, (initial b).Nonempty
  target : O → Q → Bool
  kernel : FiniteMarkovKernel Q
  silent : Setoid Q
  observe : Q → O
  consistent : Q → Prop
  target_iff : ∀ b q, target b q = true ↔ consistent q ∧ observe q = b
  initial_observe : ∀ {b q}, b ∈ P → q ∈ initial b → observe q = b
  rewrite : Q → Q → Prop
  rewrite_observe : ObservationPreserving rewrite observe
  normal_iff_consistent : ∀ q, IsNormalForm rewrite q ↔ consistent q

namespace CommonProfileData

def toProfile {P : Finset O} (D : CommonProfileData O Q P)
    (hP : P.Nonempty) : Profile O Q where
  core := {
    protectedSet := P
    initial := D.initial
    initial_nonempty := D.initial_nonempty
    target := D.target
    kernel := D.kernel
    silent := D.silent
  }
  protected_nonempty := hP
  observe := D.observe
  consistent := D.consistent
  target_iff := by
    intro b q
    simpa [FixedBehaviorModel.IsTarget] using D.target_iff b q
  initial_observe := D.initial_observe
  rewrite := D.rewrite
  rewrite_observe := D.rewrite_observe
  normal_iff_consistent := D.normal_iff_consistent

end CommonProfileData

structure VariantFamily (Variant : Type wGC) (O : Type uGC) (Q : Type vGC)
    [Fintype Variant] [DecidableEq Variant]
    [Fintype O] [DecidableEq O] [Fintype Q] [DecidableEq Q] where
  protectedCarrier : Finset O
  protected_nonempty : protectedCarrier.Nonempty
  data : Variant → CommonProfileData O Q protectedCarrier

namespace VariantFamily

variable {Variant : Type wGC}
  [Fintype Variant] [DecidableEq Variant]

def profile (F : VariantFamily Variant O Q) (i : Variant) : Profile O Q :=
  (F.data i).toProfile F.protected_nonempty

theorem profile_protected_rfl (F : VariantFamily Variant O Q) (i : Variant) :
    Eq (F.profile i).core.protectedSet F.protectedCarrier := rfl

def cutAt (F : VariantFamily Variant O Q) (j : Fin 4) (i : Variant) : Set O :=
  if j = 0 then
    BehaviorCut (F.profile i).supportedDynamics (F.protectedCarrier : Set O)
      (F.profile i).Admissible0 (F.profile i).traceObserve ()
  else if j = 1 then
    BehaviorCut (F.profile i).supportedDynamics (F.protectedCarrier : Set O)
      (F.profile i).Admissible1 (F.profile i).traceObserve ()
  else if j = 2 then
    BehaviorCut (F.profile i).supportedDynamics (F.protectedCarrier : Set O)
      (F.profile i).Admissible2 (F.profile i).traceObserve ()
  else
    BehaviorCut (F.profile i).supportedDynamics (F.protectedCarrier : Set O)
      (F.profile i).Admissible3 (F.profile i).traceObserve ()

theorem cutAt_zero (F : VariantFamily Variant O Q) (i : Variant) :
    Eq (F.cutAt 0 i) (F.profile i).cut0 := by
  simpa [cutAt, profile, CommonProfileData.toProfile] using
    (F.profile i).t7_behaviorCut0

theorem cutAt_one (F : VariantFamily Variant O Q) (i : Variant) :
    Eq (F.cutAt 1 i) (F.profile i).cut1 := by
  simpa [cutAt, profile, CommonProfileData.toProfile] using
    (F.profile i).t7_behaviorCut1

theorem cutAt_two (F : VariantFamily Variant O Q) (i : Variant) :
    Eq (F.cutAt 2 i) (F.profile i).cut2 := by
  simpa [cutAt, profile, CommonProfileData.toProfile] using
    (F.profile i).t7_behaviorCut2

theorem cutAt_three (F : VariantFamily Variant O Q) (i : Variant) :
    Eq (F.cutAt 3 i) (F.profile i).cut3 := by
  simpa [cutAt, profile, CommonProfileData.toProfile] using
    (F.profile i).t7_behaviorCut3

def LE (F : VariantFamily Variant O Q) (i k : Variant) : Prop :=
  ProfileLE F.cutAt i k

def LT (F : VariantFamily Variant O Q) (i k : Variant) : Prop :=
  ProfileLT F.cutAt i k

theorem le_refl (F : VariantFamily Variant O Q) (i : Variant) : F.LE i i :=
  profileLE_refl F.cutAt i

theorem le_trans (F : VariantFamily Variant O Q) {i j k : Variant} :
    F.LE i j → F.LE j k → F.LE i k :=
  profileLE_trans F.cutAt

end VariantFamily

/-! ## Frozen protected-carrier and observation package -/

variable {O₁ : Type uGC} {Q₁ : Type vGC} {O₂ : Type wGC} {Q₂ : Type xGC}
  {O₃ : Type yGC} {Q₃ : Type zGC}
  [Fintype O₁] [DecidableEq O₁] [Fintype Q₁] [DecidableEq Q₁]
  [Fintype O₂] [DecidableEq O₂] [Fintype Q₂] [DecidableEq Q₂]
  [Fintype O₃] [DecidableEq O₃] [Fintype Q₃] [DecidableEq Q₃]

structure FrozenMorphism (A : Profile O₁ Q₁) (B : Profile O₂ Q₂) where
  core : StructuralMorphism A.core B.core
  protectedEquiv : {b // b ∈ A.core.protectedSet} ≃
    {d // d ∈ B.core.protectedSet}
  protectedEquiv_coe : ∀ b, (protectedEquiv b).1 = core.behaviorMap b.1
  observation_commute : ∀ q, B.observe (core.stateMap q) =
    core.behaviorMap (A.observe q)

namespace FrozenMorphism

def id (A : Profile O₁ Q₁) : FrozenMorphism A A where
  core := CoreStructuralMorphism.id A.core
  protectedEquiv := Equiv.refl _
  protectedEquiv_coe := fun _ => rfl
  observation_commute := fun _ => rfl

def comp {A : Profile O₁ Q₁} {B : Profile O₂ Q₂}
    {C : Profile O₃ Q₃} (G : FrozenMorphism B C)
    (F : FrozenMorphism A B) : FrozenMorphism A C where
  core := CoreStructuralMorphism.comp G.core F.core
  protectedEquiv := F.protectedEquiv.trans G.protectedEquiv
  protectedEquiv_coe := by
    intro b
    change (G.protectedEquiv (F.protectedEquiv b)).1 =
      G.core.behaviorMap (F.core.behaviorMap b.1)
    rw [G.protectedEquiv_coe, F.protectedEquiv_coe]
  observation_commute := by
    intro q
    change C.observe (G.core.stateMap (F.core.stateMap q)) =
      G.core.behaviorMap (F.core.behaviorMap (A.observe q))
    rw [G.observation_commute, F.observation_commute]

theorem toObservationCompatibility {A : Profile O₁ Q₁}
    {B : Profile O₂ Q₂} (F : FrozenMorphism A B) :
    ObservationCompatibility F.core A.observe B.observe :=
  ⟨F.observation_commute⟩

theorem carrier_map_bijective {A : Profile O₁ Q₁}
    {B : Profile O₂ Q₂} (F : FrozenMorphism A B) :
    Function.Bijective F.protectedEquiv := F.protectedEquiv.bijective

end FrozenMorphism

end

/-! ## Finite full-support initial-law quantitative corollaries -/

open scoped BigOperators
open ObservableNormalForms
open ObservableNormalForms.ProtectedObstructions.NonidentityExactness

universe uIL vIL wIL xIL

noncomputable section

variable {O : Type uIL} {Q : Type vIL}
  [Fintype O] [DecidableEq O] [Fintype Q] [DecidableEq Q]

/-- A finite probability law whose support is exactly the active source set.
This is quantitative input, not a qualitative-layer premise. -/
structure InitialLaw (M : CoreModel O Q) (b : O) where
  mass : Q → ℝ
  mass_nonneg : ∀ q, 0 ≤ mass q
  mass_sum_one : ∑ q : Q, mass q = 1
  mass_pos_iff_initial : ∀ q, 0 < mass q ↔ q ∈ M.initial b

namespace InitialLaw

def pooledHit {M : CoreModel O Q} {b : O} (law : InitialLaw M b) : ℝ :=
  ∑ q : Q, law.mass q * M.hit b q

def pooledEndpoint {M : CoreModel O Q} {b : O} (law : InitialLaw M b)
    (c : Q) : ℝ :=
  ∑ q : Q, law.mass q * M.endpointMass b q c

/-- Pairwise singleton-support form for the pooled quotient endpoint law. -/
def PooledEndpointUnique {M : CoreModel O Q} {b : O}
    (law : InitialLaw M b) : Prop :=
  (∃ c, 0 < law.pooledEndpoint c) ∧
  ∀ c, 0 < law.pooledEndpoint c →
    ∀ d, 0 < law.pooledEndpoint d → M.silent.r c d

theorem pooledHit_nonneg {M : CoreModel O Q} {b : O}
    (law : InitialLaw M b) : 0 ≤ law.pooledHit := by
  exact Finset.sum_nonneg fun q _ =>
    mul_nonneg (law.mass_nonneg q) (hitProbability_nonneg M.kernel (M.IsTarget b) q)

theorem pooledEndpoint_nonneg {M : CoreModel O Q} {b : O}
    (law : InitialLaw M b) (c : Q) : 0 ≤ law.pooledEndpoint c := by
  exact Finset.sum_nonneg fun q _ =>
    mul_nonneg (law.mass_nonneg q) (M.endpointMass_nonneg b q c)

theorem pooledHit_pos_iff {M : CoreModel O Q} {b : O}
    (law : InitialLaw M b) :
    0 < law.pooledHit ↔ M.Positive b := by
  rw [pooledHit, Finset.sum_pos_iff_of_nonneg]
  · constructor
    · rintro ⟨q, hq, hprod⟩
      have hm : 0 < law.mass q :=
        pos_of_mul_pos_left hprod (hitProbability_nonneg M.kernel (M.IsTarget b) q)
      have hh : 0 < M.hit b q :=
        pos_of_mul_pos_right hprod (law.mass_nonneg q)
      exact ⟨q, (law.mass_pos_iff_initial q).mp hm, hh⟩
    · rintro ⟨q, hq, hh⟩
      exact ⟨q, Finset.mem_univ q,
        mul_pos ((law.mass_pos_iff_initial q).mpr hq) hh⟩
  · intro q hq
    exact mul_nonneg (law.mass_nonneg q)
      (hitProbability_nonneg M.kernel (M.IsTarget b) q)

theorem pooledEndpoint_pos_iff {M : CoreModel O Q} {b : O}
    (law : InitialLaw M b) (c : Q) :
    0 < law.pooledEndpoint c ↔
      ∃ q ∈ M.initial b, 0 < M.endpointMass b q c := by
  rw [pooledEndpoint, Finset.sum_pos_iff_of_nonneg]
  · constructor
    · rintro ⟨q, hq, hprod⟩
      have hm : 0 < law.mass q :=
        pos_of_mul_pos_left hprod (M.endpointMass_nonneg b q c)
      have he : 0 < M.endpointMass b q c :=
        pos_of_mul_pos_right hprod (law.mass_nonneg q)
      exact ⟨q, (law.mass_pos_iff_initial q).mp hm, he⟩
    · rintro ⟨q, hq, he⟩
      exact ⟨q, Finset.mem_univ q,
        mul_pos ((law.mass_pos_iff_initial q).mpr hq) he⟩
  · intro q hq
    exact mul_nonneg (law.mass_nonneg q) (M.endpointMass_nonneg b q c)

theorem pooledEndpointUnique_iff {M : CoreModel O Q} {b : O}
    (law : InitialLaw M b) :
    law.PooledEndpointUnique ↔ M.EndpointUnique b := by
  constructor
  · rintro ⟨⟨c, hc⟩, hunique⟩
    rcases (law.pooledEndpoint_pos_iff c).mp hc with ⟨q, hq, hqc⟩
    refine ⟨⟨q, hq, c, hqc⟩, ?_⟩
    intro x hx d hd y hy e he
    exact hunique d ((law.pooledEndpoint_pos_iff d).mpr ⟨x, hx, hd⟩)
      e ((law.pooledEndpoint_pos_iff e).mpr ⟨y, hy, he⟩)
  · rintro ⟨⟨q, hq, c, hc⟩, hunique⟩
    refine ⟨⟨c, (law.pooledEndpoint_pos_iff c).mpr ⟨q, hq, hc⟩⟩, ?_⟩
    intro d hd e he
    rcases (law.pooledEndpoint_pos_iff d).mp hd with ⟨x, hx, hxd⟩
    rcases (law.pooledEndpoint_pos_iff e).mp he with ⟨y, hy, hye⟩
    exact hunique x hx d hxd y hy e hye

private theorem one_sub_pooledHit {M : CoreModel O Q} {b : O}
    (law : InitialLaw M b) :
    1 - law.pooledHit = ∑ q : Q, law.mass q * (1 - M.hit b q) := by
  calc
    1 - law.pooledHit =
        (∑ q : Q, law.mass q) - law.pooledHit := by
      exact congrArg (fun z : ℝ => z - law.pooledHit) law.mass_sum_one.symm
    _ = (∑ q : Q, law.mass q) -
        ∑ q : Q, law.mass q * M.hit b q := by rfl
    _ =
        ∑ q : Q, (law.mass q - law.mass q * M.hit b q) := by
      rw [Finset.sum_sub_distrib]
    _ = ∑ q : Q, law.mass q * (1 - M.hit b q) := by
      apply Finset.sum_congr rfl
      intro q hq
      ring

theorem pooledHit_eq_one_iff {M : CoreModel O Q} {b : O}
    (law : InitialLaw M b) :
    law.pooledHit = 1 ↔ M.AlmostSure b := by
  constructor
  · intro hone q hq
    have hsum : ∑ r : Q, law.mass r * (1 - M.hit b r) = 0 := by
      rw [← one_sub_pooledHit law, hone, sub_self]
    have hnonneg : ∀ r ∈ (Finset.univ : Finset Q),
        0 ≤ law.mass r * (1 - M.hit b r) := by
      intro r hr
      exact mul_nonneg (law.mass_nonneg r)
        (sub_nonneg.mpr (hitProbability_le_one M.kernel (M.IsTarget b) r))
    have hz : law.mass q * (1 - M.hit b q) = 0 :=
      (Finset.sum_eq_zero_iff_of_nonneg hnonneg).mp hsum q (Finset.mem_univ q)
    have hm := (law.mass_pos_iff_initial q).mpr hq
    rcases mul_eq_zero.mp hz with hm0 | hdef
    · exact False.elim (ne_of_gt hm hm0)
    · linarith
  · intro hall
    unfold pooledHit
    calc
      (∑ q : Q, law.mass q * M.hit b q) = ∑ q : Q, law.mass q := by
        apply Finset.sum_congr rfl
        intro q hq
        by_cases hm : 0 < law.mass q
        · rw [hall q ((law.mass_pos_iff_initial q).mp hm), mul_one]
        · have hm0 : law.mass q = 0 :=
            le_antisymm (le_of_not_gt hm) (law.mass_nonneg q)
          simp [hm0]
      _ = 1 := law.mass_sum_one

end InitialLaw

/-! ## Quantitative transport, derived from the canonical endpoint theorem -/

variable {O₁ : Type uIL} {Q₁ : Type vIL} {O₂ : Type wIL} {Q₂ : Type xIL}
  [Fintype O₁] [DecidableEq O₁] [Fintype Q₁] [DecidableEq Q₁]
  [Fintype O₂] [DecidableEq O₂] [Fintype Q₂] [DecidableEq Q₂]

namespace Quantitative

def PushesInitial {A : CoreModel O₁ Q₁} {B : CoreModel O₂ Q₂}
    (F : StructuralMorphism A B) {b : O₁}
    (lawA : InitialLaw A b) (lawB : InitialLaw B (F.behaviorMap b)) : Prop :=
  ∀ y, lawB.mass y = F.fiberSum y lawA.mass

theorem expectation_pushforward
    {A : CoreModel O₁ Q₁} {B : CoreModel O₂ Q₂}
    (F : StructuralMorphism A B) {b : O₁}
    (lawA : InitialLaw A b) (lawB : InitialLaw B (F.behaviorMap b))
    (hpush : PushesInitial F lawA lawB) (phi : Q₂ → ℝ) :
    (∑ y : Q₂, lawB.mass y * phi y) =
      ∑ q : Q₁, lawA.mass q * phi (F.stateMap q) := by
  classical
  unfold PushesInitial at hpush
  simp_rw [hpush]
  unfold StructuralMorphism.fiberSum
  calc
    (∑ y : Q₂, (∑ q : Q₁, if F.stateMap q = y then lawA.mass q else 0) * phi y) =
        ∑ y : Q₂, ∑ q : Q₁,
          if F.stateMap q = y then lawA.mass q * phi y else 0 := by
      apply Finset.sum_congr rfl
      intro y hy
      rw [Finset.sum_mul]
      apply Finset.sum_congr rfl
      intro q hq
      by_cases h : F.stateMap q = y <;> simp [h]
    _ = ∑ q : Q₁, ∑ y : Q₂,
          if F.stateMap q = y then lawA.mass q * phi y else 0 := by
      rw [Finset.sum_comm]
    _ = ∑ q : Q₁, lawA.mass q * phi (F.stateMap q) := by
      apply Finset.sum_congr rfl
      intro q hq
      rw [Finset.sum_eq_single (F.stateMap q)]
      · simp
      · intro y hy hne
        simp [Ne.symm hne]
      · simp

theorem pooledHit_exact_derived
    {A : CoreModel O₁ Q₁} {B : CoreModel O₂ Q₂}
    (F : StructuralMorphism A B) {b : O₁} (hb : b ∈ A.protectedSet)
    (lawA : InitialLaw A b) (lawB : InitialLaw B (F.behaviorMap b))
    (hpush : PushesInitial F lawA lawB) :
    lawB.pooledHit = lawA.pooledHit := by
  unfold InitialLaw.pooledHit
  rw [expectation_pushforward F lawA lawB hpush
    (fun y => B.hit (F.behaviorMap b) y)]
  apply Finset.sum_congr rfl
  intro q hq
  rw [F.hitProbability_exact_derived hb]

theorem pooledEndpoint_pushforward_derived
    {A : CoreModel O₁ Q₁} {B : CoreModel O₂ Q₂}
    (F : StructuralMorphism A B) {b : O₁} (hb : b ∈ A.protectedSet)
    (lawA : InitialLaw A b) (lawB : InitialLaw B (F.behaviorMap b))
    (hpush : PushesInitial F lawA lawB) (d : Q₂) :
    lawB.pooledEndpoint d =
      F.fiberSum d (fun c => lawA.pooledEndpoint c) := by
  classical
  unfold InitialLaw.pooledEndpoint
  rw [expectation_pushforward F lawA lawB hpush
    (fun y => B.endpointMass (F.behaviorMap b) y d)]
  simp_rw [F.endpointHitProbability_pushforward_derived hb]
  unfold StructuralMorphism.fiberSum
  calc
    (∑ q : Q₁, lawA.mass q *
        ∑ c : Q₁, if F.stateMap c = d then A.endpointMass b q c else 0) =
        ∑ q : Q₁, ∑ c : Q₁,
          if F.stateMap c = d then
            lawA.mass q * A.endpointMass b q c else 0 := by
      apply Finset.sum_congr rfl
      intro q hq
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro c hc
      by_cases h : F.stateMap c = d <;> simp [h]
    _ = ∑ c : Q₁, ∑ q : Q₁,
          if F.stateMap c = d then
            lawA.mass q * A.endpointMass b q c else 0 := by
      rw [Finset.sum_comm]
    _ = ∑ c : Q₁, if F.stateMap c = d then
          ∑ q : Q₁, lawA.mass q * A.endpointMass b q c else 0 := by
      apply Finset.sum_congr rfl
      intro c hc
      by_cases h : F.stateMap c = d <;> simp [h]

end Quantitative

end

end GenericCompletion

end ObservableNormalForms.ProtectedObstructions.PublicAdapter
