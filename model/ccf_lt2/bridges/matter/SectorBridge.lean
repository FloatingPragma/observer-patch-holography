import ObserverPatchHolography.EinsteinBranch.PerCutCollarComposition

namespace CCF2MatterBridge

def ChargeSupport {S : Type*} (H : S → S → ℂ) (q : S → ℤ) : Prop :=
  ∀ i j, H i j ≠ 0 → q i = q j

theorem support_commutator_entry {S : Type*} (H : S → S → ℂ)
    (q : S → ℤ) (h : ChargeSupport H q) (i j : S) :
    (q i : ℂ) * H i j - H i j * (q j : ℂ) = 0 := by
  by_cases hz : H i j = 0
  · simp [hz]
  · rw [h i j hz]
    ring

theorem support_preserves_sector {S : Type*} [Fintype S]
    (H : S → S → ℂ) (q : S → ℤ) (h : ChargeSupport H q)
    (v : S → ℂ) (k : ℤ) (hv : ∀ j, q j ≠ k → v j = 0)
    (i : S) (hi : q i ≠ k) : (∑ j, H i j * v j) = 0 := by
  apply Finset.sum_eq_zero
  intro j _
  by_cases hz : H i j = 0
  · simp [hz]
  · have hj : q j ≠ k := by
      intro hk
      exact hi ((h i j hz).trans hk)
    simp [hv j hj]

noncomputable def lapseJoin {A S : Type*} [DecidableEq S]
    (aux lapse : A → A → ℂ) (H : S → S → ℂ)
    (i j : A × S) : ℂ :=
  aux i.1 j.1 * (if i.2 = j.2 then 1 else 0) + lapse i.1 j.1 * H i.2 j.2

theorem support_lapseJoin {A S : Type*} [DecidableEq S]
    (aux lapse : A → A → ℂ) (H : S → S → ℂ)
    (q : S → ℤ) (h : ChargeSupport H q) :
    ChargeSupport (lapseJoin aux lapse H) (fun i => q i.2) := by
  intro i j hn
  by_contra hq
  have hij : i.2 ≠ j.2 := by
    intro he
    exact hq (congrArg q he)
  have hz : H i.2 j.2 = 0 := by
    by_contra hh
    exact hq (h i.2 j.2 hh)
  exact hn (by simp [lapseJoin, hij, hz])

/-- The independently joined local capacity terms preserve the same sector
when each term does. Their supports may overlap and need not commute. -/
theorem support_sum {V S : Type*} (vertices : Finset V)
    (terms : V → S → S → ℂ) (q : S → ℤ)
    (h : ∀ v, ChargeSupport (terms v) q) :
    ChargeSupport (fun i j => ∑ v ∈ vertices, terms v i j) q := by
  intro i j hn
  by_contra hq
  apply hn
  apply Finset.sum_eq_zero
  intro v _
  by_contra hz
  exact hq (h v i j hz)

theorem gibbs_contrast_offset_cancels {V : Type*} [AddCommGroup V]
    [Module ℝ V] (contrast : V →ₗ[ℝ] ℝ) (H identity K : V)
    (beta logZ : ℝ) (hK : K = beta • H + logZ • identity)
    (hi : contrast identity = 0) : contrast K = beta * contrast H := by
  simp [hK, hi]

#print axioms support_commutator_entry
#print axioms support_preserves_sector
#print axioms support_lapseJoin
#print axioms support_sum
#print axioms gibbs_contrast_offset_cancels
#print axioms OPH.EinsteinBranch.geometry_responds_to_all_charge
#print axioms OPH.EinsteinBranch.exists_dark_source_of_geometric_excess
#print axioms OPH.EinsteinBranch.anomalousStress_tendsto_zero
#print axioms OPH.EinsteinBranch.deepProfile_characterization
#print axioms OPH.EinsteinBranch.perCut_a0_dictionary

end CCF2MatterBridge
