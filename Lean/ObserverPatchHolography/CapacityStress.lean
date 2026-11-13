import ObserverPatchHolography.EinsteinBranch.FixedCapacityWLaw
import ObserverPatchHolography.CapacityNonidentifiability

/-! Exact work of a declared capacity-volume coupling. These identities do not
select a physical capacity, a cosmological amplitude, or a metric dictionary. -/
namespace OPH.CapacityStress

noncomputable def energy (epsilon volume capacity : ℝ) : ℝ :=
  epsilon * volume / capacity

noncomputable def pressure (epsilon capacity : ℝ) : ℝ :=
  -epsilon / capacity

theorem volume_work (epsilon v₀ v₁ capacity : ℝ) :
    energy epsilon v₁ capacity - energy epsilon v₀ capacity =
      -pressure epsilon capacity * (v₁ - v₀) := by
  unfold energy pressure
  ring

theorem density (epsilon volume capacity : ℝ) (hv : volume ≠ 0) :
    energy epsilon volume capacity / volume = epsilon / capacity := by
  unfold energy
  field_simp

theorem capacity_work (epsilon volume k₀ k₁ : ℝ) :
    energy epsilon volume k₁ - energy epsilon volume k₀ =
      epsilon * volume * (1 / k₁ - 1 / k₀) := by
  unfold energy
  ring

theorem complete_cycle (epsilon v₀ v₁ k₀ k₁ : ℝ) :
    (energy epsilon v₁ k₀ - energy epsilon v₀ k₀) +
    (energy epsilon v₁ k₁ - energy epsilon v₁ k₀) +
    (energy epsilon v₀ k₁ - energy epsilon v₁ k₁) +
    (energy epsilon v₀ k₀ - energy epsilon v₀ k₁) = 0 := by
  ring

theorem pressure_is_negative_density (epsilon volume capacity : ℝ)
    (hv : volume ≠ 0) :
    pressure epsilon capacity = -energy epsilon volume capacity / volume := by
  rw [neg_div, density epsilon volume capacity hv]
  unfold pressure
  ring

#print axioms volume_work
#print axioms density
#print axioms capacity_work
#print axioms complete_cycle
#print axioms pressure_is_negative_density
end OPH.CapacityStress
