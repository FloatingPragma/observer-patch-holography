import ObserverPatchHolography.CapacityStress

/-! Algebraic compatibility conditions for declared source-unit and capacity
maps. These theorems supply no physical choice of stiffness or horizon. -/
namespace OPH.SourceScaleCompatibility
noncomputable section

def lightRatio (u length hc : ℝ) := u * length / hc
def newtonValue (g u length hc : ℝ) := g * length / (u * hc)

theorem gravity_light_invariant (g u length hc : ℝ)
    (hu : u ≠ 0) (hl : length ≠ 0) (hh : hc ≠ 0) :
    newtonValue g u length hc / lightRatio u length hc = g / u^2 := by
  unfold newtonValue lightRatio
  field_simp

theorem simultaneous_match_requires (g u length hc target : ℝ)
    (hu : u ≠ 0) (hl : length ≠ 0) (hh : hc ≠ 0)
    (hv : lightRatio u length hc = 1)
    (hg : newtonValue g u length hc = target) : g = target * u^2 := by
  have hi := gravity_light_invariant g u length hc hu hl hh
  rw [hg, hv, div_one] at hi
  exact (eq_div_iff (pow_ne_zero 2 hu)).mp hi |>.symm

theorem speed_matched_length (u hc : ℝ) (hu : u ≠ 0) (hh : hc ≠ 0) :
    lightRatio u (hc/u) hc = 1 := by
  unfold lightRatio
  field_simp

theorem hubble_area_dust_branch (H2 H02 om scale : ℝ) (hom : om ≠ 0)
    (friedmann : H2 = om * H02 * scale + (1-om) * H2) :
    H2 = H02 * scale := by
  have h : om * (H2 - H02*scale) = 0 := by linarith
  exact sub_eq_zero.mp ((mul_eq_zero.mp h).resolve_left hom)

theorem comoving_area_equation_of_state (N0 : ℝ) (hN0 : 0 < N0)
    {a : ℝ} (ha : a ∈ Set.Ioi (0 : ℝ)) :
    (OPH.EinsteinBranch.powerCapacity N0 hN0 2).w a = -1/3 := by
  rw [OPH.EinsteinBranch.powerCapacity_w N0 hN0 2 ha]
  norm_num

theorem inverse_spacing_energy (energy a : ℝ) (ha : a ≠ 0) :
    (energy / a) * a = energy := by
  field_simp

#print axioms gravity_light_invariant
#print axioms simultaneous_match_requires
#print axioms speed_matched_length
#print axioms hubble_area_dust_branch
#print axioms comoving_area_equation_of_state
#print axioms inverse_spacing_energy
end
end OPH.SourceScaleCompatibility
