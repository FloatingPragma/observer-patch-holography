import Geometry.LorentzOverlapCocycle

/-!
# C1: the kinematic mass shell on the declared Lorentz module

For one declared positive mass parameter along a frame direction of the
future unit-timelike hyperboloid, the scaled vector `mass • frame` satisfies
the exact kinematic mass-shell identities: its Lorentz square is `mass ^ 2`,
its scalar coordinate squared equals `mass ^ 2` plus the spatial squared
norm, the scalar coordinate dominates the mass parameter with equality
exactly at vanishing spatial coordinates, and every oriented Lorentz map
transports the construction with the same Lorentz square.  At the standard
frame the scalar coordinate equals the mass parameter exactly, and a future
null vector satisfies the zero-mass shell identity.

The mass parameter is a declared hypothesis of each statement, not a derived
quantity.  Nothing here derives a mass from observer repair data, identifies
the scalar coordinate with a physical energy readout, attaches a clock,
calibration anchor, or SI unit, or supplies the energy-clock receipt.  The
module works in the dimensionless convention of the declared Hermitian
Lorentz module, where the invariant-speed conversion is one; no decimal
magnitude and no physical rest-energy claim is made.
-/

namespace OPH.C1Lorentz

noncomputable section

open OPH.C2Soldering

/-- The Lorentz square of a scalar multiple. -/
theorem lorentzQ_smul (scale : ℝ) (vector : Herm2) :
    lorentzQ (scale • vector) = scale ^ 2 * lorentzQ vector := by
  calc
    lorentzQ (scale • vector) =
        lorentzB (scale • vector) (scale • vector) :=
      (lorentzB_self _).symm
    _ = scale * lorentzB vector (scale • vector) :=
      lorentzB_smul_left scale vector (scale • vector)
    _ = scale * (scale * lorentzB vector vector) := by
      rw [lorentzB_smul_right]
    _ = scale ^ 2 * lorentzQ vector := by
      rw [lorentzB_self]
      ring

/-- Every oriented Lorentz map preserves the Lorentz quadratic form. -/
theorem lorentzQ_oriented (L : OrientedLorentzEquiv) (v : Herm2) :
    lorentzQ (L v) = lorentzQ v := by
  rw [← lorentzB_self, L.map_lorentzB, lorentzB_self]

/-- The declared-mass kinematic momentum along a frame direction.  The mass
parameter is a hypothesis, not a derived quantity. -/
def fourMomentum (mass : ℝ) (frame : FrameHyperboloid) : Herm2 :=
  mass • (frame : Herm2)

@[simp] theorem fourMomentum_time (mass : ℝ) (frame : FrameHyperboloid) :
    (fourMomentum mass frame).1 = mass * frame.1.1 := rfl

@[simp] theorem fourMomentum_spatial (mass : ℝ) (frame : FrameHyperboloid) :
    (fourMomentum mass frame).2 = mass • frame.1.2 := rfl

/-- The Lorentz square of the declared-mass momentum is exactly the squared
mass parameter. -/
theorem lorentzQ_fourMomentum (mass : ℝ) (frame : FrameHyperboloid) :
    lorentzQ (fourMomentum mass frame) = mass ^ 2 := by
  rw [fourMomentum, lorentzQ_smul, frame.2.1, mul_one]

/-- The kinematic mass-shell identity: the squared scalar coordinate equals
the squared mass parameter plus the spatial squared norm. -/
theorem fourMomentum_time_sq_eq_mass_sq_add_spatial (mass : ℝ)
    (frame : FrameHyperboloid) :
    (fourMomentum mass frame).1 ^ 2 =
      mass ^ 2 + spatialNormSq (fourMomentum mass frame).2 := by
  have hframe := frame_time_sq_eq_one_add_spatial frame
  simp only [fourMomentum_time, fourMomentum_spatial, spatialNormSq_smul]
  calc
    (mass * frame.1.1) ^ 2 = mass ^ 2 * frame.1.1 ^ 2 := by ring
    _ = mass ^ 2 * (1 + spatialNormSq frame.1.2) := by rw [hframe]
    _ = mass ^ 2 + mass ^ 2 * spatialNormSq frame.1.2 := by ring

/-- The scalar coordinate dominates a nonnegative mass parameter. -/
theorem mass_le_fourMomentum_time {mass : ℝ} (hm : 0 ≤ mass)
    (frame : FrameHyperboloid) :
    mass ≤ (fourMomentum mass frame).1 := by
  have h1 := frame_time_ge_one frame
  have h2 : mass * 1 ≤ mass * frame.1.1 :=
    mul_le_mul_of_nonneg_left h1 hm
  simpa using h2

/-- The scalar coordinate of a positive-mass momentum is positive. -/
theorem fourMomentum_time_pos {mass : ℝ} (hm : 0 < mass)
    (frame : FrameHyperboloid) :
    0 < (fourMomentum mass frame).1 := by
  simp only [fourMomentum_time]
  exact mul_pos hm frame.2.2

/-- Rest characterization: the scalar coordinate equals the mass parameter
exactly when the spatial coordinates vanish. -/
theorem fourMomentum_time_eq_mass_iff_spatial_zero {mass : ℝ} (hm : 0 < mass)
    (frame : FrameHyperboloid) :
    (fourMomentum mass frame).1 = mass ↔ (fourMomentum mass frame).2 = 0 := by
  constructor
  · intro hE
    simp only [fourMomentum_time] at hE
    have hu : frame.1.1 = 1 := by
      apply mul_left_cancel₀ (ne_of_gt hm)
      rw [hE, mul_one]
    have hframe := frame_time_sq_eq_one_add_spatial frame
    rw [hu] at hframe
    have hs : spatialNormSq frame.1.2 = 0 := by
      have h1 : (1 : ℝ) ^ 2 = 1 := one_pow 2
      linarith
    have hx : frame.1.2 = 0 := by
      by_contra hne
      exact (spatialNormSq_pos hne).ne' hs
    simp [hx]
  · intro hp
    simp only [fourMomentum_spatial] at hp
    have hx : frame.1.2 = 0 := by
      rcases smul_eq_zero.mp hp with hm0 | hx
      · exact absurd hm0 (ne_of_gt hm)
      · exact hx
    have hs : spatialNormSq frame.1.2 = 0 := by
      simp [hx, spatialNormSq]
    have hframe := frame_time_sq_eq_one_add_spatial frame
    rw [hs] at hframe
    have hu : frame.1.1 = 1 := by
      have hfac : (frame.1.1 - 1) * (frame.1.1 + 1) =
          frame.1.1 ^ 2 - 1 := by ring
      have hzero : (frame.1.1 - 1) * (frame.1.1 + 1) = 0 := by
        rw [hfac, hframe]
        ring
      rcases mul_eq_zero.mp hzero with h1 | h2
      · linarith
      · exfalso
        linarith [frame.2.2]
    simp [hu]

/-- The declared-mass momentum at the standard frame. -/
theorem fourMomentum_standardFrame (mass : ℝ) :
    fourMomentum mass standardFrame = (mass, 0) := by
  apply Prod.ext
  · simp [fourMomentum, standardFrame]
  · funext i
    simp [fourMomentum, standardFrame]

/-- At the standard frame the scalar coordinate is the mass parameter: the
rest form of the kinematic mass-shell identity in the module's dimensionless
convention. -/
theorem standardFrame_time_eq_mass (mass : ℝ) :
    (fourMomentum mass standardFrame).1 = mass := by
  rw [fourMomentum_standardFrame]

/-- Oriented Lorentz maps transport the declared-mass momentum to the
declared-mass momentum of the transported frame. -/
theorem oriented_map_fourMomentum (L : OrientedLorentzEquiv) (mass : ℝ)
    (frame : FrameHyperboloid) :
    L (fourMomentum mass frame) = fourMomentum mass (L.mapFrame frame) := by
  simp [fourMomentum]

/-- The squared mass parameter is invariant under every oriented Lorentz
map. -/
theorem lorentzQ_oriented_fourMomentum (L : OrientedLorentzEquiv) (mass : ℝ)
    (frame : FrameHyperboloid) :
    lorentzQ (L (fourMomentum mass frame)) = mass ^ 2 := by
  rw [lorentzQ_oriented, lorentzQ_fourMomentum]

/-- The kinematic mass-shell identity holds in every oriented Lorentz
chart. -/
theorem fourMomentum_time_sq_oriented (L : OrientedLorentzEquiv) (mass : ℝ)
    (frame : FrameHyperboloid) :
    (L (fourMomentum mass frame)).1 ^ 2 =
      mass ^ 2 + spatialNormSq (L (fourMomentum mass frame)).2 := by
  rw [oriented_map_fourMomentum]
  exact fourMomentum_time_sq_eq_mass_sq_add_spatial mass (L.mapFrame frame)

/-- The massless branch: a future null vector satisfies the zero-mass shell
identity exactly. -/
theorem futureNull_time_sq_eq_spatial (v : FutureNullVector) :
    v.1.1 ^ 2 = spatialNormSq v.1.2 := by
  have h := v.2.1
  simp only [lorentzQ] at h
  linarith

#print axioms lorentzQ_fourMomentum
#print axioms fourMomentum_time_sq_eq_mass_sq_add_spatial
#print axioms fourMomentum_time_eq_mass_iff_spatial_zero
#print axioms standardFrame_time_eq_mass
#print axioms lorentzQ_oriented_fourMomentum
#print axioms fourMomentum_time_sq_oriented
#print axioms futureNull_time_sq_eq_spatial

end

end OPH.C1Lorentz
