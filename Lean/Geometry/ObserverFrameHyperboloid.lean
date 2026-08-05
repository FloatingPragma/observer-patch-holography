import Geometry.CanonicalLorentzModule

/-!
# C1: the future unit-timelike frame hyperboloid

This is the algebraic future-unit-timelike hyperboloid model underlying `H³`
in the declared Hermitian Lorentz module.  No hyperbolic metric or topology is
bundled here, and `FrameHyperboloid` is not attached to a physical observer,
clock, worldline, or continuum spacetime.
-/

namespace OPH.C1Lorentz

noncomputable section

/-- The internal algebraic future unit-timelike hyperboloid model. -/
def IsFutureUnitTimelike (u : Herm2) : Prop := lorentzQ u = 1 ∧ 0 < u.1

abbrev FrameHyperboloid := {u : Herm2 // IsFutureUnitTimelike u}

/-- The canonical origin of the frame hyperboloid. -/
def standardFrame : FrameHyperboloid :=
  ⟨(1, 0), by simp [IsFutureUnitTimelike, lorentzQ, spatialNormSq]⟩

theorem frame_time_sq_eq_one_add_spatial (u : FrameHyperboloid) :
    u.1.1 ^ 2 = 1 + spatialNormSq u.1.2 := by
  have h := u.2.1
  simp only [lorentzQ] at h
  linarith

theorem frame_time_ge_one (u : FrameHyperboloid) : 1 ≤ u.1.1 := by
  have hs : 0 ≤ spatialNormSq u.1.2 :=
    Finset.sum_nonneg fun i _ ↦ sq_nonneg (u.1.2 i)
  have ht := u.2.2
  have hsq := frame_time_sq_eq_one_add_spatial u
  nlinarith

#print axioms standardFrame
#print axioms frame_time_ge_one

end

end OPH.C1Lorentz
