import Geometry.CelestialNullCone
import ObserverPatchHolography.EinsteinBranch.Tensor

/-!
# C1: exact bridge to the Einstein-branch coordinate convention

The intrinsic Hermitian module uses signature `(+---)`, while the existing
Einstein tensor library uses `(-+++)`.  This file gives the explicit Pauli
coordinate linear equivalence and proves that the quadratic forms differ by
exactly one minus sign.  Consequently their null cones coincide on the nose.

This is a coordinate-compatibility theorem.  It does not identify either
module with physical spacetime or construct an event-frame soldering map.
-/

namespace OPH.C1Lorentz

noncomputable section

/-- Pauli coordinates in the exact `Fin 4 → ℝ` convention used by the
Einstein-branch tensor stack. -/
def pauliEinsteinChart : Herm2 ≃ₗ[ℝ] OPH.EinsteinBranch.V 3 where
  toFun v := ![v.1, v.2 0, v.2 1, v.2 2]
  invFun w := (w 0, ![w 1, w 2, w 3])
  left_inv v := by
    apply Prod.ext
    · rfl
    · funext i
      fin_cases i <;> rfl
  right_inv w := by
    funext i
    fin_cases i <;> rfl
  map_add' u v := by
    funext i
    fin_cases i <;> rfl
  map_smul' a v := by
    funext i
    fin_cases i <;> rfl

@[simp]
theorem pauliEinsteinChart_time (v : Herm2) : pauliEinsteinChart v 0 = v.1 :=
  rfl

/-- Exact sign conversion between determinant signature `(+---)` and the
Einstein branch's `(-+++)` convention. -/
theorem lorentzQ_eq_neg_einsteinQuad (v : Herm2) :
    lorentzQ v =
      -OPH.EinsteinBranch.quadOf (OPH.EinsteinBranch.eta 3)
        (pauliEinsteinChart v) := by
  simp [lorentzQ, spatialNormSq, pauliEinsteinChart,
    OPH.EinsteinBranch.quadOf, OPH.EinsteinBranch.bilinOf,
    OPH.EinsteinBranch.eta, Fin.sum_univ_succ]
  ring

/-- The two coordinate conventions have exactly the same null cone. -/
theorem lorentzQ_eq_zero_iff_einsteinQuad_eq_zero (v : Herm2) :
    lorentzQ v = 0 ↔
      OPH.EinsteinBranch.quadOf (OPH.EinsteinBranch.eta 3)
        (pauliEinsteinChart v) = 0 := by
  rw [lorentzQ_eq_neg_einsteinQuad]
  exact neg_eq_zero

/-- The future-null predicate transports exactly, including its time
orientation, into the existing Einstein coordinate type. -/
theorem isFutureNull_iff_einstein (v : Herm2) :
    IsFutureNull v ↔
      OPH.EinsteinBranch.quadOf (OPH.EinsteinBranch.eta 3)
          (pauliEinsteinChart v) = 0 ∧
        0 < pauliEinsteinChart v 0 := by
  simp only [IsFutureNull, pauliEinsteinChart_time]
  rw [lorentzQ_eq_zero_iff_einsteinQuad_eq_zero]

#print axioms pauliEinsteinChart
#print axioms lorentzQ_eq_neg_einsteinQuad
#print axioms lorentzQ_eq_zero_iff_einsteinQuad_eq_zero
#print axioms isFutureNull_iff_einstein

end

end OPH.C1Lorentz
