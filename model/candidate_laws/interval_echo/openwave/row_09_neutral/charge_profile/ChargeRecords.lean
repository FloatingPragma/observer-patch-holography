import Mathlib

namespace Screen.NeutralChargeRecords

/-- A net charge and a signed difference determine both component records. -/
theorem reconstruct {h t h' t' : ℤ} (net : h+t=h'+t') (gradient : t-h=t'-h') :
    h=h' ∧ t=t' := by omega

/-- A change of port order preserves both measured powers. -/
theorem swap_powers (h t : ℝ) : ((h+t)^2,(t-h)^2)=((t+h)^2,(h-t)^2) := by
  apply Prod.ext <;> dsimp <;> ring

/-- For two ternary ports a nontrivial neutral word has opposite unit charges. -/
theorem resolved_neutral {h t : ℤ} (hmin : -1 ≤ h) (hmax : h ≤ 1)
    (net : h+t=0) (resolved : t ≠ h) :
    ((h = -1 ∧ t = 1) ∨ (h = 1 ∧ t = -1)) := by omega

/-- Opposite dipole orientations can have zero mean signed contrast but retain full power. -/
theorem balanced_dipoles : (((2:ℚ)+(-2))/2 = 0) ∧
    (((2:ℚ)^2+(-2)^2)/2 = 4) := by norm_num

end Screen.NeutralChargeRecords
