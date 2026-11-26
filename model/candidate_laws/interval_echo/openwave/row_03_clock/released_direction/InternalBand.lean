import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith

namespace OPH.ReleasedDirection

/-- On a sphere with six load contacts including private references, the
one-writable-core energy enclosures separate every internal Bohr line from
every formation energy as soon as 3 κ > 7. -/
theorem formation_exceeds_every_gap (k a b formation : ℝ)
    (hk : 7 < 3*k) (hal : 5*k ≤ a) (hau : a ≤ 7*k+7)
    (hbl : 5*k ≤ b) (hbu : b ≤ 7*k+7)
    (hformation : 5*k ≤ formation) :
    |a-b| ≤ 2*k+7 ∧ |a-b| < formation := by
  have band : |a-b| ≤ 2*k+7 := abs_le.mpr ⟨by linarith, by linarith⟩
  exact ⟨band, by linarith⟩

#print axioms formation_exceeds_every_gap

end OPH.ReleasedDirection
