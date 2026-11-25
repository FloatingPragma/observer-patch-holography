import Mathlib

namespace OPH.ProductShell

/-- Record of one outward nearest-neighbour crossing. The first two
registers name its unique changing axis and the face on that axis.
The other two coordinates label its tangential positions. -/
abbrev FaceRead (n : ℕ) := Fin 3 × Bool × Fin (2*n) × Fin (2*n)

theorem faceRead_card (n : ℕ) :
    Fintype.card (FaceRead n) = 24*n^2 := by
  simp [FaceRead]
  ring

end OPH.ProductShell
