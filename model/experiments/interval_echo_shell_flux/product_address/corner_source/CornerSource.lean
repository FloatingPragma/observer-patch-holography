import Mathlib

namespace OPH.ProductShell.CornerSource

/-- Weight of one sign triple, including the independent first-guide skew. -/
def weight (s b x y z : ℚ) := (1 + s*x*y*z + b*x) / 8

theorem twoGuideXY (s b x y : ℚ) :
    weight s b x y 1 + weight s b x y (-1) = (1+b*x)/4 := by
  unfold weight
  ring

theorem twoGuideXZ (s b x z : ℚ) :
    weight s b x 1 z + weight s b x (-1) z = (1+b*x)/4 := by
  unfold weight
  ring

theorem twoGuideYZ (s b y z : ℚ) :
    weight s b 1 y z + weight s b (-1) y z = 1/4 := by
  unfold weight
  ring

/-- One clock tick consists of the six single-coordinate changes and two idle tags.
The coordinates are translated and doubled without changing the native transitions. -/
theorem harmonicScore (x y z : ℤ) :
    (x+2)*y*z + (x-2)*y*z + x*(y+2)*z + x*(y-2)*z +
    x*y*(z+2) + x*y*(z-2) + 2*(x*y*z) = 8*(x*y*z) := by ring

/-- Eliminate the two sharp energy branches of the finite census. -/
theorem meanEnergy (b : ℚ) :
    (65/2:ℚ)*(22*((1-b)/2)+21*((1+b)/2)) = 2795/4-65*b/4 := by ring

/-- A calibrated four-point Euclidean metric must obey Ptolemy. Min/max length
ratios between the four native ranges and such a metric obey this finite bound. -/
theorem clockDistortion (lo hi : ℝ) (ptolemy : (2*lo)^2 ≤ 2*hi^2) :
    2*lo^2 ≤ hi^2 := by linarith [show (2*lo)^2 = 4*lo^2 by ring]

/-- On cube diagonal directions, a quadratic in a null vector has these seven
sign coefficients. Squared signs contribute to the constant coefficient. -/
def quadraticRay (a ux uy uz uxy uxz uyz x y z : ℚ) :=
  a + ux*x + uy*y + uz*z + uxy*x*y + uxz*x*z + uyz*y*z

def tripleContrast (f : ℚ → ℚ → ℚ → ℚ) :=
  (f 1 1 1 - f (-1) 1 1 - f 1 (-1) 1 - f 1 1 (-1) +
   f (-1) (-1) 1 + f (-1) 1 (-1) + f 1 (-1) (-1) - f (-1) (-1) (-1)) / 8

theorem quadraticProjectionZero (a ux uy uz uxy uxz uyz : ℚ) :
    tripleContrast (quadraticRay a ux uy uz uxy uxz uyz) = 0 := by
  unfold tripleContrast quadraticRay
  ring

theorem clockSourceProjection :
    tripleContrast (fun x y z => -(x*y*z)/8) = -1/8 := by
  norm_num [tripleContrast]

#print axioms twoGuideXY
#print axioms twoGuideXZ
#print axioms twoGuideYZ
#print axioms harmonicScore
#print axioms meanEnergy
#print axioms clockDistortion
#print axioms quadraticProjectionZero
#print axioms clockSourceProjection

end OPH.ProductShell.CornerSource
