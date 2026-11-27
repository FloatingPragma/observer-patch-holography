import Mathlib.Tactic

/-! Channel thresholds and a polarization sign after dressing. Energies use
the native exchange unit in the first theorem. The interval hypotheses
are supplied by the separate integer word census. -/
namespace IntervalEcho.SignedPair

/-- The same inventory may bind relative to one partition and allow a
different partition below every joined energy. -/
theorem channel_order {positive neutral joined allExcess allVacancy : ℝ}
    (p : 112 ≤ positive) (n : 136 ≤ neutral) (jlo : 200 ≤ joined)
    (jhi : joined ≤ 209) (e : allExcess ≤ 133) (v : allVacancy ≤ 52.5) :
    joined < positive + neutral ∧ allExcess + allVacancy < joined := by
  constructor <;> linarith

/- The three physical rotations of the two sharp competing fragments. A
mask labels the twelve physical ports. A native swap can connect distinct
masks only if their symmetric difference consists of its two endpoints. -/
private def seams : List (Nat × Nat) :=
  [(0,1),(0,5),(0,7),(0,10),(0,11),(1,5),(1,7),(1,8),(1,9),
   (2,3),(2,4),(2,6),(2,10),(2,11),(3,4),(3,6),(3,8),(3,9),
   (4,5),(4,9),(4,11),(5,9),(5,11),(6,7),(6,8),(6,10),(7,8),
   (7,10),(8,9),(10,11)]

private def plusOrbit := [1155,778,1100]
private def minusOrbit := [2081,560,2068]
private def turn : Array Nat := #[9,3,0,10,11,4,7,8,6,2,1,5]
private def rotate (m : Nat) :=
  (List.range 12).foldl (fun b i => if m.testBit i then b + 2^turn[i]! else b) 0

private def diagonalTwice (m : Nat) (positive : Bool) : Int :=
  let cut : Int := ((seams.filter fun e => m.testBit e.1 != m.testBit e.2).length : Nat)
  let referenced : Int := (([0,2,4,5,6,7,8,9,11].filter m.testBit).length : Nat)
  16 * (cut + (if positive then 3 else -1) * referenced) + cut

private def isolatedOrbit (orbit : List Nat) : Bool := orbit.Nodup &&
  orbit.all (fun x => orbit.all (fun y => x == y ||
    seams.all (fun e => (x ^^^ y) != 2^e.1 + 2^e.2)))

theorem competing_product_orbit :
    plusOrbit.map rotate = [778,1100,1155] ∧
    minusOrbit.map rotate = [560,2068,2081] ∧
    plusOrbit.map (fun x => diagonalTwice x true) = [266,266,266] ∧
    minusOrbit.map (fun x => diagonalTwice x false) = [105,105,105] ∧
    isolatedOrbit plusOrbit = true ∧ isolatedOrbit minusOrbit = true := by decide +kernel

/-- The simultaneous orbit has twice the energy compression 371 I. All
three normalized Fourier columns have the same trial energy. -/
theorem projected_trial {a b c : ℝ} (normalized : a + b + c = 1) :
    185.5 * a + 185.5 * b + 185.5 * c = 185.5 := by linarith

/-- Charge and centroid use the same nonnegative material weights on the old
six-load caps. The complementary words and the cap words have disjoint
support, so a diagonal quadrupole is a convex sum on these two classes. -/
theorem circular_band_positive {p cap rest : ℝ}
    (hp : 0 ≤ p) (leak : p ≤ 10 / (4 * 11 ^ 2)) (ideal : 2/3 ≤ cap)
    (floor : -6 ≤ rest) : 0 < (1-p)*cap + p*rest := by
  have hw : 0 ≤ 1-p := by linarith
  nlinarith [mul_nonneg hw (sub_nonneg.mpr ideal),
             mul_nonneg hp (sub_nonneg.mpr floor)]

theorem axial_band_negative {p cap rest : ℝ}
    (hp : 0 ≤ p) (leak : p ≤ 10 / (4 * 11 ^ 2)) (ideal : cap ≤ -4/3)
    (ceil : rest ≤ 12) : (1-p)*cap + p*rest < 0 := by
  have hw : 0 ≤ 1-p := by linarith
  nlinarith [mul_nonneg hw (sub_nonneg.mpr ideal),
             mul_nonneg hp (sub_nonneg.mpr ceil)]

end IntervalEcho.SignedPair
