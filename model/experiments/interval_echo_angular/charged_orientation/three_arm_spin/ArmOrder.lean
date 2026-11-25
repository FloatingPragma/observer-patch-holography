import Mathlib.Data.Complex.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-! # Arm-order interference on three one-hot head amplitudes

The two operations are literal swaps of physical arm labels. The sum-zero
subspace is the standard two-dimensional arm representation. The balanced
interferometer reads the squared norm of the half-difference of its two
branches. The result concerns this permutation action on amplitudes; a
physical circuit that accumulates additional phases needs its own evidence.
-/

namespace OPH.ThreeArmSpin.ArmOrder

structure Amplitudes where
  a : ℂ
  b : ℂ
  c : ℂ

def swapAB (v : Amplitudes) : Amplitudes := ⟨v.b, v.a, v.c⟩
def swapBC (v : Amplitudes) : Amplitudes := ⟨v.a, v.c, v.b⟩

def firstOrder (v : Amplitudes) : Amplitudes := swapAB (swapBC v)
def otherOrder (v : Amplitudes) : Amplitudes := swapBC (swapAB v)
def opposite (v : Amplitudes) : Amplitudes := ⟨-v.a, -v.b, -v.c⟩

def power (z : ℂ) : ℝ := z.re ^ 2 + z.im ^ 2

def weight (v : Amplitudes) : ℝ := power v.a + power v.b + power v.c

/-- Real part of the usual complex inner product, in coordinates. -/
def realOverlap (v w : Amplitudes) : ℝ :=
  v.a.re * w.a.re + v.a.im * w.a.im +
  v.b.re * w.b.re + v.b.im * w.b.im +
  v.c.re * w.c.re + v.c.im * w.c.im

noncomputable def dark (v w : Amplitudes) : ℝ :=
  (power (v.a - w.a) + power (v.b - w.b) + power (v.c - w.c)) / 4

theorem firstOrder_cube (v : Amplitudes) :
    firstOrder (firstOrder (firstOrder v)) = v := by cases v; rfl

theorem otherOrder_cube (v : Amplitudes) :
    otherOrder (otherOrder (otherOrder v)) = v := by cases v; rfl

/-- This identity also shows exactly what the symmetric arm adds to the fringe. -/
theorem order_overlap_identity (v : Amplitudes) :
    2 * realOverlap (firstOrder v) (otherOrder v) + weight v =
    power (v.a + v.b + v.c) := by
  simp [firstOrder, otherOrder, swapAB, swapBC, realOverlap, weight, power]
  ring

theorem order_overlap (v : Amplitudes) (h : v.a + v.b + v.c = 0) :
    realOverlap (firstOrder v) (otherOrder v) = -(weight v) / 2 := by
  have hi := order_overlap_identity v
  rw [h] at hi
  simp [power] at hi
  linarith

theorem balanced_dark_identity (v w : Amplitudes) :
    4 * dark v w = weight v + weight w - 2 * realOverlap v w := by
  simp [dark, weight, power, realOverlap]
  ring

theorem literal_order_dark (v : Amplitudes) (h : v.a + v.b + v.c = 0) :
    dark (firstOrder v) (otherOrder v) = 3 * weight v / 4 := by
  have h₁ : weight (firstOrder v) = weight v := by
    dsimp [firstOrder, swapAB, swapBC, weight]
    ring
  have h₂ : weight (otherOrder v) = weight v := by
    dsimp [otherOrder, swapAB, swapBC, weight]
    ring
  have hd := balanced_dark_identity (firstOrder v) (otherOrder v)
  rw [h₁, h₂, order_overlap v h] at hd
  linarith

theorem full_minus_dark (v : Amplitudes) : dark v (opposite v) = weight v := by
  simp [dark, opposite, weight, power]
  ring

theorem normalized_signatures (v : Amplitudes)
    (h : v.a + v.b + v.c = 0) (hn : weight v = 1) :
    realOverlap (firstOrder v) (otherOrder v) = -(1 : ℝ) / 2 ∧
    dark (firstOrder v) (otherOrder v) = (3 : ℝ) / 4 ∧
    dark v (opposite v) = 1 := by
  rw [order_overlap v h, literal_order_dark v h, full_minus_dark, hn]
  norm_num

end OPH.ThreeArmSpin.ArmOrder
