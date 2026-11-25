import ExteriorSelection

/-! Neutral-current identification with a calibrated charge axis.
The group acts on the five-mode split by blocks (C,W), det C * det W = 1.
The signed strengths use a unit transition in one weak doublet.
No free-flight or source-to-exterior identification is assumed. -/

namespace OPH.NeutralObserver

open OPH.ExteriorSelection

/-- The source dictionary, in units 6Y, including both charged singlet stocks. -/
theorem exterior_slots :
    charge 9 = -3 ∧ charge 1 = 3 ∧ charge 7 = -6 ∧ charge 4 = 6 ∧
    isDoublet 9 = true ∧ isDoublet 1 = true ∧
    isDoublet 7 = false ∧ isDoublet 4 = false ∧
    conj 9 = 1 ∧ conj 7 = 4 := by decide

/-- None denotes a trivial exterior line. For a weak doublet a true second
coordinate denotes the upper slot. -/
abbrev Site := Option (Fin 10 × Bool)

def signCharge6 : Site → ℤ
  | none => 0
  | some (i, upper) => charge i + if isDoublet i then (if upper then 3 else -3) else 0

def strength (lowering : Bool) : Site → ℤ
  | none => 0
  | some (i, upper) => if isDoublet i && upper == lowering then 1 else 0

/-- The spectator has weak rank one. The last two singletons in the single-site
set have identical site records (vacuum and top). -/
def sample : Fin 7 → Site × Site := ![
  (none, some (9,true)), (none, some (1,false)), (none,none), (none,none),
  (some (7,false), some (1,true)), (some (4,false), some (9,false)),
  (some (7,false), some (4,false))]

def columns (i : Fin 7) : ℤ × ℤ × ℤ :=
  let (a,b) := sample i
  (strength true b, strength false b, signCharge6 a ^ 2)

theorem neutral_stocks (i : Fin 7) :
    signCharge6 (sample i).1 + signCharge6 (sample i).2 = 0 ∧
    strength true (sample i).1 = 0 ∧ strength false (sample i).1 = 0 := by
  fin_cases i <;> decide

/-- Requiring both signs and inner charge resolves six response classes.
Even this catalogue cannot separate the two trivial exterior lines. -/
theorem catalogue (i j : Fin 7) :
    columns i = columns j ↔
      i = j ∨ (i = 2 ∧ j = 3) ∨ (i = 3 ∧ j = 2) := by
  fin_cases i <;> fin_cases j <;> decide

/-- At every group element the negative spectator times the conjugate doublet
acts by d • W. The positive spectator times the direct doublet acts by W.
d = det C, w = det W. These are representation equalities in arbitrary frames. -/
theorem gauge_twists {d w : ℂ} (balance : d * w = 1)
    (W : Matrix (Fin 2) (Fin 2) ℂ) :
    w⁻¹ • W = d • W ∧ w • (d • W) = W := by
  have hw : w ≠ 0 := by intro hz; simp [hz] at balance
  have hd : d = w⁻¹ := by
    calc d = (d * w) * w⁻¹ := by simp [hw]
         _ = w⁻¹ := by rw [balance, one_mul]
  refine ⟨by rw [hd], ?_⟩
  rw [smul_smul, mul_comm w d, balance, one_smul]

/-- Transporting source, output charge analyzer and weak transition together
conjugates the composed effect. b is the inverse (adjoint) frame, with b*a=1.
Here k denotes J†, p denotes the final-charge projector, and j denotes J. -/
theorem transported_effect {R : Type*} [Monoid R] (a b k p j : R) (h : b * a = 1) :
    (a*k*b) * (a*p*b) * (a*j*b) = a*(k*p*j)*b := by
  simp only [mul_assoc, ← mul_assoc b a, h, one_mul]

/-- The same observation holds for every ordered word in global detectors. -/
theorem intertwined_word {R : Type*} [Monoid R] {X : Type*}
    (A B : X → R) (u : R) (h : ∀ x, A x * u = u * B x) (xs : List X) :
    (xs.map A).prod * u = u * (xs.map B).prod := by
  induction xs with
  | nil => simp
  | cons x xs ih =>
    simp only [List.map_cons, List.prod_cons]
    rw [mul_assoc, ih, ← mul_assoc, h, mul_assoc]

/-- Signed counts with independently calibrated confusion probabilities
determine both active weights exactly when the two calibration columns span
the plane (sufficiency direction). Missed launches require an independent
arrival measurement to determine the remaining dark population. -/
theorem signed_calibration {a b c d p q r s : ℝ} (det : a*d - b*c ≠ 0)
    (minus : a*p+b*q = a*r+b*s) (plus : c*p+d*q = c*r+d*s) :
    p = r ∧ q = s := by
  have hp : (a*d-b*c)*(p-r)=0 := by linear_combination d*minus - b*plus
  have hq : (a*d-b*c)*(q-s)=0 := by linear_combination a*plus - c*minus
  exact ⟨sub_eq_zero.mp ((mul_eq_zero.mp hp).resolve_left det),
    sub_eq_zero.mp ((mul_eq_zero.mp hq).resolve_left det)⟩

end OPH.NeutralObserver
