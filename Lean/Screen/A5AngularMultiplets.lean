import Mathlib

namespace OPH.A5AngularMultiplets

/-! # The A5 angular multiplet branching table, exactly

The frozen angular multiplet signature (FZ-02) declares the branching of
the rotation multiplets of dimension `2l + 1` under the icosahedral group
`A5` for `l = 0` through `6`, the first nonconstant invariant at `l = 6`,
and the face-phase multiplicity vector `(0, 1, 1, 1, 2)`.  This file checks
that table in exact integer arithmetic over `Z[sqrt 5]`, with quadrupled
cosines and doubled irreducible characters so that every value is an
integer pair `(a, b)` meaning `a + b * sqrt 5`.

INPUT DATA AND ITS PINNING.

* The five class angles of `A5` acting by rotations are `0`, `pi`,
  `2 pi / 3`, `2 pi / 5`, `4 pi / 5`, with class sizes `1, 15, 20, 12, 12`.
* The quadrupled pentagonal cosine `q1 = 4 cos(2 pi / 5) = sqrt 5 - 1` is
  pinned by its defining identity `q1^2 + 2 q1 - 4 = 0` with positive
  radical part, and `q2 = 4 cos(4 pi / 5)` by the double-angle identity
  `2 q2 = q1^2 - 8`; both are checked below.
* The rotation character at level `l` is built by the closed-form
  unrolling `chi_l = 1 + 2 (cos theta + ... + cos(l theta))`; quadrupled,
  `4 chi_l = 4 + 2 (q(theta) + ... + q(l theta))`.
* The doubled `A5` character table is verified orthogonal with the
  correct norms (`sum |class| (2 chi_V)(2 chi_W) = 240 delta`) before the
  projection is used.

CONTENT.

* `branching_table`: `sum |class| (4 chi_l)(2 chi_V) = 480 m(l, V)` with
  `m` the frozen multiplicity table, including `l = 2 -> 5`,
  `l = 3 -> 3' + 4`, `l = 4 -> 4 + 5`, `l = 5 -> 3 + 3' + 5`, and
  `l = 6 -> 1 + 3 + 4 + 5`.
* `first_nonconstant_invariant_at_six`: levels one through five carry no
  trivial summand and level six carries exactly one.
* `face_phase_multiplicities`: `6 m_omega(V) = 2 dim V - 2 chi_V(C3)` on
  `(1, 3, 3', 4, 5)` equals `(0, 6, 6, 6, 12)`, the multiplicity vector
  `(0, 1, 1, 1, 2)`: the dimension-minimal irreducible extension of a
  nontrivial face phase is three-dimensional.

BOUNDARY.  Frequencies, amplitudes, linewidths, and the screen-to-
observable coupling are dynamical inputs; this file is the exact
representation arithmetic of the frozen signature only. -/

/-- `a + b * sqrt 5` as a pair of integers. -/
abbrev Z5 := ℤ × ℤ

def zAdd (x y : Z5) : Z5 := (x.1 + y.1, x.2 + y.2)

def zMul (x y : Z5) : Z5 := (x.1 * y.1 + 5 * x.2 * y.2, x.1 * y.2 + x.2 * y.1)

def zSmul (c : ℤ) (x : Z5) : Z5 := (c * x.1, c * x.2)

def zSum (xs : List Z5) : Z5 := xs.foldl zAdd (0, 0)

/-- `q1 = 4 cos(2 pi / 5) = sqrt 5 - 1`. -/
def q1 : Z5 := (-1, 1)

/-- `q2 = 4 cos(4 pi / 5) = -(sqrt 5 + 1)`. -/
def q2 : Z5 := (-1, -1)

/-- The quadrupled pentagonal cosine satisfies its defining identity
`q1^2 + 2 q1 - 4 = 0`. -/
theorem q1_defining_identity :
    zAdd (zAdd (zMul q1 q1) (zSmul 2 q1)) (-4, 0) = (0, 0) := by decide

/-- The double-angle identity `2 q2 = q1^2 - 8` pins the second cosine. -/
theorem q2_double_angle : zSmul 2 q2 = zAdd (zMul q1 q1) (-8, 0) := by decide

/-- `4 cos(m theta)` on the five classes: identity, `pi`, `2 pi / 3`,
`2 pi / 5`, `4 pi / 5`. -/
def cos4 (m : ℕ) (k : Fin 5) : Z5 :=
  match k.val with
  | 0 => (4, 0)
  | 1 => if m % 2 = 0 then (4, 0) else (-4, 0)
  | 2 => if m % 3 = 0 then (4, 0) else (-2, 0)
  | 3 =>
      match m % 5 with
      | 0 => (4, 0)
      | 1 => q1
      | 4 => q1
      | _ => q2
  | _ =>
      match m % 5 with
      | 0 => (4, 0)
      | 2 => q1
      | 3 => q1
      | _ => q2

/-- The quadrupled rotation character at level `l` on class `k`:
`4 chi_l = 4 + 2 (4 cos theta + ... + 4 cos(l theta))`. -/
def chi4 (l : ℕ) (k : Fin 5) : Z5 :=
  zAdd (4, 0)
    (zSmul 2 (zSum ((List.range l).map fun i => cos4 (i + 1) k)))

/-- The doubled `A5` character table on classes
(identity, `C2`, `C3`, `C5`, `C5^2`), rows `(1, 3, 3', 4, 5)`. -/
def chi2 (v : Fin 5) (k : Fin 5) : Z5 :=
  match v.val, k.val with
  | 0, _ => (2, 0)
  | 1, 0 => (6, 0)
  | 1, 1 => (-2, 0)
  | 1, 2 => (0, 0)
  | 1, 3 => (1, 1)
  | 1, _ => (1, -1)
  | 2, 0 => (6, 0)
  | 2, 1 => (-2, 0)
  | 2, 2 => (0, 0)
  | 2, 3 => (1, -1)
  | 2, _ => (1, 1)
  | 3, 0 => (8, 0)
  | 3, 1 => (0, 0)
  | 3, 2 => (2, 0)
  | 3, _ => (-2, 0)
  | _, 0 => (10, 0)
  | _, 1 => (2, 0)
  | _, 2 => (-2, 0)
  | _, _ => (0, 0)

def classSize (k : Fin 5) : ℤ :=
  match k.val with
  | 0 => 1
  | 1 => 15
  | 2 => 20
  | 3 => 12
  | _ => 12

/-- `sum |class| f(class) g(class)` (all characters here are real). -/
def pairingSum (f g : Fin 5 → Z5) : Z5 :=
  zSum ((List.finRange 5).map fun k => zSmul (classSize k) (zMul (f k) (g k)))

/-- The doubled character table is orthogonal with the correct norms:
`sum |class| (2 chi_V)(2 chi_W) = 240 delta`, so the projection below is
exact. -/
theorem character_table_orthonormal :
    ∀ v w : Fin 5,
      pairingSum (chi2 v) (chi2 w) =
        (if v = w then ((240 : ℤ), (0 : ℤ)) else ((0 : ℤ), (0 : ℤ))) := by
  decide

/-- The frozen branching table for levels zero through six, rows
`(1, 3, 3', 4, 5)`. -/
def frozenTable (l : ℕ) (v : Fin 5) : ℤ :=
  match l, v.val with
  | 0, 0 => 1
  | 1, 1 => 1
  | 2, 4 => 1
  | 3, 2 => 1
  | 3, 3 => 1
  | 4, 3 => 1
  | 4, 4 => 1
  | 5, 1 => 1
  | 5, 2 => 1
  | 5, 4 => 1
  | 6, 0 => 1
  | 6, 1 => 1
  | 6, 3 => 1
  | 6, 4 => 1
  | _, _ => 0

/-- The exact branching equals the frozen table:
`sum |class| (4 chi_l)(2 chi_V) = 480 m(l, V)`.  In particular `l = 2` is
the irreducible quintet, `l = 3` is `3' + 4`, `l = 4` is `4 + 5`, `l = 5`
is `3 + 3' + 5`, and `l = 6` is `1 + 3 + 4 + 5`. -/
theorem branching_table :
    ∀ l : Fin 7, ∀ v : Fin 5,
      pairingSum (chi4 l.val) (chi2 v) = (480 * frozenTable l.val v, 0) := by
  decide

/-- Levels one through five carry no trivial summand; level six carries
exactly one: the first nonconstant invariant sits at `l = 6`. -/
theorem first_nonconstant_invariant_at_six :
    (∀ l : Fin 7, 0 < l.val → l.val < 6 →
        pairingSum (chi4 l.val) (chi2 ⟨0, by omega⟩) = (0, 0)) ∧
      pairingSum (chi4 6) (chi2 ⟨0, by omega⟩) = (480, 0) := by
  decide

/-- The sextupled face-phase multiplicity `6 m_omega(V) = 2 dim V -
2 chi_V(C3)` on `(1, 3, 3', 4, 5)` is `(0, 6, 6, 6, 12)`: the multiplicity
vector is `(0, 1, 1, 1, 2)` and the dimension-minimal irreducible
extension of a nontrivial face phase is three-dimensional. -/
theorem face_phase_multiplicities :
    (List.finRange 5).map (fun v =>
        zAdd (chi2 v ⟨0, by omega⟩) (zSmul (-1) (chi2 v ⟨2, by omega⟩))) =
      [((0 : ℤ), (0 : ℤ)), (6, 0), (6, 0), (6, 0), (12, 0)] := by
  decide

end OPH.A5AngularMultiplets

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.A5AngularMultiplets.q1_defining_identity
#print axioms OPH.A5AngularMultiplets.q2_double_angle
#print axioms OPH.A5AngularMultiplets.character_table_orthonormal
#print axioms OPH.A5AngularMultiplets.branching_table
#print axioms OPH.A5AngularMultiplets.first_nonconstant_invariant_at_six
#print axioms OPH.A5AngularMultiplets.face_phase_multiplicities
