import Mathlib

/-!
# Canonical oriented-face bracket discriminator

This is the proof-kernel companion to
`code/b14_jacobi/oriented_face_bracket_selector.certificate.json`.  It uses the
port order of the pinned echosahedral manifest, which is also the port frame of
the stage-1 Reynolds basis.  (The older `CoreAxioms.orientedFaces` literal is an
isomorphic but different port relabeling.)

The pinned incidence, together with the declared equal-weight cyclic rule,
gives a bracket with exact identity `B_face = 60 R13`.  The bracket fails Jacobi.  Exact certificate
replay finds the same unique nearest compact family `G` for total-absolute,
Hilbert--Schmidt, and worst-coordinate edits; this file kernel-checks all three
exact distance orderings.  **No one of those repair norms, minimum-repair
rules, Jacobi-repair laws, or physical source selections is derived here.**
-/

namespace OPH.OrientedFaceBracketSelector

structure Face where
  a : Fin 12
  b : Fin 12
  c : Fin 12
deriving DecidableEq, Repr

/-- The twenty source-manifest faces in the pinned Reynolds port frame. -/
def sourceFaces : List Face :=
  [⟨0, 2, 8⟩, ⟨0, 10, 2⟩, ⟨0, 4, 6⟩, ⟨0, 8, 4⟩, ⟨0, 6, 10⟩,
   ⟨1, 9, 3⟩, ⟨1, 3, 11⟩, ⟨1, 6, 4⟩, ⟨1, 4, 9⟩, ⟨1, 11, 6⟩,
   ⟨2, 7, 5⟩, ⟨2, 5, 8⟩, ⟨2, 10, 7⟩, ⟨3, 5, 7⟩, ⟨3, 9, 5⟩,
   ⟨3, 7, 11⟩, ⟨4, 8, 9⟩, ⟨5, 9, 8⟩, ⟨6, 11, 10⟩, ⟨7, 10, 11⟩]

theorem source_face_count : sourceFaces.length = 20 := by decide +kernel

/-- Positive structure constants contributed by the three cyclic orders. -/
def positiveEntries : List (Fin 12 × Fin 12 × Fin 12) :=
  sourceFaces.flatMap fun f =>
    [(f.c, f.a, f.b), (f.a, f.b, f.c), (f.b, f.c, f.a)]

/-- Integral coefficient of `e_out` in `[e_left,e_right]`. -/
def faceCoeff (out left right : Fin 12) : ℤ :=
  (positiveEntries.count (out, left, right) : ℤ) -
    (positiveEntries.count (out, right, left) : ℤ)

/-- One upper-triangular numerator in the independently serialized `R13` row. -/
structure R13Entry where
  out : Fin 12
  left : Fin 12
  right : Fin 12
  numerator : ℤ
deriving DecidableEq, Repr

/-- The literal nonzero upper-triangular numerator support of pinned `R13`. -/
def r13UpperEntries : List R13Entry :=
  [⟨0,2,8,1⟩, ⟨0,2,10,-1⟩, ⟨0,4,6,1⟩, ⟨0,4,8,-1⟩, ⟨0,6,10,1⟩,
   ⟨1,3,9,-1⟩, ⟨1,3,11,1⟩, ⟨1,4,6,-1⟩, ⟨1,4,9,1⟩, ⟨1,6,11,-1⟩,
   ⟨2,0,8,-1⟩, ⟨2,0,10,1⟩, ⟨2,5,7,-1⟩, ⟨2,5,8,1⟩, ⟨2,7,10,-1⟩,
   ⟨3,1,9,1⟩, ⟨3,1,11,-1⟩, ⟨3,5,7,1⟩, ⟨3,5,9,-1⟩, ⟨3,7,11,1⟩,
   ⟨4,0,6,-1⟩, ⟨4,0,8,1⟩, ⟨4,1,6,1⟩, ⟨4,1,9,-1⟩, ⟨4,8,9,1⟩,
   ⟨5,2,7,1⟩, ⟨5,2,8,-1⟩, ⟨5,3,7,-1⟩, ⟨5,3,9,1⟩, ⟨5,8,9,-1⟩,
   ⟨6,0,4,1⟩, ⟨6,0,10,-1⟩, ⟨6,1,4,-1⟩, ⟨6,1,11,1⟩, ⟨6,10,11,-1⟩,
   ⟨7,2,5,-1⟩, ⟨7,2,10,1⟩, ⟨7,3,5,1⟩, ⟨7,3,11,-1⟩, ⟨7,10,11,1⟩,
   ⟨8,0,2,1⟩, ⟨8,0,4,-1⟩, ⟨8,2,5,1⟩, ⟨8,4,9,-1⟩, ⟨8,5,9,1⟩,
   ⟨9,1,3,-1⟩, ⟨9,1,4,1⟩, ⟨9,3,5,-1⟩, ⟨9,4,8,1⟩, ⟨9,5,8,-1⟩,
   ⟨10,0,2,-1⟩, ⟨10,0,6,1⟩, ⟨10,2,7,-1⟩, ⟨10,6,11,1⟩, ⟨10,7,11,-1⟩,
   ⟨11,1,3,1⟩, ⟨11,1,6,-1⟩, ⟨11,3,7,1⟩, ⟨11,6,10,-1⟩, ⟨11,7,10,1⟩]

def r13UpperNumerator (out left right : Fin 12) : ℤ :=
  match r13UpperEntries.find? (fun e =>
      decide (e.out = out ∧ e.left = left ∧ e.right = right)) with
  | some e => e.numerator
  | none => 0

/-- Antisymmetric extension of the pinned `R13` row. -/
def r13Coeff (out left right : Fin 12) : ℚ :=
  if left < right then (r13UpperNumerator out left right : ℚ) / 60
  else if right < left then -(r13UpperNumerator out right left : ℚ) / 60
  else 0

/-- Exact tensor-coordinate identity, not merely equality modulo an orbit. -/
theorem face_bracket_eq_sixty_r13 :
    ∀ out left right : Fin 12,
      (faceCoeff out left right : ℚ) = 60 * r13Coeff out left right := by
  decide +kernel

/-- One coordinate of the Jacobiator on basis vectors. -/
def jacobiCoeff (out i j k : Fin 12) : ℤ :=
  ∑ m : Fin 12,
    (faceCoeff m i j * faceCoeff out m k +
     faceCoeff m j k * faceCoeff out m i +
     faceCoeff m k i * faceCoeff out m j)

/-- A concrete failure witness: this face bracket is not a Lie bracket. -/
theorem jacobi_failure_witness : jacobiCoeff 0 1 4 10 = -1 := by decide +kernel

/-
The complete 2640-coordinate census (240 nonzero entries, split 120/120 with
unit magnitude and squared norm 240) is carried by the independently replayed
Python certificate.  We intentionally do not smuggle that bounded enumeration
into Lean through a native-evaluation trust shortcut; this module kernel-checks the concrete
failure witness and the all-coordinate `60 R13` identity.
-/

/-! ## Conditional nearest-compact-locus discriminator -/

inductive CompactFamily
  | P | F | G
deriving DecidableEq, Repr

noncomputable def distanceP : ℝ := 45
noncomputable def distanceF : ℝ := (615 + 123 * Real.sqrt 5) / 22
noncomputable def distanceG : ℝ := (615 - 123 * Real.sqrt 5) / 22

noncomputable def squaredDistance : CompactFamily → ℝ
  | .P => distanceP
  | .F => distanceF
  | .G => distanceG

/-- Exact orthogonal projection parameters in `(p,m)` order. -/
noncomputable def projectionP : ℝ × ℝ :=
  (Real.sqrt 5 / 20, -Real.sqrt 5 / 20)

/-- Exact orthogonal projection parameters in `(a,b,e)` order. -/
noncomputable def projectionF : ℝ × ℝ × ℝ :=
  (-1 / 22 + Real.sqrt 5 / 220, Real.sqrt 5 / 20, -Real.sqrt 5 / 20)

/-- Exact orthogonal projection parameters in `(a,b,e)` order. -/
noncomputable def projectionG : ℝ × ℝ × ℝ :=
  (-1 / 22 - Real.sqrt 5 / 220, -Real.sqrt 5 / 20, Real.sqrt 5 / 20)

/-- The `F` projection obeys its compact-stratum sign `a*b < 0`. -/
theorem projectionF_compact_sign : projectionF.1 * projectionF.2.1 < 0 := by
  have hs : 0 < Real.sqrt 5 := by positivity
  have hs2 : Real.sqrt 5 ^ 2 = 5 := by norm_num
  simp only [projectionF]
  nlinarith

/-- The `G` projection obeys its compact-stratum sign `a*b > 0`. -/
theorem projectionG_compact_sign : 0 < projectionG.1 * projectionG.2.1 := by
  have hs : 0 < Real.sqrt 5 := by positivity
  have hs2 : Real.sqrt 5 ^ 2 = 5 := by norm_num
  simp only [projectionG]
  nlinarith

theorem exact_distance_gaps :
    distanceF - distanceG = 123 * Real.sqrt 5 / 11 ∧
      distanceP - distanceG = (375 + 123 * Real.sqrt 5) / 22 := by
  constructor <;> simp [distanceF, distanceG, distanceP] <;> ring

theorem distanceG_lt_distanceF : distanceG < distanceF := by
  have hs : 0 < Real.sqrt 5 := by positivity
  rw [show distanceF = distanceG + 123 * Real.sqrt 5 / 11 by
    linarith [exact_distance_gaps.1]]
  linarith

theorem distanceG_lt_distanceP : distanceG < distanceP := by
  have hs : 0 ≤ Real.sqrt 5 := Real.sqrt_nonneg 5
  dsimp [distanceG, distanceP]
  nlinarith

/-- `G` is the unique nearest member of `{P,F,G}` under the added HS metric. -/
theorem unique_nearest_G (family : CompactFamily) (h : family ≠ .G) :
    squaredDistance .G < squaredDistance family := by
  cases family with
  | P => exact distanceG_lt_distanceP
  | F => exact distanceG_lt_distanceF
  | G => exact False.elim (h rfl)

/-! ## Exact endpoint-norm robustness

The independently replayed certificate supplies exact primal/dual witnesses
for the `L1` and `Linfinity` distances below.  For family `F`, the displayed
`L1` value is the compact-stratum infimum: its linear-family optimizer is zero
and compact points scale to zero.  The `G` optimizer lies strictly inside its
compact sign stratum.  These declarations check the radical arithmetic and
the common family ordering; they do not replace the certificate's 792-entry
primal/dual replay.
-/

inductive EditNorm
  | totalAbsolute
  | hilbertSchmidtSquared
  | worstCoordinate
deriving DecidableEq, Repr

noncomputable def totalAbsoluteDistance : CompactFamily → ℝ
  | .P => 60
  | .F => 60
  | .G => 30 * (Real.sqrt 5 - 1)

noncomputable def worstCoordinateDistance : CompactFamily → ℝ
  | .P => 1 / 2
  | .F => Real.sqrt 5 / 5
  | .G => (5 - Real.sqrt 5) / 10

noncomputable def certifiedDistance : EditNorm → CompactFamily → ℝ
  | .totalAbsolute => totalAbsoluteDistance
  | .hilbertSchmidtSquared => squaredDistance
  | .worstCoordinate => worstCoordinateDistance

theorem sqrt_five_lt_three : Real.sqrt 5 < 3 := by
  have hs : 0 ≤ Real.sqrt 5 := Real.sqrt_nonneg 5
  have hs2 : Real.sqrt 5 ^ 2 = 5 := by norm_num
  nlinarith

theorem five_lt_three_sqrt_five : 5 < 3 * Real.sqrt 5 := by
  have hs : 0 ≤ Real.sqrt 5 := Real.sqrt_nonneg 5
  have hs2 : Real.sqrt 5 ^ 2 = 5 := by norm_num
  nlinarith

theorem two_sqrt_five_lt_five : 2 * Real.sqrt 5 < 5 := by
  have hs : 0 ≤ Real.sqrt 5 := Real.sqrt_nonneg 5
  have hs2 : Real.sqrt 5 ^ 2 = 5 := by norm_num
  nlinarith

theorem totalAbsoluteDistance_G_lt_P :
    totalAbsoluteDistance .G < totalAbsoluteDistance .P := by
  simp only [totalAbsoluteDistance]
  nlinarith [sqrt_five_lt_three]

theorem totalAbsoluteDistance_G_lt_F :
    totalAbsoluteDistance .G < totalAbsoluteDistance .F := by
  simpa [totalAbsoluteDistance] using totalAbsoluteDistance_G_lt_P

theorem worstCoordinateDistance_G_lt_F :
    worstCoordinateDistance .G < worstCoordinateDistance .F := by
  simp only [worstCoordinateDistance]
  nlinarith [five_lt_three_sqrt_five]

theorem worstCoordinateDistance_F_lt_P :
    worstCoordinateDistance .F < worstCoordinateDistance .P := by
  simp only [worstCoordinateDistance]
  nlinarith [two_sqrt_five_lt_five]

theorem worstCoordinateDistance_G_lt_P :
    worstCoordinateDistance .G < worstCoordinateDistance .P :=
  lt_trans worstCoordinateDistance_G_lt_F worstCoordinateDistance_F_lt_P

/-- All three exactly certified coordinate-edit rules distinguish the same
compact family `G`.  This is norm robustness, not derivation of a repair rule. -/
theorem three_norm_unique_nearest_G
    (norm : EditNorm) (family : CompactFamily) (h : family ≠ .G) :
    certifiedDistance norm .G < certifiedDistance norm family := by
  cases norm with
  | totalAbsolute =>
      cases family with
      | P => exact totalAbsoluteDistance_G_lt_P
      | F => exact totalAbsoluteDistance_G_lt_F
      | G => exact False.elim (h rfl)
  | hilbertSchmidtSquared =>
      cases family with
      | P => exact distanceG_lt_distanceP
      | F => exact distanceG_lt_distanceF
      | G => exact False.elim (h rfl)
  | worstCoordinate =>
      cases family with
      | P => exact worstCoordinateDistance_G_lt_P
      | F => exact worstCoordinateDistance_G_lt_F
      | G => exact False.elim (h rfl)

#print axioms face_bracket_eq_sixty_r13
#print axioms jacobi_failure_witness
#print axioms projectionF_compact_sign
#print axioms projectionG_compact_sign
#print axioms unique_nearest_G
#print axioms three_norm_unique_nearest_G

end OPH.OrientedFaceBracketSelector
