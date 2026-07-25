import Mathlib

/-!
# The W5 residual-stabilizer boundary

This module formalizes the coordinate algebra used by the conditional
`W5 = Sym²₀(ℝ³)` flavor-orbit receipt.  A form is recorded by its six
symmetric entries, together with the tracelessness equation.

For a rotation about the third axis, invariance under `A ↦ R A Rᵀ` is
equivalent to commutation with the orthogonal rotation.  The predicate
`AxialCommutes` is the coordinate expansion of that commutation relation.
If the planar sine is nonzero, as it is for the canonical order-three and
order-five generators, every commuting form is

`diag(t, t, -2t)`.

Thus the fixed locus is one-dimensional and the two planar diagonal entries
are repeated.  For the canonical half-turn, the fixed forms instead have the
three-coordinate form

`[[a,b,0],[b,d,0],[0,0,-a-d]]`.

An explicit equivalence with `ℝ × ℝ × ℝ` certifies those three free
coordinates, and `diag(1,2,-3)` is a fixed diagonal witness with pairwise
distinct entries.

These are stabilizer facts only.  This module does not construct an
`A5`-invariant potential, select a physical orbit, or state a universal
impossibility theorem for flavor models.
-/

namespace OPH.W5Stabilizer

/-- Coordinates on a real symmetric traceless three-by-three form. -/
@[ext]
structure Form where
  xx : ℝ
  yy : ℝ
  zz : ℝ
  xy : ℝ
  xz : ℝ
  yz : ℝ
  trace_zero : xx + yy + zz = 0

/-- A proper axial rotation whose angle is neither zero nor a half-turn.
Canonical generators of `C3` and `C5` satisfy this contract. -/
structure NonHalfTurnAxialRotation where
  cosine : ℝ
  sine : ℝ
  unit : cosine ^ 2 + sine ^ 2 = 1
  sine_ne_zero : sine ≠ 0

/-- Coordinate expansion of `A R = R A` for a symmetric form and a rotation
about the third axis.  Orthogonality (`rotation.unit`) makes this equivalent
to the fixed-point equation `R A Rᵀ = A`. -/
def AxialCommutes (rotation : NonHalfTurnAxialRotation) (A : Form) : Prop :=
  rotation.cosine * A.xx + rotation.sine * A.xy =
      rotation.cosine * A.xx - rotation.sine * A.xy
    ∧ -rotation.sine * A.xx + rotation.cosine * A.xy =
      rotation.cosine * A.xy - rotation.sine * A.yy
    ∧ A.xz =
      rotation.cosine * A.xz - rotation.sine * A.yz
    ∧ rotation.cosine * A.xz + rotation.sine * A.yz =
      A.xz
    ∧ A.yz =
      rotation.sine * A.xz + rotation.cosine * A.yz
    ∧ -rotation.sine * A.xz + rotation.cosine * A.yz =
      A.yz

/-- The one-parameter axial normal form.  Its first two diagonal entries are
equal, so as a diagonal symmetric matrix it has a repeated eigenvalue. -/
def axialForm (t : ℝ) : Form where
  xx := t
  yy := t
  zz := -2 * t
  xy := 0
  xz := 0
  yz := 0
  trace_zero := by ring

/-- Every axial normal form commutes with every declared proper axial
rotation. -/
theorem axialForm_commutes
    (rotation : NonHalfTurnAxialRotation)
    (t : ℝ) :
    AxialCommutes rotation (axialForm t) := by
  simp [AxialCommutes, axialForm]

/-- The fixed-form equations for any non-half-turn axial generator force all
off-diagonal coordinates to vanish and force a repeated planar diagonal
entry.  Tracelessness fixes the remaining entry. -/
theorem nonHalfTurnAxial_fixed_normal_form
    (rotation : NonHalfTurnAxialRotation)
    (A : Form)
    (hcomm : AxialCommutes rotation A) :
    A.xy = 0
      ∧ A.xz = 0
      ∧ A.yz = 0
      ∧ A.yy = A.xx
      ∧ A.zz = -2 * A.xx := by
  rcases hcomm with ⟨hxx, hxy, hxzColumn, hxzRow, hyzColumn, hyzRow⟩
  have hxyProduct : rotation.sine * A.xy = 0 := by
    linarith
  have hxyZero : A.xy = 0 :=
    (mul_eq_zero.mp hxyProduct).resolve_left rotation.sine_ne_zero
  have hdiagProduct : rotation.sine * (A.xx - A.yy) = 0 := by
    nlinarith
  have hdiag : A.xx = A.yy := by
    have : A.xx - A.yy = 0 :=
      (mul_eq_zero.mp hdiagProduct).resolve_left rotation.sine_ne_zero
    linarith
  have hyzProduct : rotation.sine * A.yz = 0 := by
    linarith
  have hyzZero : A.yz = 0 :=
    (mul_eq_zero.mp hyzProduct).resolve_left rotation.sine_ne_zero
  have hxzProduct : rotation.sine * A.xz = 0 := by
    linarith
  have hxzZero : A.xz = 0 :=
    (mul_eq_zero.mp hxzProduct).resolve_left rotation.sine_ne_zero
  constructor
  · exact hxyZero
  constructor
  · exact hxzZero
  constructor
  · exact hyzZero
  constructor
  · exact hdiag.symm
  · linarith [A.trace_zero]

/-- Hence every fixed form is exactly the one-parameter axial form.  This is
the algebraic fixed-locus statement used for the `C3` and `C5` rows of the
machine receipt. -/
theorem nonHalfTurnAxial_fixed_eq_form
    (rotation : NonHalfTurnAxialRotation)
    (A : Form)
    (hcomm : AxialCommutes rotation A) :
    A = axialForm A.xx := by
  rcases nonHalfTurnAxial_fixed_normal_form rotation A hcomm with
    ⟨hxy, hxz, hyz, hyy, hzz⟩
  apply Form.ext
  · rfl
  · simpa [axialForm] using hyy
  · simpa [axialForm] using hzz
  · simpa [axialForm] using hxy
  · simpa [axialForm] using hxz
  · simpa [axialForm] using hyz

/-- The axial fixed locus is explicitly equivalent to one real coordinate. -/
def axialFixedEquiv
    (rotation : NonHalfTurnAxialRotation) :
    {A : Form // AxialCommutes rotation A} ≃ ℝ where
  toFun A := A.1.xx
  invFun t := ⟨axialForm t, axialForm_commutes rotation t⟩
  left_inv A := by
    apply Subtype.ext
    exact (nonHalfTurnAxial_fixed_eq_form rotation A.1 A.2).symm
  right_inv _ := rfl

/-- Fixed-point equations for the canonical half-turn
`diag(-1,-1,1)`.  Exactly the two transverse entries change sign. -/
def C2Fixed (A : Form) : Prop :=
  A.xz = -A.xz ∧ A.yz = -A.yz

/-- The half-turn fixed condition is exactly vanishing of the two transverse
coordinates. -/
theorem c2_fixed_iff (A : Form) :
    C2Fixed A ↔ A.xz = 0 ∧ A.yz = 0 := by
  constructor
  · rintro ⟨hxz, hyz⟩
    constructor <;> linarith
  · rintro ⟨hxz, hyz⟩
    simp [C2Fixed, hxz, hyz]

/-- General three-coordinate half-turn fixed form. -/
def c2Form (a b d : ℝ) : Form where
  xx := a
  yy := d
  zz := -a - d
  xy := b
  xz := 0
  yz := 0
  trace_zero := by ring

/-- Every declared three-coordinate form is half-turn fixed. -/
theorem c2Form_fixed (a b d : ℝ) :
    C2Fixed (c2Form a b d) := by
  simp [C2Fixed, c2Form]

/-- Every half-turn fixed form has the declared three-coordinate
parameterization. -/
theorem c2_fixed_eq_form
    (A : Form)
    (hfixed : C2Fixed A) :
    A = c2Form A.xx A.xy A.yy := by
  have hzero := (c2_fixed_iff A).mp hfixed
  apply Form.ext
  · rfl
  · rfl
  · simp only [c2Form]
    linarith [A.trace_zero]
  · rfl
  · simpa [c2Form] using hzero.1
  · simpa [c2Form] using hzero.2

/-- The `C2` fixed locus is explicitly equivalent to three real
coordinates.  Projectivizing by nonzero overall scale therefore leaves the
two continuous ratios described in the paper receipt; that quotient is not
constructed here. -/
def c2FixedEquiv : {A : Form // C2Fixed A} ≃ ℝ × ℝ × ℝ where
  toFun A := (A.1.xx, A.1.xy, A.1.yy)
  invFun coordinates :=
    ⟨c2Form coordinates.1 coordinates.2.1 coordinates.2.2,
      c2Form_fixed coordinates.1 coordinates.2.1 coordinates.2.2⟩
  left_inv A := by
    apply Subtype.ext
    exact (c2_fixed_eq_form A.1 A.2).symm
  right_inv coordinates := by
    rcases coordinates with ⟨a, b, d⟩
    rfl

/-- The machine receipt's simple-spectrum `C2` witness. -/
def c2SimpleWitness : Form :=
  c2Form 1 0 2

/-- The explicit witness is half-turn fixed. -/
theorem c2SimpleWitness_fixed :
    C2Fixed c2SimpleWitness :=
  c2Form_fixed 1 0 2

/-- The witness is diagonal. -/
theorem c2SimpleWitness_diagonal :
    c2SimpleWitness.xy = 0
      ∧ c2SimpleWitness.xz = 0
      ∧ c2SimpleWitness.yz = 0 := by
  norm_num [c2SimpleWitness, c2Form]

/-- Its diagonal entries are exactly the receipt values. -/
theorem c2SimpleWitness_entries :
    c2SimpleWitness.xx = 1
      ∧ c2SimpleWitness.yy = 2
      ∧ c2SimpleWitness.zz = -3 := by
  norm_num [c2SimpleWitness, c2Form]

/-- The diagonal entries, hence the eigenvalues of this diagonal symmetric
matrix, are pairwise distinct. -/
theorem c2SimpleWitness_pairwise_distinct :
    c2SimpleWitness.xx ≠ c2SimpleWitness.yy
      ∧ c2SimpleWitness.xx ≠ c2SimpleWitness.zz
      ∧ c2SimpleWitness.yy ≠ c2SimpleWitness.zz := by
  norm_num [c2SimpleWitness, c2Form]

end OPH.W5Stabilizer
