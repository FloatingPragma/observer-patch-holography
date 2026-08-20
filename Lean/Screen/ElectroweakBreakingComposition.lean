import SMStructureAdequacySurface
import MatterGrammarIndexBridge

/-!
# The composed electroweak-breaking/Yukawa inhabitant (V3, issue #735)

One typed premise bundle, one composed receipt.  Register rows PR-67,
PR-68, and PR-69 declare the finite scalar grammar, quartic breaking
packet, and single Yukawa-line coefficient used here.  The broader
physical scalar, vacuum, Yukawa-matrix, and source-action attachments
(PR-48, PR-49, PR-50, PR-54) remain open.  This module inhabits only
the finite clauses on the committed coarse grammar: it declares the
premise parameters as typed fields of one bundle
`EWBreakingPremiseData` (extending the committed
`SMStructurePremiseData` of register rows PR-09, PR-10, PR-11,
PR-59), builds the scalar doublet carrier on the quantum numbers of
the committed exterior component table, and proves the composed
symmetry-breaking receipt from the one bundle.

The scalar carrier is `Fin 2 → ℂ`: its complex dimension is the
committed weak dimension `weakDim 1 = 2` of the unique colorless
weak-doublet row of the exterior table (register row PR-59), and its
declared abelian charge is the committed integer charge
`charge 1 = 3` of that row in `q = 6Y` units, so the hypercharge
reading is `Y = 1/2`, exactly the PR-48 candidate identification.
The quartic potential is `V h = lam * (|h|^2 - vev^2)^2` with `lam`
and `vev` declared positive bundle parameters under PR-68.  Here
`vev` is the model-shell radius and the chosen vector is `(0, vev)`.
If it is compared with the conventional vector
`vSM / sqrt 2 * (0, 1)`, then `vev = vSM / sqrt 2`; no physical
kinetic or mass normalization is inferred.  The proved clauses:

* the global minimizer set of `V` is exactly the nonzero shell
  `|h|^2 = vev^2`;
* the stabilizer of the chosen minimum `(0, vev)` inside the
  four-parameter electroweak direction space is exactly one
  one-parameter unbroken line, the mix of the third weak direction
  and the hypercharge direction with weights `(gY, gW)`;
* the second variation of `V` at the minimum is an exact quadratic
  coefficient of a polynomial identity, positive semidefinite, with
  kernel exactly the broken directions: the kernel coincides with the
  set of infinitesimal gauge images at the minimum and is exactly
  three-parameter (the Goldstone count);
* the gauge mass form (the exact norm-square of the infinitesimal
  action at the minimum, the coefficient that a gauge-covariant
  kinetic term assigns to constant gauge directions when register row
  PR-54 supplies such a term) vanishes exactly on the unbroken line,
  takes the value `gW^2 vev^2 / 4` on the charged direction and
  `(gW^2 + gY^2) vev^2 / 4` on the normalized neutral broken
  direction, and the tree ratio of the two is
  `gW^2 / (gW^2 + gY^2)`, the cosine-squared shape in the two
  declared coupling constants;
* one Yukawa line through the committed grammar: the epsilon-paired
  coupling of the scalar to a doublet carrier and a singlet produces
  at the minimum a bilinear with the exact coefficient `yuk * vev`,
  and that coefficient vanishes exactly when `yuk = 0` or `vev = 0`;
  the bundled selection mask supplies, in either committed parity
  sector, a colorless doublet row and a colorless singlet row of
  matching parity grading and matching Weyl column (through the
  committed index dictionary of `MatterGrammarIndexBridge`) whose
  charges balance against the scalar row charge or its conjugate.

The committed instance `committedEWBreakingPremiseData` (all declared
parameters set to one over the committed even-sector base instance)
proves the bundle nonvacuous.

Boundary and nonclaims.  This is a conditional structural composition
on the committed coarse grammar of register row PR-59.  It is not the
physical Higgs sector.  The couplings `gW`, `gY`, the quartic `lam`,
the scale `vev`, and the Yukawa coefficient `yuk` are declared bundle
parameters under PR-67--PR-69: no committed theorem selects them, and the committed
invariant-form receipts of `GaugeKineticInvariantForms` leave the
relative coupling normalization free, so the mass ratio is stated in
declared constants and is not a derived weak mixing angle.  No
measured mass, width, or angle appears.  There is no spacetime, no
covariant derivative, and no source-produced action: register row
PR-54 stays open, and the gauge mass form enters as the exact
algebraic coefficient defined here, not as a term of a constructed
Lagrangian.  The scalar carrier is built on the quantum numbers of
the committed table; its physical existence and kinetic normalization
stay on register row PR-48.  Physical vacuum selection stays on
PR-49.  The Yukawa line is one line with one
declared coefficient; no Yukawa matrix, family structure, texture, or
phase is claimed (register row PR-50 and its non-identifiability
receipts are untouched, and register row PR-36 three-family structure
is not consumed).  The parity grading is table provenance; a physical
chirality claim stays on register row PR-47.  OL-G7 may move from
owed to conditional structural partial at most.
-/

namespace OPH.ElectroweakBreakingComposition

open OPH.SMStructureAdequacySurface (SMStructurePremiseData
  committedSMStructurePremiseData)

/-! ## The scalar doublet carrier on the committed grammar -/

/-- The scalar doublet carrier: complex dimension two, matching the
committed weak dimension of the colorless weak-doublet row of the
exterior component table (register row PR-59). -/
abbrev Doublet := Fin 2 → ℂ

/-- The declared abelian charge of the scalar carrier in the
committed `q = 6Y` units: the integer charge of the colorless
weak-doublet row of the committed exterior table.  Its value is
three, so the hypercharge reading is `1/2` (the PR-48 candidate
identification). -/
def scalarChargeSix : ℤ := OPH.ExteriorSelection.charge 1

/-- The committed-table receipt behind the scalar carrier: the row
`1` of the exterior component table is the colorless weak doublet
with integer charge three, and the carrier dimension of `Doublet`
matches its committed weak dimension. -/
theorem scalar_row_receipt :
    OPH.ExteriorSelection.colorDim 1 = 1
      ∧ OPH.ExteriorSelection.weakDim 1 = 2
      ∧ OPH.ExteriorSelection.isColored 1 = false
      ∧ OPH.ExteriorSelection.isDoublet 1 = true
      ∧ scalarChargeSix = 3
      ∧ (Fintype.card (Fin 2) : ℤ) = OPH.ExteriorSelection.weakDim 1 := by
  refine ⟨by decide, by decide, by decide, by decide, by decide, by decide⟩

/-- The squared Hermitian length of a doublet: the sum of the two
componentwise norm squares. -/
def normSqSum (h : Doublet) : ℝ :=
  Complex.normSq (h 0) + Complex.normSq (h 1)

/-- The chosen minimum: the lower-component real vacuum
configuration `(0, vev)`. -/
noncomputable def vacuum (vev : ℝ) : Doublet := ![0, (vev : ℂ)]

theorem normSqSum_nonneg (h : Doublet) : 0 ≤ normSqSum h :=
  add_nonneg (Complex.normSq_nonneg _) (Complex.normSq_nonneg _)

theorem normSqSum_eq_zero_iff (u : Doublet) : normSqSum u = 0 ↔ u = 0 := by
  constructor
  · intro h
    unfold normSqSum at h
    have h0 : Complex.normSq (u 0) = 0 := by
      nlinarith [Complex.normSq_nonneg (u 0), Complex.normSq_nonneg (u 1)]
    have h1 : Complex.normSq (u 1) = 0 := by
      nlinarith [Complex.normSq_nonneg (u 0), Complex.normSq_nonneg (u 1)]
    have e0 := Complex.normSq_eq_zero.mp h0
    have e1 := Complex.normSq_eq_zero.mp h1
    funext i
    fin_cases i <;> simp [e0, e1]
  · intro h
    rw [h]
    simp [normSqSum]

theorem vacuum_normSqSum (vev : ℝ) : normSqSum (vacuum vev) = vev ^ 2 := by
  simp [normSqSum, vacuum, Complex.normSq_ofReal]
  ring

/-! ## The quartic potential and its exact minimizer set (finite PR-68 packet) -/

/-- The declared finite quartic potential `lam * (|h|^2 - vev^2)^2`
(register row PR-68).  The physical vacuum attachment remains PR-49. -/
def potential (lam vev : ℝ) (h : Doublet) : ℝ :=
  lam * (normSqSum h - vev ^ 2) ^ 2

theorem potential_nonneg (lam vev : ℝ) (hlam : 0 < lam) (h : Doublet) :
    0 ≤ potential lam vev h :=
  mul_nonneg hlam.le (sq_nonneg _)

theorem potential_vacuum (lam vev : ℝ) :
    potential lam vev (vacuum vev) = 0 := by
  rw [potential, vacuum_normSqSum]
  ring

/-- **The exact minimizer set.**  Under positive declared `lam`, a
doublet is a global minimizer of the potential exactly when its
squared length is `vev^2`. -/
theorem minimizer_iff (lam vev : ℝ) (hlam : 0 < lam) (h : Doublet) :
    (∀ h' : Doublet, potential lam vev h ≤ potential lam vev h')
      ↔ normSqSum h = vev ^ 2 := by
  constructor
  · intro hmin
    have hle := hmin (vacuum vev)
    rw [potential_vacuum] at hle
    have hzero : potential lam vev h = 0 :=
      le_antisymm hle (potential_nonneg lam vev hlam h)
    rw [potential] at hzero
    rcases mul_eq_zero.mp hzero with h' | h'
    · exact absurd h' hlam.ne'
    · have := (pow_eq_zero_iff two_ne_zero).mp h'
      linarith [sub_eq_zero.mp this]
  · intro hnorm h'
    have : potential lam vev h = 0 := by
      rw [potential, hnorm]
      ring
    rw [this]
    exact potential_nonneg lam vev hlam h'

/-- Every minimizer is nonzero under a positive declared scale: the
minimizer set is the nonzero shell. -/
theorem minimizer_nonzero (vev : ℝ) (hvev : 0 < vev) (h : Doublet)
    (hh : normSqSum h = vev ^ 2) : h ≠ 0 := by
  intro h0
  rw [h0] at hh
  simp only [normSqSum, Pi.zero_apply, Complex.normSq_zero, add_zero] at hh
  have hpos : (0 : ℝ) < vev ^ 2 := by positivity
  exact hpos.ne' hh.symm

/-! ## The electroweak direction space and the infinitesimal action -/

/-- One electroweak direction: three weak parameters and one
hypercharge parameter.  The four-parameter space matches the
committed dimension data `1 + 3` of the electroweak label pair inside
the forced `(1, {3, 8})` type of the composed structure surface. -/
@[ext]
structure EWParam where
  /-- First weak direction. -/
  a1 : ℝ
  /-- Second weak direction. -/
  a2 : ℝ
  /-- Third weak direction. -/
  a3 : ℝ
  /-- Hypercharge direction (the committed scalar charge reading
  `Y = 1/2` is folded into the coefficient convention). -/
  b : ℝ

/-- The infinitesimal electroweak action on the doublet with declared
coupling constants `g` (weak) and `g'` (hypercharge): the generator
`i (g a·sigma/2 + g' b/2)` applied to `h`, written out in
components. -/
noncomputable def infAction (g g' : ℝ) (A : EWParam) (h : Doublet) :
    Doublet :=
  ![((g * A.a2 / 2 : ℝ) : ℂ) * h 1
      + Complex.I * (((g * A.a1 / 2 : ℝ) : ℂ) * h 1)
      + Complex.I * ((((g * A.a3 + g' * A.b) / 2 : ℝ) : ℂ) * h 0),
    -(((g * A.a2 / 2 : ℝ) : ℂ) * h 0)
      + Complex.I * (((g * A.a1 / 2 : ℝ) : ℂ) * h 0)
      + Complex.I * ((((g' * A.b - g * A.a3) / 2 : ℝ) : ℂ) * h 1)]

/-- The infinitesimal action at the chosen minimum in real/imaginary
normal form. -/
theorem infAction_vacuum_coords (g g' vev : ℝ) (A : EWParam) :
    infAction g g' A (vacuum vev) =
      ![((g * A.a2 * vev / 2 : ℝ) : ℂ)
          + ((g * A.a1 * vev / 2 : ℝ) : ℂ) * Complex.I,
        ((0 : ℝ) : ℂ)
          + (((g' * A.b - g * A.a3) * vev / 2 : ℝ) : ℂ) * Complex.I] := by
  funext i
  fin_cases i <;>
    · simp [infAction, vacuum]
      ring

/-! ## The gauge mass form and the exact stabilizer -/

/-- The gauge mass form: the exact norm square of the infinitesimal
action at the chosen minimum.  When register row PR-54 supplies a
gauge-covariant kinetic term, this is the exact coefficient that term
assigns to a constant gauge direction; here it is the committed
algebraic object itself. -/
noncomputable def gaugeMassForm (g g' vev : ℝ) (A : EWParam) : ℝ :=
  normSqSum (infAction g g' A (vacuum vev))

/-- The exact closed form of the gauge mass form in the declared
constants. -/
theorem gaugeMassForm_exact (g g' vev : ℝ) (A : EWParam) :
    gaugeMassForm g g' vev A
      = vev ^ 2 / 4
          * (g ^ 2 * (A.a1 ^ 2 + A.a2 ^ 2)
            + (g' * A.b - g * A.a3) ^ 2) := by
  unfold gaugeMassForm
  rw [infAction_vacuum_coords]
  simp only [normSqSum, Matrix.cons_val_zero, Matrix.cons_val_one,
    Complex.normSq_apply, Complex.add_re,
    Complex.add_im, Complex.mul_re, Complex.mul_im, Complex.I_re,
    Complex.I_im, Complex.ofReal_re, Complex.ofReal_im]
  ring

/-- The gauge mass form vanishes exactly on the stabilizer of the
minimum: mass terms appear exactly for the broken directions. -/
theorem gaugeMassForm_zero_iff (g g' vev : ℝ) (A : EWParam) :
    gaugeMassForm g g' vev A = 0
      ↔ infAction g g' A (vacuum vev) = 0 := by
  unfold gaugeMassForm
  exact normSqSum_eq_zero_iff _

/-- The kernel equations of the infinitesimal action at the minimum:
with positive declared couplings and scale, the action vanishes
exactly when both charged weak parameters vanish and the neutral
parameters satisfy the single unbroken relation. -/
theorem infAction_vacuum_zero_iff (g g' vev : ℝ) (hg : 0 < g)
    (hvev : 0 < vev) (A : EWParam) :
    infAction g g' A (vacuum vev) = 0
      ↔ A.a1 = 0 ∧ A.a2 = 0 ∧ g' * A.b = g * A.a3 := by
  constructor
  · intro h
    have hz : gaugeMassForm g g' vev A = 0 :=
      (gaugeMassForm_zero_iff g g' vev A).mpr h
    rw [gaugeMassForm_exact] at hz
    have hv4 : (0 : ℝ) < vev ^ 2 / 4 := by positivity
    have hsum : g ^ 2 * (A.a1 ^ 2 + A.a2 ^ 2)
        + (g' * A.b - g * A.a3) ^ 2 = 0 := by
      rcases mul_eq_zero.mp hz with h' | h'
      · exact absurd h' hv4.ne'
      · exact h'
    have hg2 : (0 : ℝ) < g ^ 2 := by positivity
    have key : ∀ x : ℝ, g ^ 2 * x ^ 2 ≤ 0 → x = 0 := by
      intro x hx
      have hxz : g ^ 2 * x ^ 2 = 0 :=
        le_antisymm hx (by positivity)
      have := (mul_eq_zero.mp hxz).resolve_left hg2.ne'
      exact (pow_eq_zero_iff two_ne_zero).mp this
    have ha1 : A.a1 = 0 := by
      refine key A.a1 ?_
      nlinarith [sq_nonneg A.a2, sq_nonneg (g' * A.b - g * A.a3),
        mul_nonneg hg2.le (sq_nonneg A.a2)]
    have ha2 : A.a2 = 0 := by
      refine key A.a2 ?_
      nlinarith [sq_nonneg A.a1, sq_nonneg (g' * A.b - g * A.a3),
        mul_nonneg hg2.le (sq_nonneg A.a1)]
    have hd : (g' * A.b - g * A.a3) ^ 2 = 0 := by
      nlinarith [mul_nonneg hg2.le (sq_nonneg A.a1),
        mul_nonneg hg2.le (sq_nonneg A.a2)]
    have := (pow_eq_zero_iff two_ne_zero).mp hd
    exact ⟨ha1, ha2, by linarith [sub_eq_zero.mp this]⟩
  · rintro ⟨h1, h2, h3⟩
    rw [infAction_vacuum_coords]
    have hd : g' * A.b - g * A.a3 = 0 := by rw [h3]; ring
    funext i
    fin_cases i <;> simp [h1, h2, hd]

/-- **The exact stabilizer (unbroken direction).**  With a positive
declared weak coupling and scale, the stabilizer of the chosen minimum is
exactly the one-parameter line mixing the third weak direction and
the hypercharge direction with weights `(g', g)`: one unbroken
abelian direction inside the committed electroweak label pair, with a
unique line parameter for every stabilizing direction. -/
theorem stabilizer_exactly_one_line (g g' vev : ℝ) (hg : 0 < g)
    (hvev : 0 < vev) (A : EWParam) :
    infAction g g' A (vacuum vev) = 0
      ↔ ∃! t : ℝ, A = ⟨0, 0, g' * t, g * t⟩ := by
  rw [infAction_vacuum_zero_iff g g' vev hg hvev]
  constructor
  · rintro ⟨h1, h2, h3⟩
    refine ⟨A.b / g, ?_, ?_⟩
    · ext
      · exact h1
      · exact h2
      · field_simp
        linarith [h3]
      · field_simp
    · intro t ht
      have hb : A.b = g * t := congrArg EWParam.b ht
      rw [hb, mul_comm g t, mul_div_assoc, div_self hg.ne', mul_one]
  · rintro ⟨t, rfl, -⟩
    refine ⟨rfl, rfl, ?_⟩
    show g' * (g * t) = g * (g' * t)
    ring

/-! ## The quadratic coefficient at the minimum and the Goldstone count -/

/-- The exact coefficient of `t²` in the potential along a ray from the
chosen minimum.  The historical name `secondVariation` uses the quadratic-
coefficient convention; the conventional second derivative at `t = 0` is
twice this value. -/
def secondVariation (lam vev : ℝ) (u : Doublet) : ℝ :=
  4 * lam * (vev * (u 1).re) ^ 2

/-- Exact squared-length expansion along a ray from the minimum. -/
theorem normSqSum_vacuum_add_smul (vev t : ℝ) (u : Doublet) :
    normSqSum (vacuum vev + t • u)
      = vev ^ 2 + 2 * t * vev * (u 1).re + t ^ 2 * normSqSum u := by
  simp only [normSqSum, vacuum, Pi.add_apply, Pi.smul_apply,
    Complex.real_smul, Matrix.cons_val_zero, Matrix.cons_val_one,
    Complex.normSq_apply, Complex.add_re,
    Complex.add_im, Complex.mul_re, Complex.mul_im, Complex.ofReal_re,
    Complex.ofReal_im, zero_add]
  ring

/-- **The exact ray expansion of the potential at the minimum.**  The
potential along any ray from the chosen minimum is an exact
polynomial in the ray parameter with no constant and no linear term;
the quadratic coefficient is `secondVariation` in the file's stated
coefficient convention.  This is the
polynomial form of criticality plus second-order behavior, with no
limit taken. -/
theorem potential_ray_expansion (lam vev t : ℝ) (u : Doublet) :
    potential lam vev (vacuum vev + t • u)
      = secondVariation lam vev u * t ^ 2
        + 4 * lam * vev * (u 1).re * normSqSum u * t ^ 3
        + lam * normSqSum u ^ 2 * t ^ 4 := by
  unfold potential secondVariation
  rw [normSqSum_vacuum_add_smul]
  ring

/-- The quadratic coefficient is positive semidefinite. -/
theorem secondVariation_nonneg (lam vev : ℝ) (hlam : 0 < lam)
    (u : Doublet) : 0 ≤ secondVariation lam vev u := by
  unfold secondVariation
  positivity

/-- The kernel equation of the quadratic coefficient: with positive
declared `lam` and `vev`, it vanishes exactly on
the directions with real-free lower component. -/
theorem secondVariation_zero_iff (lam vev : ℝ) (hlam : 0 < lam)
    (hvev : 0 < vev) (u : Doublet) :
    secondVariation lam vev u = 0 ↔ (u 1).re = 0 := by
  unfold secondVariation
  constructor
  · intro h
    have h4 : (0 : ℝ) < 4 * lam := by positivity
    have hsq : (vev * (u 1).re) ^ 2 = 0 :=
      (mul_eq_zero.mp h).resolve_left h4.ne'
    have := (pow_eq_zero_iff two_ne_zero).mp hsq
    exact (mul_eq_zero.mp this).resolve_left hvev.ne'
  · intro h
    rw [h]
    ring

/-- **The broken directions are the gauge orbit directions.**  With
positive declared couplings and scale, a doublet direction is an
infinitesimal gauge image at the minimum exactly when its lower
component has no real part, which is exactly the kernel of the second
variation. -/
theorem broken_directions_exact (g g' vev : ℝ) (hg : 0 < g)
    (hg' : 0 < g') (hvev : 0 < vev) (u : Doublet) :
    (∃ A : EWParam, infAction g g' A (vacuum vev) = u)
      ↔ (u 1).re = 0 := by
  constructor
  · rintro ⟨A, rfl⟩
    rw [infAction_vacuum_coords]
    simp
  · intro hre
    have e1 : g * (2 * (u 0).re / (g * vev)) * vev / 2 = (u 0).re := by
      field_simp [hg.ne', hvev.ne']
    have e2 : g * (2 * (u 0).im / (g * vev)) * vev / 2 = (u 0).im := by
      field_simp [hg.ne', hvev.ne']
    have e3 : (g' * (2 * (u 1).im / (g' * vev)) - g * 0) * vev / 2
        = (u 1).im := by
      field_simp [hg'.ne', hvev.ne']
      ring
    refine ⟨⟨2 * (u 0).im / (g * vev), 2 * (u 0).re / (g * vev), 0,
      2 * (u 1).im / (g' * vev)⟩, ?_⟩
    rw [infAction_vacuum_coords]
    funext i
    fin_cases i
    · show ((g * (2 * (u 0).re / (g * vev)) * vev / 2 : ℝ) : ℂ)
          + ((g * (2 * (u 0).im / (g * vev)) * vev / 2 : ℝ) : ℂ)
            * Complex.I
        = u 0
      rw [e1, e2, Complex.re_add_im]
    · show ((0 : ℝ) : ℂ)
          + (((g' * (2 * (u 1).im / (g' * vev)) - g * 0) * vev / 2 : ℝ)
              : ℂ) * Complex.I
        = u 1
      rw [e3]
      apply Complex.ext
      · simpa using hre.symm
      · simp

/-- The three-parameter normal form of a broken direction. -/
noncomputable def goldstoneForm (c : ℝ × ℝ × ℝ) : Doublet :=
  ![(c.1 : ℂ) + (c.2.1 : ℂ) * Complex.I, (c.2.2 : ℂ) * Complex.I]

/-- **The exact Goldstone count.**  A doublet direction lies in the
kernel of the second variation exactly when it has the unique
three-real-parameter normal form: the broken direction space is
exactly three-parameter. -/
theorem goldstone_exactly_three_parameters (u : Doublet) :
    (u 1).re = 0 ↔ ∃! c : ℝ × ℝ × ℝ, u = goldstoneForm c := by
  constructor
  · intro hre
    refine ⟨((u 0).re, (u 0).im, (u 1).im), ?_, ?_⟩
    · funext i
      fin_cases i
      · simp only [goldstoneForm]
        exact (Complex.re_add_im (u 0)).symm
      · simp only [goldstoneForm]
        apply Complex.ext <;> simp [hre]
    · intro c hc
      have h0 : u 0 = (c.1 : ℂ) + (c.2.1 : ℂ) * Complex.I := by
        rw [hc]; simp [goldstoneForm]
      have h1 : u 1 = (c.2.2 : ℂ) * Complex.I := by
        rw [hc]; simp [goldstoneForm]
      refine Prod.ext ?_ (Prod.ext ?_ ?_) <;> simp [h0, h1]
  · rintro ⟨c, rfl, -⟩
    simp [goldstoneForm]

/-! ## The tree-level gauge mass values and their ratio in the
declared constants -/

/-- The charged (W-shaped) direction. -/
def wDirection : EWParam := ⟨1, 0, 0, 0⟩

/-- The neutral broken (Z-shaped) direction, unnormalized. -/
def zDirection (g g' : ℝ) : EWParam := ⟨0, 0, g, -g'⟩

/-- The unbroken (photon-shaped) direction, unnormalized. -/
def photonDirection (g g' : ℝ) : EWParam := ⟨0, 0, g', g⟩

/-- The squared Euclidean length of a direction parameter, used to
normalize the unnormalized neutral directions. -/
def paramNormSq (A : EWParam) : ℝ :=
  A.a1 ^ 2 + A.a2 ^ 2 + A.a3 ^ 2 + A.b ^ 2

theorem wDirection_paramNormSq : paramNormSq wDirection = 1 := by
  simp only [paramNormSq, wDirection]
  norm_num

theorem zDirection_paramNormSq (g g' : ℝ) :
    paramNormSq (zDirection g g') = g ^ 2 + g' ^ 2 := by
  simp only [paramNormSq, zDirection]
  ring

/-- The neutral broken and unbroken directions are orthogonal in the
parameter space. -/
theorem z_photon_orthogonal (g g' : ℝ) :
    (zDirection g g').a3 * (photonDirection g g').a3
      + (zDirection g g').b * (photonDirection g g').b = 0 := by
  simp only [zDirection, photonDirection]
  ring

/-- The charged-direction mass value `g^2 vev^2 / 4`. -/
theorem wMass_value (g g' vev : ℝ) :
    gaugeMassForm g g' vev wDirection = g ^ 2 * vev ^ 2 / 4 := by
  rw [gaugeMassForm_exact]
  simp only [wDirection]
  ring

/-- The unnormalized neutral broken-direction mass value. -/
theorem zMass_value (g g' vev : ℝ) :
    gaugeMassForm g g' vev (zDirection g g')
      = (g ^ 2 + g' ^ 2) ^ 2 * vev ^ 2 / 4 := by
  rw [gaugeMassForm_exact]
  simp only [zDirection]
  ring

/-- The normalized neutral broken-direction mass value
`(g^2 + g'^2) vev^2 / 4`. -/
theorem zMass_normalized (g g' vev : ℝ) (hg : 0 < g) (hg' : 0 < g') :
    gaugeMassForm g g' vev (zDirection g g')
        / paramNormSq (zDirection g g')
      = (g ^ 2 + g' ^ 2) * vev ^ 2 / 4 := by
  rw [zMass_value, zDirection_paramNormSq]
  have hgg : (0 : ℝ) < g ^ 2 + g' ^ 2 := by positivity
  field_simp

/-- The unbroken direction is exactly massless. -/
theorem photon_massless (g g' vev : ℝ) :
    gaugeMassForm g g' vev (photonDirection g g') = 0 := by
  rw [gaugeMassForm_exact]
  simp only [photonDirection]
  ring

/-- **The tree mass ratio in the declared constants.**  The ratio of
the normalized charged and neutral broken mass values is
`g^2 / (g^2 + g'^2)`: the cosine-squared shape.  The committed
grammar does not select `g'/g` (the committed invariant-form receipts
leave the relative coupling normalization free), so both constants
are declared bundle parameters and this ratio is a declared-parameter
identity, not a derived weak mixing angle. -/
theorem tree_mass_ratio (g g' vev : ℝ) (hg : 0 < g) (hg' : 0 < g')
    (hvev : 0 < vev) :
    (gaugeMassForm g g' vev wDirection / paramNormSq wDirection)
        / (gaugeMassForm g g' vev (zDirection g g')
            / paramNormSq (zDirection g g'))
      = g ^ 2 / (g ^ 2 + g' ^ 2) := by
  rw [wMass_value, wDirection_paramNormSq,
    zMass_normalized g g' vev hg hg']
  have hgg : (0 : ℝ) < g ^ 2 + g' ^ 2 := by positivity
  have hv2 : (0 : ℝ) < vev ^ 2 := by positivity
  field_simp

/-! ## One Yukawa line on the committed grammar -/

/-- The weak epsilon contraction of two doublets. -/
def epsPair (a b : Doublet) : ℂ := a 0 * b 1 - a 1 * b 0

/-- One Yukawa line with one declared coefficient: the epsilon-paired
coupling of the scalar to a doublet carrier and a singlet (register
row PR-69: one finite line; physical matrices and family structure remain PR-50). -/
noncomputable def yukawaLine (y : ℝ) (h ψ : Doublet) (χ : ℂ) : ℂ :=
  (y : ℂ) * epsPair ψ h * χ

/-- **The fermion mass coefficient at the minimum.**  At the chosen
minimum the Yukawa line is an exact bilinear in the doublet and
singlet carriers with coefficient `y * vev`. -/
theorem yukawa_mass_at_vacuum (y vev : ℝ) (ψ : Doublet) (χ : ℂ) :
    yukawaLine y (vacuum vev) ψ χ
      = ((y * vev : ℝ) : ℂ) * (ψ 0 * χ) := by
  unfold yukawaLine epsPair vacuum
  simp
  ring

/-- The mass bilinear vanishes identically exactly when the declared
Yukawa coefficient vanishes or the declared scale vanishes. -/
theorem yukawa_line_zero_iff (y vev : ℝ) :
    (∀ (ψ : Doublet) (χ : ℂ), yukawaLine y (vacuum vev) ψ χ = 0)
      ↔ y = 0 ∨ vev = 0 := by
  constructor
  · intro h
    have h1 := h ![1, 0] 1
    rw [yukawa_mass_at_vacuum] at h1
    simpa using h1
  · intro h ψ χ
    rw [yukawa_mass_at_vacuum]
    rcases h with h | h <;> simp [h]

/-! ## The typed premise bundle (hard rule: one antecedent bundle) -/

/-- The composed electroweak-breaking premise bundle: the committed
Standard-Model structure premises (register rows PR-09, PR-10, PR-11,
PR-59) extended by the declared scalar-sector, coupling, and Yukawa
parameters.  Every added field is a declared finite hypothesis under
PR-67--PR-69.  The physical scalar/action, vacuum selection, complete
Yukawa matrices, and source action remain open under PR-48--PR-50 and
PR-54.  None is derived here. -/
structure EWBreakingPremiseData where
  /-- The committed Standard-Model structure premises. -/
  base : SMStructurePremiseData
  /-- The declared finite quartic coupling (register row PR-68). -/
  lam : ℝ
  /-- The quartic coupling is positive: the declared stability clause
  of register row PR-68. -/
  lam_pos : 0 < lam
  /-- The declared model-shell radius (register row PR-68). -/
  vev : ℝ
  /-- The breaking scale is positive: the declared broken-branch
  clause of register row PR-68. -/
  vev_pos : 0 < vev
  /-- The declared finite weak coupling constant (register row PR-67). -/
  gW : ℝ
  /-- Positivity of the declared weak coupling. -/
  gW_pos : 0 < gW
  /-- The declared finite hypercharge coupling constant (register row PR-67). -/
  gY : ℝ
  /-- Positivity of the declared hypercharge coupling. -/
  gY_pos : 0 < gY
  /-- The declared coefficient of the single finite Yukawa line (PR-69). -/
  yuk : ℝ

/-- The bundled selection mask supplies, in either committed parity
sector, one colorless doublet row and one colorless singlet row of
matching parity grading and matching Weyl column (through the
committed index dictionary), whose committed charges balance against
the scalar row charge or its conjugate: the charge bookkeeping
required of one gauge-invariant Yukawa line on the committed grammar
(the abelian charge-sum balance; an invariance statement for the line
itself is not proved here). -/
theorem selected_sector_supports_yukawa_line
    (D : EWBreakingPremiseData) :
    ∃ dRow sRow : Fin 10,
      OPH.ExteriorSelection.mem D.base.selectionMask.val dRow = true
        ∧ OPH.ExteriorSelection.mem D.base.selectionMask.val sRow = true
        ∧ OPH.ExteriorSelection.isDoublet dRow = true
        ∧ OPH.ExteriorSelection.isDoublet sRow = false
        ∧ OPH.ExteriorSelection.isColored dRow = false
        ∧ OPH.ExteriorSelection.isColored sRow = false
        ∧ OPH.ExteriorSelection.odd dRow = OPH.ExteriorSelection.odd sRow
        ∧ OPH.BaryonDimensionSix.left
              (OPH.MatterGrammarIndexBridge.fieldIndex dRow)
            = OPH.BaryonDimensionSix.left
              (OPH.MatterGrammarIndexBridge.fieldIndex sRow)
        ∧ (OPH.ExteriorSelection.charge dRow
              + OPH.ExteriorSelection.charge sRow + scalarChargeSix = 0
            ∨ OPH.ExteriorSelection.charge dRow
              + OPH.ExteriorSelection.charge sRow - scalarChargeSix = 0) := by
  rcases D.base.matter_selection_is_parity_sector with h | h
  · exact ⟨9, 4, by rw [h]; decide, by rw [h]; decide, by decide,
      by decide, by decide, by decide, by decide, by decide,
      Or.inr (by decide)⟩
  · exact ⟨1, 7, by rw [h]; decide, by rw [h]; decide, by decide,
      by decide, by decide, by decide, by decide, by decide,
      Or.inl (by decide)⟩

/-! ## The composed receipt -/

/-- **The composed electroweak-breaking/Yukawa receipt (issue #735,
OL-G7 partial).**  For every premise bundle `D`: the minimizer set of
the declared quartic potential is exactly the nonzero shell
`|h|^2 = vev^2` and contains the chosen minimum; the stabilizer of
the chosen minimum is exactly one one-parameter unbroken line inside
the four-parameter electroweak direction space; the second variation
at the minimum is positive semidefinite and is the exact quadratic
coefficient of the polynomial ray expansion of the declared potential
(no limit taken), its kernel is exactly the set of infinitesimal
gauge images at the minimum, and that kernel is exactly
three-parameter (the Goldstone count); the gauge mass form
vanishes exactly on the stabilizer, takes the values
`gW^2 vev^2 / 4` and `(gW^2 + gY^2) vev^2 / 4` on the normalized
charged and neutral broken directions, keeps the unbroken direction
exactly massless, and the tree ratio of the broken values is
`gW^2 / (gW^2 + gY^2)` in the declared constants; the Yukawa line at
the minimum is an exact bilinear with coefficient `yuk * vev`,
vanishing identically exactly when `yuk = 0` or `vev = 0`; and the
bundled selection mask supports the charge bookkeeping required of
one such line on the committed grammar.  Every clause consumes the
bundle.
The theorem consumes only the finite declarations PR-67--PR-69.  It
does not discharge PR-47, PR-48, PR-49, PR-50, or PR-54: no physical
scalar, kinetic normalization, vacuum selection, complete Yukawa
matrices, spacetime action, measured mass, or family structure is claimed. -/
theorem electroweakBreakingComposition_receipt
    (D : EWBreakingPremiseData) :
    ((∀ h : Doublet,
        (∀ h' : Doublet,
          potential D.lam D.vev h ≤ potential D.lam D.vev h')
          ↔ normSqSum h = D.vev ^ 2)
      ∧ normSqSum (vacuum D.vev) = D.vev ^ 2
      ∧ (∀ h : Doublet, normSqSum h = D.vev ^ 2 → h ≠ 0))
    ∧ (∀ A : EWParam,
        infAction D.gW D.gY A (vacuum D.vev) = 0
          ↔ ∃! t : ℝ, A = ⟨0, 0, D.gY * t, D.gW * t⟩)
    ∧ ((∀ u : Doublet, 0 ≤ secondVariation D.lam D.vev u)
      ∧ (∀ (t : ℝ) (u : Doublet),
          potential D.lam D.vev (vacuum D.vev + t • u)
            = secondVariation D.lam D.vev u * t ^ 2
              + 4 * D.lam * D.vev * (u 1).re * normSqSum u * t ^ 3
              + D.lam * normSqSum u ^ 2 * t ^ 4)
      ∧ (∀ u : Doublet,
          secondVariation D.lam D.vev u = 0
            ↔ ∃ A : EWParam, infAction D.gW D.gY A (vacuum D.vev) = u)
      ∧ (∀ u : Doublet,
          secondVariation D.lam D.vev u = 0
            ↔ ∃! c : ℝ × ℝ × ℝ, u = goldstoneForm c))
    ∧ ((∀ A : EWParam,
          gaugeMassForm D.gW D.gY D.vev A = 0
            ↔ infAction D.gW D.gY A (vacuum D.vev) = 0)
      ∧ gaugeMassForm D.gW D.gY D.vev wDirection
          = D.gW ^ 2 * D.vev ^ 2 / 4
      ∧ gaugeMassForm D.gW D.gY D.vev (zDirection D.gW D.gY)
            / paramNormSq (zDirection D.gW D.gY)
          = (D.gW ^ 2 + D.gY ^ 2) * D.vev ^ 2 / 4
      ∧ gaugeMassForm D.gW D.gY D.vev (photonDirection D.gW D.gY) = 0
      ∧ (gaugeMassForm D.gW D.gY D.vev wDirection
            / paramNormSq wDirection)
          / (gaugeMassForm D.gW D.gY D.vev (zDirection D.gW D.gY)
              / paramNormSq (zDirection D.gW D.gY))
          = D.gW ^ 2 / (D.gW ^ 2 + D.gY ^ 2))
    ∧ ((∀ (ψ : Doublet) (χ : ℂ),
          yukawaLine D.yuk (vacuum D.vev) ψ χ
            = ((D.yuk * D.vev : ℝ) : ℂ) * (ψ 0 * χ))
      ∧ ((∀ (ψ : Doublet) (χ : ℂ),
            yukawaLine D.yuk (vacuum D.vev) ψ χ = 0)
          ↔ D.yuk = 0 ∨ D.vev = 0)
      ∧ ∃ dRow sRow : Fin 10,
          OPH.ExteriorSelection.mem D.base.selectionMask.val dRow = true
            ∧ OPH.ExteriorSelection.mem D.base.selectionMask.val sRow
                = true
            ∧ OPH.ExteriorSelection.isDoublet dRow = true
            ∧ OPH.ExteriorSelection.isDoublet sRow = false
            ∧ OPH.ExteriorSelection.isColored dRow = false
            ∧ OPH.ExteriorSelection.isColored sRow = false
            ∧ OPH.ExteriorSelection.odd dRow
                = OPH.ExteriorSelection.odd sRow
            ∧ OPH.BaryonDimensionSix.left
                  (OPH.MatterGrammarIndexBridge.fieldIndex dRow)
                = OPH.BaryonDimensionSix.left
                  (OPH.MatterGrammarIndexBridge.fieldIndex sRow)
            ∧ (OPH.ExteriorSelection.charge dRow
                  + OPH.ExteriorSelection.charge sRow
                  + scalarChargeSix = 0
                ∨ OPH.ExteriorSelection.charge dRow
                  + OPH.ExteriorSelection.charge sRow
                  - scalarChargeSix = 0)) := by
  refine ⟨⟨fun h => minimizer_iff D.lam D.vev D.lam_pos h,
      vacuum_normSqSum D.vev,
      fun h => minimizer_nonzero D.vev D.vev_pos h⟩,
    fun A => stabilizer_exactly_one_line D.gW D.gY D.vev D.gW_pos
      D.vev_pos A,
    ⟨fun u => secondVariation_nonneg D.lam D.vev D.lam_pos u,
      fun t u => potential_ray_expansion D.lam D.vev t u,
      fun u => ?_, fun u => ?_⟩,
    ⟨fun A => gaugeMassForm_zero_iff D.gW D.gY D.vev A,
      wMass_value D.gW D.gY D.vev,
      zMass_normalized D.gW D.gY D.vev D.gW_pos D.gY_pos,
      photon_massless D.gW D.gY D.vev,
      tree_mass_ratio D.gW D.gY D.vev D.gW_pos D.gY_pos D.vev_pos⟩,
    ⟨fun ψ χ => yukawa_mass_at_vacuum D.yuk D.vev ψ χ,
      yukawa_line_zero_iff D.yuk D.vev,
      selected_sector_supports_yukawa_line D⟩⟩
  · rw [secondVariation_zero_iff D.lam D.vev D.lam_pos D.vev_pos u]
    exact (broken_directions_exact D.gW D.gY D.vev D.gW_pos D.gY_pos
      D.vev_pos u).symm
  · rw [secondVariation_zero_iff D.lam D.vev D.lam_pos D.vev_pos u]
    exact goldstone_exactly_three_parameters u

/-! ## Nonvacuity -/

/-- The committed inhabitant: the committed even-sector base instance
with every declared parameter set to one.  The instance is a
nonvacuity control; the parameter values select nothing. -/
noncomputable def committedEWBreakingPremiseData :
    EWBreakingPremiseData where
  base := committedSMStructurePremiseData
  lam := 1
  lam_pos := one_pos
  vev := 1
  vev_pos := one_pos
  gW := 1
  gW_pos := one_pos
  gY := 1
  gY_pos := one_pos
  yuk := 1

/-- The premise bundle is satisfiable. -/
theorem ewBreakingPremiseData_nonvacuous :
    Nonempty EWBreakingPremiseData :=
  ⟨committedEWBreakingPremiseData⟩

/-- Exact tree values at the committed inhabitant (all declared
parameters one): charged value `1/4`, normalized neutral broken value
`1/2`, unbroken value `0`, and Yukawa mass coefficient `1`.  These
are witness values of the declared parameters, not predictions. -/
theorem committed_instance_values :
    gaugeMassForm 1 1 1 wDirection = 1 / 4
      ∧ gaugeMassForm 1 1 1 (zDirection 1 1)
          / paramNormSq (zDirection 1 1) = 1 / 2
      ∧ gaugeMassForm 1 1 1 (photonDirection 1 1) = 0
      ∧ ((1 : ℝ) * 1 : ℝ) = 1 := by
  refine ⟨?_, ?_, photon_massless 1 1 1, by norm_num⟩
  · rw [wMass_value]
    norm_num
  · rw [zMass_normalized 1 1 1 one_pos one_pos]
    norm_num

end OPH.ElectroweakBreakingComposition

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.ElectroweakBreakingComposition.scalar_row_receipt
#print axioms OPH.ElectroweakBreakingComposition.minimizer_iff
#print axioms OPH.ElectroweakBreakingComposition.minimizer_nonzero
#print axioms OPH.ElectroweakBreakingComposition.stabilizer_exactly_one_line
#print axioms OPH.ElectroweakBreakingComposition.potential_ray_expansion
#print axioms OPH.ElectroweakBreakingComposition.broken_directions_exact
#print axioms OPH.ElectroweakBreakingComposition.goldstone_exactly_three_parameters
#print axioms OPH.ElectroweakBreakingComposition.tree_mass_ratio
#print axioms OPH.ElectroweakBreakingComposition.yukawa_mass_at_vacuum
#print axioms OPH.ElectroweakBreakingComposition.yukawa_line_zero_iff
#print axioms OPH.ElectroweakBreakingComposition.selected_sector_supports_yukawa_line
#print axioms OPH.ElectroweakBreakingComposition.electroweakBreakingComposition_receipt
#print axioms OPH.ElectroweakBreakingComposition.ewBreakingPremiseData_nonvacuous
#print axioms OPH.ElectroweakBreakingComposition.committed_instance_values
