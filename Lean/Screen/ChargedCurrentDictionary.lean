import NeutralCurrentDictionary

/-!
# The charged-current coupling dictionary on the committed grammar
(V3, issue #745)

The charged-current companion of the committed neutral-current
dictionary: the exact raising/lowering pair built from the committed
broken directions, the per-row charged coupling dictionary on the
weak slots of the committed exterior component table, the exact
selection rule, and the charge bookkeeping, all in the declared
constants of the committed breaking bundle.

The committed inputs.  `Screen/ExteriorSelection.lean` (register row
PR-59) carries the ten-row exterior component table with integer
charges in the `q = 6Y` normalization and the weak-doublet indicator
column; `Screen/ElectroweakBreakingComposition.lean` (issue #735)
carries the four-parameter electroweak direction space `EWParam`, the
declared coupling constants `gW`, `gY` of its premise bundle
(register rows PR-48/PR-54: declared, never selected), the committed
infinitesimal action on the scalar carrier, the committed unbroken
direction `photonDirection`, the committed neutral broken direction
`zDirection`, and the gauge mass form at the chosen minimum;
`Screen/NeutralCurrentDictionary.lean` (the committed first #745
receipt) carries the declared diagonal pairing `diagPairing`, the
third-isospin column `t3Six`, the electric-charge column `q6EM`, and
the neutral normalization data `gZ`, `eCoupling`.

The charged pair.  The two charged directions are `w1Direction`
(the committed `wDirection` of the breaking module) and
`w2Direction`, the two coordinate directions of the charged plane.
The proved clauses:

* **The charged plane.**  Both charged directions are orthogonal to
  the committed photon and Z directions and to each other under the
  declared diagonal pairing, and at positive declared constants a
  direction is orthogonal to both committed neutral directions
  exactly when it is a unique combination of the two charged
  directions (`charged_plane_exact`): the charged plane is exactly
  the orthogonal complement of the committed neutral plane.  Both
  charged directions carry the same broken gauge mass value
  `gW^2 vev^2 / 4`.
* **The raising/lowering pair.**  `raiseShift` and `lowerShift` are
  the two complex combinations `w1 ± i w2` of the committed
  infinitesimal actions.  Their exact closed forms are the strict
  off-diagonal shift operators with coefficient `gW`
  (`raiseShift_exact`, `lowerShift_exact`): the raising shift moves
  the lower weak slot to the upper one and annihilates the rest, the
  lowering shift is its transpose.
* **The derived coupling constant.**  The direction length of the
  charged pair under the committed pairing is
  `chargedNorm = sqrt 2` (`chargedNorm_sq`, from
  `paramNormSq w1 + paramNormSq w2 = 2` with the two directions
  orthogonal), so the coupling of the normalized charged combination
  is exactly `gW / sqrt 2`: the constant is derived from the
  committed pairing, never inserted.  The alignment receipts
  (`committed_action_alignment_raise`,
  `committed_action_alignment_lower`) pin the dictionary value to the
  committed infinitesimal action on the scalar row, and the rigidity
  receipt (`coefficient_rigidity`) shows any drifted coefficient
  breaks the alignment identity.
* **The per-row dictionary and the exact selection rule.**  The
  charged coupling of an ordered weak slot pair is `gW / chargedNorm`
  exactly on the committed partner pairs (the two weak components of
  a committed doublet row, upper receiving from lower) and zero
  otherwise; at a positive declared weak coupling the nonvanishing is
  an exact equivalence (`charged_selection_rule`): a nonzero coupling
  exists exactly for the committed partner pairs, and committed
  singlets are annihilated (`charged_annihilates_singlets`).
* **Charge bookkeeping.**  Every nonzero charged coupling changes the
  committed electric-charge entry by exactly `6`, one unit in the
  committed `q = 6Y` normalization, proved from the committed charge
  columns per row (`charged_coupling_charge_step`); the committed
  hypercharge column is constant on each row, so the step is pure
  third isospin.  The operator form: the commutator of the committed
  photon-direction action with the raising shift is exactly
  `i (gZ * e)` times the raising shift, and with the lowering shift
  exactly the negative (`raise_carries_unit_charge`,
  `lowerShift_photon_commutator`), so per unit of the committed
  photon normalization `gZ` the raising direction carries exactly one
  declared electromagnetic unit `e`.
* **Consistency joins with the committed bracket data.**  The
  commutator of the two charged-direction actions lands exactly on
  the committed neutral directions with the exact coefficients
  `-(gW gY)/(gW^2 + gY^2)` and `-(gW^2)/(gW^2 + gY^2)` of the
  committed photon and Z directions
  (`charged_commutator_neutral_span`), and the raising/lowering
  commutator is exactly `2 gW i` times the committed third-weak
  action (`ladder_commutator_neutral`).  Grade: the committed
  electroweak receipts carry no abstract Lie bracket on `EWParam`;
  the joins are proved as exact operator-commutator identities of the
  committed infinitesimal action on the scalar-row carrier, which is
  the exact level the committed bracket receipts support.
* **The composed receipt.**  `chargedCurrentDictionary_receipt`
  states all clauses from the one committed premise bundle
  `EWBreakingPremiseData`, with the dictionary and bookkeeping
  clauses quantified over the rows of the bundled selection mask and
  a mutation witness in either committed parity sector: a selected
  row with singlet flag gets zero coupling at every slot pair, and
  every reversed or repeated slot pair gets zero coupling on every
  row.  Nonvacuity is inherited from the committed
  all-parameters-one instance.

Consumed register rows.  The bundle carries the committed
Standard-Model structure premises (register rows PR-09, PR-10, PR-11,
PR-59) through its base, and the declared constants `gW`, `gY` whose
selection stays on register rows PR-48 and PR-54; the potential
parameters of row PR-49 enter only through the committed vacuum used
in the gauge-mass receipts.  The charge and doublet columns are the
committed PR-59 table entries.

Boundary and nonclaims.  Every constant is declared: no measured
`gW`, no Fermi constant, no measured width.  The dictionary is exact
tree-level algebra on the committed finite grammar of register row
PR-59: there is no physical current operator, no spacetime, no
covariant derivative, no source-produced action (register row PR-54
open), no W propagator, no decay amplitude, no Fermi limit.  CKM and
flavor structure are absent: the dictionary is one-generation
bookkeeping on the committed table, with no mixing matrix, no phase,
and no family structure claimed (register row PR-50 untouched,
register row PR-36 three-family structure not consumed).  The weak
components of a doublet row are the committed weak-dimension
bookkeeping of the table; a physical chirality attachment stays on
register row PR-47.  This is a structural receipt on issue #745,
companion to the committed neutral dictionary, and it closes
nothing: rows OL-M1, OL-M2, OL-M3 keep their owed physical content.

Falsifier.  The dictionary fails if either charged direction is not
orthogonal to both committed neutral directions, if the charged
plane is not exactly the orthogonal complement of the neutral plane,
if the raising or lowering closed form differs from the strict shift
with coefficient `gW`, if the derived normalization differs from
`sqrt 2`, if a nonzero charged coupling exists off the committed
partner pairs or vanishes on one at positive `gW`, if a nonzero
coupling fails to step the committed electric charge by exactly `6`,
if the photon commutator eigenvalue differs from `i gZ e` up to
sign, or if either charged commutator lands off the committed
neutral span or the ladder commutator off the committed third-weak
action with the stated exact coefficients.
-/

namespace OPH.ChargedCurrentDictionary

open OPH.ElectroweakBreakingComposition (Doublet vacuum infAction
  EWParam paramNormSq photonDirection zDirection wDirection
  gaugeMassForm gaugeMassForm_exact EWBreakingPremiseData
  committedEWBreakingPremiseData)

open OPH.NeutralCurrentDictionary (t3Six q6EM diagPairing gZ
  eCoupling gZ_pos gZ_sq)

/-! ## The two charged directions and the charged plane -/

/-- The first charged direction: the committed `wDirection` of the
breaking module. -/
def w1Direction : EWParam := ⟨1, 0, 0, 0⟩

/-- The second charged direction: the second coordinate direction of
the charged plane. -/
def w2Direction : EWParam := ⟨0, 1, 0, 0⟩

/-- The first charged direction is the committed broken direction of
the breaking module. -/
theorem w1_is_committed_wDirection : w1Direction = wDirection := rfl

/-- Both charged directions are orthogonal to both committed neutral
directions and to each other under the declared diagonal pairing. -/
theorem charged_neutral_orthogonal (gW gY : ℝ) :
    diagPairing w1Direction (photonDirection gW gY) = 0
      ∧ diagPairing w1Direction (zDirection gW gY) = 0
      ∧ diagPairing w2Direction (photonDirection gW gY) = 0
      ∧ diagPairing w2Direction (zDirection gW gY) = 0
      ∧ diagPairing w1Direction w2Direction = 0 := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;>
    · simp only [diagPairing, w1Direction, w2Direction,
        photonDirection, zDirection]
      ring

/-- **The charged plane is exactly the orthogonal complement of the
committed neutral plane.**  At positive declared constants, a
direction is orthogonal to both the committed photon direction and
the committed Z direction under the declared diagonal pairing exactly
when it is a unique combination of the two charged directions. -/
theorem charged_plane_exact (gW gY : ℝ) (hgW : 0 < gW) (hgY : 0 < gY)
    (A : EWParam) :
    (diagPairing A (photonDirection gW gY) = 0
        ∧ diagPairing A (zDirection gW gY) = 0)
      ↔ ∃! p : ℝ × ℝ, A = ⟨p.1, p.2, 0, 0⟩ := by
  constructor
  · rintro ⟨hp, hz⟩
    simp only [diagPairing, photonDirection, zDirection] at hp hz
    have hgg : gW ^ 2 + gY ^ 2 ≠ 0 := by positivity
    have hkey : (gW ^ 2 + gY ^ 2) * A.a3 = 0 := by
      linear_combination gW * hz + gY * hp
    have ha3 : A.a3 = 0 := (mul_eq_zero.mp hkey).resolve_left hgg
    have hb : A.b = 0 := by
      rw [ha3] at hp
      have hb' : A.b * gW = 0 := by linear_combination hp
      exact (mul_eq_zero.mp hb').resolve_right hgW.ne'
    refine ⟨(A.a1, A.a2), ?_, ?_⟩
    · ext
      · rfl
      · rfl
      · exact ha3
      · exact hb
    · rintro ⟨x, y⟩ hxy
      have h1 : A.a1 = x := by rw [hxy]
      have h2 : A.a2 = y := by rw [hxy]
      rw [h1, h2]
  · rintro ⟨p, rfl, -⟩
    constructor <;>
      · simp only [diagPairing, photonDirection, zDirection]
        ring

/-- Both charged directions carry the same broken gauge mass value in
the declared constants. -/
theorem charged_directions_mass (gW gY vev : ℝ) :
    gaugeMassForm gW gY vev w1Direction = gW ^ 2 * vev ^ 2 / 4
      ∧ gaugeMassForm gW gY vev w2Direction
          = gW ^ 2 * vev ^ 2 / 4 := by
  constructor <;>
    · rw [gaugeMassForm_exact]
      simp only [w1Direction, w2Direction]
      ring

/-! ## The derived normalization of the charged pair -/

/-- The direction length of the charged pair under the committed
pairing: derived from the committed squared lengths of the two
charged directions, never inserted. -/
noncomputable def chargedNorm : ℝ :=
  Real.sqrt (paramNormSq w1Direction + paramNormSq w2Direction)

/-- The committed squared lengths of the charged pair total two. -/
theorem charged_pair_normSq :
    paramNormSq w1Direction + paramNormSq w2Direction = 2 := by
  simp only [paramNormSq, w1Direction, w2Direction]
  norm_num

/-- The squared derived normalization is exactly two. -/
theorem chargedNorm_sq : chargedNorm ^ 2 = 2 := by
  unfold chargedNorm
  rw [charged_pair_normSq]
  exact Real.sq_sqrt (by norm_num)

theorem chargedNorm_pos : 0 < chargedNorm := by
  unfold chargedNorm
  rw [charged_pair_normSq]
  exact Real.sqrt_pos.mpr (by norm_num)

/-! ## The raising and lowering shifts from the committed action -/

/-- The raising shift: the `w1 + i w2` complex combination of the
committed infinitesimal actions of the two charged directions. -/
noncomputable def raiseShift (gW gY : ℝ) (h : Doublet) : Doublet :=
  fun j => infAction gW gY w1Direction h j
    + Complex.I * infAction gW gY w2Direction h j

/-- The lowering shift: the `w1 - i w2` complex combination of the
committed infinitesimal actions of the two charged directions. -/
noncomputable def lowerShift (gW gY : ℝ) (h : Doublet) : Doublet :=
  fun j => infAction gW gY w1Direction h j
    - Complex.I * infAction gW gY w2Direction h j

/-- **The exact raising closed form.**  The raising shift is the
strict off-diagonal shift with coefficient `gW`: it moves the lower
weak slot to the upper one and annihilates the rest. -/
theorem raiseShift_exact (gW gY : ℝ) (h : Doublet) :
    raiseShift gW gY h
      = ![Complex.I * ((gW : ℝ) : ℂ) * h 1, 0] := by
  funext j
  fin_cases j <;>
    · simp [raiseShift, infAction, w1Direction, w2Direction]
      try ring

/-- **The exact lowering closed form.**  The lowering shift is the
transpose strict shift with the same coefficient `gW`. -/
theorem lowerShift_exact (gW gY : ℝ) (h : Doublet) :
    lowerShift gW gY h
      = ![0, Complex.I * ((gW : ℝ) : ℂ) * h 0] := by
  funext j
  fin_cases j <;>
    · simp [lowerShift, infAction, w1Direction, w2Direction]
      try ring

/-! ## The per-row charged coupling dictionary -/

/-- The charged coupling of the ordered weak slot pair `(c, d)` of
row `i`: the coefficient with which the normalized raising direction
connects slot `d` to slot `c`.  The value `gW / chargedNorm` on the
committed partner pairs is pinned to the committed infinitesimal
action by the alignment receipts; the partner condition is the
committed third-isospin step of six sixths. -/
noncomputable def chargedCoupling (gW : ℝ) (i : Fin 10) (c d : Fin 2) :
    ℝ :=
  if t3Six i c - t3Six i d = 6 then gW / chargedNorm else 0

/-- The partner condition in committed columns: the third-isospin
step is six sixths exactly on the committed doublet partner pairs,
upper component receiving from lower. -/
theorem partner_condition :
    ∀ (i : Fin 10) (c d : Fin 2),
      t3Six i c - t3Six i d = 6
        ↔ (OPH.ExteriorSelection.isDoublet i = true
            ∧ c = 0 ∧ d = 1) := by
  decide

/-- The dictionary value on a committed partner pair. -/
theorem chargedCoupling_partner (gW : ℝ) (i : Fin 10)
    (hi : OPH.ExteriorSelection.isDoublet i = true) :
    chargedCoupling gW i 0 1 = gW / chargedNorm := by
  unfold chargedCoupling
  rw [if_pos ((partner_condition i 0 1).mpr ⟨hi, rfl, rfl⟩)]

/-- **The exact selection rule.**  At a positive declared weak
coupling, a nonzero charged coupling exists exactly on the committed
partner pairs: the two weak components of a committed doublet row,
upper receiving from lower. -/
theorem charged_selection_rule (gW : ℝ) (hgW : 0 < gW) (i : Fin 10)
    (c d : Fin 2) :
    chargedCoupling gW i c d ≠ 0
      ↔ (OPH.ExteriorSelection.isDoublet i = true ∧ c = 0 ∧ d = 1) := by
  unfold chargedCoupling
  split_ifs with hcd
  · constructor
    · intro _
      exact (partner_condition i c d).mp hcd
    · intro _
      exact div_ne_zero hgW.ne' chargedNorm_pos.ne'
  · constructor
    · intro h0
      exact absurd rfl h0
    · intro hpc
      exact absurd ((partner_condition i c d).mpr hpc) hcd

/-- Committed singlet rows are annihilated: every slot pair of a
row with singlet flag gets zero coupling. -/
theorem charged_annihilates_singlets (gW : ℝ) (i : Fin 10)
    (hi : OPH.ExteriorSelection.isDoublet i = false) (c d : Fin 2) :
    chargedCoupling gW i c d = 0 := by
  unfold chargedCoupling
  rw [if_neg]
  intro hcd
  have h2 := (partner_condition i c d).mp hcd
  rw [hi] at h2
  exact Bool.false_ne_true h2.1

/-- Reversed and repeated slot pairs get zero coupling on every row:
the mismatched-pair mutation receipt. -/
theorem charged_mutation_pairs (gW : ℝ) :
    ∀ i : Fin 10,
      chargedCoupling gW i 1 0 = 0
        ∧ chargedCoupling gW i 0 0 = 0
        ∧ chargedCoupling gW i 1 1 = 0 := by
  have hno : ∀ i : Fin 10,
      t3Six i 1 - t3Six i 0 ≠ 6
        ∧ t3Six i 0 - t3Six i 0 ≠ 6
        ∧ t3Six i 1 - t3Six i 1 ≠ 6 := by decide
  intro i
  refine ⟨?_, ?_, ?_⟩ <;>
    · unfold chargedCoupling
      rw [if_neg]
      first
        | exact (hno i).1
        | exact (hno i).2.1
        | exact (hno i).2.2

/-! ## Alignment with the committed action and coefficient rigidity -/

/-- **Alignment of the raising dictionary.**  On the scalar row the
raising shift is exactly the strict shift whose coefficient is the
derived normalization times the dictionary value: the dictionary is
the committed infinitesimal action, with no freedom left. -/
theorem committed_action_alignment_raise (gW gY : ℝ) (h : Doublet) :
    raiseShift gW gY h
      = ![Complex.I
            * ((chargedNorm * chargedCoupling gW 1 0 1 : ℝ) : ℂ)
            * h 1, 0] := by
  have hc : chargedNorm * chargedCoupling gW 1 0 1 = gW := by
    rw [chargedCoupling_partner gW 1 (by decide), mul_comm]
    exact div_mul_cancel₀ gW chargedNorm_pos.ne'
  rw [raiseShift_exact, hc]

/-- **Alignment of the lowering dictionary.**  The lowering shift on
the scalar row is the transpose strict shift with the same derived
coefficient. -/
theorem committed_action_alignment_lower (gW gY : ℝ) (h : Doublet) :
    lowerShift gW gY h
      = ![0, Complex.I
            * ((chargedNorm * chargedCoupling gW 1 0 1 : ℝ) : ℂ)
            * h 0] := by
  have hc : chargedNorm * chargedCoupling gW 1 0 1 = gW := by
    rw [chargedCoupling_partner gW 1 (by decide), mul_comm]
    exact div_mul_cancel₀ gW chargedNorm_pos.ne'
  rw [lowerShift_exact, hc]

/-- **Coefficient rigidity: the drifted-coefficient mutation
receipt.**  Any coefficient satisfying the raising alignment identity
equals the dictionary value `gW / chargedNorm`: a drifted coefficient
breaks an exact identity. -/
theorem coefficient_rigidity (gW gY κ : ℝ)
    (hκ : ∀ h : Doublet,
      raiseShift gW gY h
        = ![Complex.I * ((chargedNorm * κ : ℝ) : ℂ) * h 1, 0]) :
    κ = gW / chargedNorm := by
  have h0 := congrFun (hκ ![0, 1]) 0
  rw [raiseShift_exact] at h0
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one,
    mul_one] at h0
  have h2 : ((gW : ℝ) : ℂ) = ((chargedNorm * κ : ℝ) : ℂ) :=
    mul_left_cancel₀ Complex.I_ne_zero h0
  have h3 : gW = chargedNorm * κ := by exact_mod_cast h2
  rw [h3, mul_comm, mul_div_assoc, div_self chargedNorm_pos.ne',
    mul_one]

/-! ## Charge bookkeeping -/

/-- A nonzero third-isospin step of six sixths steps the committed
electric-charge entry by exactly six: one unit in the committed
`q = 6Y` normalization, from the committed charge columns per
row. -/
theorem coupling_steps_unit_charge :
    ∀ (i : Fin 10) (c d : Fin 2),
      t3Six i c - t3Six i d = 6 → q6EM i c = q6EM i d + 6 := by
  decide

/-- **The charge step of the dictionary.**  Every nonzero charged
coupling changes the committed electric-charge entry by exactly six
sixths: one committed unit of electric charge. -/
theorem charged_coupling_charge_step (gW : ℝ) (i : Fin 10)
    (c d : Fin 2) (h : chargedCoupling gW i c d ≠ 0) :
    q6EM i c = q6EM i d + 6 := by
  unfold chargedCoupling at h
  by_cases hcd : t3Six i c - t3Six i d = 6
  · exact coupling_steps_unit_charge i c d hcd
  · rw [if_neg hcd] at h
    exact absurd rfl h

/-- The per-row doublet form of the charge step on the committed
table. -/
theorem charged_charge_step :
    ∀ i : Fin 10,
      OPH.ExteriorSelection.isDoublet i = true →
        q6EM i 0 = q6EM i 1 + 6 := by
  decide

/-- The commutator of the committed photon-direction action with the
raising shift is exactly `i (gW gY)` times the raising shift. -/
theorem raiseShift_photon_commutator (gW gY : ℝ) (h : Doublet) :
    (fun j => infAction gW gY (photonDirection gW gY)
        (raiseShift gW gY h) j
      - raiseShift gW gY (infAction gW gY (photonDirection gW gY) h) j)
      = fun j => Complex.I * ((gW * gY : ℝ) : ℂ)
          * raiseShift gW gY h j := by
  funext j
  fin_cases j <;>
    · simp [raiseShift, infAction, w1Direction, w2Direction,
        photonDirection]
      try ring

/-- The commutator of the committed photon-direction action with the
lowering shift is exactly `-i (gW gY)` times the lowering shift. -/
theorem lowerShift_photon_commutator (gW gY : ℝ) (h : Doublet) :
    (fun j => infAction gW gY (photonDirection gW gY)
        (lowerShift gW gY h) j
      - lowerShift gW gY (infAction gW gY (photonDirection gW gY) h) j)
      = fun j => -(Complex.I * ((gW * gY : ℝ) : ℂ))
          * lowerShift gW gY h j := by
  funext j
  fin_cases j <;>
    · simp [lowerShift, infAction, w1Direction, w2Direction,
        photonDirection]
      try ring

/-- The commutator eigenvalue in committed neutral units: `gW gY` is
exactly the committed photon normalization `gZ` times the declared
electromagnetic unit `e`. -/
theorem unit_charge_in_neutral_units (gW gY : ℝ) (hgW : 0 < gW) :
    gW * gY = gZ gW gY * eCoupling gW gY := by
  unfold eCoupling
  rw [mul_div_assoc', mul_comm (gZ gW gY), mul_div_assoc,
    div_self (gZ_pos gW gY hgW).ne', mul_one]

/-- **The raising direction carries exactly one committed unit of
electric charge.**  Per unit of the committed photon normalization
`gZ`, the photon commutator eigenvalue of the raising shift is
exactly the declared electromagnetic unit `e`. -/
theorem raise_carries_unit_charge (gW gY : ℝ) (hgW : 0 < gW)
    (h : Doublet) :
    (fun j => infAction gW gY (photonDirection gW gY)
        (raiseShift gW gY h) j
      - raiseShift gW gY (infAction gW gY (photonDirection gW gY) h) j)
      = fun j => Complex.I
          * ((gZ gW gY * eCoupling gW gY : ℝ) : ℂ)
          * raiseShift gW gY h j := by
  rw [raiseShift_photon_commutator gW gY h,
    unit_charge_in_neutral_units gW gY hgW]

/-! ## Consistency joins with the committed neutral directions -/

/-- The third-weak coordinate direction. -/
def a3Direction : EWParam := ⟨0, 0, 1, 0⟩

/-- The commutator of the two charged-direction actions is exactly
`-gW` times the committed third-weak action. -/
theorem charged_commutator_exact (gW gY : ℝ) (h : Doublet) :
    (fun j => infAction gW gY w1Direction
        (infAction gW gY w2Direction h) j
      - infAction gW gY w2Direction
          (infAction gW gY w1Direction h) j)
      = fun j => ((-gW : ℝ) : ℂ)
          * infAction gW gY a3Direction h j := by
  funext j
  fin_cases j <;>
    simp [infAction, w1Direction, w2Direction, a3Direction] <;>
    ring

/-- **The charged commutator lands on the committed neutral
directions.**  The commutator of the two charged-direction actions
decomposes exactly on the committed photon and Z actions, with the
exact coefficients `-(gW gY)/(gW^2 + gY^2)` and
`-(gW^2)/(gW^2 + gY^2)`.  Grade: an exact operator-commutator
identity of the committed infinitesimal action on the scalar-row
carrier; no abstract bracket on the direction space is committed. -/
theorem charged_commutator_neutral_span (gW gY : ℝ)
    (hgg : gW ^ 2 + gY ^ 2 ≠ 0) (h : Doublet) :
    (fun j => infAction gW gY w1Direction
        (infAction gW gY w2Direction h) j
      - infAction gW gY w2Direction
          (infAction gW gY w1Direction h) j)
      = fun j => ((-(gW * gY) / (gW ^ 2 + gY ^ 2) : ℝ) : ℂ)
            * infAction gW gY (photonDirection gW gY) h j
          + ((-(gW ^ 2) / (gW ^ 2 + gY ^ 2) : ℝ) : ℂ)
            * infAction gW gY (zDirection gW gY) h j := by
  have hgc : (gW : ℂ) ^ 2 + (gY : ℂ) ^ 2 ≠ 0 := by
    have hg0 : ((gW ^ 2 + gY ^ 2 : ℝ) : ℂ) ≠ 0 :=
      Complex.ofReal_ne_zero.mpr hgg
    push_cast at hg0
    exact hg0
  funext j
  fin_cases j <;>
    · simp [infAction, w1Direction, w2Direction, photonDirection,
        zDirection]
      field_simp
      try ring

/-- **The ladder commutator.**  The commutator of the raising and
lowering shifts is exactly `2 gW i` times the committed third-weak
action: the su(2) ladder algebra at the exact operator level of the
committed action. -/
theorem ladder_commutator_neutral (gW gY : ℝ) (h : Doublet) :
    (fun j => raiseShift gW gY (lowerShift gW gY h) j
      - lowerShift gW gY (raiseShift gW gY h) j)
      = fun j => ((2 * gW : ℝ) : ℂ) * Complex.I
          * infAction gW gY a3Direction h j := by
  funext j
  fin_cases j <;>
    simp [raiseShift, lowerShift, infAction, w1Direction,
      w2Direction, a3Direction] <;>
    ring

/-! ## The composed receipt -/

/-- **The charged-current dictionary receipt (issue #745, charged
companion of the committed neutral dictionary).**  For every
committed premise bundle `D` (register rows PR-09, PR-10, PR-11,
PR-59 through the base; `gW`, `gY` declared on rows PR-48/PR-54; the
scalar parameters of row PR-49 entering through the committed
vacuum): the two charged directions are orthogonal to both committed
neutral directions and to each other and span exactly the orthogonal
complement of the committed neutral plane, both at the same broken
gauge mass value; the raising and lowering shifts align exactly with
the strict shift whose coefficient is the derived normalization times
the dictionary value; on every doublet row of the bundled selection
mask the dictionary value is `gW / chargedNorm` and nonzero, the
selection rule is an exact equivalence on every row and slot pair,
and every nonzero coupling steps the committed electric charge by
exactly one committed unit; the photon commutator eigenvalues of the
two shifts are exactly `± i gZ e`; the charged commutator lands on
the committed neutral directions and the ladder commutator on the
committed third-weak action with the stated exact coefficients; and
the bundled sector contains a selected singlet row with zero coupling
at every slot pair, while reversed and repeated slot pairs get zero
coupling on every row.  Every clause consumes the bundle.  The
theorem does not consume or discharge register rows PR-36, PR-47,
PR-50, or PR-54: no physical current, no source action, no decay
amplitude, and no flavor or mixing-matrix content is claimed. -/
theorem chargedCurrentDictionary_receipt (D : EWBreakingPremiseData) :
    ((diagPairing w1Direction (photonDirection D.gW D.gY) = 0
        ∧ diagPairing w1Direction (zDirection D.gW D.gY) = 0
        ∧ diagPairing w2Direction (photonDirection D.gW D.gY) = 0
        ∧ diagPairing w2Direction (zDirection D.gW D.gY) = 0
        ∧ diagPairing w1Direction w2Direction = 0)
      ∧ (∀ A : EWParam,
          (diagPairing A (photonDirection D.gW D.gY) = 0
              ∧ diagPairing A (zDirection D.gW D.gY) = 0)
            ↔ ∃! p : ℝ × ℝ, A = ⟨p.1, p.2, 0, 0⟩)
      ∧ gaugeMassForm D.gW D.gY D.vev w1Direction
          = D.gW ^ 2 * D.vev ^ 2 / 4
      ∧ gaugeMassForm D.gW D.gY D.vev w2Direction
          = D.gW ^ 2 * D.vev ^ 2 / 4)
    ∧ ((∀ h : Doublet,
          raiseShift D.gW D.gY h
            = ![Complex.I
                  * ((chargedNorm * chargedCoupling D.gW 1 0 1 : ℝ) : ℂ)
                  * h 1, 0])
      ∧ (∀ h : Doublet,
          lowerShift D.gW D.gY h
            = ![0, Complex.I
                  * ((chargedNorm * chargedCoupling D.gW 1 0 1 : ℝ) : ℂ)
                  * h 0])
      ∧ (∀ i : Fin 10,
          OPH.ExteriorSelection.mem D.base.selectionMask.val i = true →
          OPH.ExteriorSelection.isDoublet i = true →
            chargedCoupling D.gW i 0 1 = D.gW / chargedNorm
              ∧ chargedCoupling D.gW i 0 1 ≠ 0)
      ∧ (∀ (i : Fin 10) (c d : Fin 2),
          chargedCoupling D.gW i c d ≠ 0
            ↔ (OPH.ExteriorSelection.isDoublet i = true
                ∧ c = 0 ∧ d = 1)))
    ∧ ((∀ (i : Fin 10) (c d : Fin 2),
          chargedCoupling D.gW i c d ≠ 0 → q6EM i c = q6EM i d + 6)
      ∧ (∀ h : Doublet,
          (fun j => infAction D.gW D.gY (photonDirection D.gW D.gY)
              (raiseShift D.gW D.gY h) j
            - raiseShift D.gW D.gY
                (infAction D.gW D.gY (photonDirection D.gW D.gY) h) j)
            = fun j => Complex.I
                * ((gZ D.gW D.gY * eCoupling D.gW D.gY : ℝ) : ℂ)
                * raiseShift D.gW D.gY h j)
      ∧ (∀ h : Doublet,
          (fun j => infAction D.gW D.gY (photonDirection D.gW D.gY)
              (lowerShift D.gW D.gY h) j
            - lowerShift D.gW D.gY
                (infAction D.gW D.gY (photonDirection D.gW D.gY) h) j)
            = fun j => -(Complex.I * ((D.gW * D.gY : ℝ) : ℂ))
                * lowerShift D.gW D.gY h j))
    ∧ ((∀ h : Doublet,
          (fun j => infAction D.gW D.gY w1Direction
              (infAction D.gW D.gY w2Direction h) j
            - infAction D.gW D.gY w2Direction
                (infAction D.gW D.gY w1Direction h) j)
            = fun j =>
                ((-(D.gW * D.gY) / (D.gW ^ 2 + D.gY ^ 2) : ℝ) : ℂ)
                  * infAction D.gW D.gY (photonDirection D.gW D.gY) h j
                + ((-(D.gW ^ 2) / (D.gW ^ 2 + D.gY ^ 2) : ℝ) : ℂ)
                  * infAction D.gW D.gY (zDirection D.gW D.gY) h j)
      ∧ (∀ h : Doublet,
          (fun j => raiseShift D.gW D.gY (lowerShift D.gW D.gY h) j
            - lowerShift D.gW D.gY (raiseShift D.gW D.gY h) j)
            = fun j => ((2 * D.gW : ℝ) : ℂ) * Complex.I
                * infAction D.gW D.gY a3Direction h j))
    ∧ ((∃ i : Fin 10,
          OPH.ExteriorSelection.mem D.base.selectionMask.val i = true
            ∧ OPH.ExteriorSelection.isDoublet i = false
            ∧ ∀ c d : Fin 2, chargedCoupling D.gW i c d = 0)
      ∧ (∀ i : Fin 10,
          chargedCoupling D.gW i 1 0 = 0
            ∧ chargedCoupling D.gW i 0 0 = 0
            ∧ chargedCoupling D.gW i 1 1 = 0)) := by
  have hgg : D.gW ^ 2 + D.gY ^ 2 ≠ 0 := by
    have := D.gW_pos
    positivity
  refine ⟨⟨charged_neutral_orthogonal D.gW D.gY,
      fun A => charged_plane_exact D.gW D.gY D.gW_pos D.gY_pos A,
      (charged_directions_mass D.gW D.gY D.vev).1,
      (charged_directions_mass D.gW D.gY D.vev).2⟩,
    ⟨fun h => committed_action_alignment_raise D.gW D.gY h,
      fun h => committed_action_alignment_lower D.gW D.gY h,
      fun i _ hi =>
        ⟨chargedCoupling_partner D.gW i hi,
          (charged_selection_rule D.gW D.gW_pos i 0 1).mpr
            ⟨hi, rfl, rfl⟩⟩,
      fun i c d => charged_selection_rule D.gW D.gW_pos i c d⟩,
    ⟨fun i c d => charged_coupling_charge_step D.gW i c d,
      fun h => raise_carries_unit_charge D.gW D.gY D.gW_pos h,
      fun h => lowerShift_photon_commutator D.gW D.gY h⟩,
    ⟨fun h => charged_commutator_neutral_span D.gW D.gY hgg h,
      fun h => ladder_commutator_neutral D.gW D.gY h⟩,
    ?_, fun i => charged_mutation_pairs D.gW i⟩
  rcases D.base.matter_selection_is_parity_sector with hm | hm
  · exact ⟨2, by rw [hm]; decide, by decide,
      fun c d => charged_annihilates_singlets D.gW 2 (by decide) c d⟩
  · exact ⟨0, by rw [hm]; decide, by decide,
      fun c d => charged_annihilates_singlets D.gW 0 (by decide) c d⟩

/-! ## Nonvacuity -/

/-- The premise bundle is satisfiable: the committed
all-parameters-one instance of the breaking module inhabits it. -/
theorem chargedDictionary_premises_nonvacuous :
    Nonempty EWBreakingPremiseData :=
  ⟨committedEWBreakingPremiseData⟩

/-- Exact witness values at the committed all-parameters-one
instance: the squared partner coupling is `1/2`, the coupling times
the derived normalization is one, and the singlet-row and
reversed-pair couplings are zero.  These are witness values of the
declared parameters, not predictions. -/
theorem committed_instance_values :
    chargedCoupling 1 1 0 1 ^ 2 = 1 / 2
      ∧ chargedCoupling 1 1 0 1 * chargedNorm = 1
      ∧ chargedCoupling 1 0 0 1 = 0
      ∧ chargedCoupling 1 1 1 0 = 0 := by
  have hc : chargedCoupling 1 1 0 1 = 1 / chargedNorm :=
    chargedCoupling_partner 1 1 (by decide)
  refine ⟨?_, ?_,
    charged_annihilates_singlets 1 0 (by decide) 0 1,
    (charged_mutation_pairs 1 1).1⟩
  · rw [hc, div_pow, one_pow, chargedNorm_sq]
  · rw [hc, div_mul_eq_mul_div, one_mul,
      div_self chargedNorm_pos.ne']

end OPH.ChargedCurrentDictionary

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.ChargedCurrentDictionary.charged_neutral_orthogonal
#print axioms OPH.ChargedCurrentDictionary.charged_plane_exact
#print axioms OPH.ChargedCurrentDictionary.charged_directions_mass
#print axioms OPH.ChargedCurrentDictionary.charged_pair_normSq
#print axioms OPH.ChargedCurrentDictionary.chargedNorm_sq
#print axioms OPH.ChargedCurrentDictionary.raiseShift_exact
#print axioms OPH.ChargedCurrentDictionary.lowerShift_exact
#print axioms OPH.ChargedCurrentDictionary.partner_condition
#print axioms OPH.ChargedCurrentDictionary.chargedCoupling_partner
#print axioms OPH.ChargedCurrentDictionary.charged_selection_rule
#print axioms OPH.ChargedCurrentDictionary.charged_annihilates_singlets
#print axioms OPH.ChargedCurrentDictionary.charged_mutation_pairs
#print axioms OPH.ChargedCurrentDictionary.committed_action_alignment_raise
#print axioms OPH.ChargedCurrentDictionary.committed_action_alignment_lower
#print axioms OPH.ChargedCurrentDictionary.coefficient_rigidity
#print axioms OPH.ChargedCurrentDictionary.coupling_steps_unit_charge
#print axioms OPH.ChargedCurrentDictionary.charged_coupling_charge_step
#print axioms OPH.ChargedCurrentDictionary.charged_charge_step
#print axioms OPH.ChargedCurrentDictionary.raiseShift_photon_commutator
#print axioms OPH.ChargedCurrentDictionary.lowerShift_photon_commutator
#print axioms OPH.ChargedCurrentDictionary.unit_charge_in_neutral_units
#print axioms OPH.ChargedCurrentDictionary.raise_carries_unit_charge
#print axioms OPH.ChargedCurrentDictionary.charged_commutator_exact
#print axioms OPH.ChargedCurrentDictionary.charged_commutator_neutral_span
#print axioms OPH.ChargedCurrentDictionary.ladder_commutator_neutral
#print axioms OPH.ChargedCurrentDictionary.chargedCurrentDictionary_receipt
#print axioms OPH.ChargedCurrentDictionary.chargedDictionary_premises_nonvacuous
#print axioms OPH.ChargedCurrentDictionary.committed_instance_values
