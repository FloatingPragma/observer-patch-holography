import ElectroweakBreakingComposition

/-!
# The neutral-current coupling dictionary on the committed grammar
(V3, issue #745)

One structural receipt for the electroweak neutral sector: the exact
coupling dictionary of the photon-shaped and Z-shaped directions on
the weak slots of the committed exterior component table, in the
declared constants of the committed breaking bundle.

The committed inputs.  `Screen/ExteriorSelection.lean` (register row
PR-59) carries the ten-row exterior component table with integer
charges in the `q = 6Y` normalization and the weak-doublet indicator
column; `Screen/ElectroweakBreakingComposition.lean` (issue #735)
carries the four-parameter electroweak direction space `EWParam`, the
declared coupling constants `gW`, `gY` of its premise bundle
(register rows PR-48/PR-54: declared, never selected), the committed
unbroken direction `photonDirection`, the committed neutral broken
direction `zDirection`, and the exact stabilizer and mass-form
receipts at the chosen minimum.

The dictionary.  A weak slot is a row of the committed table together
with a weak component index; its third-isospin entry in sixths is
`t3Six` (`3` and `-3` on the two components of a doublet row, `0` on
a singlet row), its hypercharge entry in sixths is the committed
integer charge column, and its electric-charge entry in sixths is the
sum `q6EM = t3Six + charge`, the `Q = T3 + Y` composition in the
committed integer units.  The neutral coupling of a slot along a
direction `A` is the diagonal action
`(gW * A.a3 * t3Six + gY * A.b * charge) / 6`.  The proved clauses:

* **Alignment.**  On the scalar row of the committed table the
  dictionary reproduces the committed infinitesimal action exactly:
  for every neutral direction (`a1 = a2 = 0`) the committed
  `infAction` is diagonal with entries `i` times the two slot
  couplings of row `1` (`committed_action_alignment`), so the `/6`
  normalization is the committed one, with no freedom left.
* **Directions.**  `photonDirection` lies on the committed unbroken
  stabilizer line and is exactly massless; `zDirection` is orthogonal
  to it under the declared diagonal pairing on the direction space,
  and every neutral direction orthogonal to `photonDirection` is a
  unique multiple of the `zDirection` shape: the Z direction is
  exactly the orthogonal complement of the unbroken line in the
  neutral plane (`z_spans_photon_complement`).
* **Mixing identities.**  With `sinSqW := gY^2 / (gW^2 + gY^2)` as a
  defined abbreviation in the declared constants (never a measured
  value), and `gZ := sqrt (gW^2 + gY^2)`, `e := gW * gY / gZ`: the
  exact identities `e = gW * sinW`, `e = gY * cosW`,
  `sinW^2 = sinSqW`, `cosW^2 = cosSqW`, `sinSqW + cosSqW = 1`, and
  the squared forms `e^2 = gW^2 * sinSqW = gY^2 * cosSqW`, stated as
  squared identities to avoid square-root case analysis.
* **The per-slot dictionary.**  Uniformly over every row and weak
  component of the committed table: the coupling along
  `photonDirection`, normalized by the direction length `gZ`, is
  exactly `e * (q6EM / 6)`, the declared unit `e` times the integer
  electric-charge entry divided by six (the exact factor in the
  committed `q = 6Y` units); the coupling along `zDirection`,
  normalized the same way, is exactly
  `gZ * (t3Six/6 - (q6EM/6) * sinSqW)`, the `T3 - Q sin^2`
  shape times the committed normalization `gZ`.
* **Anti-false-green rows.**  In either committed parity sector of
  the selection mask: two committed rows with different
  electric-charge entries have provably different photon couplings
  (rows `3` and `8` in the even sector, rows `0` and `5` in the odd
  sector), and the neutrino-type slot (the zero-electric-charge weak
  component of the colorless doublet row: row `9` upper component in
  the even sector, row `1` lower component in the odd sector) has
  photon coupling exactly zero and Z coupling exactly
  `(gW^2 + gY^2) / 2` up to sign, nonzero at positive declared
  couplings.
* **The composed receipt.**  `neutralCurrentDictionary_receipt`
  states all clauses from the one committed premise bundle
  `EWBreakingPremiseData`, with the per-slot dictionary quantified
  over the rows of the bundled selection mask.  Nonvacuity is
  inherited from the committed all-parameters-one instance.

Consumed register rows.  The bundle carries the committed
Standard-Model structure premises (register rows PR-09, PR-10, PR-11,
PR-59) through its base, and the declared constants `gW`, `gY`
whose selection stays on register rows PR-48 and PR-54; the potential
parameters of row PR-49 enter only through the committed vacuum used
in the direction receipts.  The charge and doublet columns are the
committed PR-59 table entries.

Boundary and nonclaims.  Every constant is declared: no measured
`sin^2 theta_W`, no running coupling, no scale, no scheme.  The
dictionary is exact tree-level algebra on the committed finite
grammar: there is no physical current operator, no spacetime, no
covariant derivative, no source-produced action (register row PR-54
open), no weak decay, no W-exchange amplitude, no Fermi limit.  The
third-isospin split of a doublet row into two weak components is the
committed weak-dimension bookkeeping of the table; a physical
chirality attachment stays on register row PR-47.  No flavor,
mixing-matrix, or neutrino-mass content is claimed beyond the charge
dictionary (register row PR-50 untouched, register row PR-36
three-family structure not consumed).  The
`zDirection` normalization `gZ` and the unit `e` are exact functions
of the declared constants.  This is the first structural receipt on
issue #745, and it closes nothing: rows OL-M1, OL-M2, OL-M3 keep
their owed physical content.

Falsifier.  The dictionary fails if the committed `infAction` is not
diagonal with the row-1 slot couplings on a neutral direction, if any
mixing identity fails, if any slot's normalized photon coupling
differs from `e * q6EM / 6` or its normalized Z coupling differs from
`gZ * (t3Six/6 - (q6EM/6) * sinSqW)`, if two slots with different
electric charges get equal photon couplings at positive declared
constants, or if a zero-charge slot gets a nonzero photon coupling or
a zero Z coupling.
-/

namespace OPH.NeutralCurrentDictionary

open OPH.ElectroweakBreakingComposition (Doublet vacuum infAction
  EWParam paramNormSq photonDirection zDirection gaugeMassForm
  photon_massless zDirection_paramNormSq stabilizer_exactly_one_line
  EWBreakingPremiseData committedEWBreakingPremiseData)

/-! ## Weak slots of the committed exterior table -/

/-- Six times the third weak isospin of weak component `c` of row
`i`: `3` and `-3` on the two components of a committed doublet row,
`0` on a singlet row.  The split is the committed weak-dimension
bookkeeping of the table; a physical chirality attachment stays on
register row PR-47. -/
def t3Six (i : Fin 10) (c : Fin 2) : ℤ :=
  if OPH.ExteriorSelection.isDoublet i then (if c = 0 then 3 else -3)
  else 0

/-- Six times the electric charge of a weak slot: the `Q = T3 + Y`
composition in the committed integer units, with the hypercharge
entry the committed `q = 6Y` charge column of register row PR-59. -/
def q6EM (i : Fin 10) (c : Fin 2) : ℤ :=
  t3Six i c + OPH.ExteriorSelection.charge i

/-- The electric-charge entry is definitionally the `T3 + Y` sum in
the committed integer units. -/
theorem q6EM_def :
    ∀ (i : Fin 10) (c : Fin 2),
      q6EM i c = t3Six i c + OPH.ExteriorSelection.charge i :=
  fun _ _ => rfl

/-- A weak component index is realized on a row exactly when the row
is a doublet or the index is the first one. -/
abbrev slotValid (i : Fin 10) (c : Fin 2) : Prop :=
  c = 0 ∨ OPH.ExteriorSelection.isDoublet i = true

/-- Slot validity matches the committed weak-dimension column of the
exterior table. -/
theorem slotValid_iff_weakDim :
    ∀ (i : Fin 10) (c : Fin 2),
      slotValid i c ↔ ((c : ℕ) : ℤ) < OPH.ExteriorSelection.weakDim i := by
  decide

/-- The third-isospin column receipt: `3` and `-3` on the two
components of every committed doublet row, `0` on every singlet
row. -/
theorem t3Six_receipt :
    ∀ i : Fin 10,
      (OPH.ExteriorSelection.isDoublet i = true →
        t3Six i 0 = 3 ∧ t3Six i 1 = -3)
        ∧ (OPH.ExteriorSelection.isDoublet i = false →
          ∀ c : Fin 2, t3Six i c = 0) := by
  decide

/-! ## The neutral coupling of a slot along a direction -/

/-- The neutral-direction coupling of weak slot `(i, c)` along the
direction `A`: the diagonal action of the neutral pair `(a3, b)` on a
slot with third isospin `t3Six / 6` and hypercharge `charge / 6`.
The `/6` normalization is pinned to the committed infinitesimal
action by `committed_action_alignment`. -/
noncomputable def neutralCoupling (gW gY : ℝ) (A : EWParam)
    (i : Fin 10) (c : Fin 2) : ℝ :=
  (gW * A.a3 * ((t3Six i c : ℤ) : ℝ)
    + gY * A.b * ((OPH.ExteriorSelection.charge i : ℤ) : ℝ)) / 6

/-- The two slot couplings of the scalar-carrier row (row `1`,
charge `3`) in closed form. -/
theorem neutralCoupling_scalar_row (gW gY : ℝ) (A : EWParam) :
    neutralCoupling gW gY A 1 0 = (gW * A.a3 + gY * A.b) / 2
      ∧ neutralCoupling gW gY A 1 1 = (gY * A.b - gW * A.a3) / 2 := by
  have h0 : t3Six 1 0 = 3 := by decide
  have h1 : t3Six 1 1 = -3 := by decide
  have hc : OPH.ExteriorSelection.charge 1 = 3 := by decide
  constructor
  · unfold neutralCoupling
    rw [h0, hc]
    push_cast
    ring
  · unfold neutralCoupling
    rw [h1, hc]
    push_cast
    ring

/-- **Alignment with the committed action.**  On every neutral
direction (`a1 = a2 = 0`) the committed infinitesimal electroweak
action of `ElectroweakBreakingComposition` is exactly diagonal, with
entries `i` times the two slot couplings of the scalar row of the
dictionary.  The `/6` normalization of `neutralCoupling` is the
committed one. -/
theorem committed_action_alignment (gW gY : ℝ) (A : EWParam)
    (h1 : A.a1 = 0) (h2 : A.a2 = 0) (h : Doublet) :
    infAction gW gY A h
      = ![Complex.I * ((neutralCoupling gW gY A 1 0 : ℝ) : ℂ) * h 0,
          Complex.I * ((neutralCoupling gW gY A 1 1 : ℝ) : ℂ) * h 1] := by
  rw [(neutralCoupling_scalar_row gW gY A).1,
    (neutralCoupling_scalar_row gW gY A).2]
  funext j
  fin_cases j <;>
    · simp [infAction, h1, h2]
      ring

/-! ## The two neutral directions: unbroken line and its orthogonal
complement -/

/-- The declared diagonal pairing on the electroweak direction
space. -/
def diagPairing (A B : EWParam) : ℝ :=
  A.a1 * B.a1 + A.a2 * B.a2 + A.a3 * B.a3 + A.b * B.b

/-- The diagonal pairing polarizes the committed squared direction
length. -/
theorem diagPairing_self (A : EWParam) :
    diagPairing A A = paramNormSq A := by
  simp only [diagPairing, paramNormSq]
  ring

/-- The committed photon direction lies on the committed unbroken
stabilizer line of the chosen minimum. -/
theorem photon_on_unbroken_line (gW gY vev : ℝ) (hgW : 0 < gW)
    (hvev : 0 < vev) :
    infAction gW gY (photonDirection gW gY) (vacuum vev) = 0 := by
  rw [stabilizer_exactly_one_line gW gY vev hgW hvev]
  refine ⟨1, ?_, ?_⟩
  · ext <;> simp [photonDirection]
  · intro t ht
    have hb : (photonDirection gW gY).b = gW * t := by rw [ht]
    simp only [photonDirection] at hb
    exact (mul_left_cancel₀ hgW.ne' (by rw [← hb]; ring)).symm

/-- **The exact orthogonality receipt.**  The committed Z direction
is orthogonal to the committed photon direction under the declared
diagonal pairing. -/
theorem z_photon_diag_orthogonal (gW gY : ℝ) :
    diagPairing (zDirection gW gY) (photonDirection gW gY) = 0 := by
  simp only [diagPairing, zDirection, photonDirection]
  ring

/-- **The Z direction spans the orthogonal complement.**  Every
neutral direction orthogonal to the committed photon direction under
the declared diagonal pairing is a unique multiple of the committed
`zDirection` shape `(gW, -gY)`. -/
theorem z_spans_photon_complement (gW gY : ℝ) (hgW : 0 < gW)
    (A : EWParam) (h1 : A.a1 = 0) (h2 : A.a2 = 0)
    (horth : diagPairing A (photonDirection gW gY) = 0) :
    ∃! t : ℝ, A = ⟨0, 0, gW * t, -(gY * t)⟩ := by
  simp only [diagPairing, photonDirection] at horth
  refine ⟨A.a3 / gW, ?_, ?_⟩
  · ext
    · exact h1
    · exact h2
    · field_simp
    · rw [h1, h2] at horth
      field_simp
      nlinarith [horth]
  · intro t ht
    have ha : A.a3 = gW * t := by rw [ht]
    rw [ha]
    field_simp

/-- The committed Z direction is the `t = 1` member of the complement
line, at its committed squared length. -/
theorem zDirection_on_complement (gW gY : ℝ) :
    zDirection gW gY = (⟨0, 0, gW * 1, -(gY * 1)⟩ : EWParam)
      ∧ paramNormSq (zDirection gW gY) = gW ^ 2 + gY ^ 2 := by
  refine ⟨?_, zDirection_paramNormSq gW gY⟩
  ext <;> simp [zDirection]

/-! ## The mixing identities in the declared constants -/

/-- The defined mixing abbreviation `sin^2` in the declared
constants: `gY^2 / (gW^2 + gY^2)`.  A definition, never a measured
value. -/
noncomputable def sinSqW (gW gY : ℝ) : ℝ := gY ^ 2 / (gW ^ 2 + gY ^ 2)

/-- The defined mixing abbreviation `cos^2` in the declared
constants. -/
noncomputable def cosSqW (gW gY : ℝ) : ℝ := gW ^ 2 / (gW ^ 2 + gY ^ 2)

/-- The neutral-direction normalization: the common direction length
of the committed photon and Z directions. -/
noncomputable def gZ (gW gY : ℝ) : ℝ := Real.sqrt (gW ^ 2 + gY ^ 2)

/-- The declared electromagnetic unit `e = gW * gY / gZ`. -/
noncomputable def eCoupling (gW gY : ℝ) : ℝ := gW * gY / gZ gW gY

/-- The sine of the mixing shape in the declared constants. -/
noncomputable def sinW (gW gY : ℝ) : ℝ := gY / gZ gW gY

/-- The cosine of the mixing shape in the declared constants. -/
noncomputable def cosW (gW gY : ℝ) : ℝ := gW / gZ gW gY

theorem gZ_pos (gW gY : ℝ) (hgW : 0 < gW) : 0 < gZ gW gY :=
  Real.sqrt_pos.mpr (by positivity)

theorem gZ_sq (gW gY : ℝ) : gZ gW gY ^ 2 = gW ^ 2 + gY ^ 2 :=
  Real.sq_sqrt (by positivity)

/-- The direction lengths of the committed photon and Z directions
are both exactly `gZ`. -/
theorem gZ_is_direction_length (gW gY : ℝ) :
    Real.sqrt (paramNormSq (photonDirection gW gY)) = gZ gW gY
      ∧ Real.sqrt (paramNormSq (zDirection gW gY)) = gZ gW gY := by
  constructor
  · have hp : paramNormSq (photonDirection gW gY) = gW ^ 2 + gY ^ 2 := by
      simp only [paramNormSq, photonDirection]
      ring
    rw [hp]
    rfl
  · rw [zDirection_paramNormSq]
    rfl

/-- `e` in its defining exact form. -/
theorem e_def_exact (gW gY : ℝ) :
    eCoupling gW gY = gW * gY / Real.sqrt (gW ^ 2 + gY ^ 2) := rfl

/-- The exact identity `e = gW * sinW`. -/
theorem e_eq_gW_sinW (gW gY : ℝ) :
    eCoupling gW gY = gW * sinW gW gY := by
  unfold eCoupling sinW
  ring

/-- The exact identity `e = gY * cosW`. -/
theorem e_eq_gY_cosW (gW gY : ℝ) :
    eCoupling gW gY = gY * cosW gW gY := by
  unfold eCoupling cosW
  ring

/-- The squared identity `sinW^2 = sinSqW`, stated in squared form to
avoid square-root case analysis. -/
theorem sinW_sq (gW gY : ℝ) : sinW gW gY ^ 2 = sinSqW gW gY := by
  unfold sinW sinSqW
  rw [div_pow, gZ_sq]

/-- The squared identity `cosW^2 = cosSqW`. -/
theorem cosW_sq (gW gY : ℝ) : cosW gW gY ^ 2 = cosSqW gW gY := by
  unfold cosW cosSqW
  rw [div_pow, gZ_sq]

/-- The Pythagorean identity of the defined mixing abbreviations at a
positive declared weak coupling. -/
theorem sinSq_add_cosSq (gW gY : ℝ) (hgW : 0 < gW) :
    sinSqW gW gY + cosSqW gW gY = 1 := by
  have hg2 : gW ^ 2 + gY ^ 2 ≠ 0 := by positivity
  unfold sinSqW cosSqW
  field_simp
  ring

/-- The squared identity `e^2 = gW^2 * sinSqW`, stated in squared
form. -/
theorem e_sq_sin (gW gY : ℝ) :
    eCoupling gW gY ^ 2 = gW ^ 2 * sinSqW gW gY := by
  unfold eCoupling sinSqW
  rw [div_pow, gZ_sq]
  ring

/-- The squared identity `e^2 = gY^2 * cosSqW`, stated in squared
form. -/
theorem e_sq_cos (gW gY : ℝ) :
    eCoupling gW gY ^ 2 = gY ^ 2 * cosSqW gW gY := by
  unfold eCoupling cosSqW
  rw [div_pow, gZ_sq]
  ring

/-! ## The per-slot dictionary -/

/-- The raw photon-direction coupling of any slot:
`gW * gY * (q6EM / 6)`. -/
theorem photon_coupling_exact (gW gY : ℝ) (i : Fin 10) (c : Fin 2) :
    neutralCoupling gW gY (photonDirection gW gY) i c
      = gW * gY * (((q6EM i c : ℤ) : ℝ) / 6) := by
  unfold neutralCoupling q6EM
  simp only [photonDirection]
  push_cast
  ring

/-- **The photon dictionary.**  The photon-direction coupling of any
slot, normalized by the direction length `gZ`, is exactly the
declared unit `e` times the integer electric-charge entry divided by
six: the exact factor in the committed `q = 6Y` units. -/
theorem photon_coupling_e (gW gY : ℝ) (i : Fin 10) (c : Fin 2) :
    neutralCoupling gW gY (photonDirection gW gY) i c / gZ gW gY
      = eCoupling gW gY * (((q6EM i c : ℤ) : ℝ) / 6) := by
  rw [photon_coupling_exact]
  unfold eCoupling
  ring

/-- The raw Z-direction coupling of any slot:
`(gW^2 * t3Six - gY^2 * charge) / 6`. -/
theorem z_coupling_exact (gW gY : ℝ) (i : Fin 10) (c : Fin 2) :
    neutralCoupling gW gY (zDirection gW gY) i c
      = (gW ^ 2 * ((t3Six i c : ℤ) : ℝ)
          - gY ^ 2 * ((OPH.ExteriorSelection.charge i : ℤ) : ℝ)) / 6 := by
  unfold neutralCoupling
  simp only [zDirection]
  ring

/-- The Z coupling of any slot in the `T3 - Q sin^2` shape, with the
unnormalized prefactor `gW^2 + gY^2`. -/
theorem z_coupling_shape (gW gY : ℝ) (hg2 : gW ^ 2 + gY ^ 2 ≠ 0)
    (i : Fin 10) (c : Fin 2) :
    neutralCoupling gW gY (zDirection gW gY) i c
      = (gW ^ 2 + gY ^ 2)
          * (((t3Six i c : ℤ) : ℝ) / 6
            - (((q6EM i c : ℤ) : ℝ) / 6) * sinSqW gW gY) := by
  rw [z_coupling_exact]
  unfold sinSqW q6EM
  push_cast
  field_simp
  ring

/-- **The Z dictionary.**  The Z-direction coupling of any slot,
normalized by the direction length `gZ`, is exactly
`gZ * (T3 - Q * sinSqW)` in the committed sixths entries: the
`T3 - Q sin^2` shape times the committed normalization. -/
theorem z_coupling_gZ (gW gY : ℝ) (hgW : 0 < gW) (i : Fin 10)
    (c : Fin 2) :
    neutralCoupling gW gY (zDirection gW gY) i c / gZ gW gY
      = gZ gW gY
          * (((t3Six i c : ℤ) : ℝ) / 6
            - (((q6EM i c : ℤ) : ℝ) / 6) * sinSqW gW gY) := by
  have hg2 : gW ^ 2 + gY ^ 2 ≠ 0 := by positivity
  have hz : gZ gW gY ≠ 0 := (gZ_pos gW gY hgW).ne'
  rw [z_coupling_shape gW gY hg2 i c, ← gZ_sq gW gY]
  field_simp

/-! ## Anti-false-green receipts -/

/-- A slot with zero electric-charge entry has photon coupling
exactly zero, at every value of the declared constants. -/
theorem photon_zero_of_neutral (gW gY : ℝ) (i : Fin 10) (c : Fin 2)
    (h : q6EM i c = 0) :
    neutralCoupling gW gY (photonDirection gW gY) i c = 0 := by
  rw [photon_coupling_exact, h]
  norm_num

/-- **Charge separation.**  At positive declared constants, slots
with different electric-charge entries have different photon
couplings. -/
theorem photon_separates_charges (gW gY : ℝ) (hgW : 0 < gW)
    (hgY : 0 < gY) (i j : Fin 10) (c d : Fin 2)
    (hne : q6EM i c ≠ q6EM j d) :
    neutralCoupling gW gY (photonDirection gW gY) i c
      ≠ neutralCoupling gW gY (photonDirection gW gY) j d := by
  rw [photon_coupling_exact, photon_coupling_exact]
  intro h
  apply hne
  have hgg : gW * gY ≠ 0 := by positivity
  have h6 := mul_left_cancel₀ hgg h
  have : ((q6EM i c : ℤ) : ℝ) = ((q6EM j d : ℤ) : ℝ) := by linarith
  exact_mod_cast this

/-- The neutrino-type slot of the even parity sector: row `9` (the
colorless doublet of charge `-3`), upper component.  Zero electric
charge, zero photon coupling, Z coupling exactly
`(gW^2 + gY^2) / 2`. -/
theorem neutrino_slot_even (gW gY : ℝ) :
    OPH.ExteriorSelection.mem OPH.ExteriorSelection.evenMask 9 = true
      ∧ OPH.ExteriorSelection.isColored 9 = false
      ∧ OPH.ExteriorSelection.isDoublet 9 = true
      ∧ q6EM 9 0 = 0
      ∧ neutralCoupling gW gY (photonDirection gW gY) 9 0 = 0
      ∧ neutralCoupling gW gY (zDirection gW gY) 9 0
          = (gW ^ 2 + gY ^ 2) / 2 := by
  have h3 : t3Six 9 0 = 3 := by decide
  have hc : OPH.ExteriorSelection.charge 9 = -3 := by decide
  refine ⟨by decide, by decide, by decide, by decide,
    photon_zero_of_neutral gW gY 9 0 (by decide), ?_⟩
  rw [z_coupling_exact, h3, hc]
  push_cast
  ring

/-- The neutrino-type slot of the odd parity sector: row `1` (the
colorless doublet of charge `3`), lower component.  Zero electric
charge, zero photon coupling, Z coupling exactly
`-(gW^2 + gY^2) / 2`. -/
theorem neutrino_slot_odd (gW gY : ℝ) :
    OPH.ExteriorSelection.mem OPH.ExteriorSelection.oddMask 1 = true
      ∧ OPH.ExteriorSelection.isColored 1 = false
      ∧ OPH.ExteriorSelection.isDoublet 1 = true
      ∧ q6EM 1 1 = 0
      ∧ neutralCoupling gW gY (photonDirection gW gY) 1 1 = 0
      ∧ neutralCoupling gW gY (zDirection gW gY) 1 1
          = -((gW ^ 2 + gY ^ 2) / 2) := by
  have h3 : t3Six 1 1 = -3 := by decide
  have hc : OPH.ExteriorSelection.charge 1 = 3 := by decide
  refine ⟨by decide, by decide, by decide, by decide,
    photon_zero_of_neutral gW gY 1 1 (by decide), ?_⟩
  rw [z_coupling_exact, h3, hc]
  push_cast
  ring

/-! ## The composed receipt -/

/-- **The neutral-current dictionary receipt (issue #745, first
structural receipt).**  For every committed premise bundle `D`
(register rows PR-09, PR-10, PR-11, PR-59 through the base; `gW`,
`gY` declared on rows PR-48/PR-54; the scalar parameters of row PR-49
entering through the committed vacuum): the committed photon
direction lies on the committed unbroken stabilizer line and is
exactly massless; the committed Z direction is orthogonal to it under
the declared diagonal pairing and spans the orthogonal complement in
the neutral plane; the mixing identities hold exactly in the declared
constants, with the squared forms stated as such; every weak slot of
every row of the bundled selection mask has normalized photon
coupling exactly `e * q6EM / 6` and normalized Z coupling exactly
`gZ * (t3Six/6 - (q6EM/6) * sinSqW)`; the bundled sector contains two
rows with different electric-charge entries and provably different
photon couplings; and the bundled sector contains a zero-charge weak
slot with photon coupling exactly zero and Z coupling nonzero.  Every
clause consumes the bundle.  The theorem does not consume or
discharge register rows PR-36, PR-47, PR-50, or PR-54: no physical
current, no measured mixing angle, no weak decay, and no flavor or
neutrino content beyond the charge dictionary is claimed. -/
theorem neutralCurrentDictionary_receipt (D : EWBreakingPremiseData) :
    (infAction D.gW D.gY (photonDirection D.gW D.gY) (vacuum D.vev) = 0
      ∧ gaugeMassForm D.gW D.gY D.vev (photonDirection D.gW D.gY) = 0
      ∧ diagPairing (zDirection D.gW D.gY) (photonDirection D.gW D.gY)
          = 0
      ∧ (∀ A : EWParam, A.a1 = 0 → A.a2 = 0 →
          diagPairing A (photonDirection D.gW D.gY) = 0 →
          ∃! t : ℝ, A = ⟨0, 0, D.gW * t, -(D.gY * t)⟩))
    ∧ (sinSqW D.gW D.gY + cosSqW D.gW D.gY = 1
      ∧ eCoupling D.gW D.gY = D.gW * sinW D.gW D.gY
      ∧ eCoupling D.gW D.gY = D.gY * cosW D.gW D.gY
      ∧ eCoupling D.gW D.gY ^ 2 = D.gW ^ 2 * sinSqW D.gW D.gY
      ∧ eCoupling D.gW D.gY ^ 2 = D.gY ^ 2 * cosSqW D.gW D.gY)
    ∧ (∀ i : Fin 10,
        OPH.ExteriorSelection.mem D.base.selectionMask.val i = true →
        ∀ c : Fin 2,
          neutralCoupling D.gW D.gY (photonDirection D.gW D.gY) i c
              / gZ D.gW D.gY
            = eCoupling D.gW D.gY * (((q6EM i c : ℤ) : ℝ) / 6)
          ∧ neutralCoupling D.gW D.gY (zDirection D.gW D.gY) i c
              / gZ D.gW D.gY
            = gZ D.gW D.gY
                * (((t3Six i c : ℤ) : ℝ) / 6
                  - (((q6EM i c : ℤ) : ℝ) / 6) * sinSqW D.gW D.gY))
    ∧ (∃ i j : Fin 10,
        OPH.ExteriorSelection.mem D.base.selectionMask.val i = true
          ∧ OPH.ExteriorSelection.mem D.base.selectionMask.val j = true
          ∧ q6EM i 0 ≠ q6EM j 0
          ∧ neutralCoupling D.gW D.gY (photonDirection D.gW D.gY) i 0
              ≠ neutralCoupling D.gW D.gY (photonDirection D.gW D.gY)
                  j 0)
    ∧ (∃ (i : Fin 10) (c : Fin 2),
        OPH.ExteriorSelection.mem D.base.selectionMask.val i = true
          ∧ slotValid i c
          ∧ q6EM i c = 0
          ∧ neutralCoupling D.gW D.gY (photonDirection D.gW D.gY) i c
              = 0
          ∧ neutralCoupling D.gW D.gY (zDirection D.gW D.gY) i c
              ≠ 0) := by
  have hpos : (0 : ℝ) < D.gW ^ 2 + D.gY ^ 2 := by
    have := D.gW_pos
    positivity
  refine ⟨⟨photon_on_unbroken_line D.gW D.gY D.vev D.gW_pos D.vev_pos,
      photon_massless D.gW D.gY D.vev,
      z_photon_diag_orthogonal D.gW D.gY,
      fun A h1 h2 horth =>
        z_spans_photon_complement D.gW D.gY D.gW_pos A h1 h2 horth⟩,
    ⟨sinSq_add_cosSq D.gW D.gY D.gW_pos,
      e_eq_gW_sinW D.gW D.gY, e_eq_gY_cosW D.gW D.gY,
      e_sq_sin D.gW D.gY, e_sq_cos D.gW D.gY⟩,
    fun i _ c =>
      ⟨photon_coupling_e D.gW D.gY i c,
        z_coupling_gZ D.gW D.gY D.gW_pos i c⟩,
    ?_, ?_⟩
  · rcases D.base.matter_selection_is_parity_sector with h | h
    · exact ⟨3, 8, by rw [h]; decide, by rw [h]; decide, by decide,
        photon_separates_charges D.gW D.gY D.gW_pos D.gY_pos 3 8 0 0
          (by decide)⟩
    · exact ⟨0, 5, by rw [h]; decide, by rw [h]; decide, by decide,
        photon_separates_charges D.gW D.gY D.gW_pos D.gY_pos 0 5 0 0
          (by decide)⟩
  · rcases D.base.matter_selection_is_parity_sector with h | h
    · refine ⟨9, 0, by rw [h]; decide, by decide, by decide,
        (neutrino_slot_even D.gW D.gY).2.2.2.2.1, ?_⟩
      rw [(neutrino_slot_even D.gW D.gY).2.2.2.2.2]
      have : (0 : ℝ) < (D.gW ^ 2 + D.gY ^ 2) / 2 := by positivity
      exact this.ne'
    · refine ⟨1, 1, by rw [h]; decide, by decide, by decide,
        (neutrino_slot_odd D.gW D.gY).2.2.2.2.1, ?_⟩
      rw [(neutrino_slot_odd D.gW D.gY).2.2.2.2.2]
      have : (0 : ℝ) < (D.gW ^ 2 + D.gY ^ 2) / 2 := by positivity
      exact neg_ne_zero.mpr this.ne'

/-! ## Nonvacuity -/

/-- The premise bundle is satisfiable: the committed
all-parameters-one instance of the breaking module inhabits it. -/
theorem dictionary_premises_nonvacuous :
    Nonempty EWBreakingPremiseData :=
  ⟨committedEWBreakingPremiseData⟩

/-- Exact witness values at the committed all-parameters-one
instance: the defined mixing abbreviations are both `1/2`, the
even-sector neutrino slot has photon coupling `0` and Z coupling `1`.
These are witness values of the declared parameters, not
predictions. -/
theorem committed_instance_values :
    sinSqW 1 1 = 1 / 2
      ∧ cosSqW 1 1 = 1 / 2
      ∧ neutralCoupling 1 1 (photonDirection 1 1) 9 0 = 0
      ∧ neutralCoupling 1 1 (zDirection 1 1) 9 0 = 1 := by
  refine ⟨by unfold sinSqW; norm_num, by unfold cosSqW; norm_num,
    (neutrino_slot_even 1 1).2.2.2.2.1, ?_⟩
  rw [(neutrino_slot_even 1 1).2.2.2.2.2]
  norm_num

end OPH.NeutralCurrentDictionary

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.NeutralCurrentDictionary.slotValid_iff_weakDim
#print axioms OPH.NeutralCurrentDictionary.t3Six_receipt
#print axioms OPH.NeutralCurrentDictionary.committed_action_alignment
#print axioms OPH.NeutralCurrentDictionary.photon_on_unbroken_line
#print axioms OPH.NeutralCurrentDictionary.z_photon_diag_orthogonal
#print axioms OPH.NeutralCurrentDictionary.z_spans_photon_complement
#print axioms OPH.NeutralCurrentDictionary.zDirection_on_complement
#print axioms OPH.NeutralCurrentDictionary.gZ_is_direction_length
#print axioms OPH.NeutralCurrentDictionary.e_eq_gW_sinW
#print axioms OPH.NeutralCurrentDictionary.e_eq_gY_cosW
#print axioms OPH.NeutralCurrentDictionary.sinW_sq
#print axioms OPH.NeutralCurrentDictionary.cosW_sq
#print axioms OPH.NeutralCurrentDictionary.sinSq_add_cosSq
#print axioms OPH.NeutralCurrentDictionary.e_sq_sin
#print axioms OPH.NeutralCurrentDictionary.e_sq_cos
#print axioms OPH.NeutralCurrentDictionary.photon_coupling_exact
#print axioms OPH.NeutralCurrentDictionary.photon_coupling_e
#print axioms OPH.NeutralCurrentDictionary.z_coupling_exact
#print axioms OPH.NeutralCurrentDictionary.z_coupling_shape
#print axioms OPH.NeutralCurrentDictionary.z_coupling_gZ
#print axioms OPH.NeutralCurrentDictionary.photon_zero_of_neutral
#print axioms OPH.NeutralCurrentDictionary.photon_separates_charges
#print axioms OPH.NeutralCurrentDictionary.neutrino_slot_even
#print axioms OPH.NeutralCurrentDictionary.neutrino_slot_odd
#print axioms OPH.NeutralCurrentDictionary.neutralCurrentDictionary_receipt
#print axioms OPH.NeutralCurrentDictionary.dictionary_premises_nonvacuous
#print axioms OPH.NeutralCurrentDictionary.committed_instance_values
