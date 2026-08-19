import SMStructureAdequacySurface

/-!
# Exact character descent bookkeeping for the global form (V3, issue #739)

The committed record carries the diagonal `Z6` tensor-action kernel of
`Screen/TraceBalancedKernel.lean`, its cyclic parametrization `genPow`
and the phase convention `phase (k, l, r) (q, t, d) = 2kt + 3ld + rq`
in sixths of a turn of `Screen/Z6Descent.lean`, the ten-row exterior
component table of `Screen/ExteriorSelection.lean` with its central
weights `componentWeight` of `Screen/QuantumMatterIntegration.lean`
(register row PR-59), and the six-axis loop-to-kernel identification
of register row PR-35.  This module states and kernel-checks the exact
finite character-descent piece of the register row PR-46 bundle on
those objects.

CONTENT.

* `centerChar`: the center character of a weight label `(q, t, d)` on
  the committed lattice `W = Z/6 x Z/3 x Z/2` in the committed
  normalization `q = 6Y`: the value `2t + 3d + q` modulo six, equal to
  the committed phase of the kernel generator
  `(omega_3 I_3, -I_2, e^{i pi/3})` (`centerChar_eq_phase_gen`), and
  additive on the label lattice (`centerCharHom`).
* `descent_criterion`: the descent iff.  A label is fixed by every
  element of the tensor-action kernel exactly when its center
  character vanishes, that is, exactly when the committed integer
  congruence `2t + 3d + q = 0 mod 6` holds
  (`descent_integer_congruence`); through the committed embedding
  `genPow` the whole kernel acts by multiples of the center character
  (`phase_genPow`).
* `committed_table_descends`: every row of the committed matter table
  satisfies the congruence: the ten exterior component rows through
  their committed central weights, the same rows as exact integer
  congruences on the committed charge and bidegree columns, and the
  realized weight table of the descent module including both carrier
  blocks and the scalar.
* `control_label_blocked` and `descent_count`: the criterion is not
  vacuous and not universal.  The fractionally charged singlet label
  `(1, 0, 0)` has center character one and is moved by the kernel
  generator, and on the committed 36-element label lattice the
  descending labels number exactly six, one sixth of the lattice, with
  every center-character fibre of size six.
* `masked_selection_descends`: the composed receipt on one antecedent
  bundle.  For every `SMStructurePremiseData` bundle, every component
  selected by the bundled mask satisfies the congruence and descends
  through the full kernel; the proof threads the bundled
  diagonal-generator receipt through the descent iff.

PREMISES.  The label lattice, phase convention, kernel, matter table,
and selection grammar are the committed objects of register rows
PR-35 and PR-59; the composed receipt consumes one bundle of the
committed Standard Model structure surface.  No premise is discharged
here.

BOUNDARY.  This is the finite character-lattice half of the register
row PR-46 bundle only.  No complete character category, no instanton
or topological normalization, no laboratory flux or line-operator
attachment, and no physical quotient selection is supplied; register
row PR-46 stays open with exactly this piece carved out, and the
global-form row OL-G2 stays partial.  The congruence classifies
finite labels on the committed lattice; it is not a derivation of the
physical global gauge-group form.  All statements are exact finite
arithmetic over `Int` and `ZMod`.
-/

namespace OPH.GlobalFormCharacterDescent

open OPH.TraceBalancedKernel (C gen tensorCharFun)
open OPH.Z6Descent (W phase genPow)
open OPH.QuantumMatterIntegration (componentWeight)

/-- The center character of a weight label `(q, t, d)` on the committed
lattice `W = Z/6 x Z/3 x Z/2` with `q = 6Y`: the value `2t + 3d + q`
in sixths of a turn.  This is the character by which the committed
kernel generator acts on the label. -/
def centerChar (w : W) : ZMod 6 :=
  ((2 * w.2.1.val + 3 * w.2.2.val + w.1.val : ℕ) : ZMod 6)

/-- A label descends to the quotient by the committed diagonal kernel:
every central parameter in the tensor-action kernel fixes it. -/
def Descends (w : W) : Prop :=
  ∀ c : C, tensorCharFun c = 0 → phase c w = 0

instance (w : W) : Decidable (Descends w) := by
  unfold Descends; infer_instance

/-- The center character is the committed phase of the kernel generator
`(omega_3 I_3, -I_2, e^{i pi/3})`: the definition reuses the committed
convention rather than introducing one. -/
theorem centerChar_eq_phase_gen : ∀ w : W, centerChar w = phase gen w := by
  decide

/-- The center character is additive on the committed label lattice. -/
theorem centerChar_add :
    ∀ v w : W, centerChar (v + w) = centerChar v + centerChar w := by
  decide

/-- The center character as an additive character of the label lattice. -/
def centerCharHom : W →+ ZMod 6 :=
  ⟨⟨centerChar, by decide⟩, fun v w => centerChar_add v w⟩

/-- Through the committed cyclic embedding `genPow` of `Z/6` into the
central parameters, the kernel element of exponent `n` acts on a label
by `n` times its center character: the whole kernel action on labels is
the center character. -/
theorem phase_genPow :
    ∀ (n : ZMod 6) (w : W), phase (genPow n) w = n * centerChar w := by
  decide

/-- **The descent criterion.**  A label on the committed lattice is
fixed by every element of the tensor-action kernel exactly when its
center character vanishes: the committed congruence `2t + 3d + q = 0`
modulo six is equivalent to descent. -/
theorem descent_criterion : ∀ w : W, Descends w ↔ centerChar w = 0 := by
  decide

/-- The center character in explicit integer-congruence form on the
committed value columns: `2t + 3d + q = 0 mod 6` over the natural
number entries of the label. -/
theorem descent_integer_congruence :
    ∀ w : W,
      centerChar w = 0 ↔ (2 * w.2.1.val + 3 * w.2.2.val + w.1.val) % 6 = 0 := by
  decide

/-- Every one of the ten committed exterior component rows satisfies the
congruence through its committed central weight; row by row over the
committed table. -/
theorem exterior_rows_descend :
    ∀ i : Fin 10, centerChar (componentWeight i) = 0 := by
  decide +kernel

/-- The same ten rows as exact integer congruences on the committed
charge and bidegree columns: `2t + 3d + q = 0 mod 6` with `t` the color
exterior degree, `d` the weak exterior degree, and `q` the committed
integer charge `q = 6Y` of the row. -/
theorem exterior_rows_integer_congruence :
    ∀ i : Fin 10,
      (2 * ((OPH.ExteriorComponentBridge.componentDegree i).1.val : ℤ)
        + 3 * ((OPH.ExteriorComponentBridge.componentDegree i).2.val : ℤ)
        + OPH.ExteriorSelection.charge i) % 6 = 0 := by
  decide +kernel

/-- The realized weight table of the descent module satisfies the
congruence on every row, including both carrier blocks and the
scalar. -/
theorem realized_weights_descend :
    ∀ w ∈ OPH.Z6Descent.realizedWeights, centerChar w = 0 := by
  decide

/-- **The committed matter table descends (packaged).**  The ten
exterior component rows through their committed central weights, the
same rows as integer congruences on the committed charge and bidegree
columns, and the realized weight table of the descent module all
satisfy the committed congruence; with `descent_criterion` every row
descends through the full kernel.  The fifteen-state generation of
either parity sector is carried by five of these rows, so every
fifteen-state label satisfies the congruence. -/
theorem committed_table_descends :
    (∀ i : Fin 10, centerChar (componentWeight i) = 0)
      ∧ (∀ i : Fin 10,
          (2 * ((OPH.ExteriorComponentBridge.componentDegree i).1.val : ℤ)
            + 3 * ((OPH.ExteriorComponentBridge.componentDegree i).2.val : ℤ)
            + OPH.ExteriorSelection.charge i) % 6 = 0)
      ∧ (∀ w ∈ OPH.Z6Descent.realizedWeights, centerChar w = 0) :=
  ⟨exterior_rows_descend, exterior_rows_integer_congruence,
    realized_weights_descend⟩

/-- **The control label.**  The fractionally charged singlet label
`(1, 0, 0)` violates the congruence: its center character is one, the
kernel generator is in the kernel, and the generator moves the label.
The criterion is not vacuous. -/
theorem control_label_blocked :
    centerChar ((1, 0, 0) : W) = 1
      ∧ tensorCharFun gen = 0
      ∧ phase gen ((1, 0, 0) : W) ≠ 0 := by
  decide

/-- The control label does not descend. -/
theorem control_label_not_descends : ¬ Descends ((1, 0, 0) : W) := by
  decide

/-- **The counting receipt.**  On the committed 36-element label
lattice in the committed hypercharge normalization `q = 6Y`, exactly
six labels descend: one sixth of the lattice exactly. -/
theorem descent_count :
    (Finset.univ.filter fun w : W => centerChar w = 0).card = 6
      ∧ Fintype.card W = 36
      ∧ 6 * (Finset.univ.filter fun w : W => centerChar w = 0).card
          = Fintype.card W := by
  refine ⟨by decide, by decide, by decide⟩

/-- The exact fibre count behind the one-sixth fraction: every value of
the center character is taken by exactly six labels, so the descending
fraction is one sixth because the character equidistributes the
lattice. -/
theorem centerChar_fibres_equal :
    ∀ z : ZMod 6,
      (Finset.univ.filter fun w : W => centerChar w = z).card = 6 := by
  decide

/-- **The composed receipt on one antecedent bundle.**  For every
premise bundle of the committed Standard Model structure surface,
every component selected by the bundled mask satisfies the committed
congruence and descends through the full tensor-action kernel.  The
proof consumes the bundled diagonal-generator receipt and threads it
through the descent criterion; the bundle grammar is register row
PR-59, and the physical reading of the kernel is register row
PR-35. -/
theorem masked_selection_descends
    (D : OPH.SMStructureAdequacySurface.SMStructurePremiseData) :
    ∀ i : Fin 10,
      OPH.ExteriorSelection.mem D.selectionMask.val i = true →
        centerChar (componentWeight i) = 0
          ∧ Descends (componentWeight i) := by
  intro i hi
  have hgen := D.masked_diagonal_z6_fixes i hi
  have h0 : centerChar (componentWeight i) = 0 :=
    (centerChar_eq_phase_gen (componentWeight i)).trans hgen
  exact ⟨h0, (descent_criterion (componentWeight i)).mpr h0⟩

/-- **The character-descent receipt (issue #739, OL-G2 evidence).**
The descent criterion as an iff on the committed lattice, the
committed matter table satisfying it on every row, the blocked
control label, and the exact one-sixth count.  The receipt is the
finite character-lattice half of the register row PR-46 bundle; the
complete character category, instanton normalization, laboratory flux
attachment, and physical quotient selection are not supplied, and no
physical global form is derived. -/
theorem globalFormCharacterDescent_receipt :
    (∀ w : W, Descends w ↔ centerChar w = 0)
      ∧ (∀ i : Fin 10, centerChar (componentWeight i) = 0)
      ∧ (∀ w ∈ OPH.Z6Descent.realizedWeights, centerChar w = 0)
      ∧ ¬ Descends ((1, 0, 0) : W)
      ∧ 6 * (Finset.univ.filter fun w : W => centerChar w = 0).card
          = Fintype.card W :=
  ⟨descent_criterion, exterior_rows_descend, realized_weights_descend,
    control_label_not_descends, descent_count.2.2⟩

end OPH.GlobalFormCharacterDescent

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.GlobalFormCharacterDescent.centerChar_eq_phase_gen
#print axioms OPH.GlobalFormCharacterDescent.centerChar_add
#print axioms OPH.GlobalFormCharacterDescent.centerCharHom
#print axioms OPH.GlobalFormCharacterDescent.phase_genPow
#print axioms OPH.GlobalFormCharacterDescent.descent_criterion
#print axioms OPH.GlobalFormCharacterDescent.descent_integer_congruence
#print axioms OPH.GlobalFormCharacterDescent.exterior_rows_descend
#print axioms OPH.GlobalFormCharacterDescent.exterior_rows_integer_congruence
#print axioms OPH.GlobalFormCharacterDescent.realized_weights_descend
#print axioms OPH.GlobalFormCharacterDescent.committed_table_descends
#print axioms OPH.GlobalFormCharacterDescent.control_label_blocked
#print axioms OPH.GlobalFormCharacterDescent.control_label_not_descends
#print axioms OPH.GlobalFormCharacterDescent.descent_count
#print axioms OPH.GlobalFormCharacterDescent.centerChar_fibres_equal
#print axioms OPH.GlobalFormCharacterDescent.masked_selection_descends
#print axioms OPH.GlobalFormCharacterDescent.globalFormCharacterDescent_receipt
