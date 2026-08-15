import SMStructureAdequacySurface
import BaryonDimensionSix

/-!
# The matter-grammar index dictionary (V3, issue #734)

Two committed `Fin 10` tables describe one-generation matter content
in different index orders: the exterior component table of
`Screen/ExteriorSelection.lean` (register row PR-59) and the frozen
Weyl field grammar of `Screen/BaryonDimensionSix.lean` in the order
`Q, L, u_R, d_R, e_R` followed by conjugates (register rows PR-55
and PR-56).  This module commits the explicit permutation between
the two index sets and kernel-checks column-by-column agreement:
charge against `q6`, the color and doublet indicators against
triality and the weak indicator, charge conjugation against
Hermitian conjugation, and even exterior parity against the undotted
left-Weyl rows.  A threaded corollary states the Weyl-parity reading
on the selection mask of the `SMStructurePremiseData` bundle.

The dictionary identifies index tables, not physics.  The census
columns `b3` and `ell` are declared labels under register row PR-55;
the dictionary transports declared labels and derives no baryon or
lepton number for the exterior components.  The dimension-six census
of `BaryonDimensionSix` stays a conditional finite operator census
under rows PR-55/PR-56 with rows PR-47, PR-54, and PR-57 open; it
proves neither proton stability nor proton decay, and this module
changes none of that.
-/

namespace OPH.MatterGrammarIndexBridge

/-- The committed index dictionary from the exterior component order
of `ExteriorSelection` (register row PR-59) to the frozen census
field order `Q, L, u_R, d_R, e_R, Q_bar, L_bar, u_R_bar, d_R_bar,
e_R_bar` of `BaryonDimensionSix` (register rows PR-55/PR-56). -/
def fieldIndex : Fin 10 → Fin 10 := ![3, 6, 7, 0, 9, 2, 5, 4, 8, 1]

/-- The dictionary is a bijection of the ten indices. -/
theorem fieldIndex_bijective : Function.Bijective fieldIndex := by
  decide

/-- The frozen integer charge column of the component table equals
the census hypercharge column `q6 = 6Y` through the dictionary. -/
theorem charge_column :
    ∀ i : Fin 10,
      OPH.ExteriorSelection.charge i
        = OPH.BaryonDimensionSix.q6 (fieldIndex i) := by decide

/-- The color indicator of the component table matches nonzero census
triality through the dictionary. -/
theorem color_column :
    ∀ i : Fin 10,
      OPH.ExteriorSelection.isColored i
        = (OPH.BaryonDimensionSix.triality (fieldIndex i) != 0) := by
  decide

/-- The weak-doublet indicator of the component table matches the
census weak indicator through the dictionary. -/
theorem doublet_column :
    ∀ i : Fin 10,
      OPH.ExteriorSelection.isDoublet i
        = (OPH.BaryonDimensionSix.weak (fieldIndex i) == 1) := by
  decide

/-- The dictionary intertwines exterior charge conjugation with the
census Hermitian conjugation. -/
theorem conj_column :
    ∀ i : Fin 10,
      fieldIndex (OPH.ExteriorSelection.conj i)
        = OPH.BaryonDimensionSix.conjugate (fieldIndex i) := by decide

/-- Even exterior parity lands exactly on the undotted left-Weyl
census rows through the dictionary. -/
theorem weyl_parity_column :
    ∀ i : Fin 10,
      OPH.BaryonDimensionSix.left (fieldIndex i) = 1
        ↔ OPH.ExteriorSelection.odd i = false := by decide

/-- **The matter-grammar index dictionary receipt (OL-G8
evidence).**  The committed permutation is bijective and matches the
charge, color, doublet, conjugation, and Weyl-parity columns of the
two frozen tables.  The identification is a kernel-checked dictionary
between index tables under register rows PR-59 and PR-55/PR-56; no
physical attachment is derived. -/
theorem matterGrammarIndexBridge_receipt :
    Function.Bijective fieldIndex
      ∧ (∀ i : Fin 10,
          OPH.ExteriorSelection.charge i
            = OPH.BaryonDimensionSix.q6 (fieldIndex i))
      ∧ (∀ i : Fin 10,
          OPH.ExteriorSelection.isColored i
            = (OPH.BaryonDimensionSix.triality (fieldIndex i) != 0))
      ∧ (∀ i : Fin 10,
          OPH.ExteriorSelection.isDoublet i
            = (OPH.BaryonDimensionSix.weak (fieldIndex i) == 1))
      ∧ (∀ i : Fin 10,
          fieldIndex (OPH.ExteriorSelection.conj i)
            = OPH.BaryonDimensionSix.conjugate (fieldIndex i))
      ∧ (∀ i : Fin 10,
          OPH.BaryonDimensionSix.left (fieldIndex i) = 1
            ↔ OPH.ExteriorSelection.odd i = false) :=
  ⟨fieldIndex_bijective, charge_column, color_column, doublet_column,
    conj_column, weyl_parity_column⟩

/-- **The threaded Weyl-parity corollary.**  For every premise bundle
`D` of the composed surface, every component selected by the bundled
mask maps to a census field whose Weyl-chirality column equals the
sector parity: the image row is undotted left-Weyl exactly when the
bundled selection is the even-parity sector.  Table provenance under
register row PR-59; a physical chirality claim consumes the
spacetime Spin attachment of register row PR-47, which no theorem
here supplies. -/
theorem selected_components_weyl_parity
    (D : OPH.SMStructureAdequacySurface.SMStructurePremiseData) :
    ∀ i : Fin 10,
      OPH.ExteriorSelection.mem D.selectionMask.val i = true →
        (OPH.BaryonDimensionSix.left (fieldIndex i) = 1
          ↔ D.selectionMask.val = OPH.ExteriorSelection.evenMask) := by
  rcases D.matter_selection_is_parity_sector with h | h <;> rw [h]
  · decide
  · decide

end OPH.MatterGrammarIndexBridge

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.MatterGrammarIndexBridge.fieldIndex_bijective
#print axioms OPH.MatterGrammarIndexBridge.charge_column
#print axioms OPH.MatterGrammarIndexBridge.color_column
#print axioms OPH.MatterGrammarIndexBridge.doublet_column
#print axioms OPH.MatterGrammarIndexBridge.conj_column
#print axioms OPH.MatterGrammarIndexBridge.weyl_parity_column
#print axioms OPH.MatterGrammarIndexBridge.matterGrammarIndexBridge_receipt
#print axioms OPH.MatterGrammarIndexBridge.selected_components_weyl_parity
