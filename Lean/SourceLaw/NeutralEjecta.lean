import Mathlib

/-! Inventory needed for disjoint neutral clicks on an initially filled, blank
Interval Echo receiver. The native seam certificate establishes conservation
of double loads and blank collars during free flight. These two lemmas apply
to the resulting terminal observer words. -/

namespace SourceLaw.NeutralEjecta

/-- In a sharp terminal word, zero net signed load equates the numbers of holes
and doubles. Any load departure then gives a double in that receiver. -/
theorem double_of_neutral_click {V : Type*} (holes doubles : Finset V)
    (neutral : holes.card = doubles.card) (click : holes.Nonempty ∨ doubles.Nonempty) :
    doubles.Nonempty := by
  rcases click with h | h
  · exact Finset.card_pos.mp (neutral ▸ Finset.card_pos.mpr h)
  · exact h

/-- Every disjoint neutral click needs its own double in addition to the
doubles retained by the charged daughter. Applies to any finite number of
receivers, any unobserved region, and any shape of the detector patches. -/
theorem disjoint_ejecta_inventory {V R : Type*} [DecidableEq V] [Fintype R]
    (doubles daughter : Finset V) (receiver : R → Finset V)
    (separate : ∀ i j, i ≠ j → Disjoint (receiver i) (receiver j))
    (away : ∀ i, Disjoint daughter (receiver i))
    (click : ∀ i, ((receiver i) ∩ doubles).Nonempty) :
    (daughter ∩ doubles).card + Fintype.card R ≤ doubles.card := by
  classical
  choose hit mem using click
  have inR (i : R) : hit i ∈ receiver i := (Finset.mem_inter.mp (mem i)).1
  have inD (i : R) : hit i ∈ doubles := (Finset.mem_inter.mp (mem i)).2
  let f : (↥(daughter ∩ doubles) ⊕ R) → ↥doubles :=
    Sum.elim (fun d => ⟨d.val, (Finset.mem_inter.mp d.property).2⟩)
      (fun i => ⟨hit i, inD i⟩)
  have inj : Function.Injective f := by
    intro x y eq
    have same := congrArg Subtype.val eq
    cases x with
    | inl x =>
      cases y with
      | inl y => exact congrArg Sum.inl (Subtype.ext same)
      | inr j =>
        change x.val = hit j at same
        exact False.elim (Finset.disjoint_left.mp (away j)
          ((Finset.mem_inter.mp x.property).1) (same ▸ inR j))
    | inr i =>
      cases y with
      | inl y =>
        change hit i = y.val at same
        exact False.elim (Finset.disjoint_left.mp (away i)
          ((Finset.mem_inter.mp y.property).1) (same ▸ inR i))
      | inr j =>
        change hit i = hit j at same
        apply congrArg Sum.inr
        by_contra ne
        exact Finset.disjoint_left.mp (separate i j ne) (inR i) (same ▸ inR j)
  simpa using Fintype.card_le_of_injective f inj

end SourceLaw.NeutralEjecta
