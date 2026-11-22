import ObserverPatchHolography.ObservablePriority.SharpOrder
import ObserverPatchHolography.ObservablePriority.PatchClosure
import Mathlib.CategoryTheory.Category.Preorder
import Mathlib.CategoryTheory.Types.Basic

namespace PatchComposition
open ObserverPatchSharp CategoryTheory Opposite

def parents (rank : Edge → Nat) (screen : List Edge) (e : Edge) : List Edge :=
  screen.filter fun f => meets e f && decide (rank f < rank e)

theorem physical_member {rank : Edge → Nat} {screen : List Edge} {e f : Edge} :
    f ∈ parents rank screen e ↔ f ∈ screen ∧ meets e f = true ∧ rank f < rank e := by
  simp [parents]

theorem certificate_closed (rank : Edge → Nat) (screen seeds cert : List Edge)
    (ok : PatchClosure.acceptCertificate seeds (parents rank screen) cert = true) :
    PredecessorClosed rank screen cert.contains := by
  intro a _ ha b hb hm hp
  have h := (PatchClosure.acceptCertificate_sound ok).1.present a (by simpa using ha)
    b (physical_member.mpr ⟨hb,hm,hp⟩)
  simpa using h

/-- A genuine downstream consumer of both kernel providers. -/
theorem certificate_predicts_signed_records (rank : Edge → Nat) (screen seeds cert : List Edge)
    (legal : Ordered rank screen)
    (ok : PatchClosure.acceptCertificate seeds (parents rank screen) cert = true) :
    AllSchedulesAgree rank screen cert.contains :=
  (predecessorClosed_iff_allSchedulesAgree rank screen cert.contains legal).mp
    (certificate_closed rank screen seeds cert ok)

abbrev Context (rank : Edge → Nat) (screen : List Edge) :=
  {s : Set Edge // (∀ a ∈ s, a ∈ screen) ∧
    (∀ a ∈ s, ∀ b ∈ parents rank screen a, b ∈ s)}

def records (rank : Edge → Nat) (screen : List Edge) :
    (Context rank screen)ᵒᵖ ⥤ Type where
  obj U := {e // e ∈ U.unop.val} → Option Int
  map f := PatchClosure.restrict (fun _ h => leOfHom f.unop h)
  map_id _ := by
    funext r
    exact PatchClosure.restrict_id r
  map_comp f g := by
    funext r
    exact (PatchClosure.restrict_comp (fun _ h => leOfHom g.unop h)
      (fun _ h => leOfHom f.unop h) r).symm

noncomputable def see (rank : Edge → Nat) (screen : List Edge) (P : Context rank screen)
    (L : Source) : (records rank screen).obj (op P) := by
  classical
  exact fun e => (sharpRun rank (screen.filter (fun x => decide (x ∈ P.val))) L).collar e.val

private theorem see_full {rank : Edge → Nat} {screen : List Edge}
    (ho : Ordered rank screen) (P : Context rank screen) (L : Source) (e : {e // e ∈ P.val}) :
    see rank screen P L e = (sharpRun rank screen L).collar e.val := by
  classical
  have hc : PredecessorClosed rank screen (fun x => decide (x ∈ P.val)) := by
    intro a _ ha b hb hm hp
    exact decide_eq_true (P.property.2 a (of_decide_eq_true ha) b (physical_member.mpr ⟨hb,hm,hp⟩))
  exact (closed_agrees rank screen _ ho hc L e.val (P.property.1 e.val e.property)
    (decide_eq_true e.property)).symm

theorem native_naturality {rank : Edge → Nat} {screen : List Edge}
    (ho : Ordered rank screen) {P Q : Context rank screen} (i : P ⟶ Q) (L : Source) :
    (records rank screen).map i.op (see rank screen Q L) = see rank screen P L := by
  funext a
  change see rank screen Q L ⟨a.val, leOfHom i a.property⟩ = see rank screen P L a
  rw [see_full ho,see_full ho]

#print axioms certificate_predicts_signed_records
#print axioms native_naturality
end PatchComposition
