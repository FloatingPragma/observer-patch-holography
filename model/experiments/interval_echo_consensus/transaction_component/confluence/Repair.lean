import ObservableNormalForms.ObserverConfluence

namespace OPH.MeasuredPairConfluence
open ObservableNormalForms Relation

abbrev Boundary := Bool × Bool
abbrev Address := Bool × Fin 3
abbrev Offer := Bool × Fin 4
/-- Given the two protected bits, the unresolved addresses determine every active data bit. -/
abbrev State := Boundary × Finset Address

/-- A proper pair omits its index; index 3 requests the complete triangle. -/
def mask (o : Offer) : Finset Address :=
  Finset.univ.filter fun a => a.1 = o.1 ∧ (o.2.val = 3 ∨ a.2.val ≠ o.2.val)

def repair (o : Offer) (x : State) : State := (x.1, x.2 \ mask o)

def step (x y : State) : Prop := ∃ o, y = repair o x ∧ ∃ a ∈ x.2, a ∈ mask o

def consistent : Set State := {x | x.2 = ∅}

theorem closes (o : Offer) (x : State) {a : Address} (h : a ∈ mask o) :
    a ∉ (repair o x).2 := by simp [repair, h]

theorem descent {x y : State} (h : step x y) : y.2.card < x.2.card := by
  rcases h with ⟨o, rfl, a, ha, hm⟩
  exact Finset.card_lt_card ((Finset.ssubset_iff_of_subset Finset.sdiff_subset).2
    ⟨a, ha, closes o x hm⟩)

theorem terminating : WellFounded (fun y x => step x y) :=
  Subrelation.wf (fun {_ _} h => descent h) (measure (fun x : State => x.2.card)).wf

theorem offer_path (x : State) (o : Offer) : ReflTransGen step x (repair o x) := by
  by_cases h : ∃ a ∈ x.2, a ∈ mask o
  · exact .single ⟨o, rfl, h⟩
  · have eq : x.2 \ mask o = x.2 := by
      ext a
      simp only [Finset.mem_sdiff]
      exact ⟨And.left, fun ha => ⟨ha, fun hm => h ⟨a, ha, hm⟩⟩⟩
    have same : repair o x = x := by simp [repair, eq]
    rw [same]

theorem repair_comm (x : State) (o p : Offer) :
    repair p (repair o x) = repair o (repair p x) := by
  apply Prod.ext
  · rfl
  · ext a
    simp only [repair, Finset.mem_sdiff]
    tauto

theorem diamond {x y z : State} (hy : step x y) (hz : step x z) :
    ∃ w, ReflTransGen step y w ∧ ReflTransGen step z w := by
  rcases hy with ⟨o, rfl, _⟩
  rcases hz with ⟨p, rfl, _⟩
  exact ⟨_, offer_path (repair o x) p,
    (repair_comm x o p).symm ▸ offer_path (repair p x) o⟩

theorem complete : CompleteFor step consistent := by
  intro x
  constructor
  · intro hn
    change x.2 = ∅
    apply Finset.eq_empty_iff_forall_notMem.mpr
    intro a ha
    let o : Offer := (a.1, 3)
    have hm : a ∈ mask o := by simp [mask, o]
    exact hn _ ⟨o, rfl, a, ha, hm⟩
  · rintro hx y ⟨o, rfl, a, ha, _⟩
    have he : x.2 = ∅ := hx
    simp [he] at ha

theorem normalizes (x : State) : WeaklyNormalizing step x := by
  let a : Offer := (false, 3)
  let b : Offer := (true, 3)
  refine ⟨_, (offer_path x a).trans (offer_path (repair a x) b), (complete _).2 ?_⟩
  change (x.2 \ mask a) \ mask b = ∅
  apply Finset.eq_empty_iff_forall_notMem.mpr
  rintro ⟨side, i⟩
  cases side <;> simp [mask, a, b]

theorem preserves : ObservationPreserving step Prod.fst := by
  rintro _ _ ⟨o, rfl, _⟩
  rfl

theorem separation_iff {T : Type} (read : Boundary → T) :
    BoundaryIdentifiesModulo consistent (fun x => read x.1) Eq ↔
      ObserverEndpointUniqueModulo step (fun x => read x.1) Eq :=
  boundaryIdentifiesModulo_iff_observerEndpointUniqueModulo
    (fun {_ _} h => congrArg read (preserves h)) complete

theorem full_read_unique : ObserverEndpointUniqueModulo step Prod.fst Eq := by
  apply (separation_iff id).mp
  intro x y hx hy hxy
  exact Prod.ext hxy ((show x.2 = ∅ from hx).trans (show y.2 = ∅ from hy).symm)

def parity (b : Boundary) : Bool := xor b.1 b.2

theorem parity_read_fails :
    ¬ ObserverEndpointUniqueModulo step (fun x => parity x.1) Eq := by
  intro h
  have eq := (separation_iff parity).mpr h
      (c := ((false,false),∅)) (d := ((true,true),∅)) rfl rfl rfl
  cases eq

#print axioms terminating
#print axioms diamond
#print axioms full_read_unique
#print axioms parity_read_fails
end OPH.MeasuredPairConfluence
