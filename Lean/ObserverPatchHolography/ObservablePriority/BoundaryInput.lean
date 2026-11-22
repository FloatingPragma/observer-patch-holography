import ObserverPatchHolography.ObservablePriority.SharpOrder

namespace PatchBoundaryInput
open ObserverPatchSharp

/-- Local read authority consists of its collars and the endpoint loads they may inspect. -/
def Same (menu : List Edge) (s t : Material) : Prop :=
  (∀ e ∈ menu, s.collar e = t.collar e) ∧
  (∀ e ∈ menu, s.load e.lo = t.load e.lo ∧ s.load e.hi = t.load e.hi)

theorem service_local {menu : List Edge} {s t : Material} (same : Same menu s t)
    {e : Edge} (he : e ∈ menu) : Same menu (service e s) (service e t) := by
  obtain ⟨hc,hl⟩ := same
  obtain ⟨hlo,hhi⟩ := hl e he
  have hr : ready s e = ready t e := by simp only [ready,initialReady,hlo,hhi,hc e he]
  have hs : sign s.load e = sign t.load e := by simp only [sign,hlo,hhi]
  unfold service
  rw [hr]
  split
  · constructor
    · intro f hf
      by_cases eq : f = e <;> simp [write,eq,hs,hc f hf]
    · intro f hf
      obtain ⟨la,lb⟩ := hl f hf
      simp only [write,la,lb,and_self]
  · exact ⟨hc,hl⟩

theorem window_local {menu run : List Edge} (inside : run ⊆ menu)
    {s t : Material} (same : Same menu s t) :
    Same menu (run.foldr service s) (run.foldr service t) := by
  induction run with
  | nil => exact same
  | cons a tail ih =>
    exact service_local (ih (fun _ hm => inside (List.mem_cons_of_mem a hm)))
      (inside (List.mem_cons_self))

def informed (L : Source) (boundaryWins : List Edge) : Source :=
  fun v => if touched boundaryWins v then .one else L v

/-- Once the higher external decisions are received, the local observer needs
    only its own initial endpoint loads. No outside sign or final state is read. -/
theorem boundary_flags_suffice (rank : Edge → Nat) (inner outer : List Edge)
    (ho : Ordered rank (inner ++ outer))
    (separate : ∀ e ∈ inner, e ∉ outer) (L : Source) :
    Same inner (sharpRun rank (inner ++ outer) L)
      (sharpRun rank inner (informed L (accepted (initialReady L) meets outer))) := by
  let wins := accepted (initialReady L) meets outer
  have hi : Ordered rank inner := (List.pairwise_append.mp ho).1
  have same : Same inner (material L wins) (initial (informed L wins)) := by
    constructor
    · intro e he
      have absent : e ∉ wins := fun h => separate e he (accepted_mem h)
      simp [material,initial,absent]
    · intro e he
      exact ⟨rfl,rfl⟩
  rw [sharpRun_eq_material _ _ _ ho, ← serviceFold_eq_material L (inner ++ outer), List.foldr_append,
      serviceFold_eq_material L outer,sharpRun_eq_material _ _ _ hi,← serviceFold_eq_material (informed L wins) inner]
  exact window_local (fun _ h => h) same

#print axioms window_local
#print axioms boundary_flags_suffice
end PatchBoundaryInput
