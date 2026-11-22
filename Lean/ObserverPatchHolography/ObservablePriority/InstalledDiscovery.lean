import ObserverPatchHolography.ObservablePriority.BlindOrder
import ObserverPatchHolography.ObservablePriority.ObservableCell
import ObserverPatchHolography.ObservablePriority.NetworkRows

/-! Recover a native observational class before consuming its minimal cell menu
    and exact communication theorem. If the observed classes differ, return the
    common preparation and address that distinguishes the two installations. -/
namespace InstalledDiscovery
open ObservablePriority ObservableCell ObserverPatchSharp CellBits PriorityCellData InstalledNetwork

def ExactCost (c : CellChart) (target : Source → List Bool) : Prop :=
  Communicates c target 3 ∧ ∀ b, Communicates c target b → 3 ≤ b

def Packet (s : Federation) (k : Nat) (law : Candidate (board s)) (o : HitOracle) : Prop :=
  ∃ order cert, blindOrder (board s) o = some order ∧
    checkedClosure (board s) o (chart s k).inner order = some cert ∧ FitsMenu (chart s k) cert ∧
    ExactCost (chart s k) (fun L => localWord (chart s k) (sharpRun law.rank law.schedule L))

theorem native_packet (s : Federation) (k : Nat) (hk : k < size s)
    (o : HitOracle) (real : Realizes (candidate s) o) : Packet s k (candidate s) o := by
  let c := chart s k
  have onBoard : menu c ⊆ board s := fun e he =>
    PriorityNetwork.mem_inventory.mpr ⟨k,hk,he⟩
  have seeds : c.inner ⊆ board s := fun e he => onBoard (List.mem_append_left c.outer.toList he)
  obtain ⟨order,found,build,_⟩ := blindOrder_closure real seeds
  let cert := PatchClosure.generate c.inner (recoverParents (board s) o) order
  have built : checkedClosure (board s) o c.inner order = some cert := build
  obtain ⟨_,hm,_,_⟩ := full_window s k hk
  have nodup : (menu c).Nodup := hm.imp (fun h => h.1)
  have take := PriorityNetwork.restriction (covered_rows s) (PriorityNetwork.locallyOrdered (admitted_rows s))
    k hk (full := board s) (.refl _)
  have mem : ∀ e ∈ board s, e ∈ menu c ↔ PriorityNetwork.keep k e = true := by
    intro e he
    have same : e ∈ (board s).filter (PriorityNetwork.keep k) ↔ e ∈ menu c := take.mem_iff
    simpa only [List.mem_filter,he,true_and] using same.symm
  have cl : PredecessorClosed cellRank (board s) (menu c).contains := by
    intro e he incl f hf m earlier
    have kept := (mem e he).mp (by simpa using incl)
    have hf' := PriorityNetwork.closed k (board s) e he kept f hf m earlier
    simpa using (mem f hf).mpr hf'
  have fit := fits_of_minimum c real built onBoard nodup cl (minimal_window s k hk)
  have hc : cellCheck c = true := (Bool.and_eq_true_iff.mp (admitted_rows s k hk)).1
  exact ⟨order,cert,found,built,fit,recovered_cell_cost c real built fit hm hc⟩

private theorem transfer_realizes {world : List Edge} {c d : Candidate world} {o : HitOracle}
    (same : ObservationalAgreement c d) (h : Realizes d o) : Realizes c o := by
  intro L e he
  obtain ⟨ys,perm,legal,hv⟩ := h L e he
  exact ⟨c.schedule,c.complete,c.legal,hv.trans
    (same c.schedule ys c.complete perm c.legal legal L e he).symm⟩

/-- An unrecognized law in the admitted family is subject to a single
    outcome: a cell decoding packet or an actually distinguishing experiment. -/
theorem classify_or_separate (s : Federation) (k : Nat) (hk : k < size s)
    (alt : Candidate (board s)) (reference unknown : HitOracle)
    (hr : Realizes (candidate s) reference) (hu : Realizes alt unknown) :
    match separator (board s) reference unknown with
    | none => Packet s k alt unknown
    | some (e,f) => e ∈ board s ∧ f ∈ board s ∧ meets e f = true ∧
        reference (wedge e f) f ≠ unknown (wedge e f) f := by
  cases found : separator (board s) reference unknown with
  | none =>
    have same := (separator_none_iff hr hu).mp found
    obtain ⟨order,cert,bo,bc,fit,cost⟩ := native_packet s k hk unknown (transfer_realizes same hu)
    have words : (fun L => localWord (chart s k) (sharpRun (candidate s).rank (candidate s).schedule L)) =
        (fun L => localWord (chart s k) (sharpRun alt.rank alt.schedule L)) := by
      funext L
      apply List.map_congr_left
      intro e he
      have inside : e ∈ board s := PriorityNetwork.mem_inventory.mpr
        ⟨k,hk,List.mem_append_left _ he⟩
      exact same _ _ (candidate s).complete alt.complete (candidate s).legal alt.legal L e inside
    rw [words] at cost
    exact ⟨order,cert,bo,bc,fit,cost⟩
  | some pair =>
    obtain ⟨e,f⟩ := pair
    obtain ⟨he,hf,_,hM,test⟩ := separator_some_bits hr hu found
    obtain ⟨_,_,diff⟩ := test _ _ (candidate s).complete alt.complete (candidate s).legal alt.legal
    exact ⟨he,hf,hM,by rcases diff with h|h <;> simp [h.1,h.2]⟩

#print axioms native_packet
#print axioms classify_or_separate
end InstalledDiscovery
