import ObserverPatchHolography.ObservablePriority.ObservablePriority
import Mathlib.Data.Finset.Max
import Mathlib.Data.List.Pairwise

/-! Constructing an executable service order from physical wedge observations.
    The search and screening code has no input of type Candidate or rank. -/
namespace ObservablePriority
open ObserverPatchSharp

/-- A complete deterministic wedge observation must name exactly one winner
    among the two incident edges. Equal bits cannot define a direction. -/
def directionsCovered (screen : List Edge) (o : HitOracle) : Bool :=
  screen.all fun a => screen.all fun b =>
    !eligible a b || (o (wedge a b) a != o (wedge a b) b)

private def available (parents : Edge → List Edge) (pending : List Edge) (a : Edge) : Bool :=
  (parents a).all fun b => !pending.contains b

/-- Kahn's algorithm. Menu order breaks ties between available vertices.
    Each removal consumes one unit of the original menu length. -/
private def peel (parents : Edge → List Edge) : Nat → List Edge → Option (List Edge)
  | _, [] => some []
  | 0, _ :: _ => none
  | fuel + 1, a :: tail => do
    let p := a :: tail
    let first ← p.find? (available parents p)
    let rest ← peel parents fuel (p.erase first)
    return first :: rest

def blindOrder (screen : List Edge) (o : HitOracle) : Option (List Edge) :=
  if screen.Nodup ∧ directionsCovered screen o then
    peel (recoverParents screen o) screen.length screen
  else none

private theorem apart {screen : List Edge} (c : Candidate screen) {a b : Edge}
    (ha : a ∈ screen) (hb : b ∈ screen) (ne : a ≠ b) (hm : meets a b = true) :
    c.rank a ≠ c.rank b := by
  have rel : c.schedule.Pairwise (fun x y => meets x y = true → c.rank x ≠ c.rank y) :=
    c.legal.imp (fun h hit => (Nat.ne_of_lt (h.2 hit)).symm)
  have sym : Symmetric (fun x y => meets x y = true → c.rank x ≠ c.rank y) := by
    intro x y h hit
    exact (h ((meets_symm ..).trans hit)).symm
  exact List.Pairwise.forall sym rel
    (c.complete.mem_iff.mpr ha) (c.complete.mem_iff.mpr hb) ne hm

private theorem service_covered {screen : List Edge} {c : Candidate screen} {o : HitOracle}
    (service : Realizes c o) : directionsCovered screen o = true := by
  apply List.all_eq_true.mpr
  intro a ha
  apply List.all_eq_true.mpr
  intro b hb
  by_cases h : eligible a b = true
  · obtain ⟨ne, hm⟩ : a ≠ b ∧ meets a b = true := by simpa [eligible] using h
    have lhs := oracle_wedge_hit service hb ha ne.symm ((meets_symm ..).trans hm)
    rw [wedge_symm b a] at lhs
    have rhs := oracle_wedge_hit service ha hb ne hm
    have n := apart c ha hb ne hm
    rcases Nat.lt_or_gt_of_ne n with lt | gt
    · simp [h, lhs, rhs, lt, Nat.not_lt_of_ge (Nat.le_of_lt lt)]
    · simp [h, lhs, rhs, gt, Nat.not_lt_of_ge (Nat.le_of_lt gt)]
  · simp [h]

/-- A returned order contains exactly the pending vertices and no later
    vertex is a parent of an earlier vertex, even for an arbitrary oracle. -/
private theorem peel_sound {parents : Edge → List Edge} {fuel : Nat} {xs order : List Edge}
    (ok : peel parents fuel xs = some order) :
    order.Perm xs ∧ order.Pairwise (fun a b => b ∉ parents a) := by
  induction fuel generalizing xs order with
  | zero =>
    cases xs with
    | nil => cases ok; exact ⟨.rfl, .nil⟩
    | cons a rest => simp [peel] at ok
  | succ fuel ih =>
    cases xs with
    | nil => cases ok; exact ⟨.rfl, .nil⟩
    | cons x rest =>
      let p := x :: rest
      change ((p.find? (available parents p)).bind (fun first =>
        (peel parents fuel (p.erase first)).bind (fun tail => some (first :: tail)))) = some order at ok
      cases choice : p.find? (available parents p) with
      | none => rw [choice] at ok; contradiction
      | some a =>
        have ha := List.mem_of_find?_eq_some choice
        have ready : available parents p a = true := List.find?_some choice
        cases result : peel parents fuel (p.erase a) with
        | none => rw [choice, Option.bind_some, result] at ok; contradiction
        | some tail =>
          have eq : a :: tail = order := by
            rw [choice, Option.bind_some, result] at ok
            exact Option.some.inj ok
          subst order
          obtain ⟨pt, ht⟩ := ih result
          refine ⟨(pt.cons a).trans (List.perm_cons_erase ha).symm,
            List.pairwise_cons.mpr ⟨?_, ht⟩⟩
          intro b hb parent
          have notIn : b ∉ p := by simpa using List.all_eq_true.mp ready b parent
          exact notIn (List.mem_of_mem_erase (pt.mem_iff.mp hb))

private theorem peel_complete {screen : List Edge} {c : Candidate screen} {o : HitOracle}
    (service : Realizes c o) (fuel : Nat) (p : List Edge) (inside : p ⊆ screen)
    (lengthOK : p.length ≤ fuel) :
    ∃ order, peel (recoverParents screen o) fuel p = some order := by
  induction fuel generalizing p with
  | zero =>
    have hp : p = [] := List.length_eq_zero_iff.mp (by omega)
    subst p
    exact ⟨[], rfl⟩
  | succ fuel ih =>
    cases p with
    | nil => exact ⟨[], rfl⟩
    | cons x rest =>
      let p := x :: rest
      let parents := recoverParents screen o
      have present : p.toFinset.Nonempty := ⟨x, by simp [p]⟩
      obtain ⟨m, hm, minimal⟩ := Finset.exists_min_image p.toFinset c.rank present
      have mem : m ∈ p := List.mem_toFinset.mp hm
      have ready : available parents p m = true := by
        apply List.all_eq_true.mpr
        intro b hb
        have lt := ((recoverParents_member service (inside mem)).mp hb).2.2
        have absent : b ∉ p := fun inP =>
          Nat.not_lt_of_ge (minimal b (List.mem_toFinset.mpr inP)) lt
        simpa using absent
      cases choice : p.find? (available parents p) with
      | none => exact False.elim ((List.find?_eq_none.mp choice m mem) ready)
      | some a =>
        have ha : a ∈ p := List.mem_of_find?_eq_some choice
        have short : (p.erase a).length ≤ fuel := by
          rw [List.length_erase_of_mem ha]
          change (x :: rest).length - 1 ≤ fuel
          omega
        obtain ⟨tail, ht⟩ := ih (p.erase a) (fun _ hb => inside (List.mem_of_mem_erase hb)) short
        refine ⟨a :: tail, ?_⟩
        change ((p.find? (available parents p)).bind (fun first =>
          (peel parents fuel (p.erase first)).bind (fun tail => some (first :: tail)))) = some (a :: tail)
        rw [choice, Option.bind_some, ht]
        rfl

theorem blindOrder_installation {screen order : List Edge} {o : HitOracle}
    (ok : blindOrder screen o = some order) :
    order.Perm screen ∧ PatchClosure.Installation (recoverParents screen o) order := by
  unfold blindOrder at ok
  split at ok
  next valid =>
    obtain ⟨perm, earlier⟩ := peel_sound ok
    refine ⟨perm, ⟨perm.nodup_iff.mpr valid.1, ?_, ?_, earlier⟩⟩
    · intro a _ self
      have head := (Bool.and_eq_true_iff.mp (List.mem_filter.mp self).2).1
      simp [eligible] at head
    · intro a _ b hb
      exact perm.mem_iff.mpr (List.mem_filter.mp hb).1
  next => contradiction

theorem blindOrder_some_legal {screen order : List Edge} {c : Candidate screen} {o : HitOracle}
    (service : Realizes c o) (ok : blindOrder screen o = some order) :
    order.Perm screen ∧ Ordered c.rank order.reverse := by
  obtain ⟨perm, install⟩ := blindOrder_installation ok
  refine ⟨perm, List.pairwise_reverse.mpr ?_⟩
  apply (install.unique.and install.earlier).imp_of_mem
  intro a b ha hb h
  refine ⟨h.1.symm, fun hm => ?_⟩
  have shared := (meets_symm ..).trans hm
  have ne := apart c (perm.mem_iff.mp ha) (perm.mem_iff.mp hb) h.1 shared
  have noEarlier : ¬c.rank b < c.rank a := fun lt => h.2
    ((recoverParents_member service (perm.mem_iff.mp ha)).mpr ⟨perm.mem_iff.mp hb, shared, lt⟩)
  omega

theorem blindOrder_complete {screen : List Edge} {c : Candidate screen} {o : HitOracle}
    (service : Realizes c o) :
    ∃ order, blindOrder screen o = some order ∧ order.Perm screen ∧ Ordered c.rank order.reverse := by
  have unique : screen.Nodup := c.complete.nodup (c.legal.imp fun h => h.1)
  obtain ⟨order, ok⟩ := peel_complete service screen.length screen (fun _ h => h) (Nat.le_refl _)
  have found : blindOrder screen o = some order := by
    simpa [blindOrder, unique, service_covered service] using ok
  exact ⟨order, found, blindOrder_some_legal service found⟩

theorem blindOrder_closure {screen seeds : List Edge} {c : Candidate screen} {o : HitOracle}
    (service : Realizes c o) (inside : seeds ⊆ screen) :
    ∃ order, blindOrder screen o = some order ∧
      checkedClosure screen o seeds order = some (PatchClosure.generate seeds (recoverParents screen o) order) ∧
      ∀ a, a ∈ PatchClosure.generate seeds (recoverParents screen o) order ↔
        PatchClosure.Needed seeds (PatchComposition.parents c.rank screen) a := by
  obtain ⟨order, found, perm, legal⟩ := blindOrder_complete service
  have accepted := legal_closure_accepted service (order.reverse_perm.trans perm) legal inside
  rw [List.reverse_reverse] at accepted
  exact ⟨order, found, accepted, (recovered_closure_exact service accepted).2⟩

theorem blindOrder_reject_unoriented {screen : List Edge} {o : HitOracle} {a b : Edge}
    (ha : a ∈ screen) (hb : b ∈ screen) (canTest : eligible a b = true)
    (sameBit : o (wedge a b) a = o (wedge a b) b) : blindOrder screen o = none := by
  have no : ¬ directionsCovered screen o = true := by
    intro full
    have := List.all_eq_true.mp (List.all_eq_true.mp full a ha) b hb
    simp [canTest, sameBit] at this
  simp [blindOrder, no]

private theorem no_internal_cycle {parents : Edge → List Edge} {cycle : List Edge} {start : Edge}
    (active : start ∈ cycle) (loop : ∀ a ∈ cycle, ∃ b ∈ cycle, b ∈ parents a) (xs : List Edge)
    (inside : cycle ⊆ xs) (noSelf : ∀ a ∈ xs, a ∉ parents a)
    (forward : xs.Pairwise (fun a b => b ∉ parents a)) : False := by
  induction xs with
  | nil => exact List.not_mem_nil (inside active)
  | cons a tail ih =>
    obtain ⟨head, remaining⟩ := List.pairwise_cons.mp forward
    by_cases hit : a ∈ cycle
    · obtain ⟨b, hb, parent⟩ := loop a hit
      rcases List.mem_cons.mp (inside hb) with rfl | later
      · exact noSelf _ (by simp) parent
      · exact head _ later parent
    · exact ih (fun b hb => (List.mem_cons.mp (inside hb)).resolve_left (fun eq => hit (eq ▸ hb)))
        (fun b hb => noSelf b (List.mem_cons_of_mem a hb)) remaining

/-- A nonempty closed set of stuck vertices (in particular a directed cycle)
    forces rejection even if each local pair reports a unique winner. -/
theorem blindOrder_reject_cycle {screen cycle : List Edge} {o : HitOracle} {a : Edge}
    (start : a ∈ cycle) (inside : cycle ⊆ screen)
    (loop : ∀ x ∈ cycle, ∃ p ∈ cycle, p ∈ recoverParents screen o x) :
    blindOrder screen o = none := by
  cases result : blindOrder screen o with
  | none => rfl
  | some xs =>
    obtain ⟨perm, inst⟩ := blindOrder_installation result
    exact False.elim (no_internal_cycle start loop xs
      (fun _ h => perm.mem_iff.mpr (inside h)) inst.noSelf inst.earlier)

#print axioms blindOrder_installation
#print axioms blindOrder_complete
#print axioms blindOrder_closure
#print axioms blindOrder_reject_unoriented
#print axioms blindOrder_reject_cycle
end ObservablePriority
