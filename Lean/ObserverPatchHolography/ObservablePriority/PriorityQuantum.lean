import ObserverPatchHolography.ObservablePriority.ObservablePriority
import ObserverPatchHolography.ObservablePriority.FullChannel

/-! Identification of deterministic priorities for the source instrument on a
    bounded physical vertex set. Input and output coordinates are common to
    both candidates. Each may use its own complete legal service order. -/
namespace ObservablePriority.Quantum
open ObserverPatchSharp ObserverPatchSharp.QuantumBridge ObserverPatchSharp.FullChannel

private theorem board_perm {n : Nat} {xs ys : List (BEdge n)} (p : xs.Perm ys) :
    (board xs).Perm (board ys) :=
  (List.reverse_perm _).trans ((p.map Subtype.val).trans (List.reverse_perm _).symm)

private theorem transfer_legal {n : Nat} {menu ys : List (BEdge n)}
    {c d : Candidate (board menu)} (same : IncidentAgreement c d) (py : ys.Perm menu)
    (hy : Ordered d.rank (board ys)) : Ordered c.rank (board ys) := by
  have pb := board_perm py
  exact hy.imp_of_mem fun ha hb h => ⟨h.1, fun hm =>
    (same _ (pb.mem_iff.mp ha) _ (pb.mem_iff.mp hb) h.1 hm).mpr (h.2 hm)⟩

/-- Same addressed branch operators and the corresponding output density
    matrices at every physical coordinate, including off-diagonal entries. -/
theorem incident_full_channel {n : Nat} {menu xs ys : List (BEdge n)}
    {c d : Candidate (board menu)} (same : IncidentAgreement c d)
    (px : xs.Perm menu) (py : ys.Perm menu)
    (hx : Ordered c.rank (board xs)) (hy : Ordered d.rank (board ys)) :
    (∀ (M : BEdge n → Bool) (x y : StrictQuantum.Loads (Fin n))
      (coll : BEdge n → StrictQuantum.Addressed.Collar),
      StrictQuantum.Addressed.element (screen n) xs M x y coll =
      StrictQuantum.Addressed.element (screen n) ys M x y coll) ∧
    (∀ (M : BEdge n → Bool) (ρ : StrictQuantum.Observed.density (Fin n))
      (y z : StrictQuantum.Loads (Fin n)) (cy cz : BEdge n → StrictQuantum.Addressed.Collar),
      densityElement xs M ρ y z cy cz = densityElement ys M ρ y z cy cz) :=
  ⟨actual_element_schedule_agreement c.rank (px.trans py.symm) hx (transfer_legal same py hy),
   actual_density_schedule_agreement c.rank (px.trans py.symm) hx (transfer_legal same py hy)⟩

def MarginalAgreement {n : Nat} (menu xs ys : List (BEdge n)) : Prop :=
  ∀ (wanted : List (BEdge n)), wanted ⊆ menu →
    ∀ (obs : List Bool) (ρ : StrictQuantum.Observed.density (Fin n)),
      wordMass xs wanted obs ρ = wordMass ys wanted obs ρ

private theorem incident_marginals {n : Nat} {menu xs ys : List (BEdge n)}
    {c d : Candidate (board menu)} (same : IncidentAgreement c d)
    (px : xs.Perm menu) (py : ys.Perm menu)
    (hx : Ordered c.rank (board xs)) (hy : Ordered d.rank (board ys)) :
    MarginalAgreement menu xs ys := by
  have hb : (board ys).Perm ((board xs).filter (fun _ => true)) := by
    simpa using board_perm (py.trans px.symm)
  have hc : PredecessorClosed c.rank (board xs) (fun _ => true) := by
    intro _ _ _ _ _ _ _
    rfl
  intro wanted hw obs ρ
  exact closed_patch_marginal c.rank xs ys wanted (fun _ => true) hx (transfer_legal same py hy)
    hb hc (fun _ hj => ⟨px.mem_iff.mpr (hw hj), rfl⟩) obs ρ

private theorem wedge_mass_bit {n : Nat} (r : Edge → Nat) (xs : List (BEdge n))
    (ho : Ordered r (board xs)) (e f : BEdge n) (hf : f ∈ xs) :
    wordMass xs [f] [true] (basisDensity (wedgeInput e f)) =
      if hit r (board xs) (wedge e.val f.val) f.val then 1 else 0 := by
  simp only [basis_wordMass, flags_sharp r xs _ ho, readWord, List.map_cons, List.map_nil,
    bitAt_map _ _ _ hf, source_wedge, hit, List.cons.injEq, and_true]
  by_cases h : ((sharpRun r (board xs) (wedge e.val f.val)).collar f.val).isSome = true
  · simp only [if_pos h]
  · simp only [if_neg h]

/-- A successful classical search provides two bounded physical edges, one
    basis preparation, and the same addressed event under both candidates.
    The actual marginal masses of that event are complementary zero and one. -/
theorem separator_some_wordMass {n : Nat} {menu xs ys : List (BEdge n)}
    {c d : Candidate (board menu)} {left right : HitOracle} (hl : Realizes c left)
    (hr : Realizes d right) (px : xs.Perm menu) (py : ys.Perm menu)
    (hx : Ordered c.rank (board xs)) (hy : Ordered d.rank (board ys)) {a b : Edge}
    (found : separator (board menu) left right = some (a, b)) :
    ∃ e ∈ menu, ∃ f ∈ menu, e.val = a ∧ f.val = b ∧
      ((wordMass xs [f] [true] (basisDensity (wedgeInput e f)) = 0 ∧
        wordMass ys [f] [true] (basisDensity (wedgeInput e f)) = 1) ∨
       (wordMass xs [f] [true] (basisDensity (wedgeInput e f)) = 1 ∧
        wordMass ys [f] [true] (basisDensity (wedgeInput e f)) = 0)) := by
  obtain ⟨ha, hb, _, _, outputs⟩ := separator_some_bits hl hr found
  obtain ⟨e, he, rfl⟩ := List.mem_map.mp (show a ∈ menu.map Subtype.val by simpa [board] using ha)
  obtain ⟨f, hf, rfl⟩ := List.mem_map.mp (show b ∈ menu.map Subtype.val by simpa [board] using hb)
  obtain ⟨l, r, values⟩ :=
    outputs (board xs) (board ys) (board_perm px) (board_perm py) hx hy
  refine ⟨e, he, f, hf, rfl, rfl, ?_⟩
  rw [wedge_mass_bit c.rank xs hx e f (px.mem_iff.mpr hf),
      wedge_mass_bit d.rank ys hy e f (py.mem_iff.mpr hf), l, r]
  rcases values with ⟨lzero, rone⟩ | ⟨lone, rzero⟩
  · simp [lzero, rone]
  · simp [lone, rzero]

/-- For any two supplied legal schedules on the common bounded menu, no
    separator is returned exactly when every addressed outcome marginal
    agrees on every input density. In the converse only a pure basis
    preparation and a singleton requested address are needed. -/
theorem separator_none_iff_quantum_marginals {n : Nat} {menu xs ys : List (BEdge n)}
    {c d : Candidate (board menu)} {left right : HitOracle}
    (hl : Realizes c left) (hr : Realizes d right)
    (px : xs.Perm menu) (py : ys.Perm menu)
    (hx : Ordered c.rank (board xs)) (hy : Ordered d.rank (board ys)) :
    separator (board menu) left right = none ↔ MarginalAgreement menu xs ys := by
  constructor
  · intro no
    exact incident_marginals ((observational_iff_incident c d).mp ((separator_none_iff hl hr).mp no))
      px py hx hy
  · intro agree
    cases result : separator (board menu) left right with
    | none => rfl
    | some pair =>
      obtain ⟨a, b⟩ := pair
      obtain ⟨e, _, f, hf, _, _, bits⟩ := separator_some_wordMass hl hr px py hx hy result
      have eq := agree [f] (by simpa using hf) [true] (basisDensity (wedgeInput e f))
      rcases bits with ⟨lzero, rone⟩ | ⟨lone, rzero⟩
      · rw [lzero, rone] at eq
        exact False.elim (zero_ne_one eq)
      · rw [lone, rzero] at eq
        exact False.elim (one_ne_zero eq)

private theorem bounded_legal_exists {n : Nat} {menu : List (BEdge n)}
    (c : Candidate (board menu)) : ∃ xs : List (BEdge n), xs.Perm menu ∧ Ordered c.rank (board xs) := by
  have perm : c.schedule.reverse.Perm (menu.map Subtype.val) :=
    c.schedule.reverse_perm.trans (c.complete.trans (List.reverse_perm _))
  have lifted : Relation.Comp (fun l xs => l = List.map Subtype.val xs)
      List.Perm c.schedule.reverse menu := by
    rwa [List.eq_map_comp_perm]
  obtain ⟨xs, eq, hx⟩ := lifted
  exact ⟨xs, hx, by simpa [board, ← eq] using c.legal⟩

theorem separator_none_iff_all_quantum_schedules {n : Nat} {menu : List (BEdge n)}
    {c d : Candidate (board menu)} {left right : HitOracle}
    (hl : Realizes c left) (hr : Realizes d right) :
    separator (board menu) left right = none ↔
    ∀ xs ys : List (BEdge n), xs.Perm menu → ys.Perm menu →
      Ordered c.rank (board xs) → Ordered d.rank (board ys) → MarginalAgreement menu xs ys := by
  constructor
  · intro no xs ys px py hx hy
    exact (separator_none_iff_quantum_marginals hl hr px py hx hy).mp no
  · intro all
    obtain ⟨xs, px, hx⟩ := bounded_legal_exists c
    obtain ⟨ys, py, hy⟩ := bounded_legal_exists d
    exact (separator_none_iff_quantum_marginals hl hr px py hx hy).mpr (all xs ys px py hx hy)

#print axioms incident_full_channel
#print axioms separator_some_wordMass
#print axioms separator_none_iff_quantum_marginals
#print axioms separator_none_iff_all_quantum_schedules
end ObservablePriority.Quantum
