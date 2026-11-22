import ObserverPatchHolography.ObservablePriority.PatchComposition

/-! Active identification on a fixed finite simple screen. The sole observable is
    the occupancy bit at a named address after a complete legal run. -/
namespace ObservablePriority
open ObserverPatchSharp

/-- The witness of a complete order prevents vacuous agreement for a tied or
    cyclic installation. The screen itself need not be in priority order. -/
structure Candidate (screen : List Edge) where
  rank : Edge → Nat
  schedule : List Edge
  complete : schedule.Perm screen
  legal : Ordered rank schedule

abbrev HitOracle := Source → Edge → Bool

def hit (rank : Edge → Nat) (schedule : List Edge) : HitOracle :=
  fun L e => ((sharpRun rank schedule L).collar e).isSome

/-- Each oracle answer can come from a different complete legal schedule. -/
def Realizes {screen : List Edge} (c : Candidate screen) (observe : HitOracle) : Prop :=
  ∀ L e, e ∈ screen → ∃ xs : List Edge, xs.Perm screen ∧ Ordered c.rank xs ∧
    observe L e = hit c.rank xs L e

theorem canonical_realizes {screen : List Edge} (c : Candidate screen) :
    Realizes c (hit c.rank c.schedule) :=
  fun _ _ _ => ⟨c.schedule, c.complete, c.legal, rfl⟩

private theorem observe_eq {screen : List Edge} {c : Candidate screen} {o : HitOracle}
    (h : Realizes c o) (L : Source) {e : Edge} (he : e ∈ screen) :
    o L e = hit c.rank c.schedule L e := by
  obtain ⟨xs, hp, ho, hv⟩ := h L e he
  rw [hv]
  unfold hit
  rw [sharpRun_eq_of_legal_perm c.rank L (hp.trans c.complete.symm) ho c.legal]

def IncidentAgreement {screen : List Edge} (c d : Candidate screen) : Prop :=
  ∀ e ∈ screen, ∀ f ∈ screen, e ≠ f → meets e f = true →
    (c.rank f < c.rank e ↔ d.rank f < d.rank e)

def ObservationalAgreement {screen : List Edge} (c d : Candidate screen) : Prop :=
  ∀ xs ys : List Edge, xs.Perm screen → ys.Perm screen →
    Ordered c.rank xs → Ordered d.rank ys →
    ∀ (L : Source) (e : Edge), e ∈ screen → hit c.rank xs L e = hit d.rank ys L e

private theorem incident_strict {r : Edge → Nat} {xs : List Edge} (h : Ordered r xs)
    {e f : Edge} (he : e ∈ xs) (hf : f ∈ xs) (hne : e ≠ f) (hm : meets e f = true) :
    r f < r e ∨ r e < r f := by
  induction xs with
  | nil => simp at he
  | cons a xs ih =>
    obtain ⟨ha, ht⟩ := List.pairwise_cons.mp h
    rcases List.mem_cons.mp he with rfl | het
    · exact .inl ((ha f ((List.mem_cons.mp hf).resolve_left (Ne.symm hne))).2 hm)
    rcases List.mem_cons.mp hf with rfl | hft
    · exact .inr ((ha e het).2 ((meets_symm ..).trans hm))
    · exact ih ht het hft

/-- The old separating preparation exposes its two full-screen addresses.
    The boundary witness gives the two global bit values. -/
private theorem wedge_two_bits {r : Edge → Nat} {xs : List Edge} (ho : Ordered r xs)
    {e f : Edge} (he : e ∈ xs) (hf : f ∈ xs) (hm : meets e f = true) (lt : r f < r e) :
    hit r xs (wedge e f) e = false ∧ hit r xs (wedge e f) f = true := by
  have ne : f ≠ e := by rintro rfl; omega
  obtain ⟨miss, win, _⟩ := boundary_witness r xs (fun x => decide (x = e)) ho he hf
    (by simp) (by simp [ne]) hm lt
  simp [hit, miss, win]

theorem wedge_hit {r : Edge → Nat} {xs : List Edge} (ho : Ordered r xs)
    {e f : Edge} (he : e ∈ xs) (hf : f ∈ xs) (hne : e ≠ f) (hm : meets e f = true) :
    hit r xs (wedge e f) f = decide (r f < r e) := by
  rcases incident_strict ho he hf hne hm with higher | lower
  · simp [higher, (wedge_two_bits ho he hf hm higher).2]
  · have miss := (wedge_two_bits ho hf he ((meets_symm ..).trans hm) lower).1
    rw [wedge_symm] at miss
    simp [miss, Nat.not_lt_of_ge (Nat.le_of_lt lower)]

theorem oracle_wedge_hit {screen : List Edge} {c : Candidate screen} {o : HitOracle}
    (service : Realizes c o) {e f : Edge} (he : e ∈ screen) (hf : f ∈ screen)
    (hne : e ≠ f) (hm : meets e f = true) :
    o (wedge e f) f = decide (c.rank f < c.rank e) := by
  rw [observe_eq service _ hf]
  exact wedge_hit c.legal (c.complete.mem_iff.mpr he) (c.complete.mem_iff.mpr hf) hne hm

theorem observational_iff_incident {screen : List Edge} (c d : Candidate screen) :
    ObservationalAgreement c d ↔ IncidentAgreement c d := by
  constructor
  · intro agree e he f hf hne hm
    have hv := agree c.schedule d.schedule c.complete d.complete c.legal d.legal (wedge e f) f hf
    rw [wedge_hit c.legal (c.complete.mem_iff.mpr he) (c.complete.mem_iff.mpr hf) hne hm,
        wedge_hit d.legal (d.complete.mem_iff.mpr he) (d.complete.mem_iff.mpr hf) hne hm] at hv
    simpa using hv
  · intro same xs ys px py hx hy L e _
    have hcy : Ordered c.rank ys := hy.imp_of_mem (fun ha hb h => ⟨h.1,
      fun hm => (same _ (py.mem_iff.mp ha) _ (py.mem_iff.mp hb) h.1 hm).mpr (h.2 hm)⟩)
    unfold hit
    rw [sharpRun_eq_of_legal_perm c.rank L (px.trans py.symm) hx hcy,
        sharpRun_eq_material c.rank ys L hcy, sharpRun_eq_material d.rank ys L hy]

/-- The experiment names an ordered pair for reading. Its two edges receive
    the same physical wedge preparation whichever candidate is queried. -/
def eligible (e f : Edge) : Bool := decide (e ≠ f) && meets e f

def recoverParents (screen : List Edge) (observe : HitOracle) (e : Edge) : List Edge :=
  screen.filter fun f => eligible e f && observe (wedge e f) f

theorem recoverParents_eq {screen : List Edge} {c : Candidate screen} {o : HitOracle}
    (service : Realizes c o) {e : Edge} (he : e ∈ screen) :
    recoverParents screen o e = PatchComposition.parents c.rank screen e := by
  apply List.filter_congr
  intro f hf
  by_cases hne : e = f
  · subst f
    simp [eligible]
  · cases hm : meets e f with
    | false => simp [eligible, hm]
    | true => simp [eligible, hne, hm, oracle_wedge_hit service he hf hne hm]

theorem recoverParents_member {screen : List Edge} {c : Candidate screen} {o : HitOracle}
    (service : Realizes c o) {e f : Edge} (he : e ∈ screen) :
    f ∈ recoverParents screen o e ↔ f ∈ screen ∧ meets e f = true ∧ c.rank f < c.rank e := by
  rw [recoverParents_eq service he, PatchComposition.physical_member]

private def pairs (screen : List Edge) : List (Edge × Edge) :=
  screen.flatMap fun e => (screen.filter (eligible e)).map fun f => (e, f)

private theorem pairs_spec {screen : List Edge} {e f : Edge} :
    (e, f) ∈ pairs screen ↔ e ∈ screen ∧ f ∈ screen ∧ e ≠ f ∧ meets e f = true := by
  simp [pairs, eligible]

private def differs (left right : HitOracle) (pair : Edge × Edge) : Bool :=
  left (wedge pair.1 pair.2) pair.2 != right (wedge pair.1 pair.2) pair.2

/-- Some (e, f) supplies the preparation wedge e f and the readout address f.
    Only addresses and actual readouts enter the search. -/
def separator (screen : List Edge) (left right : HitOracle) : Option (Edge × Edge) :=
  (pairs screen).find? (differs left right)

theorem separator_none_iff {screen : List Edge} {c d : Candidate screen}
    {left right : HitOracle} (hl : Realizes c left) (hr : Realizes d right) :
    separator screen left right = none ↔ ObservationalAgreement c d := by
  rw [observational_iff_incident, separator, List.find?_eq_none]
  constructor
  · intro no e he f hf hne hm
    have impossible := no (e, f) (pairs_spec.mpr ⟨he, hf, hne, hm⟩)
    simpa [differs, oracle_wedge_hit hl he hf hne hm, oracle_wedge_hit hr he hf hne hm]
      using impossible
  · intro same ⟨e, f⟩ hp
    obtain ⟨he, hf, hne, hm⟩ := pairs_spec.mp hp
    simp [differs, oracle_wedge_hit hl he hf hne hm, oracle_wedge_hit hr he hf hne hm,
      same e he f hf hne hm]

/-- Both the oracle readings and the readings on any other pair of legal
    schedules are the complementary bits of this concrete physical experiment. -/
theorem separator_some_bits {screen : List Edge} {c d : Candidate screen}
    {left right : HitOracle} (hl : Realizes c left) (hr : Realizes d right)
    {e f : Edge} (found : separator screen left right = some (e, f)) :
    e ∈ screen ∧ f ∈ screen ∧ e ≠ f ∧ meets e f = true ∧
    ∀ xs ys : List Edge, xs.Perm screen → ys.Perm screen →
      Ordered c.rank xs → Ordered d.rank ys →
      hit c.rank xs (wedge e f) f = left (wedge e f) f ∧
      hit d.rank ys (wedge e f) f = right (wedge e f) f ∧
      ((left (wedge e f) f = false ∧ right (wedge e f) f = true) ∨
       (left (wedge e f) f = true ∧ right (wedge e f) f = false)) := by
  have mem := List.mem_of_find?_eq_some found
  obtain ⟨he, hf, hne, hm⟩ := pairs_spec.mp mem
  have diff : differs left right (e, f) = true := List.find?_some found
  refine ⟨he, hf, hne, hm, fun xs ys px py hx hy => ⟨?_, ?_, ?_⟩⟩
  · rw [wedge_hit hx (px.mem_iff.mpr he) (px.mem_iff.mpr hf) hne hm,
        oracle_wedge_hit hl he hf hne hm]
  · rw [wedge_hit hy (py.mem_iff.mpr he) (py.mem_iff.mpr hf) hne hm,
        oracle_wedge_hit hr he hf hne hm]
  · cases l : left (wedge e f) f <;> cases r : right (wedge e f) f <;>
      simp_all [differs]

/-- The producer can accept an independently reconstructed service order.
    Every check here uses readouts and addresses or the supplied order. -/
def checkedClosure (screen : List Edge) (o : HitOracle) (seeds order : List Edge) :
    Option (List Edge) :=
  let parents := recoverParents screen o
  if order.Perm screen ∧ seeds ⊆ order ∧ PatchClosure.orderedCheck parents order = true then
    some (PatchClosure.generate seeds parents order)
  else none

private theorem needed_equiv {screen seeds : List Edge} {c : Candidate screen}
    {o : HitOracle} (service : Realizes c o) (hseed : seeds ⊆ screen) (a : Edge) :
    PatchClosure.Needed seeds (recoverParents screen o) a ↔
    PatchClosure.Needed seeds (PatchComposition.parents c.rank screen) a := by
  constructor
  · intro h
    suffices a ∈ screen ∧ PatchClosure.Needed seeds (PatchComposition.parents c.rank screen) a from this.2
    induction h with
    | seed x hx => exact ⟨hseed hx, .seed x hx⟩
    | parent prev hp ih =>
      rw [recoverParents_eq service ih.1] at hp
      exact ⟨(PatchComposition.physical_member.mp hp).1, .parent ih.2 hp⟩
  · intro h
    suffices a ∈ screen ∧ PatchClosure.Needed seeds (recoverParents screen o) a from this.2
    induction h with
    | seed x hx => exact ⟨hseed hx, .seed x hx⟩
    | parent prev hp ih =>
      exact ⟨(PatchComposition.physical_member.mp hp).1, .parent ih.2
        (by rw [recoverParents_eq service ih.1]; exact hp)⟩

theorem recovered_closure_exact {screen seeds order cert : List Edge} {c : Candidate screen}
    {o : HitOracle} (service : Realizes c o) (found : checkedClosure screen o seeds order = some cert) :
    PatchClosure.acceptCertificate seeds (recoverParents screen o) cert = true ∧
    ∀ a, a ∈ cert ↔ PatchClosure.Needed seeds (PatchComposition.parents c.rank screen) a := by
  dsimp only [checkedClosure] at found
  split at found
  next ok =>
    obtain ⟨perm, seedIn, legal⟩ := ok
    have seedInScreen : seeds ⊆ screen := fun _ h => perm.mem_iff.mp (seedIn h)
    obtain ⟨unique, noSelf, earlier⟩ := (PatchClosure.orderedCheck_spec _ _).mp legal
    have install : PatchClosure.Installation (recoverParents screen o) order :=
      ⟨unique, noSelf, fun a ha b hb =>
        perm.mem_iff.mpr ((recoverParents_member service (perm.mem_iff.mp ha)).mp hb).1, earlier⟩
    cases found
    exact ⟨PatchClosure.producer_accepted install seedIn, fun a =>
      (PatchClosure.generate_exact install seedIn a).trans (needed_equiv service seedInScreen a)⟩
  next => contradiction

/-- Every legal schedule supplies an order the rank-blind closure gate accepts.
    Sharp input lists are executed tail first, while the closure consumer
    expects its argument in forward service order. -/
theorem legal_closure_accepted {screen seeds xs : List Edge} {c : Candidate screen}
    {o : HitOracle} (service : Realizes c o) (p : xs.Perm screen)
    (ho : Ordered c.rank xs) (hs : seeds ⊆ screen) :
    checkedClosure screen o seeds xs.reverse =
      some (PatchClosure.generate seeds (recoverParents screen o) xs.reverse) := by
  have perm : xs.reverse.Perm screen := xs.reverse_perm.trans p
  have seedIn : seeds ⊆ xs.reverse := fun _ h => perm.mem_iff.mpr (hs h)
  have ordered : PatchClosure.orderedCheck (recoverParents screen o) xs.reverse = true := by
    apply (PatchClosure.orderedCheck_spec ..).mpr
    refine ⟨List.pairwise_reverse.mpr (ho.imp (fun h => Ne.symm h.1)), ?_, List.pairwise_reverse.mpr ?_⟩
    · intro e he h
      have := ((recoverParents_member service (perm.mem_iff.mp he)).mp h).2.2
      omega
    · apply ho.imp_of_mem
      intro a b ha hb h hpar
      obtain ⟨_, hm, lt⟩ := (recoverParents_member service (p.mem_iff.mp hb)).mp hpar
      have := h.2 ((meets_symm ..).trans hm)
      omega
  simp [checkedClosure, perm, seedIn, ordered]

#print axioms observational_iff_incident
#print axioms recoverParents_eq
#print axioms separator_none_iff
#print axioms separator_some_bits
#print axioms recovered_closure_exact
#print axioms legal_closure_accepted
end ObservablePriority
