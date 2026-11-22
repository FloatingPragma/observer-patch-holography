import ObserverPatchHolography.ObservablePriority.SharpOrder
import ObserverPatchHolography.ObservablePriority.ObservedChannel

/-! A finite numbered OPH screen. Each quantum sign bit is matched with the
    oriented integer written by the current-load sharp service. -/
namespace ObserverPatchSharp.QuantumBridge

open scoped BigOperators

abbrev BEdge (n : Nat) := {e : Edge // e.hi < n}
def screen (n : Nat) : StrictQuantum.Screen (BEdge n) (Fin n) where
  left := fun e => ⟨e.1.lo, Nat.lt_trans e.1.sorted e.2⟩
  right := fun e => ⟨e.1.hi, e.2⟩
  distinct := fun e h => (Nat.ne_of_lt e.1.sorted) (congrArg Fin.val h)

def decode : StrictQuantum.Load → Load
  | .empty => .zero | .single => .one | .double => .two
def contrast (b : StrictQuantum.Sign) : Int := if b then 2 else -2

def source {n : Nat} (x : StrictQuantum.Loads (Fin n)) : Source := fun p =>
  if h : p < n then decode (x ⟨p, h⟩) else .one

def Reads {n : Nat} (s : Material) (x : StrictQuantum.Loads (Fin n)) :=
  ∀ p : Fin n, s.load p.val = decode (x p)

def board {n : Nat} (es : List (BEdge n)) : List Edge := (es.map Subtype.val).reverse

private theorem service_collar_other (s : Material) (i j : Edge) (h : j ≠ i) :
    (service i s).collar j = s.collar j := by
  unfold service
  split <;> simp [write, h]

private theorem one_step {n : Nat} (e : BEdge n) (s : Material)
    (x : StrictQuantum.Loads (Fin n)) (hx : Reads s x) (hb : s.collar e.val = none) :
    Reads (service e.val s) (StrictQuantum.fire (screen n) e x) ∧
      (StrictQuantum.tag (screen n) e x).map contrast = (service e.val s).collar e.val := by
  let q := screen n
  have hl : s.load e.1.lo = decode (x (q.left e)) := hx (q.left e)
  have hr : s.load e.1.hi = decode (x (q.right e)) := hx (q.right e)
  have hready : ready s e.val = (StrictQuantum.tag q e x).isSome := by
    simp only [ready, initialReady, hb, hl, hr, StrictQuantum.tag]
    cases x (q.left e) <;> cases x (q.right e) <;> rfl
  constructor
  · intro p
    have hleft : p.val = e.1.lo ↔ p = q.left e :=
      ⟨fun h => Fin.ext h, fun h => congrArg Fin.val h⟩
    have hright : p.val = e.1.hi ↔ p = q.right e :=
      ⟨fun h => Fin.ext h, fun h => congrArg Fin.val h⟩
    change (service e.1 s).load p.val = decode (StrictQuantum.fire q e x p)
    cases h : (StrictQuantum.tag q e x).isSome with
    | false => simpa [service, StrictQuantum.fire, hready, h] using hx p
    | true =>
      by_cases pl : p = q.left e <;> by_cases pr : p = q.right e <;>
        simp [service, hready, h, write, incident, hleft, hright, StrictQuantum.fire,
          StrictQuantum.putPair, pl, pr, decode, hx p]
      all_goals simp [q, screen]
  · change (StrictQuantum.tag q e x).map contrast = (service e.val s).collar e.val
    simp only [service, hready, StrictQuantum.tag]
    cases xl : x (q.left e) <;> cases xr : x (q.right e) <;>
      simp [StrictQuantum.sign, write, sign, hl, hr, xl, xr, decode, contrast, Load.value, hb]

def advance {n : Nat} (es : List (BEdge n)) (s : Material) :=
  es.foldl (fun st e => service e.val st) s

private theorem unvisited {n : Nat} (es : List (BEdge n)) (s : Material) (i : Edge)
    (h : i ∉ es.map Subtype.val) : (advance es s).collar i = s.collar i := by
  induction es generalizing s with
  | nil => rfl
  | cons e rest ih =>
    have hn : i ≠ e.1 ∧ i ∉ rest.map Subtype.val := by simpa using h
    change (advance rest (service e.val s)).collar i = s.collar i
    rw [ih _ hn.2, service_collar_other _ _ _ hn.1]

private theorem scan_advance {n : Nat} (es : List (BEdge n)) (s : Material)
    (x : StrictQuantum.Loads (Fin n)) (hn : (es.map Subtype.val).Nodup)
    (hx : Reads s x) (hb : ∀ e ∈ es, s.collar e.val = none) :
    Reads (advance es s) (StrictQuantum.scan (screen n) es x).1 ∧
    (StrictQuantum.scan (screen n) es x).2.map (Option.map contrast) =
      es.map (fun e => (advance es s).collar e.val) := by
  induction es generalizing s x with
  | nil => exact ⟨hx, rfl⟩
  | cons e rest ih =>
    obtain ⟨hnHead, hnTail⟩ := List.nodup_cons.mp (by exact hn)
    obtain ⟨hx', htag⟩ := one_step e s x hx (hb e (by simp))
    have hb' : ∀ j ∈ rest, (service e.val s).collar j.val = none := by
      intro j hj
      have hne : j.val ≠ e.val := by
        intro hEq
        exact hnHead (List.mem_map.mpr ⟨j, hj, hEq⟩)
      rw [service_collar_other _ _ _ hne, hb j (by simp [hj])]
    obtain ⟨hl, hs⟩ := ih _ _ hnTail hx' hb'
    refine ⟨hl, ?_⟩
    change _ :: _ = (advance rest (service e.val s)).collar e.val :: _
    rw [unvisited _ _ _ hnHead]
    exact congrArg₂ List.cons htag hs

theorem signed_scan_sharp {n : Nat} (rank : Edge → Nat) (es : List (BEdge n))
    (x : StrictQuantum.Loads (Fin n)) (ho : Ordered rank (board es)) :
    Reads (sharpRun rank (board es) (source x)) (StrictQuantum.scan (screen n) es x).1 ∧
    (StrictQuantum.scan (screen n) es x).2.map (Option.map contrast) =
      es.map (fun e => (sharpRun rank (board es) (source x)).collar e.val) := by
  have hn : (board es).Nodup := ho.imp (fun h => h.1)
  have hn' : (es.map Subtype.val).Nodup := by simpa [board] using hn
  have ha : sharpRun rank (board es) (source x) = advance es (initial (source x)) := by
    rw [sharpRun_eq_material rank _ _ ho, ← serviceFold_eq_material]
    simp [board, advance]
  rw [ha]
  apply scan_advance _ _ _ hn'
  · intro p
    simp [initial, source, p.isLt]
  · simp [initial]

theorem flags_sharp {n : Nat} (rank : Edge → Nat) (es : List (BEdge n))
    (x : StrictQuantum.Loads (Fin n)) (ho : Ordered rank (board es)) :
    StrictQuantum.flags (StrictQuantum.scan (screen n) es x).2 =
      es.map (fun e => ((sharpRun rank (board es) (source x)).collar e.val).isSome) := by
  have h := congrArg (List.map Option.isSome) (signed_scan_sharp rank es x ho).2
  simpa [StrictQuantum.flags, List.map_map, Function.comp_def] using h

private theorem mem_board {n : Nat} {e : BEdge n} {es : List (BEdge n)} :
    e.val ∈ board es ↔ e ∈ es := by
  simp only [board, List.mem_reverse, List.mem_map]
  constructor
  · rintro ⟨i, hi, he⟩
    have : i = e := Subtype.ext he
    simpa [this] using hi
  · exact fun h => ⟨e, h, rfl⟩

def bitAt {A : Type} [DecidableEq A] (e : A) : List A → List Bool → Bool
  | a :: es, b :: bs => if e = a then b else bitAt e es bs
  | _, _ => false

theorem bitAt_map {A : Type} [DecidableEq A] (f : A → Bool) (es : List A)
    (e : A) (he : e ∈ es) : bitAt e es (es.map f) = f e := by
  induction es with
  | nil => simp at he
  | cons a tail ih =>
    by_cases h : e = a
    · subst e
      simp [bitAt]
    · have ht := (List.mem_cons.mp he).resolve_left h
      simp [bitAt, h, ih ht]

def readWord {A : Type} [DecidableEq A] (es wanted : List A) (w : List Bool) :=
  wanted.map (fun e => bitAt e es w)

theorem closed_readout {n : Nat} (rank : Edge → Nat) (full localOrder wanted : List (BEdge n))
    (keep : Edge → Bool) (hf : Ordered rank (board full)) (hl : Ordered rank (board localOrder))
    (hp : (board localOrder).Perm ((board full).filter keep))
    (hc : PredecessorClosed rank (board full) keep)
    (hw : ∀ e ∈ wanted, e ∈ full ∧ keep e.val = true) (x : StrictQuantum.Loads (Fin n)) :
    readWord full wanted (StrictQuantum.flags (StrictQuantum.scan (screen n) full x).2) =
    readWord localOrder wanted (StrictQuantum.flags (StrictQuantum.scan (screen n) localOrder x).2) := by
  rw [flags_sharp rank full x hf, flags_sharp rank localOrder x hl]
  apply List.map_congr_left
  intro e he
  obtain ⟨hef, hkeep⟩ := hw e he
  have hel : e ∈ localOrder :=
    mem_board.mp (hp.mem_iff.mpr (List.mem_filter.mpr ⟨mem_board.mpr hef, hkeep⟩))
  rw [bitAt_map _ _ _ hef, bitAt_map _ _ _ hel]
  have hh := (predecessorClosed_iff_allSchedulesAgree rank _ keep hf).mp hc
    (board full) (board localOrder) .rfl hp hf hl (source x) e.val (mem_board.mpr hef) hkeep
  rw [hh]

open scoped BigOperators in
/-- The addressed branch trace of the actual labelled Kraus columns uses the
    acceptance bits computed by the physical rank-guarded sharp execution. -/
theorem branch_trace_sharp {n : Nat} (rank : Edge → Nat) (es : List (BEdge n))
    (ho : Ordered rank (board es)) (M : BEdge n → Bool) (ρ : StrictQuantum.Observed.density (Fin n)) :
    (∑ x, ∑ c : StrictQuantum.Coins (StrictQuantum.scan (screen n) es x).2,
      StrictQuantum.Observed.entry (screen n) es M ρ x x c c) =
    ∑ x, if es.map (fun e => ((sharpRun rank (board es) (source x)).collar e.val).isSome) = es.map M
         then ρ x x else 0 := by
  rw [StrictQuantum.Observed.branch_trace]
  unfold StrictQuantum.Observed.mass
  apply Finset.sum_congr rfl
  intro x _
  rw [flags_sharp rank es x ho]

open scoped BigOperators in
/-- Sum over all unobserved outcomes of the sequential source instrument. Two
    independent legal schedules give the same marginal on any retained address
    word, for every density matrix on the common finite input screen. -/
theorem closed_patch_marginal {n : Nat} (rank : Edge → Nat)
    (full localOrder wanted : List (BEdge n)) (keep : Edge → Bool)
    (hf : Ordered rank (board full)) (hl : Ordered rank (board localOrder))
    (hp : (board localOrder).Perm ((board full).filter keep))
    (hc : PredecessorClosed rank (board full) keep)
    (hw : ∀ e ∈ wanted, e ∈ full ∧ keep e.val = true)
    (obs : List Bool) (ρ : StrictQuantum.Observed.density (Fin n)) :
    (∑ w ∈ (StrictQuantum.Observed.outcomes (screen n) full).filter
        (fun w => readWord full wanted w = obs),
      StrictQuantum.Observed.mass (screen n) full w (fun x => ρ x x)) =
    (∑ w ∈ (StrictQuantum.Observed.outcomes (screen n) localOrder).filter
        (fun w => readWord localOrder wanted w = obs),
      StrictQuantum.Observed.mass (screen n) localOrder w (fun x => ρ x x)) := by
  rw [StrictQuantum.Observed.marginal, StrictQuantum.Observed.marginal]
  apply Finset.sum_congr rfl
  intro x _
  rw [closed_readout rank full localOrder wanted keep hf hl hp hc hw x]

def encode : Load → StrictQuantum.Load
  | .zero => .empty | .one => .single | .two => .double
@[simp] theorem decode_encode (a : Load) : decode (encode a) = a := by cases a <;> rfl

def wedgeInput {n : Nat} (e f : BEdge n) : StrictQuantum.Loads (Fin n) :=
  fun p => encode (wedge e.val f.val p.val)

theorem source_wedge {n : Nat} (e f : BEdge n) : source (wedgeInput e f) = wedge e.val f.val := by
  funext p
  by_cases h : p < n
  · simp [source, h, wedgeInput]
  · have outside (j : BEdge n) : incident j.val p = false := by
      have hi := j.property
      have lo := j.val.sorted
      simp only [incident, decide_eq_false_iff_not]
      omega
    simp [source, h, wedge, outside]

noncomputable def basisDensity {n : Nat} (x : StrictQuantum.Loads (Fin n)) :
    StrictQuantum.Observed.density (Fin n) := fun y z => if y = x ∧ z = x then 1 else 0

noncomputable def wordMass {n : Nat} (es wanted : List (BEdge n)) (obs : List Bool)
    (ρ : StrictQuantum.Observed.density (Fin n)) : ℂ :=
  ∑ w ∈ (StrictQuantum.Observed.outcomes (screen n) es).filter (fun w => readWord es wanted w = obs),
    StrictQuantum.Observed.mass (screen n) es w (fun x => ρ x x)

theorem basis_wordMass {n : Nat} (es wanted : List (BEdge n)) (obs : List Bool)
    (x : StrictQuantum.Loads (Fin n)) :
    wordMass es wanted obs (basisDensity x) =
      if readWord es wanted (StrictQuantum.flags (StrictQuantum.scan (screen n) es x).2) = obs
      then 1 else 0 := by
  classical
  rw [wordMass, StrictQuantum.Observed.marginal]
  rw [Finset.sum_eq_single x]
  · simp [basisDensity]
  · intro y _ hne
    split <;> simp [basisDensity, hne]
  · simp

/-- Every omitted higher incident seam has a physical basis preparation for
    which observing the retained seam's hit has probabilities zero and one. -/
theorem boundary_patch_probability {n : Nat} (rank : Edge → Nat)
    (full localOrder : List (BEdge n)) (keep : Edge → Bool)
    (hf : Ordered rank (board full)) (hl : Ordered rank (board localOrder))
    (hp : (board localOrder).Perm ((board full).filter keep))
    (e f : BEdge n) (he : e ∈ full) (hF : f ∈ full) (hke : keep e.val = true)
    (hkf : keep f.val = false) (hm : meets e.val f.val = true) (hr : rank f.val < rank e.val) :
    wordMass full [e] [true] (basisDensity (wedgeInput e f)) = 0 ∧
    wordMass localOrder [e] [true] (basisDensity (wedgeInput e f)) = 1 := by
  have hel : e ∈ localOrder := mem_board.mp
    (hp.mem_iff.mpr (List.mem_filter.mpr ⟨mem_board.mpr he, hke⟩))
  obtain ⟨ha, _, hb⟩ := boundary_witness rank (board full) keep hf
    (mem_board.mpr he) (mem_board.mpr hF) hke hkf hm hr
  have hp' := sharpRun_eq_of_legal_perm rank (wedge e.val f.val) hp hl (List.Pairwise.filter keep hf)
  simp only [basis_wordMass, flags_sharp rank full _ hf, flags_sharp rank localOrder _ hl,
    readWord, List.map_cons, List.map_nil, bitAt_map _ full _ he, bitAt_map _ localOrder _ hel,
    source_wedge, hp', ha, hb]
  exact ⟨rfl, rfl⟩

/-- The geometric condition characterizes equality of all addressed quantum
    outcome marginals. Necessity uses only a single address and a pure basis density. -/
theorem closed_iff_quantum_marginals {n : Nat} (rank : Edge → Nat)
    (full localOrder : List (BEdge n)) (keep : Edge → Bool)
    (hf : Ordered rank (board full)) (hl : Ordered rank (board localOrder))
    (hp : (board localOrder).Perm ((board full).filter keep)) :
    PredecessorClosed rank (board full) keep ↔
      ∀ (wanted : List (BEdge n)), (∀ e ∈ wanted, e ∈ full ∧ keep e.val = true) →
        ∀ (obs : List Bool) (ρ : StrictQuantum.Observed.density (Fin n)),
          wordMass full wanted obs ρ = wordMass localOrder wanted obs ρ := by
  constructor
  · intro hc wanted hw obs ρ
    exact closed_patch_marginal rank full localOrder wanted keep hf hl hp hc hw obs ρ
  · intro h a ha hka b hb hm hr
    have ha' : a ∈ full.map Subtype.val := by simpa [board] using ha
    have hb' : b ∈ full.map Subtype.val := by simpa [board] using hb
    obtain ⟨e, he, rfl⟩ := List.mem_map.mp ha'
    obtain ⟨f, hF, rfl⟩ := List.mem_map.mp hb'
    cases hk : keep f.val with
    | true => rfl
    | false =>
      obtain ⟨hx, hy⟩ := boundary_patch_probability rank full localOrder keep hf hl hp
        e f he hF hka hk hm hr
      have hword : ∀ j ∈ [e], j ∈ full ∧ keep j.val = true := by
        intro j hj
        have : j = e := by simpa using hj
        subst j
        exact ⟨he, hka⟩
      have hh := h [e] hword [true] (basisDensity (wedgeInput e f))
      rw [hx, hy] at hh
      exact False.elim (zero_ne_one hh)

#print axioms signed_scan_sharp
#print axioms flags_sharp
#print axioms closed_readout
#print axioms branch_trace_sharp
#print axioms closed_patch_marginal
#print axioms boundary_patch_probability
#print axioms closed_iff_quantum_marginals

end ObserverPatchSharp.QuantumBridge
