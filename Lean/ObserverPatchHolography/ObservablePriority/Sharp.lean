import Std

/-! Finite priority-ordered strict repair, signed addressed collars, and observer patches.
    The input list stores the strongest seams at the tail, which is executed first. -/
namespace ObserverPatchSharp

variable {A : Type}

def accepted (ready : A → Bool) (conflict : A → A → Bool) : List A → List A
  | [] => []
  | a :: prior =>
    let selected := accepted ready conflict prior
    if ready a && !selected.any (conflict a) then a :: selected else selected

def closed (keep : A → Bool) (conflict : A → A → Bool) : List A → Bool
  | [] => true
  | a :: prior =>
    (if keep a then prior.all (fun b => !conflict a b || keep b) else true) &&
      closed keep conflict prior

theorem accepted_mem {ready : A → Bool} {conflict : A → A → Bool}
    {s : List A} {b : A} (hb : b ∈ accepted ready conflict s) : b ∈ s := by
  induction s with
  | nil => exact hb
  | cons a prior ih =>
    simp only [accepted] at hb
    split at hb
    · rcases List.mem_cons.mp hb with h | h
      · exact List.mem_cons.mpr (Or.inl h)
      · exact List.mem_cons_of_mem a (ih h)
    · exact List.mem_cons_of_mem a (ih hb)

private theorem surviving_conflict (keep test : A → Bool) (xs : List A)
    (h : ∀ b ∈ xs, test b = true → keep b = true) :
    (xs.filter keep).any test = xs.any test := by
  apply Bool.eq_iff_iff.mpr
  simp only [List.any_eq_true, List.mem_filter]
  constructor
  · rintro ⟨b, ⟨hb, _⟩, ht⟩
    exact ⟨b, hb, ht⟩
  · rintro ⟨b, hb, ht⟩
    exact ⟨b, ⟨hb, h b hb ht⟩, ht⟩

/-- An earlier conflict of every retained job remains visible. For every
    initial readiness bitmap this suffices to commute selection with projection. -/
theorem closed_projection (keep ready : A → Bool) (conflict : A → A → Bool)
    (stream : List A) (h : closed keep conflict stream = true) :
    (accepted ready conflict stream).filter keep =
      accepted ready conflict (stream.filter keep) := by
  induction stream with
  | nil => rfl
  | cons a prior ih =>
    obtain ⟨hc, ht⟩ := Bool.and_eq_true_iff.mp h
    have ip := ih ht
    by_cases hk : keep a = true
    · have priorKept : ∀ b ∈ prior, conflict a b = true → keep b = true := by
        intro b hb hconf
        have hall : prior.all (fun b => !conflict a b || keep b) = true := by
          simpa [hk] using hc
        simpa [hconf] using List.all_eq_true.mp hall b hb
      have same := surviving_conflict keep (conflict a) (accepted ready conflict prior)
        (fun b hb hconf => priorKept b (accepted_mem hb) hconf)
      rw [ip] at same
      simp only [List.filter_cons, hk, if_true, accepted, same]
      split <;> simp [hk, ip]
    · have hf : keep a = false := Bool.eq_false_iff.mpr hk
      simp only [List.filter_cons, hf, Bool.false_eq_true, if_false, accepted]
      split <;> simp [hf, ip]

/-- A sorted endpoint pair names a single undirected loop-free edge. A list
    with distinct such pairs represents any finite simple screen on port IDs. -/
structure Edge where
  lo : Nat
  hi : Nat
  sorted : lo < hi
  deriving DecidableEq, Repr

def incident (e : Edge) (p : Nat) : Bool := decide (p = e.lo ∨ p = e.hi)
def meets (a b : Edge) : Bool := incident b a.lo || incident b a.hi

theorem meets_symm (a b : Edge) : meets a b = meets b a := by
  apply Bool.eq_iff_iff.mpr
  simp only [meets, incident, Bool.or_eq_true, decide_eq_true_eq]
  grind

@[simp] theorem meets_self (a : Edge) : meets a a = true := by simp [meets, incident]

inductive Load | zero | one | two deriving DecidableEq, Repr
def Load.value : Load → Int | .zero => 0 | .one => 1 | .two => 2
def gap : Load → Load → Bool | .zero, .two | .two, .zero => true | _, _ => false

@[simp] theorem gap_one_left (x : Load) : gap .one x = false := by cases x <;> rfl
@[simp] theorem gap_one_right (x : Load) : gap x .one = false := by cases x <;> rfl

abbrev Source := Nat → Load
def initialReady (L : Source) (e : Edge) := gap (L e.lo) (L e.hi)
def sign (L : Source) (e : Edge) := (L e.lo).value - (L e.hi).value

structure Material where
  load : Source
  collar : Edge → Option Int

def initial (L : Source) : Material := ⟨L, fun _ => none⟩
def ready (s : Material) (e : Edge) := (s.collar e).isNone && initialReady s.load e

def write (s : Material) (e : Edge) : Material where
  load := fun p => if incident e p then .one else s.load p
  collar := fun f => if f = e then some (sign s.load e) else s.collar f

def enabled (rank : Edge → Nat) (screen : List Edge) (s : Material) (a : Edge) :=
  ready s a && !screen.any (fun b =>
    meets a b && decide (rank b < rank a) && ready s b)

def attempt (rank : Edge → Nat) (screen : List Edge) (s : Material) (a : Edge) :=
  if enabled rank screen s a then write s a else s

def runOn (rank : Edge → Nat) (screen : List Edge) (L : Source) : List Edge → Material
  | [] => initial L
  | a :: prior => attempt rank screen (runOn rank screen L prior) a

def sharpRun (rank : Edge → Nat) (screen : List Edge) (L : Source) :=
  runOn rank screen L screen

def touched (wins : List Edge) (p : Nat) := wins.any (fun b => incident b p)

def material (L : Source) (wins : List Edge) : Material where
  load := fun p => if touched wins p then .one else L p
  collar := fun e => if e ∈ wins then some (sign L e) else none

theorem initial_material (L : Source) : initial L = material L [] := by
  simp [initial, material, touched]

theorem any_meets (wins : List Edge) (e : Edge) :
    wins.any (meets e) = (touched wins e.lo || touched wins e.hi) := by
  apply Bool.eq_iff_iff.mpr
  simp [touched, meets]
  grind

theorem ready_material (L : Source) (wins : List Edge) (e : Edge) :
    ready (material L wins) e = (initialReady L e && !wins.any (meets e)) := by
  rw [any_meets]
  by_cases ha : touched wins e.lo = true
  · simp [ready, material, ha, initialReady]
  by_cases hb : touched wins e.hi = true
  · simp [ready, material, ha, hb, initialReady]
  have hn : e ∉ wins := by
    intro hm
    exact ha (List.any_eq_true.mpr ⟨e, hm, by simp [incident]⟩)
  simp [ready, material, ha, hb, hn, initialReady]

theorem write_material (L : Source) (wins : List Edge) (e : Edge)
    (hr : ready (material L wins) e = true) :
    write (material L wins) e = material L (e :: wins) := by
  have hn := (Bool.and_eq_true_iff.mp ((ready_material L wins e).symm.trans hr)).2
  rw [any_meets] at hn
  have hab : (touched wins e.lo || touched wins e.hi) = false := by simpa using hn
  obtain ⟨ha, hb⟩ := Bool.or_eq_false_iff.mp hab
  unfold write material
  congr 1
  · funext p
    cases hp : incident e p <;> simp [touched, hp]
  · funext f
    by_cases hf : f = e
    · subst f
      simp [sign, ha, hb]
    · simp [hf]

theorem write_energy (s : Material) (e : Edge) (hr : ready s e = true) :
    ((s.load e.lo).value ^ 2 + (s.load e.hi).value ^ 2 =
      ((write s e).load e.lo).value ^ 2 + ((write s e).load e.hi).value ^ 2 + 2) ∧
    (s.load e.lo).value + (s.load e.hi).value =
      ((write s e).load e.lo).value + ((write s e).load e.hi).value ∧
    (∀ p, incident e p = false → (write s e).load p = s.load p) := by
  have hg := (Bool.and_eq_true_iff.mp hr).2
  unfold initialReady at hg
  constructor
  · cases h : s.load e.lo <;> cases k : s.load e.hi <;>
      simp_all [write, incident, gap, Load.value]
  constructor
  · cases h : s.load e.lo <;> cases k : s.load e.hi <;>
      simp_all [write, incident, gap, Load.value]
  · intro p hp
    simp [write, hp]

/-- All originally-ready processed edges have either won or met a winner. -/
theorem covered (r : Edge → Bool) {e : Edge} {xs : List Edge}
    (he : e ∈ xs) (hr : r e = true) :
    (accepted r meets xs).any (meets e) = true := by
  induction xs with
  | nil => simp at he
  | cons a rest ih =>
    simp only [accepted]
    split
    next hy =>
      rcases List.mem_cons.mp he with rfl | hm
      · simp
      · simp [ih hm]
    next hn =>
      rcases List.mem_cons.mp he with rfl | hm
      · cases hb : (accepted r meets rest).any (meets e) <;> simp_all
      · exact ih hm

theorem spent (L : Source) {e : Edge} {xs : List Edge} (he : e ∈ xs) :
    ready (material L (accepted (initialReady L) meets xs)) e = false := by
  rw [ready_material]
  cases hr : initialReady L e with
  | false => rfl
  | true => simp [covered (initialReady L) he hr]

/-- Input orientation: for incident edges, the smaller priority occurs later
    and executes earlier. Disjoint edges can occur in either order. -/
def Ordered (rank : Edge → Nat) (xs : List Edge) :=
  xs.Pairwise (fun a b => a ≠ b ∧ (meets a b = true → rank b < rank a))

def ContainsHigher (rank : Edge → Nat) (screen xs : List Edge) :=
  ∀ a ∈ xs, ∀ b ∈ screen, meets a b = true → rank b < rank a → b ∈ xs

theorem runOn_eq_material (rank : Edge → Nat) (screen xs : List Edge) (L : Source)
    (ho : Ordered rank xs) (hc : ContainsHigher rank screen xs) :
    runOn rank screen L xs = material L (accepted (initialReady L) meets xs) := by
  induction xs with
  | nil => exact initial_material L
  | cons a prior ih =>
    obtain ⟨haf, ht⟩ := List.pairwise_cons.mp ho
    have hc' : ContainsHigher rank screen prior := by
      intro x hx b hb hm hp
      rcases List.mem_cons.mp (hc x (by simp [hx]) b hb hm hp) with h | h
      · subst b
        have := (haf x hx).2 ((meets_symm ..).trans hm)
        omega
      · exact h
    let wins := accepted (initialReady L) meets prior
    have noGuard : screen.any (fun b => meets a b && decide (rank b < rank a) &&
        ready (material L wins) b) = false := by
      apply List.any_eq_false.mpr
      intro b hb h
      obtain ⟨hp, hr⟩ := Bool.and_eq_true_iff.mp h
      obtain ⟨hm, hlt⟩ := Bool.and_eq_true_iff.mp hp
      have hlt' := of_decide_eq_true hlt
      have hbp : b ∈ prior := by
        rcases List.mem_cons.mp (hc a (by simp) b hb hm hlt') with rfl | hp'
        · omega
        · exact hp'
      rw [spent L hbp] at hr
      contradiction
    have hv : enabled rank screen (material L wins) a = ready (material L wins) a := by
      simp [enabled, noGuard]
    simp only [runOn, ih ht hc', attempt]
    change (if enabled rank screen (material L wins) a then _ else _) = _
    rw [hv]
    simp only [accepted]
    change (if ready (material L wins) a then _ else _) =
      material L (if initialReady L a && !wins.any (meets a) then a :: wins else wins)
    rw [← ready_material]
    split
    · exact write_material L wins a (by assumption)
    · rfl

theorem sharpRun_eq_material (rank : Edge → Nat) (screen : List Edge) (L : Source)
    (ho : Ordered rank screen) :
    sharpRun rank screen L = material L (accepted (initialReady L) meets screen) :=
  runOn_eq_material rank screen screen L ho (fun _ _ _ hb _ _ => hb)

theorem sharpRun_quiescent (rank : Edge → Nat) (screen : List Edge) (L : Source)
    (ho : Ordered rank screen) (e : Edge) (he : e ∈ screen) :
    ready (sharpRun rank screen L) e = false := by
  rw [sharpRun_eq_material rank screen L ho]
  exact spent L he

def PredecessorClosed (rank : Edge → Nat) (screen : List Edge) (keep : Edge → Bool) :=
  ∀ a ∈ screen, keep a = true → ∀ b ∈ screen, meets a b = true →
    rank b < rank a → keep b = true

theorem closed_iff_predecessor (rank : Edge → Nat) (screen : List Edge) (keep : Edge → Bool)
    (ho : Ordered rank screen) :
    closed keep meets screen = true ↔ PredecessorClosed rank screen keep := by
  induction screen with
  | nil => simp [closed, PredecessorClosed]
  | cons a prior ih =>
    obtain ⟨hh, ht⟩ := List.pairwise_cons.mp ho
    unfold closed
    rw [Bool.and_eq_true_iff, ih ht]
    constructor
    · rintro ⟨hhead, htail⟩ x hx hk b hb hm hp
      rcases List.mem_cons.mp hx with rfl | hxTail
      · rcases List.mem_cons.mp hb with rfl | hb
        · omega
        · have hall : prior.all (fun b => !meets x b || keep b) = true := by
            simpa [hk] using hhead
          simpa [hm] using List.all_eq_true.mp hall b hb
      · rcases List.mem_cons.mp hb with rfl | hbTail
        · have := (hh x hxTail).2 ((meets_symm ..).trans hm)
          omega
        · exact htail x hxTail hk b hbTail hm hp
    · intro h
      constructor
      · by_cases hk : keep a = true
        · simp only [hk, if_true, List.all_eq_true]
          intro b hb
          by_cases hm : meets a b = true
          · exact by simp [h a (by simp) hk b (by simp [hb]) hm ((hh b hb).2 hm)]
          · simp [hm]
        · simp [hk]
      · exact fun x hx hk b hb => h x (by simp [hx]) hk b (by simp [hb])

@[simp] theorem incident_lo (e : Edge) : incident e e.lo = true := by simp [incident]
@[simp] theorem incident_hi (e : Edge) : incident e e.hi = true := by simp [incident]

theorem eq_of_incident (a b : Edge) (hl : incident b a.lo = true)
    (hh : incident b a.hi = true) : a = b := by
  cases a with
  | mk x y hxy =>
    cases b with
    | mk u v huv =>
      simp [incident] at hl hh
      congr 1 <;> omega

/-- Adversarial load from two actual boundary edges. Their common endpoint
    gets two, the two other endpoints zero, everything else one. -/
def wedge (e f : Edge) : Source := fun p =>
  if incident e p && incident f p then .two
  else if incident e p || incident f p then .zero else .one

theorem wedge_symm (e f : Edge) : wedge e f = wedge f e := by
  funext p
  simp [wedge, Bool.and_comm, Bool.or_comm]

theorem wedge_first (e f : Edge) (hne : e ≠ f) (hm : meets e f = true) :
    initialReady (wedge e f) e = true := by
  by_cases hl : incident f e.lo = true
  · have hh : incident f e.hi = false :=
      Bool.eq_false_iff.mpr (fun h => hne (eq_of_incident e f hl h))
    simp [initialReady, wedge, hl, hh, gap]
  · have hh : incident f e.hi = true := by simpa [meets, hl] using hm
    simp [initialReady, wedge, hl, hh, gap]

theorem wedge_ready (e f : Edge) (hne : e ≠ f) (hm : meets e f = true) (a : Edge) :
    initialReady (wedge e f) a = true ↔ a = e ∨ a = f := by
  constructor
  · intro hr
    have he := eq_of_incident a e
    have hf := eq_of_incident a f
    by_cases el : incident e a.lo = true <;> by_cases eh : incident e a.hi = true <;>
      by_cases fl : incident f a.lo = true <;> by_cases fh : incident f a.hi = true <;>
      simp_all [initialReady, wedge, gap]
  · rintro (rfl | rfl)
    · exact wedge_first _ _ hne hm
    · rw [wedge_symm]
      exact wedge_first _ _ (Ne.symm hne) ((meets_symm ..).trans hm)

theorem wedge_accepted (rank : Edge → Nat) (screen : List Edge) (e f : Edge)
    (ho : Ordered rank screen) (hm : meets e f = true) (hp : rank f < rank e) :
    accepted (initialReady (wedge e f)) meets screen =
      (if f ∈ screen then [f] else if e ∈ screen then [e] else []) := by
  have hne : e ≠ f := by intro he; subst f; omega
  have hready := wedge_ready e f hne hm
  induction screen with
  | nil => simp [accepted]
  | cons a rest ih =>
    obtain ⟨hh, ht⟩ := List.pairwise_cons.mp ho
    rw [accepted, ih ht]
    by_cases hae : a = e
    · subst a
      have he : e ∉ rest := fun he => (hh e he).1 rfl
      by_cases hf : f ∈ rest <;> simp [hready, hf, he, hm, hne, Ne.symm hne]
    · by_cases haf : a = f
      · subst a
        have hf : f ∉ rest := fun hf => (hh f hf).1 rfl
        have he : e ∉ rest := by
          intro h
          have := (hh e h).2 ((meets_symm ..).trans hm)
          omega
        simp [hready, he, hf]
      · have hquiet : initialReady (wedge e f) a = false :=
          Bool.eq_false_iff.mpr (fun h => (hready a).mp h |>.elim hae haf)
        simp [hquiet, Ne.symm hae, Ne.symm haf]

def Agrees (rank : Edge → Nat) (screen : List Edge) (keep : Edge → Bool) :=
  ∀ (L : Source) (e : Edge), e ∈ screen → keep e = true →
    (sharpRun rank screen L).collar e = (sharpRun rank (screen.filter keep) L).collar e

theorem closed_agrees (rank : Edge → Nat) (screen : List Edge) (keep : Edge → Bool)
    (ho : Ordered rank screen) (hc : PredecessorClosed rank screen keep) :
    Agrees rank screen keep := by
  intro L e _ hk
  have hp := closed_projection keep (initialReady L) meets screen
    ((closed_iff_predecessor rank screen keep ho).mpr hc)
  rw [sharpRun_eq_material rank screen L ho,
      sharpRun_eq_material rank (screen.filter keep) L (List.Pairwise.filter keep ho), ← hp]
  simp [material, hk]

theorem boundary_witness (rank : Edge → Nat) (screen : List Edge) (keep : Edge → Bool)
    (ho : Ordered rank screen) {e f : Edge} (he : e ∈ screen) (hf : f ∈ screen)
    (hke : keep e = true) (hkf : keep f = false) (hm : meets e f = true)
    (hp : rank f < rank e) :
    (sharpRun rank screen (wedge e f)).collar e = none ∧
    (sharpRun rank screen (wedge e f)).collar f = some (sign (wedge e f) f) ∧
    (sharpRun rank (screen.filter keep) (wedge e f)).collar e =
      some (sign (wedge e f) e) := by
  have hne : e ≠ f := by intro h; subst f; omega
  have hpatch := List.Pairwise.filter keep ho
  rw [sharpRun_eq_material rank screen _ ho, sharpRun_eq_material rank _ _ hpatch,
      wedge_accepted rank screen e f ho hm hp,
      wedge_accepted rank (screen.filter keep) e f hpatch hm hp]
  simp [he, hf, hke, hkf, hne, material]

theorem predecessorClosed_iff_agrees (rank : Edge → Nat) (screen : List Edge)
    (keep : Edge → Bool) (ho : Ordered rank screen) :
    PredecessorClosed rank screen keep ↔ Agrees rank screen keep := by
  constructor
  · exact closed_agrees rank screen keep ho
  · intro h e he hke f hf hm hp
    cases hkf : keep f with
    | true => rfl
    | false =>
      obtain ⟨hfull, _, hpatch⟩ := boundary_witness rank screen keep ho he hf hke hkf hm hp
      have hh := h (wedge e f) e he hke
      rw [hfull, hpatch] at hh
      contradiction

def AcceptedAgrees (rank : Edge → Nat) (screen : List Edge) (keep : Edge → Bool) :=
  ∀ (L : Source) (e : Edge), e ∈ screen → keep e = true →
    ((sharpRun rank screen L).collar e).isSome =
      ((sharpRun rank (screen.filter keep) L).collar e).isSome

/-- The converse holds if the observer reads only an occupancy bit. -/
theorem predecessorClosed_iff_acceptedAgrees (rank : Edge → Nat) (screen : List Edge)
    (keep : Edge → Bool) (ho : Ordered rank screen) :
    PredecessorClosed rank screen keep ↔ AcceptedAgrees rank screen keep := by
  constructor
  · intro h L e he hk
    exact congrArg Option.isSome (closed_agrees rank screen keep ho h L e he hk)
  · intro h e he hke f hf hm hp
    cases hkf : keep f with
    | true => rfl
    | false =>
      obtain ⟨hfull, _, hpatch⟩ := boundary_witness rank screen keep ho he hf hke hkf hm hp
      have hh := h (wedge e f) e he hke
      simp [hfull, hpatch] at hh

theorem settled_attempt (rank : Edge → Nat) (screen : List Edge) (L : Source)
    (ho : Ordered rank screen) (e : Edge) (he : e ∈ screen) :
    attempt rank screen (sharpRun rank screen L) e = sharpRun rank screen L := by
  simp [attempt, enabled, sharpRun_quiescent rank screen L ho e he]

theorem ready_sign (L : Source) (e : Edge) (h : initialReady L e = true) :
    sign L e = 2 ∨ sign L e = -2 := by
  unfold initialReady at h
  cases hl : L e.lo <;> cases hh : L e.hi <;> simp_all [gap, sign, Load.value]

namespace FivePorts
def e : Edge := ⟨0, 1, by decide⟩
def f : Edge := ⟨1, 82, by decide⟩
def q : Edge := ⟨72, 73, by decide⟩
def r : Edge := ⟨72, 82, by decide⟩
def board : List Edge := [r, q, e, f]
def rank (b : Edge) : Nat := if b = f then 0 else if b = r then 4 else 1
def keepEF (b : Edge) : Bool := decide (b = e ∨ b = f)
def keepE (b : Edge) : Bool := decide (b = e)
theorem legal : Ordered rank board := by
  unfold Ordered board
  decide
theorem proper_closed : PredecessorClosed rank board keepEF :=
  (closed_iff_predecessor rank board keepEF legal).mp (by decide)
theorem open_wedge :
    (sharpRun rank board (wedge e f)).collar e = none ∧
    (sharpRun rank board (wedge e f)).collar f = some 2 ∧
    (sharpRun rank (board.filter keepE) (wedge e f)).collar e = some (-2) := by decide
end FivePorts

#print axioms write_energy
#print axioms sharpRun_eq_material
#print axioms sharpRun_quiescent
#print axioms closed_iff_predecessor
#print axioms wedge_ready
#print axioms wedge_accepted
#print axioms boundary_witness
#print axioms predecessorClosed_iff_agrees
#print axioms predecessorClosed_iff_acceptedAgrees
#print axioms settled_attempt
#print axioms ready_sign
#print axioms FivePorts.proper_closed
#print axioms FivePorts.open_wedge

end ObserverPatchSharp
