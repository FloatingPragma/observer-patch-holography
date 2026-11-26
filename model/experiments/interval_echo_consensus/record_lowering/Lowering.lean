import Geometry.SourceRecordGluing
import Mathlib.Tactic

namespace OPH.NativeRecordLowering

/-- Native detector against a zero standard on the material mean of two extremal rails. -/
def load (b : Bool) : Nat := if b then 2 else 0

def hit (a b : Bool) : Bool := decide ((load a + load b) / 2 = 2)

theorem hit_eq (a b : Bool) : hit a b = (a && b) := by
  cases a <;> cases b <;> decide

inductive Op | push | pop | toggle0 | toggle1 | toggle2
  deriving DecidableEq

/-- Live controls and the stack are views of the fixed address calendar.
Pop reverses the physical scratch exchange and releases its clean dimer. -/
structure State where
  live : List Bool
  saved : List (Bool × Bool)
  free : List Bool
  out : Bool
  deriving DecidableEq

def tick : Op → State → State
  | .push, ⟨a::b::xs, stack, t::free, y⟩ =>
    ⟨(t ^^ hit a b)::xs, (a,b)::stack, free, y⟩
  | .pop, ⟨t::xs, (a,b)::stack, free, y⟩ =>
    ⟨a::b::xs, stack, (t ^^ hit a b)::free, y⟩
  | .toggle0, q => {q with out := !q.out}
  | .toggle1, ⟨[a], stack, free, y⟩ => ⟨[a], stack, free, y ^^ a⟩
  | .toggle2, ⟨[a,b], stack, free, y⟩ => ⟨[a,b], stack, free, y ^^ hit a b⟩
  | _, q => q

def compile : Nat → List Op
  | 0 => [.toggle0]
  | 1 => [.toggle1]
  | 2 => [.toggle2]
  | n+3 => .push :: compile (n+2) ++ [.pop]

abbrev run {I S : Type*} (step : I → S → S) := SourceRecordGluing.run step

def clean (n : Nat) := List.replicate (n - 2) false

/-- Every arbitrary-arity minterm uses only the calibrated shutter, and
returns all auxiliary dimers in their original state. The unused tail of
scratch, previous enclosing gates and the input literals are arbitrary. -/
theorem minterm (xs : List Bool) (stack : List (Bool × Bool))
    (free : List Bool) (y : Bool) :
    run tick (compile xs.length) ⟨xs, stack, clean xs.length ++ free, y⟩ =
      ⟨xs, stack, clean xs.length ++ free, y ^^ xs.all id⟩ := by
  match xs with
  | [] => cases y <;> rfl
  | [a] => cases a <;> cases y <;> rfl
  | [a,b] => cases a <;> cases b <;> cases y <;> rfl
  | a::b::c::rest =>
    have inner :=  minterm (hit a b :: c :: rest) ((a,b)::stack) free y
    change run tick (.push :: compile (rest.length + 2) ++ [.pop])
        ⟨_, _, false :: (clean (rest.length+2) ++ free), _⟩ = _
    simp only [run] at inner ⊢
    rw [SourceRecordGluing.run_append]
    simp only [SourceRecordGluing.run, tick, Bool.false_xor]
    rw [show (hit a b :: c :: rest).length = rest.length+2 by rfl] at inner
    rw [inner]
    simp [clean, List.replicate_succ, hit_eq,
      Bool.and_assoc]
termination_by xs.length

/-- For fixed public input, this tape is the concrete minterm routine. -/
def nativeRow (xs : List Bool) (out : Bool) :=
  (run tick (compile xs.length) ⟨xs, [], clean xs.length, out⟩).out

theorem nativeRow_eq (xs : List Bool) (y : Bool) : nativeRow xs y = (y ^^ xs.all id) := by
  simpa [nativeRow] using congrArg State.out (minterm xs [] [] y)

def equalLiterals : List Bool → List Bool → List Bool
  | [], [] => []
  | a::as, b::bs => (a == b) :: equalLiterals as bs
  | _, _ => [false]

theorem all_equal (a b : List Bool) : (equalLiterals a b).all id = (a == b) := by
  induction a generalizing b with
  | nil => cases b <;> rfl
  | cons a as ih => cases b with
    | nil => rfl
    | cons b bs => simp [equalLiterals, ih, List.cons_beq_cons]

/-- A finite law is a list of exactly those source rows on which one
published bit is true. Each row is compiled before the source is prepared. -/
def rows (yes : List (List Bool)) (x : List Bool) (y : Bool) : Bool :=
  SourceRecordGluing.run (fun r b => nativeRow (equalLiterals r x) b) yes y

theorem absent {yes : List (List Bool)} (x : List Bool) (h : x ∉ yes) (y : Bool) :
    rows yes x y = y := by
  induction yes generalizing y with
  | nil => rfl
  | cons r rs ih =>
    have neq : r ≠ x := fun he => h (by simp [he])
    have tail : x ∉ rs := fun hx => h (by simp [hx])
    change rows rs x (nativeRow (equalLiterals r x) y) = y
    simp [nativeRow_eq, all_equal, beq_eq_false_iff_ne.mpr neq, ih tail]

theorem finite_map (yes : List (List Bool)) (distinct : yes.Nodup) (x : List Bool) (y : Bool) :
    rows yes x y = (y ^^ decide (x ∈ yes)) := by
  induction yes generalizing y with
  | nil => simp [rows, SourceRecordGluing.run]
  | cons r rs ih =>
    obtain ⟨fresh, nd⟩ := List.nodup_cons.mp distinct
    change rows rs x (nativeRow (equalLiterals r x) y) = _
    by_cases same : r = x
    · subst r
      simp [nativeRow_eq, all_equal, absent x fresh]
    · have rev := Ne.symm same
      simp [nativeRow_eq, all_equal, beq_eq_false_iff_ne.mpr same, ih nd, rev]

/-- A replacement of the current draft is itself a named instruction.
Previously sealed registers remain part of the native state at every prefix. -/
inductive Address | source (n : Nat) | version (n : Nat)
  deriving DecidableEq

inductive Inst
  | map (yes : List (List Bool)) (operands : List Address)
  | intervene (b : Bool)
  | seal
  | wait

structure Published where
  sources : List Bool
  draft : Bool
  versions : List Bool
  deriving DecidableEq

def read (p : Published) : Address → Bool
  | .source i => p.sources[i]?.getD false
  | .version i => p.versions[i]?.getD false

def native : Inst → Published → Published
  | .map yes operands, p => {p with draft := rows yes (operands.map (read p)) p.draft}
  | .intervene b, p => {p with draft := p.draft ^^ b}
  | .seal, p => {p with versions := p.versions ++ [nativeRow [p.draft] false]}
  | .wait, p => p

def ideal : Inst → Published → Published
  | .map yes operands, p => {p with draft := p.draft ^^ decide (operands.map (read p) ∈ yes)}
  | .intervene b, p => {p with draft := p.draft ^^ b}
  | .seal, p => {p with versions := p.versions ++ [p.draft]}
  | .wait, p => p

theorem isolated_intervention (b : Bool) (p : Published) :
    (native (.intervene b) p).sources = p.sources ∧
    (native (.intervene b) p).versions = p.versions := ⟨rfl,rfl⟩

def Valid : Inst → Prop
  | .map yes operands => yes.Nodup ∧ operands.Nodup
  | _ => True

theorem prefix_faithful (w : List Inst) (hv : ∀ i ∈ w, Valid i) (p : Published) :
    run native w p = run ideal w p := by
  induction w generalizing p with
  | nil => rfl
  | cons a w ih =>
    have head := hv a (by simp)
    have tail : ∀ i ∈ w, Valid i := fun i h => hv i (by simp [h])
    have eq : native a p = ideal a p := by
      cases a <;> simp_all [Valid, native, ideal, finite_map, nativeRow_eq]
    change run native w (native a p) = run ideal w (ideal a p)
    rw [eq, ih tail]

theorem separate_contact {A B : Type*} (f : A → A) (g : B → B) (p : A × B) :
    Prod.map f id (Prod.map id g p) = Prod.map id g (Prod.map f id p) := rfl

/-- Apply independent instructions in the same appointment, explicitly retaining
padding on either branch. Both native words can be read back from the joint calendar. -/
def joint {A B : Type*} : List A → List B → List (Option A × Option B)
  | [], b => b.map fun y => (none, some y)
  | a::as, [] => (a::as).map fun x => (some x, none)
  | a::as, b::bs => (some a, some b) :: joint as bs

theorem shared_elapsed {A B : Type*} (a : List A) (b : List B) :
    (joint a b).length = max a.length b.length := by
  induction a generalizing b with
  | nil => simp [joint]
  | cons x xs ih => cases b with
    | nil => simp [joint]
    | cons y ys => simp [joint, ih, Nat.add_max_add_right]

theorem retained_words {A B : Type*} (a : List A) (b : List B) :
    (joint a b).filterMap Prod.fst = a ∧ (joint a b).filterMap Prod.snd = b := by
  induction a generalizing b with
  | nil => simp [joint, List.filterMap_map]
  | cons x xs ih => cases b with
    | nil => simp [joint, List.filterMap_map]
    | cons y ys => simpa [joint] using ih ys

end OPH.NativeRecordLowering
