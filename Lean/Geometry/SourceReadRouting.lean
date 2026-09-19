import Std

/-!
# Full-family compilation under the retained M1 transport law

The generic compiler takes a transport identity as a hypothesis.  The concrete
six-event integer word below discharges it, including scratch cleanup and
preservation of the protected source.  The real-valued error bound is proved
separately in `Geometry.SourceFeedbackTransport`. No theorem selects that feedback
law, the read menu, the host placement, or a physical clock from A1--A3.

The finite theorem applies to every number of layers, every immutable initial
version, and every integer local intervention.  Its order theorem consumes
edge and read-path certificates, checked against every retained primitive in
the independent full-family replay.  Custody/control order is not erased or
identified with semantic value-consumption order.
-/

namespace OPH.SourceReadRouting

variable {Site : Type}

def sumReads (menu : List Site) (values : Site → Int) : Int :=
  (menu.map values).foldl (· + ·) 0

def logicalHistory (menu : Site → List Site) (bias : Nat → Site → Int)
    (initial : Site → Int) : Nat → Site → Int
  | 0 => initial
  | t + 1 => fun s => bias t s + sumReads (menu s) (logicalHistory menu bias initial t)

def compiledHistory (transport : Nat → Int → Int)
    (depth : Nat → Site → Site → Nat) (menu : Site → List Site)
    (bias : Nat → Site → Int) (initial : Site → Int) : Nat → Site → Int
  | 0 => initial
  | t + 1 => fun s => bias t s + sumReads (menu s)
      (fun u => transport (depth t u s) (compiledHistory transport depth menu bias initial t u))

/-- All initial-value and local-commit interventions are preserved.  There
is no inference from a finite list of sampled perturbations. -/
theorem compiledHistory_eq (transport : Nat → Int → Int)
    (exactTransport : ∀ d x, transport d x = x)
    (depth : Nat → Site → Site → Nat) (menu : Site → List Site)
    (bias : Nat → Site → Int) (initial : Site → Int) (t : Nat) :
    compiledHistory transport depth menu bias initial t = logicalHistory menu bias initial t := by
  induction t with
  | zero => rfl
  | succ t ih =>
    funext s
    simp only [compiledHistory, logicalHistory, exactTransport, ih]

/-- Register allocation, routing and schedules may vary, provided their
certified transports implement the same immutable value interface. -/
theorem schedule_independent_readouts (first second : Nat → Int → Int)
    (hf : ∀ d x, first d x = x) (hs : ∀ d x, second d x = x)
    (df ds : Nat → Site → Site → Nat) (menu : Site → List Site)
    (bias : Nat → Site → Int) (initial : Site → Int) (t : Nat) :
    compiledHistory first df menu bias initial t = compiledHistory second ds menu bias initial t := by
  rw [compiledHistory_eq first hf, compiledHistory_eq second hs]

def executeWord {State : Type} : List (State → State) → State → State
  | [] => id
  | step :: rest => fun x => executeWord rest (step x)

/-- Values in half-units, as in the native event rows. `receiverArchive`
models the newly allocated capture slot, not an overwrite of an old version.
The constant zero records and the admitted local feedback are supplied M1.
This is unbounded integer arithmetic, not a claim about arbitrary Int64 inputs. -/
structure HopState where
  sourceArchive : Int
  sourcePort : Int
  receiverPort : Int
  receiverArchive : Int
  deriving DecidableEq

inductive HopOp where
  | exportPayload | resetReceiver | pairMean | capture | resetSource
  deriving DecidableEq

def hopStep (op : HopOp) (s : HopState) : HopState :=
  match op with
  | .exportPayload => { s with sourcePort := s.sourceArchive }
  | .resetReceiver => { s with receiverPort := 0 }
  | .pairMean =>
      let mean := (s.sourcePort + s.receiverPort) / 2
      { s with sourcePort := mean, receiverPort := mean }
  | .capture => { s with receiverArchive := 2 * s.receiverPort }
  | .resetSource => { s with sourcePort := 0 }

def sixEventOps : List HopOp :=
  [.exportPayload, .resetReceiver, .pairMean, .capture, .resetSource, .resetReceiver]

def sixEventHop : HopState → HopState := executeWord (sixEventOps.map hopStep)

/-- The actual six operations transmit any integer payload and erase arbitrary
incoming scratch loads. Both the source archive and the captured version equal
the original payload; both ports are zero on exit. -/
theorem six_event_hop_exact (x a b old : Int) :
    sixEventHop ⟨2*x, a, b, old⟩ = ⟨2*x, 0, 0, 2*x⟩ := by
  simp [sixEventHop, sixEventOps, executeWord, hopStep]

def hopReads : HopOp → Nat
  | .pairMean => 2
  | _ => 1

def hopWrites : HopOp → Nat
  | .pairMean => 2
  | _ => 1

theorem six_event_hop_cost :
    sixEventOps.length = 6 ∧
    (sixEventOps.map hopReads).sum = 7 ∧ (sixEventOps.map hopWrites).sum = 7 := by
  decide

def feedbackHopValue (x : Int) : Int :=
  (sixEventHop ⟨2*x, 0, 0, 0⟩).receiverArchive / 2

theorem feedbackHopValue_eq (x : Int) : feedbackHopValue x = x := by
  simp [feedbackHopValue, six_event_hop_exact]

def feedbackRoute : Nat → Int → Int
  | 0, x => x
  | d + 1, x => feedbackHopValue (feedbackRoute d x)

theorem feedbackRoute_eq (d : Nat) (x : Int) : feedbackRoute d x = x := by
  induction d with
  | zero => rfl
  | succ d ih => simp [feedbackRoute, feedbackHopValue_eq, ih]

/-- The concrete feedback word discharges the generic compiler's transport
hypothesis. Fixed routing and local-commit interventions are arbitrary here;
the graph/writer-path certificate is separately checked by native replay. -/
theorem six_event_compiledHistory_eq
    (depth : Nat → Site → Site → Nat) (menu : Site → List Site)
    (bias : Nat → Site → Int) (initial : Site → Int) (t : Nat) :
    compiledHistory feedbackRoute depth menu bias initial t =
      logicalHistory menu bias initial t :=
  compiledHistory_eq feedbackRoute feedbackRoute_eq depth menu bias initial t

/-- Class-wide finite-word invariant. Applied analytically to total live
load, every canonical pair mean satisfies the hypothesis. -/
theorem word_preserves_invariant {State Value : Type}
    (invariant : State → Value) (word : List (State → State))
    (preserves : ∀ step ∈ word, ∀ x, invariant (step x) = invariant x) (x : State) :
    invariant (executeWord word x) = invariant x := by
  induction word generalizing x with
  | nil => rfl
  | cons step rest ih =>
    have tail : ∀ f ∈ rest, ∀ y, invariant (f y) = invariant y := by
      intro f member y
      exact preserves f (List.mem_cons_of_mem step member) y
    exact (ih tail (step x)).trans (preserves step (List.mem_cons_self) x)

/-- A reset changing total load cannot be an all-registers-closed word of
sum-preserving operations. This does not exclude routes with a changed
ancilla, protected readback, feedback, or a different operational projection. -/
theorem closed_sum_preserving_word_cannot_reset
    (word : List ((Int × Int) → (Int × Int)))
    (preserves : ∀ step ∈ word, ∀ x, (step x).1 + (step x).2 = x.1 + x.2) :
    ¬ (∀ x, executeWord word x = (0, x.2)) := by
  intro reset
  have kept := word_preserves_invariant (fun x : Int × Int => x.1+x.2) word preserves (1,0)
  rw [reset] at kept
  simp at kept

inductive Reach {α : Type} (edge : α → α → Prop) : α → α → Prop where
  | refl (a) : Reach edge a a
  | step {a b c} : edge a b → Reach edge b c → Reach edge a c

theorem reach_trans {α : Type} {edge : α → α → Prop} {a b c : α}
    (first : Reach edge a b) (second : Reach edge b c) : Reach edge a c := by
  induction first with
  | refl => exact second
  | step h _ ih => exact Reach.step h (ih second)

def Below {α : Type} (edge : α → α → Prop) : Option α → Option α → Prop
  | none, _ => True
  | some _, none => False
  | some a, some b => Reach edge a b

theorem below_trans {α : Type} {edge : α → α → Prop} {a b c : Option α}
    (ab : Below edge a b) (bc : Below edge b c) : Below edge a c := by
  cases a with
  | none => trivial
  | some a =>
    cases b with
    | none => exact False.elim ab
    | some b =>
      cases c with
      | none => exact False.elim bc
      | some c => exact reach_trans ab bc

/-- Every semantic edge is checked, including initialization, resets and
intermediate relays.  `none` means no logical-version payload was consumed. -/
theorem projection_sound {Logical Event : Type}
    (logical : Logical → Logical → Prop) (operational : Event → Event → Prop)
    (label : Event → Option Logical)
    (edges : ∀ a b, operational a b → Below logical (label a) (label b))
    {a b : Event} (path : Reach operational a b) : Below logical (label a) (label b) := by
  induction path with
  | refl a =>
    cases label a with
    | none => trivial
    | some x => exact Reach.refl x
  | step h _ ih => exact below_trans (edges _ _ h) ih

/-- The independent replay checks a consumed writer chain for each requested
read.  This lifts the entire generated order, not just isolated pairs. -/
theorem projection_complete {Logical Event : Type}
    (logical : Logical → Logical → Prop) (operational : Event → Event → Prop)
    (embed : Logical → Event)
    (reads : ∀ a b, logical a b → Reach operational (embed a) (embed b))
    {a b : Logical} (path : Reach logical a b) : Reach operational (embed a) (embed b) := by
  induction path with
  | refl a => exact Reach.refl (embed a)
  | step h _ ih => exact reach_trans (reads _ _ h) ih

theorem exact_induced_order {Logical Event : Type}
    (logical : Logical → Logical → Prop) (operational : Event → Event → Prop)
    (embed : Logical → Event) (label : Event → Option Logical)
    (hsection : ∀ a, label (embed a) = some a)
    (edges : ∀ a b, operational a b → Below logical (label a) (label b))
    (reads : ∀ a b, logical a b → Reach operational (embed a) (embed b))
    (a b : Logical) :
    Reach operational (embed a) (embed b) ↔ Reach logical a b := by
  constructor
  · intro path
    have result := projection_sound logical operational label edges path
    simpa only [hsection, Below] using result
  · exact projection_complete logical operational embed reads

end OPH.SourceReadRouting

#print axioms OPH.SourceReadRouting.compiledHistory_eq
#print axioms OPH.SourceReadRouting.six_event_hop_exact
#print axioms OPH.SourceReadRouting.six_event_hop_cost
#print axioms OPH.SourceReadRouting.six_event_compiledHistory_eq
#print axioms OPH.SourceReadRouting.schedule_independent_readouts
#print axioms OPH.SourceReadRouting.exact_induced_order
#print axioms OPH.SourceReadRouting.closed_sum_preserving_word_cannot_reset
