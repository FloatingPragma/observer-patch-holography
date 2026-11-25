import Computation.FixedFederationComplexity

/-! # A strict echo as the local instrument of a Boolean repair node

The executable integer matrix certificate checks which outgoing columns
light the balanced-pair pointer on sharp extremal rails. Each probe is
unwritten before the next. Here its two-query decoder is connected to
the formula compiler and its deterministic accepted step.
-/

namespace OPH.RepairUniversality.NativeBooleanNodes

def rail : Bool → Nat
  | false => 0
  | true => 2

/-- The registered pointer on sharp rails after the strict echo. -/
def balancedPointer (x y : Bool) : Bool :=
  decide (rail x + rail y = 2)

/-- The second shot, when needed, pairs the first operand with a known double. -/
def firstRead (x y : Bool) : Nat :=
  if balancedPointer x y then 1 else if balancedPointer x true then 0 else 2

theorem pointer_decoder_eq_half_sum (x y : Bool) :
    firstRead x y = (rail x + rail y) / 2 := by
  cases x <;> cases y <;> decide

inductive Service where
  | copy (i : Nat)
  | const (b : Bool)
  | nand (i j : Nat)

def Service.bit (op : Service) (s : State) : Bool :=
  match op with
  | .copy i => decide (firstRead (s i) true = 2)
  | .const b => decide (firstRead b true = 2)
  | .nand i j => !(decide (firstRead (s i) (s j) = 2))

@[simp] theorem copy_bit (x : Bool) : decide (firstRead x true = 2) = x := by
  cases x <;> decide

@[simp] theorem nand_bit (x y : Bool) :
    (!(decide (firstRead x y = 2))) = (!(x && y)) := by
  cases x <;> cases y <;> decide

/-- The second read, on a scratch rail and the old output, clocks the swap. -/
def attempt (reg : Nat) (op : Service) (s : State) : State :=
  if firstRead (op.bit s) (s reg) = 1 then write s reg (!(s reg)) else s

theorem attempt_eq_write (reg : Nat) (op : Service) (s : State) :
    attempt reg op s = write s reg (op.bit s) := by
  funext j
  by_cases hj : j = reg
  · subst j
    cases h : op.bit s <;> cases k : s reg <;> simp [attempt, h, k, write, firstRead, balancedPointer, rail]
  · unfold attempt
    split <;> simp [write, hj]

theorem emitted_has_native_service {k : Nat} (φ : Formula k) (c : Nat) :
    ∀ n ∈ (go φ c).nodes, ∃ op : Service, ∀ s : State, op.bit s = n.val s := by
  induction φ generalizing c with
  | var i =>
    intro n hn
    simp only [go, List.mem_singleton] at hn
    subst n
    exact ⟨.copy i.val, fun s => copy_bit (s i.val)⟩
  | const b =>
    intro n hn
    simp only [go, List.mem_singleton] at hn
    subst n
    exact ⟨.const b, fun _ => copy_bit b⟩
  | nand f g hf hg =>
    intro n hn
    simp only [go, List.mem_append, List.mem_singleton] at hn
    rcases hn with h | h | h
    · exact hf c n h
    · exact hg (go f c).next n h
    · subst n
      exact ⟨.nand (go f c).out (go g (go f c).next).out, fun s => nand_bit _ _⟩

theorem native_attempt_of_compiled_node {k : Nat} (φ : Formula k)
    (n : Node) (hn : n ∈ FixedFederation.fixedProgram φ) :
    ∃ op : Service, ∀ s, attempt n.reg op s = FixedFederation.repairNode n s := by
  obtain ⟨op, h⟩ := emitted_has_native_service φ k n hn
  refine ⟨op, fun s => ?_⟩
  change attempt n.reg op s = write s n.reg (n.val s)
  rw [attempt_eq_write, h s]

/-- The scheduled address is itself a member of the compiled list. -/
def serviceRun (L : List Node) (op : {n : Node // n ∈ L} → Service)
    (calendar : FixedFederation.NodeScheduler L) (s : State) : Nat → State
  | 0 => s
  | t + 1 =>
    let before := serviceRun L op calendar s t
    let node := calendar t before
    attempt node.val.reg (op node) before

/-- The one-time assignment of a service to each emitted node works for
all initial registers, adaptive address calendars and finite prefixes. -/
theorem compiled_histories {k : Nat} (φ : Formula k) :
    ∃ op : {n : Node // n ∈ FixedFederation.fixedProgram φ} → Service,
      ∀ calendar s t,
        serviceRun (FixedFederation.fixedProgram φ) op calendar s t =
          FixedFederation.attemptRun (FixedFederation.fixedProgram φ) t calendar s := by
  classical
  let op := fun n : {n : Node // n ∈ FixedFederation.fixedProgram φ} =>
    Classical.choose (native_attempt_of_compiled_node φ n.val n.property)
  refine ⟨op, ?_⟩
  intro calendar s t
  induction t with
  | zero => rfl
  | succ t ih =>
    simp only [serviceRun, FixedFederation.attemptRun, ih]
    let state := FixedFederation.attemptRun (FixedFederation.fixedProgram φ) t calendar s
    let node := calendar t state
    exact Classical.choose_spec (native_attempt_of_compiled_node φ node.val node.property) state

/-- The state-blind calendar consumed by the native services revisits every
member and bounds the delay to a genuine write from any inconsistent state. -/
theorem compiled_calendar {k : Nat} (φ : Formula k) (s : State) :
    FixedFederation.NodePathwiseWeakFair (FixedFederation.fixedProgram φ)
        (FixedFederation.fixedRoundRobinScheduler φ) s ∧
    FixedFederation.NodeBoundedWaste (FixedFederation.fixedProgram φ)
        ((FixedFederation.fixedProgram φ).length - 1)
        (FixedFederation.fixedRoundRobinScheduler φ) :=
  ⟨FixedFederation.fixedRoundRobin_pathwiseWeakFair φ s,
    FixedFederation.fixedRoundRobin_nodeBoundedWaste φ⟩

#print axioms pointer_decoder_eq_half_sum
#print axioms compiled_histories
#print axioms compiled_calendar
#print axioms native_attempt_of_compiled_node

end OPH.RepairUniversality.NativeBooleanNodes
