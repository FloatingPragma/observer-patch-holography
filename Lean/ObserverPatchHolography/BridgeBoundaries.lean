import Mathlib

/-!
# Bridge boundaries for repair confluence and union-collar gluing

This module records two finite counterexamples used by
`reality_as_consensus_protocol.tex`.

* Disjoint writes that separately preserve a nonlinear protected observable
  at one state need not compose to another protected state.
* All one-bit and two-bit marginals can agree while two tripartite laws differ.

The positive confluence and coherent-gluing statements therefore consume
separate receipts. No physical current, equilibrium law, or continuum claim
is formalized here.
-/

namespace OPH.BridgeBoundaries

/-! ## Protected-observable counterexample -/

abbrev RepairState := Bool × Bool

def origin : RepairState := (false, false)

def firstWrite (x : RepairState) : RepairState := (true, x.2)

def secondWrite (x : RepairState) : RepairState := (x.1, true)

def protectedObservable (x : RepairState) : Bool := x.1 && x.2

def repairMeasure : RepairState → Nat
  | (false, false) => 2
  | (true, false) => 1
  | (false, true) => 1
  | (true, true) => 0

theorem first_write_descends_at_origin :
    repairMeasure (firstWrite origin) < repairMeasure origin := by
  decide

theorem second_write_descends_at_origin :
    repairMeasure (secondWrite origin) < repairMeasure origin := by
  decide

theorem first_write_preserves_protected_at_origin :
    protectedObservable (firstWrite origin) = protectedObservable origin := rfl

theorem second_write_preserves_protected_at_origin :
    protectedObservable (secondWrite origin) = protectedObservable origin := rfl

theorem writes_are_coordinate_disjoint :
    (firstWrite origin).2 = origin.2 ∧
      (secondWrite origin).1 = origin.1 := by
  decide

theorem composed_writes_change_protected :
    protectedObservable (secondWrite (firstWrite origin)) ≠ protectedObservable origin := by
  decide

theorem peak_endpoints_differ :
    firstWrite origin ≠ secondWrite origin := by
  decide

/-! ## Pair-marginal counterexample -/

structure Bit3Law where
  p000 : ℚ
  p001 : ℚ
  p010 : ℚ
  p011 : ℚ
  p100 : ℚ
  p101 : ℚ
  p110 : ℚ
  p111 : ℚ
  deriving DecidableEq

def evenParity : Bit3Law :=
  ⟨1 / 4, 0, 0, 1 / 4, 0, 1 / 4, 1 / 4, 0⟩

def oddParity : Bit3Law :=
  ⟨0, 1 / 4, 1 / 4, 0, 1 / 4, 0, 0, 1 / 4⟩

def marginalAB (p : Bit3Law) : Bool × Bool → ℚ
  | (false, false) => p.p000 + p.p001
  | (false, true) => p.p010 + p.p011
  | (true, false) => p.p100 + p.p101
  | (true, true) => p.p110 + p.p111

def marginalAC (p : Bit3Law) : Bool × Bool → ℚ
  | (false, false) => p.p000 + p.p010
  | (false, true) => p.p001 + p.p011
  | (true, false) => p.p100 + p.p110
  | (true, true) => p.p101 + p.p111

def marginalBC (p : Bit3Law) : Bool × Bool → ℚ
  | (false, false) => p.p000 + p.p100
  | (false, true) => p.p001 + p.p101
  | (true, false) => p.p010 + p.p110
  | (true, true) => p.p011 + p.p111

theorem even_odd_are_distinct : evenParity ≠ oddParity := by
  intro h
  have h000 := congrArg Bit3Law.p000 h
  norm_num [evenParity, oddParity] at h000

theorem even_odd_marginalAB : marginalAB evenParity = marginalAB oddParity := by
  funext x
  rcases x with ⟨a, b⟩
  cases a <;> cases b <;> norm_num [marginalAB, evenParity, oddParity]

theorem even_odd_marginalAC : marginalAC evenParity = marginalAC oddParity := by
  funext x
  rcases x with ⟨a, c⟩
  cases a <;> cases c <;> norm_num [marginalAC, evenParity, oddParity]

theorem even_odd_marginalBC : marginalBC evenParity = marginalBC oddParity := by
  funext x
  rcases x with ⟨b, c⟩
  cases b <;> cases c <;> norm_num [marginalBC, evenParity, oddParity]

theorem pair_marginals_do_not_determine_global_law :
    evenParity ≠ oddParity ∧
      marginalAB evenParity = marginalAB oddParity ∧
      marginalAC evenParity = marginalAC oddParity ∧
      marginalBC evenParity = marginalBC oddParity :=
  ⟨even_odd_are_distinct, even_odd_marginalAB,
    even_odd_marginalAC, even_odd_marginalBC⟩

#print axioms composed_writes_change_protected
#print axioms pair_marginals_do_not_determine_global_law

end OPH.BridgeBoundaries
