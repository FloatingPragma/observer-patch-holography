import Mathlib

/-!
# OPH Primitives (placeholders, sorry-bearing)

These are the primitives Proposition 4.2 depends on. They are **not**
formalised here — they are TeX macros in *Paradise as Fixed-Point
Consensus* (lines 28–31, 285–303) whose structural content is outsourced
to the companion paper *Reality as a Consensus Protocol* (cited as
`OPHConsensus`; reference at *Paradise* line 1615).

This file makes the cross-paper dependency visible at the **type level**:
any downstream definition that wants to claim Prop-4.2-relevance must
instantiate these primitives, not paper over them. The admitted declarations
below are `sorry`-bearing on purpose: `lake build` warns on every one, and
CI checks that the count stays fixed until the primitives are discharged.

## Filling in (from the paper)

* `Records`, `Patch`, `Obs` — line 28–31 TeX macros; structural content in
  OPHConsensus.
* `Repair : Records → Records` — line 30 macro, "built from local recovery
  moves" (line 297). The companion consensus paper expands this into
  local accepted repair steps executed under asynchronous schedules.
* `Φ : Records → NNReal` — line 300 concrete formula:
  `Φ(x) = Σ_e w_e · d_e(π_{i,e}(x_i), π_{j,e}(x_j))`.
* `gaugeEquiv` (`∼_gauge`) — line 311: identifies hidden local
  presentations with the same declared observable overlap data.
* `repair_respects_gauge` — `∼_gauge` is a `Repair`-congruence; this is
  the load-bearing obligation Prop 4.2 sentence 2 ("on the physical
  quotient") imposes.
* `acceptedStep`, `LyapunovDescent`, `Termination`, `Confluence`, and
  `Completeness` — the asynchronous-repair obligations from OPHConsensus:
  descent gives termination on the finite patch-net branch; confluence gives
  schedule independence; completeness says normal forms are exactly
  consistent states.
-/

namespace OPH

def Records : Type := sorry
def Patch : Type := sorry
def Obs : Type := sorry
def Site : Type := sorry

noncomputable def Repair : Records → Records := sorry

noncomputable def localRepair : Site → Records → Records := sorry

/-- One accepted asynchronous repair step. This is the relation the generic
    abstract-rewriting skeleton must eventually instantiate. -/
def acceptedStep (x y : Records) : Prop :=
  ∃ i : Site, y = localRepair i x ∧ localRepair i x ≠ x

noncomputable def Φ : Records → NNReal := sorry

def NormalForm (x : Records) : Prop :=
  ∀ y : Records, ¬ acceptedStep x y

def Consistent (x : Records) : Prop :=
  Φ x = 0

def LyapunovDescent : Prop :=
  ∀ x y : Records, acceptedStep x y → Φ y < Φ x

def Termination : Prop :=
  WellFounded (fun y x : Records => acceptedStep x y)

def gaugeEquiv : Records → Records → Prop := sorry

theorem gaugeEquiv_equivalence : Equivalence gaugeEquiv := sorry

/-- `∼_gauge` is a `Repair`-congruence. Required by Prop 4.2 sentence 2
    (independence on the physical quotient). -/
theorem repair_respects_gauge :
    ∀ x y : Records, gaugeEquiv x y → gaugeEquiv (Repair x) (Repair y) :=
  sorry

/-- OPH confluence condition for accepted asynchronous repair steps
    (Prop 4.2 hypothesis; defined per OPHConsensus). -/
def Confluence : Prop :=
  ∀ x y z : Records, ReflTransGen acceptedStep x y → ReflTransGen acceptedStep x z →
    ∃ w : Records, ReflTransGen acceptedStep y w ∧ ReflTransGen acceptedStep z w

/-- OPH repair completeness: normal forms are exactly consistent states.
    Termination is a separate Lyapunov/finite-state obligation. -/
def Completeness : Prop :=
  ∀ x : Records, NormalForm x ↔ Consistent x

end OPH
