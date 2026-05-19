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
instantiate these primitives, not paper over them. Every signature below
is `sorry`-bearing on purpose: `lake build` warns on every one, so a
stale primitive cannot hide behind a green CI run.

## Filling in (from the paper)

* `Records`, `Patch`, `Obs` — line 28–31 TeX macros; structural content in
  OPHConsensus.
* `Repair : Records → Records` — line 30 macro, "built from local recovery
  moves" (line 297). Paper says **local**, not *asynchronous*.
* `Φ : Records → NNReal` — line 300 concrete formula:
  `Φ(x) = Σ_e w_e · d_e(π_{i,e}(x_i), π_{j,e}(x_j))`.
* `gaugeEquiv` (`∼_gauge`) — line 311: identifies hidden local
  presentations with the same declared observable overlap data.
* `repair_respects_gauge` — `∼_gauge` is a `Repair`-congruence; this is
  the load-bearing obligation Prop 4.2 sentence 2 ("on the physical
  quotient") imposes.
* `Confluence`, `Completeness` predicates — referenced in Prop 4.2
  hypothesis (line 326) but not structurally defined in *Paradise*; see
  OPHConsensus.
-/

namespace OPH

def Records : Type := sorry
def Patch : Type := sorry
def Obs : Type := sorry

noncomputable def Repair : Records → Records := sorry

noncomputable def Φ : Records → NNReal := sorry

def gaugeEquiv : Records → Records → Prop := sorry

theorem gaugeEquiv_equivalence : Equivalence gaugeEquiv := sorry

/-- `∼_gauge` is a `Repair`-congruence. Required by Prop 4.2 sentence 2
    (independence on the physical quotient). -/
theorem repair_respects_gauge :
    ∀ x y : Records, gaugeEquiv x y → gaugeEquiv (Repair x) (Repair y) :=
  sorry

/-- OPH confluence condition (Prop 4.2 hypothesis; defined per
    OPHConsensus). -/
def Confluence : Prop := sorry

/-- OPH completeness condition (Prop 4.2 hypothesis; defined per
    OPHConsensus). -/
def Completeness : Prop := sorry

end OPH
