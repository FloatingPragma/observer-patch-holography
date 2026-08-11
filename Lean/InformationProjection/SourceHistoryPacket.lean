import InformationProjection.HistoryLaw
import MixingChainRealization

namespace OPH.InformationProjection

open OPH.Thermodynamics

/-!
# Source-attached history-law packet (B7, issue #683)

This module mirrors, as exact literals, the history-law packet extracted
from the locally hash-pinned bounded source run `runs/b12_prereg_16k_20260806`
of the simulator (seed `20260806`, engine commit
`b39b78faf894894ebe573571e0902ccfaaeac32a`).  The payload with the full
extraction and its exact verification lives at
`oph-physics-sim/docs/B7_HISTORY_PACKET_PAYLOAD.json`
(`scripts/extract_b7_history_packet.py`), computed over the committed
mixing-chain payload `docs/B12_MIXING_CHAIN_PAYLOAD.json`
(sha256 `01111b9d2338296ba4e181ec1875c1dfbbfd7d8ce1d76775ddbfba1604a76cb2`)
and the pinned run inputs `observer_views.jsonl`
(sha256 `3da6f04d770b81b49a082ad1e6ecb93c611517dd85ecaccd7cc89077fcd6c0ca`)
and `finite_repair_transition_matrix_report.json`
(sha256 `e38bb28475fb7b34127b864abd20414efc920d8ef9ddbdf5c78612874ae37f2d`).

The path space is the eight length-3 paths (three states, two
transitions) over the restricted two-state recurrent class `{18, 19}` of
that run, encoded `g = 4*s0 + 2*s1 + s2`.  Each maximal in-class segment
of an observer transition history contributes every consecutive state
triple as one window; the 1754 windows give the empirical path law by
exact normalization.  The reference path law is the stationary Markov
path law `pi(s0) * P(s0,s1) * P(s1,s2)` of the committed chain of
`MixingChainRealization`, and the action is the repair-move count of the
path, computed from the path alone.

Receipt structure.  The integer layer (window counts and totals, the
path encoding, the action as the state-change count, the nonconstant
action, and the cleared-denominator Markov-product identity of the
reference law) is kernel-decided on the literals.  The rational layer
(both laws with their positivity, normalization, and exact mean
actions, and the Markov-product display) is elaborated from the same
literals by `norm_num`; rational division is kept out of kernel
reduction.  The real layer casts the same literals and instantiates the
`HistoryLaw` interface at them: the tilt of the reference law by the
action is strictly positive and normalized at every multiplier, it
equals the reference exactly at multiplier zero, and at every multiplier
it is the unique relative-entropy minimizer on its mean-action surface.
An intermediate-value receipt locates a multiplier matching the declared
empirical mean action inside the rational bracket
`(31901/10000, 3190101/1000000)` carried by the payload, and the
declared multiplier point `lambda0 = -log (3190100906251/10^12)` has an
exactly rational tilted law whose mean action lies within `10^-12` of
the declared mean.

Claim boundary.  Representation-level mirror over literals extracted
from one retained locally hash-pinned bounded source run.  The declared
postprocessing constructs the Markov reference, repair-count action, and
empirical path law from those data; these objects were not all fixed by
the run's local pre-execution contract.  The multiplier is matched to the
DECLARED alternative mean, the empirical mean action of the same run.
No source mechanism selects the multiplier: the multiplier selection
principle from source data is the open core of issue B7.
-/

section TiltZero

variable {Γ : Type*} [Fintype Γ]

/-- At multiplier zero the exponential tilt of a normalized reference
law is the reference law itself. -/
theorem tilt_zero (τ S : Γ → ℝ) (hτ1 : ∑ g, τ g = 1) :
    tilt τ S 0 = τ := by
  funext g
  unfold tilt tiltZ
  simp only [zero_mul, neg_zero, Real.exp_zero, mul_one]
  rw [hτ1, div_one]

end TiltZero

/-! ## The literal layer: extracted integers and rationals -/

/-- First state of path `g` under the encoding `g = 4*s0 + 2*s1 + s2`
(local `0` = global 18, local `1` = global 19). -/
def sourceState0 : Fin 8 → Fin 2 := ![0, 0, 0, 0, 1, 1, 1, 1]

/-- Middle state of path `g`. -/
def sourceState1 : Fin 8 → Fin 2 := ![0, 0, 1, 1, 0, 0, 1, 1]

/-- Last state of path `g`. -/
def sourceState2 : Fin 8 → Fin 2 := ![0, 1, 0, 1, 0, 1, 0, 1]

/-- Integer window counts of the eight paths: every maximal in-class
segment of an observer transition history contributes each consecutive
state triple as one window. -/
def sourceWindowCount : Fin 8 → ℕ := ![1149, 96, 2, 87, 3, 2, 3, 412]

/-- The declared repair-count action: the number of state
changes) of the length-3 path. -/
def sourceAction : Fin 8 → ℕ := ![0, 1, 2, 1, 1, 2, 1, 0]

/-- The empirical path law: the exact normalization of the integer
window counts. -/
def sourceTauEmpQ : Fin 8 → ℚ :=
  ![1149/1754, 48/877, 1/877, 87/1754, 3/1754, 1/877, 3/1754, 206/877]

/-- The reference path law: the stationary Markov path law of the
committed mixing chain. -/
def sourceTauChainQ : Fin 8 → ℚ :=
  ![8764880/88022241, 708340/88022241, 2675/31247588, 269105/31247588,
    708340/88022241, 57245/88022241, 269105/31247588, 27071963/31247588]

/-- The path encoding is `g = 4*s0 + 2*s1 + s2`. -/
theorem sourceState_encoding :
    ∀ g : Fin 8, 4 * (sourceState0 g).val + 2 * (sourceState1 g).val
      + (sourceState2 g).val = g.val := by decide

/-- The action is the state-change count of the path. -/
theorem sourceAction_eq_changes :
    ∀ g : Fin 8, sourceAction g
      = (if sourceState0 g = sourceState1 g then 0 else 1)
        + (if sourceState1 g = sourceState2 g then 0 else 1) := by decide

/-- Every length-3 path over the recurrent class is realized: all eight
window counts are strictly positive. -/
theorem sourceWindowCount_pos : ∀ g, 0 < sourceWindowCount g := by decide

/-- The window counts total 1754. -/
theorem sourceWindowCount_total :
    sourceWindowCount 0 + sourceWindowCount 1 + sourceWindowCount 2
      + sourceWindowCount 3 + sourceWindowCount 4 + sourceWindowCount 5
      + sourceWindowCount 6 + sourceWindowCount 7 = 1754 := by decide

/-- The count-weighted action total is 197: the numerator of the
empirical mean action `197/1754`. -/
theorem sourceWindowCount_action_total :
    sourceWindowCount 0 * sourceAction 0 + sourceWindowCount 1 * sourceAction 1
      + sourceWindowCount 2 * sourceAction 2 + sourceWindowCount 3 * sourceAction 3
      + sourceWindowCount 4 * sourceAction 4 + sourceWindowCount 5 * sourceAction 5
      + sourceWindowCount 6 * sourceAction 6 + sourceWindowCount 7 * sourceAction 7
      = 197 := by decide

/-- Numerators of the reference law over the per-path denominators. -/
def sourceTauChainNum : Fin 8 → ℕ :=
  ![8764880, 708340, 2675, 269105, 708340, 57245, 269105, 27071963]

/-- Per-path denominators of the reference law, in reduced form. -/
def sourceTauChainDen : Fin 8 → ℕ :=
  ![88022241, 88022241, 31247588, 31247588,
    88022241, 88022241, 31247588, 31247588]

/-- **Cleared-denominator Markov product, kernel-decided.**  Over ℕ,
every path satisfies
`num_g * (61511 * rowSum(s0) * rowSum(s1))
  = piNum(s0) * count(s0,s1) * count(s1,s2) * den_g`:
the identity `tau_chain(g) = pi(s0) P(s0,s1) P(s1,s2)` with every
denominator cleared against the committed integer literals of
`MixingChainRealization`. -/
theorem sourceTauChainNum_markov_integer_identity :
    ∀ g : Fin 8, sourceTauChainNum g
        * (61511 * mixingChainRowSums (sourceState0 g)
          * mixingChainRowSums (sourceState1 g))
      = mixingChainStationaryNum (sourceState0 g)
        * mixingChainCounts (sourceState0 g) (sourceState1 g)
        * mixingChainCounts (sourceState1 g) (sourceState2 g)
        * sourceTauChainDen g := by decide

/-- The action is nonconstant across paths, kernel-decided. -/
theorem sourceAction_nonconstant :
    ¬ ∀ g g' : Fin 8, sourceAction g = sourceAction g' := by decide

theorem sourceTauEmpQ_apply_0 : sourceTauEmpQ 0 = 1149/1754 := rfl
theorem sourceTauEmpQ_apply_1 : sourceTauEmpQ 1 = 48/877 := rfl
theorem sourceTauEmpQ_apply_2 : sourceTauEmpQ 2 = 1/877 := rfl
theorem sourceTauEmpQ_apply_3 : sourceTauEmpQ 3 = 87/1754 := rfl
theorem sourceTauEmpQ_apply_4 : sourceTauEmpQ 4 = 3/1754 := rfl
theorem sourceTauEmpQ_apply_5 : sourceTauEmpQ 5 = 1/877 := rfl
theorem sourceTauEmpQ_apply_6 : sourceTauEmpQ 6 = 3/1754 := rfl
theorem sourceTauEmpQ_apply_7 : sourceTauEmpQ 7 = 206/877 := rfl

theorem sourceTauChainQ_apply_0 : sourceTauChainQ 0 = 8764880/88022241 := rfl
theorem sourceTauChainQ_apply_1 : sourceTauChainQ 1 = 708340/88022241 := rfl
theorem sourceTauChainQ_apply_2 : sourceTauChainQ 2 = 2675/31247588 := rfl
theorem sourceTauChainQ_apply_3 : sourceTauChainQ 3 = 269105/31247588 := rfl
theorem sourceTauChainQ_apply_4 : sourceTauChainQ 4 = 708340/88022241 := rfl
theorem sourceTauChainQ_apply_5 : sourceTauChainQ 5 = 57245/88022241 := rfl
theorem sourceTauChainQ_apply_6 : sourceTauChainQ 6 = 269105/31247588 := rfl
theorem sourceTauChainQ_apply_7 : sourceTauChainQ 7 = 27071963/31247588 := rfl

theorem sourceAction_apply_0 : sourceAction 0 = 0 := rfl
theorem sourceAction_apply_1 : sourceAction 1 = 1 := rfl
theorem sourceAction_apply_2 : sourceAction 2 = 2 := rfl
theorem sourceAction_apply_3 : sourceAction 3 = 1 := rfl
theorem sourceAction_apply_4 : sourceAction 4 = 1 := rfl
theorem sourceAction_apply_5 : sourceAction 5 = 2 := rfl
theorem sourceAction_apply_6 : sourceAction 6 = 1 := rfl
theorem sourceAction_apply_7 : sourceAction 7 = 0 := rfl

theorem sourceWindowCount_apply_0 : sourceWindowCount 0 = 1149 := rfl
theorem sourceWindowCount_apply_1 : sourceWindowCount 1 = 96 := rfl
theorem sourceWindowCount_apply_2 : sourceWindowCount 2 = 2 := rfl
theorem sourceWindowCount_apply_3 : sourceWindowCount 3 = 87 := rfl
theorem sourceWindowCount_apply_4 : sourceWindowCount 4 = 3 := rfl
theorem sourceWindowCount_apply_5 : sourceWindowCount 5 = 2 := rfl
theorem sourceWindowCount_apply_6 : sourceWindowCount 6 = 3 := rfl
theorem sourceWindowCount_apply_7 : sourceWindowCount 7 = 412 := rfl

/-- The empirical law is the exact normalization of the window counts,
entry by entry. -/
theorem sourceTauEmpQ_eq_counts :
    ∀ g, sourceTauEmpQ g = (sourceWindowCount g : ℚ) / 1754 := by
  have h0 : sourceTauEmpQ 0 = (sourceWindowCount 0 : ℚ) / 1754 := by
    rw [sourceTauEmpQ_apply_0, sourceWindowCount_apply_0]; norm_num
  have h1 : sourceTauEmpQ 1 = (sourceWindowCount 1 : ℚ) / 1754 := by
    rw [sourceTauEmpQ_apply_1, sourceWindowCount_apply_1]; norm_num
  have h2 : sourceTauEmpQ 2 = (sourceWindowCount 2 : ℚ) / 1754 := by
    rw [sourceTauEmpQ_apply_2, sourceWindowCount_apply_2]; norm_num
  have h3 : sourceTauEmpQ 3 = (sourceWindowCount 3 : ℚ) / 1754 := by
    rw [sourceTauEmpQ_apply_3, sourceWindowCount_apply_3]; norm_num
  have h4 : sourceTauEmpQ 4 = (sourceWindowCount 4 : ℚ) / 1754 := by
    rw [sourceTauEmpQ_apply_4, sourceWindowCount_apply_4]; norm_num
  have h5 : sourceTauEmpQ 5 = (sourceWindowCount 5 : ℚ) / 1754 := by
    rw [sourceTauEmpQ_apply_5, sourceWindowCount_apply_5]; norm_num
  have h6 : sourceTauEmpQ 6 = (sourceWindowCount 6 : ℚ) / 1754 := by
    rw [sourceTauEmpQ_apply_6, sourceWindowCount_apply_6]; norm_num
  have h7 : sourceTauEmpQ 7 = (sourceWindowCount 7 : ℚ) / 1754 := by
    rw [sourceTauEmpQ_apply_7, sourceWindowCount_apply_7]; norm_num
  intro g
  fin_cases g
  · exact h0
  · exact h1
  · exact h2
  · exact h3
  · exact h4
  · exact h5
  · exact h6
  · exact h7

/-- **The chain path law is the Markov product.**  Every entry of the
reference law equals `pi(s0) * P(s0,s1) * P(s1,s2)` with the committed
exact stationary law and transition matrix of `MixingChainRealization`;
the cleared-denominator form of the same identity is kernel-decided at
`sourceTauChainNum_markov_integer_identity`. -/
theorem sourceTauChainQ_eq_markov_product :
    ∀ g, sourceTauChainQ g
      = mixingChainStationary (sourceState0 g)
        * mixingChain (sourceState0 g) (sourceState1 g)
        * mixingChain (sourceState1 g) (sourceState2 g) := by
  have h0 : sourceTauChainQ 0
      = mixingChainStationary 0 * mixingChain 0 0 * mixingChain 0 0 := by
    rw [sourceTauChainQ_apply_0]
    norm_num [mixingChainStationary, mixingChain]
  have h1 : sourceTauChainQ 1
      = mixingChainStationary 0 * mixingChain 0 0 * mixingChain 0 1 := by
    rw [sourceTauChainQ_apply_1]
    norm_num [mixingChainStationary, mixingChain]
  have h2 : sourceTauChainQ 2
      = mixingChainStationary 0 * mixingChain 0 1 * mixingChain 1 0 := by
    rw [sourceTauChainQ_apply_2]
    norm_num [mixingChainStationary, mixingChain]
  have h3 : sourceTauChainQ 3
      = mixingChainStationary 0 * mixingChain 0 1 * mixingChain 1 1 := by
    rw [sourceTauChainQ_apply_3]
    norm_num [mixingChainStationary, mixingChain]
  have h4 : sourceTauChainQ 4
      = mixingChainStationary 1 * mixingChain 1 0 * mixingChain 0 0 := by
    rw [sourceTauChainQ_apply_4]
    norm_num [mixingChainStationary, mixingChain]
  have h5 : sourceTauChainQ 5
      = mixingChainStationary 1 * mixingChain 1 0 * mixingChain 0 1 := by
    rw [sourceTauChainQ_apply_5]
    norm_num [mixingChainStationary, mixingChain]
  have h6 : sourceTauChainQ 6
      = mixingChainStationary 1 * mixingChain 1 1 * mixingChain 1 0 := by
    rw [sourceTauChainQ_apply_6]
    norm_num [mixingChainStationary, mixingChain]
  have h7 : sourceTauChainQ 7
      = mixingChainStationary 1 * mixingChain 1 1 * mixingChain 1 1 := by
    rw [sourceTauChainQ_apply_7]
    norm_num [mixingChainStationary, mixingChain]
  intro g
  fin_cases g
  · exact h0
  · exact h1
  · exact h2
  · exact h3
  · exact h4
  · exact h5
  · exact h6
  · exact h7

/-- The empirical path law is strictly positive. -/
theorem sourceTauEmpQ_pos : ∀ g, 0 < sourceTauEmpQ g := by
  intro g
  rw [sourceTauEmpQ_eq_counts g]
  have hc : (0 : ℚ) < (sourceWindowCount g : ℚ) := by
    exact_mod_cast sourceWindowCount_pos g
  positivity

/-- The reference path law is strictly positive: the Markov product of
the strictly positive committed literals. -/
theorem sourceTauChainQ_pos : ∀ g, 0 < sourceTauChainQ g := by
  intro g
  rw [sourceTauChainQ_eq_markov_product g]
  exact mul_pos
    (mul_pos (mixingChainStationary_pos _) (mixingChain_entries_pos _ _))
    (mixingChain_entries_pos _ _)

/-- The two laws disagree (at the all-`18` path). -/
theorem sourceTauEmpQ_ne_chain : sourceTauEmpQ 0 ≠ sourceTauChainQ 0 := by
  rw [sourceTauEmpQ_apply_0, sourceTauChainQ_apply_0]
  norm_num

/-- The empirical law is normalized over the eight paths. -/
theorem sourceTauEmpQ_sum : ∑ g, sourceTauEmpQ g = 1 := by
  rw [Fin.sum_univ_eight, sourceTauEmpQ_apply_0, sourceTauEmpQ_apply_1,
    sourceTauEmpQ_apply_2, sourceTauEmpQ_apply_3, sourceTauEmpQ_apply_4,
    sourceTauEmpQ_apply_5, sourceTauEmpQ_apply_6, sourceTauEmpQ_apply_7]
  norm_num

/-- The reference law is normalized over the eight paths. -/
theorem sourceTauChainQ_sum : ∑ g, sourceTauChainQ g = 1 := by
  rw [Fin.sum_univ_eight, sourceTauChainQ_apply_0, sourceTauChainQ_apply_1,
    sourceTauChainQ_apply_2, sourceTauChainQ_apply_3,
    sourceTauChainQ_apply_4, sourceTauChainQ_apply_5,
    sourceTauChainQ_apply_6, sourceTauChainQ_apply_7]
  norm_num

/-- Exact empirical mean action `197/1754`. -/
theorem sourceTauEmpQ_meanAction :
    ∑ g, sourceTauEmpQ g * (sourceAction g : ℚ) = 197 / 1754 := by
  rw [Fin.sum_univ_eight, sourceTauEmpQ_apply_0, sourceTauEmpQ_apply_1,
    sourceTauEmpQ_apply_2, sourceTauEmpQ_apply_3, sourceTauEmpQ_apply_4,
    sourceTauEmpQ_apply_5, sourceTauEmpQ_apply_6, sourceTauEmpQ_apply_7,
    sourceAction_apply_0, sourceAction_apply_1, sourceAction_apply_2,
    sourceAction_apply_3, sourceAction_apply_4, sourceAction_apply_5,
    sourceAction_apply_6, sourceAction_apply_7]
  norm_num

/-- Exact reference mean action `2140/61511`. -/
theorem sourceTauChainQ_meanAction :
    ∑ g, sourceTauChainQ g * (sourceAction g : ℚ) = 2140 / 61511 := by
  rw [Fin.sum_univ_eight, sourceTauChainQ_apply_0, sourceTauChainQ_apply_1,
    sourceTauChainQ_apply_2, sourceTauChainQ_apply_3,
    sourceTauChainQ_apply_4, sourceTauChainQ_apply_5,
    sourceTauChainQ_apply_6, sourceTauChainQ_apply_7,
    sourceAction_apply_0, sourceAction_apply_1, sourceAction_apply_2,
    sourceAction_apply_3, sourceAction_apply_4, sourceAction_apply_5,
    sourceAction_apply_6, sourceAction_apply_7]
  norm_num

/-! ## The real mirror -/

/-- The empirical path law over ℝ. -/
noncomputable def sourceTauEmpR : Fin 8 → ℝ := fun g => (sourceTauEmpQ g : ℝ)

/-- The reference path law over ℝ. -/
noncomputable def sourceTauChainR : Fin 8 → ℝ :=
  fun g => (sourceTauChainQ g : ℝ)

/-- The action over ℝ. -/
noncomputable def sourceActionR : Fin 8 → ℝ := fun g => (sourceAction g : ℝ)

theorem sourceTauEmpR_pos : ∀ g, 0 < sourceTauEmpR g := fun g => by
  simp only [sourceTauEmpR]
  exact_mod_cast sourceTauEmpQ_pos g

theorem sourceTauChainR_pos : ∀ g, 0 < sourceTauChainR g := fun g => by
  simp only [sourceTauChainR]
  exact_mod_cast sourceTauChainQ_pos g

theorem sourceTauEmpR_sum : ∑ g, sourceTauEmpR g = 1 := by
  simp only [sourceTauEmpR]
  rw [← Rat.cast_sum, sourceTauEmpQ_sum, Rat.cast_one]

theorem sourceTauChainR_sum : ∑ g, sourceTauChainR g = 1 := by
  simp only [sourceTauChainR]
  rw [← Rat.cast_sum, sourceTauChainQ_sum, Rat.cast_one]

/-- Exact empirical mean action over ℝ. -/
theorem sourceTauEmpR_meanAction :
    meanAction sourceTauEmpR sourceActionR = 197 / 1754 := by
  unfold meanAction
  rw [Finset.sum_congr rfl fun g _ =>
    show sourceTauEmpR g * sourceActionR g
        = ((sourceTauEmpQ g * (sourceAction g : ℚ) : ℚ) : ℝ) by
      simp only [sourceTauEmpR, sourceActionR]
      push_cast
      ring]
  rw [← Rat.cast_sum, sourceTauEmpQ_meanAction]
  norm_num

/-- Exact reference mean action over ℝ. -/
theorem sourceTauChainR_meanAction :
    meanAction sourceTauChainR sourceActionR = 2140 / 61511 := by
  unfold meanAction
  rw [Finset.sum_congr rfl fun g _ =>
    show sourceTauChainR g * sourceActionR g
        = ((sourceTauChainQ g * (sourceAction g : ℚ) : ℚ) : ℝ) by
      simp only [sourceTauChainR, sourceActionR]
      push_cast
      ring]
  rw [← Rat.cast_sum, sourceTauChainQ_meanAction]
  norm_num

/-- The two laws differ as real families. -/
theorem sourceTauEmpR_ne_chain : sourceTauEmpR ≠ sourceTauChainR := by
  intro h
  have h0 := congrFun h 0
  simp only [sourceTauEmpR, sourceTauChainR] at h0
  exact absurd (by exact_mod_cast h0 :
    sourceTauEmpQ 0 = sourceTauChainQ 0) sourceTauEmpQ_ne_chain

/-- **Strict divergence of the two produced laws.**  The empirical path
law scores strictly positive relative entropy against the reference:
the Lean side of the payload's exact-enclosure receipt
`kl_emp_chain_strictly_positive`. -/
theorem sourceKl_emp_chain_pos : 0 < kl sourceTauEmpR sourceTauChainR :=
  kl_pos_of_ne sourceTauEmpR sourceTauChainR
    (fun g => le_of_lt (sourceTauEmpR_pos g))
    (fun g => le_of_lt (sourceTauChainR_pos g))
    (fun g hg => absurd hg (ne_of_gt (sourceTauChainR_pos g)))
    (by rw [sourceTauEmpR_sum, sourceTauChainR_sum])
    sourceTauEmpR_ne_chain

/-- At multiplier zero the tilt of the reference law by the action is
the reference law exactly. -/
theorem sourceTilt_zero :
    tilt sourceTauChainR sourceActionR 0 = sourceTauChainR :=
  tilt_zero sourceTauChainR sourceActionR sourceTauChainR_sum

/-- **Exact packet inhabitation of the history-law interface.**  The
`HistoryLaw` projection theorems hold at the literals of this packet:
for every multiplier, the tilt of the declared extracted reference law
by the declared repair-count action is a strictly positive normalized path law,
and it is the unique minimizer of relative entropy to the reference
among laws sharing its mean action.  The multiplier is the one slot the
packet does not fill: `lam` is a free input here, the payload matches it
to the declared empirical mean action of the same run, and the
multiplier selection principle from source data is the open core of
issue B7. -/
theorem sourceHistoryLaw_packet_inhabitation (lam : ℝ) :
    (∀ g, 0 < tilt sourceTauChainR sourceActionR lam g)
      ∧ (∑ g, tilt sourceTauChainR sourceActionR lam g = 1)
      ∧ (∀ ρ : Fin 8 → ℝ, (∀ g, 0 ≤ ρ g) → (∑ g, ρ g = 1) →
          meanAction ρ sourceActionR
            = meanAction (tilt sourceTauChainR sourceActionR lam)
                sourceActionR →
          kl (tilt sourceTauChainR sourceActionR lam) sourceTauChainR
            ≤ kl ρ sourceTauChainR)
      ∧ (∀ ρ : Fin 8 → ℝ, (∀ g, 0 ≤ ρ g) → (∑ g, ρ g = 1) →
          meanAction ρ sourceActionR
            = meanAction (tilt sourceTauChainR sourceActionR lam)
                sourceActionR →
          kl ρ sourceTauChainR
            = kl (tilt sourceTauChainR sourceActionR lam) sourceTauChainR →
          ρ = tilt sourceTauChainR sourceActionR lam) :=
  ⟨tilt_pos sourceTauChainR sourceActionR lam sourceTauChainR_pos,
    tilt_total sourceTauChainR sourceActionR lam sourceTauChainR_pos,
    fun ρ h0 h1 hm => tilt_constrained_minimizer sourceTauChainR
      sourceActionR lam ρ h0 h1 sourceTauChainR_pos hm,
    fun ρ h0 h1 hm he => tilt_constrained_unique sourceTauChainR
      sourceActionR lam ρ h0 h1 sourceTauChainR_pos hm he⟩

/-- The empirical law is admissible for the constrained comparison: it
is nonnegative and normalized, so at any multiplier whose tilted mean
action equals the empirical mean action, the tilt scores at most the
empirical law's relative entropy to the reference. -/
theorem sourceTauEmpR_admissible :
    (∀ g, 0 ≤ sourceTauEmpR g) ∧ ∑ g, sourceTauEmpR g = 1 :=
  ⟨fun g => le_of_lt (sourceTauEmpR_pos g), sourceTauEmpR_sum⟩

/-! ## The matching multiplier: intermediate-value receipt -/

/-- The mean-action matching quadratic of the payload, in
`x = exp (-lambda)`: its positive root is the multiplier whose tilted
mean action equals the declared empirical mean `197/1754`. -/
noncomputable def sourceMatchQuad (x : ℝ) : ℝ :=
  1771385 / 1275066792 * x ^ 2 + 128872193675 / 4357257413484 * x
    - 8508930004321 / 78430633442712

theorem sourceTauChainR_apply_0 :
    sourceTauChainR 0 = 8764880 / 88022241 := by
  simp only [sourceTauChainR, sourceTauChainQ_apply_0]
  norm_num

theorem sourceTauChainR_apply_1 :
    sourceTauChainR 1 = 708340 / 88022241 := by
  simp only [sourceTauChainR, sourceTauChainQ_apply_1]
  norm_num

theorem sourceTauChainR_apply_2 :
    sourceTauChainR 2 = 2675 / 31247588 := by
  simp only [sourceTauChainR, sourceTauChainQ_apply_2]
  norm_num

theorem sourceTauChainR_apply_3 :
    sourceTauChainR 3 = 269105 / 31247588 := by
  simp only [sourceTauChainR, sourceTauChainQ_apply_3]
  norm_num

theorem sourceTauChainR_apply_4 :
    sourceTauChainR 4 = 708340 / 88022241 := by
  simp only [sourceTauChainR, sourceTauChainQ_apply_4]
  norm_num

theorem sourceTauChainR_apply_5 :
    sourceTauChainR 5 = 57245 / 88022241 := by
  simp only [sourceTauChainR, sourceTauChainQ_apply_5]
  norm_num

theorem sourceTauChainR_apply_6 :
    sourceTauChainR 6 = 269105 / 31247588 := by
  simp only [sourceTauChainR, sourceTauChainQ_apply_6]
  norm_num

theorem sourceTauChainR_apply_7 :
    sourceTauChainR 7 = 27071963 / 31247588 := by
  simp only [sourceTauChainR, sourceTauChainQ_apply_7]
  norm_num

theorem sourceActionR_apply_0 : sourceActionR 0 = 0 := by
  simp only [sourceActionR, sourceAction_apply_0]
  norm_num

theorem sourceActionR_apply_1 : sourceActionR 1 = 1 := by
  simp only [sourceActionR, sourceAction_apply_1]
  norm_num

theorem sourceActionR_apply_2 : sourceActionR 2 = 2 := by
  simp only [sourceActionR, sourceAction_apply_2]
  norm_num

theorem sourceActionR_apply_3 : sourceActionR 3 = 1 := by
  simp only [sourceActionR, sourceAction_apply_3]
  norm_num

theorem sourceActionR_apply_4 : sourceActionR 4 = 1 := by
  simp only [sourceActionR, sourceAction_apply_4]
  norm_num

theorem sourceActionR_apply_5 : sourceActionR 5 = 2 := by
  simp only [sourceActionR, sourceAction_apply_5]
  norm_num

theorem sourceActionR_apply_6 : sourceActionR 6 = 1 := by
  simp only [sourceActionR, sourceAction_apply_6]
  norm_num

theorem sourceActionR_apply_7 : sourceActionR 7 = 0 := by
  simp only [sourceActionR, sourceAction_apply_7]
  norm_num

/-- Value of the exponential tilt at multiplier `-log x` for `x > 0`:
the reference times `x` to the action, normalized. -/
theorem sourceTilt_neg_log_apply {x : ℝ} (hx : 0 < x) (g : Fin 8) :
    tilt sourceTauChainR sourceActionR (-Real.log x) g
      = sourceTauChainR g * x ^ sourceAction g
        / ∑ g', sourceTauChainR g' * x ^ sourceAction g' := by
  have hexp : ∀ g' : Fin 8,
      Real.exp (-(-Real.log x * sourceActionR g'))
        = x ^ sourceAction g' := by
    intro g'
    have harg : -(-Real.log x * sourceActionR g')
        = (sourceAction g' : ℝ) * Real.log x := by
      simp only [sourceActionR]
      ring
    rw [harg, Real.exp_nat_mul, Real.exp_log hx]
  unfold tilt tiltZ
  rw [hexp g, Finset.sum_congr rfl fun g' _ => by rw [hexp g']]

/-- The tilt normalizer at `-log x` is strictly positive for `x > 0`. -/
theorem sourceTiltDenom_pos {x : ℝ} (hx : 0 < x) :
    0 < ∑ g, sourceTauChainR g * x ^ sourceAction g :=
  Finset.sum_pos
    (fun g _ => mul_pos (sourceTauChainR_pos g) (pow_pos hx _))
    Finset.univ_nonempty

/-- Cleared-denominator mean action of the tilt at `-log x`. -/
theorem sourceTilt_neg_log_meanAction_mul {x : ℝ} (hx : 0 < x) :
    meanAction (tilt sourceTauChainR sourceActionR (-Real.log x))
        sourceActionR
        * ∑ g, sourceTauChainR g * x ^ sourceAction g
      = ∑ g, sourceTauChainR g * sourceActionR g * x ^ sourceAction g := by
  unfold meanAction
  rw [Finset.sum_congr rfl fun g _ =>
    show tilt sourceTauChainR sourceActionR (-Real.log x) g
        * sourceActionR g
      = sourceTauChainR g * sourceActionR g * x ^ sourceAction g
        / ∑ g', sourceTauChainR g' * x ^ sourceAction g' by
      rw [sourceTilt_neg_log_apply hx g]
      ring]
  rw [← Finset.sum_div, div_mul_cancel₀ _ (ne_of_gt (sourceTiltDenom_pos hx))]

/-- The matching quadratic measures the gap between the tilted mean
action and the declared empirical mean, with the normalizer cleared. -/
theorem sourceMatchQuad_key (x : ℝ) :
    sourceMatchQuad x
      = (∑ g, sourceTauChainR g * sourceActionR g * x ^ sourceAction g)
        - 197 / 1754 * ∑ g, sourceTauChainR g * x ^ sourceAction g := by
  rw [Fin.sum_univ_eight, Fin.sum_univ_eight, sourceTauChainR_apply_0,
    sourceTauChainR_apply_1, sourceTauChainR_apply_2,
    sourceTauChainR_apply_3, sourceTauChainR_apply_4,
    sourceTauChainR_apply_5, sourceTauChainR_apply_6,
    sourceTauChainR_apply_7, sourceActionR_apply_0, sourceActionR_apply_1,
    sourceActionR_apply_2, sourceActionR_apply_3, sourceActionR_apply_4,
    sourceActionR_apply_5, sourceActionR_apply_6, sourceActionR_apply_7,
    sourceAction_apply_0, sourceAction_apply_1, sourceAction_apply_2,
    sourceAction_apply_3, sourceAction_apply_4, sourceAction_apply_5,
    sourceAction_apply_6, sourceAction_apply_7]
  unfold sourceMatchQuad
  ring

/-- Sign of the matching quadratic at the lower bracket endpoint:
`F(31901/10000) = -18214945532951/522870889618080000000 < 0`. -/
theorem sourceMatchQuad_lo_neg : sourceMatchQuad (31901 / 10000) < 0 := by
  unfold sourceMatchQuad
  norm_num

/-- Sign of the matching quadratic at the upper bracket endpoint:
`F(3190101/1000000) = 56528441069771947/15686126688542400000000000 > 0`. -/
theorem sourceMatchQuad_hi_pos : 0 < sourceMatchQuad (3190101 / 1000000) := by
  unfold sourceMatchQuad
  norm_num

/-- **Intermediate-value receipt.**  A multiplier matching the declared
empirical mean action exists inside the payload's rational bracket:
some `x` with `31901/10000 < x < 3190101/1000000` makes the tilted mean
action at multiplier `-log x` equal `197/1754` exactly.  The bracket
endpoints and their quadratic signs are the committed payload literals;
existence is the one nonconstructive step, and the payload carries an
exact-rational bisection enclosure of width at most `10^-40` for the
same root. -/
theorem sourceMatchingMultiplier_exists :
    ∃ x : ℝ, 31901 / 10000 < x ∧ x < 3190101 / 1000000
      ∧ meanAction (tilt sourceTauChainR sourceActionR (-Real.log x))
          sourceActionR = 197 / 1754 := by
  have hle : (31901 : ℝ) / 10000 ≤ 3190101 / 1000000 := by norm_num
  have hcont : ContinuousOn sourceMatchQuad
      (Set.Icc (31901 / 10000) (3190101 / 1000000)) := by
    apply Continuous.continuousOn
    unfold sourceMatchQuad
    fun_prop
  have hmem : (0 : ℝ) ∈ Set.Ioo (sourceMatchQuad (31901 / 10000))
      (sourceMatchQuad (3190101 / 1000000)) :=
    ⟨sourceMatchQuad_lo_neg, sourceMatchQuad_hi_pos⟩
  obtain ⟨x, hx, hfx⟩ := intermediate_value_Ioo hle hcont hmem
  have hxpos : (0 : ℝ) < x := lt_trans (by norm_num) hx.1
  refine ⟨x, hx.1, hx.2, ?_⟩
  have hkey := sourceMatchQuad_key x
  rw [hfx] at hkey
  have hmean := sourceTilt_neg_log_meanAction_mul hxpos
  have hZ := sourceTiltDenom_pos hxpos
  have hprod : meanAction (tilt sourceTauChainR sourceActionR
        (-Real.log x)) sourceActionR
        * ∑ g, sourceTauChainR g * x ^ sourceAction g
      = 197 / 1754 * ∑ g, sourceTauChainR g * x ^ sourceAction g := by
    rw [hmean]
    linarith
  exact mul_right_cancel₀ (ne_of_gt hZ) hprod

/-! ## The declared multiplier point -/

/-- The declared multiplier of the payload:
`lambda0 = -log x0` with `x0 = 3190100906251/10^12`, the payload's
truncation of the matching root. -/
noncomputable def sourceLambdaZero : ℝ :=
  -Real.log (3190100906251 / 1000000000000)

/-- The exactly rational tilted law at the declared multiplier point. -/
def sourceTiltX0Q : Fin 8 → ℚ :=
  ![890511808000000000000000000000000/9656043270326080709108107703956677,
    229583089314877467344000000000000/9656043270326080709108107703956677,
    288561570223959167023760653355/357631232234299285522522507553951,
    9099804306392758763000000000000/357631232234299285522522507553951,
    229583089314877467344000000000000/9656043270326080709108107703956677,
    59188878155069903708466166316092/9656043270326080709108107703956677,
    9099804306392758763000000000000/357631232234299285522522507553951,
    286962807800000000000000000000000/357631232234299285522522507553951]

theorem sourceTiltX0Q_apply_0 : sourceTiltX0Q 0
    = 890511808000000000000000000000000/9656043270326080709108107703956677 := rfl
theorem sourceTiltX0Q_apply_1 : sourceTiltX0Q 1
    = 229583089314877467344000000000000/9656043270326080709108107703956677 := rfl
theorem sourceTiltX0Q_apply_2 : sourceTiltX0Q 2
    = 288561570223959167023760653355/357631232234299285522522507553951 := rfl
theorem sourceTiltX0Q_apply_3 : sourceTiltX0Q 3
    = 9099804306392758763000000000000/357631232234299285522522507553951 := rfl
theorem sourceTiltX0Q_apply_4 : sourceTiltX0Q 4
    = 229583089314877467344000000000000/9656043270326080709108107703956677 := rfl
theorem sourceTiltX0Q_apply_5 : sourceTiltX0Q 5
    = 59188878155069903708466166316092/9656043270326080709108107703956677 := rfl
theorem sourceTiltX0Q_apply_6 : sourceTiltX0Q 6
    = 9099804306392758763000000000000/357631232234299285522522507553951 := rfl
theorem sourceTiltX0Q_apply_7 : sourceTiltX0Q 7
    = 286962807800000000000000000000000/357631232234299285522522507553951 := rfl

/-- The declared tilted law is normalized. -/
theorem sourceTiltX0Q_sum : ∑ g, sourceTiltX0Q g = 1 := by
  rw [Fin.sum_univ_eight, sourceTiltX0Q_apply_0, sourceTiltX0Q_apply_1,
    sourceTiltX0Q_apply_2, sourceTiltX0Q_apply_3, sourceTiltX0Q_apply_4,
    sourceTiltX0Q_apply_5, sourceTiltX0Q_apply_6, sourceTiltX0Q_apply_7]
  norm_num

/-- The tilt normalizer at the declared point evaluates to the exact
rational `z0` of the payload. -/
theorem sourceTiltDenom_x0_eval :
    ∑ g, sourceTauChainR g
        * ((3190100906251 : ℝ) / 1000000000000) ^ sourceAction g
      = 357631232234299285522522507553951
        / 331224432800000000000000000000000 := by
  rw [Fin.sum_univ_eight, sourceTauChainR_apply_0, sourceTauChainR_apply_1,
    sourceTauChainR_apply_2, sourceTauChainR_apply_3,
    sourceTauChainR_apply_4, sourceTauChainR_apply_5,
    sourceTauChainR_apply_6, sourceTauChainR_apply_7,
    sourceAction_apply_0, sourceAction_apply_1, sourceAction_apply_2,
    sourceAction_apply_3, sourceAction_apply_4, sourceAction_apply_5,
    sourceAction_apply_6, sourceAction_apply_7]
  norm_num

/-- **The tilt at the declared multiplier is exactly rational.**  At
`lambda0 = -log (3190100906251/10^12)` the tilt of the reference law by
the action equals the declared rational law of the payload, entry by
entry. -/
theorem sourceTilt_lambdaZero_apply :
    ∀ g, tilt sourceTauChainR sourceActionR sourceLambdaZero g
      = ((sourceTiltX0Q g : ℚ) : ℝ) := by
  have hx : (0 : ℝ) < 3190100906251 / 1000000000000 := by norm_num
  have key : ∀ g : Fin 8,
      tilt sourceTauChainR sourceActionR sourceLambdaZero g
        = sourceTauChainR g
          * ((3190100906251 : ℝ) / 1000000000000) ^ sourceAction g
          / (357631232234299285522522507553951
            / 331224432800000000000000000000000) := by
    intro g
    rw [sourceLambdaZero, sourceTilt_neg_log_apply hx g,
      sourceTiltDenom_x0_eval]
  have h0 := key 0
  have h1 := key 1
  have h2 := key 2
  have h3 := key 3
  have h4 := key 4
  have h5 := key 5
  have h6 := key 6
  have h7 := key 7
  rw [sourceTauChainR_apply_0, sourceAction_apply_0] at h0
  rw [sourceTauChainR_apply_1, sourceAction_apply_1] at h1
  rw [sourceTauChainR_apply_2, sourceAction_apply_2] at h2
  rw [sourceTauChainR_apply_3, sourceAction_apply_3] at h3
  rw [sourceTauChainR_apply_4, sourceAction_apply_4] at h4
  rw [sourceTauChainR_apply_5, sourceAction_apply_5] at h5
  rw [sourceTauChainR_apply_6, sourceAction_apply_6] at h6
  rw [sourceTauChainR_apply_7, sourceAction_apply_7] at h7
  have e0 : tilt sourceTauChainR sourceActionR sourceLambdaZero 0
      = ((sourceTiltX0Q 0 : ℚ) : ℝ) := by
    rw [h0, sourceTiltX0Q_apply_0]
    push_cast
    norm_num
  have e1 : tilt sourceTauChainR sourceActionR sourceLambdaZero 1
      = ((sourceTiltX0Q 1 : ℚ) : ℝ) := by
    rw [h1, sourceTiltX0Q_apply_1]
    push_cast
    norm_num
  have e2 : tilt sourceTauChainR sourceActionR sourceLambdaZero 2
      = ((sourceTiltX0Q 2 : ℚ) : ℝ) := by
    rw [h2, sourceTiltX0Q_apply_2]
    push_cast
    norm_num
  have e3 : tilt sourceTauChainR sourceActionR sourceLambdaZero 3
      = ((sourceTiltX0Q 3 : ℚ) : ℝ) := by
    rw [h3, sourceTiltX0Q_apply_3]
    push_cast
    norm_num
  have e4 : tilt sourceTauChainR sourceActionR sourceLambdaZero 4
      = ((sourceTiltX0Q 4 : ℚ) : ℝ) := by
    rw [h4, sourceTiltX0Q_apply_4]
    push_cast
    norm_num
  have e5 : tilt sourceTauChainR sourceActionR sourceLambdaZero 5
      = ((sourceTiltX0Q 5 : ℚ) : ℝ) := by
    rw [h5, sourceTiltX0Q_apply_5]
    push_cast
    norm_num
  have e6 : tilt sourceTauChainR sourceActionR sourceLambdaZero 6
      = ((sourceTiltX0Q 6 : ℚ) : ℝ) := by
    rw [h6, sourceTiltX0Q_apply_6]
    push_cast
    norm_num
  have e7 : tilt sourceTauChainR sourceActionR sourceLambdaZero 7
      = ((sourceTiltX0Q 7 : ℚ) : ℝ) := by
    rw [h7, sourceTiltX0Q_apply_7]
    push_cast
    norm_num
  intro g
  fin_cases g
  · exact e0
  · exact e1
  · exact e2
  · exact e3
  · exact e4
  · exact e5
  · exact e6
  · exact e7

/-- Exact mean action at the declared multiplier point. -/
theorem sourceTilt_lambdaZero_meanAction :
    meanAction (tilt sourceTauChainR sourceActionR sourceLambdaZero)
        sourceActionR
      = 120501743586355278925135045323706
        / 1072893696702897856567567522661853 := by
  unfold meanAction
  rw [Finset.sum_congr rfl fun g _ => by
    rw [sourceTilt_lambdaZero_apply g]]
  rw [Fin.sum_univ_eight, sourceTiltX0Q_apply_0, sourceTiltX0Q_apply_1,
    sourceTiltX0Q_apply_2, sourceTiltX0Q_apply_3, sourceTiltX0Q_apply_4,
    sourceTiltX0Q_apply_5, sourceTiltX0Q_apply_6, sourceTiltX0Q_apply_7,
    sourceActionR_apply_0, sourceActionR_apply_1, sourceActionR_apply_2,
    sourceActionR_apply_3, sourceActionR_apply_4, sourceActionR_apply_5,
    sourceActionR_apply_6, sourceActionR_apply_7]
  push_cast
  norm_num

/-- **Declared-point mismatch bound.**  The mean action of the tilt at
the declared multiplier lies within `10^-12` of the declared empirical
mean `197/1754`: the Lean side of the payload's
`mean_action_mismatch_bound` receipt. -/
theorem sourceTilt_lambdaZero_mean_close :
    |meanAction (tilt sourceTauChainR sourceActionR sourceLambdaZero)
        sourceActionR - 197 / 1754| ≤ 1 / 1000000000000 := by
  rw [sourceTilt_lambdaZero_meanAction, abs_le]
  constructor <;> norm_num

/-! ## Packet receipt -/

/-- **Source-attached history-law packet receipt (B7, issue #683).**
Both declared extracted path laws are strictly positive and normalized, the
action is nonconstant, the reference law is the Markov product of the
committed mixing chain (`sourceTauChainQ_eq_markov_product`, with the
cleared-denominator identity kernel-decided at
`sourceTauChainNum_markov_integer_identity`), the tilt at multiplier
zero is the reference law exactly, the empirical law diverges strictly
from the reference, and a multiplier matching the declared empirical
mean action exists inside the payload's rational bracket.  The multiplier slot
itself remains a declared input: producing it from source data is the
open core of B7. -/
theorem sourceHistoryPacket_receipt :
    (∀ g, 0 < sourceTauEmpR g) ∧ (∀ g, 0 < sourceTauChainR g)
      ∧ (∑ g, sourceTauEmpR g = 1) ∧ (∑ g, sourceTauChainR g = 1)
      ∧ (¬ ∀ g g' : Fin 8, sourceAction g = sourceAction g')
      ∧ tilt sourceTauChainR sourceActionR 0 = sourceTauChainR
      ∧ 0 < kl sourceTauEmpR sourceTauChainR
      ∧ (∃ x : ℝ, 31901 / 10000 < x ∧ x < 3190101 / 1000000
          ∧ meanAction (tilt sourceTauChainR sourceActionR (-Real.log x))
              sourceActionR = 197 / 1754) :=
  ⟨sourceTauEmpR_pos, sourceTauChainR_pos, sourceTauEmpR_sum,
    sourceTauChainR_sum, sourceAction_nonconstant, sourceTilt_zero,
    sourceKl_emp_chain_pos, sourceMatchingMultiplier_exists⟩

end OPH.InformationProjection

#print axioms OPH.InformationProjection.sourceState_encoding
#print axioms OPH.InformationProjection.sourceAction_eq_changes
#print axioms OPH.InformationProjection.sourceTauChainNum_markov_integer_identity
#print axioms OPH.InformationProjection.sourceTauEmpQ_eq_counts
#print axioms OPH.InformationProjection.sourceTauChainQ_eq_markov_product
#print axioms OPH.InformationProjection.sourceTauEmpQ_pos
#print axioms OPH.InformationProjection.sourceTauChainQ_pos
#print axioms OPH.InformationProjection.sourceAction_nonconstant
#print axioms OPH.InformationProjection.sourceTauEmpR_meanAction
#print axioms OPH.InformationProjection.sourceTauChainR_meanAction
#print axioms OPH.InformationProjection.sourceKl_emp_chain_pos
#print axioms OPH.InformationProjection.sourceTilt_zero
#print axioms OPH.InformationProjection.sourceHistoryLaw_packet_inhabitation
#print axioms OPH.InformationProjection.sourceMatchingMultiplier_exists
#print axioms OPH.InformationProjection.sourceTilt_lambdaZero_apply
#print axioms OPH.InformationProjection.sourceTilt_lambdaZero_mean_close
#print axioms OPH.InformationProjection.sourceHistoryPacket_receipt
