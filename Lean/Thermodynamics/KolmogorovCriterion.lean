import Mathlib

/-!
# Finite rate reconstruction on the reversible layer (Kolmogorov cycle criterion)

This file formalizes the theorem "Finite rate reconstruction on the reversible layer"
of `paper/screen_microphysics_and_observer_synchronization.tex` (the finite-graph
potential-theory core of the tracked claim `OPH-SCREEN-MODULAR-GEARING`).

Paper statement.  Let `G = (X,E)` be a connected finite state graph with positive rates
`q_{u→v}` in both directions of every edge.  Define the log-rate edge field
`a_{u→v} := log (q_{v→u} / q_{u→v})`.  The following are equivalent:
  (i)   every oriented cycle has zero affinity;
  (ii)  there is a potential `K : X → ℝ`, unique up to an additive constant, with
        `a_{u→v} = K_v − K_u`;
  (iii) `π_x ∝ e^{−K_x}` satisfies detailed balance `π_u q_{u→v} = π_v q_{v→u}`.

Formalization notes.
* Graph objects are Mathlib's `SimpleGraph`, `SimpleGraph.Walk`, `SimpleGraph.Connected`.
  This is the faithful rendering of "connected finite state graph with both directions of
  every edge": `SimpleGraph.Adj` is symmetric and irreflexive.
* The coboundary is written in the same `head − tail` convention as
  `Lean/Thermodynamics/GraphDiffusion.lean` (`edgeDifference G u e = u (head e) − u (tail e)`):
  here a dart `d` plays the role of an oriented edge and the potential condition reads
  `affinity q d.fst d.snd = K d.snd − K d.fst`.  Nothing from `GraphDiffusion` is imported;
  only its coboundary convention is matched, so the statement could later slot in.
* Condition (i) is rendered as: every closed walk has zero affinity.  This is exactly the
  paper's proof object ("zero cycle affinity gives path independence") and the cochain
  reading "the log-rate 1-cochain has zero integral on every cycle" (a cycle = a closed
  1-chain).  A simple cycle is a special closed walk, so this (i) implies the narrow
  simple-cycle form; the converse (simple cycles ⇒ all closed walks) needs a closed-walk
  decomposition not used here.
* `[Fintype V]` (paper: "finite state graph") is carried on the main theorem as a hypothesis;
  the equivalence itself does not consume it, so the result is in fact slightly stronger than
  the paper's (it holds for any connected graph).  A stronger theorem is not a weaker sibling.
-/

namespace OPH.Thermodynamics.KolmogorovCriterion

open SimpleGraph

variable {V : Type*}

/-- The log-rate edge field (affinity 1-cochain): `a_{u→v} = log (q_{v→u} / q_{u→v})`. -/
noncomputable def affinity (q : V → V → ℝ) (u v : V) : ℝ :=
  Real.log (q v u / q u v)

/-- The affinity 1-cochain is antisymmetric, with no positivity hypothesis, because
`log (x⁻¹) = − log x` and `(a/b)⁻¹ = b/a` hold unconditionally in `ℝ`. -/
theorem affinity_antisymm (q : V → V → ℝ) (u v : V) :
    affinity q u v = - affinity q v u := by
  unfold affinity
  rw [← Real.log_inv, inv_div]

/-- The affinity accumulated along a walk: the sum of `affinity` over the walk's darts.
This is the discrete line integral of the 1-cochain along the walk. -/
noncomputable def walkAffinity {G : SimpleGraph V} (q : V → V → ℝ) {u v : V}
    (p : G.Walk u v) : ℝ :=
  (p.darts.map (fun d => affinity q d.fst d.snd)).sum

@[simp] theorem walkAffinity_nil {G : SimpleGraph V} (q : V → V → ℝ) (u : V) :
    walkAffinity q (Walk.nil : G.Walk u u) = 0 := by
  simp [walkAffinity]

theorem walkAffinity_cons {G : SimpleGraph V} (q : V → V → ℝ) {u v w : V}
    (h : G.Adj u v) (p : G.Walk v w) :
    walkAffinity q (Walk.cons h p) = affinity q u v + walkAffinity q p := by
  simp [walkAffinity, Walk.darts_cons]

theorem walkAffinity_append {G : SimpleGraph V} (q : V → V → ℝ) {u v w : V}
    (p : G.Walk u v) (p' : G.Walk v w) :
    walkAffinity q (p.append p') = walkAffinity q p + walkAffinity q p' := by
  simp [walkAffinity, Walk.darts_append]

/-- Reversing a walk negates its affinity (uses antisymmetry of the cochain). -/
theorem walkAffinity_reverse {G : SimpleGraph V} (q : V → V → ℝ) {u v : V}
    (p : G.Walk u v) :
    walkAffinity q p.reverse = - walkAffinity q p := by
  induction p with
  | nil => simp
  | @cons a b c h p ih =>
      rw [Walk.reverse_cons, walkAffinity_append, ih]
      simp only [walkAffinity_cons, walkAffinity_nil]
      rw [affinity_antisymm q b a]
      ring

/-! ## The three conditions -/

/-- (i) The Kolmogorov cycle condition: every closed walk has zero affinity. -/
def CycleCondition (q : V → V → ℝ) (G : SimpleGraph V) : Prop :=
  ∀ (u : V) (p : G.Walk u u), walkAffinity q p = 0

/-- `K` is a potential for the affinity cochain: `a_{u→v} = K_v − K_u` on every edge
(the coboundary / `head − tail` condition of `GraphDiffusion`). -/
def IsPotential (q : V → V → ℝ) (G : SimpleGraph V) (K : V → ℝ) : Prop :=
  ∀ u v, G.Adj u v → affinity q u v = K v - K u

/-- (ii) A potential exists. -/
def PotentialExists (q : V → V → ℝ) (G : SimpleGraph V) : Prop :=
  ∃ K : V → ℝ, IsPotential q G K

/-- (iii) A strictly positive detailed-balance stationary weight exists:
`π_u q_{u→v} = π_v q_{v→u}` on every edge. -/
def DetailedBalanceExists (q : V → V → ℝ) (G : SimpleGraph V) : Prop :=
  ∃ π : V → ℝ, (∀ x, 0 < π x) ∧ ∀ u v, G.Adj u v → π u * q u v = π v * q v u

/-! ## Telescoping: a potential makes walk affinity a difference of endpoint values -/

/-- If `K` is a potential then the affinity along any walk telescopes to `K(end) − K(start)`. -/
theorem walkAffinity_of_isPotential {G : SimpleGraph V} (q : V → V → ℝ) {K : V → ℝ}
    (hK : IsPotential q G K) {u v : V} (p : G.Walk u v) :
    walkAffinity q p = K v - K u := by
  induction p with
  | nil => simp
  | @cons a b c h p ih =>
      rw [walkAffinity_cons, ih, hK a b h]
      ring

/-! ## (ii) ⇒ (i): potential ⇒ cycle condition (telescoping) -/

theorem cycle_of_potential (q : V → V → ℝ) (G : SimpleGraph V) :
    PotentialExists q G → CycleCondition q G := by
  rintro ⟨K, hK⟩ u p
  rw [walkAffinity_of_isPotential q hK p, sub_self]

/-! ## (ii) ⇒ (iii): potential ⇒ detailed balance (log/exp transform) -/

/-- The literal `π_x = e^{−K_x}` detailed-balance identity of the paper's (iii). -/
theorem detailedBalance_exp_neg (q : V → V → ℝ) (G : SimpleGraph V) {K : V → ℝ}
    (hq : ∀ u v, G.Adj u v → 0 < q u v) (hK : IsPotential q G K)
    (u v : V) (huv : G.Adj u v) :
    Real.exp (-K u) * q u v = Real.exp (-K v) * q v u := by
  have hqp : 0 < q u v := hq u v huv
  have hqm : 0 < q v u := hq v u huv.symm
  have hqvu : q v u = Real.exp (K v - K u) * q u v := by
    have hlog : affinity q u v = K v - K u := hK u v huv
    unfold affinity at hlog
    have hexp : Real.exp (K v - K u) = q v u / q u v := by
      rw [← hlog, Real.exp_log (div_pos hqm hqp)]
    rw [hexp, div_mul_cancel₀ _ hqp.ne']
  rw [hqvu, ← mul_assoc, ← Real.exp_add]
  have harg : -K v + (K v - K u) = -K u := by ring
  rw [harg]

theorem detailedBalance_of_potential (q : V → V → ℝ) (G : SimpleGraph V)
    (hq : ∀ u v, G.Adj u v → 0 < q u v) :
    PotentialExists q G → DetailedBalanceExists q G := by
  rintro ⟨K, hK⟩
  exact ⟨fun x => Real.exp (-K x), fun x => Real.exp_pos _,
    fun u v huv => detailedBalance_exp_neg q G hq hK u v huv⟩

/-! ## (iii) ⇒ (ii): detailed balance ⇒ potential (log transform) -/

theorem potential_of_detailedBalance (q : V → V → ℝ) (G : SimpleGraph V)
    (hq : ∀ u v, G.Adj u v → 0 < q u v) :
    DetailedBalanceExists q G → PotentialExists q G := by
  rintro ⟨π, hπ, hdb⟩
  refine ⟨fun x => - Real.log (π x), ?_⟩
  intro u v huv
  have hqp : 0 < q u v := hq u v huv
  have hpu : 0 < π u := hπ u
  have hpv : 0 < π v := hπ v
  have hratio : q v u / q u v = π u / π v := by
    rw [div_eq_div_iff hqp.ne' hpv.ne']
    linear_combination - hdb u v huv
  show affinity q u v = - Real.log (π v) - - Real.log (π u)
  unfold affinity
  rw [hratio, Real.log_div hpu.ne' hpv.ne']
  ring

/-! ## (i) ⇒ (ii): cycle condition ⇒ potential (discrete Hodge / exactness) -/

theorem potential_of_cycle (q : V → V → ℝ) (G : SimpleGraph V) (hconn : G.Connected) :
    CycleCondition q G → PotentialExists q G := by
  intro hcyc
  obtain ⟨r⟩ := hconn.nonempty
  -- A chosen walk from the root `r` to every vertex (Classical.choice via connectivity).
  set W : (v : V) → G.Walk r v := fun v => Nonempty.some (hconn.preconnected r v) with hW
  -- Path independence: any two walks `r → v` carry equal affinity, because their
  -- concatenation with one reversed is a closed walk, hence has zero affinity.
  have hpi : ∀ (v : V) (p p' : G.Walk r v), walkAffinity q p = walkAffinity q p' := by
    intro v p p'
    have h0 : walkAffinity q (p.append p'.reverse) = 0 := hcyc r _
    rw [walkAffinity_append, walkAffinity_reverse] at h0
    linarith
  refine ⟨fun v => walkAffinity q (W v), ?_⟩
  intro u v huv
  -- The walk `r → u → v` and the chosen walk `r → v` carry equal affinity.
  have key : walkAffinity q ((W u).append (Walk.cons huv Walk.nil))
      = walkAffinity q (W v) := hpi v _ (W v)
  rw [walkAffinity_append, walkAffinity_cons, walkAffinity_nil] at key
  show affinity q u v = walkAffinity q (W v) - walkAffinity q (W u)
  linarith

/-! ## Uniqueness of the potential up to an additive constant (part of the paper's (ii)) -/

/-- A real vertex field that is constant across every edge is constant along every walk. -/
theorem constant_of_adj_constant {G : SimpleGraph V} {f : V → ℝ}
    (hf : ∀ u v, G.Adj u v → f u = f v) {a b : V} (p : G.Walk a b) :
    f a = f b := by
  induction p with
  | nil => rfl
  | cons h p ih => exact (hf _ _ h).trans ih

/-- The potential of the paper's (ii) is unique up to an additive constant. -/
theorem potential_unique (q : V → V → ℝ) (G : SimpleGraph V) (hconn : G.Connected)
    {K K' : V → ℝ} (hK : IsPotential q G K) (hK' : IsPotential q G K') :
    ∃ c : ℝ, ∀ v, K' v = K v + c := by
  obtain ⟨r⟩ := hconn.nonempty
  refine ⟨K' r - K r, fun v => ?_⟩
  have hedge : ∀ a b, G.Adj a b → K' a - K a = K' b - K b := by
    intro a b hab
    have h1 := hK a b hab
    have h2 := hK' a b hab
    linarith
  have hconst : K' r - K r = K' v - K v :=
    constant_of_adj_constant (f := fun w => K' w - K w) hedge
      (Nonempty.some (hconn.preconnected r v))
  linarith

/-! ## The three-way equivalence -/

/-- **Finite rate reconstruction on the reversible layer.**
For a finite connected graph with positive rates in both directions of every edge, the
Kolmogorov cycle condition, existence of a Gibbs potential for the log-rate cochain, and
existence of a positive detailed-balance weight are equivalent. -/
theorem kolmogorov_criterion [Fintype V] (q : V → V → ℝ) (G : SimpleGraph V)
    (hconn : G.Connected) (hq : ∀ u v, G.Adj u v → 0 < q u v) :
    [CycleCondition q G, PotentialExists q G, DetailedBalanceExists q G].TFAE := by
  tfae_have 1 → 2 := potential_of_cycle q G hconn
  tfae_have 2 → 3 := detailedBalance_of_potential q G hq
  tfae_have 3 → 1 := fun h => cycle_of_potential q G (potential_of_detailedBalance q G hq h)
  tfae_finish

/-! ## Non-vacuity witness: the three conditions are genuinely refutable

A concrete finite connected graph (the complete graph `⊤` on `Fin 3`, i.e. a triangle) with
strictly positive rates whose oriented triangle `0 → 1 → 2 → 0` has affinity `log 2 ≠ 0`.
Hence `CycleCondition` fails here — it is not vacuously true — and by `kolmogorov_criterion`
this `q` admits no potential and no detailed-balance weight either.  So all three sides of the
equivalence have real content (none is always true), and the `TFAE` is not degenerate. -/

namespace Witness

/-- Positive rates on the triangle: `1` everywhere except `q 1 0 = 2`, which breaks the cycle. -/
noncomputable def qEx : Fin 3 → Fin 3 → ℝ :=
  fun i j => if (i, j) = ((1 : Fin 3), (0 : Fin 3)) then 2 else 1

theorem qEx_pos (u v : Fin 3) (_ : (⊤ : SimpleGraph (Fin 3)).Adj u v) : 0 < qEx u v := by
  unfold qEx; split <;> norm_num

/-- The oriented triangle `0 → 1 → 2 → 0` as a closed walk.  The intermediate vertices are
pinned via `show` so the `Adj` goals contain no metavariables. -/
def triangle : (⊤ : SimpleGraph (Fin 3)).Walk 0 0 :=
  Walk.cons (show (⊤ : SimpleGraph (Fin 3)).Adj 0 1 by decide)
    (Walk.cons (show (⊤ : SimpleGraph (Fin 3)).Adj 1 2 by decide)
      (Walk.cons (show (⊤ : SimpleGraph (Fin 3)).Adj 2 0 by decide) Walk.nil))

theorem triangle_affinity : walkAffinity qEx triangle = Real.log 2 := by
  unfold triangle
  rw [walkAffinity_cons, walkAffinity_cons, walkAffinity_cons, walkAffinity_nil]
  have e1 : qEx (1 : Fin 3) 0 = 2 := rfl
  have e2 : qEx (0 : Fin 3) 1 = 1 := rfl
  have e3 : qEx (2 : Fin 3) 1 = 1 := rfl
  have e4 : qEx (1 : Fin 3) 2 = 1 := rfl
  have e5 : qEx (0 : Fin 3) 2 = 1 := rfl
  have e6 : qEx (2 : Fin 3) 0 = 1 := rfl
  simp [affinity, e1, e2, e3, e4, e5, e6, div_one, Real.log_one]

theorem cycleCondition_fails : ¬ CycleCondition qEx (⊤ : SimpleGraph (Fin 3)) := by
  intro h
  have h0 : walkAffinity qEx triangle = 0 := h 0 triangle
  rw [triangle_affinity] at h0
  exact (Real.log_pos (by norm_num)).ne' h0

/-- Via the main theorem, this `q` refutes all three conditions at once. -/
theorem all_three_fail :
    ¬ CycleCondition qEx (⊤ : SimpleGraph (Fin 3)) ∧
      ¬ PotentialExists qEx (⊤ : SimpleGraph (Fin 3)) ∧
        ¬ DetailedBalanceExists qEx (⊤ : SimpleGraph (Fin 3)) := by
  have htfae := kolmogorov_criterion qEx (⊤ : SimpleGraph (Fin 3)) connected_top qEx_pos
  have hc := cycleCondition_fails
  exact ⟨hc, fun hp => hc ((htfae.out 1 0).mp hp), fun hd => hc ((htfae.out 2 0).mp hd)⟩

end Witness

#print axioms kolmogorov_criterion
#print axioms potential_unique
#print axioms detailedBalance_exp_neg
#print axioms Witness.cycleCondition_fails
#print axioms Witness.all_three_fail

end OPH.Thermodynamics.KolmogorovCriterion
