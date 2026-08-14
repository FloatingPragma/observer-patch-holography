import Variational.LegendreBridge
import Variational.RealizedHistoryLegendreNoGo
import Variational.StationarySaddleCoverage
import InformationProjection.LogTransitionAction

namespace OPH.Variational

/-!
# The source-to-Hamiltonian composed theorem (V3, issue #731)

One typed antecedent bundle and one conclusion package on the same
objects.  `ComposedMechanicsData` carries a supplied positive normalized
source law, the record embedding with a real two-point enrichment that
reproduces the derived corner action of the same transition law
(register row PR-06), and the derivative, regularity, and solver data of
that enrichment (register row PR-45).  The registered step-uniform path
reference (register row PR-05) enters the conclusion as the reference of
the exponential-tilt clause, matching its register statement.
`ComposedSymmetryData` extends the bundle by the PR-45 symmetry block: a
one-parameter flow with generator that leaves the enrichment invariant.

`source_to_hamiltonian_composed` states, for one bundle `D` and one
history `s` whose embedded path minimizes the enriched local action
under every real fixed-endpoint single-site variation at interior
junctions: the enriched local action of the embedded path equals the
derived log-transition action of the same kernel; the path law of the
same `(pi, P, n)` is the exponential tilt of the registered reference by
that action at multiplier one; the same history is a most probable
history of that path law among interior single-site alphabet variations;
and the embedded path satisfies the discrete Hamilton equations of the
same enrichment at every interior junction.  The boundary conjuncts
carry the committed limits: the realized binary history law does not
select the enrichment (PR-06 stays a register row), and stationarity
alone yields neither a mode nor a minimizer.
`source_to_hamiltonian_composed_noether` adds, for the symmetry-extended
bundle, one constant Noether segment momentum along the whole chain.
`composed_two_faces` is the two-sided clause: fixed-endpoint alphabet
minimality of the enriched action on embedded variations holds exactly
when the same history is a fixed-endpoint single-site mode of the same
path law, with no real-extremality hypothesis.

The real-variation hypothesis is quantified over interior junctions
only.  The committed obstruction of `Variational.FiniteHistoryBridge`
separates alphabet variations from real variations, so real
fixed-endpoint minimality is a hypothesis here and mode-ness is a
consequence by restriction; no clause states or implies that a mode of
the path law satisfies the Hamilton equations, which the committed
stationary-maximum scoping and the enrichment no-go both forbid.
Endpoint sites are excluded because the endpoint variation functional of
the committed curvature family is affine and unbounded below.

A committed instance is included: `chainComposedData` discharges every
field at the two-state mixing chain of the retained locally hash-pinned
run, at every supplied positive curvature and every path length, with
the enrichment family of `RealizedHistoryLegendreNoGo` and its exact
solver.  Its derivative packet `chainCurvedLagrangian_hasFDerivAt` is
proved from the general polynomial packet `polyTwoPoint_hasFDerivAt`.

Boundary.  Finite exact statements under the registered premises.  No
theorem here claims that the source selects the reference, the
enrichment, or the derivative data; physical units, clocks, and
amplitudes stay register rows, and no continuum variational principle is
claimed.
-/

/-- The composed mechanics antecedent bundle.  The first block carries
the supplied source law: a strictly positive normalized initial law, a
strictly positive row-stochastic transition law, and a path length.  The
second block is register row PR-06: the record embedding and a real
two-point enrichment that reproduces the derived corner action
`-log P` at embedded records.  The third block is register row PR-45:
the supplied derivative packet of the enrichment, the strict
monotonicity of its momentum map, and a declared velocity-solver
section.  The registered step-uniform reference (register row PR-05)
appears in the conclusion clauses, not as a field.  The conclusions of
`source_to_hamiltonian_composed` consume: the source block in every
clause; PR-06 in the action identity, the mode clause, and the Hamilton
clause; PR-45 in the Hamilton clause. -/
structure ComposedMechanicsData (Ω : Type*) [Fintype Ω] [Nonempty Ω] where
  /-- The initial law of the supplied transition system. -/
  pi : Ω → ℝ
  /-- The supplied Markov transition law. -/
  P : Ω → Ω → ℝ
  /-- The path length: histories carry `n` transitions. -/
  n : ℕ
  /-- The initial law is strictly positive. -/
  pi_pos : ∀ x, 0 < pi x
  /-- The transition law is strictly positive. -/
  P_pos : ∀ x y, 0 < P x y
  /-- The initial law is normalized. -/
  pi_sum : ∑ x, pi x = 1
  /-- The transition law is row stochastic. -/
  P_rows : ∀ x, ∑ y, P x y = 1
  /-- The record embedding of register row PR-06: alphabet records sit
  at declared real values. -/
  emb : Ω → ℝ
  /-- The real two-point enrichment of register row PR-06. -/
  Lenr : ℝ → ℝ → ℝ
  /-- The PR-06 corner identity: the enrichment reproduces the derived
  corner action of the same supplied kernel at embedded records. -/
  corner : ∀ x y, Lenr (emb x) (emb y) = -Real.log (P x y)
  /-- The first-slot derivative of the enrichment (register row
  PR-45). -/
  d1 : ℝ → ℝ → ℝ
  /-- The second-slot (momentum) derivative of the enrichment (register
  row PR-45). -/
  d2 : ℝ → ℝ → ℝ
  /-- The declared velocity solver of the enrichment (register row
  PR-45). -/
  vel : ℝ → ℝ → ℝ
  /-- The supplied uncurried derivative packet of the enrichment
  (register row PR-45). -/
  deriv_packet : ∀ p : ℝ × ℝ, HasFDerivAt (Function.uncurry Lenr)
    ((d1 p.1 p.2) • (ContinuousLinearMap.fst ℝ ℝ ℝ)
      + (d2 p.1 p.2) • (ContinuousLinearMap.snd ℝ ℝ ℝ)) p
  /-- The momentum map is strictly monotone in the velocity slot: the
  enrichment is regular (register row PR-45). -/
  mono : ∀ x, StrictMono fun v => d2 x v
  /-- The velocity solver is a section of the momentum map (register
  row PR-45). -/
  sec : ∀ x p, d2 x (vel x p) = p

/-- The symmetry-extended bundle: the PR-45 symmetry block.  A
one-parameter transformation family with generator `xi` that fixes every
point at parameter zero, is differentiable in the parameter, and leaves
the PR-06 enrichment exactly invariant. -/
structure ComposedSymmetryData (Ω : Type*) [Fintype Ω] [Nonempty Ω]
    extends ComposedMechanicsData Ω where
  /-- The one-parameter transformation family (register row PR-45). -/
  T : ℝ → ℝ → ℝ
  /-- The generator of the family (register row PR-45). -/
  xi : ℝ → ℝ
  /-- The family fixes every point at parameter zero. -/
  T0 : ∀ x, T 0 x = x
  /-- The family is differentiable in the parameter with derivative the
  generator. -/
  Td : ∀ x, HasDerivAt (fun s => T s x) (xi x) 0
  /-- The family leaves the PR-06 enrichment exactly invariant. -/
  inv : ∀ s x y, Lenr (T s x) (T s y) = Lenr x y

/-- The two committed spellings of the log-transition action agree: the
negated sum of `InformationProjection.LogTransitionAction` is the sum of
negations of `Variational.LegendreBridge`, on every kernel and path. -/
theorem informationProjection_logTransitionAction_eq {Ω : Type*} {n : ℕ}
    (P : Ω → Ω → ℝ) (t : Fin (n + 1) → Ω) :
    InformationProjection.logTransitionAction P n t
      = logTransitionAction P t := by
  unfold InformationProjection.logTransitionAction logTransitionAction
  rw [Finset.sum_neg_distrib]

/-- The Markov path law factors as initial mass times the conditional
transition weight of the committed variational stack. -/
theorem markovPathLaw_eq_pi_mul_pathTransitionWeight {Ω : Type*}
    [Fintype Ω] (pi : Ω → ℝ) (P : Ω → Ω → ℝ) (n : ℕ)
    (t : Fin (n + 1) → Ω) :
    InformationProjection.markovPathLaw pi P n t
      = pi (t 0) * pathTransitionWeight P t := rfl

namespace ComposedMechanicsData

variable {Ω : Type*} [Fintype Ω] [Nonempty Ω] (D : ComposedMechanicsData Ω)

/-- **The glue identity (PR-06 consumed).**  On every embedded history
of the bundle, the enriched local action equals the derived
log-transition action of the same supplied kernel. -/
theorem localAction_emb (t : Fin (D.n + 1) → Ω) :
    localAction D.Lenr (D.emb ∘ t)
      = InformationProjection.logTransitionAction D.P D.n t := by
  rw [informationProjection_logTransitionAction_eq]
  unfold localAction logTransitionAction
  refine Finset.sum_congr rfl fun i _ => ?_
  simp only [Function.comp_apply]
  exact D.corner (t i.castSucc) (t i.succ)

/-- **Interior mode from real minimality.**  A history whose embedded
path minimizes the enriched local action under every real single-site
variation at interior junctions is a most probable history of the same
path law among interior single-site alphabet variations.  The alphabet
variations are restrictions of the real variations through the PR-06
embedding, and interior updates never touch the initial record. -/
theorem interior_mode_of_min (s : Fin (D.n + 1) → Ω)
    (hmin_real : ∀ (k m : Fin D.n), k.succ = m.castSucc → ∀ x : ℝ,
      localAction D.Lenr (D.emb ∘ s)
        ≤ localAction D.Lenr (Function.update (D.emb ∘ s) k.succ x)) :
    ∀ (k m : Fin D.n), k.succ = m.castSucc → ∀ x : Ω,
      InformationProjection.markovPathLaw D.pi D.P D.n
          (Function.update s k.succ x)
        ≤ InformationProjection.markovPathLaw D.pi D.P D.n s := by
  intro k m hkm x
  have h1 := hmin_real k m hkm (D.emb x)
  rw [← Function.comp_update] at h1
  rw [D.localAction_emb s, D.localAction_emb (Function.update s k.succ x),
    informationProjection_logTransitionAction_eq,
    informationProjection_logTransitionAction_eq] at h1
  have h2 := (logTransitionAction_le_iff D.P D.P_pos s
    (Function.update s k.succ x)).mp h1
  have h0 : Function.update s k.succ x 0 = s 0 :=
    Function.update_of_ne (Fin.succ_ne_zero k).symm _ _
  rw [markovPathLaw_eq_pi_mul_pathTransitionWeight,
    markovPathLaw_eq_pi_mul_pathTransitionWeight, h0]
  exact mul_le_mul_of_nonneg_left h2 (le_of_lt (D.pi_pos (s 0)))

/-- **The two-faces clause on the bundle.**  Fixed-endpoint alphabet
minimality of the enriched local action on embedded variations holds
exactly when the same history is a fixed-endpoint single-site mode of
the path law of the same supplied system.  Both directions; no
real-extremality hypothesis enters. -/
theorem composed_two_faces (s : Fin (D.n + 1) → Ω) :
    (∀ (k m : Fin D.n), k.succ = m.castSucc → ∀ x : Ω,
        localAction D.Lenr (D.emb ∘ s)
          ≤ localAction D.Lenr (D.emb ∘ Function.update s k.succ x))
      ↔ ∀ (k m : Fin D.n), k.succ = m.castSucc → ∀ x : Ω,
          InformationProjection.markovPathLaw D.pi D.P D.n
              (Function.update s k.succ x)
            ≤ InformationProjection.markovPathLaw D.pi D.P D.n s := by
  refine forall_congr' fun k => forall_congr' fun m => ?_
  refine imp_congr_right fun hkm => forall_congr' fun x => ?_
  rw [D.localAction_emb s, D.localAction_emb (Function.update s k.succ x),
    informationProjection_logTransitionAction_eq,
    informationProjection_logTransitionAction_eq,
    markovPathLaw_eq_pi_mul_pathTransitionWeight,
    markovPathLaw_eq_pi_mul_pathTransitionWeight,
    Function.update_of_ne (Fin.succ_ne_zero k).symm x s]
  constructor
  · intro h
    exact mul_le_mul_of_nonneg_left
      ((logTransitionAction_le_iff D.P D.P_pos s
        (Function.update s k.succ x)).mp h)
      (le_of_lt (D.pi_pos (s 0)))
  · intro h
    exact (logTransitionAction_le_iff D.P D.P_pos s
      (Function.update s k.succ x)).mpr
      (le_of_mul_le_mul_left h (D.pi_pos (s 0)))

end ComposedMechanicsData

/-- **The source-to-Hamiltonian composed theorem (issue #731, rows
PR-05, PR-06, PR-45).**  For one bundle `D` and one history `s` whose
embedded path minimizes the enriched local action under every real
fixed-endpoint single-site variation at interior junctions: the
enriched local action of the embedded path equals the derived
log-transition action of the same kernel (PR-06 consumed); the path law
of the same `(pi, P, n)` is the exponential tilt of the registered
step-uniform reference (row PR-05) by that action at multiplier one;
the same history is a most probable history of the same path law among
interior single-site alphabet variations; and the embedded path
satisfies the discrete Hamilton equations of the same enrichment at
every interior junction (PR-45 consumed).  The boundary conjuncts carry
the committed limits: the canonical bilinear extension of the realized
action has no velocity solver and two realized-history-equivalent
enrichments differ, so PR-06 stays a register row; and a stationary
history can fail minimality while the Gibbs weight prefers a variation
at every positive multiplier, so stationarity alone yields neither a
mode nor a minimizer.  The real-minimality hypothesis quantifies over
interior junctions only, and no clause states the converse direction
from mode to Hamilton equations. -/
theorem source_to_hamiltonian_composed {Ω : Type*} [Fintype Ω]
    [Nonempty Ω] (D : ComposedMechanicsData Ω) (s : Fin (D.n + 1) → Ω)
    (hmin_real : ∀ (k m : Fin D.n), k.succ = m.castSucc → ∀ x : ℝ,
      localAction D.Lenr (D.emb ∘ s)
        ≤ localAction D.Lenr (Function.update (D.emb ∘ s) k.succ x)) :
    localAction D.Lenr (D.emb ∘ s)
        = InformationProjection.logTransitionAction D.P D.n s
      ∧ InformationProjection.markovPathLaw D.pi D.P D.n
          = InformationProjection.tilt
              (InformationProjection.stepUniformRef D.pi D.n)
              (InformationProjection.logTransitionAction D.P D.n) 1
      ∧ (∀ (k m : Fin D.n), k.succ = m.castSucc → ∀ x : Ω,
          InformationProjection.markovPathLaw D.pi D.P D.n
              (Function.update s k.succ x)
            ≤ InformationProjection.markovPathLaw D.pi D.P D.n s)
      ∧ (∀ (k m : Fin D.n), k.succ = m.castSucc →
          (D.emb ∘ s) m.castSucc
              = D.vel ((D.emb ∘ s) k.castSucc)
                  (D.d2 ((D.emb ∘ s) k.castSucc) ((D.emb ∘ s) k.succ))
            ∧ D.d2 ((D.emb ∘ s) k.castSucc) ((D.emb ∘ s) k.succ)
                + D.d1 ((D.emb ∘ s) m.castSucc)
                    (D.vel ((D.emb ∘ s) m.castSucc)
                      (D.d2 ((D.emb ∘ s) m.castSucc)
                        ((D.emb ∘ s) m.succ))) = 0)
      ∧ (¬ ∃ vel : ℝ → ℝ → ℝ, SolvesMomentum chainLogLagrangian vel)
      ∧ chainCurvedLagrangian 1 ≠ chainCurvedLagrangian 2
      ∧ (¬ ∀ x, localAction concaveControlL stationaryMaximumHistory
          ≤ localAction concaveControlL
              (Function.update stationaryMaximumHistory
                ((0 : Fin 2).succ) x))
      ∧ ∀ lam : ℝ, 0 < lam →
          Real.exp (-lam
              * localAction concaveControlL stationaryMaximumHistory)
            < Real.exp (-lam * localAction concaveControlL
                (Function.update stationaryMaximumHistory
                  ((0 : Fin 2).succ) 1)) :=
  ⟨D.localAction_emb s,
    InformationProjection.markov_path_law_eq_gibbs D.pi D.P D.n
      D.P_pos D.pi_sum D.P_rows,
    D.interior_mode_of_min s hmin_real,
    fun k m hkm => stationary_hamilton_strictConvex D.Lenr D.d1 D.d2
      D.vel D.deriv_packet D.mono D.sec (D.emb ∘ s) hkm
      (hmin_real k m hkm),
    chainLogLagrangian_no_velocity_solver,
    chainCurvedLagrangian_one_ne_two,
    stationaryMaximumHistory_not_minimal,
    fun _ hlam => gibbs_prefers_nonstationary hlam⟩

/-- Noether constancy from junction-quantified minimality: the
chain-wide theorem of `DiscreteNoether` restated over the interior
junction quantifier of the composed bundle.  At length zero the segment
set is empty and any scalar serves. -/
theorem noether_constant_of_junction_min {n : ℕ}
    (L : ℝ → ℝ → ℝ) (T : ℝ → ℝ → ℝ) (ξ : ℝ → ℝ)
    (d₁ d₂ : ℝ → ℝ → ℝ)
    (hT0 : ∀ x, T 0 x = x)
    (hTd : ∀ x, HasDerivAt (fun s => T s x) (ξ x) 0)
    (hL : ∀ p : ℝ × ℝ, HasFDerivAt (Function.uncurry L)
      ((d₁ p.1 p.2) • (ContinuousLinearMap.fst ℝ ℝ ℝ)
        + (d₂ p.1 p.2) • (ContinuousLinearMap.snd ℝ ℝ ℝ)) p)
    (hinv : ∀ s x y, L (T s x) (T s y) = L x y)
    (γ : Fin (n + 1) → ℝ)
    (hmin : ∀ (k m : Fin n), k.succ = m.castSucc → ∀ x : ℝ,
      localAction L γ ≤ localAction L (Function.update γ k.succ x)) :
    ∃ J : ℝ, ∀ seg : Fin n, segmentMomentum d₂ ξ γ seg = J := by
  cases n with
  | zero => exact ⟨0, fun seg => seg.elim0⟩
  | succ M =>
      exact noether_current_constant_on_finite_chain L T ξ d₁ d₂
        hT0 hTd hL hinv γ (fun i x => hmin i.castSucc i.succ rfl x)

/-- **The composed Noether clause (row PR-45 symmetry block).**  For the
symmetry-extended bundle and a history whose embedded path is a real
fixed-endpoint extremizer at every interior junction, the segment
momenta of the same embedded path are one constant along the whole
finite chain.  The scalar is conditional on the supplied symmetry data;
no physical charge, clock, or continuum current is claimed. -/
theorem source_to_hamiltonian_composed_noether {Ω : Type*} [Fintype Ω]
    [Nonempty Ω] (E : ComposedSymmetryData Ω) (s : Fin (E.n + 1) → Ω)
    (hmin_real : ∀ (k m : Fin E.n), k.succ = m.castSucc → ∀ x : ℝ,
      localAction E.Lenr (E.emb ∘ s)
        ≤ localAction E.Lenr (Function.update (E.emb ∘ s) k.succ x)) :
    ∃ J : ℝ, ∀ seg : Fin E.n,
      segmentMomentum E.d2 E.xi (E.emb ∘ s) seg = J :=
  noether_constant_of_junction_min E.Lenr E.T E.xi E.d1 E.d2 E.T0 E.Td
    E.deriv_packet E.inv (E.emb ∘ s) hmin_real

/-! ## The general polynomial derivative packet -/

/-- The uncurried derivative packet of a two-point polynomial: for the
function `α + β x + γ y + δ x y + ε x² + ζ y²` on the product plane, the
Fréchet derivative at `p` is the stated combination of the coordinate
projections.  This is the one analysis lemma the composed instances
consume; the committed enrichment families are all of this shape. -/
theorem polyTwoPoint_hasFDerivAt (α β γ δ ε ζ : ℝ) (p : ℝ × ℝ) :
    HasFDerivAt
      (fun q : ℝ × ℝ => α + β * q.1 + γ * q.2 + δ * (q.1 * q.2)
        + ε * (q.1 * q.1) + ζ * (q.2 * q.2))
      ((β + δ * p.2 + 2 * ε * p.1) • ContinuousLinearMap.fst ℝ ℝ ℝ
        + (γ + δ * p.1 + 2 * ζ * p.2) • ContinuousLinearMap.snd ℝ ℝ ℝ)
      p := by
  have hx : HasFDerivAt (fun q : ℝ × ℝ => q.1)
      (ContinuousLinearMap.fst ℝ ℝ ℝ) p := hasFDerivAt_fst
  have hy : HasFDerivAt (fun q : ℝ × ℝ => q.2)
      (ContinuousLinearMap.snd ℝ ℝ ℝ) p := hasFDerivAt_snd
  have h : HasFDerivAt
      (fun q : ℝ × ℝ => α + β * q.1 + γ * q.2 + δ * (q.1 * q.2)
        + ε * (q.1 * q.1) + ζ * (q.2 * q.2))
      ((0 : ℝ × ℝ →L[ℝ] ℝ)
        + β • ContinuousLinearMap.fst ℝ ℝ ℝ
        + γ • ContinuousLinearMap.snd ℝ ℝ ℝ
        + δ • (p.1 • ContinuousLinearMap.snd ℝ ℝ ℝ
            + p.2 • ContinuousLinearMap.fst ℝ ℝ ℝ)
        + ε • (p.1 • ContinuousLinearMap.fst ℝ ℝ ℝ
            + p.1 • ContinuousLinearMap.fst ℝ ℝ ℝ)
        + ζ • (p.2 • ContinuousLinearMap.snd ℝ ℝ ℝ
            + p.2 • ContinuousLinearMap.snd ℝ ℝ ℝ)) p :=
    (((((hasFDerivAt_const α p).add (hx.const_mul β)).add
      (hy.const_mul γ)).add ((hx.mul hy).const_mul δ)).add
      ((hx.mul hx).const_mul ε)).add ((hy.mul hy).const_mul ζ)
  have hD : (β + δ * p.2 + 2 * ε * p.1) • ContinuousLinearMap.fst ℝ ℝ ℝ
        + (γ + δ * p.1 + 2 * ζ * p.2) • ContinuousLinearMap.snd ℝ ℝ ℝ
      = (0 : ℝ × ℝ →L[ℝ] ℝ)
        + β • ContinuousLinearMap.fst ℝ ℝ ℝ
        + γ • ContinuousLinearMap.snd ℝ ℝ ℝ
        + δ • (p.1 • ContinuousLinearMap.snd ℝ ℝ ℝ
            + p.2 • ContinuousLinearMap.fst ℝ ℝ ℝ)
        + ε • (p.1 • ContinuousLinearMap.fst ℝ ℝ ℝ
            + p.1 • ContinuousLinearMap.fst ℝ ℝ ℝ)
        + ζ • (p.2 • ContinuousLinearMap.snd ℝ ℝ ℝ
            + p.2 • ContinuousLinearMap.snd ℝ ℝ ℝ) := by
    refine ContinuousLinearMap.ext fun v => ?_
    simp only [ContinuousLinearMap.add_apply,
      ContinuousLinearMap.smul_apply, ContinuousLinearMap.coe_fst',
      ContinuousLinearMap.coe_snd', ContinuousLinearMap.zero_apply,
      smul_eq_mul]
    ring
  rw [hD]
  exact h

/-! ## The committed instance at the mixing chain -/

/-- The first-slot slope of the committed curvature family: the
derivative of `chainLogBase`. -/
noncomputable def chainBaseSlope : ℝ :=
  Real.log (chainWeight 0 0) - Real.log (chainWeight 1 0)

/-- The rate of the fiber slope of the committed curvature family: the
derivative of `chainFiberSlope`. -/
noncomputable def chainSlopeRate : ℝ :=
  (Real.log (chainWeight 1 0) - Real.log (chainWeight 1 1))
    - (Real.log (chainWeight 0 0) - Real.log (chainWeight 0 1))

/-- The first-slot derivative of the committed curvature family.  The
family is affine in its first slot, so the value does not depend on the
first argument or on the curvature. -/
noncomputable def chainCurvedD1 (_x y : ℝ) : ℝ :=
  chainBaseSlope + chainSlopeRate * y

/-- The second-slot (momentum) derivative of the committed curvature
family: the value certified by `chainCurvedLagrangian_momentum`. -/
noncomputable def chainCurvedD2 (a x y : ℝ) : ℝ :=
  chainFiberSlope x + a * y - a / 2

/-- The uncurried derivative packet of the committed curvature family,
from the general polynomial packet. -/
theorem chainCurvedLagrangian_hasFDerivAt (a : ℝ) (p : ℝ × ℝ) :
    HasFDerivAt (Function.uncurry (chainCurvedLagrangian a))
      (chainCurvedD1 p.1 p.2 • ContinuousLinearMap.fst ℝ ℝ ℝ
        + chainCurvedD2 a p.1 p.2 • ContinuousLinearMap.snd ℝ ℝ ℝ)
      p := by
  have hfun : Function.uncurry (chainCurvedLagrangian a)
      = fun q : ℝ × ℝ =>
          -Real.log (chainWeight 0 0)
            + chainBaseSlope * q.1
            + (Real.log (chainWeight 0 0) - Real.log (chainWeight 0 1)
                - a / 2) * q.2
            + chainSlopeRate * (q.1 * q.2)
            + 0 * (q.1 * q.1)
            + (a / 2) * (q.2 * q.2) := by
    funext q
    simp only [Function.uncurry, chainCurvedLagrangian_eq, chainLogBase,
      chainFiberSlope, chainBaseSlope, chainSlopeRate]
    ring
  have h := polyTwoPoint_hasFDerivAt (-Real.log (chainWeight 0 0))
    chainBaseSlope
    (Real.log (chainWeight 0 0) - Real.log (chainWeight 0 1) - a / 2)
    chainSlopeRate 0 (a / 2) p
  have h1 : chainCurvedD1 p.1 p.2
      = chainBaseSlope + chainSlopeRate * p.2 + 2 * 0 * p.1 := by
    simp only [chainCurvedD1]
    ring
  have h2 : chainCurvedD2 a p.1 p.2
      = Real.log (chainWeight 0 0) - Real.log (chainWeight 0 1) - a / 2
        + chainSlopeRate * p.1 + 2 * (a / 2) * p.2 := by
    simp only [chainCurvedD2, chainFiberSlope, chainSlopeRate]
    ring
  rw [hfun, h1, h2]
  exact h

/-- The momentum map of the committed curvature family is strictly
monotone in the velocity slot at every positive curvature. -/
theorem chainCurvedD2_strictMono (a : ℝ) (ha : 0 < a) (x : ℝ) :
    StrictMono fun v => chainCurvedD2 a x v := by
  intro u v huv
  simp only [chainCurvedD2]
  have h := mul_lt_mul_of_pos_left huv ha
  linarith

/-- The committed exact solver is a section of the momentum map of the
committed curvature family at every positive curvature. -/
theorem chainCurvedD2_solver_section (a : ℝ) (ha : 0 < a) (x p : ℝ) :
    chainCurvedD2 a x (chainCurvedVelocitySolver a x p) = p := by
  have ha' : a ≠ 0 := ne_of_gt ha
  simp only [chainCurvedD2, chainCurvedVelocitySolver]
  field_simp
  ring

/-- The cast committed chain is the real chain weight table. -/
theorem sourcePR_eq_chainWeight :
    InformationProjection.sourcePR = chainWeight := rfl

/-- **The committed bundle instance.**  The two-state mixing chain of
the retained locally hash-pinned run, at every supplied positive
curvature and every path length, discharges every field of the composed
bundle: the cast stationary law and chain, the record embedding at the
alphabet values, the curvature-family enrichment of
`RealizedHistoryLegendreNoGo` with its corner identity, and the
derivative, monotonicity, and solver data of that family.  The
curvature is a supplied field, matching the PR-06 no-go: the realized
history law does not select it. -/
noncomputable def chainComposedData (a : ℝ) (ha : 0 < a) (M : ℕ) :
    ComposedMechanicsData (Fin 2) where
  pi := InformationProjection.sourcePiR
  P := InformationProjection.sourcePR
  n := M
  pi_pos := fun i => by
    simp only [InformationProjection.sourcePiR]
    exact_mod_cast OPH.Thermodynamics.mixingChainStationary_pos i
  P_pos := InformationProjection.sourcePR_pos
  pi_sum := InformationProjection.sourcePiR_sum
  P_rows := InformationProjection.sourcePR_rows
  emb := fun i => ((i : ℕ) : ℝ)
  Lenr := chainCurvedLagrangian a
  corner := fun i j => by
    rw [sourcePR_eq_chainWeight]
    exact chainCurvedLagrangian_corner a i j
  d1 := chainCurvedD1
  d2 := chainCurvedD2 a
  vel := chainCurvedVelocitySolver a
  deriv_packet := chainCurvedLagrangian_hasFDerivAt a
  mono := chainCurvedD2_strictMono a ha
  sec := chainCurvedD2_solver_section a ha

end OPH.Variational

#print axioms OPH.Variational.informationProjection_logTransitionAction_eq
#print axioms OPH.Variational.markovPathLaw_eq_pi_mul_pathTransitionWeight
#print axioms OPH.Variational.ComposedMechanicsData.localAction_emb
#print axioms OPH.Variational.ComposedMechanicsData.interior_mode_of_min
#print axioms OPH.Variational.ComposedMechanicsData.composed_two_faces
#print axioms OPH.Variational.source_to_hamiltonian_composed
#print axioms OPH.Variational.noether_constant_of_junction_min
#print axioms OPH.Variational.source_to_hamiltonian_composed_noether
#print axioms OPH.Variational.polyTwoPoint_hasFDerivAt
#print axioms OPH.Variational.chainCurvedLagrangian_hasFDerivAt
#print axioms OPH.Variational.chainCurvedD2_strictMono
#print axioms OPH.Variational.chainCurvedD2_solver_section
#print axioms OPH.Variational.sourcePR_eq_chainWeight
