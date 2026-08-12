import Variational.LegendreBridge
import Variational.RealizedHistoryLegendreNoGo
import Variational.StationarySaddleCoverage
import InformationProjection.LogTransitionAction

namespace OPH.Variational

/-!
# The composed mechanics adequacy surface (V3, issue #731)

One citable package for the observation-ledger rows OL-D1, OL-D2, and
OL-D3: least/stationary action, the discrete Euler-Lagrange and Hamilton
equivalence, and Noether conservation, composed from the committed
mechanics modules with every supplied object named.  This module adds no
mathematical content; every theorem re-exports a committed result
through the premise structure below.

`MechanicsPremiseData` bundles the supplied objects the committed
theorems consume: a normalized strictly positive initial law, a strictly
positive row-stochastic transition law, and a path length (the inputs of
the derived-action theorem of
`InformationProjection.LogTransitionAction`), together with a supplied
regular quadratic Lagrangian given by a strictly positive stiffness, a
potential, and its derivative (the inputs of the Legendre bridge of
`Variational.LegendreBridge`).  Two of the composed objects are premise
register rows.  The step-uniform path reference of the derived-action
theorem is register row PR-05, and the real velocity enrichment carried
by the supplied stiffness is register row PR-06.  The surface is the
composition of the committed theorems over those rows, and no theorem
here claims that the source selects either row.

The two premise families do not compose into a single derivation from
the transition law to the Hamiltonian face, and the surface records the
obstruction as its own boundary.  By the committed no-go of
`Variational.RealizedHistoryLegendreNoGo`, every real curvature
enrichment reproduces the realized binary action, so the realized
history law does not determine the regular Lagrangian that the Legendre
bridge consumes.  The derived-action theorems and the Legendre-bridge
theorems therefore stand side by side inside the
structure-parameterized namespace, with the enrichment a declared
register row.

Contents.  On the derived-action side: the Markov path law is the
exponential tilt of the step-uniform reference by the log-transition
action at multiplier one, the reproducing action-multiplier pairs form
exactly the additive-constant gauge class, the partition constant is
computed, and the bare multiplier is forced to one under nonconstancy
of the action.  On the Legendre side: the finite Legendre transform of
the supplied Lagrangian with involutivity and Fenchel-Young, the
junction equivalence of the discrete Euler-Lagrange equation with one
step of the discrete Hamilton flow, and the composed forms in which
single-site minimality yields the discrete Euler-Lagrange equation and
the Hamilton step.  On the Noether side: the discrete Hamilton flow of
the free class at the supplied stiffness conserves the momentum and
Hamiltonian readouts, and the committed translation witness carries the
constant chain current with exact conserved values `2` and `1`.  The
realized two-faces theorem identifies single-site minimizers of the
committed-form log-transition action with the most probable paths among
single-site variations of the realized chain at every length.  Boundary corollaries carry the
committed limits: the realized binary history law does not select the
enrichment, and a stationary history can fail minimality while the
Gibbs weight prefers a variation at every positive multiplier, so the
mode/minimizer reading of stationarity stays scoped.

A committed instance of the premise data is included: the two-state
mixing chain of the retained locally hash-pinned run, with stiffness
`2` and zero potential, discharges every field with exact literals.

Boundary.  Finite exact statements under the registered premises.
Physical units, clocks, and amplitudes stay register rows, and no
continuum variational principle is claimed.
-/

/-- The supplied objects of the composed mechanics surface.  The first
block carries the derived-action inputs: a normalized strictly positive
initial law `pi`, a strictly positive row-stochastic transition law `P`,
and a path length `n`.  The second block carries the Legendre-bridge
inputs: a supplied regular quadratic Lagrangian with strictly positive
stiffness `a`, potential `V`, and derivative `V'`.  The stiffness is the
velocity enrichment of register row PR-06; the theorems below consume it
as a declared object, and the no-go corollary states why. -/
structure MechanicsPremiseData (Ω : Type*) [Fintype Ω] where
  /-- The initial law of the supplied transition system. -/
  pi : Ω → ℝ
  /-- The supplied Markov transition law. -/
  P : Ω → Ω → ℝ
  /-- The path length: histories carry `n` transitions. -/
  n : ℕ
  /-- The stiffness of the supplied regular Lagrangian: the velocity
  enrichment of register row PR-06. -/
  a : ℝ
  /-- The potential of the supplied regular Lagrangian. -/
  V : ℝ → ℝ
  /-- The supplied derivative of the potential. -/
  V' : ℝ → ℝ
  /-- The initial law is strictly positive. -/
  pi_pos : ∀ x, 0 < pi x
  /-- The transition law is strictly positive. -/
  P_pos : ∀ x y, 0 < P x y
  /-- The initial law is normalized. -/
  pi_sum : ∑ x, pi x = 1
  /-- The transition law is row stochastic. -/
  P_rows : ∀ x, ∑ y, P x y = 1
  /-- The stiffness is strictly positive: the supplied Lagrangian is
  regular. -/
  a_pos : 0 < a
  /-- The potential derivative packet. -/
  potential_hasDerivAt : ∀ t, HasDerivAt V (V' t) t

namespace MechanicsPremiseData

variable {Ω : Type*} [Fintype Ω] (D : MechanicsPremiseData Ω)

/-! ## The derived action over the registered reference (row PR-05) -/

/-- **The derived action (OL-D1 input).**  The path law of the supplied
transition system is the exponential tilt of the step-uniform reference
(register row PR-05) by the log-transition action at multiplier one.
Re-export of `InformationProjection.markov_path_law_eq_gibbs` at the
bundled premises. -/
theorem path_law_eq_exponential_tilt [Nonempty Ω] :
    InformationProjection.markovPathLaw D.pi D.P D.n
      = InformationProjection.tilt
          (InformationProjection.stepUniformRef D.pi D.n)
          (InformationProjection.logTransitionAction D.P D.n) 1 :=
  InformationProjection.markov_path_law_eq_gibbs D.pi D.P D.n
    D.P_pos D.pi_sum D.P_rows

/-- The partition constant of the derived action is computed exactly:
the uniform-step counting weight of the path space.  Re-export of
`InformationProjection.markov_tiltZ_eq`. -/
theorem partition_constant [Nonempty Ω] :
    InformationProjection.tiltZ
        (InformationProjection.stepUniformRef D.pi D.n)
        (InformationProjection.logTransitionAction D.P D.n) 1
      = ((Fintype.card Ω : ℝ))⁻¹ ^ D.n :=
  InformationProjection.markov_tiltZ_eq D.pi D.P D.n
    D.P_pos D.pi_sum D.P_rows

/-- The derived action is unique up to gauge: an action-multiplier pair
reproduces the path law over the registered reference exactly when its
multiplier-weighted action is the log-transition action plus a constant.
Re-export of `InformationProjection.action_unique_up_to_gauge`. -/
theorem derived_action_gauge_class [Nonempty Ω]
    (S' : InformationProjection.PathSpace Ω D.n → ℝ) (lam' : ℝ) :
    InformationProjection.tilt
          (InformationProjection.stepUniformRef D.pi D.n) S' lam'
        = InformationProjection.markovPathLaw D.pi D.P D.n
      ↔ ∃ c : ℝ, ∀ γ : InformationProjection.PathSpace Ω D.n,
          lam' * S' γ
            = InformationProjection.logTransitionAction D.P D.n γ + c :=
  InformationProjection.action_unique_up_to_gauge D.pi D.P D.n
    D.pi_pos D.P_pos D.pi_sum D.P_rows S' lam'

/-- The bare log-transition action forces multiplier one whenever the
action is nonconstant on the path space.  Re-export of
`InformationProjection.bare_log_action_multiplier_unique_of_nonconstant`. -/
theorem bare_multiplier_forced [Nonempty Ω]
    (hnonconst : ∃ γ δ : InformationProjection.PathSpace Ω D.n,
      InformationProjection.logTransitionAction D.P D.n γ
        ≠ InformationProjection.logTransitionAction D.P D.n δ)
    (lam : ℝ)
    (hfit : InformationProjection.tilt
        (InformationProjection.stepUniformRef D.pi D.n)
        (InformationProjection.logTransitionAction D.P D.n) lam
      = InformationProjection.markovPathLaw D.pi D.P D.n) :
    lam = 1 :=
  InformationProjection.bare_log_action_multiplier_unique_of_nonconstant
    D.pi D.P D.n D.pi_pos D.P_pos D.pi_sum D.P_rows hnonconst lam hfit

/-! ## The Legendre bridge for the supplied regular Lagrangian (row PR-06) -/

/-- The finite Legendre transform of the supplied regular Lagrangian is
the quadratic Hamiltonian.  Re-export of `quad_legendreTransform`. -/
theorem legendre_transform_eq :
    legendreTransform (quadLagrangian D.a D.V) (quadVelocitySolver D.a)
      = quadHamiltonian D.a D.V :=
  quad_legendreTransform D.a D.V (ne_of_gt D.a_pos)

/-- Applying the transform twice returns the supplied Lagrangian.
Re-export of `quad_legendre_involutive`. -/
theorem legendre_involutive :
    legendreTransform (quadHamiltonian D.a D.V) (quadVelocitySolver D.a⁻¹)
      = quadLagrangian D.a D.V :=
  quad_legendre_involutive D.a D.V (ne_of_gt D.a_pos)

/-- The exact velocity solver of the supplied Lagrangian solves its
momentum relation.  Re-export of `quad_solver_solves`. -/
theorem velocity_solver_solves :
    SolvesMomentum (quadLagrangian D.a D.V) (quadVelocitySolver D.a) :=
  quad_solver_solves D.a D.V (ne_of_gt D.a_pos)

/-- The supplied Lagrangian is strictly convex in the velocity slot.
Re-export of `quadLagrangian_strictConvexOn`. -/
theorem lagrangian_strictConvexOn (x : ℝ) :
    StrictConvexOn ℝ Set.univ fun v : ℝ => quadLagrangian D.a D.V x v :=
  quadLagrangian_strictConvexOn D.a D.V D.a_pos x

/-- Fenchel-Young for the supplied Lagrangian.  Re-export of
`quad_fenchel_young`. -/
theorem fenchel_young (x p w : ℝ) :
    p * w - quadLagrangian D.a D.V x w ≤ quadHamiltonian D.a D.V x p :=
  quad_fenchel_young D.a D.V D.a_pos x p w

/-- **The Euler-Lagrange/Hamilton equivalence (OL-D2 input).**  At a
junction of the committed stack, the discrete Euler-Lagrange sum of the
supplied Lagrangian vanishes exactly when one step of the discrete
Hamilton flow carries the incoming junction state to the outgoing one.
Re-export of `quad_discreteEulerLagrange_iff_hamiltonStep`. -/
theorem eulerLagrange_iff_hamiltonStep {N : ℕ} (γ : Fin (N + 1) → ℝ)
    {k m : Fin N} (hkm : k.succ = m.castSucc) :
    quadSegmentMomentum D.a γ k
        + (-(quadSegmentMomentum D.a γ m) - D.V' (γ m.castSucc)) = 0
      ↔ quadHamiltonStep D.a D.V'
            (γ k.castSucc, quadSegmentMomentum D.a γ k)
          = (γ m.castSucc, quadSegmentMomentum D.a γ m) :=
  quad_discreteEulerLagrange_iff_hamiltonStep D.a (ne_of_gt D.a_pos)
    D.V' γ hkm

/-- **Stationarity yields the Hamilton step.**  A history that minimizes
the local action of the supplied Lagrangian among single-site variations
at a junction satisfies the discrete Hamilton equations there.
Re-export of `quad_stationary_hamiltonStep`. -/
theorem stationary_hamiltonStep {N : ℕ} (γ : Fin (N + 1) → ℝ)
    {k m : Fin N} (hkm : k.succ = m.castSucc)
    (hmin : ∀ x, localAction (quadTwoPoint D.a D.V) γ
      ≤ localAction (quadTwoPoint D.a D.V) (Function.update γ k.succ x)) :
    quadHamiltonStep D.a D.V'
        (γ k.castSucc, quadSegmentMomentum D.a γ k)
      = (γ m.castSucc, quadSegmentMomentum D.a γ m) :=
  quad_stationary_hamiltonStep D.a (ne_of_gt D.a_pos) D.V D.V'
    D.potential_hasDerivAt γ hkm hmin

/-- **Stationarity yields the discrete Euler-Lagrange equation.**  The
composition of the two re-exports above: single-site minimality at a
junction gives the vanishing Euler-Lagrange sum of the supplied
Lagrangian. -/
theorem stationary_eulerLagrange {N : ℕ} (γ : Fin (N + 1) → ℝ)
    {k m : Fin N} (hkm : k.succ = m.castSucc)
    (hmin : ∀ x, localAction (quadTwoPoint D.a D.V) γ
      ≤ localAction (quadTwoPoint D.a D.V) (Function.update γ k.succ x)) :
    quadSegmentMomentum D.a γ k
        + (-(quadSegmentMomentum D.a γ m) - D.V' (γ m.castSucc)) = 0 :=
  (D.eulerLagrange_iff_hamiltonStep γ hkm).mpr
    (D.stationary_hamiltonStep γ hkm hmin)

/-! ## Noether conservation at the supplied stiffness (OL-D3 input) -/

/-- The discrete Hamilton flow of the free class at the supplied
stiffness conserves the momentum readout: the Legendre image of the
translation Noether current.  Re-export of
`free_hamiltonStep_conserves_momentum`. -/
theorem free_flow_conserves_momentum (q : ℝ × ℝ) :
    (quadHamiltonStep D.a (fun _ => 0) q).2 = q.2 :=
  free_hamiltonStep_conserves_momentum D.a q

/-- The discrete Hamilton flow of the free class at the supplied
stiffness conserves the transformed Hamiltonian readout.  Re-export of
`free_hamiltonStep_conserves_hamiltonian`. -/
theorem free_flow_conserves_hamiltonian (q : ℝ × ℝ) :
    quadHamiltonian D.a (fun _ => 0)
        (quadHamiltonStep D.a (fun _ => 0) q).1
        (quadHamiltonStep D.a (fun _ => 0) q).2
      = quadHamiltonian D.a (fun _ => 0) q.1 q.2 :=
  free_hamiltonStep_conserves_hamiltonian D.a q

end MechanicsPremiseData

/-! ## The committed instance of the premise data

The two-state mixing chain of the retained locally hash-pinned run
discharges every field of `MechanicsPremiseData` with exact literals:
the cast stationary law, the cast chain, path length two, stiffness `2`,
and zero potential.  The stiffness value matches the committed free
witness Lagrangian (`witnessL_eq_quadTwoPoint`). -/

/-- The committed premise instance: the two-state mixing chain with the
free witness stiffness. -/
noncomputable def committedMechanicsPremiseData :
    MechanicsPremiseData (Fin 2) where
  pi := InformationProjection.sourcePiR
  P := InformationProjection.sourcePR
  n := 2
  a := 2
  V := fun _ => 0
  V' := fun _ => 0
  pi_pos := fun i => by
    simp only [InformationProjection.sourcePiR]
    exact_mod_cast OPH.Thermodynamics.mixingChainStationary_pos i
  P_pos := InformationProjection.sourcePR_pos
  pi_sum := InformationProjection.sourcePiR_sum
  P_rows := InformationProjection.sourcePR_rows
  a_pos := two_pos
  potential_hasDerivAt := fun t => hasDerivAt_const t 0

/-- The derived-action theorem instantiated at the committed premise
data: the committed chain's two-transition path law is the exponential
tilt of the registered reference by its log-transition action at
multiplier one. -/
theorem committedMechanicsPremiseData_path_law :
    InformationProjection.markovPathLaw
        committedMechanicsPremiseData.pi committedMechanicsPremiseData.P
        committedMechanicsPremiseData.n
      = InformationProjection.tilt
          (InformationProjection.stepUniformRef
            committedMechanicsPremiseData.pi
            committedMechanicsPremiseData.n)
          (InformationProjection.logTransitionAction
            committedMechanicsPremiseData.P
            committedMechanicsPremiseData.n) 1 :=
  committedMechanicsPremiseData.path_law_eq_exponential_tilt

/-- The exact eight-path packet reading of the committed instance: the
tilt of the declared reference by the derived action at multiplier one
is the committed chain path law, with partition constant `1/4`.
Re-export of `InformationProjection.sourceLogTilt_eq_chain` and
`InformationProjection.sourceLogTiltZ_eq`. -/
theorem mechanicsSurface_source_instance :
    InformationProjection.tilt InformationProjection.sourceLogRefR
        InformationProjection.sourceLogActionR 1
      = InformationProjection.sourceTauChainR
    ∧ InformationProjection.tiltZ InformationProjection.sourceLogRefR
        InformationProjection.sourceLogActionR 1 = 1 / 4 :=
  ⟨InformationProjection.sourceLogTilt_eq_chain,
    InformationProjection.sourceLogTiltZ_eq⟩

/-! ## The realized least-action face -/

/-- **Least action on the realized chain (OL-D1 input).**  A chain path
minimizes the committed-form local action of the log-transition
Lagrangian among single-site variations exactly when it is a most
probable path among them, at every path length.  Re-export of
`two_faces_mixing_chain`. -/
theorem mechanicsSurface_least_action_two_faces (M : ℕ)
    (s : Fin (M + 1) → Fin 2) :
    (∀ (i : Fin (M + 1)) (x : Fin 2),
        localAction chainLogLagrangian (chainEmb s)
          ≤ localAction chainLogLagrangian
              (chainEmb (Function.update s i x)))
      ↔ ∀ (i : Fin (M + 1)) (x : Fin 2),
          pathTransitionWeight chainWeight (Function.update s i x)
            ≤ pathTransitionWeight chainWeight s :=
  two_faces_mixing_chain s

/-! ## The committed Noether witnesses -/

/-- The committed translation Noether witness: the two adjacent
free-particle segment momenta of the affine path `0,1,2` agree and are
nonzero.  Re-export of `free_translation_nonzero_noether_witness`. -/
theorem mechanicsSurface_noether_witness :
    segmentMomentum freeDTwo unitTranslationGenerator affineFreePath3
        (0 : Fin 2)
      = segmentMomentum freeDTwo unitTranslationGenerator affineFreePath3
        (1 : Fin 2)
    ∧ segmentMomentum freeDTwo unitTranslationGenerator affineFreePath3
        (0 : Fin 2) ≠ 0 :=
  free_translation_nonzero_noether_witness

/-- The committed Noether correspondence: on the affine free path the
constant chain current equals the Legendre segment momentum on every
segment, the discrete Hamilton flow reproduces the witness chain, and
the conserved momentum and Hamiltonian readouts take the exact values
`2` and `1`.  Re-export of `witness_noether_correspondence`. -/
theorem mechanicsSurface_noether_correspondence :
    (∀ n : Fin 2,
        segmentMomentum freeDTwo unitTranslationGenerator
            affineFreePath3 n
          = quadSegmentMomentum 2 affineFreePath3 n)
      ∧ quadHamiltonStep 2 (fun _ => 0)
          (affineFreePath3 0, quadSegmentMomentum 2 affineFreePath3 0)
          = (affineFreePath3 1, quadSegmentMomentum 2 affineFreePath3 1)
      ∧ quadHamiltonStep 2 (fun _ => 0)
          (affineFreePath3 1, quadSegmentMomentum 2 affineFreePath3 1)
          = (affineFreePath3 2, quadSegmentMomentum 2 affineFreePath3 2)
      ∧ quadSegmentMomentum 2 affineFreePath3 0 = 2
      ∧ quadHamiltonian 2 (fun _ => 0) (affineFreePath3 0)
          (quadSegmentMomentum 2 affineFreePath3 0) = 1
      ∧ quadHamiltonian 2 (fun _ => 0) (affineFreePath3 1)
          (quadSegmentMomentum 2 affineFreePath3 1) = 1 :=
  witness_noether_correspondence

/-! ## Boundary corollaries: the surface carries its own limits -/

/-- **The enrichment is a register row, with proof (row PR-06
boundary).**  Every real curvature enrichment reproduces the realized
binary action at every path length; the canonical bilinear member has
no velocity solver; and the enrichments `1` and `2` differ on both the
Lagrangian and the Hamiltonian face.  The realized history law
therefore does not select the enrichment that the Legendre bridge
consumes, which is why `MechanicsPremiseData` carries the stiffness as
a supplied field.  Re-export of the committed no-go of
`RealizedHistoryLegendreNoGo`. -/
theorem mechanicsSurface_enrichment_no_go :
    (∀ (a : ℝ) (M : ℕ) (s : Fin (M + 1) → Fin 2),
        localAction (chainCurvedLagrangian a) (chainEmb s)
          = logTransitionAction chainWeight s)
      ∧ (¬ ∃ vel : ℝ → ℝ → ℝ, SolvesMomentum chainLogLagrangian vel)
      ∧ chainCurvedLagrangian 1 ≠ chainCurvedLagrangian 2
      ∧ chainCurvedHamiltonian 1 ≠ chainCurvedHamiltonian 2 :=
  ⟨fun a _ s => localAction_chainCurvedLagrangian a s,
    chainLogLagrangian_no_velocity_solver,
    chainCurvedLagrangian_one_ne_two,
    chainCurvedHamiltonian_one_ne_two⟩

/-- **The mode/minimizer boundary (OL-D1 scoping).**  A history can
satisfy the discrete Euler-Lagrange stationarity condition with explicit
derivative witnesses, fail single-site minimality outright, and carry
strictly smaller Gibbs weight than a non-stationary variation at every
positive multiplier.  Stationarity alone does not make a history a mode
or a minimizer.  Re-export of the committed stationary-maximum
counterexample of `StationarySaddleCoverage`. -/
theorem mechanicsSurface_stationary_maximum_scoping :
    (∃ d₂ d₁ : ℝ,
        HasDerivAt (fun x => concaveControlL
            (stationaryMaximumHistory (0 : Fin 2).castSucc) x)
          d₂ (stationaryMaximumHistory (0 : Fin 2).succ)
        ∧ HasDerivAt (fun x => concaveControlL x
            (stationaryMaximumHistory (1 : Fin 2).succ))
          d₁ (stationaryMaximumHistory (1 : Fin 2).castSucc)
        ∧ d₂ + d₁ = 0)
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
  ⟨stationaryMaximumHistory_stationary,
    stationaryMaximumHistory_not_minimal,
    fun _ hlam => gibbs_prefers_nonstationary hlam⟩

/-! ## Composite receipt -/

/-- **The mechanics adequacy surface receipt (issue #731).**  For every
premise bundle: the path law is the exponential tilt of the registered
reference by the log-transition action at multiplier one; single-site
minimizers of the realized log-transition action are exactly the most
probable paths among single-site variations of the realized chain; the finite Legendre transform of
the supplied regular Lagrangian is its quadratic Hamiltonian;
single-site minimality at a junction yields one step of the discrete
Hamilton flow; the free flow at the supplied stiffness conserves the
Noether momentum readout; the canonical bilinear extension of the
realized action has no velocity solver; two realized-history-equivalent
enrichments differ; and the committed stationary history fails
minimality.  The first five conjuncts are the positive surface, the
last three are its boundary. -/
theorem mechanicsAdequacySurface_receipt {Ω : Type*} [Fintype Ω]
    [Nonempty Ω] (D : MechanicsPremiseData Ω) :
    (InformationProjection.markovPathLaw D.pi D.P D.n
        = InformationProjection.tilt
            (InformationProjection.stepUniformRef D.pi D.n)
            (InformationProjection.logTransitionAction D.P D.n) 1)
      ∧ (∀ (M : ℕ) (s : Fin (M + 1) → Fin 2),
          (∀ (i : Fin (M + 1)) (x : Fin 2),
              localAction chainLogLagrangian (chainEmb s)
                ≤ localAction chainLogLagrangian
                    (chainEmb (Function.update s i x)))
            ↔ ∀ (i : Fin (M + 1)) (x : Fin 2),
                pathTransitionWeight chainWeight (Function.update s i x)
                  ≤ pathTransitionWeight chainWeight s)
      ∧ (legendreTransform (quadLagrangian D.a D.V)
            (quadVelocitySolver D.a)
          = quadHamiltonian D.a D.V)
      ∧ (∀ (N : ℕ) (γ : Fin (N + 1) → ℝ) (k m : Fin N),
          k.succ = m.castSucc
            → (∀ x, localAction (quadTwoPoint D.a D.V) γ
                ≤ localAction (quadTwoPoint D.a D.V)
                    (Function.update γ k.succ x))
            → quadHamiltonStep D.a D.V'
                  (γ k.castSucc, quadSegmentMomentum D.a γ k)
                = (γ m.castSucc, quadSegmentMomentum D.a γ m))
      ∧ (∀ q : ℝ × ℝ, (quadHamiltonStep D.a (fun _ => 0) q).2 = q.2)
      ∧ (¬ ∃ vel : ℝ → ℝ → ℝ, SolvesMomentum chainLogLagrangian vel)
      ∧ chainCurvedLagrangian 1 ≠ chainCurvedLagrangian 2
      ∧ (¬ ∀ x, localAction concaveControlL stationaryMaximumHistory
          ≤ localAction concaveControlL
              (Function.update stationaryMaximumHistory
                ((0 : Fin 2).succ) x)) :=
  ⟨D.path_law_eq_exponential_tilt,
    fun _ s => two_faces_mixing_chain s,
    D.legendre_transform_eq,
    fun _ γ _ _ hkm hmin => D.stationary_hamiltonStep γ hkm hmin,
    D.free_flow_conserves_momentum,
    chainLogLagrangian_no_velocity_solver,
    chainCurvedLagrangian_one_ne_two,
    stationaryMaximumHistory_not_minimal⟩

end OPH.Variational

#print axioms OPH.Variational.MechanicsPremiseData.path_law_eq_exponential_tilt
#print axioms OPH.Variational.MechanicsPremiseData.partition_constant
#print axioms OPH.Variational.MechanicsPremiseData.derived_action_gauge_class
#print axioms OPH.Variational.MechanicsPremiseData.bare_multiplier_forced
#print axioms OPH.Variational.MechanicsPremiseData.legendre_transform_eq
#print axioms OPH.Variational.MechanicsPremiseData.legendre_involutive
#print axioms OPH.Variational.MechanicsPremiseData.velocity_solver_solves
#print axioms OPH.Variational.MechanicsPremiseData.lagrangian_strictConvexOn
#print axioms OPH.Variational.MechanicsPremiseData.fenchel_young
#print axioms OPH.Variational.MechanicsPremiseData.eulerLagrange_iff_hamiltonStep
#print axioms OPH.Variational.MechanicsPremiseData.stationary_hamiltonStep
#print axioms OPH.Variational.MechanicsPremiseData.stationary_eulerLagrange
#print axioms OPH.Variational.MechanicsPremiseData.free_flow_conserves_momentum
#print axioms OPH.Variational.MechanicsPremiseData.free_flow_conserves_hamiltonian
#print axioms OPH.Variational.committedMechanicsPremiseData_path_law
#print axioms OPH.Variational.mechanicsSurface_source_instance
#print axioms OPH.Variational.mechanicsSurface_least_action_two_faces
#print axioms OPH.Variational.mechanicsSurface_noether_witness
#print axioms OPH.Variational.mechanicsSurface_noether_correspondence
#print axioms OPH.Variational.mechanicsSurface_enrichment_no_go
#print axioms OPH.Variational.mechanicsSurface_stationary_maximum_scoping
#print axioms OPH.Variational.mechanicsAdequacySurface_receipt
