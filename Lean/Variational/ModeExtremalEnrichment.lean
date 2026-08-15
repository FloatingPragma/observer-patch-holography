import Variational.SourceToHamiltonianComposed

namespace OPH.Variational

/-!
# Mode-extremal representative in a registered enrichment family
(V3, issue #731)

By the committed no-go of `RealizedHistoryLegendreNoGo`, the realized
history law of the committed mixing chain does not select a member of
the registered curvature family: every curvature reproduces the same
corner action on every realized path.  This module imposes a declared
principle on one deliberately restricted one-parameter ansatz and proves
that the principle fixes its coefficient inside that ansatz.

The principle: the transition-weight global mode must be a real
fixed-endpoint extremizer of the enrichment.  The constant state-1
history has globally maximal conditional transition weight at every
length, because every kernel entry is at most the stay weight `W11 = 503/508`
(`chainWeight_le_stay`, `constOne_global_mode`).  Interior-junction
stationarity of the curvature member `a` at the embedded constant-1
history is a linear equation in `a` with exactly one solution
(`modeExtremal_forced`): the mode-extremal curvature
`a* = 2 * log (W11 ^ 2 / (W10 * W01)) = 2 * log (362055879 / 271780)`.
The value is strictly positive exactly because the interior
mode-dominance inequality `W10 * W01 < W11 ^ 2` is strict
(`modeExtremalCurvature_pos_iff`), so the forced member is strictly
convex and regular.  At the forced member, the embedded mode is a real
fixed-endpoint single-site minimizer at every junction and length, by
an exact quadratic gap identity (`chainCurved_modeExtremal_gap`,
`constOne_realMin`).  The receipt discharges the real-minimality
hypothesis of `source_to_hamiltonian_composed` at the committed bundle
instance, so the composed conclusion holds of the committed kernel at
every length with the forced curvature
(`modeExtremal_composed_receipt`).

Boundary.  The principle is a declared rule, and it references
off-alphabet real variations; the realized history law itself selects
nothing, and the committed no-go stays in force
(`mechanicsSurface_enrichment_no_go`).  The uniqueness theorem applies
only after restricting the enrichment to the one-parameter family that
adds curvature in the second slot and forbids independent first-slot
off-alphabet curvature.  The two-parameter family below proves that this
restriction is load-bearing: a continuum of distinct positive velocity
curvatures has the same binary corners and the same constant-history
stationarity and minimality once a compensating first-slot term is
allowed (`modeExtremal_not_unique_beyond_oneParameter`).  Thus the
declared principle does not select a real enrichment or velocity
curvature beyond the registered ansatz.  The extremal
hypothesis and conclusion are fixed-endpoint: endpoint sites of the
constant-1 path carry an affine variation functional and are excluded.
The derivative data consumed here is register row PR-45 through the
committed bundle instance.
-/

/-- The constant state-1 history at every length: the global maximizer
of the committed chain's conditional transition weight. -/
def constOneHistory (M : ℕ) : Fin (M + 1) → Fin 2 := fun _ => 1

/-- The mode-extremal curvature: twice the log of the interior
mode-dominance ratio of the committed chain.  This is the unique
coefficient compatible with the declared mode-extremality principle
inside the one-parameter family `chainCurvedLagrangian`; it is not a
unique real enrichment (`modeExtremal_not_unique_beyond_oneParameter`). -/
noncomputable def modeExtremalCurvature : ℝ :=
  2 * Real.log (chainWeight 1 1 ^ 2 / (chainWeight 1 0 * chainWeight 0 1))

/-- Every entry of the committed chain is at most the stay weight
`W11 = 503/508`. -/
theorem chainWeight_le_stay (i j : Fin 2) :
    chainWeight i j ≤ chainWeight 1 1 := by
  fin_cases i <;> fin_cases j <;>
    norm_num [chainWeight, OPH.Thermodynamics.mixingChain]

/-- **The constant state-1 history is the global mode at every
length.**  Its transition weight dominates the transition weight of
every path of the same length, over all paths, with no variation
restriction. -/
theorem constOne_global_mode (M : ℕ) (t : Fin (M + 1) → Fin 2) :
    pathTransitionWeight chainWeight t
      ≤ pathTransitionWeight chainWeight (constOneHistory M) := by
  have heq : pathTransitionWeight chainWeight (constOneHistory M)
      = ∏ _n : Fin M, chainWeight 1 1 := rfl
  rw [heq]
  unfold pathTransitionWeight
  exact Finset.prod_le_prod (fun i _ => le_of_lt (chainWeight_pos _ _))
    (fun i _ => chainWeight_le_stay _ _)

/-- The embedded constant-1 history sits at the real value `1` on every
record. -/
theorem chainEmb_constOne (M : ℕ) (n : Fin (M + 1)) :
    chainEmb (constOneHistory M) n = 1 := by
  simp [chainEmb, constOneHistory]

/-- The mode-extremal curvature in log-combination form. -/
theorem modeExtremalCurvature_eq :
    modeExtremalCurvature
      = 2 * (2 * Real.log (chainWeight 1 1)
          - (Real.log (chainWeight 1 0) + Real.log (chainWeight 0 1))) := by
  unfold modeExtremalCurvature
  rw [Real.log_div (pow_ne_zero 2 (ne_of_gt (chainWeight_pos 1 1)))
      (ne_of_gt (mul_pos (chainWeight_pos 1 0) (chainWeight_pos 0 1))),
    Real.log_mul (ne_of_gt (chainWeight_pos 1 0))
      (ne_of_gt (chainWeight_pos 0 1)),
    Real.log_pow]
  push_cast
  ring

/-- The mode-extremal curvature as an exact rational-ratio literal:
`2 * log (362055879 / 271780)`. -/
theorem modeExtremalCurvature_value :
    modeExtremalCurvature = 2 * Real.log (362055879 / 271780) := by
  have h11 : chainWeight 1 1 = 503 / 508 := by
    norm_num [chainWeight, OPH.Thermodynamics.mixingChain]
  have h10 : chainWeight 1 0 = 5 / 508 := by
    norm_num [chainWeight, OPH.Thermodynamics.mixingChain]
  have h01 : chainWeight 0 1 = 107 / 1431 := by
    norm_num [chainWeight, OPH.Thermodynamics.mixingChain]
  unfold modeExtremalCurvature
  rw [h11, h10, h01]
  norm_num

/-- **The forcing equation.**  Interior stationarity of the curvature
member `a` at the record pair `(1, 1)`, in the exact form of the
discrete Euler-Lagrange sum of the committed derivative packet, is a
linear equation in `a` with exactly one solution: the mode-extremal
curvature. -/
theorem modeExtremal_forced (a : ℝ) :
    chainCurvedD2 a 1 1 + chainCurvedD1 1 1 = 0
      ↔ a = modeExtremalCurvature := by
  have he := modeExtremalCurvature_eq
  simp only [chainCurvedD2, chainCurvedD1, chainFiberSlope,
    chainBaseSlope, chainSlopeRate]
  constructor
  · intro h
    linarith [he, h]
  · intro h
    rw [h, he]
    ring

/-- The forcing equation on the embedded chain: at every length and
interior junction of the embedded constant-1 history, the discrete
Euler-Lagrange sum of the curvature member `a` vanishes exactly when
`a` is the mode-extremal curvature. -/
theorem modeExtremal_forced_on_chain (M : ℕ) {k m : Fin M}
    (_hkm : k.succ = m.castSucc) (a : ℝ) :
    chainCurvedD2 a (chainEmb (constOneHistory M) k.castSucc)
          (chainEmb (constOneHistory M) k.succ)
        + chainCurvedD1 (chainEmb (constOneHistory M) m.castSucc)
            (chainEmb (constOneHistory M) m.succ) = 0
      ↔ a = modeExtremalCurvature := by
  simp only [chainEmb_constOne]
  exact modeExtremal_forced a

/-- **The forced curvature is strictly positive**, so the forced member
is a strictly convex regular member of the registered family. -/
theorem modeExtremalCurvature_pos : 0 < modeExtremalCurvature := by
  have h11 : chainWeight 1 1 = 503 / 508 := by
    norm_num [chainWeight, OPH.Thermodynamics.mixingChain]
  have h10 : chainWeight 1 0 = 5 / 508 := by
    norm_num [chainWeight, OPH.Thermodynamics.mixingChain]
  have h01 : chainWeight 0 1 = 107 / 1431 := by
    norm_num [chainWeight, OPH.Thermodynamics.mixingChain]
  have hratio : (1 : ℝ)
      < chainWeight 1 1 ^ 2 / (chainWeight 1 0 * chainWeight 0 1) := by
    rw [h11, h10, h01]
    norm_num
  exact mul_pos two_pos (Real.log_pos hratio)

/-- Positivity of the forced curvature is interior mode dominance: the
curvature is positive exactly when the cross product `W10 * W01` is
strictly below the stay weight squared. -/
theorem modeExtremalCurvature_pos_iff :
    0 < modeExtremalCurvature
      ↔ chainWeight 1 0 * chainWeight 0 1 < chainWeight 1 1 ^ 2 := by
  have hprod : 0 < chainWeight 1 0 * chainWeight 0 1 :=
    mul_pos (chainWeight_pos 1 0) (chainWeight_pos 0 1)
  unfold modeExtremalCurvature
  constructor
  · intro h
    have hlog : 0 < Real.log
        (chainWeight 1 1 ^ 2 / (chainWeight 1 0 * chainWeight 0 1)) := by
      linarith
    have h1 := (Real.log_pos_iff
      (le_of_lt (div_pos (pow_pos (chainWeight_pos 1 1) 2) hprod))).mp hlog
    exact (one_lt_div hprod).mp h1
  · intro h
    have h1 : (1 : ℝ)
        < chainWeight 1 1 ^ 2 / (chainWeight 1 0 * chainWeight 0 1) :=
      (one_lt_div hprod).mpr h
    have := Real.log_pos h1
    linarith

/-- **The exact quadratic gap identity at the forced member.**  The
single-site variation gap of the forced member at the embedded
constant-1 junction is exactly `(a* / 2) * (x - 1) ^ 2`. -/
theorem chainCurved_modeExtremal_gap (x : ℝ) :
    (chainCurvedLagrangian modeExtremalCurvature 1 x
          - chainCurvedLagrangian modeExtremalCurvature 1 1)
        + (chainCurvedLagrangian modeExtremalCurvature x 1
          - chainCurvedLagrangian modeExtremalCurvature 1 1)
      = modeExtremalCurvature / 2 * (x - 1) ^ 2 := by
  have he := modeExtremalCurvature_eq
  simp only [chainCurvedLagrangian, chainLogLagrangian]
  rw [he]
  ring

/-- **The embedded mode is a real fixed-endpoint single-site minimizer at the
forced member**, at every length and interior junction, for every real
single-site replacement value: the variation gap is the nonnegative exact quadratic
of `chainCurved_modeExtremal_gap`. -/
theorem constOne_realMin (M : ℕ) {k m : Fin M}
    (hkm : k.succ = m.castSucc) (x : ℝ) :
    localAction (chainCurvedLagrangian modeExtremalCurvature)
        (chainEmb (constOneHistory M))
      ≤ localAction (chainCurvedLagrangian modeExtremalCurvature)
          (Function.update (chainEmb (constOneHistory M)) k.succ x) := by
  have hdiff := localAction_update_diff
    (chainCurvedLagrangian modeExtremalCurvature)
    (chainEmb (constOneHistory M)) hkm x
  simp only [chainEmb_constOne] at hdiff
  have hgap := chainCurved_modeExtremal_gap x
  have hnn : 0 ≤ modeExtremalCurvature / 2 * (x - 1) ^ 2 :=
    mul_nonneg (le_of_lt (half_pos modeExtremalCurvature_pos))
      (sq_nonneg _)
  linarith [hdiff, hgap, hnn]

/-! ## Boundary: mode extremality does not select the full enrichment

The one-parameter family above permits a corner-invisible quadratic in the
second slot only.  The equally natural two-parameter family below also permits
an independent corner-invisible quadratic in the first slot.  The second-slot
coefficient `a` remains the velocity curvature, while `c` changes the
off-alphabet first-slot dependence.  At the constant-one junction the declared
mode-extremality condition fixes only `a + c = a*`.  Consequently every
`0 < a < a*`, paired with `c = a* - a`, gives a distinct regular velocity
curvature with the same source corners and the same real fixed-endpoint
single-site variation minimum. -/

/-- A two-parameter corner-invisible enrichment.  The parameter `a` is the
second-slot (velocity) curvature and `c` is an independent first-slot
off-alphabet curvature. -/
noncomputable def chainTwoSlotCurvedLagrangian (a c x y : ℝ) : ℝ :=
  chainLogLagrangian x y
    + (a / 2) * y * (y - 1)
    + (c / 2) * x * (x - 1)

/-- The first-slot derivative of the two-parameter family. -/
noncomputable def chainTwoSlotD1 (c x y : ℝ) : ℝ :=
  chainCurvedD1 x y + c * x - c / 2

/-- The second-slot derivative of the two-parameter family.  It is unchanged
by the first-slot parameter `c`, so `a` remains the velocity curvature. -/
noncomputable def chainTwoSlotD2 (a x y : ℝ) : ℝ :=
  chainCurvedD2 a x y

/-- Both new quadratic terms vanish on the binary alphabet, so every member
of the two-parameter family reproduces the exact source corner action. -/
theorem chainTwoSlotCurvedLagrangian_corner (a c : ℝ) (i j : Fin 2) :
    chainTwoSlotCurvedLagrangian a c ((i : ℕ) : ℝ) ((j : ℕ) : ℝ)
      = -Real.log (chainWeight i j) := by
  rw [chainTwoSlotCurvedLagrangian, chainLogLagrangian_corner]
  fin_cases i <;> fin_cases j <;> norm_num

/-- Every two-parameter enrichment has the same local action as the exact
log-transition action on every realized binary history. -/
theorem localAction_chainTwoSlotCurvedLagrangian (a c : ℝ) {M : ℕ}
    (s : Fin (M + 1) → Fin 2) :
    localAction (chainTwoSlotCurvedLagrangian a c) (chainEmb s)
      = logTransitionAction chainWeight s := by
  unfold localAction logTransitionAction chainEmb
  exact Finset.sum_congr rfl fun n _ =>
    chainTwoSlotCurvedLagrangian_corner a c (s n.castSucc) (s n.succ)

/-- The uncurried derivative packet of the two-parameter family.  This makes
the compensating first-slot freedom part of the same exact regular-mechanics
interface as the one-parameter representative. -/
theorem chainTwoSlotCurvedLagrangian_hasFDerivAt (a c : ℝ)
    (p : ℝ × ℝ) :
    HasFDerivAt (Function.uncurry (chainTwoSlotCurvedLagrangian a c))
      (chainTwoSlotD1 c p.1 p.2 • ContinuousLinearMap.fst ℝ ℝ ℝ
        + chainTwoSlotD2 a p.1 p.2 • ContinuousLinearMap.snd ℝ ℝ ℝ)
      p := by
  have hfun : Function.uncurry (chainTwoSlotCurvedLagrangian a c)
      = fun q : ℝ × ℝ =>
          -Real.log (chainWeight 0 0)
            + (chainBaseSlope - c / 2) * q.1
            + (Real.log (chainWeight 0 0) - Real.log (chainWeight 0 1)
                - a / 2) * q.2
            + chainSlopeRate * (q.1 * q.2)
            + (c / 2) * (q.1 * q.1)
            + (a / 2) * (q.2 * q.2) := by
    funext q
    simp only [Function.uncurry, chainTwoSlotCurvedLagrangian,
      chainLogLagrangian, chainBaseSlope, chainSlopeRate]
    ring
  have h := polyTwoPoint_hasFDerivAt (-Real.log (chainWeight 0 0))
    (chainBaseSlope - c / 2)
    (Real.log (chainWeight 0 0) - Real.log (chainWeight 0 1) - a / 2)
    chainSlopeRate (c / 2) (a / 2) p
  have h1 : chainTwoSlotD1 c p.1 p.2
      = chainBaseSlope - c / 2 + chainSlopeRate * p.2
        + 2 * (c / 2) * p.1 := by
    simp only [chainTwoSlotD1, chainCurvedD1]
    ring
  have h2 : chainTwoSlotD2 a p.1 p.2
      = Real.log (chainWeight 0 0) - Real.log (chainWeight 0 1) - a / 2
        + chainSlopeRate * p.1 + 2 * (a / 2) * p.2 := by
    simp only [chainTwoSlotD2, chainCurvedD2, chainFiberSlope,
      chainSlopeRate]
    ring
  rw [hfun, h1, h2]
  exact h

/-- Positive `a` makes the second-slot momentum map of the two-parameter
family strictly monotone, independently of `c`. -/
theorem chainTwoSlotD2_strictMono (a : ℝ) (ha : 0 < a) (x : ℝ) :
    StrictMono fun y => chainTwoSlotD2 a x y := by
  simpa only [chainTwoSlotD2] using chainCurvedD2_strictMono a ha x

/-- The committed velocity solver remains a section of the two-parameter
momentum map at every positive velocity curvature. -/
theorem chainTwoSlotD2_solver_section (a : ℝ) (ha : 0 < a) (x p : ℝ) :
    chainTwoSlotD2 a x (chainCurvedVelocitySolver a x p) = p := by
  simpa only [chainTwoSlotD2] using
    chainCurvedD2_solver_section a ha x p

/-- Distinct velocity-curvature coefficients give distinct second-slot
momentum maps in the two-parameter family. -/
theorem chainTwoSlotD2_curvature_injective :
    Function.Injective (fun a : ℝ => chainTwoSlotD2 a) := by
  intro a b h
  have h00 := congrFun (congrFun h 0) 0
  simp only [chainTwoSlotD2, chainCurvedD2] at h00
  linarith

/-- In the enlarged quadratic family, constant-one stationarity fixes only
the sum of the two off-alphabet curvatures. -/
theorem chainTwoSlot_modeExtremal_forced (a c : ℝ) :
    chainTwoSlotD2 a 1 1 + chainTwoSlotD1 c 1 1 = 0
      ↔ a + c = modeExtremalCurvature := by
  have he := modeExtremalCurvature_eq
  simp only [chainTwoSlotD2, chainTwoSlotD1, chainCurvedD2,
    chainCurvedD1, chainFiberSlope, chainBaseSlope, chainSlopeRate]
  constructor <;> intro h <;> linarith [he]

/-- When `a + c = a*`, the two-parameter family has exactly the same
nonnegative fixed-endpoint variation gap as the one-parameter representative,
even though its velocity curvature may differ. -/
theorem chainTwoSlot_modeExtremal_gap (a c x : ℝ)
    (hsum : a + c = modeExtremalCurvature) :
    (chainTwoSlotCurvedLagrangian a c 1 x
          - chainTwoSlotCurvedLagrangian a c 1 1)
        + (chainTwoSlotCurvedLagrangian a c x 1
          - chainTwoSlotCurvedLagrangian a c 1 1)
      = modeExtremalCurvature / 2 * (x - 1) ^ 2 := by
  have he := modeExtremalCurvature_eq
  have hc : c = modeExtremalCurvature - a := by linarith
  rw [hc]
  simp only [chainTwoSlotCurvedLagrangian, chainLogLagrangian]
  rw [he]
  ring

/-- Every positive split `a + c = a*` makes the embedded constant-one path a
real fixed-endpoint single-site minimizer, not only the split `a=a*`, `c=0`
chosen by the registered one-parameter ansatz. -/
theorem constOne_twoSlot_realMin (a c : ℝ)
    (hsum : a + c = modeExtremalCurvature) (M : ℕ) {k m : Fin M}
    (hkm : k.succ = m.castSucc) (x : ℝ) :
    localAction (chainTwoSlotCurvedLagrangian a c)
        (chainEmb (constOneHistory M))
      ≤ localAction (chainTwoSlotCurvedLagrangian a c)
          (Function.update (chainEmb (constOneHistory M)) k.succ x) := by
  have hdiff := localAction_update_diff
    (chainTwoSlotCurvedLagrangian a c)
    (chainEmb (constOneHistory M)) hkm x
  simp only [chainEmb_constOne] at hdiff
  have hgap := chainTwoSlot_modeExtremal_gap a c x hsum
  have hnn : 0 ≤ modeExtremalCurvature / 2 * (x - 1) ^ 2 :=
    mul_nonneg (le_of_lt (half_pos modeExtremalCurvature_pos))
      (sq_nonneg _)
  linarith [hdiff, hgap, hnn]

/-- **Counterfamily to unrestricted mode-extremal uniqueness.**  Splitting
`a*` equally between the second- and first-slot corner-invisible terms gives a
positive regular velocity curvature `a*/2`, different from `a*`, while keeping
the exact source corners, constant-one stationarity, and real fixed-endpoint
single-site minimality at every length.  The derivative packet, strict monotonicity, and
solver section are included so this is a regular Legendre enrichment, not only
an action-value counterexample. -/
theorem modeExtremal_not_unique_beyond_oneParameter :
    ∃ a c : ℝ,
      0 < a
        ∧ a ≠ modeExtremalCurvature
        ∧ a + c = modeExtremalCurvature
        ∧ chainTwoSlotD2 a ≠ chainTwoSlotD2 modeExtremalCurvature
        ∧ (∀ i j : Fin 2,
            chainTwoSlotCurvedLagrangian a c ((i : ℕ) : ℝ) ((j : ℕ) : ℝ)
              = -Real.log (chainWeight i j))
        ∧ chainTwoSlotD2 a 1 1 + chainTwoSlotD1 c 1 1 = 0
        ∧ (∀ (M : ℕ) (k m : Fin M), k.succ = m.castSucc → ∀ x : ℝ,
            localAction (chainTwoSlotCurvedLagrangian a c)
                (chainEmb (constOneHistory M))
              ≤ localAction (chainTwoSlotCurvedLagrangian a c)
                  (Function.update (chainEmb (constOneHistory M)) k.succ x))
        ∧ (∀ p : ℝ × ℝ,
            HasFDerivAt (Function.uncurry (chainTwoSlotCurvedLagrangian a c))
              (chainTwoSlotD1 c p.1 p.2 • ContinuousLinearMap.fst ℝ ℝ ℝ
                + chainTwoSlotD2 a p.1 p.2 •
                    ContinuousLinearMap.snd ℝ ℝ ℝ) p)
        ∧ (∀ x : ℝ, StrictMono fun y => chainTwoSlotD2 a x y)
        ∧ ∀ x p : ℝ,
            chainTwoSlotD2 a x (chainCurvedVelocitySolver a x p) = p := by
  let a : ℝ := modeExtremalCurvature / 2
  let c : ℝ := modeExtremalCurvature / 2
  have ha : 0 < a := by
    dsimp [a]
    exact half_pos modeExtremalCurvature_pos
  have hne : a ≠ modeExtremalCurvature := by
    dsimp [a]
    linarith [modeExtremalCurvature_pos]
  have hsum : a + c = modeExtremalCurvature := by
    dsimp [a, c]
    ring
  refine ⟨a, c, ha, hne, hsum, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact fun h => hne (chainTwoSlotD2_curvature_injective h)
  · exact chainTwoSlotCurvedLagrangian_corner a c
  · exact (chainTwoSlot_modeExtremal_forced a c).2 hsum
  · exact fun M k m hkm x => constOne_twoSlot_realMin a c hsum M hkm x
  · exact chainTwoSlotCurvedLagrangian_hasFDerivAt a c
  · exact chainTwoSlotD2_strictMono a ha
  · exact chainTwoSlotD2_solver_section a ha

/-- **The composed receipt at the selected one-parameter representative
(issue #731).**  The committed bundle instance at the mode-extremal
coefficient inside `chainCurvedLagrangian`, with the transition-weight
global mode as history, discharges the real-minimality
hypothesis of `source_to_hamiltonian_composed` at every length.  The
conclusion holds of the committed kernel: the enriched local action of
the embedded mode equals its derived log-transition action, the
committed path law is the exponential tilt of the registered reference
(row PR-05) by that action at multiplier one, the mode is a most
probable history among interior single-site alphabet variations, and
the embedded mode satisfies the discrete Hamilton equations of that
representative at every interior junction.  The counterfamily above
shows that this receipt does not select the curvature beyond the chosen
one-parameter ansatz. -/
theorem modeExtremal_composed_receipt (M : ℕ) :
    localAction (chainCurvedLagrangian modeExtremalCurvature)
        (chainEmb (constOneHistory M))
      = InformationProjection.logTransitionAction
          InformationProjection.sourcePR M (constOneHistory M)
    ∧ InformationProjection.markovPathLaw InformationProjection.sourcePiR
          InformationProjection.sourcePR M
        = InformationProjection.tilt
            (InformationProjection.stepUniformRef
              InformationProjection.sourcePiR M)
            (InformationProjection.logTransitionAction
              InformationProjection.sourcePR M) 1
    ∧ (∀ (k m : Fin M), k.succ = m.castSucc → ∀ x : Fin 2,
        InformationProjection.markovPathLaw
            InformationProjection.sourcePiR
            InformationProjection.sourcePR M
            (Function.update (constOneHistory M) k.succ x)
          ≤ InformationProjection.markovPathLaw
              InformationProjection.sourcePiR
              InformationProjection.sourcePR M (constOneHistory M))
    ∧ (∀ (k m : Fin M), k.succ = m.castSucc →
        chainEmb (constOneHistory M) m.castSucc
            = chainCurvedVelocitySolver modeExtremalCurvature
                (chainEmb (constOneHistory M) k.castSucc)
                (chainCurvedD2 modeExtremalCurvature
                  (chainEmb (constOneHistory M) k.castSucc)
                  (chainEmb (constOneHistory M) k.succ))
          ∧ chainCurvedD2 modeExtremalCurvature
                (chainEmb (constOneHistory M) k.castSucc)
                (chainEmb (constOneHistory M) k.succ)
              + chainCurvedD1 (chainEmb (constOneHistory M) m.castSucc)
                  (chainCurvedVelocitySolver modeExtremalCurvature
                    (chainEmb (constOneHistory M) m.castSucc)
                    (chainCurvedD2 modeExtremalCurvature
                      (chainEmb (constOneHistory M) m.castSucc)
                      (chainEmb (constOneHistory M) m.succ))) = 0) := by
  have hmin : ∀ (k m : Fin M), k.succ = m.castSucc → ∀ x : ℝ,
      localAction (chainCurvedLagrangian modeExtremalCurvature)
          (chainEmb (constOneHistory M))
        ≤ localAction (chainCurvedLagrangian modeExtremalCurvature)
            (Function.update (chainEmb (constOneHistory M)) k.succ x) :=
    fun k m hkm x => constOne_realMin M hkm x
  have h := source_to_hamiltonian_composed
    (chainComposedData modeExtremalCurvature modeExtremalCurvature_pos M)
    (constOneHistory M) hmin
  exact ⟨h.1, h.2.1, h.2.2.1, h.2.2.2.1⟩

end OPH.Variational

#print axioms OPH.Variational.chainWeight_le_stay
#print axioms OPH.Variational.constOne_global_mode
#print axioms OPH.Variational.chainEmb_constOne
#print axioms OPH.Variational.modeExtremalCurvature_eq
#print axioms OPH.Variational.modeExtremalCurvature_value
#print axioms OPH.Variational.modeExtremal_forced
#print axioms OPH.Variational.modeExtremal_forced_on_chain
#print axioms OPH.Variational.modeExtremalCurvature_pos
#print axioms OPH.Variational.modeExtremalCurvature_pos_iff
#print axioms OPH.Variational.chainCurved_modeExtremal_gap
#print axioms OPH.Variational.constOne_realMin
#print axioms OPH.Variational.chainTwoSlotCurvedLagrangian_corner
#print axioms OPH.Variational.localAction_chainTwoSlotCurvedLagrangian
#print axioms OPH.Variational.chainTwoSlotCurvedLagrangian_hasFDerivAt
#print axioms OPH.Variational.chainTwoSlotD2_strictMono
#print axioms OPH.Variational.chainTwoSlotD2_solver_section
#print axioms OPH.Variational.chainTwoSlotD2_curvature_injective
#print axioms OPH.Variational.chainTwoSlot_modeExtremal_forced
#print axioms OPH.Variational.chainTwoSlot_modeExtremal_gap
#print axioms OPH.Variational.constOne_twoSlot_realMin
#print axioms OPH.Variational.modeExtremal_not_unique_beyond_oneParameter
#print axioms OPH.Variational.modeExtremal_composed_receipt
