import Variational.SourceToHamiltonianComposed

namespace OPH.Variational

/-!
# Mode-extremal forcing of the enrichment member on the committed chain
(V3, issue #731)

By the committed no-go of `RealizedHistoryLegendreNoGo`, the realized
history law of the committed mixing chain does not select a member of
the registered curvature family: every curvature reproduces the same
corner action on every realized path.  This module selects a member by a
declared principle instead, and proves the selection unique within the
registered family.

The principle: the realized global mode of the chain path law must be a
real fixed-endpoint extremizer of the enrichment.  The constant state-1
history is the global most probable path at every length, because every
kernel entry is at most the stay weight `W11 = 503/508`
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

Boundary.  The selection principle is a declared rule, and it
references off-alphabet real variations; the realized history law
itself selects nothing, and the committed no-go stays in force
(`mechanicsSurface_enrichment_no_go`).  Uniqueness holds within the
registered one-parameter curvature family (register row PR-06), and
among all real enrichments no uniqueness is claimed.  The extremal
hypothesis and conclusion are fixed-endpoint: endpoint sites of the
constant-1 path carry an affine variation functional and are excluded.
The derivative data consumed here is register row PR-45 through the
committed bundle instance.
-/

/-- The constant state-1 history at every length: the realized global
mode of the committed chain path law. -/
def constOneHistory (M : ℕ) : Fin (M + 1) → Fin 2 := fun _ => 1

/-- The mode-extremal curvature: twice the log of the interior
mode-dominance ratio of the committed chain.  A distinguished member of
the register row PR-06 curvature family, selected by the declared
mode-extremality principle, and provably the unique member compatible
with it (`modeExtremal_forced`). -/
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

/-- **The embedded mode is a real fixed-endpoint minimizer at the
forced member**, at every length, interior junction, and real
replacement value: the variation gap is the nonnegative exact quadratic
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

/-- **The composed receipt at the forced member (issue #731).**  The
committed bundle instance at the mode-extremal curvature, with the
realized global mode as history, discharges the real-minimality
hypothesis of `source_to_hamiltonian_composed` at every length.  The
conclusion holds of the committed kernel: the enriched local action of
the embedded mode equals its derived log-transition action, the
committed path law is the exponential tilt of the registered reference
(row PR-05) by that action at multiplier one, the mode is a most
probable history among interior single-site alphabet variations, and
the embedded mode satisfies the discrete Hamilton equations of the
forced member at every interior junction. -/
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
#print axioms OPH.Variational.modeExtremal_composed_receipt
