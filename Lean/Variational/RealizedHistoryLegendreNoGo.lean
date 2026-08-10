import Variational.LegendreBridge

namespace OPH.Variational

variable {N : ℕ}

/-!
# Realized-history to Legendre non-identifiability (B7, issue #683)

`InformationProjection.LogTransitionAction` produces the finite action of the
committed two-state source chain from its strictly positive transition kernel.
`Variational.LegendreBridge` embeds that source-produced four-corner table in
the reals by the bilinear function `chainLogLagrangian`.  This module audits
whether that realized-history object itself determines a regular
Legendre--Hamilton system.

The answer is a sharp finite no-go.

* The bilinear extension is affine in its second (velocity) slot.  Its
  momentum is independent of velocity, so no global velocity solver exists.
* For every supplied real number `a`, adding
  `(a / 2) * y * (y - 1)` changes no binary corner.  Consequently it changes
  no local action on any realized chain path and reproduces exactly the same
  source log-transition action at every path length.
* For every supplied `a > 0`, the enriched extension is strictly convex in
  the velocity slot, has an exact velocity solver, and has an explicit
  Legendre transform.  The choices `a = 1` and `a = 2` give different
  Lagrangians and different Hamiltonians while agreeing on every realized
  source history.

Thus the finite source dynamics **does produce** the corner action, but it
does not produce the real velocity curvature, regular Legendre map, or
Hamiltonian.  The curvature parameter and off-alphabet extension below are
constructed enrichments.  Selecting one requires additional source data or a
new principle; it cannot be inferred from the current realized-history law.
No continuum, physical clock, units, amplitudes, or laboratory current are
claimed.
-/

/-! ## Affine-fiber form of the constructed corner extension -/

/-- The value on the `y = 0` edge of the bilinear extension constructed from
the source-produced corner table. -/
noncomputable def chainLogBase (x : ℝ) : ℝ :=
  (1 - x) * -Real.log (chainWeight 0 0)
    + x * -Real.log (chainWeight 1 0)

/-- The second-slot slope of the constructed bilinear extension.  It depends
on the first record but not on the second record. -/
noncomputable def chainFiberSlope (x : ℝ) : ℝ :=
  (1 - x) * (-Real.log (chainWeight 0 1)
      + Real.log (chainWeight 0 0))
    + x * (-Real.log (chainWeight 1 1)
      + Real.log (chainWeight 1 0))

/-- The committed bilinear extension is affine in its second slot. -/
theorem chainLogLagrangian_affine_fiber (x y : ℝ) :
    chainLogLagrangian x y = chainLogBase x + chainFiberSlope x * y := by
  unfold chainLogLagrangian chainLogBase chainFiberSlope
  ring

/-- The momentum of the bilinear extension of the source-produced corners is
the same at every velocity. -/
theorem chainLogLagrangian_momentum (x y : ℝ) :
    MomentumRelation chainLogLagrangian x y (chainFiberSlope x) := by
  have h : HasDerivAt
      (fun z : ℝ => chainLogBase x + chainFiberSlope x * z)
      (chainFiberSlope x) y := by
    convert (hasDerivAt_const y (chainLogBase x)).add
      ((hasDerivAt_id y).const_mul (chainFiberSlope x)) using 1
    all_goals ring_nf
  show HasDerivAt (fun z => chainLogLagrangian x z)
    (chainFiberSlope x) y
  rw [show (fun z => chainLogLagrangian x z)
      = (fun z => chainLogBase x + chainFiberSlope x * z) by
        funext z
        exact chainLogLagrangian_affine_fiber x z]
  exact h

/-- Any claimed momentum packet for the bilinear extension is forced to be
the velocity-independent source slope. -/
theorem chainLogLagrangian_momentum_iff (x y p : ℝ) :
    MomentumRelation chainLogLagrangian x y p
      ↔ p = chainFiberSlope x := by
  constructor
  · intro h
    exact HasDerivAt.unique h (chainLogLagrangian_momentum x y)
  · rintro rfl
    exact chainLogLagrangian_momentum x y

/-- **Canonical-extension no-go.**  The bilinear real extension of the
source-produced corner table has no global velocity solver: at a fixed base
record its momentum image is a singleton, not all of `ℝ`. -/
theorem chainLogLagrangian_no_velocity_solver :
    ¬ ∃ vel : ℝ → ℝ → ℝ, SolvesMomentum chainLogLagrangian vel := by
  rintro ⟨vel, hvel⟩
  have h := (chainLogLagrangian_momentum_iff
    0 (vel 0 (chainFiberSlope 0 + 1)) (chainFiberSlope 0 + 1)).mp
      (hvel 0 (chainFiberSlope 0 + 1))
  linarith

/-! ## A corner-invisible family of regular enrichments -/

/-- A constructed real extension of the source corner action.  The supplied
parameter `a` is invisible on the binary alphabet because `y * (y - 1)`
vanishes at `y = 0, 1`. -/
noncomputable def chainCurvedLagrangian (a x y : ℝ) : ℝ :=
  chainLogLagrangian x y + (a / 2) * y * (y - 1)

/-- The zero-curvature member is exactly the canonical bilinear extension. -/
theorem chainCurvedLagrangian_zero :
    chainCurvedLagrangian 0 = chainLogLagrangian := by
  funext x y
  simp [chainCurvedLagrangian]

/-- Polynomial form of the constructed extension. -/
theorem chainCurvedLagrangian_eq (a x y : ℝ) :
    chainCurvedLagrangian a x y
      = chainLogBase x + chainFiberSlope x * y
        + (a / 2) * y * (y - 1) := by
  rw [chainCurvedLagrangian, chainLogLagrangian_affine_fiber]

/-- Every member of the enrichment family reads the source-produced corner
table exactly. -/
theorem chainCurvedLagrangian_corner (a : ℝ) (i j : Fin 2) :
    chainCurvedLagrangian a ((i : ℕ) : ℝ) ((j : ℕ) : ℝ)
      = -Real.log (chainWeight i j) := by
  rw [chainCurvedLagrangian, chainLogLagrangian_corner]
  fin_cases j <;> norm_num

/-- **Same-system realized-history identity.**  At every path length and for
every supplied curvature, the enriched real local action equals the exact
log-transition action of the same committed source chain. -/
theorem localAction_chainCurvedLagrangian (a : ℝ)
    (s : Fin (N + 1) → Fin 2) :
    localAction (chainCurvedLagrangian a) (chainEmb s)
      = logTransitionAction chainWeight s := by
  unfold localAction logTransitionAction chainEmb
  exact Finset.sum_congr rfl fun n _ =>
    chainCurvedLagrangian_corner a (s n.castSucc) (s n.succ)

/-- All constructed enrichments are observationally identical on the full
finite realized-history carrier, not merely on one selected path. -/
theorem chainCurvedLagrangian_realized_indistinguishable (a b : ℝ)
    (s : Fin (N + 1) → Fin 2) :
    localAction (chainCurvedLagrangian a) (chainEmb s)
      = localAction (chainCurvedLagrangian b) (chainEmb s) := by
  rw [localAction_chainCurvedLagrangian,
    localAction_chainCurvedLagrangian]

/-- Every enriched local action gives the same conditional path probability
as the committed source kernel on the full realized carrier. -/
theorem pathTransitionWeight_chainCurvedLagrangian (a : ℝ)
    (s : Fin (N + 1) → Fin 2) :
    pathTransitionWeight chainWeight s
      = Real.exp
          (-localAction (chainCurvedLagrangian a) (chainEmb s)) := by
  rw [localAction_chainCurvedLagrangian,
    pathTransitionWeight_eq_exp_neg chainWeight chainWeight_pos]

/-- The corner-invisible curvature factor vanishes at every embedded source
record. -/
theorem fin2_curvatureFactor_zero (i : Fin 2) :
    (((i : ℕ) : ℝ)) * (((i : ℕ) : ℝ) - 1) = 0 := by
  fin_cases i <;> norm_num

theorem chainEmb_curvatureFactor_zero (s : Fin (N + 1) → Fin 2)
    (n : Fin (N + 1)) :
    chainEmb s n * (chainEmb s n - 1) = 0 := by
  unfold chainEmb
  exact fin2_curvatureFactor_zero (s n)

/-- **Real-variation non-identifiability.**  Although curvatures `a` and `b`
agree on every realized source path, replacing the middle record of any
three-record path by an arbitrary real `z` separates their local actions by
exactly `(a-b) z(z-1)/2`.  The source carrier therefore does not determine
the action seen by the universal real variations used in the Euler--Lagrange
bridge. -/
theorem chainCurvedLagrangian_middleVariation_gap (a b z : ℝ)
    (s : Fin 3 → Fin 2) :
    localAction (chainCurvedLagrangian a)
          (Function.update (chainEmb s) 1 z)
        - localAction (chainCurvedLagrangian b)
          (Function.update (chainEmb s) 1 z)
      = ((a - b) / 2) * z * (z - 1) := by
  rw [localAction_update_one, localAction_update_one]
  simp only [chainCurvedLagrangian]
  have hterminal := chainEmb_curvatureFactor_zero s (2 : Fin 3)
  have ha0 : (a / 2) * chainEmb s 2 * (chainEmb s 2 - 1) = 0 := by
    calc
      _ = (a / 2) * (chainEmb s 2 * (chainEmb s 2 - 1)) := by ring
      _ = 0 := by rw [hterminal, mul_zero]
  have hb0 : (b / 2) * chainEmb s 2 * (chainEmb s 2 - 1) = 0 := by
    calc
      _ = (b / 2) * (chainEmb s 2 * (chainEmb s 2 - 1)) := by ring
      _ = 0 := by rw [hterminal, mul_zero]
  rw [ha0, hb0]
  ring

/-- At the first off-alphabet midpoint, the two regular enrichments `a = 1`
and `a = 2` disagree by the exact amount `1/8` for every realized path. -/
theorem chainCurvedLagrangian_one_two_midpoint_gap (s : Fin 3 → Fin 2) :
    localAction (chainCurvedLagrangian 1)
          (Function.update (chainEmb s) 1 (1 / 2 : ℝ))
        - localAction (chainCurvedLagrangian 2)
          (Function.update (chainEmb s) 1 (1 / 2 : ℝ))
      = 1 / 8 := by
  rw [chainCurvedLagrangian_middleVariation_gap]
  norm_num

/-- The momentum of the enriched extension.  The supplied curvature is
exactly the coefficient that restores velocity dependence. -/
theorem chainCurvedLagrangian_momentum (a x y : ℝ) :
    MomentumRelation (chainCurvedLagrangian a) x y
      (chainFiberSlope x + a * y - a / 2) := by
  have hbase : HasDerivAt
      (fun z : ℝ => chainLogBase x + chainFiberSlope x * z)
      (chainFiberSlope x) y := by
    convert (hasDerivAt_const y (chainLogBase x)).add
      ((hasDerivAt_id y).const_mul (chainFiberSlope x)) using 1
    all_goals ring_nf
  have hbump : HasDerivAt (fun z : ℝ => (a / 2) * z * (z - 1))
      (a * y - a / 2) y := by
    convert (((hasDerivAt_id y).const_mul (a / 2)).mul
      ((hasDerivAt_id y).sub_const 1)) using 1
    all_goals simp only [id_eq]
    all_goals ring
  have h := hbase.add hbump
  show HasDerivAt (fun z => chainCurvedLagrangian a x z)
    (chainFiberSlope x + a * y - a / 2) y
  rw [show (fun z => chainCurvedLagrangian a x z)
      = (fun z => chainLogBase x + chainFiberSlope x * z
          + (a / 2) * z * (z - 1)) by
        funext z
        exact chainCurvedLagrangian_eq a x z]
  convert h using 1
  all_goals ring_nf

/-- The exact velocity solver of a nondegenerate enrichment. -/
noncomputable def chainCurvedVelocitySolver (a : ℝ) : ℝ → ℝ → ℝ :=
  fun x p => (p - chainFiberSlope x + a / 2) / a

/-- Every nonzero supplied curvature gives a global Legendre velocity
solver. -/
theorem chainCurvedVelocitySolver_solves (a : ℝ) (ha : a ≠ 0) :
    SolvesMomentum (chainCurvedLagrangian a)
      (chainCurvedVelocitySolver a) := by
  intro x p
  have h := chainCurvedLagrangian_momentum a x
    (chainCurvedVelocitySolver a x p)
  have hp : chainFiberSlope x
        + a * chainCurvedVelocitySolver a x p - a / 2 = p := by
    unfold chainCurvedVelocitySolver
    field_simp [ha]
    ring
  rwa [hp] at h

/-- Positive supplied curvature makes the enriched velocity fiber strictly
convex on the whole real line. -/
theorem chainCurvedLagrangian_strictConvexOn (a : ℝ) (ha : 0 < a)
    (x : ℝ) :
    StrictConvexOn ℝ Set.univ (fun y => chainCurvedLagrangian a x y) := by
  refine ⟨convex_univ, fun y _ z _ hyz b c hb hc hbc => ?_⟩
  have hc' : c = 1 - b := by linarith
  subst c
  have hsq : 0 < b * (1 - b) * (y - z) ^ 2 :=
    mul_pos (mul_pos hb hc)
      (pow_two_pos_of_ne_zero (sub_ne_zero.mpr hyz))
  have hgap :
      b * chainCurvedLagrangian a x y
          + (1 - b) * chainCurvedLagrangian a x z
          - chainCurvedLagrangian a x (b * y + (1 - b) * z)
        = (a / 2) * (b * (1 - b) * (y - z) ^ 2) := by
    simp only [chainCurvedLagrangian_eq]
    ring
  have hpos : 0 < (a / 2) * (b * (1 - b) * (y - z) ^ 2) := by
    exact mul_pos (half_pos ha) hsq
  simp only [smul_eq_mul]
  nlinarith [hgap, hpos]

/-- The explicit Hamiltonian produced after one supplies a nonzero
curvature enrichment. -/
noncomputable def chainCurvedHamiltonian (a : ℝ) : ℝ → ℝ → ℝ :=
  fun x p =>
    (p - chainFiberSlope x + a / 2) ^ 2 / (2 * a) - chainLogBase x

/-- Exact Legendre transform of the constructed regular enrichment. -/
theorem chainCurved_legendreTransform (a : ℝ) (ha : a ≠ 0) :
    legendreTransform (chainCurvedLagrangian a)
        (chainCurvedVelocitySolver a)
      = chainCurvedHamiltonian a := by
  funext x p
  simp only [legendreTransform, chainCurvedVelocitySolver,
    chainCurvedHamiltonian, chainCurvedLagrangian_eq]
  field_simp [ha]
  ring

/-! ## Constructive non-identifiability witnesses -/

/-- The source histories cannot distinguish the regular enrichments `a = 1`
and `a = 2`, but the real Lagrangians are different off the source alphabet. -/
theorem chainCurvedLagrangian_one_ne_two :
    chainCurvedLagrangian 1 ≠ chainCurvedLagrangian 2 := by
  intro h
  have hval := congrFun (congrFun h (0 : ℝ)) (1 / 2 : ℝ)
  simp only [chainCurvedLagrangian] at hval
  norm_num at hval

/-- The two source-indistinguishable regular enrichments also produce
different Hamiltonians. -/
theorem chainCurvedHamiltonian_one_ne_two :
    chainCurvedHamiltonian 1 ≠ chainCurvedHamiltonian 2 := by
  intro h
  have hval := congrFun (congrFun h (0 : ℝ)) (chainFiberSlope 0)
  simp only [chainCurvedHamiltonian] at hval
  norm_num at hval

/-- **B7 same-system non-identifiability receipt.**  The exact source
log-action is shared by the entire real-extension family; its canonical
bilinear member has no Legendre solver; every positive supplied curvature
gives a strictly convex regular member with an explicit Hamiltonian; and two
such members are provably distinct on both faces.  Therefore the current
realized-history law does not select a real variation/Legendre--Hamilton
bridge. -/
theorem realizedHistory_legendre_nonidentifiability_receipt :
    (∀ (a : ℝ) (M : ℕ) (s : Fin (M + 1) → Fin 2),
        localAction (chainCurvedLagrangian a) (chainEmb s)
          = logTransitionAction chainWeight s)
      ∧ (¬ ∃ vel : ℝ → ℝ → ℝ,
          SolvesMomentum chainLogLagrangian vel)
      ∧ (∀ a : ℝ, 0 < a →
          (∀ x : ℝ,
            StrictConvexOn ℝ Set.univ
              (fun y => chainCurvedLagrangian a x y))
            ∧ SolvesMomentum (chainCurvedLagrangian a)
                (chainCurvedVelocitySolver a)
            ∧ legendreTransform (chainCurvedLagrangian a)
                (chainCurvedVelocitySolver a)
                = chainCurvedHamiltonian a)
      ∧ chainCurvedLagrangian 1 ≠ chainCurvedLagrangian 2
      ∧ chainCurvedHamiltonian 1 ≠ chainCurvedHamiltonian 2
      ∧ (∀ s : Fin 3 → Fin 2,
          localAction (chainCurvedLagrangian 1)
                (Function.update (chainEmb s) 1 (1 / 2 : ℝ))
              - localAction (chainCurvedLagrangian 2)
                (Function.update (chainEmb s) 1 (1 / 2 : ℝ))
            = 1 / 8) := by
  refine ⟨fun a _ s => localAction_chainCurvedLagrangian a s,
    chainLogLagrangian_no_velocity_solver, ?_,
    chainCurvedLagrangian_one_ne_two,
    chainCurvedHamiltonian_one_ne_two,
    chainCurvedLagrangian_one_two_midpoint_gap⟩
  intro a ha
  exact ⟨chainCurvedLagrangian_strictConvexOn a ha,
    chainCurvedVelocitySolver_solves a (ne_of_gt ha),
    chainCurved_legendreTransform a (ne_of_gt ha)⟩

end OPH.Variational

#print axioms OPH.Variational.chainLogLagrangian_no_velocity_solver
#print axioms OPH.Variational.localAction_chainCurvedLagrangian
#print axioms OPH.Variational.chainCurvedLagrangian_middleVariation_gap
#print axioms OPH.Variational.chainCurvedLagrangian_strictConvexOn
#print axioms OPH.Variational.chainCurvedVelocitySolver_solves
#print axioms OPH.Variational.chainCurved_legendreTransform
#print axioms OPH.Variational.realizedHistory_legendre_nonidentifiability_receipt
