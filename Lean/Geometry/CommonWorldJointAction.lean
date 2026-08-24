import Geometry.CommonWorldMaxwellClockJoin

set_option autoImplicit false

open scoped BigOperators Matrix

/-!
# CW1 joint action on the Maxwell-clock joined record (issue #740; PR-54 stays open)

WHAT IS PROVED.  The joined step record `MaxwellClockJoinedArchitecture` of
`Geometry/CommonWorldMaxwellClockJoin.lean` carries the scaled Maxwell
island and the clock worldline on one record with the dynamics stated
separately: the scaled Ampere update and the Gauss constraint are the
Euler-Lagrange system of the committed window action `windowAction` of
`Screen/ScaledMaxwellStability.lean` through the committed equivalences
`action_stationary_A_iff_ampere` and `action_stationary_phi_iff_gauss`,
while the worldline is defined, not varied.  This module builds the missing
variational principle for the kinematics island and composes one action on
the joined record, at the attainable first rung: a direct sum, not a
coupled action.

(1) The discrete worldline action.  For a declared step window of length
`M` and a path `x : ℕ → Herm2` read at the nodes `0, …, M`, the action
`clockAction M x` is the sum over the `M` increments of the committed
Lorentz quadratic form `lorentzQ` of `Geometry/CanonicalLorentzModule.lean`,
used as committed in its `(+---)` signature convention.  The action depends
on the window values alone (`clockAction_congr`), so the `ℕ`-indexed path
is a finite path on the window.  The expansion `clockAction_expansion` is
exact: the action of a displaced path is the action plus the first
variation plus the action of the displacement, with no truncation, so the
quadratic remainder of the clock sector is the clock action of the
variation itself.  `ClockStationary M x` states that every variation
vanishing at the two declared window endpoints moves the action by exactly
that remainder.  The discrete geodesic characterization
`clockStationary_iff_uniform` proves stationarity equivalent to the
equal-increment clause: the first variation is the sum over interior nodes
`k + 1` of `2 B(Δ_k - Δ_{k+1}, η (k + 1))` for the polarizing bilinear form
`B` (`clockFirstVariation_interior`), vanishing for every variation forces
equal consecutive increments through the nondegeneracy of `B`
(`lorentzB_nondegenerate`, proved in this module; the imported modules
state bilinearity and signature), and equal increments annihilate every
variation by telescoping.  The corollary: the record's clock worldline,
read as the path `joinedPath` of joined events, has constant increment
`stepDuration • frame` (`joinedPath_increment`), is stationary on every
window (`worldline_clockStationary`), and is the unique stationary path
with its committed endpoints (`worldline_is_stationary_path`), so the
kinematics island's worldline is variational, not merely declared.

(2) The joint action.  `jointAction J N A φ x` is the sum of the committed
window action of the Maxwell configuration `(A, φ)` at the record's
declared sources over the window `0, …, N` and the clock action of the
path `x` over the window `0, …, N + 1`, the two windows sharing the
endpoint-fixing structure: both sectors' variations vanish at the nodes
`0` and `N + 1`.  The composed theorem `jointAction_one_principle` states,
on one typed record, both equivalences as conjuncts: stationarity of the
joint action under Maxwell-sector variations at the record's committed
bundle is equivalent to the conjunction of the scaled Ampere residual
clause at the interior steps and the Gauss constraint clause at every
window step (`jointAction_maxwell_sector`, through the committed
equivalences, which are reused and not re-proved), and stationarity under
worldline variations is equivalent to the equal-increment clause
(`jointAction_worldline_sector`), hence, with the committed endpoints, to
the record's clock line (`worldline_sector_solution`).  The committed
configuration of every inhabitant is jointly stationary
(`committed_configuration_stationary`): one variational principle carries
both dynamics.  The receipt `jointAction_receipt` composes the three.

(3) The exact delimitation.  `jointAction_sector_decoupling` proves the
joint action is a direct sum: the Maxwell variation of the joint action is
independent of the worldline argument and the worldline variation is
independent of the Maxwell arguments, as exact partial-variation
independence identities.  The direct-sum action carries both dynamics on
one carrier but couples nothing: a genuine interaction term, a
worldline-field coupling producing a Lorentz-force-shaped clause, is
exactly what the open source-action premise PR-54 still owes.  After this
module the common-world missing joins are exactly the four labels of
`MissingJoin`: the coupling term of the common action, the carrier map,
the matter dynamics, and the observer readout.

(4) Receipts.  The committed joined inhabitant's worldline satisfies the
equal-increment clause with the computed increment `(1, 0)`
(`committedJoinedWitness_increment`).  The comparison path `kinkPath` is
not stationary: at the interior node `1` the explicit variation
`kinkVariation` has first variation `2`, not `0`
(`kink_firstVariation_eq_two`, `kinkPath_not_stationary`), so the
characterization is nonvacuous in both directions
(`stationarity_characterization_nonvacuous`).  Because the form is
Lorentzian, stationarity is not minimality and no minimality is claimed:
the endpoint-fixed spatial variation `spatialVariation` has strictly
negative clock action (`spatialVariation_negative_action`), so the
stationary worldline is strictly beaten by a displaced comparison history
(`stationarity_not_minimality`).  The non-forcing receipt
`jointAction_not_forcing_calibration` exhibits two joint inhabitants with
one certified bundle whose step durations differ and whose joined paths
differ while both carry both stationarity properties: the action selects
the shape of the worldline, not the calibration (PR-15 open).

DECLARED DATA.  The window length, the endpoint-fixing of the variation
classes, the step duration, and every field of the joined record are
declared exactly as recorded in the underlying modules; the joint action
adds no datum and selects none of them.

FALSIFIER.  The module fails if the exact expansion misses a term, if some
stationary path on a positive window has unequal increments or some
equal-increment path fails stationarity, if the record's worldline path
fails the equal-increment clause or fails to be the unique stationary path
with its committed endpoints, if the Maxwell-sector stationarity of the
joint action differs from the committed Ampere and Gauss clauses, if one
of the partial-variation independence identities fails, if the kink path
is stationary, if the spatial variation has nonnegative clock action, or
if the two calibration-distinct inhabitants fail a stationarity property.

WHAT IS NOT PROVED HERE.  No physical units, no coupling, and no continuum
limit: the joint action is dimensionless, its two summands share only the
step index, and no clause relates the Maxwell field at a step to the
worldline event at that step beyond the shared window.  The window and the
endpoints are declared; the step duration stays a declared positive real
with no unit, and the worldline action does not fix it (PR-15 open).  No
port, seam, or face is mapped to a point, interval, or cone of the Lorentz
module (PR-53 open).  The direct sum contains no worldline-field coupling
and produces no Lorentz-force-shaped clause; the source gauge-field,
current, and action attachment is open (PR-54), and issue #740 stays open.
The missing joins after this module are the coupling term of the common
action, the carrier map, the matter dynamics, and the observer readout.

Axiom audit.  Every proof composes committed receipts with exact
mathematics; the module adds no project axiom and uses no native decision
procedure.  The guard lines at the end of the file show at most `propext`,
`Classical.choice`, and `Quot.sound`.
-/

namespace OPH.CommonWorldJointAction

open OPH.C1Lorentz OPH.CommonWorld OPH.CommonWorldIslandBridge
open OPH.CausalComposition
open OPH.TemporalMaxwellEvolution OPH.ScaledMaxwellStability
open OPH.LocalFaceMaxwellAction OPH.DiscreteCoulombGreen
open OPH.CommonWorldInstrumentJoin OPH.CertifiedScaledStepInstrument
open OPH.CommonWorldMaxwellClockJoin

noncomputable section

/-! ## The committed bilinear form: missing algebraic clauses -/

/-- The Lorentz pairing against the zero vector vanishes on the right. -/
theorem lorentzB_zero_right (u : Herm2) : lorentzB u 0 = 0 := by
  have h := lorentzB_smul_right 0 u 0
  rw [zero_smul, zero_mul] at h
  exact h

/-- The Lorentz pairing against the zero vector vanishes on the left. -/
theorem lorentzB_zero_left (v : Herm2) : lorentzB 0 v = 0 := by
  rw [lorentzB_symm]
  exact lorentzB_zero_right v

/-- The Lorentz pairing is odd in its right argument. -/
theorem lorentzB_neg_right (u v : Herm2) : lorentzB u (-v) = -lorentzB u v := by
  have h := lorentzB_smul_right (-1) u v
  rwa [neg_one_smul, neg_one_mul] at h

/-- The Lorentz pairing is subtractive in its right argument. -/
theorem lorentzB_sub_right (u v w : Herm2) :
    lorentzB u (v - w) = lorentzB u v - lorentzB u w := by
  rw [sub_eq_add_neg, lorentzB_add_right, lorentzB_neg_right, ← sub_eq_add_neg]

/-- Exact polarization: the quadratic form of a sum is the two forms plus
twice the pairing. -/
theorem lorentzQ_add (u v : Herm2) :
    lorentzQ (u + v) = lorentzQ u + 2 * lorentzB u v + lorentzQ v := by
  simp only [← lorentzB_self, lorentzB_add_left, lorentzB_add_right]
  rw [lorentzB_symm v u]
  ring

/-- **Nondegeneracy of the committed Lorentz pairing.**  A vector pairing
to zero against every vector is zero: the time axis and the three spatial
coordinate directions separate it. -/
theorem lorentzB_nondegenerate {u : Herm2}
    (h : ∀ w : Herm2, lorentzB u w = 0) : u = 0 := by
  have ht : u.1 = 0 := by
    have h1 : u.1 * 1 - spatialDot u.2 (0 : Spatial) = 0 :=
      h ((1 : ℝ), (0 : Spatial))
    have hdot : spatialDot u.2 (0 : Spatial) = 0 := by
      simp [spatialDot]
    rw [hdot, mul_one, sub_zero] at h1
    exact h1
  have hs : ∀ i : Fin 3, u.2 i = 0 := by
    intro i
    have h2 : u.1 * 0 - spatialDot u.2 (Pi.single i (1 : ℝ)) = 0 :=
      h ((0 : ℝ), Pi.single i (1 : ℝ))
    have hdot : spatialDot u.2 (Pi.single i (1 : ℝ)) = u.2 i := by
      unfold spatialDot
      rw [Finset.sum_eq_single i]
      · simp
      · intro j _ hji
        simp [hji]
      · intro hi
        exact absurd (Finset.mem_univ i) hi
    rw [hdot, mul_zero, zero_sub, neg_eq_zero] at h2
    exact h2
  refine Prod.ext ?_ ?_
  · simpa using ht
  · funext i
    simpa using hs i

/-! ## Deliverable 1: the discrete worldline action -/

/-- The discrete worldline action on the declared window `0, …, M`: the sum
over the `M` step increments of the committed Lorentz quadratic form, in
the committed `(+---)` signature convention. -/
def clockAction (M : ℕ) (x : ℕ → Herm2) : ℝ :=
  ∑ k ∈ Finset.range M, lorentzQ (x (k + 1) - x k)

/-- The exact first variation of the worldline action: the sum over the
window steps of twice the pairing of the path increment against the
variation increment. -/
def clockFirstVariation (M : ℕ) (x η : ℕ → Herm2) : ℝ :=
  ∑ k ∈ Finset.range M, 2 * lorentzB (x (k + 1) - x k) (η (k + 1) - η k)

/-- The worldline action depends on the window values alone: two paths
agreeing at the nodes `0, …, M` have one action.  The `ℕ`-indexed path is
a finite path on the window. -/
theorem clockAction_congr {M : ℕ} {x y : ℕ → Herm2}
    (h : ∀ n, n ≤ M → x n = y n) : clockAction M x = clockAction M y := by
  unfold clockAction
  refine Finset.sum_congr rfl fun k hk => ?_
  have hk' := Finset.mem_range.mp hk
  rw [h k (by omega), h (k + 1) (by omega)]

/-- **Exact expansion of the worldline action.**  The action of a displaced
path is the action, plus the first variation, plus the action of the
displacement; the expansion is exact with no truncation, so the quadratic
remainder is the clock action of the variation itself. -/
theorem clockAction_expansion (M : ℕ) (x η : ℕ → Herm2) :
    clockAction M (x + η) =
      clockAction M x + clockFirstVariation M x η + clockAction M η := by
  unfold clockAction clockFirstVariation
  rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun k _ => ?_
  have hstep : (x + η) (k + 1) - (x + η) k =
      (x (k + 1) - x k) + (η (k + 1) - η k) := by
    simp only [Pi.add_apply]
    abel
  rw [hstep, lorentzQ_add]

/-- Stationarity of the worldline action on the window `0, …, M`: every
variation vanishing at the two declared endpoints moves the action by
exactly the quadratic remainder, the clock action of the variation. -/
def ClockStationary (M : ℕ) (x : ℕ → Herm2) : Prop :=
  ∀ η : ℕ → Herm2, η 0 = 0 → η M = 0 →
    clockAction M (x + η) = clockAction M x + clockAction M η

/-- A stationary path has vanishing first variation at every endpoint-fixed
variation: immediate from the exact expansion. -/
theorem clockStationary_firstVariation_zero {M : ℕ} {x : ℕ → Herm2}
    (hs : ClockStationary M x) {η : ℕ → Herm2} (h0 : η 0 = 0)
    (hM : η M = 0) : clockFirstVariation M x η = 0 := by
  have hexp := clockAction_expansion M x η
  have hst := hs η h0 hM
  linarith

/-- **Interior form of the first variation.**  For a variation vanishing at
the endpoints of the window `0, …, N + 1`, the first variation is the sum
over the interior nodes `k + 1` of `2 B(Δ_k - Δ_{k+1}, η (k + 1))`: the
boundary terms telescope away. -/
theorem clockFirstVariation_interior (N : ℕ) (x η : ℕ → Herm2)
    (h0 : η 0 = 0) (hN : η (N + 1) = 0) :
    clockFirstVariation (N + 1) x η =
      ∑ k ∈ Finset.range N,
        2 * lorentzB ((x (k + 1) - x k) - (x (k + 1 + 1) - x (k + 1)))
          (η (k + 1)) := by
  unfold clockFirstVariation
  have hsplit : ∀ k ∈ Finset.range (N + 1),
      2 * lorentzB (x (k + 1) - x k) (η (k + 1) - η k) =
        2 * lorentzB (x (k + 1) - x k) (η (k + 1)) -
          2 * lorentzB (x (k + 1) - x k) (η k) := by
    intro k _
    rw [lorentzB_sub_right]
    ring
  rw [Finset.sum_congr rfl hsplit, Finset.sum_sub_distrib]
  have hS1 : (∑ k ∈ Finset.range (N + 1),
      2 * lorentzB (x (k + 1) - x k) (η (k + 1))) =
      ∑ k ∈ Finset.range N, 2 * lorentzB (x (k + 1) - x k) (η (k + 1)) := by
    rw [Finset.sum_range_succ, hN, lorentzB_zero_right, mul_zero, add_zero]
  have hS2 : (∑ k ∈ Finset.range (N + 1),
      2 * lorentzB (x (k + 1) - x k) (η k)) =
      ∑ k ∈ Finset.range N,
        2 * lorentzB (x (k + 1 + 1) - x (k + 1)) (η (k + 1)) := by
    rw [Finset.sum_range_succ', h0, lorentzB_zero_right, mul_zero, add_zero]
  rw [hS1, hS2, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [lorentzB_sub_left (x (k + 1) - x k) (x (k + 1 + 1) - x (k + 1))
    (η (k + 1))]
  ring

/-- Equal increments annihilate every endpoint-fixed variation: the first
variation telescopes to the boundary values, which vanish. -/
theorem clockStationary_of_uniform (M : ℕ) (x : ℕ → Herm2)
    (hu : ∀ k, k < M → x (k + 1) - x k = x 1 - x 0) :
    ClockStationary M x := by
  unfold ClockStationary
  intro η h0 hM
  have hfv : clockFirstVariation M x η = 0 := by
    unfold clockFirstVariation
    have hterm : ∀ k ∈ Finset.range M,
        2 * lorentzB (x (k + 1) - x k) (η (k + 1) - η k) =
          2 * lorentzB (x 1 - x 0) (η (k + 1)) -
            2 * lorentzB (x 1 - x 0) (η k) := by
      intro k hk
      rw [hu k (Finset.mem_range.mp hk), lorentzB_sub_right]
      ring
    rw [Finset.sum_congr rfl hterm,
      Finset.sum_range_sub (fun m => 2 * lorentzB (x 1 - x 0) (η m)) M,
      h0, hM, lorentzB_zero_right, mul_zero, sub_zero]
  have hexp := clockAction_expansion M x η
  rw [hexp, hfv, add_zero]

/-- A stationary path has equal consecutive increments at every interior
node: the variation supported at the node `j + 1` with value `w` has first
variation `2 B(Δ_j - Δ_{j+1}, w)`, and nondegeneracy forces the difference
to vanish. -/
theorem clockStationary_consecutive {N : ℕ} {x : ℕ → Herm2}
    (hs : ClockStationary (N + 1) x) {j : ℕ} (hj : j < N) :
    x (j + 1 + 1) - x (j + 1) = x (j + 1) - x j := by
  have key : ∀ w : Herm2,
      lorentzB ((x (j + 1) - x j) - (x (j + 1 + 1) - x (j + 1))) w = 0 := by
    intro w
    let η : ℕ → Herm2 := fun m => if m = j + 1 then w else 0
    have h0 : η 0 = 0 := if_neg (by omega)
    have hN1 : η (N + 1) = 0 := if_neg (by omega)
    have hfv := clockStationary_firstVariation_zero hs h0 hN1
    rw [clockFirstVariation_interior N x η h0 hN1] at hfv
    have hterm : ∀ k ∈ Finset.range N,
        2 * lorentzB ((x (k + 1) - x k) - (x (k + 1 + 1) - x (k + 1)))
            (η (k + 1)) =
          if k = j then
            2 * lorentzB ((x (j + 1) - x j) - (x (j + 1 + 1) - x (j + 1))) w
          else 0 := by
      intro k _
      by_cases hkj : k = j
      · subst hkj
        rw [if_pos rfl]
        have hv : η (k + 1) = w := if_pos rfl
        rw [hv]
      · rw [if_neg hkj]
        have hv : η (k + 1) = 0 := if_neg (by omega)
        rw [hv, lorentzB_zero_right, mul_zero]
    rw [Finset.sum_congr rfl hterm, Finset.sum_ite_eq',
      if_pos (Finset.mem_range.mpr hj)] at hfv
    linarith
  have hD : (x (j + 1) - x j) - (x (j + 1 + 1) - x (j + 1)) = 0 :=
    lorentzB_nondegenerate key
  exact (sub_eq_zero.mp hD).symm

/-- **Discrete geodesic characterization.**  On a positive window, a path
is stationary for the worldline action under all endpoint-fixed variations
exactly when all its window increments are equal. -/
theorem clockStationary_iff_uniform (M : ℕ) (x : ℕ → Herm2) (hM : 0 < M) :
    ClockStationary M x ↔
      ∀ k, k < M → x (k + 1) - x k = x 1 - x 0 := by
  obtain ⟨N, rfl⟩ : ∃ N', M = N' + 1 := ⟨M - 1, by omega⟩
  constructor
  · intro hs k
    induction k with
    | zero => exact fun _ => rfl
    | succ j ih =>
      intro hk
      rw [clockStationary_consecutive hs (show j < N by omega)]
      exact ih (by omega)
  · exact clockStationary_of_uniform (N + 1) x

/-- A path with equal window increments is affine in the step index on the
window. -/
theorem uniform_path_formula {M : ℕ} {x : ℕ → Herm2}
    (hu : ∀ k, k < M → x (k + 1) - x k = x 1 - x 0) :
    ∀ n, n ≤ M → x n = x 0 + (n : ℝ) • (x 1 - x 0) := by
  intro n
  induction n with
  | zero =>
    intro _
    simp
  | succ j ih =>
    intro hj
    have hxj := ih (by omega)
    have hinc := hu j (by omega)
    have hstep : x (j + 1) = x j + (x 1 - x 0) := by
      rw [← hinc]
      abel
    rw [hstep, hxj]
    push_cast
    rw [add_smul, one_smul]
    abel

/-! ## The record's worldline path -/

/-- The worldline path of the joined record: the event component of the
join at each step index. -/
def joinedPath (J : MaxwellClockJoinedArchitecture) (n : ℕ) : Herm2 :=
  (J.join n).2

/-- The joined path in coordinates: the declared step time along the
record's frame direction. -/
theorem joinedPath_eq (J : MaxwellClockJoinedArchitecture) (n : ℕ) :
    joinedPath J n = ((n : ℝ) * J.stepDuration) • (J.frame : Herm2) :=
  J.join_event n

/-- Projection lemma: the joined path is the record's step event at every
index. -/
theorem joinedPath_eq_stepEvent (J : MaxwellClockJoinedArchitecture) (n : ℕ) :
    joinedPath J n = stepEvent J.toInstrumentedCommonWorldArchitecture n :=
  join_stepEvent J n

/-- The joined path starts at the origin. -/
theorem joinedPath_zero (J : MaxwellClockJoinedArchitecture) :
    joinedPath J 0 = 0 := by
  rw [joinedPath_eq]
  simp

/-- The joined path has the constant increment `stepDuration • frame`. -/
theorem joinedPath_increment (J : MaxwellClockJoinedArchitecture) (n : ℕ) :
    joinedPath J (n + 1) - joinedPath J n =
      J.stepDuration • (J.frame : Herm2) := by
  rw [joinedPath_eq, joinedPath_eq, ← sub_smul]
  congr 1
  push_cast
  ring

/-- The equal-increment clause of the joined path: every increment is the
first increment. -/
theorem joinedPath_uniform (J : MaxwellClockJoinedArchitecture) (k : ℕ) :
    joinedPath J (k + 1) - joinedPath J k = joinedPath J 1 - joinedPath J 0 := by
  rw [joinedPath_increment, joinedPath_increment]

/-- The record's worldline path is stationary on every window. -/
theorem worldline_clockStationary (J : MaxwellClockJoinedArchitecture)
    (M : ℕ) : ClockStationary M (joinedPath J) :=
  clockStationary_of_uniform M (joinedPath J) fun k _ => joinedPath_uniform J k

/-- A stationary path with the record's committed endpoints is the record's
worldline path on the whole window. -/
theorem clockStationary_fixed_endpoints (J : MaxwellClockJoinedArchitecture)
    {M : ℕ} (hM : 0 < M) {x : ℕ → Herm2}
    (h0 : x 0 = joinedPath J 0) (hMx : x M = joinedPath J M)
    (hs : ClockStationary M x) :
    ∀ n, n ≤ M → x n = joinedPath J n := by
  have hu := (clockStationary_iff_uniform M x hM).mp hs
  have hform := uniform_path_formula hu
  set c : Herm2 := x 1 - x 0 with hc
  have hx0 : x 0 = 0 := by rw [h0, joinedPath_zero]
  have hMr : (M : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hM.ne'
  have hMform := hform M le_rfl
  rw [hx0, zero_add] at hMform
  have hMc : (M : ℝ) • c =
      ((M : ℝ) * J.stepDuration) • (J.frame : Herm2) := by
    rw [← hMform, hMx, joinedPath_eq]
  have hcval : c = J.stepDuration • (J.frame : Herm2) := by
    have h2 := congrArg (fun z : Herm2 => (M : ℝ)⁻¹ • z) hMc
    simp only [smul_smul] at h2
    rw [inv_mul_cancel₀ hMr, one_smul, ← mul_assoc, inv_mul_cancel₀ hMr,
      one_mul] at h2
    exact h2
  intro n hn
  rw [hform n hn, hx0, hcval, zero_add, smul_smul, joinedPath_eq]

/-- **The worldline is exactly the stationary path with its committed
endpoints.**  On a positive window, a path carries the committed endpoints
and stationarity exactly when it agrees with the record's worldline path at
every window node: the kinematics island's worldline is variational, not
merely declared. -/
theorem worldline_is_stationary_path (J : MaxwellClockJoinedArchitecture)
    {M : ℕ} (hM : 0 < M) (x : ℕ → Herm2) :
    (x 0 = joinedPath J 0 ∧ x M = joinedPath J M ∧ ClockStationary M x) ↔
      ∀ n, n ≤ M → x n = joinedPath J n := by
  constructor
  · rintro ⟨h0, hMx, hs⟩
    exact clockStationary_fixed_endpoints J hM h0 hMx hs
  · intro hagree
    refine ⟨hagree 0 (Nat.zero_le M), hagree M le_rfl,
      clockStationary_of_uniform M x ?_⟩
    intro k hk
    rw [hagree (k + 1) (by omega), hagree k (by omega), hagree 1 (by omega),
      hagree 0 (by omega)]
    exact joinedPath_uniform J k

/-! ## Deliverable 2: the joint action on the joined record -/

/-- The joint action of the joined record on the declared window: the
committed window action of the Maxwell configuration at the record's
declared sources over the steps `0, …, N`, plus the worldline action of the
path over the nodes `0, …, N + 1`.  The two summands share the
endpoint-fixing structure: both sectors' variations vanish at the nodes `0`
and `N + 1`. -/
def jointAction (J : MaxwellClockJoinedArchitecture) (N : ℕ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (x : ℕ → Herm2) : ℝ :=
  windowAction J.scaled.h N A φ J.scaled.rho J.scaled.J +
    clockAction (N + 1) x

/-- A Maxwell-sector displacement of the joint action moves it exactly as
it moves the window action: the clock summand cancels. -/
theorem jointAction_maxwell_shift (J : MaxwellClockJoinedArchitecture)
    (N : ℕ) (A A' : ℕ → Fin 30 → ℝ) (φ φ' : ℕ → Fin 12 → ℝ)
    (x : ℕ → Herm2) (r : ℝ) :
    jointAction J N A' φ' x = jointAction J N A φ x + r ↔
      windowAction J.scaled.h N A' φ' J.scaled.rho J.scaled.J =
        windowAction J.scaled.h N A φ J.scaled.rho J.scaled.J + r := by
  unfold jointAction
  constructor <;> intro h <;> linarith

/-- A worldline displacement of the joint action moves it exactly as it
moves the clock action: the Maxwell summand cancels. -/
theorem jointAction_worldline_shift (J : MaxwellClockJoinedArchitecture)
    (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (x x' : ℕ → Herm2)
    (r : ℝ) :
    jointAction J N A φ x' = jointAction J N A φ x + r ↔
      clockAction (N + 1) x' = clockAction (N + 1) x + r := by
  unfold jointAction
  constructor <;> intro h <;> linarith

/-- Stationarity of the joint action under Maxwell-sector variations at the
record's committed bundle: seam-potential variations vanishing at the
window endpoints and unrestricted port-potential variations move the joint
action by exactly the committed quadratic remainder. -/
def MaxwellSectorStationary (J : MaxwellClockJoinedArchitecture) (N : ℕ)
    (x : ℕ → Herm2) : Prop :=
  (∀ a : ℕ → Fin 30 → ℝ, a 0 = 0 → a (N + 1) = 0 →
    jointAction J N (J.scaled.A + a) J.scaled.phi x =
      jointAction J N J.scaled.A J.scaled.phi x +
        quadraticRemainder J.scaled.h N a 0) ∧
  (∀ f : ℕ → Fin 12 → ℝ,
    jointAction J N J.scaled.A (J.scaled.phi + f) x =
      jointAction J N J.scaled.A J.scaled.phi x +
        quadraticRemainder J.scaled.h N 0 f)

/-- Stationarity of the joint action under worldline variations: variations
of the path vanishing at the nodes `0` and `N + 1` move the joint action by
exactly the clock action of the variation. -/
def WorldlineSectorStationary (J : MaxwellClockJoinedArchitecture) (N : ℕ)
    (x : ℕ → Herm2) : Prop :=
  ∀ η : ℕ → Herm2, η 0 = 0 → η (N + 1) = 0 →
    jointAction J N J.scaled.A J.scaled.phi (x + η) =
      jointAction J N J.scaled.A J.scaled.phi x + clockAction (N + 1) η

/-- The Maxwell-sector stationarity of the joint action is the committed
window-action stationarity, clause by clause. -/
theorem maxwellSector_iff (J : MaxwellClockJoinedArchitecture) (N : ℕ)
    (x : ℕ → Herm2) :
    MaxwellSectorStationary J N x ↔
      ((∀ a : ℕ → Fin 30 → ℝ, a 0 = 0 → a (N + 1) = 0 →
        windowAction J.scaled.h N (J.scaled.A + a) J.scaled.phi J.scaled.rho
            J.scaled.J =
          windowAction J.scaled.h N J.scaled.A J.scaled.phi J.scaled.rho
              J.scaled.J +
            quadraticRemainder J.scaled.h N a 0) ∧
      (∀ f : ℕ → Fin 12 → ℝ,
        windowAction J.scaled.h N J.scaled.A (J.scaled.phi + f) J.scaled.rho
            J.scaled.J =
          windowAction J.scaled.h N J.scaled.A J.scaled.phi J.scaled.rho
              J.scaled.J +
            quadraticRemainder J.scaled.h N 0 f)) := by
  unfold MaxwellSectorStationary
  refine and_congr ?_ ?_
  · exact forall_congr' fun a => imp_congr Iff.rfl (imp_congr Iff.rfl
      (jointAction_maxwell_shift J N J.scaled.A (J.scaled.A + a) J.scaled.phi
        J.scaled.phi x (quadraticRemainder J.scaled.h N a 0)))
  · exact forall_congr' fun f =>
      jointAction_maxwell_shift J N J.scaled.A J.scaled.A J.scaled.phi
        (J.scaled.phi + f) x (quadraticRemainder J.scaled.h N 0 f)

/-- **Maxwell sector of the one principle.**  Stationarity of the joint
action under Maxwell-sector variations is equivalent to the conjunction of
the scaled Ampere residual clause at the interior steps and the Gauss
constraint clause at every window step, through the committed equivalences
of `Screen/ScaledMaxwellStability.lean`. -/
theorem jointAction_maxwell_sector (J : MaxwellClockJoinedArchitecture)
    (N : ℕ) (x : ℕ → Herm2) :
    MaxwellSectorStationary J N x ↔
      ((∀ m, m < N →
          ampereResidual J.scaled.h J.scaled.A J.scaled.phi J.scaled.J m = 0) ∧
        (∀ n, n < N + 1 →
          realBoundary (electricFieldScaled J.scaled.h J.scaled.A
              J.scaled.phi n) =
            J.scaled.rho n)) := by
  rw [maxwellSector_iff]
  exact and_congr
    (action_stationary_A_iff_ampere J.scaled.h J.scaled.h_pos.ne' N
      J.scaled.A J.scaled.phi J.scaled.rho J.scaled.J)
    (action_stationary_phi_iff_gauss J.scaled.h J.scaled.h_pos.ne' N
      J.scaled.A J.scaled.phi J.scaled.rho J.scaled.J)

/-- The worldline-sector stationarity of the joint action is the clock
stationarity of the path on the window `0, …, N + 1`. -/
theorem worldlineSector_iff_clockStationary
    (J : MaxwellClockJoinedArchitecture) (N : ℕ) (x : ℕ → Herm2) :
    WorldlineSectorStationary J N x ↔ ClockStationary (N + 1) x := by
  unfold WorldlineSectorStationary ClockStationary
  refine forall_congr' fun η => imp_congr Iff.rfl (imp_congr Iff.rfl ?_)
  exact jointAction_worldline_shift J N J.scaled.A J.scaled.phi x (x + η)
    (clockAction (N + 1) η)

/-- **Worldline sector of the one principle.**  Stationarity of the joint
action under worldline variations is equivalent to the equal-increment
clause of the path on the window. -/
theorem jointAction_worldline_sector (J : MaxwellClockJoinedArchitecture)
    (N : ℕ) (x : ℕ → Herm2) :
    WorldlineSectorStationary J N x ↔
      ∀ k, k < N + 1 → x (k + 1) - x k = x 1 - x 0 :=
  (worldlineSector_iff_clockStationary J N x).trans
    (clockStationary_iff_uniform (N + 1) x (Nat.succ_pos N))

/-- With the committed endpoints, worldline-sector stationarity of the
joint action is equivalent to being the record's clock line on the whole
window. -/
theorem worldline_sector_solution (J : MaxwellClockJoinedArchitecture)
    (N : ℕ) (x : ℕ → Herm2) (h0 : x 0 = joinedPath J 0)
    (hend : x (N + 1) = joinedPath J (N + 1)) :
    WorldlineSectorStationary J N x ↔
      ∀ n, n ≤ N + 1 → x n = joinedPath J n := by
  rw [worldlineSector_iff_clockStationary J N x]
  constructor
  · intro hs
    exact clockStationary_fixed_endpoints J (Nat.succ_pos N) h0 hend hs
  · intro hagree
    exact ((worldline_is_stationary_path J (Nat.succ_pos N) x).mpr hagree).2.2

/-- **The committed configuration is jointly stationary.**  On every
inhabitant, the record's certified bundle and the record's worldline path
are a stationary point of the one joint action in both sectors: the Ampere
residuals vanish by the bundle's committed law, the Gauss clause holds by
the committed propagation, and the worldline has equal increments. -/
theorem committed_configuration_stationary
    (J : MaxwellClockJoinedArchitecture) (N : ℕ) :
    MaxwellSectorStationary J N (joinedPath J) ∧
      WorldlineSectorStationary J N (joinedPath J) := by
  constructor
  · refine (jointAction_maxwell_sector J N (joinedPath J)).mpr ⟨?_, ?_⟩
    · intro m _
      exact (ampereEvolutionScaled_iff_residual J.scaled.h J.scaled.A
        J.scaled.phi J.scaled.J).mp J.scaled.ampere m
    · intro n _
      exact gauss_propagation_scaled J.scaled.h J.scaled.A J.scaled.phi
        J.scaled.J J.scaled.rho J.scaled.ampere J.scaled.gauss_init
        J.scaled.continuity n
  · refine (jointAction_worldline_sector J N (joinedPath J)).mpr ?_
    intro k _
    exact joinedPath_uniform J k

/-- **One variational principle, both dynamics (issue #740; PR-54 stays
open).**  On one typed record, the two sector equivalences as conjuncts:
Maxwell-sector stationarity of the joint action is the scaled Ampere and
Gauss system, and worldline-sector stationarity is the equal-increment
clause. -/
theorem jointAction_one_principle (J : MaxwellClockJoinedArchitecture)
    (N : ℕ) (x : ℕ → Herm2) :
    (MaxwellSectorStationary J N x ↔
      ((∀ m, m < N →
          ampereResidual J.scaled.h J.scaled.A J.scaled.phi J.scaled.J m = 0) ∧
        (∀ n, n < N + 1 →
          realBoundary (electricFieldScaled J.scaled.h J.scaled.A
              J.scaled.phi n) =
            J.scaled.rho n))) ∧
    (WorldlineSectorStationary J N x ↔
      ∀ k, k < N + 1 → x (k + 1) - x k = x 1 - x 0) :=
  ⟨jointAction_maxwell_sector J N x, jointAction_worldline_sector J N x⟩

/-- **The joint-action receipt.**  Every joined record carries: the joint
stationarity of its committed configuration; the two sector equivalences at
every path; and, at the committed endpoints, the identification of
worldline-sector stationarity with the record's clock line. -/
theorem jointAction_receipt (J : MaxwellClockJoinedArchitecture) (N : ℕ) :
    (MaxwellSectorStationary J N (joinedPath J) ∧
      WorldlineSectorStationary J N (joinedPath J)) ∧
    (∀ x : ℕ → Herm2,
      (MaxwellSectorStationary J N x ↔
        ((∀ m, m < N →
            ampereResidual J.scaled.h J.scaled.A J.scaled.phi J.scaled.J
              m = 0) ∧
          (∀ n, n < N + 1 →
            realBoundary (electricFieldScaled J.scaled.h J.scaled.A
                J.scaled.phi n) =
              J.scaled.rho n))) ∧
      (WorldlineSectorStationary J N x ↔
        ∀ k, k < N + 1 → x (k + 1) - x k = x 1 - x 0)) ∧
    (∀ x : ℕ → Herm2, x 0 = joinedPath J 0 →
      x (N + 1) = joinedPath J (N + 1) →
      (WorldlineSectorStationary J N x ↔
        ∀ n, n ≤ N + 1 → x n = joinedPath J n)) :=
  ⟨committed_configuration_stationary J N,
    fun x => jointAction_one_principle J N x,
    fun x h0 hend => worldline_sector_solution J N x h0 hend⟩

/-! ## Deliverable 3: the exact delimitation -/

/-- **Sector decoupling.**  The joint action is a direct sum: it splits as
window action plus clock action, a Maxwell displacement moves it
independently of the worldline argument in both slots, and a worldline
displacement moves it independently of both Maxwell arguments.  The
direct-sum action carries both dynamics on one carrier but couples
nothing. -/
theorem jointAction_sector_decoupling (J : MaxwellClockJoinedArchitecture)
    (N : ℕ) :
    (∀ (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (x : ℕ → Herm2),
      jointAction J N A φ x =
        windowAction J.scaled.h N A φ J.scaled.rho J.scaled.J +
          clockAction (N + 1) x) ∧
    (∀ (A a : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (x y : ℕ → Herm2),
      jointAction J N (A + a) φ x - jointAction J N A φ x =
        jointAction J N (A + a) φ y - jointAction J N A φ y) ∧
    (∀ (A : ℕ → Fin 30 → ℝ) (φ f : ℕ → Fin 12 → ℝ) (x y : ℕ → Herm2),
      jointAction J N A (φ + f) x - jointAction J N A φ x =
        jointAction J N A (φ + f) y - jointAction J N A φ y) ∧
    (∀ (A A' : ℕ → Fin 30 → ℝ) (φ φ' : ℕ → Fin 12 → ℝ) (x η : ℕ → Herm2),
      jointAction J N A φ (x + η) - jointAction J N A φ x =
        jointAction J N A' φ' (x + η) - jointAction J N A' φ' x) := by
  refine ⟨fun _ _ _ => rfl, ?_, ?_, ?_⟩
  · intro A a φ x y
    unfold jointAction
    ring
  · intro A φ f x y
    unfold jointAction
    ring
  · intro A A' φ φ' x η
    unfold jointAction
    ring

/-- The common-world joins that stay missing after this module, as named
labels: the coupling term of the common action, that is a worldline-field
interaction producing a Lorentz-force-shaped clause (PR-54 open); the
carrier map from ports, seams, and faces to points, intervals, and cones of
the Lorentz module (PR-53 open); the matter dynamics beyond the committed
structure premises; and the observer readout. -/
inductive MissingJoin : Type
  /-- The interaction term of the common action: a worldline-field coupling
  producing a Lorentz-force-shaped clause. -/
  | couplingTerm
  /-- Ports, seams, and faces sent to points, intervals, or cones of the
  declared module. -/
  | carrierMap
  /-- Matter dynamics beyond the committed structure premises. -/
  | matterDynamics
  /-- Public outcome, readback, and provenance of the instrument. -/
  | observerReadout

/-- The explicit inhabitant of the missing-join label type: the coupling
term heads the list. -/
def firstMissingJoin : MissingJoin := MissingJoin.couplingTerm

/-! ## Deliverable 4: inhabitants and receipts -/

/-- The committed joined inhabitant's worldline increment, computed: at
step duration `1` along the standard frame, every increment is `(1, 0)`. -/
theorem committedJoinedWitness_increment (n : ℕ) :
    joinedPath committedJoinedWitness (n + 1) -
        joinedPath committedJoinedWitness n =
      ((1 : ℝ), (0 : Spatial)) := by
  have hδ : committedJoinedWitness.stepDuration = 1 := rfl
  have hf : (committedJoinedWitness.frame : Herm2) =
      ((1 : ℝ), (0 : Spatial)) := rfl
  rw [joinedPath_increment, hδ, hf, one_smul]

/-- The committed joined inhabitant's worldline satisfies the
equal-increment clause. -/
theorem committedJoinedWitness_equal_increments (k : ℕ) :
    joinedPath committedJoinedWitness (k + 1) -
        joinedPath committedJoinedWitness k =
      joinedPath committedJoinedWitness 1 -
        joinedPath committedJoinedWitness 0 :=
  joinedPath_uniform committedJoinedWitness k

/-- The comparison path: at the origin at node `0` and at the fixed event
`(1, 0)` from node `1` onward, so the first increment is `(1, 0)` and every
later increment is `0`. -/
def kinkPath : ℕ → Herm2
  | 0 => 0
  | _ + 1 => ((1 : ℝ), (0 : Spatial))

theorem kinkPath_zero : kinkPath 0 = 0 := rfl

theorem kinkPath_succ (n : ℕ) :
    kinkPath (n + 1) = ((1 : ℝ), (0 : Spatial)) := rfl

/-- The two increments of the comparison path on the window `0, …, 2`. -/
theorem kinkPath_increments :
    kinkPath 1 - kinkPath 0 = ((1 : ℝ), (0 : Spatial)) ∧
      kinkPath 2 - kinkPath 1 = 0 := by
  constructor
  · rw [kinkPath_zero,
      show kinkPath 1 = ((1 : ℝ), (0 : Spatial)) from kinkPath_succ 0,
      sub_zero]
  · rw [show kinkPath 2 = ((1 : ℝ), (0 : Spatial)) from kinkPath_succ 1,
      show kinkPath 1 = ((1 : ℝ), (0 : Spatial)) from kinkPath_succ 0,
      sub_self]

/-- The comparison path has unequal increments. -/
theorem kinkPath_nonuniform :
    kinkPath 2 - kinkPath 1 ≠ kinkPath 1 - kinkPath 0 := by
  rw [kinkPath_increments.2, kinkPath_increments.1]
  intro h
  have h1 : (0 : ℝ) = 1 := congrArg Prod.fst h
  norm_num at h1

/-- The explicit variation supported at the interior node `1` of the window
`0, …, 2`, with the timelike value `(1, 0)`. -/
def kinkVariation : ℕ → Herm2 :=
  fun n => if n = 1 then ((1 : ℝ), (0 : Spatial)) else 0

/-- The first variation of the comparison path at the explicit interior
variation is `2`, not `0`: the interior node witnesses the failure of
stationarity. -/
theorem kink_firstVariation_eq_two :
    clockFirstVariation 2 kinkPath kinkVariation = 2 := by
  unfold clockFirstVariation
  rw [Finset.sum_range_succ, Finset.sum_range_succ, Finset.sum_range_zero]
  have hv0 : kinkVariation 0 = 0 := by simp [kinkVariation]
  have hv1 : kinkVariation 1 = ((1 : ℝ), (0 : Spatial)) := by
    simp [kinkVariation]
  have hv2 : kinkVariation 2 = 0 := by simp [kinkVariation]
  have hp1 : kinkPath 1 - kinkPath 0 = ((1 : ℝ), (0 : Spatial)) :=
    kinkPath_increments.1
  have hp2 : kinkPath 2 - kinkPath 1 = 0 := kinkPath_increments.2
  rw [zero_add, hv0, hv1, hv2, hp1, hp2, sub_zero, zero_sub,
    lorentzB_zero_left, mul_zero, add_zero]
  have hB : lorentzB (((1 : ℝ), (0 : Spatial)) : Herm2)
      ((1 : ℝ), (0 : Spatial)) = 1 := by
    simp [lorentzB, spatialDot]
  rw [hB]
  norm_num

/-- The comparison path is not stationary: the explicit interior variation
has nonzero first variation. -/
theorem kinkPath_not_stationary : ¬ ClockStationary 2 kinkPath := by
  intro hs
  have h0 : kinkVariation 0 = 0 := by simp [kinkVariation]
  have h2 : kinkVariation 2 = 0 := by simp [kinkVariation]
  have hfv := clockStationary_firstVariation_zero hs h0 h2
  rw [kink_firstVariation_eq_two] at hfv
  norm_num at hfv

/-- The characterization is nonvacuous in both directions: the committed
inhabitant's worldline is stationary on every window, and the comparison
path is not stationary. -/
theorem stationarity_characterization_nonvacuous :
    (∀ M : ℕ, ClockStationary M (joinedPath committedJoinedWitness)) ∧
      ¬ ClockStationary 2 kinkPath :=
  ⟨fun M => worldline_clockStationary committedJoinedWitness M,
    kinkPath_not_stationary⟩

/-- The endpoint-fixed spatial variation at the interior node `1` of the
window `0, …, 2`. -/
def spatialVariation : ℕ → Herm2 :=
  fun n => if n = 1 then ((0 : ℝ), Pi.single 0 (1 : ℝ)) else 0

/-- The spatial variation vanishes at the window endpoints and has clock
action `-2`: the quadratic remainder of the Lorentzian form is
indefinite. -/
theorem spatialVariation_negative_action :
    spatialVariation 0 = 0 ∧ spatialVariation 2 = 0 ∧
      clockAction 2 spatialVariation = -2 := by
  refine ⟨by simp [spatialVariation], by simp [spatialVariation], ?_⟩
  unfold clockAction
  rw [Finset.sum_range_succ, Finset.sum_range_succ, Finset.sum_range_zero]
  have hv0 : spatialVariation 0 = 0 := by simp [spatialVariation]
  have hv1 : spatialVariation 1 = ((0 : ℝ), Pi.single 0 (1 : ℝ)) := by
    simp [spatialVariation]
  have hv2 : spatialVariation 2 = 0 := by simp [spatialVariation]
  rw [zero_add, hv0, hv1, hv2, sub_zero, zero_sub]
  have hQ : lorentzQ (((0 : ℝ), Pi.single 0 (1 : ℝ)) : Herm2) = -1 := by
    simp [lorentzQ, spatialNormSq, Pi.single_apply]
  have hQneg : lorentzQ (-(((0 : ℝ), Pi.single 0 (1 : ℝ)) : Herm2)) = -1 := by
    rw [show (-(((0 : ℝ), Pi.single 0 (1 : ℝ)) : Herm2)) =
        ((-1 : ℝ) • (((0 : ℝ), Pi.single 0 (1 : ℝ)) : Herm2)) from
      (neg_one_smul ℝ _).symm, lorentzQ_smul, hQ]
    norm_num
  rw [hQ, hQneg]
  norm_num

/-- **Stationarity is not minimality.**  On the Lorentzian form the
stationary worldline is strictly beaten by the spatially displaced
comparison history: the displaced clock action is strictly smaller.  No
minimality is claimed anywhere in this module. -/
theorem stationarity_not_minimality (J : MaxwellClockJoinedArchitecture) :
    clockAction 2 (joinedPath J + spatialVariation) <
      clockAction 2 (joinedPath J) := by
  have hexp := clockAction_expansion 2 (joinedPath J) spatialVariation
  have hfv := clockStationary_firstVariation_zero
    (worldline_clockStationary J 2) spatialVariation_negative_action.1
    spatialVariation_negative_action.2.1
  have hneg := spatialVariation_negative_action.2.2
  linarith

/-- **The joint action does not fix the clock calibration.**  Two joint
inhabitants with one certified bundle differ in the step duration and in
the joined path while both carry both stationarity properties of the joint
action: the action selects the shape of the worldline, not the calibration
(PR-15 open). -/
theorem jointAction_not_forcing_calibration :
    ∃ J₁ J₂ : MaxwellClockJoinedArchitecture,
      J₁.scaled = J₂.scaled ∧
        J₁.stepDuration ≠ J₂.stepDuration ∧
        joinedPath J₁ 1 ≠ joinedPath J₂ 1 ∧
        (∀ N : ℕ, MaxwellSectorStationary J₁ N (joinedPath J₁) ∧
          WorldlineSectorStationary J₁ N (joinedPath J₁)) ∧
        (∀ N : ℕ, MaxwellSectorStationary J₂ N (joinedPath J₂) ∧
          WorldlineSectorStationary J₂ N (joinedPath J₂)) := by
  refine ⟨doubleStepJoinedWitness, committedJoinedWitness, rfl, ?_, ?_,
    fun N => committed_configuration_stationary doubleStepJoinedWitness N,
    fun N => committed_configuration_stationary committedJoinedWitness N⟩
  · show (2 : ℝ) ≠ 1
    norm_num
  · exact joined_event_ne_of_stepDuration_ne committedJoinedWitness
      doubleStepJoinedWitness rfl (by show (2 : ℝ) ≠ 1; norm_num)

end

end OPH.CommonWorldJointAction

/- Axiom audit: committed receipts and exact mathematics only.  Expected
axioms per line: at most `propext`, `Classical.choice`, `Quot.sound`.  No
native decision procedure is used. -/

#print axioms OPH.CommonWorldJointAction.lorentzB_zero_right
#print axioms OPH.CommonWorldJointAction.lorentzB_zero_left
#print axioms OPH.CommonWorldJointAction.lorentzB_neg_right
#print axioms OPH.CommonWorldJointAction.lorentzB_sub_right
#print axioms OPH.CommonWorldJointAction.lorentzQ_add
#print axioms OPH.CommonWorldJointAction.lorentzB_nondegenerate
#print axioms OPH.CommonWorldJointAction.clockAction_congr
#print axioms OPH.CommonWorldJointAction.clockAction_expansion
#print axioms OPH.CommonWorldJointAction.clockStationary_firstVariation_zero
#print axioms OPH.CommonWorldJointAction.clockFirstVariation_interior
#print axioms OPH.CommonWorldJointAction.clockStationary_of_uniform
#print axioms OPH.CommonWorldJointAction.clockStationary_consecutive
#print axioms OPH.CommonWorldJointAction.clockStationary_iff_uniform
#print axioms OPH.CommonWorldJointAction.uniform_path_formula
#print axioms OPH.CommonWorldJointAction.joinedPath_eq
#print axioms OPH.CommonWorldJointAction.joinedPath_eq_stepEvent
#print axioms OPH.CommonWorldJointAction.joinedPath_zero
#print axioms OPH.CommonWorldJointAction.joinedPath_increment
#print axioms OPH.CommonWorldJointAction.joinedPath_uniform
#print axioms OPH.CommonWorldJointAction.worldline_clockStationary
#print axioms OPH.CommonWorldJointAction.clockStationary_fixed_endpoints
#print axioms OPH.CommonWorldJointAction.worldline_is_stationary_path
#print axioms OPH.CommonWorldJointAction.jointAction_maxwell_shift
#print axioms OPH.CommonWorldJointAction.jointAction_worldline_shift
#print axioms OPH.CommonWorldJointAction.maxwellSector_iff
#print axioms OPH.CommonWorldJointAction.jointAction_maxwell_sector
#print axioms OPH.CommonWorldJointAction.worldlineSector_iff_clockStationary
#print axioms OPH.CommonWorldJointAction.jointAction_worldline_sector
#print axioms OPH.CommonWorldJointAction.worldline_sector_solution
#print axioms OPH.CommonWorldJointAction.committed_configuration_stationary
#print axioms OPH.CommonWorldJointAction.jointAction_one_principle
#print axioms OPH.CommonWorldJointAction.jointAction_receipt
#print axioms OPH.CommonWorldJointAction.jointAction_sector_decoupling
#print axioms OPH.CommonWorldJointAction.committedJoinedWitness_increment
#print axioms OPH.CommonWorldJointAction.committedJoinedWitness_equal_increments
#print axioms OPH.CommonWorldJointAction.kinkPath_increments
#print axioms OPH.CommonWorldJointAction.kinkPath_nonuniform
#print axioms OPH.CommonWorldJointAction.kink_firstVariation_eq_two
#print axioms OPH.CommonWorldJointAction.kinkPath_not_stationary
#print axioms OPH.CommonWorldJointAction.stationarity_characterization_nonvacuous
#print axioms OPH.CommonWorldJointAction.spatialVariation_negative_action
#print axioms OPH.CommonWorldJointAction.stationarity_not_minimality
#print axioms OPH.CommonWorldJointAction.jointAction_not_forcing_calibration
#print axioms OPH.CommonWorldJointAction.kinkPath_zero
#print axioms OPH.CommonWorldJointAction.kinkPath_succ
