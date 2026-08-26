import ScaledMaxwellStability
import Geometry.WorldlineHopTransport

set_option autoImplicit false

open scoped BigOperators Matrix

namespace OPH.SeamChargeContinuity

open OPH.SeamCurrentCarrierQuotient
open OPH.DiscreteCoulombGreen
open OPH.PositionSpaceMaxwellAction
open OPH.LocalFaceMaxwellAction
open OPH.TemporalMaxwellEvolution
open OPH.ScaledMaxwellStability
open OPH.PortChargeMinimalCoupling
open OPH.WorldlineHopTransport

/-!
# Consistency of the sourced field equations on the committed carrier
(issues #740, #733)

WHAT IS PROVED.  Exact finite real linear algebra on the committed
twelve-port, thirty-seam complex, reusing the committed port boundary
`realBoundary` (∂ = dᵀ), the committed scaled staggered packet
`electricFieldScaled`, `magneticField`, `AmpereEvolutionScaled`, the
committed hopping sources `hoppingLoad`, `hoppingCurrent` of a declared
`HoppingPath`, the committed continuity theorem `hopping_continuity`, the
committed Gauss step equivalence `gauss_step_iff_scaled` and propagation
`gauss_propagation_scaled`, and the committed canonical Coulomb field
`realCoulombField` with its Gauss identity `realCoulombField_gauss`.

(A) Continuity of the hopping current, pointwise
(`hoppingCurrent_continuity`).  At every port `p` and step `n`,
`ρ (n+1) p - ρ n p = -(h * (∂ J n) p)` for the committed hopping load
`ρ = hoppingLoad q γ` and hopping current `J = hoppingCurrent q h γ`,
`h ≠ 0`.  This is the pointwise reading of the committed
`hopping_continuity`; nothing new is assumed.  A finite family of hopping
charges (`familyLoad`, `familyCurrent`) satisfies the same continuity
equation (`family_continuity`) and has total load `∑ i, q i` at every step
(`familyLoad_total`).

(B) Gauss propagation under the sourced Ampere step
(`gauss_propagates_under_sourced_ampere`): the forward direction of the
committed `gauss_step_iff_scaled`, stated as a one-step propagation with
the continuity residual as hypothesis.  Composed with (A) this gives the
hopping-family propagation `family_gauss_propagation` and, for the join of
a declared background source pair `(ρ_b, J_b)` obeying continuity with one
hopping charge, `hopping_gauss_propagation_with_background`, in the exact
shape `J_b + hoppingCurrent`, `ρ_b + hoppingLoad` of the committed
`transported_field_equations`.

(C) HEADLINE, solvability of the sourced system
(`sourced_maxwell_solvable_iff`).  For declared sources `(ρ, J)` and
`h ≠ 0`: a field history `(A, φ)` solving the committed scaled Ampere
evolution with current `J` and the committed Gauss constraint
`∂ (E n) = ρ n` at EVERY step exists if and only if the initial load is
neutral, `∑ p, ρ 0 p = 0`, and the sources obey the committed continuity
equation at every step.  Direction (→): neutrality is the committed
`realBoundary_total`, continuity is the committed converse direction of
`gauss_step_iff_scaled`.  Direction (←): the explicit history
`coulombSourcedHistory h (ρ 0) J` in the gauge `φ = 0`, built by the
sourced second-order recursion `sourcedScaledA` from the initial pair
`A 0 = 0`, `A 1 = -(h • realCoulombField (ρ 0))`; its Ampere law is
`sourcedScaledA_ampere`, its initial electric field is the committed
Coulomb field (`coulombSourcedHistory_electric_zero`), and Gauss at every
step follows by (B).  The committed zero-current demonstration history
`demoScaledA` is the `J = 0`, `A 0 = 0`, `A 1 = demoInitial` member of the
same recursion (`demoScaledA_eq_sourced`).

(D) The closed-carrier obstruction (`gauss_consistent_load_neutral`,
`single_hopping_charge_obstruction`,
`single_hopping_charge_solvable_iff_zero`).  Because `∂` conserves total
load, a Gauss-consistent seam field exists for the hopping load of ONE
charge only if `q = 0`; the single-charge Gauss propagation
(`single_hopping_charge_gauss_propagation`) is therefore correct but its
hypothesis bundle is inhabited only at `q = 0`.  The non-vacuous
statement on the closed committed carrier is the neutral family or the
neutralised background: `hopping_family_solvable_iff` states that the
sourced system with a family of hopping charges is solvable exactly when
`∑ i, q i = 0`.  A consumer that joins one worldline to the field sector
must carry a neutral partner or a neutralising background load; no
theorem here supplies one.

(E) Non-vacuity and numeric verification.
  * Sign-consistency check of continuity: the integer load and current of
    the unit forward crossing of every seam satisfy the continuity equation
    with the committed `incidenceZ`, by kernel `decide` over all thirty
    seams and twelve ports (`crossing_continuity_table`), and the real
    statement for `q = h = 1` is recovered by casting
    (`crossing_continuity_from_table`), independently of
    `hopping_continuity`.  This check is DEFINITIONAL: unfolding
    `crossingLoadZ`, `crossingCurrentZ` and `incidenceZ` reduces it to
    `([p = R e] - [p = L e]) - ([p = R e] - [p = L e]) = 0`, which holds
    for any seam table; it fixes the sign convention between
    `crossingCurrentZ` and `incidenceZ` and carries no information about
    the committed table entries.  The checks that depend on the committed
    table content are `pairPotentialZ_green` and `pairPotentialZ_gauss`
    below.
  * The neutral pair: charge `1` crossing seam `0` (ports `0 → 1`) and
    charge `-1` crossing seam `29` (ports `10 → 11`), both through the
    committed `crossingPath`.  Its step-`0` load is the integer vector
    `pairLoadZ`, the committed Green matrix maps it to the integer potential
    `pairPotentialZ` over `180` (`pairPotentialZ_green`, kernel `decide`),
    and the coboundary of that potential has boundary `180` times the load
    on the committed table (`pairPotentialZ_gauss`, kernel `decide`): the
    Coulomb datum is Gauss-consistent by table arithmetic as well as by the
    committed identity.  Its value on seam `0` is `-1/6`
    (`pairCoulomb_seam_zero`).
  * `neutralPairBundle` is an explicit inhabitant of the committed
    `ScaledMaxwellBundle` at `h = 1/2`, `Λ = 6` whose current is the
    nonzero hopping current of the pair (`neutralPairBundle_nonvacuous`:
    `J 0 0 = -2`, `ρ 0 0 = 1`, `E 0 0 = -1/6`).  The committed demonstration
    bundle has zero current; this one is sourced.  The unconditional
    clauses of the committed receipts on `ScaledMaxwellBundle` (Gauss at
    every step, the Gauss step equivalence, Faraday, the energy balance
    with its work term, gauge invariance, the action characterisation, and
    the face-local balance) apply to it; the clauses guarded by `S.J = 0`
    (energy conservation, the wave law, the uniform bounds) are inert on
    this sourced inhabitant, and no boundedness is asserted for it.
  * `neutralPairBundle_with_background` reads the same inhabitant in the
    source shape of `hopping_gauss_propagation_with_background`: the charge
    `1` crossing seam `0` is the hopping charge and the charge `-1`
    crossing seam `29` is the declared background `(ρ_b, J_b)`
    (`pairLoad_eq_background`, `pairCurrent_eq_background`), so the
    with-background hypothesis bundle has an explicit inhabitant.

SIGN AND ORIENTATION CONVENTION, stated once.  The committed
`realBoundary c p = ∑ e, ([p = seamRight e] - [p = seamLeft e]) * c e` is
the net INFLOW at `p` of a seam value counted positive in the direction
`seamLeft e → seamRight e` (`incidenceZ` credits `seamRight`, debits
`seamLeft`; `seamLeft e < seamRight e` by `seam_table_sound`).  The
committed continuity form `ρ (n+1) - ρ n + h • ∂ (J n) = 0` therefore says
that the load at a port decreases by `h` times the net inflow of `J` in
that orientation; consistently, the committed `hoppingCurrent` of a
forward hop `seamLeft e → seamRight e` carries the value `-(q / h)` on `e`
(`hoppingCurrent_forward`).  Equivalently: `realBoundary J` is the net
OUTFLOW (the divergence) of a charge flow whose seam values are counted
positive in the direction `seamRight e → seamLeft e`, and with that
reading the committed equations take the textbook shape
`ρ (n+1) - ρ n = -(h • div J)`.  On the field side the same operator
gives, for a static field `E = -(d φ)`, `∂ E = ρ ↔ L φ = -ρ` with
`L = ∂ d` the committed nonnegative Laplacian
(`static_gauss_iff_laplacian`); so the committed scalar potential of a
positive committed load is `-(G ρ)` plus a constant, while the committed
canonical seam field of the same load is `+ d (G ρ)`
(`realCoulombField`).  Both facts are consequences of the committed
definitions.  Which sign of the committed load is the physical charge
sign is part of the physical attachment row and is not decided here.

DIRECTIONS OF THE EQUIVALENCES.  `sourced_maxwell_solvable_iff`,
`hopping_family_solvable_iff`, `single_hopping_charge_solvable_iff_zero`,
and `static_gauss_iff_laplacian` are stated as `↔` and both directions are
proved; the propagation theorems are one-directional implications whose
hypotheses are listed in the statements.

DECLARED INPUTS.  The step `h ≠ 0` is declared and is not a clock; the
charge `q` and the hopping paths are declared; the background sources
`(ρ_b, J_b)` are declared; the kinetic term behind `AmpereEvolutionScaled`
is the committed declared enrichment.  No physical attachment of ports,
seams, loads, or currents is made or used.

ROWS CITED.  OL-N1 (the physical attachment of the committed carrier and
its sources) is cited at scope and NOT discharged.  The source clock and
duration row, the coupled-action row, and PR-53/PR-54/PR-15 stay open and
are not consumed.  No frozen prediction is added.

WHAT IS NOT PROVED.  No equation of motion for the hopping path; no
selection of the path by the field; no continuum limit; no identification
of the committed load sign with a physical charge sign; no clock.  The
solvability theorem produces one history in the gauge `φ = 0`; uniqueness
up to gauge and initial data is not stated here.

Axiom audit.  Every proof composes the committed receipts with exact real
linear algebra and kernel `decide` on the committed integer tables; the
module adds no project axiom and uses no native decision procedure.  The
audit lines at the end of the file show at most `propext`,
`Classical.choice`, and `Quot.sound`.
-/

noncomputable section

/-! ## (A) Continuity of the hopping current, pointwise -/

/-- **Continuity of the hopping current, pointwise.**  For `h ≠ 0`, at every
port `p` and step `n`, the change of the committed hopping load across the
step `n → n+1` equals minus `h` times the committed port boundary
(incidence transpose, net inflow in the `seamLeft → seamRight` orientation)
of the committed hopping current.  Pointwise reading of the committed
`hopping_continuity`. -/
theorem hoppingCurrent_continuity (q h : ℝ) (hh : h ≠ 0) (γ : HoppingPath) (n : ℕ)
    (p : Fin 12) :
    hoppingLoad q γ (n + 1) p - hoppingLoad q γ n p =
      -(h * realBoundary (hoppingCurrent q h γ n) p) := by
  have hc := congrFun (hopping_continuity q h hh γ n) p
  simp only [Pi.add_apply, Pi.sub_apply, Pi.smul_apply, smul_eq_mul, Pi.zero_apply] at hc
  linear_combination hc

/-- Port load of a finite family of declared hopping charges: the sum of
the committed hopping loads. -/
def familyLoad {k : ℕ} (q : Fin k → ℝ) (γ : Fin k → HoppingPath) (n : ℕ) :
    Fin 12 → ℝ :=
  ∑ i : Fin k, hoppingLoad (q i) (γ i) n

/-- Seam current of a finite family of declared hopping charges: the sum of
the committed hopping currents. -/
def familyCurrent {k : ℕ} (q : Fin k → ℝ) (h : ℝ) (γ : Fin k → HoppingPath) (n : ℕ) :
    Fin 30 → ℝ :=
  ∑ i : Fin k, hoppingCurrent (q i) h (γ i) n

theorem familyLoad_apply {k : ℕ} (q : Fin k → ℝ) (γ : Fin k → HoppingPath) (n : ℕ)
    (p : Fin 12) :
    familyLoad q γ n p = ∑ i : Fin k, (if p = (γ i).port n then q i else 0) := by
  unfold familyLoad
  rw [Finset.sum_apply]
  exact Finset.sum_congr rfl fun i _ ↦ hoppingLoad_apply (q i) (γ i) n p

theorem familyCurrent_apply {k : ℕ} (q : Fin k → ℝ) (h : ℝ) (γ : Fin k → HoppingPath)
    (n : ℕ) (e : Fin 30) :
    familyCurrent q h γ n e = ∑ i : Fin k, hoppingCurrent (q i) h (γ i) n e := by
  unfold familyCurrent
  rw [Finset.sum_apply]

/-- The total load of a family is the sum of its charges at every step. -/
theorem familyLoad_total {k : ℕ} (q : Fin k → ℝ) (γ : Fin k → HoppingPath) (n : ℕ) :
    (∑ p : Fin 12, familyLoad q γ n p) = ∑ i : Fin k, q i := by
  unfold familyLoad
  simp only [Finset.sum_apply]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun i _ ↦ hoppingLoad_total (q i) (γ i) n

/-- **Continuity of a family of hopping charges.**  For `h ≠ 0`, the family
sources satisfy the committed scaled continuity equation identically. -/
theorem family_continuity {k : ℕ} (q : Fin k → ℝ) (h : ℝ) (hh : h ≠ 0)
    (γ : Fin k → HoppingPath) (n : ℕ) :
    familyLoad q γ (n + 1) - familyLoad q γ n +
      h • realBoundary (familyCurrent q h γ n) = 0 := by
  unfold familyLoad familyCurrent
  rw [map_sum, Finset.smul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
  exact Finset.sum_eq_zero fun i _ ↦ hopping_continuity (q i) h hh (γ i) n

/-! ## Sign convention on the field side -/

/-- For a static seam field `E = -(d φ)`, the committed Gauss constraint
`∂ E = ρ` is the Poisson equation `L φ = -ρ` for the committed nonnegative
Laplacian `L = ∂ d`.  Both directions. -/
theorem static_gauss_iff_laplacian (φ ρ : Fin 12 → ℝ) :
    realBoundary (-realCoboundary φ) = ρ ↔ realLaplacian φ = -ρ := by
  rw [map_neg, realLaplacian_apply, neg_eq_iff_eq_neg]

/-! ## (B) Gauss propagation under the sourced Ampere step -/

/-- **Gauss propagates under the sourced Ampere step.**  If the Gauss
constraint `∂ (E n) = ρ n` holds at step `n`, the committed scaled Ampere
step with current `J` is applied, and the sources obey the committed
continuity equation at step `n`, then the Gauss constraint holds at step
`n+1`.  Forward direction of the committed `gauss_step_iff_scaled`; the
converse (Gauss at `n` and `n+1` forces continuity at `n`) is the other
direction of that committed equivalence. -/
theorem gauss_propagates_under_sourced_ampere (h : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) (ρ : ℕ → Fin 12 → ℝ)
    (hAmp : AmpereEvolutionScaled h A φ J) (n : ℕ)
    (hn : realBoundary (electricFieldScaled h A φ n) = ρ n)
    (hcont : ρ (n + 1) - ρ n + h • realBoundary (J n) = 0) :
    realBoundary (electricFieldScaled h A φ (n + 1)) = ρ (n + 1) :=
  (gauss_step_iff_scaled h A φ J ρ hAmp n hn).mpr hcont

/-- One Gauss step for a family of hopping charges. -/
theorem family_gauss_step {k : ℕ} (q : Fin k → ℝ) (h : ℝ) (hh : h ≠ 0)
    (γ : Fin k → HoppingPath) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (hAmp : AmpereEvolutionScaled h A φ (familyCurrent q h γ)) (n : ℕ)
    (hn : realBoundary (electricFieldScaled h A φ n) = familyLoad q γ n) :
    realBoundary (electricFieldScaled h A φ (n + 1)) = familyLoad q γ (n + 1) :=
  gauss_propagates_under_sourced_ampere h A φ (familyCurrent q h γ) (familyLoad q γ)
    hAmp n hn (family_continuity q h hh γ n)

/-- **Gauss propagation for a family of hopping charges.**  Any field
history solving the committed scaled Ampere evolution with the family
current, starting from data satisfying the Gauss constraint at the family
load, satisfies the Gauss constraint at every step. -/
theorem family_gauss_propagation {k : ℕ} (q : Fin k → ℝ) (h : ℝ) (hh : h ≠ 0)
    (γ : Fin k → HoppingPath) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (hAmp : AmpereEvolutionScaled h A φ (familyCurrent q h γ))
    (h0 : realBoundary (electricFieldScaled h A φ 0) = familyLoad q γ 0) :
    ∀ n : ℕ, realBoundary (electricFieldScaled h A φ n) = familyLoad q γ n :=
  gauss_propagation_scaled h A φ (familyCurrent q h γ) (familyLoad q γ) hAmp h0
    (family_continuity q h hh γ)

/-- **Gauss propagation for one hopping charge joined to a declared
background.**  With background sources `(ρ_b, J_b)` obeying the committed
continuity equation, any field history solving the committed scaled Ampere
evolution with the current `J_b + hoppingCurrent q h γ`, starting from
data with `∂ (E 0) = ρ_b 0 + hoppingLoad q γ 0`, satisfies
`∂ (E n) = ρ_b n + hoppingLoad q γ n` at every step.  This is the source
shape of the committed `transported_field_equations`; the hopping path
`hoppingPath w` of any seam-step worldline `w` is an instance of `γ`. -/
theorem hopping_gauss_propagation_with_background (q h : ℝ) (hh : h ≠ 0)
    (γ : HoppingPath) (ρb : ℕ → Fin 12 → ℝ) (Jb : ℕ → Fin 30 → ℝ)
    (hcont : ∀ n, ρb (n + 1) - ρb n + h • realBoundary (Jb n) = 0)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (hAmp : AmpereEvolutionScaled h A φ (fun m ↦ Jb m + hoppingCurrent q h γ m))
    (h0 : realBoundary (electricFieldScaled h A φ 0) = ρb 0 + hoppingLoad q γ 0) :
    ∀ n : ℕ, realBoundary (electricFieldScaled h A φ n) = ρb n + hoppingLoad q γ n := by
  refine gauss_propagation_scaled h A φ _ (fun n ↦ ρb n + hoppingLoad q γ n) hAmp h0 ?_
  intro n
  have h1 := hcont n
  have h2 := hopping_continuity q h hh γ n
  rw [map_add, smul_add]
  calc ρb (n + 1) + hoppingLoad q γ (n + 1) - (ρb n + hoppingLoad q γ n) +
        (h • realBoundary (Jb n) + h • realBoundary (hoppingCurrent q h γ n))
      = (ρb (n + 1) - ρb n + h • realBoundary (Jb n)) +
        (hoppingLoad q γ (n + 1) - hoppingLoad q γ n +
          h • realBoundary (hoppingCurrent q h γ n)) := by abel
    _ = 0 := by rw [h1, h2, add_zero]

/-- Gauss propagation for a single hopping charge with no background.  The
statement is correct as an implication; by
`single_hopping_charge_obstruction` its hypothesis `h0` is satisfiable only
at `q = 0`, so on the closed committed carrier this corollary has no
nonzero-charge inhabitant. -/
theorem single_hopping_charge_gauss_propagation (q h : ℝ) (hh : h ≠ 0)
    (γ : HoppingPath) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (hAmp : AmpereEvolutionScaled h A φ (hoppingCurrent q h γ))
    (h0 : realBoundary (electricFieldScaled h A φ 0) = hoppingLoad q γ 0) :
    ∀ n : ℕ, realBoundary (electricFieldScaled h A φ n) = hoppingLoad q γ n :=
  gauss_propagation_scaled h A φ (hoppingCurrent q h γ) (hoppingLoad q γ) hAmp h0
    (hopping_continuity q h hh γ)

/-! ## (D) The closed-carrier obstruction -/

/-- A Gauss-consistent load on the closed committed carrier is neutral:
the committed boundary conserves total load. -/
theorem gauss_consistent_load_neutral (E : Fin 30 → ℝ) (ρ : Fin 12 → ℝ)
    (hG : realBoundary E = ρ) : (∑ p : Fin 12, ρ p) = 0 := by
  rw [← hG]
  exact realBoundary_total E

/-- **Single-charge obstruction.**  No seam field has boundary equal to the
hopping load of one charge unless the charge vanishes. -/
theorem single_hopping_charge_obstruction (q : ℝ) (γ : HoppingPath) (n : ℕ)
    (E : Fin 30 → ℝ) (hG : realBoundary E = hoppingLoad q γ n) : q = 0 := by
  have h := gauss_consistent_load_neutral E _ hG
  rw [hoppingLoad_total] at h
  exact h

/-- Family obstruction: a Gauss-consistent family load has zero total
charge. -/
theorem family_obstruction {k : ℕ} (q : Fin k → ℝ) (γ : Fin k → HoppingPath) (n : ℕ)
    (E : Fin 30 → ℝ) (hG : realBoundary E = familyLoad q γ n) :
    (∑ i : Fin k, q i) = 0 := by
  have h := gauss_consistent_load_neutral E _ hG
  rw [familyLoad_total] at h
  exact h

/-- Background obstruction: with a background load, Gauss-consistency forces
the background total to be `-q`. -/
theorem background_obstruction (q : ℝ) (γ : HoppingPath) (ρb : Fin 12 → ℝ) (n : ℕ)
    (E : Fin 30 → ℝ) (hG : realBoundary E = ρb + hoppingLoad q γ n) :
    (∑ p : Fin 12, ρb p) + q = 0 := by
  have h := gauss_consistent_load_neutral E _ hG
  simp only [Pi.add_apply] at h
  rw [Finset.sum_add_distrib, hoppingLoad_total] at h
  exact h

/-! ## (C) Existence: the sourced recursion and the Coulomb initial datum -/

/-- The sourced second-order recursion in the gauge `φ = 0`:
`A (n+2) = 2 A (n+1) - A n - h² CᵀC (A (n+1)) + h² J n` from a declared
initial pair.  The committed zero-current `demoScaledA` is its `J = 0`,
`A 0 = 0`, `A 1 = demoInitial` member (`demoScaledA_eq_sourced`). -/
def sourcedScaledA (h : ℝ) (J : ℕ → Fin 30 → ℝ) (A0 A1 : Fin 30 → ℝ) : ℕ → Fin 30 → ℝ
  | 0 => A0
  | 1 => A1
  | n + 2 => (2 : ℝ) • sourcedScaledA h J A0 A1 (n + 1) - sourcedScaledA h J A0 A1 n -
      (h ^ 2) • localMaxwellOperator (sourcedScaledA h J A0 A1 (n + 1)) + (h ^ 2) • J n

theorem sourcedScaledA_zero (h : ℝ) (J : ℕ → Fin 30 → ℝ) (A0 A1 : Fin 30 → ℝ) :
    sourcedScaledA h J A0 A1 0 = A0 := by
  simp only [sourcedScaledA]

theorem sourcedScaledA_one (h : ℝ) (J : ℕ → Fin 30 → ℝ) (A0 A1 : Fin 30 → ℝ) :
    sourcedScaledA h J A0 A1 1 = A1 := by
  simp only [sourcedScaledA]

/-- **The sourced recursion solves the committed scaled Ampere evolution**
with current `J` in the gauge `φ = 0`, for every `h ≠ 0` and every
initial pair. -/
theorem sourcedScaledA_ampere (h : ℝ) (hh : h ≠ 0) (J : ℕ → Fin 30 → ℝ)
    (A0 A1 : Fin 30 → ℝ) :
    AmpereEvolutionScaled h (sourcedScaledA h J A0 A1) (fun _ ↦ 0) J := by
  intro n
  rw [electricFieldScaled_temporal_gauge, electricFieldScaled_temporal_gauge]
  have hlm : faceCodifferential (magneticField (sourcedScaledA h J A0 A1) (n + 1)) =
      localMaxwellOperator (sourcedScaledA h J A0 A1 (n + 1)) := rfl
  rw [hlm]
  have hidx : sourcedScaledA h J A0 A1 (n + 1 + 1) =
      (2 : ℝ) • sourcedScaledA h J A0 A1 (n + 1) - sourcedScaledA h J A0 A1 n -
        (h ^ 2) • localMaxwellOperator (sourcedScaledA h J A0 A1 (n + 1)) +
          (h ^ 2) • J n := by
    simp only [sourcedScaledA]
  rw [hidx]
  have hinv : h * h⁻¹ = 1 := mul_inv_cancel₀ hh
  funext e
  simp only [Pi.sub_apply, Pi.neg_apply, Pi.smul_apply, Pi.add_apply, smul_eq_mul]
  linear_combination (h * (localMaxwellOperator (sourcedScaledA h J A0 A1 (n + 1)) e
    - J n e)) * hinv

/-- The committed zero-current demonstration history is a member of the
sourced recursion. -/
theorem demoScaledA_eq_sourced (h : ℝ) :
    demoScaledA h = sourcedScaledA h (fun _ ↦ 0) 0 demoInitial := by
  funext n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    match n with
    | 0 => simp only [demoScaledA, sourcedScaledA]
    | 1 => simp only [demoScaledA, sourcedScaledA]
    | m + 2 =>
      have h1 := ih (m + 1) (by omega)
      have h0 := ih m (by omega)
      simp only [demoScaledA, sourcedScaledA, h1, h0, smul_zero, add_zero]

/-- The Coulomb-started sourced history: `A 0 = 0`,
`A 1 = -(h • realCoulombField ρ0)`, then the sourced recursion with
current `J`, in the gauge `φ = 0`. -/
def coulombSourcedHistory (h : ℝ) (ρ0 : Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) :
    ℕ → Fin 30 → ℝ :=
  sourcedScaledA h J 0 (-(h • realCoulombField ρ0))

theorem coulombSourcedHistory_ampere (h : ℝ) (hh : h ≠ 0) (ρ0 : Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) :
    AmpereEvolutionScaled h (coulombSourcedHistory h ρ0 J) (fun _ ↦ 0) J :=
  sourcedScaledA_ampere h hh J 0 _

/-- The initial electric field of the Coulomb-started history is the
committed canonical Coulomb field of `ρ0`. -/
theorem coulombSourcedHistory_electric_zero (h : ℝ) (hh : h ≠ 0) (ρ0 : Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) :
    electricFieldScaled h (coulombSourcedHistory h ρ0 J) (fun _ ↦ 0) 0 =
      realCoulombField ρ0 := by
  rw [electricFieldScaled_temporal_gauge]
  unfold coulombSourcedHistory
  rw [sourcedScaledA_one, sourcedScaledA_zero, sub_zero, smul_neg, smul_smul,
    inv_mul_cancel₀ hh, one_smul, neg_neg]

/-- For a neutral initial load the Coulomb-started history satisfies the
Gauss constraint at step `0`. -/
theorem coulombSourcedHistory_gauss_zero (h : ℝ) (hh : h ≠ 0) (ρ0 : Fin 12 → ℝ)
    (hρ : (∑ p : Fin 12, ρ0 p) = 0) (J : ℕ → Fin 30 → ℝ) :
    realBoundary (electricFieldScaled h (coulombSourcedHistory h ρ0 J) (fun _ ↦ 0) 0) =
      ρ0 := by
  rw [coulombSourcedHistory_electric_zero h hh ρ0 J]
  exact realCoulombField_gauss ρ0 hρ

/-- **HEADLINE: solvability of the sourced field equations.**  For declared
sources `(ρ, J)` and `h ≠ 0`, a field history solving the committed scaled
Ampere evolution with current `J` and the committed Gauss constraint at
every step exists if and only if the initial load is neutral and the
sources obey the committed continuity equation at every step.  (→) uses
the committed `realBoundary_total` and the converse direction of the
committed `gauss_step_iff_scaled`; (←) is the explicit Coulomb-started
sourced history in the gauge `φ = 0`. -/
theorem sourced_maxwell_solvable_iff (h : ℝ) (hh : h ≠ 0) (ρ : ℕ → Fin 12 → ℝ)
    (J : ℕ → Fin 30 → ℝ) :
    (∃ (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ), AmpereEvolutionScaled h A φ J ∧
        ∀ n, realBoundary (electricFieldScaled h A φ n) = ρ n) ↔
      ((∑ p : Fin 12, ρ 0 p) = 0 ∧
        ∀ n, ρ (n + 1) - ρ n + h • realBoundary (J n) = 0) := by
  constructor
  · rintro ⟨A, φ, hAmp, hG⟩
    refine ⟨gauss_consistent_load_neutral _ _ (hG 0), fun n ↦ ?_⟩
    exact (gauss_step_iff_scaled h A φ J ρ hAmp n (hG n)).mp (hG (n + 1))
  · rintro ⟨hneutral, hcont⟩
    refine ⟨coulombSourcedHistory h (ρ 0) J, fun _ ↦ 0,
      coulombSourcedHistory_ampere h hh (ρ 0) J, ?_⟩
    exact gauss_propagation_scaled h _ _ J ρ (coulombSourcedHistory_ampere h hh (ρ 0) J)
      (coulombSourcedHistory_gauss_zero h hh (ρ 0) hneutral J) hcont

/-- **Solvability for a family of hopping charges.**  The sourced system
with the family sources is solvable exactly when the total charge
vanishes. -/
theorem hopping_family_solvable_iff {k : ℕ} (q : Fin k → ℝ) (h : ℝ) (hh : h ≠ 0)
    (γ : Fin k → HoppingPath) :
    (∃ (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ),
        AmpereEvolutionScaled h A φ (familyCurrent q h γ) ∧
        ∀ n, realBoundary (electricFieldScaled h A φ n) = familyLoad q γ n) ↔
      (∑ i : Fin k, q i) = 0 := by
  rw [sourced_maxwell_solvable_iff h hh (familyLoad q γ) (familyCurrent q h γ),
    familyLoad_total]
  exact ⟨fun hs ↦ hs.1, fun hq ↦ ⟨hq, family_continuity q h hh γ⟩⟩

/-- **Solvability for one hopping charge.**  The sourced system with the
hopping sources of a single charge is solvable exactly when the charge is
zero: on the closed committed carrier one monopole has no Gauss-consistent
field history. -/
theorem single_hopping_charge_solvable_iff_zero (q h : ℝ) (hh : h ≠ 0) (γ : HoppingPath) :
    (∃ (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ),
        AmpereEvolutionScaled h A φ (hoppingCurrent q h γ) ∧
        ∀ n, realBoundary (electricFieldScaled h A φ n) = hoppingLoad q γ n) ↔
      q = 0 := by
  rw [sourced_maxwell_solvable_iff h hh (hoppingLoad q γ) (hoppingCurrent q h γ),
    hoppingLoad_total]
  exact ⟨fun hs ↦ hs.1, fun hq ↦ ⟨hq, hopping_continuity q h hh γ⟩⟩

/-! ## (E) Non-vacuity: table verification of the unit crossing -/

/-- Integer load of the unit forward crossing of seam `e`: the indicator of
`seamLeft e` at step `0`, of `seamRight e` afterwards. -/
def crossingLoadZ (e : Fin 30) (n : ℕ) (p : Fin 12) : ℤ :=
  if p = (if n = 0 then seamLeft e else seamRight e) then 1 else 0

/-- Integer current of the unit forward crossing of seam `e` at `q = h = 1`:
`-1` on `e`, zero elsewhere, in the committed sign convention. -/
def crossingCurrentZ (e e' : Fin 30) : ℤ :=
  if e' = e then -1 else 0

set_option maxRecDepth 32768 in
/-- **Sign-consistency check of continuity.**  For every seam and every
port, the integer crossing load and current satisfy
`ρ 1 p - ρ 0 p + ∑ e', incidenceZ e' p * J e' = 0` with the committed
`incidenceZ`.  Kernel `decide` over thirty seams and twelve ports.  The
identity is definitional (it reduces to
`([p = R e] - [p = L e]) - ([p = R e] - [p = L e]) = 0` for any seam
table), so it verifies the sign convention between `crossingCurrentZ` and
`incidenceZ` and not the entries of the committed table; the table-content
checks are `pairPotentialZ_green` and `pairPotentialZ_gauss`. -/
theorem crossing_continuity_table :
    ∀ (e : Fin 30) (p : Fin 12),
      crossingLoadZ e 1 p - crossingLoadZ e 0 p +
        ∑ e' : Fin 30, incidenceZ e' p * crossingCurrentZ e e' = 0 := by
  decide

/-- The committed real boundary of an integer-cast seam field is the cast
of the integer incidence sum. -/
theorem realBoundary_intCast (c : Fin 30 → ℤ) (p : Fin 12) :
    realBoundary (fun e ↦ (c e : ℝ)) p = ((∑ e : Fin 30, incidenceZ e p * c e : ℤ) : ℝ) := by
  rw [realBoundary_apply]
  push_cast
  refine Finset.sum_congr rfl fun e _ ↦ ?_
  unfold incidenceZ
  push_cast
  split_ifs <;> ring

theorem hoppingLoad_crossingPath_cast (e : Fin 30) (n : ℕ) (p : Fin 12) :
    hoppingLoad 1 (crossingPath e) n p = (crossingLoadZ e n p : ℝ) := by
  rw [hoppingLoad_apply]
  unfold crossingLoadZ
  have hport : (crossingPath e).port n = if n = 0 then seamLeft e else seamRight e := rfl
  rw [hport]
  split_ifs <;> simp

theorem hoppingCurrent_crossingPath_cast (e : Fin 30) :
    hoppingCurrent 1 1 (crossingPath e) 0 = fun e' ↦ (crossingCurrentZ e e' : ℝ) := by
  rw [hoppingCurrent_forward 1 1 (crossingPath e) 0 e (crossingPath_zero e)
    (crossingPath_one e)]
  funext e'
  unfold crossingCurrentZ
  split_ifs <;> simp

/-- **Continuity of the unit crossing from the table.**  The real
continuity equation at `q = h = 1` for the crossing of any seam, recovered
from the integer table check by casting, independently of the committed
`hopping_continuity`. -/
theorem crossing_continuity_from_table (e : Fin 30) :
    hoppingLoad 1 (crossingPath e) 1 - hoppingLoad 1 (crossingPath e) 0 +
      (1 : ℝ) • realBoundary (hoppingCurrent 1 1 (crossingPath e) 0) = 0 := by
  funext p
  simp only [Pi.add_apply, Pi.sub_apply, Pi.smul_apply, smul_eq_mul, Pi.zero_apply, one_mul]
  rw [hoppingLoad_crossingPath_cast, hoppingLoad_crossingPath_cast,
    hoppingCurrent_crossingPath_cast, realBoundary_intCast]
  have h := crossing_continuity_table e p
  exact_mod_cast h

/-! ## (E) Non-vacuity: the neutral pair and its Coulomb datum -/

/-- Charges of the neutral pair: `1` and `-1`. -/
def pairCharge : Fin 2 → ℝ := ![1, -1]

/-- Paths of the neutral pair: the committed crossing paths of seam `0`
(ports `0 → 1`, committed `seam_zero_endpoints`) and seam `29`
(ports `10 → 11`, `seam_29_endpoints`). -/
def pairPath : Fin 2 → HoppingPath := ![crossingPath 0, crossingPath 29]

/-- The pair load at every step. -/
def pairLoad (n : ℕ) : Fin 12 → ℝ := familyLoad pairCharge pairPath n

/-- The pair current at step size `h`, step `n`. -/
def pairCurrent (h : ℝ) (n : ℕ) : Fin 30 → ℝ := familyCurrent pairCharge h pairPath n

theorem pairCharge_total : (∑ i : Fin 2, pairCharge i) = 0 := by
  simp [pairCharge, Fin.sum_univ_two]

/-- The pair load is neutral at every step. -/
theorem pairLoad_total (n : ℕ) : (∑ p : Fin 12, pairLoad n p) = 0 := by
  unfold pairLoad
  rw [familyLoad_total]
  exact pairCharge_total

/-- Endpoints of seam `29`: ports `10 → 11`.  Seam `0` is the committed
`seam_zero_endpoints`. -/
theorem seam_29_endpoints : seamLeft 29 = 10 ∧ seamRight 29 = 11 := by
  decide

/-- Integer step-`0` load of the pair: `+1` at port `0`, `-1` at port `10`. -/
def pairLoadZ : Fin 12 → ℤ := ![1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0]

/-- One hundred eighty times the committed Green potential of the pair
load. -/
def pairPotentialZ : Fin 12 → ℤ := ![42, 12, 9, 9, 0, 0, 0, 0, -9, -9, -42, -12]

/-- The committed integer Green table maps the pair load to the pair
potential.  Kernel `decide`. -/
theorem pairPotentialZ_green :
    ∀ p : Fin 12, (∑ r : Fin 12, greenNumZ p r * pairLoadZ r) = pairPotentialZ p := by
  decide

set_option maxRecDepth 32768 in
/-- **Table check of the Gauss datum.**  The incidence boundary of the
coboundary of the pair potential is `180` times the pair load on the
committed table.  Kernel `decide`. -/
theorem pairPotentialZ_gauss :
    ∀ p : Fin 12, (∑ e : Fin 30, incidenceZ e p *
      (pairPotentialZ (seamRight e) - pairPotentialZ (seamLeft e))) = 180 * pairLoadZ p := by
  decide

/-- The real step-`0` pair load is the cast of the integer pair load. -/
theorem pairLoad_zero_cast : pairLoad 0 = fun p ↦ (pairLoadZ p : ℝ) := by
  funext p
  unfold pairLoad
  rw [familyLoad_apply, Fin.sum_univ_two]
  have h0 : (pairPath 0).port 0 = seamLeft 0 := rfl
  have h1 : (pairPath 1).port 0 = seamLeft 29 := rfl
  rw [h0, h1, seam_zero_endpoints.1, seam_29_endpoints.1]
  simp only [pairCharge, pairLoadZ, Matrix.cons_val_zero, Matrix.cons_val_one]
  fin_cases p <;> simp

/-- The committed Coulomb field of the pair load, seam by seam, is the
cast coboundary of the integer pair potential over `180`. -/
theorem pairCoulomb_eq_table (e : Fin 30) :
    realCoulombField (pairLoad 0) e =
      ((pairPotentialZ (seamRight e) - pairPotentialZ (seamLeft e) : ℤ) : ℝ) / 180 := by
  rw [pairLoad_zero_cast]
  unfold realCoulombField
  rw [realCoboundary_apply]
  have hpot : ∀ p : Fin 12, greenMatrixR.mulVec (fun r ↦ (pairLoadZ r : ℝ)) p =
      (pairPotentialZ p : ℝ) / 180 := by
    intro p
    rw [← pairPotentialZ_green p]
    simp only [Matrix.mulVec, dotProduct, greenMatrixR_apply]
    push_cast
    rw [Finset.sum_div]
    exact Finset.sum_congr rfl fun r _ ↦ by ring
  rw [hpot, hpot]
  push_cast
  ring

/-- The Coulomb datum of the pair has boundary equal to the pair load. -/
theorem pairCoulomb_gauss : realBoundary (realCoulombField (pairLoad 0)) = pairLoad 0 :=
  realCoulombField_gauss _ (pairLoad_total 0)

/-- The pair Coulomb field on seam `0` is `-1/6`. -/
theorem pairCoulomb_seam_zero : realCoulombField (pairLoad 0) 0 = -1 / 6 := by
  rw [pairCoulomb_eq_table, seam_zero_endpoints.1, seam_zero_endpoints.2]
  simp only [pairPotentialZ, Matrix.cons_val_zero, Matrix.cons_val_one]
  norm_num

/-- The pair current at `h = 1/2`, step `0`, seam `0` is `-2`. -/
theorem pairCurrent_seam_zero : pairCurrent (1 / 2) 0 0 = -2 := by
  unfold pairCurrent
  rw [familyCurrent_apply, Fin.sum_univ_two]
  have h0 : pairPath 0 = crossingPath 0 := rfl
  have h1 : pairPath 1 = crossingPath 29 := rfl
  have hq0 : pairCharge 0 = 1 := rfl
  have hq1 : pairCharge 1 = -1 := rfl
  rw [h0, h1, hq0, hq1,
    hoppingCurrent_forward 1 (1 / 2) (crossingPath 0) 0 0 (crossingPath_zero 0)
      (crossingPath_one 0),
    hoppingCurrent_forward (-1) (1 / 2) (crossingPath 29) 0 29 (crossingPath_zero 29)
      (crossingPath_one 29)]
  have h29 : (0 : Fin 30) ≠ 29 := by decide
  norm_num [h29]

/-- **The sourced inhabitant.**  The committed `ScaledMaxwellBundle` at
`h = 1/2`, `Λ = 6` with the neutral pair as its sources and the
Coulomb-started sourced history as its field history. -/
def neutralPairBundle : ScaledMaxwellBundle where
  h := 1 / 2
  Λ := 6
  h_pos := by norm_num
  Λ_nonneg := by norm_num
  courant := committed_courant
  courant_strict := by norm_num
  A := coulombSourcedHistory (1 / 2) (pairLoad 0) (pairCurrent (1 / 2))
  phi := fun _ ↦ 0
  rho := pairLoad
  J := pairCurrent (1 / 2)
  ampere := coulombSourcedHistory_ampere (1 / 2) (by norm_num) _ _
  gauss_init := coulombSourcedHistory_gauss_zero (1 / 2) (by norm_num) _
    (pairLoad_total 0) _
  continuity := family_continuity pairCharge (1 / 2) (by norm_num) pairPath

/-- The sourced inhabitant carries a nonzero current, a nonzero load, and
the Coulomb initial field. -/
theorem neutralPairBundle_nonvacuous :
    neutralPairBundle.J 0 0 = -2 ∧ neutralPairBundle.rho 0 0 = 1 ∧
      electricFieldScaled neutralPairBundle.h neutralPairBundle.A neutralPairBundle.phi 0 0 =
        -1 / 6 := by
  refine ⟨pairCurrent_seam_zero, ?_, ?_⟩
  · show pairLoad 0 0 = 1
    rw [pairLoad_zero_cast]
    simp [pairLoadZ]
  · show electricFieldScaled (1 / 2) (coulombSourcedHistory (1 / 2) (pairLoad 0)
      (pairCurrent (1 / 2))) (fun _ ↦ 0) 0 0 = -1 / 6
    rw [coulombSourcedHistory_electric_zero (1 / 2) (by norm_num), pairCoulomb_seam_zero]

/-- Gauss holds at every step on the sourced inhabitant: the committed
receipt applied to the sourced bundle. -/
theorem neutralPairBundle_gauss :
    ∀ n, realBoundary (electricFieldScaled neutralPairBundle.h neutralPairBundle.A
      neutralPairBundle.phi n) = pairLoad n :=
  (scaledMaxwellStability_receipt neutralPairBundle).1

/-- The neutral pair is solvable, as an instance of the family criterion. -/
theorem pair_solvable :
    ∃ (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ),
      AmpereEvolutionScaled (1 / 2) A φ (pairCurrent (1 / 2)) ∧
      ∀ n, realBoundary (electricFieldScaled (1 / 2) A φ n) = pairLoad n :=
  (hopping_family_solvable_iff pairCharge (1 / 2) (by norm_num) pairPath).mpr pairCharge_total

/-- The pair load is the hopping load of the charge `1` on seam `0` on top
of the declared background load of the charge `-1` on seam `29`. -/
theorem pairLoad_eq_background (n : ℕ) :
    pairLoad n = hoppingLoad (-1) (crossingPath 29) n + hoppingLoad 1 (crossingPath 0) n := by
  unfold pairLoad familyLoad
  rw [Fin.sum_univ_two]
  have h0 : pairPath 0 = crossingPath 0 := rfl
  have h1 : pairPath 1 = crossingPath 29 := rfl
  have hq0 : pairCharge 0 = 1 := rfl
  have hq1 : pairCharge 1 = -1 := rfl
  rw [h0, h1, hq0, hq1, add_comm]

/-- The pair current is the hopping current of the charge `1` on seam `0`
on top of the declared background current of the charge `-1` on seam
`29`. -/
theorem pairCurrent_eq_background (h : ℝ) (n : ℕ) :
    pairCurrent h n =
      hoppingCurrent (-1) h (crossingPath 29) n + hoppingCurrent 1 h (crossingPath 0) n := by
  unfold pairCurrent familyCurrent
  rw [Fin.sum_univ_two]
  have h0 : pairPath 0 = crossingPath 0 := rfl
  have h1 : pairPath 1 = crossingPath 29 := rfl
  have hq0 : pairCharge 0 = 1 := rfl
  have hq1 : pairCharge 1 = -1 := rfl
  rw [h0, h1, hq0, hq1, add_comm]

/-- **Explicit inhabitant of the with-background hypothesis bundle.**  The
field history of `neutralPairBundle`, read with the charge `1` crossing
seam `0` as the hopping charge and the charge `-1` crossing seam `29` as
the declared background `(ρ_b, J_b)`, satisfies every hypothesis of
`hopping_gauss_propagation_with_background` at `h = 1/2`, and the
conclusion is Gauss at every step in the background-plus-hopping source
shape. -/
theorem neutralPairBundle_with_background :
    ∀ n : ℕ, realBoundary (electricFieldScaled neutralPairBundle.h neutralPairBundle.A
        neutralPairBundle.phi n) =
      hoppingLoad (-1) (crossingPath 29) n + hoppingLoad 1 (crossingPath 0) n := by
  have hJ : (fun m ↦ hoppingCurrent (-1) (1 / 2) (crossingPath 29) m +
      hoppingCurrent 1 (1 / 2) (crossingPath 0) m) = pairCurrent (1 / 2) := by
    funext m
    rw [pairCurrent_eq_background]
  have hAmp : AmpereEvolutionScaled (1 / 2) neutralPairBundle.A neutralPairBundle.phi
      (fun m ↦ hoppingCurrent (-1) (1 / 2) (crossingPath 29) m +
        hoppingCurrent 1 (1 / 2) (crossingPath 0) m) := by
    rw [hJ]
    exact neutralPairBundle.ampere
  have h0 : realBoundary (electricFieldScaled (1 / 2) neutralPairBundle.A
      neutralPairBundle.phi 0) =
      hoppingLoad (-1) (crossingPath 29) 0 + hoppingLoad 1 (crossingPath 0) 0 := by
    rw [← pairLoad_eq_background]
    exact neutralPairBundle.gauss_init
  exact hopping_gauss_propagation_with_background 1 (1 / 2) (by norm_num) (crossingPath 0)
    (hoppingLoad (-1) (crossingPath 29)) (hoppingCurrent (-1) (1 / 2) (crossingPath 29))
    (hopping_continuity (-1) (1 / 2) (by norm_num) (crossingPath 29))
    neutralPairBundle.A neutralPairBundle.phi hAmp h0

/-! ## Seam-step worldlines -/

/-- The crossing worldline of `WorldlineHopTransport` has the hopping path
of the committed crossing path, so the table verification above covers
it. -/
theorem crossingWorldline_continuity_from_table (e : Fin 30) :
    hoppingLoad 1 (hoppingPath (crossingWorldline e)) 1 -
        hoppingLoad 1 (hoppingPath (crossingWorldline e)) 0 +
      (1 : ℝ) • realBoundary (hoppingCurrent 1 1 (hoppingPath (crossingWorldline e)) 0) = 0 := by
  have hport := (crossingWorldline_generated e).2
  have hload : ∀ n, hoppingLoad 1 (hoppingPath (crossingWorldline e)) n =
      hoppingLoad 1 (crossingPath e) n := by
    intro n
    unfold hoppingLoad
    rw [hport]
  have hcur : hoppingCurrent 1 1 (hoppingPath (crossingWorldline e)) 0 =
      hoppingCurrent 1 1 (crossingPath e) 0 := by
    unfold hoppingCurrent
    rw [hport]
  rw [hload, hload, hcur]
  exact crossing_continuity_from_table e

end

end OPH.SeamChargeContinuity

/- Axiom audit: the committed receipts, exact real linear algebra, and
kernel `decide` on the committed integer tables only.  Expected axioms per
line: at most `propext`, `Classical.choice`, `Quot.sound`.  No native
decision procedure is used. -/

#print axioms OPH.SeamChargeContinuity.hoppingCurrent_continuity
#print axioms OPH.SeamChargeContinuity.family_continuity
#print axioms OPH.SeamChargeContinuity.familyLoad_total
#print axioms OPH.SeamChargeContinuity.static_gauss_iff_laplacian
#print axioms OPH.SeamChargeContinuity.gauss_propagates_under_sourced_ampere
#print axioms OPH.SeamChargeContinuity.family_gauss_propagation
#print axioms OPH.SeamChargeContinuity.hopping_gauss_propagation_with_background
#print axioms OPH.SeamChargeContinuity.single_hopping_charge_gauss_propagation
#print axioms OPH.SeamChargeContinuity.single_hopping_charge_obstruction
#print axioms OPH.SeamChargeContinuity.family_obstruction
#print axioms OPH.SeamChargeContinuity.background_obstruction
#print axioms OPH.SeamChargeContinuity.sourcedScaledA_ampere
#print axioms OPH.SeamChargeContinuity.demoScaledA_eq_sourced
#print axioms OPH.SeamChargeContinuity.coulombSourcedHistory_gauss_zero
#print axioms OPH.SeamChargeContinuity.sourced_maxwell_solvable_iff
#print axioms OPH.SeamChargeContinuity.hopping_family_solvable_iff
#print axioms OPH.SeamChargeContinuity.single_hopping_charge_solvable_iff_zero
#print axioms OPH.SeamChargeContinuity.crossing_continuity_table
#print axioms OPH.SeamChargeContinuity.crossing_continuity_from_table
#print axioms OPH.SeamChargeContinuity.pairPotentialZ_green
#print axioms OPH.SeamChargeContinuity.pairPotentialZ_gauss
#print axioms OPH.SeamChargeContinuity.pairCoulomb_eq_table
#print axioms OPH.SeamChargeContinuity.pairCoulomb_gauss
#print axioms OPH.SeamChargeContinuity.pairCoulomb_seam_zero
#print axioms OPH.SeamChargeContinuity.neutralPairBundle_nonvacuous
#print axioms OPH.SeamChargeContinuity.neutralPairBundle_gauss
#print axioms OPH.SeamChargeContinuity.pair_solvable
#print axioms OPH.SeamChargeContinuity.neutralPairBundle_with_background
#print axioms OPH.SeamChargeContinuity.crossingWorldline_continuity_from_table
