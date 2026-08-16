import Variational.ModeExtremalEnrichment

set_option autoImplicit false

namespace OPH.Variational

/-!
# Complete characterization of corner-invisible enrichments in the
quadratic grammar (V3, issue #739, register row PR-06)

The committed counterfamily of `ModeExtremalEnrichment` shows that the
declared mode-extremality principle does not select a unique enrichment
beyond the registered one-parameter ansatz.  This module upgrades that
audit finding to an exact classification inside one declared grammar.

The declared grammar is the real quadratic two-slot polynomial
`quadPoly α β γ δ ε ζ = α + β x + γ y + δ x y + ε x² + ζ y²` added to
the committed bilinear extension `chainLogLagrangian`.  What IS proved,
all of it exact finite real algebra with no new premises:

* **Completeness (E1).**  A quadratic sentence is corner-invisible
  (vanishes at the four binary corners) exactly when `α = 0`, `δ = 0`,
  `β = -ε`, `γ = -ζ` (`quadPoly_cornerInvisible_iff`), i.e. exactly
  when it is `ε x(x-1) + ζ y(y-1)`
  (`quadPoly_cornerInvisible_normal_form`), i.e. exactly when the
  enriched Lagrangian is the committed two-parameter counterfamily
  member `(a, c) = (2ζ, 2ε)` (`quadPoly_cornerInvisible_iff_twoSlot`).
  The parametrization is injective
  (`chainTwoSlot_parameters_injective`): the corner-invisible quadratic
  enrichment space is exactly two-dimensional and the committed
  counterfamily is complete — there is no third corner-invisible
  quadratic direction.

* **Residual freedom (E2).**  With the derivative packet
  (`enrichedQuad_hasFDerivAt`, from the committed polynomial packet),
  the committed interior-junction stationarity condition on a
  corner-invisible quadratic enrichment holds exactly when
  `2ζ + 2ε = a*` (`quadPoly_stationary_iff`, composed with the
  committed `chainTwoSlot_modeExtremal_forced`).  The solution set is
  exactly the one-dimensional affine line `c = a* - a`
  (`chainTwoSlot_stationary_line`,
  `chainTwoSlot_stationary_existsUnique`).

* **Selection (E3).**  The named clause `VelocityOnly` (the added
  sentence is a function of the second slot alone; on corner-invisible
  quadratics this is exactly `ε = 0`, `quadPoly_velocityOnly_iff`)
  restores uniqueness: quadratic + corner-invisible + `VelocityOnly` +
  stationarity forces `(a, c) = (a*, 0)`, the committed one-parameter
  rule, with all six coefficients pinned
  (`quadPoly_selection_unique`, `quadPoly_selection_eq_committed`,
  `chainTwoSlot_velocityOnly_selection`).

* **Necessity of `VelocityOnly` (E4-i).**  The committed counterfamily
  point `(a*/2, a*/2)` satisfies every clause except `VelocityOnly`
  and carries a different momentum map and a different Lagrangian
  (`velocityOnly_clause_necessary`, re-exporting the committed
  injectivity `chainTwoSlotD2_curvature_injective`).

* **Necessity of quadraticity inside the stationarity-only test
  (E4-ii).**  For velocity-only cubic
  corner-invisible sentences `y(y-1)(μy + ν)` the committed
  stationarity equation fixes only `μ + ν = a*/2`
  (`chainCubic_stationary_iff`): the explicit cubic
  `μ = a*/2, ν = 0` is corner-invisible, velocity-only, stationary at
  the constant-one junction with a certified derivative packet, and is
  not equal to any quadratic-grammar sentence
  (`quadraticity_clause_necessary`, `chainCubic_ne_twoSlot`,
  `chainCubic_ne_quadGrammar`).  Thus quadraticity is needed for the
  stated selection theorem if only corner invisibility, velocity-only,
  and stationarity are assumed.  The stronger committed receipt does
  exclude this declared cubic family: real fixed-endpoint
  single-site minimality of the embedded constant-one history forces
  `μ = 0` and returns the committed rule exactly
  (`chainCubic_gap`, `chainCubic_realMin_forces_quadratic`).

* **Quartic continuation and global-selection boundary (E5).**  The
  velocity-only quartic increment `λ y²(y-1)²` is corner-invisible and
  leaves the committed stationarity equation unchanged.  For every
  `0 < λ < a*`, its exact single-site gap is
  `(x-1)²(a*/2 + λx²)`, so the embedded constant-one history remains a
  global real single-site minimizer.  Its momentum derivative has the
  exact positive Hessian
  `(a* - λ) + 12λ(y-1/2)²`; consequently the momentum map is strictly
  increasing and surjective, with a unique global velocity solution.
  Nevertheless the Lagrangian is neither the committed rule nor any
  member of the quadratic grammar
  (`quartic_continuation_blocks_global_selection`).  Therefore the
  cubic exclusion is not a global uniqueness theorem.

What is NOT proved.  Nothing here states or implies that the committed
source selects the enrichment; the realized-history no-go of
`RealizedHistoryLegendreNoGo` stays in force, register row PR-06 stays
a register row, and OL-D1 stays partial.  `VelocityOnly` is a named
clause of a declared axiomatization, not a source-produced fact.  In
particular, no theorem in this module proves that corner invisibility,
velocity-only, stationarity, single-site minimality, or regular
momentum inversion is a globally selecting or minimal axiom system;
the quartic continuation proves that those tested clauses still admit
nonquadratic alternatives.  All
statements are exact real algebra about the declared grammar around
the committed two-state mixing chain; no continuum, units, clock, or
amplitude is claimed.  The derivative data reference is register row
PR-45 through the committed packets.
-/

/-! ## The declared quadratic grammar and its clauses -/

/-- A sentence of the declared quadratic two-slot grammar: the general
real quadratic polynomial in the two record slots. -/
noncomputable def quadPoly (α β γ δ ε ζ : ℝ) (x y : ℝ) : ℝ :=
  α + β * x + γ * y + δ * (x * y) + ε * (x * x) + ζ * (y * y)

/-- Corner invisibility: the added sentence vanishes at all four
binary corners, so the enriched Lagrangian reads the committed corner
table exactly. -/
def CornerInvisible (q : ℝ → ℝ → ℝ) : Prop :=
  q 0 0 = 0 ∧ q 1 0 = 0 ∧ q 0 1 = 0 ∧ q 1 1 = 0

/-- The named selection clause: the added sentence is a function of
the second (velocity) slot alone. -/
def VelocityOnly (q : ℝ → ℝ → ℝ) : Prop :=
  ∀ x x' y : ℝ, q x y = q x' y

/-- The corner-invisible increment of the committed two-parameter
counterfamily: `chainTwoSlotCurvedLagrangian a c` minus the committed
bilinear extension. -/
noncomputable def twoSlotIncrement (a c : ℝ) (x y : ℝ) : ℝ :=
  (a / 2) * y * (y - 1) + (c / 2) * x * (x - 1)

/-! ## E1: completeness of the two-parameter counterfamily -/

/-- **Corner invisibility solved in coefficients.**  A quadratic
sentence vanishes at the four binary corners exactly when `α = 0`,
`δ = 0`, `β = -ε`, `γ = -ζ`.  Plain linear algebra on the four corner
evaluations. -/
theorem quadPoly_cornerInvisible_iff (α β γ δ ε ζ : ℝ) :
    CornerInvisible (quadPoly α β γ δ ε ζ)
      ↔ α = 0 ∧ δ = 0 ∧ β = -ε ∧ γ = -ζ := by
  unfold CornerInvisible quadPoly
  constructor
  · rintro ⟨h00, h10, h01, h11⟩
    norm_num at h00 h10 h01 h11
    refine ⟨h00, ?_, ?_, ?_⟩ <;> linarith
  · rintro ⟨hα, hδ, hβ, hγ⟩
    subst hα; subst hδ; subst hβ; subst hγ
    norm_num

/-- **The corner-invisible normal form.**  A quadratic sentence is
corner-invisible exactly when it is `ε x(x-1) + ζ y(y-1)`. -/
theorem quadPoly_cornerInvisible_normal_form (α β γ δ ε ζ : ℝ) :
    CornerInvisible (quadPoly α β γ δ ε ζ)
      ↔ ∀ x y : ℝ, quadPoly α β γ δ ε ζ x y
          = ε * (x * (x - 1)) + ζ * (y * (y - 1)) := by
  rw [quadPoly_cornerInvisible_iff]
  constructor
  · rintro ⟨hα, hδ, hβ, hγ⟩ x y
    subst hα; subst hδ; subst hβ; subst hγ
    unfold quadPoly
    ring
  · intro h
    have h00 := h 0 0
    have h10 := h 1 0
    have h01 := h 0 1
    have h11 := h 1 1
    unfold quadPoly at h00 h10 h01 h11
    norm_num at h00 h10 h01 h11
    refine ⟨h00, ?_, ?_, ?_⟩ <;> linarith

/-- The committed counterfamily member is the committed bilinear
extension plus the two-parameter increment. -/
theorem chainTwoSlot_eq_increment (a c x y : ℝ) :
    chainTwoSlotCurvedLagrangian a c x y
      = chainLogLagrangian x y + twoSlotIncrement a c x y := by
  unfold chainTwoSlotCurvedLagrangian twoSlotIncrement
  ring

/-- The two-parameter increment is corner-invisible at every
parameter pair: the counterfamily sits inside the classified space. -/
theorem twoSlotIncrement_cornerInvisible (a c : ℝ) :
    CornerInvisible (twoSlotIncrement a c) := by
  unfold CornerInvisible twoSlotIncrement
  norm_num

/-- **Completeness of the committed counterfamily (E1).**  A quadratic
sentence is corner-invisible exactly when the enriched Lagrangian is
the committed two-parameter member `(a, c) = (2ζ, 2ε)`.  There is no
corner-invisible quadratic direction outside the committed family. -/
theorem quadPoly_cornerInvisible_iff_twoSlot (α β γ δ ε ζ : ℝ) :
    CornerInvisible (quadPoly α β γ δ ε ζ)
      ↔ ∀ x y : ℝ,
          chainLogLagrangian x y + quadPoly α β γ δ ε ζ x y
            = chainTwoSlotCurvedLagrangian (2 * ζ) (2 * ε) x y := by
  constructor
  · intro h x y
    have hn := (quadPoly_cornerInvisible_normal_form α β γ δ ε ζ).mp h x y
    rw [hn, chainTwoSlot_eq_increment]
    unfold twoSlotIncrement
    ring
  · intro h
    refine (quadPoly_cornerInvisible_normal_form α β γ δ ε ζ).mpr ?_
    intro x y
    have hx := h x y
    rw [chainTwoSlot_eq_increment] at hx
    rw [add_left_cancel hx]
    unfold twoSlotIncrement
    ring

/-- **The parametrization is injective (E1 dimension count).**
Distinct parameter pairs give distinct enrichments: the
corner-invisible quadratic enrichment space is exactly
two-dimensional. -/
theorem chainTwoSlot_parameters_injective :
    Function.Injective
      (fun p : ℝ × ℝ => chainTwoSlotCurvedLagrangian p.1 p.2) := by
  rintro ⟨a, c⟩ ⟨a', c'⟩ h
  have h02 := congrFun (congrFun h 0) 2
  have h20 := congrFun (congrFun h 2) 0
  simp only [chainTwoSlotCurvedLagrangian] at h02 h20
  norm_num at h02 h20
  have ha : a = a' := by linarith
  have hc : c = c' := by linarith
  rw [ha, hc]

/-! ## E2: residual freedom under the committed stationarity -/

/-- The first-slot derivative of a quadratic sentence. -/
noncomputable def quadD1 (β δ ε : ℝ) (x y : ℝ) : ℝ :=
  β + δ * y + 2 * ε * x

/-- The second-slot derivative of a quadratic sentence. -/
noncomputable def quadD2 (γ δ ζ : ℝ) (x y : ℝ) : ℝ :=
  γ + δ * x + 2 * ζ * y

/-- The uncurried derivative packet of a quadratic-grammar enrichment,
from the committed polynomial packet: the slot derivatives of the
enriched Lagrangian are the committed slot derivatives of the bilinear
extension plus the quadratic sentence's slot derivatives. -/
theorem enrichedQuad_hasFDerivAt (α β γ δ ε ζ : ℝ) (p : ℝ × ℝ) :
    HasFDerivAt
      (Function.uncurry
        (fun x y => chainLogLagrangian x y + quadPoly α β γ δ ε ζ x y))
      ((chainCurvedD1 p.1 p.2 + quadD1 β δ ε p.1 p.2)
          • ContinuousLinearMap.fst ℝ ℝ ℝ
        + (chainFiberSlope p.1 + quadD2 γ δ ζ p.1 p.2)
          • ContinuousLinearMap.snd ℝ ℝ ℝ) p := by
  have hfun : Function.uncurry
      (fun x y => chainLogLagrangian x y + quadPoly α β γ δ ε ζ x y)
      = fun q : ℝ × ℝ =>
          (-Real.log (chainWeight 0 0) + α)
            + (chainBaseSlope + β) * q.1
            + (Real.log (chainWeight 0 0) - Real.log (chainWeight 0 1)
                + γ) * q.2
            + (chainSlopeRate + δ) * (q.1 * q.2)
            + ε * (q.1 * q.1)
            + ζ * (q.2 * q.2) := by
    funext q
    simp only [Function.uncurry, quadPoly, chainLogLagrangian,
      chainBaseSlope, chainSlopeRate]
    ring
  have h := polyTwoPoint_hasFDerivAt (-Real.log (chainWeight 0 0) + α)
    (chainBaseSlope + β)
    (Real.log (chainWeight 0 0) - Real.log (chainWeight 0 1) + γ)
    (chainSlopeRate + δ) ε ζ p
  have h1 : chainCurvedD1 p.1 p.2 + quadD1 β δ ε p.1 p.2
      = chainBaseSlope + β + (chainSlopeRate + δ) * p.2
        + 2 * ε * p.1 := by
    simp only [chainCurvedD1, quadD1]
    ring
  have h2 : chainFiberSlope p.1 + quadD2 γ δ ζ p.1 p.2
      = Real.log (chainWeight 0 0) - Real.log (chainWeight 0 1) + γ
        + (chainSlopeRate + δ) * p.1 + 2 * ζ * p.2 := by
    simp only [quadD2, chainFiberSlope, chainSlopeRate]
    ring
  rw [hfun, h1, h2]
  exact h

/-- **Residual freedom (E2).**  On a corner-invisible quadratic
enrichment, the committed interior-junction stationarity condition at
the constant-one junction holds exactly when `2ζ + 2ε = a*`.  Composed
with the committed `chainTwoSlot_modeExtremal_forced`; nothing is
reproved. -/
theorem quadPoly_stationary_iff (α β γ δ ε ζ : ℝ)
    (hcorner : CornerInvisible (quadPoly α β γ δ ε ζ)) :
    (chainFiberSlope 1 + quadD2 γ δ ζ 1 1)
        + (chainCurvedD1 1 1 + quadD1 β δ ε 1 1) = 0
      ↔ 2 * ζ + 2 * ε = modeExtremalCurvature := by
  obtain ⟨hα, hδ, hβ, hγ⟩ :=
    (quadPoly_cornerInvisible_iff α β γ δ ε ζ).mp hcorner
  have hEq : (chainFiberSlope 1 + quadD2 γ δ ζ 1 1)
        + (chainCurvedD1 1 1 + quadD1 β δ ε 1 1)
      = chainTwoSlotD2 (2 * ζ) 1 1 + chainTwoSlotD1 (2 * ε) 1 1 := by
    subst hδ; subst hβ; subst hγ
    simp only [quadD1, quadD2, chainTwoSlotD1, chainTwoSlotD2,
      chainCurvedD2]
    ring
  rw [hEq, chainTwoSlot_modeExtremal_forced]

/-- The stationary solution set in family coordinates is the affine
line `c = a* - a`. -/
theorem chainTwoSlot_stationary_line (a c : ℝ) :
    chainTwoSlotD2 a 1 1 + chainTwoSlotD1 c 1 1 = 0
      ↔ c = modeExtremalCurvature - a := by
  rw [chainTwoSlot_modeExtremal_forced]
  constructor <;> intro h <;> linarith

/-- For every velocity curvature `a` there is exactly one compensating
first-slot curvature meeting the committed stationarity condition: the
solution set is a one-dimensional affine line, not a point. -/
theorem chainTwoSlot_stationary_existsUnique (a : ℝ) :
    ∃! c : ℝ, chainTwoSlotD2 a 1 1 + chainTwoSlotD1 c 1 1 = 0 := by
  refine ⟨modeExtremalCurvature - a,
    (chainTwoSlot_stationary_line a _).mpr rfl, fun c hc => ?_⟩
  exact (chainTwoSlot_stationary_line a c).mp hc

/-! ## E3: the velocity-only clause restores uniqueness -/

/-- On the committed counterfamily increment, the velocity-only clause
is exactly the vanishing of the first-slot parameter. -/
theorem twoSlotIncrement_velocityOnly_iff (a c : ℝ) :
    VelocityOnly (twoSlotIncrement a c) ↔ c = 0 := by
  constructor
  · intro h
    have h20 := h 2 0 0
    unfold twoSlotIncrement at h20
    norm_num at h20
    linarith
  · intro hc
    subst hc
    intro x x' y
    unfold twoSlotIncrement
    ring

/-- On a corner-invisible quadratic sentence, the velocity-only clause
is exactly `ε = 0`: the `x(x-1)` coefficient vanishes. -/
theorem quadPoly_velocityOnly_iff (α β γ δ ε ζ : ℝ)
    (hcorner : CornerInvisible (quadPoly α β γ δ ε ζ)) :
    VelocityOnly (quadPoly α β γ δ ε ζ) ↔ ε = 0 := by
  have hn := (quadPoly_cornerInvisible_normal_form α β γ δ ε ζ).mp hcorner
  constructor
  · intro h
    have h20 := h 2 0 0
    rw [hn 2 0, hn 0 0] at h20
    norm_num at h20
    linarith
  · intro hε x x' y
    rw [hn x y, hn x' y, hε]
    ring

/-- The `c = 0` axis of the counterfamily is the committed
one-parameter family. -/
theorem chainTwoSlot_zero_eq_curved (a : ℝ) :
    chainTwoSlotCurvedLagrangian a 0 = chainCurvedLagrangian a := by
  funext x y
  unfold chainTwoSlotCurvedLagrangian chainCurvedLagrangian
  ring

/-- **Selection in family coordinates (E3).**  Committed stationarity
plus the velocity-only clause forces `(a, c) = (a*, 0)`: the committed
one-parameter rule, uniquely. -/
theorem chainTwoSlot_velocityOnly_selection (a c : ℝ)
    (hstat : chainTwoSlotD2 a 1 1 + chainTwoSlotD1 c 1 1 = 0)
    (hvel : VelocityOnly (twoSlotIncrement a c)) :
    a = modeExtremalCurvature ∧ c = 0
      ∧ chainTwoSlotCurvedLagrangian a c
          = chainCurvedLagrangian modeExtremalCurvature := by
  have hc : c = 0 := (twoSlotIncrement_velocityOnly_iff a c).mp hvel
  have hsum := (chainTwoSlot_modeExtremal_forced a c).mp hstat
  have ha : a = modeExtremalCurvature := by
    rw [hc] at hsum
    linarith
  refine ⟨ha, hc, ?_⟩
  rw [ha, hc, chainTwoSlot_zero_eq_curved]

/-- **The two-clause axiomatization pins all six coefficients (E3).**
A quadratic sentence satisfies corner invisibility, the velocity-only
clause, and the committed stationarity condition exactly when it is
the committed one-parameter rule's increment: `α = β = δ = ε = 0`,
`γ = -a*/2`, `ζ = a*/2`. -/
theorem quadPoly_selection_unique (α β γ δ ε ζ : ℝ) :
    (CornerInvisible (quadPoly α β γ δ ε ζ)
        ∧ VelocityOnly (quadPoly α β γ δ ε ζ)
        ∧ (chainFiberSlope 1 + quadD2 γ δ ζ 1 1)
            + (chainCurvedD1 1 1 + quadD1 β δ ε 1 1) = 0)
      ↔ (α = 0 ∧ β = 0 ∧ δ = 0 ∧ ε = 0
          ∧ γ = -(modeExtremalCurvature / 2)
          ∧ ζ = modeExtremalCurvature / 2) := by
  constructor
  · rintro ⟨hcorner, hvel, hstat⟩
    obtain ⟨hα, hδ, hβ, hγ⟩ :=
      (quadPoly_cornerInvisible_iff α β γ δ ε ζ).mp hcorner
    have hε : ε = 0 :=
      (quadPoly_velocityOnly_iff α β γ δ ε ζ hcorner).mp hvel
    have hsum :=
      (quadPoly_stationary_iff α β γ δ ε ζ hcorner).mp hstat
    refine ⟨hα, by linarith, hδ, hε, by linarith, by linarith⟩
  · rintro ⟨hα, hβ, hδ, hε, hγ, hζ⟩
    have hcorner : CornerInvisible (quadPoly α β γ δ ε ζ) :=
      (quadPoly_cornerInvisible_iff α β γ δ ε ζ).mpr
        ⟨hα, hδ, by rw [hβ, hε]; norm_num, by rw [hγ, hζ]⟩
    refine ⟨hcorner,
      (quadPoly_velocityOnly_iff α β γ δ ε ζ hcorner).mpr hε, ?_⟩
    refine (quadPoly_stationary_iff α β γ δ ε ζ hcorner).mpr ?_
    rw [hζ, hε]
    ring

/-- **The selected enrichment is the committed rule (E3).**  Under the
two clauses and stationarity, the enriched Lagrangian equals the
committed mode-extremal representative pointwise. -/
theorem quadPoly_selection_eq_committed (α β γ δ ε ζ : ℝ)
    (hcorner : CornerInvisible (quadPoly α β γ δ ε ζ))
    (hvel : VelocityOnly (quadPoly α β γ δ ε ζ))
    (hstat : (chainFiberSlope 1 + quadD2 γ δ ζ 1 1)
        + (chainCurvedD1 1 1 + quadD1 β δ ε 1 1) = 0) :
    ∀ x y : ℝ, chainLogLagrangian x y + quadPoly α β γ δ ε ζ x y
      = chainCurvedLagrangian modeExtremalCurvature x y := by
  intro x y
  have hmem :=
    (quadPoly_cornerInvisible_iff_twoSlot α β γ δ ε ζ).mp hcorner x y
  have hε : ε = 0 :=
    (quadPoly_velocityOnly_iff α β γ δ ε ζ hcorner).mp hvel
  have hsum := (quadPoly_stationary_iff α β γ δ ε ζ hcorner).mp hstat
  have hζ : 2 * ζ = modeExtremalCurvature := by linarith
  rw [hmem, hε, hζ, mul_zero, chainTwoSlot_zero_eq_curved]

/-! ## E4-i: the velocity-only clause is necessary -/

/-- **Necessity of the velocity-only clause (E4-i).**  The committed
counterfamily point `(a*/2, a*/2)` is corner-invisible and satisfies
the committed stationarity condition, is not velocity-only, and has a
momentum map and a Lagrangian different from the committed rule's.  So
dropping `VelocityOnly` breaks uniqueness: the clause is necessary. -/
theorem velocityOnly_clause_necessary :
    ∃ a c : ℝ,
      0 < a
        ∧ chainTwoSlotD2 a 1 1 + chainTwoSlotD1 c 1 1 = 0
        ∧ CornerInvisible (twoSlotIncrement a c)
        ∧ ¬ VelocityOnly (twoSlotIncrement a c)
        ∧ chainTwoSlotD2 a ≠ chainTwoSlotD2 modeExtremalCurvature
        ∧ chainTwoSlotCurvedLagrangian a c
            ≠ chainCurvedLagrangian modeExtremalCurvature := by
  have hpos := modeExtremalCurvature_pos
  refine ⟨modeExtremalCurvature / 2, modeExtremalCurvature / 2,
    half_pos hpos, (chainTwoSlot_modeExtremal_forced _ _).mpr (by ring),
    twoSlotIncrement_cornerInvisible _ _, ?_, ?_, ?_⟩
  · intro hvel
    have hc := (twoSlotIncrement_velocityOnly_iff _ _).mp hvel
    linarith
  · intro h
    have := chainTwoSlotD2_curvature_injective h
    linarith
  · intro h
    have hval := congrFun (congrFun h (1 / 2 : ℝ)) 0
    simp only [chainTwoSlotCurvedLagrangian, chainCurvedLagrangian]
      at hval
    ring_nf at hval
    linarith

/-! ## E4-ii: quadraticity is necessary — the cubic counterexample -/

/-- A velocity-only cubic corner-invisible sentence:
`y (y-1) (μ y + ν)`. -/
noncomputable def cubicIncrement (μ ν : ℝ) (y : ℝ) : ℝ :=
  y * (y - 1) * (μ * y + ν)

/-- The derivative of the cubic sentence in the velocity slot. -/
noncomputable def cubicIncrementDeriv (μ ν : ℝ) (y : ℝ) : ℝ :=
  3 * μ * y ^ 2 + 2 * (ν - μ) * y - ν

/-- The cubic-enriched Lagrangian: the committed bilinear extension
plus a velocity-only cubic corner-invisible sentence. -/
noncomputable def chainCubicLagrangian (μ ν : ℝ) (x y : ℝ) : ℝ :=
  chainLogLagrangian x y + cubicIncrement μ ν y

/-- The second-slot derivative of the cubic-enriched Lagrangian. -/
noncomputable def chainCubicD2 (μ ν : ℝ) (x y : ℝ) : ℝ :=
  chainFiberSlope x + cubicIncrementDeriv μ ν y

/-- Every cubic sentence is corner-invisible: `y (y-1)` divides it. -/
theorem cubicIncrement_cornerInvisible (μ ν : ℝ) :
    CornerInvisible (fun _ y => cubicIncrement μ ν y) := by
  unfold CornerInvisible cubicIncrement
  norm_num

/-- Every cubic sentence is velocity-only by construction. -/
theorem cubicIncrement_velocityOnly (μ ν : ℝ) :
    VelocityOnly (fun _ y => cubicIncrement μ ν y) :=
  fun _ _ _ => rfl

/-- The cubic-enriched Lagrangian reads the committed corner table
exactly, at every `(μ, ν)`. -/
theorem chainCubic_corner (μ ν : ℝ) (i j : Fin 2) :
    chainCubicLagrangian μ ν ((i : ℕ) : ℝ) ((j : ℕ) : ℝ)
      = -Real.log (chainWeight i j) := by
  rw [chainCubicLagrangian, chainLogLagrangian_corner]
  fin_cases j <;> norm_num [cubicIncrement]

/-- The velocity-slot derivative of the cubic sentence. -/
theorem cubicIncrement_hasDerivAt (μ ν y : ℝ) :
    HasDerivAt (cubicIncrement μ ν) (cubicIncrementDeriv μ ν y) y := by
  unfold cubicIncrement cubicIncrementDeriv
  convert (((hasDerivAt_id y).mul ((hasDerivAt_id y).sub_const 1)).mul
    (((hasDerivAt_id y).const_mul μ).add_const ν)) using 1
  simp only [id_eq, Pi.mul_apply]
  ring

/-- The uncurried derivative packet of the cubic-enriched Lagrangian:
the committed bilinear packet plus the cubic velocity-slot
derivative. -/
theorem chainCubic_hasFDerivAt (μ ν : ℝ) (p : ℝ × ℝ) :
    HasFDerivAt (Function.uncurry (chainCubicLagrangian μ ν))
      (chainCurvedD1 p.1 p.2 • ContinuousLinearMap.fst ℝ ℝ ℝ
        + chainCubicD2 μ ν p.1 p.2 • ContinuousLinearMap.snd ℝ ℝ ℝ)
      p := by
  have hL0 : HasFDerivAt (Function.uncurry chainLogLagrangian)
      (chainCurvedD1 p.1 p.2 • ContinuousLinearMap.fst ℝ ℝ ℝ
        + chainCurvedD2 0 p.1 p.2 • ContinuousLinearMap.snd ℝ ℝ ℝ)
      p := by
    have h := chainCurvedLagrangian_hasFDerivAt 0 p
    rwa [chainCurvedLagrangian_zero] at h
  have hg : HasFDerivAt (fun q : ℝ × ℝ => cubicIncrement μ ν q.2)
      (cubicIncrementDeriv μ ν p.2 • ContinuousLinearMap.snd ℝ ℝ ℝ)
      p := by
    have h2 := ((cubicIncrement_hasDerivAt μ ν p.2).hasFDerivAt).comp p
      hasFDerivAt_snd
    have hclm : (ContinuousLinearMap.toSpanSingleton ℝ
          (cubicIncrementDeriv μ ν p.2)).comp
          (ContinuousLinearMap.snd ℝ ℝ ℝ)
        = cubicIncrementDeriv μ ν p.2
            • ContinuousLinearMap.snd ℝ ℝ ℝ := by
      refine ContinuousLinearMap.ext fun v => ?_
      simp only [ContinuousLinearMap.coe_comp', Function.comp_apply,
        ContinuousLinearMap.toSpanSingleton_apply,
        ContinuousLinearMap.coe_snd', ContinuousLinearMap.smul_apply,
        smul_eq_mul]
      ring
    rwa [hclm] at h2
  have hsum := hL0.add hg
  have hfun : Function.uncurry (chainCubicLagrangian μ ν)
      = fun q : ℝ × ℝ =>
          Function.uncurry chainLogLagrangian q
            + cubicIncrement μ ν q.2 := by
    funext q
    simp only [Function.uncurry, chainCubicLagrangian]
  have hD : chainCurvedD1 p.1 p.2 • ContinuousLinearMap.fst ℝ ℝ ℝ
        + chainCubicD2 μ ν p.1 p.2 • ContinuousLinearMap.snd ℝ ℝ ℝ
      = (chainCurvedD1 p.1 p.2 • ContinuousLinearMap.fst ℝ ℝ ℝ
          + chainCurvedD2 0 p.1 p.2 • ContinuousLinearMap.snd ℝ ℝ ℝ)
        + cubicIncrementDeriv μ ν p.2
            • ContinuousLinearMap.snd ℝ ℝ ℝ := by
    refine ContinuousLinearMap.ext fun v => ?_
    simp only [ContinuousLinearMap.add_apply,
      ContinuousLinearMap.smul_apply, ContinuousLinearMap.coe_fst',
      ContinuousLinearMap.coe_snd', smul_eq_mul, chainCubicD2,
      chainCurvedD2]
    ring
  rw [hfun, hD]
  exact hsum

/-- **The committed stationarity equation on the cubics.**  At the
constant-one junction it fixes only the sum `μ + ν = a*/2`: a
one-parameter family of cubics survives the committed condition. -/
theorem chainCubic_stationary_iff (μ ν : ℝ) :
    chainCubicD2 μ ν 1 1 + chainCurvedD1 1 1 = 0
      ↔ μ + ν = modeExtremalCurvature / 2 := by
  have he := modeExtremalCurvature_eq
  simp only [chainCubicD2, cubicIncrementDeriv, chainCurvedD1,
    chainFiberSlope, chainBaseSlope, chainSlopeRate]
  constructor <;> intro h <;> nlinarith [he, h]

/-- A cubic sentence with `μ ≠ 0` is not any member of the committed
two-parameter counterfamily. -/
theorem chainCubic_ne_twoSlot (μ ν : ℝ) (hμ : μ ≠ 0) (a c : ℝ) :
    chainCubicLagrangian μ ν ≠ chainTwoSlotCurvedLagrangian a c := by
  intro h
  have h1 := congrFun (congrFun h 0) (-1)
  have h2 := congrFun (congrFun h 0) 2
  simp only [chainCubicLagrangian, chainTwoSlotCurvedLagrangian,
    cubicIncrement] at h1 h2
  ring_nf at h1 h2
  apply hμ
  linarith

/-- A cubic sentence with `μ ≠ 0` is not any sentence of the quadratic
grammar at all, corner-invisible or not. -/
theorem chainCubic_ne_quadGrammar (μ ν : ℝ) (hμ : μ ≠ 0)
    (α β γ δ ε ζ : ℝ) :
    chainCubicLagrangian μ ν
      ≠ fun x y => chainLogLagrangian x y + quadPoly α β γ δ ε ζ x y := by
  intro h
  have h0 := congrFun (congrFun h 0) 0
  have h1 := congrFun (congrFun h 0) 1
  have hm := congrFun (congrFun h 0) (-1)
  have h2 := congrFun (congrFun h 0) 2
  simp only [chainCubicLagrangian, cubicIncrement, quadPoly]
    at h0 h1 hm h2
  ring_nf at h0 h1 hm h2
  apply hμ
  linarith

/-- **Necessity of quadraticity (E4-ii): the explicit counterexample.**
The cubic `μ = a*/2, ν = 0` is corner-invisible, velocity-only, and
satisfies the committed interior-junction stationarity equation with a
certified derivative packet, yet equals no sentence of the quadratic
grammar.  Corner invisibility + velocity-only + stationarity WITHOUT
quadraticity does not return the committed rule: the quadraticity
clause is load-bearing. -/
theorem quadraticity_clause_necessary :
    ∃ μ ν : ℝ,
      μ ≠ 0
        ∧ CornerInvisible (fun _ y => cubicIncrement μ ν y)
        ∧ VelocityOnly (fun _ y => cubicIncrement μ ν y)
        ∧ chainCubicD2 μ ν 1 1 + chainCurvedD1 1 1 = 0
        ∧ (∀ i j : Fin 2,
            chainCubicLagrangian μ ν ((i : ℕ) : ℝ) ((j : ℕ) : ℝ)
              = -Real.log (chainWeight i j))
        ∧ (∀ p : ℝ × ℝ,
            HasFDerivAt (Function.uncurry (chainCubicLagrangian μ ν))
              (chainCurvedD1 p.1 p.2 • ContinuousLinearMap.fst ℝ ℝ ℝ
                + chainCubicD2 μ ν p.1 p.2
                    • ContinuousLinearMap.snd ℝ ℝ ℝ) p)
        ∧ (∀ a c : ℝ,
            chainCubicLagrangian μ ν ≠ chainTwoSlotCurvedLagrangian a c)
        ∧ ∀ α β γ δ ε ζ : ℝ,
            chainCubicLagrangian μ ν
              ≠ fun x y =>
                  chainLogLagrangian x y + quadPoly α β γ δ ε ζ x y := by
  have hne : modeExtremalCurvature / 2 ≠ 0 :=
    ne_of_gt (half_pos modeExtremalCurvature_pos)
  exact ⟨modeExtremalCurvature / 2, 0, hne,
    cubicIncrement_cornerInvisible _ _,
    cubicIncrement_velocityOnly _ _,
    (chainCubic_stationary_iff _ _).mpr (by ring),
    chainCubic_corner _ _,
    chainCubic_hasFDerivAt _ _,
    chainCubic_ne_twoSlot _ _ hne,
    fun α β γ δ ε ζ => chainCubic_ne_quadGrammar _ _ hne α β γ δ ε ζ⟩

/-- **The exact single-site variation gap of a stationary cubic.**
Under `μ + ν = a*/2`, the gap at the embedded constant-one junction is
exactly `(x-1)² (μ x + a*/2)`: nonnegative for all real `x` only when
`μ = 0`. -/
theorem chainCubic_gap (μ ν x : ℝ)
    (hsum : μ + ν = modeExtremalCurvature / 2) :
    (chainCubicLagrangian μ ν 1 x - chainCubicLagrangian μ ν 1 1)
        + (chainCubicLagrangian μ ν x 1 - chainCubicLagrangian μ ν 1 1)
      = (x - 1) ^ 2 * (μ * x + modeExtremalCurvature / 2) := by
  have he := modeExtremalCurvature_eq
  have hν : ν = modeExtremalCurvature / 2 - μ := by linarith
  subst hν
  simp only [chainCubicLagrangian, cubicIncrement, chainLogLagrangian]
  rw [he]
  ring

/-- **The committed minimality receipt excludes the cubics
(E4-ii, strengthening).**  A stationary cubic whose embedded
constant-one history is a real fixed-endpoint single-site minimizer at
one interior junction has `μ = 0`, `ν = a*/2`, and equals the
committed one-parameter rule exactly.  The stationarity equation alone
does not exclude the cubics (`quadraticity_clause_necessary`); the
minimality receipt does. -/
theorem chainCubic_realMin_forces_quadratic (μ ν : ℝ)
    (hsum : μ + ν = modeExtremalCurvature / 2)
    (M : ℕ) {k m : Fin M} (hkm : k.succ = m.castSucc)
    (hmin : ∀ x : ℝ,
      localAction (chainCubicLagrangian μ ν)
          (chainEmb (constOneHistory M))
        ≤ localAction (chainCubicLagrangian μ ν)
            (Function.update (chainEmb (constOneHistory M)) k.succ x)) :
    μ = 0 ∧ ν = modeExtremalCurvature / 2
      ∧ chainCubicLagrangian μ ν
          = chainCurvedLagrangian modeExtremalCurvature := by
  have hgap : ∀ x : ℝ,
      0 ≤ (x - 1) ^ 2 * (μ * x + modeExtremalCurvature / 2) := by
    intro x
    have hdiff := localAction_update_diff (chainCubicLagrangian μ ν)
      (chainEmb (constOneHistory M)) hkm x
    simp only [chainEmb_constOne] at hdiff
    have hid := chainCubic_gap μ ν x hsum
    have h := hmin x
    linarith
  have hμ : μ = 0 := by
    by_contra hμ
    have key : ∀ t : ℝ, 0 < t →
        (-(modeExtremalCurvature / 2) - t) / μ = 1 := by
      intro t ht
      set x : ℝ := (-(modeExtremalCurvature / 2) - t) / μ with hxdef
      have hx := hgap x
      have hv : μ * x + modeExtremalCurvature / 2 = -t := by
        rw [hxdef]
        field_simp
        ring
      rw [hv] at hx
      have hsq : (x - 1) ^ 2 ≤ 0 := by
        nlinarith [hx, ht, sq_nonneg (x - 1)]
      have h0 : (x - 1) ^ 2 = 0 := le_antisymm hsq (sq_nonneg _)
      have hx1 : x - 1 = 0 := by
        exact sq_eq_zero_iff.mp h0
      linarith
    have h1 := key 1 one_pos
    have h2 := key 2 two_pos
    rw [div_eq_one_iff_eq hμ] at h1 h2
    linarith
  have hν : ν = modeExtremalCurvature / 2 := by
    rw [hμ] at hsum
    linarith
  refine ⟨hμ, hν, ?_⟩
  subst hμ; subst hν
  funext x y
  simp only [chainCubicLagrangian, cubicIncrement, chainCurvedLagrangian]
  ring

/-! ## E5: a quartic continuation beyond the cubic boundary -/

/-- The quartic velocity-only increment used to test whether the
single-site minimality receipt globally selects the committed rule. -/
noncomputable def quarticIncrement (lam y : ℝ) : ℝ :=
  lam * y ^ 2 * (y - 1) ^ 2

/-- The first velocity derivative of `quarticIncrement`. -/
noncomputable def quarticIncrementDeriv (lam y : ℝ) : ℝ :=
  lam * (4 * y ^ 3 - 6 * y ^ 2 + 2 * y)

/-- The second velocity derivative of `quarticIncrement`. -/
noncomputable def quarticIncrementSecondDeriv (lam y : ℝ) : ℝ :=
  lam * (12 * y ^ 2 - 12 * y + 2)

/-- The committed curved rule plus the quartic velocity increment. -/
noncomputable def chainQuarticLagrangian (lam : ℝ) (x y : ℝ) : ℝ :=
  chainCurvedLagrangian modeExtremalCurvature x y
    + quarticIncrement lam y

/-- The exact second-slot derivative of `chainQuarticLagrangian`. -/
noncomputable def chainQuarticD2 (lam : ℝ) (x y : ℝ) : ℝ :=
  chainCurvedD2 modeExtremalCurvature x y
    + quarticIncrementDeriv lam y

/-- The quartic increment vanishes at every binary corner. -/
theorem quarticIncrement_cornerInvisible (lam : ℝ) :
    CornerInvisible (fun _ y => quarticIncrement lam y) := by
  unfold CornerInvisible quarticIncrement
  norm_num

/-- The quartic increment is velocity-only by construction. -/
theorem quarticIncrement_velocityOnly (lam : ℝ) :
    VelocityOnly (fun _ y => quarticIncrement lam y) :=
  fun _ _ _ => rfl

/-- The quartic continuation reads the committed corner table exactly. -/
theorem chainQuartic_corner (lam : ℝ) (i j : Fin 2) :
    chainQuarticLagrangian lam ((i : ℕ) : ℝ) ((j : ℕ) : ℝ)
      = -Real.log (chainWeight i j) := by
  rw [chainQuarticLagrangian, chainCurvedLagrangian_corner]
  fin_cases j <;> norm_num [quarticIncrement]

/-- Exact first derivative of the quartic increment. -/
theorem quarticIncrement_hasDerivAt (lam y : ℝ) :
    HasDerivAt (quarticIncrement lam) (quarticIncrementDeriv lam y) y := by
  have h := ((((hasDerivAt_id y).pow 4).const_mul lam).sub
    (((hasDerivAt_id y).pow 3).const_mul (2 * lam))).add
      (((hasDerivAt_id y).pow 2).const_mul lam)
  have hfun : quarticIncrement lam
      = fun t : ℝ => lam * t ^ 4 - (2 * lam) * t ^ 3 + lam * t ^ 2 := by
    funext t
    unfold quarticIncrement
    ring
  have hval : quarticIncrementDeriv lam y
      = 4 * lam * y ^ 3 - 6 * lam * y ^ 2 + 2 * lam * y := by
    unfold quarticIncrementDeriv
    ring
  rw [hfun, hval]
  convert h using 1
  all_goals (try simp only [id_eq])
  all_goals ring

/-- Exact derivative of the first-derivative polynomial. -/
theorem quarticIncrementDeriv_hasDerivAt (lam y : ℝ) :
    HasDerivAt (quarticIncrementDeriv lam)
      (quarticIncrementSecondDeriv lam y) y := by
  unfold quarticIncrementDeriv quarticIncrementSecondDeriv
  convert ((((((hasDerivAt_id y).pow 3).const_mul 4).sub
    (((hasDerivAt_id y).pow 2).const_mul 6)).add
      ((hasDerivAt_id y).const_mul 2)).const_mul lam) using 1
  all_goals simp only [id_eq]
  all_goals ring

/-- The exact momentum relation of the quartic continuation. -/
theorem chainQuartic_momentum (lam x y : ℝ) :
    MomentumRelation (chainQuarticLagrangian lam) x y
      (chainQuarticD2 lam x y) := by
  have hbase :=
    chainCurvedLagrangian_momentum modeExtremalCurvature x y
  have hquartic := quarticIncrement_hasDerivAt lam y
  have h := hbase.add hquartic
  simpa only [chainQuarticLagrangian, chainQuarticD2,
    MomentumRelation] using h

/-- The uncurried derivative packet of the quartic continuation. -/
theorem chainQuartic_hasFDerivAt (lam : ℝ) (p : ℝ × ℝ) :
    HasFDerivAt (Function.uncurry (chainQuarticLagrangian lam))
      (chainCurvedD1 p.1 p.2 • ContinuousLinearMap.fst ℝ ℝ ℝ
        + chainQuarticD2 lam p.1 p.2
            • ContinuousLinearMap.snd ℝ ℝ ℝ) p := by
  have hL := chainCurvedLagrangian_hasFDerivAt
    modeExtremalCurvature p
  have hq : HasFDerivAt (fun z : ℝ × ℝ => quarticIncrement lam z.2)
      (quarticIncrementDeriv lam p.2
        • ContinuousLinearMap.snd ℝ ℝ ℝ) p := by
    have h2 := ((quarticIncrement_hasDerivAt lam p.2).hasFDerivAt).comp p
      hasFDerivAt_snd
    have hclm : (ContinuousLinearMap.toSpanSingleton ℝ
          (quarticIncrementDeriv lam p.2)).comp
          (ContinuousLinearMap.snd ℝ ℝ ℝ)
        = quarticIncrementDeriv lam p.2
            • ContinuousLinearMap.snd ℝ ℝ ℝ := by
      refine ContinuousLinearMap.ext fun v => ?_
      simp only [ContinuousLinearMap.coe_comp', Function.comp_apply,
        ContinuousLinearMap.toSpanSingleton_apply,
        ContinuousLinearMap.coe_snd', ContinuousLinearMap.smul_apply,
        smul_eq_mul]
      ring
    rwa [hclm] at h2
  have hsum := hL.add hq
  have hfun : Function.uncurry (chainQuarticLagrangian lam)
      = fun z : ℝ × ℝ =>
          Function.uncurry
            (chainCurvedLagrangian modeExtremalCurvature) z
              + quarticIncrement lam z.2 := by
    funext z
    simp only [Function.uncurry, chainQuarticLagrangian]
  have hD : chainCurvedD1 p.1 p.2 • ContinuousLinearMap.fst ℝ ℝ ℝ
        + chainQuarticD2 lam p.1 p.2
            • ContinuousLinearMap.snd ℝ ℝ ℝ
      = (chainCurvedD1 p.1 p.2 • ContinuousLinearMap.fst ℝ ℝ ℝ
          + chainCurvedD2 modeExtremalCurvature p.1 p.2
              • ContinuousLinearMap.snd ℝ ℝ ℝ)
        + quarticIncrementDeriv lam p.2
            • ContinuousLinearMap.snd ℝ ℝ ℝ := by
    refine ContinuousLinearMap.ext fun v => ?_
    simp only [ContinuousLinearMap.add_apply,
      ContinuousLinearMap.smul_apply, ContinuousLinearMap.coe_fst',
      ContinuousLinearMap.coe_snd', smul_eq_mul, chainQuarticD2]
    ring
  rw [hfun, hD]
  exact hsum

/-- The quartic increment and its first derivative vanish at the
constant-one junction, so the committed stationarity equation is
unchanged for every `lam`. -/
theorem chainQuartic_stationary (lam : ℝ) :
    chainQuarticD2 lam 1 1 + chainCurvedD1 1 1 = 0 := by
  have h := (modeExtremal_forced modeExtremalCurvature).mpr rfl
  unfold chainQuarticD2 quarticIncrementDeriv
  norm_num
  exact h

/-- Exact fixed-endpoint single-site gap of the quartic continuation. -/
theorem chainQuartic_gap (lam x : ℝ) :
    (chainQuarticLagrangian lam 1 x
        - chainQuarticLagrangian lam 1 1)
      + (chainQuarticLagrangian lam x 1
        - chainQuarticLagrangian lam 1 1)
      = (x - 1) ^ 2
          * (modeExtremalCurvature / 2 + lam * x ^ 2) := by
  have h := chainCurved_modeExtremal_gap x
  unfold chainQuarticLagrangian quarticIncrement
  norm_num
  linarith

/-- For nonnegative `lam`, the embedded constant-one history remains
a global real fixed-endpoint single-site minimizer at every admissible
interior junction. -/
theorem chainQuartic_constOne_realMin (lam : ℝ) (hlam : 0 ≤ lam)
    (M : ℕ) {k m : Fin M} (hkm : k.succ = m.castSucc) (x : ℝ) :
    localAction (chainQuarticLagrangian lam)
        (chainEmb (constOneHistory M))
      ≤ localAction (chainQuarticLagrangian lam)
          (Function.update (chainEmb (constOneHistory M)) k.succ x) := by
  have hdiff := localAction_update_diff (chainQuarticLagrangian lam)
    (chainEmb (constOneHistory M)) hkm x
  simp only [chainEmb_constOne] at hdiff
  have hgap := chainQuartic_gap lam x
  have hcoef :
      0 ≤ modeExtremalCurvature / 2 + lam * x ^ 2 :=
    add_nonneg (le_of_lt (half_pos modeExtremalCurvature_pos))
      (mul_nonneg hlam (sq_nonneg x))
  have hnonneg :
      0 ≤ (x - 1) ^ 2
        * (modeExtremalCurvature / 2 + lam * x ^ 2) :=
    mul_nonneg (sq_nonneg _) hcoef
  linarith

/-- Exact square-completed form of the velocity Hessian. -/
theorem chainQuartic_hessian_identity (lam y : ℝ) :
    modeExtremalCurvature + quarticIncrementSecondDeriv lam y
      = (modeExtremalCurvature - lam)
          + 12 * lam * (y - (1 / 2 : ℝ)) ^ 2 := by
  unfold quarticIncrementSecondDeriv
  ring

/-- If `0 < lam < a*`, the velocity Hessian is positive everywhere. -/
theorem chainQuartic_hessian_pos (lam y : ℝ)
    (hlam0 : 0 < lam) (hlam : lam < modeExtremalCurvature) :
    0 < modeExtremalCurvature + quarticIncrementSecondDeriv lam y := by
  rw [chainQuartic_hessian_identity]
  have hs : 0 ≤ (y - (1 / 2 : ℝ)) ^ 2 := sq_nonneg _
  nlinarith

/-- The quartic momentum map has the displayed positive derivative. -/
theorem chainQuarticD2_hasDerivAt (lam x y : ℝ) :
    HasDerivAt (chainQuarticD2 lam x)
      (modeExtremalCurvature + quarticIncrementSecondDeriv lam y) y := by
  have hbase : HasDerivAt
      (chainCurvedD2 modeExtremalCurvature x)
      modeExtremalCurvature y := by
    unfold chainCurvedD2
    convert ((hasDerivAt_const y (chainFiberSlope x)).add
      ((hasDerivAt_id y).const_mul modeExtremalCurvature)).sub_const
        (modeExtremalCurvature / 2) using 1
    all_goals ring_nf
  have hquartic := quarticIncrementDeriv_hasDerivAt lam y
  simpa only [chainQuarticD2] using hbase.add hquartic

/-- For `0 < lam < a*`, the momentum map is strictly increasing. -/
theorem chainQuarticD2_strictMono (lam : ℝ)
    (hlam0 : 0 < lam) (hlam : lam < modeExtremalCurvature) (x : ℝ) :
    StrictMono (chainQuarticD2 lam x) := by
  apply strictMono_of_deriv_pos
  intro y
  rw [(chainQuarticD2_hasDerivAt lam x y).deriv]
  exact chainQuartic_hessian_pos lam y hlam0 hlam

/-- The positive-Hessian quartic momentum map is onto.  The proof uses
an explicit linear comparison and the intermediate value theorem. -/
theorem chainQuarticD2_surjective (lam : ℝ)
    (hlam0 : 0 < lam) (hlam : lam < modeExtremalCurvature) (x : ℝ) :
    Function.Surjective (chainQuarticD2 lam x) := by
  intro p
  let c := modeExtremalCurvature - lam
  have hc : 0 < c := by dsimp [c]; linarith
  have hcne : c ≠ 0 := ne_of_gt hc
  let z := (p - chainQuarticD2 lam x 0) / c
  have hz_eq : chainQuarticD2 lam x 0 + c * z = p := by
    dsimp [z]
    field_simp [hcne]
    ring
  have haux : 0 < 4 * z ^ 2 - 6 * z + 3 := by
    nlinarith [sq_nonneg (2 * z - (3 / 2 : ℝ))]
  have hfactor : chainQuarticD2 lam x z
      = chainQuarticD2 lam x 0 + c * z
          + lam * z * (4 * z ^ 2 - 6 * z + 3) := by
    dsimp [c]
    unfold chainQuarticD2 quarticIncrementDeriv chainCurvedD2
    ring
  have hcont : Continuous (chainQuarticD2 lam x) := by
    unfold chainQuarticD2 quarticIncrementDeriv chainCurvedD2
    fun_prop
  by_cases hz : 0 ≤ z
  · have h0p : chainQuarticD2 lam x 0 ≤ p := by
      rw [← hz_eq]
      exact le_add_of_nonneg_right (mul_nonneg hc.le hz)
    have hpfz : p ≤ chainQuarticD2 lam x z := by
      rw [hfactor, hz_eq]
      exact le_add_of_nonneg_right
        (mul_nonneg (mul_nonneg hlam0.le hz) haux.le)
    exact intermediate_value_univ 0 z hcont ⟨h0p, hpfz⟩
  · have hzle : z ≤ 0 := le_of_not_ge hz
    have hcz : c * z ≤ 0 :=
      mul_nonpos_of_nonneg_of_nonpos hc.le hzle
    have hp0 : p ≤ chainQuarticD2 lam x 0 := by
      linarith [hz_eq, hcz]
    have hrest : lam * z * (4 * z ^ 2 - 6 * z + 3) ≤ 0 :=
      mul_nonpos_of_nonpos_of_nonneg
        (mul_nonpos_of_nonneg_of_nonpos hlam0.le hzle) haux.le
    have hfzp : chainQuarticD2 lam x z ≤ p := by
      rw [hfactor, hz_eq]
      linarith
    exact intermediate_value_univ z 0 hcont ⟨hfzp, hp0⟩

/-- A choice of velocity solving the quartic momentum relation. -/
noncomputable def chainQuarticVelocitySolver (lam : ℝ)
    (hlam0 : 0 < lam) (hlam : lam < modeExtremalCurvature)
    (x p : ℝ) : ℝ :=
  Classical.choose (chainQuarticD2_surjective lam hlam0 hlam x p)

theorem chainQuarticVelocitySolver_spec (lam : ℝ)
    (hlam0 : 0 < lam) (hlam : lam < modeExtremalCurvature)
    (x p : ℝ) :
    chainQuarticD2 lam x
      (chainQuarticVelocitySolver lam hlam0 hlam x p) = p :=
  Classical.choose_spec (chainQuarticD2_surjective lam hlam0 hlam x p)

/-- The chosen velocity solves every quartic momentum equation. -/
theorem chainQuartic_velocitySolver_solves (lam : ℝ)
    (hlam0 : 0 < lam) (hlam : lam < modeExtremalCurvature) :
    SolvesMomentum (chainQuarticLagrangian lam)
      (chainQuarticVelocitySolver lam hlam0 hlam) := by
  intro x p
  have h := chainQuartic_momentum lam x
    (chainQuarticVelocitySolver lam hlam0 hlam x p)
  rw [chainQuarticVelocitySolver_spec lam hlam0 hlam x p] at h
  exact h

/-- Every quartic momentum equation has at most one velocity solution. -/
theorem chainQuartic_momentum_unique (lam : ℝ)
    (hlam0 : 0 < lam) (hlam : lam < modeExtremalCurvature)
    (x p v₁ v₂ : ℝ)
    (h₁ : MomentumRelation (chainQuarticLagrangian lam) x v₁ p)
    (h₂ : MomentumRelation (chainQuarticLagrangian lam) x v₂ p) :
    v₁ = v₂ := by
  have hc₁ := (chainQuartic_momentum lam x v₁).deriv
  have hc₂ := (chainQuartic_momentum lam x v₂).deriv
  have hp₁ := h₁.deriv
  have hp₂ := h₂.deriv
  apply (chainQuarticD2_strictMono lam hlam0 hlam x).injective
  calc
    chainQuarticD2 lam x v₁ = p := hc₁.symm.trans hp₁
    _ = chainQuarticD2 lam x v₂ := (hc₂.symm.trans hp₂).symm

/-- A nonzero quartic continuation differs from the committed rule. -/
theorem chainQuartic_ne_committed (lam : ℝ) (hlam : lam ≠ 0) :
    chainQuarticLagrangian lam
      ≠ chainCurvedLagrangian modeExtremalCurvature := by
  intro h
  have h2 := congrFun (congrFun h 0) 2
  simp only [chainQuarticLagrangian, quarticIncrement] at h2
  ring_nf at h2
  apply hlam
  linarith

/-- A nonzero quartic continuation is not any sentence of the declared
quadratic grammar. -/
theorem chainQuartic_ne_quadGrammar (lam : ℝ) (hlam : lam ≠ 0)
    (α β γ δ ε ζ : ℝ) :
    chainQuarticLagrangian lam
      ≠ fun x y =>
          chainLogLagrangian x y + quadPoly α β γ δ ε ζ x y := by
  intro h
  have h0 := congrFun (congrFun h 0) 0
  have h1 := congrFun (congrFun h 0) 1
  have h2 := congrFun (congrFun h 0) 2
  have h3 := congrFun (congrFun h 0) 3
  simp only [chainQuarticLagrangian, chainCurvedLagrangian,
    quarticIncrement, quadPoly] at h0 h1 h2 h3
  ring_nf at h0 h1 h2 h3
  apply hlam
  linarith

/-- **Explicit global-selection boundary.**  The parameter
`lam = a*/2` gives a nonquadratic, corner-invisible, velocity-only,
stationary continuation which retains the exact corner table, global
real single-site minimality, a positive Hessian, and a unique global
momentum inverse.  Thus the current tested clauses do not select the
committed rule beyond the declared polynomial subfamilies. -/
theorem quartic_continuation_blocks_global_selection :
    ∃ lam : ℝ,
      0 < lam ∧ lam < modeExtremalCurvature
        ∧ CornerInvisible (fun _ y => quarticIncrement lam y)
        ∧ VelocityOnly (fun _ y => quarticIncrement lam y)
        ∧ chainQuarticD2 lam 1 1 + chainCurvedD1 1 1 = 0
        ∧ (∀ i j : Fin 2,
            chainQuarticLagrangian lam ((i : ℕ) : ℝ) ((j : ℕ) : ℝ)
              = -Real.log (chainWeight i j))
        ∧ (∀ M : ℕ, ∀ {k m : Fin M}, k.succ = m.castSucc →
            ∀ x : ℝ,
              localAction (chainQuarticLagrangian lam)
                  (chainEmb (constOneHistory M))
                ≤ localAction (chainQuarticLagrangian lam)
                    (Function.update (chainEmb (constOneHistory M))
                      k.succ x))
        ∧ (∀ _x y : ℝ,
            0 < modeExtremalCurvature
              + quarticIncrementSecondDeriv lam y)
        ∧ (∃ vel : ℝ → ℝ → ℝ,
            SolvesMomentum (chainQuarticLagrangian lam) vel)
        ∧ chainQuarticLagrangian lam
            ≠ chainCurvedLagrangian modeExtremalCurvature
        ∧ ∀ α β γ δ ε ζ : ℝ,
            chainQuarticLagrangian lam
              ≠ fun x y =>
                  chainLogLagrangian x y
                    + quadPoly α β γ δ ε ζ x y := by
  let lam := modeExtremalCurvature / 2
  have hlam0 : 0 < lam := by
    dsimp [lam]
    exact half_pos modeExtremalCurvature_pos
  have hlam : lam < modeExtremalCurvature := by
    dsimp [lam]
    exact half_lt_self modeExtremalCurvature_pos
  have hlamne : lam ≠ 0 := ne_of_gt hlam0
  refine ⟨lam, hlam0, hlam,
    quarticIncrement_cornerInvisible lam,
    quarticIncrement_velocityOnly lam,
    chainQuartic_stationary lam,
    chainQuartic_corner lam,
    ?_, ?_, ⟨chainQuarticVelocitySolver lam hlam0 hlam,
      chainQuartic_velocitySolver_solves lam hlam0 hlam⟩,
    chainQuartic_ne_committed lam hlamne, ?_⟩
  · intro M k m hkm x
    exact chainQuartic_constOne_realMin lam hlam0.le M hkm x
  · intro x y
    exact chainQuartic_hessian_pos lam y hlam0 hlam
  · intro α β γ δ ε ζ
    exact chainQuartic_ne_quadGrammar lam hlamne α β γ δ ε ζ

end OPH.Variational

#print axioms OPH.Variational.quadPoly_cornerInvisible_iff
#print axioms OPH.Variational.quadPoly_cornerInvisible_normal_form
#print axioms OPH.Variational.chainTwoSlot_eq_increment
#print axioms OPH.Variational.twoSlotIncrement_cornerInvisible
#print axioms OPH.Variational.quadPoly_cornerInvisible_iff_twoSlot
#print axioms OPH.Variational.chainTwoSlot_parameters_injective
#print axioms OPH.Variational.enrichedQuad_hasFDerivAt
#print axioms OPH.Variational.quadPoly_stationary_iff
#print axioms OPH.Variational.chainTwoSlot_stationary_line
#print axioms OPH.Variational.chainTwoSlot_stationary_existsUnique
#print axioms OPH.Variational.twoSlotIncrement_velocityOnly_iff
#print axioms OPH.Variational.quadPoly_velocityOnly_iff
#print axioms OPH.Variational.chainTwoSlot_zero_eq_curved
#print axioms OPH.Variational.chainTwoSlot_velocityOnly_selection
#print axioms OPH.Variational.quadPoly_selection_unique
#print axioms OPH.Variational.quadPoly_selection_eq_committed
#print axioms OPH.Variational.velocityOnly_clause_necessary
#print axioms OPH.Variational.cubicIncrement_cornerInvisible
#print axioms OPH.Variational.cubicIncrement_velocityOnly
#print axioms OPH.Variational.chainCubic_corner
#print axioms OPH.Variational.cubicIncrement_hasDerivAt
#print axioms OPH.Variational.chainCubic_hasFDerivAt
#print axioms OPH.Variational.chainCubic_stationary_iff
#print axioms OPH.Variational.chainCubic_ne_twoSlot
#print axioms OPH.Variational.chainCubic_ne_quadGrammar
#print axioms OPH.Variational.quadraticity_clause_necessary
#print axioms OPH.Variational.chainCubic_gap
#print axioms OPH.Variational.chainCubic_realMin_forces_quadratic
#print axioms OPH.Variational.quarticIncrement_cornerInvisible
#print axioms OPH.Variational.quarticIncrement_velocityOnly
#print axioms OPH.Variational.chainQuartic_corner
#print axioms OPH.Variational.quarticIncrement_hasDerivAt
#print axioms OPH.Variational.quarticIncrementDeriv_hasDerivAt
#print axioms OPH.Variational.chainQuartic_momentum
#print axioms OPH.Variational.chainQuartic_hasFDerivAt
#print axioms OPH.Variational.chainQuartic_stationary
#print axioms OPH.Variational.chainQuartic_gap
#print axioms OPH.Variational.chainQuartic_constOne_realMin
#print axioms OPH.Variational.chainQuartic_hessian_identity
#print axioms OPH.Variational.chainQuartic_hessian_pos
#print axioms OPH.Variational.chainQuarticD2_hasDerivAt
#print axioms OPH.Variational.chainQuarticD2_strictMono
#print axioms OPH.Variational.chainQuarticD2_surjective
#print axioms OPH.Variational.chainQuarticVelocitySolver_spec
#print axioms OPH.Variational.chainQuartic_velocitySolver_solves
#print axioms OPH.Variational.chainQuartic_momentum_unique
#print axioms OPH.Variational.chainQuartic_ne_committed
#print axioms OPH.Variational.chainQuartic_ne_quadGrammar
#print axioms OPH.Variational.quartic_continuation_blocks_global_selection

-- Expected axioms for every theorem above: propext, Classical.choice,
-- Quot.sound (real analysis via Mathlib).  No native_decide, no decide.
