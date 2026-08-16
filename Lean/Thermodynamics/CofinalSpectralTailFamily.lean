import CoherentRefinementFamily

set_option autoImplicit false

namespace OPH.Thermodynamics

/-!
# Genuine cofinal regulator families for the PR-08 spectral-tail route

This module strengthens `CoherentRefinementFamily` by adding the ambient
regulator semantics that directedness alone does not provide.  The finite
stage family is embedded order-reflectingly in an explicit ambient preorder;
its image is cofinal; and the ambient preorder has no maximal regulator.

The interface is deliberately stronger than a renamed directed family:

* `AmbientCofinality` contains an explicit ambient regulator type, an order
  embedding of the finite-stage index, a cofinality receipt, and a no-maximal
  regulator receipt;
* `CofinalSpectralTailFamily` combines that interface with the already checked
  coherent spectral-tail family and an explicit strict carrier-growth receipt;
* `cofinalLadderFamily` is a nonconstant, unbounded-cardinality inhabitant;
* a generic incompatibility theorem and the two counterexamples at the end
  prove that the old uniformly bounded-cardinality structure cannot coexist
  with this strict-growth route on the same family;
* `FourLawCofinalConclusions` is a replacement conclusion record.  Its common
  finite-law core contains no refinement field, while every regulator-family
  third-law field is stated directly on the supplied cofinal spectral-tail family.
  Thus the new theorem does not expose a degenerate uniform-gap family as the
  thermodynamic refinement conclusion.

Nothing here produces a physical regulator, energy, inverse temperature, or
repair law from the OPH source axioms.  The ambient/cofinal data and the
spectral envelope are explicit hypotheses, so register row PR-08 is narrowed
but not discharged by this formal interface alone.

Strict finite-carrier growth is one sufficient, deliberately concrete route
to genuine refinement.  This module does not claim that it is necessary for
every possible PR-08 route; a different physical carrier may need a different
data-level nondegeneracy receipt rather than this particular cardinal ladder.
-/

universe u

/-! ## Ambient regulator semantics -/

/-- An explicit genuinely cofinal embedding of finite-stage indices into an
ambient regulator preorder.  `stage` is an order embedding rather than merely
a monotone map, so ambient progress reflects an actual refinement relation in
the index.  `cofinal` says every ambient regulator is reached or exceeded;
`ambient_noMax` rules out a terminal regulator and, together with cofinality,
  rules out singleton or terminal index systems.  Data-level progress is the
  separate `strict_carrier_progress` field below. -/
structure AmbientCofinality (I : Type u) [Preorder I] where
  /-- The ambient regulator type. -/
  Regulator : Type u
  /-- The ambient refinement preorder. -/
  [preorderRegulator : Preorder Regulator]
  /-- A declared ambient regulator exists; this prevents the entire package
  from being inhabited vacuously on empty index and regulator types. -/
  [nonemptyRegulator : Nonempty Regulator]
  /-- Order-reflecting placement of finite stages in the ambient system. -/
  stage : I ↪o Regulator
  /-- The finite-stage image is cofinal in the ambient regulator system. -/
  cofinal : ∀ a : Regulator, ∃ i : I, a ≤ stage i
  /-- The ambient regulator system has no maximal stage. -/
  ambient_noMax : ∀ a : Regulator, ∃ b : Regulator, a < b

attribute [instance] AmbientCofinality.preorderRegulator
attribute [instance] AmbientCofinality.nonemptyRegulator

namespace AmbientCofinality

variable {I : Type u} [Preorder I]

/-- Every represented finite stage has a strictly later represented stage.
This conclusion genuinely consumes the ambient no-maximal-stage receipt,
cofinality, and order reflection. -/
theorem exists_strict_later (C : AmbientCofinality I) (i : I) :
    ∃ j : I, i < j := by
  obtain ⟨a, hia⟩ := C.ambient_noMax (C.stage i)
  obtain ⟨j, haj⟩ := C.cofinal a
  refine ⟨j, C.stage.lt_iff_lt.mp ?_⟩
  exact lt_of_lt_of_le hia haj

/-- A genuine cofinal no-max regulator attachment forces infinitely many
finite-stage indices.  This is the nonvacuity theorem that excludes both
empty and every finite stage system. -/
theorem index_infinite (C : AmbientCofinality I) : Infinite I := by
  let a : C.Regulator := Classical.choice inferInstance
  obtain ⟨i0, _⟩ := C.cofinal a
  letI : Nonempty I := ⟨i0⟩
  rw [← Set.infinite_univ_iff]
  apply Set.infinite_of_forall_exists_gt
  intro i
  obtain ⟨j, hij⟩ := C.exists_strict_later i
  exact ⟨j, Set.mem_univ j, hij⟩

/-- The ambient preorder is directed whenever the finite-stage index is
directed.  The proof first lifts the ambient pair to cofinal finite stages and
then uses the finite family's common refinement. -/
theorem ambient_directed (C : AmbientCofinality I)
    (hdir : ∀ i j : I, ∃ k : I, i ≤ k ∧ j ≤ k) :
    ∀ a b : C.Regulator, ∃ c : C.Regulator, a ≤ c ∧ b ≤ c := by
  intro a b
  obtain ⟨i, hai⟩ := C.cofinal a
  obtain ⟨j, hbj⟩ := C.cofinal b
  obtain ⟨k, hik, hjk⟩ := hdir i j
  exact ⟨C.stage k,
    le_trans hai (C.stage.monotone hik),
    le_trans hbj (C.stage.monotone hjk)⟩

end AmbientCofinality

/-! ## Cofinal spectral-tail families -/

/-- A coherent finite spectral-tail family equipped with genuine ambient
cofinal regulator semantics and a data-level nondegeneracy receipt.  The
strict carrier-growth field prevents a constant finite carrier from passing
merely because its index is embedded cofinally in a larger preorder.  It is a
sufficient genuine-refinement route, not a necessity claim about all possible
PR-08 regulator models. -/
structure CofinalSpectralTailFamily {I : Type u} [Preorder I]
    (X : I → Type*) [∀ r, Fintype (X r)] [∀ r, DecidableEq (X r)]
    [∀ r, Nonempty (X r)] (E : ∀ r, X r → ℝ) where
  /-- The finite coherent family and its uniform spectral-tail envelope. -/
  family : CoherentRefinementFamily X E
  /-- The explicit ambient cofinal regulator attachment. -/
  ambient : AmbientCofinality I
  /-- Every represented stage has a strictly later stage with strictly more
  finite carrier data. -/
  strict_carrier_progress : ∀ r : I, ∃ s : I,
    r < s ∧ Fintype.card (X r) < Fintype.card (X s)

namespace CofinalSpectralTailFamily

variable {I : Type u} [Preorder I] {X : I → Type*}
  [∀ r, Fintype (X r)] [∀ r, DecidableEq (X r)] [∀ r, Nonempty (X r)]
  {E : ∀ r, X r → ℝ}

/-- Honest promotion from a coherent family requires both the entire ambient
cofinality package and a separate data-level carrier-growth receipt. -/
def ofCoherent (F : CoherentRefinementFamily X E)
    (C : AmbientCofinality I)
    (hprogress : ∀ r : I, ∃ s : I,
      r < s ∧ Fintype.card (X r) < Fintype.card (X s)) :
    CofinalSpectralTailFamily X E :=
  ⟨F, C, hprogress⟩

/-- Forgetting ambient semantics recovers the coherent directed family. -/
def toCoherent (F : CofinalSpectralTailFamily X E) :
    CoherentRefinementFamily X E :=
  F.family

/-- The uniform spectral envelope still gives one temperature threshold for
all stages. -/
theorem uniform_concentration (F : CofinalSpectralTailFamily X E) :
    ∀ eps : ℝ, 0 < eps → ∃ beta1 : ℝ, ∀ beta : ℝ, beta1 ≤ beta →
      ∀ r : I, offMinMass (E r) beta < eps :=
  F.family.uniform_concentration

/-- **Cofinal-tail concentration.** Beyond one inverse-temperature threshold,
every ambient regulator admits a represented stage above it such that that
stage and every finer finite stage satisfy the same concentration bound.
This is stronger than merely appending a cofinality proposition beside an
unrelated concentration theorem: the quantified energy is the energy of the
cofinally embedded family itself. -/
theorem cofinal_tail_concentration (F : CofinalSpectralTailFamily X E) :
    ∀ eps : ℝ, 0 < eps → ∃ beta1 : ℝ, ∀ beta : ℝ, beta1 ≤ beta →
      ∀ a : F.ambient.Regulator, ∃ r : I,
        a ≤ F.ambient.stage r ∧
        ∀ s : I, r ≤ s →
          a ≤ F.ambient.stage s ∧ offMinMass (E s) beta < eps := by
  intro eps heps
  obtain ⟨beta1, hbeta1⟩ := F.family.uniform_concentration eps heps
  refine ⟨beta1, fun beta hbeta a ↦ ?_⟩
  obtain ⟨r, har⟩ := F.ambient.cofinal a
  exact ⟨r, har, fun s hrs ↦
    ⟨le_trans har (F.ambient.stage.monotone hrs), hbeta1 beta hbeta s⟩⟩

/-- Genuine cofinality supplies strict finite-stage progress. -/
theorem exists_strict_later (F : CofinalSpectralTailFamily X E) (r : I) :
    ∃ s : I, r < s :=
  F.ambient.exists_strict_later r

/-- The supplied data-level receipt gives a later stage with strictly larger
finite carrier.  Unlike `exists_strict_later`, this cannot be witnessed by a
constant family living over a nonconstant index. -/
theorem exists_strict_carrier_later (F : CofinalSpectralTailFamily X E)
    (r : I) : ∃ s : I,
      r < s ∧ Fintype.card (X r) < Fintype.card (X s) :=
  F.strict_carrier_progress r

/-- Genuine cofinality forces an infinite finite-stage index. -/
theorem index_infinite (F : CofinalSpectralTailFamily X E) : Infinite I :=
  F.ambient.index_infinite

/-- The ambient regulator system is directed, derived from the coherent
family's directed finite-stage index and the cofinal embedding. -/
theorem ambient_directed (F : CofinalSpectralTailFamily X E) :
    ∀ a b : F.ambient.Regulator, ∃ c : F.ambient.Regulator,
      a ≤ c ∧ b ≤ c :=
  F.ambient.ambient_directed F.family.directed

/-- A uniformly bounded-cardinality refinement cannot also have strict
carrier growth from every stage of a nonempty index.  Iterating the supplied
progress choice produces a strictly increasing sequence of natural-number
cardinalities, while `UniformGapRefinement.card_le` bounds every term by one
fixed natural number.  This is why there is no promotion constructor from the
old interface to this module's strict-growth interface. -/
theorem uniformGap_strictCarrierProgress_false
    (G : UniformGapRefinement X E) [Nonempty I]
    (hprogress : ∀ r : I, ∃ s : I,
      r < s ∧ Fintype.card (X r) < Fintype.card (X s)) : False := by
  let next : I → I := fun r ↦ Classical.choose (hprogress r)
  have hnext (r : I) :
      Fintype.card (X r) < Fintype.card (X (next r)) :=
    (Classical.choose_spec (hprogress r)).2
  let r0 : I := Classical.choice inferInstance
  let cards : ℕ → ℕ := fun n ↦ Fintype.card (X ((next^[n]) r0))
  have hcards_succ (n : ℕ) : cards n < cards (n + 1) := by
    change Fintype.card (X ((next^[n]) r0)) <
      Fintype.card (X ((next^[n + 1]) r0))
    rw [Function.iterate_succ_apply']
    exact hnext _
  have hcards : StrictMono cards := strictMono_nat_of_lt_succ hcards_succ
  have hlower : G.cardBound + 1 ≤ cards (G.cardBound + 1) :=
    hcards.id_le (G.cardBound + 1)
  have hupper : cards (G.cardBound + 1) ≤ G.cardBound :=
    G.card_le _
  exact Nat.not_succ_le_self G.cardBound (le_trans hlower hupper)

/-- In particular, an old uniform-gap object and a cofinal spectral-tail
object with this module's strict carrier progress cannot inhabit the same
indexed carriers and energies.  Ambient cofinality supplies the nonempty
index needed by the generic incompatibility theorem. -/
theorem uniformGap_cofinalFamily_isEmpty
    (G : UniformGapRefinement X E) :
    IsEmpty (CofinalSpectralTailFamily X E) := by
  refine ⟨fun F ↦ ?_⟩
  let a : F.ambient.Regulator := Classical.choice inferInstance
  obtain ⟨r, _⟩ := F.ambient.cofinal a
  letI : Nonempty I := ⟨r⟩
  exact uniformGap_strictCarrierProgress_false G F.strict_carrier_progress

end CofinalSpectralTailFamily

/-! ## A nonconstant, unbounded-cardinality cofinal inhabitant -/

/-- Inclusion of ladder stage `n` into stage `n + 1`, preserving both the
energy level and its degeneracy label. -/
def ladderSuccInclusion (n : ℕ) (x : LadderState n) :
    LadderState (n + 1) :=
  ⟨⟨(x.1 : ℕ), lt_trans x.1.isLt (Nat.lt_succ_self (n + 1))⟩, x.2⟩

/-- Truncating the successor-stage inclusion recovers the original ladder
state. -/
theorem ladderTrunc_succInclusion (n : ℕ) (x : LadderState n) :
    ladderTrunc n (ladderSuccInclusion n x) = x := by
  unfold ladderSuccInclusion ladderTrunc
  rw [dif_pos (Nat.lt_succ_iff.mp x.1.isLt)]

/-- The successor-stage inclusion is injective. -/
theorem ladderSuccInclusion_injective (n : ℕ) :
    Function.Injective (ladderSuccInclusion n) :=
  Function.LeftInverse.injective (ladderTrunc_succInclusion n)

/-- One state at the new top level of ladder stage `n + 1`. -/
def ladderNewTop (n : ℕ) : LadderState (n + 1) :=
  ⟨⟨n + 1, Nat.lt_succ_self (n + 1)⟩,
    ⟨0, pow_pos (by norm_num) _⟩⟩

/-- The new top-level state is not in the image of the previous stage. -/
theorem ladderNewTop_notMem_range (n : ℕ) :
    ladderNewTop n ∉ Set.range (ladderSuccInclusion n) := by
  rintro ⟨x, hx⟩
  have hlevel := congrArg
    (fun z : LadderState (n + 1) ↦ (z.1 : ℕ)) hx
  exact (Nat.ne_of_lt x.1.isLt) hlevel

/-- Every successor ladder stage has strictly more carrier states. -/
theorem ladder_card_strict_succ (n : ℕ) :
    Fintype.card (LadderState n) <
      Fintype.card (LadderState (n + 1)) :=
  Fintype.card_lt_of_injective_of_notMem
    (ladderSuccInclusion n) (ladderSuccInclusion_injective n)
    (ladderNewTop_notMem_range n)

/-- The exponential ladder with its natural-number stages embedded identically
in the ambient natural-number regulator.  This is genuinely cofinal and has no
terminal stage. -/
noncomputable def cofinalLadderFamily (delta : ℝ) (hdelta : 0 < delta) :
    CofinalSpectralTailFamily LadderState (ladderEnergy delta) where
  family := ladderFamily delta hdelta
  ambient :=
    { Regulator := ℕ
      stage := (OrderIso.refl ℕ).toOrderEmbedding
      cofinal := fun a ↦ ⟨a, le_rfl⟩
      ambient_noMax := fun a ↦ ⟨a + 1, Nat.lt_succ_self a⟩ }
  strict_carrier_progress := fun n ↦
    ⟨n + 1, Nat.lt_succ_self n, ladder_card_strict_succ n⟩

/-- The cofinal ladder's carriers are unbounded in cardinality. -/
theorem cofinalLadder_card_unbounded (m : ℕ) :
    m < Fintype.card (LadderState m) :=
  ladder_card_unbounded m

/-- The cofinal ladder has an explicitly later represented stage. -/
theorem cofinalLadder_strict_progress (delta : ℝ) (hdelta : 0 < delta)
    (n : ℕ) : ∃ m : ℕ, n < m :=
  (cofinalLadderFamily delta hdelta).exists_strict_later n

/-- The cofinal ladder's explicit data-level progress receipt. -/
theorem cofinalLadder_strict_carrier_progress (delta : ℝ)
    (hdelta : 0 < delta) (n : ℕ) : ∃ m : ℕ,
      n < m ∧ Fintype.card (LadderState n) <
        Fintype.card (LadderState m) :=
  (cofinalLadderFamily delta hdelta).exists_strict_carrier_later n

/-! ## A genuinely rethreaded four-law conclusion -/

open FourLawSurface

/-- The refinement-independent finite four-law core.  This record is stated
only on one repair law and one calibration; it contains every conclusion of
`FourLawConclusions` except the two fields tied to its old uniform-gap
refinement member. -/
structure FiniteFourLawCore {Omega : Type u} [Fintype Omega]
    [DecidableEq Omega] [Nonempty Omega] {B : Type*} [DecidableEq B]
    (D : RepairLawData Omega B) (C : CalibrationData D) : Prop where
  kernel_nonneg : ∀ x y, 0 ≤ repairKernel D x y
  kernel_row_sum : ∀ x, ∑ y, repairKernel D x y = 1
  kernel_stationary : ∀ y, push D.ref (repairKernel D) y = D.ref y
  kernel_detailed_balance : ∀ x y,
    D.ref x * repairKernel D x y = D.ref y * repairKernel D y x
  kernel_idempotent : ∀ x y,
    ∑ z, repairKernel D x z * repairKernel D z y = repairKernel D x y
  kernel_row_optimal : ∀ (x : Omega) (r : Omega → ℝ), (∀ y, 0 ≤ r y) →
    (∑ y, r y = 1) →
    (∀ y, D.visible y ≠ D.visible x → r y = 0) →
    kl (repairKernel D x) D.ref ≤ kl r D.ref
  zeroth_ref_eq_gibbs : ∀ x, D.ref x = gibbs C.energy C.beta x
  zeroth_thermometer : ∀ (beta1 beta2 : ℝ) (i j : Omega),
    C.energy i ≠ C.energy j →
    gibbs C.energy beta1 = gibbs C.energy beta2 → beta1 = beta2
  zeroth_multiplier_unique : ∀ C2 : CalibrationData D,
    C2.energy = C.energy → ∀ i j : Omega, C.energy i ≠ C.energy j →
      C.beta = C2.beta
  first_cap_exact : ∀ p : Omega → ℝ,
    shannon p - shannon D.ref
      = (∑ x, p x * (-Real.log (D.ref x)))
        - (∑ x, D.ref x * (-Real.log (D.ref x))) - kl p D.ref
  first_repair_conserves_fibre_mean : ∀ Q : Omega → ℝ,
    (∀ x y, D.visible x = D.visible y → Q x = Q y) →
    ∀ p : Omega → ℝ,
      ∑ y, push p (repairKernel D) y * Q y = ∑ x, p x * Q x
  first_heat_channel : ∀ p : Omega → ℝ,
    (internalEnergy (diagR (push p (repairKernel D))) (diagR C.energy)
        - internalEnergy (diagR p) (diagR C.energy)
      = heatIncrement (diagR (push p (repairKernel D)) - diagR p)
          (diagR C.energy))
    ∧ heatIncrement (diagR (push p (repairKernel D)) - diagR p)
        (diagR C.energy)
      = ∑ x, (push p (repairKernel D) x - p x) * C.energy x
  second_contraction : ∀ p : Omega → ℝ, (∀ x, 0 ≤ p x) →
    kl (push p (repairKernel D)) D.ref ≤ kl p D.ref
  second_contraction_of_stationary :
    ∀ (p : Omega → ℝ) (K : Omega → Omega → ℝ), (∀ x, 0 ≤ p x) →
    (∀ x y, 0 ≤ K x y) → (∀ x, ∑ y, K x y = 1) →
    push D.ref K = D.ref → kl (push p K) D.ref ≤ kl p D.ref
  second_clausius : ∀ p : Omega → ℝ, (∀ x, 0 ≤ p x) →
    (∑ x, push p (repairKernel D) x * (-Real.log (D.ref x)))
      - (∑ x, p x * (-Real.log (D.ref x)))
      ≤ shannon (push p (repairKernel D)) - shannon p
  second_integral_fluctuation : ∀ p : Omega → ℝ, (∀ x, 0 < p x) →
    (∑ x, p x = 1) →
    ∑ x, ∑ y, p x * repairKernel D x y
      * Real.exp (-(sigmaEP D.ref p (push p (repairKernel D)) x y)) = 1
  second_crooks : ∀ p : Omega → ℝ, (∀ x, 0 < p x) → ∀ x y,
    p x * repairKernel D x y
      = Real.exp (sigmaEP D.ref p (push p (repairKernel D)) x y)
        * (push p (repairKernel D) y * repairKernel D y x)
  arrow_mean_entropy_production_nonneg : ∀ p : Omega → ℝ,
    (∀ x, 0 < p x) →
    0 ≤ ∑ x, ∑ y, p x * repairKernel D x y
      * sigmaEP D.ref p (push p (repairKernel D)) x y
  arrow_mean_ep_eq_kl_to_repaired : ∀ p : Omega → ℝ,
    (∀ x, 0 < p x) →
    ∑ x, ∑ y, p x * repairKernel D x y
      * sigmaEP D.ref p (push p (repairKernel D)) x y
      = kl p (push p (repairKernel D))
  arrow_strict : ∀ p : Omega → ℝ, (∀ x, 0 < p x) →
    (push p (repairKernel D) ≠ p ↔
      0 < ∑ x, ∑ y, p x * repairKernel D x y
        * sigmaEP D.ref p (push p (repairKernel D)) x y)
  repair_fixed_iff : ∀ p : Omega → ℝ,
    (push p (repairKernel D) = p ↔
      ∀ y, p y = D.ref y * fiberMass p D.visible y
        / fiberMass D.ref D.visible y)
  third_excited_mass_bound : ∀ beta E0 gap : ℝ, 0 ≤ beta → 0 < gap →
    (∀ x, C.energy x = E0 ∨ E0 + gap ≤ C.energy x) →
    (Finset.univ.filter (fun x ↦ C.energy x = E0)).Nonempty →
    excitedMass C.energy beta E0
      ≤ ((Finset.univ.filter (fun x ↦ C.energy x ≠ E0)).card : ℝ)
        / ((Finset.univ.filter (fun x ↦ C.energy x = E0)).card : ℝ)
        * Real.exp (-beta * gap)
  third_excited_mass_threshold : ∀ beta E0 gap eps : ℝ, 0 ≤ beta →
    0 < gap → 0 < eps →
    (∀ x, C.energy x = E0 ∨ E0 + gap ≤ C.energy x) →
    (Finset.univ.filter (fun x ↦ C.energy x = E0)).Nonempty →
    Real.log
        (((Finset.univ.filter (fun x ↦ C.energy x ≠ E0)).card : ℝ)
          / ((Finset.univ.filter (fun x ↦ C.energy x = E0)).card : ℝ)
          / eps) < beta * gap →
    excitedMass C.energy beta E0 < eps
  third_no_step_extinguishes : ∀ p : Omega → ℝ, (∀ x, 0 < p x) →
    ∀ y, 0 < push p (repairKernel D) y
  landauer_bound : ∀ (p : Omega → ℝ) (c : ℝ), (∀ x, 0 ≤ p x) →
    shannon (push p (repairKernel D)) - shannon p ≤ -c →
    (∑ x, p x * C.energy x)
      - (∑ x, push p (repairKernel D) x * C.energy x) ≥ c / C.beta

/-- Extract the refinement-independent core from the committed conclusion
record.  Its old PR-08 fields are intentionally not projected. -/
def FiniteFourLawCore.ofFourLawConclusions {Omega : Type u}
    [Fintype Omega] [DecidableEq Omega] [Nonempty Omega]
    {B : Type*} [DecidableEq B] (A : FourLawAntecedent Omega B)
    (H : FourLawConclusions A) : FiniteFourLawCore A.repair A.calib where
  kernel_nonneg := H.kernel_nonneg
  kernel_row_sum := H.kernel_row_sum
  kernel_stationary := H.kernel_stationary
  kernel_detailed_balance := H.kernel_detailed_balance
  kernel_idempotent := H.kernel_idempotent
  kernel_row_optimal := H.kernel_row_optimal
  zeroth_ref_eq_gibbs := H.zeroth_ref_eq_gibbs
  zeroth_thermometer := H.zeroth_thermometer
  zeroth_multiplier_unique := H.zeroth_multiplier_unique
  first_cap_exact := H.first_cap_exact
  first_repair_conserves_fibre_mean := H.first_repair_conserves_fibre_mean
  first_heat_channel := H.first_heat_channel
  second_contraction := H.second_contraction
  second_contraction_of_stationary := H.second_contraction_of_stationary
  second_clausius := H.second_clausius
  second_integral_fluctuation := H.second_integral_fluctuation
  second_crooks := H.second_crooks
  arrow_mean_entropy_production_nonneg := H.arrow_mean_entropy_production_nonneg
  arrow_mean_ep_eq_kl_to_repaired := H.arrow_mean_ep_eq_kl_to_repaired
  arrow_strict := H.arrow_strict
  repair_fixed_iff := H.repair_fixed_iff
  third_excited_mass_bound := H.third_excited_mass_bound
  third_excited_mass_threshold := H.third_excited_mass_threshold
  third_no_step_extinguishes := H.third_no_step_extinguishes
  landauer_bound := H.landauer_bound

/-- Direct construction of the refinement-independent finite-law core from
one repair law and one calibration.  Unlike the compatibility projection
above, this proof never constructs a uniform-gap or degenerate refinement
family. -/
theorem finiteFourLawCore {Omega : Type u} [Fintype Omega]
    [DecidableEq Omega] [Nonempty Omega] {B : Type*} [DecidableEq B]
    (D : RepairLawData Omega B) (C : CalibrationData D) :
    FiniteFourLawCore D C where
  kernel_nonneg := repairKernel_nonneg D
  kernel_row_sum := repairKernel_row_sum D
  kernel_stationary := repairKernel_stationary D
  kernel_detailed_balance := repairKernel_detailedBalance D
  kernel_idempotent := repairKernel_idempotent D
  kernel_row_optimal := repairKernel_row_optimal D
  zeroth_ref_eq_gibbs := ref_eq_gibbs D C
  zeroth_thermometer := zeroth_thermometer D C
  zeroth_multiplier_unique := fun C2 hE i j hij ↦
    zeroth_multiplier_unique D C C2 hE.symm i j hij
  first_cap_exact := first_cap_exact D
  first_repair_conserves_fibre_mean := first_repair_conserves_fibre_mean D
  first_heat_channel := first_heat_channel D C
  second_contraction := second_contraction D
  second_contraction_of_stationary := second_contraction_of_stationary D
  second_clausius := second_clausius D
  second_integral_fluctuation := second_integral_fluctuation D
  second_crooks := second_crooks D
  arrow_mean_entropy_production_nonneg :=
    arrow_mean_entropy_production_nonneg D
  arrow_mean_ep_eq_kl_to_repaired :=
    arrow_mean_entropy_production_eq_kl_to_repaired D
  arrow_strict := arrow_strict D
  repair_fixed_iff := repair_fixed_iff D
  third_excited_mass_bound := third_excited_mass_bound D C
  third_excited_mass_threshold := third_excited_mass_threshold D C
  third_no_step_extinguishes := third_no_step_extinguishes D
  landauer_bound := landauer_bound D C

/-- PR-08 attachment with a genuinely cofinal ambient regulator family tied
to the calibrated energy at one distinguished finite stage. -/
structure CofinalRefinementAttachment {Omega : Type u} [Fintype Omega]
    [DecidableEq Omega] {B : Type*} [DecidableEq B]
    (D : RepairLawData Omega B) (C : CalibrationData D) where
  I : Type
  [preorderI : Preorder I]
  X : I → Type u
  [fintypeX : ∀ r, Fintype (X r)]
  [decEqX : ∀ r, DecidableEq (X r)]
  [nonemptyX : ∀ r, Nonempty (X r)]
  E : ∀ r, X r → ℝ
  family : CofinalSpectralTailFamily X E
  r0 : I
  e : X r0 ≃ Omega
  energy_eq : ∀ x, E r0 x = C.energy (e x)

attribute [instance] CofinalRefinementAttachment.preorderI
attribute [instance] CofinalRefinementAttachment.fintypeX
attribute [instance] CofinalRefinementAttachment.decEqX
attribute [instance] CofinalRefinementAttachment.nonemptyX

/-- One typed four-law antecedent whose PR-08 field is the genuinely cofinal
spectral-tail family. -/
structure FourLawCofinalAntecedent (Omega : Type u) [Fintype Omega]
    [DecidableEq Omega] (B : Type*) [DecidableEq B] where
  repair : RepairLawData Omega B
  calib : CalibrationData repair
  refinement : CofinalRefinementAttachment repair calib

/-- The cofinal attachment makes the surface carrier inhabited. -/
theorem FourLawCofinalAntecedent.nonempty {Omega : Type u} [Fintype Omega]
    [DecidableEq Omega] {B : Type*} [DecidableEq B]
    (A : FourLawCofinalAntecedent Omega B) : Nonempty Omega :=
  A.refinement.e.nonempty_congr.mp
    (A.refinement.nonemptyX A.refinement.r0)

/-- Compatibility bridge that forgets only ambient cofinality and retains the
same coherent spectral-tail family.  The rethreaded theorem below does not
convert this object to the old uniform-gap antecedent. -/
def FourLawCofinalAntecedent.toCoherent {Omega : Type u} [Fintype Omega]
    [DecidableEq Omega] {B : Type*} [DecidableEq B]
    (A : FourLawCofinalAntecedent Omega B) :
    FourLawCoherentAntecedent Omega B where
  repair := A.repair
  calib := A.calib
  refinement :=
    { I := A.refinement.I
      X := A.refinement.X
      E := A.refinement.E
      family := A.refinement.family.family
      r0 := A.refinement.r0
      e := A.refinement.e
      energy_eq := A.refinement.energy_eq }

/-- The final cofinal thermodynamic conclusion type.  Its finite core is
independent of any refinement representation.  All refinement/continuum
third-law clauses refer directly to the supplied cofinal family. -/
structure FourLawCofinalConclusions {Omega : Type u} [Fintype Omega]
    [DecidableEq Omega] [Nonempty Omega] {B : Type*} [DecidableEq B]
    (A : FourLawCofinalAntecedent Omega B) : Prop where
  /-- All exact finite zeroth/first/second/third and Landauer conclusions,
  with no refinement member hidden in the type. -/
  finite : FiniteFourLawCore A.repair A.calib
  /-- The one envelope bounds every represented regulator stage. -/
  third_refinement_envelope : ∀ beta : ℝ,
    A.refinement.family.family.beta0 ≤ beta →
    ∀ r : A.refinement.I,
      offMinMass (A.refinement.E r) beta
        ≤ A.refinement.family.family.Env beta
  /-- One inverse-temperature threshold controls all represented stages. -/
  third_uniform_concentration : ∀ eps : ℝ, 0 < eps →
    ∃ beta1 : ℝ, ∀ beta : ℝ, beta1 ≤ beta →
      ∀ r : A.refinement.I,
        offMinMass (A.refinement.E r) beta < eps
  /-- The same bound holds on a genuine cofinal tail of the ambient regulator
  system, not merely on an unrelated directed index. -/
  third_cofinal_tail : ∀ eps : ℝ, 0 < eps →
    ∃ beta1 : ℝ, ∀ beta : ℝ, beta1 ≤ beta →
      ∀ a : A.refinement.family.ambient.Regulator,
        ∃ r : A.refinement.I,
          a ≤ A.refinement.family.ambient.stage r ∧
          ∀ s : A.refinement.I, r ≤ s →
            a ≤ A.refinement.family.ambient.stage s
            ∧ offMinMass (A.refinement.E s) beta < eps
  /-- The calibrated energy is genuinely one member of this same family and
  inherits its envelope. -/
  third_calibrated_member : ∀ beta : ℝ,
    offMinMass A.calib.energy beta
      = offMinMass (A.refinement.E A.refinement.r0) beta
    ∧ (A.refinement.family.family.beta0 ≤ beta →
      offMinMass A.calib.energy beta
        ≤ A.refinement.family.family.Env beta)
  /-- The finite-stage system has no terminal represented stage. -/
  strict_refinement_progress : ∀ r : A.refinement.I,
    ∃ s : A.refinement.I, r < s
  /-- The supplied family also makes strict progress in its finite carrier
  data, rather than only in the names of its represented stages. -/
  strict_carrier_refinement_progress : ∀ r : A.refinement.I,
    ∃ s : A.refinement.I, r < s ∧
      Fintype.card (A.refinement.X r) <
        Fintype.card (A.refinement.X s)
  /-- In particular, the represented refinement-stage type is infinite. -/
  infinite_refinement_stages : Infinite A.refinement.I

/-- **Rethreaded four-law theorem.** The exact finite-law core is extracted
without the old PR-08 projections.  The spectral envelope, uniform
concentration, cofinal-tail concentration, calibrated-member conclusion, and
strict regulator and carrier progress are all proved for the supplied
genuinely cofinal family itself. -/
theorem fourLaws_composed_cofinal {Omega : Type u} [Fintype Omega]
    [DecidableEq Omega] [Nonempty Omega] {B : Type*} [DecidableEq B]
    (A : FourLawCofinalAntecedent Omega B) :
    FourLawCofinalConclusions A := by
  have hmember : ∀ beta : ℝ,
      offMinMass A.calib.energy beta
        = offMinMass (A.refinement.E A.refinement.r0) beta := by
    intro beta
    have hE : A.refinement.E A.refinement.r0
        = fun x ↦ A.calib.energy (A.refinement.e x) :=
      funext A.refinement.energy_eq
    rw [hE]
    exact (offMinMass_equiv A.refinement.e A.calib.energy beta).symm
  refine
    { finite := finiteFourLawCore A.repair A.calib
      third_refinement_envelope := A.refinement.family.family.env_bound
      third_uniform_concentration :=
        A.refinement.family.family.uniform_concentration
      third_cofinal_tail := A.refinement.family.cofinal_tail_concentration
      third_calibrated_member := fun beta ↦ ⟨hmember beta, fun hbeta ↦ ?_⟩
      strict_refinement_progress := A.refinement.family.exists_strict_later
      strict_carrier_refinement_progress :=
        A.refinement.family.strict_carrier_progress
      infinite_refinement_stages := A.refinement.family.index_infinite }
  rw [hmember beta]
  exact A.refinement.family.family.env_bound beta hbeta A.refinement.r0

/-! ## Non-nesting and load-bearing negative controls -/

/-- A genuine ambient cofinality package cannot exist on a singleton index.
This is the false-green guard against silently relabeling a constant family as
cofinal. -/
theorem ambientCofinality_pUnit_isEmpty :
    IsEmpty (AmbientCofinality (PUnit.{1})) := by
  refine ⟨fun C ↦ ?_⟩
  obtain ⟨j, hj⟩ := C.exists_strict_later (PUnit.unit : PUnit.{1})
  have hju : j = (PUnit.unit : PUnit.{1}) := Subsingleton.elim _ _
  rw [hju] at hj
  exact (lt_irrefl (PUnit.unit : PUnit.{1})) hj

/-- Consequently, no cofinal spectral-tail family can use a singleton stage
index, regardless of its finite carrier and energy. -/
theorem cofinalFamily_pUnit_isEmpty {X : PUnit.{1} → Type*}
    [∀ r, Fintype (X r)] [∀ r, DecidableEq (X r)] [∀ r, Nonempty (X r)]
    (E : ∀ r, X r → ℝ) :
    IsEmpty (CofinalSpectralTailFamily X E) :=
  ⟨fun F ↦ ambientCofinality_pUnit_isEmpty.false F.ambient⟩

/-- **Constant-carrier countercontrol.** No family with one fixed finite
carrier at every stage can satisfy this module's strict carrier-growth route,
even if its index and ambient regulator have arbitrarily rich order
structure.  This rules out the synthetic constant-family countermodel while
making no necessity claim about other possible PR-08 progress notions. -/
theorem constantCarrier_cofinalFamily_isEmpty {I : Type u} [Preorder I]
    (Omega : Type*) [Fintype Omega] [DecidableEq Omega] [Nonempty Omega]
    (E : I → Omega → ℝ) :
    IsEmpty (CofinalSpectralTailFamily (fun _ : I ↦ Omega) E) := by
  refine ⟨fun F ↦ ?_⟩
  let a : F.ambient.Regulator := Classical.choice inferInstance
  obtain ⟨r, _⟩ := F.ambient.cofinal a
  obtain ⟨s, _, hcard⟩ := F.strict_carrier_progress r
  exact (lt_irrefl (Fintype.card Omega)) hcard

/-- **Old does not imply new.** The old uniform-gap structure is inhabited on
a singleton constant family, while a genuine cofinal family on that same
index is impossible. -/
theorem old_uniformGap_not_imply_cofinal (delta : ℝ) :
    Nonempty (UniformGapRefinement (fun _ : PUnit.{1} ↦ Bool)
      (fun _ ↦ twoLevelEnergy delta))
    ∧ IsEmpty (CofinalSpectralTailFamily (fun _ : PUnit.{1} ↦ Bool)
      (fun _ ↦ twoLevelEnergy delta)) := by
  constructor
  · exact ⟨constUniformGap (PUnit.{1}) Bool (twoLevelEnergy delta)⟩
  · exact cofinalFamily_pUnit_isEmpty _

/-- **New does not imply old.** The genuinely cofinal ladder is inhabited,
while its unbounded stage cardinalities make the old uniform-cardinality
structure impossible on the same carriers and energies. -/
theorem cofinal_not_imply_old_uniformGap (delta : ℝ) (hdelta : 0 < delta) :
    Nonempty (CofinalSpectralTailFamily LadderState (ladderEnergy delta))
    ∧ IsEmpty (UniformGapRefinement LadderState (ladderEnergy delta)) := by
  exact ⟨⟨cofinalLadderFamily delta hdelta⟩,
    ladder_uniformGap_isEmpty delta⟩

end OPH.Thermodynamics

#print axioms OPH.Thermodynamics.AmbientCofinality.exists_strict_later
#print axioms OPH.Thermodynamics.AmbientCofinality.index_infinite
#print axioms OPH.Thermodynamics.AmbientCofinality.ambient_directed
#print axioms OPH.Thermodynamics.CofinalSpectralTailFamily.uniform_concentration
#print axioms OPH.Thermodynamics.CofinalSpectralTailFamily.cofinal_tail_concentration
#print axioms OPH.Thermodynamics.CofinalSpectralTailFamily.exists_strict_carrier_later
#print axioms OPH.Thermodynamics.CofinalSpectralTailFamily.uniformGap_strictCarrierProgress_false
#print axioms OPH.Thermodynamics.CofinalSpectralTailFamily.uniformGap_cofinalFamily_isEmpty
#print axioms OPH.Thermodynamics.cofinalLadderFamily
#print axioms OPH.Thermodynamics.cofinalLadder_card_unbounded
#print axioms OPH.Thermodynamics.ladder_card_strict_succ
#print axioms OPH.Thermodynamics.cofinalLadder_strict_carrier_progress
#print axioms OPH.Thermodynamics.FiniteFourLawCore.ofFourLawConclusions
#print axioms OPH.Thermodynamics.finiteFourLawCore
#print axioms OPH.Thermodynamics.FourLawCofinalAntecedent.toCoherent
#print axioms OPH.Thermodynamics.fourLaws_composed_cofinal
#print axioms OPH.Thermodynamics.ambientCofinality_pUnit_isEmpty
#print axioms OPH.Thermodynamics.constantCarrier_cofinalFamily_isEmpty
#print axioms OPH.Thermodynamics.old_uniformGap_not_imply_cofinal
#print axioms OPH.Thermodynamics.cofinal_not_imply_old_uniformGap
-- Expected axioms: propext, Classical.choice, Quot.sound only.
