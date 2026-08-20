import QFT.SectorInheritance
import QFT.PathTimeSliceInterface
import QFT.TripleCarrierJoin
import LocalFaceMaxwellAction

/-!
# The E4 structural-inheritance status matrix

This module records the current partial E4 status of completion-plan issue
`#701` by representing all seven inheritance targets as one first-class
artifact.  It does not close that issue.  Issue `#701` names seven targets —
CPT,
spin-statistics, DHR/sector reconstruction, field/action
reconstruction, KMS/thermality, relative-Cauchy stress response, and
particle/scattering — and licenses exactly three verdict shapes per
target: a theorem, a countermodel, or a typed non-evaluable exit.  No
blanket QFT inheritance claim is allowed, and none is made here.

**Summary, stated so it cannot be quoted as a positive inheritance
result: the constructed net does not inherit QFT structure.  Of the
seven targets, five are not statable against the committed carriers,
one has a refutable restricted surrogate (every bijective localized
endomorphism is inner), and one
carries a finite provable core whose net-level reading is degenerate
and whose public-level reading is refuted.**  The five non-statable
rows are the point of this document, not a gap in it: each names the
precise structure that would have to exist first, with the citation.

The matrix is data (`e4StatusMatrix`), the counts behind the summary
sentence are kernel-checked (`e4_theorem_shaped_rows`,
`e4_exit_or_gated_rows`), and each verdict of the two built rows is
pinned to its carrying theorem by a restatement below, so a change to
any cited declaration breaks this module rather than silently
orphaning the matrix.

## The seven rows

| # | Target | Verdict shape | Carrier / reason |
|---|--------|---------------|------------------|
| 1 | CPT | TYPED NON-EVALUABLE EXIT | No carrier.  `EventRegion` is order-theoretic data of the E1 net with no involution, region reflection, or distinguished vacuum; the claim boundary of `QFT/LocallyCovariantLimit.lean` attaches no Lorentzian reading.  In pinned Mathlib the Tomita conjugation is an explicit TODO (`Mathlib/Analysis/InnerProductSpace/StandardSubspace.lean`, "Define the Tomita conjugation, prove Tomita's theorem, prove the KMS condition"). |
| 2 | Spin-statistics | TYPED NON-EVALUABLE EXIT | No carrier.  The hom-type `CausalEmbedding` is proof-only and `Subsingleton` (receipt `causalEmbedding_proof_irrelevant` below), so the event-region category is thin and carries no symmetry-group action data, hence no rotation double cover; no spin structure exists in the tree.  `Mathlib/LinearAlgebra/CliffordAlgebra/SpinGroup.lean` exists but nothing connects it to the net; the statistics side is additionally gated on row 3's exit. |
| 3 | DHR / sector reconstruction | RESTRICTED SURROGATE COLLAPSES + exit (full DHR) | `QFT/SectorInheritance.lean`: every bijective localized endomorphism of the finite constructed net is inner (`privateSector_localized_inner`) and hence equivalent to the identity under the file's declared finite relation; nonvacuity by `sectorWitness_inner`.  This does not construct a DHR category or perform Doplicher–Roberts reconstruction.  Full DHR (quasi-local algebra, vacuum representation, transportability) is that module's typed exit. |
| 4 | Field / action reconstruction | GATED BY DESIGN | Exact finite precursors now include the response-word algebra, the dense global Hodge comparison, the genuinely local face-curvature/static scalar-vector action of `Screen/LocalFaceMaxwellAction.lean`, and the declared-step leapfrog evolution of `Screen/TemporalMaxwellEvolution.lean` with its proved Gauss-continuity equivalence, exact energy balance, and discrete wave law.  The step index is a declared evolution parameter, and the sources are declared histories: no source-selected physical current, physical matter/Spin action, spacetime carrier, physical time, or continuum limit is supplied.  Thus the full field/action reconstruction packet E4 would consume is still not supplied. |
| 5 | KMS / thermality | PROVABLE (finite core) + REFUTABLE (net and public) + exit (horizon) | `QFT/StructuralInheritance.lean`: the finite algebraic Gibbs-KMS theorem on the witness private block (`privateGibbsState_kms`); the degeneracy countermodel — on the support-graded net every functional is KMS for every grading-preserving flow at every complex time (`supportGradedNet_regionalExpectation_kms_degenerate`); the cited B3 rigidity refuting public thermality (`publicRecord_no_thermal_flow`).  Horizon thermality is that module's typed exit. |
| 6 | Relative-Cauchy stress response | TYPED NON-EVALUABLE EXIT | A finite one-step interface exists in `QFT/PathTimeSliceInterface.lean`.  `QFT/TripleCarrierJoin.lean` additionally gives one common state/path/checkpoint carrier and exact marginal intertwining over the 31 actual source transitions.  `QFT/LocallyCovariantLimit.lean` now supplies the order-theoretic E3 Cauchy embedding class `IsCauchyEmbedding`, defined through the net's contravariant restriction maps, with the limit-level time-slice property; its claim boundary attaches no Lorentzian reading.  None of these supplies a regional-net morphism, coherent physical evolution, Lorentzian metric perturbations, or a stress-response derivative. |
| 7 | Particle / scattering | TYPED NON-EVALUABLE EXIT | No carrier.  Haag–Ruelle theory needs a translation group action on the net, a spectrum condition, and asymptotic limits; `EventRegion` carries no ℝ-action and no energy-momentum object exists in the tree (B9 `#685`, the optional spectral adapter, has no module).  Pinned Mathlib has zero Haag/scattering material. |

Four of the five non-statable rows (1, 2, 6, 7) still lack structures
required by the full target in the tree and in pinned Mathlib; row 4 is
gated by design on three other bricks.  This is a statement about what
the committed carriers presently type, not a no-go against extending
them.  In particular row 6 now has the finite transport-aware precursor and
the order-theoretic E3 Cauchy embedding class recorded above, while
Lorentzian Cauchy geometry, coherent physical evolution and a
stress-response derivative remain possible construction routes rather
than being ruled out by this matrix.

## Downstream gates

### F1 (`#694`, Einstein-branch continuation)

F1's readout completion adds modular, stress, entropy, and scale
readouts.  What E4 supplies to F1:

* the finite algebraic Gibbs-KMS theorem on the private block
  (`e4Row5_provable_anchor`), a finite-block algebra fact F1's modular
  readout may consume as algebra — with no physical temperature, no
  spacetime, and no vacuum attached;
* two constraints F1 must respect rather than results it may cite: the
  public record algebra admits no nontrivial pointwise-continuous star
  flow (`e4Row5_public_refutation_anchor`), so no F1 modular readout
  can be carried by a public flow; and on the constructed example net
  the KMS condition distinguishes no state and no temperature
  (`e4Row5_degeneracy_anchor`), so a modular readout must bring its own
  selection structure and may not point at E4 for one.

What E4 withholds from F1: everything in rows 1, 6, and the row-5
horizon exit.  Row 6 now acknowledges a finite time-indexed
net-evolution interface and the order-theoretic E3 Cauchy embedding
class of `QFT/LocallyCovariantLimit.lean`, but the stress readout
still receives no relative-Cauchy response: the embedding class is
order-theoretic with no Lorentzian reading, and there is no Lorentzian
metric or metric perturbation, relative-Cauchy evolution, or stress
derivative.  No statement of this module or of the two row
modules may be cited by F1 as "the net is thermal" or "the net carries
stress response".  Where thermality was testable, E4 refuted it.

### H2 (`#697`, particle masses and mixings)

H2's promotion rule requires "the E3/E4 quantum, spin-statistics,
pole, and continuum receipts".  What E4 supplies to H2:

* exactly two finite receipts: the finite Gibbs-KMS core (row 5) and
  innerness of every bijective localized endomorphism in the restricted
  finite surrogate (row 3);
* one constraint: no H2 candidate may infer a nontrivial charge or flavor
  quantum number from that restricted automorphism class.  Non-bijective
  transportable endomorphisms and full DHR reconstruction are not available.
  Flavor carriers must come from completed source-attached field and matter
  packets, which row 4 records as not yet delivered despite exact finite
  algebraic and global-Hodge precursors.

What E4 withholds from H2: the spin-statistics receipt (row 2 is a
typed exit — there is no spin structure and the thin category carries
no action data), the pole/particle receipt (row 7 is a typed exit —
there is no ℝ-action and no spectrum condition), and the physical or
norm-completed continuum receipt.  `QFT/LocallyCovariantLimit.lean`
constructs the algebraic filtered colimit and its order-theoretic local,
public, locality, refinement-independence, and time-slice properties;
`QFT/ColimitNormedCompletion.lean` supplies the inductive-limit C*-seminorm,
its exact kernel, and the C*-completion of that colimit; neither provides a
physical Lorentzian attachment.
The E4 clause of H2's promotion rule is therefore NOT dischargeable
today, and H2 may not quote E4 as having supplied its quantum receipts
in general.  H2 physical promotion remains blocked on structure the
committed carriers do not presently supply.

## What this module does not license

No statement here may be quoted as "the constructed net inherits QFT
structure", as a CPT, spin-statistics, stress-response, or scattering
result of any polarity, or as a thermality or superselection claim
beyond the two finite theorems cited by name above.  The matrix is the
result: two finite theorems, five precisely cited absences.
-/

namespace OPH.QFT

open OPH.Tower

/-! ## The matrix as data -/

/-- A verdict shape licensed by issue `#701`: a theorem
(`provable`), a countermodel (`refutable`), a typed non-evaluable
exit, or (for row 4 only) an explicit design gate on other bricks. -/
inductive E4Verdict
  | provable
  | refutable
  | typedNonEvaluableExit
  | gated
  deriving DecidableEq

/-- One row of the E4 status matrix: the target, its verdict shape
(possibly composite, e.g. row 5 is provable core + countermodel +
exit), the carrying module if one exists, and the citation for the
verdict — for exit rows, the precise missing structure. -/
structure E4Row where
  index : Nat
  target : String
  shape : List E4Verdict
  carrier : Option String
  citation : String

/-- The seven-row E4 status matrix, in issue order.  The docstring of
this module is the prose form; the anchors below pin the built rows'
verdicts to their carrying theorems. -/
def e4StatusMatrix : List E4Row :=
  [{ index := 1
     target := "CPT"
     shape := [.typedNonEvaluableExit]
     carrier := none
     citation := "No antiunitary conjugation, region reflection, or \
       distinguished vacuum on EventRegion (order-theoretic only, \
       QFT/LocallyCovariantLimit.lean claim boundary); Tomita \
       conjugation is an explicit TODO in pinned Mathlib \
       (Analysis/InnerProductSpace/StandardSubspace.lean)." },
   { index := 2
     target := "spin-statistics"
     shape := [.typedNonEvaluableExit]
     carrier := none
     citation := "CausalEmbedding is proof-only and Subsingleton \
       (causalEmbedding_proof_irrelevant), so the category is thin and \
       carries no symmetry-group action data; no spin structure in \
       tree; Mathlib SpinGroup.lean unconnected to the net." },
   { index := 3
     target := "DHR / sector reconstruction"
     shape := [.refutable, .typedNonEvaluableExit]
     carrier := some "QFT/SectorInheritance.lean"
     citation := "Restricted finite surrogate: every bijective localized \
       endomorphism is inner (privateSector_localized_inner) and equivalent \
       to the identity under the declared relation; no DHR category or \
       reconstructed gauge group is constructed. Full DHR is that module's \
       typed exit (no quasi-local algebra, no vacuum representation)." },
   { index := 4
     target := "field / action reconstruction"
     shape := [.gated]
     carrier := none
     citation := "Screen/LocalFaceMaxwellAction.lean supplies an exact local \
       finite static face-curvature action and separate scalar Coulomb \
       action, and Screen/TemporalMaxwellEvolution.lean supplies a \
       declared-step leapfrog evolution with the proved Gauss-continuity \
       equivalence on the same carrier.  Physical \
       source/matter/Spin attachment, physical time, spacetime covariance, \
       and a continuum reconstruction packet are absent." },
   { index := 5
     target := "KMS / thermality"
     shape := [.provable, .refutable, .typedNonEvaluableExit]
     carrier := some "QFT/StructuralInheritance.lean"
     citation := "Finite Gibbs-KMS theorem (privateGibbsState_kms); \
       net-level degeneracy countermodel \
       (supportGradedNet_regionalExpectation_kms_degenerate) and \
       public-flow rigidity (publicRecord_no_thermal_flow); horizon \
       thermality is that module's typed exit." },
   { index := 6
     target := "relative-Cauchy stress response"
     shape := [.typedNonEvaluableExit]
     carrier := none
     citation := "QFT/PathTimeSliceInterface.lean supplies a finite \
       one-step net-evolution interface; QFT/TripleCarrierJoin.lean supplies \
       a common state/path/checkpoint carrier and 31-step marginal \
       intertwining; QFT/LocallyCovariantLimit.lean supplies the \
       order-theoretic E3 Cauchy embedding class with the limit time-slice \
       property and no Lorentzian reading; but no net morphism, Lorentzian \
       metric or metric perturbation, relative-Cauchy evolution, or \
       stress-response derivative exists; the required \
       Lorentzian/globally-hyperbolic carrier is also absent from pinned \
       Mathlib." },
   { index := 7
     target := "particle / scattering"
     shape := [.typedNonEvaluableExit]
     carrier := none
     citation := "No translation action on EventRegion, no spectrum \
       condition object, B9 (#685) absent from tree; zero \
       Haag/scattering material in pinned Mathlib." }]

/-! ## The kernel-checked summary counts

The summary sentence of this module quotes three numbers: seven rows,
two rows carrying a theorem-shaped verdict, five rows carrying no
theorem of either polarity.  Each is checked by the kernel against the
matrix data, so the summary cannot drift from the matrix. -/

theorem e4StatusMatrix_complete : e4StatusMatrix.map (·.index) = [1, 2, 3, 4, 5, 6, 7] := rfl

/-- Exactly one row carries a provable component: row 5. -/
theorem e4_provable_rows :
    (e4StatusMatrix.filter fun r =>
      r.shape.any fun v => v matches .provable).map (·.index) = [5] := rfl

/-- Exactly two rows carry a refutable component: rows 3 and 5. -/
theorem e4_refutable_rows :
    (e4StatusMatrix.filter fun r =>
      r.shape.any fun v => v matches .refutable).map (·.index) = [3, 5] := rfl

/-- The rows carrying any theorem-shaped verdict, of either polarity,
are exactly rows 3 and 5 — two of seven. -/
theorem e4_theorem_shaped_rows :
    (e4StatusMatrix.filter fun r =>
      r.shape.any fun v => v matches .provable | .refutable).map (·.index) =
      [3, 5] := rfl

/-- The rows carrying no theorem of either polarity — typed exits and
the design gate — are exactly rows 1, 2, 4, 6, 7: five of seven. -/
theorem e4_exit_or_gated_rows :
    (e4StatusMatrix.filter fun r =>
      r.shape.all fun v => v matches .typedNonEvaluableExit | .gated).map (·.index) =
      [1, 2, 4, 6, 7] := rfl

/-- Every non-statable row names its missing structure: no row has an
empty citation, and a row lacks a carrier module exactly when it
carries no theorem-shaped verdict. -/
theorem e4_no_silent_blanks :
    (e4StatusMatrix.filter fun r =>
      r.carrier.isNone != (r.shape.all fun v =>
        v matches .typedNonEvaluableExit | .gated)).map (·.index) = [] := rfl

/-! ## Row anchors

Each verdict of the two built rows, and the thinness reason of row 2's
exit, restated against its carrying theorem.  These are citations that
the kernel checks: renaming or weakening any cited declaration breaks
this module. -/

/-- Row 2 reason receipt: parallel causal embeddings are equal — the
event-region category is thin, so it carries no nontrivial
symmetry-group action data, which is why spin-statistics is not
statable against it.  Cites the `Subsingleton` instance of
`QFT/LocallyCovariantLimit.lean`. -/
theorem causalEmbedding_proof_irrelevant {ι : Type*} [Preorder ι]
    {T : ConsensusTower ι} {N : FiniteCausalObserverNet T}
    {U V : EventRegion N} (f g : CausalEmbedding U V) : f = g :=
  Subsingleton.elim f g

/-- Row 3 anchor, restricted-surrogate half: every bijective localized
endomorphism of the witness private block is inner.  Cites
`privateSector_localized_inner` of `QFT/SectorInheritance.lean`. -/
theorem e4Row3_refutable_anchor
    (ρ : ConsensusTower.PrivateAlgebra witnessTower () →⋆ₐ[ℂ]
      ConsensusTower.PrivateAlgebra witnessTower ())
    {U : Finset (Fin 2)} (hloc : LocalizedIn ρ U)
    (hbij : Function.Bijective ρ) :
    ∃ u ∈ unitary (ConsensusTower.PrivateAlgebra witnessTower ()),
      ∀ x, ρ x = u * x * star u :=
  privateSector_localized_inner ρ hloc hbij

/-- Row 3 anchor, nonvacuity: the sector theorem applies to a
demonstrably nontrivial localized endomorphism.  Cites
`sectorWitness_inner` of `QFT/SectorInheritance.lean`. -/
theorem e4Row3_nonvacuity_anchor :
    ∃ u ∈ unitary (Matrix (Fin 2) (Fin 2) ℂ),
      ∀ x, sectorWitness x = u * x * star u :=
  sectorWitness_inner

/-- Row 6 presence anchor: the finite time-indexed net-evolution
interface is genuinely inhabited on the committed 86/247 carrier.
This is only a presence receipt.  It supplies no Cauchy embedding,
Lorentzian perturbation, relative-Cauchy evolution, or stress
derivative, so the relative-Cauchy stress-response row remains a typed
non-evaluable exit. -/
theorem e4Row6_interface_anchor :
    Nonempty (TimeIndexedNetEvolution PairIndex247 twoSlotNet247) :=
  ⟨sourceGeneratorEvolution⟩

/-- Row 5 anchor, provable half: the finite algebraic Gibbs-KMS
theorem on the witness private block.  Cites `privateGibbsState_kms`
of `QFT/StructuralInheritance.lean`. -/
theorem e4Row5_provable_anchor
    (H : ConsensusTower.PrivateAlgebra witnessTower ()) (β : ℝ)
    (a b : ConsensusTower.PrivateAlgebra witnessTower ()) :
    gibbsState H β (a * heisenbergFlow H b (Complex.I * β)) =
      gibbsState H β (b * a) :=
  privateGibbsState_kms H β a b

/-- Row 5 anchor, net-level refutable half: the A3 selected-state
regional pairing of the constructed example net satisfies the KMS
identity for every observer, region, support-graded Hamiltonian, and
complex time — thermality that distinguishes nothing.  Cites
`supportGradedNet_regionalExpectation_kms_degenerate` of
`QFT/StructuralInheritance.lean`. -/
theorem e4Row5_degeneracy_anchor
    (o : witnessTower.Observer ()) {U : Finset (Fin 2)}
    {H a b : Matrix (Fin 2) (Fin 2) ℂ}
    (hH : H ∈ supportDiagonalAlgebra totalSupportRegion)
    (ha : a ∈ supportGradedNet.localAlgebra () U)
    (hb : b ∈ supportGradedNet.localAlgebra () U) (z : ℂ) :
    (witnessSelectedState o * (a * heisenbergFlow H b z)).trace =
      (witnessSelectedState o * (b * a)).trace :=
  supportGradedNet_regionalExpectation_kms_degenerate o hH ha hb z

/-- Row 5 anchor, public refutable half: the public record algebra
admits no nontrivial pointwise-continuous star flow, so public
thermality is refuted rather than unproved.  Cites
`publicRecord_no_thermal_flow` of `QFT/StructuralInheritance.lean`,
itself citing B3's `ContinuousPublicStarFlow.toAut_eq_refl`. -/
theorem e4Row5_public_refutation_anchor {κ : Type*} [Fintype κ]
    [DecidableEq κ] (A : OPH.Dynamics.ContinuousPublicStarFlow κ) (t : ℝ) :
    A.toAut t = StarAlgEquiv.refl :=
  publicRecord_no_thermal_flow A t

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.e4StatusMatrix
#print axioms OPH.QFT.e4StatusMatrix_complete
#print axioms OPH.QFT.e4_provable_rows
#print axioms OPH.QFT.e4_refutable_rows
#print axioms OPH.QFT.e4_theorem_shaped_rows
#print axioms OPH.QFT.e4_exit_or_gated_rows
#print axioms OPH.QFT.e4_no_silent_blanks
#print axioms OPH.QFT.causalEmbedding_proof_irrelevant
#print axioms OPH.QFT.e4Row3_refutable_anchor
#print axioms OPH.QFT.e4Row3_nonvacuity_anchor
#print axioms OPH.QFT.e4Row6_interface_anchor
#print axioms OPH.QFT.e4Row5_provable_anchor
#print axioms OPH.QFT.e4Row5_degeneracy_anchor
#print axioms OPH.QFT.e4Row5_public_refutation_anchor
