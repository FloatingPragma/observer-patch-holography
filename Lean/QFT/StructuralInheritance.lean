import Mathlib.Analysis.Normed.Algebra.MatrixExponential
import Mathlib.Analysis.SpecialFunctions.Exponential
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.LinearAlgebra.Matrix.PosDef
import QFT.LocallyCovariantLimit
import Dynamics.PublicAutomorphism

/-!
# The E4 structural-inheritance matrix: the KMS and thermality row

This module opens the E4 inheritance matrix of completion-plan issue
`#701` with its one row that is fully statable against the committed
E1–E3 carriers: KMS and thermality.  Per the issue, the row is treated
as a separate theorem target with explicit hypotheses, a nonvacuous
example, a countermodel, and a typed non-evaluable exit; no blanket QFT
inheritance claim is made or licensed by anything below.  The row has
three explicitly labelled parts, and part 2 is half the result, not a
caveat on part 1.

## Part 1 — PROVABLE: the finite algebraic KMS theorem

On a finite full matrix block the modular flow at imaginary time is a
plain product, so no analytic continuation is involved and no KMS
predicate beyond the boundary identity itself is needed.  For a
Hamiltonian `H` the Heisenberg flow `heisenbergFlow H b z` is
conjugation by `NormedSpace.exp (z • (Complex.I • H))`, the matrix
mirror of the B3 flow `OPH.Dynamics.innerExponentialFlow`.  The
normalised Gibbs functional `gibbsState H β` divides the trace against
`gibbsDensity H β = exp (-(β • H))` by the partition trace
`gibbsPartition H β`.  The receipts are:

* `gibbsPartition_ne_zero` / `gibbsPartition_pos`: for Hermitian `H`
  the partition trace is a nonzero, indeed positive, complex number, by
  the Hermitian square root `gibbsSqrt` of the density — no matrix
  inverse and no nonsingularity side condition is ever taken;
* `gibbsState_one`: the Gibbs functional is genuinely normalised;
* `gibbsState_star_mul_self_nonneg`: the Gibbs functional is positive,
  so it is a state and not merely a normalised functional;
* `gibbsState_kms`: the finite KMS boundary identity
  `ω (a * σ_{iβ} b) = ω (b * a)` at the analytic-continuation point
  `z = I * β`, by trace cyclicity and the exponential unit
  `exp (β • H) * exp (-(β • H)) = 1` obtained from `exp_add_of_commute`
  and `exp_zero` — the `exp_neg`-based-unit route, never `⁻¹`.

The identity itself is algebraic and is stated without a Hermiticity
hypothesis; Hermiticity is exactly what makes the normalisation and
positivity receipts true, and the private-algebra restatements bundle
it back in.  The `privateGibbs*` declarations restate the receipts
verbatim on `ConsensusTower.PrivateAlgebra witnessTower ()`, the E6
witness tower's private block, which is definitionally
`Matrix (Fin 2) (Fin 2) ℂ` (`witnessPrivateAlgebra_eq`).

## Part 2 — REFUTABLE: the degeneracy countermodel, net level and public level

Part 1 does NOT transfer to the constructed net as genuine thermal
structure, and this part proves it.

Net level: on the E3 example `supportGradedNet`, a support-graded
Hamiltonian (one lying in the total region's graded algebra, i.e. a
diagonal one) generates a flow that fixes every element of every
regional algebra pointwise (`supportGraded_flow_fixes_local`), hence
preserves the grading (`supportGraded_flow_preserves`), and the KMS
boundary identity for such a flow holds at the level of algebra
elements before any state is applied (`supportGraded_kms_degenerate`):
`a * σ_z b = b * a` for all regional `a`, `b` and ALL complex `z`, not
only at `z = I * β`.  Consequently EVERY functional whatsoever — not
just Gibbs states — satisfies the KMS identity for every support-graded
Hamiltonian at every complex time
(`supportGraded_kms_degenerate_functional`,
`supportGradedNet_regionalExpectation_kms_degenerate`).  The
constructed net therefore inherits KMS structure only DEGENERATELY: the
KMS condition distinguishes no regional state and no inverse
temperature.  The hypothesis class is inhabited and nontrivial: the
two-site number operator `witnessHamiltonian = diagonal ![0, 1]` is
Hermitian, support-graded, not a scalar, and its flow genuinely moves
the off-diagonal matrix unit (`witnessHamiltonian_flow_moves`), so the
degeneracy is a property of the graded algebras, not of a flow that was
trivial anyway.

Public level: by the cited B3 rigidity theorem
`OPH.Dynamics.ContinuousPublicStarFlow.toAut_eq_refl`, every
pointwise-continuous real-parameter group of star automorphisms of the
finite public record algebra is the identity
(`publicRecord_no_thermal_flow`, `publicRecord_flow_fixes`).  KMS
relative to any pointwise-continuous public star flow is therefore
degenerate, because every such flow is the identity.  The theorem
excludes continuous automorphism groups only, stated over an abstract
finite label type: classical detailed-balance thermality of public
records relative to continuous row-stochastic semigroups, the class
`OPH.Dynamics.positive_unital_iff_stochastic` keeps live, is not
evaluated here.

## Part 3 — TYPED NON-EVALUABLE EXIT: horizon thermality

The "horizon" half of the KMS/thermality target — Bisognano–Wichmann or
Hawking–Unruh thermality of a wedge or horizon state — is not
evaluable against the committed carriers, and this module declines it
rather than approximating it.  The missing types are: no wedge type, no
horizon type, no bifurcation-surface type, and no boost action exist on
`EventRegion`, which is order-theoretic data of the E1 net only; its
own claim boundary (`LocallyCovariantLimit.lean`) states that no
Lorentzian or spacetime reading is attached to causal embeddings.
Until a declared wedge/boost structure exists on the event-region
category, horizon thermality is a typed exit of this row, exactly as
time-slice and the limit algebra are typed exits of E3.

No statement in this module may be quoted as "the observer net
satisfies the KMS condition".  What is proved is: the finite private
block carries the algebraic Gibbs-KMS identity (part 1), and the
constructed net inherits it only in the degenerate form in which every
functional is KMS for every grading-preserving flow (part 2).
-/

namespace OPH.QFT

open OPH.Tower
open scoped ComplexOrder
open scoped Matrix

/-! ## Part 1 — the finite algebraic KMS theorem -/

section FiniteKMS

variable {d : Type*} [Fintype d] [DecidableEq d]

/-- The inner exponential flow on a finite matrix block, the matrix
mirror of `OPH.Dynamics.innerExponentialFlow`: conjugation of `b` by
the exponential of `z • A`, with the inverse side written through the
exponential of the negated generator rather than through `⁻¹`. -/
noncomputable def matrixInnerFlow (A b : Matrix d d ℂ) (z : ℂ) : Matrix d d ℂ :=
  NormedSpace.exp (z • A) * b * NormedSpace.exp (z • (-A))

/-- The Heisenberg flow of the Hamiltonian `H` at complex time `z`:
the inner exponential flow of the generator `Complex.I • H`.  Its real
slice is the usual `e^{itH} b e^{-itH}`; no analytic continuation is
needed to evaluate it at imaginary time in finite dimension. -/
noncomputable def heisenbergFlow (H b : Matrix d d ℂ) (z : ℂ) : Matrix d d ℂ :=
  matrixInnerFlow (Complex.I • H) b z

/-- The unnormalised Gibbs density `e^{-βH}`. -/
noncomputable def gibbsDensity (H : Matrix d d ℂ) (β : ℝ) : Matrix d d ℂ :=
  NormedSpace.exp (-((β : ℂ) • H))

/-- The partition trace `Tr e^{-βH}`. -/
noncomputable def gibbsPartition (H : Matrix d d ℂ) (β : ℝ) : ℂ :=
  (gibbsDensity H β).trace

/-- The normalised Gibbs functional `a ↦ Tr(e^{-βH} a) / Tr e^{-βH}`.
The normalised state, not the raw trace functional, is used, because
the raw trace functional is not a state and the nonvanishing of the
denominator is itself one of the theorems of this row
(`gibbsPartition_ne_zero`). -/
noncomputable def gibbsState (H : Matrix d d ℂ) (β : ℝ) (a : Matrix d d ℂ) : ℂ :=
  (gibbsDensity H β * a).trace / gibbsPartition H β

/-- The exponential unit, right half: `e^A e^{-A} = 1`.  This is the
`exp_neg`-based unit of the scope correction — invertibility of the
Gibbs density by construction, with no `⁻¹` and no nonsingularity side
condition. -/
theorem matrix_exp_mul_exp_neg (A : Matrix d d ℂ) :
    NormedSpace.exp A * NormedSpace.exp (-A) = 1 := by
  rw [← Matrix.exp_add_of_commute A (-A) ((Commute.refl A).neg_right),
    add_neg_cancel, NormedSpace.exp_zero]

/-- The exponential unit, left half: `e^{-A} e^A = 1`. -/
theorem matrix_exp_neg_mul_exp (A : Matrix d d ℂ) :
    NormedSpace.exp (-A) * NormedSpace.exp A = 1 := by
  rw [← Matrix.exp_add_of_commute (-A) A ((Commute.refl A).neg_left),
    neg_add_cancel, NormedSpace.exp_zero]

/-- The Hermitian square root `e^{-(β/2)H}` of the Gibbs density. -/
noncomputable def gibbsSqrt (H : Matrix d d ℂ) (β : ℝ) : Matrix d d ℂ :=
  NormedSpace.exp (-(((β / 2 : ℝ) : ℂ) • H))

theorem gibbsSqrt_mul_self (H : Matrix d d ℂ) (β : ℝ) :
    gibbsSqrt H β * gibbsSqrt H β = gibbsDensity H β := by
  rw [gibbsSqrt, gibbsDensity,
    ← Matrix.exp_add_of_commute _ _ (Commute.refl _)]
  congr 1
  rw [← neg_add, ← add_smul]
  congr 2
  push_cast
  ring

theorem gibbsSqrt_isHermitian {H : Matrix d d ℂ} (hH : H.IsHermitian)
    (β : ℝ) : (gibbsSqrt H β).IsHermitian := by
  have hc : IsSelfAdjoint ((β / 2 : ℝ) : ℂ) := Complex.conj_ofReal (β / 2)
  exact ((hH.smul hc).neg).exp

/-- For Hermitian `H` the partition trace is nonzero: the density is
the square of its Hermitian square root, so a vanishing trace would
force the square root — an exponential, hence a unit — to vanish. -/
theorem gibbsPartition_ne_zero [Nonempty d] {H : Matrix d d ℂ}
    (hH : H.IsHermitian) (β : ℝ) : gibbsPartition H β ≠ 0 := by
  intro h0
  have hM : gibbsSqrt H β = 0 := by
    refine Matrix.trace_conjTranspose_mul_self_eq_zero_iff.mp ?_
    rw [(gibbsSqrt_isHermitian hH β).eq, gibbsSqrt_mul_self]
    exact h0
  exact (Matrix.isUnit_exp (-(((β / 2 : ℝ) : ℂ) • H))).ne_zero hM

/-- For Hermitian `H` the partition trace is positive in the complex
order: it is the trace of `Mᴴ M` for the Hermitian square root `M`. -/
theorem gibbsPartition_pos [Nonempty d] {H : Matrix d d ℂ}
    (hH : H.IsHermitian) (β : ℝ) : 0 < gibbsPartition H β := by
  have hnn : 0 ≤ gibbsPartition H β := by
    have h := (Matrix.posSemidef_conjTranspose_mul_self
      (gibbsSqrt H β)).trace_nonneg
    rwa [(gibbsSqrt_isHermitian hH β).eq, gibbsSqrt_mul_self] at h
  exact lt_of_le_of_ne hnn (Ne.symm (gibbsPartition_ne_zero hH β))

/-- The Gibbs functional is normalised: `ω(1) = 1`. -/
theorem gibbsState_one [Nonempty d] {H : Matrix d d ℂ}
    (hH : H.IsHermitian) (β : ℝ) : gibbsState H β 1 = 1 := by
  unfold gibbsState
  rw [mul_one]
  exact div_self (gibbsPartition_ne_zero hH β)

/-- The unnormalised Gibbs functional is positive on `aᴴ a`: cycling
the Hermitian square root turns the pairing into the trace of
`(a M)ᴴ (a M)`. -/
theorem gibbsExpectation_star_mul_self_nonneg {H : Matrix d d ℂ}
    (hH : H.IsHermitian) (β : ℝ) (a : Matrix d d ℂ) :
    0 ≤ (gibbsDensity H β * (aᴴ * a)).trace := by
  have h1 : (gibbsDensity H β * (aᴴ * a)).trace =
      ((a * gibbsSqrt H β)ᴴ * (a * gibbsSqrt H β)).trace := by
    rw [Matrix.conjTranspose_mul, (gibbsSqrt_isHermitian hH β).eq,
      ← gibbsSqrt_mul_self H β]
    calc (gibbsSqrt H β * gibbsSqrt H β * (aᴴ * a)).trace
        = (gibbsSqrt H β * (gibbsSqrt H β * aᴴ * a)).trace := by
          congr 1; noncomm_ring
      _ = ((gibbsSqrt H β * aᴴ * a) * gibbsSqrt H β).trace :=
          Matrix.trace_mul_comm _ _
      _ = (gibbsSqrt H β * aᴴ * (a * gibbsSqrt H β)).trace := by
          congr 1; noncomm_ring
  rw [h1]
  exact (Matrix.posSemidef_conjTranspose_mul_self _).trace_nonneg

/-- The normalised Gibbs functional is positive, hence a state. -/
theorem gibbsState_star_mul_self_nonneg [Nonempty d] {H : Matrix d d ℂ}
    (hH : H.IsHermitian) (β : ℝ) (a : Matrix d d ℂ) :
    0 ≤ gibbsState H β (aᴴ * a) := by
  unfold gibbsState
  have hnum := gibbsExpectation_star_mul_self_nonneg hH β a
  have hden := gibbsPartition_pos hH β
  exact div_nonneg hnum hden.le

/-- At the analytic-continuation point `z = I β` the Heisenberg flow is
the plain product `e^{-βH} b e^{βH}` — the finite-dimensional modular
flow at imaginary time, with no analytic continuation performed. -/
theorem heisenbergFlow_kms_point (H b : Matrix d d ℂ) (β : ℝ) :
    heisenbergFlow H b (Complex.I * β) =
      gibbsDensity H β * b * NormedSpace.exp ((β : ℂ) • H) := by
  have hscal : (Complex.I * (β : ℂ)) • (Complex.I • H) = -((β : ℂ) • H) := by
    rw [smul_smul, ← neg_smul]
    congr 1
    rw [mul_right_comm, Complex.I_mul_I]
    ring
  have hscal' : (Complex.I * (β : ℂ)) • (-(Complex.I • H)) = (β : ℂ) • H := by
    rw [smul_neg, hscal]
    exact neg_neg _
  unfold heisenbergFlow matrixInnerFlow gibbsDensity
  rw [hscal, hscal']

/-- **The finite algebraic KMS theorem.**  The normalised Gibbs
functional satisfies the KMS boundary identity
`ω (a σ_{iβ}(b)) = ω (b a)` for the inner Heisenberg flow of `H` at
inverse temperature `β`.  The identity is algebraic — trace cyclicity
plus the exponential unit — and holds for every `H`; Hermiticity of `H`
is what additionally makes `gibbsState` a genuine state
(`gibbsState_one`, `gibbsState_star_mul_self_nonneg`). -/
theorem gibbsState_kms (H : Matrix d d ℂ) (β : ℝ) (a b : Matrix d d ℂ) :
    gibbsState H β (a * heisenbergFlow H b (Complex.I * β)) =
      gibbsState H β (b * a) := by
  unfold gibbsState
  congr 1
  rw [heisenbergFlow_kms_point]
  have hED : NormedSpace.exp ((β : ℂ) • H) * gibbsDensity H β = 1 :=
    matrix_exp_mul_exp_neg ((β : ℂ) • H)
  calc (gibbsDensity H β *
        (a * (gibbsDensity H β * b * NormedSpace.exp ((β : ℂ) • H)))).trace
      = ((gibbsDensity H β * a * gibbsDensity H β * b) *
          NormedSpace.exp ((β : ℂ) • H)).trace := by
        congr 1; noncomm_ring
    _ = (NormedSpace.exp ((β : ℂ) • H) *
          (gibbsDensity H β * a * gibbsDensity H β * b)).trace :=
        Matrix.trace_mul_comm _ _
    _ = ((NormedSpace.exp ((β : ℂ) • H) * gibbsDensity H β) *
          (a * (gibbsDensity H β * b))).trace := by
        congr 1; noncomm_ring
    _ = (a * (gibbsDensity H β * b)).trace := by rw [hED, one_mul]
    _ = ((gibbsDensity H β * b) * a).trace := Matrix.trace_mul_comm _ _
    _ = (gibbsDensity H β * (b * a)).trace := by congr 1; noncomm_ring

end FiniteKMS

/-! ### Part 1 restated on the witness tower's private block -/

/-- Carrier receipt: the private algebra of the E6 witness tower at its
sole regulator is definitionally the two-by-two complex matrix
algebra. -/
theorem witnessPrivateAlgebra_eq :
    ConsensusTower.PrivateAlgebra witnessTower () =
      Matrix (Fin 2) (Fin 2) ℂ := rfl

instance : Nonempty (Fin (witnessTower.dim ())) := ⟨(0 : Fin 2)⟩

theorem privateGibbsPartition_ne_zero
    {H : ConsensusTower.PrivateAlgebra witnessTower ()}
    (hH : H.IsHermitian) (β : ℝ) : gibbsPartition H β ≠ 0 :=
  gibbsPartition_ne_zero hH β

theorem privateGibbsState_one
    {H : ConsensusTower.PrivateAlgebra witnessTower ()}
    (hH : H.IsHermitian) (β : ℝ) : gibbsState H β 1 = 1 :=
  gibbsState_one hH β

/-- The finite algebraic KMS theorem on the witness tower's private
block, the carrier named by the E4 issue. -/
theorem privateGibbsState_kms
    (H : ConsensusTower.PrivateAlgebra witnessTower ()) (β : ℝ)
    (a b : ConsensusTower.PrivateAlgebra witnessTower ()) :
    gibbsState H β (a * heisenbergFlow H b (Complex.I * β)) =
      gibbsState H β (b * a) :=
  gibbsState_kms H β a b

/-! ## Part 2 — the degeneracy countermodel on the E3 net -/

section DegenerateInheritance

/-- The total region of the two-site support grading. -/
def totalSupportRegion : Finset (Fin 2) := {0, 1}

theorem mem_totalSupportRegion : ∀ i : Fin 2, i ∈ totalSupportRegion := by
  decide

/-- Entry formula for the Heisenberg flow of a diagonal Hamiltonian:
conjugation by diagonal exponentials scales the `(i, j)` entry by
`e^{z i h_i} e^{-z i h_j}`. -/
theorem heisenbergFlow_diagonal_apply (h : Fin 2 → ℂ)
    (b : Matrix (Fin 2) (Fin 2) ℂ) (z : ℂ) (i j : Fin 2) :
    heisenbergFlow (Matrix.diagonal h) b z i j =
      NormedSpace.exp (z * (Complex.I * h i)) * b i j *
        NormedSpace.exp (-(z * (Complex.I * h j))) := by
  have hdiag : ∀ w : ℂ, w • (Complex.I • Matrix.diagonal h) =
      Matrix.diagonal (w • (Complex.I • h)) := fun w => by
    rw [Matrix.diagonal_smul, Matrix.diagonal_smul]
  have hneg : z • (-(Complex.I • Matrix.diagonal h)) =
      (-z) • (Complex.I • Matrix.diagonal h) := by
    rw [smul_neg, neg_smul]
  unfold heisenbergFlow matrixInnerFlow
  rw [hneg, hdiag z, hdiag (-z), Matrix.exp_diagonal, Matrix.exp_diagonal,
    Matrix.mul_diagonal, Matrix.diagonal_mul, Pi.coe_exp, Pi.coe_exp]
  simp [Pi.smul_apply, smul_eq_mul, neg_mul]

/-- Net-level degeneracy, elementwise: the flow of a support-graded
Hamiltonian fixes every element of every regional algebra pointwise,
at every complex time. -/
theorem supportGraded_flow_fixes_local {H b : Matrix (Fin 2) (Fin 2) ℂ}
    {U : Finset (Fin 2)}
    (hH : H ∈ supportDiagonalAlgebra totalSupportRegion)
    (hb : b ∈ supportDiagonalAlgebra U) (z : ℂ) :
    heisenbergFlow H b z = b := by
  obtain ⟨h, rfl, -⟩ := hH
  obtain ⟨f, rfl, -⟩ := hb
  ext i j
  rw [heisenbergFlow_diagonal_apply]
  by_cases hij : i = j
  · subst hij
    rw [Matrix.diagonal_apply_eq, ← Complex.exp_eq_exp_ℂ]
    rw [Complex.exp_neg, mul_comm (Complex.exp _) (f i), mul_assoc,
      mul_inv_cancel₀ (Complex.exp_ne_zero _), mul_one]
  · rw [Matrix.diagonal_apply_ne _ hij, mul_zero, zero_mul]

/-- Support-graded Hamiltonians are grading preserving: their flow
maps each regional algebra into itself (indeed fixes it pointwise, by
`supportGraded_flow_fixes_local`). -/
theorem supportGraded_flow_preserves {H : Matrix (Fin 2) (Fin 2) ℂ}
    (hH : H ∈ supportDiagonalAlgebra totalSupportRegion)
    {U : Finset (Fin 2)} {b : Matrix (Fin 2) (Fin 2) ℂ}
    (hb : b ∈ supportDiagonalAlgebra U) (z : ℂ) :
    heisenbergFlow H b z ∈ supportDiagonalAlgebra U := by
  rw [supportGraded_flow_fixes_local hH hb z]
  exact hb

/-- **The degeneracy theorem.**  On the E3 support-graded net the KMS
boundary identity for a support-graded Hamiltonian holds at the level
of algebra elements, for every pair of regional observables and every
complex time — before any state is applied and with no distinguished
inverse temperature.  This is the precise sense in which the
constructed net inherits KMS structure only degenerately. -/
theorem supportGraded_kms_degenerate {H a b : Matrix (Fin 2) (Fin 2) ℂ}
    {U : Finset (Fin 2)}
    (hH : H ∈ supportDiagonalAlgebra totalSupportRegion)
    (ha : a ∈ supportDiagonalAlgebra U)
    (hb : b ∈ supportDiagonalAlgebra U) (z : ℂ) :
    a * heisenbergFlow H b z = b * a := by
  rw [supportGraded_flow_fixes_local hH hb z]
  obtain ⟨f, rfl, -⟩ := ha
  obtain ⟨g, rfl, -⟩ := hb
  exact diagonal_commute f g

/-- Every functional whatsoever — linear or not, state or not —
satisfies the KMS identity for every support-graded Hamiltonian at
every complex time.  A thermality notion satisfied by every functional
distinguishes nothing; this is the countermodel to any blanket
KMS-inheritance reading of part 1. -/
theorem supportGraded_kms_degenerate_functional
    {H a b : Matrix (Fin 2) (Fin 2) ℂ} {U : Finset (Fin 2)}
    (hH : H ∈ supportDiagonalAlgebra totalSupportRegion)
    (ha : a ∈ supportDiagonalAlgebra U)
    (hb : b ∈ supportDiagonalAlgebra U) (z : ℂ)
    (φ : Matrix (Fin 2) (Fin 2) ℂ → ℂ) :
    φ (a * heisenbergFlow H b z) = φ (b * a) :=
  congrArg φ (supportGraded_kms_degenerate hH ha hb z)

/-- The A3 selected state of the witness tower, read on the matrix
carrier of the private block. -/
noncomputable def witnessSelectedState (o : witnessTower.Observer ()) :
    Matrix (Fin 2) (Fin 2) ℂ :=
  witnessTower.state () o

/-- The degeneracy read through the E1 pairing: the A3 selected-state
regional pairing of `supportGradedNet` (the trace against the tower
state, which is definitionally the value of
`FiniteCausalObserverNet.regionalExpectation`) satisfies the KMS
identity for every observer, every region, and every support-graded
Hamiltonian at every complex time. -/
theorem supportGradedNet_regionalExpectation_kms_degenerate
    (o : witnessTower.Observer ()) {U : Finset (Fin 2)}
    {H a b : Matrix (Fin 2) (Fin 2) ℂ}
    (hH : H ∈ supportDiagonalAlgebra totalSupportRegion)
    (ha : a ∈ supportGradedNet.localAlgebra () U)
    (hb : b ∈ supportGradedNet.localAlgebra () U) (z : ℂ) :
    (witnessSelectedState o * (a * heisenbergFlow H b z)).trace =
      (witnessSelectedState o * (b * a)).trace := by
  have ha' : a ∈ supportDiagonalAlgebra U := ha
  have hb' : b ∈ supportDiagonalAlgebra U := hb
  rw [supportGraded_kms_degenerate hH ha' hb' z]

/-! ### Nonvacuity of the degeneracy hypothesis -/

/-- The concrete grading-preserving Hamiltonian: the two-site number
operator `diag(0, 1)`. -/
noncomputable def witnessHamiltonian : Matrix (Fin 2) (Fin 2) ℂ :=
  Matrix.diagonal ![0, 1]

theorem witnessHamiltonian_mem :
    witnessHamiltonian ∈ supportDiagonalAlgebra totalSupportRegion :=
  ⟨![0, 1], rfl, fun i hi _ _ => absurd (mem_totalSupportRegion i) hi⟩

theorem witnessHamiltonian_isHermitian : witnessHamiltonian.IsHermitian := by
  refine Matrix.isHermitian_diagonal_of_self_adjoint _ ?_
  show star ![(0 : ℂ), 1] = ![0, 1]
  funext i
  fin_cases i <;> simp

/-- The witness Hamiltonian is not a scalar: the hypothesis class of
the degeneracy theorem contains more than multiples of the identity. -/
theorem witnessHamiltonian_not_mem_empty :
    witnessHamiltonian ∉ supportDiagonalAlgebra (∅ : Finset (Fin 2)) := by
  rintro ⟨f, hf, hc⟩
  have h0 : (0 : ℂ) = f 0 := by
    have h := congrArg (fun M => M 0 0) hf
    simpa [witnessHamiltonian, Matrix.diagonal_apply_eq] using h
  have h1 : (1 : ℂ) = f 1 := by
    have h := congrArg (fun M => M 1 1) hf
    simpa [witnessHamiltonian, Matrix.diagonal_apply_eq] using h
  have h01 : (0 : ℂ) = 1 := by
    rw [h0, h1]
    exact hc 0 (by simp) 1 (by simp)
  exact zero_ne_one h01

/-- The witness Hamiltonian's flow is genuinely nontrivial on the
ambient block: at time `π` it negates the off-diagonal matrix unit.
So the degeneracy theorem is about the graded algebras, not about a
flow that was trivial anyway. -/
theorem witnessHamiltonian_flow_moves :
    heisenbergFlow witnessHamiltonian (Matrix.single 0 1 (1 : ℂ))
      ((Real.pi : ℝ) : ℂ) ≠ Matrix.single 0 1 (1 : ℂ) := by
  intro heq
  have hentry := congrArg (fun M => M 0 1) heq
  rw [witnessHamiltonian] at hentry
  simp only [heisenbergFlow_diagonal_apply, Matrix.single_apply_same] at hentry
  rw [← Complex.exp_eq_exp_ℂ] at hentry
  have h0 : ((Real.pi : ℝ) : ℂ) * (Complex.I * ![(0 : ℂ), 1] 0) = 0 := by
    rw [show (![(0 : ℂ), 1] 0) = 0 from rfl, mul_zero, mul_zero]
  have h1 : -(((Real.pi : ℝ) : ℂ) * (Complex.I * ![(0 : ℂ), 1] 1)) =
      -(((Real.pi : ℝ) : ℂ) * Complex.I) := by
    rw [show (![(0 : ℂ), 1] 1) = 1 from rfl, mul_one]
  rw [h0, h1, Complex.exp_zero, Complex.exp_neg, Complex.exp_pi_mul_I] at hentry
  norm_num at hentry

end DegenerateInheritance

/-! ## Part 2, public side: continuous-automorphism KMS degenerate by the B3 rigidity theorem -/

/-- Citation of B3: every pointwise-continuous real-parameter group of
star automorphisms of the finite public record algebra is the identity.
KMS relative to any such public star flow is degenerate because the flow
is the identity.  Excluded: continuous automorphism groups only, over an
abstract finite label type.  Live: continuous row-stochastic public
semigroups (`OPH.Dynamics.positive_unital_iff_stochastic`), whose
detailed-balance thermality is not evaluated here.  The proof is
`OPH.Dynamics.ContinuousPublicStarFlow.toAut_eq_refl`, cited, not
re-proved. -/
theorem publicRecord_no_thermal_flow {κ : Type*} [Fintype κ] [DecidableEq κ]
    (A : OPH.Dynamics.ContinuousPublicStarFlow κ) (t : ℝ) :
    A.toAut t = StarAlgEquiv.refl :=
  OPH.Dynamics.ContinuousPublicStarFlow.toAut_eq_refl A t

/-- Pointwise form of the public rigidity citation: every continuous
public star flow fixes every public record function at every time. -/
theorem publicRecord_flow_fixes {κ : Type*} [Fintype κ] [DecidableEq κ]
    (A : OPH.Dynamics.ContinuousPublicStarFlow κ) (t : ℝ) (f : κ → ℂ) :
    A.toAut t f = f := by
  rw [publicRecord_no_thermal_flow A t]
  rfl

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.matrix_exp_mul_exp_neg
#print axioms OPH.QFT.matrix_exp_neg_mul_exp
#print axioms OPH.QFT.gibbsSqrt_mul_self
#print axioms OPH.QFT.gibbsSqrt_isHermitian
#print axioms OPH.QFT.gibbsPartition_ne_zero
#print axioms OPH.QFT.gibbsPartition_pos
#print axioms OPH.QFT.gibbsState_one
#print axioms OPH.QFT.gibbsExpectation_star_mul_self_nonneg
#print axioms OPH.QFT.gibbsState_star_mul_self_nonneg
#print axioms OPH.QFT.heisenbergFlow_kms_point
#print axioms OPH.QFT.gibbsState_kms
#print axioms OPH.QFT.witnessPrivateAlgebra_eq
#print axioms OPH.QFT.privateGibbsPartition_ne_zero
#print axioms OPH.QFT.privateGibbsState_one
#print axioms OPH.QFT.privateGibbsState_kms
#print axioms OPH.QFT.mem_totalSupportRegion
#print axioms OPH.QFT.heisenbergFlow_diagonal_apply
#print axioms OPH.QFT.supportGraded_flow_fixes_local
#print axioms OPH.QFT.supportGraded_flow_preserves
#print axioms OPH.QFT.supportGraded_kms_degenerate
#print axioms OPH.QFT.supportGraded_kms_degenerate_functional
#print axioms OPH.QFT.supportGradedNet_regionalExpectation_kms_degenerate
#print axioms OPH.QFT.witnessHamiltonian_mem
#print axioms OPH.QFT.witnessHamiltonian_isHermitian
#print axioms OPH.QFT.witnessHamiltonian_not_mem_empty
#print axioms OPH.QFT.witnessHamiltonian_flow_moves
#print axioms OPH.QFT.publicRecord_no_thermal_flow
#print axioms OPH.QFT.publicRecord_flow_fixes
