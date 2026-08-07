import Mathlib.Analysis.CStarAlgebra.Matrix
import Mathlib.Analysis.InnerProductSpace.Adjoint
import QFT.StructuralInheritance
import Dynamics.PrivateInner
import EventAlgebra.Superselection
import Tower.OperationalObserver

/-!
# The E4 structural-inheritance matrix: the DHR / superselection-sector row

This module delivers the second row of the E4 inheritance matrix of
completion-plan issue `#701`: Doplicher–Haag–Roberts sector structure
and gauge-group reconstruction.  Per the issue, the row is a separate
theorem target with explicit hypotheses, a nonvacuous example, and a
typed non-evaluable exit; no blanket QFT inheritance claim is made or
licensed by anything below.  The verdict shape of this row is
REFUTABLE, and that is the result, not a caveat: Doplicher–Roberts
reconstruction over a finite type-I net reconstructs the TRIVIAL gauge
group, so the constructed net inherits DHR structure only degenerately.

## Part 1 — PROVED: the transported inner-implementation engine

The engine is B3's `OPH.Dynamics.finitePrivateStarAutomorphism_inner`:
every star-algebra automorphism of the continuous-endomorphism algebra
of a finite-dimensional complex Hilbert space is implemented by a
unitary.  The committed carriers speak the `Matrix (Fin 2) (Fin 2) ℂ`
language, so the theorem is transported across the star-algebra
equivalence `Matrix.toEuclideanCLM` between square matrices and
continuous endomorphisms of Euclidean space
(`matrixStarAlgAutomorphism_inner`): every star-algebra automorphism
of a finite matrix algebra is conjugation by a unitary matrix, with
the implementing unitary produced through
`Unitary.linearIsometryEquiv` and the conjugation read back through
`Unitary.conjStarAlgAut`.

## Part 2 — REFUTABLE: sector triviality of the finite constructed net

`LocalizedIn ρ U` is the DHR-type localization predicate against the
E3 example net: the endomorphism `ρ` of the private algebra acts as
the identity on the regional algebra of every region declared disjoint
from `U` by `supportGradedNet.disjoint`, the net's own declared
disjointness.  The sector theorem
(`supportGraded_localized_inner`) states that every bijective
localized endomorphism is inner, hence sector-equivalent to the
identity (`supportGraded_localized_sector_trivial`): the DHR sector
category of THIS FINITE CONSTRUCTED NET has one object up to
equivalence, and the reconstructed gauge group is trivial.  The proof
consumes part 1 and discards the localization hypothesis — that
discard is the content, not a defect: over a finite type-I block every
automorphism whatsoever is already inner, so localization can create
no charge.  The row therefore refutes, for this net, the hope that
DHR superselection structure is inherited nontrivially.

Nonvacuity: the hypothesis class contains more than the identity.
`sectorWitness`, conjugation by the diagonal involution
`sectorFlip = diagonal (1, -1)`, is a bijective endomorphism localized
in the nonempty region `{0}` (`sectorWitness_localizedIn`) that
genuinely moves the off-diagonal matrix unit
(`sectorWitness_moves`, `sectorWitness_ne_id_on`), and the main
theorem applies to it (`sectorWitness_inner`).  The definitional
degeneracy at the empty region is stated, not hidden:
`localizedIn_empty` records that `LocalizedIn ρ ∅` holds vacuously for
every `ρ` because the net's declared disjointness requires both
regions nonempty; the witness therefore uses a nonempty region.

The row is tied to the finite superselection module: relative to the
declared two-sector partition `witnessPartition`, the nontrivial
localized witness is operationally invisible — `sectorWitness X` and
`X` have identical trace statistics against every matrix of the
sector-preserving commutant
(`sectorWitness_operationally_invisible`, via
`EventAlgebra.PartitionOperationallyEquivalent`).  The localized
endomorphism carries no charge visible to any sector-preserving
readout, which is the operational face of sector triviality.  It is
visible to sector-off-diagonal tests (it moves the matrix unit), so
the commutant restriction is load-bearing, not decorative.

## Part 3 — TYPED NON-EVALUABLE EXIT: full DHR reconstruction

Full DHR sector theory — localized transportable endomorphisms of the
quasi-local algebra relative to a vacuum representation, their
intertwiner categories, and Doplicher–Roberts reconstruction of a
compact gauge group — is not evaluable against the committed carriers,
and this module declines it rather than approximating it.  The missing
types are: no quasi-local (inductive-limit) observable algebra exists
— E3's own claim boundary (`LocallyCovariantLimit.lean`) records the
filtered-colimit algebra as a declared later deliverable, explicitly
not constructed; no Hilbert-space representation of the net, hence no
vacuum representation and no superselection criterion relative to one;
no transportability statement is even statable, since transporting a
localized morphism requires an unbounded supply of mutually disjoint
regions and the two-site example has exactly one declared disjoint
pair; and the pinned Mathlib has no representation theory of
C*-algebras and no sector category.  Until the limit algebra and a
represented net exist, full DHR is a typed exit of this row, exactly
as horizon thermality is a typed exit of the KMS row and time-slice is
a typed exit of E3.

No statement in this module may be quoted as "the observer net has
trivial superselection structure" or as any general DHR result.  What
is proved is: every star-automorphism of the finite private block is
unitarily inner (part 1), and FOR THIS FINITE CONSTRUCTED NET every
bijective localized endomorphism is inner, so the finite DHR-type
sector structure the net inherits is the trivial one (part 2).
-/

namespace OPH.QFT

open OPH.Tower

/-! ## Part 1 — the transported inner-implementation engine -/

set_option maxHeartbeats 800000 in
/-- **The transported engine.**  Every star-algebra automorphism of a
finite complex matrix algebra is conjugation by a unitary matrix.
This is B3's `OPH.Dynamics.finitePrivateStarAutomorphism_inner`
transported across `Matrix.toEuclideanCLM` from the
continuous-endomorphism spelling to the matrix spelling: the
automorphism is conjugated into an automorphism of
`EuclideanSpace ℂ n →L[ℂ] EuclideanSpace ℂ n`, implemented there by a
unitary via B3, and the implementing unitary is carried back through
the same star-algebra equivalence.

The heartbeat bump is for the instance-level defeq work of the
`EuclideanSpace` (`PiLp`) spelling; Mathlib's own bridge lemma
`Unitary.conjStarAlgAut_symm_unitaryLinearIsometryEquiv` carries a
`maxHeartbeats 400000` bump for the same reason. -/
theorem matrixStarAlgAutomorphism_inner {n : Type*} [Fintype n] [DecidableEq n]
    (F : Matrix n n ℂ ≃⋆ₐ[ℂ] Matrix n n ℂ) :
    ∃ u ∈ unitary (Matrix n n ℂ), ∀ x, F x = u * x * star u := by
  have e : Matrix n n ℂ ≃⋆ₐ[ℂ] (EuclideanSpace ℂ n →L[ℂ] EuclideanSpace ℂ n) :=
    Matrix.toEuclideanCLM
  obtain ⟨U, hU⟩ :=
    OPH.Dynamics.finitePrivateStarAutomorphism_inner ((e.symm.trans F).trans e)
  have hconj : (e.symm.trans F).trans e =
      Unitary.conjStarAlgAut ℂ (EuclideanSpace ℂ n →L[ℂ] EuclideanSpace ℂ n)
        (Unitary.linearIsometryEquiv.symm U) := by
    rw [hU]
    exact (Unitary.conjStarAlgAut_symm_unitaryLinearIsometryEquiv U).symm
  have hv1 := (Unitary.mem_iff.mp (Unitary.linearIsometryEquiv.symm U).2).1
  have hv2 := (Unitary.mem_iff.mp (Unitary.linearIsometryEquiv.symm U).2).2
  refine ⟨e.symm ↑(Unitary.linearIsometryEquiv.symm U),
    Unitary.mem_iff.mpr ⟨?_, ?_⟩, fun x => ?_⟩
  · rw [← map_star, ← map_mul, hv1, map_one]
  · rw [← map_star, ← map_mul, hv2, map_one]
  · have hx : F x = e.symm (((e.symm.trans F).trans e) (e x)) := by
      rw [StarAlgEquiv.trans_apply, StarAlgEquiv.trans_apply,
        StarAlgEquiv.symm_apply_apply, StarAlgEquiv.symm_apply_apply]
    rw [hx, hconj, Unitary.conjStarAlgAut_apply, map_mul, map_mul,
      StarAlgEquiv.symm_apply_apply, map_star]

/-- The transported engine restated on the witness tower's private
block, the carrier named by the E4 issue (definitionally
`Matrix (Fin 2) (Fin 2) ℂ`, receipt `witnessPrivateAlgebra_eq`). -/
theorem privateMatrixAutomorphism_inner
    (F : ConsensusTower.PrivateAlgebra witnessTower () ≃⋆ₐ[ℂ]
      ConsensusTower.PrivateAlgebra witnessTower ()) :
    ∃ u ∈ unitary (ConsensusTower.PrivateAlgebra witnessTower ()),
      ∀ x, F x = u * x * star u :=
  matrixStarAlgAutomorphism_inner F

/-! ## Part 2 — DHR-type localization against the E3 net -/

/-- DHR-type localization in the region `U`, relative to the E3
example net: the endomorphism acts as the identity on the regional
algebra of every region declared disjoint from `U` by the net's own
declared disjointness `supportGradedNet.disjoint`.  Note that the
declared disjointness requires both regions nonempty, so localization
in the empty region is vacuous (`localizedIn_empty`); the nonvacuity
witness below is localized in a nonempty region. -/
def LocalizedIn
    (ρ : Matrix (Fin 2) (Fin 2) ℂ →⋆ₐ[ℂ] Matrix (Fin 2) (Fin 2) ℂ)
    (U : Finset (Fin 2)) : Prop :=
  ∀ V : Finset (Fin 2), supportGradedNet.disjoint () U V →
    ∀ b : Matrix (Fin 2) (Fin 2) ℂ,
      b ∈ supportGradedNet.localAlgebra () V → ρ b = b

/-- Honesty receipt for the definitional boundary: because the net's
declared disjointness requires both regions nonempty, EVERY
endomorphism is (vacuously) localized in the empty region.  A sector
statement quantified over `LocalizedIn ρ ∅` would therefore be about
all endomorphisms; the nonvacuity witness uses `{0}` instead. -/
theorem localizedIn_empty
    (ρ : Matrix (Fin 2) (Fin 2) ℂ →⋆ₐ[ℂ] Matrix (Fin 2) (Fin 2) ℂ) :
    LocalizedIn ρ (∅ : Finset (Fin 2)) := by
  intro V hV b _
  have hV' : (∅ : Finset (Fin 2)).Nonempty ∧ V.Nonempty ∧
      Disjoint (∅ : Finset (Fin 2)) V := hV
  exact absurd hV'.1 Finset.not_nonempty_empty

/-- DHR sector equivalence of two endomorphisms of the private
algebra: one is a unitary conjugate of the other.  Sector triviality
is `SectorEquivalent ρ (StarAlgHom.id ℂ _)`. -/
def SectorEquivalent
    (ρ σ : Matrix (Fin 2) (Fin 2) ℂ →⋆ₐ[ℂ] Matrix (Fin 2) (Fin 2) ℂ) : Prop :=
  ∃ u ∈ unitary (Matrix (Fin 2) (Fin 2) ℂ), ∀ x, ρ x = u * σ x * star u

/-- **The sector theorem for the finite constructed net.**  Every
bijective endomorphism of the private algebra that is DHR-localized in
a region of the E3 example net is inner: it is conjugation by a
unitary matrix.  The proof consumes the transported engine and
DISCARDS the localization hypothesis — over a finite type-I block
every automorphism whatsoever is inner, which is exactly why
Doplicher–Roberts reconstruction over this net returns the trivial
gauge group.  This is a statement about THIS finite constructed net;
no general DHR claim is made. -/
theorem supportGraded_localized_inner
    (ρ : Matrix (Fin 2) (Fin 2) ℂ →⋆ₐ[ℂ] Matrix (Fin 2) (Fin 2) ℂ)
    {U : Finset (Fin 2)} (_hloc : LocalizedIn ρ U)
    (hbij : Function.Bijective ρ) :
    ∃ u ∈ unitary (Matrix (Fin 2) (Fin 2) ℂ), ∀ x, ρ x = u * x * star u := by
  obtain ⟨u, hu, hconj⟩ :=
    matrixStarAlgAutomorphism_inner (StarAlgEquiv.ofBijective ρ hbij)
  exact ⟨u, hu, fun x =>
    ((StarAlgEquiv.ofBijective_apply hbij x).symm).trans (hconj x)⟩

/-- Sector triviality, stated as sector equivalence: every bijective
localized endomorphism of the finite constructed net is
sector-equivalent to the identity endomorphism.  The DHR-type sector
category of this net has a single sector. -/
theorem supportGraded_localized_sector_trivial
    (ρ : Matrix (Fin 2) (Fin 2) ℂ →⋆ₐ[ℂ] Matrix (Fin 2) (Fin 2) ℂ)
    {U : Finset (Fin 2)} (hloc : LocalizedIn ρ U)
    (hbij : Function.Bijective ρ) :
    SectorEquivalent ρ (StarAlgHom.id ℂ (Matrix (Fin 2) (Fin 2) ℂ)) := by
  obtain ⟨u, hu, h⟩ := supportGraded_localized_inner ρ hloc hbij
  exact ⟨u, hu, fun x => h x⟩

/-- The sector theorem restated with the endomorphism typed on the
witness tower's private block, the carrier named by the E4 issue. -/
theorem privateSector_localized_inner
    (ρ : ConsensusTower.PrivateAlgebra witnessTower () →⋆ₐ[ℂ]
      ConsensusTower.PrivateAlgebra witnessTower ())
    {U : Finset (Fin 2)} (hloc : LocalizedIn ρ U)
    (hbij : Function.Bijective ρ) :
    ∃ u ∈ unitary (ConsensusTower.PrivateAlgebra witnessTower ()),
      ∀ x, ρ x = u * x * star u :=
  supportGraded_localized_inner ρ hloc hbij

/-! ### Nonvacuity: a concrete nontrivial localized endomorphism -/

/-- The diagonal involution `diag(1, -1)`, the implementing unitary of
the nonvacuity witness. -/
noncomputable def sectorFlip : Matrix (Fin 2) (Fin 2) ℂ :=
  Matrix.diagonal ![1, -1]

theorem sectorFlip_mul_self : sectorFlip * sectorFlip = 1 := by
  show Matrix.diagonal ![1, -1] * Matrix.diagonal ![1, -1] = 1
  rw [Matrix.diagonal_mul_diagonal, ← Matrix.diagonal_one]
  congr 1
  funext i
  fin_cases i <;> simp

theorem sectorFlip_star : star sectorFlip = sectorFlip := by
  show star (Matrix.diagonal ![1, -1]) = Matrix.diagonal ![1, -1]
  rw [Matrix.star_eq_conjTranspose, Matrix.diagonal_conjTranspose]
  congr 1
  funext i
  fin_cases i <;> simp

theorem sectorFlip_mem_unitary :
    sectorFlip ∈ unitary (Matrix (Fin 2) (Fin 2) ℂ) :=
  Unitary.mem_iff.mpr
    ⟨by rw [sectorFlip_star, sectorFlip_mul_self],
     by rw [sectorFlip_star, sectorFlip_mul_self]⟩

/-- The flip commutes with every diagonal matrix; proved in the plain
matrix language via the E3 locality receipt `diagonal_commute`. -/
theorem sectorFlip_comm_diagonal (f : Fin 2 → ℂ) :
    sectorFlip * Matrix.diagonal f = Matrix.diagonal f * sectorFlip :=
  diagonal_commute ![1, -1] f

/-- The nonvacuity witness: conjugation by the diagonal involution
`sectorFlip`, bundled as a star-algebra endomorphism of the private
algebra. -/
noncomputable def sectorWitness :
    Matrix (Fin 2) (Fin 2) ℂ →⋆ₐ[ℂ] Matrix (Fin 2) (Fin 2) ℂ where
  toFun x := sectorFlip * x * sectorFlip
  map_one' := by rw [mul_one, sectorFlip_mul_self]
  map_mul' x y := by
    calc sectorFlip * (x * y) * sectorFlip
        = sectorFlip * x * (sectorFlip * sectorFlip) * y * sectorFlip := by
          rw [sectorFlip_mul_self]; noncomm_ring
      _ = sectorFlip * x * sectorFlip * (sectorFlip * y * sectorFlip) := by
          noncomm_ring
  map_zero' := by rw [mul_zero, zero_mul]
  map_add' x y := by rw [mul_add, add_mul]
  commutes' c := by
    rw [Algebra.algebraMap_eq_smul_one, mul_smul_comm, smul_mul_assoc,
      mul_one, sectorFlip_mul_self]
  map_star' x := by
    rw [star_mul, star_mul, sectorFlip_star, mul_assoc]

theorem sectorWitness_apply (x : Matrix (Fin 2) (Fin 2) ℂ) :
    sectorWitness x = sectorFlip * x * sectorFlip :=
  rfl

/-- The witness fixes every diagonal matrix, hence every element of
every regional algebra of the support-graded net. -/
theorem sectorWitness_fixes_diagonal (f : Fin 2 → ℂ) :
    sectorWitness (Matrix.diagonal f) = Matrix.diagonal f := by
  rw [sectorWitness_apply, sectorFlip_comm_diagonal, mul_assoc,
    sectorFlip_mul_self, mul_one]

/-- The witness is localized in the nonempty region `{0}`: it acts as
the identity on the regional algebra of every region declared disjoint
from `{0}` (indeed on every regional algebra, since all of them are
diagonal). -/
theorem sectorWitness_localizedIn :
    LocalizedIn sectorWitness ({0} : Finset (Fin 2)) := by
  intro V _ b hb
  have hb' : b ∈ supportDiagonalAlgebra V := hb
  obtain ⟨f, rfl, -⟩ := hb'
  exact sectorWitness_fixes_diagonal f

theorem sectorWitness_involutive : Function.Involutive sectorWitness := by
  intro x
  rw [sectorWitness_apply, sectorWitness_apply]
  calc sectorFlip * (sectorFlip * x * sectorFlip) * sectorFlip
      = sectorFlip * sectorFlip * x * (sectorFlip * sectorFlip) := by
        noncomm_ring
    _ = x := by rw [sectorFlip_mul_self, one_mul, mul_one]

theorem sectorWitness_bijective : Function.Bijective sectorWitness :=
  sectorWitness_involutive.bijective

/-- The witness genuinely moves the off-diagonal matrix unit: it
negates it.  The hypothesis class of the sector theorem contains more
than the identity endomorphism. -/
theorem sectorWitness_moves :
    sectorWitness (Matrix.single 0 1 (1 : ℂ)) =
      -Matrix.single 0 1 (1 : ℂ) := by
  rw [sectorWitness_apply]
  show Matrix.diagonal ![1, -1] * Matrix.single 0 1 (1 : ℂ) *
      Matrix.diagonal ![1, -1] = -Matrix.single 0 1 (1 : ℂ)
  ext i j
  rw [Matrix.mul_diagonal, Matrix.diagonal_mul, Matrix.neg_apply]
  fin_cases i <;> fin_cases j <;> simp

theorem sectorWitness_ne_id_on :
    sectorWitness (Matrix.single 0 1 (1 : ℂ)) ≠
      Matrix.single 0 1 (1 : ℂ) := by
  rw [sectorWitness_moves]
  intro heq
  have h := congrArg (fun M => M 0 1) heq
  simp only [Matrix.neg_apply, Matrix.single_apply_same] at h
  norm_num at h

/-- The sector theorem applied to the nonvacuity witness: the
nontrivial localized endomorphism is inner, as the theorem asserts of
every member of its (now demonstrably inhabited) hypothesis class. -/
theorem sectorWitness_inner :
    ∃ u ∈ unitary (Matrix (Fin 2) (Fin 2) ℂ),
      ∀ x, sectorWitness x = u * x * star u :=
  supportGraded_localized_inner sectorWitness sectorWitness_localizedIn
    sectorWitness_bijective

/-! ### The superselection cross-receipt

Relative to the declared two-sector partition `witnessPartition`, the
nontrivial localized witness is operationally invisible: it changes no
trace statistic against any matrix of the sector-preserving
commutant.  This ties the sector triviality of this row to the finite
superselection module: the localized endomorphism carries no charge
any sector-preserving readout can see.  It is NOT invisible to
sector-off-diagonal tests — it negates the off-diagonal matrix unit —
so the commutant restriction in
`EventAlgebra.PartitionOperationallyEquivalent` is load-bearing. -/

/-- For the rank-one diagonal partition, pinching is exactly the
diagonal part. -/
theorem witnessPartition_pinching_eq (X : Matrix (Fin 2) (Fin 2) ℂ) :
    EventAlgebra.partitionPinching witnessPartition X =
      Matrix.diagonal fun i => X i i := by
  have hproj : ∀ i : Fin 2, witnessPartition.proj i =
      Matrix.diagonal (Pi.single i (1 : ℂ)) := fun _ => rfl
  show ∑ i, witnessPartition.proj i * X * witnessPartition.proj i =
    Matrix.diagonal fun i => X i i
  rw [Fin.sum_univ_two, hproj 0, hproj 1]
  ext i j
  rw [Matrix.add_apply, Matrix.mul_diagonal, Matrix.diagonal_mul,
    Matrix.mul_diagonal, Matrix.diagonal_mul]
  by_cases hij : i = j
  · subst hij
    rw [Matrix.diagonal_apply_eq]
    fin_cases i <;> simp
  · rw [Matrix.diagonal_apply_ne _ hij]
    fin_cases i <;> fin_cases j <;> simp_all

/-- The witness preserves every diagonal entry. -/
theorem sectorWitness_diag_apply (X : Matrix (Fin 2) (Fin 2) ℂ) (i : Fin 2) :
    sectorWitness X i i = X i i := by
  have h : sectorWitness X =
      Matrix.diagonal ![1, -1] * X * Matrix.diagonal ![1, -1] := rfl
  rw [h, Matrix.mul_diagonal, Matrix.diagonal_mul]
  fin_cases i <;> simp

/-- **The charge-invisibility receipt.**  The nontrivial localized
witness is operationally equivalent to the identity's action on every
input, relative to the declared two-sector partition: no
sector-preserving trace readout distinguishes `sectorWitness X` from
`X`. -/
theorem sectorWitness_operationally_invisible (X : Matrix (Fin 2) (Fin 2) ℂ) :
    EventAlgebra.PartitionOperationallyEquivalent witnessPartition
      (sectorWitness X) X := by
  rw [EventAlgebra.partitionOperationallyEquivalent_iff_pinching_eq,
    witnessPartition_pinching_eq, witnessPartition_pinching_eq]
  congr 1
  funext i
  exact sectorWitness_diag_apply X i

end OPH.QFT

-- Axiom audit: no project-specific axiom or admission is permitted here.
#print axioms OPH.QFT.matrixStarAlgAutomorphism_inner
#print axioms OPH.QFT.privateMatrixAutomorphism_inner
#print axioms OPH.QFT.LocalizedIn
#print axioms OPH.QFT.localizedIn_empty
#print axioms OPH.QFT.SectorEquivalent
#print axioms OPH.QFT.supportGraded_localized_inner
#print axioms OPH.QFT.supportGraded_localized_sector_trivial
#print axioms OPH.QFT.privateSector_localized_inner
#print axioms OPH.QFT.sectorFlip
#print axioms OPH.QFT.sectorFlip_mul_self
#print axioms OPH.QFT.sectorFlip_star
#print axioms OPH.QFT.sectorFlip_mem_unitary
#print axioms OPH.QFT.sectorFlip_comm_diagonal
#print axioms OPH.QFT.sectorWitness
#print axioms OPH.QFT.sectorWitness_apply
#print axioms OPH.QFT.sectorWitness_fixes_diagonal
#print axioms OPH.QFT.sectorWitness_localizedIn
#print axioms OPH.QFT.sectorWitness_involutive
#print axioms OPH.QFT.sectorWitness_bijective
#print axioms OPH.QFT.sectorWitness_moves
#print axioms OPH.QFT.sectorWitness_ne_id_on
#print axioms OPH.QFT.sectorWitness_inner
#print axioms OPH.QFT.witnessPartition_pinching_eq
#print axioms OPH.QFT.sectorWitness_diag_apply
#print axioms OPH.QFT.sectorWitness_operationally_invisible
