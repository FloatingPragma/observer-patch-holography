import Mathlib.Analysis.Matrix.Spectrum
import Mathlib.LinearAlgebra.Lagrange
import Mathlib.LinearAlgebra.Dimension.Finite
import Mathlib.LinearAlgebra.FiniteDimensional.Lemmas
import Mathlib.RingTheory.Adjoin.Polynomial.Basic
import Dynamics.CentralBlocks

/-!
# Spectral decomposition of the center of a star subalgebra

`CentralBlocks.lean` controls the two-block direct sum.  This module carries
step 1 of the remaining gap of issue 679: the center of a star subalgebra of
a finite complex matrix algebra is spanned by finitely many minimal
orthogonal central projections summing to one.

1. Relative-center vocabulary: membership in the center of a star
   subalgebra `S` (`IsCentralIn`), central star projections of `S`
   (`IsCentralProjectionIn`), and minimal central projections of `S`
   (`IsMinimalCentralProjectionIn`).  The center is a submodule
   (`centerSubmodule`) and each central projection cuts a corner submodule
   of the center (`centerCorner`).
2. Finite descent: every central projection of `S` splits into a finite
   orthogonal family of minimal central projections of `S`
   (`IsCentralProjectionIn.exists_minimal_decomposition`), by strict
   descent of the corner dimension.  Applied to the unit this yields a
   minimal orthogonal resolution of the identity
   (`exists_minimal_central_resolution`).  This part is finite linear
   algebra over any finite-dimensional complex star algebra.
3. Spectral projections: on `Matrix n n ℂ`, a self-adjoint element of the
   center of `S` has a family of spectral projections obtained by Lagrange
   polynomial functional calculus on its eigenvalue support
   (`spectralProjection`).  They are central star projections of `S`,
   pairwise orthogonal, nonzero, sum to one, and reconstruct the element
   (`IsCentralIn.exists_spectral_resolution`).  A non-scalar self-adjoint
   central element therefore splits off a nontrivial central projection
   (`IsCentralIn.exists_nontrivial_central_projection`).
4. Spanning: against a minimal central projection every self-adjoint
   central element acts as a scalar (`IsCentralIn.mul_minimal_eq_smul`),
   so the center is exactly the span of any minimal orthogonal resolution
   of the identity (`exists_minimal_central_resolution_spanning`).  The
   minimal resolution is unique up to permutation
   (`minimal_central_resolution_unique`).
5. Connection to the central cut: a two-projection resolution of the
   identity by central star projections factors the algebra through
   `centralCutEquiv` (`pairResolutionCutEquiv`); the second projection is
   forced to be the complement (`pair_resolution_eq_one_sub`).

The iterated `k`-fold product equivalence over a full minimal resolution
and the induction over block counts remain separate obligations of issue
679.  No physical clock, source-selected Hamiltonian, or physical time
identification is claimed.
-/

namespace OPH.Dynamics

universe u

open Polynomial

/-! ## The relative center of a star subalgebra -/

section RelativeCenter

variable {A : Type u} [Ring A] [StarRing A] [Algebra ℂ A] [StarModule ℂ A]

/-- Membership in the center of a star subalgebra: an element of the
subalgebra commuting with every element of the subalgebra. -/
def IsCentralIn (S : StarSubalgebra ℂ A) (z : A) : Prop :=
  z ∈ S ∧ ∀ x ∈ S, z * x = x * z

namespace IsCentralIn

variable {S : StarSubalgebra ℂ A} {z w : A}

theorem mem (hz : IsCentralIn S z) : z ∈ S := hz.1

theorem commutes (hz : IsCentralIn S z) {x : A} (hx : x ∈ S) :
    z * x = x * z := hz.2 x hx

protected theorem one : IsCentralIn S 1 :=
  ⟨one_mem S, fun x _ => by rw [one_mul, mul_one]⟩

protected theorem zero : IsCentralIn S 0 :=
  ⟨zero_mem S, fun x _ => by rw [zero_mul, mul_zero]⟩

protected theorem add (hz : IsCentralIn S z) (hw : IsCentralIn S w) :
    IsCentralIn S (z + w) :=
  ⟨add_mem hz.mem hw.mem, fun x hx => by
    rw [add_mul, mul_add, hz.commutes hx, hw.commutes hx]⟩

protected theorem sub (hz : IsCentralIn S z) (hw : IsCentralIn S w) :
    IsCentralIn S (z - w) :=
  ⟨sub_mem hz.mem hw.mem, fun x hx => by
    rw [sub_mul, mul_sub, hz.commutes hx, hw.commutes hx]⟩

protected theorem smul (c : ℂ) (hz : IsCentralIn S z) :
    IsCentralIn S (c • z) :=
  ⟨SMulMemClass.smul_mem c hz.mem, fun x hx => by
    rw [smul_mul_assoc, mul_smul_comm, hz.commutes hx]⟩

protected theorem mul (hz : IsCentralIn S z) (hw : IsCentralIn S w) :
    IsCentralIn S (z * w) :=
  ⟨mul_mem hz.mem hw.mem, fun x hx => by
    rw [mul_assoc, hw.commutes hx, ← mul_assoc, hz.commutes hx, mul_assoc]⟩

protected theorem star (hz : IsCentralIn S z) : IsCentralIn S (star z) := by
  refine ⟨star_mem hz.mem, fun x hx => ?_⟩
  calc star z * x = star (star x * z) := by rw [star_mul, star_star]
    _ = star (z * star x) := by rw [hz.commutes (star_mem hx)]
    _ = x * star z := by rw [star_mul, star_star]

/-- Two members of the center commute with each other. -/
theorem commute (hz : IsCentralIn S z) (hw : IsCentralIn S w) :
    z * w = w * z := hz.commutes hw.mem

end IsCentralIn

/-- A central star projection of a star subalgebra: a self-adjoint
idempotent lying in the subalgebra and commuting with all of it. -/
structure IsCentralProjectionIn (S : StarSubalgebra ℂ A) (p : A) : Prop
    extends IsStarProjection p where
  centralIn : IsCentralIn S p

/-- A minimal central projection of a star subalgebra: a nonzero central
star projection whose only central subprojections are zero and itself. -/
structure IsMinimalCentralProjectionIn (S : StarSubalgebra ℂ A) (p : A) :
    Prop extends IsCentralProjectionIn S p where
  ne_zero : p ≠ 0
  minimal : ∀ q, IsCentralProjectionIn S q → q * p = q → q = 0 ∨ q = p

variable {S : StarSubalgebra ℂ A}

/-- The unit is a central star projection of every star subalgebra. -/
theorem isCentralProjectionIn_one : IsCentralProjectionIn S 1 where
  isIdempotentElem := IsIdempotentElem.one
  isSelfAdjoint := IsSelfAdjoint.one A
  centralIn := IsCentralIn.one

/-- The product of two central star projections is a central star
projection. -/
theorem IsCentralProjectionIn.mul {p q : A} (hp : IsCentralProjectionIn S p)
    (hq : IsCentralProjectionIn S q) : IsCentralProjectionIn S (p * q) where
  isIdempotentElem := by
    have hcomm : q * p = p * q := hq.centralIn.commute hp.centralIn
    show p * q * (p * q) = p * q
    rw [mul_assoc, ← mul_assoc q p q, hcomm, mul_assoc p q q,
      hq.isIdempotentElem.eq, ← mul_assoc, hp.isIdempotentElem.eq]
  isSelfAdjoint := by
    show star (p * q) = p * q
    rw [star_mul, hp.isSelfAdjoint.star_eq, hq.isSelfAdjoint.star_eq]
    exact hq.centralIn.commute hp.centralIn
  centralIn := hp.centralIn.mul hq.centralIn

/-- The difference of a central star projection and a central
subprojection is a central star projection. -/
theorem IsCentralProjectionIn.sub {p q : A} (hp : IsCentralProjectionIn S p)
    (hq : IsCentralProjectionIn S q) (hle : q * p = q) :
    IsCentralProjectionIn S (p - q) where
  isIdempotentElem := by
    have hpq : p * q = q := by
      rw [← hq.centralIn.commutes hp.centralIn.mem]; exact hle
    show (p - q) * (p - q) = p - q
    rw [sub_mul, mul_sub, mul_sub, hp.isIdempotentElem.eq, hpq, hle,
      hq.isIdempotentElem.eq]
    abel
  isSelfAdjoint := hp.isSelfAdjoint.sub hq.isSelfAdjoint
  centralIn := hp.centralIn.sub hq.centralIn

/-- The center of a star subalgebra as a complex submodule of the ambient
algebra. -/
def centerSubmodule (S : StarSubalgebra ℂ A) : Submodule ℂ A where
  carrier := {z : A | IsCentralIn S z}
  add_mem' := fun ha hb => ha.add hb
  zero_mem' := IsCentralIn.zero
  smul_mem' := fun c _ ha => ha.smul c

@[simp] theorem mem_centerSubmodule_iff {z : A} :
    z ∈ centerSubmodule S ↔ IsCentralIn S z := Iff.rfl

/-- The corner of the center cut by an element: the image of the center
under left multiplication. -/
def centerCorner (S : StarSubalgebra ℂ A) (p : A) : Submodule ℂ A :=
  Submodule.map (LinearMap.mulLeft ℂ p) (centerSubmodule S)

theorem mem_centerCorner_iff' {p x : A} :
    x ∈ centerCorner S p ↔ ∃ z : A, IsCentralIn S z ∧ p * z = x := by
  constructor
  · rintro ⟨z, hz, hx⟩
    exact ⟨z, hz, hx⟩
  · rintro ⟨z, hz, hx⟩
    exact ⟨z, hz, hx⟩

/-- Every element lies in its own corner of the center, cut against the
central unit. -/
theorem self_mem_centerCorner {p : A} : p ∈ centerCorner S p :=
  mem_centerCorner_iff'.mpr ⟨1, IsCentralIn.one, mul_one p⟩

/-- The corner of a central subprojection is contained in the corner of
the projection. -/
theorem centerCorner_le {p q : A} (hp : IsCentralProjectionIn S p)
    (hq : IsCentralProjectionIn S q) (hle : q * p = q) :
    centerCorner S q ≤ centerCorner S p := by
  intro x hx
  obtain ⟨z, hz, hxz⟩ := mem_centerCorner_iff'.mp hx
  have hpq : p * q = q := by
    rw [← hq.centralIn.commutes hp.centralIn.mem]; exact hle
  refine mem_centerCorner_iff'.mpr ⟨q * z, hq.centralIn.mul hz, ?_⟩
  rw [← mul_assoc, hpq, hxz]

/-- A central projection is outside the corner of any strictly smaller
central subprojection. -/
theorem notMem_centerCorner {p q : A} (hq : IsCentralProjectionIn S q)
    (hle : q * p = q) (hne : q ≠ p) : p ∉ centerCorner S q := by
  intro hmem
  obtain ⟨z, _, hz⟩ := mem_centerCorner_iff'.mp hmem
  have h1 : q * p = p := by
    rw [← hz, ← mul_assoc, hq.isIdempotentElem.eq, hz]
  exact hne (hle.symm.trans h1)

/-- Strict descent of corners: a proper nonzero central subprojection has
a strictly smaller corner of the center. -/
theorem centerCorner_lt {p q : A} (hp : IsCentralProjectionIn S p)
    (hq : IsCentralProjectionIn S q) (hle : q * p = q) (hne : q ≠ p) :
    centerCorner S q < centerCorner S p := by
  refine lt_of_le_of_ne (centerCorner_le hp hq hle) (fun h => ?_)
  exact notMem_centerCorner hq hle hne (h ▸ self_mem_centerCorner)

end RelativeCenter

/-! ## Existence of a minimal decomposition by finite descent -/

section Descent

variable {A : Type u} [Ring A] [StarRing A] [Algebra ℂ A] [StarModule ℂ A]
  {S : StarSubalgebra ℂ A}

/-- Every central star projection of a star subalgebra of a
finite-dimensional complex star algebra decomposes into a finite
orthogonal family of minimal central projections summing to it.  The
proof is a strict descent on the dimension of the corner of the center
cut by the projection. -/
theorem IsCentralProjectionIn.exists_minimal_decomposition
    [FiniteDimensional ℂ A] {p : A} (hp : IsCentralProjectionIn S p) :
    ∃ l : List A, (∀ e ∈ l, IsMinimalCentralProjectionIn S e) ∧
      l.Pairwise (fun a b => a * b = 0) ∧ (∀ e ∈ l, e * p = e) ∧
      l.sum = p := by
  suffices H : ∀ N : ℕ, ∀ p : A, IsCentralProjectionIn S p →
      Module.finrank ℂ (centerCorner S p) ≤ N →
      ∃ l : List A, (∀ e ∈ l, IsMinimalCentralProjectionIn S e) ∧
        l.Pairwise (fun a b => a * b = 0) ∧ (∀ e ∈ l, e * p = e) ∧
        l.sum = p by
    exact H (Module.finrank ℂ (centerCorner S p)) p hp le_rfl
  intro N
  induction N with
  | zero =>
    intro p hp hrank
    have hbot : centerCorner S p = ⊥ :=
      Submodule.finrank_eq_zero.mp (Nat.le_zero.mp hrank)
    have hp0 : p = 0 := by
      have := self_mem_centerCorner (S := S) (p := p)
      rw [hbot, Submodule.mem_bot] at this
      exact this
    exact ⟨[], by simp, by simp, by simp, by rw [List.sum_nil, hp0]⟩
  | succ N IH =>
    intro p hp hrank
    by_cases hp0 : p = 0
    · exact ⟨[], by simp, by simp, by simp, by rw [List.sum_nil, hp0]⟩
    by_cases hmin : ∀ q, IsCentralProjectionIn S q → q * p = q →
        q = 0 ∨ q = p
    · refine ⟨[p], ?_, ?_, ?_, ?_⟩
      · intro e he
        rw [List.mem_singleton] at he
        subst he
        exact ⟨hp, hp0, hmin⟩
      · exact List.pairwise_singleton _ p
      · intro e he
        rw [List.mem_singleton] at he
        subst he
        exact hp.isIdempotentElem.eq
      · exact List.sum_singleton
    · push Not at hmin
      obtain ⟨q, hq, hqle, hq0, hqne⟩ := hmin
      have hpq : IsCentralProjectionIn S (p - q) := hp.sub hq hqle
      have hsuble : (p - q) * p = p - q := by
        rw [sub_mul, hp.isIdempotentElem.eq, hqle]
      have hsubne : p - q ≠ p := by
        intro h1
        apply hq0
        have h2 : p - (p - q) = p - p := by rw [h1]
        simpa using h2
      have hqrank : Module.finrank ℂ (centerCorner S q) ≤ N :=
        Nat.lt_succ_iff.mp (lt_of_lt_of_le
          (Submodule.finrank_lt_finrank_of_lt
            (centerCorner_lt hp hq hqle hqne)) hrank)
      have hsubrank : Module.finrank ℂ (centerCorner S (p - q)) ≤ N :=
        Nat.lt_succ_iff.mp (lt_of_lt_of_le
          (Submodule.finrank_lt_finrank_of_lt
            (centerCorner_lt hp hpq hsuble hsubne)) hrank)
      obtain ⟨l₁, h₁min, h₁orth, h₁le, h₁sum⟩ := IH q hq hqrank
      obtain ⟨l₂, h₂min, h₂orth, h₂le, h₂sum⟩ := IH (p - q) hpq hsubrank
      have hqsub : q * (p - q) = 0 := by
        rw [mul_sub, hqle, hq.isIdempotentElem.eq, sub_self]
      have hcross : ∀ a ∈ l₁, ∀ b ∈ l₂, a * b = 0 := by
        intro a ha b hb
        have hb' : (p - q) * b = b := by
          rw [← (h₂min b hb).centralIn.commutes hpq.centralIn.mem]
          exact h₂le b hb
        calc a * b = (a * q) * ((p - q) * b) := by rw [h₁le a ha, hb']
          _ = a * (q * ((p - q) * b)) := by rw [mul_assoc]
          _ = a * ((q * (p - q)) * b) := by rw [← mul_assoc q _ b]
          _ = 0 := by rw [hqsub, zero_mul, mul_zero]
      refine ⟨l₁ ++ l₂, ?_, ?_, ?_, ?_⟩
      · intro e he
        rcases List.mem_append.mp he with h1 | h2
        · exact h₁min e h1
        · exact h₂min e h2
      · exact List.pairwise_append.mpr ⟨h₁orth, h₂orth, hcross⟩
      · intro e he
        rcases List.mem_append.mp he with h1 | h2
        · rw [← h₁le e h1, mul_assoc, hqle]
        · rw [← h₂le e h2, mul_assoc, hsuble]
      · rw [List.sum_append, h₁sum, h₂sum, add_comm, sub_add_cancel]

/-- Minimal orthogonal resolution of the identity: the unit of a star
subalgebra of a finite-dimensional complex star algebra splits into a
finite orthogonal family of minimal central projections. -/
theorem exists_minimal_central_resolution [FiniteDimensional ℂ A]
    (S : StarSubalgebra ℂ A) :
    ∃ l : List A, (∀ e ∈ l, IsMinimalCentralProjectionIn S e) ∧
      l.Pairwise (fun a b => a * b = 0) ∧ l.sum = 1 := by
  obtain ⟨l, hmin, horth, _, hsum⟩ :=
    (isCentralProjectionIn_one (S := S)).exists_minimal_decomposition
  exact ⟨l, hmin, horth, hsum⟩

end Descent

/-! ## Spectral projections of a central self-adjoint matrix -/

section Spectral

variable {n : Type u} [Fintype n] [DecidableEq n]

/-- The eigenvalue support of a Hermitian matrix inside `ℂ`. -/
noncomputable def spectralSupport {A : Matrix n n ℂ} (hA : A.IsHermitian) :
    Finset ℂ :=
  Finset.image (fun i => (RCLike.ofReal (hA.eigenvalues i) : ℂ)) Finset.univ

/-- The spectral projection of a Hermitian matrix at a point of its
eigenvalue support, realized by Lagrange polynomial functional calculus:
the Lagrange basis polynomial of the support evaluated on the matrix. -/
noncomputable def spectralProjection {A : Matrix n n ℂ} (hA : A.IsHermitian)
    (a : ℂ) : Matrix n n ℂ :=
  Polynomial.aeval A (Lagrange.basis (spectralSupport hA) id a)

private theorem aeval_conjStarAlgAut (u : unitary (Matrix n n ℂ))
    (D : Matrix n n ℂ) (p : Polynomial ℂ) :
    Polynomial.aeval (Unitary.conjStarAlgAut ℂ (Matrix n n ℂ) u D) p =
      Unitary.conjStarAlgAut ℂ (Matrix n n ℂ) u (Polynomial.aeval D p) := by
  have h1 := Polynomial.aeval_algHom_apply
    (Unitary.conjStarAlgAut ℂ (Matrix n n ℂ) u).toAlgEquiv D p
  simpa using h1

private theorem aeval_diagonal (d : n → ℂ) (p : Polynomial ℂ) :
    Polynomial.aeval (Matrix.diagonal d) p =
      Matrix.diagonal (fun i => p.eval (d i)) := by
  have h1 := Polynomial.aeval_algHom_apply (Matrix.diagonalAlgHom ℂ) d p
  have h2 : Matrix.diagonalAlgHom ℂ d = Matrix.diagonal d := rfl
  have h3 : Matrix.diagonalAlgHom ℂ (Polynomial.aeval d p) =
      Matrix.diagonal (Polynomial.aeval d p) := rfl
  rw [h2, h3] at h1
  rw [h1]
  congr 1
  funext i
  rw [Polynomial.aeval_fn_apply]
  exact congrFun (Polynomial.coe_aeval_eq_eval (d i)) p

/-- Polynomial functional calculus through the spectral theorem: the
polynomial acts on the diagonal of eigenvalues inside the eigenvector
frame. -/
private theorem aeval_spectral {A : Matrix n n ℂ} (hA : A.IsHermitian)
    (p : Polynomial ℂ) :
    Polynomial.aeval A p =
      (hA.eigenvectorUnitary : Matrix n n ℂ) *
        Matrix.diagonal
          (fun i => p.eval ((RCLike.ofReal (hA.eigenvalues i) : ℂ))) *
        (star hA.eigenvectorUnitary : Matrix n n ℂ) := by
  conv_lhs => rw [hA.spectral_theorem]
  rw [aeval_conjStarAlgAut, aeval_diagonal, Unitary.conjStarAlgAut_apply]
  rfl

/-- Sandwiching by the eigenvector unitary is injective. -/
private theorem unitary_sandwich_cancel (u : unitary (Matrix n n ℂ))
    (X : Matrix n n ℂ) :
    (star u : Matrix n n ℂ) *
      ((u : Matrix n n ℂ) * X * (star u : Matrix n n ℂ)) *
      (u : Matrix n n ℂ) = X := by
  have h1 : (star u : Matrix n n ℂ) * (u : Matrix n n ℂ) = 1 :=
    Unitary.coe_star_mul_self u
  calc (star u : Matrix n n ℂ) *
      ((u : Matrix n n ℂ) * X * (star u : Matrix n n ℂ)) *
      (u : Matrix n n ℂ)
      = ((star u : Matrix n n ℂ) * (u : Matrix n n ℂ)) * X *
        ((star u : Matrix n n ℂ) * (u : Matrix n n ℂ)) := by
        simp only [mul_assoc]
    _ = X := by rw [h1, one_mul, mul_one]

private theorem basis_eval_eig {A : Matrix n n ℂ} (hA : A.IsHermitian)
    {a : ℂ} (ha : a ∈ spectralSupport hA) (i : n) :
    (Lagrange.basis (spectralSupport hA) id a).eval
        ((RCLike.ofReal (hA.eigenvalues i) : ℂ)) =
      if (RCLike.ofReal (hA.eigenvalues i) : ℂ) = a then 1 else 0 := by
  by_cases hcase : (RCLike.ofReal (hA.eigenvalues i) : ℂ) = a
  · rw [if_pos hcase, hcase]
    simpa using Lagrange.eval_basis_self (Set.injOn_id _) ha
  · rw [if_neg hcase]
    have hmem : (RCLike.ofReal (hA.eigenvalues i) : ℂ) ∈ spectralSupport hA :=
      Finset.mem_image_of_mem _ (Finset.mem_univ i)
    simpa using Lagrange.eval_basis_of_ne (s := spectralSupport hA)
      (v := id) (i := a) (j := (RCLike.ofReal (hA.eigenvalues i) : ℂ))
      (fun h => hcase (h.symm)) hmem

/-- Two polynomials agreeing on the eigenvalues act identically. -/
private theorem aeval_eq_of_eval_eq {A : Matrix n n ℂ} (hA : A.IsHermitian)
    {p q : Polynomial ℂ}
    (hpq : ∀ i : n, p.eval ((RCLike.ofReal (hA.eigenvalues i) : ℂ)) =
      q.eval ((RCLike.ofReal (hA.eigenvalues i) : ℂ))) :
    Polynomial.aeval A p = Polynomial.aeval A q := by
  rw [aeval_spectral hA p, aeval_spectral hA q]
  have hd : (fun i => p.eval ((RCLike.ofReal (hA.eigenvalues i) : ℂ))) =
      (fun i => q.eval ((RCLike.ofReal (hA.eigenvalues i) : ℂ))) :=
    funext hpq
  rw [hd]

/-- Spectral projections at distinct support points are orthogonal. -/
theorem spectralProjection_mul_of_ne {A : Matrix n n ℂ}
    (hA : A.IsHermitian) {a b : ℂ} (ha : a ∈ spectralSupport hA)
    (hb : b ∈ spectralSupport hA) (hab : a ≠ b) :
    spectralProjection hA a * spectralProjection hA b = 0 := by
  rw [spectralProjection, spectralProjection, ← map_mul]
  have h0 : Polynomial.aeval A
      (Lagrange.basis (spectralSupport hA) id a *
        Lagrange.basis (spectralSupport hA) id b) =
      Polynomial.aeval A (0 : Polynomial ℂ) := by
    refine aeval_eq_of_eval_eq hA (fun i => ?_)
    rw [Polynomial.eval_mul, basis_eval_eig hA ha i, basis_eval_eig hA hb i,
      Polynomial.eval_zero]
    by_cases h1 : (RCLike.ofReal (hA.eigenvalues i) : ℂ) = a
    · rw [if_pos h1, if_neg (fun h2 => hab (h1.symm.trans h2)), mul_zero]
    · rw [if_neg h1, zero_mul]
  rw [h0, map_zero]

/-- Each spectral projection is idempotent. -/
theorem isIdempotentElem_spectralProjection {A : Matrix n n ℂ}
    (hA : A.IsHermitian) {a : ℂ} (ha : a ∈ spectralSupport hA) :
    IsIdempotentElem (spectralProjection hA a) := by
  show spectralProjection hA a * spectralProjection hA a =
    spectralProjection hA a
  rw [spectralProjection, ← map_mul]
  refine aeval_eq_of_eval_eq hA (fun i => ?_)
  rw [Polynomial.eval_mul, basis_eval_eig hA ha i]
  by_cases h1 : (RCLike.ofReal (hA.eigenvalues i) : ℂ) = a
  · rw [if_pos h1, one_mul]
  · rw [if_neg h1, zero_mul]

/-- The spectral projection in the eigenvector frame: a zero-one
diagonal conjugated by the eigenvector unitary. -/
private theorem spectralProjection_eq_diag {A : Matrix n n ℂ}
    (hA : A.IsHermitian) {a : ℂ} (ha : a ∈ spectralSupport hA) :
    spectralProjection hA a =
      (hA.eigenvectorUnitary : Matrix n n ℂ) *
        Matrix.diagonal
          (fun i =>
            if (RCLike.ofReal (hA.eigenvalues i) : ℂ) = a then 1 else 0) *
        (star hA.eigenvectorUnitary : Matrix n n ℂ) := by
  rw [spectralProjection, aeval_spectral hA]
  have hd : (fun i => (Lagrange.basis (spectralSupport hA) id a).eval
      ((RCLike.ofReal (hA.eigenvalues i) : ℂ))) =
      (fun i =>
        if (RCLike.ofReal (hA.eigenvalues i) : ℂ) = a then (1 : ℂ) else 0) :=
    funext (basis_eval_eig hA ha)
  rw [hd]

/-- Each spectral projection is self-adjoint: its diagonal values on the
eigenvalue frame are zero and one. -/
theorem isSelfAdjoint_spectralProjection {A : Matrix n n ℂ}
    (hA : A.IsHermitian) {a : ℂ} (ha : a ∈ spectralSupport hA) :
    IsSelfAdjoint (spectralProjection hA a) := by
  have hdiag : star (Matrix.diagonal
      (fun i =>
        if (RCLike.ofReal (hA.eigenvalues i) : ℂ) = a then (1 : ℂ) else 0)) =
      Matrix.diagonal
        (fun i =>
          if (RCLike.ofReal (hA.eigenvalues i) : ℂ) = a then 1 else 0) := by
    rw [Matrix.star_eq_conjTranspose, Matrix.diagonal_conjTranspose]
    congr 1
    funext i
    rw [Pi.star_apply]
    by_cases h1 : (RCLike.ofReal (hA.eigenvalues i) : ℂ) = a
    · rw [if_pos h1, star_one]
    · rw [if_neg h1, star_zero]
  show star (spectralProjection hA a) = spectralProjection hA a
  rw [spectralProjection_eq_diag hA ha, star_mul, star_star, star_mul,
    hdiag, ← mul_assoc]

/-- The spectral projections sum to the identity. -/
theorem sum_spectralProjection [Nonempty n] {A : Matrix n n ℂ}
    (hA : A.IsHermitian) :
    ∑ a ∈ spectralSupport hA, spectralProjection hA a = 1 := by
  have hne : (spectralSupport hA).Nonempty :=
    Finset.univ_nonempty.image _
  calc ∑ a ∈ spectralSupport hA, spectralProjection hA a
      = Polynomial.aeval A
        (∑ a ∈ spectralSupport hA,
          Lagrange.basis (spectralSupport hA) id a) :=
        (map_sum (Polynomial.aeval A) _ _).symm
    _ = Polynomial.aeval A (1 : Polynomial ℂ) := by
        rw [Lagrange.sum_basis (Set.injOn_id _) hne]
    _ = 1 := map_one _

/-- The matrix is reconstructed from its spectral projections. -/
theorem eq_sum_smul_spectralProjection [Nonempty n] {A : Matrix n n ℂ}
    (hA : A.IsHermitian) :
    A = ∑ a ∈ spectralSupport hA, a • spectralProjection hA a := by
  have h1 : ∑ a ∈ spectralSupport hA, a • spectralProjection hA a =
      Polynomial.aeval A
        (∑ a ∈ spectralSupport hA,
          Polynomial.C a * Lagrange.basis (spectralSupport hA) id a) := by
    rw [map_sum]
    refine Finset.sum_congr rfl (fun a _ => ?_)
    rw [map_mul, Polynomial.aeval_C, ← Algebra.smul_def, spectralProjection]
  have h2 : Polynomial.aeval A
      (∑ a ∈ spectralSupport hA,
        Polynomial.C a * Lagrange.basis (spectralSupport hA) id a) =
      Polynomial.aeval A (Polynomial.X : Polynomial ℂ) := by
    refine aeval_eq_of_eval_eq hA (fun i => ?_)
    rw [Polynomial.eval_finset_sum, Polynomial.eval_X]
    have hmem : (RCLike.ofReal (hA.eigenvalues i) : ℂ) ∈ spectralSupport hA :=
      Finset.mem_image_of_mem _ (Finset.mem_univ i)
    have hterm : ∀ a ∈ spectralSupport hA,
        (Polynomial.C a *
          Lagrange.basis (spectralSupport hA) id a).eval
            ((RCLike.ofReal (hA.eigenvalues i) : ℂ)) =
        if (RCLike.ofReal (hA.eigenvalues i) : ℂ) = a then a else 0 := by
      intro a ha
      rw [Polynomial.eval_mul, Polynomial.eval_C, basis_eval_eig hA ha i,
        mul_ite, mul_one, mul_zero]
    rw [Finset.sum_congr rfl hterm, Finset.sum_ite_eq, if_pos hmem]
  rw [h1, h2, Polynomial.aeval_X]

/-- Spectral projections at support points are nonzero. -/
theorem spectralProjection_ne_zero {A : Matrix n n ℂ}
    (hA : A.IsHermitian) {a : ℂ} (ha : a ∈ spectralSupport hA) :
    spectralProjection hA a ≠ 0 := by
  obtain ⟨i₀, _, hi₀⟩ := Finset.mem_image.mp ha
  intro hzero
  have hπ := spectralProjection_eq_diag hA ha
  have hD : Matrix.diagonal
      (fun i =>
        if (RCLike.ofReal (hA.eigenvalues i) : ℂ) = a then 1 else 0) =
      (0 : Matrix n n ℂ) := by
    have h1 := unitary_sandwich_cancel hA.eigenvectorUnitary
      (Matrix.diagonal
        (fun i =>
          if (RCLike.ofReal (hA.eigenvalues i) : ℂ) = a then 1 else 0))
    rw [← hπ, hzero, mul_zero, zero_mul] at h1
    exact h1.symm
  have hentry : Matrix.diagonal
      (fun i =>
        if (RCLike.ofReal (hA.eigenvalues i) : ℂ) = a then 1 else 0)
        i₀ i₀ = (0 : Matrix n n ℂ) i₀ i₀ := by
    rw [hD]
  rw [Matrix.diagonal_apply_eq, Matrix.zero_apply, if_pos hi₀] at hentry
  exact one_ne_zero hentry

variable {S : StarSubalgebra ℂ (Matrix n n ℂ)}

/-- Polynomial functional calculus stays inside the center of a star
subalgebra containing the element. -/
theorem IsCentralIn.aeval {h : Matrix n n ℂ} (hh : IsCentralIn S h)
    (p : Polynomial ℂ) : IsCentralIn S (Polynomial.aeval h p) := by
  constructor
  · have hle : Algebra.adjoin ℂ ({h} : Set (Matrix n n ℂ)) ≤
        S.toSubalgebra :=
      Algebra.adjoin_le (Set.singleton_subset_iff.mpr
        (StarSubalgebra.mem_toSubalgebra.mpr hh.mem))
    exact StarSubalgebra.mem_toSubalgebra.mp
      (hle (Polynomial.aeval_mem_adjoin_singleton ℂ h))
  · intro x hx
    have hmem : h ∈ Subalgebra.centralizer ℂ (S : Set (Matrix n n ℂ)) :=
      (Subalgebra.mem_centralizer_iff ℂ).mpr
        (fun g hg => (hh.commutes hg).symm)
    have hle : Algebra.adjoin ℂ ({h} : Set (Matrix n n ℂ)) ≤
        Subalgebra.centralizer ℂ (S : Set (Matrix n n ℂ)) :=
      Algebra.adjoin_le (Set.singleton_subset_iff.mpr hmem)
    have h2 := (Subalgebra.mem_centralizer_iff ℂ).mp
      (hle (Polynomial.aeval_mem_adjoin_singleton ℂ h (p := p))) x hx
    exact h2.symm

/-- The spectral projections of a self-adjoint central element are
central star projections of the subalgebra. -/
theorem isCentralProjectionIn_spectralProjection {h : Matrix n n ℂ}
    (hh : IsCentralIn S h) (hsa : IsSelfAdjoint h) {a : ℂ}
    (ha : a ∈ spectralSupport hsa.isHermitian) :
    IsCentralProjectionIn S (spectralProjection hsa.isHermitian a) where
  isIdempotentElem := isIdempotentElem_spectralProjection _ ha
  isSelfAdjoint := isSelfAdjoint_spectralProjection _ ha
  centralIn := hh.aeval _

/-- Spectral resolution of a self-adjoint central element: a finite
family of central star projections of the subalgebra, pairwise
orthogonal, nonzero, summing to one, and reconstructing the element.
This is the functional-calculus splitting step of the center
decomposition. -/
theorem IsCentralIn.exists_spectral_resolution [Nonempty n]
    {h : Matrix n n ℂ} (hh : IsCentralIn S h) (hsa : IsSelfAdjoint h) :
    ∃ (T : Finset ℂ) (π : ℂ → Matrix n n ℂ),
      (∀ a ∈ T, IsCentralProjectionIn S (π a)) ∧
      (∀ a ∈ T, π a ≠ 0) ∧
      (∀ a ∈ T, ∀ b ∈ T, a ≠ b → π a * π b = 0) ∧
      ∑ a ∈ T, π a = 1 ∧
      h = ∑ a ∈ T, a • π a :=
  ⟨spectralSupport hsa.isHermitian, spectralProjection hsa.isHermitian,
    fun _ ha => isCentralProjectionIn_spectralProjection hh hsa ha,
    fun _ ha => spectralProjection_ne_zero _ ha,
    fun _ ha _ hb hab => spectralProjection_mul_of_ne _ ha hb hab,
    sum_spectralProjection _,
    eq_sum_smul_spectralProjection _⟩

/-- A non-scalar self-adjoint central element splits off a nontrivial
central projection.  This is the induction driver for the block
decomposition of issue 679. -/
theorem IsCentralIn.exists_nontrivial_central_projection [Nonempty n]
    {h : Matrix n n ℂ} (hh : IsCentralIn S h) (hsa : IsSelfAdjoint h)
    (hns : ∀ c : ℂ, h ≠ c • 1) :
    ∃ p : Matrix n n ℂ,
      IsCentralProjectionIn S p ∧ p ≠ 0 ∧ p ≠ 1 := by
  obtain ⟨T, π, hproj, hne, horth, hsum, hrec⟩ :=
    hh.exists_spectral_resolution hsa
  have hTne : T.Nonempty := by
    by_contra hT
    rw [Finset.not_nonempty_iff_eq_empty] at hT
    rw [hT, Finset.sum_empty] at hsum
    exact one_ne_zero hsum.symm
  obtain ⟨a₀, ha₀⟩ := hTne
  by_cases hsingle : ∀ b ∈ T, b = a₀
  · exfalso
    have hT : T = {a₀} := Finset.eq_singleton_iff_unique_mem.mpr
      ⟨ha₀, hsingle⟩
    rw [hT, Finset.sum_singleton] at hsum hrec
    exact hns a₀ (by rw [hrec, hsum])
  · push Not at hsingle
    obtain ⟨b, hbT, hba⟩ := hsingle
    refine ⟨π a₀, hproj a₀ ha₀, hne a₀ ha₀, fun hone => ?_⟩
    have h0 : π b * π a₀ = 0 := horth b hbT a₀ ha₀ hba
    rw [hone, mul_one] at h0
    exact hne b hbT h0

end Spectral

/-! ## The center is spanned by the minimal central projections -/

section Span

variable {n : Type u} [Fintype n] [DecidableEq n]
  {S : StarSubalgebra ℂ (Matrix n n ℂ)}

/-- Against a minimal central projection, every self-adjoint central
element acts as a scalar: exactly one spectral projection survives the
cut. -/
theorem IsCentralIn.mul_minimal_eq_smul [Nonempty n]
    {h e : Matrix n n ℂ} (hh : IsCentralIn S h) (hsa : IsSelfAdjoint h)
    (he : IsMinimalCentralProjectionIn S e) :
    ∃ c : ℂ, h * e = c • e := by
  obtain ⟨T, π, hproj, hne, horth, hsum, hrec⟩ :=
    hh.exists_spectral_resolution hsa
  have heproj : IsCentralProjectionIn S e := he.toIsCentralProjectionIn
  have hq : ∀ a ∈ T, IsCentralProjectionIn S (e * π a) :=
    fun a ha => heproj.mul (hproj a ha)
  have hcomm : ∀ a ∈ T, π a * e = e * π a :=
    fun a ha => (hproj a ha).centralIn.commutes heproj.centralIn.mem
  have hqle : ∀ a ∈ T, (e * π a) * e = e * π a := by
    intro a ha
    rw [mul_assoc, hcomm a ha, ← mul_assoc, heproj.isIdempotentElem.eq]
  have hcases : ∀ a ∈ T, e * π a = 0 ∨ e * π a = e :=
    fun a ha => he.minimal _ (hq a ha) (hqle a ha)
  have hesum : e = ∑ a ∈ T, e * π a := by
    rw [← Finset.mul_sum, hsum, mul_one]
  have hex : ∃ a₀ ∈ T, e * π a₀ ≠ 0 := by
    by_contra hall
    push Not at hall
    exact he.ne_zero (hesum.trans (Finset.sum_eq_zero hall))
  obtain ⟨a₀, ha₀T, ha₀⟩ := hex
  have ha₀e : e * π a₀ = e := (hcases a₀ ha₀T).resolve_left ha₀
  have hzero : ∀ b ∈ T, b ≠ a₀ → e * π b = 0 := by
    intro b hbT hba
    rcases hcases b hbT with h0 | hbe
    · exact h0
    · exfalso
      apply he.ne_zero
      have h1 : π a₀ * π b = 0 := horth a₀ ha₀T b hbT (fun hc => hba hc.symm)
      have h2 : π a₀ * (e * π b) = 0 := by
        rw [← mul_assoc, hcomm a₀ ha₀T, mul_assoc, h1, mul_zero]
      calc e = e * e := heproj.isIdempotentElem.eq.symm
        _ = (e * π a₀) * (e * π b) := by rw [ha₀e, hbe]
        _ = e * (π a₀ * (e * π b)) := by rw [mul_assoc]
        _ = 0 := by rw [h2, mul_zero]
  refine ⟨a₀, ?_⟩
  conv_lhs => rw [hrec]
  rw [Finset.sum_mul]
  have hterm : ∀ b ∈ T, (b • π b) * e = if b = a₀ then a₀ • e else 0 := by
    intro b hbT
    by_cases hb : b = a₀
    · subst hb
      rw [if_pos rfl, smul_mul_assoc, hcomm b hbT, ha₀e]
    · rw [if_neg hb, smul_mul_assoc, hcomm b hbT, hzero b hbT hb,
        smul_zero]
  rw [Finset.sum_congr rfl hterm, Finset.sum_ite_eq', if_pos ha₀T]

private theorem mul_list_sum {M : Type u} [NonUnitalNonAssocSemiring M]
    (a : M) (l : List M) : a * l.sum = (l.map (fun x => a * x)).sum := by
  induction l with
  | nil => simp
  | cons b l ih => simp [List.sum_cons, mul_add, ih]

/-- Every central element lies in the span of a minimal orthogonal
resolution of the identity: split into self-adjoint parts and cut each
part against the minimal projections. -/
theorem IsCentralIn.mem_span_of_resolution [Nonempty n]
    {l : List (Matrix n n ℂ)}
    (hl : ∀ e ∈ l, IsMinimalCentralProjectionIn S e) (hsum : l.sum = 1)
    {z : Matrix n n ℂ} (hz : IsCentralIn S z) :
    z ∈ Submodule.span ℂ {x : Matrix n n ℂ | x ∈ l} := by
  have key : ∀ h : Matrix n n ℂ, IsCentralIn S h → IsSelfAdjoint h →
      h ∈ Submodule.span ℂ {x : Matrix n n ℂ | x ∈ l} := by
    intro h hc hsa
    have hrepr : h = (l.map (fun e => h * e)).sum := by
      rw [← mul_list_sum, hsum, mul_one]
    rw [hrepr]
    refine list_sum_mem (fun x hx => ?_)
    rw [List.mem_map] at hx
    obtain ⟨e, hel, rfl⟩ := hx
    obtain ⟨c, hce⟩ := hc.mul_minimal_eq_smul hsa (hl e hel)
    rw [hce]
    exact Submodule.smul_mem _ c (Submodule.subset_span hel)
  have hsa₁ : IsSelfAdjoint (((2 : ℂ)⁻¹) • (z + star z)) := by
    show star _ = _
    rw [star_smul, star_add, star_star, star_inv₀, star_ofNat,
      add_comm (star z) z]
  have hsa₂ : IsSelfAdjoint (((2 : ℂ)⁻¹ * Complex.I) • (star z - z)) := by
    show star _ = _
    rw [star_smul, star_sub, star_star, star_mul', star_inv₀, star_ofNat,
      Complex.star_def, Complex.conj_I]
    rw [show (2 : ℂ)⁻¹ * -Complex.I = -((2 : ℂ)⁻¹ * Complex.I) by ring]
    rw [neg_smul, ← smul_neg, neg_sub]
  have hz₁ : IsCentralIn S (((2 : ℂ)⁻¹) • (z + star z)) :=
    (hz.add hz.star).smul _
  have hz₂ : IsCentralIn S (((2 : ℂ)⁻¹ * Complex.I) • (star z - z)) :=
    (hz.star.sub hz).smul _
  have hdecomp : z = ((2 : ℂ)⁻¹) • (z + star z) +
      Complex.I • (((2 : ℂ)⁻¹ * Complex.I) • (star z - z)) := by
    rw [smul_smul, smul_add, smul_sub]
    rw [show Complex.I * ((2 : ℂ)⁻¹ * Complex.I) = -(2 : ℂ)⁻¹ by
      rw [mul_comm (2 : ℂ)⁻¹ Complex.I, ← mul_assoc, Complex.I_mul_I]; ring]
    rw [neg_smul, neg_smul, neg_sub_neg]
    rw [show ((2 : ℂ)⁻¹ • z + (2 : ℂ)⁻¹ • star z +
        ((2 : ℂ)⁻¹ • z - (2 : ℂ)⁻¹ • star z)) =
        ((2 : ℂ)⁻¹ + (2 : ℂ)⁻¹) • z by
      rw [add_smul]; abel]
    rw [show ((2 : ℂ)⁻¹ + (2 : ℂ)⁻¹ : ℂ) = 1 by norm_num, one_smul]
  rw [hdecomp]
  exact Submodule.add_mem _ (key _ hz₁ hsa₁)
    (Submodule.smul_mem _ _ (key _ hz₂ hsa₂))

/-- Step 1 of the block decomposition: the center of a star subalgebra of
a finite complex matrix algebra is spanned by a finite orthogonal family
of minimal central projections summing to one.  Existence comes from the
finite descent, spanning from the spectral splitting of self-adjoint
central elements against the minimal projections. -/
theorem exists_minimal_central_resolution_spanning [Nonempty n]
    (S : StarSubalgebra ℂ (Matrix n n ℂ)) :
    ∃ l : List (Matrix n n ℂ),
      (∀ e ∈ l, IsMinimalCentralProjectionIn S e) ∧
      l.Pairwise (fun a b => a * b = 0) ∧
      l.sum = 1 ∧
      ∀ z : Matrix n n ℂ, IsCentralIn S z ↔
        z ∈ Submodule.span ℂ {x : Matrix n n ℂ | x ∈ l} := by
  obtain ⟨l, hmin, horth, hsum⟩ := exists_minimal_central_resolution S
  refine ⟨l, hmin, horth, hsum, fun z => ⟨fun hz =>
    hz.mem_span_of_resolution hmin hsum, fun hz => ?_⟩⟩
  have hsub : {x : Matrix n n ℂ | x ∈ l} ⊆
      ((centerSubmodule S : Submodule ℂ (Matrix n n ℂ)) :
        Set (Matrix n n ℂ)) :=
    fun x hx => (hmin x hx).centralIn
  exact Submodule.span_le.mpr hsub hz

/-- Uniqueness of the minimal resolution up to permutation: two minimal
orthogonal resolutions of the identity list the same projections. -/
theorem minimal_central_resolution_unique
    {l₁ l₂ : List (Matrix n n ℂ)}
    (h₁ : ∀ e ∈ l₁, IsMinimalCentralProjectionIn S e)
    (o₁ : l₁.Pairwise (fun a b => a * b = 0)) (s₁ : l₁.sum = 1)
    (h₂ : ∀ e ∈ l₂, IsMinimalCentralProjectionIn S e)
    (o₂ : l₂.Pairwise (fun a b => a * b = 0)) (s₂ : l₂.sum = 1) :
    l₁.Perm l₂ := by
  have hmem : ∀ (la lb : List (Matrix n n ℂ)),
      (∀ e ∈ la, IsMinimalCentralProjectionIn S e) →
      (∀ e ∈ lb, IsMinimalCentralProjectionIn S e) →
      lb.sum = 1 → ∀ e ∈ la, e ∈ lb := by
    intro la lb ha hb hbsum e hea
    have hee := ha e hea
    have hesum : e = (lb.map (fun f => e * f)).sum := by
      rw [← mul_list_sum, hbsum, mul_one]
    have hex : ∃ f ∈ lb, e * f ≠ 0 := by
      by_contra hall
      push Not at hall
      refine hee.ne_zero (hesum.trans (List.sum_eq_zero ?_))
      intro x hx
      rw [List.mem_map] at hx
      obtain ⟨f, hf, rfl⟩ := hx
      exact hall f hf
    obtain ⟨f, hfb, hef⟩ := hex
    have hff := hb f hfb
    have hprod : IsCentralProjectionIn S (e * f) :=
      hee.toIsCentralProjectionIn.mul hff.toIsCentralProjectionIn
    have hcomm : e * f = f * e :=
      hee.centralIn.commute hff.centralIn
    have hle_e : (e * f) * e = e * f := by
      rw [hcomm, mul_assoc, hee.isIdempotentElem.eq, ← hcomm]
    have hle_f : (e * f) * f = e * f := by
      rw [mul_assoc, hff.isIdempotentElem.eq]
    have he_eq : e * f = e :=
      (hee.minimal _ hprod hle_e).resolve_left hef
    have hf_eq : e * f = f :=
      (hff.minimal _ hprod hle_f).resolve_left hef
    rw [he_eq.symm.trans hf_eq]
    exact hfb
  have hnodup : ∀ (la : List (Matrix n n ℂ)),
      (∀ e ∈ la, IsMinimalCentralProjectionIn S e) →
      la.Pairwise (fun a b => a * b = 0) → la.Nodup := by
    intro la ha horth
    induction la with
    | nil => exact List.nodup_nil
    | cons e la ih =>
      rw [List.pairwise_cons] at horth
      rw [List.nodup_cons]
      refine ⟨fun hmem => ?_, ih (fun x hx => ha x (List.mem_cons_of_mem e hx))
        horth.2⟩
      have h0 : e * e = 0 := horth.1 e hmem
      rw [(ha e List.mem_cons_self).isIdempotentElem.eq] at h0
      exact (ha e List.mem_cons_self).ne_zero h0
  exact ((hnodup l₁ h₁ o₁).subperm
      (fun e he => hmem l₁ l₂ h₁ h₂ s₂ e he)).antisymm
    ((hnodup l₂ h₂ o₂).subperm (fun e he => hmem l₂ l₁ h₂ h₁ s₁ e he))

end Span

/-! ## Two-projection factorization through the central cut -/

section PairCut

variable {A : Type u} [Ring A] [StarRing A] [Algebra ℂ A] [StarModule ℂ A]

/-- Centrality relative to the full subalgebra is ambient centrality. -/
theorem isCentralProjectionIn_top_iff {p : A} :
    IsCentralProjectionIn (⊤ : StarSubalgebra ℂ A) p ↔
      IsCentralStarProjection p := by
  constructor
  · intro hp
    exact { toIsStarProjection := hp.toIsStarProjection
            central := fun x => hp.centralIn.2 x (by trivial) }
  · intro hp
    exact { toIsStarProjection := hp.toIsStarProjection
            centralIn := ⟨by trivial, fun x _ => hp.central x⟩ }

omit [StarRing A] [Algebra ℂ A] [StarModule ℂ A] in
/-- In a two-projection resolution of the identity the second projection
is the complement of the first. -/
theorem pair_resolution_eq_one_sub {p q : A} (hsum : p + q = 1) :
    q = 1 - p :=
  eq_sub_of_add_eq' hsum

omit [Algebra ℂ A] [StarModule ℂ A] in
/-- A two-projection resolution is automatically orthogonal. -/
theorem pair_resolution_orthogonal {p q : A}
    (hp : IsCentralStarProjection p) (hsum : p + q = 1) : q * p = 0 := by
  rw [pair_resolution_eq_one_sub hsum, sub_mul, one_mul,
    hp.isIdempotentElem.eq, sub_self]

/-- Two-projection case of the iterated factorization: a resolution of
the identity by two central star projections factors the algebra as the
direct sum of the two corners, through the central cut equivalence of
`CentralBlocks.lean`.  The `k`-fold iterated product over a full minimal
resolution is a separate obligation of issue 679. -/
noncomputable def pairResolutionCutEquiv {p q : A}
    (hp : IsCentralStarProjection p) (hq : IsCentralStarProjection q)
    (hsum : p + q = 1) :
    A ≃⋆ₐ[ℂ] (centralCorner p hp × centralCorner q hq) := by
  have hq1 : q = 1 - p := pair_resolution_eq_one_sub hsum
  subst hq1
  exact centralCutEquiv p hp

end PairCut

#print axioms OPH.Dynamics.IsCentralProjectionIn.exists_minimal_decomposition
#print axioms OPH.Dynamics.exists_minimal_central_resolution
#print axioms OPH.Dynamics.IsCentralIn.exists_spectral_resolution
#print axioms OPH.Dynamics.IsCentralIn.exists_nontrivial_central_projection
#print axioms OPH.Dynamics.IsCentralIn.mul_minimal_eq_smul
#print axioms OPH.Dynamics.exists_minimal_central_resolution_spanning
#print axioms OPH.Dynamics.minimal_central_resolution_unique
#print axioms OPH.Dynamics.pairResolutionCutEquiv

end OPH.Dynamics
