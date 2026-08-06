import EventAlgebra.Basic

/-!
# Finite frequency concentration for the Born weight

Route 4 of the Born packet: the finite Finkelstein-Hartle argument, stated
with its circularity boundary.

For a state `ρ` on `Fin d` and a projection event `E`, the `N`-copy product
state is the iterated Kronecker power `kronPow ρ N`. The frequency operator
`freqOp E N` is the average of `E` acting on each tensor factor, the sum of
the position embeddings `embedAt E N i` divided by `N`. The file proves,
exactly and for every `N`:

* the Born mean of `freqOp E N` in the product state is `Tr(ρ E)`
  (`concentrationMass_freqOp`);
* the Born second moment is `Tr(ρE)/N + (1 - 1/N) Tr(ρE)^2`, by the exact
  cross-term factorization for the commuting position embeddings
  (`concentrationMass_freqOp_sq`, `embedAt_mul_embedAt_comm`);
* the variance is `(Tr(ρE) - Tr(ρE)^2)/N`, exactly (`variance_freqOp`);
* the Chebyshev bound: every event dominated by the squared deviation
  operator at scale `ε` carries mass at most `(Tr(ρE) - Tr(ρE)^2)/(N ε^2)`
  (`chebyshev_mass_bound`), a bound tending to zero in `N`
  (`varianceRatio_tendsto_zero`);
* uniqueness of the stable frequency: for `p ≠ Tr(ρ E)` and any threshold
  `θ > 0`, every `N` beyond an explicit bound forces the mass of every band
  event around `p` below `θ` (`born_unique_stable_frequency`,
  `exists_concentration_refuting_N`).

The spectral reading: `freqOp E N` is Hermitian (`freqOp_isHermitian`), and
for each `ε`-band around a value `p` with `|p - Tr(ρE)| ≥ 2ε` the spectral
projection outside the band around `Tr(ρE)`, in particular the band
projection at `p`, satisfies the domination hypothesis entrywise in the
eigenbasis; the witness section instantiates this with an explicit rational
band projection at `d = 2`, `N = 2`.

## Circularity boundary

The mass functional in every statement above is `concentrationMass`, the
trace pairing of the product state against the event: the Born pairing
itself, one tensor level up. The named Prop `ConcentrationMassIsBornPairing`
records this exactly, with a `rfl` proof. The packet is therefore a
consistency-and-uniqueness result: the Born weight is the unique frequency
value consistent with the concentration receipts at every `N`, while the
receipts are themselves stated in the Born measure of the product state.
This route never derives the Born rule from frequency data alone; it
complements the noncontextuality and envariance routes by showing that Born
is the unique self-consistent frequency law.
-/

namespace EventAlgebra.FrequencyConcentration

open Matrix
open scoped ComplexOrder Kronecker

variable {d : ℕ}

/-! ## The iterated tensor index and the Kronecker power -/

/-- Index type of the `N`-fold tensor power of `Fin d`: nested pairs with a
`PUnit` terminator. Reducible so that the successor case unfolds to the
product type during unification. -/
@[reducible] def PowIdx (d : ℕ) : ℕ → Type
  | 0 => PUnit
  | N + 1 => Fin d × PowIdx d N

@[reducible] instance instFintypePowIdx (d : ℕ) :
    (N : ℕ) → Fintype (PowIdx d N)
  | 0 => inferInstanceAs (Fintype PUnit)
  | N + 1 =>
    letI := instFintypePowIdx d N
    inferInstanceAs (Fintype (Fin d × PowIdx d N))

@[reducible] instance instDecidableEqPowIdx (d : ℕ) :
    (N : ℕ) → DecidableEq (PowIdx d N)
  | 0 => inferInstanceAs (DecidableEq PUnit)
  | N + 1 =>
    letI := instDecidableEqPowIdx d N
    inferInstanceAs (DecidableEq (Fin d × PowIdx d N))

/-- The `N`-fold Kronecker power `A ⊗ ⋯ ⊗ A`, the `N`-copy product matrix. -/
def kronPow (A : Matrix (Fin d) (Fin d) ℂ) :
    (N : ℕ) → Matrix (PowIdx d N) (PowIdx d N) ℂ
  | 0 => 1
  | N + 1 => A ⊗ₖ kronPow A N

@[simp]
theorem kronPow_zero (A : Matrix (Fin d) (Fin d) ℂ) :
    kronPow A 0 = 1 := rfl

@[simp]
theorem kronPow_succ (A : Matrix (Fin d) (Fin d) ℂ) (N : ℕ) :
    kronPow A (N + 1) = A ⊗ₖ kronPow A N := rfl

/-- **Trace-dependent.** The trace of the Kronecker power is the power of
the trace; for a state the product matrix therefore has trace one. -/
theorem trace_kronPow (A : Matrix (Fin d) (Fin d) ℂ) :
    ∀ N : ℕ, (kronPow A N).trace = A.trace ^ N
  | 0 => by
    rw [kronPow_zero, pow_zero, Matrix.trace_one,
      show Fintype.card (PowIdx d 0) = 1 from rfl, Nat.cast_one]
  | N + 1 => by
    rw [kronPow_succ]
    show (A ⊗ₖ kronPow A N).trace = A.trace ^ (N + 1)
    rw [Matrix.trace_kronecker, trace_kronPow A N, pow_succ']

/-- **Algebra-only.** The Kronecker power of a positive-semidefinite matrix
is positive semidefinite. -/
theorem kronPow_posSemidef {A : Matrix (Fin d) (Fin d) ℂ}
    (hA : A.PosSemidef) : ∀ N : ℕ, (kronPow A N).PosSemidef
  | 0 => by rw [kronPow_zero]; exact Matrix.PosSemidef.one
  | N + 1 => by
    rw [kronPow_succ]
    exact hA.kronecker (kronPow_posSemidef hA N)

/-- **Trace-dependent.** Positivity of the trace pairing of two
positive-semidefinite matrices over an arbitrary finite index type. This
generalizes `EventAlgebra.expectation_nonneg` from `Fin n` to the tensor
index types of this file. -/
theorem trace_mul_nonneg_of_posSemidef {ι : Type*} [Fintype ι]
    [DecidableEq ι] {A B : Matrix ι ι ℂ} (hA : A.PosSemidef)
    (hB : B.PosSemidef) : 0 ≤ (A * B).trace := by
  classical
  set V : Matrix ι ι ℂ := (hB.1.eigenvectorUnitary : Matrix ι ι ℂ) with hV
  set dg : ι → ℂ := RCLike.ofReal ∘ hB.1.eigenvalues with hd
  have hspec : B = V * diagonal dg * star V := by
    rw [hV, hd]
    conv_lhs => rw [hB.1.spectral_theorem, Unitary.conjStarAlgAut_apply]
  have hsandwich : (star V * A * V).PosSemidef := by
    have := hA.conjTranspose_mul_mul_same V
    rwa [← star_eq_conjTranspose] at this
  have hcycle : (A * B).trace = ((star V * A * V) * diagonal dg).trace := by
    rw [hspec,
      show A * (V * diagonal dg * star V) = (A * V * diagonal dg) * star V by
        simp only [mul_assoc],
      trace_mul_comm]
    simp only [mul_assoc]
  rw [hcycle]
  simp only [Matrix.trace, Matrix.diag, Matrix.mul_diagonal]
  refine Finset.sum_nonneg fun i _ => mul_nonneg hsandwich.diag_nonneg ?_
  rw [hd]
  simp only [Function.comp_apply]
  exact RCLike.ofReal_nonneg.mpr (hB.eigenvalues_nonneg i)

/-! ## Position embeddings and the frequency operator -/

/-- The embedding of a single-copy observable `E` at tensor position `i` of
the `N`-fold product: `1 ⊗ ⋯ ⊗ E ⊗ ⋯ ⊗ 1` with `E` in slot `i`. -/
def embedAt (E : Matrix (Fin d) (Fin d) ℂ) :
    (N : ℕ) → Fin N → Matrix (PowIdx d N) (PowIdx d N) ℂ
  | 0, i => i.elim0
  | N + 1, i =>
    Fin.cases (E ⊗ₖ (1 : Matrix (PowIdx d N) (PowIdx d N) ℂ))
      (fun j => (1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt E N j) i

@[simp]
theorem embedAt_zero (E : Matrix (Fin d) (Fin d) ℂ) (N : ℕ) :
    embedAt E (N + 1) 0 = E ⊗ₖ (1 : Matrix (PowIdx d N) (PowIdx d N) ℂ) := by
  simp [embedAt]

@[simp]
theorem embedAt_succ (E : Matrix (Fin d) (Fin d) ℂ) (N : ℕ) (j : Fin N) :
    embedAt E (N + 1) j.succ
      = (1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt E N j := by
  simp [embedAt]

/-- **Algebra-only.** Conjugate transpose passes through the position
embedding. -/
theorem embedAt_conjTranspose (E : Matrix (Fin d) (Fin d) ℂ) :
    ∀ (N : ℕ) (i : Fin N), (embedAt E N i)ᴴ = embedAt Eᴴ N i
  | N + 1, i => by
    rcases Fin.eq_zero_or_eq_succ i with rfl | ⟨j, rfl⟩
    · rw [embedAt_zero, embedAt_zero]
      show (E ⊗ₖ (1 : Matrix (PowIdx d N) (PowIdx d N) ℂ))ᴴ
          = Eᴴ ⊗ₖ (1 : Matrix (PowIdx d N) (PowIdx d N) ℂ)
      rw [Matrix.conjTranspose_kronecker, Matrix.conjTranspose_one]
    · rw [embedAt_succ, embedAt_succ]
      show ((1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt E N j)ᴴ
          = (1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt Eᴴ N j
      rw [Matrix.conjTranspose_kronecker, Matrix.conjTranspose_one,
        embedAt_conjTranspose E N j]

/-- **Algebra-only.** Same-position embeddings multiply as the embedded
product; for an event this collapses the diagonal terms of the second
moment. -/
theorem embedAt_mul_embedAt_same (E F : Matrix (Fin d) (Fin d) ℂ) :
    ∀ (N : ℕ) (i : Fin N),
      embedAt E N i * embedAt F N i = embedAt (E * F) N i
  | N + 1, i => by
    rcases Fin.eq_zero_or_eq_succ i with rfl | ⟨j, rfl⟩
    · rw [embedAt_zero, embedAt_zero, embedAt_zero]
      show (E ⊗ₖ (1 : Matrix (PowIdx d N) (PowIdx d N) ℂ))
            * (F ⊗ₖ (1 : Matrix (PowIdx d N) (PowIdx d N) ℂ))
          = (E * F) ⊗ₖ (1 : Matrix (PowIdx d N) (PowIdx d N) ℂ)
      rw [← Matrix.mul_kronecker_mul, one_mul]
    · rw [embedAt_succ, embedAt_succ, embedAt_succ]
      show ((1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt E N j)
            * ((1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt F N j)
          = (1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt (E * F) N j
      rw [← Matrix.mul_kronecker_mul, one_mul,
        embedAt_mul_embedAt_same E F N j]

/-- **Algebra-only.** Embeddings at distinct positions commute: the
frequency operator is an average of commuting projections. -/
theorem embedAt_mul_embedAt_comm (E F : Matrix (Fin d) (Fin d) ℂ) :
    ∀ (N : ℕ) (i j : Fin N), i ≠ j →
      embedAt E N i * embedAt F N j = embedAt F N j * embedAt E N i
  | N + 1, i, j, hij => by
    rcases Fin.eq_zero_or_eq_succ i with rfl | ⟨i', rfl⟩ <;>
      rcases Fin.eq_zero_or_eq_succ j with rfl | ⟨j', rfl⟩
    · exact absurd rfl hij
    · rw [embedAt_zero, embedAt_succ]
      show (E ⊗ₖ (1 : Matrix (PowIdx d N) (PowIdx d N) ℂ))
            * ((1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt F N j')
          = ((1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt F N j')
            * (E ⊗ₖ (1 : Matrix (PowIdx d N) (PowIdx d N) ℂ))
      rw [← Matrix.mul_kronecker_mul, ← Matrix.mul_kronecker_mul, one_mul,
        one_mul, mul_one, mul_one]
    · rw [embedAt_zero, embedAt_succ]
      show ((1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt E N i')
            * (F ⊗ₖ (1 : Matrix (PowIdx d N) (PowIdx d N) ℂ))
          = (F ⊗ₖ (1 : Matrix (PowIdx d N) (PowIdx d N) ℂ))
            * ((1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt E N i')
      rw [← Matrix.mul_kronecker_mul, ← Matrix.mul_kronecker_mul, one_mul,
        one_mul, mul_one, mul_one]
    · rw [embedAt_succ, embedAt_succ]
      show ((1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt E N i')
            * ((1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt F N j')
          = ((1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt F N j')
            * ((1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt E N i')
      rw [← Matrix.mul_kronecker_mul, ← Matrix.mul_kronecker_mul, one_mul,
        embedAt_mul_embedAt_comm E F N i' j'
          (fun h => hij (congrArg Fin.succ h))]

/-- **Trace-dependent.** The single-position mean: the product state pairs
with one embedded copy of `E` to the single-copy pairing `Tr(ρ E)`. -/
theorem trace_kronPow_mul_embedAt (ρ E : Matrix (Fin d) (Fin d) ℂ)
    (hρ : ρ.trace = 1) :
    ∀ (N : ℕ) (i : Fin N),
      (kronPow ρ N * embedAt E N i).trace = (ρ * E).trace
  | N + 1, i => by
    rcases Fin.eq_zero_or_eq_succ i with rfl | ⟨j, rfl⟩
    · rw [kronPow_succ, embedAt_zero]
      show ((ρ ⊗ₖ kronPow ρ N)
            * (E ⊗ₖ (1 : Matrix (PowIdx d N) (PowIdx d N) ℂ))).trace
          = (ρ * E).trace
      rw [← Matrix.mul_kronecker_mul, mul_one, Matrix.trace_kronecker,
        trace_kronPow, hρ, one_pow, mul_one]
    · rw [kronPow_succ, embedAt_succ]
      show ((ρ ⊗ₖ kronPow ρ N)
            * ((1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt E N j)).trace
          = (ρ * E).trace
      rw [← Matrix.mul_kronecker_mul, mul_one, Matrix.trace_kronecker, hρ,
        one_mul, trace_kronPow_mul_embedAt ρ E hρ N j]

/-- **Trace-dependent.** The exact cross-term factorization: the product
state pairs with two embedded copies of `E` at distinct positions to
`Tr(ρ E)^2`. This is the cancellation that makes the variance scale as
`1/N`. -/
theorem trace_kronPow_mul_embedAt_pair (ρ E : Matrix (Fin d) (Fin d) ℂ)
    (hρ : ρ.trace = 1) :
    ∀ (N : ℕ) (i j : Fin N), i ≠ j →
      (kronPow ρ N * (embedAt E N i * embedAt E N j)).trace
        = (ρ * E).trace ^ 2
  | N + 1, i, j, hij => by
    rcases Fin.eq_zero_or_eq_succ i with rfl | ⟨i', rfl⟩ <;>
      rcases Fin.eq_zero_or_eq_succ j with rfl | ⟨j', rfl⟩
    · exact absurd rfl hij
    · rw [kronPow_succ, embedAt_zero, embedAt_succ]
      show ((ρ ⊗ₖ kronPow ρ N)
            * ((E ⊗ₖ (1 : Matrix (PowIdx d N) (PowIdx d N) ℂ))
              * ((1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt E N j'))).trace
          = (ρ * E).trace ^ 2
      rw [← Matrix.mul_kronecker_mul, one_mul, mul_one,
        ← Matrix.mul_kronecker_mul, Matrix.trace_kronecker,
        trace_kronPow_mul_embedAt ρ E hρ N j', pow_two]
    · rw [kronPow_succ, embedAt_zero, embedAt_succ]
      show ((ρ ⊗ₖ kronPow ρ N)
            * (((1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt E N i')
              * (E ⊗ₖ (1 : Matrix (PowIdx d N) (PowIdx d N) ℂ)))).trace
          = (ρ * E).trace ^ 2
      rw [← Matrix.mul_kronecker_mul, one_mul, mul_one,
        ← Matrix.mul_kronecker_mul, Matrix.trace_kronecker,
        trace_kronPow_mul_embedAt ρ E hρ N i', pow_two]
    · rw [kronPow_succ, embedAt_succ, embedAt_succ]
      show ((ρ ⊗ₖ kronPow ρ N)
            * (((1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt E N i')
              * ((1 : Matrix (Fin d) (Fin d) ℂ) ⊗ₖ embedAt E N j'))).trace
          = (ρ * E).trace ^ 2
      rw [← Matrix.mul_kronecker_mul, one_mul, ← Matrix.mul_kronecker_mul,
        mul_one, Matrix.trace_kronecker, hρ, one_mul,
        trace_kronPow_mul_embedAt_pair ρ E hρ N i' j'
          (fun h => hij (congrArg Fin.succ h))]

/-- The frequency operator `F_N`: the average of `E` over the `N` tensor
positions. Its spectrum for a projection `E` is `{k/N : 0 ≤ k ≤ N}`, the
possible relative frequencies. -/
noncomputable def freqOp (E : Matrix (Fin d) (Fin d) ℂ) (N : ℕ) :
    Matrix (PowIdx d N) (PowIdx d N) ℂ :=
  (N : ℂ)⁻¹ • ∑ i : Fin N, embedAt E N i

/-- **Algebra-only.** The frequency operator of a Hermitian observable is
Hermitian, so its spectral calculus is available. -/
theorem freqOp_isHermitian {E : Matrix (Fin d) (Fin d) ℂ}
    (hE : E.IsHermitian) (N : ℕ) : (freqOp E N).IsHermitian := by
  show (freqOp E N)ᴴ = freqOp E N
  rw [freqOp, Matrix.conjTranspose_smul, Matrix.conjTranspose_sum]
  congr 1
  · simp
  · exact Finset.sum_congr rfl fun i _ => by
      rw [embedAt_conjTranspose, hE.eq]

/-! ## The mass functional and the circularity boundary -/

/-- The mass functional used in every concentration statement of this file:
the trace pairing of the `N`-copy product state against an event of the
product algebra. -/
noncomputable def concentrationMass (ρ : Matrix (Fin d) (Fin d) ℂ) (N : ℕ)
    (Q : Matrix (PowIdx d N) (PowIdx d N) ℂ) : ℂ :=
  (kronPow ρ N * Q).trace

@[simp]
theorem concentrationMass_def (ρ : Matrix (Fin d) (Fin d) ℂ) (N : ℕ)
    (Q : Matrix (PowIdx d N) (PowIdx d N) ℂ) :
    concentrationMass ρ N Q = (kronPow ρ N * Q).trace := rfl

/-- **The circularity boundary, as a named Prop.** The mass functional of
the concentration statements is the Born pairing of the product state: the
same trace pairing `Tr(σ Q)` that `EventAlgebra.bornWeight` uses one tensor
level down. Concentration therefore quantifies Born weights of frequency
bands inside the Born measure of `kronPow ρ N`; it never defines a
probability notion outside it. This is why the packet is a
consistency-and-uniqueness result and never a standalone derivation of
the Born rule. -/
def ConcentrationMassIsBornPairing (d : ℕ) : Prop :=
  ∀ (ρ : Matrix (Fin d) (Fin d) ℂ) (N : ℕ)
    (Q : Matrix (PowIdx d N) (PowIdx d N) ℂ),
    concentrationMass ρ N Q = (kronPow ρ N * Q).trace

/-- **The circularity boundary holds definitionally.** The proof is `rfl`:
`concentrationMass` is the Born trace pairing of the product state by
definition, so every concentration receipt below is a statement inside the
Born measure, and the packet is a consistency-and-uniqueness result. -/
theorem concentrationMass_isBornPairing (d : ℕ) :
    ConcentrationMassIsBornPairing d := fun _ _ _ => rfl

/-! ## The exact moment theorems -/

/-- **Trace-dependent. The Born mean.** For every `N ≥ 1` the product state
assigns the frequency operator the mean `Tr(ρ E)`: the Born weight of the
single-copy event, exactly, with no large-`N` limit. -/
theorem concentrationMass_freqOp (ρ E : Matrix (Fin d) (Fin d) ℂ)
    (hρ : ρ.trace = 1) {N : ℕ} (hN : 0 < N) :
    concentrationMass ρ N (freqOp E N) = bornWeight ρ E := by
  have hNc : (N : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hN.ne'
  rw [concentrationMass_def, freqOp, Matrix.mul_smul, Matrix.trace_smul,
    Finset.mul_sum, Matrix.trace_sum,
    Finset.sum_congr rfl fun i _ => trace_kronPow_mul_embedAt ρ E hρ N i,
    Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul,
    smul_eq_mul, bornWeight, inv_mul_cancel_left₀ hNc]

/-- **Trace-dependent. The Born second moment, exact.** For a projection
`E` and every `N ≥ 1`,
`⟨F_N²⟩ = Tr(ρE)/N + (1 - 1/N) Tr(ρE)²`: the diagonal terms contribute the
first summand, the `N(N-1)` commuting cross terms factor exactly into the
second. -/
theorem concentrationMass_freqOp_sq (ρ E : Matrix (Fin d) (Fin d) ℂ)
    (hρ : ρ.trace = 1) (hE : E * E = E) {N : ℕ} (hN : 0 < N) :
    concentrationMass ρ N (freqOp E N * freqOp E N)
      = (N : ℂ)⁻¹ * bornWeight ρ E
        + (1 - (N : ℂ)⁻¹) * bornWeight ρ E ^ 2 := by
  have hNc : (N : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hN.ne'
  have hinner : ∀ i : Fin N,
      ∑ j : Fin N, (kronPow ρ N * (embedAt E N i * embedAt E N j)).trace
        = (ρ * E).trace + ((N : ℂ) - 1) * (ρ * E).trace ^ 2 := by
    intro i
    rw [← Finset.add_sum_erase Finset.univ _ (Finset.mem_univ i)]
    congr 1
    · rw [embedAt_mul_embedAt_same, hE, trace_kronPow_mul_embedAt ρ E hρ N i]
    · rw [Finset.sum_congr rfl fun j hj =>
        trace_kronPow_mul_embedAt_pair ρ E hρ N i j
          (Ne.symm (Finset.ne_of_mem_erase hj)),
        Finset.sum_const, Finset.card_erase_of_mem (Finset.mem_univ i),
        Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, Nat.cast_sub hN,
        Nat.cast_one]
  have hFF : freqOp E N * freqOp E N
      = ((N : ℂ)⁻¹ * (N : ℂ)⁻¹) • ∑ i : Fin N, ∑ j : Fin N,
          embedAt E N i * embedAt E N j := by
    rw [freqOp, Matrix.smul_mul, Matrix.mul_smul, smul_smul,
      Finset.sum_mul_sum]
  rw [concentrationMass_def, hFF, Matrix.mul_smul, Matrix.trace_smul,
    smul_eq_mul]
  have hexp : (kronPow ρ N * ∑ i : Fin N, ∑ j : Fin N,
      embedAt E N i * embedAt E N j).trace
      = ∑ i : Fin N, ∑ j : Fin N,
          (kronPow ρ N * (embedAt E N i * embedAt E N j)).trace := by
    simp only [Finset.mul_sum, Matrix.trace_sum]
  rw [hexp, Finset.sum_congr rfl fun i _ => hinner i, Finset.sum_const,
    Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, bornWeight]
  field_simp

/-- The deviation of the frequency operator from the Born weight of the
single-copy event. -/
noncomputable def deviationOp (ρ E : Matrix (Fin d) (Fin d) ℂ) (N : ℕ) :
    Matrix (PowIdx d N) (PowIdx d N) ℂ :=
  freqOp E N - bornWeight ρ E • (1 : Matrix (PowIdx d N) (PowIdx d N) ℂ)

/-- **Trace-dependent. The exact variance.** For a projection `E` and every
`N ≥ 1` the product state assigns the squared deviation the value
`(Tr(ρE) - Tr(ρE)²)/N`: the exact finite variance of the frequency
operator, with the `1/N` decay carried by the cross-term factorization. -/
theorem variance_freqOp (ρ E : Matrix (Fin d) (Fin d) ℂ)
    (hρ : ρ.trace = 1) (hE : E * E = E) {N : ℕ} (hN : 0 < N) :
    concentrationMass ρ N (deviationOp ρ E N * deviationOp ρ E N)
      = (bornWeight ρ E - bornWeight ρ E ^ 2) / N := by
  have hNc : (N : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hN.ne'
  have hexpand : deviationOp ρ E N * deviationOp ρ E N
      = freqOp E N * freqOp E N - bornWeight ρ E • freqOp E N
        - bornWeight ρ E • freqOp E N
        + (bornWeight ρ E * bornWeight ρ E)
            • (1 : Matrix (PowIdx d N) (PowIdx d N) ℂ) := by
    simp only [deviationOp, sub_mul, mul_sub, Matrix.smul_mul,
      Matrix.mul_smul, Matrix.one_mul, Matrix.mul_one, smul_smul]
    abel
  rw [concentrationMass_def, hexpand]
  simp only [Matrix.mul_sub, Matrix.mul_add, Matrix.mul_smul,
    Matrix.trace_sub, Matrix.trace_add, Matrix.trace_smul, smul_eq_mul,
    Matrix.mul_one]
  have hM1 : (kronPow ρ N * freqOp E N).trace = bornWeight ρ E := by
    have := concentrationMass_freqOp ρ E hρ hN
    rwa [concentrationMass_def] at this
  have hM2 : (kronPow ρ N * (freqOp E N * freqOp E N)).trace
      = (N : ℂ)⁻¹ * bornWeight ρ E
        + (1 - (N : ℂ)⁻¹) * bornWeight ρ E ^ 2 := by
    have := concentrationMass_freqOp_sq ρ E hρ hE hN
    rwa [concentrationMass_def] at this
  rw [hM1, hM2, trace_kronPow, hρ, one_pow]
  field_simp
  ring

/-! ## The exact Chebyshev concentration bound

The bound is stated for every event `Q` of the product algebra dominated by
the squared deviation operator at scale `ε`: matrices with
`ε² • Q ≤ (F_N - Tr(ρE))²` in the positive-semidefinite order. For the
Hermitian `F_N` this hypothesis holds in particular for every spectral
projection of `F_N` onto eigenvalues at distance at least `ε` from
`Tr(ρ E)`, since both sides are then diagonal in one eigenbasis and the
entrywise inequality is `(λ - Tr(ρE))² ≥ ε²` on the selected eigenvalues.
The witness section instantiates the hypothesis with an explicit band
projection. -/

/-- **Trace-dependent. Chebyshev, multiplied form, exact.** For every event
`Q` dominated by the squared deviation at scale `ε`, the mass of `Q` under
the product state obeys `ε² · mass ≤ (Tr(ρE) - Tr(ρE)²)/N` in the order of
`ℂ`. -/
theorem chebyshev_smul_bound (ρ E : Matrix (Fin d) (Fin d) ℂ)
    (hρ : IsState ρ) (hEev : IsEvent E) {N : ℕ} (hN : 0 < N) (ε : ℝ)
    {Q : Matrix (PowIdx d N) (PowIdx d N) ℂ}
    (hdom : ((deviationOp ρ E N * deviationOp ρ E N)
        - ((ε : ℂ) ^ 2) • Q).PosSemidef) :
    (ε : ℂ) ^ 2 * concentrationMass ρ N Q
      ≤ (bornWeight ρ E - bornWeight ρ E ^ 2) / N := by
  have h0 : 0 ≤ (kronPow ρ N * ((deviationOp ρ E N * deviationOp ρ E N)
      - ((ε : ℂ) ^ 2) • Q)).trace :=
    trace_mul_nonneg_of_posSemidef (kronPow_posSemidef hρ.1 N) hdom
  rw [Matrix.mul_sub, Matrix.trace_sub, Matrix.mul_smul, Matrix.trace_smul,
    smul_eq_mul, sub_nonneg] at h0
  calc (ε : ℂ) ^ 2 * concentrationMass ρ N Q
      = (ε : ℂ) ^ 2 * (kronPow ρ N * Q).trace := by
        rw [concentrationMass_def]
    _ ≤ (kronPow ρ N * (deviationOp ρ E N * deviationOp ρ E N)).trace := h0
    _ = (bornWeight ρ E - bornWeight ρ E ^ 2) / N := by
        have hv := variance_freqOp ρ E hρ.2 hEev.2 hN
        rwa [concentrationMass_def] at hv

/-- **Trace-dependent. Chebyshev, divided real form, exact.** The real mass
of every event dominated by the squared deviation at scale `ε ≠ 0` is at
most `(Tr(ρE) - Tr(ρE)²)/(N ε²)`: the exact variance ratio, decaying as
`1/N`. -/
theorem chebyshev_mass_bound (ρ E : Matrix (Fin d) (Fin d) ℂ)
    (hρ : IsState ρ) (hEev : IsEvent E) {N : ℕ} (hN : 0 < N) {ε : ℝ}
    (hε : ε ≠ 0) {Q : Matrix (PowIdx d N) (PowIdx d N) ℂ}
    (hdom : ((deviationOp ρ E N * deviationOp ρ E N)
        - ((ε : ℂ) ^ 2) • Q).PosSemidef) :
    (concentrationMass ρ N Q).re
      ≤ ((bornWeight ρ E).re - (bornWeight ρ E).re ^ 2) / (N * ε ^ 2) := by
  have hb := chebyshev_smul_bound ρ E hρ hEev hN ε hdom
  have hreal : (((bornWeight ρ E).re : ℝ) : ℂ) = bornWeight ρ E :=
    bornWeight_eq_re hρ.1.1 hEev.1
  rw [← hreal] at hb
  have hRHS : (((bornWeight ρ E).re : ℂ)
        - ((bornWeight ρ E).re : ℂ) ^ 2) / (N : ℂ)
      = ((((bornWeight ρ E).re - (bornWeight ρ E).re ^ 2)
          / (N : ℝ) : ℝ) : ℂ) := by
    push_cast
    ring
  have hLHS : (ε : ℂ) ^ 2 * concentrationMass ρ N Q
      = ((ε ^ 2 : ℝ) : ℂ) * concentrationMass ρ N Q := by
    push_cast
    ring
  rw [hRHS, hLHS] at hb
  have hre := (Complex.le_def.mp hb).1
  rw [Complex.re_ofReal_mul, Complex.ofReal_re] at hre
  have hε2 : (0 : ℝ) < ε ^ 2 :=
    (sq_nonneg ε).lt_of_ne' (pow_ne_zero 2 hε)
  have hdiv : (concentrationMass ρ N Q).re
      ≤ (((bornWeight ρ E).re - (bornWeight ρ E).re ^ 2) / (N : ℝ))
        / ε ^ 2 :=
    (le_div_iff₀ hε2).mpr (by rw [mul_comm]; exact hre)
  rwa [div_div] at hdiv

/-- The variance ratio of the Chebyshev bound tends to zero: the exact
finite bounds `(Tr(ρE) - Tr(ρE)²)/(N ε²)` form a null sequence in `N` for
every fixed band scale `ε`. -/
theorem varianceRatio_tendsto_zero (c ε : ℝ) :
    Filter.Tendsto (fun N : ℕ => c / (N * ε ^ 2)) Filter.atTop (nhds 0) := by
  have h : ∀ N : ℕ, c / ((N : ℝ) * ε ^ 2) = (c / ε ^ 2) / (N : ℝ) := by
    intro N
    rw [div_div, mul_comm]
  simp only [h]
  exact tendsto_const_div_atTop_nhds_zero_nat (c / ε ^ 2)

/-! ## Uniqueness of the stable frequency value -/

/-- **Trace-dependent. Uniqueness of the stable point, exact finite form.**
If a frequency assignment `p` differs from the Born weight `Tr(ρ E)`, then
for every declared threshold `θ > 0` and every
`N > (Tr(ρE) - Tr(ρE)²) / (θ ((p - Tr(ρE))/2)²)`, the explicit bound
displayed in the hypothesis `hbig`, every band event around `p` of
half-width `(p - Tr(ρE))/2`, encoded by the domination hypothesis, has
product-state mass below `θ`. A frequency assignment concentrated near `p`
is therefore inconsistent with the concentration receipts at this explicit
finite `N`: the Born weight is the unique value stable under them.

Circularity boundary: the mass in the conclusion is `concentrationMass`,
the Born pairing of the product state (`ConcentrationMassIsBornPairing`),
so this uniqueness is consistency inside the Born measure, never a
state-independent derivation of it. -/
theorem born_unique_stable_frequency (ρ E : Matrix (Fin d) (Fin d) ℂ)
    (hρ : IsState ρ) (hEev : IsEvent E) {p θ : ℝ}
    (hp : p ≠ (bornWeight ρ E).re) (hθ : 0 < θ) {N : ℕ}
    (hbig : ((bornWeight ρ E).re - (bornWeight ρ E).re ^ 2)
        / (θ * ((p - (bornWeight ρ E).re) / 2) ^ 2) < N)
    {Q : Matrix (PowIdx d N) (PowIdx d N) ℂ}
    (hdom : ((deviationOp ρ E N * deviationOp ρ E N)
        - ((((p - (bornWeight ρ E).re) / 2 : ℝ) : ℂ) ^ 2) • Q).PosSemidef) :
    (concentrationMass ρ N Q).re < θ := by
  set t : ℝ := (bornWeight ρ E).re with htdef
  set ε : ℝ := (p - t) / 2 with hεdef
  have hεne : ε ≠ 0 := div_ne_zero (sub_ne_zero.mpr hp) two_ne_zero
  have hε2 : (0 : ℝ) < ε ^ 2 :=
    (sq_nonneg ε).lt_of_ne' (pow_ne_zero 2 hεne)
  have ht0 : 0 ≤ t := bornWeight_re_nonneg hρ.1 hEev
  have ht1 : t ≤ 1 := bornWeight_re_le_one hρ hEev
  have hvar0 : 0 ≤ t - t ^ 2 := by nlinarith
  have hNR : (0 : ℝ) < N :=
    lt_of_le_of_lt (div_nonneg hvar0 (le_of_lt (mul_pos hθ hε2))) hbig
  have hN : 0 < N := by exact_mod_cast hNR
  have hmass := chebyshev_mass_bound ρ E hρ hEev hN hεne hdom
  have hstep : (t - t ^ 2) / ((N : ℝ) * ε ^ 2) < θ := by
    rw [div_lt_iff₀ (mul_pos hNR hε2)]
    have h2 : t - t ^ 2 < (N : ℝ) * (θ * ε ^ 2) :=
      (div_lt_iff₀ (mul_pos hθ hε2)).mp hbig
    linarith [h2,
      (by ring : θ * ((N : ℝ) * ε ^ 2) = (N : ℝ) * (θ * ε ^ 2))]
  exact lt_of_le_of_lt hmass hstep

/-- **The explicit refuting copy number.** For every `p ≠ Tr(ρE)` and every
threshold `θ > 0` there is an explicit `N`, the floor of the variance ratio
plus one, past which every band event around `p` has mass below `θ`. This
is the finite Finkelstein-Hartle statement: no frequency value other than
the Born weight survives the concentration receipts. -/
theorem exists_concentration_refuting_N (ρ E : Matrix (Fin d) (Fin d) ℂ)
    (hρ : IsState ρ) (hEev : IsEvent E) {p θ : ℝ}
    (hp : p ≠ (bornWeight ρ E).re) (hθ : 0 < θ) :
    ∃ N : ℕ, 0 < N ∧ ∀ Q : Matrix (PowIdx d N) (PowIdx d N) ℂ,
      ((deviationOp ρ E N * deviationOp ρ E N)
        - ((((p - (bornWeight ρ E).re) / 2 : ℝ) : ℂ) ^ 2) • Q).PosSemidef →
      (concentrationMass ρ N Q).re < θ := by
  refine ⟨⌊((bornWeight ρ E).re - (bornWeight ρ E).re ^ 2)
      / (θ * ((p - (bornWeight ρ E).re) / 2) ^ 2)⌋₊ + 1, Nat.succ_pos _,
    fun Q hdom => born_unique_stable_frequency ρ E hρ hEev hp hθ ?_ hdom⟩
  have h := Nat.lt_floor_add_one
    (((bornWeight ρ E).re - (bornWeight ρ E).re ^ 2)
      / (θ * ((p - (bornWeight ρ E).re) / 2) ^ 2))
  rwa [← Nat.cast_add_one] at h

/-! ## Witnesses at `d = 2`, `N = 2, 3`, with exact rational spectra

The witness state is the rational qubit density matrix `diag(1/3, 2/3)` and
the witness event is the projection onto the first basis vector, with Born
weight `1/3`. The frequency operators at `N = 2, 3` are computed entrywise
as explicit diagonal matrices: their spectra are the exact rational
frequency lists `{0, 1/2, 1}` and `{0, 1/3, 2/3, 1}`. The band projection
at frequency `p = 1` receives an explicit domination receipt at scale
`ε = 1/3 = (p - 1/3)/2`, and its exact masses `1/9` at `N = 2` and `1/27`
at `N = 3` display the concentration away from any value other than the
Born weight `1/3`. -/

/-- Witness state: the rational qubit density matrix `diag(1/3, 2/3)`. -/
noncomputable def stateW : Matrix (Fin 2) (Fin 2) ℂ :=
  !![((1 / 3 : ℝ) : ℂ), 0; 0, ((2 / 3 : ℝ) : ℂ)]

/-- Witness event: the projection onto the first basis vector. -/
def eventW : Matrix (Fin 2) (Fin 2) ℂ := !![1, 0; 0, 0]

theorem stateW_isState : IsState stateW := by
  constructor
  · have hdiag : stateW
        = Matrix.diagonal ![((1 / 3 : ℝ) : ℂ), ((2 / 3 : ℝ) : ℂ)] := by
      ext i j
      fin_cases i <;> fin_cases j <;> simp [stateW]
    rw [hdiag]
    refine Matrix.PosSemidef.diagonal ?_
    intro i
    fin_cases i
    · exact Complex.zero_le_real.mpr (by norm_num)
    · exact Complex.zero_le_real.mpr (by norm_num)
  · rw [Matrix.trace_fin_two]
    have h00 : stateW 0 0 = ((1 / 3 : ℝ) : ℂ) := rfl
    have h11 : stateW 1 1 = ((2 / 3 : ℝ) : ℂ) := rfl
    rw [h00, h11, ← Complex.ofReal_add]
    norm_num

theorem eventW_isEvent : IsEvent eventW := by
  constructor
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [eventW, Matrix.conjTranspose_apply]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [eventW, Matrix.mul_apply, Fin.sum_univ_two]

/-- The exact witness Born weight `Tr(ρ E) = 1/3`. -/
theorem bornWeight_stateW_eventW :
    bornWeight stateW eventW = ((1 / 3 : ℝ) : ℂ) := by
  rw [bornWeight, Matrix.trace_fin_two]
  simp [stateW, eventW, Matrix.mul_apply, Fin.sum_univ_two]

/-- **Exact mean at `N = 2`.** -/
theorem witness_mean_two :
    concentrationMass stateW 2 (freqOp eventW 2) = ((1 / 3 : ℝ) : ℂ) := by
  rw [concentrationMass_freqOp stateW eventW stateW_isState.2 (by norm_num),
    bornWeight_stateW_eventW]

/-- **Exact mean at `N = 3`.** -/
theorem witness_mean_three :
    concentrationMass stateW 3 (freqOp eventW 3) = ((1 / 3 : ℝ) : ℂ) := by
  rw [concentrationMass_freqOp stateW eventW stateW_isState.2 (by norm_num),
    bornWeight_stateW_eventW]

/-- **Exact variance at `N = 2`:** `(1/3 - 1/9)/2 = 1/9`. -/
theorem witness_variance_two :
    concentrationMass stateW 2
      (deviationOp stateW eventW 2 * deviationOp stateW eventW 2)
      = ((1 / 9 : ℝ) : ℂ) := by
  rw [variance_freqOp stateW eventW stateW_isState.2 eventW_isEvent.2
    (by norm_num), bornWeight_stateW_eventW]
  push_cast
  norm_num

/-- **Exact variance at `N = 3`:** `(1/3 - 1/9)/3 = 2/27`. -/
theorem witness_variance_three :
    concentrationMass stateW 3
      (deviationOp stateW eventW 3 * deviationOp stateW eventW 3)
      = ((2 / 27 : ℝ) : ℂ) := by
  rw [variance_freqOp stateW eventW stateW_isState.2 eventW_isEvent.2
    (by norm_num), bornWeight_stateW_eventW]
  push_cast
  norm_num

/-- **The exact rational spectrum at `N = 2`.** The frequency operator is
the explicit diagonal matrix of relative frequencies: eigenvalue `1` on the
all-zeros index, `1/2` on the two mixed indices, `0` on the all-ones index.
-/
theorem freqOp_eventW_two :
    freqOp eventW 2
      = Matrix.diagonal (fun x : Fin 2 × Fin 2 × PUnit =>
          ((((if x.1 = 0 then 1 else 0)
            + (if x.2.1 = 0 then 1 else 0)) / 2 : ℝ) : ℂ)) := by
  have h1 : (1 : Fin 2) = (0 : Fin 1).succ := rfl
  rw [freqOp, Fin.sum_univ_two, h1, embedAt_succ, embedAt_zero, embedAt_zero]
  ext ⟨a, b, ⟨⟩⟩ ⟨a', b', ⟨⟩⟩
  fin_cases a <;> fin_cases a' <;> fin_cases b <;> fin_cases b'
  all_goals simp [Matrix.smul_apply, Matrix.add_apply,
    Matrix.kroneckerMap_apply, eventW, Prod.ext_iff]
  all_goals norm_num

/-- **The exact rational spectrum at `N = 3`:** eigenvalues `k/3` with `k`
the number of tensor slots reading `0`. -/
theorem freqOp_eventW_three :
    freqOp eventW 3
      = Matrix.diagonal (fun x : Fin 2 × Fin 2 × Fin 2 × PUnit =>
          ((((if x.1 = 0 then 1 else 0) + (if x.2.1 = 0 then 1 else 0)
            + (if x.2.2.1 = 0 then 1 else 0)) / 3 : ℝ) : ℂ)) := by
  have h1 : (1 : Fin 3) = (0 : Fin 2).succ := rfl
  have h2 : (2 : Fin 3) = (1 : Fin 2).succ := rfl
  have h1' : (1 : Fin 2) = (0 : Fin 1).succ := rfl
  rw [freqOp, Fin.sum_univ_three, h1, h2, embedAt_succ, embedAt_succ, h1',
    embedAt_succ, embedAt_zero, embedAt_zero, embedAt_zero]
  ext ⟨a, b, c, ⟨⟩⟩ ⟨a', b', c', ⟨⟩⟩
  fin_cases a <;> fin_cases a' <;> fin_cases b <;> fin_cases b' <;>
    fin_cases c <;> fin_cases c'
  all_goals simp [Matrix.smul_apply, Matrix.add_apply,
    Matrix.kroneckerMap_apply, eventW, Prod.ext_iff]
  all_goals norm_num

/-- The `2`-copy witness product state, as an explicit rational diagonal. -/
theorem kronPow_stateW_two :
    kronPow stateW 2
      = Matrix.diagonal (fun x : Fin 2 × Fin 2 × PUnit =>
          (((if x.1 = 0 then (1 : ℝ) / 3 else 2 / 3)
            * (if x.2.1 = 0 then (1 : ℝ) / 3 else 2 / 3) : ℝ) : ℂ)) := by
  ext ⟨a, b, ⟨⟩⟩ ⟨a', b', ⟨⟩⟩
  fin_cases a <;> fin_cases a' <;> fin_cases b <;> fin_cases b' <;>
    simp [stateW, Matrix.kroneckerMap_apply, Prod.ext_iff]

/-- The `3`-copy witness product state, as an explicit rational diagonal. -/
theorem kronPow_stateW_three :
    kronPow stateW 3
      = Matrix.diagonal (fun x : Fin 2 × Fin 2 × Fin 2 × PUnit =>
          (((if x.1 = 0 then (1 : ℝ) / 3 else 2 / 3)
            * ((if x.2.1 = 0 then (1 : ℝ) / 3 else 2 / 3)
              * (if x.2.2.1 = 0 then (1 : ℝ) / 3 else 2 / 3)) : ℝ) : ℂ)) := by
  ext ⟨a, b, c, ⟨⟩⟩ ⟨a', b', c', ⟨⟩⟩
  fin_cases a <;> fin_cases a' <;> fin_cases b <;> fin_cases b' <;>
    fin_cases c <;> fin_cases c' <;>
    simp [stateW, Matrix.kroneckerMap_apply, Prod.ext_iff]

/-- The band event at frequency `p = 1` for `N = 2`: the spectral projection
of `freqOp eventW 2` onto its eigenvalue `1`. -/
def bandW2 : Matrix (PowIdx 2 2) (PowIdx 2 2) ℂ :=
  Matrix.diagonal (fun x : Fin 2 × Fin 2 × PUnit =>
    if x.1 = 0 ∧ x.2.1 = 0 then 1 else 0)

/-- The band event at frequency `p = 1` for `N = 3`. -/
def bandW3 : Matrix (PowIdx 2 3) (PowIdx 2 3) ℂ :=
  Matrix.diagonal (fun x : Fin 2 × Fin 2 × Fin 2 × PUnit =>
    if x.1 = 0 ∧ x.2.1 = 0 ∧ x.2.2.1 = 0 then 1 else 0)

/-- **Exact band mass at `N = 2`:** the product state gives the frequency-1
band the mass `(1/3)² = 1/9`, against the Born weight `1/3`. -/
theorem bandW2_mass : concentrationMass stateW 2 bandW2 = ((1 / 9 : ℝ) : ℂ) := by
  rw [concentrationMass_def, kronPow_stateW_two, bandW2,
    Matrix.diagonal_mul_diagonal, Matrix.trace_diagonal]
  simp [Fintype.sum_prod_type, Fin.sum_univ_two]
  norm_num

/-- **Exact band mass at `N = 3`:** the frequency-1 band mass decays to
`(1/3)³ = 1/27`: concentration away from the value `p = 1 ≠ 1/3`. -/
theorem bandW3_mass : concentrationMass stateW 3 bandW3 = ((1 / 27 : ℝ) : ℂ) := by
  rw [concentrationMass_def, kronPow_stateW_three, bandW3,
    Matrix.diagonal_mul_diagonal, Matrix.trace_diagonal]
  simp [Fintype.sum_prod_type, Fin.sum_univ_two]
  norm_num

/-- The deviation operator of the witness pair at `N = 2`, as an explicit
rational diagonal. -/
theorem deviationOp_W_two :
    deviationOp stateW eventW 2
      = Matrix.diagonal (fun x : Fin 2 × Fin 2 × PUnit =>
          ((((if x.1 = 0 then 1 else 0) + (if x.2.1 = 0 then 1 else 0)) / 2
            - 1 / 3 : ℝ) : ℂ)) := by
  rw [deviationOp, freqOp_eventW_two, bornWeight_stateW_eventW]
  ext ⟨a, b, ⟨⟩⟩ ⟨a', b', ⟨⟩⟩
  fin_cases a <;> fin_cases a' <;> fin_cases b <;> fin_cases b'
  all_goals simp [Matrix.sub_apply, Matrix.smul_apply, Prod.ext_iff]

/-- **The domination receipt at `N = 2`.** The frequency-1 band projection
is dominated by the squared deviation at scale `ε = 1/3 = (1 - 1/3)/2`:
the entrywise inequalities are `(2/3)² ≥ 1/9` on the band and
`(1/6)², (1/6)², (1/3)² ≥ 0` off it. This instantiates the Chebyshev
hypothesis with an exact rational witness. -/
theorem bandW2_dominated :
    ((deviationOp stateW eventW 2 * deviationOp stateW eventW 2)
      - ((((1 : ℝ) / 3 : ℝ) : ℂ) ^ 2) • bandW2).PosSemidef := by
  have hδ : (deviationOp stateW eventW 2 * deviationOp stateW eventW 2)
      - ((((1 : ℝ) / 3 : ℝ) : ℂ) ^ 2) • bandW2
      = Matrix.diagonal (fun x : Fin 2 × Fin 2 × PUnit =>
          (((((if x.1 = 0 then 1 else 0) + (if x.2.1 = 0 then 1 else 0)) / 2
              - 1 / 3) ^ 2
            - (1 / 9) * (if x.1 = 0 ∧ x.2.1 = 0 then 1 else 0) : ℝ) : ℂ)) := by
    rw [deviationOp_W_two, bandW2, Matrix.diagonal_mul_diagonal,
      ← Matrix.diagonal_smul, Matrix.diagonal_sub]
    congr 1
    funext x
    rcases x with ⟨a, b, ⟨⟩⟩
    fin_cases a <;> fin_cases b <;>
      · simp only [Pi.smul_apply, smul_eq_mul]
        push_cast
        norm_num
  rw [hδ]
  refine Matrix.PosSemidef.diagonal ?_
  intro x
  rcases x with ⟨a, b, ⟨⟩⟩
  fin_cases a <;> fin_cases b <;>
    exact Complex.zero_le_real.mpr (by norm_num)

-- Axiom audit: each must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms trace_kronPow
#print axioms kronPow_posSemidef
#print axioms trace_mul_nonneg_of_posSemidef
#print axioms embedAt_mul_embedAt_comm
#print axioms trace_kronPow_mul_embedAt
#print axioms trace_kronPow_mul_embedAt_pair
#print axioms freqOp_isHermitian
#print axioms concentrationMass_isBornPairing
#print axioms concentrationMass_freqOp
#print axioms concentrationMass_freqOp_sq
#print axioms variance_freqOp
#print axioms chebyshev_smul_bound
#print axioms chebyshev_mass_bound
#print axioms varianceRatio_tendsto_zero
#print axioms born_unique_stable_frequency
#print axioms exists_concentration_refuting_N
#print axioms freqOp_eventW_two
#print axioms freqOp_eventW_three
#print axioms witness_variance_two
#print axioms witness_variance_three
#print axioms bandW2_mass
#print axioms bandW3_mass
#print axioms bandW2_dominated

end EventAlgebra.FrequencyConcentration
