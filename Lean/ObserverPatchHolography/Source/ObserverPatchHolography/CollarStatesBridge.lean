import Mathlib
import ObserverPatchHolography.CollarStatesT1
import EventAlgebra.CenterExpectation

/-!
# Issue #544 bridge receipts: the flux expectation against the
# `EventAlgebra` pinching library

This module records, admission-free, the precise relationship between the
T1 flux conditional expectation `EfluxL` (from `CollarStatesT1.lean`) and
the `EventAlgebra` library's commutant-valued pinching
`EventAlgebra.pinchingExpectation` (from `EventAlgebra/CenterExpectation.lean`).

## Why a bridge and not a replacement

Two independent obstructions rule out re-expressing `EfluxL` *as* a
pinching, so the T0/T1/T2 payloads are kept verbatim and the relationship
is carried here as additive receipts:

1. **Commutant vs centre.** For any `ProjectivePartition` the pinching
   `X ↦ ∑ i, Pᵢ X Pᵢ` lands in the *commutant* `⊕ᵢ M_{rank Pᵢ}(ℂ)`, of
   dimension `∑ rankᵢ² ≥ 4` on a four-dimensional carrier. `EfluxL` lands
   in the two-dimensional flux *centre* `span{1, uuC}`. No partition makes
   the two maps equal; `EfluxL` is the blockwise normalised-trace
   compression of the pinching for the partition `{quP 0, quP 1}`, one
   step finer than anything the library defines.
2. **Index types.** The library is stated over `Matrix (Fin n) (Fin n) ℂ`
   while the collar algebra is indexed by `Fin 2 × Fin 2`, so every
   contact point must pass through the reindexing algebra isomorphism
   `collarReindex` along `finProdFinEquiv`.

## The receipts

* `fluxPartition` (**R1**): the `uuC` half-projections `quP` genuinely
  instantiate the library's `ProjectivePartition` after reindexing — the
  fields are discharged by the already-frisked T0/T1 lemmas
  (`quP_isHermitian`, `quP_idem`, `quP_mul`, `sum_quP`).
* `collarReindex_EfluxL_eq_compressed_pinching` (**R2a**): transported
  along `collarReindex`, `EfluxL` *is* the blockwise normalised-trace
  compression of `pinchingExpectation fluxPartition`; the block traces
  are computed with the library's `pinchingExpectation_mul_proj` and
  `trace_sandwich`.
* `EfluxL_absorbs_pinching` (**R2b**): `EfluxL` factors through the
  library pinching — pinching first changes nothing. Proved with the
  library's defining conditional-expectation property
  `trace_pinchingExpectation_mul_central`.

Nothing in this file is load-bearing for the T0/T1/T2 payloads; deleting
it changes no payload statement. It exists so the state-side no-gos cite
the audited `EventAlgebra` development instead of appearing to reinvent
it.
-/

namespace OPH

open Matrix

/-- The collar algebra reindexed onto `Fin (2 * 2)`, as a ℂ-algebra
isomorphism. Every contact point with the `EventAlgebra` library (which is
stated over `Matrix (Fin n) (Fin n) ℂ`) passes through this map. -/
noncomputable def collarReindex :
    CollarC ≃ₐ[ℂ] Matrix (Fin (2 * 2)) (Fin (2 * 2)) ℂ :=
  Matrix.reindexAlgEquiv ℂ ℂ finProdFinEquiv

@[simp] theorem collarReindex_apply (m : CollarC) :
    collarReindex m = Matrix.reindex finProdFinEquiv finProdFinEquiv m :=
  rfl

/-- Reindexing preserves the trace. -/
theorem trace_collarReindex (m : CollarC) :
    (collarReindex m).trace = m.trace := by
  rw [collarReindex_apply, Matrix.reindex_apply, Matrix.trace, Matrix.trace]
  exact Equiv.sum_comp finProdFinEquiv.symm fun p => m.diag p

/-- Reindexing preserves hermiticity. -/
theorem isHermitian_collarReindex {m : CollarC} (hm : m.IsHermitian) :
    (collarReindex m).IsHermitian := by
  rw [collarReindex_apply, Matrix.reindex_apply]
  exact hm.submatrix _

/-- **R1.** The `uuC` half-projections `quP`, transported along
`collarReindex`, form a `ProjectivePartition` in the sense of the
`EventAlgebra` library: each is an event (`quP_isHermitian`, `quP_idem`),
distinct members are orthogonal (`quP_mul`), and they resolve the identity
(`sum_quP`). -/
noncomputable def fluxPartition : EventAlgebra.ProjectivePartition (2 * 2) 2 where
  proj s := collarReindex (quP s)
  isEvent s :=
    ⟨isHermitian_collarReindex (quP_isHermitian s), by
      rw [← map_mul, quP_idem]⟩
  orthogonal s t hst := by
    rw [← map_mul, quP_mul, if_neg hst, map_zero]
  complete := by
    rw [← map_sum, sum_quP, map_one]

@[simp] theorem fluxPartition_proj (s : Fin 2) :
    fluxPartition.proj s = collarReindex (quP s) :=
  rfl

/-- The reindexed generator `uuC` is the signed sum of the partition
members — hence central for the partition. -/
theorem collarReindex_uuC_eq :
    collarReindex uuC = fluxPartition.proj 0 - fluxPartition.proj 1 := by
  rw [fluxPartition_proj, fluxPartition_proj, ← map_sub]
  congr 1
  rw [quP, quP, sgnR_zero, sgnR_one, Complex.ofReal_one, Complex.ofReal_neg]
  match_scalars <;> ring

/-- Members of the flux partition commute with each other, so each is
central for the partition in the library's sense. -/
theorem fluxPartition_proj_comm (s t : Fin 2) :
    fluxPartition.proj s * fluxPartition.proj t
      = fluxPartition.proj t * fluxPartition.proj s :=
  fluxPartition.proj_commute s t

/-- **R2a.** Transported along `collarReindex`, the flux expectation
`EfluxL` is exactly the blockwise normalised-trace compression of the
library pinching for `fluxPartition`: each block trace of the pinching
(computed with the library's `pinchingExpectation_mul_proj` and
`trace_sandwich`) is divided by the block rank `2` and spread back over
the block projector. This is the precise sense in which `EfluxL` is one
step finer than — not an instance of — `pinchingExpectation`. -/
theorem collarReindex_EfluxL_eq_compressed_pinching (m : CollarC) :
    collarReindex (EfluxL m) =
      ∑ s : Fin 2,
        ((EventAlgebra.pinchingExpectation fluxPartition (collarReindex m)
            * fluxPartition.proj s).trace / 2) • fluxPartition.proj s := by
  have hblock : ∀ s : Fin 2,
      (EventAlgebra.pinchingExpectation fluxPartition (collarReindex m)
          * fluxPartition.proj s).trace = (quP s * m).trace := by
    intro s
    rw [EventAlgebra.pinchingExpectation_mul_proj,
      EventAlgebra.trace_sandwich (fluxPartition.isEvent s).2,
      EventAlgebra.bornWeight, fluxPartition_proj, ← map_mul,
      trace_collarReindex, Matrix.trace_mul_comm]
  calc collarReindex (EfluxL m)
      = collarReindex (∑ s : Fin 2, ((quP s * m).trace / 2) • quP s) := by
        rw [EfluxL_spectral]
    _ = ∑ s : Fin 2, ((quP s * m).trace / 2) • collarReindex (quP s) := by
        rw [map_sum]
        exact Finset.sum_congr rfl fun s _ => map_smul _ _ _
    _ = ∑ s : Fin 2,
          ((EventAlgebra.pinchingExpectation fluxPartition (collarReindex m)
              * fluxPartition.proj s).trace / 2) • fluxPartition.proj s := by
        exact Finset.sum_congr rfl fun s _ => by rw [hblock, fluxPartition_proj]

/-- The identity is central for any partition — trace compatibility
specialised to `C = 1`. -/
private theorem pinching_trace_eq (X : Matrix (Fin (2 * 2)) (Fin (2 * 2)) ℂ) :
    (EventAlgebra.pinchingExpectation fluxPartition X).trace = X.trace := by
  have h := EventAlgebra.trace_pinchingExpectation_mul_central fluxPartition X 1
    (fun i => by rw [one_mul, mul_one])
  rwa [mul_one, mul_one] at h

/-- **R2b.** `EfluxL` absorbs the library pinching: applying
`pinchingExpectation fluxPartition` first (transported back along
`collarReindex.symm`) changes nothing. `EfluxL` reads only the two trace
pairings `Tr(m)` and `Tr(uuC · m)`, and both are invariant under the
pinching by the library's defining conditional-expectation property
`trace_pinchingExpectation_mul_central` (tested against the central
elements `1` and `collarReindex uuC`). -/
theorem EfluxL_absorbs_pinching (m : CollarC) :
    EfluxL (collarReindex.symm
      (EventAlgebra.pinchingExpectation fluxPartition (collarReindex m)))
      = EfluxL m := by
  set x : CollarC := collarReindex.symm
    (EventAlgebra.pinchingExpectation fluxPartition (collarReindex m)) with hx
  have hxr : collarReindex x
      = EventAlgebra.pinchingExpectation fluxPartition (collarReindex m) := by
    rw [hx, AlgEquiv.apply_symm_apply]
  have htrace : x.trace = m.trace := by
    rw [← trace_collarReindex x, hxr, pinching_trace_eq, trace_collarReindex]
  have huuC_central : ∀ i, collarReindex uuC * fluxPartition.proj i
      = fluxPartition.proj i * collarReindex uuC := by
    intro i
    rw [collarReindex_uuC_eq, sub_mul, mul_sub,
      fluxPartition_proj_comm 0 i, fluxPartition_proj_comm 1 i]
  have huu : (uuC * x).trace = (uuC * m).trace := by
    have h := EventAlgebra.trace_pinchingExpectation_mul_central fluxPartition
      (collarReindex m) (collarReindex uuC) huuC_central
    rw [← hxr] at h
    calc (uuC * x).trace
        = (x * uuC).trace := Matrix.trace_mul_comm uuC x
      _ = (collarReindex (x * uuC)).trace := (trace_collarReindex _).symm
      _ = (collarReindex x * collarReindex uuC).trace := by rw [map_mul]
      _ = (collarReindex m * collarReindex uuC).trace := h
      _ = (collarReindex (m * uuC)).trace := by rw [map_mul]
      _ = (m * uuC).trace := trace_collarReindex _
      _ = (uuC * m).trace := Matrix.trace_mul_comm m uuC
  rw [EfluxL_apply, EfluxL_apply, htrace, huu]

/-! ## Axiom audit

Expected footprint for every entry: `[propext, Classical.choice,
Quot.sound]`.  No `sorry`, no `native_decide`, no extra axioms. -/

#print axioms fluxPartition
#print axioms collarReindex_EfluxL_eq_compressed_pinching
#print axioms EfluxL_absorbs_pinching

end OPH
