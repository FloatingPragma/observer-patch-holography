import Mathlib

namespace OPH.Thermodynamics

/-!
# Collar-matrix realization on a mixing recurrent chain (B12 receipt 3)

This module mirrors, as exact rational literals, the mixing recurrent chain
extracted from the locally hash-pinned bounded source run
`runs/b12_prereg_16k_20260806` of the simulator (seed `20260806`, engine
commit `b39b78faf894894ebe573571e0902ccfaaeac32a`).  The payload with the full
extraction and its exact verification lives at
`oph-physics-sim/docs/B12_MIXING_CHAIN_PAYLOAD.json`
(`scripts/extract_b12_mixing_chain.py`); the pinned transition-clock report is
`finite_repair_transition_matrix_report.json`
(sha256 `e38bb28475fb7b34127b864abd20414efc920d8ef9ddbdf5c78612874ae37f2d`)
over `observer_views.jsonl`
(sha256 `3da6f04d770b81b49a082ad1e6ecb93c611517dd85ecaccd7cc89077fcd6c0ca`).

The raw 26-state packet-quotient chain of that run (fields
`checkpoint_class`, `stable_flag`, `s3_sector_class`, `repair_load_bucket`)
is reducible; its unique closed recurrent class consists of the two committed
record-checkpoint states

* global state 18: `checkpoint_class = 3` (committed, stable commit),
* global state 19: `checkpoint_class = 2` (committed, unstable commit),

with unweighted integer transition counts `[[1324, 107], [5, 503]]` and row
sums `[1431, 508]`.  Issue #688 admits an explicitly restricted recurrent
chain, and the restriction is certified exactly on the simulator side:
closure, shared positivity support with the pinned weighted count matrix,
irreducibility, aperiodicity, and exact stationarity all hold over exact
rationals.

Receipt structure.  The integer layer (counts, row sums, strict positivity,
the cleared-denominator stationarity identity, the cleared-denominator trace
identity, and the nonconstant protected labelling) is kernel-decided on the
literals.  The rational layer (row stochasticity, entrywise positivity, the
stationary law, and the exact second eigenvalue) is elaborated from the same
literals by `norm_num`; rational division is kept out of kernel reduction.

Claim boundary.  Realization at the representation level over extracted
literals.  Every entry of the chain is strictly positive, so the chain is a
primitive stochastic matrix: irreducible, aperiodic, and mixing to its unique
stationary law.  The protected labelling carried by this chain is the
committed record-checkpoint class; the `record_family` projection of the
state-side protected record (`record_signature`, 32 classes, preserved by the
idempotent conditional-resampling kernel of the same run) is constant on the
transition alphabet of this run.  The state-side and transition-side objects
share one pinned source bundle and provenance record (run id, seed, config
hash, commit, and sha256 digests recorded in the payload); that does not
identify a common state/transition probability law. Source-side common-law and
physical-collar construction remains under #739; energy-clock calibration is
the empirical import PR-15.
-/

/-- The exact restricted recurrent chain: row-normalized unweighted integer
transition counts on the unique closed recurrent class `{18, 19}` of the
locally hash-pinned run's 26-state quotient chain. -/
def mixingChain : Fin 2 → Fin 2 → ℚ := fun i =>
  ![![1324/1431, 107/1431], ![5/508, 503/508]] i

/-- The unweighted integer transition counts on the recurrent class. -/
def mixingChainCounts : Fin 2 → Fin 2 → ℕ := fun i =>
  ![![1324, 107], ![5, 503]] i

/-- Row sums of the integer counts. -/
def mixingChainRowSums : Fin 2 → ℕ := ![1431, 508]

/-- The exact stationary law of `mixingChain`. -/
def mixingChainStationary : Fin 2 → ℚ := ![7155/61511, 54356/61511]

/-- Numerators of the stationary law over the common denominator `61511`. -/
def mixingChainStationaryNum : Fin 2 → ℕ := ![7155, 54356]

/-- Cofactors clearing the row denominators: `726948 / rowSum i`. -/
def rowSumCofactor : Fin 2 → ℕ := ![508, 1431]

/-- Protected record labelling of the recurrent states: the committed
record-checkpoint class (`checkpoint_class` of the packet quotient;
`3` = committed stable, `2` = committed unstable). -/
def protectedRecordLabel : Fin 2 → ℕ := ![3, 2]

/-! ## Kernel-decided integer receipts -/

/-- Row sums of the integer counts are as declared. -/
theorem mixingChainCounts_row_sums :
    ∀ i, (∑ j, mixingChainCounts i j) = mixingChainRowSums i := by
  decide

/-- Every integer transition count on the recurrent class is strictly
positive: a one-step primitivity witness.  The positive-entry digraph is
complete, hence strongly connected (irreducibility), and the positive
diagonal forces cycle-length gcd one (aperiodicity). -/
theorem mixingChainCounts_pos : ∀ i j, 0 < mixingChainCounts i j := by
  decide

/-- The declared cofactors clear the row denominators against the common
denominator `726948 = 1431 * 508`. -/
theorem rowSumCofactor_spec :
    ∀ i, mixingChainRowSums i * rowSumCofactor i = 726948 := by
  decide

/-- Cleared-denominator stationarity: the identity
`∑ i, π_num i * C i j * cofactor i = π_num j * 726948` over ℕ is the
stationarity equation of the normalized chain with every denominator
cleared. -/
theorem mixingChainStationary_integer_identity :
    ∀ j, (∑ i, mixingChainStationaryNum i * mixingChainCounts i j *
      rowSumCofactor i) = mixingChainStationaryNum j * 726948 := by
  decide

/-- The stationary numerators sum to the common denominator `61511`. -/
theorem mixingChainStationaryNum_sum :
    (∑ i, mixingChainStationaryNum i) = 61511 := by
  decide

/-- Cleared-denominator trace identity for the second eigenvalue: over ℕ,
`1324 * 508 + 503 * 1431 = 726948 + 665437`, the integer form of
`trace - 1 = 665437 / 726948`. -/
theorem mixingChain_trace_integer_identity :
    mixingChainCounts 0 0 * rowSumCofactor 0 +
      mixingChainCounts 1 1 * rowSumCofactor 1 = 726948 + 665437 := by
  decide

/-- The protected record labelling is nonconstant across the recurrent
states. -/
theorem protectedRecordLabel_nonconstant :
    protectedRecordLabel 0 ≠ protectedRecordLabel 1 := by
  decide

/-! ## Rational receipts over the same literals -/

/-- The chain is the row normalization of the integer counts. -/
theorem mixingChain_eq_normalized_counts :
    ∀ i j, mixingChain i j =
      (mixingChainCounts i j : ℚ) / (mixingChainRowSums i : ℚ) := by
  intro i j
  fin_cases i <;> fin_cases j <;>
    norm_num [mixingChain, mixingChainCounts, mixingChainRowSums]

/-- Every entry of the chain is strictly positive. -/
theorem mixingChain_entries_pos : ∀ i j, 0 < mixingChain i j := by
  intro i j
  fin_cases i <;> fin_cases j <;> norm_num [mixingChain]

/-- The chain is row stochastic. -/
theorem mixingChain_row_stochastic : ∀ i, (∑ j, mixingChain i j) = 1 := by
  intro i
  fin_cases i <;> norm_num [mixingChain, Fin.sum_univ_succ]

/-- The diagonal is strictly positive, an explicit aperiodicity witness. -/
theorem mixingChain_diagonal_pos : ∀ i, 0 < mixingChain i i := by
  intro i
  exact mixingChain_entries_pos i i

/-- The stationary law is strictly positive. -/
theorem mixingChainStationary_pos : ∀ i, 0 < mixingChainStationary i := by
  intro i
  fin_cases i <;> norm_num [mixingChainStationary]

/-- The stationary law is normalized. -/
theorem mixingChainStationary_sum : (∑ i, mixingChainStationary i) = 1 := by
  norm_num [mixingChainStationary, Fin.sum_univ_succ]

/-- Exact stationarity: the stationary law is fixed by the chain. -/
theorem mixingChainStationary_stationary :
    ∀ j, (∑ i, mixingChainStationary i * mixingChain i j) =
      mixingChainStationary j := by
  intro j
  fin_cases j <;>
    norm_num [mixingChainStationary, mixingChain, Fin.sum_univ_succ]

/-- Two-state trace identity: the second eigenvalue of a `2 × 2` row
stochastic matrix is `trace - 1`, here exactly `665437/726948`, strictly
inside the open unit interval, so the exact spectral gap is
`61511/726948`. -/
theorem mixingChain_second_eigenvalue :
    mixingChain 0 0 + mixingChain 1 1 - 1 = 665437/726948 ∧
      (0 : ℚ) < 665437/726948 ∧ (665437 : ℚ)/726948 < 1 := by
  refine ⟨?_, by norm_num, by norm_num⟩
  norm_num [mixingChain]

/-- **Transition-side matrix realization receipt.**  The exact chain
extracted from the locally hash-pinned source run is row stochastic with every
entry strictly positive (a one-step primitivity witness, so the chain is an
irreducible aperiodic recurrent chain mixing to its stationary law), it
fixes a strictly positive exact stationary law, and its states carry a
nonconstant protected record labelling.  The realization is at the
representation level over committed extracted literals; it does not identify the
state-side probability reference. Source-side physical-collar and
common-reference construction remains under #739; energy-clock calibration is
the empirical import PR-15. -/
theorem mixingChainRealization_receipt :
    (∀ i, (∑ j, mixingChain i j) = 1) ∧
      (∀ i j, 0 < mixingChain i j) ∧
      (∀ j, (∑ i, mixingChainStationary i * mixingChain i j) =
        mixingChainStationary j) ∧
      (∀ i, 0 < mixingChainStationary i) ∧
      (∑ i, mixingChainStationary i) = 1 ∧
      protectedRecordLabel 0 ≠ protectedRecordLabel 1 :=
  ⟨mixingChain_row_stochastic, mixingChain_entries_pos,
    mixingChainStationary_stationary, mixingChainStationary_pos,
    mixingChainStationary_sum, protectedRecordLabel_nonconstant⟩

#print axioms mixingChainRealization_receipt
#print axioms OPH.Thermodynamics.mixingChainCounts_row_sums
#print axioms OPH.Thermodynamics.mixingChainStationary_integer_identity
#print axioms OPH.Thermodynamics.mixingChain_trace_integer_identity

end OPH.Thermodynamics
