import Mathlib
import A5Commutant

open scoped BigOperators

namespace OPH.A5ResponseWordAlgebra

open OPH.PortFrameGram (adj antipode)
open OPH.A5Commutant (dOrb aOrb nOrb pOrb Equivariant commutant_decomposition
  orbitals_independent equivariant_dOrb equivariant_aOrb equivariant_nOrb
  equivariant_pOrb)

/-! # The response-word algebra of the registered adjacency, exactly

Issue #705, boundary item one.  The registered adjacency recurrence and its
polynomial antipode readback generate an exactly four-dimensional commutative
response-word algebra over `ℚ`.  This file carries the arithmetic out on the
explicit twelve-port matrices rather than assuming it:

* `antipode_readback`: the antipode orbital is the cubic polynomial
  `P = (A³ - 4·A² - 5·A + 10·I) / 10` in the adjacency matrix, entry by
  entry.  The integral form of this identity is already on record as
  `Screen/A5IncidenceResponse.lean` (`antipode_polynomial`); here it is
  rederived over `ℚ` from the orbital product table as the input to the
  closure argument (`mulZ_*` kernel checks feed every identity through an
  exact `ℤ → ℚ` cast; no floating point and no `native_decide` anywhere).
* `adjacency_quartic`: the registered adjacency satisfies the quartic
  `A⁴ = 4·A³ + 10·A² - 20·A - 25·I`, so every response word reduces to a
  rational combination of `I, A, A², A³`.
* `word_in_power_span`, `powers_independent`: every word in the generators
  `{A, P}` lies in the span of `I, A, A², A³`, and those four matrices are
  `ℚ`-linearly independent — the response-word algebra has dimension exactly
  four, both directions.
* `word_mul_comm`, `no_nonzero_word_commutator`: the algebra is commutative;
  every commutator of response words vanishes identically.
* `equivariant_iff_word_span`: the word span is exactly the equivariant
  commutant of `Screen/A5Commutant.lean` — nothing outside the four-term
  span survives equivariance, and every word is equivariant.
* `no_twelve_independent_span_members`: no twelve members of the word span
  are `ℚ`-linearly independent, so the registered response cannot supply
  twelve independent current generators (issue #705, boundary item two,
  generator half).

NEGATIVE CONTROLS.

* `broken_readback_fails`, `broken_quartic_fails`: deleting the single edge
  `{0,1}` from the registered incidence breaks both load-bearing polynomial
  identities at explicit entries, so the identities test the registered
  incidence rather than holding vacuously.
* `elementary_commutator_nonzero`, `comparison_commutator_nonzero`,
  `comparison_not_in_span`: the elementary matrix `E₀₁` has a nonzero
  commutator with the adjacency and lies outside the word span — the zero
  commutator on the word algebra is a property of the algebra, not of all
  twelve-port matrices.
* `cube_not_in_three_span`, `antipode_not_quadratic`: `A³` is outside the
  span of `I, A, A²` and the antipode readback admits no polynomial of
  degree at most two — dimension four and the cubic readback are sharp.

BOUNDARY.  This file is exact coefficient algebra on the registered
twelve-port incidence.  It does not construct a compact current bracket,
does not touch the `u(12)` comparison rank or the fourteen-dimensional
alternating equivariant space (issue #705, boundary items three and four),
and carries no physical current, holonomy, or laboratory attachment. -/

/-! ## Integer layer: the registered matrices and their kernel-checked
products

All matrix identities are established first over `ℤ` by kernel `decide` and
then transported to `ℚ` through the exact cast; `ℚ` arithmetic is never
kernel-evaluated. -/

/-- The identity matrix over `ℤ`. -/
def dZ (i j : Fin 12) : ℤ := if i = j then 1 else 0

/-- The registered adjacency matrix over `ℤ`. -/
def aZ (i j : Fin 12) : ℤ := if adj i j then 1 else 0

/-- The distance-two orbital matrix over `ℤ`. -/
def nZ (i j : Fin 12) : ℤ :=
  if i ≠ j ∧ adj i j = false ∧ j ≠ antipode i then 1 else 0

/-- The antipode permutation matrix over `ℤ`. -/
def pZ (i j : Fin 12) : ℤ := if j = antipode i then 1 else 0

/-- One entry of an integer matrix product. -/
def mulZ (M N : Fin 12 → Fin 12 → ℤ) (i j : Fin 12) : ℤ :=
  ∑ k, M i k * N k j

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- `D` is a two-sided unit on the orbital matrices, kernel-checked. -/
theorem mulZ_dd : ∀ i j : Fin 12, mulZ dZ dZ i j = dZ i j := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
theorem mulZ_da : ∀ i j : Fin 12, mulZ dZ aZ i j = aZ i j := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
theorem mulZ_dn : ∀ i j : Fin 12, mulZ dZ nZ i j = nZ i j := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
theorem mulZ_dp : ∀ i j : Fin 12, mulZ dZ pZ i j = pZ i j := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
theorem mulZ_ad : ∀ i j : Fin 12, mulZ aZ dZ i j = aZ i j := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
theorem mulZ_nd : ∀ i j : Fin 12, mulZ nZ dZ i j = nZ i j := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
theorem mulZ_pd : ∀ i j : Fin 12, mulZ pZ dZ i j = pZ i j := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- `A² = 5·D + 2·A + 2·N`: the registered adjacency squared, kernel-checked. -/
theorem mulZ_aa :
    ∀ i j : Fin 12, mulZ aZ aZ i j = 5 * dZ i j + 2 * aZ i j + 2 * nZ i j := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- `A·N = 2·A + 2·N + 5·P`, kernel-checked. -/
theorem mulZ_an :
    ∀ i j : Fin 12, mulZ aZ nZ i j = 2 * aZ i j + 2 * nZ i j + 5 * pZ i j := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- `A·P = N`, kernel-checked. -/
theorem mulZ_ap : ∀ i j : Fin 12, mulZ aZ pZ i j = nZ i j := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- `N·A = 2·A + 2·N + 5·P`, kernel-checked. -/
theorem mulZ_na :
    ∀ i j : Fin 12, mulZ nZ aZ i j = 2 * aZ i j + 2 * nZ i j + 5 * pZ i j := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- `N² = 5·D + 2·A + 2·N`: the distance-two graph is again icosahedral,
kernel-checked. -/
theorem mulZ_nn :
    ∀ i j : Fin 12, mulZ nZ nZ i j = 5 * dZ i j + 2 * aZ i j + 2 * nZ i j := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- `N·P = A`, kernel-checked. -/
theorem mulZ_np : ∀ i j : Fin 12, mulZ nZ pZ i j = aZ i j := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- `P·A = N`, kernel-checked. -/
theorem mulZ_pa : ∀ i j : Fin 12, mulZ pZ aZ i j = nZ i j := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- `P·N = A`, kernel-checked. -/
theorem mulZ_pn : ∀ i j : Fin 12, mulZ pZ nZ i j = aZ i j := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- `P² = D`: the antipode is an involution, kernel-checked at matrix level. -/
theorem mulZ_pp : ∀ i j : Fin 12, mulZ pZ pZ i j = dZ i j := by decide

/-! ## The exact cast to `ℚ` -/

/-- Entrywise cast of an integer matrix to `ℚ`. -/
def castM (M : Fin 12 → Fin 12 → ℤ) : Fin 12 → Fin 12 → ℚ :=
  fun i j => (M i j : ℚ)

/-- One entry of a rational matrix product. -/
def mulQ (M N : Fin 12 → Fin 12 → ℚ) (i j : Fin 12) : ℚ :=
  ∑ k, M i k * N k j

/-- The cast is multiplicative entry by entry: rational products of cast
matrices are casts of the kernel-checked integer products. -/
theorem castM_mul (M N : Fin 12 → Fin 12 → ℤ) (i j : Fin 12) :
    mulQ (castM M) (castM N) i j = castM (mulZ M N) i j := by
  unfold mulQ mulZ castM
  push_cast
  rfl

theorem dOrb_cast : dOrb = castM dZ := by
  funext i j
  by_cases h : i = j <;> simp [OPH.A5Commutant.dOrb, dZ, castM, h]

theorem aOrb_cast : aOrb = castM aZ := by
  funext i j
  by_cases h : adj i j = true <;> simp [OPH.A5Commutant.aOrb, aZ, castM, h]

theorem nOrb_cast : nOrb = castM nZ := by
  funext i j
  by_cases h : i ≠ j ∧ adj i j = false ∧ j ≠ antipode i <;>
    simp [OPH.A5Commutant.nOrb, nZ, castM, h]

theorem pOrb_cast : pOrb = castM pZ := by
  funext i j
  by_cases h : j = antipode i <;> simp [OPH.A5Commutant.pOrb, pZ, castM, h]

/-! ## The rational orbital product table -/

theorem mulQ_dd (i j : Fin 12) : mulQ dOrb dOrb i j = dOrb i j := by
  rw [dOrb_cast, castM_mul]; unfold castM; rw [mulZ_dd i j]

theorem mulQ_da (i j : Fin 12) : mulQ dOrb aOrb i j = aOrb i j := by
  rw [dOrb_cast, aOrb_cast, castM_mul]; unfold castM; rw [mulZ_da i j]

theorem mulQ_dn (i j : Fin 12) : mulQ dOrb nOrb i j = nOrb i j := by
  rw [dOrb_cast, nOrb_cast, castM_mul]; unfold castM; rw [mulZ_dn i j]

theorem mulQ_dp (i j : Fin 12) : mulQ dOrb pOrb i j = pOrb i j := by
  rw [dOrb_cast, pOrb_cast, castM_mul]; unfold castM; rw [mulZ_dp i j]

theorem mulQ_ad (i j : Fin 12) : mulQ aOrb dOrb i j = aOrb i j := by
  rw [dOrb_cast, aOrb_cast, castM_mul]; unfold castM; rw [mulZ_ad i j]

theorem mulQ_nd (i j : Fin 12) : mulQ nOrb dOrb i j = nOrb i j := by
  rw [dOrb_cast, nOrb_cast, castM_mul]; unfold castM; rw [mulZ_nd i j]

theorem mulQ_pd (i j : Fin 12) : mulQ pOrb dOrb i j = pOrb i j := by
  rw [dOrb_cast, pOrb_cast, castM_mul]; unfold castM; rw [mulZ_pd i j]

theorem mulQ_aa (i j : Fin 12) :
    mulQ aOrb aOrb i j = 5 * dOrb i j + 2 * aOrb i j + 2 * nOrb i j := by
  rw [dOrb_cast, aOrb_cast, nOrb_cast, castM_mul]
  unfold castM
  rw [mulZ_aa i j]; push_cast; ring

theorem mulQ_an (i j : Fin 12) :
    mulQ aOrb nOrb i j = 2 * aOrb i j + 2 * nOrb i j + 5 * pOrb i j := by
  rw [aOrb_cast, nOrb_cast, pOrb_cast, castM_mul]
  unfold castM
  rw [mulZ_an i j]; push_cast; ring

theorem mulQ_ap (i j : Fin 12) : mulQ aOrb pOrb i j = nOrb i j := by
  rw [aOrb_cast, nOrb_cast, pOrb_cast, castM_mul]
  unfold castM
  rw [mulZ_ap i j]

theorem mulQ_na (i j : Fin 12) :
    mulQ nOrb aOrb i j = 2 * aOrb i j + 2 * nOrb i j + 5 * pOrb i j := by
  rw [aOrb_cast, nOrb_cast, pOrb_cast, castM_mul]
  unfold castM
  rw [mulZ_na i j]; push_cast; ring

theorem mulQ_nn (i j : Fin 12) :
    mulQ nOrb nOrb i j = 5 * dOrb i j + 2 * aOrb i j + 2 * nOrb i j := by
  rw [dOrb_cast, aOrb_cast, nOrb_cast, castM_mul]
  unfold castM
  rw [mulZ_nn i j]; push_cast; ring

theorem mulQ_np (i j : Fin 12) : mulQ nOrb pOrb i j = aOrb i j := by
  rw [aOrb_cast, nOrb_cast, pOrb_cast, castM_mul]
  unfold castM
  rw [mulZ_np i j]

theorem mulQ_pa (i j : Fin 12) : mulQ pOrb aOrb i j = nOrb i j := by
  rw [aOrb_cast, nOrb_cast, pOrb_cast, castM_mul]
  unfold castM
  rw [mulZ_pa i j]

theorem mulQ_pn (i j : Fin 12) : mulQ pOrb nOrb i j = aOrb i j := by
  rw [aOrb_cast, nOrb_cast, pOrb_cast, castM_mul]
  unfold castM
  rw [mulZ_pn i j]

theorem mulQ_pp (i j : Fin 12) : mulQ pOrb pOrb i j = dOrb i j := by
  rw [dOrb_cast, pOrb_cast, castM_mul]
  unfold castM
  rw [mulZ_pp i j]

/-! ## The span and its multiplicative closure -/

/-- Membership in the four-dimensional orbital span. -/
def InSpan (M : Fin 12 → Fin 12 → ℚ) : Prop :=
  ∃ a b c d : ℚ, ∀ i j : Fin 12,
    M i j = a * dOrb i j + b * aOrb i j + c * nOrb i j + d * pOrb i j

/-- Products of span members expand along the orbital product table; the
coefficient formula is symmetric in the two factors. -/
theorem mulQ_span_expand (a b c d a' b' c' d' : ℚ)
    (M N : Fin 12 → Fin 12 → ℚ)
    (hM : ∀ i j : Fin 12,
      M i j = a * dOrb i j + b * aOrb i j + c * nOrb i j + d * pOrb i j)
    (hN : ∀ i j : Fin 12,
      N i j = a' * dOrb i j + b' * aOrb i j + c' * nOrb i j + d' * pOrb i j)
    (i j : Fin 12) :
    mulQ M N i j =
      (a * a' + 5 * (b * b') + 5 * (c * c') + d * d') * dOrb i j +
      (a * b' + b * a' + 2 * (b * b') + 2 * (b * c') + 2 * (c * b') +
        2 * (c * c') + c * d' + d * c') * aOrb i j +
      (a * c' + c * a' + 2 * (b * b') + 2 * (b * c') + 2 * (c * b') +
        2 * (c * c') + b * d' + d * b') * nOrb i j +
      (a * d' + d * a' + 5 * (b * c') + 5 * (c * b')) * pOrb i j := by
  have step1 : mulQ M N i j =
      a * a' * mulQ dOrb dOrb i j + a * b' * mulQ dOrb aOrb i j +
      a * c' * mulQ dOrb nOrb i j + a * d' * mulQ dOrb pOrb i j +
      b * a' * mulQ aOrb dOrb i j + b * b' * mulQ aOrb aOrb i j +
      b * c' * mulQ aOrb nOrb i j + b * d' * mulQ aOrb pOrb i j +
      c * a' * mulQ nOrb dOrb i j + c * b' * mulQ nOrb aOrb i j +
      c * c' * mulQ nOrb nOrb i j + c * d' * mulQ nOrb pOrb i j +
      d * a' * mulQ pOrb dOrb i j + d * b' * mulQ pOrb aOrb i j +
      d * c' * mulQ pOrb nOrb i j + d * d' * mulQ pOrb pOrb i j := by
    unfold mulQ
    rw [Finset.sum_congr rfl fun k _ =>
      (by rw [hM i k, hN k j]; ring :
        M i k * N k j =
          a * a' * (dOrb i k * dOrb k j) + a * b' * (dOrb i k * aOrb k j) +
          a * c' * (dOrb i k * nOrb k j) + a * d' * (dOrb i k * pOrb k j) +
          b * a' * (aOrb i k * dOrb k j) + b * b' * (aOrb i k * aOrb k j) +
          b * c' * (aOrb i k * nOrb k j) + b * d' * (aOrb i k * pOrb k j) +
          c * a' * (nOrb i k * dOrb k j) + c * b' * (nOrb i k * aOrb k j) +
          c * c' * (nOrb i k * nOrb k j) + c * d' * (nOrb i k * pOrb k j) +
          d * a' * (pOrb i k * dOrb k j) + d * b' * (pOrb i k * aOrb k j) +
          d * c' * (pOrb i k * nOrb k j) + d * d' * (pOrb i k * pOrb k j))]
    simp only [Finset.sum_add_distrib, ← Finset.mul_sum]
  rw [step1, mulQ_dd i j, mulQ_da i j, mulQ_dn i j, mulQ_dp i j,
    mulQ_ad i j, mulQ_aa i j, mulQ_an i j, mulQ_ap i j,
    mulQ_nd i j, mulQ_na i j, mulQ_nn i j, mulQ_np i j,
    mulQ_pd i j, mulQ_pa i j, mulQ_pn i j, mulQ_pp i j]
  ring

/-- The orbital span is closed under matrix multiplication. -/
theorem inSpan_mul {M N : Fin 12 → Fin 12 → ℚ}
    (hM : InSpan M) (hN : InSpan N) : InSpan (mulQ M N) := by
  obtain ⟨a, b, c, d, hM⟩ := hM
  obtain ⟨a', b', c', d', hN⟩ := hN
  exact ⟨_, _, _, _, fun i j => mulQ_span_expand a b c d a' b' c' d' M N hM hN i j⟩

/-- Multiplication on the orbital span is commutative: the expansion
coefficients are symmetric in the two factors. -/
theorem inSpan_mul_comm {M N : Fin 12 → Fin 12 → ℚ}
    (hM : InSpan M) (hN : InSpan N) (i j : Fin 12) :
    mulQ M N i j = mulQ N M i j := by
  obtain ⟨a, b, c, d, hM⟩ := hM
  obtain ⟨a', b', c', d', hN⟩ := hN
  rw [mulQ_span_expand a b c d a' b' c' d' M N hM hN i j,
    mulQ_span_expand a' b' c' d' a b c d N M hN hM i j]
  ring

theorem inSpan_dOrb : InSpan dOrb := ⟨1, 0, 0, 0, fun i j => by ring⟩

theorem inSpan_aOrb : InSpan aOrb := ⟨0, 1, 0, 0, fun i j => by ring⟩

theorem inSpan_pOrb : InSpan pOrb := ⟨0, 0, 0, 1, fun i j => by ring⟩

/-! ## Adjacency powers, the antipode readback, and the quartic -/

/-- The adjacency squared, as an actual matrix product. -/
def sqPow : Fin 12 → Fin 12 → ℚ := fun i j => mulQ aOrb aOrb i j

/-- The adjacency cubed. -/
def cubePow : Fin 12 → Fin 12 → ℚ := fun i j => mulQ aOrb sqPow i j

/-- The adjacency to the fourth. -/
def fourthPow : Fin 12 → Fin 12 → ℚ := fun i j => mulQ aOrb cubePow i j

theorem sq_orbital (i j : Fin 12) :
    sqPow i j = 5 * dOrb i j + 2 * aOrb i j + 2 * nOrb i j := mulQ_aa i j

theorem cube_orbital (i j : Fin 12) :
    cubePow i j =
      10 * dOrb i j + 13 * aOrb i j + 8 * nOrb i j + 10 * pOrb i j := by
  have h := mulQ_span_expand 0 1 0 0 5 2 2 0 aOrb sqPow
    (fun i j => by ring) (fun i j => by rw [sq_orbital i j]; ring) i j
  calc cubePow i j = mulQ aOrb sqPow i j := rfl
    _ = _ := by rw [h]; ring

theorem fourth_orbital (i j : Fin 12) :
    fourthPow i j =
      65 * dOrb i j + 52 * aOrb i j + 52 * nOrb i j + 40 * pOrb i j := by
  have h := mulQ_span_expand 0 1 0 0 10 13 8 10 aOrb cubePow
    (fun i j => by ring) (fun i j => cube_orbital i j) i j
  calc fourthPow i j = mulQ aOrb cubePow i j := rfl
    _ = _ := by rw [h]; ring

/-- The polynomial antipode readback, with the arithmetic carried out: the
antipode permutation is the exact cubic
`P = (A³ - 4·A² - 5·A + 10·I) / 10` in the registered adjacency. -/
theorem antipode_readback (i j : Fin 12) :
    pOrb i j =
      (cubePow i j - 4 * sqPow i j - 5 * aOrb i j + 10 * dOrb i j) / 10 := by
  rw [cube_orbital i j, sq_orbital i j]
  ring

/-- The registered adjacency satisfies the exact quartic
`A⁴ = 4·A³ + 10·A² - 20·A - 25·I`: the adjacency recurrence closes at
length four. -/
theorem adjacency_quartic (i j : Fin 12) :
    fourthPow i j =
      4 * cubePow i j + 10 * sqPow i j - 20 * aOrb i j - 25 * dOrb i j := by
  rw [fourth_orbital i j, cube_orbital i j, sq_orbital i j]
  ring

/-! ## The response-word algebra -/

/-- Generators of the response-word algebra: one adjacency recurrence step,
one polynomial antipode readback. -/
inductive Gen
  | step
  | readback

/-- The matrix of one generator. -/
def genMat : Gen → (Fin 12 → Fin 12 → ℚ)
  | .step => aOrb
  | .readback => pOrb

/-- The matrix of a response word; the empty word is the identity. -/
def wordMat : List Gen → (Fin 12 → Fin 12 → ℚ)
  | [] => dOrb
  | g :: w => fun i j => mulQ (genMat g) (wordMat w) i j

theorem inSpan_genMat : ∀ g : Gen, InSpan (genMat g)
  | .step => inSpan_aOrb
  | .readback => inSpan_pOrb

/-- Every response word lies in the four-dimensional orbital span. -/
theorem inSpan_wordMat : ∀ w : List Gen, InSpan (wordMat w)
  | [] => inSpan_dOrb
  | g :: w => inSpan_mul (inSpan_genMat g) (inSpan_wordMat w)

/-- The response-word algebra is commutative: any two words multiply the
same way in either order. -/
theorem word_mul_comm (v w : List Gen) (i j : Fin 12) :
    mulQ (wordMat v) (wordMat w) i j = mulQ (wordMat w) (wordMat v) i j :=
  inSpan_mul_comm (inSpan_wordMat v) (inSpan_wordMat w) i j

/-- No response word has a nonzero commutator with any other: the registered
response cannot supply a nonzero commutator (issue #705, boundary item two,
commutator half). -/
theorem no_nonzero_word_commutator (v w : List Gen) (i j : Fin 12) :
    mulQ (wordMat v) (wordMat w) i j - mulQ (wordMat w) (wordMat v) i j = 0 := by
  rw [word_mul_comm v w i j]
  ring

/-! ## Dimension exactly four, both directions -/

/-- The orbital span and the power span `⟨I, A, A², A³⟩` coincide: the
change of basis is invertible because the readback denominator `10` and the
distance-two coefficient `2` are nonzero. -/
theorem inSpan_iff_power_combo (M : Fin 12 → Fin 12 → ℚ) :
    InSpan M ↔
      ∃ x y z w : ℚ, ∀ i j : Fin 12,
        M i j = x * dOrb i j + y * aOrb i j + z * sqPow i j + w * cubePow i j := by
  constructor
  · rintro ⟨a, b, c, d, h⟩
    refine ⟨a - 5 * ((c - 8 * (d / 10)) / 2) - 10 * (d / 10),
      b - 2 * ((c - 8 * (d / 10)) / 2) - 13 * (d / 10),
      (c - 8 * (d / 10)) / 2, d / 10, fun i j => ?_⟩
    rw [h i j, sq_orbital i j, cube_orbital i j]
    ring
  · rintro ⟨x, y, z, w, h⟩
    exact ⟨x + 5 * z + 10 * w, y + 2 * z + 13 * w, 2 * z + 8 * w, 10 * w,
      fun i j => by rw [h i j, sq_orbital i j, cube_orbital i j]; ring⟩

/-- Spanning direction: every response word is a rational combination of
`I, A, A², A³`. -/
theorem word_in_power_span (w : List Gen) :
    ∃ x y z t : ℚ, ∀ i j : Fin 12,
      wordMat w i j =
        x * dOrb i j + y * aOrb i j + z * sqPow i j + t * cubePow i j :=
  (inSpan_iff_power_combo (wordMat w)).mp (inSpan_wordMat w)

/-- Independence direction: `I, A, A², A³` are `ℚ`-linearly independent, so
the response-word algebra has dimension exactly four. -/
theorem powers_independent (x y z w : ℚ)
    (h : ∀ i j : Fin 12,
      x * dOrb i j + y * aOrb i j + z * sqPow i j + w * cubePow i j = 0) :
    x = 0 ∧ y = 0 ∧ z = 0 ∧ w = 0 := by
  have horb := orbitals_independent (x + 5 * z + 10 * w) (y + 2 * z + 13 * w)
    (2 * z + 8 * w) (10 * w) (fun i j => by
      have hh := h i j
      rw [sq_orbital i j, cube_orbital i j] at hh
      linear_combination hh)
  obtain ⟨h1, h2, h3, h4⟩ := horb
  have hw : w = 0 := by linarith
  have hz : z = 0 := by linarith
  have hy : y = 0 := by linarith
  have hx : x = 0 := by linarith
  exact ⟨hx, hy, hz, hw⟩

/-! ## The word span is exactly the commutant -/

/-- Every equivariant matrix is a rational combination of `I, A, A², A³`,
and conversely every such combination is equivariant: the response-word
algebra is exactly the four-dimensional commutant of
`Screen/A5Commutant.lean`.  Nothing outside the span survives equivariance. -/
theorem equivariant_iff_word_span (M : Fin 12 → Fin 12 → ℚ) :
    Equivariant M ↔
      ∃ x y z w : ℚ, ∀ i j : Fin 12,
        M i j = x * dOrb i j + y * aOrb i j + z * sqPow i j + w * cubePow i j := by
  constructor
  · intro hM
    exact (inSpan_iff_power_combo M).mp
      ⟨M 0 0, M 0 1, M 0 5, M 0 11, commutant_decomposition M hM⟩
  · intro hM
    obtain ⟨a, b, c, d, h⟩ := (inSpan_iff_power_combo M).mpr hM
    intro g hg i j
    rw [h (g i) (g j), h i j,
      equivariant_dOrb g hg i j, equivariant_aOrb g hg i j,
      equivariant_nOrb g hg i j, equivariant_pOrb g hg i j]

/-- Every response word is equivariant. -/
theorem wordMat_equivariant (w : List Gen) : Equivariant (wordMat w) :=
  (equivariant_iff_word_span (wordMat w)).mpr (word_in_power_span w)

/-! ## No twelve independent generators -/

/-- No twelve members of the word span are `ℚ`-linearly independent: any
twelve carry a nontrivial vanishing rational combination, because the span
has dimension four.  In particular the registered response cannot supply
twelve independent current generators. -/
theorem no_twelve_independent_span_members
    (f : Fin 12 → (Fin 12 → Fin 12 → ℚ)) (hf : ∀ t, InSpan (f t)) :
    ∃ c : Fin 12 → ℚ, (∃ t, c t ≠ 0) ∧
      ∀ i j : Fin 12, ∑ t, c t * f t i j = 0 := by
  choose a b c d hrep using hf
  set v : Fin 12 → (Fin 4 → ℚ) := fun t => ![a t, b t, c t, d t] with hv
  have hnli : ¬ LinearIndependent ℚ v := by
    intro hli
    have hcard := hli.fintype_card_le_finrank
    rw [Module.finrank_fin_fun] at hcard
    simp at hcard
  obtain ⟨g, hsum, t0, hg⟩ := Fintype.not_linearIndependent_iff.mp hnli
  have hcoord : ∀ m : Fin 4, (∑ t, g t * v t m) = 0 := by
    intro m
    have hm := congrFun hsum m
    simpa [Finset.sum_apply, Pi.smul_apply, smul_eq_mul] using hm
  have ha0 : (∑ t, g t * a t) = 0 := by simpa [hv] using hcoord 0
  have hb0 : (∑ t, g t * b t) = 0 := by simpa [hv] using hcoord 1
  have hc0 : (∑ t, g t * c t) = 0 := by simpa [hv] using hcoord 2
  have hd0 : (∑ t, g t * d t) = 0 := by simpa [hv] using hcoord 3
  refine ⟨g, ⟨t0, hg⟩, fun i j => ?_⟩
  calc ∑ t, g t * f t i j
      = ∑ t, (g t * a t * dOrb i j + g t * b t * aOrb i j +
          g t * c t * nOrb i j + g t * d t * pOrb i j) :=
        Finset.sum_congr rfl fun t _ => by rw [hrep t i j]; ring
    _ = (∑ t, g t * a t) * dOrb i j + (∑ t, g t * b t) * aOrb i j +
        (∑ t, g t * c t) * nOrb i j + (∑ t, g t * d t) * pOrb i j := by
        simp only [Finset.sum_add_distrib, ← Finset.sum_mul]
    _ = 0 := by rw [ha0, hb0, hc0, hd0]; ring

/-- Twelve response words are never `ℚ`-linearly independent. -/
theorem no_twelve_independent_words (f : Fin 12 → List Gen) :
    ∃ c : Fin 12 → ℚ, (∃ t, c t ≠ 0) ∧
      ∀ i j : Fin 12, ∑ t, c t * wordMat (f t) i j = 0 :=
  no_twelve_independent_span_members (fun t => wordMat (f t))
    (fun t => inSpan_wordMat (f t))

/-! ## Negative controls -/

/-- Broken incidence: the registered adjacency with the single edge `{0,1}`
deleted. -/
def adjB (i j : Fin 12) : Bool :=
  adj i j && !((i == 0 && j == 1) || (i == 1 && j == 0))

/-- The broken adjacency matrix over `ℤ`. -/
def aB (i j : Fin 12) : ℤ := if adjB i j then 1 else 0

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- CONTROL.  On the broken incidence the cubic antipode readback fails at
entry `(0,1)`: the readback identity tests the registered incidence rather
than holding for arbitrary graphs. -/
theorem broken_readback_fails :
    ¬ (mulZ aB (mulZ aB aB) 0 1 - 4 * mulZ aB aB 0 1 - 5 * aB 0 1 +
        10 * dZ 0 1 = 10 * pZ 0 1) := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- CONTROL.  On the broken incidence the quartic fails at entry `(0,0)`:
deleting one edge already destroys the length-four closure of the adjacency
recurrence, so `adjacency_quartic` is not an artifact of the proof method. -/
theorem broken_quartic_fails :
    ¬ (mulZ aB (mulZ aB (mulZ aB aB)) 0 0 =
        4 * mulZ aB (mulZ aB aB) 0 0 + 10 * mulZ aB aB 0 0 -
        20 * aB 0 0 - 25 * dZ 0 0) := by decide

/-- The elementary comparison matrix `E₀₁` over `ℤ`. -/
def eZ (i j : Fin 12) : ℤ := if i = 0 ∧ j = 1 then 1 else 0

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- CONTROL.  The elementary matrix `E₀₁` has a nonzero commutator with the
registered adjacency at entry `(2,1)`: the vanishing of all word
commutators is a property of the response-word algebra, not of twelve-port
matrices in general. -/
theorem elementary_commutator_nonzero :
    ¬ (mulZ aZ eZ 2 1 = mulZ eZ aZ 2 1) := by decide

/-- The same commutator control at the `ℚ` level. -/
theorem comparison_commutator_nonzero :
    ¬ (∀ i j : Fin 12,
        mulQ aOrb (castM eZ) i j = mulQ (castM eZ) aOrb i j) := by
  intro h
  have h21 := h 2 1
  rw [aOrb_cast, castM_mul, castM_mul] at h21
  unfold castM at h21
  exact elementary_commutator_nonzero (by exact_mod_cast h21)

/-- Entry values of the four orbitals at the probe entries. -/
theorem orbital_entries :
    dOrb 0 1 = 0 ∧ aOrb 0 1 = 1 ∧ nOrb 0 1 = 0 ∧ pOrb 0 1 = 0 ∧
    dOrb 0 2 = 0 ∧ aOrb 0 2 = 1 ∧ nOrb 0 2 = 0 ∧ pOrb 0 2 = 0 ∧
    dOrb 0 11 = 0 ∧ aOrb 0 11 = 0 ∧ nOrb 0 11 = 0 ∧ pOrb 0 11 = 1 := by
  have e1 : adj (0 : Fin 12) 1 = true := by decide
  have e2 : adj (0 : Fin 12) 2 = true := by decide
  have e3 : adj (0 : Fin 12) 11 = false := by decide
  have e4 : antipode (0 : Fin 12) = 11 := by decide
  have f1 : (0 : Fin 12) ≠ 1 := by decide
  have f2 : (0 : Fin 12) ≠ 2 := by decide
  have f3 : (0 : Fin 12) ≠ 11 := by decide
  have g1 : (1 : Fin 12) ≠ 11 := by decide
  have g2 : (2 : Fin 12) ≠ 11 := by decide
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
    simp [OPH.A5Commutant.dOrb, OPH.A5Commutant.aOrb, OPH.A5Commutant.nOrb,
      OPH.A5Commutant.pOrb, e1, e2, e3, e4, f1, f2, f3, g1, g2]

/-- CONTROL.  The comparison matrix `E₀₁` is outside the word span: it takes
different values on the equivalent adjacent entries `(0,1)` and `(0,2)`,
which no span member does. -/
theorem comparison_not_in_span : ¬ InSpan (castM eZ) := by
  rintro ⟨a, b, c, d, h⟩
  obtain ⟨d1, a1, n1, p1, d2, a2, n2, p2, _, _, _, _⟩ := orbital_entries
  have h1 := h 0 1
  have h2 := h 0 2
  rw [d1, a1, n1, p1] at h1
  rw [d2, a2, n2, p2] at h2
  have v1 : castM eZ 0 1 = 1 := by unfold castM eZ; norm_num
  have v2 : castM eZ 0 2 = 0 := by
    unfold castM eZ
    rw [if_neg (by decide : ¬ ((0 : Fin 12) = 0 ∧ (2 : Fin 12) = 1))]
    norm_num
  rw [v1] at h1
  rw [v2] at h2
  linarith

/-- CONTROL.  `A³` is outside the span of `I, A, A²`: the response-word
algebra is not three-dimensional, so `powers_independent` is sharp. -/
theorem cube_not_in_three_span :
    ¬ ∃ x y z : ℚ, ∀ i j : Fin 12,
      cubePow i j = x * dOrb i j + y * aOrb i j + z * sqPow i j := by
  rintro ⟨x, y, z, h⟩
  obtain ⟨_, _, _, _, _, _, _, _, d3, a3, n3, p3⟩ := orbital_entries
  have h011 := h 0 11
  rw [cube_orbital 0 11, sq_orbital 0 11, d3, a3, n3, p3] at h011
  norm_num at h011

/-- CONTROL.  The antipode readback admits no polynomial of degree at most
two in the adjacency: the cubic in `antipode_readback` is sharp. -/
theorem antipode_not_quadratic :
    ¬ ∃ x y z : ℚ, ∀ i j : Fin 12,
      pOrb i j = x * dOrb i j + y * aOrb i j + z * sqPow i j := by
  rintro ⟨x, y, z, h⟩
  obtain ⟨_, _, _, _, _, _, _, _, d3, a3, n3, p3⟩ := orbital_entries
  have h011 := h 0 11
  rw [sq_orbital 0 11, d3, a3, n3, p3] at h011
  norm_num at h011

end OPH.A5ResponseWordAlgebra

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.A5ResponseWordAlgebra.antipode_readback
#print axioms OPH.A5ResponseWordAlgebra.adjacency_quartic
#print axioms OPH.A5ResponseWordAlgebra.word_in_power_span
#print axioms OPH.A5ResponseWordAlgebra.powers_independent
#print axioms OPH.A5ResponseWordAlgebra.word_mul_comm
#print axioms OPH.A5ResponseWordAlgebra.no_nonzero_word_commutator
#print axioms OPH.A5ResponseWordAlgebra.equivariant_iff_word_span
#print axioms OPH.A5ResponseWordAlgebra.wordMat_equivariant
#print axioms OPH.A5ResponseWordAlgebra.no_twelve_independent_words
#print axioms OPH.A5ResponseWordAlgebra.broken_readback_fails
#print axioms OPH.A5ResponseWordAlgebra.broken_quartic_fails
#print axioms OPH.A5ResponseWordAlgebra.comparison_commutator_nonzero
#print axioms OPH.A5ResponseWordAlgebra.comparison_not_in_span
#print axioms OPH.A5ResponseWordAlgebra.cube_not_in_three_span
#print axioms OPH.A5ResponseWordAlgebra.antipode_not_quadratic
