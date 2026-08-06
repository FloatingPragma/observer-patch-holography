import QFT.ConvexAffinityBridge

/-!
# Gauge irreducibility and the Born weight of the invariant state

This module is the gauge route to the Born rule on the realized split
fibre: the third of three independent derivation routes, next to the
Busch/Gleason effect-additivity route and the realized-frequency route.
The route derives the Born weight for the gauge-invariant state alone,
from symmetry and convex affinity, with no additivity postulate over
effect decompositions.

## The exact representation

`gaugePair` records the exact two-dimensional irreducible representation
of the symmetric group `S3 = Equiv.Perm (Fin 3)` through doubled integer
components: `gaugePair g = (A, B)` encodes the orthogonal matrix
`gaugeRep g = (A + sqrt 3 * B) / 2`, so every entry of `2 * gaugeRep g`
lies in the ring of integers of `ℚ(sqrt 3)`.  The component tables are
kernel-decided over `ℤ` (`gaugePair_mulTable`, `gaugePair_orth`), and a
single `Real.sq_sqrt` receipt (`sqrt3_sq`, `sqrt3C_mul_self`) lifts both
to the complex representation: the homomorphism law `gaugeRep_mul` on
the full six-element multiplication table and the unitarity receipts
`gaugeRep_conjTranspose_mul_self`, `gaugeRep_mul_conjTranspose_self`.
The six matrices are stated entrywise in
`gaugeRep_one_eq` through `gaugeRep_trr_eq`.

The matrices agree entry for entry with the `irrep.elements` block of
the committed payload `docs/BORN_CONTEXT_WEB_PAYLOAD.json` of the
`oph-physics-sim` repository (schema
`oph.sim.born_context_web_payload.v1`): scalar encoding `a + b*sqrt(3)`
with rational pairs, composition convention `(g*h)[i] = g[h[i]]`, and
generator assignment `rho((1,2,0)) = rotation by 120 degrees`,
`rho((0,2,1)) = diag(1, -1)`.  Here `rGen` is the permutation
`(1,2,0)` and `tGen` is `(0,2,1)`.

## Irreducibility, the invariant state, and the derivation

* `commutant_pair_scalar` / `commutant_scalar`: finite Schur form.
  Every complex `2 x 2` matrix commuting with all six images is scalar;
  the proof is direct linear algebra from the order-3 rotation and one
  reflection alone.
* `invariant_state_eq_tracialState`: every density matrix invariant
  under conjugation by all six images equals `tracialState`, the state
  `(1/2) • 1`, by Schur plus the trace normalization.
* `born_weight_forced_at_invariant_state`: the derivation statement.
  Any outcome assignment `v` on states and effects that (a) mixes
  affinely in the effect slot at the invariant state, in the exact
  binary-mixing grammar of the convex-affinity bridge, (b) is covariant
  under the realized gauge action, assigning to a rotated effect in the
  rotated state the original value, and (c) is normalized at the sure
  and impossible effects, takes on every projection event `E` the value
  `v tracialState E = Re (Tr ((1/2) • 1 * E))`
  (`bornWeight_tracialState_def` displays the right side in that shape;
  the trace is real for events, so the real part is a display choice).
  The proof averages the effect over the six realized gauge images
  through five binary mixes, identifies the average as a scalar by
  Schur, and evaluates the scalar by the trace.

## The convex-affinity bridge, consumed

The mixing grammar of hypothesis (a) is the mixing law of
`ConvexPreparationFamily` from `QFT/ConvexAffinityBridge.lean`,
transplanted from the state slot to the effect slot.  The bridge itself
is consumed by instantiation in `bornAssignment_family_mix`: for every
convex preparation family, the trace-pairing assignment satisfies
state-slot affinity as the bridge theorem
`ConvexPreparationFamily.outcome_re_affine`, with `segmentFamily` and
`towerSegmentFamily` as the committed producers.  The negative control
`traceSquareAssignment` proves that state-slot affinity, covariance,
and normalization together cannot force the conclusion: the effect-slot
mixing hypothesis is load-bearing, and the theorem states the sharpest
clean form along this route.

## The source tie

The committed run `b12_prereg_16k_20260806` (seed `20260806`, commit
`b39b78faf894894ebe573571e0902ccfaaeac32a`, gauge state sha256
`5f2cd276b84b1917c4d321a532db6efd12442931fbc80fd764cbcf1918dadca2`)
realizes the six gauge elements on the edges incident to the designated
split fibre of observer 64, record class 18, nodes 64 and 43: of the 21
incident edge slots the payload counts 8, 4, 3, 2, 3, 1 across the six
elements, every element with positive count.  The literals are
transcribed in `designatedFibreEdgeCount`; the receipt-shaped Prop
`RealizesAllGaugeElements` names exactly what the run supplies, and
`born_weight_forced_with_realized_covariance` consumes it: covariance
is assumed only at gauge elements with positive realized edge count,
and the census receipt `designatedFibre_realizes_all` discharges the
rest.  The payload also names its own capability gap: no realized
gauge-conjugated outcome frequencies exist on the committed fibre, and
Born predictions reported as rotated-context frequencies are a
prohibited fill.  This module consumes only the element realization;
frequency receipts belong to the realized-frequency route and are out
of scope.  The interface is stated here without importing any payload
module, so no dependency on `QFT/SourceContextWeb.lean` exists.

## Witnesses and negative controls

`tracial_born_vertexZero` and `tracial_born_vertexPlus` compute the
Born weights of the tracial state on the two noncommuting projection
events of the bridge witness layer (`vertexZero_vertexPlus_noncomm`):
both are the exact half, and `admissible_values_on_noncommuting_pair`
shows every admissible assignment takes that value on both.  The
negative controls: `cornerAssignment` satisfies effect-slot mixing and
normalization, breaks covariance on a realized rotation
(`cornerAssignment_breaks_covariance`), and evades the conclusion
(`cornerAssignment_evades_born`); `traceSquareAssignment` satisfies
state-slot affinity, covariance, and normalization, violates
effect-slot mixing (`traceSquareAssignment_not_effect_mix`), and evades
the conclusion (`traceSquareAssignment_evades_born`).

## Claim boundary

The derivation fixes the outcome assignment at the gauge-invariant
state only; states that break the symmetry are untouched by this route
and need the effect-additivity route or the realized-frequency route.
The three routes are independent lanes converging on the one trace
rule, and that convergence is the point of the program: each route
carries its own hypotheses, and the shared conclusion is the Born
weight.  The gauge elements and the fibre census are realized data of
one committed run; the outcome assignments quantified over here are
mathematical objects, and no physical measurement attachment is
claimed.
-/

namespace OPH.QFT.GaugeBorn

open EventAlgebra
open Matrix
open scoped ComplexOrder

/-- The gauge group of the split fibre: permutations of three letters. -/
abbrev S3 := Equiv.Perm (Fin 3)

/-- The order-3 generator, the payload permutation `(1,2,0)`. -/
def rGen : S3 := Equiv.swap 0 1 * Equiv.swap 1 2

/-- The order-2 generator, the payload permutation `(0,2,1)`. -/
def tGen : S3 := Equiv.swap 1 2

theorem rGen_apply : rGen 0 = 1 ∧ rGen 1 = 2 ∧ rGen 2 = 0 := by decide

theorem tGen_apply : tGen 0 = 0 ∧ tGen 1 = 2 ∧ tGen 2 = 1 := by decide

/-! ## The exact component tables over `ℤ` -/

/-- Doubled integer components of the standard two-dimensional
irreducible representation: `gaugePair g = (A, B)` encodes the matrix
`(A + sqrt 3 * B) / 2`.  The six values mirror the `irrep.elements`
block of the committed payload, entry for entry after doubling. -/
def gaugePair (g : S3) : Matrix (Fin 2) (Fin 2) ℤ × Matrix (Fin 2) (Fin 2) ℤ :=
  if g = 1 then (!![2, 0; 0, 2], 0)
  else if g = rGen then (!![-1, 0; 0, -1], !![0, -1; 1, 0])
  else if g = rGen * rGen then (!![-1, 0; 0, -1], !![0, 1; -1, 0])
  else if g = tGen then (!![2, 0; 0, -2], 0)
  else if g = tGen * rGen then (!![-1, 0; 0, 1], !![0, -1; -1, 0])
  else (!![-1, 0; 0, 1], !![0, 1; 1, 0])

set_option maxRecDepth 8000 in
/-- Kernel-decided homomorphism law on the full six-element
multiplication table, at the level of exact integer components: the
rational parts and the `sqrt 3` parts of `4 * (gaugeRep g * gaugeRep h)`
against `4 * gaugeRep (g * h)`. -/
theorem gaugePair_mulTable : ∀ g h : S3,
    (gaugePair g).1 * (gaugePair h).1
        + ((gaugePair g).2 * (gaugePair h).2 + (gaugePair g).2 * (gaugePair h).2
          + (gaugePair g).2 * (gaugePair h).2)
      = (gaugePair (g * h)).1 + (gaugePair (g * h)).1
    ∧ (gaugePair g).1 * (gaugePair h).2 + (gaugePair g).2 * (gaugePair h).1
      = (gaugePair (g * h)).2 + (gaugePair (g * h)).2 := by
  decide

set_option maxRecDepth 8000 in
/-- Kernel-decided unitarity components: the exact integer form of
`(2 * gaugeRep g)ᵀ * (2 * gaugeRep g) = 4 * 1`. -/
theorem gaugePair_orth : ∀ g : S3,
    (gaugePair g).1ᵀ * (gaugePair g).1
        + ((gaugePair g).2ᵀ * (gaugePair g).2 + (gaugePair g).2ᵀ * (gaugePair g).2
          + (gaugePair g).2ᵀ * (gaugePair g).2)
      = !![4, 0; 0, 4]
    ∧ (gaugePair g).1ᵀ * (gaugePair g).2 + (gaugePair g).2ᵀ * (gaugePair g).1 = 0 := by
  decide

theorem gaugePair_one : gaugePair 1 = (!![2, 0; 0, 2], 0) := by decide

theorem gaugePair_r : gaugePair rGen = (!![-1, 0; 0, -1], !![0, -1; 1, 0]) := by decide

theorem gaugePair_rr :
    gaugePair (rGen * rGen) = (!![-1, 0; 0, -1], !![0, 1; -1, 0]) := by decide

theorem gaugePair_t : gaugePair tGen = (!![2, 0; 0, -2], 0) := by decide

theorem gaugePair_tr :
    gaugePair (tGen * rGen) = (!![-1, 0; 0, 1], !![0, -1; -1, 0]) := by decide

theorem gaugePair_trr :
    gaugePair (tGen * (rGen * rGen)) = (!![-1, 0; 0, 1], !![0, 1; 1, 0]) := by decide

/-! ## The `sqrt 3` receipts and the lift to `ℂ` -/

/-- The `Real.sq_sqrt` receipt for the adjoined scalar. -/
theorem sqrt3_sq : Real.sqrt 3 ^ 2 = 3 := Real.sq_sqrt (by norm_num)

/-- The adjoined scalar, encoded as `Real.sqrt 3` inside `ℂ`. -/
noncomputable def sqrt3C : ℂ := ((Real.sqrt 3 : ℝ) : ℂ)

theorem sqrt3C_mul_self : sqrt3C * sqrt3C = 3 := by
  have h : Real.sqrt 3 * Real.sqrt 3 = 3 := Real.mul_self_sqrt (by norm_num)
  rw [sqrt3C, ← Complex.ofReal_mul, h]
  norm_num

theorem sqrt3C_ne_zero : sqrt3C ≠ 0 := by
  rw [sqrt3C, Complex.ofReal_ne_zero]
  positivity

theorem star_sqrt3C : star sqrt3C = sqrt3C := by
  rw [sqrt3C, Complex.star_def, Complex.conj_ofReal]

/-- Entrywise integer cast into complex matrices, as a ring
homomorphism. -/
noncomputable def intLift : Matrix (Fin 2) (Fin 2) ℤ →+* Matrix (Fin 2) (Fin 2) ℂ :=
  (Int.castRingHom ℂ).mapMatrix

theorem intLift_conjTranspose (M : Matrix (Fin 2) (Fin 2) ℤ) :
    (intLift M)ᴴ = intLift Mᵀ := by
  ext i j
  simp [intLift, Matrix.conjTranspose_apply, Matrix.map_apply]

/-- The complex representation: `gaugeRep g = (A + sqrt 3 * B) / 2` for
the decided integer components `(A, B) = gaugePair g`. -/
noncomputable def gaugeRep (g : S3) : Matrix (Fin 2) (Fin 2) ℂ :=
  (2⁻¹ : ℂ) • (intLift (gaugePair g).1 + sqrt3C • intLift (gaugePair g).2)

/-- Product expansion of two half-scaled `ℤ[sqrt 3]` matrix pairs; the
one place the `sqrt 3` receipt enters the multiplication law. -/
private theorem intLift_pair_mul (A B C D : Matrix (Fin 2) (Fin 2) ℤ) :
    ((2⁻¹ : ℂ) • (intLift A + sqrt3C • intLift B))
        * ((2⁻¹ : ℂ) • (intLift C + sqrt3C • intLift D))
      = (4⁻¹ : ℂ) • (intLift (A * C + (B * D + B * D + B * D))
          + sqrt3C • intLift (A * D + B * C)) := by
  simp only [map_add, map_mul]
  rw [Matrix.smul_mul, Matrix.mul_smul, smul_smul]
  rw [add_mul, mul_add, mul_add]
  simp only [smul_mul_assoc, mul_smul_comm, smul_smul, sqrt3C_mul_self]
  norm_num
  module

/-- The homomorphism law on the six-element multiplication table,
lifted from the kernel-decided integer components. -/
theorem gaugeRep_mul (g h : S3) : gaugeRep (g * h) = gaugeRep g * gaugeRep h := by
  obtain ⟨h1, h2⟩ := gaugePair_mulTable g h
  rw [gaugeRep, gaugeRep, gaugeRep, intLift_pair_mul, h1, h2]
  simp only [map_add]
  module

theorem gaugeRep_one : gaugeRep 1 = 1 := by
  rw [gaugeRep, gaugePair_one]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [intLift, Matrix.map_apply]

/-- Conjugate transpose acts on the integer components as the
transpose: all entries are real. -/
theorem gaugeRep_conjTranspose (g : S3) :
    (gaugeRep g)ᴴ
      = (2⁻¹ : ℂ) • (intLift (gaugePair g).1ᵀ + sqrt3C • intLift (gaugePair g).2ᵀ) := by
  rw [gaugeRep, Matrix.conjTranspose_smul, Matrix.conjTranspose_add,
    Matrix.conjTranspose_smul, intLift_conjTranspose, intLift_conjTranspose,
    star_sqrt3C]
  norm_num

/-- Unitarity receipt: `Uᴴ * U = 1` for every image, lifted from the
kernel-decided orthogonality components. -/
theorem gaugeRep_conjTranspose_mul_self (g : S3) :
    (gaugeRep g)ᴴ * gaugeRep g = 1 := by
  obtain ⟨h1, h2⟩ := gaugePair_orth g
  rw [gaugeRep_conjTranspose, gaugeRep, intLift_pair_mul, h1, h2]
  rw [map_zero, smul_zero, add_zero]
  have h4 : intLift !![4, 0; 0, 4] = (4 : ℂ) • 1 := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [intLift, Matrix.map_apply]
  rw [h4, smul_smul]
  norm_num

/-- Unitarity receipt, right-handed form. -/
theorem gaugeRep_mul_conjTranspose_self (g : S3) :
    gaugeRep g * (gaugeRep g)ᴴ = 1 :=
  mul_eq_one_comm.mp (gaugeRep_conjTranspose_mul_self g)

theorem gaugeRep_cancel_left (g : S3) (Z : Matrix (Fin 2) (Fin 2) ℂ) :
    (gaugeRep g)ᴴ * (gaugeRep g * Z) = Z := by
  rw [← Matrix.mul_assoc, gaugeRep_conjTranspose_mul_self, one_mul]

/-! ## The six matrices, entrywise -/

theorem gaugeRep_one_eq : gaugeRep 1 = 1 := gaugeRep_one

theorem gaugeRep_r_eq :
    gaugeRep rGen
      = !![-(2⁻¹ : ℂ), -(2⁻¹ * sqrt3C); 2⁻¹ * sqrt3C, -(2⁻¹ : ℂ)] := by
  rw [gaugeRep, gaugePair_r]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [intLift, Matrix.map_apply]

theorem gaugeRep_rr_eq :
    gaugeRep (rGen * rGen)
      = !![-(2⁻¹ : ℂ), 2⁻¹ * sqrt3C; -(2⁻¹ * sqrt3C), -(2⁻¹ : ℂ)] := by
  rw [gaugeRep, gaugePair_rr]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [intLift, Matrix.map_apply]

theorem gaugeRep_t_eq : gaugeRep tGen = !![(1 : ℂ), 0; 0, -1] := by
  rw [gaugeRep, gaugePair_t]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [intLift, Matrix.map_apply]

theorem gaugeRep_tr_eq :
    gaugeRep (tGen * rGen)
      = !![-(2⁻¹ : ℂ), -(2⁻¹ * sqrt3C); -(2⁻¹ * sqrt3C), (2⁻¹ : ℂ)] := by
  rw [gaugeRep, gaugePair_tr]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [intLift, Matrix.map_apply]

theorem gaugeRep_trr_eq :
    gaugeRep (tGen * (rGen * rGen))
      = !![-(2⁻¹ : ℂ), 2⁻¹ * sqrt3C; 2⁻¹ * sqrt3C, (2⁻¹ : ℂ)] := by
  rw [gaugeRep, gaugePair_trr]
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [intLift, Matrix.map_apply]

/-! ## Irreducibility: the finite Schur form -/

/-- Finite Schur form, generator version: a matrix commuting with the
order-3 rotation image and one reflection image is scalar.  Direct
linear algebra: the reflection kills the off-diagonal entries, the
rotation identifies the diagonal ones through `sqrt 3 ≠ 0`. -/
theorem commutant_pair_scalar {N : Matrix (Fin 2) (Fin 2) ℂ}
    (hr : N * gaugeRep rGen = gaugeRep rGen * N)
    (ht : N * gaugeRep tGen = gaugeRep tGen * N) :
    N = N 0 0 • 1 := by
  rw [gaugeRep_t_eq] at ht
  rw [gaugeRep_r_eq] at hr
  have h01 : N 0 1 = 0 := by
    have h := congrFun (congrFun ht 0) 1
    simp [Matrix.mul_apply, Fin.sum_univ_two] at h
    linear_combination (-(1 : ℂ) / 2) * h
  have h10 : N 1 0 = 0 := by
    have h := congrFun (congrFun ht 1) 0
    simp [Matrix.mul_apply, Fin.sum_univ_two] at h
    linear_combination ((1 : ℂ) / 2) * h
  have hdiag : N 1 1 = N 0 0 := by
    have h := congrFun (congrFun hr 0) 1
    simp [Matrix.mul_apply, Fin.sum_univ_two, h01] at h
    have h' : sqrt3C * N 1 1 = sqrt3C * N 0 0 := by
      linear_combination (-2 : ℂ) * h
    exact mul_left_cancel₀ sqrt3C_ne_zero h'
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [h01, h10, hdiag]

/-- Finite Schur form: every complex `2 x 2` matrix commuting with all
six images of the representation is scalar. -/
theorem commutant_scalar {N : Matrix (Fin 2) (Fin 2) ℂ}
    (h : ∀ g : S3, N * gaugeRep g = gaugeRep g * N) :
    ∃ c : ℂ, N = c • 1 :=
  ⟨N 0 0, commutant_pair_scalar (h rGen) (h tGen)⟩

/-! ## The gauge action on states and effects -/

/-- Conjugation by a gauge image. -/
noncomputable def conjBy (g : S3) (M : Matrix (Fin 2) (Fin 2) ℂ) :
    Matrix (Fin 2) (Fin 2) ℂ :=
  gaugeRep g * M * (gaugeRep g)ᴴ

theorem conjBy_conjBy (a b : S3) (M : Matrix (Fin 2) (Fin 2) ℂ) :
    conjBy a (conjBy b M) = conjBy (a * b) M := by
  simp only [conjBy, gaugeRep_mul, Matrix.conjTranspose_mul, Matrix.mul_assoc]

theorem conjBy_smul (g : S3) (c : ℂ) (M : Matrix (Fin 2) (Fin 2) ℂ) :
    conjBy g (c • M) = c • conjBy g M := by
  simp [conjBy]

theorem conjBy_add (g : S3) (M N : Matrix (Fin 2) (Fin 2) ℂ) :
    conjBy g (M + N) = conjBy g M + conjBy g N := by
  simp [conjBy, Matrix.mul_add, Matrix.add_mul]

theorem conjBy_mul (g : S3) (M N : Matrix (Fin 2) (Fin 2) ℂ) :
    conjBy g M * conjBy g N = conjBy g (M * N) := by
  simp only [conjBy, Matrix.mul_assoc, gaugeRep_cancel_left]

theorem trace_conjBy (g : S3) (M : Matrix (Fin 2) (Fin 2) ℂ) :
    (conjBy g M).trace = M.trace := by
  rw [conjBy, Matrix.trace_mul_cycle, gaugeRep_conjTranspose_mul_self, one_mul]

/-- A matrix fixed by conjugation commutes with the conjugating
image. -/
theorem comm_of_conjBy_fixed {g : S3} {M : Matrix (Fin 2) (Fin 2) ℂ}
    (h : conjBy g M = M) : M * gaugeRep g = gaugeRep g * M := by
  conv_lhs => rw [← h]
  show gaugeRep g * M * (gaugeRep g)ᴴ * gaugeRep g = gaugeRep g * M
  rw [Matrix.mul_assoc (gaugeRep g * M), gaugeRep_conjTranspose_mul_self, mul_one]

/-! ## The invariant state -/

/-- The tracial state on `M₂`: one half of the identity. -/
noncomputable def tracialState : Matrix (Fin 2) (Fin 2) ℂ := (2⁻¹ : ℂ) • 1

/-- The tracial state is a certified state. -/
theorem tracialState_isState : IsState tracialState := by
  constructor
  · refine isEvent_one.posSemidef.smul (a := (2⁻¹ : ℂ)) ?_
    rw [show (2⁻¹ : ℂ) = ((2⁻¹ : ℝ) : ℂ) by norm_num]
    exact_mod_cast (by norm_num : (0 : ℝ) ≤ 2⁻¹)
  · rw [tracialState, Matrix.trace_smul, Matrix.trace_one]
    simp [Fintype.card_fin]

/-- The tracial state is invariant under the whole gauge action. -/
theorem conjBy_tracialState (g : S3) : conjBy g tracialState = tracialState := by
  rw [conjBy, tracialState, Matrix.mul_smul, mul_one, Matrix.smul_mul,
    gaugeRep_mul_conjTranspose_self]

/-- Uniqueness: every density matrix invariant under conjugation by all
six images is the tracial state.  Schur makes it scalar; the trace
makes the scalar one half. -/
theorem invariant_state_eq_tracialState {ρ : Matrix (Fin 2) (Fin 2) ℂ}
    (hρ : IsState ρ) (hinv : ∀ g : S3, conjBy g ρ = ρ) : ρ = tracialState := by
  have hsc : ρ = ρ 0 0 • 1 := commutant_pair_scalar
    (comm_of_conjBy_fixed (hinv rGen)) (comm_of_conjBy_fixed (hinv tGen))
  have htr := hρ.2
  rw [hsc, Matrix.trace_smul, Matrix.trace_one, Fintype.card_fin] at htr
  have hval : ρ 0 0 = 2⁻¹ := by
    have h : ρ 0 0 * 2 = 1 := by simpa [smul_eq_mul] using htr
    linear_combination (2⁻¹ : ℂ) * h
  rw [hsc, hval, tracialState]

/-- The Born weight of the tracial state is half the trace. -/
theorem bornWeight_tracialState (E : Matrix (Fin 2) (Fin 2) ℂ) :
    bornWeight tracialState E = ((2⁻¹ : ℝ) : ℂ) * E.trace := by
  rw [bornWeight, tracialState, Matrix.smul_mul, one_mul, Matrix.trace_smul,
    smul_eq_mul]
  norm_num

/-- Display form: the Born weight of the tracial state is the trace
pairing against `(1/2) • 1`, definitionally. -/
theorem bornWeight_tracialState_def (E : Matrix (Fin 2) (Fin 2) ℂ) :
    bornWeight tracialState E
      = ((2⁻¹ : ℂ) • (1 : Matrix (Fin 2) (Fin 2) ℂ) * E).trace := rfl

/-! ## The twirl over the six realized gauge elements -/

/-- The equal-weight average of an effect over the six explicit gauge
conjugates.  The six group words enumerate `S3`. -/
noncomputable def twirl (E : Matrix (Fin 2) (Fin 2) ℂ) :
    Matrix (Fin 2) (Fin 2) ℂ :=
  (6 : ℂ)⁻¹ • (conjBy 1 E + conjBy rGen E + conjBy (rGen * rGen) E
    + conjBy tGen E + conjBy (tGen * rGen) E + conjBy (tGen * (rGen * rGen)) E)

/-- Left multiplication by the rotation permutes the six twirl words. -/
theorem conjBy_r_twirl (E : Matrix (Fin 2) (Fin 2) ℂ) :
    conjBy rGen (twirl E) = twirl E := by
  have g3 : rGen * (rGen * rGen) = 1 := by decide
  have g4 : rGen * tGen = tGen * (rGen * rGen) := by decide
  have g5 : rGen * (tGen * rGen) = tGen := by decide
  have g6 : rGen * (tGen * (rGen * rGen)) = tGen * rGen := by decide
  simp only [twirl, conjBy_smul, conjBy_add, conjBy_conjBy, mul_one,
    g3, g4, g5, g6]
  congr 1
  abel

/-- Left multiplication by the reflection permutes the six twirl
words. -/
theorem conjBy_t_twirl (E : Matrix (Fin 2) (Fin 2) ℂ) :
    conjBy tGen (twirl E) = twirl E := by
  have g4 : tGen * tGen = 1 := by decide
  have g5 : tGen * (tGen * rGen) = rGen := by decide
  have g6 : tGen * (tGen * (rGen * rGen)) = rGen * rGen := by decide
  simp only [twirl, conjBy_smul, conjBy_add, conjBy_conjBy, mul_one,
    g4, g5, g6]
  congr 1
  abel

theorem twirl_comm_r (E : Matrix (Fin 2) (Fin 2) ℂ) :
    twirl E * gaugeRep rGen = gaugeRep rGen * twirl E :=
  comm_of_conjBy_fixed (conjBy_r_twirl E)

theorem twirl_comm_t (E : Matrix (Fin 2) (Fin 2) ℂ) :
    twirl E * gaugeRep tGen = gaugeRep tGen * twirl E :=
  comm_of_conjBy_fixed (conjBy_t_twirl E)

/-- The twirl preserves the trace: each conjugate has the trace of the
original by cyclicity and unitarity. -/
theorem twirl_trace (E : Matrix (Fin 2) (Fin 2) ℂ) :
    (twirl E).trace = E.trace := by
  simp only [twirl, Matrix.trace_smul, Matrix.trace_add, trace_conjBy,
    smul_eq_mul]
  ring

/-- On projection events the twirl is the Born-weight multiple of the
identity: scalar by Schur, evaluated by the trace, real by
hermiticity. -/
theorem twirl_event_eq {E : Matrix (Fin 2) (Fin 2) ℂ} (hE : IsEvent E) :
    twirl E = (((bornWeight tracialState E).re : ℝ) : ℂ) • 1 := by
  have hsc : twirl E = twirl E 0 0 • 1 :=
    commutant_pair_scalar (twirl_comm_r E) (twirl_comm_t E)
  have htr : twirl E 0 0 * 2 = E.trace := by
    have h := twirl_trace E
    rw [hsc, Matrix.trace_smul, Matrix.trace_one, Fintype.card_fin] at h
    simpa [smul_eq_mul] using h
  have hreal : ((E.trace.re : ℝ) : ℂ) = E.trace := by
    rw [← Complex.conj_eq_iff_re]
    show star E.trace = E.trace
    rw [← Matrix.trace_conjTranspose, hE.1.eq]
  have hborn : (((bornWeight tracialState E).re : ℝ) : ℂ) = twirl E 0 0 := by
    rw [bornWeight_tracialState, Complex.re_ofReal_mul]
    push_cast
    rw [hreal]
    linear_combination (-(2⁻¹ : ℂ)) * htr
  rw [hsc, hborn]

/-! ## Binary mixes generate the six-average -/

/-- Five binary mixes in the bridge's mixing grammar produce the
equal-weight six-average of outcome values. -/
private theorem outcome_six_average (w : Matrix (Fin 2) (Fin 2) ℂ → ℝ)
    (hmix : ∀ (l : ℝ), 0 ≤ l → l ≤ 1 → ∀ X Y : Matrix (Fin 2) (Fin 2) ℂ,
      w ((l : ℂ) • X + ((1 - l : ℝ) : ℂ) • Y) = l * w X + (1 - l) * w Y)
    (X₁ X₂ X₃ X₄ X₅ X₆ : Matrix (Fin 2) (Fin 2) ℂ) :
    w ((6 : ℂ)⁻¹ • (X₁ + X₂ + X₃ + X₄ + X₅ + X₆))
      = 6⁻¹ * (w X₁ + w X₂ + w X₃ + w X₄ + w X₅ + w X₆) := by
  have h2 : ∀ Y Z : Matrix (Fin 2) (Fin 2) ℂ,
      w ((2 : ℂ)⁻¹ • (Y + Z)) = 2⁻¹ * w Y + 2⁻¹ * w Z := by
    intro Y Z
    have e : ((2⁻¹ : ℝ) : ℂ) • Y + ((1 - 2⁻¹ : ℝ) : ℂ) • Z
        = (2 : ℂ)⁻¹ • (Y + Z) := by
      push_cast
      module
    rw [← e, hmix 2⁻¹ (by norm_num) (by norm_num)]
    norm_num
  have h3 : ∀ Y Z W : Matrix (Fin 2) (Fin 2) ℂ,
      w ((3 : ℂ)⁻¹ • (Y + Z + W)) = 3⁻¹ * w Y + 3⁻¹ * w Z + 3⁻¹ * w W := by
    intro Y Z W
    have e : ((3⁻¹ : ℝ) : ℂ) • Y + ((1 - 3⁻¹ : ℝ) : ℂ) • ((2 : ℂ)⁻¹ • (Z + W))
        = (3 : ℂ)⁻¹ • (Y + Z + W) := by
      push_cast
      module
    rw [← e, hmix 3⁻¹ (by norm_num) (by norm_num), h2]
    ring
  have e : ((2⁻¹ : ℝ) : ℂ) • ((3 : ℂ)⁻¹ • (X₁ + X₂ + X₃))
        + ((1 - 2⁻¹ : ℝ) : ℂ) • ((3 : ℂ)⁻¹ • (X₄ + X₅ + X₆))
      = (6 : ℂ)⁻¹ • (X₁ + X₂ + X₃ + X₄ + X₅ + X₆) := by
    push_cast
    module
  rw [← e, hmix 2⁻¹ (by norm_num) (by norm_num), h3, h3]
  ring

/-! ## The derivation statement -/

/-- **The gauge route to the Born weight.**  Any outcome assignment `v`
on pairs of a state and an effect that

* (a) mixes affinely in the effect slot at the invariant state, in the
  binary-mixing grammar of the convex-affinity bridge,
* (b) is covariant under the realized gauge action, and
* (c) is normalized at the sure and impossible effects,

takes on every projection event the Born weight of the tracial state.
The twirl of the effect is reached through five binary mixes and fixed
values, is scalar by the finite Schur form, and the scalar is evaluated
by hypothesis (a) at the endpoints `1` and `0`. -/
theorem born_weight_forced_at_invariant_state
    (v : Matrix (Fin 2) (Fin 2) ℂ → Matrix (Fin 2) (Fin 2) ℂ → ℝ)
    (hmix : ∀ (l : ℝ), 0 ≤ l → l ≤ 1 → ∀ X Y : Matrix (Fin 2) (Fin 2) ℂ,
      v tracialState ((l : ℂ) • X + ((1 - l : ℝ) : ℂ) • Y)
        = l * v tracialState X + (1 - l) * v tracialState Y)
    (hcov : ∀ (g : S3) (ρ E : Matrix (Fin 2) (Fin 2) ℂ),
      v (conjBy g ρ) (conjBy g E) = v ρ E)
    (hone : v tracialState 1 = 1) (hzero : v tracialState 0 = 0)
    {E : Matrix (Fin 2) (Fin 2) ℂ} (hE : IsEvent E) :
    v tracialState E = (bornWeight tracialState E).re := by
  have hfix : ∀ (g : S3) (F : Matrix (Fin 2) (Fin 2) ℂ),
      v tracialState (conjBy g F) = v tracialState F := by
    intro g F
    have h := hcov g tracialState F
    rwa [conjBy_tracialState] at h
  have h1 : v tracialState (twirl E) = v tracialState E := by
    simp only [twirl]
    rw [outcome_six_average (v tracialState) hmix]
    simp only [hfix]
    ring
  have hc0 : 0 ≤ (bornWeight tracialState E).re :=
    bornWeight_re_nonneg tracialState_isState.1 hE
  have hc1 : (bornWeight tracialState E).re ≤ 1 :=
    bornWeight_re_le_one tracialState_isState hE
  have h2 : v tracialState ((((bornWeight tracialState E).re : ℝ) : ℂ) • 1)
      = (bornWeight tracialState E).re := by
    have h := hmix ((bornWeight tracialState E).re) hc0 hc1 1 0
    rw [smul_zero, add_zero] at h
    rw [h, hone, hzero]
    ring
  rw [← h1, twirl_event_eq hE, h2]

/-! ## The source tie: the realized edge census -/

/-- The per-element incident edge counts of the designated split fibre
(observer 64, record class 18, nodes 64 and 43), transcribed from the
`edge_gauge_extraction.per_split_fibre` block of the committed payload.
Order of the literals follows the payload index order under the
permutation dictionary: index 0 is the identity, 1 is `tGen`, 2 is
`tGen * (rGen * rGen)`, 3 is `rGen`, 4 is `rGen * rGen`, 5 is
`tGen * rGen`. -/
def designatedFibreEdgeCount (g : S3) : ℕ :=
  if g = 1 then 8
  else if g = tGen then 4
  else if g = tGen * (rGen * rGen) then 3
  else if g = rGen then 2
  else if g = rGen * rGen then 3
  else 1

/-- Receipt-shaped interface Prop for the source tie: an edge census
realizes the gauge action when every group element has positive count.
The committed payload is the intended producer; any module binding the
payload literals can discharge this Prop without this file importing
it. -/
def RealizesAllGaugeElements (count : S3 → ℕ) : Prop :=
  ∀ g : S3, 0 < count g

/-- The designated-fibre census realizes all six gauge elements. -/
theorem designatedFibre_realizes_all :
    RealizesAllGaugeElements designatedFibreEdgeCount := by
  show ∀ g : S3, 0 < designatedFibreEdgeCount g
  decide

set_option maxRecDepth 8000 in
/-- The census totals the payload's incident edge count. -/
theorem designatedFibre_edge_total :
    ∑ g : S3, designatedFibreEdgeCount g = 21 := by
  decide

/-- The derivation with covariance assumed only at gauge elements of
positive realized edge count: any census realizing all six elements
discharges the full covariance need, so the run's census is exactly
what the route consumes from the source. -/
theorem born_weight_forced_with_realized_covariance
    (v : Matrix (Fin 2) (Fin 2) ℂ → Matrix (Fin 2) (Fin 2) ℂ → ℝ)
    (hmix : ∀ (l : ℝ), 0 ≤ l → l ≤ 1 → ∀ X Y : Matrix (Fin 2) (Fin 2) ℂ,
      v tracialState ((l : ℂ) • X + ((1 - l : ℝ) : ℂ) • Y)
        = l * v tracialState X + (1 - l) * v tracialState Y)
    {count : S3 → ℕ} (hcount : RealizesAllGaugeElements count)
    (hcov : ∀ g : S3, 0 < count g → ∀ ρ E : Matrix (Fin 2) (Fin 2) ℂ,
      v (conjBy g ρ) (conjBy g E) = v ρ E)
    (hone : v tracialState 1 = 1) (hzero : v tracialState 0 = 0)
    {E : Matrix (Fin 2) (Fin 2) ℂ} (hE : IsEvent E) :
    v tracialState E = (bornWeight tracialState E).re :=
  born_weight_forced_at_invariant_state v hmix
    (fun g => hcov g (hcount g)) hone hzero hE

/-! ## The Born assignment inhabits the hypotheses -/

/-- The trace-pairing assignment. -/
noncomputable def bornAssignment (ρ E : Matrix (Fin 2) (Fin 2) ℂ) : ℝ :=
  (bornWeight ρ E).re

/-- Bridge consumption by instantiation: on every convex preparation
family the Born assignment satisfies the state-slot mixing law, and the
receipt is the bridge theorem `ConvexPreparationFamily.outcome_re_affine`
itself; `segmentFamily` and `towerSegmentFamily` supply the committed
producers. -/
theorem bornAssignment_family_mix {P : Type} (F : ConvexPreparationFamily 2 P)
    (E : Matrix (Fin 2) (Fin 2) ℂ) (l : ℝ) (hl : l ∈ Set.Icc (0 : ℝ) 1)
    (p q : P) :
    bornAssignment (F.state (F.mix l hl p q)) E
      = l * bornAssignment (F.state p) E + (1 - l) * bornAssignment (F.state q) E :=
  F.outcome_re_affine E l hl p q

/-- The Born assignment mixes affinely in the effect slot. -/
theorem bornAssignment_effect_mix (l : ℝ) (X Y : Matrix (Fin 2) (Fin 2) ℂ) :
    bornAssignment tracialState ((l : ℂ) • X + ((1 - l : ℝ) : ℂ) • Y)
      = l * bornAssignment tracialState X + (1 - l) * bornAssignment tracialState Y := by
  simp only [bornAssignment, bornWeight, Matrix.mul_add, Matrix.mul_smul,
    Matrix.trace_add, Matrix.trace_smul, smul_eq_mul]
  simp [Complex.add_re]

/-- The Born assignment is covariant under the gauge action. -/
theorem bornAssignment_covariant (g : S3) (ρ E : Matrix (Fin 2) (Fin 2) ℂ) :
    bornAssignment (conjBy g ρ) (conjBy g E) = bornAssignment ρ E := by
  simp only [bornAssignment, bornWeight, conjBy_mul, trace_conjBy]

theorem bornAssignment_norm_one : bornAssignment tracialState 1 = 1 := by
  simp [bornAssignment, bornWeight_one tracialState_isState]

theorem bornAssignment_norm_zero : bornAssignment tracialState 0 = 0 := by
  simp [bornAssignment, bornWeight]

/-- Non-vacuity: the Born assignment satisfies every hypothesis of the
derivation, and the theorem returns its own value on it. -/
theorem bornAssignment_forced {E : Matrix (Fin 2) (Fin 2) ℂ} (hE : IsEvent E) :
    bornAssignment tracialState E = (bornWeight tracialState E).re :=
  born_weight_forced_at_invariant_state bornAssignment
    (fun l _ _ X Y => bornAssignment_effect_mix l X Y)
    bornAssignment_covariant bornAssignment_norm_one bornAssignment_norm_zero hE

/-! ## Witness computations: exact halves on a noncommuting pair -/

/-- The two bridge witness projectors fail to commute. -/
theorem vertexZero_vertexPlus_noncomm :
    vertexZero * vertexPlus ≠ vertexPlus * vertexZero := by
  intro h
  have h01 := congrFun (congrFun h 0) 1
  simp [vertexZero, vertexPlus, Matrix.mul_apply, Fin.sum_univ_two] at h01

/-- Exact half: the Born weight of the tracial state on the first
witness projector. -/
theorem tracial_born_vertexZero :
    (bornWeight tracialState vertexZero).re = 2⁻¹ := by
  rw [bornWeight_tracialState, trace_vertexZero, mul_one, Complex.ofReal_re]

/-- Exact half: the Born weight of the tracial state on the second
witness projector. -/
theorem tracial_born_vertexPlus :
    (bornWeight tracialState vertexPlus).re = 2⁻¹ := by
  rw [bornWeight_tracialState, trace_vertexPlus, mul_one, Complex.ofReal_re]

/-- Every admissible assignment takes the exact half on both
noncommuting witness projectors. -/
theorem admissible_values_on_noncommuting_pair
    (v : Matrix (Fin 2) (Fin 2) ℂ → Matrix (Fin 2) (Fin 2) ℂ → ℝ)
    (hmix : ∀ (l : ℝ), 0 ≤ l → l ≤ 1 → ∀ X Y : Matrix (Fin 2) (Fin 2) ℂ,
      v tracialState ((l : ℂ) • X + ((1 - l : ℝ) : ℂ) • Y)
        = l * v tracialState X + (1 - l) * v tracialState Y)
    (hcov : ∀ (g : S3) (ρ E : Matrix (Fin 2) (Fin 2) ℂ),
      v (conjBy g ρ) (conjBy g E) = v ρ E)
    (hone : v tracialState 1 = 1) (hzero : v tracialState 0 = 0) :
    v tracialState vertexZero = 2⁻¹ ∧ v tracialState vertexPlus = 2⁻¹ := by
  constructor
  · rw [born_weight_forced_at_invariant_state v hmix hcov hone hzero
      isEvent_vertexZero, tracial_born_vertexZero]
  · rw [born_weight_forced_at_invariant_state v hmix hcov hone hzero
      isEvent_vertexPlus, tracial_born_vertexPlus]

/-! ## Negative control: covariance is load-bearing -/

/-- The corner assignment reads the `(0,0)` entry of the effect and
ignores the state. -/
noncomputable def cornerAssignment
    (_ρ E : Matrix (Fin 2) (Fin 2) ℂ) : ℝ :=
  (E 0 0).re

theorem cornerAssignment_effect_mix (l : ℝ) (X Y : Matrix (Fin 2) (Fin 2) ℂ) :
    cornerAssignment tracialState ((l : ℂ) • X + ((1 - l : ℝ) : ℂ) • Y)
      = l * cornerAssignment tracialState X
        + (1 - l) * cornerAssignment tracialState Y := by
  simp only [cornerAssignment, Matrix.add_apply, Matrix.smul_apply,
    smul_eq_mul, Complex.add_re, Complex.re_ofReal_mul]

theorem cornerAssignment_norm_one : cornerAssignment tracialState 1 = 1 := by
  simp [cornerAssignment, Matrix.one_apply_eq]

theorem cornerAssignment_norm_zero : cornerAssignment tracialState 0 = 0 := by
  simp [cornerAssignment]

/-- The corner assignment breaks covariance on the realized rotation:
the rotated first witness projector carries corner value `1/4` against
the unrotated value `1`. -/
theorem cornerAssignment_breaks_covariance :
    cornerAssignment (conjBy rGen tracialState) (conjBy rGen vertexZero)
      ≠ cornerAssignment tracialState vertexZero := by
  have hval : ((conjBy rGen vertexZero) 0 0).re = 4⁻¹ := by
    rw [conjBy, gaugeRep_r_eq]
    simp [vertexZero, Matrix.mul_apply, Fin.sum_univ_two,
      Matrix.conjTranspose_apply]
    norm_num
  have hcorner : ((vertexZero : Matrix (Fin 2) (Fin 2) ℂ) 0 0).re = 1 := by
    simp [vertexZero]
  simp only [cornerAssignment, hval, hcorner]
  norm_num

/-- The corner assignment evades the conclusion on the first witness
projector: value `1` against the Born weight `1/2`. -/
theorem cornerAssignment_evades_born :
    cornerAssignment tracialState vertexZero
      ≠ (bornWeight tracialState vertexZero).re := by
  rw [tracial_born_vertexZero]
  simp [cornerAssignment, vertexZero]

/-! ## Negative control: effect-slot mixing is load-bearing -/

/-- The squared-trace assignment: state-affine (constant in the state),
covariant, and normalized, with a nonlinear effect response. -/
noncomputable def traceSquareAssignment
    (_ρ E : Matrix (Fin 2) (Fin 2) ℂ) : ℝ :=
  (E.trace.re * 2⁻¹) ^ 2

/-- State-slot mixing holds for the squared-trace assignment. -/
theorem traceSquareAssignment_state_mix (l : ℝ)
    (ρ σ E : Matrix (Fin 2) (Fin 2) ℂ) :
    traceSquareAssignment ((l : ℂ) • ρ + ((1 - l : ℝ) : ℂ) • σ) E
      = l * traceSquareAssignment ρ E + (1 - l) * traceSquareAssignment σ E := by
  simp only [traceSquareAssignment]
  ring

/-- Covariance holds for the squared-trace assignment. -/
theorem traceSquareAssignment_covariant (g : S3)
    (ρ E : Matrix (Fin 2) (Fin 2) ℂ) :
    traceSquareAssignment (conjBy g ρ) (conjBy g E)
      = traceSquareAssignment ρ E := by
  simp [traceSquareAssignment, trace_conjBy]

theorem traceSquareAssignment_norm_one :
    traceSquareAssignment tracialState 1 = 1 := by
  simp [traceSquareAssignment, Matrix.trace_one, Fintype.card_fin]

theorem traceSquareAssignment_norm_zero :
    traceSquareAssignment tracialState 0 = 0 := by
  simp [traceSquareAssignment]

/-- Effect-slot mixing fails for the squared-trace assignment at the
midpoint of the sure and impossible effects. -/
theorem traceSquareAssignment_not_effect_mix :
    traceSquareAssignment tracialState
        (((2⁻¹ : ℝ) : ℂ) • 1 + ((1 - 2⁻¹ : ℝ) : ℂ) • 0)
      ≠ 2⁻¹ * traceSquareAssignment tracialState 1
        + (1 - 2⁻¹) * traceSquareAssignment tracialState 0 := by
  simp [traceSquareAssignment, Matrix.trace_smul, Matrix.trace_one,
    Fintype.card_fin]
  norm_num

/-- The squared-trace assignment evades the conclusion on the first
witness projector: value `1/4` against the Born weight `1/2`.  With
`traceSquareAssignment_state_mix` and
`traceSquareAssignment_covariant` this shows state-slot affinity plus
covariance plus normalization cannot force the Born weight; the
effect-slot hypothesis of the derivation is sharp. -/
theorem traceSquareAssignment_evades_born :
    traceSquareAssignment tracialState vertexZero
      ≠ (bornWeight tracialState vertexZero).re := by
  rw [tracial_born_vertexZero]
  simp [traceSquareAssignment, trace_vertexZero]
  norm_num

end OPH.QFT.GaugeBorn

-- Axiom audit: each must report nothing beyond
-- `[propext, Classical.choice, Quot.sound]`.
#print axioms OPH.QFT.GaugeBorn.gaugePair_mulTable
#print axioms OPH.QFT.GaugeBorn.gaugePair_orth
#print axioms OPH.QFT.GaugeBorn.gaugeRep_mul
#print axioms OPH.QFT.GaugeBorn.gaugeRep_conjTranspose_mul_self
#print axioms OPH.QFT.GaugeBorn.gaugeRep_mul_conjTranspose_self
#print axioms OPH.QFT.GaugeBorn.gaugeRep_r_eq
#print axioms OPH.QFT.GaugeBorn.gaugeRep_t_eq
#print axioms OPH.QFT.GaugeBorn.commutant_pair_scalar
#print axioms OPH.QFT.GaugeBorn.commutant_scalar
#print axioms OPH.QFT.GaugeBorn.invariant_state_eq_tracialState
#print axioms OPH.QFT.GaugeBorn.twirl_event_eq
#print axioms OPH.QFT.GaugeBorn.born_weight_forced_at_invariant_state
#print axioms OPH.QFT.GaugeBorn.born_weight_forced_with_realized_covariance
#print axioms OPH.QFT.GaugeBorn.designatedFibre_realizes_all
#print axioms OPH.QFT.GaugeBorn.designatedFibre_edge_total
#print axioms OPH.QFT.GaugeBorn.bornAssignment_family_mix
#print axioms OPH.QFT.GaugeBorn.bornAssignment_forced
#print axioms OPH.QFT.GaugeBorn.vertexZero_vertexPlus_noncomm
#print axioms OPH.QFT.GaugeBorn.tracial_born_vertexZero
#print axioms OPH.QFT.GaugeBorn.tracial_born_vertexPlus
#print axioms OPH.QFT.GaugeBorn.admissible_values_on_noncommuting_pair
#print axioms OPH.QFT.GaugeBorn.cornerAssignment_breaks_covariance
#print axioms OPH.QFT.GaugeBorn.cornerAssignment_evades_born
#print axioms OPH.QFT.GaugeBorn.traceSquareAssignment_not_effect_mix
#print axioms OPH.QFT.GaugeBorn.traceSquareAssignment_evades_born
