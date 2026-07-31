import Mathlib
import A5IncidenceResponse

namespace OPH.A5AngularKernels

open OPH.A5IncidenceResponse
  (identityMatrix adjacencyMatrix antipodeMatrix matMul)

/-!
# The four Legendre Gram kernels on the twelve ports

The twelve icosahedral ports embed in the sphere with pairwise dot products
in `{1, 1/√5, -1/√5, -1}`, so the Legendre Gram kernel
`B_ℓ(i, j) = P_ℓ(v_i · v_j)` lives in the adjacency algebra
`span{I, A, A′, Π}` over `ℚ(√5)`, where `A′ = AΠ`.  Writing
`5·B_ℓ = X_ℓ + Y_ℓ·√5` with integral `X_ℓ, Y_ℓ`, the exact components are

* `X₀ = 5J`, `Y₀ = 0` (with `J` the all-ones matrix),
* `X₁ = 5I - 5Π`, `Y₁ = A - A′`,
* `X₂ = 5I - A - A′ + 5Π`, `Y₂ = 0`,
* `X₃ = 5I - 5Π`, `Y₃ = A′ - A`.

This file proves, by kernel computation over `ℤ`:

* the dot-product partition `I + A + A′ + Π = J`;
* the four kernels are mutually orthogonal scaled projectors with scales
  `(12, 4, 12/5, 4)` — in integral form `(5B_ℓ)² = 5·c_ℓ·(5B_ℓ)` with
  `5·c_ℓ ∈ {60, 20, 12, 20}`;
* the traces `tr X_ℓ = 60`, `tr Y_ℓ = 0`, which force the ranks
  `(1, 3, 5, 3)` given the scales;
* the resolution of the identity
  `B₀/12 + B₁/4 + 5B₂/12 + B₃/4 = I`;
* the band binding `A B₁ = √5 B₁`, `A B₂ = -B₂`, `A B₃ = -√5 B₃`,
  `A B₀ = 5 B₀`;
* the antipodal parity `Π B_ℓ = (-1)^ℓ B_ℓ`, hence the alternating
  readback parity `(-Π) B_ℓ = (-1)^(ℓ+1) B_ℓ`.

A product in `ℤ[√5]` splits as
`(X + Y√5)(X′ + Y′√5) = (XX′ + 5YY′) + (XY′ + YX′)√5`, so every identity
below is a pair of integral matrix identities checked by `decide`.
-/

/-- The all-ones matrix on the twelve ports. -/
def allOnes (_ _ : Fin 12) : ℤ := 1

/-- The second-neighbor permutation kernel `A′ = AΠ`. -/
def secondMatrix : Fin 12 → Fin 12 → ℤ :=
  matMul adjacencyMatrix antipodeMatrix

/-- Rational component of `5·B₀`. -/
def x0 (i j : Fin 12) : ℤ := 5 * allOnes i j

/-- Rational component of `5·B₁` (and of `5·B₃`). -/
def x1 (i j : Fin 12) : ℤ :=
  5 * identityMatrix i j - 5 * antipodeMatrix i j

/-- Rational component of `5·B₂`. -/
def x2 (i j : Fin 12) : ℤ :=
  5 * identityMatrix i j - adjacencyMatrix i j - secondMatrix i j +
    5 * antipodeMatrix i j

/-- `√5` component of `5·B₁`. -/
def y1 (i j : Fin 12) : ℤ := adjacencyMatrix i j - secondMatrix i j

/-- The dot-product partition: `I + A + A′ + Π = J`. -/
theorem partition :
    ∀ i j : Fin 12,
      identityMatrix i j + adjacencyMatrix i j + secondMatrix i j +
        antipodeMatrix i j = allOnes i j := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- `B₀` is a scaled projector: `(5B₀)² = 60·(5B₀)`, scale `12`. -/
theorem projector_zero :
    ∀ i j : Fin 12, matMul x0 x0 i j = 60 * x0 i j := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- `B₁` is a scaled projector, scale `4`: rational part of
`(5B₁)² = 20·(5B₁)`. -/
theorem projector_one_rat :
    ∀ i j : Fin 12,
      matMul x1 x1 i j + 5 * matMul y1 y1 i j = 20 * x1 i j := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- `B₁` is a scaled projector, scale `4`: `√5` part of
`(5B₁)² = 20·(5B₁)`. -/
theorem projector_one_irr :
    ∀ i j : Fin 12,
      matMul x1 y1 i j + matMul y1 x1 i j = 20 * y1 i j := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- `B₂` is a scaled projector: `(5B₂)² = 12·(5B₂)`, scale `12/5`. -/
theorem projector_two :
    ∀ i j : Fin 12, matMul x2 x2 i j = 12 * x2 i j := by
  decide

/-- `5B₃ = X₁ - Y₁√5`, so `(5B₃)² = 20·(5B₃)` follows from the `B₁`
projector identities with the sign of the `√5` part flipped. -/
theorem projector_three_rat :
    ∀ i j : Fin 12,
      matMul x1 x1 i j + 5 * matMul y1 y1 i j = 20 * x1 i j :=
  projector_one_rat

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Orthogonality `B₀B₁ = 0`: rational part `X₀X₁ + 5·Y₀Y₁ = 0`. -/
theorem orthogonal_zero_one_rat :
    ∀ i j : Fin 12, matMul x0 x1 i j = 0 := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Orthogonality `B₀B₁ = 0`: `√5` part `X₀Y₁ + Y₀X₁ = 0`. -/
theorem orthogonal_zero_one_irr :
    ∀ i j : Fin 12, matMul x0 y1 i j = 0 := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Orthogonality `B₀B₂ = 0`. -/
theorem orthogonal_zero_two :
    ∀ i j : Fin 12, matMul x0 x2 i j = 0 := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Orthogonality `B₁B₂ = 0`: rational part. -/
theorem orthogonal_one_two_rat :
    ∀ i j : Fin 12, matMul x1 x2 i j = 0 := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Orthogonality `B₁B₂ = 0`: `√5` part `Y₁X₂ = 0`. -/
theorem orthogonal_one_two_irr :
    ∀ i j : Fin 12, matMul y1 x2 i j = 0 := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Orthogonality `B₁B₃ = 0`: with `5B₃ = X₁ - Y₁√5` the rational part
is `X₁X₁ - 5·Y₁Y₁ = 0`. -/
theorem orthogonal_one_three_rat :
    ∀ i j : Fin 12,
      matMul x1 x1 i j - 5 * matMul y1 y1 i j = 0 := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Orthogonality `B₁B₃ = 0`: `√5` part `Y₁X₁ - X₁Y₁ = 0`. -/
theorem orthogonal_one_three_irr :
    ∀ i j : Fin 12,
      matMul y1 x1 i j - matMul x1 y1 i j = 0 := by
  decide

/-- Traces: `tr X_ℓ = 60` and `tr Y_ℓ = 0` for every level, which with
the scales `(12, 4, 12/5, 4)` forces `tr B_ℓ = 12` and the ranks
`(1, 3, 5, 3)`. -/
theorem traces :
    ((List.finRange 12).map fun i => x0 i i).sum = 60 ∧
    ((List.finRange 12).map fun i => x1 i i).sum = 60 ∧
    ((List.finRange 12).map fun i => x2 i i).sum = 60 ∧
    ((List.finRange 12).map fun i => y1 i i).sum = 0 := by
  decide

/-- Resolution of the identity `B₀/12 + B₁/4 + 5B₂/12 + B₃/4 = I`,
cleared to the integral form `5X₀ + 30X₁ + 25X₂ = 300·I` (the `√5`
parts cancel as `15Y₁ + 15Y₃ = 0` by construction). -/
theorem resolution_of_identity :
    ∀ i j : Fin 12,
      5 * x0 i j + 30 * x1 i j + 25 * x2 i j =
        300 * identityMatrix i j := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Band binding `A B₀ = 5 B₀`. -/
theorem binding_zero :
    ∀ i j : Fin 12, matMul adjacencyMatrix x0 i j = 5 * x0 i j := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Band binding `A B₁ = √5 B₁`: rational part `A X₁ = 5 Y₁`. -/
theorem binding_one_rat :
    ∀ i j : Fin 12, matMul adjacencyMatrix x1 i j = 5 * y1 i j := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Band binding `A B₁ = √5 B₁`: `√5` part `A Y₁ = X₁`. -/
theorem binding_one_irr :
    ∀ i j : Fin 12, matMul adjacencyMatrix y1 i j = x1 i j := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Band binding `A B₂ = -B₂`. -/
theorem binding_two :
    ∀ i j : Fin 12, matMul adjacencyMatrix x2 i j = -x2 i j := by
  decide

/-- Band binding `A B₃ = -√5 B₃`: with `5B₃ = X₁ - Y₁√5` the rational
part `A X₁ = -5·(-Y₁)` and the `√5` part `A(-Y₁) = -X₁` are the `B₁`
binding identities with flipped signs. -/
theorem binding_three_rat :
    ∀ i j : Fin 12, matMul adjacencyMatrix x1 i j = -(5 * -y1 i j) := by
  intro i j
  have h := binding_one_rat i j
  omega

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Antipodal parity `Π B_ℓ = (-1)^ℓ B_ℓ`, hence the alternating
readback parity `(-Π) B_ℓ = (-1)^(ℓ+1) B_ℓ`: even levels. -/
theorem parity_even :
    ∀ i j : Fin 12,
      matMul antipodeMatrix x0 i j = x0 i j ∧
      matMul antipodeMatrix x2 i j = x2 i j := by
  decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- Antipodal parity, odd levels: `Π X₁ = -X₁` and `Π Y₁ = -Y₁`. -/
theorem parity_odd :
    ∀ i j : Fin 12,
      matMul antipodeMatrix x1 i j = -x1 i j ∧
      matMul antipodeMatrix y1 i j = -y1 i j := by
  decide

end OPH.A5AngularKernels

/- Axiom audit: kernel computation with standard axioms only. -/

#print axioms OPH.A5AngularKernels.partition
#print axioms OPH.A5AngularKernels.projector_zero
#print axioms OPH.A5AngularKernels.projector_one_rat
#print axioms OPH.A5AngularKernels.projector_one_irr
#print axioms OPH.A5AngularKernels.projector_two
#print axioms OPH.A5AngularKernels.orthogonal_zero_one_rat
#print axioms OPH.A5AngularKernels.orthogonal_zero_one_irr
#print axioms OPH.A5AngularKernels.orthogonal_zero_two
#print axioms OPH.A5AngularKernels.orthogonal_one_two_rat
#print axioms OPH.A5AngularKernels.orthogonal_one_two_irr
#print axioms OPH.A5AngularKernels.orthogonal_one_three_rat
#print axioms OPH.A5AngularKernels.orthogonal_one_three_irr
#print axioms OPH.A5AngularKernels.traces
#print axioms OPH.A5AngularKernels.resolution_of_identity
#print axioms OPH.A5AngularKernels.binding_zero
#print axioms OPH.A5AngularKernels.binding_one_rat
#print axioms OPH.A5AngularKernels.binding_one_irr
#print axioms OPH.A5AngularKernels.binding_two
#print axioms OPH.A5AngularKernels.parity_even
#print axioms OPH.A5AngularKernels.parity_odd
