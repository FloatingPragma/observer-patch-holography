import Mathlib
import A5IncidenceResponse

namespace OPH.A5AngularKernels

open OPH.A5IncidenceResponse
  (identityMatrix adjacencyMatrix antipodeMatrix matMul)

/-!
# Integral component identities for the four port kernels

This file defines integral 12 x 12 matrices on the port set: `allOnes`
(written `J` below), `secondMatrix` (`A' = A * Pi`, with `A` the
adjacency matrix and `Pi` the antipode from `A5IncidenceResponse`), and

* `x0 = 5*J`,
* `x1 = 5*I - 5*Pi`, `y1 = A - A'`,
* `x2 = 5*I - A - A' + 5*Pi`.

The declared reading is `5*B_l = X_l + Y_l*sqrt5` with
`(X_0, Y_0) = (x0, 0)`, `(X_1, Y_1) = (x1, y1)`, `(X_2, Y_2) = (x2, 0)`,
and `5*B_3 = x1 - y1*sqrt5`. A product in `Z[sqrt5]` splits as
`(X + Y*sqrt5)(X' + Y'*sqrt5) = (XX' + 5YY') + (XY' + YX')*sqrt5`, and
each theorem below is one integral matrix identity checked by kernel
computation (`decide`).

Proved as theorems in this file:

* `partition`: the dot-product partition `I + A + A' + Pi = J`;
* `projector_zero`: `x0*x0 = 60*x0`; `projector_one_rat` and
  `projector_one_irr`: the two components of `(5*B_1)^2 = 20*(5*B_1)`;
  `projector_two`: `x2*x2 = 12*x2`; `projector_three_rat`: a
  restatement of `projector_one_rat` (the `B_3` reading flips the sign
  of the `sqrt5` component);
* `orthogonal_zero_one_rat`, `orthogonal_zero_one_irr`,
  `orthogonal_zero_two`, `orthogonal_one_two_rat`,
  `orthogonal_one_two_irr`, `orthogonal_one_three_rat`,
  `orthogonal_one_three_irr`: all components of `B_0 B_1 = 0`,
  `B_0 B_2 = 0`, `B_1 B_2 = 0`, `B_1 B_3 = 0`;
* `traces`: `tr x0 = tr x1 = tr x2 = 60` and `tr y1 = 0`;
* `resolution_of_identity`: `5*x0 + 30*x1 + 25*x2 = 300*I`, the
  integral form of `B_0/12 + B_1/4 + 5*B_2/12 + B_3/4 = I` (the
  `sqrt5` parts cancel by the definitional reading of `B_3`);
* `binding_zero`, `binding_one_rat`, `binding_one_irr`, `binding_two`:
  the components of `A B_0 = 5 B_0`, `A B_1 = sqrt5 B_1`,
  `A B_2 = -B_2`; `binding_three_rat`: a sign-flip restatement of
  `binding_one_rat` under the `B_3` reading;
* `parity_even`, `parity_odd`: `Pi x0 = x0`, `Pi x2 = x2`,
  `Pi x1 = -x1`, `Pi y1 = -y1`, the component form of
  `Pi B_l = (-1)^l B_l` (the level-3 reading reuses the level-1
  identities with the flipped `sqrt5` sign).

Not proved in this file. Each companion claim lives in the Python
certificate `code/angular_sprint/angular_interpolant_certificate.py`
unless a different home is named:

* the sphere embedding of the twelve ports, the dot-product classes
  `{1, 1/sqrt5, -1/sqrt5, -1}`, and the identification of the reading
  above with the Legendre Gram kernels `B_l(i, j) = P_l(v_i . v_j)`
  (`build_kernels`, exact `Q(sqrt5)` arithmetic);
* the rank values `(1, 3, 5, 3)`: only the traces are proved here; the
  trace-over-scale rank readout is computed in `projector_certificate`;
* the remaining `B_3` table entries: the `sqrt5` component of the
  `B_3` projector square, `B_0 B_3 = 0`, `B_2 B_3 = 0`, and the
  `sqrt5` component of `A B_3 = -sqrt5 B_3` are sign-flip and
  transpose corollaries of the identities above, not stated as
  theorems here; the full four-by-four tables are checked in exact
  arithmetic in the Python certificate;
* A5 equivariance of the kernels and the order-120 generator closure
  (`equivariance_certificate`);
* the readback response `J = (A^3 - 4A^2 - 5A + 10I)/10`, `R = -J`,
  and the alternating readback parity `R B_l = (-1)^(l+1) B_l`
  (`parity_response_certificate`); the parity theorems here state the
  `Pi` sign table only;
* refinement naturality: recorded as a premise-typed conditional in
  `band_binding_certificate`, with the premise discharged in no
  registered producer;
* the equal-port sequence `I_l` and its even initial vector: the
  companion Lean module `Lean/Screen/A5AngularBands.lean`
  (`OPH.AngularBands`).
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

/-- Traces: `tr x0 = tr x1 = tr x2 = 60` and `tr y1 = 0`. The rank
readout `(1, 3, 5, 3)` from trace over scale is not stated here; it is
computed in the Python certificate (`projector_certificate`). -/
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
/-- Antipodal parity, even levels: `Pi x0 = x0` and `Pi x2 = x2`. The
readback parity through `R = -J` is a Python-side claim
(`parity_response_certificate`); this file states the `Pi` signs only. -/
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
