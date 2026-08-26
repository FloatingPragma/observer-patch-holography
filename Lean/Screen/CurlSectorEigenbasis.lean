import CarrierEvolutionFlow
import FieldSectorEnergyInnerProduct
import GoldenSectorCharacters

set_option autoImplicit false

open scoped BigOperators Matrix

namespace OPH.CurlSectorEigenbasis

open OPH.DiscreteCoulombGreen
open OPH.PositionSpaceMaxwellAction
open OPH.LocalFaceMaxwellAction
open OPH.TemporalMaxwellEvolution
open OPH.ScaledMaxwellStability
open OPH.CarrierModeOscillators
open OPH.CarrierEvolutionFlow
open OPH.FieldSectorEnergyInnerProduct
open OPH.GoldenSectorCharacters
open OPH.ScreenCarrierMapCandidate
open OPH.SeamCurrentCarrierQuotient

/-!
# An explicit orthogonal eigenbasis of the curl sector of the committed seam
operator, and the carrier evolution flow on the whole seam space (issue #733)

STATUS.  Kernel `decide` checks on committed finite tables (the twenty
oriented faces and thirty seams of `LocalFaceMaxwellAction`, the face
incidence `faceIncidenceZ`, the seam table `seamLeft`/`seamRight`) and on
the new `ℤ[φ]` tables listed below, transported to exact real linear
algebra through the committed `evalPhi : Zphi → ℝ`, `(a, b) ↦ a + bφ`.
The carrier, its face orientation, the operator `CᵀC =
localMaxwellOperator`, the zero-current scaled Ampere evolution
`AmpereEvolutionScaled h A 0 0` in the temporal gauge, and the step `h`
are declared in the imported files.  No register row is discharged.  No
`native_decide` is used.

WHAT IS PROVED.
1. Spectrum and the nineteen-vector selection.  `curlZ : Fin 19 → Fin 30 →
   Zphi` lists nineteen seam vectors with entries in `ℤ[φ]`, grouped by
   eigenvalue: rows `0..4` for `2`, `5..8` for `3`, `9..12` for `5`,
   `13..15` for `λ₊ = 2 + 2φ = 3 + √5`, `16..18` for `λ₋ = 4 - 2φ = 3 - √5`
   (`curlLamZ`).  Each is transported by the committed incidence to the
   listed face vector `curlFaceZ` and back to the eigenvalue times itself
   (`curl_curvature_Z`, `curl_codifferential_Z`), so each real cast `curlR i`
   is an eigenvector of `CᵀC` with eigenvalue `curlLamR i` (`curl_eigen`).
   The `ℤ[φ]` Gram matrix of the nineteen is diagonal with the listed
   nonzero norms `curlNormZ` (`curl_gram_Z`); hence the real vectors are
   pairwise orthogonal (`curl_orth`) and nonzero (`curl_ne_zero`), so
   linearly independent (`curl_linearIndependent`).  Each lies in the
   range of the face codifferential, exhibited as `Cᵀ` of `λ⁻¹` times its
   listed face vector (`curl_mem_range_codifferential`), and is orthogonal
   to every gradient seam field (`curl_orth_gradient`).  The span of the
   nineteen is exactly the range of `Cᵀ`, which equals the range of `CᵀC`
   and the kernel of the committed boundary (`curl_span`,
   `curl_sector_characterizations`), all of dimension `19`.  Row `0` is
   the committed `twoMode`, row `5` is `threeMode`, row `9` is `fiveMode`,
   and row `13` is `goldenMode = √5 • goldenA + goldenB` written as
   `(B - A) + 2A φ` (`curl_committed_members`; the five families of
   `carrier_flow_five_modes` are placed in their sectors by
   `fiveFamily_mem_sector`); rows `16..18` are the
   entrywise Galois conjugates `a + bφ ↦ (a + b) - bφ` of rows `13..15`
   (`curl_minus_eq_conj_plus`).  The multiplicities `5, 4, 4, 3, 3` sum
   to `19 = 30 - 11`; the exact characteristic polynomial
   `x¹¹ (x-2)⁵ (x-3)⁴ (x-5)⁴ (x² - 6x + 4)³` is an inference outside Lean
   (verified in the lane script, exact in `ℚ(√5)`), and inside Lean the
   count is carried by the dimension identities above.
2. Gradient sector.  `gradPortZ : Fin 11 → Fin 12 → Zphi` lists eleven
   port loads (eigenvectors of the committed vertex Laplacian, an
   observation not used); their coboundaries `gradR i = realCoboundary
   (gradPortR i)` have diagonal `ℤ[φ]` Gram matrix with nonzero norms
   (`grad_gram_Z`), so they are pairwise orthogonal, nonzero, linearly
   independent, and span exactly the range of the coboundary
   (`grad_span`), the kernel of `CᵀC`.
3. Exact complement.  The range of `Cᵀ` and the range of the coboundary
   are complementary: disjoint and with sum the whole seam space
   (`curl_grad_isCompl`, `curl_sup_grad_eq_top`), by orthogonality and the
   dimension count `19 + 11 = 30`.
4. The thirty-vector family.  `fullR = Sum.elim curlR gradR` on
   `Fin 19 ⊕ Fin 11` with eigenvalues `fullLam = Sum.elim curlLamR 0` is a
   family of pairwise orthogonal nonzero eigenvectors (`full_eigen`,
   `full_orth`, `full_ne_zero`) spanning `Fin 30 → ℝ` (`full_span_top`);
   every pair of seam fields `(A₀, E₀)` is the potential and electric
   readout of some coefficient state (`exists_state_of_fields`).
5. Flow on the whole seam space.  For a declared step `h ≠ 0` inside the
   window `h² (3 + √5) < 4`, every eigenvalue of the family is admissible
   in the sense of `CarrierEvolutionFlow.Admissible`
   (`fullLam_admissible`), and conversely admissibility of every eigenvalue
   of the family forces `h² (3 + √5) < 4` (`fullLam_admissible_iff`, through
   the row `13` eigenvalue `3 + √5`, the largest listed,
   `curlLamR_le_golden`); so `h² (3 + √5) < 4` is the admissibility window
   of the family in the committed sense of `Admissible`, no further
   boundedness statement being made here.  `carrier_flow_full_curl` assembles
   the committed `assembled_flow` packet on the thirty-vector family
   together with: the span is the whole seam space, every initial pair
   `(A₀, E₀)` is realised by a coefficient state, and every zero-current
   temporal-gauge solution of the committed evolution is the assembled
   history of some coefficient state (`ampere_history_eq_assembled`,
   through the two-step form of the evolution `ampere_step_eq` and its
   uniqueness `ampere_unique`).  So the `h`-step of the assembled flow is
   the committed scaled Ampere step on every zero-current history, not
   only on the five listed families of `carrier_flow_five_modes`.  On the
   nineteen nonzero curl modes, the energy-derived positive definite real
   form and Hermitian form of `FieldSectorEnergyInnerProduct` are preserved,
   the complex coordinates have the displayed phase and generator, the
   field energy is their diagonal, and the coefficients are recovered from
   the potential and electric readouts (`curl_sector_energy_hilbert_packet`).

WHAT IS NOT PROVED.  No physical attachment: the step `h` and the flow
parameter `t` are declared, the seam vectors are not identified with
physical modes, no propagation speed, mass, or laboratory frequency is
attached, and `t` is not identified with physical time (source clock and
duration row).  The identification of the unique continuous flow of the
OL-C2 row surface (`EventAlgebra/QuantumAdequacySurface.lean`,
`unique_continuous_flow`) with this coefficient flow is not stated: OL-C2
is cited, not discharged.  Multiplicities are read from dimension
identities of subspaces, not from a characteristic polynomial computed in
Lean.  The family `Fin 19 ⊕ Fin 11` is an orthogonal basis, not
orthonormal; the flow is stated on coefficient states and transported to
field histories through `potentialOf`/`electricOf` as in the imported
files.  Its finite coefficient-space Hermitian isometry is proved, but it is
not identified with a physical photon Hilbert space or with the
private-algebra Stone flow.

PRIOR WORK.  `ScaledMaxwellStability` (`fiveMode`, `goldenMode`, the face
projectors, `localMaxwellOperator_cast`); `CarrierModeOscillators`
(`twoMode`, `threeMode`, `codifferential_eigen`); `GoldenSectorCharacters`
(`Zphi` projector tables, `zint`, `evalPhi_sum`, `evalPhi_zint`,
`lamPlus`, `lamMinus`, `evalPhi_lamPlus`, `evalPhi_lamMinus`);
`Geometry/ScreenCarrierMapCandidate` (`Zphi`, `zmul`, `zsub`, `evalPhi`,
`evalPhi_ne_zero`, the irrationality-of-`φ` step used for nonvanishing);
`CarrierEvolutionFlow` (`Admissible`, `assembledFlow`, `assembled_flow`,
`eigen_orthogonal`, `gradient_eigen`, `carrier_flow_five_modes`, whose
docstring names the nineteen-vector selection as the missing input);
`LocalFaceMaxwellAction` (`faceCurvature_coboundary`,
`faceCurvature_codifferential_adjoint`, `range_localMaxwellOperator`);
`PositionSpaceMaxwellAction` (`gauge_finrank`, `cycleSpace_finrank`).

ROWS TOUCHED (none discharged).  OL-C2 (unique continuous flow on the
private algebra): cited as the target reading; this module supplies the
finite field-sector flow on the whole seam space and joins nothing to the
private algebra.  Source clock and duration row: `h`, `t` declared, no
unit.  Physical spacetime attachment row: the carrier is the committed
combinatorial complex.  Coupled-action row: the kinetic term is declared
in `ScaledMaxwellStability`.

CONVENTIONS.  `Zphi = ℤ × ℤ`, `(a, b) = a + bφ`, `φ² = φ + 1`,
`√5 = 2φ - 1`; products `zmul`, differences `zsub`, integer scaling
`zint`; `evalPhi (a, b) = a + bφ`.  Forward differences in the temporal
gauge as in the imported files: `E n = -(h⁻¹ • (A (n+1) - A n))`,
`B n = C (A n)`; the two-step form of the zero-current evolution is
`A (n+2) = 2 A (n+1) - A n - h² CᵀC A (n+1)`.

FALSIFIER.  A face-vector entry off the listed `curlFaceZ`, a
codifferential entry off `λ` times the listed seam entry, an off-diagonal
Gram entry that is nonzero, a listed norm that is zero, or a solution of
the committed evolution not equal to the assembled history of its initial
state would make the corresponding theorem fail.

Axiom audit.  The `#print axioms` lines at the end of the file show at
most `propext`, `Classical.choice`, and `Quot.sound`; no `native_decide`.
-/

/-! ## 1. The tables -/

/-- The nineteen curl vectors, entries `(a, b) = a + bφ`.  Rows `0..4`:
eigenvalue `2`; `5..8`: `3`; `9..12`: `5`; `13..15`: `2 + 2φ = 3 + √5`;
`16..18`: `4 - 2φ = 3 - √5`.  Row `0` is `twoModeZ`, row `5` is
`threeModeZ`, row `9` is `fiveModeZ`, row `13` is `goldenMode` as
`(B - A) + 2A φ`; rows `16..18` are the Galois conjugates of rows `13..15`.
Within each eigenvalue the rows were produced by exact Gram--Schmidt over
`ℚ(φ)` seeded with the committed mode, then scaled into `ℤ[φ]`; the
receipts below are the only facts used. -/
def curlZ : Fin 19 → Fin 30 → Zphi :=
  ![
    ![(1, 0), (-1, 0), (1, 0), (-1, 0), (0, 0), (1, 0), (-1, 0), (1, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (-1, 0), (1, 0), (1, 0), (-1, 0), (-1, 0), (1, 0), (1, 0), (-1, 0), (1, 0)],
    ![(0, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (-1, 0), (1, 0), (1, 0), (-1, 0), (0, 0), (-1, 0), (1, 0), (0, 0), (0, 0), (1, 0), (1, 0), (0, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (1, 0), (1, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0)],
    ![(0, 0), (0, 0), (2, 0), (-1, 0), (-1, 0), (0, 0), (2, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (2, 0), (1, 0), (1, 0), (2, 0), (-2, 0), (1, 0), (-1, 0), (-2, 0), (1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (2, 0), (2, 0), (0, 0), (0, 0), (0, 0)],
    ![(-2, 0), (2, 0), (1, 0), (-1, 0), (0, 0), (4, 0), (-1, 0), (-2, 0), (-3, 0), (1, 0), (2, 0), (3, 0), (-3, 0), (0, 0), (3, 0), (3, 0), (0, 0), (-3, 0), (3, 0), (-3, 0), (0, 0), (2, 0), (-2, 0), (1, 0), (-1, 0), (-1, 0), (1, 0), (4, 0), (2, 0), (-2, 0)],
    ![(2, 0), (2, 0), (-1, 0), (-1, 0), (-2, 0), (0, 0), (1, 0), (0, 0), (1, 0), (1, 0), (0, 0), (1, 0), (1, 0), (-2, 0), (1, 0), (1, 0), (-2, 0), (1, 0), (1, 0), (1, 0), (-2, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (1, 0), (-1, 0), (0, 0), (2, 0), (2, 0)],
    ![(2, 0), (-2, 0), (1, 0), (-1, 0), (0, 0), (2, 0), (-1, 0), (1, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (2, 0), (-2, 0), (0, 0), (-2, 0), (2, 0), (0, 0), (2, 0), (-2, 0), (0, 0), (1, 0), (-1, 0), (-1, 0), (1, 0), (1, 0), (-1, 0), (-2, 0), (2, 0), (-2, 0)],
    ![(0, 0), (0, 0), (-1, 0), (-1, 0), (2, 0), (0, 0), (-1, 0), (-1, 0), (2, 0), (-1, 0), (-1, 0), (2, 0), (0, 0), (0, 0), (-2, 0), (0, 0), (0, 0), (-2, 0), (0, 0), (0, 0), (-2, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0), (0, 0), (0, 0), (0, 0)],
    ![(1, 0), (-1, 0), (-2, 0), (2, 0), (0, 0), (-2, 0), (0, 0), (2, 0), (1, 0), (0, 0), (-2, 0), (-1, 0), (-1, 0), (-2, 0), (1, 0), (1, 0), (2, 0), (-1, 0), (-1, 0), (1, 0), (0, 0), (2, 0), (-2, 0), (0, 0), (-2, 0), (0, 0), (2, 0), (2, 0), (1, 0), (-1, 0)],
    ![(-3, 0), (-3, 0), (2, 0), (2, 0), (2, 0), (0, 0), (-4, 0), (2, 0), (-1, 0), (-4, 0), (2, 0), (-1, 0), (-3, 0), (0, 0), (1, 0), (-3, 0), (0, 0), (1, 0), (3, 0), (3, 0), (-2, 0), (-2, 0), (-2, 0), (4, 0), (-2, 0), (4, 0), (-2, 0), (0, 0), (3, 0), (3, 0)],
    ![(0, 0), (0, 0), (0, 0), (-1, 0), (1, 0), (1, 0), (-1, 0), (-2, 0), (2, 0), (0, 0), (1, 0), (0, 0), (0, 0), (-1, 0), (0, 0), (-2, 0), (-1, 0), (2, 0), (-2, 0), (0, 0), (1, 0), (1, 0), (-2, 0), (0, 0), (-1, 0), (-1, 0), (0, 0), (1, 0), (0, 0), (0, 0)],
    ![(0, 0), (0, 0), (6, 0), (1, 0), (-7, 0), (5, 0), (1, 0), (-4, 0), (-2, 0), (6, 0), (11, 0), (-12, 0), (12, 0), (7, 0), (-12, 0), (2, 0), (7, 0), (-2, 0), (2, 0), (12, 0), (-7, 0), (11, 0), (-4, 0), (6, 0), (1, 0), (1, 0), (6, 0), (5, 0), (0, 0), (0, 0)],
    ![(0, 0), (7, 0), (-5, 0), (-9, 0), (7, 0), (-3, 0), (5, 0), (1, 0), (-3, 0), (9, 0), (-1, 0), (-4, 0), (-3, 0), (7, 0), (-4, 0), (-4, 0), (7, 0), (-3, 0), (-4, 0), (-3, 0), (7, 0), (-1, 0), (1, 0), (9, 0), (-9, 0), (5, 0), (-5, 0), (-3, 0), (7, 0), (0, 0)],
    ![(4, 0), (-3, 0), (-3, 0), (1, 0), (1, 0), (3, 0), (3, 0), (-1, 0), (-1, 0), (-1, 0), (1, 0), (0, 0), (-1, 0), (1, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (-1, 0), (1, 0), (1, 0), (-1, 0), (-1, 0), (1, 0), (3, 0), (-3, 0), (3, 0), (-3, 0), (4, 0)],
    ![(0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (2, 0), (-2, 0), (0, -2), (0, 2), (2, 0), (0, 2), (0, -2), (-2, 0), (0, -2), (0, 2), (2, 0), (0, 2), (0, -2), (-2, 0), (2, 0), (0, 0), (0, -2), (0, 2), (-2, 0), (0, 0), (2, 0), (0, 0), (-2, 0), (0, 0), (0, 0)],
    ![(0, 0), (-5, 0), (5, 0), (0, 5), (0, -5), (-1, 2), (1, -2), (-2, -1), (2, 1), (-1, -3), (-3, 1), (-2, 4), (1, 3), (3, -1), (2, -4), (4, 2), (-3, 1), (-2, -1), (-4, -2), (-1, -3), (0, 5), (3, -1), (2, 1), (1, 3), (0, -5), (-1, 2), (-5, 0), (1, -2), (5, 0), (0, 0)],
    ![(-2, 0), (0, 1), (0, 1), (1, -1), (1, -1), (0, -1), (0, -1), (-1, 1), (-1, 1), (1, 0), (-1, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (1, -1), (0, 0), (-1, 0), (-1, 1), (1, 0), (1, -1), (-1, 0), (-1, 1), (0, 1), (0, -1), (0, 1), (0, -1), (2, 0)],
    ![(0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (2, 0), (-2, 0), (-2, 2), (2, -2), (2, 0), (2, -2), (-2, 2), (-2, 0), (-2, 2), (2, -2), (2, 0), (2, -2), (-2, 2), (-2, 0), (2, 0), (0, 0), (-2, 2), (2, -2), (-2, 0), (0, 0), (2, 0), (0, 0), (-2, 0), (0, 0), (0, 0)],
    ![(0, 0), (-5, 0), (5, 0), (5, -5), (-5, 5), (1, -2), (-1, 2), (-3, 1), (3, -1), (-4, 3), (-2, -1), (2, -4), (4, -3), (2, 1), (-2, 4), (6, -2), (-2, -1), (-3, 1), (-6, 2), (-4, 3), (5, -5), (2, 1), (3, -1), (4, -3), (-5, 5), (1, -2), (-5, 0), (-1, 2), (5, 0), (0, 0)],
    ![(-2, 0), (1, -1), (1, -1), (0, 1), (0, 1), (-1, 1), (-1, 1), (0, -1), (0, -1), (1, 0), (-1, 0), (0, 0), (1, 0), (-1, 0), (0, 0), (0, 0), (1, 0), (0, 1), (0, 0), (-1, 0), (0, -1), (1, 0), (0, 1), (-1, 0), (0, -1), (1, -1), (-1, 1), (1, -1), (-1, 1), (2, 0)]]

/-- The face vectors `C (curlZ i)`. -/
def curlFaceZ : Fin 19 → Fin 20 → Zphi :=
  ![
    ![(3, 0), (1, 0), (1, 0), (-1, 0), (-1, 0), (1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (1, 0), (1, 0), (1, 0), (3, 0)],
    ![(0, 0), (0, 0), (0, 0), (0, 0), (2, 0), (0, 0), (0, 0), (2, 0), (-2, 0), (-2, 0), (-2, 0), (-2, 0), (2, 0), (0, 0), (2, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0)],
    ![(0, 0), (0, 0), (0, 0), (-4, 0), (-2, 0), (0, 0), (4, 0), (2, 0), (2, 0), (-2, 0), (-2, 0), (2, 0), (2, 0), (4, 0), (-2, 0), (-4, 0), (0, 0), (0, 0), (0, 0), (0, 0)],
    ![(0, 0), (4, 0), (4, 0), (2, 0), (2, 0), (-8, 0), (2, 0), (-4, 0), (2, 0), (-4, 0), (-4, 0), (2, 0), (-4, 0), (2, 0), (2, 0), (2, 0), (-8, 0), (4, 0), (4, 0), (0, 0)],
    ![(0, 0), (-4, 0), (4, 0), (-2, 0), (2, 0), (0, 0), (-2, 0), (0, 0), (2, 0), (0, 0), (0, 0), (2, 0), (0, 0), (-2, 0), (2, 0), (-2, 0), (0, 0), (4, 0), (-4, 0), (0, 0)],
    ![(6, 0), (0, 0), (0, 0), (-3, 0), (-3, 0), (0, 0), (-3, 0), (-3, 0), (-3, 0), (-3, 0), (3, 0), (3, 0), (3, 0), (3, 0), (3, 0), (3, 0), (0, 0), (0, 0), (0, 0), (-6, 0)],
    ![(0, 0), (0, 0), (0, 0), (3, 0), (-3, 0), (0, 0), (-3, 0), (3, 0), (3, 0), (-3, 0), (3, 0), (-3, 0), (-3, 0), (3, 0), (3, 0), (-3, 0), (0, 0), (0, 0), (0, 0), (0, 0)],
    ![(0, 0), (-3, 0), (-3, 0), (3, 0), (3, 0), (6, 0), (-3, 0), (0, 0), (-3, 0), (0, 0), (0, 0), (3, 0), (0, 0), (3, 0), (-3, 0), (-3, 0), (-6, 0), (3, 0), (3, 0), (0, 0)],
    ![(0, 0), (9, 0), (-9, 0), (3, 0), (-3, 0), (0, 0), (-3, 0), (-6, 0), (3, 0), (6, 0), (-6, 0), (-3, 0), (6, 0), (3, 0), (3, 0), (-3, 0), (0, 0), (9, 0), (-9, 0), (0, 0)],
    ![(1, 0), (1, 0), (1, 0), (1, 0), (-4, 0), (-4, 0), (-4, 0), (6, 0), (1, 0), (1, 0), (1, 0), (1, 0), (6, 0), (-4, 0), (-4, 0), (1, 0), (-4, 0), (1, 0), (1, 0), (1, 0)],
    ![(5, 0), (5, 0), (5, 0), (-25, 0), (10, 0), (-20, 0), (10, 0), (0, 0), (-25, 0), (35, 0), (35, 0), (-25, 0), (0, 0), (10, 0), (10, 0), (-25, 0), (-20, 0), (5, 0), (5, 0), (5, 0)],
    ![(-10, 0), (-10, 0), (25, 0), (15, 0), (-20, 0), (5, 0), (15, 0), (0, 0), (-20, 0), (0, 0), (0, 0), (-20, 0), (0, 0), (15, 0), (-20, 0), (15, 0), (5, 0), (25, 0), (-10, 0), (-10, 0)],
    ![(10, 0), (-10, 0), (-5, 0), (5, 0), (0, 0), (-5, 0), (5, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (5, 0), (0, 0), (5, 0), (-5, 0), (-5, 0), (-10, 0), (10, 0)],
    ![(2, 0), (2, 0), (2, 0), (2, 0), (2, 0), (-2, -4), (-2, -4), (2, 4), (-2, -4), (2, 4), (-2, -4), (2, 4), (-2, -4), (2, 4), (-2, 0), (-2, 0), (2, 4), (-2, 0), (-2, 0), (-2, 0)],
    ![(4, 2), (4, 2), (-6, -8), (-6, -8), (4, 12), (2, -4), (2, -4), (8, 4), (2, 6), (-2, -6), (2, 6), (-2, -6), (-8, -4), (-2, 4), (-4, -12), (6, 8), (-2, 4), (6, 8), (-4, -2), (-4, -2)],
    ![(-2, -2), (2, 2), (0, 2), (0, -2), (0, 0), (0, 2), (0, -2), (0, 0), (-2, 0), (-2, 0), (2, 0), (2, 0), (0, 0), (0, 2), (0, 0), (0, 2), (0, -2), (0, -2), (-2, -2), (2, 2)],
    ![(2, 0), (2, 0), (2, 0), (2, 0), (2, 0), (-6, 4), (-6, 4), (6, -4), (-6, 4), (6, -4), (-6, 4), (6, -4), (-6, 4), (6, -4), (-2, 0), (-2, 0), (6, -4), (-2, 0), (-2, 0), (-2, 0)],
    ![(6, -2), (6, -2), (-14, 8), (-14, 8), (16, -12), (-2, 4), (-2, 4), (12, -4), (8, -6), (-8, 6), (8, -6), (-8, 6), (-12, 4), (2, -4), (-16, 12), (14, -8), (2, -4), (14, -8), (-6, 2), (-6, 2)],
    ![(-4, 2), (4, -2), (2, -2), (-2, 2), (0, 0), (2, -2), (-2, 2), (0, 0), (-2, 0), (-2, 0), (2, 0), (2, 0), (0, 0), (2, -2), (0, 0), (2, -2), (-2, 2), (-2, 2), (-4, 2), (4, -2)]]

/-- The `ℤ[φ]` eigenvalues in row order. -/
def curlLamZ : Fin 19 → Zphi :=
  ![(2, 0), (2, 0), (2, 0), (2, 0), (2, 0), (3, 0), (3, 0), (3, 0), (3, 0), (5, 0), (5, 0), (5, 0), (5, 0), lamPlus, lamPlus, lamPlus, lamMinus, lamMinus, lamMinus]

/-- The `ℤ[φ]` norms `∑ e, (curlZ i e)²`. -/
def curlNormZ : Fin 19 → Zphi :=
  ![(18, 0), (16, 0), (48, 0), (144, 0), (48, 0), (60, 0), (36, 0), (60, 0), (180, 0), (36, 0), (1260, 0), (840, 0), (120, 0), (80, 40), (400, 200), (40, 0), (120, -40), (600, -200), (40, 0)]

/-- Eleven port loads whose coboundaries form an orthogonal basis of the
gradient sector (they are eigenvectors of the committed vertex Laplacian
for `5 - √5`, `6`, `5 + √5`; that reading is not used). -/
def gradPortZ : Fin 11 → Fin 12 → Zphi :=
  ![
    ![(0, 0), (0, 0), (-1, 0), (1, 0), (1, -1), (1, -1), (-1, 1), (-1, 1), (-1, 0), (1, 0), (0, 0), (0, 0)],
    ![(0, 0), (-2, 0), (1, -1), (1, -1), (0, 1), (0, -1), (0, 1), (0, -1), (-1, 1), (-1, 1), (2, 0), (0, 0)],
    ![(-5, 0), (1, -2), (1, -2), (1, -2), (1, -2), (-1, 2), (1, -2), (-1, 2), (-1, 2), (-1, 2), (-1, 2), (5, 0)],
    ![(0, 0), (0, 0), (0, 0), (0, 0), (1, 0), (-1, 0), (-1, 0), (1, 0), (0, 0), (0, 0), (0, 0), (0, 0)],
    ![(0, 0), (0, 0), (0, 0), (2, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (2, 0), (0, 0), (0, 0), (0, 0)],
    ![(0, 0), (0, 0), (3, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (3, 0), (0, 0), (0, 0)],
    ![(0, 0), (4, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (4, 0), (0, 0)],
    ![(5, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (-1, 0), (5, 0)],
    ![(0, 0), (0, 0), (-1, 0), (1, 0), (0, 1), (0, 1), (0, -1), (0, -1), (-1, 0), (1, 0), (0, 0), (0, 0)],
    ![(0, 0), (-2, 0), (0, 1), (0, 1), (1, -1), (-1, 1), (1, -1), (-1, 1), (0, -1), (0, -1), (2, 0), (0, 0)],
    ![(-5, 0), (-1, 2), (-1, 2), (-1, 2), (-1, 2), (1, -2), (-1, 2), (1, -2), (1, -2), (1, -2), (1, -2), (5, 0)]]

/-- The `ℤ[φ]` norms `∑ e, (gradZ i e)²`. -/
def gradNormZ : Fin 11 → Zphi :=
  ![(80, -40), (120, -40), (600, -200), (24, 0), (72, 0), (144, 0), (240, 0), (360, 0), (40, 40), (80, 40), (400, 200)]

/-- Seam vector of the `ℤ[φ]` gradient family: the coboundary of the
listed port load, `gradZ i e = gradPortZ i (seamRight e) - gradPortZ i (seamLeft e)`. -/
def gradZ (i : Fin 11) (e : Fin 30) : Zphi :=
  zsub (gradPortZ i (seamRight e)) (gradPortZ i (seamLeft e))

/-! ## 2. Kernel checks on the tables -/

set_option maxRecDepth 32768 in
/-- **Curvature receipt.**  `C` carries each listed seam vector to its listed
face vector. -/
theorem curl_curvature_Z :
    ∀ (i : Fin 19) (f : Fin 20),
      (∑ e : Fin 30, zint (faceIncidenceZ f e) (curlZ i e)) = curlFaceZ i f := by
  decide +kernel

set_option maxRecDepth 32768 in
/-- **Codifferential receipt.**  `Cᵀ` carries each listed face vector to the
listed eigenvalue times the listed seam vector. -/
theorem curl_codifferential_Z :
    ∀ (i : Fin 19) (e : Fin 30),
      (∑ f : Fin 20, zint (faceIncidenceZ f e) (curlFaceZ i f)) =
        zmul (curlLamZ i) (curlZ i e) := by
  decide +kernel

set_option maxRecDepth 32768 in
/-- **Diagonal Gram matrix** of the nineteen curl vectors over `ℤ[φ]`. -/
theorem curl_gram_Z :
    ∀ i j : Fin 19,
      (∑ e : Fin 30, zmul (curlZ i e) (curlZ j e)) =
        if i = j then curlNormZ i else 0 := by
  decide +kernel

set_option maxRecDepth 32768 in
/-- **Diagonal Gram matrix** of the eleven gradient vectors over `ℤ[φ]`. -/
theorem grad_gram_Z :
    ∀ i j : Fin 11,
      (∑ e : Fin 30, zmul (gradZ i e) (gradZ j e)) =
        if i = j then gradNormZ i else 0 := by
  decide +kernel

theorem curlNormZ_ne_zero : ∀ i : Fin 19, curlNormZ i ≠ 0 := by decide
theorem gradNormZ_ne_zero : ∀ i : Fin 11, gradNormZ i ≠ 0 := by decide

/-- Rows `16..18` are the entrywise Galois conjugates of rows `13..15`. -/
theorem curl_minus_eq_conj_plus :
    ∀ e : Fin 30, curlZ 16 e = zconj (curlZ 13 e) ∧ curlZ 17 e = zconj (curlZ 14 e) ∧
      curlZ 18 e = zconj (curlZ 15 e) := by
  decide

/-- Row `0` is `twoModeZ`, row `5` is `threeModeZ`, row `9` is `fiveModeZ`,
row `13` is `(goldenBZ - goldenAZ) + 2 goldenAZ φ`. -/
theorem curl_committed_rows_Z :
    ∀ e : Fin 30, curlZ 0 e = (twoModeZ e, 0) ∧ curlZ 5 e = (threeModeZ e, 0) ∧
      curlZ 9 e = (fiveModeZ e, 0) ∧
      curlZ 13 e = (goldenBZ e - goldenAZ e, 2 * goldenAZ e) := by
  decide

/-! ## 3. Real transport through `evalPhi` -/

noncomputable section

theorem faceCurvature_evalPhi (vZ : Fin 30 → Zphi) (wZ : Fin 20 → Zphi)
    (h : ∀ f : Fin 20, (∑ e : Fin 30, zint (faceIncidenceZ f e) (vZ e)) = wZ f) :
    faceCurvature (fun e ↦ evalPhi (vZ e)) = fun f ↦ evalPhi (wZ f) := by
  funext f
  rw [faceCurvature_apply, ← h f, evalPhi_sum]
  refine Finset.sum_congr rfl fun e _ ↦ ?_
  rw [evalPhi_zint]
  rfl

theorem faceCodifferential_evalPhi (wZ : Fin 20 → Zphi) (uZ : Fin 30 → Zphi)
    (h : ∀ e : Fin 30, (∑ f : Fin 20, zint (faceIncidenceZ f e) (wZ f)) = uZ e) :
    faceCodifferential (fun f ↦ evalPhi (wZ f)) = fun e ↦ evalPhi (uZ e) := by
  funext e
  rw [faceCodifferential_apply, ← h e, evalPhi_sum]
  refine Finset.sum_congr rfl fun f _ ↦ ?_
  rw [evalPhi_zint]
  rfl

/-- Transport of an eigen receipt over `ℤ[φ]` to the real operator. -/
theorem localMaxwellOperator_evalPhi (vZ : Fin 30 → Zphi) (wZ : Fin 20 → Zphi)
    (lam : Zphi)
    (h1 : ∀ f : Fin 20, (∑ e : Fin 30, zint (faceIncidenceZ f e) (vZ e)) = wZ f)
    (h2 : ∀ e : Fin 30, (∑ f : Fin 20, zint (faceIncidenceZ f e) (wZ f)) = zmul lam (vZ e)) :
    localMaxwellOperator (fun e ↦ evalPhi (vZ e)) =
      evalPhi lam • fun e ↦ evalPhi (vZ e) := by
  show faceCodifferential (faceCurvature _) = _
  rw [faceCurvature_evalPhi vZ wZ h1, faceCodifferential_evalPhi wZ _ h2]
  funext e
  simp only [Pi.smul_apply, smul_eq_mul, evalPhi_zmul]

/-- The seam pairing of two `ℤ[φ]` vectors is the evaluation of their `ℤ[φ]`
pairing. -/
theorem realSeamInner_evalPhi (uZ vZ : Fin 30 → Zphi) :
    realSeamInner (fun e ↦ evalPhi (uZ e)) (fun e ↦ evalPhi (vZ e)) =
      evalPhi (∑ e : Fin 30, zmul (uZ e) (vZ e)) := by
  unfold realSeamInner
  rw [evalPhi_sum]
  refine Finset.sum_congr rfl fun e _ ↦ ?_
  rw [evalPhi_zmul]

/-! ## 4. The nineteen real curl vectors -/

/-- The real cast of the `i`-th curl vector. -/
def curlR (i : Fin 19) : Fin 30 → ℝ := fun e ↦ evalPhi (curlZ i e)

/-- The real eigenvalues `2, 3, 5, 3 + √5, 3 - √5` in row order. -/
def curlLamR : Fin 19 → ℝ :=
  ![2, 2, 2, 2, 2, 3, 3, 3, 3, 5, 5, 5, 5,
    3 + Real.sqrt 5, 3 + Real.sqrt 5, 3 + Real.sqrt 5,
    3 - Real.sqrt 5, 3 - Real.sqrt 5, 3 - Real.sqrt 5]

theorem evalPhi_curlLam (i : Fin 19) : evalPhi (curlLamZ i) = curlLamR i := by
  fin_cases i <;> simp [curlLamZ, curlLamR, lamPlus, lamMinus, evalPhi, Real.goldenRatio] <;>
    ring

theorem curlLamR_pos (i : Fin 19) : 0 < curlLamR i := by
  have h3 := OPH.CarrierModeOscillators.sqrt5_lt_three
  have h0 := sqrt5_pos
  fin_cases i <;> simp [curlLamR] <;> linarith

theorem curlLamR_ne_zero (i : Fin 19) : curlLamR i ≠ 0 := (curlLamR_pos i).ne'

/-- The largest listed eigenvalue is `3 + √5`. -/
theorem curlLamR_le_golden (i : Fin 19) : curlLamR i ≤ 3 + Real.sqrt 5 := by
  have h2 := two_lt_sqrt5
  have h0 := sqrt5_pos
  fin_cases i <;> simp [curlLamR] <;> linarith

/-- **(1) Eigen receipt.**  `CᵀC (curlR i) = curlLamR i • curlR i`. -/
theorem curl_eigen (i : Fin 19) : localMaxwellOperator (curlR i) = curlLamR i • curlR i := by
  rw [← evalPhi_curlLam]
  exact localMaxwellOperator_evalPhi (curlZ i) (curlFaceZ i) (curlLamZ i)
    (curl_curvature_Z i) (curl_codifferential_Z i)

theorem curl_inner (i j : Fin 19) :
    realSeamInner (curlR i) (curlR j) = evalPhi (if i = j then curlNormZ i else 0) := by
  unfold curlR
  rw [realSeamInner_evalPhi, curl_gram_Z i j]

/-- **(1) Pairwise orthogonality.** -/
theorem curl_orth (i j : Fin 19) (hij : i ≠ j) : realSeamInner (curlR i) (curlR j) = 0 := by
  rw [curl_inner, if_neg hij, evalPhi_zero]

theorem curl_energy (i : Fin 19) : realSeamEnergy (curlR i) = evalPhi (curlNormZ i) := by
  rw [← realSeamInner_self_eq_energy, curl_inner, if_pos rfl]

/-- **(1) Nonzero.** -/
theorem curl_ne_zero (i : Fin 19) : curlR i ≠ 0 := by
  intro h
  have h1 := curl_energy i
  rw [h, realSeamEnergy_zero] at h1
  exact evalPhi_ne_zero (curlNormZ_ne_zero i) h1.symm

/-- Pairwise orthogonal nonzero seam vectors are linearly independent. -/
theorem linearIndependent_of_orth {ι : Type} [Fintype ι] (v : ι → Fin 30 → ℝ)
    (hne : ∀ i, v i ≠ 0) (horth : ∀ i j, i ≠ j → realSeamInner (v i) (v j) = 0) :
    LinearIndependent ℝ v := by
  rw [Fintype.linearIndependent_iff]
  intro g hg j
  have h1 : realSeamInner (∑ i, g i • v i) (v j) = g j * realSeamEnergy (v j) := by
    rw [realSeamInner_sum_left, Finset.sum_eq_single j]
    · rw [seamInner_smul_left, realSeamInner_self_eq_energy]
    · intro i _ hij
      rw [seamInner_smul_left, horth i j hij, mul_zero]
    · intro h
      exact absurd (Finset.mem_univ j) h
  rw [hg, realSeamInner_zero_left] at h1
  have hpos := seamEnergy_pos_of_ne_zero (v j) (hne j)
  rcases mul_eq_zero.mp h1.symm with h | h
  · exact h
  · exact absurd h hpos.ne'

/-- **(1) Linear independence** of the nineteen. -/
theorem curl_linearIndependent : LinearIndependent ℝ curlR :=
  linearIndependent_of_orth curlR curl_ne_zero curl_orth

/-- **(1) Membership in the curl sector**, exhibited: `curlR i = Cᵀ (λ⁻¹ • F_i)`
with `F_i` the listed face vector. -/
theorem curl_eq_codifferential (i : Fin 19) :
    curlR i = faceCodifferential ((curlLamR i)⁻¹ • fun f ↦ evalPhi (curlFaceZ i f)) := by
  rw [map_smul, faceCodifferential_evalPhi _ _ (curl_codifferential_Z i)]
  funext e
  simp only [curlR, Pi.smul_apply, smul_eq_mul, evalPhi_zmul, evalPhi_curlLam]
  rw [← mul_assoc, inv_mul_cancel₀ (curlLamR_ne_zero i), one_mul]

theorem curl_mem_range_codifferential (i : Fin 19) :
    curlR i ∈ LinearMap.range faceCodifferential :=
  ⟨_, (curl_eq_codifferential i).symm⟩

/-- **(1) Orthogonal to every gradient seam field.** -/
theorem curl_orth_gradient (i : Fin 19) (χ : Fin 12 → ℝ) :
    realSeamInner (curlR i) (realCoboundary χ) = 0 :=
  eigen_orthogonal _ _ _ _ (curl_eigen i) (gradient_eigen χ) (curlLamR_ne_zero i)

/-- The range of `Cᵀ` lies in the kernel of the committed boundary. -/
theorem range_codifferential_le_ker_boundary :
    LinearMap.range faceCodifferential ≤ LinearMap.ker realBoundary := by
  rintro _ ⟨w, rfl⟩
  rw [LinearMap.mem_ker]
  apply (realPortInner_self_eq_zero_iff _).mp
  rw [← realCoboundary_boundary_adjoint, ← faceCurvature_codifferential_adjoint,
    faceCurvature_coboundary]
  simp [faceInner]

theorem span_curl_le : Submodule.span ℝ (Set.range curlR) ≤ LinearMap.range faceCodifferential := by
  rw [Submodule.span_le]
  rintro _ ⟨i, rfl⟩
  exact curl_mem_range_codifferential i

theorem span_curl_finrank :
    Module.finrank ℝ (Submodule.span ℝ (Set.range curlR)) = 19 := by
  rw [finrank_span_eq_card curl_linearIndependent, Fintype.card_fin]

/-- The range of `Cᵀ` has dimension `19`. -/
theorem range_codifferential_finrank :
    Module.finrank ℝ (LinearMap.range faceCodifferential) = 19 := by
  apply le_antisymm
  · have h := Submodule.finrank_mono range_codifferential_le_ker_boundary
    rwa [cycleSpace_finrank] at h
  · have h := Submodule.finrank_mono span_curl_le
    rwa [span_curl_finrank] at h

/-- **(1) The nineteen span exactly the curl sector**, the range of `Cᵀ`. -/
theorem curl_span : Submodule.span ℝ (Set.range curlR) = LinearMap.range faceCodifferential :=
  Submodule.eq_of_le_of_finrank_eq span_curl_le
    (by rw [span_curl_finrank, range_codifferential_finrank])

/-- **(1) Three descriptions of the curl sector coincide**: the range of `Cᵀ`,
the range of `CᵀC`, and the kernel of the committed boundary. -/
theorem curl_sector_characterizations :
    LinearMap.range faceCodifferential = LinearMap.ker realBoundary ∧
      LinearMap.range faceCodifferential = LinearMap.range localMaxwellOperator := by
  have h1 : LinearMap.range faceCodifferential = LinearMap.ker realBoundary :=
    Submodule.eq_of_le_of_finrank_eq range_codifferential_le_ker_boundary
      (by rw [range_codifferential_finrank, cycleSpace_finrank])
  exact ⟨h1, by rw [h1, range_localMaxwellOperator]⟩

/-- **(1) The committed modes are members**: rows `0, 5, 9, 13` are
`twoMode, threeMode, fiveMode, goldenMode`. -/
theorem curl_committed_members :
    curlR 0 = twoMode ∧ curlR 5 = threeMode ∧ curlR 9 = fiveMode ∧ curlR 13 = goldenMode := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · funext e
    simp only [curlR, (curl_committed_rows_Z e).1, twoMode, evalPhi]
    push_cast; ring
  · funext e
    simp only [curlR, (curl_committed_rows_Z e).2.1, threeMode, evalPhi]
    push_cast; ring
  · funext e
    simp only [curlR, (curl_committed_rows_Z e).2.2.1, fiveMode, evalPhi]
    push_cast; ring
  · funext e
    simp only [curlR, (curl_committed_rows_Z e).2.2.2, goldenMode, goldenA, goldenB,
      Pi.add_apply, Pi.smul_apply, smul_eq_mul, evalPhi]
    unfold Real.goldenRatio
    push_cast; ring

/-! ## 5. The eleven real gradient vectors and the exact complement -/

/-- The real cast of the `i`-th listed port load. -/
def gradPortR (i : Fin 11) : Fin 12 → ℝ := fun p ↦ evalPhi (gradPortZ i p)

/-- The `i`-th gradient seam vector: the coboundary of the listed port load. -/
def gradR (i : Fin 11) : Fin 30 → ℝ := realCoboundary (gradPortR i)

theorem gradR_eq (i : Fin 11) : gradR i = fun e ↦ evalPhi (gradZ i e) := by
  funext e
  show gradPortR i (seamRight e) - gradPortR i (seamLeft e) = _
  unfold gradZ gradPortR
  rw [evalPhi_zsub]

theorem grad_inner (i j : Fin 11) :
    realSeamInner (gradR i) (gradR j) = evalPhi (if i = j then gradNormZ i else 0) := by
  rw [gradR_eq, gradR_eq, realSeamInner_evalPhi, grad_gram_Z i j]

theorem grad_orth (i j : Fin 11) (hij : i ≠ j) : realSeamInner (gradR i) (gradR j) = 0 := by
  rw [grad_inner, if_neg hij, evalPhi_zero]

theorem grad_energy (i : Fin 11) : realSeamEnergy (gradR i) = evalPhi (gradNormZ i) := by
  rw [← realSeamInner_self_eq_energy, grad_inner, if_pos rfl]

theorem grad_ne_zero (i : Fin 11) : gradR i ≠ 0 := by
  intro h
  have h1 := grad_energy i
  rw [h, realSeamEnergy_zero] at h1
  exact evalPhi_ne_zero (gradNormZ_ne_zero i) h1.symm

theorem grad_linearIndependent : LinearIndependent ℝ gradR :=
  linearIndependent_of_orth gradR grad_ne_zero grad_orth

theorem grad_mem_range (i : Fin 11) : gradR i ∈ LinearMap.range realCoboundary :=
  ⟨gradPortR i, rfl⟩

/-- **(2) The eleven span exactly the gradient sector.** -/
theorem grad_span : Submodule.span ℝ (Set.range gradR) = LinearMap.range realCoboundary := by
  apply Submodule.eq_of_le_of_finrank_eq
  · rw [Submodule.span_le]
    rintro _ ⟨i, rfl⟩
    exact grad_mem_range i
  · rw [finrank_span_eq_card grad_linearIndependent, Fintype.card_fin, gauge_finrank]

/-- **(3) Disjointness** of the curl and gradient sectors, by orthogonality. -/
theorem curl_grad_disjoint :
    Disjoint (LinearMap.range faceCodifferential) (LinearMap.range realCoboundary) := by
  rw [Submodule.disjoint_def]
  rintro v ⟨w, rfl⟩ ⟨χ, hχ⟩
  apply (realSeamEnergy_eq_zero_iff _).mp
  rw [← realSeamInner_self_eq_energy]
  have h : realSeamInner (faceCodifferential w) (faceCodifferential w) =
      realSeamInner (realCoboundary χ) (faceCodifferential w) := by rw [hχ]
  rw [h, ← faceCurvature_codifferential_adjoint, faceCurvature_coboundary]
  simp [faceInner]

/-- **(3) The two sectors sum to the whole seam space**, by `19 + 11 = 30`. -/
theorem curl_sup_grad_eq_top :
    LinearMap.range faceCodifferential ⊔ LinearMap.range realCoboundary = ⊤ := by
  apply Submodule.eq_top_of_finrank_eq
  have h := Submodule.finrank_sup_add_finrank_inf_eq (LinearMap.range faceCodifferential)
    (LinearMap.range realCoboundary)
  rw [curl_grad_disjoint.eq_bot, finrank_bot, range_codifferential_finrank, gauge_finrank] at h
  rw [Module.finrank_fin_fun]
  omega

/-- **(3) Exact complement.** -/
theorem curl_grad_isCompl :
    IsCompl (LinearMap.range faceCodifferential) (LinearMap.range realCoboundary) :=
  ⟨curl_grad_disjoint, codisjoint_iff.mpr curl_sup_grad_eq_top⟩

/-! ## 6. The thirty-vector family and the flow on the whole seam space -/

/-- The thirty-vector family: nineteen curl vectors and eleven gradient vectors. -/
def fullR : Fin 19 ⊕ Fin 11 → Fin 30 → ℝ := Sum.elim curlR gradR

/-- Its eigenvalues: the listed curl eigenvalues and `0` on the gradient part. -/
def fullLam : Fin 19 ⊕ Fin 11 → ℝ := Sum.elim curlLamR fun _ ↦ 0

theorem full_eigen (i : Fin 19 ⊕ Fin 11) :
    localMaxwellOperator (fullR i) = fullLam i • fullR i := by
  rcases i with i | i
  · exact curl_eigen i
  · exact gradient_eigen (gradPortR i)

theorem full_orth (i j : Fin 19 ⊕ Fin 11) (hij : i ≠ j) :
    realSeamInner (fullR i) (fullR j) = 0 := by
  rcases i with i | i <;> rcases j with j | j
  · exact curl_orth i j fun h ↦ hij (congrArg Sum.inl h)
  · exact curl_orth_gradient i (gradPortR j)
  · rw [realSeamInner_comm]
    exact curl_orth_gradient j (gradPortR i)
  · exact grad_orth i j fun h ↦ hij (congrArg Sum.inr h)

theorem full_ne_zero (i : Fin 19 ⊕ Fin 11) : fullR i ≠ 0 := by
  rcases i with i | i
  · exact curl_ne_zero i
  · exact grad_ne_zero i

theorem full_linearIndependent : LinearIndependent ℝ fullR :=
  linearIndependent_of_orth fullR full_ne_zero full_orth

/-- **(4) The thirty vectors span the whole seam space.** -/
theorem full_span_top : Submodule.span ℝ (Set.range fullR) = ⊤ := by
  unfold fullR
  rw [Set.Sum.elim_range, Submodule.span_union, curl_span, grad_span, curl_sup_grad_eq_top]

/-- **(4) Every pair of seam fields is the potential and electric readout of a
coefficient state.** -/
theorem exists_state_of_fields (A₀ E₀ : Fin 30 → ℝ) :
    ∃ x : Fin 19 ⊕ Fin 11 → Fin 2 → ℝ, potentialOf fullR x = A₀ ∧ electricOf fullR x = E₀ := by
  have hA : A₀ ∈ Submodule.span ℝ (Set.range fullR) := by
    rw [full_span_top]; exact Submodule.mem_top
  have hE : -E₀ ∈ Submodule.span ℝ (Set.range fullR) := by
    rw [full_span_top]; exact Submodule.mem_top
  obtain ⟨a, ha⟩ := (Submodule.mem_span_range_iff_exists_fun ℝ).mp hA
  obtain ⟨b, hb⟩ := (Submodule.mem_span_range_iff_exists_fun ℝ).mp hE
  refine ⟨fun i ↦ ![a i, b i], ?_, ?_⟩
  · unfold potentialOf
    simpa using ha
  · unfold electricOf
    simp only [Matrix.cons_val_one, Matrix.cons_val_fin_one]
    rw [hb, neg_neg]

/-- **(5) Admissibility of every eigenvalue** inside the window `h² (3 + √5) < 4`. -/
theorem fullLam_admissible (h : ℝ) (hh : h ≠ 0) (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4)
    (i : Fin 19 ⊕ Fin 11) : Admissible h (fullLam i) := by
  rcases i with i | i
  · right
    have hp : 0 < h ^ 2 := by positivity
    refine ⟨mul_pos hp (curlLamR_pos i), lt_of_le_of_lt ?_ h4⟩
    exact mul_le_mul_of_nonneg_left (curlLamR_le_golden i) hp.le
  · left
    rfl

/-- **(5) The window is two-sided for the family**: with `h ≠ 0`, every
eigenvalue of the thirty-vector family is admissible (in the committed sense
`Admissible`) iff `h² (3 + √5) < 4`; the forward direction is read off the
row `13` eigenvalue `3 + √5`. -/
theorem fullLam_admissible_iff (h : ℝ) (hh : h ≠ 0) :
    (∀ i, Admissible h (fullLam i)) ↔ h ^ 2 * (3 + Real.sqrt 5) < 4 := by
  refine ⟨fun hall ↦ ?_, fullLam_admissible h hh⟩
  have h13 := hall (Sum.inl 13)
  have hl : fullLam (Sum.inl 13) = 3 + Real.sqrt 5 := by
    simp [fullLam, curlLamR]
  rw [hl] at h13
  rcases h13 with h0 | ⟨_, h4⟩
  · have := sqrt5_pos; linarith
  · exact h4

/-- Every one of the nineteen nonzero curl eigenvalues lies in the strict
energy window whenever the largest eigenvalue does. -/
theorem curlLamR_window (h : ℝ) (hh : h ≠ 0)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) (i : Fin 19) :
    0 < h ^ 2 * curlLamR i ∧ h ^ 2 * curlLamR i < 4 := by
  have hp : 0 < h ^ 2 := by positivity
  exact ⟨mul_pos hp (curlLamR_pos i),
    lt_of_le_of_lt (mul_le_mul_of_nonneg_left (curlLamR_le_golden i) hp.le) h4⟩

/-- **Energy/Hermitian bridge for the full curl sector.**  Inside the exact
global window, the nineteen listed nonzero curl modes instantiate the
energy-derived Hilbert reading.  This is a theorem about the finite
coefficient space and its field readouts; it makes no physical-time,
private-algebra, or photon-Hilbert-space identification. -/
theorem curl_sector_energy_hilbert_packet (h : ℝ) (hh : h ≠ 0)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) :
    (∀ x, x ≠ 0 → 0 < assembledInner h curlLamR curlR x x) ∧
    (∀ t x y, assembledInner h curlLamR curlR (assembledFlow h curlLamR t x)
      (assembledFlow h curlLamR t y) = assembledInner h curlLamR curlR x y) ∧
    (∀ x n, fieldEnergyScaled h (assembledHistory h curlLamR curlR x) (fun _ ↦ 0) n =
      assembledInner h curlLamR curlR x x) ∧
    (∀ t x i, assembledCoordinate h curlLamR (assembledFlow h curlLamR t x) i =
      Complex.exp ((modeAngle h (curlLamR i) * t / h : ℝ) * Complex.I) *
        assembledCoordinate h curlLamR x i) ∧
    (∀ x i t, HasDerivAt
      (fun s : ℝ ↦ assembledCoordinate h curlLamR (assembledFlow h curlLamR s x) i)
      (Complex.I * (modeAngle h (curlLamR i) / h : ℝ) *
        assembledCoordinate h curlLamR (assembledFlow h curlLamR t x) i) t) ∧
    (∀ t x y, assembledHermitian h curlLamR curlR (assembledFlow h curlLamR t x)
      (assembledFlow h curlLamR t y) = assembledHermitian h curlLamR curlR x y) ∧
    (∀ x, assembledHermitian h curlLamR curlR x x =
      (assembledInner h curlLamR curlR x x : ℂ)) ∧
    (∀ x y, potentialOf curlR x = potentialOf curlR y →
      electricOf curlR x = electricOf curlR y → x = y) := by
  obtain ⟨hpos, _, _, _, hflow, henergy, hphase, hderiv, hherm, hself, hrecover⟩ :=
    orthogonal_family_hilbert_reading h curlLamR curlR hh (curlLamR_window h hh h4)
      curl_eigen curl_ne_zero curl_orth
  exact ⟨hpos, hflow, henergy, hphase, hderiv, hherm, hself, hrecover⟩

/-- **Two-step form of the zero-current temporal-gauge evolution:**
`A (n+2) = 2 A (n+1) - A n - h² CᵀC A (n+1)`. -/
theorem ampere_step_eq (h : ℝ) (hh : h ≠ 0) (A : ℕ → Fin 30 → ℝ)
    (hA : AmpereEvolutionScaled h A (fun _ ↦ 0) (fun _ ↦ 0)) (n : ℕ) :
    A (n + 2) = (2 : ℝ) • A (n + 1) - A n - (h ^ 2) • localMaxwellOperator (A (n + 1)) := by
  have hn := hA n
  have hL : faceCodifferential (magneticField A (n + 1)) = localMaxwellOperator (A (n + 1)) := rfl
  unfold electricFieldScaled at hn
  rw [map_zero, sub_zero, sub_zero, hL, sub_zero] at hn
  funext e
  have he := congrFun hn e
  simp only [Pi.sub_apply, Pi.neg_apply, Pi.smul_apply, smul_eq_mul] at he ⊢
  have hinv : h * h⁻¹ = 1 := mul_inv_cancel₀ hh
  rw [show n + 1 + 1 = n + 2 from rfl] at he
  linear_combination (-h) * he -
    ((A (n + 2) e - A (n + 1) e) - (A (n + 1) e - A n e)) * hinv

/-- **Uniqueness:** two zero-current temporal-gauge solutions with the same
first two values coincide. -/
theorem ampere_unique (h : ℝ) (hh : h ≠ 0) (A B : ℕ → Fin 30 → ℝ)
    (hA : AmpereEvolutionScaled h A (fun _ ↦ 0) (fun _ ↦ 0))
    (hB : AmpereEvolutionScaled h B (fun _ ↦ 0) (fun _ ↦ 0))
    (h0 : A 0 = B 0) (h1 : A 1 = B 1) : A = B := by
  have key : ∀ n, A n = B n ∧ A (n + 1) = B (n + 1) := by
    intro n
    induction n with
    | zero => exact ⟨h0, h1⟩
    | succ n ih =>
      refine ⟨ih.2, ?_⟩
      show A (n + 2) = B (n + 2)
      rw [ampere_step_eq h hh A hA n, ampere_step_eq h hh B hB n, ih.1, ih.2]
  funext n
  exact (key n).1

/-- **(5) Every zero-current temporal-gauge history is an assembled history**
of the thirty-vector family: the assembled flow reproduces every solution
of the committed evolution at the step values. -/
theorem ampere_history_eq_assembled (h : ℝ) (hh : h ≠ 0)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) (A : ℕ → Fin 30 → ℝ)
    (hA : AmpereEvolutionScaled h A (fun _ ↦ 0) (fun _ ↦ 0)) :
    ∃ x : Fin 19 ⊕ Fin 11 → Fin 2 → ℝ, A = assembledHistory h fullLam fullR x := by
  have hadm := fullLam_admissible h hh h4
  obtain ⟨x, hx0, hxE⟩ := exists_state_of_fields (A 0) (electricFieldScaled h A (fun _ ↦ 0) 0)
  have hX0 : assembledHistory h fullLam fullR x 0 = A 0 := by
    unfold assembledHistory
    rw [Nat.cast_zero, zero_mul, assembledFlow_zero h fullLam hh hadm, hx0]
  refine ⟨x, ampere_unique h hh A _ hA
    (assembledHistory_ampere h fullLam fullR hh hadm full_eigen x) hX0.symm ?_⟩
  have hE := assembledHistory_electricField h fullLam fullR hh hadm x 0
  rw [Nat.cast_zero, zero_mul, assembledFlow_zero h fullLam hh hadm, hxE] at hE
  unfold electricFieldScaled at hE
  rw [map_zero, sub_zero, sub_zero, zero_add, hX0] at hE
  have h1 := smul_right_injective (Fin 30 → ℝ) (inv_ne_zero hh) (neg_injective hE)
  exact (sub_left_inj.mp h1).symm

/-- **(5) Carrier flow on the whole seam space.**  On the thirty-vector
orthogonal eigenbasis (nineteen curl vectors spanning the range of `Cᵀ`,
eleven gradient vectors spanning its exact complement), at any declared
step `h ≠ 0` inside the global window `h² (3 + √5) < 4`, the assembled
flow packet of `CarrierEvolutionFlow.assembled_flow` holds; the family
spans `Fin 30 → ℝ`; every initial pair of seam fields is realised by a
coefficient state; and every zero-current temporal-gauge solution of the
committed evolution is the assembled history of some coefficient state,
so the `h`-step of the flow is the committed scaled Ampere step on every
zero-current history. -/
theorem carrier_flow_full_curl (h : ℝ) (hh : h ≠ 0)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) :
    (∀ x, assembledFlow h fullLam 0 x = x) ∧
    (∀ s t x, assembledFlow h fullLam s (assembledFlow h fullLam t x) =
      assembledFlow h fullLam (s + t) x) ∧
    (∀ t x y, assembledFlow h fullLam t (x + y) =
      assembledFlow h fullLam t x + assembledFlow h fullLam t y) ∧
    (∀ (t c : ℝ) x, assembledFlow h fullLam t (c • x) = c • assembledFlow h fullLam t x) ∧
    (∀ x, Continuous fun t ↦ assembledFlow h fullLam t x) ∧
    (∀ x i, assembledFlow h fullLam h x i = (stepMatrix h (fullLam i)).mulVec (x i)) ∧
    (∀ t x, assembledEnergy h fullLam fullR (assembledFlow h fullLam t x) =
      assembledEnergy h fullLam fullR x) ∧
    (∀ x, AmpereEvolutionScaled h (assembledHistory h fullLam fullR x)
      (fun _ ↦ 0) (fun _ ↦ 0)) ∧
    (∀ x n, electricFieldScaled h (assembledHistory h fullLam fullR x) (fun _ ↦ 0) n =
      electricOf fullR (assembledFlow h fullLam (n * h) x)) ∧
    (∀ x n, fieldEnergyScaled h (assembledHistory h fullLam fullR x) (fun _ ↦ 0) n =
      assembledEnergy h fullLam fullR x) ∧
    Submodule.span ℝ (Set.range fullR) = ⊤ ∧
    (∀ A₀ E₀ : Fin 30 → ℝ, ∃ x, potentialOf fullR x = A₀ ∧ electricOf fullR x = E₀) ∧
    (∀ A : ℕ → Fin 30 → ℝ, AmpereEvolutionScaled h A (fun _ ↦ 0) (fun _ ↦ 0) →
      ∃ x, A = assembledHistory h fullLam fullR x) := by
  obtain ⟨p1, p2, p3, p4, p5, p6, p7, p8, p9, p10⟩ :=
    assembled_flow h fullLam fullR hh (fullLam_admissible h hh h4) full_eigen full_orth
  exact ⟨p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, full_span_top, exists_state_of_fields,
    fun A hA ↦ ampere_history_eq_assembled h hh h4 A hA⟩

/-- The five committed families of `carrier_flow_five_modes` lie in the span
of the thirty-vector family (which is everything), with the four named
modes literally among the rows. -/
theorem fiveFamily_mem_span (χ : Fin 12 → ℝ) (i : Fin 5) :
    fiveFamily χ i ∈ Submodule.span ℝ (Set.range fullR) := by
  rw [full_span_top]; exact Submodule.mem_top

/-- The five committed families placed in their sectors: `twoMode, threeMode,
fiveMode, goldenMode` (rows `0, 5, 9, 13`) lie in the span of the nineteen
curl vectors, and the gradient `d χ` lies in the span of the eleven gradient
vectors. -/
theorem fiveFamily_mem_sector (χ : Fin 12 → ℝ) :
    fiveFamily χ 0 ∈ Submodule.span ℝ (Set.range curlR) ∧
    fiveFamily χ 1 ∈ Submodule.span ℝ (Set.range curlR) ∧
    fiveFamily χ 2 ∈ Submodule.span ℝ (Set.range curlR) ∧
    fiveFamily χ 3 ∈ Submodule.span ℝ (Set.range curlR) ∧
    fiveFamily χ 4 ∈ Submodule.span ℝ (Set.range gradR) := by
  obtain ⟨h0, h5, h9, h13⟩ := curl_committed_members
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · show twoMode ∈ _
    rw [← h0]; exact Submodule.subset_span ⟨0, rfl⟩
  · show threeMode ∈ _
    rw [← h5]; exact Submodule.subset_span ⟨5, rfl⟩
  · show fiveMode ∈ _
    rw [← h9]; exact Submodule.subset_span ⟨9, rfl⟩
  · show goldenMode ∈ _
    rw [← h13]; exact Submodule.subset_span ⟨13, rfl⟩
  · show realCoboundary χ ∈ _
    rw [grad_span]; exact ⟨χ, rfl⟩

end

end OPH.CurlSectorEigenbasis

#print axioms OPH.CurlSectorEigenbasis.curl_curvature_Z
#print axioms OPH.CurlSectorEigenbasis.curl_codifferential_Z
#print axioms OPH.CurlSectorEigenbasis.curl_gram_Z
#print axioms OPH.CurlSectorEigenbasis.grad_gram_Z
#print axioms OPH.CurlSectorEigenbasis.curl_eigen
#print axioms OPH.CurlSectorEigenbasis.curl_orth
#print axioms OPH.CurlSectorEigenbasis.curl_ne_zero
#print axioms OPH.CurlSectorEigenbasis.curl_linearIndependent
#print axioms OPH.CurlSectorEigenbasis.curl_eq_codifferential
#print axioms OPH.CurlSectorEigenbasis.curl_orth_gradient
#print axioms OPH.CurlSectorEigenbasis.curl_span
#print axioms OPH.CurlSectorEigenbasis.curl_sector_characterizations
#print axioms OPH.CurlSectorEigenbasis.curl_committed_members
#print axioms OPH.CurlSectorEigenbasis.curl_minus_eq_conj_plus
#print axioms OPH.CurlSectorEigenbasis.grad_span
#print axioms OPH.CurlSectorEigenbasis.curl_grad_isCompl
#print axioms OPH.CurlSectorEigenbasis.full_span_top
#print axioms OPH.CurlSectorEigenbasis.exists_state_of_fields
#print axioms OPH.CurlSectorEigenbasis.fullLam_admissible
#print axioms OPH.CurlSectorEigenbasis.fullLam_admissible_iff
#print axioms OPH.CurlSectorEigenbasis.curlLamR_window
#print axioms OPH.CurlSectorEigenbasis.curl_sector_energy_hilbert_packet
#print axioms OPH.CurlSectorEigenbasis.ampere_step_eq
#print axioms OPH.CurlSectorEigenbasis.ampere_unique
#print axioms OPH.CurlSectorEigenbasis.ampere_history_eq_assembled
#print axioms OPH.CurlSectorEigenbasis.carrier_flow_full_curl
#print axioms OPH.CurlSectorEigenbasis.fiveFamily_mem_span
#print axioms OPH.CurlSectorEigenbasis.fiveFamily_mem_sector
