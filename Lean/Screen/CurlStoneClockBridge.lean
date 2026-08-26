import CurlSectorEigenbasis
import Dynamics.StoneConverse
import Geometry.CommonWorldMaxwellClockJoin

set_option autoImplicit false

open scoped BigOperators Matrix

namespace OPH.CurlStoneClockBridge

open OPH.CarrierModeOscillators
open OPH.CarrierEvolutionFlow
open OPH.FieldSectorEnergyInnerProduct
open OPH.CurlSectorEigenbasis
open OPH.Dynamics
open OPH.CommonWorldInstrumentJoin

noncomputable section

/-!
# Conditional curl-sector Stone and clock bridge

This module compares three already-constructed one-parameter objects without
identifying them merely because their parameters are all written `t`:

* the energy-derived Hermitian coefficient flow on the nineteen exact curl
  modes of `CurlSectorEigenbasis`;
* a finite full matrix algebra and its Stone propagator from
  `Dynamics/StoneConverse`;
* the declared internal-clock readout on an instrumented common-world record.

The comparison has an exact algebraic part and a deliberately conditional
clock part.  Complex conjugation selects the Schrödinger orientation of the
coefficient coordinates.  Their nineteen frequencies form a positive
diagonal self-adjoint Hamiltonian.  Its Stone conjugation carries the
coefficient outer product exactly to the outer product of the curl flow.  Thus
the vector flow and its rank-one projective matrix image are genuinely
intertwined, not matched by dimension or notation.  The outer product loses
global phase and is neither onto the full matrix algebra nor an algebra map.

For a selected mode, the common-world clock drives that coefficient at the
discrete Maxwell steps provided two explicit sufficient calibration hypotheses
are supplied: the record's step duration equals the Courant step, and its one
positive mass/rate equals that mode's selected principal-branch frequency.  A
sharp obstruction shows only that one identical scalar rate cannot equal all
of those principal-branch frequencies: already the eigenvalue-2 and
eigenvalue-3 sectors differ at every positive admissible step.

No physical identification is proved.  The diagonal matrix algebra is a new
mathematical realization of the coefficient space, not the committed or a
source-selected private algebra; the principal logarithm branch and calibration
equalities are choices/hypotheses, not consequences;
the clock remains dimensionless and uncalibrated; and no Born readout,
laboratory frequency, photon identification, or source provenance is supplied.
-/

/-- The selected principal-angle interpolation frequency of curl mode `i`.
The discrete `h`-step also admits logarithm branches shifted by integer
multiples of `2π / h`; this definition does not prove a unique physical
frequency. -/
def curlFrequency (h : ℝ) (i : Fin 19) : ℝ :=
  modeAngle h (curlLamR i) / h

/-- The conjugate coefficient coordinate, choosing the Schrödinger phase
`exp (-i ωᵢ t)` rather than the equally valid conjugate orientation
`exp (+i ωᵢ t)` used by `assembledCoordinate`. -/
def schrodingerCoordinate (h : ℝ) (x : Fin 19 → Fin 2 → ℝ) (i : Fin 19) : ℂ :=
  (starRingEnd ℂ) (assembledCoordinate h curlLamR x i)

/-- The nineteen Schrödinger-oriented coordinates as one complex vector. -/
def schrodingerVector (h : ℝ) (x : Fin 19 → Fin 2 → ℝ) : Fin 19 → ℂ :=
  fun i ↦ schrodingerCoordinate h x i

/-- Real scale that converts the raw coefficient coordinate into the standard
complex coordinate whose squared norm is the corresponding term of the
energy-derived assembled Hermitian form. -/
def energyAmplitudeScale (h : ℝ) (i : Fin 19) : ℝ :=
  Real.sqrt (OPH.DiscreteCoulombGreen.realSeamEnergy (curlR i) / 2) *
    Real.sin (modeAngle h (curlLamR i))

/-- The energy-normalized (but not unit-normalized) Schrödinger coefficient
vector.  Its standard squared norm is the committed assembled field energy. -/
def energySchrodingerVector
    (h : ℝ) (x : Fin 19 → Fin 2 → ℝ) : Fin 19 → ℂ :=
  fun i ↦ (energyAmplitudeScale h i : ℂ) * schrodingerCoordinate h x i

/-- The exact Schrödinger phase of mode `i`. -/
def curlPhase (h t : ℝ) (i : Fin 19) : ℂ :=
  Complex.exp (-Complex.I * ((curlFrequency h i * t : ℝ) : ℂ))

/-- The nineteen-mode coefficient Hamiltonian on one newly constructed full
matrix block.  It is diagonal in the exact curl eigenbasis. -/
def curlHamiltonian (h : ℝ) : Matrix (Fin 19) (Fin 19) ℂ :=
  Matrix.diagonal fun i ↦ (curlFrequency h i : ℂ)

/-- The coefficient outer product `z z†` in the full matrix block.
Normalization is not imposed, so this is not called a density matrix; the
identity below is homogeneous and applies also to the zero state. -/
def coefficientOuter (z : Fin 19 → ℂ) : Matrix (Fin 19) (Fin 19) ℂ :=
  fun i j ↦ z i * (starRingEnd ℂ) (z j)

/-- The coefficient outer product forgets global phase.  It therefore gives a
projective rank-one matrix image, not an injective vector encoding or an
identification with every observable of the full matrix block. -/
theorem coefficientOuter_global_phase_invariant
    (alpha : ℝ) (z : Fin 19 → ℂ) :
    coefficientOuter
        (fun i ↦ Complex.exp (Complex.I * (alpha : ℂ)) * z i) =
      coefficientOuter z := by
  ext i j
  unfold coefficientOuter
  rw [map_mul, ← Complex.exp_conj]
  have hconj :
      starRingEnd ℂ (Complex.I * (alpha : ℂ)) =
        -(Complex.I * (alpha : ℂ)) := by
    simp
  rw [hconj]
  calc
    _ = (Complex.exp (Complex.I * (alpha : ℂ)) *
          Complex.exp (-(Complex.I * (alpha : ℂ)))) *
        (z i * (starRingEnd ℂ) (z j)) := by ring
    _ = _ := by rw [← Complex.exp_add]; simp

/-- The coefficient Hamiltonian is self-adjoint because every diagonal entry
is real. -/
theorem curlHamiltonian_isSelfAdjoint (h : ℝ) :
    IsSelfAdjoint (curlHamiltonian h) := by
  rw [isSelfAdjoint_iff]
  ext i j
  rw [Matrix.star_apply]
  by_cases hij : i = j
  · subst j
    simp [curlHamiltonian]
  · simp [curlHamiltonian, hij, Ne.symm hij]

/-- Inside the strict Courant window and for positive `h`, every diagonal
frequency is strictly positive. -/
theorem curlFrequency_pos (h : ℝ) (hh : 0 < h)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) (i : Fin 19) :
    0 < curlFrequency h i := by
  have hw := curlLamR_window h hh.ne' h4 i
  have hs := sin_modeAngle_pos h (curlLamR i) hw.1 hw.2
  have hn := modeAngle_nonneg h (curlLamR i)
  have hne : modeAngle h (curlLamR i) ≠ 0 := by
    intro hz
    rw [hz] at hs
    simp at hs
  have htheta : 0 < modeAngle h (curlLamR i) := lt_of_le_of_ne hn (Ne.symm hne)
  exact div_pos htheta hh

/-- Conjugating the energy coordinate gives the exact Schrödinger-sign phase
law for each of the nineteen curl modes. -/
theorem schrodingerCoordinate_flow (h : ℝ) (hh : h ≠ 0)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) (t : ℝ)
    (x : Fin 19 → Fin 2 → ℝ) (i : Fin 19) :
    schrodingerCoordinate h (assembledFlow h curlLamR t x) i =
      curlPhase h t i * schrodingerCoordinate h x i := by
  rw [schrodingerCoordinate, schrodingerCoordinate,
    assembledCoordinate_flow h curlLamR hh (curlLamR_window h hh h4) t x i,
    map_mul, ← Complex.exp_conj]
  congr 1
  simp [curlPhase, curlFrequency]
  ring_nf

/-- The energy-normalized coefficient vector obeys the same diagonal phase
law; the energy weights are fixed by the exact curl eigenbasis. -/
theorem energySchrodingerVector_flow (h : ℝ) (hh : h ≠ 0)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) (t : ℝ)
    (x : Fin 19 → Fin 2 → ℝ) (i : Fin 19) :
    energySchrodingerVector h (assembledFlow h curlLamR t x) i =
      curlPhase h t i * energySchrodingerVector h x i := by
  unfold energySchrodingerVector
  rw [schrodingerCoordinate_flow h hh h4]
  ring_nf

/-- The trace of the energy-weighted coefficient outer product is exactly the
assembled Hermitian diagonal, hence the committed assembled field energy.
This is the metric contact point between the field-sector Hilbert packet and
the standard matrix-block norm. -/
theorem energyCoefficientOuter_trace (h : ℝ) (hh : h ≠ 0)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4)
    (x : Fin 19 → Fin 2 → ℝ) :
    (coefficientOuter (energySchrodingerVector h x)).trace =
      assembledHermitian h curlLamR curlR x x := by
  rw [assembledHermitian_eq_coordinate h curlLamR curlR hh
    (curlLamR_window h hh h4)]
  unfold Matrix.trace Matrix.diag coefficientOuter energySchrodingerVector
    energyAmplitudeScale schrodingerCoordinate
  apply Finset.sum_congr rfl
  intro i _
  have henergy : 0 ≤ OPH.DiscreteCoulombGreen.realSeamEnergy (curlR i) / 2 := by
    exact div_nonneg (OPH.DiscreteCoulombGreen.realSeamEnergy_nonneg _) (by norm_num)
  have hsqrt :
      ((Real.sqrt (OPH.DiscreteCoulombGreen.realSeamEnergy (curlR i) / 2) : ℝ) : ℂ) ^ 2 =
        ((OPH.DiscreteCoulombGreen.realSeamEnergy (curlR i) / 2 : ℝ) : ℂ) := by
    exact_mod_cast Real.sq_sqrt henergy
  rw [map_mul, Complex.conj_ofReal, Complex.conj_conj]
  simp_rw [Complex.ofReal_mul]
  calc
    _ = ((Real.sqrt
            (OPH.DiscreteCoulombGreen.realSeamEnergy (curlR i) / 2) : ℝ) : ℂ) ^ 2 *
          (Real.sin (modeAngle h (curlLamR i)) : ℂ) ^ 2 *
          (starRingEnd ℂ) (assembledCoordinate h curlLamR x i) *
          assembledCoordinate h curlLamR x i := by ring_nf
    _ = _ := by rw [hsqrt]

/-- Equivalent real-energy form of `energyCoefficientOuter_trace`: the matrix
trace is the complex cast of the committed assembled field energy. -/
theorem energyCoefficientOuter_trace_eq_energy (h : ℝ) (hh : h ≠ 0)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4)
    (x : Fin 19 → Fin 2 → ℝ) :
    (coefficientOuter (energySchrodingerVector h x)).trace =
      (assembledEnergy h curlLamR curlR x : ℂ) := by
  rw [energyCoefficientOuter_trace h hh h4,
    assembledHermitian_self h curlLamR curlR hh (curlLamR_window h hh h4)]

/-- The explicit diagonal phase matrix is exactly the Stone propagator of the
coefficient Hamiltonian. -/
theorem curlStonePropagator_eq_diagonal (h t : ℝ) :
    stonePropagator (curlHamiltonian h) t = Matrix.diagonal (curlPhase h t) := by
  unfold stonePropagator curlHamiltonian curlPhase curlFrequency
  have harg :
      t • ((-Complex.I) • Matrix.diagonal
          (fun i : Fin 19 ↦ ((modeAngle h (curlLamR i) / h : ℝ) : ℂ))) =
        Matrix.diagonal (fun i : Fin 19 ↦
          -Complex.I * (((modeAngle h (curlLamR i) / h) * t : ℝ) : ℂ)) := by
    ext i j
    by_cases hij : i = j
    · subst j
      simp
      ring_nf
    · simp [hij]
  rw [harg, Matrix.exp_diagonal]
  ext i j
  by_cases hij : i = j
  · subst j
    simp only [Matrix.diagonal_apply_eq]
    rw [Pi.coe_exp]
    rw [← Complex.exp_eq_exp_ℂ]
  · simp [hij]

/-- Reversing the flow parameter conjugates each diagonal phase. -/
theorem curlPhase_neg (h t : ℝ) (i : Fin 19) :
    curlPhase h (-t) i = (starRingEnd ℂ) (curlPhase h t i) := by
  unfold curlPhase
  rw [← Complex.exp_conj]
  congr 1
  simp [curlFrequency]

/-- The diagonal Stone family is unitary.  This uses the generic private
matrix-algebra theorem with the explicitly proved self-adjoint Hamiltonian. -/
theorem curlStonePropagator_unitary (h t : ℝ) :
    stonePropagator (curlHamiltonian h) t ∈
      unitary (Matrix (Fin 19) (Fin 19) ℂ) :=
  stonePropagator_mem_unitary (curlHamiltonian h)
    (curlHamiltonian_isSelfAdjoint h) t

/-- **Exact Hilbert-vector intertwiner.**  The Stone propagator acting on the
energy-normalized complex coefficient vector is exactly the nineteen-mode
curl flow in those coordinates.  This retains the vector phase; the outer
product theorem below is its private-algebra conjugation counterpart. -/
theorem energySchrodingerVector_stone_flow (h : ℝ) (hh : h ≠ 0)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) (t : ℝ)
    (x : Fin 19 → Fin 2 → ℝ) :
    (stonePropagator (curlHamiltonian h) t).mulVec
        (energySchrodingerVector h x) =
      energySchrodingerVector h (assembledFlow h curlLamR t x) := by
  ext i
  rw [curlStonePropagator_eq_diagonal, Matrix.mulVec_diagonal,
    energySchrodingerVector_flow h hh h4]

/-- **Exact rank-one matrix-flow intertwiner.**  The Stone conjugation of the
coefficient outer product is exactly the outer product of the evolved curl coefficient
state.  This identifies every entry in the rank-one image and retains pairwise
relative phases.  It loses global phase; the nineteen coordinates carry five
distinct frequencies with the proved degeneracies, not nineteen distinct
frequencies. -/
theorem coefficientOuter_stone_flow (h : ℝ) (hh : h ≠ 0)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) (t : ℝ)
    (x : Fin 19 → Fin 2 → ℝ) :
    stonePropagator (curlHamiltonian h) t * coefficientOuter (schrodingerVector h x) *
        stonePropagator (curlHamiltonian h) (-t) =
      coefficientOuter (schrodingerVector h (assembledFlow h curlLamR t x)) := by
  ext i j
  rw [curlStonePropagator_eq_diagonal, curlStonePropagator_eq_diagonal,
    Matrix.mul_diagonal, Matrix.diagonal_mul]
  unfold coefficientOuter schrodingerVector
  rw [schrodingerCoordinate_flow h hh h4 t x i,
    schrodingerCoordinate_flow h hh h4 t x j, map_mul, ← curlPhase_neg]
  ring_nf

/-- The exact Stone intertwiner also holds for the energy-normalized vector,
so its matrix trace is the energy-derived Hermitian norm, not an arbitrary
Euclidean normalization. -/
theorem energyCoefficientOuter_stone_flow (h : ℝ) (hh : h ≠ 0)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) (t : ℝ)
    (x : Fin 19 → Fin 2 → ℝ) :
    stonePropagator (curlHamiltonian h) t *
        coefficientOuter (energySchrodingerVector h x) *
        stonePropagator (curlHamiltonian h) (-t) =
      coefficientOuter
        (energySchrodingerVector h (assembledFlow h curlLamR t x)) := by
  ext i j
  rw [curlStonePropagator_eq_diagonal, curlStonePropagator_eq_diagonal,
    Matrix.mul_diagonal, Matrix.diagonal_mul]
  unfold coefficientOuter
  rw [energySchrodingerVector_flow h hh h4 t x i,
    energySchrodingerVector_flow h hh h4 t x j, map_mul, ← curlPhase_neg]
  ring_nf

/-! ## The conditional clock attachment -/

/-- Two sufficient calibration equalities for reading one curl mode with the
common-world scalar clock.  Neither equality is inferred from the record, and
discrete phase agreement does not make them necessary because rates are
defined only modulo the relevant `2π` phase aliases. -/
structure ModeClockCalibration
    (W : InstrumentedCommonWorldArchitecture) (i : Fin 19) : Prop where
  duration_eq_step : W.stepDuration = W.scaled.h
  mass_eq_frequency : W.mass = curlFrequency W.scaled.h i

/-- The two scalar calibration equations are algebraically consistent for
every selected mode in the positive strict window.  This constructs only the
positive rate and duration data; it does not construct a source or a physical
clock record carrying them. -/
theorem positive_mode_calibration_scalars_exist (h : ℝ) (hh : 0 < h)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) (i : Fin 19) :
    ∃ mass duration : ℝ,
      0 < mass ∧ 0 < duration ∧ duration = h ∧ mass = curlFrequency h i :=
  ⟨curlFrequency h i, h, curlFrequency_pos h hh h4 i, hh, rfl, rfl⟩

/-- Under the explicit duration and rate calibration, the record's clock
multiplies the selected Schrödinger coefficient by exactly the same phase as
one discrete Maxwell step count. -/
theorem stepClock_drives_mode
    (W : InstrumentedCommonWorldArchitecture) (i : Fin 19)
    (cal : ModeClockCalibration W i) (x : Fin 19 → Fin 2 → ℝ) (n : ℕ) :
    schrodingerCoordinate W.scaled.h
        (assembledFlow W.scaled.h curlLamR
          (stepTime W.scaled.h n) x) i =
      stepClock W n * schrodingerCoordinate W.scaled.h x i := by
  rw [schrodingerCoordinate_flow W.scaled.h W.scaled.h_pos.ne'
      (OPH.CommonWorldMaxwellClockJoin.joined_step_certified
        (OPH.CommonWorldMaxwellClockJoin.joinOf W))
      (stepTime W.scaled.h n) x i,
    stepClock_eq_exp]
  congr 1
  unfold curlPhase
  rw [cal.duration_eq_step, cal.mass_eq_frequency]

/-! ## Sharp limits of the attachment -/

/-- The `M₁₉(ℂ)` block constructed here is not the committed observer record's
`M₂(ℂ)` carrier by a relabelling of matrix indices.  A future physical bridge
therefore needs a specified compression/channel or lower-dimensional sector,
or an enlarged private carrier; the generic Stone theorem alone provides none
of these. -/
theorem curl_index_not_equiv_committed_private_index :
    ¬ Nonempty (Fin 19 ≃ Fin 2) := by
  rintro ⟨e⟩
  have hc := Fintype.card_congr e
  norm_num at hc

/-- The currently committed common-world witness does not carry the new
calibration: its declared clock duration is `1`, whereas its certified
Maxwell step is `1/2`. -/
theorem committedWitness_not_mode_clock_calibrated (i : Fin 19) :
    ¬ ModeClockCalibration instrumentedCommittedWitness i := by
  intro cal
  have hd := cal.duration_eq_step
  norm_num [instrumentedCommittedWitness,
    OPH.ScaledMaxwellStability.demoScaledBundle] at hd

/-- Already the eigenvalue-2 row and eigenvalue-3 row have different lattice
frequencies for every positive step in the full curl window. -/
theorem curlFrequency_eigen_two_ne_three (h : ℝ) (hh : 0 < h)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) :
    curlFrequency h (0 : Fin 19) ≠ curlFrequency h (5 : Fin 19) := by
  change modeAngle h 2 / h ≠ modeAngle h 3 / h
  intro heq
  have hangle : modeAngle h 2 = modeAngle h 3 := by
    have hmul := (div_eq_div_iff hh.ne' hh.ne').mp heq
    exact mul_right_cancel₀ hh.ne' hmul
  have hcos := congrArg Real.cos hangle
  have hw2 := curlLamR_window h hh.ne' h4 (0 : Fin 19)
  have hw3 := curlLamR_window h hh.ne' h4 (5 : Fin 19)
  change 0 < h ^ 2 * 2 ∧ h ^ 2 * 2 < 4 at hw2
  change 0 < h ^ 2 * 3 ∧ h ^ 2 * 3 < 4 at hw3
  rw [cos_modeAngle h 2 hw2.1.le hw2.2.le,
    cos_modeAngle h 3 hw3.1.le hw3.2.le] at hcos
  nlinarith [sq_pos_of_pos hh]

/-- Consequently no one identical scalar rate can equal every selected
principal-branch curl frequency.  This narrow result does not rule out one
shared time parameter, phase aliases, harmonics, nonlinear/finite-horizon
decoding, mode-resolved rates, or an operator-valued readout. -/
theorem no_single_scalar_clock_calibrates_all_curl_modes (h : ℝ) (hh : 0 < h)
    (h4 : h ^ 2 * (3 + Real.sqrt 5) < 4) :
    ¬ ∃ mass : ℝ, ∀ i : Fin 19, mass = curlFrequency h i := by
  rintro ⟨mass, hall⟩
  apply curlFrequency_eigen_two_ne_three h hh h4
  rw [← hall 0, ← hall 5]

#print axioms curlHamiltonian_isSelfAdjoint
#print axioms coefficientOuter_global_phase_invariant
#print axioms energyCoefficientOuter_trace
#print axioms energyCoefficientOuter_trace_eq_energy
#print axioms energySchrodingerVector_stone_flow
#print axioms energyCoefficientOuter_stone_flow
#print axioms stepClock_drives_mode
#print axioms curl_index_not_equiv_committed_private_index
#print axioms committedWitness_not_mode_clock_calibrated
#print axioms no_single_scalar_clock_calibrates_all_curl_modes

end

end OPH.CurlStoneClockBridge
