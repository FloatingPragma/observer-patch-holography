import EventAlgebra.OperationalPhaseInstrument

set_option autoImplicit false

namespace EventAlgebra

/-!
# Exact determination boundary of the ideal phase model (PR-04 delimitation)

This module sharpens the PR-04 delimitation carried by
`EventAlgebra.OperationalPhaseInstrument`.  It proves, on the committed
`Fin 2` record algebra, that the mathematical content of any inhabitant
of `IdealPhasePOVMCountModel` over the committed algebraic core is
exhausted by one complex off-diagonal coordinate of its model state,
and that the count literals of such an inhabitant satisfy exact
integer constraints.  Four groups of results:

* **Exact phase-sensitivity characterization.**  For a Hermitian
  effect, the phase clause of the model holds exactly when the
  imaginary part of the off-diagonal entry is nonzero
  (`phase_sensitivity_iff_offdiag_im`), with the exact separation
  identity `bornWeight rhoYPlus E - bornWeight rhoYMinus E =
  -2 (E 0 1).im`.  The diagonal and complexified-real insufficiency
  theorems of the base module are special cases, and every conjugation
  of a real effect by a real matrix (gauge images, permutations,
  reflections, and their compositions) is excluded from the phase slot
  in one statement (`phase_effect_ne_real_conjugated_web`).
* **State-coordinate pinning.**  Every inhabitant's model state has
  diagonal entries pinned to the committed run frequencies `111/179`
  and `68/179`, and its lower off-diagonal entry is the conjugate of
  the upper one (`prep_coordinates`).  One complex number is the whole
  remaining state freedom.
* **Exact determination.**  Over the committed core, the phase-count
  frequency equals `1/2 - Im(prep 0 1)` and the rotated-context
  frequency equals `315/716 - (sqrt 3 / 2) Re(prep 0 1)`; more
  generally, every count frequency of every inhabitant, in every
  instrument context, is an exact affine function of the one free
  coordinate (`frequency_affine`), with the phase and rotated-3
  closed forms as instances.  Two committed-core inhabitants have
  equal states exactly when all their count frequencies agree
  (`prep_eq_iff_frequencies_eq`).
* **Exact receipt window.**  Positivity of the pinned state bounds the
  off-diagonal coordinate, which yields a decidable integer inequality
  that the phase counts of every committed-core inhabitant satisfy:
  `32041 (a - b)^2 <= 30192 (a + b)^2`
  (`phase_counts_receipt_window`).  Both phase outcomes carry positive
  count mass (`phase_counts_both_pos`), and the minimal phase count
  mass over committed-core inhabitants is exactly two
  (`minimal_phase_count_mass`), attained by the synthetic inhabitant.

**Boundary.**  Nothing here produces counts, operations, preparations,
or receipts, and nothing here discharges register row PR-04.  The
theorems reduce the residue of PR-04 to its non-mathematical part: for
a committed-core inhabitant, the state diagonal, the affine dependence
of every count frequency on the one free coordinate, and the
admissible count window are forced by the
committed run literals `(111, 68)` and the committed effects, so what
the row keeps open is exactly count provenance (a custody condition
the register holds outside the type system) together with the
source production of the phase operation and the common preparation.
The receipt window is a necessary condition on any
producer-bound phase receipt, never a certificate of one; the
synthetic inhabitant of the base module shows the window is
inhabitable without any producer, so membership in the window carries
no provenance information.
-/

open Matrix
open OPH.QFT
open scoped ComplexOrder

noncomputable section

/-! ## Born weights on the two-by-two algebra in entry coordinates -/

/-- The Born pairing on the committed dimension-two algebra, written in
matrix entries. -/
theorem bornWeight_fin_two (rho P : Matrix (Fin 2) (Fin 2) ℂ) :
    bornWeight rho P =
      rho 0 0 * P 0 0 + rho 0 1 * P 1 0 + rho 1 0 * P 0 1 +
        rho 1 1 * P 1 1 := by
  simp only [bornWeight, Matrix.trace_fin_two, Matrix.mul_apply,
    Fin.sum_univ_two]
  ring

/-- **Exact separation identity.**  The Born-weight gap between the two
committed Pauli-Y states on any matrix is `i (E 0 1 - E 1 0)`. -/
theorem bornWeight_rhoY_sub (E : Matrix (Fin 2) (Fin 2) ℂ) :
    bornWeight rhoYPlus E - bornWeight rhoYMinus E =
      Complex.I * (E 0 1 - E 1 0) := by
  rw [bornWeight_fin_two, bornWeight_fin_two]
  simp [rhoYPlus, rhoYMinus]
  ring

/-- On a Hermitian matrix the separation identity closes to twice the
off-diagonal imaginary part. -/
theorem bornWeight_rhoY_sub_of_isHermitian {E : Matrix (Fin 2) (Fin 2) ℂ}
    (hE : E.IsHermitian) :
    bornWeight rhoYPlus E - bornWeight rhoYMinus E =
      -(2 * ((E 0 1).im : ℂ)) := by
  have h10 : E 1 0 = (starRingEnd ℂ) (E 0 1) := by
    have h := congrFun (congrFun hE.eq 1) 0
    rw [Matrix.conjTranspose_apply] at h
    simpa using h.symm
  rw [bornWeight_rhoY_sub, h10, Complex.sub_conj]
  push_cast
  linear_combination (2 * ((E 0 1).im : ℂ)) * Complex.I_sq

/-- **Exact phase-sensitivity characterization.**  A Hermitian effect
satisfies the phase clause of the ideal model exactly when the
imaginary part of its off-diagonal entry is nonzero.  The diagonal and
complexified-real insufficiency theorems of the base module are the
`im = 0` instances, and the positive direction states that no further
obstruction exists at the level of single Hermitian effects. -/
theorem phase_sensitivity_iff_offdiag_im {E : Matrix (Fin 2) (Fin 2) ℂ}
    (hE : E.IsHermitian) :
    bornWeight rhoYPlus E ≠ bornWeight rhoYMinus E ↔ (E 0 1).im ≠ 0 := by
  have hd := bornWeight_rhoY_sub_of_isHermitian hE
  constructor
  · intro hne him
    apply hne
    have hzero : bornWeight rhoYPlus E - bornWeight rhoYMinus E = 0 := by
      rw [hd, him]
      simp
    exact sub_eq_zero.mp hzero
  · intro him heq
    apply him
    rw [heq, sub_self] at hd
    have h2 : ((E 0 1).im : ℂ) = 0 := by linear_combination (1 / 2 : ℂ) * hd
    exact_mod_cast h2

/-- The phase outcome effect of every inhabitant has nonzero
off-diagonal imaginary part: the exact coordinate form of the phase
clause. -/
theorem phase_effect_offdiag_im_ne_zero (I : IdealPhasePOVMCountModel) :
    ((I.effect InstrumentContext.phase 0) 0 1).im ≠ 0 :=
  (phase_sensitivity_iff_offdiag_im
    (I.effect_posSemidef InstrumentContext.phase 0).isHermitian).mp
    I.phase_sensitive

/-- No complexified real matrix occupies the phase slot of any
inhabitant: the coordinate form of the real-closure exclusion. -/
theorem phase_effect_ne_complexified (I : IdealPhasePOVMCountModel)
    (M : Matrix (Fin 2) (Fin 2) ℝ) :
    I.effect InstrumentContext.phase 0 ≠ complexifyRealMatrix M := by
  intro h
  apply phase_effect_offdiag_im_ne_zero I
  rw [h]
  simp [complexifyRealMatrix]

/-- **Closed escape route.**  No conjugation of a committed web-context
projector by any real matrix (gauge images, permutation matrices,
reflections, and all their compositions are real) satisfies the phase
clause of any inhabitant. -/
theorem phase_effect_ne_real_conjugated_web (I : IdealPhasePOVMCountModel)
    (R : Matrix (Fin 2) (Fin 2) ℝ) (c : WebContext) :
    I.effect InstrumentContext.phase 0 ≠
      complexifyRealMatrix (R * webContextProjector c * Rᵀ) :=
  phase_effect_ne_complexified I _

/-! ## State-coordinate pinning by the committed run literals -/

/-- The lower off-diagonal entry of every inhabitant's model state is the
conjugate of the upper one. -/
theorem prep_offdiag_conj (I : IdealPhasePOVMCountModel) :
    I.prep 1 0 = (starRingEnd ℂ) (I.prep 0 1) := by
  have h := congrFun (congrFun I.prep_isState.1.isHermitian.eq 1) 0
  rw [Matrix.conjTranspose_apply] at h
  simpa using h.symm

/-- The upper diagonal entry of every inhabitant's model state is the
committed run frequency `111/179`. -/
theorem prep_00 (I : IdealPhasePOVMCountModel) :
    I.prep 0 0 = 111 / 179 := by
  have h := I.born_matches (InstrumentContext.web WebContext.diagonal)
  rw [I.web_effect WebContext.diagonal, webContextProjector_diagonal,
    I.diagonal_counts_run, binaryFrequency_diagonal_run,
    bornWeight_fin_two] at h
  simpa [complexifyRealMatrix, recordProjector] using h

/-- The lower diagonal entry of every inhabitant's model state is the
committed run frequency `68/179`. -/
theorem prep_11 (I : IdealPhasePOVMCountModel) :
    I.prep 1 1 = 68 / 179 := by
  have htr := I.prep_isState.2
  rw [Matrix.trace_fin_two, prep_00 I] at htr
  linear_combination htr

/-- **State-coordinate pinning.**  Every inhabitant's model state is the
committed diagonal matrix plus one free off-diagonal complex
coordinate.  The committed run literals `(111, 68)` pin the entire
diagonal; the single complex number `I.prep 0 1` is the whole remaining
state content of the model. -/
theorem prep_coordinates (I : IdealPhasePOVMCountModel) :
    I.prep = !![(111 / 179 : ℂ), I.prep 0 1;
      (starRingEnd ℂ) (I.prep 0 1), 68 / 179] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [prep_00 I, prep_11 I, prep_offdiag_conj I]

/-! ## Exact determination of all count frequencies over the committed core -/

/-- The committed algebraic core carries the committed phase lift. -/
theorem committed_core_lift :
    QuantumSurface.committedAlgebraicPhaseCompletion.lift = sourcePhaseLift :=
  rfl

/-- The phase-count frequency of a committed-core inhabitant, in the
commutator coordinates of the state. -/
theorem phase_frequency_commutator_form (I : IdealPhasePOVMCountModel)
    (hcore : I.core = QuantumSurface.committedAlgebraicPhaseCompletion) :
    binaryFrequency (I.counts InstrumentContext.phase) =
      1 / 2 + Complex.I / 2 * (I.prep 0 1 - I.prep 1 0) := by
  have h := I.born_matches InstrumentContext.phase
  rw [I.phase_effect, hcore, committed_core_lift,
    sourcePhaseLift_eq_rhoYPlus, bornWeight_fin_two] at h
  rw [← h, prep_00 I, prep_11 I]
  simp [rhoYPlus]
  ring

/-- **Exact phase determination.**  Over the committed core, the
phase-count frequency of every inhabitant equals `1/2` minus the
imaginary part of the free state coordinate.  The count literals of
the phase context carry exactly one real number of information, and
that number is already a coordinate of the model state. -/
theorem phase_frequency_eq (I : IdealPhasePOVMCountModel)
    (hcore : I.core = QuantumSurface.committedAlgebraicPhaseCompletion) :
    binaryFrequency (I.counts InstrumentContext.phase) =
      1 / 2 - ((I.prep 0 1).im : ℂ) := by
  rw [phase_frequency_commutator_form I hcore, prep_offdiag_conj I,
    Complex.sub_conj]
  push_cast
  linear_combination ((I.prep 0 1).im : ℂ) * Complex.I_sq

/-- The rotated-context frequency of any inhabitant, in the
anticommutator coordinates of the state. -/
theorem rotated_frequency_commutator_form (I : IdealPhasePOVMCountModel) :
    binaryFrequency
        (I.counts (InstrumentContext.web (WebContext.conjugated 3))) =
      315 / 716 - ((sqrt3 : ℝ) : ℂ) / 4 * (I.prep 0 1 + I.prep 1 0) := by
  have h := I.born_matches (InstrumentContext.web (WebContext.conjugated 3))
  rw [I.web_effect (WebContext.conjugated 3), webContextProjector_conjugated,
    conjProjector_three, bornWeight_fin_two] at h
  rw [← h, prep_00 I, prep_11 I]
  simp [complexifyRealMatrix]
  ring

/-- **Exact rotated determination.**  The rotated-context frequency of
every inhabitant equals `315/716` minus `sqrt 3 / 2` times the real
part of the free state coordinate. -/
theorem rotated_frequency_eq (I : IdealPhasePOVMCountModel) :
    binaryFrequency
        (I.counts (InstrumentContext.web (WebContext.conjugated 3))) =
      315 / 716 - ((sqrt3 : ℝ) : ℂ) / 2 * ((I.prep 0 1).re : ℂ) := by
  rw [rotated_frequency_commutator_form I, prep_offdiag_conj I,
    Complex.add_conj]
  push_cast
  ring

/-- **Every count frequency is affine in the free coordinate.**  For
every inhabitant and every instrument context there are real
coefficients `u`, `v`, `w`, determined by the entries of that
context's effect, with `binaryFrequency (I.counts c) =
u + v Re(prep 0 1) + w Im(prep 0 1)`.  The diagonal context is the
constant instance, and `phase_frequency_eq` and `rotated_frequency_eq`
are the closed-form instances for the phase and rotated-3 contexts. -/
theorem frequency_affine (I : IdealPhasePOVMCountModel)
    (c : InstrumentContext) :
    ∃ u v w : ℝ, binaryFrequency (I.counts c) =
      (u : ℂ) + (v : ℂ) * ((I.prep 0 1).re : ℂ) +
        (w : ℂ) * ((I.prep 0 1).im : ℂ) := by
  have hHerm : (I.effect c 0).IsHermitian :=
    (I.effect_posSemidef c 0).isHermitian
  have h10 : I.effect c 0 1 0 = (starRingEnd ℂ) (I.effect c 0 0 1) := by
    have h := congrFun (congrFun hHerm.eq 1) 0
    rw [Matrix.conjTranspose_apply] at h
    simpa using h.symm
  have h00 : ((I.effect c 0 0 0).re : ℂ) = I.effect c 0 0 0 := by
    apply Complex.conj_eq_iff_re.mp
    have h := congrFun (congrFun hHerm.eq 0) 0
    rw [Matrix.conjTranspose_apply] at h
    simpa using h
  have h11 : ((I.effect c 0 1 1).re : ℂ) = I.effect c 0 1 1 := by
    apply Complex.conj_eq_iff_re.mp
    have h := congrFun (congrFun hHerm.eq 1) 1
    rw [Matrix.conjTranspose_apply] at h
    simpa using h
  have hzw : I.prep 0 1 * (starRingEnd ℂ) (I.effect c 0 0 1) +
      (starRingEnd ℂ) (I.prep 0 1) * I.effect c 0 0 1 =
      2 * (((I.prep 0 1).re : ℂ) * ((I.effect c 0 0 1).re : ℂ) +
        ((I.prep 0 1).im : ℂ) * ((I.effect c 0 0 1).im : ℂ)) := by
    have h := Complex.add_conj
      (I.prep 0 1 * (starRingEnd ℂ) (I.effect c 0 0 1))
    rw [map_mul, Complex.conj_conj] at h
    rw [h, Complex.mul_re, Complex.conj_re, Complex.conj_im]
    push_cast
    ring
  refine ⟨111 / 179 * (I.effect c 0 0 0).re + 68 / 179 * (I.effect c 0 1 1).re,
    2 * (I.effect c 0 0 1).re, 2 * (I.effect c 0 0 1).im, ?_⟩
  rw [← I.born_matches c, bornWeight_fin_two, prep_00 I, prep_11 I,
    prep_offdiag_conj I, h10]
  push_cast
  linear_combination (-(111 : ℂ) / 179) * h00 + (-(68 : ℂ) / 179) * h11 + hzw

/-- Two committed-core inhabitants carry the same effect family: every
context effect is pinned by the structure fields and the shared core. -/
theorem committed_core_effect_eq {I J : IdealPhasePOVMCountModel}
    (hI : I.core = QuantumSurface.committedAlgebraicPhaseCompletion)
    (hJ : J.core = QuantumSurface.committedAlgebraicPhaseCompletion) :
    I.effect = J.effect := by
  funext c i
  have h0 : I.effect c 0 = J.effect c 0 := by
    cases c with
    | web wc => rw [I.web_effect wc, J.web_effect wc]
    | phase => rw [I.phase_effect, J.phase_effect, hI, hJ]
  fin_cases i
  · exact h0
  · show I.effect c 1 = J.effect c 1
    rw [I.effect_one c, J.effect_one c, h0]

/-- **Exact determination.**  Two committed-core inhabitants have equal
model states exactly when all their count frequencies agree.  Together
with `prep_coordinates`, `phase_frequency_eq`, and
`rotated_frequency_eq`, this reduces the entire mathematical content of
a committed-core inhabitant to the one free state coordinate; the
residue of register row PR-04 is count provenance and preparation
production, neither of which is expressible inside the model. -/
theorem prep_eq_iff_frequencies_eq {I J : IdealPhasePOVMCountModel}
    (hI : I.core = QuantumSurface.committedAlgebraicPhaseCompletion)
    (hJ : J.core = QuantumSurface.committedAlgebraicPhaseCompletion) :
    I.prep = J.prep ↔
      ∀ c : InstrumentContext,
        binaryFrequency (I.counts c) = binaryFrequency (J.counts c) := by
  have heff := committed_core_effect_eq hI hJ
  constructor
  · intro hp c
    rw [← I.born_matches c, ← J.born_matches c, hp, heff]
  · intro hf
    refine (preparation_identified_by_exact_fit_data I J.prep_isState ?_).symm
    intro c
    rw [heff, J.born_matches c]
    exact (hf c).symm

/-! ## The exact receipt window forced by positivity -/

/-- Positivity of the pinned state bounds the free coordinate: the
squared modulus of the off-diagonal entry is at most
`(111/179)(68/179) = 7548/32041`. -/
theorem prep_offdiag_normSq_le (I : IdealPhasePOVMCountModel) :
    Complex.normSq (I.prep 0 1) ≤ 7548 / 32041 := by
  have hx := I.prep_isState.1.dotProduct_mulVec_nonneg
    ![I.prep 1 1, -(I.prep 1 0)]
  have hval : star (![I.prep 1 1, -(I.prep 1 0)]) ⬝ᵥ
      (I.prep *ᵥ ![I.prep 1 1, -(I.prep 1 0)]) =
      (((68 : ℝ) / 179 *
        (7548 / 32041 - Complex.normSq (I.prep 0 1)) : ℝ) : ℂ) := by
    have h00 := prep_00 I
    have h11 := prep_11 I
    have h10 := prep_offdiag_conj I
    have hmc := Complex.mul_conj (I.prep 0 1)
    simp only [Matrix.mulVec, dotProduct, Fin.sum_univ_two, Pi.star_apply,
      Matrix.cons_val_zero, Matrix.cons_val_one,
      h00, h10, h11, Complex.star_def,
      map_div₀, map_ofNat]
    push_cast
    linear_combination (-(68 : ℂ) / 179) * hmc
  rw [hval, Complex.zero_le_real] at hx
  nlinarith [hx]

/-- **The exact integer receipt window.**  The phase counts `(a, b)` of
every committed-core inhabitant satisfy
`32041 (a - b)^2 <= 30192 (a + b)^2`.  This is a decidable necessary
condition on any producer-bound phase receipt over the committed
diagonal run; it is not a provenance certificate, since the synthetic
inhabitant meets it without any producer. -/
theorem phase_counts_receipt_window (I : IdealPhasePOVMCountModel)
    (hcore : I.core = QuantumSurface.committedAlgebraicPhaseCompletion) :
    (32041 : ℤ) * (((I.counts InstrumentContext.phase).1 : ℤ) -
        ((I.counts InstrumentContext.phase).2 : ℤ)) ^ 2 ≤
      30192 * (((I.counts InstrumentContext.phase).1 : ℤ) +
        ((I.counts InstrumentContext.phase).2 : ℤ)) ^ 2 := by
  set a := (I.counts InstrumentContext.phase).1 with ha
  set b := (I.counts InstrumentContext.phase).2 with hb
  have hpos : (0 : ℝ) < (a : ℝ) + (b : ℝ) := by
    have h := I.counts_pos InstrumentContext.phase
    rw [← ha, ← hb] at h
    exact_mod_cast h
  have hfreq := phase_frequency_eq I hcore
  simp only [binaryFrequency, ← ha, ← hb] at hfreq
  push_cast at hfreq
  have him : 1 / 2 - (a : ℝ) / ((a : ℝ) + (b : ℝ)) = (I.prep 0 1).im := by
    have h' : ((1 / 2 - (a : ℝ) / ((a : ℝ) + (b : ℝ)) : ℝ) : ℂ) =
        (((I.prep 0 1).im : ℝ) : ℂ) := by
      push_cast
      linear_combination -hfreq
    exact_mod_cast h'
  have hb2 := prep_offdiag_normSq_le I
  rw [Complex.normSq_apply] at hb2
  have him2 : (1 / 2 - (a : ℝ) / ((a : ℝ) + (b : ℝ))) ^ 2 ≤ 7548 / 32041 := by
    rw [him]
    nlinarith [mul_self_nonneg (I.prep 0 1).re]
  have hne : ((a : ℝ) + (b : ℝ)) ≠ 0 := ne_of_gt hpos
  have hstep : 1 / 2 - (a : ℝ) / ((a : ℝ) + (b : ℝ)) =
      ((b : ℝ) - (a : ℝ)) / (2 * ((a : ℝ) + (b : ℝ))) := by
    field_simp
    ring
  rw [hstep] at him2
  have h4 : (0 : ℝ) < (2 * ((a : ℝ) + (b : ℝ))) ^ 2 := by positivity
  rw [div_pow, div_le_iff₀ h4] at him2
  have hreal : (32041 : ℝ) * ((a : ℝ) - (b : ℝ)) ^ 2 ≤
      30192 * ((a : ℝ) + (b : ℝ)) ^ 2 := by nlinarith [him2]
  exact_mod_cast hreal

/-- Both phase outcomes of every committed-core inhabitant carry
positive count mass: the receipt window excludes every deterministic
phase count pair. -/
theorem phase_counts_both_pos (I : IdealPhasePOVMCountModel)
    (hcore : I.core = QuantumSurface.committedAlgebraicPhaseCompletion) :
    0 < (I.counts InstrumentContext.phase).1 ∧
      0 < (I.counts InstrumentContext.phase).2 := by
  have hw := phase_counts_receipt_window I hcore
  have hpos := I.counts_pos InstrumentContext.phase
  constructor
  · by_contra h
    have ha0 : (I.counts InstrumentContext.phase).1 = 0 := by omega
    have hbpos : 0 < (I.counts InstrumentContext.phase).2 := by omega
    have hB1 : (1 : ℤ) ≤ ((I.counts InstrumentContext.phase).2 : ℤ) := by
      exact_mod_cast hbpos
    rw [ha0] at hw
    push_cast at hw
    have hB2 : (1 : ℤ) ≤ ((I.counts InstrumentContext.phase).2 : ℤ) ^ 2 := by
      nlinarith [hB1]
    nlinarith [hw, hB2]
  · by_contra h
    have hb0 : (I.counts InstrumentContext.phase).2 = 0 := by omega
    have hapos : 0 < (I.counts InstrumentContext.phase).1 := by omega
    have hA1 : (1 : ℤ) ≤ ((I.counts InstrumentContext.phase).1 : ℤ) := by
      exact_mod_cast hapos
    rw [hb0] at hw
    push_cast at hw
    have hA2 : (1 : ℤ) ≤ ((I.counts InstrumentContext.phase).1 : ℤ) ^ 2 := by
      nlinarith [hA1]
    nlinarith [hw, hA2]

/-- The phase count mass of every committed-core inhabitant is at least
two. -/
theorem phase_count_mass_ge_two (I : IdealPhasePOVMCountModel)
    (hcore : I.core = QuantumSurface.committedAlgebraicPhaseCompletion) :
    2 ≤ (I.counts InstrumentContext.phase).1 +
      (I.counts InstrumentContext.phase).2 := by
  obtain ⟨h1, h2⟩ := phase_counts_both_pos I hcore
  omega

/-- The set of phase count masses realized by committed-core
inhabitants of the ideal model. -/
def committedCorePhaseMass : Set ℕ :=
  {n | ∃ I : IdealPhasePOVMCountModel,
    I.core = QuantumSurface.committedAlgebraicPhaseCompletion ∧
    (I.counts InstrumentContext.phase).1 +
      (I.counts InstrumentContext.phase).2 = n}

/-- **Exact minimal witness.**  The minimal phase count mass over
committed-core inhabitants is exactly two: the window forces both
outcomes to occur, and the synthetic inhabitant of the base module
attains mass two.  The lower bound is a constraint on every
producer receipt; the attainment witness carries no provenance. -/
theorem minimal_phase_count_mass : IsLeast committedCorePhaseMass 2 := by
  constructor
  · exact ⟨modelOfIdealPhaseFitData syntheticIdealPhaseFitData, rfl, rfl⟩
  · rintro n ⟨I, hcore, hn⟩
    rw [← hn]
    exact phase_count_mass_ge_two I hcore

end

-- Axiom audit: each must report only a subset of
-- `[propext, Classical.choice, Quot.sound]`.  No `native_decide` is used.
#print axioms bornWeight_fin_two
#print axioms bornWeight_rhoY_sub
#print axioms bornWeight_rhoY_sub_of_isHermitian
#print axioms phase_sensitivity_iff_offdiag_im
#print axioms phase_effect_offdiag_im_ne_zero
#print axioms phase_effect_ne_complexified
#print axioms phase_effect_ne_real_conjugated_web
#print axioms prep_offdiag_conj
#print axioms prep_00
#print axioms prep_11
#print axioms prep_coordinates
#print axioms committed_core_lift
#print axioms phase_frequency_commutator_form
#print axioms phase_frequency_eq
#print axioms rotated_frequency_commutator_form
#print axioms rotated_frequency_eq
#print axioms frequency_affine
#print axioms committed_core_effect_eq
#print axioms prep_eq_iff_frequencies_eq
#print axioms prep_offdiag_normSq_le
#print axioms phase_counts_receipt_window
#print axioms phase_counts_both_pos
#print axioms phase_count_mass_ge_two
#print axioms minimal_phase_count_mass

end EventAlgebra
