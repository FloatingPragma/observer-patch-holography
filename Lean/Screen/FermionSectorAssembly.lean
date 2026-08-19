import AssembledActionComposition

/-!
# Fermion-sector completion of the assembled structural action (V3, issue #735, OL-G9)

One typed fermion sector on the committed matter grammar, and one
extended assembled action that sums it with the committed assembled
structural action of `AssembledActionComposition` on one extended
field tuple.  The committed assembled action covers three sectors
(gauge kinetic, scalar, realified Yukawa line); this module adds the
fermion sector as a fourth summand and proves the composition,
restriction, invariance, census, and mutation receipts on one typed
premise bundle.

The fermion carrier is one complex amplitude per row of the committed
ten-component exterior table (register row PR-59): the committed
selection mask of the bundled base restricts the sector functional to
the committed mask rows, and each row carries the committed integer
charge column (`q = 6Y`) and the committed chirality grading through
the index dictionary of `MatterGrammarIndexBridge`.  The fermion
bilinear is the row-pair form `conj (psi i) * psi j`: the
kinetic/covariant pairing shape at the exact structural level of the
committed assembled action.  The committed action's sectors are
algebraic functionals with no derivative and no spacetime carrier
(register row PR-54 open), so the fermion sector is one too: the
sector functional is the masked diagonal sum of declared coefficients
against the squared row amplitudes, the finite stand-in for a kinetic
term, exactly as the committed gauge sector pairs declared weights
against squared strengths.

Proved on the extended bundle `FermionExtendedPremiseData` (the
committed `AssembledActionPremiseData` plus one declared coefficient
per row):

1. the abelian phase action through the committed charge column
   leaves the row-pair form invariant exactly when the committed row
   charges balance (`rowPairForm_invariant_iff_charge_balance`): the
   forward direction extracts divisibility of the charge difference
   from two explicit unit phases (a fourth and a third root of
   unity), and a kernel check on the committed charge column closes
   the gap; on the committed mask rows the balanced pairs are exactly
   the diagonal ones (`committed_mask_invariant_pairs_diagonal`), so
   the charge column forces the kinetic-shaped pairing to be
   diagonal;
2. the extended action is one typed sum: the committed assembled
   action plus the fermion sector functional at the bundled mask,
   with exact restriction receipts: on the committed sub-tuple (zero
   fermion slot) and at zero fermion coefficients the extended action
   is exactly the committed assembled action, on the fermion slice it
   is exactly the fermion sector functional, and the three committed
   sector restrictions persist through the extension;
3. the Yukawa line of the committed electroweak composition is
   related to the fermion tuple by a typed identification of the
   shared mask rows: `matterSlotsOfFermion` embeds the amplitudes at
   the committed doublet and singlet rows (`yukawaRows` of the
   bundled mask) into the committed matter slots, the embedded
   evaluation carries the exact committed mass bilinear coefficient,
   and the embedding intertwines the fermion phase action with the
   committed matter-slot phase action at the committed row charges;
4. the extended parameter census: bundles agreeing on the committed
   declared coefficients (the two gauge weights, the quartic, the
   scale, the two couplings, the Yukawa coefficient), the committed
   mask, and the named fermion coefficients have equal extended
   actions (`extended_parameter_census`);
5. an inhabitant extending the committed inhabitant, with one exact
   witness evaluation;
6. mutation receipts: a charge-imbalanced control pair of committed
   mask rows breaks the invariance receipt at an exhibited unit
   phase and tuple (`charge_imbalanced_control`), and the
   zero-coefficient restriction pins the extension against silently
   changing the committed sectors;
7. one composed receipt (`fermionSectorAssembly_receipt`) with
   `#print axioms` guards.

Consumed register rows.  The base field carries the committed rows of
the assembled action: PR-09, PR-10, PR-11, PR-59 through the
committed base, declared parameter rows PR-48, PR-49, PR-50, PR-54.
The charge and chirality columns are the committed table content of
row PR-59 through the committed index dictionary.  The fermion
coefficients are declared here on the boundary of rows PR-47 and
PR-54: the covariant kinetic term those rows name is absent, and the
coefficients are its finite stand-in weights.

Boundary and nonclaims.  The fermion coefficients are declared and
select nothing.  There is no source action: register row PR-54 stays
open, and no derivative, spacetime carrier, or field-strength
structure appears beyond what the committed assembled action itself
carries; every sector, the fermion sector included, is an algebraic
functional on a finite tuple.  There is no physical Spin or chirality
attachment: register row PR-47 stays open, and the chirality grading
is table provenance through the committed dictionary.  There is no
quantization, no measured coupling, no three-family structure (row
PR-36 untouched), no Yukawa matrix (row PR-50 untouched beyond the
committed one-line coefficient), and no global-form attachment (rows
PR-35, PR-46 untouched).  OL-G9 stays partial, with the
covered-sector list extended exactly by the fermion sector.

Falsifier.  The composition fails if the row-pair invariance holds at
a committed-mask pair with unbalanced committed charges or fails at a
balanced pair, if any restriction receipt differs from the committed
functional, if the shared-row identification misses the committed
mass bilinear coefficient, if the census equality fails on agreeing
parameter lists, or if the exhibited control pair fails to break
invariance.
-/

namespace OPH.FermionSectorAssembly

open OPH.AssembledActionComposition (FieldTuple SectorStrength
  AssembledActionPremiseData assembledAction baseTuple scaleTuple
  committedAssembledActionPremiseData assembledAction_base
  assembledAction_explicit restrict_gauge restrict_scalar
  restrict_yukawa restrict_yukawa_mass_coefficient
  assembled_abelian_invariance committed_charge_abelian_invariance
  committed_inhabitant_value gaugeKineticEnergy)
open OPH.ElectroweakBreakingComposition (Doublet vacuum yukawaLine
  epsPair scalarChargeSix wDirection normSqSum)
open OPH.GaugeKineticInvariantForms (gNormalForm)
open OPH.ExteriorSelection (charge mem evenMask oddMask isDoublet)

/-! ## The fermion carrier on the committed matter grammar -/

/-- The fermion field tuple: one complex amplitude per row of the
committed ten-component exterior table (register row PR-59).  The
committed selection mask of the bundled base restricts the sector
functional to the committed mask rows below. -/
abbrev FermionTuple := Fin 10 → ℂ

/-- The abelian phase action through the committed charge column:
each row amplitude is rescaled by the phase raised to the committed
integer charge (`q = 6Y`) of its row. -/
noncomputable def chargeScale (p : ℂ) (ψ : FermionTuple) :
    FermionTuple :=
  fun i => p ^ charge i * ψ i

/-- The fermion row-pair form `conj (psi i) * psi j`: the
kinetic/covariant pairing shape on one pair of committed rows, at the
algebraic level of the committed assembled action (no derivative,
register row PR-54 open). -/
def rowPairForm (i j : Fin 10) (ψ : FermionTuple) : ℂ :=
  starRingEnd ℂ (ψ i) * ψ j

/-- The diagonal row-pair form is the squared row amplitude. -/
theorem rowPairForm_diag (i : Fin 10) (ψ : FermionTuple) :
    rowPairForm i i ψ = (Complex.normSq (ψ i) : ℂ) := by
  rw [rowPairForm, Complex.normSq_eq_conj_mul_self]

/-- Transformation law of the row-pair form under the charge phase
action. -/
theorem rowPairForm_chargeScale (p : ℂ) (i j : Fin 10)
    (ψ : FermionTuple) :
    rowPairForm i j (chargeScale p ψ)
      = starRingEnd ℂ p ^ charge i * p ^ charge j
          * rowPairForm i j ψ := by
  simp only [rowPairForm, chargeScale, map_mul, map_zpow₀]
  ring

/-- At a unit phase the transformation factor is one power of the
phase at the charge difference: the exact balance bookkeeping of the
conjugated pairing. -/
theorem rowPairForm_chargeScale_unit (p : ℂ)
    (hp : Complex.normSq p = 1) (i j : Fin 10) (ψ : FermionTuple) :
    rowPairForm i j (chargeScale p ψ)
      = p ^ (charge j - charge i) * rowPairForm i j ψ := by
  have hp0 : p ≠ 0 := by
    intro h
    rw [h] at hp
    simp at hp
  have hconj : starRingEnd ℂ p = p⁻¹ := by
    have hmc : p * starRingEnd ℂ p = 1 := by
      rw [Complex.mul_conj, hp, Complex.ofReal_one]
    exact (inv_eq_of_mul_eq_one_right hmc).symm
  rw [rowPairForm_chargeScale, hconj, inv_zpow', ← zpow_add₀ hp0,
    neg_add_eq_sub]

/-- Balanced committed charges give exact invariance at every unit
phase: the positive half of the invariance receipt, proved from the
committed charge column. -/
theorem rowPairForm_invariant_of_charge_balance {i j : Fin 10}
    (hbal : charge i = charge j) (p : ℂ)
    (hp : Complex.normSq p = 1) (ψ : FermionTuple) :
    rowPairForm i j (chargeScale p ψ) = rowPairForm i j ψ := by
  rw [rowPairForm_chargeScale_unit p hp, hbal, sub_self, zpow_zero,
    one_mul]

/-! ## Root-of-unity divisibility toolkit for the converse -/

/-- Divisibility extraction for integer powers: a nonzero base with
exact period `n` on natural powers divides every integer exponent
sent to one. -/
theorem zpow_dvd_of_prim_period {p : ℂ} (hp0 : p ≠ 0) {n : ℕ}
    (hn : 0 < n) (hpn : p ^ n = 1)
    (hmin : ∀ r : ℕ, 0 < r → r < n → p ^ r ≠ 1)
    {d : ℤ} (hd : p ^ d = 1) : (n : ℤ) ∣ d := by
  have hnz : (n : ℤ) ≠ 0 := by exact_mod_cast hn.ne'
  have hsplit : (n : ℤ) * (d / (n : ℤ)) + d % (n : ℤ) = d :=
    Int.mul_ediv_add_emod d (n : ℤ)
  have hrem : p ^ (d % (n : ℤ)) = 1 := by
    have hpow : p ^ ((n : ℤ) * (d / (n : ℤ)) + d % (n : ℤ)) = 1 := by
      rw [hsplit]
      exact hd
    rw [zpow_add₀ hp0, zpow_mul, zpow_natCast, hpn, one_zpow,
      one_mul] at hpow
    exact hpow
  have h0 : 0 ≤ d % (n : ℤ) := Int.emod_nonneg d hnz
  have hlt : d % (n : ℤ) < (n : ℤ) :=
    Int.emod_lt_of_pos d (by exact_mod_cast hn)
  have hcast : (((d % (n : ℤ)).toNat : ℤ)) = d % (n : ℤ) :=
    Int.toNat_of_nonneg h0
  have hrn : (d % (n : ℤ)).toNat < n := by omega
  have hr1 : p ^ (d % (n : ℤ)).toNat = 1 := by
    rw [← zpow_natCast, hcast]
    exact hrem
  have hrz : (d % (n : ℤ)).toNat = 0 := by
    by_contra hne
    exact hmin _ (Nat.pos_of_ne_zero hne) hrn hr1
  have hz : d % (n : ℤ) = 0 := by omega
  exact Int.dvd_of_emod_eq_zero hz

/-- Integer powers of `I` equal to one have exponent divisible by
four. -/
theorem four_dvd_of_I_zpow_eq_one {d : ℤ} (hd : Complex.I ^ d = 1) :
    (4 : ℤ) ∣ d := by
  have hnI : (-Complex.I : ℂ) ≠ 1 := by
    simp [Complex.ext_iff]
  refine zpow_dvd_of_prim_period Complex.I_ne_zero (by norm_num)
    ?_ ?_ hd
  · exact Complex.I_pow_four
  · intro r hr0 hr4
    interval_cases r
    · rw [pow_one]
      simp [Complex.ext_iff]
    · rw [Complex.I_sq]
      norm_num
    · rw [pow_succ, Complex.I_sq, neg_one_mul]
      exact hnI

/-- One explicit third root of unity: the unit phase at angle two
thirds of a turn. -/
noncomputable def omegaPhase : ℂ :=
  ((-1 / 2 : ℝ) : ℂ) + ((Real.sqrt 3 / 2 : ℝ) : ℂ) * Complex.I

theorem omegaPhase_normSq : Complex.normSq omegaPhase = 1 := by
  have h3 : Real.sqrt 3 ^ 2 = 3 := Real.sq_sqrt (by norm_num)
  rw [omegaPhase, Complex.normSq_add_mul_I]
  nlinarith [h3]

theorem omegaPhase_ne_zero : omegaPhase ≠ 0 := by
  intro h
  have hn := omegaPhase_normSq
  rw [h] at hn
  simp at hn

theorem omegaPhase_quadratic :
    omegaPhase ^ 2 + omegaPhase + 1 = 0 := by
  have h3 : ((Real.sqrt 3 : ℝ) : ℂ) ^ 2 = 3 := by
    rw [← Complex.ofReal_pow, Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 3)]
    norm_num
  rw [omegaPhase]
  push_cast
  linear_combination (3 / 4 : ℂ) * Complex.I_sq
    + (Complex.I ^ 2 / 4) * h3

theorem omegaPhase_pow_three : omegaPhase ^ 3 = 1 := by
  linear_combination (omegaPhase - 1) * omegaPhase_quadratic

theorem omegaPhase_ne_one : omegaPhase ≠ 1 := by
  intro h
  have him := congrArg Complex.im h
  simp only [omegaPhase, Complex.add_im, Complex.ofReal_im,
    Complex.mul_im, Complex.ofReal_re, Complex.I_im, Complex.I_re,
    Complex.one_im, mul_one, mul_zero, add_zero, zero_add] at him
  have h3 : (0 : ℝ) < Real.sqrt 3 := Real.sqrt_pos.mpr (by norm_num)
  linarith

theorem omegaPhase_sq_ne_one : omegaPhase ^ 2 ≠ 1 := by
  intro h
  have hq := omegaPhase_quadratic
  rw [h] at hq
  have hval : omegaPhase = -2 := by linear_combination hq
  have hn := omegaPhase_normSq
  rw [hval] at hn
  norm_num [Complex.normSq_apply] at hn

/-- Integer powers of the third root of unity equal to one have
exponent divisible by three. -/
theorem three_dvd_of_omegaPhase_zpow_eq_one {d : ℤ}
    (hd : omegaPhase ^ d = 1) : (3 : ℤ) ∣ d := by
  refine zpow_dvd_of_prim_period omegaPhase_ne_zero (by norm_num)
    omegaPhase_pow_three ?_ hd
  intro r hr0 hr3
  interval_cases r
  · rw [pow_one]
    exact omegaPhase_ne_one
  · exact omegaPhase_sq_ne_one

/-! ## Kernel checks on the committed charge column -/

/-- On either committed parity mask, a charge difference divisible by
twelve is zero: the committed charge column has all its mask-internal
differences inside the window that the two exhibited roots of unity
resolve.  Kernel-checked from the committed table. -/
theorem mask_twelve_dvd_charge_balance :
    ∀ i j : Fin 10,
      (mem evenMask i = true → mem evenMask j = true →
        (charge j - charge i) % 12 = 0 → charge i = charge j)
      ∧ (mem oddMask i = true → mem oddMask j = true →
        (charge j - charge i) % 12 = 0 → charge i = charge j) := by
  decide

/-- The committed charge column separates the rows of either
committed parity mask: equal charges force equal rows.
Kernel-checked from the committed table. -/
theorem mask_charge_injective :
    ∀ i j : Fin 10,
      (mem evenMask i = true → mem evenMask j = true →
        charge i = charge j → i = j)
      ∧ (mem oddMask i = true → mem oddMask j = true →
        charge i = charge j → i = j) := by
  decide

/-! ## The invariance receipt: balance exactly, diagonal exactly -/

/-- **Charge-balance invariance, exact.**  On committed mask rows,
the row-pair form is invariant under the abelian phase action at
every unit phase exactly when the committed row charges balance.  The
forward direction instantiates the hypothesis at `I` and at the third
root of unity, extracts divisibility of the charge difference by
twelve, and closes with the kernel check on the committed charge
column. -/
theorem rowPairForm_invariant_iff_charge_balance (m : ℕ)
    (hm : m = evenMask ∨ m = oddMask) {i j : Fin 10}
    (hi : mem m i = true) (hj : mem m j = true) :
    (∀ p : ℂ, Complex.normSq p = 1 →
        ∀ ψ : FermionTuple,
          rowPairForm i j (chargeScale p ψ) = rowPairForm i j ψ)
      ↔ charge i = charge j := by
  constructor
  · intro hinv
    have hone : rowPairForm i j (fun _ => 1) = 1 := by
      simp [rowPairForm]
    have hphase : ∀ p : ℂ, Complex.normSq p = 1 →
        p ^ (charge j - charge i) = 1 := by
      intro p hp
      have hcall := hinv p hp (fun _ => 1)
      rw [rowPairForm_chargeScale_unit p hp, hone, mul_one] at hcall
      exact hcall
    have h4 : (4 : ℤ) ∣ (charge j - charge i) :=
      four_dvd_of_I_zpow_eq_one
        (hphase Complex.I Complex.normSq_I)
    have h3 : (3 : ℤ) ∣ (charge j - charge i) :=
      three_dvd_of_omegaPhase_zpow_eq_one
        (hphase omegaPhase omegaPhase_normSq)
    have h12 : (charge j - charge i) % 12 = 0 := by omega
    rcases hm with rfl | rfl
    · exact (mask_twelve_dvd_charge_balance i j).1 hi hj h12
    · exact (mask_twelve_dvd_charge_balance i j).2 hi hj h12
  · intro hbal p hp ψ
    exact rowPairForm_invariant_of_charge_balance hbal p hp ψ

/-- **The committed charge column forces the diagonal.**  On the
committed mask rows the invariant row-pair forms are exactly the
diagonal ones: the kinetic/covariant pairing shape compatible with
the abelian phase action at the committed charges pairs each selected
row with itself. -/
theorem committed_mask_invariant_pairs_diagonal (m : ℕ)
    (hm : m = evenMask ∨ m = oddMask) {i j : Fin 10}
    (hi : mem m i = true) (hj : mem m j = true) :
    (∀ p : ℂ, Complex.normSq p = 1 →
        ∀ ψ : FermionTuple,
          rowPairForm i j (chargeScale p ψ) = rowPairForm i j ψ)
      ↔ i = j := by
  rw [rowPairForm_invariant_iff_charge_balance m hm hi hj]
  constructor
  · intro hc
    rcases hm with rfl | rfl
    · exact (mask_charge_injective i j).1 hi hj hc
    · exact (mask_charge_injective i j).2 hi hj hc
  · rintro rfl
    rfl

/-- **The charge-imbalanced control.**  Rows `3` and `8` of the
committed even mask carry distinct committed charges, and the
exhibited unit phase `-1` with the all-ones tuple changes the value
of their row-pair form: the invariance receipt is load-bearing on the
committed charge balance. -/
theorem charge_imbalanced_control :
    mem evenMask 3 = true ∧ mem evenMask 8 = true
      ∧ charge (3 : Fin 10) ≠ charge 8
      ∧ Complex.normSq (-1 : ℂ) = 1
      ∧ rowPairForm 3 8 (chargeScale (-1) fun _ => (1 : ℂ))
          ≠ rowPairForm 3 8 (fun _ => (1 : ℂ)) := by
  refine ⟨by decide, by decide, by decide, by simp, ?_⟩
  have hc3 : charge (3 : Fin 10) = 1 := by decide
  have hc8 : charge (8 : Fin 10) = 2 := by decide
  have hval : rowPairForm 3 8 (chargeScale (-1) fun _ => (1 : ℂ))
      = -1 := by
    simp only [rowPairForm, chargeScale, hc3, hc8]
    norm_num
  have hone : rowPairForm 3 8 (fun _ => (1 : ℂ)) = 1 := by
    simp [rowPairForm]
  rw [hval, hone]
  norm_num

/-! ## The fermion sector functional -/

/-- The fermion sector functional: the masked diagonal sum of
declared per-row coefficients against the squared row amplitudes.
The diagonal shape is the one the committed charge column permits
(`committed_mask_invariant_pairs_diagonal`), and the pairing level
matches the committed gauge kinetic sector: declared weights against
squared magnitudes, no derivative (register row PR-54 open). -/
def fermionSector (κ : Fin 10 → ℝ) (m : ℕ) (ψ : FermionTuple) : ℝ :=
  ∑ i : Fin 10, if mem m i then κ i * Complex.normSq (ψ i) else 0

/-- The sector functional is the masked sum of the diagonal row-pair
forms with the declared coefficients: the typed link between the
sector and the fermion bilinear. -/
theorem fermionSector_diagonal_display (κ : Fin 10 → ℝ) (m : ℕ)
    (ψ : FermionTuple) :
    fermionSector κ m ψ
      = ∑ i : Fin 10,
          if mem m i then κ i * (rowPairForm i i ψ).re else 0 := by
  unfold fermionSector
  refine Finset.sum_congr rfl fun i _ => ?_
  by_cases h : mem m i = true
  · simp only [h, if_true]
    rw [rowPairForm_diag]
    simp
  · simp [h]

/-- The fermion sector is exactly invariant under the abelian phase
action through the committed charge column, at every unit phase: each
diagonal term is balanced by the committed bookkeeping. -/
theorem fermionSector_chargeScale (κ : Fin 10 → ℝ) (m : ℕ) (p : ℂ)
    (hp : Complex.normSq p = 1) (ψ : FermionTuple) :
    fermionSector κ m (chargeScale p ψ) = fermionSector κ m ψ := by
  unfold fermionSector
  refine Finset.sum_congr rfl fun i _ => ?_
  have hnorm : Complex.normSq (p ^ charge i * ψ i)
      = Complex.normSq (ψ i) := by
    rw [map_mul, map_zpow₀, hp, one_zpow, one_mul]
  simp only [chargeScale]
  rw [hnorm]

theorem fermionSector_zero_tuple (κ : Fin 10 → ℝ) (m : ℕ) :
    fermionSector κ m 0 = 0 := by
  unfold fermionSector
  simp

theorem fermionSector_zero_coefficients (m : ℕ) (ψ : FermionTuple) :
    fermionSector (fun _ => 0) m ψ = 0 := by
  unfold fermionSector
  simp

/-- Nonnegativity of the fermion sector under declared nonnegative
coefficients: the kinetic-positivity clause at the committed level. -/
theorem fermionSector_nonneg (κ : Fin 10 → ℝ) (hκ : ∀ i, 0 ≤ κ i)
    (m : ℕ) (ψ : FermionTuple) : 0 ≤ fermionSector κ m ψ := by
  unfold fermionSector
  refine Finset.sum_nonneg fun i _ => ?_
  by_cases h : mem m i = true
  · simp only [h, if_true]
    exact mul_nonneg (hκ i) (Complex.normSq_nonneg _)
  · simp [h]

/-! ## The extended premise bundle and the extended field tuple -/

/-- The extended premise bundle: the committed assembled-action
bundle plus one declared fermion coefficient per committed row.  The
coefficients are declared on the boundary of register rows PR-47 and
PR-54: the covariant kinetic term those rows name is absent, and the
coefficients are its finite stand-in weights.  No coefficient is
derived or selected. -/
structure FermionExtendedPremiseData where
  /-- The committed assembled-action premise bundle. -/
  base : AssembledActionPremiseData
  /-- The declared fermion coefficients, one per committed row. -/
  kappa : Fin 10 → ℝ

/-- The extended field tuple: the committed field tuple plus the
fermion tuple. -/
structure ExtendedFieldTuple where
  /-- The committed boson field tuple of the assembled action. -/
  boson : FieldTuple
  /-- The fermion field tuple on the committed rows. -/
  fermion : FermionTuple

/-- **The extended assembled action.**  The committed assembled
structural action plus the fermion sector at the bundled committed
mask, as one real function of one extended field tuple. -/
noncomputable def extendedAction (E : FermionExtendedPremiseData)
    (Ψ : ExtendedFieldTuple) : ℝ :=
  assembledAction E.base Ψ.boson
    + fermionSector E.kappa E.base.ew.base.selectionMask.val Ψ.fermion

/-- The extended base tuple: the committed base tuple with the zero
fermion tuple. -/
noncomputable def extendedBaseTuple (E : FermionExtendedPremiseData) :
    ExtendedFieldTuple :=
  ⟨baseTuple E.base, 0⟩

/-! ## Restriction receipts: the extension pins the committed
sectors -/

/-- **Committed sub-tuple restriction.**  On the zero fermion slot
the extended action is exactly the committed assembled action: the
extension changes no committed sector. -/
theorem extended_restricts_on_committed_subtuple
    (E : FermionExtendedPremiseData) (Φ : FieldTuple) :
    extendedAction E ⟨Φ, 0⟩ = assembledAction E.base Φ := by
  unfold extendedAction
  rw [fermionSector_zero_tuple, add_zero]

/-- **Zero-coefficient restriction.**  At the zero fermion
coefficient set the extended action is exactly the committed
assembled action on every extended tuple: the pinning receipt against
silent changes to committed sectors. -/
theorem extended_restricts_zero_coefficients
    (E : FermionExtendedPremiseData) (hκ : ∀ i, E.kappa i = 0)
    (Ψ : ExtendedFieldTuple) :
    extendedAction E Ψ = assembledAction E.base Ψ.boson := by
  unfold extendedAction
  have hfun : E.kappa = fun _ => 0 := funext hκ
  rw [hfun, fermionSector_zero_coefficients, add_zero]

/-- **Fermion restriction.**  On the fermion slice through the
committed base tuple the extended action is exactly the fermion
sector functional. -/
theorem restrict_fermion (E : FermionExtendedPremiseData)
    (ψ : FermionTuple) :
    extendedAction E ⟨baseTuple E.base, ψ⟩
      = fermionSector E.kappa E.base.ew.base.selectionMask.val ψ := by
  unfold extendedAction
  rw [assembledAction_base, zero_add]

/-- The committed gauge restriction persists through the
extension. -/
theorem extended_restrict_gauge (E : FermionExtendedPremiseData)
    (s : SectorStrength) :
    extendedAction E ⟨{ baseTuple E.base with strength := s }, 0⟩
      = gaugeKineticEnergy
          (gNormalForm E.base.gaugeFreePlus E.base.gaugeFive) s := by
  rw [extended_restricts_on_committed_subtuple]
  exact restrict_gauge E.base s

/-- The committed scalar restriction persists through the
extension. -/
theorem extended_restrict_scalar (E : FermionExtendedPremiseData)
    (h : Doublet) (A : OPH.ElectroweakBreakingComposition.EWParam) :
    extendedAction E
        ⟨{ baseTuple E.base with scalar := h, direction := A }, 0⟩
      = OPH.ElectroweakBreakingComposition.potential E.base.ew.lam
          E.base.ew.vev h
        + OPH.ElectroweakBreakingComposition.gaugeMassForm
            E.base.ew.gW E.base.ew.gY E.base.ew.vev A := by
  rw [extended_restricts_on_committed_subtuple]
  exact restrict_scalar E.base h A

/-- The committed Yukawa restriction persists through the
extension. -/
theorem extended_restrict_yukawa (E : FermionExtendedPremiseData)
    (ψ : Doublet) (χ : ℂ) :
    extendedAction E
        ⟨{ baseTuple E.base with
            matterDoublet := ψ, matterSinglet := χ }, 0⟩
      = (yukawaLine E.base.ew.yuk (vacuum E.base.ew.vev) ψ χ).re := by
  rw [extended_restricts_on_committed_subtuple]
  exact restrict_yukawa E.base ψ χ

/-- The extended action vanishes at the extended base tuple. -/
theorem extendedAction_base (E : FermionExtendedPremiseData) :
    extendedAction E (extendedBaseTuple E) = 0 := by
  unfold extendedBaseTuple
  rw [extended_restricts_on_committed_subtuple, assembledAction_base]

/-! ## The explicit display and the extended parameter census -/

/-- **The explicit extended sum.**  The extended action written out
with every coefficient a declared bundle field: the committed
polynomial-plus-realified-line display of the assembled action, plus
the masked fermion sum with the declared per-row coefficients. -/
theorem extendedAction_explicit (E : FermionExtendedPremiseData)
    (Ψ : ExtendedFieldTuple) :
    extendedAction E Ψ
      = E.base.gaugeFreePlus * Ψ.boson.strength.threePlus ^ 2
        + Real.sqrt 5 * E.base.gaugeFive
            * Ψ.boson.strength.threeMinus ^ 2
        + E.base.gaugeFive * Ψ.boson.strength.five ^ 2
        + E.base.ew.lam
            * (normSqSum Ψ.boson.scalar - E.base.ew.vev ^ 2) ^ 2
        + E.base.ew.vev ^ 2 / 4
            * (E.base.ew.gW ^ 2
                * (Ψ.boson.direction.a1 ^ 2
                  + Ψ.boson.direction.a2 ^ 2)
              + (E.base.ew.gY * Ψ.boson.direction.b
                  - E.base.ew.gW * Ψ.boson.direction.a3) ^ 2)
        + ((E.base.ew.yuk : ℂ)
            * epsPair Ψ.boson.matterDoublet Ψ.boson.scalar
            * Ψ.boson.matterSinglet).re
        + ∑ i : Fin 10,
            (if mem E.base.ew.base.selectionMask.val i then
              E.kappa i * Complex.normSq (Ψ.fermion i) else 0) := by
  unfold extendedAction fermionSector
  rw [assembledAction_explicit]

/-- **The extended parameter census.**  Two extended bundles agreeing
on the committed declared coefficient list of the assembled action
(the two gauge weights, the quartic, the scale, the two couplings,
the Yukawa coefficient), on the committed selection mask, and on the
named fermion coefficients have equal extended actions: the declared
parameter list of the extension is exactly the committed census plus
the fermion coefficients. -/
theorem extended_parameter_census (E E' : FermionExtendedPremiseData)
    (hplus : E.base.gaugeFreePlus = E'.base.gaugeFreePlus)
    (hfive : E.base.gaugeFive = E'.base.gaugeFive)
    (hlam : E.base.ew.lam = E'.base.ew.lam)
    (hvev : E.base.ew.vev = E'.base.ew.vev)
    (hgW : E.base.ew.gW = E'.base.ew.gW)
    (hgY : E.base.ew.gY = E'.base.ew.gY)
    (hyuk : E.base.ew.yuk = E'.base.ew.yuk)
    (hmask : E.base.ew.base.selectionMask.val
        = E'.base.ew.base.selectionMask.val)
    (hkappa : ∀ i : Fin 10, E.kappa i = E'.kappa i) :
    ∀ Ψ : ExtendedFieldTuple,
      extendedAction E Ψ = extendedAction E' Ψ := by
  intro Ψ
  rw [extendedAction_explicit, extendedAction_explicit, hplus, hfive,
    hlam, hvev, hgW, hgY, hyuk, hmask]
  congr 1
  exact Finset.sum_congr rfl fun i _ => by rw [hkappa]

/-! ## The typed Yukawa identification on the shared mask rows -/

/-- The committed Yukawa rows of a parity mask: the colorless doublet
row and the balancing singlet row that the committed
electroweak-breaking receipt selects in each sector. -/
def yukawaRows (m : ℕ) : Fin 10 × Fin 10 :=
  if m = evenMask then (9, 4) else (1, 7)

/-- The typed embedding of the fermion tuple into the committed
matter slots at the shared mask rows: the doublet slot carries the
doublet-row amplitude in its mass-bilinear component, the singlet
slot carries the singlet-row amplitude. -/
def matterSlotsOfFermion (m : ℕ) (ψ : FermionTuple) : Doublet × ℂ :=
  (![ψ (yukawaRows m).1, 0], ψ (yukawaRows m).2)

/-- The committed-table receipt of the shared rows: on the bundled
committed mask, the two Yukawa rows are selected, carry the doublet
and singlet indicators, share the committed Weyl-chirality column
through the committed index dictionary, and their committed charges
balance the scalar row charge or its negative. -/
theorem yukawaRows_committed_receipt (E : FermionExtendedPremiseData) :
    mem E.base.ew.base.selectionMask.val
        (yukawaRows E.base.ew.base.selectionMask.val).1 = true
      ∧ mem E.base.ew.base.selectionMask.val
          (yukawaRows E.base.ew.base.selectionMask.val).2 = true
      ∧ isDoublet (yukawaRows E.base.ew.base.selectionMask.val).1
          = true
      ∧ isDoublet (yukawaRows E.base.ew.base.selectionMask.val).2
          = false
      ∧ OPH.BaryonDimensionSix.left
            (OPH.MatterGrammarIndexBridge.fieldIndex
              (yukawaRows E.base.ew.base.selectionMask.val).1)
          = OPH.BaryonDimensionSix.left
            (OPH.MatterGrammarIndexBridge.fieldIndex
              (yukawaRows E.base.ew.base.selectionMask.val).2)
      ∧ (charge (yukawaRows E.base.ew.base.selectionMask.val).1
            + charge (yukawaRows E.base.ew.base.selectionMask.val).2
            + scalarChargeSix = 0
          ∨ charge (yukawaRows E.base.ew.base.selectionMask.val).1
            + charge (yukawaRows E.base.ew.base.selectionMask.val).2
            - scalarChargeSix = 0) := by
  rcases E.base.ew.base.matter_selection_is_parity_sector with h | h <;>
    rw [h] <;> decide

/-- **The typed Yukawa identification.**  For every extended bundle
and fermion tuple, evaluating the extended action on the embedded
matter slots (fermion slot zero) yields exactly the committed mass
bilinear with coefficient `yuk * vev` in the shared-row amplitudes:
the Yukawa line of the committed electroweak composition and the
fermion tuple are related by the typed embedding, with no prose
identification. -/
theorem fermion_matter_identification (E : FermionExtendedPremiseData)
    (ψ : FermionTuple) :
    extendedAction E
        ⟨{ baseTuple E.base with
            matterDoublet :=
              (matterSlotsOfFermion
                E.base.ew.base.selectionMask.val ψ).1,
            matterSinglet :=
              (matterSlotsOfFermion
                E.base.ew.base.selectionMask.val ψ).2 }, 0⟩
      = E.base.ew.yuk * E.base.ew.vev
          * (ψ (yukawaRows E.base.ew.base.selectionMask.val).1
              * ψ (yukawaRows E.base.ew.base.selectionMask.val).2).re := by
  rw [extended_restricts_on_committed_subtuple,
    restrict_yukawa_mass_coefficient]
  simp [matterSlotsOfFermion]

/-- The embedding intertwines the fermion charge phase action with
the committed matter-slot phase action at the committed row charges:
the doublet slot scales by the phase at the doublet-row charge, the
singlet slot by the phase at the singlet-row charge, exactly the
powers of the committed charge-balance invariance receipt. -/
theorem matterSlots_phase_equivariant (m : ℕ) (p : ℂ)
    (ψ : FermionTuple) :
    ((matterSlotsOfFermion m (chargeScale p ψ)).1
        = fun k => p ^ charge (yukawaRows m).1
            * (matterSlotsOfFermion m ψ).1 k)
      ∧ (matterSlotsOfFermion m (chargeScale p ψ)).2
        = p ^ charge (yukawaRows m).2
            * (matterSlotsOfFermion m ψ).2 := by
  constructor
  · funext k
    fin_cases k <;> simp [matterSlotsOfFermion, chargeScale]
  · rfl

/-! ## Invariance of the extended sum -/

/-- The extended phase rescaling: the committed abelian scaling on
the boson slots and the charge phase action on the fermion slot. -/
noncomputable def extendedScaleTuple (u v w p : ℂ)
    (Ψ : ExtendedFieldTuple) : ExtendedFieldTuple :=
  ⟨scaleTuple u v w Ψ.boson, chargeScale p Ψ.fermion⟩

/-- **Exact abelian invariance of the extended sum.**  For unit
scalars with product one on the boson slots and any unit phase on the
fermion slot, the extended action is exactly invariant: the committed
invariance of the assembled action and the diagonal invariance of the
fermion sector compose. -/
theorem extended_abelian_invariance (E : FermionExtendedPremiseData)
    (u v w p : ℂ) (hu : Complex.normSq u = 1) (hprod : u * v * w = 1)
    (hp : Complex.normSq p = 1) (Ψ : ExtendedFieldTuple) :
    extendedAction E (extendedScaleTuple u v w p Ψ)
      = extendedAction E Ψ := by
  unfold extendedAction extendedScaleTuple
  rw [assembled_abelian_invariance E.base u v w hu hprod Ψ.boson,
    fermionSector_chargeScale E.kappa _ p hp Ψ.fermion]

/-- **Committed-charge abelian invariance of the extended sum.**  The
bundled selection mask supplies the committed doublet and singlet
rows whose integer charges balance the scalar row charge or its
negative; at that committed charge assignment every unit phase,
acting through the integer charges on the boson matter slots and
through the committed charge column on the fermion slot, leaves the
extended action exactly fixed. -/
theorem extended_committed_charge_abelian_invariance
    (E : FermionExtendedPremiseData) :
    ∃ (cH : ℤ) (dRow sRow : Fin 10),
      (cH = scalarChargeSix ∨ cH = -scalarChargeSix)
        ∧ mem E.base.ew.base.selectionMask.val dRow = true
        ∧ mem E.base.ew.base.selectionMask.val sRow = true
        ∧ cH + charge dRow + charge sRow = 0
        ∧ ∀ p : ℂ, Complex.normSq p = 1 →
            ∀ Ψ : ExtendedFieldTuple,
              extendedAction E
                  (extendedScaleTuple (p ^ cH) (p ^ charge dRow)
                    (p ^ charge sRow) p Ψ)
                = extendedAction E Ψ := by
  obtain ⟨cH, dRow, sRow, hcH, hd, hs, hbal, hinv⟩ :=
    committed_charge_abelian_invariance E.base
  refine ⟨cH, dRow, sRow, hcH, hd, hs, hbal, fun p hp Ψ => ?_⟩
  unfold extendedAction extendedScaleTuple
  rw [hinv p hp Ψ.boson,
    fermionSector_chargeScale E.kappa _ p hp Ψ.fermion]

/-! ## Nonvacuity: the extended committed inhabitant -/

/-- The extended committed inhabitant: the committed assembled-action
instance with every fermion coefficient set to one.  The instance is
a nonvacuity control; the coefficient values select nothing. -/
noncomputable def committedFermionExtendedPremiseData :
    FermionExtendedPremiseData where
  base := committedAssembledActionPremiseData
  kappa := fun _ => 1

/-- The extended premise bundle is satisfiable. -/
theorem fermionExtendedPremiseData_nonvacuous :
    Nonempty FermionExtendedPremiseData :=
  ⟨committedFermionExtendedPremiseData⟩

/-- Exact witness evaluation at the extended committed inhabitant:
the committed boson witness value `9/4` plus one unit per selected
row of the committed even mask, five in total.  A witness value of
the declared parameters, not a prediction. -/
theorem committed_extended_inhabitant_value :
    extendedAction committedFermionExtendedPremiseData
        ⟨⟨⟨1, 0, 0⟩, vacuum 1, wDirection, ![1, 0], 1⟩, fun _ => 1⟩
      = 29 / 4 := by
  have hferm : fermionSector (fun _ => (1 : ℝ)) evenMask
      (fun _ => (1 : ℂ)) = 5 := by
    unfold fermionSector
    simp
    norm_cast
  show assembledAction committedAssembledActionPremiseData
        ⟨⟨1, 0, 0⟩, vacuum 1, wDirection, ![1, 0], 1⟩
      + fermionSector (fun _ => (1 : ℝ)) evenMask (fun _ => (1 : ℂ))
      = 29 / 4
  rw [committed_inhabitant_value, hferm]
  norm_num

/-! ## The chirality-grading re-export at the fermion bundle -/

/-- Every fermion row selected by the bundled committed mask carries
the committed Weyl-chirality column through the committed index
dictionary: the image census row is undotted left-Weyl exactly when
the bundled selection is the even-parity sector.  Table provenance
under register row PR-59; a physical chirality claim stays on
register row PR-47. -/
theorem fermion_rows_weyl_graded (E : FermionExtendedPremiseData) :
    ∀ i : Fin 10,
      mem E.base.ew.base.selectionMask.val i = true →
        (OPH.BaryonDimensionSix.left
            (OPH.MatterGrammarIndexBridge.fieldIndex i) = 1
          ↔ E.base.ew.base.selectionMask.val = evenMask) :=
  OPH.MatterGrammarIndexBridge.selected_components_weyl_parity
    E.base.ew.base

/-! ## The composed receipt -/

/-- **The composed fermion-sector assembly receipt (issue #735,
OL-G9 partial).**  For every extended premise bundle `E`: the
extended action is the committed assembled action plus the fermion
sector at the bundled committed mask; it restricts exactly to the
committed assembled action on the committed sub-tuple and at zero
fermion coefficients, exactly to the fermion sector on the fermion
slice, and vanishes at the extended base tuple; the fermion sector is
exactly invariant under the charge phase action at every unit phase,
the invariant row-pair forms on the bundled mask rows are exactly the
diagonal ones, and the extended action is exactly invariant at the
committed charge balance of the bundled selection mask; the shared
Yukawa rows carry the committed table receipt and the typed embedding
evaluates to the exact committed mass bilinear; the exhibited
charge-imbalanced control pair breaks the invariance receipt.  Every
bundle-dependent clause consumes the bundle.  Consumed premise rows
through the bundle: PR-09, PR-10, PR-11, PR-59 (committed base),
PR-48, PR-49, PR-50, PR-54 (declared parameter fields; the rows stay
open); the fermion coefficients are declared on the boundary of rows
PR-47 and PR-54.  Rows PR-35, PR-36, PR-46, PR-47 are untouched. -/
theorem fermionSectorAssembly_receipt (E : FermionExtendedPremiseData) :
    (∀ Ψ : ExtendedFieldTuple,
        extendedAction E Ψ
          = assembledAction E.base Ψ.boson
            + fermionSector E.kappa
                E.base.ew.base.selectionMask.val Ψ.fermion)
    ∧ ((∀ Φ : FieldTuple,
          extendedAction E ⟨Φ, 0⟩ = assembledAction E.base Φ)
      ∧ ((∀ i : Fin 10, E.kappa i = 0) →
          ∀ Ψ : ExtendedFieldTuple,
            extendedAction E Ψ = assembledAction E.base Ψ.boson)
      ∧ (∀ ψ : FermionTuple,
          extendedAction E ⟨baseTuple E.base, ψ⟩
            = fermionSector E.kappa
                E.base.ew.base.selectionMask.val ψ)
      ∧ extendedAction E (extendedBaseTuple E) = 0)
    ∧ ((∀ p : ℂ, Complex.normSq p = 1 →
          ∀ ψ : FermionTuple,
            fermionSector E.kappa E.base.ew.base.selectionMask.val
                (chargeScale p ψ)
              = fermionSector E.kappa
                  E.base.ew.base.selectionMask.val ψ)
      ∧ (∀ i j : Fin 10,
          mem E.base.ew.base.selectionMask.val i = true →
          mem E.base.ew.base.selectionMask.val j = true →
            ((∀ p : ℂ, Complex.normSq p = 1 →
                ∀ ψ : FermionTuple,
                  rowPairForm i j (chargeScale p ψ)
                    = rowPairForm i j ψ)
              ↔ i = j))
      ∧ (∃ (cH : ℤ) (dRow sRow : Fin 10),
          (cH = scalarChargeSix ∨ cH = -scalarChargeSix)
            ∧ mem E.base.ew.base.selectionMask.val dRow = true
            ∧ mem E.base.ew.base.selectionMask.val sRow = true
            ∧ cH + charge dRow + charge sRow = 0
            ∧ ∀ p : ℂ, Complex.normSq p = 1 →
                ∀ Ψ : ExtendedFieldTuple,
                  extendedAction E
                      (extendedScaleTuple (p ^ cH) (p ^ charge dRow)
                        (p ^ charge sRow) p Ψ)
                    = extendedAction E Ψ))
    ∧ ((mem E.base.ew.base.selectionMask.val
            (yukawaRows E.base.ew.base.selectionMask.val).1 = true
        ∧ mem E.base.ew.base.selectionMask.val
            (yukawaRows E.base.ew.base.selectionMask.val).2 = true
        ∧ (charge (yukawaRows E.base.ew.base.selectionMask.val).1
              + charge (yukawaRows E.base.ew.base.selectionMask.val).2
              + scalarChargeSix = 0
            ∨ charge (yukawaRows E.base.ew.base.selectionMask.val).1
              + charge (yukawaRows E.base.ew.base.selectionMask.val).2
              - scalarChargeSix = 0))
      ∧ (∀ ψ : FermionTuple,
          extendedAction E
              ⟨{ baseTuple E.base with
                  matterDoublet :=
                    (matterSlotsOfFermion
                      E.base.ew.base.selectionMask.val ψ).1,
                  matterSinglet :=
                    (matterSlotsOfFermion
                      E.base.ew.base.selectionMask.val ψ).2 }, 0⟩
            = E.base.ew.yuk * E.base.ew.vev
                * (ψ (yukawaRows E.base.ew.base.selectionMask.val).1
                    * ψ (yukawaRows
                        E.base.ew.base.selectionMask.val).2).re))
    ∧ (charge (3 : Fin 10) ≠ charge 8
      ∧ rowPairForm 3 8 (chargeScale (-1) fun _ => (1 : ℂ))
          ≠ rowPairForm 3 8 (fun _ => (1 : ℂ))) := by
  obtain ⟨hd, hs, -, -, -, hbal⟩ := yukawaRows_committed_receipt E
  refine ⟨fun Ψ => rfl,
    ⟨fun Φ => extended_restricts_on_committed_subtuple E Φ,
      fun hκ Ψ => extended_restricts_zero_coefficients E hκ Ψ,
      fun ψ => restrict_fermion E ψ,
      extendedAction_base E⟩,
    ⟨fun p hp ψ => fermionSector_chargeScale E.kappa _ p hp ψ,
      fun i j hi hj =>
        committed_mask_invariant_pairs_diagonal
          E.base.ew.base.selectionMask.val
          E.base.ew.base.matter_selection_is_parity_sector hi hj,
      extended_committed_charge_abelian_invariance E⟩,
    ⟨⟨hd, hs, hbal⟩,
      fun ψ => fermion_matter_identification E ψ⟩,
    ⟨charge_imbalanced_control.2.2.1,
      charge_imbalanced_control.2.2.2.2⟩⟩

end OPH.FermionSectorAssembly

/- Axiom audit: standard axioms only; no native_decide. -/

#print axioms OPH.FermionSectorAssembly.rowPairForm_diag
#print axioms OPH.FermionSectorAssembly.rowPairForm_chargeScale_unit
#print axioms OPH.FermionSectorAssembly.rowPairForm_invariant_of_charge_balance
#print axioms OPH.FermionSectorAssembly.zpow_dvd_of_prim_period
#print axioms OPH.FermionSectorAssembly.four_dvd_of_I_zpow_eq_one
#print axioms OPH.FermionSectorAssembly.omegaPhase_pow_three
#print axioms OPH.FermionSectorAssembly.three_dvd_of_omegaPhase_zpow_eq_one
#print axioms OPH.FermionSectorAssembly.mask_twelve_dvd_charge_balance
#print axioms OPH.FermionSectorAssembly.mask_charge_injective
#print axioms OPH.FermionSectorAssembly.rowPairForm_invariant_iff_charge_balance
#print axioms OPH.FermionSectorAssembly.committed_mask_invariant_pairs_diagonal
#print axioms OPH.FermionSectorAssembly.charge_imbalanced_control
#print axioms OPH.FermionSectorAssembly.fermionSector_diagonal_display
#print axioms OPH.FermionSectorAssembly.fermionSector_chargeScale
#print axioms OPH.FermionSectorAssembly.fermionSector_nonneg
#print axioms OPH.FermionSectorAssembly.extended_restricts_on_committed_subtuple
#print axioms OPH.FermionSectorAssembly.extended_restricts_zero_coefficients
#print axioms OPH.FermionSectorAssembly.restrict_fermion
#print axioms OPH.FermionSectorAssembly.extended_restrict_gauge
#print axioms OPH.FermionSectorAssembly.extended_restrict_scalar
#print axioms OPH.FermionSectorAssembly.extended_restrict_yukawa
#print axioms OPH.FermionSectorAssembly.extendedAction_base
#print axioms OPH.FermionSectorAssembly.extendedAction_explicit
#print axioms OPH.FermionSectorAssembly.extended_parameter_census
#print axioms OPH.FermionSectorAssembly.yukawaRows_committed_receipt
#print axioms OPH.FermionSectorAssembly.fermion_matter_identification
#print axioms OPH.FermionSectorAssembly.matterSlots_phase_equivariant
#print axioms OPH.FermionSectorAssembly.extended_abelian_invariance
#print axioms OPH.FermionSectorAssembly.extended_committed_charge_abelian_invariance
#print axioms OPH.FermionSectorAssembly.fermionExtendedPremiseData_nonvacuous
#print axioms OPH.FermionSectorAssembly.committed_extended_inhabitant_value
#print axioms OPH.FermionSectorAssembly.fermion_rows_weyl_graded
#print axioms OPH.FermionSectorAssembly.fermionSectorAssembly_receipt
