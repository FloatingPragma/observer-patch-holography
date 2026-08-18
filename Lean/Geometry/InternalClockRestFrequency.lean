import Geometry.MassShellKinematics

/-!
# C1: internal clock and rest frequency on the declared Lorentz module

On the declared Hermitian Lorentz module, this file attaches the declared
plane-wave phase form to the kinematic mass-shell packet of
`Geometry/MassShellKinematics.lean` and proves the internal-clock identities.
The phase of the declared-mass momentum along the declared frame worldline is
exactly `mass * τ`, so the internal rest rotation rate is the declared mass
parameter in the module's dimensionless convention, in every oriented Lorentz
chart.  The declared internal clock satisfies the exact rotation law
`I * derivative = mass * clock`, its periods are exactly the integer multiples
of `2π / mass` for nonzero mass, and on the declared plane-wave form the
spatial phase gradient equals the momentum coordinate while the temporal
phase rate equals the scalar (energy) coordinate, with the committed shell
identity tying the two readouts.  On the zero-mass branch the internal clock
is constant, every real number is a period, so no fundamental period exists,
and no oriented Lorentz chart brings a future null vector to vanishing
spatial coordinates.  A completeness lemma shows every vector on the positive
shell arises as a declared-mass momentum of some frame.

The evolution is a declared plane-wave form on the declared module, and the
mass parameter is a declared hypothesis of each statement.  Nothing here
constructs a physical clock, a measured frequency, an electron or a Compton
value, a source-produced state, or an energy readout in laboratory units; the
clock and energy calibration anchors and the observer-to-physical-spacetime
attachment remain open register premises.  This is a formal
precursor for the internal-particle-clock observation row, not a physics
attainment.
-/

namespace OPH.C1Lorentz

noncomputable section

open OPH.C2Soldering

/-! ## Declared phase form, worldline, and internal clock -/

/-- The declared plane-wave phase form: the Lorentz pairing of a momentum
with an event coordinate, `E t - p · x` in the `(+---)` convention. -/
def planeWavePhase (p x : Herm2) : ℝ := lorentzB p x

/-- The declared inertial worldline through the origin with unit-timelike
velocity given by a frame, parameterized by the dimensionless proper-time
coordinate. -/
def frameWorldline (frame : FrameHyperboloid) (τ : ℝ) : Herm2 :=
  τ • (frame : Herm2)

/-- The declared plane wave `exp(i (p · x - E t)) = exp(-i lorentzB p x)`. -/
def planeWave (p x : Herm2) : ℂ :=
  Complex.exp (-Complex.I * (planeWavePhase p x : ℝ))

/-- The declared internal clock: the plane wave of the declared-mass momentum
read along the frame's own worldline. -/
def internalClock (mass : ℝ) (frame : FrameHyperboloid) (τ : ℝ) : ℂ :=
  planeWave (fourMomentum mass frame) (frameWorldline frame τ)

/-! ## The internal rest phase and its rate -/

/-- The phase of the declared-mass momentum along its own frame worldline is
exactly `mass * τ`: the internal rest rotation advances at the declared mass
parameter in the module's dimensionless convention. -/
theorem planeWavePhase_frameWorldline (mass : ℝ) (frame : FrameHyperboloid)
    (τ : ℝ) :
    planeWavePhase (fourMomentum mass frame) (frameWorldline frame τ) =
      mass * τ := by
  simp only [planeWavePhase, fourMomentum, frameWorldline, lorentzB_smul_left,
    lorentzB_smul_right, lorentzB_self, frame.2.1]
  ring

/-- Affine real functions have the expected derivative. -/
theorem hasDerivAt_const_add_mul (c0 c1 s : ℝ) :
    HasDerivAt (fun t : ℝ => c0 + t * c1) c1 s := by
  simpa using ((hasDerivAt_id s).mul_const c1).const_add c0

/-- The internal rest rotation rate is exactly the declared mass parameter:
the proper-time derivative of the rest phase equals `mass` at every
proper-time value. -/
theorem hasDerivAt_restPhase (mass : ℝ) (frame : FrameHyperboloid) (τ : ℝ) :
    HasDerivAt
      (fun σ : ℝ =>
        planeWavePhase (fourMomentum mass frame) (frameWorldline frame σ))
      mass τ := by
  have hfun :
      (fun σ : ℝ =>
        planeWavePhase (fourMomentum mass frame) (frameWorldline frame σ)) =
      fun σ : ℝ => 0 + σ * mass := by
    funext σ
    rw [planeWavePhase_frameWorldline]
    ring
  rw [hfun]
  exact hasDerivAt_const_add_mul 0 mass τ

/-- The internal clock is the exact unit-circle rotation at the declared mass
rate. -/
theorem internalClock_eq_exp (mass : ℝ) (frame : FrameHyperboloid) (τ : ℝ) :
    internalClock mass frame τ =
      Complex.exp (-Complex.I * ((mass * τ : ℝ) : ℂ)) := by
  rw [internalClock, planeWave, planeWavePhase_frameWorldline]

/-- The internal clock satisfies the exact rotation law of a phase evolving
at the declared mass rate. -/
theorem internalClock_hasDerivAt (mass : ℝ) (frame : FrameHyperboloid)
    (τ : ℝ) :
    HasDerivAt (internalClock mass frame)
      (-Complex.I * mass * internalClock mass frame τ) τ := by
  have hfun : internalClock mass frame =
      fun σ : ℝ => Complex.exp (-Complex.I * (mass : ℂ) * (σ : ℂ)) := by
    funext σ
    rw [internalClock_eq_exp]
    push_cast
    ring_nf
  have hlin : HasDerivAt (fun σ : ℝ => -Complex.I * (mass : ℂ) * (σ : ℂ))
      (-Complex.I * (mass : ℂ)) τ := by
    simpa using
      (Complex.ofRealCLM.hasDerivAt (x := τ)).const_mul
        (-Complex.I * (mass : ℂ))
  have hexp := hlin.cexp
  rw [hfun]
  simpa [mul_comm, mul_left_comm] using hexp

/-- The Zitterbewegung-shaped internal-clock identity: `I` times the clock
derivative equals the declared mass parameter times the clock. -/
theorem internalClock_rotation_law (mass : ℝ) (frame : FrameHyperboloid)
    (τ : ℝ) :
    Complex.I * deriv (internalClock mass frame) τ =
      (mass : ℂ) * internalClock mass frame τ := by
  rw [(internalClock_hasDerivAt mass frame τ).deriv]
  have hI : Complex.I * -Complex.I = 1 := by
    rw [mul_neg, Complex.I_mul_I, neg_neg]
  calc
    Complex.I * (-Complex.I * (mass : ℂ) * internalClock mass frame τ) =
        Complex.I * -Complex.I * ((mass : ℂ) * internalClock mass frame τ) :=
      by ring
    _ = (mass : ℂ) * internalClock mass frame τ := by rw [hI, one_mul]

/-! ## Oriented-chart invariance -/

/-- Oriented Lorentz maps transport frame worldlines to the transported
frame's worldlines. -/
theorem frameWorldline_oriented (L : OrientedLorentzEquiv)
    (frame : FrameHyperboloid) (τ : ℝ) :
    L (frameWorldline frame τ) = frameWorldline (L.mapFrame frame) τ := by
  simp [frameWorldline]

/-- The plane-wave phase is invariant under joint oriented transport of
momentum and event coordinate. -/
theorem planeWavePhase_oriented (L : OrientedLorentzEquiv) (p x : Herm2) :
    planeWavePhase (L p) (L x) = planeWavePhase p x :=
  L.map_lorentzB p x

/-- The rest phase read in any oriented chart is the same function
`mass * τ`. -/
theorem restPhase_oriented (L : OrientedLorentzEquiv) (mass : ℝ)
    (frame : FrameHyperboloid) (τ : ℝ) :
    planeWavePhase (L (fourMomentum mass frame))
      (L (frameWorldline frame τ)) = mass * τ := by
  rw [planeWavePhase_oriented, planeWavePhase_frameWorldline]

/-- Oriented-chart invariance of the internal clock: the transported frame
carries the identical clock function, hence the identical internal rest
frequency. -/
theorem internalClock_oriented (L : OrientedLorentzEquiv) (mass : ℝ)
    (frame : FrameHyperboloid) :
    internalClock mass (L.mapFrame frame) = internalClock mass frame := by
  funext τ
  rw [internalClock_eq_exp, internalClock_eq_exp]

/-- The internal rest rotation rate equals the declared mass parameter in
every oriented chart. -/
theorem hasDerivAt_restPhase_oriented (L : OrientedLorentzEquiv) (mass : ℝ)
    (frame : FrameHyperboloid) (τ : ℝ) :
    HasDerivAt
      (fun σ : ℝ =>
        planeWavePhase (fourMomentum mass (L.mapFrame frame))
          (frameWorldline (L.mapFrame frame) σ))
      mass τ :=
  hasDerivAt_restPhase mass (L.mapFrame frame) τ

/-! ## Exact period structure and frequency rigidity -/

/-- `exp(-I a) = 1` exactly on the integer multiples of `2π`. -/
theorem exp_neg_I_mul_real_eq_one_iff (a : ℝ) :
    Complex.exp (-Complex.I * (a : ℂ)) = 1 ↔
      ∃ n : ℤ, a = n * (2 * Real.pi) := by
  rw [Complex.exp_eq_one_iff]
  constructor
  · rintro ⟨n, hn⟩
    refine ⟨-n, ?_⟩
    have hmul := congrArg (fun z : ℂ => Complex.I * z) hn
    have hleft : Complex.I * (-Complex.I * (a : ℂ)) = (a : ℂ) := by
      rw [show Complex.I * (-Complex.I * (a : ℂ)) =
          -(Complex.I * Complex.I) * (a : ℂ) by ring, Complex.I_mul_I]
      ring
    have hright : Complex.I * ((n : ℂ) * (2 * (Real.pi : ℂ) * Complex.I)) =
        -((n : ℂ) * (2 * (Real.pi : ℂ))) := by
      rw [show Complex.I * ((n : ℂ) * (2 * (Real.pi : ℂ) * Complex.I)) =
          (n : ℂ) * (2 * (Real.pi : ℂ)) * (Complex.I * Complex.I) by ring,
        Complex.I_mul_I]
      ring
    simp only at hmul
    rw [hleft, hright] at hmul
    have hreal : a = -((n : ℝ) * (2 * Real.pi)) := by
      exact_mod_cast hmul
    push_cast
    linarith
  · rintro ⟨n, hn⟩
    refine ⟨-n, ?_⟩
    rw [hn]
    push_cast
    ring

/-- Exact period classification: a real number is a period of the internal
clock exactly when the declared mass times it is an integer multiple of
`2π`.  This is the frequency-rigidity statement for the declared clock. -/
theorem internalClock_periodic_iff (mass : ℝ) (frame : FrameHyperboloid)
    (T : ℝ) :
    Function.Periodic (internalClock mass frame) T ↔
      ∃ n : ℤ, mass * T = n * (2 * Real.pi) := by
  constructor
  · intro hper
    have h0 := hper 0
    rw [internalClock_eq_exp, internalClock_eq_exp] at h0
    have h1 : Complex.exp (-Complex.I * ((mass * T : ℝ) : ℂ)) = 1 := by
      have hz : ((mass * (0 + T) : ℝ) : ℂ) = ((mass * T : ℝ) : ℂ) := by
        push_cast
        ring
      have hz0 : ((mass * 0 : ℝ) : ℂ) = 0 := by
        push_cast
        ring
      rw [hz, hz0] at h0
      simpa using h0
    exact (exp_neg_I_mul_real_eq_one_iff (mass * T)).mp h1
  · rintro ⟨n, hn⟩
    intro τ
    rw [internalClock_eq_exp, internalClock_eq_exp]
    have hsplit : ((mass * (τ + T) : ℝ) : ℂ) =
        ((mass * τ : ℝ) : ℂ) + ((mass * T : ℝ) : ℂ) := by
      push_cast
      ring
    rw [hsplit, mul_add, Complex.exp_add,
      (exp_neg_I_mul_real_eq_one_iff (mass * T)).mpr ⟨n, hn⟩, mul_one]

/-- For nonzero declared mass, `2π / mass` is a period of the internal
clock: the de Broglie internal period in the module's dimensionless
convention. -/
theorem internalClock_periodic {mass : ℝ} (hm : mass ≠ 0)
    (frame : FrameHyperboloid) :
    Function.Periodic (internalClock mass frame) (2 * Real.pi / mass) := by
  refine (internalClock_periodic_iff mass frame _).mpr ⟨1, ?_⟩
  push_cast
  field_simp

/-- For nonzero declared mass the period set is exactly the integer lattice
generated by `2π / mass`. -/
theorem internalClock_period_rigidity {mass : ℝ} (hm : mass ≠ 0)
    (frame : FrameHyperboloid) (T : ℝ) :
    Function.Periodic (internalClock mass frame) T ↔
      ∃ n : ℤ, T = n * (2 * Real.pi / mass) := by
  rw [internalClock_periodic_iff]
  constructor
  · rintro ⟨n, hn⟩
    refine ⟨n, ?_⟩
    field_simp
    linarith [hn]
  · rintro ⟨n, hn⟩
    refine ⟨n, ?_⟩
    rw [hn]
    field_simp

/-! ## The de Broglie relations on the declared plane-wave form -/

/-- The de Broglie orientation of the phase, `p · x - E t`: the negative of
the invariant Lorentz phase. -/
def wavePhase (p x : Herm2) : ℝ := spatialDot p.2 x.2 - p.1 * x.1

theorem wavePhase_eq_neg_planeWavePhase (p x : Herm2) :
    wavePhase p x = -planeWavePhase p x := by
  simp only [wavePhase, planeWavePhase, lorentzB]
  ring

/-- Bilinearity of the spatial pairing along an affine spatial line. -/
theorem spatialDot_add_smul (a xs w : Spatial) (s : ℝ) :
    spatialDot a (xs + s • w) = spatialDot a xs + s * spatialDot a w := by
  simp only [spatialDot, Pi.add_apply, Pi.smul_apply, smul_eq_mul,
    Finset.mul_sum]
  rw [← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro i _
  ring

/-- The spatial pairing against a coordinate axis reads off the momentum
coordinate. -/
theorem spatialDot_single (x : Spatial) (i : Fin 3) :
    spatialDot x (Pi.single i 1) = x i := by
  simp [spatialDot, Pi.single_apply, mul_ite]

/-- The spatial phase gradient of the de Broglie phase along any spatial
direction is the spatial momentum pairing with that direction. -/
theorem hasDerivAt_wavePhase_spatial (p : Herm2) (t : ℝ) (xs w : Spatial)
    (s : ℝ) :
    HasDerivAt (fun σ : ℝ => wavePhase p (t, xs + σ • w))
      (spatialDot p.2 w) s := by
  have hfun : (fun σ : ℝ => wavePhase p (t, xs + σ • w)) =
      fun σ : ℝ => (spatialDot p.2 xs - p.1 * t) + σ * spatialDot p.2 w := by
    funext σ
    simp only [wavePhase, spatialDot_add_smul]
    ring
  rw [hfun]
  exact hasDerivAt_const_add_mul _ _ s

/-- The exact de Broglie relation, spatial half: the phase gradient along
the `i`-th spatial axis equals the `i`-th momentum coordinate. -/
theorem hasDerivAt_wavePhase_spatial_axis (p : Herm2) (t : ℝ) (xs : Spatial)
    (i : Fin 3) (s : ℝ) :
    HasDerivAt
      (fun σ : ℝ => wavePhase p (t, xs + σ • (Pi.single i 1 : Spatial)))
      (p.2 i) s := by
  have h := hasDerivAt_wavePhase_spatial p t xs (Pi.single i 1 : Spatial) s
  rwa [spatialDot_single] at h

/-- The exact de Broglie relation, temporal half: the temporal rate of the
invariant Lorentz phase equals the scalar (energy) coordinate. -/
theorem hasDerivAt_planeWavePhase_time (p : Herm2) (xs : Spatial) (t : ℝ) :
    HasDerivAt (fun σ : ℝ => planeWavePhase p (σ, xs)) p.1 t := by
  have hfun : (fun σ : ℝ => planeWavePhase p (σ, xs)) =
      fun σ : ℝ => -spatialDot p.2 xs + σ * p.1 := by
    funext σ
    simp only [planeWavePhase, lorentzB]
    ring
  rw [hfun]
  exact hasDerivAt_const_add_mul _ _ t

/-- On the de Broglie orientation the temporal phase rate is minus the
scalar coordinate: the wave rotates as `exp(-i E t)` in time. -/
theorem hasDerivAt_wavePhase_time (p : Herm2) (xs : Spatial) (t : ℝ) :
    HasDerivAt (fun σ : ℝ => wavePhase p (σ, xs)) (-p.1) t := by
  have hfun : (fun σ : ℝ => wavePhase p (σ, xs)) =
      fun σ : ℝ => spatialDot p.2 xs + σ * -p.1 := by
    funext σ
    simp only [wavePhase]
    ring
  rw [hfun]
  exact hasDerivAt_const_add_mul _ _ t

/-- The plane wave factorizes into the temporal rotation at the scalar
(energy) coordinate rate and the spatial wave at the momentum coordinates. -/
theorem planeWave_time_factorization (p : Herm2) (t : ℝ) (xs : Spatial) :
    planeWave p (t, xs) =
      Complex.exp (-Complex.I * ((p.1 * t : ℝ) : ℂ)) *
        Complex.exp (Complex.I * ((spatialDot p.2 xs : ℝ) : ℂ)) := by
  rw [planeWave, ← Complex.exp_add]
  have harg : -Complex.I * ((planeWavePhase p (t, xs) : ℝ) : ℂ) =
      -Complex.I * ((p.1 * t : ℝ) : ℂ) +
        Complex.I * ((spatialDot p.2 xs : ℝ) : ℂ) := by
    have hphase : planeWavePhase p (t, xs) =
        p.1 * t - spatialDot p.2 xs := rfl
    rw [hphase]
    push_cast
    ring
  rw [harg]

/-! ## The zero-mass branch -/

/-- The zero-mass momentum degenerates to the zero vector: the massless
branch is not reached by the massive parameterization. -/
theorem fourMomentum_zero_mass (frame : FrameHyperboloid) :
    fourMomentum 0 frame = 0 := by
  simp [fourMomentum]

/-- At zero declared mass the internal clock is constant: it never
advances. -/
theorem internalClock_zero_mass (frame : FrameHyperboloid) (τ : ℝ) :
    internalClock 0 frame τ = 1 := by
  rw [internalClock_eq_exp]
  simp

/-- At zero declared mass every real number is a period, so the internal
period degenerates: no fundamental period exists. -/
theorem internalClock_zero_mass_every_period (frame : FrameHyperboloid)
    (T : ℝ) :
    Function.Periodic (internalClock 0 frame) T := by
  intro τ
  rw [internalClock_zero_mass, internalClock_zero_mass]

/-- A future null vector has nonvanishing spatial coordinates. -/
theorem futureNull_spatial_ne_zero (v : FutureNullVector) : v.1.2 ≠ 0 := by
  intro hzero
  have hq := v.2.1
  have ht := v.2.2
  have hs : spatialNormSq (0 : Spatial) = 0 := by
    simp [spatialNormSq, pow_two]
  rw [lorentzQ, hzero, hs] at hq
  nlinarith

/-- No rest chart exists on the zero-mass branch: no oriented Lorentz map
brings a future null vector to vanishing spatial coordinates. -/
theorem futureNull_no_rest_chart (L : OrientedLorentzEquiv)
    (v : FutureNullVector) : (L v.1).2 ≠ 0 :=
  futureNull_spatial_ne_zero (L.mapFutureNullVector v)

/-! ## Shell completeness -/

/-- Every vector on the positive-scalar shell of a positive declared mass is
the declared-mass momentum of exactly the rescaled frame: the massive
parameterization covers the shell. -/
theorem exists_frame_fourMomentum_eq {mass : ℝ} (hm : 0 < mass) {v : Herm2}
    (hshell : lorentzQ v = mass ^ 2) (hpos : 0 < v.1) :
    ∃ u : FrameHyperboloid, fourMomentum mass u = v := by
  have hmne : mass ≠ 0 := ne_of_gt hm
  refine ⟨⟨mass⁻¹ • v, ?_, ?_⟩, ?_⟩
  · rw [lorentzQ_smul, hshell]
    field_simp
  · have : (mass⁻¹ • v).1 = mass⁻¹ * v.1 := rfl
    rw [this]
    exact mul_pos (inv_pos.mpr hm) hpos
  · show mass • (mass⁻¹ • v) = v
    rw [smul_smul, mul_inv_cancel₀ hmne, one_smul]

/-! ## Composed receipt -/

/-- The single typed antecedent bundle of the internal-clock packet: one
declared positive mass parameter and one declared frame on the declared
Hermitian Lorentz module. -/
structure InternalClockAntecedents where
  /-- The declared mass parameter, a hypothesis and never a derived
  quantity. -/
  mass : ℝ
  /-- Declared positivity of the mass parameter. -/
  mass_pos : 0 < mass
  /-- The declared frame on the future unit-timelike hyperboloid. -/
  frame : FrameHyperboloid

/-- The composed internal-clock receipt.  From the single declared
antecedent bundle: the rest phase is exactly `mass * τ` with proper-time
rate exactly `mass`; the internal clock is oriented-chart invariant; the
period set is exactly the `2π / mass` integer lattice; the de Broglie
spatial gradient and temporal rate readouts equal the momentum and scalar
coordinates; and the committed shell identity ties those two readouts.  All
clauses are conditional on the declared module, the declared plane-wave
evolution, and the declared antecedents; none of them attaches a physical
clock, frequency, or laboratory unit. -/
theorem internalClock_receipt (P : InternalClockAntecedents) :
    (∀ τ : ℝ,
      planeWavePhase (fourMomentum P.mass P.frame)
        (frameWorldline P.frame τ) = P.mass * τ) ∧
    (∀ τ : ℝ,
      HasDerivAt
        (fun σ : ℝ =>
          planeWavePhase (fourMomentum P.mass P.frame)
            (frameWorldline P.frame σ)) P.mass τ) ∧
    (∀ L : OrientedLorentzEquiv,
      internalClock P.mass (L.mapFrame P.frame) =
        internalClock P.mass P.frame) ∧
    (∀ T : ℝ,
      Function.Periodic (internalClock P.mass P.frame) T ↔
        ∃ n : ℤ, T = n * (2 * Real.pi / P.mass)) ∧
    (∀ (t : ℝ) (xs : Spatial) (i : Fin 3) (s : ℝ),
      HasDerivAt
        (fun σ : ℝ =>
          wavePhase (fourMomentum P.mass P.frame)
            (t, xs + σ • (Pi.single i 1 : Spatial)))
        ((fourMomentum P.mass P.frame).2 i) s) ∧
    (∀ (xs : Spatial) (t : ℝ),
      HasDerivAt
        (fun σ : ℝ => planeWavePhase (fourMomentum P.mass P.frame) (σ, xs))
        (fourMomentum P.mass P.frame).1 t) ∧
    (fourMomentum P.mass P.frame).1 ^ 2 =
      P.mass ^ 2 + spatialNormSq (fourMomentum P.mass P.frame).2 := by
  refine ⟨fun τ => planeWavePhase_frameWorldline P.mass P.frame τ,
    fun τ => hasDerivAt_restPhase P.mass P.frame τ,
    fun L => internalClock_oriented L P.mass P.frame,
    fun T => internalClock_period_rigidity (ne_of_gt P.mass_pos) P.frame T,
    fun t xs i s => hasDerivAt_wavePhase_spatial_axis _ t xs i s,
    fun xs t => hasDerivAt_planeWavePhase_time _ xs t,
    fourMomentum_time_sq_eq_mass_sq_add_spatial P.mass P.frame⟩

#print axioms planeWavePhase_frameWorldline
#print axioms hasDerivAt_restPhase
#print axioms internalClock_rotation_law
#print axioms internalClock_oriented
#print axioms internalClock_periodic_iff
#print axioms internalClock_period_rigidity
#print axioms hasDerivAt_wavePhase_spatial_axis
#print axioms hasDerivAt_planeWavePhase_time
#print axioms internalClock_zero_mass_every_period
#print axioms futureNull_no_rest_chart
#print axioms exists_frame_fourMomentum_eq
#print axioms internalClock_receipt

end

end OPH.C1Lorentz
