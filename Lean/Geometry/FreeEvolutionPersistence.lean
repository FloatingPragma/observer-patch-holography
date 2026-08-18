import Geometry.InternalClockRestFrequency

/-!
# C1: free-evolution persistence on the declared Lorentz module

On the declared Hermitian Lorentz module, this file declares the free
evolution as the one-parameter translation flow along a declared frame
worldline, acting on wave functions by precomposition, and proves the
persistence packet for the declared plane-wave states of
`Geometry/InternalClockRestFrequency.lean`.  The flow is an exact
one-parameter group.  Under the flow a declared plane-wave state changes by
a global unit-modulus phase only: its momentum label is unchanged, so
mass-shell membership is preserved exactly at every parameter value for
every declared mass; the temporal and spatial de Broglie phase readouts of
the evolved state equal the unchanged scalar and momentum coordinates at
every parameter, so the labels are constants of the flow; the evolved state
has modulus one at every parameter and every event, so it never decays and
never vanishes; and the flow is injective on declared positive-mass labels,
so two states with distinct mass or momentum labels stay distinct at every
parameter, with no spontaneous transition.  Along the state's own frame the
global flow phase is exactly the internal clock of the wave-1 module.  On
the zero-mass branch the massive parameterization degenerates to the
constant state, which the flow fixes and whose frame label is not
recoverable, while a future null vector carries the same persistence packet
on the zero-mass shell.  A composed receipt derives all clauses from a
single typed antecedent bundle.

This packet is an exact invariance property of a declared kinematic flow.
It is distinct from record-repair stability: nothing here concerns the
declared repair law, conditional expectations, entropy production, or
thermodynamic persistence of records.  The flow, the module, and the mass
and frame labels are declared hypotheses, never derived quantities.
Nothing here constructs a physical particle, a source-produced state, an
interaction, a decay rate, a localized wave packet, or a laboratory
persistence time; the clock and energy calibration anchors, the
observer-to-physical-spacetime attachment, and the source field and action
attachment remain open register premises.  This is a formal precursor for
the particle-persistence observation row, not a physics attainment.
-/

namespace OPH.C1Lorentz

noncomputable section

open OPH.C2Soldering

/-! ## The declared free-evolution flow -/

/-- The declared free evolution: translation of the wave function along the
declared frame worldline by the flow parameter.  The flow direction is a
declared hypothesis, never a derived quantity. -/
def freeEvolution (u : FrameHyperboloid) (τ : ℝ) (ψ : Herm2 → ℂ) :
    Herm2 → ℂ :=
  fun x => ψ (x + frameWorldline u τ)

/-- At parameter zero the declared flow is the identity. -/
theorem freeEvolution_zero (u : FrameHyperboloid) (ψ : Herm2 → ℂ) :
    freeEvolution u 0 ψ = ψ := by
  funext x
  simp [freeEvolution, frameWorldline]

/-- The declared flow is an exact one-parameter group: composing the flow at
two parameters is the flow at their sum. -/
theorem freeEvolution_comp (u : FrameHyperboloid) (τ σ : ℝ) (ψ : Herm2 → ℂ) :
    freeEvolution u σ (freeEvolution u τ ψ) = freeEvolution u (τ + σ) ψ := by
  funext x
  simp only [freeEvolution, frameWorldline]
  congr 1
  rw [add_smul]
  abel

/-! ## Global-phase covariance of plane-wave states -/

/-- The plane-wave phase is additive in the event coordinate. -/
theorem planeWavePhase_add_event (p x y : Herm2) :
    planeWavePhase p (x + y) = planeWavePhase p x + planeWavePhase p y := by
  simp only [planeWavePhase]
  exact lorentzB_add_right p x y

/-- The plane-wave phase along the flow displacement. -/
theorem planeWavePhase_flow_displacement (u : FrameHyperboloid) (p : Herm2)
    (τ : ℝ) :
    planeWavePhase p (frameWorldline u τ) = τ * lorentzB p (u : Herm2) := by
  simp only [planeWavePhase, frameWorldline, lorentzB_smul_right]

/-- The global phase acquired by a plane-wave state under the declared flow. -/
def flowPhase (u : FrameHyperboloid) (p : Herm2) (τ : ℝ) : ℂ :=
  Complex.exp (-Complex.I * ((τ * lorentzB p (u : Herm2) : ℝ) : ℂ))

/-- Under the declared flow a plane-wave state changes by the global flow
phase only: the momentum label is unchanged at every parameter value. -/
theorem freeEvolution_planeWave (u : FrameHyperboloid) (p : Herm2) (τ : ℝ) :
    freeEvolution u τ (planeWave p) =
      fun x => flowPhase u p τ * planeWave p x := by
  funext x
  simp only [freeEvolution, planeWave, flowPhase]
  rw [← Complex.exp_add]
  congr 1
  rw [planeWavePhase_add_event, planeWavePhase_flow_displacement]
  push_cast
  ring

/-- Along the state's own frame the global flow phase is exactly the
internal clock of the wave-1 module. -/
theorem flowPhase_self_eq_internalClock (mass : ℝ) (frame : FrameHyperboloid)
    (τ : ℝ) :
    flowPhase frame (fourMomentum mass frame) τ = internalClock mass frame τ := by
  have hB : lorentzB (fourMomentum mass frame) (frame : Herm2) = mass := by
    rw [fourMomentum, lorentzB_smul_left, lorentzB_self, frame.2.1, mul_one]
  rw [internalClock_eq_exp, flowPhase, hB, mul_comm τ mass]

/-! ## Clause 1: exact shell invariance -/

/-- Shell invariance: at every parameter value and for every declared mass,
the evolved state is a global-phase multiple of the plane wave of the same
momentum, and that momentum satisfies the exact mass-shell identity. -/
theorem freeEvolution_shell_invariance (u : FrameHyperboloid) (mass : ℝ)
    (frame : FrameHyperboloid) (τ : ℝ) :
    lorentzQ (fourMomentum mass frame) = mass ^ 2 ∧
      freeEvolution u τ (planeWave (fourMomentum mass frame)) =
        fun x => flowPhase u (fourMomentum mass frame) τ *
          planeWave (fourMomentum mass frame) x :=
  ⟨lorentzQ_fourMomentum mass frame, freeEvolution_planeWave u _ τ⟩

/-! ## Clause 2: label conservation along the flow -/

/-- Label conservation, temporal readout: at every flow parameter the
temporal rate of the evolved phase equals the scalar coordinate of the
unchanged momentum label. -/
theorem hasDerivAt_evolved_phase_time (u : FrameHyperboloid) (p : Herm2)
    (τ : ℝ) (xs : Spatial) (t : ℝ) :
    HasDerivAt
      (fun σ : ℝ => planeWavePhase p ((σ, xs) + frameWorldline u τ))
      p.1 t := by
  have hfun :
      (fun σ : ℝ => planeWavePhase p ((σ, xs) + frameWorldline u τ)) =
      fun σ : ℝ =>
        (τ * lorentzB p (u : Herm2) - spatialDot p.2 xs) + σ * p.1 := by
    funext σ
    rw [planeWavePhase_add_event, planeWavePhase_flow_displacement]
    simp only [planeWavePhase, lorentzB]
    ring
  rw [hfun]
  exact hasDerivAt_const_add_mul _ _ t

/-- Label conservation, spatial readout: at every flow parameter the spatial
phase gradient of the evolved state along the `i`-th axis equals the `i`-th
coordinate of the unchanged momentum label. -/
theorem hasDerivAt_evolved_phase_spatial_axis (u : FrameHyperboloid)
    (p : Herm2) (τ : ℝ) (t : ℝ) (xs : Spatial) (i : Fin 3) (s : ℝ) :
    HasDerivAt
      (fun σ : ℝ =>
        wavePhase p
          ((t, xs + σ • (Pi.single i 1 : Spatial)) + frameWorldline u τ))
      (p.2 i) s := by
  have hfun :
      (fun σ : ℝ =>
        wavePhase p
          ((t, xs + σ • (Pi.single i 1 : Spatial)) + frameWorldline u τ)) =
      fun σ : ℝ =>
        (spatialDot p.2 (xs + (frameWorldline u τ).2) -
          p.1 * (t + (frameWorldline u τ).1)) + σ * p.2 i := by
    funext σ
    simp only [wavePhase, Prod.fst_add, Prod.snd_add]
    rw [add_right_comm xs (σ • (Pi.single i 1 : Spatial))
        (frameWorldline u τ).2,
      spatialDot_add_smul, spatialDot_single]
    ring
  rw [hfun]
  exact hasDerivAt_const_add_mul _ _ s

/-! ## Clause 3: norm persistence -/

/-- The modulus of a unit-circle exponential with real angle. -/
theorem norm_exp_neg_I_mul_real (a : ℝ) :
    ‖Complex.exp (-Complex.I * (a : ℂ))‖ = 1 := by
  have h : -Complex.I * (a : ℂ) = ((-a : ℝ) : ℂ) * Complex.I := by
    push_cast
    ring
  rw [h, Complex.norm_exp_ofReal_mul_I]

/-- A plane-wave state has modulus one at every event. -/
theorem planeWave_norm (p x : Herm2) : ‖planeWave p x‖ = 1 := by
  rw [planeWave]
  exact norm_exp_neg_I_mul_real _

/-- The global flow phase has modulus one at every parameter. -/
theorem flowPhase_norm (u : FrameHyperboloid) (p : Herm2) (τ : ℝ) :
    ‖flowPhase u p τ‖ = 1 := by
  rw [flowPhase]
  exact norm_exp_neg_I_mul_real _

/-- Norm persistence: the evolved state has modulus one at every parameter
value and every event, so the state never decays under the declared flow. -/
theorem freeEvolution_planeWave_norm (u : FrameHyperboloid) (p : Herm2)
    (τ : ℝ) (x : Herm2) :
    ‖freeEvolution u τ (planeWave p) x‖ = 1 := by
  simp only [freeEvolution]
  exact planeWave_norm p _

/-- The evolved state never vanishes at any parameter value or event. -/
theorem freeEvolution_planeWave_ne_zero (u : FrameHyperboloid) (p : Herm2)
    (τ : ℝ) (x : Herm2) :
    freeEvolution u τ (planeWave p) x ≠ 0 := by
  intro h0
  have h1 := freeEvolution_planeWave_norm u p τ x
  rw [h0] at h1
  simp at h1

/-- The wave-1 internal clock has modulus one at every proper time. -/
theorem internalClock_norm (mass : ℝ) (frame : FrameHyperboloid) (τ : ℝ) :
    ‖internalClock mass frame τ‖ = 1 := by
  rw [internalClock_eq_exp]
  exact norm_exp_neg_I_mul_real _

/-! ## Clause 4: no spontaneous transition -/

/-- The Lorentz pairing is subtractive in its left argument. -/
theorem lorentzB_sub_left (p q x : Herm2) :
    lorentzB (p - q) x = lorentzB p x - lorentzB q x := by
  have h : p - q = p + (-1 : ℝ) • q := by
    rw [neg_one_smul, sub_eq_add_neg]
  rw [h, lorentzB_add_left, lorentzB_smul_left]
  ring

/-- A real number whose every real multiple lies on the `2π` integer lattice
is zero. -/
theorem eq_zero_of_forall_mul_mem_two_pi_lattice {c : ℝ}
    (h : ∀ s : ℝ, ∃ n : ℤ, s * c = n * (2 * Real.pi)) : c = 0 := by
  by_contra hc
  obtain ⟨n, hn⟩ := h (Real.pi / c)
  have hcancel : Real.pi / c * c = Real.pi := by
    field_simp
  rw [hcancel] at hn
  have h2 : (2 * (n : ℝ)) * Real.pi = 1 * Real.pi := by
    rw [one_mul]
    calc (2 * (n : ℝ)) * Real.pi = (n : ℝ) * (2 * Real.pi) := by ring
      _ = Real.pi := hn.symm
  have h3 : (2 * (n : ℝ)) = 1 := mul_right_cancel₀ Real.pi_ne_zero h2
  have h4 : (2 * n : ℤ) = 1 := by exact_mod_cast h3
  omega

/-- Plane-wave states determine their momentum labels: equal wave functions
force equal momenta. -/
theorem planeWave_injective : Function.Injective planeWave := by
  intro p q h
  have hone : ∀ x : Herm2, ∃ n : ℤ,
      planeWavePhase p x - planeWavePhase q x = n * (2 * Real.pi) := by
    intro x
    have hx : Complex.exp (-Complex.I * ((planeWavePhase p x : ℝ) : ℂ)) =
        Complex.exp (-Complex.I * ((planeWavePhase q x : ℝ) : ℂ)) :=
      congrFun h x
    have hdiff : Complex.exp
        (-Complex.I *
          ((planeWavePhase p x - planeWavePhase q x : ℝ) : ℂ)) = 1 := by
      rw [show -Complex.I *
            ((planeWavePhase p x - planeWavePhase q x : ℝ) : ℂ) =
          -Complex.I * ((planeWavePhase p x : ℝ) : ℂ) -
            -Complex.I * ((planeWavePhase q x : ℝ) : ℂ) by push_cast; ring,
        Complex.exp_sub, hx, div_self (Complex.exp_ne_zero _)]
    exact (exp_neg_I_mul_real_eq_one_iff _).mp hdiff
  have hlin : ∀ x : Herm2, lorentzB (p - q) x = 0 := by
    intro x
    apply eq_zero_of_forall_mul_mem_two_pi_lattice
    intro s
    obtain ⟨n, hn⟩ := hone (s • x)
    refine ⟨n, ?_⟩
    have hphase : planeWavePhase p (s • x) - planeWavePhase q (s • x) =
        s * lorentzB (p - q) x := by
      simp only [planeWavePhase]
      rw [lorentzB_sub_left]
      rw [lorentzB_smul_right, lorentzB_smul_right]
      ring
    rw [← hphase, hn]
  have hs0 : spatialDot (p - q).2 (0 : Spatial) = 0 := by
    simp [spatialDot]
  have ht : (p - q).1 = 0 := by
    have h1 := hlin ((1 : ℝ), (0 : Spatial))
    simp only [lorentzB] at h1
    rw [hs0] at h1
    linarith
  have hsp : ∀ i : Fin 3, (p - q).2 i = 0 := by
    intro i
    have h1 := hlin ((0 : ℝ), (Pi.single i 1 : Spatial))
    simp only [lorentzB] at h1
    rw [spatialDot_single] at h1
    linarith
  have hpq : p - q = 0 := by
    apply Prod.ext
    · exact ht
    · funext i
      exact hsp i
  exact sub_eq_zero.mp hpq

/-- The declared flow is injective on plane-wave states at every parameter:
evolved states that agree at one parameter carry equal momentum labels. -/
theorem freeEvolution_planeWave_injective (u : FrameHyperboloid) (τ : ℝ)
    {p q : Herm2}
    (h : freeEvolution u τ (planeWave p) = freeEvolution u τ (planeWave q)) :
    p = q := by
  apply planeWave_injective
  funext y
  have hpt : y - frameWorldline u τ + frameWorldline u τ = y := by abel
  have hy := congrFun h (y - frameWorldline u τ)
  simp only [freeEvolution] at hy
  rwa [hpt] at hy

/-- Positive-mass labels are determined by their momenta: equal declared
momenta force equal mass parameters and equal frames. -/
theorem fourMomentum_label_injective {m₁ m₂ : ℝ} (h₁ : 0 < m₁)
    (h₂ : 0 < m₂) {f₁ f₂ : FrameHyperboloid}
    (h : fourMomentum m₁ f₁ = fourMomentum m₂ f₂) : m₁ = m₂ ∧ f₁ = f₂ := by
  have hq : m₁ ^ 2 = m₂ ^ 2 := by
    have hQ := congrArg lorentzQ h
    rwa [lorentzQ_fourMomentum, lorentzQ_fourMomentum] at hQ
  have hfac : (m₁ - m₂) * (m₁ + m₂) = 0 := by
    calc (m₁ - m₂) * (m₁ + m₂) = m₁ ^ 2 - m₂ ^ 2 := by ring
      _ = 0 := by rw [hq]; ring
  have hm : m₁ = m₂ := by
    rcases mul_eq_zero.mp hfac with h0 | h0
    · linarith
    · linarith
  refine ⟨hm, ?_⟩
  apply Subtype.ext
  have hsm : m₂ • (f₁ : Herm2) = m₂ • (f₂ : Herm2) := by
    have h' := h
    rw [fourMomentum, fourMomentum, hm] at h'
    exact h'
  have hcancel := congrArg (fun v : Herm2 => m₂⁻¹ • v) hsm
  simpa [smul_smul, inv_mul_cancel₀ (ne_of_gt h₂)] using hcancel

/-- No spontaneous transition: if the evolved states of two declared
positive-mass labels agree at any single parameter value, the mass
parameters and the frames agree. -/
theorem freeEvolution_no_transition (u : FrameHyperboloid) {m₁ m₂ : ℝ}
    (h₁ : 0 < m₁) (h₂ : 0 < m₂) (f₁ f₂ : FrameHyperboloid) (τ : ℝ)
    (h : freeEvolution u τ (planeWave (fourMomentum m₁ f₁)) =
      freeEvolution u τ (planeWave (fourMomentum m₂ f₂))) :
    m₁ = m₂ ∧ f₁ = f₂ :=
  fourMomentum_label_injective h₁ h₂
    (freeEvolution_planeWave_injective u τ h)

/-- Distinct declared positive-mass labels stay distinct at every parameter
value of the declared flow. -/
theorem freeEvolution_distinct_labels_stay_distinct (u : FrameHyperboloid)
    {m₁ m₂ : ℝ} (h₁ : 0 < m₁) (h₂ : 0 < m₂) {f₁ f₂ : FrameHyperboloid}
    (hne : ¬(m₁ = m₂ ∧ f₁ = f₂)) (τ : ℝ) :
    freeEvolution u τ (planeWave (fourMomentum m₁ f₁)) ≠
      freeEvolution u τ (planeWave (fourMomentum m₂ f₂)) :=
  fun h => hne (freeEvolution_no_transition u h₁ h₂ f₁ f₂ τ h)

/-! ## Clause 5: the zero-mass branch -/

/-- The zero-momentum plane wave is the constant state. -/
theorem planeWave_zero_momentum (x : Herm2) : planeWave 0 x = 1 := by
  have h0 : planeWavePhase 0 x = 0 := by
    simp [planeWavePhase, lorentzB, spatialDot]
  rw [planeWave, h0]
  simp

/-- On the zero-mass branch the declared flow fixes the state exactly: the
massive parameterization degenerates to the constant state. -/
theorem freeEvolution_zero_mass (u frame : FrameHyperboloid) (τ : ℝ) :
    freeEvolution u τ (planeWave (fourMomentum 0 frame)) =
      planeWave (fourMomentum 0 frame) := by
  funext x
  rw [fourMomentum_zero_mass]
  show planeWave 0 (x + frameWorldline u τ) = planeWave 0 x
  rw [planeWave_zero_momentum, planeWave_zero_momentum]

/-- At zero declared mass the frame label is not recoverable: every frame
yields the same momentum, so the positive-mass hypothesis of the
no-transition clause is sharp. -/
theorem fourMomentum_zero_mass_label_degenerate (f₁ f₂ : FrameHyperboloid) :
    fourMomentum 0 f₁ = fourMomentum 0 f₂ := by
  rw [fourMomentum_zero_mass, fourMomentum_zero_mass]

/-- A future null vector carries the persistence packet on the zero-mass
shell: exact shell invariance with vanishing Lorentz square, global-phase
covariance, and norm persistence at every parameter and event. -/
theorem freeEvolution_null_persistence (u : FrameHyperboloid)
    (v : FutureNullVector) (τ : ℝ) :
    lorentzQ v.1 = 0 ∧
      freeEvolution u τ (planeWave v.1) =
        (fun x => flowPhase u v.1 τ * planeWave v.1 x) ∧
      ∀ x : Herm2, ‖freeEvolution u τ (planeWave v.1) x‖ = 1 :=
  ⟨v.2.1, freeEvolution_planeWave u v.1 τ,
    fun x => freeEvolution_planeWave_norm u v.1 τ x⟩

/-! ## Clause 6: the composed receipt -/

/-- The single typed antecedent bundle of the persistence packet: one
declared positive mass parameter, one declared state frame, and one declared
flow frame on the declared Hermitian Lorentz module.  Every field is a
hypothesis, never a derived quantity. -/
structure FreeEvolutionAntecedents where
  /-- The declared mass parameter. -/
  mass : ℝ
  /-- Declared positivity of the mass parameter. -/
  mass_pos : 0 < mass
  /-- The declared frame carrying the state's momentum label. -/
  frame : FrameHyperboloid
  /-- The declared frame whose worldline directs the flow. -/
  flowFrame : FrameHyperboloid

/-- The composed persistence receipt.  From the single declared antecedent
bundle: the momentum label satisfies the exact mass-shell identity and the
evolved state is a global-phase multiple of the unchanged plane wave at
every parameter (shell invariance); the temporal and spatial phase readouts
of the evolved state equal the unchanged scalar and momentum coordinates at
every parameter (label conservation); the evolved state has modulus one at
every parameter and event (norm persistence); agreement of evolved states at
any parameter forces equal mass and frame labels (no spontaneous
transition); the flow is an exact one-parameter group; and along the state's
own frame the global flow phase equals the internal clock of the wave-1
module (clock identity), which pins the declared flow and its phase as
nontrivial.  All clauses are conditional on the declared module, the
declared flow, and the declared antecedents; none of them attaches a
physical particle, a decay rate, or a laboratory persistence time. -/
theorem freeEvolution_receipt (P : FreeEvolutionAntecedents) :
    lorentzQ (fourMomentum P.mass P.frame) = P.mass ^ 2 ∧
    (∀ τ : ℝ,
      freeEvolution P.flowFrame τ (planeWave (fourMomentum P.mass P.frame)) =
        fun x => flowPhase P.flowFrame (fourMomentum P.mass P.frame) τ *
          planeWave (fourMomentum P.mass P.frame) x) ∧
    (∀ (τ t : ℝ) (xs : Spatial),
      HasDerivAt
        (fun σ : ℝ =>
          planeWavePhase (fourMomentum P.mass P.frame)
            ((σ, xs) + frameWorldline P.flowFrame τ))
        (fourMomentum P.mass P.frame).1 t) ∧
    (∀ (τ t : ℝ) (xs : Spatial) (i : Fin 3) (s : ℝ),
      HasDerivAt
        (fun σ : ℝ =>
          wavePhase (fourMomentum P.mass P.frame)
            ((t, xs + σ • (Pi.single i 1 : Spatial)) +
              frameWorldline P.flowFrame τ))
        ((fourMomentum P.mass P.frame).2 i) s) ∧
    (∀ (τ : ℝ) (x : Herm2),
      ‖freeEvolution P.flowFrame τ
        (planeWave (fourMomentum P.mass P.frame)) x‖ = 1) ∧
    (∀ m' : ℝ, 0 < m' → ∀ (f' : FrameHyperboloid) (τ : ℝ),
      freeEvolution P.flowFrame τ (planeWave (fourMomentum P.mass P.frame)) =
        freeEvolution P.flowFrame τ (planeWave (fourMomentum m' f')) →
      P.mass = m' ∧ P.frame = f') ∧
    (∀ ψ : Herm2 → ℂ, freeEvolution P.flowFrame 0 ψ = ψ) ∧
    (∀ (τ σ : ℝ) (ψ : Herm2 → ℂ),
      freeEvolution P.flowFrame σ (freeEvolution P.flowFrame τ ψ) =
        freeEvolution P.flowFrame (τ + σ) ψ) ∧
    (∀ τ : ℝ,
      flowPhase P.frame (fourMomentum P.mass P.frame) τ =
        internalClock P.mass P.frame τ) := by
  refine ⟨lorentzQ_fourMomentum P.mass P.frame,
    fun τ => freeEvolution_planeWave P.flowFrame _ τ,
    fun τ t xs => hasDerivAt_evolved_phase_time P.flowFrame _ τ xs t,
    fun τ t xs i s =>
      hasDerivAt_evolved_phase_spatial_axis P.flowFrame _ τ t xs i s,
    fun τ x => freeEvolution_planeWave_norm P.flowFrame _ τ x,
    fun m' hm' f' τ h =>
      freeEvolution_no_transition P.flowFrame P.mass_pos hm' P.frame f' τ h,
    fun ψ => freeEvolution_zero P.flowFrame ψ,
    fun τ σ ψ => freeEvolution_comp P.flowFrame τ σ ψ,
    fun τ => flowPhase_self_eq_internalClock P.mass P.frame τ⟩

/-- The antecedent bundle is inhabited: unit declared mass at the standard
frame, with the standard frame directing the flow.  The receipt is therefore
not vacuous. -/
def unitRestAntecedents : FreeEvolutionAntecedents :=
  ⟨1, one_pos, standardFrame, standardFrame⟩

example := freeEvolution_receipt unitRestAntecedents

#print axioms freeEvolution_zero
#print axioms freeEvolution_comp
#print axioms freeEvolution_planeWave
#print axioms flowPhase_self_eq_internalClock
#print axioms freeEvolution_shell_invariance
#print axioms hasDerivAt_evolved_phase_time
#print axioms hasDerivAt_evolved_phase_spatial_axis
#print axioms freeEvolution_planeWave_norm
#print axioms freeEvolution_planeWave_ne_zero
#print axioms internalClock_norm
#print axioms planeWave_injective
#print axioms freeEvolution_planeWave_injective
#print axioms fourMomentum_label_injective
#print axioms freeEvolution_no_transition
#print axioms freeEvolution_distinct_labels_stay_distinct
#print axioms freeEvolution_zero_mass
#print axioms fourMomentum_zero_mass_label_degenerate
#print axioms freeEvolution_null_persistence
#print axioms freeEvolution_receipt

end

end OPH.C1Lorentz
