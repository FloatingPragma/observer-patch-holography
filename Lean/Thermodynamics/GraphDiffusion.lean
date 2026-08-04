import Mathlib

namespace OPH.Thermodynamics.GraphDiffusion

/-!
# Exact transport identities on a finite oriented graph

This file isolates the algebraic core of finite graph diffusion.  An edge
orientation is bookkeeping: positive edge flux points from `tail` to `head`.
The divergence convention is net outward flux, so a source-free continuity
update is `qNext - q = -dt * divergence`.

The quantity wrappers below keep concentration, temperature, particle amount,
energy, cell volume, heat capacity, particle flux, heat flux, distance, and
clock increment distinct at the Lean type level.  The file does not derive any
physical unit or identify a graph field with a laboratory observable.  Such an
identification is an external realization receipt.  Its explicit one-step maps
prove balance and constitutive identities; they do not assert state positivity,
time-step stability, convergence, or an entropy-production theorem.
-/

variable {V E : Type*}

/-- A finite graph with a chosen bookkeeping orientation on every edge. -/
structure OrientedGraph (V E : Type*) where
  tail : E → V
  head : E → V

/-- Positive edge lengths.  Their physical interpretation is external. -/
structure EdgeDistance (E : Type*) where
  value : E → ℝ
  positive : ∀ e, 0 < value e

/-- A positive clock increment.  Calibration to physical time is external. -/
structure ClockIncrement where
  value : ℝ
  positive : 0 < value

/-- Positive vertex cell volumes relating concentration to particle amount. -/
structure CellVolume (V : Type*) where
  value : V → ℝ
  positive : ∀ v, 0 < value v

/-- Positive vertex heat capacities relating temperature to stored energy. -/
structure HeatCapacity (V : Type*) where
  value : V → ℝ
  positive : ∀ v, 0 < value v

/-- A vertex concentration field. -/
structure Concentration (V : Type*) where
  value : V → ℝ

/-- A vertex temperature field. -/
structure Temperature (V : Type*) where
  value : V → ℝ

/-- A vertex particle-amount field. -/
structure ParticleAmount (V : Type*) where
  value : V → ℝ

/-- A vertex energy field. -/
structure Energy (V : Type*) where
  value : V → ℝ

/-- A particle source rate at each vertex. -/
structure ParticleSourceRate (V : Type*) where
  value : V → ℝ

/-- An energy source rate at each vertex. -/
structure EnergySourceRate (V : Type*) where
  value : V → ℝ

/-- Oriented particle flux on edges. -/
structure ParticleFlux (E : Type*) where
  value : E → ℝ

/-- Oriented heat flux on edges. -/
structure HeatFlux (E : Type*) where
  value : E → ℝ

/-- Diffusive conductance.  Nonnegativity is an explicit theorem premise. -/
structure DiffusiveConductance (E : Type*) where
  value : E → ℝ

/-- Thermal conductance.  Nonnegativity is an explicit theorem premise. -/
structure ThermalConductance (E : Type*) where
  value : E → ℝ

/-- Particle amount represented by a concentration and positive cell volume. -/
def particleAmountOfConcentration (volume : CellVolume V)
    (c : Concentration V) : ParticleAmount V where
  value v := volume.value v * c.value v

/-- Concentration represented by particle amount in positive-volume cells. -/
noncomputable def concentrationOfParticleAmount (volume : CellVolume V)
    (q : ParticleAmount V) : Concentration V where
  value v := q.value v / volume.value v

/-- Typed pointwise relation between particle amount and concentration. -/
def AmountConcentrationBridge (volume : CellVolume V)
    (q : ParticleAmount V) (c : Concentration V) : Prop :=
  ∀ v, q.value v = volume.value v * c.value v

/-- Energy represented by temperature and positive heat capacity. -/
def energyOfTemperature (capacity : HeatCapacity V)
    (T : Temperature V) : Energy V where
  value v := capacity.value v * T.value v

/-- Temperature represented by energy and positive heat capacity. -/
noncomputable def temperatureOfEnergy (capacity : HeatCapacity V)
    (q : Energy V) : Temperature V where
  value v := q.value v / capacity.value v

/-- Typed pointwise relation between stored energy and temperature. -/
def EnergyTemperatureBridge (capacity : HeatCapacity V)
    (q : Energy V) (T : Temperature V) : Prop :=
  ∀ v, q.value v = capacity.value v * T.value v

/-- The amount constructed from concentration satisfies its typed bridge. -/
theorem particleAmountOfConcentration_bridge (volume : CellVolume V)
    (c : Concentration V) :
    AmountConcentrationBridge volume (particleAmountOfConcentration volume c) c := by
  intro v
  rfl

/-- Positive volume makes amount-to-concentration conversion a pointwise inverse. -/
theorem concentrationOfParticleAmount_bridge (volume : CellVolume V)
    (q : ParticleAmount V) :
    AmountConcentrationBridge volume q (concentrationOfParticleAmount volume q) := by
  intro v
  simp only [concentrationOfParticleAmount]
  field_simp [ne_of_gt (volume.positive v)]

/-- The energy constructed from temperature satisfies its typed bridge. -/
theorem energyOfTemperature_bridge (capacity : HeatCapacity V)
    (T : Temperature V) :
    EnergyTemperatureBridge capacity (energyOfTemperature capacity T) T := by
  intro v
  rfl

/-- Positive heat capacity makes energy-to-temperature conversion a pointwise inverse. -/
theorem temperatureOfEnergy_bridge (capacity : HeatCapacity V) (q : Energy V) :
    EnergyTemperatureBridge capacity q (temperatureOfEnergy capacity q) := by
  intro v
  simp only [temperatureOfEnergy]
  field_simp [ne_of_gt (capacity.positive v)]

/-- Difference of a scalar vertex field in the chosen edge orientation. -/
def edgeDifference (G : OrientedGraph V E) (u : V → ℝ) (e : E) : ℝ :=
  u (G.head e) - u (G.tail e)

/-- Difference per positive edge length. -/
noncomputable def edgeGradient (G : OrientedGraph V E) (d : EdgeDistance E)
    (u : V → ℝ) (e : E) : ℝ :=
  edgeDifference G u e / d.value e

/-- Net outward oriented flux: outgoing minus incoming flux at a vertex. -/
noncomputable def divergence [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (j : E → ℝ) (v : V) : ℝ :=
  (∑ e, if G.tail e = v then j e else 0)
    - ∑ e, if G.head e = v then j e else 0

/-- Every internal oriented edge cancels in the global divergence sum. -/
theorem sum_divergence_zero [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (j : E → ℝ) :
    ∑ v, divergence G j v = 0 := by
  classical
  unfold divergence
  rw [Finset.sum_sub_distrib]
  have hout : (∑ v, (∑ e, if G.tail e = v then j e else 0)) = ∑ e, j e := by
    rw [Finset.sum_comm]
    apply Finset.sum_congr rfl
    intro e _
    simp
  have hin : (∑ v, (∑ e, if G.head e = v then j e else 0)) = ∑ e, j e := by
    rw [Finset.sum_comm]
    apply Finset.sum_congr rfl
    intro e _
    simp
  rw [hout, hin, sub_self]

/-- Finite-graph summation by parts for the outward-divergence convention. -/
theorem summation_by_parts [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (u : V → ℝ) (j : E → ℝ) :
    ∑ v, u v * divergence G j v
      = -∑ e, edgeDifference G u e * j e := by
  classical
  unfold divergence edgeDifference
  have hsplit : (∑ v, u v *
      ((∑ e, if G.tail e = v then j e else 0)
        - ∑ e, if G.head e = v then j e else 0))
      = (∑ v, u v * (∑ e, if G.tail e = v then j e else 0))
        - ∑ v, u v * (∑ e, if G.head e = v then j e else 0) := by
    rw [← Finset.sum_sub_distrib]
    apply Finset.sum_congr rfl
    intro v _
    ring
  rw [hsplit]
  have hout : (∑ v, u v * (∑ e, if G.tail e = v then j e else 0))
      = ∑ e, u (G.tail e) * j e := by
    simp_rw [Finset.mul_sum]
    rw [Finset.sum_comm]
    apply Finset.sum_congr rfl
    intro e _
    simp
  have hin : (∑ v, u v * (∑ e, if G.head e = v then j e else 0))
      = ∑ e, u (G.head e) * j e := by
    simp_rw [Finset.mul_sum]
    rw [Finset.sum_comm]
    apply Finset.sum_congr rfl
    intro e _
    simp
  rw [hout, hin, ← Finset.sum_sub_distrib, ← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro e _
  ring

/-- A generic one-step local continuity equation. -/
def LocalContinuity [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q qNext source : V → ℝ) (j : E → ℝ) : Prop :=
  ∀ v, qNext v - q v = dt.value * (source v - divergence G j v)

/-- Local continuity implies the exact global balance law. -/
theorem total_balance_of_local_continuity
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q qNext source : V → ℝ) (j : E → ℝ)
    (hlocal : LocalContinuity G dt q qNext source j) :
    (∑ v, qNext v) - ∑ v, q v = dt.value * ∑ v, source v := by
  classical
  have hsum : (∑ v, (qNext v - q v))
      = ∑ v, (dt.value * (source v - divergence G j v)) :=
    Finset.sum_congr rfl (fun v _ => hlocal v)
  rw [Finset.sum_sub_distrib] at hsum
  calc
    (∑ v, qNext v) - ∑ v, q v
        = ∑ v, dt.value * (source v - divergence G j v) := hsum
    _ = dt.value * ((∑ v, source v) - ∑ v, divergence G j v) := by
      rw [← Finset.sum_sub_distrib, Finset.mul_sum]
    _ = dt.value * ∑ v, source v := by rw [sum_divergence_zero, sub_zero]

/-- With zero total source, every local continuity update conserves total amount. -/
theorem total_conservation_of_zero_source
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q qNext source : V → ℝ) (j : E → ℝ)
    (hlocal : LocalContinuity G dt q qNext source j)
    (hsource : ∑ v, source v = 0) :
    ∑ v, qNext v = ∑ v, q v := by
  have hbalance := total_balance_of_local_continuity
    G dt q qNext source j hlocal
  rw [hsource, mul_zero] at hbalance
  linarith

/-- Fick flux: minus diffusive conductance times concentration gradient. -/
noncomputable def fickParticleFlux (G : OrientedGraph V E) (d : EdgeDistance E)
    (κ : DiffusiveConductance E) (c : Concentration V) : ParticleFlux E where
  value e := -κ.value e * edgeGradient G d c.value e

/-- Fourier heat flux is the same typed graph law for temperature and heat. -/
noncomputable def fourierHeatFlux (G : OrientedGraph V E) (d : EdgeDistance E)
    (κ : ThermalConductance E) (T : Temperature V) : HeatFlux E where
  value e := -κ.value e * edgeGradient G d T.value e

/-- Particle continuity, with amount, source rate, clock, and flux kept distinct. -/
def ParticleLocalContinuity [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q qNext : ParticleAmount V) (source : ParticleSourceRate V)
    (j : ParticleFlux E) : Prop :=
  LocalContinuity G dt q.value qNext.value source.value j.value

/-- Energy continuity, with energy, source rate, clock, and heat flux distinct. -/
def EnergyLocalContinuity [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q qNext : Energy V) (source : EnergySourceRate V)
    (j : HeatFlux E) : Prop :=
  LocalContinuity G dt q.value qNext.value source.value j.value

/-- Canonical explicit continuity update for particle amount.  This definition
does not assert positivity or time-step stability. -/
noncomputable def particleContinuityStep [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q : ParticleAmount V) (source : ParticleSourceRate V)
    (j : ParticleFlux E) : ParticleAmount V where
  value v := q.value v + dt.value * (source.value v - divergence G j.value v)

/-- The canonical particle update satisfies local continuity by construction. -/
theorem particleContinuityStep_local [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q : ParticleAmount V) (source : ParticleSourceRate V)
    (j : ParticleFlux E) :
    ParticleLocalContinuity G dt q
      (particleContinuityStep G dt q source j) source j := by
  intro v
  simp only [particleContinuityStep]
  ring

/-- The canonical particle update has exact global source balance. -/
theorem particleContinuityStep_total_balance
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q : ParticleAmount V) (source : ParticleSourceRate V)
    (j : ParticleFlux E) :
    (∑ v, (particleContinuityStep G dt q source j).value v) - ∑ v, q.value v
      = dt.value * ∑ v, source.value v :=
  total_balance_of_local_continuity G dt q.value
    (particleContinuityStep G dt q source j).value source.value j.value
    (particleContinuityStep_local G dt q source j)

/-- The source-free canonical particle update conserves total amount. -/
theorem particleContinuityStep_total_conservation
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q : ParticleAmount V) (j : ParticleFlux E) :
    ∑ v, (particleContinuityStep G dt q ⟨fun _ => 0⟩ j).value v
      = ∑ v, q.value v := by
  have h := particleContinuityStep_total_balance G dt q
    (⟨fun _ => 0⟩ : ParticleSourceRate V) j
  simp only [Finset.sum_const_zero, mul_zero] at h
  linarith

/-- Canonical explicit continuity update for stored energy.  This definition
does not assert positivity or time-step stability. -/
noncomputable def energyContinuityStep [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q : Energy V) (source : EnergySourceRate V)
    (j : HeatFlux E) : Energy V where
  value v := q.value v + dt.value * (source.value v - divergence G j.value v)

/-- The canonical energy update satisfies local continuity by construction. -/
theorem energyContinuityStep_local [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q : Energy V) (source : EnergySourceRate V) (j : HeatFlux E) :
    EnergyLocalContinuity G dt q (energyContinuityStep G dt q source j) source j := by
  intro v
  simp only [energyContinuityStep]
  ring

/-- The canonical energy update has exact global source balance. -/
theorem energyContinuityStep_total_balance
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q : Energy V) (source : EnergySourceRate V) (j : HeatFlux E) :
    (∑ v, (energyContinuityStep G dt q source j).value v) - ∑ v, q.value v
      = dt.value * ∑ v, source.value v :=
  total_balance_of_local_continuity G dt q.value
    (energyContinuityStep G dt q source j).value source.value j.value
    (energyContinuityStep_local G dt q source j)

/-- The source-free canonical energy update conserves total energy. -/
theorem energyContinuityStep_total_conservation
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q : Energy V) (j : HeatFlux E) :
    ∑ v, (energyContinuityStep G dt q ⟨fun _ => 0⟩ j).value v
      = ∑ v, q.value v := by
  have h := energyContinuityStep_total_balance G dt q
    (⟨fun _ => 0⟩ : EnergySourceRate V) j
  simp only [Finset.sum_const_zero, mul_zero] at h
  linarith

/-- Typed particle continuity gives exact global particle balance. -/
theorem particle_total_balance
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q qNext : ParticleAmount V) (source : ParticleSourceRate V)
    (j : ParticleFlux E)
    (hlocal : ParticleLocalContinuity G dt q qNext source j) :
    (∑ v, qNext.value v) - ∑ v, q.value v
      = dt.value * ∑ v, source.value v :=
  total_balance_of_local_continuity G dt q.value qNext.value
    source.value j.value hlocal

/-- Typed source-free particle continuity conserves total particle amount. -/
theorem particle_total_conservation
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q qNext : ParticleAmount V) (j : ParticleFlux E)
    (hlocal : ParticleLocalContinuity G dt q qNext
      ⟨fun _ => 0⟩ j) :
    ∑ v, qNext.value v = ∑ v, q.value v := by
  apply total_conservation_of_zero_source G dt q.value qNext.value
    (fun _ => 0) j.value hlocal
  simp

/-- Typed energy continuity gives exact global energy balance. -/
theorem energy_total_balance
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q qNext : Energy V) (source : EnergySourceRate V) (j : HeatFlux E)
    (hlocal : EnergyLocalContinuity G dt q qNext source j) :
    (∑ v, qNext.value v) - ∑ v, q.value v
      = dt.value * ∑ v, source.value v :=
  total_balance_of_local_continuity G dt q.value qNext.value
    source.value j.value hlocal

/-- Typed source-free energy continuity conserves total energy. -/
theorem energy_total_conservation
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (dt : ClockIncrement)
    (q qNext : Energy V) (j : HeatFlux E)
    (hlocal : EnergyLocalContinuity G dt q qNext
      ⟨fun _ => 0⟩ j) :
    ∑ v, qNext.value v = ∑ v, q.value v := by
  apply total_conservation_of_zero_source G dt q.value qNext.value
    (fun _ => 0) j.value hlocal
  simp

/-- Fick transport is a conserved-particle specialization once its local
continuity equation is supplied. -/
theorem fick_particle_total_conservation
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (d : EdgeDistance E)
    (dt : ClockIncrement) (κ : DiffusiveConductance E)
    (c : Concentration V) (q qNext : ParticleAmount V)
    (hlocal : ParticleLocalContinuity G dt q qNext
      ⟨fun _ => 0⟩ (fickParticleFlux G d κ c)) :
    ∑ v, qNext.value v = ∑ v, q.value v :=
  particle_total_conservation G dt q qNext
    (fickParticleFlux G d κ c) hlocal

/-- Fourier transport is a conserved-energy specialization once its local
continuity equation is supplied. -/
theorem fourier_energy_total_conservation
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (d : EdgeDistance E)
    (dt : ClockIncrement) (κ : ThermalConductance E)
    (T : Temperature V) (q qNext : Energy V)
    (hlocal : EnergyLocalContinuity G dt q qNext
      ⟨fun _ => 0⟩ (fourierHeatFlux G d κ T)) :
    ∑ v, qNext.value v = ∑ v, q.value v :=
  energy_total_conservation G dt q qNext
    (fourierHeatFlux G d κ T) hlocal

/-- Canonical Fick update, starting from the amount represented by the supplied
concentration and positive cell volumes. -/
noncomputable def fickParticleAmountStep [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (d : EdgeDistance E) (dt : ClockIncrement)
    (volume : CellVolume V) (κ : DiffusiveConductance E)
    (c : Concentration V) (source : ParticleSourceRate V) : ParticleAmount V :=
  particleContinuityStep G dt (particleAmountOfConcentration volume c) source
    (fickParticleFlux G d κ c)

/-- Concentration readback of the canonical Fick amount update. -/
noncomputable def fickConcentrationStep [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (d : EdgeDistance E) (dt : ClockIncrement)
    (volume : CellVolume V) (κ : DiffusiveConductance E)
    (c : Concentration V) (source : ParticleSourceRate V) : Concentration V :=
  concentrationOfParticleAmount volume
    (fickParticleAmountStep G d dt volume κ c source)

/-- The canonical Fick amount update satisfies its local continuity equation. -/
theorem fickParticleAmountStep_local [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (d : EdgeDistance E) (dt : ClockIncrement)
    (volume : CellVolume V) (κ : DiffusiveConductance E)
    (c : Concentration V) (source : ParticleSourceRate V) :
    ParticleLocalContinuity G dt (particleAmountOfConcentration volume c)
      (fickParticleAmountStep G d dt volume κ c source) source
      (fickParticleFlux G d κ c) :=
  particleContinuityStep_local G dt
    (particleAmountOfConcentration volume c) source (fickParticleFlux G d κ c)

/-- Positive cell volume relates the canonical next amount to its concentration
readback, so the next state is not an unconnected amount field. -/
theorem fickParticleAmountStep_bridge [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (d : EdgeDistance E) (dt : ClockIncrement)
    (volume : CellVolume V) (κ : DiffusiveConductance E)
    (c : Concentration V) (source : ParticleSourceRate V) :
    AmountConcentrationBridge volume
      (fickParticleAmountStep G d dt volume κ c source)
      (fickConcentrationStep G d dt volume κ c source) :=
  concentrationOfParticleAmount_bridge volume
    (fickParticleAmountStep G d dt volume κ c source)

/-- The canonical Fick update has exact global source balance without a supplied
local-continuity witness. -/
theorem fickParticleAmountStep_total_balance
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (d : EdgeDistance E) (dt : ClockIncrement)
    (volume : CellVolume V) (κ : DiffusiveConductance E)
    (c : Concentration V) (source : ParticleSourceRate V) :
    (∑ v, (fickParticleAmountStep G d dt volume κ c source).value v)
        - ∑ v, (particleAmountOfConcentration volume c).value v
      = dt.value * ∑ v, source.value v :=
  particleContinuityStep_total_balance G dt
    (particleAmountOfConcentration volume c) source (fickParticleFlux G d κ c)

/-- The source-free canonical Fick update conserves total particle amount. -/
theorem fickParticleAmountStep_total_conservation
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (d : EdgeDistance E) (dt : ClockIncrement)
    (volume : CellVolume V) (κ : DiffusiveConductance E) (c : Concentration V) :
    ∑ v, (fickParticleAmountStep G d dt volume κ c ⟨fun _ => 0⟩).value v
      = ∑ v, (particleAmountOfConcentration volume c).value v :=
  particleContinuityStep_total_conservation G dt
    (particleAmountOfConcentration volume c) (fickParticleFlux G d κ c)

/-- Canonical Fourier update, starting from the energy represented by the
supplied temperature and positive heat capacities. -/
noncomputable def fourierEnergyStep [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (d : EdgeDistance E) (dt : ClockIncrement)
    (capacity : HeatCapacity V) (κ : ThermalConductance E)
    (T : Temperature V) (source : EnergySourceRate V) : Energy V :=
  energyContinuityStep G dt (energyOfTemperature capacity T) source
    (fourierHeatFlux G d κ T)

/-- Temperature readback of the canonical Fourier energy update. -/
noncomputable def fourierTemperatureStep [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (d : EdgeDistance E) (dt : ClockIncrement)
    (capacity : HeatCapacity V) (κ : ThermalConductance E)
    (T : Temperature V) (source : EnergySourceRate V) : Temperature V :=
  temperatureOfEnergy capacity (fourierEnergyStep G d dt capacity κ T source)

/-- The canonical Fourier energy update satisfies its local continuity equation. -/
theorem fourierEnergyStep_local [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (d : EdgeDistance E) (dt : ClockIncrement)
    (capacity : HeatCapacity V) (κ : ThermalConductance E)
    (T : Temperature V) (source : EnergySourceRate V) :
    EnergyLocalContinuity G dt (energyOfTemperature capacity T)
      (fourierEnergyStep G d dt capacity κ T source) source
      (fourierHeatFlux G d κ T) :=
  energyContinuityStep_local G dt (energyOfTemperature capacity T) source
    (fourierHeatFlux G d κ T)

/-- Positive heat capacity relates the canonical next energy to its temperature
readback, so the next state is not an unconnected energy field. -/
theorem fourierEnergyStep_bridge [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (d : EdgeDistance E) (dt : ClockIncrement)
    (capacity : HeatCapacity V) (κ : ThermalConductance E)
    (T : Temperature V) (source : EnergySourceRate V) :
    EnergyTemperatureBridge capacity
      (fourierEnergyStep G d dt capacity κ T source)
      (fourierTemperatureStep G d dt capacity κ T source) :=
  temperatureOfEnergy_bridge capacity
    (fourierEnergyStep G d dt capacity κ T source)

/-- The canonical Fourier update has exact global source balance without a
supplied local-continuity witness. -/
theorem fourierEnergyStep_total_balance
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (d : EdgeDistance E) (dt : ClockIncrement)
    (capacity : HeatCapacity V) (κ : ThermalConductance E)
    (T : Temperature V) (source : EnergySourceRate V) :
    (∑ v, (fourierEnergyStep G d dt capacity κ T source).value v)
        - ∑ v, (energyOfTemperature capacity T).value v
      = dt.value * ∑ v, source.value v :=
  energyContinuityStep_total_balance G dt (energyOfTemperature capacity T)
    source (fourierHeatFlux G d κ T)

/-- The source-free canonical Fourier update conserves total energy. -/
theorem fourierEnergyStep_total_conservation
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (d : EdgeDistance E) (dt : ClockIncrement)
    (capacity : HeatCapacity V) (κ : ThermalConductance E) (T : Temperature V) :
    ∑ v, (fourierEnergyStep G d dt capacity κ T ⟨fun _ => 0⟩).value v
      = ∑ v, (energyOfTemperature capacity T).value v :=
  energyContinuityStep_total_conservation G dt (energyOfTemperature capacity T)
    (fourierHeatFlux G d κ T)

/-- Length-weighted Dirichlet dissipation of a scalar vertex field. -/
noncomputable def dirichletDissipation [Fintype E]
    (G : OrientedGraph V E) (d : EdgeDistance E)
    (κ : E → ℝ) (u : V → ℝ) : ℝ :=
  ∑ e, κ e * d.value e * (edgeGradient G d u e) ^ 2

/-- Fick's law turns flux-gradient power into minus Dirichlet dissipation. -/
theorem fick_flux_gradient_power [Fintype E]
    (G : OrientedGraph V E) (d : EdgeDistance E)
    (κ : DiffusiveConductance E) (c : Concentration V) :
    ∑ e, d.value e * (fickParticleFlux G d κ c).value e
        * edgeGradient G d c.value e
      = -dirichletDissipation G d κ.value c.value := by
  unfold fickParticleFlux dirichletDissipation
  rw [← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro e _
  ring

/-- Fourier's law turns heat-flux-gradient power into minus Dirichlet
dissipation. -/
theorem fourier_flux_gradient_power [Fintype E]
    (G : OrientedGraph V E) (d : EdgeDistance E)
    (κ : ThermalConductance E) (T : Temperature V) :
    ∑ e, d.value e * (fourierHeatFlux G d κ T).value e
        * edgeGradient G d T.value e
      = -dirichletDissipation G d κ.value T.value := by
  unfold fourierHeatFlux dirichletDissipation
  rw [← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro e _
  ring

/-- Nonnegative conductances give nonnegative Dirichlet dissipation. -/
theorem dirichletDissipation_nonnegative [Fintype E]
    (G : OrientedGraph V E) (d : EdgeDistance E)
    (κ : E → ℝ) (u : V → ℝ) (hκ : ∀ e, 0 ≤ κ e) :
    0 ≤ dirichletDissipation G d κ u := by
  unfold dirichletDissipation
  apply Finset.sum_nonneg
  intro e _
  exact mul_nonneg (mul_nonneg (hκ e) (le_of_lt (d.positive e))) (sq_nonneg _)

/-- The Fick flux-gradient power is nonpositive for nonnegative conductance. -/
theorem fick_flux_gradient_power_nonpositive [Fintype E]
    (G : OrientedGraph V E) (d : EdgeDistance E)
    (κ : DiffusiveConductance E) (c : Concentration V)
    (hκ : ∀ e, 0 ≤ κ.value e) :
    ∑ e, d.value e * (fickParticleFlux G d κ c).value e
        * edgeGradient G d c.value e ≤ 0 := by
  rw [fick_flux_gradient_power]
  exact neg_nonpos.mpr (dirichletDissipation_nonnegative G d κ.value c.value hκ)

/-- The Fourier heat-flux-gradient power is nonpositive for nonnegative
thermal conductance. -/
theorem fourier_flux_gradient_power_nonpositive [Fintype E]
    (G : OrientedGraph V E) (d : EdgeDistance E)
    (κ : ThermalConductance E) (T : Temperature V)
    (hκ : ∀ e, 0 ≤ κ.value e) :
    ∑ e, d.value e * (fourierHeatFlux G d κ T).value e
        * edgeGradient G d T.value e ≤ 0 := by
  rw [fourier_flux_gradient_power]
  exact neg_nonpos.mpr (dirichletDissipation_nonnegative G d κ.value T.value hκ)

/-- Divergence pairing of a Fick flux equals the Dirichlet form. -/
theorem fick_divergence_pairing
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (d : EdgeDistance E)
    (κ : DiffusiveConductance E) (c : Concentration V) :
    ∑ v, c.value v * divergence G (fickParticleFlux G d κ c).value v
      = dirichletDissipation G d κ.value c.value := by
  rw [summation_by_parts]
  unfold dirichletDissipation fickParticleFlux edgeGradient
  rw [← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro e _
  have hd : d.value e ≠ 0 := ne_of_gt (d.positive e)
  field_simp

/-- Fourier flux obeys the same exact dissipation identity, with distinct types. -/
theorem fourier_divergence_pairing
    [Fintype V] [Fintype E] [DecidableEq V]
    (G : OrientedGraph V E) (d : EdgeDistance E)
    (κ : ThermalConductance E) (T : Temperature V) :
    ∑ v, T.value v * divergence G (fourierHeatFlux G d κ T).value v
      = dirichletDissipation G d κ.value T.value := by
  rw [summation_by_parts]
  unfold dirichletDissipation fourierHeatFlux edgeGradient
  rw [← Finset.sum_neg_distrib]
  apply Finset.sum_congr rfl
  intro e _
  have hd : d.value e ≠ 0 := ne_of_gt (d.positive e)
  field_simp

section TwoVertexExample

/-- One edge oriented from vertex `0` to vertex `1`. -/
def twoVertexGraph : OrientedGraph (Fin 2) Unit where
  tail _ := 0
  head _ := 1

def unitEdgeDistance : EdgeDistance Unit where
  value _ := 1
  positive _ := by norm_num

def unitDiffusiveConductance : DiffusiveConductance Unit where
  value _ := 1

def negativeDiffusiveConductance : DiffusiveConductance Unit where
  value _ := -1

def unitThermalConductance : ThermalConductance Unit where
  value _ := 1

def negativeThermalConductance : ThermalConductance Unit where
  value _ := -1

def unitCellVolume : CellVolume (Fin 2) where
  value _ := 1
  positive _ := by norm_num

def unitHeatCapacity : HeatCapacity (Fin 2) where
  value _ := 1
  positive _ := by norm_num

noncomputable def quarterClockIncrement : ClockIncrement where
  value := 1 / 4
  positive := by norm_num

def highAtTail : Concentration (Fin 2) where
  value v := if v = 0 then 2 else 0

def hotAtTail : Temperature (Fin 2) where
  value v := if v = 0 then 2 else 0

/-- In the two-vertex test, Fick flux points from high to low concentration. -/
theorem twoVertex_fick_flux :
    (fickParticleFlux twoVertexGraph unitEdgeDistance
      unitDiffusiveConductance highAtTail).value () = 2 := by
  norm_num [fickParticleFlux, edgeGradient, edgeDifference,
    twoVertexGraph, unitEdgeDistance, unitDiffusiveConductance, highAtTail]

/-- The same test has strictly positive Dirichlet dissipation. -/
theorem twoVertex_positive_dissipation :
    dirichletDissipation twoVertexGraph unitEdgeDistance
      unitDiffusiveConductance.value highAtTail.value = 4 := by
  norm_num [dirichletDissipation, edgeGradient, edgeDifference,
    twoVertexGraph, unitEdgeDistance, unitDiffusiveConductance, highAtTail]

/-- Negative conductance is a sharp negative control for nonnegativity. -/
theorem negative_conductance_counterexample :
    dirichletDissipation twoVertexGraph unitEdgeDistance
      negativeDiffusiveConductance.value highAtTail.value = -4 := by
  norm_num [dirichletDissipation, edgeGradient, edgeDifference,
    twoVertexGraph, unitEdgeDistance, negativeDiffusiveConductance, highAtTail]

/-- Negative thermal conductance reverses the downhill heat-flux sign. -/
theorem negative_thermal_conductance_uphill :
    (fourierHeatFlux twoVertexGraph unitEdgeDistance
      negativeThermalConductance hotAtTail).value () = -2 := by
  norm_num [fourierHeatFlux, edgeGradient, edgeDifference,
    twoVertexGraph, unitEdgeDistance, negativeThermalConductance, hotAtTail]

/-- Negative thermal conductance is a sharp control for thermal dissipation. -/
theorem negative_thermal_conductance_counterexample :
    dirichletDissipation twoVertexGraph unitEdgeDistance
      negativeThermalConductance.value hotAtTail.value = -4 := by
  norm_num [dirichletDissipation, edgeGradient, edgeDifference,
    twoVertexGraph, unitEdgeDistance, negativeThermalConductance, hotAtTail]

/-- A source-free quarter-step gives a nonvacuous closed Fick update: amount
moves from the high-concentration tail to the low-concentration head. -/
theorem twoVertex_fick_closed_step :
    (fickParticleAmountStep twoVertexGraph unitEdgeDistance quarterClockIncrement
        unitCellVolume unitDiffusiveConductance highAtTail
        (⟨fun _ => 0⟩ : ParticleSourceRate (Fin 2))).value 0 = 3 / 2 ∧
      (fickParticleAmountStep twoVertexGraph unitEdgeDistance quarterClockIncrement
        unitCellVolume unitDiffusiveConductance highAtTail
        (⟨fun _ => 0⟩ : ParticleSourceRate (Fin 2))).value 1 = 1 / 2 := by
  norm_num [fickParticleAmountStep, particleContinuityStep,
    particleAmountOfConcentration, fickParticleFlux, divergence, edgeGradient,
    edgeDifference, twoVertexGraph, unitEdgeDistance, quarterClockIncrement,
    unitCellVolume, unitDiffusiveConductance, highAtTail]

/-- The analogous source-free quarter-step transfers stored energy from the
hot tail to the cold head. -/
theorem twoVertex_fourier_closed_step :
    (fourierEnergyStep twoVertexGraph unitEdgeDistance quarterClockIncrement
        unitHeatCapacity unitThermalConductance hotAtTail
        (⟨fun _ => 0⟩ : EnergySourceRate (Fin 2))).value 0 = 3 / 2 ∧
      (fourierEnergyStep twoVertexGraph unitEdgeDistance quarterClockIncrement
        unitHeatCapacity unitThermalConductance hotAtTail
        (⟨fun _ => 0⟩ : EnergySourceRate (Fin 2))).value 1 = 1 / 2 := by
  norm_num [fourierEnergyStep, energyContinuityStep, energyOfTemperature,
    fourierHeatFlux, divergence, edgeGradient, edgeDifference, twoVertexGraph,
    unitEdgeDistance, quarterClockIncrement, unitHeatCapacity,
    unitThermalConductance, hotAtTail]

end TwoVertexExample

end OPH.Thermodynamics.GraphDiffusion

#print axioms OPH.Thermodynamics.GraphDiffusion.sum_divergence_zero
#print axioms OPH.Thermodynamics.GraphDiffusion.summation_by_parts
#print axioms OPH.Thermodynamics.GraphDiffusion.total_balance_of_local_continuity
#print axioms OPH.Thermodynamics.GraphDiffusion.particle_total_conservation
#print axioms OPH.Thermodynamics.GraphDiffusion.concentrationOfParticleAmount_bridge
#print axioms OPH.Thermodynamics.GraphDiffusion.temperatureOfEnergy_bridge
#print axioms OPH.Thermodynamics.GraphDiffusion.particleContinuityStep_total_balance
#print axioms OPH.Thermodynamics.GraphDiffusion.energyContinuityStep_total_balance
#print axioms OPH.Thermodynamics.GraphDiffusion.fickParticleAmountStep_bridge
#print axioms OPH.Thermodynamics.GraphDiffusion.fickParticleAmountStep_total_conservation
#print axioms OPH.Thermodynamics.GraphDiffusion.fourierEnergyStep_bridge
#print axioms OPH.Thermodynamics.GraphDiffusion.fourierEnergyStep_total_conservation
#print axioms OPH.Thermodynamics.GraphDiffusion.dirichletDissipation_nonnegative
#print axioms OPH.Thermodynamics.GraphDiffusion.fourier_flux_gradient_power_nonpositive
#print axioms OPH.Thermodynamics.GraphDiffusion.fick_divergence_pairing
#print axioms OPH.Thermodynamics.GraphDiffusion.negative_conductance_counterexample
#print axioms OPH.Thermodynamics.GraphDiffusion.negative_thermal_conductance_counterexample
#print axioms OPH.Thermodynamics.GraphDiffusion.twoVertex_fick_closed_step
#print axioms OPH.Thermodynamics.GraphDiffusion.twoVertex_fourier_closed_step
