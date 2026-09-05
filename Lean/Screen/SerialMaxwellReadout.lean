import LabeledEventReadout
import NeutralPairJointStationaryWitness

set_option autoImplicit false

namespace OPH.SerialMaxwellReadout

open OPH.LabeledEventReadout
open ObserverPatchHolography.ScalarSeamRepair
open OPH.ScaledMaxwellStability OPH.TemporalMaxwellEvolution
open OPH.NeutralPairCoupledAction OPH.NeutralPairJointStationaryWitness

noncomputable section

/-! Exact classical records and writable ports allow a serial readout of one
state: retain a parent baseline, pair-average, read the response, then restore
the pair using those two records. This is an instrument construction with
declared read/write access, not a quantum measurement or physical selection
theorem. Potential coordinates remain typed separately from electric and
magnetic fields. -/

abbrev Probe := (Fin 12) × (Fin 30)

/-- Feedback consumes only the baseline and response records at the two
written slots; all other coordinates are carried through. -/
def feedback (p : Probe) (b r : ℝ) (y : Field) : Field :=
  fun s => if s = .inl p.1 then b else if s = .inr p.2 then 2*r-b else y s

def response (p : Probe) (x : Field) : ℝ :=
  pairAverage (.inl p.1) (.inr p.2) x (.inl p.1)

def cycle (base : Fin 12 → ℝ) (p : Probe) (x : Field) : Field :=
  feedback p (base p.1) (response p x) (pairAverage (.inl p.1) (.inr p.2) x)

theorem cycle_restores (p : Probe) (x : Field) :
    cycle (baseline x) p x = x := by
  funext s
  by_cases h₁ : s = .inl p.1
  · subst s; simp [cycle, feedback, baseline]
  · by_cases h₂ : s = .inr p.2
    · subst s; simp [cycle, feedback, response, baseline]; ring
    · simp [cycle, feedback, pairAverage, h₁, h₂]

/-- Any finite probe word, including repeated or overlapping probes. -/
def serialState (base : Fin 12 → ℝ) : List Probe → Field → Field
  | [], x => x
  | p :: ps, x => serialState base ps (cycle base p x)

theorem serial_restores (ps : List Probe) (x : Field) :
    serialState (baseline x) ps x = x := by
  induction ps with
  | nil => rfl
  | cons p ps ih => simpa [serialState, cycle_restores] using ih

theorem response_after_prefix (ps : List Probe) (p : Probe) (x : Field) :
    response p (serialState (baseline x) ps x) =
      (x (.inl p.1) + x (.inr p.2))/2 := by
  rw [serial_restores]; simp [response]

def serialMeasure (parent : Fin 30 → Fin 12) (x : Field) : Readout :=
  (baseline x, fun m => response (parent m, m)
    (serialState (baseline x)
      (((List.finRange 30).take m.val).map (fun k => (parent k, k))) x))

theorem serialMeasure_eq (parent : Fin 30 → Fin 12) (x : Field) :
    serialMeasure parent x = measure parent x := by
  apply Prod.ext
  · rfl
  · funext m; simp [serialMeasure, OPH.LabeledEventReadout.measure, response_after_prefix,
      eventResponse_eq_average]

theorem decode_serial (parent : Fin 30 → Fin 12) (x : Field) :
    reconstruct parent (serialMeasure parent x) = x := by
  rw [serialMeasure_eq, reconstruct_measure]

/-- Exact local error law, not a global noisy repeated-probe stability claim. -/
theorem feedback_error (p : Probe) (x : Field) (eb er : ℝ) :
    feedback p (x (.inl p.1)+eb) (response p x+er)
        (pairAverage (.inl p.1) (.inr p.2) x) (.inr p.2) - x (.inr p.2)
      = 2*er-eb := by
  simp [feedback, response]; ring

def potentialPacket (A : Fin 30 → ℝ) (φ : Fin 12 → ℝ) : Field
  | .inl u => φ u
  | .inr e => A e

def decodedPacket (parent : Fin 30 → Fin 12) (A : Fin 30 → ℝ)
    (φ : Fin 12 → ℝ) : Field :=
  reconstruct parent (serialMeasure parent (potentialPacket A φ))

def decodedA (parent : Fin 30 → Fin 12) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) : ℕ → Fin 30 → ℝ :=
  fun n e => decodedPacket parent (A n) (φ n) (.inr e)

def decodedPhi (parent : Fin 30 → Fin 12) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) : ℕ → Fin 12 → ℝ :=
  fun n u => decodedPacket parent (A n) (φ n) (.inl u)

@[simp] theorem decodedA_eq (parent : Fin 30 → Fin 12)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) : decodedA parent A φ = A := by
  funext n e; simp [decodedA, decodedPacket, decode_serial, potentialPacket]

@[simp] theorem decodedPhi_eq (parent : Fin 30 → Fin 12)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) : decodedPhi parent A φ = φ := by
  funext n u; simp [decodedPhi, decodedPacket, decode_serial, potentialPacket]

theorem decoded_fields (parent : Fin 30 → Fin 12) (h : ℝ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    electricFieldScaled h (decodedA parent A φ) (decodedPhi parent A φ) n =
        electricFieldScaled h A φ n ∧
      magneticField (decodedA parent A φ) n = magneticField A n := by simp

/-- Equality holds for all fields and paths, so also for their variations. -/
theorem decoded_coupled_action (parent : Fin 30 → Fin 12) (q h tp tn : ℝ)
    (N : ℕ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (wp wn : OPH.WorldlineHopTransport.SeamStepWorldline) :
    neutralPairCoupledAction q h tp tn N (decodedA parent A φ)
      (decodedPhi parent A φ) wp wn = neutralPairCoupledAction q h tp tn N A φ wp wn := by
  simp

theorem decoded_joint_field_stationary (parent : Fin 30 → Fin 12) (tp tn : ℝ) :
    NeutralPairFieldStationary 1 (1/2) tp tn 1
      (decodedA parent jointPotential jointScalarPotential)
      (decodedPhi parent jointPotential jointScalarPotential)
      (OPH.WorldlineHopTransport.crossingWorldline 0)
      (OPH.WorldlineHopTransport.crossingWorldline 29) := by
  simpa using joint_field_stationary tp tn

end
end OPH.SerialMaxwellReadout

#print axioms OPH.SerialMaxwellReadout.serial_restores
#print axioms OPH.SerialMaxwellReadout.decode_serial
#print axioms OPH.SerialMaxwellReadout.feedback_error
#print axioms OPH.SerialMaxwellReadout.decoded_coupled_action
#print axioms OPH.SerialMaxwellReadout.decoded_joint_field_stationary
