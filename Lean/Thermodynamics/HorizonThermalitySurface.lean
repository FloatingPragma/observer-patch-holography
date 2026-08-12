import CapFirstLaw
import EinsteinPremiseLink
import FourLawAdequacySurface
import ObserverPatchHolography.DeSitterCapacityShock

set_option autoImplicit false

namespace OPH.Thermodynamics

/-!
# The horizon thermality surface

One parameterized package for the observation-ledger row OL-B3 (horizon
entropy and temperature), V3 issue #729. Every theorem below is a
committed result of the thermodynamics modules or of the finite-capacity
module, re-exported through one named premise bundle so that the ledger
row cites a single surface.

`HorizonRecordData` is the horizon-facing declared premise of the V3
program: the central-charge identification of the horizon record. Over
the declared repair law (register row PR-07, the structure
`RepairLawData` of `FourLawAdequacySurface`) it declares a split of the
modular weight of the declared reference, `-log ref = 2π B + Z`, with
the central charge `Z` measurable through the repaired visible datum.
The identification of `B` with a physical boost energy and of `Z` with
the edge record of a physical horizon is content of this declared row;
no theorem here derives it.

The surface exports three groups of statements.

* The cap Clausius inequality `2π Δ⟨B⟩ ≤ ΔS` for one step of the
  declared repair kernel: the horizon-thermality-shaped statement on
  the finite model. The coefficient `2π` on the boost charge is the
  modular normalization declared in the split, and the repair kernel
  conserves the central charge with no further hypothesis, so no
  independent temperature input enters the inequality.
* The entanglement first law: the exact split
  `ΔS = 2π Δ⟨B⟩ + Δ⟨Z⟩ - D(p ‖ ref)` whose first-order truncation is
  `δS = 2π δ⟨B⟩ + δ⟨Z⟩`, and the differential form of that truncation
  on the simplex tangent space through the Einstein-branch premise
  discharge of `EinsteinPremiseLink`.
* The finite-capacity and de Sitter identities of
  `ObserverPatchHolography.DeSitterCapacityShock` at their committed
  level. Their hypotheses do not factor through the declared repair law
  or the horizon record, so they are re-exported side by side inside
  the namespace; the attachment of the analytic area functional to a
  physical horizon stays with the paper-level shock-spectrum
  diagnostic, as recorded in that module.

Boundary. This surface carries horizon thermality in shape, and it
carries no horizon temperature claim. A temperature in physical units
requires the clock-and-energy calibration anchors row (the structure
`CalibrationData` of `FourLawAdequacySurface`, shared with the
constants lane) together with the horizon dictionary identifying the
abstract split with physical horizon data; both are declared rows, and
no continuum horizon exists on the finite model.
`modular_weight_double_identification` records that under the horizon
record and the calibration row the boost-central split and the thermal
calibration identify one and the same modular weight, so every
temperature statement reachable from this surface is a composite of two
declared identifications with the declared repair law. That composition
replaces any independent temperature input.
-/

/-- **The horizon record: central-charge identification.** A declared
premise of the V3 program over the declared repair law: the modular
weight of the declared reference splits as `-log ref = 2π B + Z`, and
the central charge `Z` is measurable through the repaired visible
datum. The identification of `B` with a physical boost energy and of
`Z` with the edge record of a physical horizon is content of this
declared row; the theorems below consume only the split and the
measurability. -/
structure HorizonRecordData {Ω : Type*} [Fintype Ω] [DecidableEq Ω]
    {B : Type*} [DecidableEq B] (D : RepairLawData Ω B) where
  /-- The bulk boost charge of the split. -/
  boostCharge : Ω → ℝ
  /-- The central charge of the split. -/
  centralCharge : Ω → ℝ
  /-- The central-charge identification of the modular weight. -/
  split : ∀ x, -Real.log (D.ref x)
    = 2 * Real.pi * boostCharge x + centralCharge x
  /-- The central charge is measurable through the repaired visible
  datum. -/
  central_measurable : ∀ x y, D.visible x = D.visible y →
    centralCharge x = centralCharge y

namespace HorizonSurface

open FourLawSurface

variable {Ω : Type*} [Fintype Ω] [DecidableEq Ω]
variable {B : Type*} [DecidableEq B]
variable (D : RepairLawData Ω B) (H : HorizonRecordData D)

/-! ## Entanglement first law -/

/-- **Entanglement first law, exact split form.** The entropy
difference to the declared reference splits through the horizon record
as `ΔS = 2π Δ⟨B⟩ + Δ⟨Z⟩ - D(p ‖ ref)` exactly; the first-order
entanglement first law `δS = 2π δ⟨B⟩ + δ⟨Z⟩` is this identity with the
quadratic-order deficit dropped. Re-export of `cap_firstLaw_split` on
the declared reference. -/
theorem entanglement_first_law_exact (p : Ω → ℝ) :
    shannon p - shannon D.ref
      = 2 * Real.pi * ((∑ x, p x * H.boostCharge x)
          - (∑ x, D.ref x * H.boostCharge x))
        + ((∑ x, p x * H.centralCharge x)
          - (∑ x, D.ref x * H.centralCharge x))
        - kl p D.ref :=
  cap_firstLaw_split p D.ref H.boostCharge H.centralCharge D.ref_pos
    H.split

/-- **Entanglement first law, first-order differential form.** On every
mass-preserving variation of the simplex tangent space, the entropy
differential at the declared reference equals `2π` times the bulk
pairing plus the edge differential of the horizon record. Re-export of
the Einstein-branch premise discharge
`thermo_first_law_on_simplex_tangent`, instantiated on the declared
reference and the horizon record. -/
theorem entanglement_first_law_first_order (δρ : MassZero Ω) :
    (thermoFirstLawData Ω D.ref H.boostCharge
        H.centralCharge).entropyDifferential δρ
      = 2 * Real.pi
          * (thermoFirstLawData Ω D.ref H.boostCharge
              H.centralCharge).bulkPairing δρ
        + (thermoFirstLawData Ω D.ref H.boostCharge
            H.centralCharge).edgeDifferential δρ :=
  thermo_first_law_on_simplex_tangent D.ref H.boostCharge
    H.centralCharge H.split δρ

/-! ## Horizon-thermality-shaped Clausius inequality -/

/-- **Cap Clausius inequality, the horizon-thermality-shaped statement
on the finite model.** One step of the declared repair kernel obeys
`2π Δ⟨B⟩ ≤ ΔS`: the entropy change dominates `2π` times the
boost-charge change. The coefficient `2π` is the declared modular
normalization of the split; the composition with the declared repair
law supplies the conservation of the central charge, so no independent
temperature input enters. Re-export of `heatBath_cap_clausius` through
the premise bundle. -/
theorem horizon_cap_clausius (p : Ω → ℝ) (hp : ∀ x, 0 ≤ p x) :
    2 * Real.pi * ((∑ x, push p (repairKernel D) x * H.boostCharge x)
        - (∑ x, p x * H.boostCharge x))
      ≤ shannon (push p (repairKernel D)) - shannon p := by
  rw [repairKernel_def]
  exact heatBath_cap_clausius D.ref_pos H.boostCharge H.centralCharge p
    hp H.split H.central_measurable

/-- **Central-charge conservation under repair.** One repair step
preserves the mean of the central charge of the horizon record, since
the charge is measurable through the repaired visible datum. Re-export
of `first_repair_conserves_fibre_mean` at the central charge. -/
theorem repair_conserves_central_charge (p : Ω → ℝ) :
    ∑ y, push p (repairKernel D) y * H.centralCharge y
      = ∑ x, p x * H.centralCharge x :=
  first_repair_conserves_fibre_mean D H.centralCharge
    H.central_measurable p

/-- **The repair displacement is central-neutral.** The central pairing
of the horizon record vanishes on the displacement of one repair step,
so the first-order entropy change of repair is purely bulk-modular.
Re-export of `repair_variation_central_pairing_zero` through the
premise bundle. -/
theorem repair_displacement_central_pairing_zero (p : Ω → ℝ) :
    pairing Ω H.centralCharge
      (fun y => push p (repairKernel D) y - p y) = 0 := by
  rw [repairKernel_def]
  exact repair_variation_central_pairing_zero D.ref D.visible D.ref_pos
    H.centralCharge H.central_measurable p

/-! ## Boundary to horizon temperature -/

/-- **One modular weight, two declared identifications.** Under the
horizon record and the clock-and-energy calibration anchors row, the
boost-central split and the thermal calibration describe the same
modular weight of the declared reference:
`2π B + Z = β E + log Z_partition` pointwise. Every temperature
statement reachable from this surface is a composite of these two
declared identifications with the declared repair law; no independent
temperature input exists on the surface. -/
theorem modular_weight_double_identification (C : CalibrationData D)
    (x : Ω) :
    2 * Real.pi * H.boostCharge x + H.centralCharge x
      = C.beta * C.energy x + C.logZ := by
  have h1 := H.split x
  have h2 := C.thermal x
  linarith

/-! ## Finite-capacity and de Sitter identities at the committed level

The four statements below re-export
`ObserverPatchHolography.DeSitterCapacityShock`. Their hypotheses do
not factor through the declared repair law or the horizon record, so
they stand side by side inside this namespace instead of consuming the
premise bundle. The attachment of the analytic area functional to a
physical horizon, of the screen kinetic operator to a graph Laplacian,
and of the icosahedral rotation triplet to a gauge mode stays with the
paper-level shock-spectrum diagnostic, as recorded in the source
module. -/

/-- **Capacity entropy identity.** Dimension-proportional sector
weights give the exact finite entropy `log M` at total capacity `M`.
Re-export of `generalizedSectorEntropy_dimension_weights`. -/
theorem capacity_dimension_weight_entropy {α : Type*} [Fintype α]
    (dimension : α → ℝ) {M : ℝ} (hM : 0 < M)
    (hdimension : ∀ i, 0 < dimension i)
    (hsum : ∑ i, dimension i = M) :
    OPH.DeSitterCapacityShock.generalizedSectorEntropy
        (fun i => dimension i / M) dimension
      = Real.log M :=
  OPH.DeSitterCapacityShock.generalizedSectorEntropy_dimension_weights
    dimension hM hdimension hsum

/-- **Capacity-transfer area loss.** Transferring the fraction `f` of a
fixed total capacity from the horizon sectors to an observer ledger
drops the symmetric-ray area by exactly `log (1 - f)` along the
positive-real interpolation. Re-export of `capacityTransferArea_sub`. -/
theorem capacity_transfer_area_drop {d f : ℝ} (hd : 0 < d)
    (hf : f < 1) :
    OPH.DeSitterCapacityShock.capacityTransferArea d f
        - OPH.DeSitterCapacityShock.symmetricArea d
      = Real.log (1 - f) :=
  OPH.DeSitterCapacityShock.capacityTransferArea_sub hd hf

/-- **Baseline maximality on the transfer path.** Positive transfer
fractions strictly decrease the symmetric-ray area, so the untransferred
configuration is a one-sided boundary maximum on `0 ≤ f < 1`. Re-export
of `capacityTransferArea_lt_baseline`. -/
theorem capacity_transfer_below_baseline {d f : ℝ} (hd : 0 < d)
    (hf0 : 0 < f) (hf1 : f < 1) :
    OPH.DeSitterCapacityShock.capacityTransferArea d f
      < OPH.DeSitterCapacityShock.symmetricArea d :=
  OPH.DeSitterCapacityShock.capacityTransferArea_lt_baseline hd hf0 hf1

/-- **De Sitter radius cancellation.** The pure-de-Sitter shock mass
parameter at spacetime dimension `Dim` equals the `ℓ = 1` spherical
Laplacian eigenvalue; the radius cancels exactly. Re-export of
`pureDeSitterMuSquared_eq_lOne`. -/
theorem deSitter_radius_cancellation (Dim : ℕ) {L : ℝ} (hL : L ≠ 0) :
    OPH.DeSitterCapacityShock.pureDeSitterMuSquared Dim L
      = OPH.DeSitterCapacityShock.lOneSphereEigenvalue Dim :=
  OPH.DeSitterCapacityShock.pureDeSitterMuSquared_eq_lOne Dim hL

end HorizonSurface

end OPH.Thermodynamics

#print axioms OPH.Thermodynamics.HorizonSurface.entanglement_first_law_exact
#print axioms OPH.Thermodynamics.HorizonSurface.entanglement_first_law_first_order
#print axioms OPH.Thermodynamics.HorizonSurface.horizon_cap_clausius
#print axioms OPH.Thermodynamics.HorizonSurface.repair_conserves_central_charge
#print axioms OPH.Thermodynamics.HorizonSurface.repair_displacement_central_pairing_zero
#print axioms OPH.Thermodynamics.HorizonSurface.modular_weight_double_identification
#print axioms OPH.Thermodynamics.HorizonSurface.capacity_dimension_weight_entropy
#print axioms OPH.Thermodynamics.HorizonSurface.capacity_transfer_area_drop
#print axioms OPH.Thermodynamics.HorizonSurface.capacity_transfer_below_baseline
#print axioms OPH.Thermodynamics.HorizonSurface.deSitter_radius_cancellation
