import Mathlib
import A5PrimitivePortPrediction

namespace OPH.A5OrbitRaySeparation

/-!
# Exact separation of three normalized A5 orbit rays

This file is a finite arithmetic certificate for the normalized coefficient
table associated with the vertex-12, face-20, and edge-30 carrier orbits.  All
three rows use the same `C4` and `B0` normalization; their `B6` coefficients,
and hence their two registered rank-six ratios, are distinct.

Boundary: the coefficient rows are inputs to this arithmetic certificate.  No
theorem here derives a row from OPH repair data, selects a propagation support,
or identifies any row with a physical field or sector.
-/

/-- Labels for the three candidate carrier-orbit rows. -/
inductive OrbitRow where
  | vertex12
  | face20
  | edge30
  deriving DecidableEq, Repr

/-- Common normalized fourth-order coefficient. -/
def C4 : ℚ := OPH.A5PrimitivePortPrediction.C4

/-- Common normalized isotropic sixth-order coefficient. -/
def B0 : ℚ := OPH.A5PrimitivePortPrediction.B0

/-- Orbit-dependent normalized rank-six coefficient. -/
def B6 : OrbitRow → ℚ
  | .vertex12 => OPH.A5PrimitivePortPrediction.B6
  | .face20 => -2 / 14175
  | .edge30 => -1 / 12600

/-- The ratio using the square of the common fourth-order coefficient. -/
def b6OverC4Squared (row : OrbitRow) : ℚ := B6 row / C4 ^ 2

/-- The ratio using the common isotropic sixth-order coefficient. -/
def b6OverB0 (row : OrbitRow) : ℚ := B6 row / B0

/-- The shared coefficients have the registered normalization. -/
theorem common_normalization :
    C4 = -1 / 20 ∧ B0 = 1 / 840 ∧ B0 / C4 ^ 2 = 10 / 21 := by
  norm_num [C4, B0, OPH.A5PrimitivePortPrediction.C4,
    OPH.A5PrimitivePortPrediction.B0]

/-- Exact `B6 / C4^2` value for the vertex-12 row. -/
theorem vertex12_b6_over_c4_squared :
    b6OverC4Squared .vertex12 = 32 / 315 := by
  norm_num [b6OverC4Squared, B6, C4,
    OPH.A5PrimitivePortPrediction.B6,
    OPH.A5PrimitivePortPrediction.C4]

/-- Exact `B6 / C4^2` value for the face-20 row. -/
theorem face20_b6_over_c4_squared :
    b6OverC4Squared .face20 = -32 / 567 := by
  norm_num [b6OverC4Squared, B6, C4,
    OPH.A5PrimitivePortPrediction.C4]

/-- Exact `B6 / C4^2` value for the edge-30 row. -/
theorem edge30_b6_over_c4_squared :
    b6OverC4Squared .edge30 = -2 / 63 := by
  norm_num [b6OverC4Squared, B6, C4,
    OPH.A5PrimitivePortPrediction.C4]

/-- Exact `B6 / B0` value for the vertex-12 row. -/
theorem vertex12_b6_over_b0 : b6OverB0 .vertex12 = 16 / 75 := by
  norm_num [b6OverB0, B6, B0,
    OPH.A5PrimitivePortPrediction.B6,
    OPH.A5PrimitivePortPrediction.B0]

/-- Exact `B6 / B0` value for the face-20 row. -/
theorem face20_b6_over_b0 : b6OverB0 .face20 = -16 / 135 := by
  norm_num [b6OverB0, B6, B0,
    OPH.A5PrimitivePortPrediction.B0]

/-- Exact `B6 / B0` value for the edge-30 row. -/
theorem edge30_b6_over_b0 : b6OverB0 .edge30 = -1 / 15 := by
  norm_num [b6OverB0, B6, B0,
    OPH.A5PrimitivePortPrediction.B0]

/-- The three `B6 / C4^2` ratios are pairwise distinct. -/
theorem b6_over_c4_squared_pairwise_distinct :
    b6OverC4Squared .vertex12 ≠ b6OverC4Squared .face20 ∧
      b6OverC4Squared .vertex12 ≠ b6OverC4Squared .edge30 ∧
      b6OverC4Squared .face20 ≠ b6OverC4Squared .edge30 := by
  norm_num [b6OverC4Squared, B6, C4,
    OPH.A5PrimitivePortPrediction.B6,
    OPH.A5PrimitivePortPrediction.C4]

/-- The three `B6 / B0` ratios are pairwise distinct. -/
theorem b6_over_b0_pairwise_distinct :
    b6OverB0 .vertex12 ≠ b6OverB0 .face20 ∧
      b6OverB0 .vertex12 ≠ b6OverB0 .edge30 ∧
      b6OverB0 .face20 ≠ b6OverB0 .edge30 := by
  norm_num [b6OverB0, B6, B0,
    OPH.A5PrimitivePortPrediction.B6,
    OPH.A5PrimitivePortPrediction.B0]

/-- The `B6 / C4^2` ratio uniquely identifies a row within this table. -/
theorem b6_over_c4_squared_injective : Function.Injective b6OverC4Squared := by
  intro row₁ row₂ h
  cases row₁ <;> cases row₂ <;>
    norm_num [b6OverC4Squared, B6, C4,
      OPH.A5PrimitivePortPrediction.B6,
      OPH.A5PrimitivePortPrediction.C4] at h <;> rfl

/-- The `B6 / B0` ratio uniquely identifies a row within this table. -/
theorem b6_over_b0_injective : Function.Injective b6OverB0 := by
  intro row₁ row₂ h
  cases row₁ <;> cases row₂ <;>
    norm_num [b6OverB0, B6, B0,
      OPH.A5PrimitivePortPrediction.B6,
      OPH.A5PrimitivePortPrediction.B0] at h <;> rfl

/-- Both ratio presentations agree through the common normalization. -/
theorem ratio_compatibility (row : OrbitRow) :
    b6OverC4Squared row = b6OverB0 row * (B0 / C4 ^ 2) := by
  cases row <;>
    norm_num [b6OverC4Squared, b6OverB0, B6, B0, C4,
      OPH.A5PrimitivePortPrediction.B6,
      OPH.A5PrimitivePortPrediction.B0,
      OPH.A5PrimitivePortPrediction.C4]

#print axioms common_normalization
#print axioms vertex12_b6_over_c4_squared
#print axioms face20_b6_over_c4_squared
#print axioms edge30_b6_over_c4_squared
#print axioms vertex12_b6_over_b0
#print axioms face20_b6_over_b0
#print axioms edge30_b6_over_b0
#print axioms b6_over_c4_squared_pairwise_distinct
#print axioms b6_over_b0_pairwise_distinct
#print axioms b6_over_c4_squared_injective
#print axioms b6_over_b0_injective
#print axioms ratio_compatibility

end OPH.A5OrbitRaySeparation
