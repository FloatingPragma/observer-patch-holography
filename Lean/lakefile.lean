import Lake
open Lake DSL

require "leanprover-community" / "mathlib" @ git "v4.29.1"

package «ObserverPatchHolography» where
  leanOptions := #[
    ⟨`autoImplicit, false⟩,
    ⟨`relaxedAutoImplicit, false⟩,
    ⟨`pp.unicode.fun, true⟩
  ]

@[default_target]
lean_lib «ObservableNormalForms» where
  srcDir := "ObservableNormalForms"

@[default_target]
lean_lib «ObserverPatchHolography» where
  srcDir := "."

@[default_target]
lean_lib «EventAlgebra» where
  srcDir := "."

@[default_target]
lean_lib «OPHThermodynamics» where
  srcDir := "Thermodynamics"
  roots := #[`FiniteConditionalRepair, `FirstLawIdentity,
    `FluctuationTheorems, `CapFirstLaw, `EinsteinPremiseLink]

@[default_target]
lean_lib «OPHScreen» where
  srcDir := "Screen"
  roots := #[`OPHScreen, `Compact12, `S2DesignSignature,
    `TopThreeKernelFix, `UnitSplit12, `Z6Exact, `PhysicalA5ForcingNoGo,
    `PortFrameGram, `A5PortAction, `A5CouplingSymmetry, `A5OPH,
    `A2HolonomyBridge, `A5CharacterField, `A5SixAxes, `A5PortModule, `A5Commutant,
    `A5IncidenceResponse, `TraceBalancedKernel, `TrichotomyCases, `Z6Descent,
    `A5AngularMultiplets, `A5AngularBands, `A5AngularKernels,
    `ExteriorSelection, `EqualStateWeights,
    `A5FamilyBand, `RGRepresentationFrontier, `CommonEWOrderUnit,
    `KineticFormDichotomy, `A5PrimitivePortPrediction,
    `A5CarrierClassBand,
    `A5OrbitRaySeparation, `DiscreteRefinement, `BaryonDimensionSix,
    `BipoSHTransferInvariant, `BipoSHInverseBoundary, `BipoSHFrameInvariant,
    `VolumeReadoutBridge, `PrimitivePortTranslationBridge,
    `PrimitivePortScaleBoundary, `PrimitivePortMetricAttachment,
    `PrimitivePortOperatorSelectionBoundary, `KineticFamilyCancellation,
    `ElectroweakPoleScaleQuotient,
    `PrimitivePortDualMeasure, `PrimitiveHopSelection,
    `PrimitivePortFrameQuotient, `PortGramRepairBand,
    `PortGramRepairCovariance, `PortGramA5Isometry,
    `RepairWordCarrierReadout, `SeamCurrentCarrierQuotient,
    `RegionalContinuity, `DiscreteGauss,
    `SeamCurrentEdge30Moment, `SeamCurrentEdge30Remainder,
    `SeamCurrentHomogeneousAction,
    `SeamCurrentDirichletGenerator, `SeamCurrentAuxiliaryOscillatorLift,
    `SeamCurrentFreePhotonLift, `SeamCurrentPhysicalMetricAttachment,
    `SeamCurrentPhotonLeptonThreshold]

@[default_target]
lean_lib «OPHConstruction» where
  srcDir := "."
  roots := #[`Dynamics]

lean_exe «oph» where
  root := `Main
