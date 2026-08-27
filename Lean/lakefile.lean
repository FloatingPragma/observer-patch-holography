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
lean_lib «OPHTime» where
  srcDir := "."
  roots := #[`Time.TimeOrderLedger, `Time.ObserverHistory,
    `Time.ClockReadout, `Time.WorldlineRealization,
    `Time.ProperTimeCalibration, `Time.ClockComparison]

@[default_target]
lean_lib «OPHComputation» where
  srcDir := "."
  roots := #[`Computation.RepairUniversality]

@[default_target]
lean_lib «OPHThermodynamics» where
  srcDir := "Thermodynamics"
  roots := #[`FiniteConditionalRepair, `FirstLawIdentity,
    `GibbsReferenceEnergyIdentification, `ModularEnergyAdditivity,
    `CollarTemperatureReading,
    `FluctuationTheorems, `CapFirstLaw, `EinsteinPremiseLink,
    `GreenKubo, `GraphDiffusion, `StationaryRealization,
    `PoissonizedRepair, `PoissonizedRepairOperatorExp,
    `LowTemperatureControl, `MixingChainRealization, `CommonObjectBinding,
    `CommonReferenceObstruction, `RepairCurrentOrientation,
    `FourLawAdequacySurface, `HorizonThermalitySurface,
    `CoherentRefinementFamily, `CofinalSpectralTailFamily,
    `PhysicalCalibrationImport]

@[default_target]
lean_lib «OPHScreen» where
  srcDir := "Screen"
  roots := #[`OPHScreen, `Compact12, `S2DesignSignature,
    `TopThreeKernelFix, `UnitSplit12, `Z6Exact, `PhysicalA5ForcingNoGo,
    `PortFrameGram, `A5PortAction, `A5CouplingSymmetry, `A5OPH,
    `A2HolonomyBridge, `HolonomyInterference, `A5CharacterField, `A5SixAxes,
    `A5PortModule, `A5Commutant, `A5ResponseWordAlgebra,
    `A5IncidenceResponse, `TraceBalancedKernel, `TrichotomyCases, `Z6Descent,
    `A5AngularMultiplets, `A5AngularBands, `A5AngularKernels,
    `ExteriorSelection, `EqualStateWeights,
    `A5FamilyBand, `RGRepresentationFrontier, `CommonEWOrderUnit,
    `KineticFormDichotomy, `A5PrimitivePortPrediction,
    `A5CarrierClassBand,
    `A5OrbitRaySeparation, `DiscreteRefinement, `LabeledEventReadout,
    `BaryonDimensionSix,
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
    `SeamCurrentPhotonLeptonThreshold, `CarrierFrequencySpeed,
    `GaugeKineticInvariantForms, `OrientedFaceBracketSelector,
    `OrientedFaceInvariantMetric, `PortDualMetricSelection,
    `LightSignalAdequacySurface, `ModalMaxwellFactorizationBoundary,
    `LightSignalMaxwellComposition, `SeamU1HolonomyClassification,
    `LayeredDiscreteGauss, `DiscreteCoulombGreen, `PositionSpaceMaxwellAction,
    `LocalFaceMaxwellAction, `TemporalMaxwellEvolution,
    `ScaledMaxwellStability, `CertifiedScaledStepInstrument,
    `DispersionArmingInterface, `KogutSusskindFiberRateComparison,
    `GaugeOrbitQuotientGap, `LocalEnergyBalance, `CarrierModeOscillators, `CarrierModeEquivariance, `GoldenSectorCharacters, `GoldenSectorIrreducibility, `CarrierEvolutionFlow, `SeamChargeContinuity, `NeutralPairCoupledAction, `NeutralPairJointStationaryWitness, `GoldenSectorComplexIrreducibility, `FieldSectorEnergyInnerProduct, `CurlSectorEigenbasis, `CurlStoneClockBridge,
    `SMStructureAdequacySurface, `SMStructureComposition,
    `MatterGrammarIndexBridge, `ElectroweakBreakingComposition,
    `AssembledActionComposition, `GlobalFormCharacterDescent,
    `FermionSectorAssembly, `NeutralCurrentDictionary,
    `ChargedCurrentDictionary, `GaugeSectorBracketCompletion,
    `ExteriorComponentBridge, `QuantumMatterIntegration,
    `B10EdgeCenterAction,
    `B10QuantumLimitations, `CarrierUniqueness]

@[default_target]
lean_lib «OPHConstruction» where
  srcDir := "."
  roots := #[`Tower, `Dynamics, `Geometry, `InformationProjection, `Locality,
    `QFT, `Variational]

lean_exe «oph» where
  root := `Main
