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

-- The frozen 8,662-carrier certificate is deliberately opt-in.  Lean CI
-- builds it in a dedicated, change-aware job so unrelated default builds do
-- not inherit its sequential kernel-evaluation cost.
lean_lib «OphGap» where
  srcDir := "."
  -- The folded proof predates this repository's project-wide explicit-binder
  -- policy. Keep its elaboration environment unchanged during the Mathlib bump.
  leanOptions := #[
    ⟨`autoImplicit, true⟩,
    ⟨`relaxedAutoImplicit, true⟩
  ]

@[default_target]
lean_lib «EventAlgebra» where
  srcDir := "."

@[default_target]
lean_lib «OPHTime» where
  srcDir := "."
  roots := #[`Time.TimeOrderLedger, `Time.ObserverHistory,
    `Time.ClockReadout, `Time.SourceCountClock, `Time.SourceCountClockEnclosure,
    `Time.WorldlineRealization,
    `Time.ProperTimeCalibration, `Time.ClockComparison]

@[default_target]
lean_lib «OPHComputation» where
  srcDir := "."
  roots := #[`Computation.RepairUniversality,
    `Computation.FixedFederationProgress,
    `Computation.FixedFederationCounterexamples,
    `Computation.FixedFederationExecution,
    `Computation.FixedFederationExecutionExamples,
    `Computation.FixedFederationComplexity,
    `Computation.FixedFederationFanoutControl]

@[default_target]
lean_lib «OPHThermodynamics» where
  srcDir := "Thermodynamics"
  roots := #[`FiniteConditionalRepair, `FirstLawIdentity, `EquationOfState,
    `GibbsReferenceEnergyIdentification, `ModularEnergyAdditivity,
    `CollarTemperatureReading,
    `FluctuationTheorems, `CapFirstLaw, `EinsteinPremiseLink,
    `GreenKubo, `NativeRepairFourCycle, `TemporalFilterMonotonicity, `GraphDiffusion, `StationaryRealization,
    `PoissonizedRepair, `PoissonizedRepairOperatorExp,
    `LowTemperatureControl, `MixingChainRealization, `CommonObjectBinding,
    `CommonReferenceObstruction, `RepairCurrentOrientation,
    `FourLawAdequacySurface, `GibbsEntropyLimit, `HorizonThermalitySurface,
    `CoherentRefinementFamily, `CofinalSpectralTailFamily,
    `ProtectedGaugeCorner,
    `PhysicalCalibrationImport,
    `KolmogorovCriterion]

@[default_target]
lean_lib «OPHScreen» where
  srcDir := "Screen"
  roots := #[`OPHScreen, `Compact12, `S2DesignSignature,
    `TopThreeKernelFix, `UnitSplit12, `Z6Exact, `PhysicalA5ForcingNoGo,
    `PortFrameGram, `A5PortAction, `A5CouplingSymmetry, `A5OPH,
    `A2HolonomyBridge, `HolonomyInterference, `A5CharacterField, `A5SixAxes,
    `PSL2F5SixAxesBridge, `GoldenSectorPSL2F5Representations, `SL2F5PortCover,
    `A5PortModule, `A5PortSixAxesBridge, `A5Commutant, `A5ResponseWordAlgebra,
    `A5IncidenceResponse, `TraceBalancedKernel, `TrichotomyCases, `Z6Descent,
    `A5AngularMultiplets, `A5AngularBands, `A5AngularKernels,
    `ExteriorSelection, `WeylYukawaConventions, `ObserverFieldInterfaces, `QuarkMeanDomain, `EqualStateWeights,
    `A5FamilyBand, `RGRepresentationFrontier, `CommonEWOrderUnit,
    `KineticFormDichotomy, `A5PrimitivePortPrediction,
    `A5CarrierClassBand,
    `A5OrbitRaySeparation, `DiscreteRefinement, `LabeledEventReadout, `SerialMaxwellReadout, `CommonHistoryFeedback,
    `ConeCochainBridge, `WhitneyTimeBridge, `WhitneyMaxwellDynamics, `WhitneyQuantumBridge,
    `WhitneyNormalModeConstruction, `WhitneyConeModes, `WhitneyFrameNaturality,
    `WhitneyFrameMorphism, `WhitneyConeNaturality, `WhitneyConeMass,
    `WhitneyConeMassNaturality, `WhitneyMassPremiseInstance, `WhitneyRotationWitness,
    `WhitneyChargedMatter, `WhitneySpatialConsistency,
    `WhitneyFiniteCertificate, `WhitneyAlgebraicLDL,
    `WhitneyCertificatePipeline, `WhitneyCertificateCLI,
    `WhitneyExactSourceProblem, `WhitneySourceProblemCLI,
    `WhitneyGeneratedCertificate,
    `WhitneyActualCertificates, `WhitneySourceGeometry, `WhitneySourceAssembly,
    `WhitneyCurlRow0, `WhitneyCurlRow1, `WhitneyCurlRow2, `WhitneyCurlRow3,
    `WhitneyCertifiedConsumers, `WhitneySourceNaturality,
    `WhitneyOmittedCellCounterexample,
    `MagneticMatterStability,
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
    `SeamCurrentPhotonLeptonThreshold, `CarrierFrequencySpeed, `SeamMaxwellContinuum,
    `GaugeKineticInvariantForms, `OrientedFaceBracketSelector,
    `OrientedFaceInvariantMetric, `PortDualMetricSelection,
    `LightSignalAdequacySurface, `ModalMaxwellFactorizationBoundary,
    `LightSignalMaxwellComposition, `SeamU1HolonomyClassification,
    `LayeredDiscreteGauss, `DiscreteCoulombGreen, `PositionSpaceMaxwellAction,
    `LocalFaceMaxwellAction, `TemporalMaxwellEvolution,
    `ScaledMaxwellStability, `CertifiedScaledStepInstrument,
    `DispersionArmingInterface, `KogutSusskindFiberRateComparison,
    `GaugeOrbitQuotientGap, `LocalEnergyBalance, `CarrierModeOscillators, `CarrierModeEquivariance, `GoldenSectorCharacters, `GoldenSectorIrreducibility, `CarrierEvolutionFlow, `SeamChargeContinuity, `NeutralPairCoupledAction, `NeutralPairJointStationaryWitness, `GoldenSectorComplexIrreducibility, `FieldSectorEnergyInnerProduct, `CurlSectorEigenbasis, `CurlStoneClockBridge, `ActionTimeGram, `SourceActionTime,
    `SMStructureAdequacySurface, `SMStructureComposition,
    `MatterGrammarIndexBridge, `ElectroweakBreakingComposition,
    `AssembledActionComposition, `LocalGaugeJetAction, `CartanScalarReduction, `GlobalFormCharacterDescent,
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

@[default_target]
lean_lib «SourceLaw» where
  srcDir := "."
  roots := #[`SourceLaw.WeightedSeams, `SourceLaw.HistoryRegister,
    `SourceLaw.CounterfamilyLocality, `SourceLaw.IdentificationLadder,
    `SourceLaw.LangevinStationarity, `SourceLaw.TransferStructure,
    `SourceLaw.ReflectionPositivity, `SourceLaw.ContinuumTimeLimit,
    `SourceLaw.TwoPointFlux, `SourceLaw.MatrixSeams,
    `SourceLaw.ScalarPacket, `SourceLaw.TransferMatrixLift,
    `SourceLaw.RecordRealization, `SourceLaw.SpinMarginalForce,
    `SourceLaw.RepairAdmissibility, `SourceLaw.CounterfamilyReader,
    `SourceLaw.StationaryWindow, `SourceLaw.GoldenCarrier,
    `SourceLaw.MatrixLimits, `SourceLaw.DynamicalResponse,
    `SourceLaw.GoldenFamily, `SourceLaw.ContinuumConnection,
    `SourceLaw.RepairConvergence, `SourceLaw.OpenWindow,
    `SourceLaw.VacuumIdentity]

lean_exe «oph» where
  root := `Main
