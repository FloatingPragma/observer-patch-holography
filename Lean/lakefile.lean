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
lean_lib «OPHScreen» where
  srcDir := "Screen"
  roots := #[`OPHScreen, `Compact12, `S2DesignSignature,
    `TopThreeKernelFix, `UnitSplit12, `Z6Exact, `PhysicalA5ForcingNoGo,
    `PortFrameGram, `A5PortAction, `A5CouplingSymmetry, `A5OPH,
    `A5CharacterField, `A5SixAxes, `A5PortModule, `A5Commutant,
    `A5IncidenceResponse, `TraceBalancedKernel, `TrichotomyCases, `Z6Descent,
    `A5AngularMultiplets, `A5AngularBands, `ExteriorSelection, `EqualStateWeights,
    `A5FamilyBand, `RGRepresentationFrontier, `CommonEWOrderUnit]

lean_exe «oph» where
  root := `Main
