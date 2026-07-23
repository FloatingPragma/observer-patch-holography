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
  srcDir := "ObserverPatchHolography/Proofs/ObservableNormalForms"

@[default_target]
lean_lib «ObserverPatchHolography» where
  srcDir := "ObserverPatchHolography/Source"

@[default_target]
lean_lib «EventAlgebra» where
  srcDir := "ObserverPatchHolography/Source"

@[default_target]
lean_lib «OPHScreen» where
  srcDir := "ObserverPatchHolography/Screen"
  roots := #[`OPHScreen, `Compact12, `S2DesignSignature,
    `TopThreeKernelFix, `UnitSplit12, `Z6Exact, `PhysicalA5ForcingNoGo,
    `PortFrameGram, `A5PortAction, `A5CouplingSymmetry, `A5OPH,
    `A5CharacterField, `A5SixAxes]

lean_exe «oph» where
  root := `Main
