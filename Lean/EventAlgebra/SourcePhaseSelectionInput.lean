import Mathlib

set_option autoImplicit false

/-!
# Generated source phase-selection input

This module is mechanically generated from the committed target-free source
packet.  It contains source states, source operations, the complete
state/event cell table, the enabled-domain split, and custody hashes only.
-/

namespace EventAlgebra.SourcePhaseSelectionInput

inductive PhaseState where
  | orbit00
  | orbit02
  | orbit04
  | phaseNegative
  | phasePositive
  deriving DecidableEq, Fintype, Repr

inductive PhaseEvent where
  | pair0002
  | pair0003
  | pair0004
  | pair0005
  | pair0102
  | pair0103
  | pair0104
  | pair0105
  | pair0204
  | pair0205
  | pair0304
  | pair0305
  deriving DecidableEq, Fintype, Repr

abbrev StateBytes := List Nat
abbrev OperationBytes := List Nat

def stateBytes : PhaseState → StateBytes
  | .orbit00 => [91, 91, 91, 34, 49, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]
  | .orbit02 => [91, 91, 91, 34, 49, 47, 52, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 45, 49, 47, 52, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 45, 49, 47, 52, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 51, 47, 52, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]
  | .orbit04 => [91, 91, 91, 34, 49, 47, 52, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 49, 47, 52, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 49, 47, 52, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 51, 47, 52, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]
  | .phaseNegative => [91, 91, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 49, 47, 50, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 45, 49, 47, 50, 34, 44, 34, 48, 34, 93, 44, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]
  | .phasePositive => [91, 91, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 45, 49, 47, 50, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 49, 47, 50, 34, 44, 34, 48, 34, 93, 44, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]

def operationBytes : PhaseEvent → OperationBytes
  | .pair0002 => [91, 91, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 45, 49, 47, 50, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 49, 47, 50, 34, 44, 34, 48, 34, 93, 44, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]
  | .pair0003 => [91, 91, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 45, 49, 47, 50, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 49, 47, 50, 34, 44, 34, 48, 34, 93, 44, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]
  | .pair0004 => [91, 91, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 49, 47, 50, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 45, 49, 47, 50, 34, 44, 34, 48, 34, 93, 44, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]
  | .pair0005 => [91, 91, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 49, 47, 50, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 45, 49, 47, 50, 34, 44, 34, 48, 34, 93, 44, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]
  | .pair0102 => [91, 91, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 45, 49, 47, 50, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 49, 47, 50, 34, 44, 34, 48, 34, 93, 44, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]
  | .pair0103 => [91, 91, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 45, 49, 47, 50, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 49, 47, 50, 34, 44, 34, 48, 34, 93, 44, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]
  | .pair0104 => [91, 91, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 49, 47, 50, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 45, 49, 47, 50, 34, 44, 34, 48, 34, 93, 44, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]
  | .pair0105 => [91, 91, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 49, 47, 50, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 45, 49, 47, 50, 34, 44, 34, 48, 34, 93, 44, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]
  | .pair0204 => [91, 91, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 45, 49, 47, 50, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 49, 47, 50, 34, 44, 34, 48, 34, 93, 44, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]
  | .pair0205 => [91, 91, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 45, 49, 47, 50, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 49, 47, 50, 34, 44, 34, 48, 34, 93, 44, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]
  | .pair0304 => [91, 91, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 45, 49, 47, 50, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 49, 47, 50, 34, 44, 34, 48, 34, 93, 44, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]
  | .pair0305 => [91, 91, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 44, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 45, 49, 47, 50, 34, 44, 34, 48, 34, 93, 93, 44, 91, 91, 34, 48, 34, 44, 34, 48, 34, 44, 34, 49, 47, 50, 34, 44, 34, 48, 34, 93, 44, 91, 34, 49, 47, 50, 34, 44, 34, 48, 34, 44, 34, 48, 34, 44, 34, 48, 34, 93, 93, 93, 10]

inductive CellStatus where
  | enabled
  | disabled
  deriving DecidableEq, Repr

abbrev SourceCell := PhaseState × PhaseEvent

def cellStatus : PhaseState → PhaseEvent → CellStatus
  | .orbit00, .pair0002 => .enabled
  | .orbit00, .pair0003 => .enabled
  | .orbit00, .pair0004 => .enabled
  | .orbit00, .pair0005 => .enabled
  | .orbit00, .pair0102 => .enabled
  | .orbit00, .pair0103 => .enabled
  | .orbit00, .pair0104 => .enabled
  | .orbit00, .pair0105 => .enabled
  | .orbit00, .pair0204 => .enabled
  | .orbit00, .pair0205 => .enabled
  | .orbit00, .pair0304 => .enabled
  | .orbit00, .pair0305 => .enabled
  | .orbit02, .pair0002 => .enabled
  | .orbit02, .pair0003 => .enabled
  | .orbit02, .pair0004 => .enabled
  | .orbit02, .pair0005 => .enabled
  | .orbit02, .pair0102 => .enabled
  | .orbit02, .pair0103 => .enabled
  | .orbit02, .pair0104 => .enabled
  | .orbit02, .pair0105 => .enabled
  | .orbit02, .pair0204 => .enabled
  | .orbit02, .pair0205 => .enabled
  | .orbit02, .pair0304 => .enabled
  | .orbit02, .pair0305 => .enabled
  | .orbit04, .pair0002 => .enabled
  | .orbit04, .pair0003 => .enabled
  | .orbit04, .pair0004 => .enabled
  | .orbit04, .pair0005 => .enabled
  | .orbit04, .pair0102 => .enabled
  | .orbit04, .pair0103 => .enabled
  | .orbit04, .pair0104 => .enabled
  | .orbit04, .pair0105 => .enabled
  | .orbit04, .pair0204 => .enabled
  | .orbit04, .pair0205 => .enabled
  | .orbit04, .pair0304 => .enabled
  | .orbit04, .pair0305 => .enabled
  | .phaseNegative, .pair0002 => .disabled
  | .phaseNegative, .pair0003 => .disabled
  | .phaseNegative, .pair0004 => .enabled
  | .phaseNegative, .pair0005 => .enabled
  | .phaseNegative, .pair0102 => .disabled
  | .phaseNegative, .pair0103 => .disabled
  | .phaseNegative, .pair0104 => .enabled
  | .phaseNegative, .pair0105 => .enabled
  | .phaseNegative, .pair0204 => .disabled
  | .phaseNegative, .pair0205 => .disabled
  | .phaseNegative, .pair0304 => .disabled
  | .phaseNegative, .pair0305 => .disabled
  | .phasePositive, .pair0002 => .enabled
  | .phasePositive, .pair0003 => .enabled
  | .phasePositive, .pair0004 => .disabled
  | .phasePositive, .pair0005 => .disabled
  | .phasePositive, .pair0102 => .enabled
  | .phasePositive, .pair0103 => .enabled
  | .phasePositive, .pair0104 => .disabled
  | .phasePositive, .pair0105 => .disabled
  | .phasePositive, .pair0204 => .enabled
  | .phasePositive, .pair0205 => .enabled
  | .phasePositive, .pair0304 => .enabled
  | .phasePositive, .pair0305 => .enabled

def allSourceCells : List SourceCell :=
  [(.orbit00, .pair0002), (.orbit00, .pair0003), (.orbit00, .pair0004), (.orbit00, .pair0005), (.orbit00, .pair0102), (.orbit00, .pair0103), (.orbit00, .pair0104), (.orbit00, .pair0105), (.orbit00, .pair0204), (.orbit00, .pair0205), (.orbit00, .pair0304), (.orbit00, .pair0305), (.orbit02, .pair0002), (.orbit02, .pair0003), (.orbit02, .pair0004), (.orbit02, .pair0005), (.orbit02, .pair0102), (.orbit02, .pair0103), (.orbit02, .pair0104), (.orbit02, .pair0105), (.orbit02, .pair0204), (.orbit02, .pair0205), (.orbit02, .pair0304), (.orbit02, .pair0305), (.orbit04, .pair0002), (.orbit04, .pair0003), (.orbit04, .pair0004), (.orbit04, .pair0005), (.orbit04, .pair0102), (.orbit04, .pair0103), (.orbit04, .pair0104), (.orbit04, .pair0105), (.orbit04, .pair0204), (.orbit04, .pair0205), (.orbit04, .pair0304), (.orbit04, .pair0305), (.phaseNegative, .pair0002), (.phaseNegative, .pair0003), (.phaseNegative, .pair0004), (.phaseNegative, .pair0005), (.phaseNegative, .pair0102), (.phaseNegative, .pair0103), (.phaseNegative, .pair0104), (.phaseNegative, .pair0105), (.phaseNegative, .pair0204), (.phaseNegative, .pair0205), (.phaseNegative, .pair0304), (.phaseNegative, .pair0305), (.phasePositive, .pair0002), (.phasePositive, .pair0003), (.phasePositive, .pair0004), (.phasePositive, .pair0005), (.phasePositive, .pair0102), (.phasePositive, .pair0103), (.phasePositive, .pair0104), (.phasePositive, .pair0105), (.phasePositive, .pair0204), (.phasePositive, .pair0205), (.phasePositive, .pair0304), (.phasePositive, .pair0305)]

abbrev EnabledCell :=
  {cell : SourceCell // cellStatus cell.1 cell.2 = .enabled}

abbrev DisabledCell :=
  {cell : SourceCell // cellStatus cell.1 cell.2 = .disabled}

def enabledCells : List EnabledCell :=
  [⟨(.orbit00, .pair0002), by rfl⟩, ⟨(.orbit00, .pair0003), by rfl⟩, ⟨(.orbit00, .pair0004), by rfl⟩, ⟨(.orbit00, .pair0005), by rfl⟩, ⟨(.orbit00, .pair0102), by rfl⟩, ⟨(.orbit00, .pair0103), by rfl⟩, ⟨(.orbit00, .pair0104), by rfl⟩, ⟨(.orbit00, .pair0105), by rfl⟩, ⟨(.orbit00, .pair0204), by rfl⟩, ⟨(.orbit00, .pair0205), by rfl⟩, ⟨(.orbit00, .pair0304), by rfl⟩, ⟨(.orbit00, .pair0305), by rfl⟩, ⟨(.orbit02, .pair0002), by rfl⟩, ⟨(.orbit02, .pair0003), by rfl⟩, ⟨(.orbit02, .pair0004), by rfl⟩, ⟨(.orbit02, .pair0005), by rfl⟩, ⟨(.orbit02, .pair0102), by rfl⟩, ⟨(.orbit02, .pair0103), by rfl⟩, ⟨(.orbit02, .pair0104), by rfl⟩, ⟨(.orbit02, .pair0105), by rfl⟩, ⟨(.orbit02, .pair0204), by rfl⟩, ⟨(.orbit02, .pair0205), by rfl⟩, ⟨(.orbit02, .pair0304), by rfl⟩, ⟨(.orbit02, .pair0305), by rfl⟩, ⟨(.orbit04, .pair0002), by rfl⟩, ⟨(.orbit04, .pair0003), by rfl⟩, ⟨(.orbit04, .pair0004), by rfl⟩, ⟨(.orbit04, .pair0005), by rfl⟩, ⟨(.orbit04, .pair0102), by rfl⟩, ⟨(.orbit04, .pair0103), by rfl⟩, ⟨(.orbit04, .pair0104), by rfl⟩, ⟨(.orbit04, .pair0105), by rfl⟩, ⟨(.orbit04, .pair0204), by rfl⟩, ⟨(.orbit04, .pair0205), by rfl⟩, ⟨(.orbit04, .pair0304), by rfl⟩, ⟨(.orbit04, .pair0305), by rfl⟩, ⟨(.phaseNegative, .pair0004), by rfl⟩, ⟨(.phaseNegative, .pair0005), by rfl⟩, ⟨(.phaseNegative, .pair0104), by rfl⟩, ⟨(.phaseNegative, .pair0105), by rfl⟩, ⟨(.phasePositive, .pair0002), by rfl⟩, ⟨(.phasePositive, .pair0003), by rfl⟩, ⟨(.phasePositive, .pair0102), by rfl⟩, ⟨(.phasePositive, .pair0103), by rfl⟩, ⟨(.phasePositive, .pair0204), by rfl⟩, ⟨(.phasePositive, .pair0205), by rfl⟩, ⟨(.phasePositive, .pair0304), by rfl⟩, ⟨(.phasePositive, .pair0305), by rfl⟩]

def disabledCells : List DisabledCell :=
  [⟨(.phaseNegative, .pair0002), by rfl⟩, ⟨(.phaseNegative, .pair0003), by rfl⟩, ⟨(.phaseNegative, .pair0102), by rfl⟩, ⟨(.phaseNegative, .pair0103), by rfl⟩, ⟨(.phaseNegative, .pair0204), by rfl⟩, ⟨(.phaseNegative, .pair0205), by rfl⟩, ⟨(.phaseNegative, .pair0304), by rfl⟩, ⟨(.phaseNegative, .pair0305), by rfl⟩, ⟨(.phasePositive, .pair0004), by rfl⟩, ⟨(.phasePositive, .pair0005), by rfl⟩, ⟨(.phasePositive, .pair0104), by rfl⟩, ⟨(.phasePositive, .pair0105), by rfl⟩]

theorem sourceCellCensus : allSourceCells.length = 60 := by decide
theorem enabledCellCensus : enabledCells.length = 48 := by decide
theorem disabledCellCensus : disabledCells.length = 12 := by decide

def sourceBoundEnabledCell : EnabledCell := enabledCells[0]
def sourceBoundDisabledCell : DisabledCell := disabledCells[0]

def sourcePayloadSha256 : String := "71a06f1c15192123cd09feb2386da702b572c8ac57c9b7633f5aa60c5d404e22"
def sourceSelectionPacketSha256 : String := "3412b8fa528635bf670e8d3ba7a1a68558d63524a9795e4603ca5148b6617970"
def sourceViewSha256 : String := "06970219d5df9fc0783cf7576020f9796945a44eb433b31918587dd81a09a6ac"

end EventAlgebra.SourcePhaseSelectionInput
