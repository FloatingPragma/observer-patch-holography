import QFT.SourcePhaseLiftBridge
import EventAlgebra.SourcePhaseSelectionInput

set_option autoImplicit false

/-! Mechanically generated from the committed target-free source packet. -/
namespace EventAlgebra.SourcePhaseSelectionSemantics

open Matrix
open EventAlgebra.SourcePhaseSelectionInput

noncomputable section

inductive EffectValue where
  | positive
  | negative
  deriving DecidableEq, Fintype, Repr

def stateMatrix : PhaseState → Matrix (Fin 2) (Fin 2) ℂ
  | .orbit00 => !![(((((1 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)), (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)); (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)), (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I))]
  | .orbit02 => !![((((((1 / 4 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)), (((((0 : ℝ) + ((-1 / 4 : ℝ) : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)); (((((0 : ℝ) + ((-1 / 4 : ℝ) : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)), ((((((3 / 4 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I))]
  | .orbit04 => !![((((((1 / 4 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)), (((((0 : ℝ) + ((1 / 4 : ℝ) : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)); (((((0 : ℝ) + ((1 / 4 : ℝ) : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)), ((((((3 / 4 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I))]
  | .phaseNegative => !![((((((1 / 2 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)), (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + ((((((1 / 2 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)); (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + ((((((-1 / 2 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)), ((((((1 / 2 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I))]
  | .phasePositive => !![((((((1 / 2 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)), (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + ((((((-1 / 2 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)); (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + ((((((1 / 2 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)), ((((((1 / 2 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I))]

def effectMatrix : EffectValue → Matrix (Fin 2) (Fin 2) ℂ
  | .positive => !![((((((1 / 2 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)), (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + ((((((-1 / 2 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)); (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + ((((((1 / 2 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)), ((((((1 / 2 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I))]
  | .negative => !![((((((1 / 2 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)), (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + ((((((1 / 2 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)); (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + ((((((-1 / 2 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I)), ((((((1 / 2 : ℝ) : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) + (((((0 : ℝ) + (0 : ℝ) * OPH.QFT.sqrt3)) : ℂ) * Complex.I))]

def effectValueOf : PhaseEvent → EffectValue
  | .pair0002 => .positive
  | .pair0003 => .positive
  | .pair0004 => .negative
  | .pair0005 => .negative
  | .pair0102 => .positive
  | .pair0103 => .positive
  | .pair0104 => .negative
  | .pair0105 => .negative
  | .pair0204 => .positive
  | .pair0205 => .positive
  | .pair0304 => .positive
  | .pair0305 => .positive

def generatedEffect (event : PhaseEvent) : Matrix (Fin 2) (Fin 2) ℂ :=
  effectMatrix (effectValueOf event)

def sourcePairOrder : List PhaseEvent :=
  [.pair0002, .pair0003, .pair0004, .pair0005, .pair0102, .pair0103, .pair0104, .pair0105, .pair0204, .pair0205, .pair0304, .pair0305]

def sourceSelectedEvent : PhaseEvent := sourcePairOrder[0]

def sourceSelectedGeneratedEffect : Matrix (Fin 2) (Fin 2) ℂ :=
  generatedEffect sourceSelectedEvent

def generatedEffectValues : List EffectValue :=
  (sourcePairOrder.map effectValueOf).eraseDups

def upperImaginaryNegative (value : EffectValue) : Prop :=
  Complex.im (effectMatrix value 0 1) < 0

def sourcePayloadSha256 : String :=
  "71a06f1c15192123cd09feb2386da702b572c8ac57c9b7633f5aa60c5d404e22"
def sourceSelectionPacketSha256 : String :=
  "3412b8fa528635bf670e8d3ba7a1a68558d63524a9795e4603ca5148b6617970"

#check stateMatrix
#check effectMatrix
#check generatedEffect
#check sourceSelectedGeneratedEffect

end

end EventAlgebra.SourcePhaseSelectionSemantics
