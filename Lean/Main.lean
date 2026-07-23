import ObserverPatchHolography
import ObserverPatchHolography.BridgeEquivalence
import ObserverPatchHolography.CapacityFixedPoint
import ObserverPatchHolography.SeedPi
import EventAlgebra

def main : IO Unit := do
  IO.println "ObserverPatchHolography Lean library loaded."
  IO.println "EventAlgebra library loaded (finite event algebras, sorry-free)."
  IO.println "See docs/PROOF_INDEX.md for current formalisation status."
