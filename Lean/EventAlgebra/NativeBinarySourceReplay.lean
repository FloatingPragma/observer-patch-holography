import EventAlgebra.NativeBinarySourceSemantics

open EventAlgebra.NativeBinarySource Lean

namespace EventAlgebra.NativeBinarySource.Replay

def requireIO (condition : Bool) (message : String) : IO Unit :=
  unless condition do throw (IO.userError message)

def fromExcept {α : Type} (result : Except String α) : IO α :=
  match result with
  | .ok value => pure value
  | .error error => throw (IO.userError error)

def main (args : List String) : IO Unit := do
  let [sourceFile, carrierFile] := args |
    throw (IO.userError "expected source.py carrier.json")
  let source ← IO.FS.readFile sourceFile
  let text := extractNative source
  let some program := parseProgram text |
    throw (IO.userError "native function is outside the supported source grammar")
  let data ← fromExcept (Json.parse (← IO.FS.readFile carrierFile))
  let states : Array (Array Int) ← fromExcept (data.getObjValAs? (Array (Array Int)) "states")
  let seams : Array (Array Nat) ← fromExcept (data.getObjValAs? (Array (Array Nat)) "seams")
  requireIO (states.size == 66 && seams.size == 30) "carrier/schedule shape mismatch"
  let mut counts := Array.replicate 66 (Array.replicate 66 (0 : Nat))
  let mut branches : Array Json := #[]
  for i in [:states.size] do
    let loads := states[i]!.toList
    for e in [:seams.size] do
      let seam := seams[e]!
      requireIO (seam.size == 2) "bad seam"
      for coin in [0, 1] do
        let some (final, (recordCoin, current)) :=
          evalProgram program loads seam[0]! seam[1]! coin |
          throw (IO.userError "source evaluator failed")
        let target := states.findIdx? (·.toList == final)
        let some j := target | throw (IO.userError "source left binary carrier")
        requireIO (recordCoin == coin) "recorded coin drift"
        let row := counts[i]!
        counts := counts.set! i (row.set! j (row[j]! + 1))
        branches := branches.push (toJson #[toJson i, toJson e, toJson coin,
                                           toJson j, toJson recordCoin, toJson current])
  let result := Json.mkObj [("native_text", toJson text), ("states", toJson states),
                            ("seams", toJson seams), ("counts", toJson counts),
                            ("branches", toJson branches)]
  IO.println result.compress

end EventAlgebra.NativeBinarySource.Replay
