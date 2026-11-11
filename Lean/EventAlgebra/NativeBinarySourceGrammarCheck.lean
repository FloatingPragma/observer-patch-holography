import EventAlgebra.NativeBinarySourceSemantics

open EventAlgebra.NativeBinarySource Lean

namespace EventAlgebra.NativeBinarySource.GrammarCheck

def main (args : List String) : IO Unit := do
  let [path] := args | throw (IO.userError "expected grammar-cases.json")
  let .ok json := Json.parse (← IO.FS.readFile path) |
    throw (IO.userError "invalid input JSON")
  let .ok (cases : Array (Array String)) := fromJson? json |
    throw (IO.userError "expected string pairs")
  let mut results : Array Json := #[]
  for row in cases do
    if row.size != 2 then throw (IO.userError "expected name/source pair")
    results := results.push (Json.mkObj [("name",toJson row[0]!),
      ("source_text",toJson row[1]!),
      ("grammar_admitted",toJson (parseProgram row[1]!).isSome)])
  IO.println (toJson results).compress

end EventAlgebra.NativeBinarySource.GrammarCheck
