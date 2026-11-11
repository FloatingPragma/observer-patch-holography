import EventAlgebra.NativeBinarySourceSemantics

set_option maxRecDepth 4096
set_option maxHeartbeats 2000000

namespace EventAlgebra.NativeBinarySource.Endpoints

-- The exact source function, with expressions decoded inside Lean.
def sourceText : String := r#"def native(loads, seam, coin):
    a, b = seam
    d = loads[a] - loads[b]
    current = (d - (2 * coin - 1) * (abs(d) % 2)) // 2
    out = list(loads)
    out[a] -= current
    out[b] += current
    return tuple(out), (coin, current)"#

theorem source_decodes : (parseProgram sourceText).isSome = true := by decide

def runSource (loads : List Int) (a b : Nat) (coin : Int) := do
  let program ← parseProgram sourceText
  evalProgram program loads a b coin

def endpointA (x y coin : Int) : Int := if x = y then x else coin
def endpointB (x y coin : Int) : Int := if x = y then y else 1-coin

-- All binary endpoint/coin cases, with an arbitrary unchanged list tail.
-- Negative d/current and zero-current branches occur among these cases.
theorem source_binary_endpoint (x y coin : Int) (tail : List Int)
    (hx : x = 0 ∨ x = 1) (hy : y = 0 ∨ y = 1) (hc : coin = 0 ∨ coin = 1) :
    runSource (x::y::tail) 0 1 coin =
      some (endpointA x y coin :: endpointB x y coin :: tail,
            (coin, x-endpointA x y coin)) := by
  rcases hx with rfl | rfl <;> rcases hy with rfl | rfl <;>
    rcases hc with rfl | rfl <;> rfl

-- Python // agrees with Euclidean division at the supported positive divisors.
theorem negative_floor_case : (-3 : Int) / 2 = -2 ∧ (-3 : Int) % 2 = 1 := by decide

#print axioms source_decodes
#print axioms source_binary_endpoint
#print axioms negative_floor_case

end EventAlgebra.NativeBinarySource.Endpoints
