import Init

/-!
# Material exchange on the two four-face gateways

The addressed experiment uses four contacts, then repeats them in the same
order. An enabled contact at exposure π/8 swaps its two material occupations,
reverses the written seam head, and contributes unit amplitude. The two origin
labels below follow the paths of the loads. The observer reads occupations and
orbital records at the end. The four head positions distinguish every active
contact spoke. All other contacts are isolated during the sequence.

The experimental generator checks that all sixteen sharp source words restrict
to the two indicated head patterns on the active contacts. The theorem covers
both of these restrictions, both choices of initial bridge head and both
material signs. It establishes the exchange and restoration of the addressed
heads. A physical spin assignment or a general self-exchange phase for a bound
composite needs further hypotheses.
-/

namespace SourceLaw.DressedMaterialExchange

/-- The stored bit is the physical end carrying the head, represented by its port address. -/
structure Contact where
  a : Nat
  b : Nat
  head : Nat
  deriving DecidableEq

structure State where
  first : Nat
  second : Nat
  seams : List Contact
  deriving DecidableEq

/-- Occupations distinguish the two origins for the path audit. At a vacancy the
eligible sending end is opposite the seam head. -/
def pulse (excess : Bool) (st : State) (index : Nat) : State :=
  match st.seams[index]? with
  | none => st
  | some link =>
    let other := if link.head = link.a then link.b else link.a
    let src := if excess then link.head else other
    let dst := if excess then other else link.head
    if st.first = dst ∨ st.second = dst then st
    else if st.first = src then
      { st with first := dst, seams := st.seams.set index { link with head := other } }
    else if st.second = src then
      { st with second := dst, seams := st.seams.set index { link with head := other } }
    else st

def prepared (leftBridgeHead firstBuffer excess : Bool) : State :=
  let centre := if leftBridgeHead then 24 else 14
  let peer := if leftBridgeHead then 14 else 24
  let x := if leftBridgeHead then 35 else 23
  let y := if leftBridgeHead then 34 else 22
  let z := if leftBridgeHead then 29 else 16
  let head (a b : Nat) (eligibleAtA : Bool) := if eligibleAtA = excess then a else b
  { first := 14, second := 24, seams := [
    ⟨centre, x, head centre x firstBuffer⟩,
    ⟨centre, y, head centre y (!firstBuffer)⟩,
    ⟨centre, peer, head centre peer false⟩,
    ⟨centre, z, head centre z true⟩] }

def exchange (excess : Bool) (st : State) : State :=
  [0, 1, 2, 3, 0, 1, 2, 3].foldl (pulse excess) st

/-- The same address sequence for a bridge orientation works for both source
buffer words and their charge conjugates. The final field is exactly the initial field. -/
theorem exchanged_loads_restored_heads (leftBridgeHead firstBuffer excess : Bool) :
    exchange excess (prepared leftBridgeHead firstBuffer excess) =
      { prepared leftBridgeHead firstBuffer excess with first := 24, second := 14 } := by
  cases leftBridgeHead <;> cases firstBuffer <;> cases excess <;> decide

end SourceLaw.DressedMaterialExchange
