import Std

/-!
# Connected neutral words of the guarded Interval Echo cell

Two excess ports and two vacant ports in the executable cell's labeling.
The table entries are doubled energies relative to the filled common host.
A row sum of the connected compression subtracts one for every allowed
off-diagonal exchange that stays inside the connected subspace.
-/
namespace OPH.NeutralClusterEnergy

def seams : List (Nat × Nat) :=
  [(0,1),(0,5),(0,7),(0,10),(0,11),(1,5),(1,7),(1,8),(1,9),
   (2,3),(2,4),(2,6),(2,10),(2,11),(3,4),(3,6),(3,8),(3,9),
   (4,5),(4,9),(4,11),(5,9),(5,11),(6,7),(6,8),(6,10),(7,8),
   (7,10),(8,9),(10,11)]

/-- Two disjoint, independently sorted positive and negative port pairs. -/
def pairs : List (Nat × Nat) := (List.range 12).flatMap fun a =>
  (List.range a).map fun b => (a,b)

def mask (p : Nat × Nat) : Nat := 2^p.1 + 2^p.2

/-- Bit s caches reachability on the four-site support s. The row proof checks
each legal bit against the independent four-vertex adjacency criterion. -/
def onCluster : Nat := 0x10001000000000000000000000000000000000001000000000000000100000000000000000000000000000001000000010001001600000001000101120001001601121028000000000000000000000000000000000000000000000000000000010001000000000000000000000000000000000000000000000000000000010116011610080000000000000000000000000000000200000000000000100000000401001008000000000000001200000006000200280000001001101020011610681228208000000000000000000000000000000001000000000000000100000000000000000000000000000000000000000000000400000000000001000000000200101008000000000000000100000000000001160000000000000116000000000000100800000000000001160000000600101068000000020110122800120028102020800000000000000001000000010001010600000001000101140001010601141408000000000000010400000006000404080000010001101400011616081448480000000000000001160000000601001608000000040110144801040408140048000000000601101668000600080000000001100000100000001668000000000000

/-- Adjacency masks in the same executable port order; checked against seams. -/
def neighbors : Array Nat := #[3234,929,3160,852,2604,2579,1420,1347,714,314,2245,1077]

def joined (a b : Nat) : Bool := neighbors[a]!.testBit b

theorem graph_agreement : (List.range 12).all (fun a => (List.range 12).all (fun b =>
    joined a b == seams.contains (min a b,max a b))) = true := by decide +kernel

theorem graph_degree : (List.range 12).all (fun a =>
    ((List.range 12).filter (joined a)).length == 5) = true := by decide +kernel

def reference (p : Nat) : Int := if p == 1 || p == 3 || p == 10 then 0 else 1

def bit (b : Bool) : Int := if b then 1 else 0

/-- Twice the diagonal potential. On a degree-five screen, expanding the
squared differences gives 20 minus twice each like bond plus twice each
unlike bond. A private site contributes three for an excess and minus one
for a vacancy. -/
def potential (p m : Nat × Nat) : Int :=
  let like := bit (joined p.1 p.2) + bit (joined m.1 m.2)
  let unlike := bit (joined p.1 m.1) + bit (joined p.1 m.2) +
    bit (joined p.2 m.1) + bit (joined p.2 m.2)
  16*(20 - 2*like + 2*unlike + 3*(reference p.1 + reference p.2) -
    reference m.1 - reference m.2)

/-- For four sites connectivity means no isolated site and at least three
internal edges: a disconnected graph with no isolated site is two pairs. -/
def connectedFour (p m : Nat × Nat) : Bool :=
  let a := joined p.1 p.2
  let b := joined m.1 m.2
  let c := joined p.1 m.1
  let d := joined p.1 m.2
  let e := joined p.2 m.1
  let f := joined p.2 m.2
  (a || c || d) && (a || e || f) && (b || c || e) && (b || d || f) &&
    bit a + bit b + bit c + bit d + bit e + bit f ≥ 3

/-- Swapping opposite signs leaves the support fixed. Each swap into a blank
that leaves the connected subspace adds one to the doubled row sum. -/
def escaping (s : Nat) : Int :=
  seams.foldl (fun n e =>
    if s.testBit e.1 != s.testBit e.2 &&
        !onCluster.testBit (s ^^^ (2^e.1 + 2^e.2)) then n+1 else n) 0

def checkPair (p m : Nat × Nat) : Bool := Id.run do
  if mask p &&& mask m != 0 then return true
  let s := mask p ||| mask m
  let connected := onCluster.testBit s
  if connected != connectedFour p m then return false
  if !connected then return true
  let base := potential p m
  return base ≥ 310 || base + escaping s ≥ 310

def batch (a : Nat) : Bool :=
  (List.range a).all fun b => pairs.all (checkPair (a,b))

/-- All connected row sums are at least 310 in doubled units. -/
theorem connected_row_floor : (List.range 12).all batch = true := by
  change (batch 0 && (batch 1 && (batch 2 && (batch 3 && (batch 4 && (batch 5 &&
    (batch 6 && (batch 7 && (batch 8 && (batch 9 && (batch 10 && (batch 11 && true)))))))))))) = true
  simp only [Bool.and_eq_true]
  repeat' constructor <;> decide +kernel

def counts (a : Nat) : Nat × Nat :=
  (List.range a).foldl (fun total b =>
    pairs.foldl (fun t m =>
      if mask (a,b) &&& mask m != 0 then t else
        (t.1+1,t.2+(if onCluster.testBit (mask (a,b) ||| mask m) then 1 else 0))) total) (0,0)

set_option maxRecDepth 8192 in
theorem counted_words :
    ((List.range 12).map counts).foldl (fun a b => (a.1+b.1,a.2+b.2)) (0,0) =
      (2970,1260) := by decide +kernel

def differing (p m : Nat × Nat) : Nat :=
  let charge := fun i => (if (mask p).testBit i then 1 else if (mask m).testBit i then 2 else 0 : Nat)
  (seams.filter fun e => charge e.1 != charge e.2).length

/-- Excesses at 0,1 and vacancies at 2,4 are an unrestricted legal trial. -/
theorem unrestricted_trial :
    (1,0) ∈ pairs ∧ (4,2) ∈ pairs ∧ mask (1,0) &&& mask (4,2) = 0 ∧
    potential (1,0) (4,2) + differing (1,0) (4,2) = 288 := by decide +kernel

end OPH.NeutralClusterEnergy
