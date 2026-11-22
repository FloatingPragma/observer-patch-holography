import ObserverPatchHolography.ObservablePriority.CellBits

namespace PriorityCellData
open ObserverPatchSharp CellBits
structure Row where
  cell : Nat
  slots : Three Nat
  remote : Three Nat
  spare : Three Nat

def originalChart : List (Nat × Nat) := [(0,1),(0,5),(0,7),(0,10),(0,11),(1,5),(1,7),(1,8),(1,9),(2,3),(2,4),(2,6),(2,10),(2,11),(3,4),(3,6),(3,8),(3,9),(4,5),(4,9),(4,11),(5,9),(5,11),(6,7),(6,8),(6,10),(7,8),(7,10),(8,9),(10,11)]

def ed (a b : Nat) : Edge := if h : a < b then ⟨a,b,h⟩ else ⟨b,max (b+1) a, by omega⟩
def cellRank (e : Edge) := if e.lo / 12 = e.hi / 12 then
  1 + originalChart.idxOf (e.lo % 12,e.hi % 12) else 0

def make (n : Nat) (row : Row) : CellChart :=
  let o := row.cell*12
  let v : Three Nat := ⟨o+row.slots.a,o+row.slots.b,o+row.slots.c⟩
  ⟨n,(List.range 12).map (o+·),
   originalChart.reverse.map (fun (a,b) => ed (o+a) (o+b)),
   ⟨ed v.a row.remote.a,ed v.b row.remote.b,ed v.c row.remote.c⟩,
   v,row.remote,row.spare⟩

def legalB (rank : Edge → Nat) : List Edge → Bool
  | [] => true
  | a::xs => xs.all (fun b => decide (a≠b) && (!meets a b || decide (rank b<rank a))) && legalB rank xs

theorem legalB_spec (rank : Edge → Nat) (xs : List Edge) :
    legalB rank xs = true ↔ Ordered rank xs := by
  induction xs with
  | nil => simp [legalB,Ordered]
  | cons a xs ih =>
    simp only [legalB,Bool.and_eq_true_iff,ih,Ordered,List.pairwise_cons,List.all_eq_true]
    simp
    grind

def admitted (n : Nat) (r : Row) : Bool :=
  cellCheck (make n r) && legalB cellRank (menu (make n r))

theorem row_exact (n : Nat) (r : Row) (cert : admitted n r = true) :
    let c := make n r
    Communicates c (fun L => localWord c (sharpRun cellRank (menu c) L)) 3 ∧
    ∀ b, Communicates c (fun L => localWord c (sharpRun cellRank (menu c) L)) b → 3 ≤ b :=
  exactly_three _ _ ((legalB_spec _ _).mp (Bool.and_eq_true_iff.mp cert).2) (Bool.and_eq_true_iff.mp cert).1

end PriorityCellData
