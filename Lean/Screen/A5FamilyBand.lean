import Mathlib
import A5PortAction
import PortFrameGram

namespace OPH.A5FamilyBand

open OPH.A5PortAction (perms)

/-! # The canonical rank-three family band, exactly (issue #569)

Integer certificate for the spectral selection theorem behind the
family-band attachment.  Over `Z[sqrt5]`, encoded as integer pairs
`(x, y) = x + y*sqrt5`, the twelve-port adjacency `A` of
`PortFrameGram.neighbors` has the exact resolution

  `120*P1 = M1`, `40*P3 = X + sqrt5*Y`, `40*P3' = X - sqrt5*Y`,
  `24*P5 = M5`,

with `M1`, `X`, `Y`, `M5` the integer tables below.  Kernel `decide`
verifies, with no floating point and no axioms beyond the standard three:

* the tables bind to the listed incidence (`table_matches_incidence`);
* each scaled projector is idempotent, self-adjoint, an exact
  eigenband (`A*M1 = 5*M1`, `A*(X + sqrt5 Y) = sqrt5*(X + sqrt5 Y)`
  componentwise, `A*M5 = -M5`), and the four resolve the identity
  (`resolution_complete`), so the multiplicities are exactly
  `1 + 3 + 3' + 5`;
* every listed rotation fixes each band (`band_equivariant`), the
  trivial band absorbs all sixty rotations while the `3`, `3'`, and `5`
  band actions have trivial kernel (`band_kernels`);
* the exact seam-cost order is strict,
  `5 - sqrt5 < 6 < 5 + sqrt5` (`cost_order_strict`), so among the
  faithful in-window bands the Laplacian cost has a unique strict
  minimizer: the rank-three band (`family_band_selected`).

BOUNDARY.  This file is coefficient algebra on the listed carrier.  The
two interface clauses that make the selection physical (realization of
the multiplicity object inside the screen coefficient space, and
comparison by the operational cost order) are premises carried by the
executable certificate `family_band_attachment_certificate.py` and
issue #569; they are not proved here, and the #617 copy-count
invisibility theorem for external completions holds unchanged. -/

/-- Integer entry lookup in a listed table. -/
def entry (m : List (List Int)) (i j : Nat) : Int :=
  (m.getD i []).getD j 0

/-- Twelve-port matrix product entry. -/
def mulEntry (m n : List (List Int)) (i j : Nat) : Int :=
  (List.range 12).foldl (fun s k => s + entry m i k * entry n k j) 0

/-- Entrywise equality of two matrix expressions on the twelve ports. -/
def matEq (f g : Nat → Nat → Int) : Bool :=
  (List.range 12).all fun i => (List.range 12).all fun j => f i j == g i j

/-- The listed adjacency table of the twelve-port screen. -/
def At : List (List Int) := [
    [0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0],
    [1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0],
    [1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0],
    [1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0],
    [1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0],
    [0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1],
    [1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0],
    [0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1],
    [0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1],
    [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1],
    [0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1],
    [0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0]]

/-- `120 * P1`: the trivial band. -/
def M1t : List (List Int) :=
  List.replicate 12 (List.replicate 12 10)

/-- Real part of `40 * P3`: ten times identity minus antipode. -/
def Xt : List (List Int) := [
    [10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -10],
    [0, 10, 0, 0, 0, 0, 0, 0, 0, 0, -10, 0],
    [0, 0, 10, 0, 0, 0, 0, 0, 0, -10, 0, 0],
    [0, 0, 0, 10, 0, 0, 0, 0, -10, 0, 0, 0],
    [0, 0, 0, 0, 10, 0, 0, -10, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 10, -10, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, -10, 10, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, -10, 0, 0, 10, 0, 0, 0, 0],
    [0, 0, 0, -10, 0, 0, 0, 0, 10, 0, 0, 0],
    [0, 0, -10, 0, 0, 0, 0, 0, 0, 10, 0, 0],
    [0, -10, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0],
    [-10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10]]

/-- `sqrt5` part of `40 * P3`: twice adjacency minus distance-two. -/
def Yt : List (List Int) := [
    [0, 2, 2, 2, 2, -2, 2, -2, -2, -2, -2, 0],
    [2, 0, 2, 2, -2, 2, -2, 2, -2, -2, 0, -2],
    [2, 2, 0, -2, 2, 2, -2, -2, 2, 0, -2, -2],
    [2, 2, -2, 0, -2, -2, 2, 2, 0, 2, -2, -2],
    [2, -2, 2, -2, 0, -2, 2, 0, 2, -2, 2, -2],
    [-2, 2, 2, -2, -2, 0, 0, 2, 2, -2, -2, 2],
    [2, -2, -2, 2, 2, 0, 0, -2, -2, 2, 2, -2],
    [-2, 2, -2, 2, 0, 2, -2, 0, -2, 2, -2, 2],
    [-2, -2, 2, 0, 2, 2, -2, -2, 0, -2, 2, 2],
    [-2, -2, 0, 2, -2, -2, 2, 2, -2, 0, 2, 2],
    [-2, 0, -2, -2, 2, -2, 2, -2, 2, 2, 0, 2],
    [0, -2, -2, -2, -2, 2, -2, 2, 2, 2, 2, 0]]

/-- `24 * P5`: the five band. -/
def M5t : List (List Int) := [
    [10, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, 10],
    [-2, 10, -2, -2, -2, -2, -2, -2, -2, -2, 10, -2],
    [-2, -2, 10, -2, -2, -2, -2, -2, -2, 10, -2, -2],
    [-2, -2, -2, 10, -2, -2, -2, -2, 10, -2, -2, -2],
    [-2, -2, -2, -2, 10, -2, -2, 10, -2, -2, -2, -2],
    [-2, -2, -2, -2, -2, 10, 10, -2, -2, -2, -2, -2],
    [-2, -2, -2, -2, -2, 10, 10, -2, -2, -2, -2, -2],
    [-2, -2, -2, -2, 10, -2, -2, 10, -2, -2, -2, -2],
    [-2, -2, -2, 10, -2, -2, -2, -2, 10, -2, -2, -2],
    [-2, -2, 10, -2, -2, -2, -2, -2, -2, 10, -2, -2],
    [-2, 10, -2, -2, -2, -2, -2, -2, -2, -2, 10, -2],
    [10, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, 10]]

set_option maxHeartbeats 4000000 in
/-- The listed adjacency table equals the corpus incidence
`PortFrameGram.adj` entry for entry. -/
theorem table_matches_incidence :
    ∀ i j : Fin 12,
      entry At i j = (if OPH.PortFrameGram.adj i j then 1 else 0) := by
  decide

set_option maxHeartbeats 4000000 in
/-- Scaled idempotence of the `3` band over `Z[sqrt5]`:
`(X + sqrt5 Y)^2 = 40 (X + sqrt5 Y)` componentwise. -/
theorem band3_idempotent :
    matEq (fun i j => mulEntry Xt Xt i j + 5 * mulEntry Yt Yt i j)
        (fun i j => 40 * entry Xt i j) = true ∧
      matEq (fun i j => mulEntry Xt Yt i j + mulEntry Yt Xt i j)
        (fun i j => 40 * entry Yt i j) = true := by
  constructor <;> decide

set_option maxHeartbeats 4000000 in
/-- The `3` band is the exact `sqrt5` eigenband of the adjacency:
`A (X + sqrt5 Y) = sqrt5 (X + sqrt5 Y)` componentwise. -/
theorem band3_eigen :
    matEq (mulEntry At Xt) (fun i j => 5 * entry Yt i j) = true ∧
      matEq (mulEntry At Yt) (entry Xt) = true := by
  constructor <;> decide

set_option maxHeartbeats 4000000 in
/-- Trivial and five bands: scaled idempotence and exact eigenvalues
`5` and `-1`. -/
theorem outer_bands :
    matEq (mulEntry M1t M1t) (fun i j => 120 * entry M1t i j) = true ∧
      matEq (mulEntry At M1t) (fun i j => 5 * entry M1t i j) = true ∧
      matEq (mulEntry M5t M5t) (fun i j => 24 * entry M5t i j) = true ∧
      matEq (mulEntry At M5t) (fun i j => - entry M5t i j) = true := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> decide

set_option maxHeartbeats 4000000 in
/-- Self-adjointness of every table. -/
theorem bands_symmetric :
    matEq (entry Xt) (fun i j => entry Xt j i) = true ∧
      matEq (entry Yt) (fun i j => entry Yt j i) = true ∧
      matEq (entry M5t) (fun i j => entry M5t j i) = true := by
  refine ⟨?_, ?_, ?_⟩ <;> decide

set_option maxHeartbeats 4000000 in
/-- The four bands resolve the identity:
`M1 + 3*(40 P3) + 3*(40 P3') + 5*M5 = 120 * I` reduces to
`M1 + 6X + 5*M5 = 120 I` because the `sqrt5` parts cancel in the
Galois pair.  With the traces below this pins the multiplicities
`1 + 3 + 3' + 5`. -/
theorem resolution_complete :
    matEq (fun i j => entry M1t i j + 6 * entry Xt i j + 5 * entry M5t i j)
      (fun i j => if i == j then 120 else 0) = true := by
  decide

set_option maxHeartbeats 4000000 in
/-- Exact traces: `tr X = 120 = 40 * 3`, `tr Y = 0`, `tr M1 = 120 * 1`,
`tr M5 = 24 * 5`, so the band ranks are `3`, `1`, and `5`. -/
theorem band_traces :
    (List.range 12).foldl (fun s i => s + entry Xt i i) 0 = 120 ∧
      (List.range 12).foldl (fun s i => s + entry Yt i i) 0 = 0 ∧
      (List.range 12).foldl (fun s i => s + entry M1t i i) 0 = 120 ∧
      (List.range 12).foldl (fun s i => s + entry M5t i i) 0 = 120 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> decide

/-- Row-and-column permutation invariance of a table under a listed
rotation. -/
def fixedBy (m : List (List Int)) (p : List Nat) : Bool :=
  (List.range 12).all fun i =>
    (List.range 12).all fun j =>
      entry m (p.getD i 0) (p.getD j 0) == entry m i j

set_option maxHeartbeats 16000000 in
set_option maxRecDepth 16384 in
/-- Every listed rotation fixes every band. -/
theorem band_equivariant :
    perms.all (fun p =>
      fixedBy Xt p && fixedBy Yt p && fixedBy M5t p) = true := by
  decide

/-- A rotation acts as the identity on a band exactly when permuting the
rows of its symmetric projector table changes nothing. -/
def actsTrivially (m : List (List Int)) (p : List Nat) : Bool :=
  (List.range 12).all fun i =>
    (List.range 12).all fun j =>
      entry m (p.getD i 0) j == entry m i j

set_option maxHeartbeats 16000000 in
set_option maxRecDepth 16384 in
/-- Band action kernels: the trivial band absorbs all sixty rotations;
the `3` and `5` band actions are faithful (kernel order one).  The `3'`
band shares the kernel of the `3` band because its tables are `X` and
`-Y`. -/
theorem band_kernels :
    ((perms.filter (actsTrivially M1t)).length = 60) ∧
      ((perms.filter (fun p => actsTrivially Xt p && actsTrivially Yt p)).length = 1) ∧
      ((perms.filter (actsTrivially M5t)).length = 1) := by
  refine ⟨?_, ?_, ?_⟩ <;> decide

/-- Strict positivity of `x + y*sqrt5` on integer pairs. -/
def z5pos (x y : Int) : Bool :=
  if y == 0 then decide (0 < x)
  else if x == 0 then decide (0 < y)
  else if 0 < x ∧ 0 < y then true
  else if x < 0 ∧ y < 0 then false
  else if 0 < x then decide (5 * y * y < x * x)
  else decide (x * x < 5 * y * y)

/-- Strict order on `Z[sqrt5]` pairs. -/
def z5lt (a : Int × Int) (b : Int × Int) : Bool :=
  z5pos (b.1 - a.1) (b.2 - a.2)

/-- The exact seam-cost order of the three admissible bands is strict:
`5 - sqrt5 < 6 < 5 + sqrt5`. -/
theorem cost_order_strict :
    z5lt (5, -1) (6, 0) = true ∧ z5lt (6, 0) (5, 1) = true := by
  constructor <;> decide

/-- The candidate table: label, complex dimension, seam cost as a
`Z[sqrt5]` pair, and kernel order of the band action from
`band_kernels`. -/
def candidates : List (String × Nat × (Int × Int) × Nat) :=
  [("1", 1, (0, 0), 60), ("3", 3, (5, -1), 1),
   ("3p", 3, (5, 1), 1), ("5", 5, (6, 0), 1)]

/-- Admissibility: faithful band action and dimension inside the pinned
physical window `[3, 5]`. -/
def admissible (c : String × Nat × (Int × Int) × Nat) : Bool :=
  c.2.2.2 == 1 && decide (3 ≤ c.2.1) && decide (c.2.1 ≤ 5)

/-- **The family-band selection.**  Among the admissible candidates the
seam cost has a unique strict minimizer, and it is the rank-three band:
every other admissible candidate costs strictly more. -/
theorem family_band_selected :
    (candidates.filter admissible).length = 3 ∧
      ((candidates.filter admissible).all fun c =>
        c.1 == "3" || z5lt (5, -1) c.2.2.1) = true := by
  constructor <;> decide

/-- Without the faithfulness clause the cost minimizer over all four
bands is the trivial band at cost zero: the clause is load-bearing and
cost minimization alone selects nothing physical. -/
theorem trivial_band_without_faithfulness :
    (candidates.all fun c => c.1 == "1" || z5lt (0, 0) c.2.2.1) = true := by
  decide

#print axioms table_matches_incidence
#print axioms band3_idempotent
#print axioms band3_eigen
#print axioms outer_bands
#print axioms resolution_complete
#print axioms band_traces
#print axioms band_equivariant
#print axioms band_kernels
#print axioms cost_order_strict
#print axioms family_band_selected
#print axioms trivial_band_without_faithfulness

end OPH.A5FamilyBand
