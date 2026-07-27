import Mathlib
import PortFrameGram

namespace OPH.CoreAxioms

open OPH.PortFrameGram (adj neighbors antipode)

/-! # Typed shadow of the three-axiom basis

Machine-checked structures shadowing the canonical basis in
`claims/axiom_registry.yaml` and `docs/AXIOM_REFERENCE.md`. No global Lean
axiom is introduced: every statement is a definition, a structure whose
fields are ordinary theorem arguments, or a `decide`-checked fact about the
concrete boundary template.

1. The exact A1 boundary packet: the oriented icosahedral template on
   `Fin 12` with its 30 edges and 20 coherently oriented faces, five-cycle
   links, and two-face edge incidence, each checked by `decide` against the
   same adjacency the Screen modules use.
2. Structure types for the A1 carrier tuple, the A2 interpretation
   naturality data, and the A3 information-projection specification, so
   downstream theorems can take the basis as typed hypotheses.

BOUNDARY. The packet facts are combinatorial statements about the template.
They do not assert semantic agreement, state selection, a response law, a
rotation-group symbol inside A1, or any physical identification. -/

/-! ## The oriented boundary template -/

/-- The twenty coherently oriented faces of the boundary template, over the
`PortFrameGram` adjacency. Each face lists its ports in positive cyclic
order; every undirected edge occurs in exactly two faces with opposite
induced directions. -/
def orientedFaces : List (Fin 12 × Fin 12 × Fin 12) :=
  [(0, 1, 2), (0, 3, 1), (0, 2, 4), (0, 6, 3), (0, 4, 6),
   (1, 5, 2), (1, 3, 7), (1, 7, 5), (2, 8, 4), (2, 5, 8),
   (3, 6, 9), (3, 9, 7), (4, 10, 6), (4, 8, 10), (5, 7, 11),
   (5, 11, 8), (6, 10, 9), (7, 9, 11), (8, 11, 10), (9, 10, 11)]

/-- The directed boundary edges of one oriented face. -/
def faceEdges (f : Fin 12 × Fin 12 × Fin 12) :
    List (Fin 12 × Fin 12) :=
  [(f.1, f.2.1), (f.2.1, f.2.2), (f.2.2, f.1)]

/-- All sixty directed edges induced by the oriented faces. -/
def directedEdges : List (Fin 12 × Fin 12) :=
  orientedFaces.flatMap faceEdges

/-- There are exactly twenty faces. -/
theorem face_count : orientedFaces.length = 20 := by decide

/-- There are exactly thirty undirected edges: every port pair is adjacent
exactly when it occurs as a directed face edge, and each undirected edge
carries exactly two directed occurrences. -/
theorem edge_count :
    (Finset.univ.filter
      (fun p : Fin 12 × Fin 12 => p.1 < p.2 ∧ adj p.1 p.2)).card = 30 := by
  decide

set_option maxRecDepth 8192 in
/-- Every face is an incidence triangle. -/
theorem faces_are_triangles :
    ∀ f ∈ orientedFaces,
      adj f.1 f.2.1 = true ∧ adj f.2.1 f.2.2 = true ∧
        adj f.2.2 f.1 = true := by decide +kernel

/-- Coherent orientation: every ordered adjacent pair occurs exactly once as
a directed face edge, so each undirected edge lies in exactly two faces with
opposite induced directions. -/
theorem edges_two_faces_opposite :
    ∀ i j : Fin 12, adj i j = true →
      (directedEdges.count (i, j) = 1 ∧ directedEdges.count (j, i) = 1) := by
  decide

/-- Every port has exactly five neighbors (restated from the Screen spine
for the packet interface). -/
theorem port_degree_five :
    ∀ i : Fin 12, (Finset.univ.filter (fun j => adj i j)).card = 5 := by
  decide

/-- Five-cycle links: within the link of each port, every neighbor is
adjacent to exactly two other neighbors. A two-regular graph on five
vertices is a single five-cycle, so each link is a pentagon. -/
theorem link_two_regular :
    ∀ i : Fin 12, ∀ j ∈ neighbors i,
      ((neighbors i).filter (fun k => adj j k)).length = 2 := by decide

/-- Adjacency closure of the first link vertex inside the link of `i`,
iterated five times. -/
def linkClosure (i : Fin 12) : List (Fin 12) :=
  (List.range 5).foldl
    (fun acc _ =>
      (acc ++ (neighbors i).filter (fun v => acc.any (fun a => adj a v))).dedup)
    ((neighbors i).take 1)

/-- Each link is connected: every neighbor of `i` lies in the adjacency
closure of the first link vertex, so the two-regular link graph is one
five-cycle rather than a disjoint union. -/
theorem link_connected :
    ∀ i : Fin 12, ∀ j ∈ neighbors i, j ∈ linkClosure i := by decide

/-! ## Typed basis structures -/

/-- A1 boundary packet over an abstract port type: the incidence and
oriented-face data of one local carrier, required to be isomorphic to the
concrete template by an orientation-preserving relabeling. -/
structure BoundaryPacket (V : Type*) [Fintype V] [DecidableEq V] where
  adjacent : V → V → Bool
  oriented : List (V × V × V)
  relabel : V ≃ Fin 12
  adjacent_iff : ∀ i j, adjacent i j = adj (relabel i) (relabel j)
  oriented_matches :
    List.Perm
      (oriented.map (fun f => (relabel f.1, relabel f.2.1, relabel f.2.2)))
      orientedFaces

/-- A1 carrier tuple shadow: the quotient-visible operational fields of one
local carrier. Algebraic components are abstract types with the maps the
axiom requires; the boundary component is the exact packet. -/
structure CarrierShadow (V State Record Interface : Type*)
    [Fintype V] [DecidableEq V] where
  boundary : BoundaryPacket V
  readback : V → State → Interface
  record : State → Record
  repair : State → State
  checkpoint : State → Record

/-- A2 shadow: an interpretation map together with the naturality squares it
must close, given as explicit families of data-transport and
meaning-transport maps. -/
structure MeaningNaturality (Data Meaning Index : Type*) where
  interpret : Index → Data → Meaning
  dataTransport : Index → Index → Data → Data
  meaningTransport : Index → Index → Meaning → Meaning
  natural :
    ∀ r s d, interpret r (dataTransport s r d)
      = meaningTransport s r (interpret s d)

/-- A3 shadow: an information-projection specification over a finite index
of local state coordinates, with a reference family, positive weights, a
feasible predicate, and the projection characterized as the unique feasible
minimizer of the weighted divergence. -/
structure InfoProjectionSpec (StateFam : Type*) where
  feasible : StateFam → Prop
  divergence : StateFam → ℝ
  realized : StateFam
  realized_feasible : feasible realized
  realized_minimal :
    ∀ s, feasible s → divergence realized ≤ divergence s
  realized_unique :
    ∀ s, feasible s → divergence s = divergence realized → s = realized

#print axioms face_count
#print axioms edge_count
#print axioms faces_are_triangles
#print axioms edges_two_faces_opposite
#print axioms port_degree_five
#print axioms link_two_regular
#print axioms link_connected

end OPH.CoreAxioms
