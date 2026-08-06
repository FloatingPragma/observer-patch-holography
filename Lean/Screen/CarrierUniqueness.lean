import Mathlib
import PortFrameGram

namespace OPH.CarrierUniqueness

open scoped BigOperators
open OPH.PortFrameGram (adj)

set_option maxRecDepth 8000
set_option maxHeartbeats 4000000

/-! # Carrier counting and characterization for the twelve-port boundary packet

Issue #676 (A2, optional).  The A1 boundary packet `K = (P, E, F, o)` declares
twelve ports, thirty edges, twenty faces, degree-five ports, five-cycle links,
and coherently oriented two-face edges, combinatorially the icosahedral
boundary complex.  This file proves the counting core of that declaration
together with an exact characterization of the port count by the Euler
characteristic, a partial rigidity theorem, a fully certified icosahedral
witness, and three exact controls.

Main statements.

* `BoundaryComplex`: a finite triangulated boundary complex over a vertex
  type: two-element edges inside the port set, three-element faces whose
  two-element subsets are edges, and every edge contained in exactly two
  faces.
* `handshake_five` and `triangulation_count`: the double-counting identities
  `2E = 5P` (degree five) and `3F = 2E` (triangle faces, two per edge).  Both
  are proved abstractly by counting one incidence set two ways; neither
  planarity machinery nor an assumed Euler relation enters.
* `six_mul_euler_eq_ports`: for every degree-five complex,
  `6 (P - E + F) = P`.  The Euler characteristic is derived rather than
  assumed, and it determines the port count exactly: `P = 12` if and only if
  `P - E + F = 2`, the spherical value (`twelve_ports_iff_euler_two`).
* `counts_forced`: twelve ports force `E = 30` and `F = 20`
  (handshake plus Euler counting, with the Euler relation as output).
* `Icosa.carrier` and `Icosa.criteria`: the icosahedral witness built on the
  committed `PortFrameGram` incidence, with kernel-checked port, edge, and
  face counts, degree five, link five-cycles listed explicitly, a coherent
  orientation on the twenty faces (`Icosa.coherent_orientation`), and edge
  connectivity.
* `faces_determined`: any boundary complex on the icosahedral edge set
  carries exactly the icosahedral face set.  The edge set forces the face
  set inside this axiom class.
* Controls.  `Octa.carrier` satisfies every complex law with degree four and
  six ports, so the degree criterion is load-bearing.  `Cube.edges` admits
  no boundary complex at all (`Cube.no_boundary_complex`), so the
  triangulation criterion is load-bearing.  `Hemi.carrier`, the
  hemi-icosahedron, satisfies every law including degree five with six ports
  and Euler characteristic one, so the twelve-port criterion is load-bearing
  and is exactly the choice of the spherical Euler value.

BOUNDARY.  Everything in this file is combinatorial characterization under
the declared finite criteria.  No physical selection claim is made, and A1
stands as an axiom independently of this file; issue #676 is optional
explanatory mathematics.  Full rigidity of the combinatorial type (any two
complexes satisfying the criteria are isomorphic) is out of scope: it needs
the rooted-flag classification argument, and the labelled edge-set search
space (66 choose 30) is beyond kernel enumeration.  The delivered level is
the forced counting with derived Euler relation, the Euler characterization
of the port count, edge-set-forces-face-set rigidity, the certified witness,
and the three controls. -/

variable {V : Type} [DecidableEq V]

/-- A finite triangulated boundary complex: ports, two-element edges inside
the port set, three-element faces whose two-element subsets are edges, and
every edge contained in exactly two faces.  Orientation data is handled
separately on the witness (`Icosa.orientedFaces`); the counting theorems do
not need it. -/
structure BoundaryComplex (V : Type) [DecidableEq V] where
  /-- The port (vertex) set. -/
  ports : Finset V
  /-- The edge set, as two-element port sets. -/
  edges : Finset (Finset V)
  /-- The face set, as three-element port sets. -/
  faces : Finset (Finset V)
  /-- Every edge has exactly two ports. -/
  edge_card : ∀ e ∈ edges, e.card = 2
  /-- Every edge lies inside the port set. -/
  edge_sub_ports : ∀ e ∈ edges, e ⊆ ports
  /-- Every face has exactly three ports. -/
  face_card : ∀ f ∈ faces, f.card = 3
  /-- Every two-element subset of a face is an edge. -/
  face_edges : ∀ f ∈ faces, ∀ e ∈ Finset.powersetCard 2 f, e ∈ edges
  /-- Every edge lies in exactly two faces. -/
  edge_two_faces : ∀ e ∈ edges, (faces.filter (fun f => e ⊆ f)).card = 2

/-! ## Double counting -/

/-- Counting an incidence set two ways: summing fiber sizes over either side
of a decidable relation gives the same total. -/
theorem count_swap {α β : Type} (s : Finset α) (t : Finset β)
    (r : α → β → Prop) [∀ a b, Decidable (r a b)] :
    ∑ a ∈ s, (t.filter (fun b => r a b)).card
      = ∑ b ∈ t, (s.filter (fun a => r a b)).card := by
  simp only [Finset.card_filter]
  exact Finset.sum_comm

/-- Handshake: the degrees over all ports sum to twice the edge count. -/
theorem handshake (K : BoundaryComplex V) :
    ∑ p ∈ K.ports, (K.edges.filter (fun e => p ∈ e)).card = 2 * K.edges.card := by
  refine (count_swap K.ports K.edges (fun p e => p ∈ e)).trans ?_
  have h : ∀ e ∈ K.edges, (K.ports.filter (fun p => p ∈ e)).card = 2 := by
    intro e he
    have hsub : e ⊆ K.ports := K.edge_sub_ports e he
    have hfe : K.ports.filter (fun p => p ∈ e) = e := by
      ext x
      simp only [Finset.mem_filter]
      exact ⟨fun hx => hx.2, fun hx => ⟨hsub hx, hx⟩⟩
    rw [hfe]
    exact K.edge_card e he
  rw [Finset.sum_congr rfl h, Finset.sum_const, smul_eq_mul, mul_comm]

/-- Degree-five handshake: if every port has degree five then `2E = 5P`. -/
theorem handshake_five (K : BoundaryComplex V)
    (hdeg : ∀ p ∈ K.ports, (K.edges.filter (fun e => p ∈ e)).card = 5) :
    2 * K.edges.card = 5 * K.ports.card := by
  have h := handshake K
  rw [Finset.sum_congr rfl hdeg, Finset.sum_const, smul_eq_mul] at h
  omega

/-- Triangulation count: `3F = 2E`.  Each face carries exactly three edges
(its two-element subsets), each edge lies in exactly two faces, and the
face-edge incidence set is counted two ways.  No degree or port-count
hypothesis enters. -/
theorem triangulation_count (K : BoundaryComplex V) :
    3 * K.faces.card = 2 * K.edges.card := by
  have hswap := count_swap K.faces K.edges (fun f e => e ⊆ f)
  have hface : ∀ f ∈ K.faces, (K.edges.filter (fun e => e ⊆ f)).card = 3 := by
    intro f hf
    have hpc : K.edges.filter (fun e => e ⊆ f) = Finset.powersetCard 2 f := by
      ext e
      simp only [Finset.mem_filter, Finset.mem_powersetCard]
      constructor
      · rintro ⟨he, hef⟩
        exact ⟨hef, K.edge_card e he⟩
      · rintro ⟨hef, h2⟩
        exact ⟨K.face_edges f hf e (Finset.mem_powersetCard.mpr ⟨hef, h2⟩), hef⟩
    rw [hpc, Finset.card_powersetCard, K.face_card f hf]
    decide
  rw [Finset.sum_congr rfl hface, Finset.sum_const, smul_eq_mul,
    Finset.sum_congr rfl K.edge_two_faces, Finset.sum_const, smul_eq_mul] at hswap
  omega

/-! ## The Euler characteristic is derived and pins the port count -/

/-- For every degree-five boundary complex, six times the Euler
characteristic equals the port count: `6 (P - E + F) = P`.  The Euler
relation is an output of the two counting identities, never a hypothesis. -/
theorem six_mul_euler_eq_ports (K : BoundaryComplex V)
    (hdeg : ∀ p ∈ K.ports, (K.edges.filter (fun e => p ∈ e)).card = 5) :
    6 * ((K.ports.card : ℤ) - K.edges.card + K.faces.card) = K.ports.card := by
  have h1 := handshake_five K hdeg
  have h2 := triangulation_count K
  omega

/-- Twelve ports and the spherical Euler value determine each other in the
degree-five class: `P = 12` if and only if `P - E + F = 2`. -/
theorem twelve_ports_iff_euler_two (K : BoundaryComplex V)
    (hdeg : ∀ p ∈ K.ports, (K.edges.filter (fun e => p ∈ e)).card = 5) :
    K.ports.card = 12 ↔ (K.ports.card : ℤ) - K.edges.card + K.faces.card = 2 := by
  have h := six_mul_euler_eq_ports K hdeg
  omega

/-- Twelve degree-five ports force thirty edges and twenty faces. -/
theorem counts_forced (K : BoundaryComplex V)
    (h12 : K.ports.card = 12)
    (hdeg : ∀ p ∈ K.ports, (K.edges.filter (fun e => p ∈ e)).card = 5) :
    K.edges.card = 30 ∧ K.faces.card = 20 := by
  have h1 := handshake_five K hdeg
  have h2 := triangulation_count K
  omega

/-- The converse direction packaged: a degree-five complex with the spherical
Euler value has exactly the packet counts `P = 12`, `E = 30`, `F = 20`. -/
theorem spherical_counts (K : BoundaryComplex V)
    (hdeg : ∀ p ∈ K.ports, (K.edges.filter (fun e => p ∈ e)).card = 5)
    (hsph : (K.ports.card : ℤ) - K.edges.card + K.faces.card = 2) :
    K.ports.card = 12 ∧ K.edges.card = 30 ∧ K.faces.card = 20 := by
  have h1 := handshake_five K hdeg
  have h2 := triangulation_count K
  omega

/-! ## The declared criteria as a structure interface -/

/-- The declared finite criteria of the A1 boundary packet, as a structure
interface: a triangulated boundary complex with twelve ports, degree-five
ports, and edge connectivity.  Five-cycle links are a consequence rather
than an input: `edge_two_faces` makes the link of each port a two-regular
graph on its five neighbors, and a two-regular simple graph on five vertices
is a single five-cycle (five admits no decomposition into shorter simple
cycles).  The witness `Icosa.criteria` lists its link five-cycles
explicitly. -/
structure EchosahedralCriteria (V : Type) [DecidableEq V]
    extends BoundaryComplex V where
  /-- Twelve ports. -/
  ports_card : ports.card = 12
  /-- Every port has degree five. -/
  degree_five : ∀ p ∈ ports, (edges.filter (fun e => p ∈ e)).card = 5
  /-- Any two ports are joined through edges. -/
  connected : ∀ p ∈ ports, ∀ q ∈ ports,
    Relation.ReflTransGen (fun a b => ({a, b} : Finset V) ∈ edges) p q

/-- Under the declared criteria the edge count is thirty. -/
theorem EchosahedralCriteria.edges_thirty (C : EchosahedralCriteria V) :
    C.edges.card = 30 :=
  (counts_forced C.toBoundaryComplex C.ports_card C.degree_five).1

/-- Under the declared criteria the face count is twenty. -/
theorem EchosahedralCriteria.faces_twenty (C : EchosahedralCriteria V) :
    C.faces.card = 20 :=
  (counts_forced C.toBoundaryComplex C.ports_card C.degree_five).2

/-- Under the declared criteria the Euler relation `P - E + F = 2` holds. -/
theorem EchosahedralCriteria.euler_two (C : EchosahedralCriteria V) :
    (C.ports.card : ℤ) - C.edges.card + C.faces.card = 2 :=
  (twelve_ports_iff_euler_two C.toBoundaryComplex C.degree_five).mp C.ports_card

/-! ## Building complexes from Boolean adjacency data -/

/-- The edge set generated by a Boolean adjacency on `Fin n`: the
two-element subsets whose distinct members are adjacent. -/
def edgesOf {n : ℕ} (adjB : Fin n → Fin n → Bool) : Finset (Finset (Fin n)) :=
  (Finset.powersetCard 2 Finset.univ).filter
    (fun e => ∀ i ∈ e, ∀ j ∈ e, i ≠ j → adjB i j = true)

/-- A boundary complex from Boolean adjacency data and a declared face set:
faces must be three-element cliques and every generated edge must lie in
exactly two faces. -/
def ofData {n : ℕ} (adjB : Fin n → Fin n → Bool) (F : Finset (Finset (Fin n)))
    (hface3 : ∀ f ∈ F, f.card = 3)
    (hclique : ∀ f ∈ F, ∀ i ∈ f, ∀ j ∈ f, i ≠ j → adjB i j = true)
    (htwo : ∀ e ∈ edgesOf adjB, (F.filter (fun f => e ⊆ f)).card = 2) :
    BoundaryComplex (Fin n) where
  ports := Finset.univ
  edges := edgesOf adjB
  faces := F
  edge_card := fun e he => (Finset.mem_powersetCard.mp (Finset.mem_filter.mp he).1).2
  edge_sub_ports := fun e _ => Finset.subset_univ e
  face_card := hface3
  face_edges := by
    intro f hf e he
    obtain ⟨hef, h2⟩ := Finset.mem_powersetCard.mp he
    refine Finset.mem_filter.mpr
      ⟨Finset.mem_powersetCard.mpr ⟨Finset.subset_univ e, h2⟩, ?_⟩
    intro i hi j hj hij
    exact hclique f hf i (hef hi) j (hef hj) hij
  edge_two_faces := htwo

/-- A two-element subset of a set from two distinct members. -/
theorem pair_mem_powersetCard {f : Finset V} {x y : V}
    (hx : x ∈ f) (hy : y ∈ f) (hxy : x ≠ y) :
    ({x, y} : Finset V) ∈ Finset.powersetCard 2 f :=
  Finset.mem_powersetCard.mpr
    ⟨Finset.insert_subset_iff.mpr ⟨hx, Finset.singleton_subset_iff.mpr hy⟩,
      Finset.card_pair hxy⟩

/-! ### Kernel-evaluation bridges

Kernel evaluation of bounded quantifiers over a `Finset` of `Finset`s is
prohibitively expensive at this scale, while Boolean list folds evaluate
directly.  The three bridges below let every witness obligation be a
kernel-checked `List.all` computation while the theorem statements stay at
the `Finset` level. -/

/-- Bounded quantification over a list-backed `Finset` from a Boolean fold. -/
theorem ball_of_all {α : Type} {l : List α} {h : Multiset.Nodup (l : Multiset α)}
    {p : α → Prop} [DecidablePred p]
    (hall : l.all (fun a => decide (p a)) = true) :
    ∀ a ∈ (⟨(l : Multiset α), h⟩ : Finset α), p a := by
  intro a ha
  exact of_decide_eq_true (List.all_eq_true.mp hall a (by simpa using ha))

/-- Bounded quantification over a list from a Boolean fold. -/
theorem ball_list_of_all {α : Type} {l : List α} {p : α → Prop} [DecidablePred p]
    (hall : l.all (fun a => decide (p a)) = true) :
    ∀ a ∈ l, p a :=
  fun a ha => of_decide_eq_true (List.all_eq_true.mp hall a ha)

/-- Filter cardinality on a list-backed `Finset` as a list count. -/
theorem filter_card_coe {α : Type} [DecidableEq α] {l : List α}
    {h : Multiset.Nodup (l : Multiset α)} (p : α → Prop) [DecidablePred p] :
    ((⟨(l : Multiset α), h⟩ : Finset α).filter p).card
      = l.countP (fun a => decide (p a)) := by
  simp [Finset.filter, Finset.card, Multiset.filter_coe, List.countP_eq_length_filter]

/-! ## The icosahedral witness

Built on the committed `PortFrameGram` incidence (`neighbors`, `adj`), the
vetted twelve-port encoding shared by the `A5` files. -/

namespace Icosa

/-- The thirty edges of the icosahedral boundary, as a duplicate-free list
of port pairs. -/
def edgesList : List (Finset (Fin 12)) :=
  [{0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 6}, {1, 2},
   {1, 3}, {1, 5}, {1, 7}, {2, 4}, {2, 5}, {2, 8},
   {3, 6}, {3, 7}, {3, 9}, {4, 6}, {4, 8}, {4, 10},
   {5, 7}, {5, 8}, {5, 11}, {6, 9}, {6, 10}, {7, 9},
   {7, 11}, {8, 10}, {8, 11}, {9, 10}, {9, 11}, {10, 11}]

/-- The twenty faces of the icosahedral boundary, as a duplicate-free list
of port triples. -/
def facesList : List (Finset (Fin 12)) :=
  [{0, 1, 2}, {0, 1, 3}, {0, 2, 4}, {0, 3, 6}, {0, 4, 6},
   {1, 2, 5}, {1, 3, 7}, {1, 5, 7}, {2, 4, 8}, {2, 5, 8},
   {3, 6, 9}, {3, 7, 9}, {4, 6, 10}, {4, 8, 10}, {5, 7, 11},
   {5, 8, 11}, {6, 9, 10}, {7, 9, 11}, {8, 10, 11}, {9, 10, 11}]

/-- The edge set. -/
def edges : Finset (Finset (Fin 12)) :=
  ⟨(edgesList : Multiset (Finset (Fin 12))), by decide⟩

/-- The face set. -/
def faces : Finset (Finset (Fin 12)) :=
  ⟨(facesList : Multiset (Finset (Fin 12))), by decide⟩

/-- The declared edge set is exactly the edge set generated by the committed
`PortFrameGram` incidence. -/
theorem edges_eq_edgesOf : edges = edgesOf adj := by decide

/-- The declared face set is exactly the set of three-element cliques of the
committed incidence: faces are determined by adjacency alone. -/
theorem faces_eq_triangles :
    faces = (Finset.powersetCard 3 Finset.univ).filter
      (fun f => ∀ i ∈ f, ∀ j ∈ f, i ≠ j → adj i j = true) := by
  decide

/-- Every edge lies in exactly two faces, as a kernel-checked list count. -/
theorem two_faces_receipt :
    edgesList.all
      (fun e => facesList.countP (fun f => decide (e ⊆ f)) == 2) = true := by
  decide

/-- Filter cardinality on the face set as a list count. -/
theorem faces_filter_card (p : Finset (Fin 12) → Prop) [DecidablePred p] :
    (faces.filter p).card = facesList.countP (fun f => decide (p f)) :=
  filter_card_coe p

/-- Every distinct pair inside a face forms an edge, as a kernel-checked
list fold. -/
theorem face_pairs_receipt :
    facesList.all (fun f =>
      decide (∀ i : Fin 12, ∀ j : Fin 12, i ∈ f → j ∈ f → i ≠ j →
        ({i, j} : Finset (Fin 12)) ∈ edges)) = true := by
  decide

/-- Every two-element subset of a face is an edge. -/
theorem face_edges_all :
    ∀ f ∈ faces, ∀ e ∈ Finset.powersetCard 2 f, e ∈ edges := by
  intro f hf e he
  obtain ⟨hef, h2⟩ := Finset.mem_powersetCard.mp he
  obtain ⟨x, y, hxy, rfl⟩ := Finset.card_eq_two.mp h2
  have h := of_decide_eq_true
    (List.all_eq_true.mp face_pairs_receipt f (by simpa [faces, Finset.mem_mk] using hf))
  exact h x y (hef (Finset.mem_insert_self _ _))
    (hef (Finset.mem_insert_of_mem (Finset.mem_singleton_self _))) hxy

/-- The icosahedral boundary complex on the committed incidence. -/
def carrier : BoundaryComplex (Fin 12) where
  ports := Finset.univ
  edges := edges
  faces := faces
  edge_card := ball_of_all (by decide)
  edge_sub_ports := fun e _ => Finset.subset_univ e
  face_card := ball_of_all (by decide)
  face_edges := face_edges_all
  edge_two_faces := by
    intro e he
    rw [faces_filter_card]
    exact beq_iff_eq.mp (List.all_eq_true.mp two_faces_receipt e
      (by simpa [edges, Finset.mem_mk] using he))

/-- Twelve ports, kernel-checked. -/
theorem ports_card : carrier.ports.card = 12 := by decide

/-- Thirty edges, kernel-checked. -/
theorem edges_card : carrier.edges.card = 30 := by decide

/-- Twenty faces, kernel-checked. -/
theorem faces_card : carrier.faces.card = 20 := by decide

/-- Degree five at every port, kernel-checked. -/
theorem degree_five_at :
    ∀ p : Fin 12, (carrier.edges.filter (fun e => p ∈ e)).card = 5 := by decide

/-- Degree five, in the bounded form used by the criteria interface. -/
theorem degree_five :
    ∀ p ∈ carrier.ports, (carrier.edges.filter (fun e => p ∈ e)).card = 5 :=
  fun p _ => degree_five_at p

/-- The Euler relation on the witness, kernel-checked. -/
theorem euler_two :
    (carrier.ports.card : ℤ) - carrier.edges.card + carrier.faces.card = 2 := by
  decide

/-- Membership in the witness edge set is exactly adjacency of a distinct
pair. -/
theorem pair_mem_edges_iff :
    ∀ i j : Fin 12,
      (({i, j} : Finset (Fin 12)) ∈ carrier.edges) ↔ i ≠ j ∧ adj i j = true := by
  decide

/-! ### Link five-cycles

`linkCycle i` lists the five neighbors of port `i` in cyclic order.  The
three theorems below certify the A1 clause that the neighbors of each port
form a five-cycle, and `linkCycle_faces` adds that consecutive link members
span the five faces at the port. -/

/-- A cyclic listing of the five neighbors of each port. -/
def linkCycle : Fin 12 → List (Fin 12)
  | 0 => [1, 2, 4, 6, 3]
  | 1 => [0, 2, 5, 7, 3]
  | 2 => [0, 1, 5, 8, 4]
  | 3 => [0, 1, 7, 9, 6]
  | 4 => [0, 2, 8, 10, 6]
  | 5 => [1, 2, 8, 11, 7]
  | 6 => [0, 3, 9, 10, 4]
  | 7 => [1, 3, 9, 11, 5]
  | 8 => [2, 4, 10, 11, 5]
  | 9 => [3, 6, 10, 11, 7]
  | 10 => [4, 6, 9, 11, 8]
  | 11 => [5, 7, 9, 10, 8]

/-- Each link listing has five distinct members. -/
theorem linkCycle_nodup : ∀ i : Fin 12, (linkCycle i).Nodup ∧ (linkCycle i).length = 5 := by
  decide

/-- The link listing of a port contains exactly its neighbors. -/
theorem linkCycle_covers : ∀ i j : Fin 12, j ∈ linkCycle i ↔ adj i j = true := by
  decide

/-- Consecutive link members are adjacent, including the wrap-around pair:
each link is a five-cycle. -/
theorem linkCycle_cyclic :
    ∀ i : Fin 12, ∀ p ∈ (linkCycle i).zip ((linkCycle i).rotate 1),
      adj p.1 p.2 = true := by
  decide

/-- Consecutive link members span a face with the port: the five faces at
each port, read off the link five-cycle. -/
theorem linkCycle_faces :
    ∀ i : Fin 12, ∀ p ∈ (linkCycle i).zip ((linkCycle i).rotate 1),
      ({i, p.1, p.2} : Finset (Fin 12)) ∈ carrier.faces := by
  decide

/-! ### Coherent orientation

The twenty faces carry an explicit orientation.  `coherent_orientation`
states the A1 clause exactly: every directed pair of adjacent ports occurs
exactly once among the sixty directed boundary pairs, so every edge lies in
two oppositely oriented faces, and non-adjacent pairs occur in no face. -/

/-- The twenty faces with a coherent orientation. -/
def orientedFaces : List (Fin 12 × Fin 12 × Fin 12) :=
  [(0, 1, 2), (0, 2, 4), (0, 4, 6), (1, 0, 3), (1, 3, 7),
   (2, 1, 5), (2, 5, 8), (3, 0, 6), (3, 6, 9), (4, 2, 8),
   (4, 8, 10), (5, 1, 7), (5, 7, 11), (6, 4, 10), (7, 3, 9),
   (7, 9, 11), (8, 5, 11), (8, 11, 10), (9, 6, 10), (9, 10, 11)]

/-- The three directed boundary pairs of an oriented face. -/
def boundaryPairs (t : Fin 12 × Fin 12 × Fin 12) : List (Fin 12 × Fin 12) :=
  [(t.1, t.2.1), (t.2.1, t.2.2), (t.2.2, t.1)]

/-- Twenty oriented faces, and each underlying set is a face of the
witness. -/
theorem orientedFaces_sub :
    orientedFaces.length = 20 ∧
      ∀ t ∈ orientedFaces,
        ({t.1, t.2.1, t.2.2} : Finset (Fin 12)) ∈ carrier.faces :=
  ⟨by decide, ball_list_of_all (by decide)⟩

/-- Every face of the witness is the underlying set of an oriented face. -/
theorem orientedFaces_complete :
    ∀ f ∈ carrier.faces, ∃ t ∈ orientedFaces,
      ({t.1, t.2.1, t.2.2} : Finset (Fin 12)) = f := by
  refine ball_of_all (l := facesList) ?_
  decide

/-- Coherence: a directed port pair occurs exactly once among the sixty
directed boundary pairs when the ports are adjacent, and never otherwise.
Each edge therefore lies in exactly two faces which traverse it in opposite
directions. -/
theorem coherent_orientation :
    ∀ i j : Fin 12,
      (orientedFaces.map boundaryPairs).flatten.count (i, j)
        = if adj i j then 1 else 0 := by
  decide

/-! ### Connectivity -/

/-- Every port is reachable from port `0` through edges. -/
theorem reach_zero :
    ∀ j : Fin 12,
      Relation.ReflTransGen
        (fun a b => ({a, b} : Finset (Fin 12)) ∈ carrier.edges) 0 j := by
  intro j
  fin_cases j
  · exact .refl
  · exact .single (by decide)
  · exact .single (by decide)
  · exact .single (by decide)
  · exact .single (by decide)
  · exact .head (b := 1) (by decide) (.single (by decide))
  · exact .single (by decide)
  · exact .head (b := 1) (by decide) (.single (by decide))
  · exact .head (b := 2) (by decide) (.single (by decide))
  · exact .head (b := 3) (by decide) (.single (by decide))
  · exact .head (b := 4) (by decide) (.single (by decide))
  · exact .head (b := 2) (by decide)
      (.head (b := 8) (by decide) (.single (by decide)))

/-- Any two ports are joined through edges. -/
theorem connected :
    ∀ p ∈ carrier.ports, ∀ q ∈ carrier.ports,
      Relation.ReflTransGen
        (fun a b => ({a, b} : Finset (Fin 12)) ∈ carrier.edges) p q := by
  intro p _ q _
  have hsym :
      Symmetric (fun a b : Fin 12 =>
        ({a, b} : Finset (Fin 12)) ∈ carrier.edges) := by
    intro a b h
    rwa [Finset.pair_comm]
  exact ((Relation.ReflTransGen.symmetric hsym) (reach_zero p)).trans
    (reach_zero q)

/-- The icosahedral witness satisfies the declared criteria. -/
def criteria : EchosahedralCriteria (Fin 12) where
  toBoundaryComplex := carrier
  ports_card := ports_card
  degree_five := degree_five
  connected := connected

/-- The abstract counting theorems and the kernel receipts agree on the
witness. -/
example : criteria.edges.card = 30 ∧ criteria.faces.card = 20 :=
  ⟨criteria.edges_thirty, criteria.faces_twenty⟩

end Icosa

/-! ## Partial rigidity: the edge set forces the face set -/

/-- Any boundary complex on the icosahedral edge set carries exactly the
icosahedral face set.  Faces are three-element sets whose two-element
subsets are edges, so every face is a triangle of the incidence; the
triangulation count forces twenty faces, and the incidence has exactly
twenty triangles, so the containment is an equality. -/
theorem faces_determined (K : BoundaryComplex (Fin 12))
    (hE : K.edges = Icosa.carrier.edges) :
    K.faces = Icosa.carrier.faces := by
  have hEcard : K.edges.card = 30 := by rw [hE]; exact Icosa.edges_card
  have h3 : 3 * K.faces.card = 2 * K.edges.card := triangulation_count K
  have hFcard : K.faces.card = 20 := by omega
  have hsub : K.faces ⊆ Icosa.carrier.faces := by
    intro f hf
    show f ∈ Icosa.faces
    rw [Icosa.faces_eq_triangles]
    refine Finset.mem_filter.mpr
      ⟨Finset.mem_powersetCard.mpr ⟨Finset.subset_univ f, K.face_card f hf⟩, ?_⟩
    intro i hi j hj hij
    have hmem : ({i, j} : Finset (Fin 12)) ∈ K.edges :=
      K.face_edges f hf _ (pair_mem_powersetCard hi hj hij)
    rw [hE] at hmem
    exact ((Icosa.pair_mem_edges_iff i j).mp hmem).2
  refine Finset.eq_of_subset_of_card_le hsub ?_
  rw [hFcard, Icosa.faces_card]

/-! ## Controls

Three exact controls locate the load-bearing criteria.  The octahedron
satisfies every complex law with degree four, the cube edge set admits no
boundary complex, and the hemi-icosahedron satisfies every law including
degree five with six ports and Euler characteristic one. -/

namespace Octa

/-- Octahedron adjacency on six ports: distinct and non-antipodal
(`i + j = 5` pairs antipodes). -/
def adjB : Fin 6 → Fin 6 → Bool :=
  fun i j => decide (i ≠ j) && decide ((i : ℕ) + (j : ℕ) ≠ 5)

/-- The eight octahedron faces: the three-element cliques of `adjB`. -/
def faces : Finset (Finset (Fin 6)) :=
  (Finset.powersetCard 3 Finset.univ).filter
    (fun f => ∀ i ∈ f, ∀ j ∈ f, i ≠ j → adjB i j = true)

/-- The octahedral boundary complex: a valid triangulated boundary. -/
def carrier : BoundaryComplex (Fin 6) :=
  ofData adjB faces (by decide) (by decide) (by decide)

/-- Six ports, twelve edges, eight faces, kernel-checked. -/
theorem counts :
    carrier.ports.card = 6 ∧ carrier.edges.card = 12 ∧ carrier.faces.card = 8 := by
  decide

/-- Every port has degree four, so the degree-five criterion fails while
every complex law holds: the degree criterion is load-bearing. -/
theorem degree_four :
    ∀ p ∈ carrier.ports, (carrier.edges.filter (fun e => p ∈ e)).card = 4 := by
  decide

/-- The octahedron is spherical: its Euler relation gives two. -/
theorem euler_two :
    (carrier.ports.card : ℤ) - carrier.edges.card + carrier.faces.card = 2 := by
  decide

end Octa

namespace Cube

/-- The twelve cube edges on eight ports (binary-coordinate labelling). -/
def edges : Finset (Finset (Fin 8)) :=
  {{0, 1}, {0, 2}, {1, 3}, {2, 3}, {4, 5}, {4, 6}, {5, 7}, {6, 7},
   {0, 4}, {1, 5}, {2, 6}, {3, 7}}

/-- Twelve edges and degree three at every port, kernel-checked. -/
theorem counts :
    edges.card = 12 ∧ ∀ p : Fin 8, (edges.filter (fun e => p ∈ e)).card = 3 := by
  decide

/-- The cube graph has no triangle, kernel-checked over all port triples. -/
theorem triangle_free :
    ∀ a b c : Fin 8, ({a, b} : Finset (Fin 8)) ∈ edges →
      ({a, c} : Finset (Fin 8)) ∈ edges →
      ({b, c} : Finset (Fin 8)) ∈ edges → False := by
  decide

/-- No boundary complex exists on the cube edge set: any face would be a
triangle of a triangle-free graph.  The triangulation criterion is
load-bearing, and the cube fails it before any counting starts. -/
theorem no_boundary_complex :
    ¬∃ K : BoundaryComplex (Fin 8), K.edges = edges := by
  rintro ⟨K, hE⟩
  have h01 : ({0, 1} : Finset (Fin 8)) ∈ K.edges := by rw [hE]; decide
  have h2 := K.edge_two_faces _ h01
  have hpos : 0 < (K.faces.filter
      (fun f => ({0, 1} : Finset (Fin 8)) ⊆ f)).card := by omega
  obtain ⟨f, hf⟩ := Finset.card_pos.mp hpos
  obtain ⟨hfF, -⟩ := Finset.mem_filter.mp hf
  obtain ⟨a, b, c, hab, hac, hbc, hfeq⟩ :=
    Finset.card_eq_three.mp (K.face_card f hfF)
  have ha : a ∈ f := by rw [hfeq]; exact Finset.mem_insert_self _ _
  have hb : b ∈ f := by
    rw [hfeq]; exact Finset.mem_insert_of_mem (Finset.mem_insert_self _ _)
  have hc : c ∈ f := by
    rw [hfeq]
    exact Finset.mem_insert_of_mem
      (Finset.mem_insert_of_mem (Finset.mem_singleton_self _))
  have e1 := K.face_edges f hfF _ (pair_mem_powersetCard ha hb hab)
  have e2 := K.face_edges f hfF _ (pair_mem_powersetCard ha hc hac)
  have e3 := K.face_edges f hfF _ (pair_mem_powersetCard hb hc hbc)
  rw [hE] at e1 e2 e3
  exact triangle_free a b c e1 e2 e3

end Cube

namespace Hemi

/-- Hemi-icosahedron adjacency: the complete graph on six ports. -/
def adjB : Fin 6 → Fin 6 → Bool := fun i j => decide (i ≠ j)

/-- The ten faces of the hemi-icosahedron, the antipodal quotient of the
icosahedral face set. -/
def faces : Finset (Finset (Fin 6)) :=
  {{0, 1, 2}, {0, 1, 3}, {0, 2, 4}, {0, 3, 5}, {0, 4, 5},
   {1, 2, 5}, {1, 3, 4}, {1, 4, 5}, {2, 3, 4}, {2, 3, 5}}

/-- The hemi-icosahedral boundary complex: a valid triangulated boundary. -/
def carrier : BoundaryComplex (Fin 6) :=
  ofData adjB faces (by decide) (by decide) (by decide)

/-- Degree five at every port: the hemi-icosahedron satisfies the degree
criterion. -/
theorem degree_five :
    ∀ p ∈ carrier.ports, (carrier.edges.filter (fun e => p ∈ e)).card = 5 := by
  decide

/-- Six ports, fifteen edges, ten faces, kernel-checked. -/
theorem counts :
    carrier.ports.card = 6 ∧ carrier.edges.card = 15 ∧ carrier.faces.card = 10 := by
  decide

/-- The hemi-icosahedron has Euler characteristic one, the projective-plane
value: every complex law and the degree criterion hold, and only the
twelve-port criterion excludes it.  The port-count criterion is load-bearing
and coincides with the choice of the spherical Euler value. -/
theorem euler_one :
    (carrier.ports.card : ℤ) - carrier.edges.card + carrier.faces.card = 1 := by
  decide

/-- Consistency with the abstract theorem: six ports give Euler
characteristic one via `six_mul_euler_eq_ports`. -/
example :
    6 * ((carrier.ports.card : ℤ) - carrier.edges.card + carrier.faces.card)
      = carrier.ports.card :=
  six_mul_euler_eq_ports carrier degree_five

end Hemi

/-! ## Axiom audit -/

#print axioms handshake_five
#print axioms triangulation_count
#print axioms six_mul_euler_eq_ports
#print axioms twelve_ports_iff_euler_two
#print axioms counts_forced
#print axioms spherical_counts
#print axioms EchosahedralCriteria.edges_thirty
#print axioms EchosahedralCriteria.faces_twenty
#print axioms EchosahedralCriteria.euler_two
#print axioms Icosa.edges_eq_edgesOf
#print axioms Icosa.faces_eq_triangles
#print axioms Icosa.degree_five
#print axioms Icosa.linkCycle_cyclic
#print axioms Icosa.linkCycle_faces
#print axioms Icosa.coherent_orientation
#print axioms Icosa.connected
#print axioms faces_determined
#print axioms Octa.degree_four
#print axioms Cube.no_boundary_complex
#print axioms Hemi.euler_one

end OPH.CarrierUniqueness
