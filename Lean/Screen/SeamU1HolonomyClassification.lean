import DiscreteGauss
import HolonomyInterference

namespace OPH.SeamU1HolonomyClassification

open OPH.SeamCurrentCarrierQuotient
open OPH.DiscreteGauss
open OPH.HolonomyInterference

/-!
# Finite U(1) holonomy classification of the seam graph (V3, issue #733)

A U(1)-valued seam connection assigns one unit-circle element to each of
the thirty canonically oriented icosahedral seams of the committed table
`seamLeft`/`seamRight`.  A port gauge transformation, one circle element
per port, acts on connections by the abelian rechart pattern of
`HolonomyInterference`: the label of a seam is conjugated by the two
endpoint factors, and `gaugeAct_eq_rechart` states this identification on
the seam alphabet exactly.

The theorems classify connections modulo port gauge at the stated finite
level.  Every connection is gauge-equivalent to exactly one
representative that is trivial on the eleven seams of the committed
spanning tree, the tree whose seam support carries the declared Gauss
section `rationalBoundarySection`
(`rationalBoundarySection_vanishes_off_tree`,
`existsUnique_treeTrivial_rep`).  Two connections are gauge-equivalent
exactly when their nineteen chord holonomies agree
(`gaugeEquiv_iff_chordHolonomy_eq`), every assignment of nineteen circle
elements is realized (`chordHolonomy_surjective`), and the induced map is
a group isomorphism from connections modulo the coboundary subgroup onto
the nineteen-fold power of the circle
(`seamConnectionQuotientMulEquiv`).  Each chord holonomy is the literal
closed-loop holonomy, in the sense of `HolonomyInterference.holonomy`, of
the based loop that runs out along the tree, across the chord, and back
(`chordHolonomy_eq_pathHolonomy`); for the abelian structure group these
loop holonomies are rechart-invariant
(`HolonomyInterference.holonomy_rechart_invariant`), and
`chordHolonomy_gaugeAct` records the invariance for the classification
map itself.

This is the discrete degree-one cohomology statement
`H^1(seam graph, U(1)) = U(1)^19` for the committed table: the moduli
count nineteen is the number of non-tree seams, thirty seams minus the
eleven tree seams, and it equals the committed source-free cycle rank of
the discrete Gauss receipts, the first Betti number `30 - 12 + 1 = 19`
(`chord_count_eq_gauss_cycle_rank`, `gauss_cycle_space_finrank`).  The
same committed number counts the source-free ambiguity of the Gauss
solution fibre and the finite U(1) moduli.

Premises.  Every theorem in this file is premise-free exact mathematics
on the committed seam table.  No field of `LightSignalPremiseData` and no
premise-register row is consumed; in particular the physical-attachment
rows PR-53 and PR-54 are not consumed and stay open.  A composed lane
receipt that carries these statements as clauses consumes zero fields of
its premise bundle through them, and any conditionality displayed there
is not contributed by this module.

Boundary.  The seam connection is a finite kinematic object.  No theorem
states curvature dynamics, an action, a source coupling, an
inhomogeneous equation, a Coulomb readout, a continuum limit, or an
identification of the connection with a laboratory electromagnetic
potential; those attachments are exactly the open rows named above.  The
gauge quotient and the tree-gauge representative are statements about
the incidence combinatorics of the thirty-seam table with unit-circle
coefficients and nothing else.

Axiom audit.  Every proof composes committed results with finite
decidable checks; the module adds no project axiom and uses no native
decision procedure.  The audit lines at the end of the file show at most
`propext`, `Classical.choice`, and `Quot.sound`.
-/

/-! ## Connections, port gauge, and the gauge action -/

/-- A U(1)-valued connection on the thirty canonically oriented seams. -/
abbrev SeamConnection := Fin 30 → Circle

/-- A port gauge transformation: one unit-circle element per port. -/
abbrev PortGauge := Fin 12 → Circle

/-- Port gauge action on a seam connection.  On the canonical orientation
of seam `e` the label is conjugated by the endpoint factors, which is the
`HolonomyInterference.rechart` pattern restricted to the seam alphabet
(`gaugeAct_eq_rechart`). -/
noncomputable def gaugeAct (g : PortGauge) (A : SeamConnection) : SeamConnection :=
  fun e => (g (seamLeft e))⁻¹ * A e * g (seamRight e)

/-- The identity gauge fixes every connection. -/
theorem gaugeAct_one (A : SeamConnection) : gaugeAct 1 A = A := by
  funext e
  simp [gaugeAct]

/-- Two gauge steps compose to the pointwise product gauge. -/
theorem gaugeAct_gaugeAct (g₂ g₁ : PortGauge) (A : SeamConnection) :
    gaugeAct g₂ (gaugeAct g₁ A) = gaugeAct (g₁ * g₂) A := by
  funext e
  simp [gaugeAct, Pi.mul_apply, mul_inv_rev, mul_assoc]

/-- Constant gauges act trivially: the global U(1) phase is in the
stabilizer of every connection. -/
theorem gaugeAct_const (c : Circle) (A : SeamConnection) :
    gaugeAct (fun _ => c) A = A := by
  funext e
  exact inv_mul_cancel_comm c (A e)

/-- Gauge equivalence of seam connections: one port gauge carries the
first connection to the second. -/
def GaugeEquiv (A B : SeamConnection) : Prop :=
  ∃ g : PortGauge, gaugeAct g A = B

theorem gaugeEquiv_refl (A : SeamConnection) : GaugeEquiv A A :=
  ⟨1, gaugeAct_one A⟩

theorem gaugeEquiv_symm {A B : SeamConnection} (h : GaugeEquiv A B) :
    GaugeEquiv B A := by
  obtain ⟨g, rfl⟩ := h
  refine ⟨g⁻¹, ?_⟩
  rw [gaugeAct_gaugeAct, mul_inv_cancel, gaugeAct_one]

theorem gaugeEquiv_trans {A B C : SeamConnection}
    (h₁ : GaugeEquiv A B) (h₂ : GaugeEquiv B C) : GaugeEquiv A C := by
  obtain ⟨g₁, rfl⟩ := h₁
  obtain ⟨g₂, rfl⟩ := h₂
  exact ⟨g₁ * g₂, (gaugeAct_gaugeAct g₂ g₁ A).symm⟩

/-! ## The committed spanning tree and its chords

The eleven tree seams are exactly the seam support of the declared Gauss
section `rationalBoundarySection`: rooted at port `0`, the tree edges are
seam `0 = (0,1)`, `1 = (0,2)`, `2 = (0,3)`, `3 = (0,4)`, `4 = (0,6)`,
`7 = (1,5)`, `8 = (1,7)`, `11 = (2,8)`, `14 = (3,9)`, `17 = (4,10)`, and
`20 = (5,11)`.  The nineteen remaining seams are the chords. -/

/-- The eleven seams of the committed spanning tree. -/
def treeSeams : Finset (Fin 30) := {0, 1, 2, 3, 4, 7, 8, 11, 14, 17, 20}

/-- The nineteen chords, the non-tree seams in increasing order. -/
def chordIndex : Fin 19 → Fin 30 :=
  ![5, 6, 9, 10, 12, 13, 15, 16, 18, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29]

/-- The committed tree has exactly eleven seams. -/
theorem treeSeams_card : treeSeams.card = 11 := by decide

/-- Exactly nineteen seams lie outside the committed tree. -/
theorem treeSeams_compl_card : (treeSeamsᶜ : Finset (Fin 30)).card = 19 := by
  decide

/-- Every seam is a tree seam or a chord. -/
theorem seam_tree_or_chord :
    ∀ e : Fin 30, e ∈ treeSeams ∨ ∃ i : Fin 19, chordIndex i = e := by decide

/-- The declared Gauss section `rationalBoundarySection` is supported on
the tree seams: its current vanishes on every chord.  This pins
`treeSeams` to the committed section rather than to a fresh choice. -/
theorem rationalBoundarySection_vanishes_off_tree
    (x : RationalPortLoad) (e : Fin 30) (he : e ∉ treeSeams) :
    rationalBoundarySection x e = 0 := by
  fin_cases e <;>
    first
      | exact absurd (by decide) he
      | simp [rationalBoundarySection]

/-- A connection in tree gauge: trivial on every tree seam. -/
def TreeTrivial (B : SeamConnection) : Prop :=
  ∀ e ∈ treeSeams, B e = 1

/-! ## Tree transport and the tree-gauge representative -/

/-- Transport of a connection from port `0` along the tree path to each
port: the ordered product of the connection labels on the tree seams of
that path.  Port by port this is the explicit committed path table. -/
noncomputable def treePotential (A : SeamConnection) : Fin 12 → Circle :=
  ![1, A 0, A 1, A 2, A 3, A 0 * A 7, A 4, A 0 * A 8,
    A 1 * A 11, A 2 * A 14, A 3 * A 17, A 0 * A 7 * A 20]

/-- The gauge built from inverse tree transport.  Applying it puts any
connection into tree gauge (`treeGaugeRep_treeTrivial`). -/
noncomputable def treeGauge (A : SeamConnection) : PortGauge :=
  fun p => (treePotential A p)⁻¹

/-- The tree-gauge representative of a connection. -/
noncomputable def treeGaugeRep (A : SeamConnection) : SeamConnection :=
  gaugeAct (treeGauge A) A

/-- The tree-gauge representative is trivial on every tree seam. -/
theorem treeGaugeRep_treeTrivial (A : SeamConnection) :
    TreeTrivial (treeGaugeRep A) := by
  intro e he
  fin_cases he <;>
    simp [treeGaugeRep, gaugeAct, treeGauge, treePotential,
      seamLeft, seamRight, mul_assoc]

/-- Every connection is gauge-equivalent to its tree-gauge
representative. -/
theorem gaugeEquiv_treeGaugeRep (A : SeamConnection) :
    GaugeEquiv A (treeGaugeRep A) :=
  ⟨treeGauge A, rfl⟩

/-! ## Chord holonomies -/

/-- The nineteen chord holonomies of a connection: for chord `i` with
endpoints `u < v`, the based loop transport out along the tree to `u`,
across the chord, and back along the tree from `v`.
`chordHolonomy_eq_pathHolonomy` identifies this with the literal
closed-loop holonomy of `HolonomyInterference`. -/
noncomputable def chordHolonomy (A : SeamConnection) : Fin 19 → Circle :=
  fun i =>
    treePotential A (seamLeft (chordIndex i)) * A (chordIndex i) *
      (treePotential A (seamRight (chordIndex i)))⁻¹

/-- On a chord, the tree-gauge representative displays the chord
holonomy itself. -/
theorem treeGaugeRep_chordIndex (A : SeamConnection) (i : Fin 19) :
    treeGaugeRep A (chordIndex i) = chordHolonomy A i := by
  simp [treeGaugeRep, gaugeAct, treeGauge, chordHolonomy, mul_assoc]

private theorem gaugeAct_based_shift (g : PortGauge) (A : SeamConnection) :
    gaugeAct g A = gaugeAct (fun p => (g 0)⁻¹ * g p) A := by
  have h := gaugeAct_gaugeAct g (fun _ : Fin 12 => (g 0)⁻¹) A
  rw [gaugeAct_const] at h
  exact h

private theorem treePotential_gaugeAct_based (g : PortGauge) (hg : g 0 = 1)
    (A : SeamConnection) (p : Fin 12) :
    treePotential (gaugeAct g A) p = treePotential A p * g p := by
  fin_cases p <;>
    simp [treePotential, gaugeAct, seamLeft, seamRight, hg, mul_assoc]

private theorem chordHolonomy_gaugeAct_based (g : PortGauge) (hg : g 0 = 1)
    (A : SeamConnection) :
    chordHolonomy (gaugeAct g A) = chordHolonomy A := by
  funext i
  unfold chordHolonomy
  rw [treePotential_gaugeAct_based g hg, treePotential_gaugeAct_based g hg]
  show _ * ((g (seamLeft (chordIndex i)))⁻¹ * A (chordIndex i) *
      g (seamRight (chordIndex i))) * _ = _
  simp [mul_inv_rev, mul_assoc]

/-- Chord holonomies are invariant under every port gauge
transformation.  This is the finite abelian holonomy invariance of
`HolonomyInterference.holonomy_rechart_invariant` instantiated on the
seam alphabet, proved here by tree-transport telescoping. -/
theorem chordHolonomy_gaugeAct (g : PortGauge) (A : SeamConnection) :
    chordHolonomy (gaugeAct g A) = chordHolonomy A := by
  rw [gaugeAct_based_shift g A]
  exact chordHolonomy_gaugeAct_based _ (by simp) A

private theorem treePotential_of_treeTrivial {B : SeamConnection}
    (hB : TreeTrivial B) (p : Fin 12) : treePotential B p = 1 := by
  have h0 := hB 0 (by decide)
  have h1 := hB 1 (by decide)
  have h2 := hB 2 (by decide)
  have h3 := hB 3 (by decide)
  have h4 := hB 4 (by decide)
  have h7 := hB 7 (by decide)
  have h8 := hB 8 (by decide)
  have h11 := hB 11 (by decide)
  have h14 := hB 14 (by decide)
  have h17 := hB 17 (by decide)
  have h20 := hB 20 (by decide)
  fin_cases p <;>
    simp [treePotential, h0, h1, h2, h3, h4, h7, h8, h11, h14, h17, h20]

/-- In tree gauge the chord holonomies are the chord labels
themselves. -/
theorem chordHolonomy_of_treeTrivial {B : SeamConnection}
    (hB : TreeTrivial B) (i : Fin 19) :
    chordHolonomy B i = B (chordIndex i) := by
  unfold chordHolonomy
  rw [treePotential_of_treeTrivial hB, treePotential_of_treeTrivial hB]
  simp

/-- A tree-gauge connection is determined by its chord holonomies. -/
theorem treeTrivial_eq_of_chordHolonomy_eq {B B' : SeamConnection}
    (hB : TreeTrivial B) (hB' : TreeTrivial B')
    (h : chordHolonomy B = chordHolonomy B') : B = B' := by
  funext e
  rcases seam_tree_or_chord e with he | ⟨i, rfl⟩
  · rw [hB e he, hB' e he]
  · rw [← chordHolonomy_of_treeTrivial hB,
      ← chordHolonomy_of_treeTrivial hB', h]

/-! ## The classification -/

/-- Two seam connections are gauge-equivalent exactly when their
nineteen chord holonomies agree. -/
theorem gaugeEquiv_iff_chordHolonomy_eq (A B : SeamConnection) :
    GaugeEquiv A B ↔ chordHolonomy A = chordHolonomy B := by
  constructor
  · rintro ⟨g, rfl⟩
    exact (chordHolonomy_gaugeAct g A).symm
  · intro h
    have hrep : treeGaugeRep A = treeGaugeRep B := by
      apply treeTrivial_eq_of_chordHolonomy_eq
        (treeGaugeRep_treeTrivial A) (treeGaugeRep_treeTrivial B)
      rw [show chordHolonomy (treeGaugeRep A) = chordHolonomy A from
          chordHolonomy_gaugeAct _ A,
        show chordHolonomy (treeGaugeRep B) = chordHolonomy B from
          chordHolonomy_gaugeAct _ B, h]
    have hA := gaugeEquiv_treeGaugeRep A
    rw [hrep] at hA
    exact gaugeEquiv_trans hA (gaugeEquiv_symm (gaugeEquiv_treeGaugeRep B))

/-- Every connection is gauge-equivalent to exactly one tree-gauge
representative. -/
theorem existsUnique_treeTrivial_rep (A : SeamConnection) :
    ∃! B : SeamConnection, TreeTrivial B ∧ GaugeEquiv A B := by
  refine ⟨treeGaugeRep A,
    ⟨treeGaugeRep_treeTrivial A, gaugeEquiv_treeGaugeRep A⟩, ?_⟩
  rintro B ⟨hBtriv, hAB⟩
  apply treeTrivial_eq_of_chordHolonomy_eq hBtriv (treeGaugeRep_treeTrivial A)
  have h₁ : chordHolonomy A = chordHolonomy B :=
    (gaugeEquiv_iff_chordHolonomy_eq A B).1 hAB
  have h₂ : chordHolonomy (treeGaugeRep A) = chordHolonomy A :=
    chordHolonomy_gaugeAct _ A
  rw [h₂, h₁]

/-- The tree-gauge connection with a prescribed chord assignment. -/
noncomputable def chordExtend (h : Fin 19 → Circle) : SeamConnection :=
  ![1, 1, 1, 1, 1, h 0, h 1, 1, 1, h 2, h 3, 1, h 4, h 5, 1, h 6, h 7,
    1, h 8, h 9, 1, h 10, h 11, h 12, h 13, h 14, h 15, h 16, h 17, h 18]

/-- The prescribed-chord connection is in tree gauge. -/
theorem chordExtend_treeTrivial (h : Fin 19 → Circle) :
    TreeTrivial (chordExtend h) := by
  intro e he
  fin_cases he <;> rfl

/-- The prescribed-chord connection realizes its prescription. -/
theorem chordHolonomy_chordExtend (h : Fin 19 → Circle) :
    chordHolonomy (chordExtend h) = h := by
  funext i
  rw [chordHolonomy_of_treeTrivial (chordExtend_treeTrivial h)]
  fin_cases i <;> rfl

/-- Every assignment of nineteen circle elements is the chord holonomy
vector of some connection. -/
theorem chordHolonomy_surjective : Function.Surjective chordHolonomy :=
  fun h => ⟨chordExtend h, chordHolonomy_chordExtend h⟩

/-! ## The rechart identification and typed loop holonomy -/

/-- Row lookup of the canonical seam between two ports, if one exists. -/
def seamBetween (u v : Fin 12) : Option (Fin 30) :=
  (List.finRange 30).find? fun e => decide (seamLeft e = u ∧ seamRight e = v)

set_option maxRecDepth 8192 in
private theorem seamBetween_canonical :
    ∀ e : Fin 30, seamBetween (seamLeft e) (seamRight e) = some e := by decide

set_option maxRecDepth 8192 in
private theorem seamBetween_reverse_none :
    ∀ e : Fin 30, seamBetween (seamRight e) (seamLeft e) = none := by decide

set_option maxRecDepth 8192 in
private theorem seamBetween_asymm :
    ∀ u v : Fin 12, seamBetween u v = none ∨ seamBetween v u = none := by decide

/-- The reversal-compatible edge-label extension of a seam connection to
ordered port pairs: the label on the canonical orientation, its inverse
on the reversed orientation, and the identity off the seam table.  This
is the interface consumed by the typed paths of `HolonomyInterference`. -/
noncomputable def connectionEdgeLabel (A : SeamConnection) (u v : Fin 12) : Circle :=
  match seamBetween u v with
  | some e => A e
  | none =>
    match seamBetween v u with
    | some e => (A e)⁻¹
    | none => 1

/-- The edge-label extension represents edge reversal by group inverse. -/
theorem connectionEdgeLabel_reversalCompatible (A : SeamConnection) :
    ReversalCompatible (connectionEdgeLabel A) := by
  intro u v
  rcases huv : seamBetween u v with _ | e₁
  · rcases hvu : seamBetween v u with _ | e₂
    · simp [connectionEdgeLabel, huv, hvu]
    · simp [connectionEdgeLabel, huv, hvu]
  · have hvu : seamBetween v u = none := by
      rcases seamBetween_asymm u v with h | h
      · rw [h] at huv
        exact absurd huv (by simp)
      · exact h
    simp [connectionEdgeLabel, huv, hvu]

/-- On the canonical orientation the extension returns the connection
label itself. -/
theorem connectionEdgeLabel_canonical (A : SeamConnection) (e : Fin 30) :
    connectionEdgeLabel A (seamLeft e) (seamRight e) = A e := by
  simp [connectionEdgeLabel, seamBetween_canonical]

/-- On the seam alphabet, the port gauge action is exactly the abelian
`rechart` of `HolonomyInterference` applied to the edge-label
extension. -/
theorem gaugeAct_eq_rechart (g : PortGauge) (A : SeamConnection) (e : Fin 30) :
    gaugeAct g A e = rechart g (connectionEdgeLabel A) (seamLeft e) (seamRight e) := by
  rw [rechart, connectionEdgeLabel_canonical]
  rfl

/-- The typed tree path from port `0` to each port, edge by edge along
the committed spanning tree. -/
def treePath : (p : Fin 12) → Path (Fin 12) 0 p
  | ⟨0, _⟩ => .nil 0
  | ⟨1, _⟩ => .snoc (.nil 0) 1
  | ⟨2, _⟩ => .snoc (.nil 0) 2
  | ⟨3, _⟩ => .snoc (.nil 0) 3
  | ⟨4, _⟩ => .snoc (.nil 0) 4
  | ⟨5, _⟩ => .snoc (.snoc (.nil 0) 1) 5
  | ⟨6, _⟩ => .snoc (.nil 0) 6
  | ⟨7, _⟩ => .snoc (.snoc (.nil 0) 1) 7
  | ⟨8, _⟩ => .snoc (.snoc (.nil 0) 2) 8
  | ⟨9, _⟩ => .snoc (.snoc (.nil 0) 3) 9
  | ⟨10, _⟩ => .snoc (.snoc (.nil 0) 4) 10
  | ⟨11, _⟩ => .snoc (.snoc (.snoc (.nil 0) 1) 5) 11

private theorem transport_treePath (A : SeamConnection) (p : Fin 12) :
    transport (connectionEdgeLabel A) (treePath p) = treePotential A p := by
  fin_cases p <;>
    simp [treePath, treePotential, connectionEdgeLabel,
      show seamBetween 0 1 = some 0 from by decide,
      show seamBetween 0 2 = some 1 from by decide,
      show seamBetween 0 3 = some 2 from by decide,
      show seamBetween 0 4 = some 3 from by decide,
      show seamBetween 0 6 = some 4 from by decide,
      show seamBetween 1 5 = some 7 from by decide,
      show seamBetween 1 7 = some 8 from by decide,
      show seamBetween 2 8 = some 11 from by decide,
      show seamBetween 3 9 = some 14 from by decide,
      show seamBetween 4 10 = some 17 from by decide,
      show seamBetween 5 11 = some 20 from by decide]

/-- The based loop of chord `i`: out along the tree to the left
endpoint, across the chord, and back along the tree from the right
endpoint. -/
def chordLoop (i : Fin 19) : Path (Fin 12) 0 0 :=
  Path.closeRatio
    (Path.snoc (treePath (seamLeft (chordIndex i))) (seamRight (chordIndex i)))
    (treePath (seamRight (chordIndex i)))

/-- Each chord holonomy is the literal closed-loop holonomy of the
edge-label extension along the chord loop, in the sense of
`HolonomyInterference.holonomy`. -/
theorem chordHolonomy_eq_pathHolonomy (A : SeamConnection) (i : Fin 19) :
    chordHolonomy A i = holonomy (connectionEdgeLabel A) (chordLoop i) := by
  rw [chordLoop, ← transportRatio_eq_closedLoopHolonomy (connectionEdgeLabel A)
    (connectionEdgeLabel_reversalCompatible A)]
  rw [transport_snoc, transport_treePath, transport_treePath,
    connectionEdgeLabel_canonical]
  rfl

/-! ## The moduli count and the Gauss cycle rank -/

/-- The finite U(1) moduli count equals the committed source-free cycle
rank of the discrete Gauss receipts: nineteen, the number of chords,
thirty seams minus the eleven tree seams. -/
theorem chord_count_eq_gauss_cycle_rank :
    ((treeSeamsᶜ : Finset (Fin 30)).card : ℕ)
      = Module.finrank ℚ (LinearMap.ker rationalSeamBoundary) := by
  rw [gauss_cycle_space_finrank, treeSeams_compl_card]

/-! ## The group isomorphism form -/

private theorem treePotential_one (p : Fin 12) :
    treePotential (1 : SeamConnection) p = 1 := by
  fin_cases p <;> simp [treePotential]

private theorem treePotential_mul (A B : SeamConnection) (p : Fin 12) :
    treePotential (A * B) p = treePotential A p * treePotential B p := by
  fin_cases p <;>
    simp [treePotential, Pi.mul_apply, mul_comm, mul_left_comm, mul_assoc]

private theorem mul_mul_mul_inv_comm (a₁ b₁ a₂ b₂ a₃ b₃ : Circle) :
    (a₁ * b₁) * (a₂ * b₂) * (a₃ * b₃)⁻¹ =
      (a₁ * a₂ * a₃⁻¹) * (b₁ * b₂ * b₃⁻¹) := by
  rw [mul_inv, mul_mul_mul_comm a₁ b₁ a₂ b₂, mul_mul_mul_comm (a₁ * a₂) (b₁ * b₂)]

/-- The coboundary homomorphism from port gauges to seam connections.
Its range is the subgroup of pure-gauge connections, and the gauge
orbits are its cosets (`gaugeAct_eq_mul_coboundary`). -/
noncomputable def coboundary : PortGauge →* (Fin 30 → Circle) where
  toFun g := fun e => (g (seamLeft e))⁻¹ * g (seamRight e)
  map_one' := by
    funext e
    simp
  map_mul' g h := by
    funext e
    simp only [Pi.mul_apply]
    rw [mul_inv]
    exact mul_mul_mul_comm _ _ _ _

/-- The gauge action is translation by the coboundary: orbits of the
port gauge action are exactly the cosets of `coboundary.range`. -/
theorem gaugeAct_eq_mul_coboundary (g : PortGauge) (A : SeamConnection) :
    gaugeAct g A = A * coboundary g := by
  funext e
  show (g (seamLeft e))⁻¹ * A e * g (seamRight e)
      = A e * ((g (seamLeft e))⁻¹ * g (seamRight e))
  rw [mul_comm ((g (seamLeft e))⁻¹) (A e), mul_assoc]

/-- The chord holonomy map as a group homomorphism into the
nineteen-fold power of the circle. -/
noncomputable def chordHolonomyHom : (Fin 30 → Circle) →* (Fin 19 → Circle) where
  toFun := chordHolonomy
  map_one' := by
    funext i
    simp [chordHolonomy, treePotential_one]
  map_mul' A B := by
    funext i
    simp only [Pi.mul_apply, chordHolonomy, treePotential_mul]
    exact mul_mul_mul_inv_comm _ _ _ _ _ _

/-- The pure-gauge subgroup is exactly the kernel of the chord holonomy
homomorphism. -/
theorem coboundary_range_eq_chordHolonomyHom_ker :
    coboundary.range = chordHolonomyHom.ker := by
  ext A
  rw [MonoidHom.mem_range, MonoidHom.mem_ker]
  constructor
  · rintro ⟨g, rfl⟩
    have hg : (coboundary g : Fin 30 → Circle) = gaugeAct g 1 := by
      rw [gaugeAct_eq_mul_coboundary, one_mul]
    show chordHolonomy (coboundary g) = 1
    rw [hg, chordHolonomy_gaugeAct]
    exact chordHolonomyHom.map_one
  · intro hA
    refine ⟨(treeGauge A)⁻¹, ?_⟩
    have hrep : A * coboundary (treeGauge A) = 1 := by
      rw [← gaugeAct_eq_mul_coboundary]
      funext e
      rcases seam_tree_or_chord e with he | ⟨i, rfl⟩
      · exact treeGaugeRep_treeTrivial A e he
      · rw [show gaugeAct (treeGauge A) A (chordIndex i) = chordHolonomy A i from
          treeGaugeRep_chordIndex A i]
        exact congrFun hA i
    rw [map_inv]
    exact inv_eq_of_mul_eq_one_left hrep

/-- The finite `H^1` statement in group form: seam connections modulo
the pure-gauge subgroup are isomorphic, as a group, to the nineteen-fold
power of the circle, indexed by the chords of the committed tree. -/
noncomputable def seamConnectionQuotientMulEquiv :
    ((Fin 30 → Circle) ⧸ coboundary.range) ≃* (Fin 19 → Circle) :=
  (QuotientGroup.quotientMulEquivOfEq coboundary_range_eq_chordHolonomyHom_ker).trans
    (QuotientGroup.quotientKerEquivOfSurjective chordHolonomyHom
      chordHolonomy_surjective)

/-- The isomorphism sends the class of a connection to its chord
holonomy vector. -/
theorem seamConnectionQuotientMulEquiv_mk (A : SeamConnection) :
    seamConnectionQuotientMulEquiv (QuotientGroup.mk A) = chordHolonomy A := rfl

/-- Classes in the quotient coincide exactly on gauge orbits, so the
quotient by `coboundary.range` is the set of connections modulo port
gauge. -/
theorem quotient_mk_eq_iff_gaugeEquiv (A B : SeamConnection) :
    (QuotientGroup.mk A : (Fin 30 → Circle) ⧸ coboundary.range)
        = QuotientGroup.mk B ↔ GaugeEquiv A B := by
  rw [QuotientGroup.eq]
  constructor
  · rintro ⟨g, hg⟩
    refine ⟨g, ?_⟩
    rw [gaugeAct_eq_mul_coboundary, hg, mul_inv_cancel_left]
  · rintro ⟨g, rfl⟩
    refine ⟨g, ?_⟩
    rw [gaugeAct_eq_mul_coboundary, inv_mul_cancel_left]

/-! ## The classification package -/

/-- One citable statement of the finite U(1) holonomy classification:
unique tree-gauge representatives, gauge equivalence as chord holonomy
agreement, realization of every chord assignment, and the identity of
the moduli count with the committed Gauss cycle rank. -/
theorem seam_u1_holonomy_classification :
    (∀ A : SeamConnection, ∃! B : SeamConnection, TreeTrivial B ∧ GaugeEquiv A B) ∧
    (∀ A B : SeamConnection, GaugeEquiv A B ↔ chordHolonomy A = chordHolonomy B) ∧
    Function.Surjective chordHolonomy ∧
    ((treeSeamsᶜ : Finset (Fin 30)).card : ℕ)
      = Module.finrank ℚ (LinearMap.ker rationalSeamBoundary) :=
  ⟨existsUnique_treeTrivial_rep, gaugeEquiv_iff_chordHolonomy_eq,
    chordHolonomy_surjective, chord_count_eq_gauss_cycle_rank⟩

end OPH.SeamU1HolonomyClassification

/- Axiom audit: finite decidable table checks, the committed incidence
and Gauss receipts, and standard Mathlib group and quotient lemmas only. -/

#print axioms OPH.SeamU1HolonomyClassification.treeSeams_card
#print axioms OPH.SeamU1HolonomyClassification.treeSeams_compl_card
#print axioms OPH.SeamU1HolonomyClassification.seam_tree_or_chord
#print axioms OPH.SeamU1HolonomyClassification.rationalBoundarySection_vanishes_off_tree
#print axioms OPH.SeamU1HolonomyClassification.gaugeAct_one
#print axioms OPH.SeamU1HolonomyClassification.gaugeAct_gaugeAct
#print axioms OPH.SeamU1HolonomyClassification.gaugeAct_const
#print axioms OPH.SeamU1HolonomyClassification.gaugeEquiv_refl
#print axioms OPH.SeamU1HolonomyClassification.gaugeEquiv_symm
#print axioms OPH.SeamU1HolonomyClassification.gaugeEquiv_trans
#print axioms OPH.SeamU1HolonomyClassification.treeGaugeRep_treeTrivial
#print axioms OPH.SeamU1HolonomyClassification.gaugeEquiv_treeGaugeRep
#print axioms OPH.SeamU1HolonomyClassification.treeGaugeRep_chordIndex
#print axioms OPH.SeamU1HolonomyClassification.chordHolonomy_gaugeAct
#print axioms OPH.SeamU1HolonomyClassification.chordHolonomy_of_treeTrivial
#print axioms OPH.SeamU1HolonomyClassification.treeTrivial_eq_of_chordHolonomy_eq
#print axioms OPH.SeamU1HolonomyClassification.gaugeEquiv_iff_chordHolonomy_eq
#print axioms OPH.SeamU1HolonomyClassification.existsUnique_treeTrivial_rep
#print axioms OPH.SeamU1HolonomyClassification.chordExtend_treeTrivial
#print axioms OPH.SeamU1HolonomyClassification.chordHolonomy_chordExtend
#print axioms OPH.SeamU1HolonomyClassification.chordHolonomy_surjective
#print axioms OPH.SeamU1HolonomyClassification.connectionEdgeLabel_reversalCompatible
#print axioms OPH.SeamU1HolonomyClassification.connectionEdgeLabel_canonical
#print axioms OPH.SeamU1HolonomyClassification.gaugeAct_eq_rechart
#print axioms OPH.SeamU1HolonomyClassification.chordHolonomy_eq_pathHolonomy
#print axioms OPH.SeamU1HolonomyClassification.chord_count_eq_gauss_cycle_rank
#print axioms OPH.SeamU1HolonomyClassification.gaugeAct_eq_mul_coboundary
#print axioms OPH.SeamU1HolonomyClassification.coboundary_range_eq_chordHolonomyHom_ker
#print axioms OPH.SeamU1HolonomyClassification.seamConnectionQuotientMulEquiv
#print axioms OPH.SeamU1HolonomyClassification.seamConnectionQuotientMulEquiv_mk
#print axioms OPH.SeamU1HolonomyClassification.quotient_mk_eq_iff_gaugeEquiv
#print axioms OPH.SeamU1HolonomyClassification.seam_u1_holonomy_classification
