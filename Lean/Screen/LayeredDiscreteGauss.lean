import RegionalContinuity

open scoped BigOperators

namespace OPH.LayeredDiscreteGauss

open OPH.RepairWordCarrierReadout
open OPH.SeamCurrentCarrierQuotient

/-!
# Layered discrete Gauss bookkeeping on a generic finite graph

V3 gravitation lane (issue #729), discrete companion row to the
continuum shell law of `Geometry/InverseSquareShellLaw.lean`. On a
generic finite graph with directed edge orientations and an integer
layer function, the module proves one exact rational identity. A flux
field drains a declared source supported in the layer-zero ball
(`steadyGauss` together with `sourceSupport`, so divergence-freedom
holds outside the source region on every constrained layer). A declared
shell-cardinality law places `shellCoeff * n ^ 2` boundary edges at
shell `n`, and a declared shell-equidistribution premise puts the same
outward flux on every boundary edge of one shell. Under these premises
the outward flux per shell boundary edge at layer `n` equals
`charge / (shellCoeff * n ^ 2)` exactly over the rationals.

The premises constrain shells up to a declared resolution `depth`, and
the truncation is itself forced, with the obstruction recorded at proof
level: a steady drained source on the whole closed carrier has zero
total load (`steady_total_source_eq_zero`), and an `n ^ 2` shell law
imposed at every positive layer forces a zero coefficient on a finite
carrier (`unbounded_shellCard_forces_zero`). Layers beyond `depth`
absorb the outward flux, standing in for the escape to spatial infinity
that a finite carrier does not contain. A three-vertex chain witness
(`chainWitness`) records that the depth-bounded premise set is
satisfiable with a nonzero coefficient, so the shell theorems have
instances.

The content of the main theorem is exact bookkeeping from the declared
premises to the shell exponent; the bridge theorems below record where
the committed conservation receipts make contact. The equidistribution and cardinality premises do
nearly all the work, and the `n ^ 2` law itself is not derived from the
carrier: the committed `20 * n ^ 2` count (`faceCount` in
`DiscreteRefinement`) is a subdivision count on one fixed screen, and
identifying subdivision frequency with radial layer distance is exactly
the open premise of the discrete route. The premises enter as named
fields of `LayeredFluxData`, declared premises of the V3 program,
parallel to the radial readout and shell-flux normalization fields of
`RadialFluxData` in the continuum module. This module exists to make
the discrete route's assumptions inspectable, and each theorem below
names the premises it consumes.

The committed regional divergence theorem
(`RegionalContinuity.seamBoundary_region_eq_neg_flux`) is stated on the
fixed twelve-port, thirty-seam carrier and does not transfer to a
generic layered graph, so the generic divergence theorem here is
developed standalone. Two bridge theorems record where the committed
receipts genuinely apply: `regionFlux_twelvePort` identifies the
generic outward flux with the committed `regionFlux` on the twelve-port
carrier after the integer-to-rational cast, and
`steady_witness_regionFlux_eq_source` derives the fixed-carrier
instance of the `steadyGauss` premise from the committed
`regional_continuity` theorem alone.

No theorem in this file identifies the layer index with a physical
radius, the rational flux with a gravitational field, or the source
load with a mass density.
-/

section GenericGraph

variable {V E : Type*} [Fintype V] [DecidableEq V] [Fintype E]

/-- Stored rational load in a finite set of vertices. Generic analogue
of `RegionalContinuity.regionLoad`, which is fixed to the twelve-port
carrier. -/
def regionLoad (region : Finset V) (q : V → ℚ) : ℚ :=
  ∑ p ∈ region, q p

/-- Outward flux carried by one oriented edge across a region boundary.
Positive flux points from `left e` to `right e`; edges with both
endpoints inside or both endpoints outside contribute zero. -/
def edgeOutwardFlux (left right : E → V) (region : Finset V) (e : E)
    (n : ℚ) : ℚ :=
  (if left e ∈ region then n else 0) - (if right e ∈ region then n else 0)

/-- Total signed outward flux through the boundary of a vertex region. -/
def regionFlux (left right : E → V) (region : Finset V) (flux : E → ℚ) : ℚ :=
  ∑ e : E, edgeOutwardFlux left right region e (flux e)

/-- Net outflow at one vertex: flux on edges leaving the vertex minus
flux on edges arriving at it. -/
def vertexOutflow (left right : E → V) (flux : E → ℚ) (v : V) : ℚ :=
  (∑ e : E, if left e = v then flux e else 0) -
    (∑ e : E, if right e = v then flux e else 0)

omit [Fintype V] in
/-- Regional divergence theorem on a generic finite graph: the sum of
vertex outflows over a region is its total outward boundary flux.
Internal edge contributions cancel exactly. Developed standalone
because the committed `RegionalContinuity` version is fixed to the
twelve-port carrier. -/
theorem sum_vertexOutflow_eq_regionFlux (left right : E → V)
    (region : Finset V) (flux : E → ℚ) :
    (∑ v ∈ region, vertexOutflow left right flux v) =
      regionFlux left right region flux := by
  classical
  unfold vertexOutflow regionFlux edgeOutwardFlux
  calc
    (∑ v ∈ region,
        ((∑ e : E, if left e = v then flux e else 0) -
          (∑ e : E, if right e = v then flux e else 0))) =
        ∑ v ∈ region, ∑ e : E,
          ((if left e = v then flux e else 0) -
            (if right e = v then flux e else 0)) := by
      apply Finset.sum_congr rfl
      intro v _
      rw [Finset.sum_sub_distrib]
    _ = ∑ e : E, ∑ v ∈ region,
          ((if left e = v then flux e else 0) -
            (if right e = v then flux e else 0)) := Finset.sum_comm
    _ = ∑ e : E,
          ((if left e ∈ region then flux e else 0) -
            (if right e ∈ region then flux e else 0)) := by
      apply Finset.sum_congr rfl
      intro e _
      rw [Finset.sum_sub_distrib]
      simp [Finset.sum_ite_eq]

/-- No edge crosses the boundary of the complete vertex region. -/
theorem regionFlux_univ (left right : E → V) (flux : E → ℚ) :
    regionFlux left right Finset.univ flux = 0 := by
  simp [regionFlux, edgeOutwardFlux]

/-- The layered ball: vertices with layer index at most `n`. -/
def layerBall (layer : V → ℤ) (n : ℤ) : Finset V :=
  Finset.univ.filter fun v => layer v ≤ n

omit [DecidableEq V] in
/-- Ball membership is the layer inequality. -/
theorem mem_layerBall {layer : V → ℤ} {n : ℤ} {v : V} :
    v ∈ layerBall layer n ↔ layer v ≤ n := by
  simp [layerBall]

/-- Boundary edges of shell `n`: edges with one endpoint in the layered
ball at `n` and the other endpoint outside, in either orientation. -/
def shellEdges (left right : E → V) (layer : V → ℤ) (n : ℤ) : Finset E :=
  Finset.univ.filter fun e =>
    (layer (left e) ≤ n ∧ n < layer (right e)) ∨
      (layer (right e) ≤ n ∧ n < layer (left e))

/-- The outward flux through a layered ball is carried entirely by its
shell boundary edges: every other edge contributes zero. -/
theorem regionFlux_layerBall_eq_shell_sum (left right : E → V)
    (layer : V → ℤ) (n : ℤ) (flux : E → ℚ) :
    regionFlux left right (layerBall layer n) flux =
      ∑ e ∈ shellEdges left right layer n,
        edgeOutwardFlux left right (layerBall layer n) e (flux e) := by
  classical
  unfold regionFlux shellEdges
  refine (Finset.sum_filter_of_ne ?_).symm
  intro e _ hne
  by_contra hcross
  push Not at hcross
  obtain ⟨himp₁, himp₂⟩ := hcross
  apply hne
  unfold edgeOutwardFlux
  by_cases hl : layer (left e) ≤ n
  · have hr : layer (right e) ≤ n := himp₁ hl
    simp [layerBall, hl, hr]
  · have hr : ¬layer (right e) ≤ n := fun hr => hl (himp₂ hr)
    simp [layerBall, hl, hr]

/-- Closed-carrier balance: total vertex outflow over the whole finite
graph vanishes. -/
theorem sum_vertexOutflow_univ (left right : E → V) (flux : E → ℚ) :
    (∑ v : V, vertexOutflow left right flux v) = 0 := by
  rw [sum_vertexOutflow_eq_regionFlux]
  exact regionFlux_univ left right flux

/-- Obstruction one: a steady drained source on the whole closed
carrier has zero total load. This identity forces the `depth` bound on
the `steadyGauss` premise of `LayeredFluxData`: a finite carrier
contains no spatial infinity for the flux to escape to, so the layers
beyond `depth` are left unconstrained and absorb it. -/
theorem steady_total_source_eq_zero (left right : E → V) (flux : E → ℚ)
    (source : V → ℚ)
    (hsteady : ∀ v : V, vertexOutflow left right flux v = source v) :
    (∑ v : V, source v) = 0 := by
  calc
    (∑ v : V, source v) = ∑ v : V, vertexOutflow left right flux v :=
      Finset.sum_congr rfl fun v _ => (hsteady v).symm
    _ = 0 := sum_vertexOutflow_univ left right flux

omit [DecidableEq V] in
/-- Obstruction two: an `n ^ 2` shell-cardinality law imposed at every
positive layer forces a zero coefficient on a finite carrier, because
the layer function is bounded and large shells are empty. This identity
forces the `depth` bound on the `shellCard` premise of
`LayeredFluxData`. -/
theorem unbounded_shellCard_forces_zero (left right : E → V)
    (layer : V → ℤ) (c : ℕ)
    (h : ∀ n : ℕ, 0 < n →
      (shellEdges left right layer (n : ℤ)).card = c * n ^ 2) :
    c = 0 := by
  classical
  obtain ⟨m, hmpos, hm⟩ : ∃ m : ℕ, 0 < m ∧ ∀ v : V, layer v < (m : ℤ) := by
    rcases isEmpty_or_nonempty V with hV | hV
    · exact ⟨1, one_pos, fun v => hV.elim v⟩
    · obtain ⟨M, hM⟩ := Finset.exists_le (Finset.univ.image layer)
      refine ⟨M.toNat + 1, Nat.succ_pos _, fun v => ?_⟩
      calc
        layer v ≤ M :=
          hM _ (Finset.mem_image_of_mem layer (Finset.mem_univ v))
        _ ≤ (M.toNat : ℤ) := Int.self_le_toNat M
        _ < ((M.toNat + 1 : ℕ) : ℤ) := by push_cast; omega
  have hempty : shellEdges left right layer (m : ℤ) = ∅ := by
    unfold shellEdges
    apply Finset.filter_false_of_mem
    intro e _
    rintro (⟨_, hgt⟩ | ⟨_, hgt⟩) <;> exact lt_asymm hgt (hm _)
  have hzero : c * m ^ 2 = 0 := by
    have hcard := h m hmpos
    rw [hempty] at hcard
    simpa using hcard.symm
  rcases Nat.mul_eq_zero.mp hzero with hc | hmsq
  · exact hc
  · exact absurd hmsq (pow_ne_zero 2 hmpos.ne')

end GenericGraph

/-- The declared premises of the discrete companion row, bundled with
the graph data they constrain. The four premise fields are declared
premises of the V3 program; no source theorem produces any of them on a
generic layered graph. -/
structure LayeredFluxData (V E : Type*) [Fintype V] [DecidableEq V]
    [Fintype E] where
  /-- Orientation tail of each edge. -/
  left : E → V
  /-- Orientation head of each edge. -/
  right : E → V
  /-- Integer layer index of each vertex. -/
  layer : V → ℤ
  /-- Rational flux on each oriented edge. -/
  flux : E → ℚ
  /-- Rational source load on each vertex. -/
  source : V → ℚ
  /-- Declared shell-cardinality coefficient. -/
  shellCoeff : ℕ
  /-- Declared per-edge outward flux of each shell. -/
  perEdgeFlux : ℕ → ℚ
  /-- Declared shell resolution: the premises constrain layers up to
  `depth` only. The truncation is forced on a finite carrier by
  `steady_total_source_eq_zero` and `unbounded_shellCard_forces_zero`;
  layers beyond `depth` absorb the outward flux. -/
  depth : ℕ
  /-- Declared premise: the source is supported in the layer-zero
  ball. -/
  sourceSupport : ∀ v : V, 0 < layer v → source v = 0
  /-- Declared premise: the flux field drains the source in steady
  state on every constrained layer, so vertex outflow equals the source
  there and the flux is divergence-free outside the source region up to
  `depth`. The twelve-port instance of this premise is a theorem,
  `steady_witness_regionFlux_eq_source` below. -/
  steadyGauss : ∀ v : V, layer v ≤ (depth : ℤ) →
    vertexOutflow left right flux v = source v
  /-- Declared premise: shell `n` has `shellCoeff * n ^ 2` boundary
  edges for `n` up to `depth`. The exponent is declared here, with no
  derivation from the committed carrier. -/
  shellCard : ∀ n : ℕ, 0 < n → n ≤ depth →
    (shellEdges left right layer (n : ℤ)).card = shellCoeff * n ^ 2
  /-- Declared premise: within shell `n`, for `n` up to `depth`, every
  boundary edge carries the same outward flux `perEdgeFlux n`, signed
  outward through `edgeOutwardFlux`. -/
  equidistributed : ∀ n : ℕ, 0 < n → n ≤ depth →
    ∀ e ∈ shellEdges left right layer (n : ℤ),
      edgeOutwardFlux left right (layerBall layer (n : ℤ)) e (flux e) =
        perEdgeFlux n

namespace LayeredFluxData

variable {V E : Type*} [Fintype V] [DecidableEq V] [Fintype E]

/-- Total declared load: the source summed over the layer-zero ball. -/
def charge (D : LayeredFluxData V E) : ℚ :=
  regionLoad (layerBall D.layer 0) D.source

/-- Discrete Gauss step: every layered ball at nonnegative layer index
up to the declared depth carries the full outward flux `charge`. This
is a composition of the generic divergence theorem with the
`steadyGauss` and `sourceSupport` premises. -/
theorem regionFlux_eq_charge (D : LayeredFluxData V E) {n : ℤ}
    (hn : 0 ≤ n) (hd : n ≤ (D.depth : ℤ)) :
    regionFlux D.left D.right (layerBall D.layer n) D.flux = D.charge := by
  classical
  calc
    regionFlux D.left D.right (layerBall D.layer n) D.flux =
        ∑ v ∈ layerBall D.layer n, vertexOutflow D.left D.right D.flux v :=
      (sum_vertexOutflow_eq_regionFlux D.left D.right _ D.flux).symm
    _ = ∑ v ∈ layerBall D.layer n, D.source v := by
      apply Finset.sum_congr rfl
      intro v hv
      rw [mem_layerBall] at hv
      exact D.steadyGauss v (le_trans hv hd)
    _ = ∑ v ∈ layerBall D.layer 0, D.source v := by
      symm
      apply Finset.sum_subset
      · intro v hv
        rw [mem_layerBall] at hv ⊢
        exact le_trans hv hn
      · intro v _ hvnot
        rw [mem_layerBall] at hvnot
        exact D.sourceSupport v (not_le.mp hvnot)
    _ = D.charge := rfl

/-- Under the `equidistributed` premise the shell boundary sum
collapses to a cardinality multiple of the per-edge flux. -/
theorem shell_sum_eq_card_mul (D : LayeredFluxData V E) {n : ℕ}
    (hn : 0 < n) (hd : n ≤ D.depth) :
    (∑ e ∈ shellEdges D.left D.right D.layer (n : ℤ),
        edgeOutwardFlux D.left D.right (layerBall D.layer (n : ℤ)) e
          (D.flux e)) =
      ((shellEdges D.left D.right D.layer (n : ℤ)).card : ℚ) *
        D.perEdgeFlux n := by
  calc
    (∑ e ∈ shellEdges D.left D.right D.layer (n : ℤ),
        edgeOutwardFlux D.left D.right (layerBall D.layer (n : ℤ)) e
          (D.flux e)) =
        ∑ _e ∈ shellEdges D.left D.right D.layer (n : ℤ), D.perEdgeFlux n :=
      Finset.sum_congr rfl fun e he => D.equidistributed n hn hd e he
    _ = ((shellEdges D.left D.right D.layer (n : ℤ)).card : ℚ) *
        D.perEdgeFlux n := by
      rw [Finset.sum_const, nsmul_eq_mul]

/-- Conservation across shells: shell cardinality times per-edge flux
is the total load on every positive shell. This is the composition of
the Gauss step, the shell localization, the `shellCard` premise, and
the `equidistributed` premise. -/
theorem shell_total_eq_charge (D : LayeredFluxData V E) {n : ℕ}
    (hn : 0 < n) (hd : n ≤ D.depth) :
    ((D.shellCoeff : ℚ) * (n : ℚ) ^ 2) * D.perEdgeFlux n = D.charge := by
  have hflux := D.regionFlux_eq_charge (n := (n : ℤ)) (Int.natCast_nonneg n)
    (by exact_mod_cast hd)
  rw [regionFlux_layerBall_eq_shell_sum, D.shell_sum_eq_card_mul hn hd,
    D.shellCard n hn hd] at hflux
  push_cast at hflux
  linear_combination hflux

/-- The discrete shell law: under the declared premises the outward
flux per shell boundary edge at layer `n` is exactly the total load
over `shellCoeff * n ^ 2`, as rational numbers. The division step is
pure bookkeeping; the exponent comes from the `shellCard` premise and
the uniform value from the `equidistributed` premise. -/
theorem perEdgeFlux_eq (D : LayeredFluxData V E) {n : ℕ} (hn : 0 < n)
    (hd : n ≤ D.depth) (hc : D.shellCoeff ≠ 0) :
    D.perEdgeFlux n = D.charge / ((D.shellCoeff : ℚ) * (n : ℚ) ^ 2) := by
  have htotal := D.shell_total_eq_charge hn hd
  have hden : ((D.shellCoeff : ℚ) * (n : ℚ) ^ 2) ≠ 0 :=
    mul_ne_zero (Nat.cast_ne_zero.mpr hc)
      (pow_ne_zero 2 (Nat.cast_ne_zero.mpr hn.ne'))
  rw [eq_div_iff hden]
  linear_combination htotal

/-- Scale-free form: per-edge flux times the squared layer index is
shell-independent under the declared premises. -/
theorem scale_free (D : LayeredFluxData V E) {n₁ n₂ : ℕ} (h₁ : 0 < n₁)
    (h₂ : 0 < n₂) (hd₁ : n₁ ≤ D.depth) (hd₂ : n₂ ≤ D.depth)
    (hc : D.shellCoeff ≠ 0) :
    D.perEdgeFlux n₁ * (n₁ : ℚ) ^ 2 = D.perEdgeFlux n₂ * (n₂ : ℚ) ^ 2 := by
  have e₁ := D.shell_total_eq_charge h₁ hd₁
  have e₂ := D.shell_total_eq_charge h₂ hd₂
  have hcq : (D.shellCoeff : ℚ) ≠ 0 := Nat.cast_ne_zero.mpr hc
  have hkey : (D.shellCoeff : ℚ) * (D.perEdgeFlux n₁ * (n₁ : ℚ) ^ 2) =
      (D.shellCoeff : ℚ) * (D.perEdgeFlux n₂ * (n₂ : ℚ) ^ 2) := by
    linear_combination e₁ - e₂
  exact mul_left_cancel₀ hcq hkey

end LayeredFluxData

/-- Joint satisfiability witness: a three-vertex chain with a unit
source at layer zero realizes every declared premise at depth one with
coefficient one, so the depth-bounded premise set is satisfiable with a
nonzero coefficient and the shell theorems are not vacuous. -/
def chainWitness : LayeredFluxData (Fin 3) (Fin 2) where
  left := ![0, 1]
  right := ![1, 2]
  layer := ![0, 1, 2]
  flux _ := 1
  source := ![1, 0, 0]
  shellCoeff := 1
  perEdgeFlux _ := 1
  depth := 1
  sourceSupport := by
    intro v hv
    fin_cases v <;> simp_all
  steadyGauss := by
    intro v hv
    fin_cases v
    · simp only [vertexOutflow, Fin.sum_univ_two, Matrix.cons_val_zero,
        Matrix.cons_val_one]
      norm_num [Fin.ext_iff]
    · simp only [vertexOutflow, Fin.sum_univ_two, Matrix.cons_val_zero,
        Matrix.cons_val_one]
      norm_num [Fin.ext_iff]
    · norm_num at hv
  shellCard := by
    intro n hn hd
    interval_cases n
    decide
  equidistributed := by
    intro n hn hd
    interval_cases n
    intro e he
    fin_cases e
    · exfalso
      simp [shellEdges] at he
    · norm_num [edgeOutwardFlux, layerBall, Matrix.cons_val_one,
        Matrix.head_cons, Matrix.cons_val_two, Matrix.tail_cons]

/-- The witness carries unit load. -/
theorem chainWitness_charge : chainWitness.charge = 1 := by
  simp only [LayeredFluxData.charge, chainWitness, regionLoad, layerBall,
    Finset.sum_filter, Fin.sum_univ_three, Matrix.cons_val_zero,
    Matrix.cons_val_one, Matrix.head_cons, Matrix.cons_val_two,
    Matrix.tail_cons]
  norm_num

/-! ## Bridges to the committed twelve-port receipts

The two theorems below record where the committed `RegionalContinuity`
receipts genuinely apply. They are re-exported side by side with the
generic development because the committed receipts are stated on the
fixed twelve-port, thirty-seam carrier with integer currents, and the
layered development on generic vertex and edge types cannot consume
them directly. The hypotheses do not literally compose; the bridges
state the exact points of contact. -/

/-- On the twelve-port carrier the generic outward flux coincides with
the committed `RegionalContinuity.regionFlux` after the
integer-to-rational cast. The generic definitions restrict to the
committed ones; nothing beyond a cast is involved. -/
theorem regionFlux_twelvePort (region : Finset (Fin 12))
    (current : SeamCurrent) :
    regionFlux seamLeft seamRight region (fun e => (current e : ℚ)) =
      ((OPH.RegionalContinuity.regionFlux region current : ℤ) : ℚ) := by
  unfold regionFlux edgeOutwardFlux OPH.RegionalContinuity.regionFlux
    OPH.RegionalContinuity.seamOutwardFlux
  push_cast [apply_ite (fun z : ℤ => (z : ℚ))]
  rfl

/-- A steady continuity witness (no stored-load change) forces the
committed regional outward flux to equal the regional source, through
the committed `regional_continuity` theorem alone. This is the
fixed-carrier instance of the `steadyGauss` premise: on the twelve-port
carrier it is a theorem, and on the generic layered graph it is a
declared premise. -/
theorem steady_witness_regionFlux_eq_source (region : Finset (Fin 12))
    (q source : PortLoad) (current : SeamCurrent)
    (hsteady :
      OPH.RegionalContinuity.SatisfiesContinuityUpdate q q source current) :
    OPH.RegionalContinuity.regionFlux region current =
      OPH.RegionalContinuity.regionLoad region source := by
  have h := OPH.RegionalContinuity.regional_continuity region q q source
    current hsteady
  rw [sub_self] at h
  linarith

/-! ## Axiom audit -/

#print axioms sum_vertexOutflow_eq_regionFlux
#print axioms regionFlux_univ
#print axioms mem_layerBall
#print axioms regionFlux_layerBall_eq_shell_sum
#print axioms sum_vertexOutflow_univ
#print axioms steady_total_source_eq_zero
#print axioms unbounded_shellCard_forces_zero
#print axioms LayeredFluxData.regionFlux_eq_charge
#print axioms LayeredFluxData.shell_sum_eq_card_mul
#print axioms LayeredFluxData.shell_total_eq_charge
#print axioms LayeredFluxData.perEdgeFlux_eq
#print axioms LayeredFluxData.scale_free
#print axioms chainWitness
#print axioms chainWitness_charge
#print axioms regionFlux_twelvePort
#print axioms steady_witness_regionFlux_eq_source

end OPH.LayeredDiscreteGauss
