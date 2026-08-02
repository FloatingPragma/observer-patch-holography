import Mathlib
import PrimitivePortTranslationBridge
import PortFrameGram

open scoped BigOperators goldenRatio

namespace OPH.PrimitivePortFrameQuotient

/-!
# Real scalar extension of the primitive six-axis port adapter

The simulator's conditional FZ-11 adapter sends the six positive signed-record
generators, in source-port order `(0, 1, 4, 5, 8, 9)`, to the six unoriented
icosahedral axes

`(-1, φ, 0)`, `(1, φ, 0)`, `(0, -1, φ)`, `(0, 1, φ)`,
`(φ, 0, -1)`, `(φ, 0, 1)`.

This file proves the exact real-linear statement behind that adapter.  Its
scalar extension from six real controls onto the declared Cartesian frame is
surjective, has a three-dimensional kernel, and becomes linearly equivalent
to Cartesian three-space after quotienting by that kernel.  The exact
source-to-RER relabel identifies the normalized coordinate Gram with the
finite `PortFrameGram` table.  The quotient is defined by the radical of that
source form, and the coordinate kernel is proved equal to it.  Irrationality
of the golden ratio proves that the operational integer image is injective and
dense.  The normalized source-Gram metric on the signed cumulative
port-record/control module completes to an abstract three-dimensional
Euclidean carrier, presented in coordinates by `Vec3`.

This is a source-frame quotient and completion, not a derivation of physical
spacetime.  It does not select the adapter from A1--A3, make the dense integer
image locally finite, or prove that observer readback uses the Gram topology
rather than record-count distance.  Long cumulative record controls can have
arbitrarily small nonzero Gram displacement.  Physical space, scale, rest
frame, Lorentzian structure, and particle-sector attachment remain separate.
The completion is canonical only up to port-label-preserving isometry; its
Cartesian axes are a witness, not preferred directions.  A faithful
isometric extension of the finite carrier action and overlap/refinement
gluing are not proved here.  Record order and repair cost are retained as
separate data rather than erased into carrier position.
-/

abbrev Vec3 := OPH.PrimitivePortTranslationBridge.Vec3
abbrev AxisControl := Fin 6 → ℝ

/-- Exact relabeling from the source coordinate packet to the registered
`PortFrameGram` incidence labels. -/
def sourceToRERPort : Fin 12 → Fin 12 :=
  ![8, 10, 1, 3, 7, 11, 0, 4, 6, 9, 2, 5]

/-- The six positive representatives in the source coordinate labeling. -/
def positiveSourcePort : Fin 6 → Fin 12 := ![0, 1, 4, 5, 8, 9]

/-- The same six positive representatives in the registered Gram-table
labeling. -/
def positivePort : Fin 6 → Fin 12 :=
  fun i ↦ sourceToRERPort (positiveSourcePort i)

/-- Real evaluation of the exact pair representation `a + b * sqrt(5)` used
by the finite port Gram table. -/
noncomputable def evalZ5 (x : ℤ × ℤ) : ℝ :=
  (x.1 : ℝ) + (x.2 : ℝ) * Real.sqrt 5

/-- Six-axis Gram matrix read directly from the finite twelve-port source
table.  The source stores `5 G`, hence the division by five. -/
noncomputable def sourceGramEntry (i j : Fin 6) : ℝ :=
  evalZ5 (OPH.PortFrameGram.g5 (positivePort i) (positivePort j)) / 5

/-- The exact six positive source axes used by the pinned simulator adapter.
The common normalization `1 / sqrt (2 + φ)` is omitted because it is a
nonzero common scale and does not change the kernel or the quotient dimension.
-/
noncomputable def rawAxis : Fin 6 → Vec3 :=
  ![
    ![-1, φ, 0],
    ![1, φ, 0],
    ![0, -1, φ],
    ![0, 1, φ],
    ![φ, 0, -1],
    ![φ, 0, 1]
  ]

set_option maxHeartbeats 4000000 in
/-- The coordinate realization is a witness for the source Gram table rather
than an independent metric choice.  Its common squared radius is `φ + 2`.
Dividing the raw Cartesian dot product by that radius reproduces every entry
of the exact finite port Gram table. -/
theorem rawAxis_dot_eq_sourceGram (i j : Fin 6) :
    OPH.PrimitivePortTranslationBridge.dot (rawAxis i) (rawAxis j) =
      (φ + 2) * sourceGramEntry i j := by
  have hsqrt : Real.sqrt 5 ^ 2 = 5 := by norm_num
  fin_cases i <;> fin_cases j <;>
    simp [OPH.PrimitivePortTranslationBridge.dot, rawAxis, sourceGramEntry,
      evalZ5, positivePort, positiveSourcePort, sourceToRERPort,
      OPH.PortFrameGram.g5, OPH.PortFrameGram.adj,
      OPH.PortFrameGram.neighbors, OPH.PortFrameGram.antipode,
      Real.goldenRatio, Fin.sum_univ_succ] <;>
    nlinarith

/-- Real scalar extension of the signed-record-to-frame adapter. -/
noncomputable def frameMap : AxisControl →ₗ[ℝ] Vec3 where
  toFun c := fun d ↦ ∑ i : Fin 6, c i * rawAxis i d
  map_add' c e := by
    ext d
    simp only [Pi.add_apply, add_mul, Finset.sum_add_distrib]
  map_smul' a c := by
    ext d
    simp only [Pi.smul_apply, smul_eq_mul, RingHom.id_apply]
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro i _
    ring

/-- A concrete right inverse.  Only four of the six axes are needed: the
differences of the first pair and the last pair give the Cartesian `x` and
`z` directions, while the sum of the first pair gives `y`. -/
noncomputable def framePreimage (v : Vec3) : AxisControl :=
  ![
    -(v 0) / 2 + (v 1) / (2 * φ),
    (v 0) / 2 + (v 1) / (2 * φ),
    0,
    0,
    -(v 2) / 2,
    (v 2) / 2
  ]

theorem frameMap_preimage (v : Vec3) : frameMap (framePreimage v) = v := by
  ext d
  fin_cases d <;>
    simp [frameMap, framePreimage, rawAxis, Fin.sum_univ_succ] <;>
    field_simp [Real.goldenRatio_ne_zero] <;> ring

/-- The real scalar extension reaches every vector in the declared frame. -/
theorem frameMap_surjective : Function.Surjective frameMap := by
  intro v
  exact ⟨framePreimage v, frameMap_preimage v⟩

/-! ## Operational integer controls

The constructive source emits integer signed-record loads.  The next map is
their exact inclusion in the real scalar extension.  Injectivity is a useful
boundary fact: the real kernel consists of genuinely irrational control
combinations and identifies no two integer records.
-/

/-- Coordinatewise cast of the operational integer control module. -/
def castIntegerControl (n : Fin 6 → ℤ) : AxisControl :=
  fun i ↦ (n i : ℝ)

/-- The integer signed-record image in the declared source frame. -/
noncomputable def integerFrame (n : Fin 6 → ℤ) : Vec3 :=
  frameMap (castIntegerControl n)

private theorem int_phi_split {a b : ℤ}
    (h : (a : ℝ) + (b : ℝ) * φ = 0) : a = 0 ∧ b = 0 := by
  have hb : b = 0 := by
    by_contra hb0
    have hirr : Irrational ((b : ℝ) * φ) :=
      Real.goldenRatio_irrational.intCast_mul hb0
    have heq : (b : ℝ) * φ = ((-a : ℤ) : ℝ) := by
      push_cast
      linarith
    rw [heq] at hirr
    exact Int.not_irrational (-a) hirr
  subst b
  simp only [Int.cast_zero, zero_mul, add_zero] at h
  exact ⟨by exact_mod_cast h, rfl⟩

/-- The six operational integer record coordinates remain distinct in the
three-real-coordinate adapter.  This is compatible with a three-dimensional
real quotient because the real kernel contains irrational, not integer,
control directions. -/
theorem integerFrame_injective : Function.Injective integerFrame := by
  intro n m hnm
  let z : Fin 6 → ℤ := fun i ↦ n i - m i
  have hz : integerFrame z = 0 := by
    funext d
    have hd := congrFun hnm d
    simp [integerFrame, frameMap, castIntegerControl, z,
      Fin.sum_univ_succ] at hd ⊢
    linarith
  have hx := congrFun hz (0 : Fin 3)
  have hy := congrFun hz (1 : Fin 3)
  have hzcoord := congrFun hz (2 : Fin 3)
  simp [integerFrame, frameMap, castIntegerControl, rawAxis,
    Fin.sum_univ_succ] at hx hy hzcoord
  have sx := int_phi_split
    (a := -z 0 + z 1) (b := z 4 + z 5) (by
      push_cast
      linarith [hx])
  have sy := int_phi_split
    (a := -z 2 + z 3) (b := z 0 + z 1) (by
      push_cast
      linarith [hy])
  have sz := int_phi_split
    (a := -z 4 + z 5) (b := z 2 + z 3) (by
      push_cast
      linarith [hzcoord])
  funext i
  fin_cases i <;> simp [z] at sx sy sz ⊢ <;> omega

/-! ## Dense integer image and completion

The density statement is proved rather than imported from the simulator
receipt.  Its arithmetic input is the irrationality of the golden ratio.
-/

/-- One real coordinate generated by an integer pair in `ℤ + φ ℤ`. -/
noncomputable def goldenLine (p : ℤ × ℤ) : ℝ :=
  (p.1 : ℝ) + (p.2 : ℝ) * φ

/-- `ℤ + φ ℤ` is dense in the real line. -/
theorem goldenLine_denseRange : DenseRange goldenLine := by
  have hd : Dense
      (AddSubgroup.closure ({φ, (1 : ℝ)} : Set ℝ) : Set ℝ) := by
    rw [dense_addSubgroupClosure_pair_iff]
    simpa using Real.goldenRatio_irrational
  rw [DenseRange]
  apply hd.mono
  intro x hx
  change x ∈ AddSubgroup.closure ({φ, (1 : ℝ)} : Set ℝ) at hx
  rw [AddSubgroup.mem_closure_pair] at hx
  rcases hx with ⟨m, n, rfl⟩
  refine ⟨(n, m), ?_⟩
  simp [goldenLine, add_comm]

/-- Integer controls that realize twice three independently chosen golden-line
coordinates.  The factor two is the finite-index parity residue of this
particularly convenient section and has no effect on density. -/
def denseIntegerControl (p : Fin 3 → ℤ × ℤ) : Fin 6 → ℤ :=
  ![
    (p 1).2 - (p 0).1,
    (p 1).2 + (p 0).1,
    (p 2).2 - (p 1).1,
    (p 2).2 + (p 1).1,
    (p 0).2 - (p 2).1,
    (p 0).2 + (p 2).1
  ]

/-- Exact section identity for the dense three-coordinate subset. -/
theorem integerFrame_denseIntegerControl (p : Fin 3 → ℤ × ℤ) :
    integerFrame (denseIntegerControl p) =
      fun d ↦ 2 * goldenLine (p d) := by
  ext d
  fin_cases d <;>
    simp [integerFrame, denseIntegerControl, goldenLine, frameMap,
      castIntegerControl, rawAxis, Fin.sum_univ_succ] <;>
    ring

/-- The exact operational integer image is dense in the declared Cartesian
frame.  This does not make the image locally finite; in fact density rules out
that interpretation in the inherited Euclidean topology. -/
theorem integerFrame_denseRange : DenseRange integerFrame := by
  have hscaled : DenseRange (fun p : ℤ × ℤ ↦ 2 * goldenLine p) := by
    have htwo : DenseRange (fun x : ℝ ↦ 2 * x) :=
      (show Function.Surjective (fun x : ℝ ↦ 2 * x) by
        intro y
        exact ⟨y / 2, by ring⟩).denseRange
    exact htwo.comp goldenLine_denseRange (by fun_prop)
  have hproduct : DenseRange
      (fun p : Fin 3 → ℤ × ℤ ↦ fun d ↦ 2 * goldenLine (p d)) := by
    simpa only [Pi.map_apply] using
      (DenseRange.piMap (fun _ : Fin 3 ↦ hscaled))
  apply hproduct.mono
  rintro _ ⟨p, rfl⟩
  exact ⟨denseIntegerControl p, integerFrame_denseIntegerControl p⟩

/-- Wrapper carrying the operational integer controls with the metric induced
by their exact frame image.  A wrapper is used so this metric cannot be
confused with the ordinary discrete/product metric on integer tuples. -/
@[ext]
structure IntegerFramePoint where
  control : Fin 6 → ℤ

/-- Exact frame map on the metric wrapper. -/
noncomputable def pointFrame (p : IntegerFramePoint) : Vec3 :=
  integerFrame p.control

theorem pointFrame_injective : Function.Injective pointFrame := by
  intro p q h
  ext i
  exact congrFun (integerFrame_injective h) i

/-- Euclidean realization of the exact three-coordinate frame. -/
abbrev EuclideanVec3 := EuclideanSpace ℝ (Fin 3)

/-- The canonical uniform equivalence between Euclidean `L²` coordinates and
the finite product presentation used by the translation calculations. -/
noncomputable def euclideanEquivVec3 : EuclideanVec3 ≃ᵤ Vec3 :=
  PiLp.uniformEquiv 2 (fun _ : Fin 3 ↦ ℝ)

/-- Common raw-axis radius used only to normalize the coordinate witness to
the unit-diagonal finite source Gram. -/
noncomputable def rawRadius : ℝ := Real.sqrt (φ + 2)

theorem rawRadius_pos : 0 < rawRadius := by
  unfold rawRadius
  exact Real.sqrt_pos.2 (by linarith [Real.goldenRatio_pos])

theorem rawRadius_ne_zero : rawRadius ≠ 0 := ne_of_gt rawRadius_pos

@[simp]
theorem euclideanEquivVec3_symm_apply (v : Vec3) (d : Fin 3) :
    (euclideanEquivVec3.symm v) d = v d := by
  rfl

/-- Exact integer frame map into Euclidean three-space. -/
noncomputable def pointEuclideanFrame (p : IntegerFramePoint) : EuclideanVec3 :=
  rawRadius⁻¹ • euclideanEquivVec3.symm (pointFrame p)

theorem pointEuclideanFrame_injective : Function.Injective pointEuclideanFrame := by
  intro p q h
  have hraw : euclideanEquivVec3.symm (pointFrame p) =
      euclideanEquivVec3.symm (pointFrame q) :=
    smul_right_injective EuclideanVec3 (inv_ne_zero rawRadius_ne_zero) h
  apply pointFrame_injective
  exact euclideanEquivVec3.symm.injective hraw

/-- Euclidean metric pulled back along the exact frame witness.  The theorem
`gram_eq_radius_mul_sourceGram` below proves that its square is the finite
source Gram metric times one common positive radius factor. -/
noncomputable instance integerFramePointMetric : MetricSpace IntegerFramePoint :=
  MetricSpace.induced pointEuclideanFrame pointEuclideanFrame_injective inferInstance

/-- The wrapped operational map is an isometry by construction of the pulled
back Euclidean metric. -/
theorem pointEuclideanFrame_isometry : Isometry pointEuclideanFrame := by
  exact MetricSpace.isometry_induced pointEuclideanFrame
    pointEuclideanFrame_injective

/-- The wrapped operational map retains the proved dense range. -/
theorem pointFrame_denseRange : DenseRange pointFrame := by
  apply integerFrame_denseRange.mono
  rintro _ ⟨n, rfl⟩
  exact ⟨⟨n⟩, rfl⟩

/-- The Euclidean realization retains the proved dense range. -/
theorem pointEuclideanFrame_denseRange : DenseRange pointEuclideanFrame := by
  have hbase : DenseRange
      (fun p : IntegerFramePoint ↦ euclideanEquivVec3.symm (pointFrame p)) :=
    euclideanEquivVec3.symm.surjective.denseRange.comp
      pointFrame_denseRange euclideanEquivVec3.symm.continuous
  have hscale : Function.Surjective
      (fun x : EuclideanVec3 ↦ rawRadius⁻¹ • x) := by
    intro y
    refine ⟨rawRadius • y, ?_⟩
    change rawRadius⁻¹ • (rawRadius • y) = y
    rw [← mul_smul]
    simp [rawRadius_ne_zero]
  simpa only [pointEuclideanFrame, Function.comp_apply] using
    hscale.denseRange.comp hbase (by fun_prop)

/-- Euclidean three-space is an abstract completion of the operational integer
adapter equipped with its source-Gram metric. -/
noncomputable def integerFrameCompletion :
    AbstractCompletion IntegerFramePoint where
  space := EuclideanVec3
  coe := pointEuclideanFrame
  uniformStruct := inferInstance
  complete := inferInstance
  separation := inferInstance
  isUniformInducing := pointEuclideanFrame_isometry.isUniformInducing
  dense := pointEuclideanFrame_denseRange

/-- Canonical uniform equivalence between the completion of the exact integer
adapter and Euclidean three-space. -/
noncomputable def integerCompletionEquivEuclidean3 :
    UniformSpace.Completion IntegerFramePoint ≃ᵤ EuclideanVec3 :=
  AbstractCompletion.compareEquiv UniformSpace.Completion.cPkg
    integerFrameCompletion

/-- Canonical uniform equivalence between the metric completion of the exact
integer adapter and the product-coordinate presentation.  Its
operational-to-physical interpretation, scale, rest frame, Lorentzian
structure, and spacetime attachment remain separate.
-/
noncomputable def integerCompletionEquivVec3 :
    UniformSpace.Completion IntegerFramePoint ≃ᵤ Vec3 :=
  integerCompletionEquivEuclidean3.trans euclideanEquivVec3

/-- The Gram form pulled back from the declared Cartesian dot product. -/
noncomputable def gram (c e : AxisControl) : ℝ :=
  OPH.PrimitivePortTranslationBridge.dot (frameMap c) (frameMap e)

/-- Bilinear form constructed directly from the finite source Gram table. -/
noncomputable def sourceGram (c e : AxisControl) : ℝ :=
  ∑ i : Fin 6, ∑ j : Fin 6, c i * e j * sourceGramEntry i j

/-- Expanding the coordinate witness expresses the pulled-back form as the
pairwise Gram sum over the six atomic records. -/
theorem gram_eq_axis_sum (c e : AxisControl) :
    gram c e = ∑ i : Fin 6, ∑ j : Fin 6,
      c i * e j * OPH.PrimitivePortTranslationBridge.dot (rawAxis i) (rawAxis j) := by
  unfold gram OPH.PrimitivePortTranslationBridge.dot frameMap
  simp only [LinearMap.coe_mk, AddHom.coe_mk]
  simp_rw [Finset.sum_mul, Finset.mul_sum]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro i _
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro j _
  apply Finset.sum_congr rfl
  intro d _
  ring

/-- The finite source table determines the entire pulled-back metric up to its
single common positive radius factor.  Cartesian coordinates only witness
this identity; they do not supply an additional Gram choice. -/
theorem gram_eq_radius_mul_sourceGram (c e : AxisControl) :
    gram c e = (φ + 2) * sourceGram c e := by
  rw [gram_eq_axis_sum]
  unfold sourceGram
  simp_rw [rawAxis_dot_eq_sourceGram]
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i _
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro j _
  ring

/-- Real control difference associated with two cumulative integer record
controls. -/
def pointControlDifference (p q : IntegerFramePoint) : AxisControl :=
  fun i ↦ ((p.control i - q.control i : ℤ) : ℝ)

theorem frameMap_pointControlDifference (p q : IntegerFramePoint) :
    frameMap (pointControlDifference p q) = pointFrame p - pointFrame q := by
  ext d
  simp [frameMap, pointControlDifference, pointFrame, integerFrame,
    castIntegerControl]
  rw [← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro i _
  ring

/-- The square of the normalized operational metric is the raw coordinate
Gram quadratic form divided by the common squared radius. -/
theorem point_dist_sq_eq_scaled_gram (p q : IntegerFramePoint) :
    dist p q ^ 2 = rawRadius⁻¹ ^ 2 *
      gram (pointControlDifference p q) (pointControlDifference p q) := by
  change dist (pointEuclideanFrame p) (pointEuclideanFrame q) ^ 2 = _
  rw [EuclideanSpace.dist_sq_eq]
  unfold gram
  rw [frameMap_pointControlDifference]
  unfold OPH.PrimitivePortTranslationBridge.dot
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro d _
  simp [pointEuclideanFrame, Real.dist_eq, euclideanEquivVec3_symm_apply]
  ring

/-- Consequently the operational metric is exactly the unit-diagonal finite
source Gram metric.  This rules out any hidden coordinate metric choice in
the completion statement. -/
theorem point_dist_sq_eq_sourceGram (p q : IntegerFramePoint) :
    dist p q ^ 2 =
      sourceGram (pointControlDifference p q) (pointControlDifference p q) := by
  rw [point_dist_sq_eq_scaled_gram, gram_eq_radius_mul_sourceGram]
  have hradius : rawRadius ^ 2 = φ + 2 := by
    unfold rawRadius
    exact Real.sq_sqrt (by linarith [Real.goldenRatio_pos])
  rw [← hradius]
  field_simp [rawRadius_ne_zero]

/-- The quadratic Gram nullspace equals the kernel of the frame map. -/
theorem gram_self_eq_zero_iff (c : AxisControl) :
    gram c c = 0 ↔ c ∈ LinearMap.ker frameMap := by
  change (∑ d : Fin 3, frameMap c d * frameMap c d) = 0 ↔ frameMap c = 0
  simpa [dotProduct] using
    (dotProduct_self_eq_zero (v := frameMap c))

/-- The bilinear radical of the pulled-back Gram form is the same kernel. -/
theorem gram_radical_iff (c : AxisControl) :
    (∀ e : AxisControl, gram c e = 0) ↔ c ∈ LinearMap.ker frameMap := by
  constructor
  · intro h
    have hself := h (framePreimage (frameMap c))
    have hdot : OPH.PrimitivePortTranslationBridge.dot
        (frameMap c) (frameMap c) = 0 := by
      simpa [gram, frameMap_preimage] using hself
    change frameMap c = 0
    apply (dotProduct_self_eq_zero (v := frameMap c)).mp
    simpa [dotProduct, OPH.PrimitivePortTranslationBridge.dot] using hdot
  · intro hc e
    change frameMap c = 0 at hc
    simp [gram, hc, OPH.PrimitivePortTranslationBridge.dot]

/-- Source-Gram bilinear form as a linear map into the algebraic dual. -/
theorem sourceGram_add_right (c e f : AxisControl) :
    sourceGram c (e + f) = sourceGram c e + sourceGram c f := by
  unfold sourceGram
  simp_rw [Pi.add_apply, mul_add, add_mul, Finset.sum_add_distrib]

theorem sourceGram_smul_right (a : ℝ) (c e : AxisControl) :
    sourceGram c (a • e) = a * sourceGram c e := by
  unfold sourceGram
  simp only [Pi.smul_apply, smul_eq_mul]
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i _
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro j _
  ring

theorem sourceGram_add_left (c d e : AxisControl) :
    sourceGram (c + d) e = sourceGram c e + sourceGram d e := by
  unfold sourceGram
  simp_rw [Pi.add_apply, add_mul, Finset.sum_add_distrib]

theorem sourceGram_smul_left (a : ℝ) (c e : AxisControl) :
    sourceGram (a • c) e = a * sourceGram c e := by
  unfold sourceGram
  simp only [Pi.smul_apply, smul_eq_mul]
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro i _
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro j _
  ring

noncomputable def sourceGramLinear :
    AxisControl →ₗ[ℝ] (AxisControl →ₗ[ℝ] ℝ) where
  toFun c :=
    { toFun := fun e ↦ sourceGram c e
      map_add' := sourceGram_add_right c
      map_smul' := by
        intro a e
        simpa only [RingHom.id_apply, smul_eq_mul] using
          sourceGram_smul_right a c e }
  map_add' := by
    intro c d
    apply LinearMap.ext
    intro e
    change sourceGram (c + d) e = sourceGram c e + sourceGram d e
    exact sourceGram_add_left c d e
  map_smul' := by
    intro a c
    apply LinearMap.ext
    intro e
    change sourceGram (a • c) e = a • sourceGram c e
    simpa only [smul_eq_mul] using sourceGram_smul_left a c e

/-- Radical defined without coordinates, directly as the kernel of the finite
source Gram map. -/
noncomputable def sourceGramRadical : Submodule ℝ AxisControl :=
  LinearMap.ker sourceGramLinear

theorem mem_sourceGramRadical_iff (c : AxisControl) :
    c ∈ sourceGramRadical ↔ ∀ e : AxisControl, sourceGram c e = 0 := by
  rw [sourceGramRadical, LinearMap.mem_ker]
  constructor
  · intro h e
    have he := LinearMap.congr_fun h e
    simpa [sourceGramLinear] using he
  · intro h
    apply LinearMap.ext
    intro e
    simpa [sourceGramLinear] using h e

/-- The coordinate witness has exactly the radical intrinsic to the finite
source Gram table. -/
theorem sourceGramRadical_eq_frameKernel :
    sourceGramRadical = LinearMap.ker frameMap := by
  ext c
  rw [mem_sourceGramRadical_iff, ← gram_radical_iff]
  constructor
  · intro h e
    rw [gram_eq_radius_mul_sourceGram, h e, mul_zero]
  · intro h e
    have he := h e
    rw [gram_eq_radius_mul_sourceGram] at he
    exact (mul_eq_zero.mp he).resolve_left
      (by linarith [Real.goldenRatio_pos])

abbrev FrameQuotient := AxisControl ⧸ sourceGramRadical

/-- First-isomorphism equivalence from the Gram quotient to the declared
three-dimensional frame. -/
noncomputable def sourceQuotientFrameMap : FrameQuotient →ₗ[ℝ] Vec3 :=
  Submodule.liftQ sourceGramRadical frameMap (by
    intro c hc
    rw [sourceGramRadical_eq_frameKernel] at hc
    exact hc)

@[simp]
theorem sourceQuotientFrameMap_mk (c : AxisControl) :
    sourceQuotientFrameMap (Submodule.Quotient.mk c) = frameMap c := by
  rfl

theorem sourceQuotientFrameMap_injective :
    Function.Injective sourceQuotientFrameMap := by
  apply LinearMap.ker_eq_bot.mp
  ext x
  obtain ⟨c, rfl⟩ := Submodule.Quotient.mk_surjective sourceGramRadical x
  simp [sourceGramRadical_eq_frameKernel]

theorem sourceQuotientFrameMap_surjective :
    Function.Surjective sourceQuotientFrameMap := by
  intro v
  exact ⟨Submodule.Quotient.mk (framePreimage v), by
    simp [frameMap_preimage]⟩

noncomputable def quotientEquivVec3 : FrameQuotient ≃ₗ[ℝ] Vec3 :=
  LinearEquiv.ofBijective sourceQuotientFrameMap
    ⟨sourceQuotientFrameMap_injective, sourceQuotientFrameMap_surjective⟩

@[simp]
theorem quotientEquivVec3_mk (c : AxisControl) :
    quotientEquivVec3 (Submodule.Quotient.mk c) = frameMap c := by
  exact sourceQuotientFrameMap_mk c

/-- The descended inner-product candidate.  Its metric content is precisely
the Cartesian form declared on `Vec3`; this theorem does not select that form
as the physical metric. -/
noncomputable def quotientGram (q r : FrameQuotient) : ℝ :=
  OPH.PrimitivePortTranslationBridge.dot (quotientEquivVec3 q) (quotientEquivVec3 r)

/-- The quotient equivalence preserves the transported Gram form exactly. -/
theorem quotientEquivVec3_preserves_gram (q r : FrameQuotient) :
    quotientGram q r =
      OPH.PrimitivePortTranslationBridge.dot
        (quotientEquivVec3 q) (quotientEquivVec3 r) := by
  rfl

/-- The quotient has real dimension three. -/
theorem frameQuotient_finrank : Module.finrank ℝ FrameQuotient = 3 := by
  rw [quotientEquivVec3.finrank_eq]
  simp

/-- The discarded Gram-null directions have real dimension three. -/
theorem frameKernel_finrank : Module.finrank ℝ (LinearMap.ker frameMap) = 3 := by
  have hrank :
      Module.finrank ℝ FrameQuotient +
        Module.finrank ℝ sourceGramRadical = 6 := by
    simpa [FrameQuotient, AxisControl] using
      sourceGramRadical.finrank_quotient_add_finrank
  rw [frameQuotient_finrank] at hrank
  have hradical : Module.finrank ℝ sourceGramRadical = 3 := by omega
  rwa [sourceGramRadical_eq_frameKernel] at hradical

/-- The descended Gram form is positive definite. -/
theorem quotientGram_self_eq_zero_iff (q : FrameQuotient) :
    quotientGram q q = 0 ↔ q = 0 := by
  change
    (∑ d : Fin 3, quotientEquivVec3 q d * quotientEquivVec3 q d) = 0 ↔ q = 0
  rw [show q = 0 ↔ quotientEquivVec3 q = 0 by
    exact quotientEquivVec3.map_eq_zero_iff.symm]
  simpa [dotProduct] using
    (dotProduct_self_eq_zero (v := quotientEquivVec3 q))

/-! ## Positive-definite Gram cores

`Vec3` is represented by a finite function type, whose ambient Mathlib norm is
not used to smuggle in a Euclidean claim.  The following two explicit cores
record the Cartesian dot product and its pullback to the quotient.  A
positive-definite `InnerProductSpace.Core` is precisely the algebraic datum
from which Mathlib constructs the associated normed inner-product structure.
-/

/-- Cartesian positive-definite inner-product core on the declared frame. -/
@[reducible] noncomputable def cartesianInnerCore : InnerProductSpace.Core ℝ Vec3 where
  inner := OPH.PrimitivePortTranslationBridge.dot
  conj_inner_symm x y := by
    simp [OPH.PrimitivePortTranslationBridge.dot, mul_comm]
  re_inner_nonneg x := by
    change 0 ≤ ∑ d : Fin 3, x d * x d
    exact Finset.sum_nonneg fun d _ ↦ mul_self_nonneg (x d)
  add_left x y z := by
    change
      (∑ d : Fin 3, (x d + y d) * z d) =
        (∑ d : Fin 3, x d * z d) + ∑ d : Fin 3, y d * z d
    rw [← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro d _
    ring
  smul_left x y a := by
    change
      (∑ d : Fin 3, (a * x d) * y d) =
        a * ∑ d : Fin 3, x d * y d
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro d _
    ring
  definite x hx := by
    apply (dotProduct_self_eq_zero (v := x)).mp
    simpa [dotProduct, OPH.PrimitivePortTranslationBridge.dot] using hx

/-- Positive-definite inner-product core on the Gram quotient. -/
@[reducible] noncomputable def quotientInnerCore : InnerProductSpace.Core ℝ FrameQuotient where
  inner := quotientGram
  conj_inner_symm q r := by
    simp [quotientGram, OPH.PrimitivePortTranslationBridge.dot, mul_comm]
  re_inner_nonneg q := by
    change 0 ≤ ∑ d : Fin 3,
      quotientEquivVec3 q d * quotientEquivVec3 q d
    exact Finset.sum_nonneg fun d _ ↦ mul_self_nonneg (quotientEquivVec3 q d)
  add_left q r s := by
    change
      (∑ d : Fin 3,
        quotientEquivVec3 (q + r) d * quotientEquivVec3 s d) =
      (∑ d : Fin 3, quotientEquivVec3 q d * quotientEquivVec3 s d) +
        ∑ d : Fin 3, quotientEquivVec3 r d * quotientEquivVec3 s d
    simp only [map_add, Pi.add_apply, add_mul]
    exact Finset.sum_add_distrib
  smul_left q r a := by
    change
      (∑ d : Fin 3,
        quotientEquivVec3 (a • q) d * quotientEquivVec3 r d) =
      a * ∑ d : Fin 3, quotientEquivVec3 q d * quotientEquivVec3 r d
    simp only [map_smul, Pi.smul_apply, smul_eq_mul]
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro d _
    ring
  definite q h := (quotientGram_self_eq_zero_iff q).mp h

/-- The quotient equivalence is an exact isometry of the two explicit
positive-definite Gram cores.  The statement is relative to the declared
Cartesian frame; it is not a physical metric-selection theorem. -/
theorem quotientEquivVec3_core_isometry (q r : FrameQuotient) :
    cartesianInnerCore.inner (quotientEquivVec3 q) (quotientEquivVec3 r) =
      quotientInnerCore.inner q r := by
  rfl

#print axioms frameMap_surjective
#print axioms rawAxis_dot_eq_sourceGram
#print axioms integerFrame_injective
#print axioms goldenLine_denseRange
#print axioms integerFrame_denseRange
#print axioms pointFrame_denseRange
#print axioms pointEuclideanFrame_denseRange
#print axioms point_dist_sq_eq_sourceGram
#print axioms gram_self_eq_zero_iff
#print axioms gram_radical_iff
#print axioms sourceGramRadical_eq_frameKernel
#print axioms frameQuotient_finrank
#print axioms frameKernel_finrank
#print axioms quotientGram_self_eq_zero_iff
#print axioms quotientEquivVec3_core_isometry

end OPH.PrimitivePortFrameQuotient
