import RepairWordCarrierReadout

open scoped BigOperators Topology

namespace OPH.SeamCurrentCarrierQuotient

open OPH.PortFrameGram
open OPH.PrimitivePortFrameQuotient
open OPH.RepairWordCarrierReadout
open ObserverPatchHolography.A2EndpointCommutator

/-!
# The thirty-seam current quotient

The finite source has thirty undirected seams.  Giving every seam its
smaller endpoint first fixes a bookkeeping orientation, and an integer seam
current has the usual signed incidence boundary on the twelve ports.  This
file computes the antipodal-odd readback of that boundary exactly.

Its image is the even-sum root lattice

`D6 = {z : Fin 6 → ℤ | Even (∑ i, z i)}`.

The proof includes an explicit six-seam section.  The residual cokernel is
therefore precisely one parity bit.  The `D6` controls still have dense image
in the response-selected three-dimensional Gram carrier, and that carrier is
an abstract metric completion of the even-sum records.

These are incidence, quotient, and metric-completion statements.  A seam is
not identified with a primitive six-axis translation.  The theorems do not
select physical position, a physical seam current, scale, time, or global
observer gluing.
-/

/-! ## Exact seam table and orientation convention -/

/-- The smaller endpoint of each of the thirty undirected icosahedral seams.
The index is bookkeeping only. -/
def seamLeft : Fin 30 → Fin 12 :=
  ![0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3,
    4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10]

/-- The larger endpoint of each seam, in the same exact table. -/
def seamRight : Fin 30 → Fin 12 :=
  ![1, 2, 3, 4, 6, 2, 3, 5, 7, 4, 5, 8, 6, 7, 9,
    6, 8, 10, 7, 8, 11, 9, 10, 9, 11, 10, 11, 10, 11, 11]

/-- Every table row is a canonically oriented incidence seam. -/
theorem seam_table_sound :
    ∀ e : Fin 30,
      seamLeft e < seamRight e ∧ adj (seamLeft e) (seamRight e) = true := by
  decide

/-- No two seam indices encode the same endpoint pair. -/
theorem seam_table_injective :
    Function.Injective (fun e : Fin 30 ↦ (seamLeft e, seamRight e)) := by
  decide

set_option maxRecDepth 8192 in
/-- The table exhausts the undirected incidence relation exactly. -/
theorem seam_table_complete :
    ∀ u v : Fin 12,
      (u < v ∧ adj u v = true) ↔
        ∃ e : Fin 30, seamLeft e = u ∧ seamRight e = v := by
  decide +kernel

/-- Completeness and the injective table give a unique row for every
canonically oriented seam. -/
theorem seam_table_unique
    {u v : Fin 12} (huv : u < v) (hadj : adj u v = true) :
    ∃! e : Fin 30, seamLeft e = u ∧ seamRight e = v := by
  obtain ⟨e, he⟩ := (seam_table_complete u v).1 ⟨huv, hadj⟩
  refine ⟨e, he, ?_⟩
  intro b hb
  apply seam_table_injective
  simp [he.1, he.2, hb.1, hb.2]

/-- Signed incidence boundary of an oriented port pair: current enters `v`
and leaves `u`. -/
def directedBoundary (u v : Fin 12) (n : ℤ) : PortLoad :=
  fun p ↦ (if p = v then n else 0) - (if p = u then n else 0)

/-- Reversing a seam and reversing its current leaves the signed boundary
unchanged.  Thus the bookkeeping orientation has no physical content. -/
theorem directedBoundary_reverse_neg (u v : Fin 12) (n : ℤ) :
    directedBoundary v u (-n) = directedBoundary u v n := by
  funext p
  simp only [directedBoundary]
  split_ifs <;> ring

/-- Integer current on the exact thirty-seam alphabet. -/
abbrev SeamCurrent := Fin 30 → ℤ

/-- Boundary load induced on the twelve ports by a complete seam current. -/
def seamBoundary (c : SeamCurrent) : PortLoad :=
  fun p ↦ ∑ e : Fin 30,
    directedBoundary (seamLeft e) (seamRight e) (c e) p

/-- Current supported on one canonical seam. -/
def seamAtom (e : Fin 30) (n : ℤ) : SeamCurrent :=
  fun k ↦ if k = e then n else 0

/-- One canonical seam atom has exactly the orientation-independent directed
boundary declared above. -/
theorem seamBoundary_atom (e : Fin 30) (n : ℤ) :
    seamBoundary (seamAtom e n) =
      directedBoundary (seamLeft e) (seamRight e) n := by
  funext p
  unfold seamBoundary seamAtom
  rw [Finset.sum_eq_single e]
  · simp
  · intro b _ hbe
    simp [hbe, directedBoundary]
  · simp

/-- One oriented incidence boundary conserves load. -/
theorem directedBoundary_total_zero (u v : Fin 12) (n : ℤ) :
    ∑ p : Fin 12, directedBoundary u v n p = 0 := by
  simp [directedBoundary, Finset.sum_sub_distrib]

/-- Every seam boundary conserves the total twelve-port load. -/
theorem seamBoundary_total_zero (c : SeamCurrent) :
    ∑ p : Fin 12, seamBoundary c p = 0 := by
  unfold seamBoundary
  rw [Finset.sum_comm]
  apply Finset.sum_eq_zero
  intro e _
  exact directedBoundary_total_zero (seamLeft e) (seamRight e) (c e)

/-! ## Boundary rank and cycle rank

The integral image statement below is proved later after the antipodal
readback.  For the ordinary incidence boundary itself, rank is most cleanly
stated after scalar extension to `ℚ`.  The extension changes no finite table
entry and permits exact use of rank-nullity.
-/

abbrev RationalSeamCurrent := Fin 30 → ℚ
abbrev RationalPortLoad := Fin 12 → ℚ

/-- Rational scalar extension of one signed incidence boundary. -/
def rationalDirectedBoundary (u v : Fin 12) (n : ℚ) : RationalPortLoad :=
  fun p ↦ (if p = v then n else 0) - (if p = u then n else 0)

theorem rationalDirectedBoundary_add (u v : Fin 12) (a b : ℚ) :
    rationalDirectedBoundary u v (a + b) =
      rationalDirectedBoundary u v a + rationalDirectedBoundary u v b := by
  funext p
  simp only [rationalDirectedBoundary, Pi.add_apply]
  split_ifs <;> ring

theorem rationalDirectedBoundary_smul (u v : Fin 12) (a n : ℚ) :
    rationalDirectedBoundary u v (a * n) =
      a • rationalDirectedBoundary u v n := by
  funext p
  simp only [rationalDirectedBoundary, Pi.smul_apply, smul_eq_mul]
  split_ifs <;> ring

/-- Rational scalar extension of the exact thirty-seam boundary map. -/
def rationalSeamBoundary : RationalSeamCurrent →ₗ[ℚ] RationalPortLoad where
  toFun c := fun p ↦ ∑ e : Fin 30,
    rationalDirectedBoundary (seamLeft e) (seamRight e) (c e) p
  map_add' c d := by
    funext p
    change (∑ e : Fin 30, rationalDirectedBoundary (seamLeft e)
      (seamRight e) (c e + d e) p) = _
    simp_rw [rationalDirectedBoundary_add]
    exact Finset.sum_add_distrib
  map_smul' a c := by
    funext p
    change (∑ e : Fin 30, rationalDirectedBoundary (seamLeft e)
      (seamRight e) (a * c e) p) = _
    simp_rw [rationalDirectedBoundary_smul]
    simp only [Pi.smul_apply, RingHom.id_apply, smul_eq_mul]
    rw [Finset.mul_sum]

/-- Total port load as a rational linear functional. -/
def rationalPortTotal : RationalPortLoad →ₗ[ℚ] ℚ where
  toFun x := ∑ p : Fin 12, x p
  map_add' x y := by simp [Finset.sum_add_distrib]
  map_smul' a x := by simp [Finset.mul_sum]

/-- Explicit spanning-tree section of the boundary on a zero-total load.
The eleven used seams form a rooted tree at port zero. -/
def rationalBoundarySection (x : RationalPortLoad) : RationalSeamCurrent :=
  ![
    x 1 + x 5 + x 7 + x 11,
    x 2 + x 8,
    x 3 + x 9,
    x 4 + x 10,
    x 6,
    0, 0,
    x 5 + x 11,
    x 7,
    0, 0,
    x 8,
    0, 0,
    x 9,
    0, 0,
    x 10,
    0, 0,
    x 11,
    0, 0, 0, 0, 0, 0, 0, 0, 0
  ]

set_option maxHeartbeats 2000000 in
/-- The tree section reconstructs every rational load of total zero. -/
theorem rationalSeamBoundary_section
    (x : RationalPortLoad) (hx : rationalPortTotal x = 0) :
    rationalSeamBoundary (rationalBoundarySection x) = x := by
  funext p
  fin_cases p <;>
    simp [rationalPortTotal, rationalSeamBoundary, rationalBoundarySection,
      rationalDirectedBoundary, seamLeft, seamRight, Fin.sum_univ_succ] at hx ⊢ <;>
    linarith

/-- The rational boundary image is exactly the zero-total hyperplane. -/
theorem rationalSeamBoundary_range :
    LinearMap.range rationalSeamBoundary = LinearMap.ker rationalPortTotal := by
  ext x
  constructor
  · rintro ⟨c, rfl⟩
    rw [LinearMap.mem_ker]
    change ∑ p : Fin 12, ∑ e : Fin 30,
      rationalDirectedBoundary (seamLeft e) (seamRight e) (c e) p = 0
    rw [Finset.sum_comm]
    apply Finset.sum_eq_zero
    intro e _
    simp [rationalDirectedBoundary, Finset.sum_sub_distrib]
  · intro hx
    rw [LinearMap.mem_ker] at hx
    exact ⟨rationalBoundarySection x, rationalSeamBoundary_section x hx⟩

theorem rationalPortTotal_surjective :
    Function.Surjective rationalPortTotal := by
  intro q
  let x : RationalPortLoad := fun p ↦ if p = 0 then q else 0
  refine ⟨x, ?_⟩
  simp [rationalPortTotal, x]

/-- The ordinary incidence boundary has exact rank eleven. -/
theorem rationalSeamBoundary_range_finrank :
    Module.finrank ℚ (LinearMap.range rationalSeamBoundary) = 11 := by
  rw [rationalSeamBoundary_range]
  have htotal := LinearMap.finrank_range_add_finrank_ker rationalPortTotal
  have hrange : LinearMap.range rationalPortTotal = ⊤ :=
    LinearMap.range_eq_top.mpr rationalPortTotal_surjective
  rw [hrange] at htotal
  simp at htotal ⊢
  omega

/-- The thirty-edge incidence kernel, the integral graph's cycle space after
scalar extension, has exact rank nineteen. -/
theorem rationalSeamBoundary_kernel_finrank :
    Module.finrank ℚ (LinearMap.ker rationalSeamBoundary) = 19 := by
  have hrank := LinearMap.finrank_range_add_finrank_ker rationalSeamBoundary
  rw [rationalSeamBoundary_range_finrank] at hrank
  simp at hrank ⊢
  omega

/-! ## The exact `D6` image -/

/-- Six-axis antipodal-odd current readback.  This is the composition of the
thirty-edge incidence boundary with the independently defined source-load
quotient; it does not rename seams as axis translations. -/
def seamAxisCurrent (c : SeamCurrent) : Axis → ℤ :=
  signedLoadControl (seamBoundary c)

/-- Expanded finite table for the incidence/readback composite.  Keeping this
normal form named makes the subsequent integer proofs small while the next
theorem checks that it is exactly the boundary construction. -/
def seamAxisTable (c : SeamCurrent) : Axis → ℤ :=
  ![
    -c 2 - c 6 + c 11 + c 12 + c 13 + c 14 + c 16 + c 19 - c 25 - c 26,
    -c 0 + c 5 + c 6 + c 7 + c 8 + c 17 + c 22 + c 25 + c 27 - c 29,
    -c 3 + c 8 - c 9 + c 13 + c 15 + c 16 + c 17 + c 18 - c 23 - c 24,
    c 0 + c 1 + c 2 + c 3 + c 4 + c 20 + c 24 + c 26 + c 28 + c 29,
    c 4 - c 7 - c 10 + c 12 + c 15 + c 18 + c 19 + c 20 - c 21 - c 22,
    -c 1 - c 5 + c 9 + c 10 + c 11 + c 14 + c 21 + c 23 - c 27 - c 28
  ]

set_option maxHeartbeats 2000000 in
/-- The displayed six-by-thirty integer table is the exact signed incidence
map, not an independent declaration. -/
theorem seamAxisCurrent_eq_table (c : SeamCurrent) :
    seamAxisCurrent c = seamAxisTable c := by
  funext i
  fin_cases i <;>
    simp [seamAxisCurrent, seamAxisTable, signedLoadControl, seamBoundary,
      directedBoundary, seamLeft, seamRight, positivePort,
      positiveSourcePort, sourceToRERPort, OPH.PortFrameGram.antipode,
      Fin.sum_univ_succ] <;>
    ring

/-- The even-sum root-lattice predicate. -/
def EvenAxisTotal (z : Axis → ℤ) : Prop :=
  Even (∑ i : Axis, z i)

/-- Exact `D6` control lattice. -/
abbrev D6Control := {z : Axis → ℤ // EvenAxisTotal z}

/-- The seam image always has even coordinate sum. -/
theorem seamAxisCurrent_even (c : SeamCurrent) :
    EvenAxisTotal (seamAxisCurrent c) := by
  rw [seamAxisCurrent_eq_table]
  refine ⟨c 4 + c 8 + c 11 + c 12 + c 13 + c 14 + c 15 + c 16 +
    c 17 + c 18 + c 19 + c 20, ?_⟩
  simp [seamAxisTable, Fin.sum_univ_succ]
  ring

/-- Half the total is integral on `D6`. -/
def halfTotal (z : D6Control) : ℤ :=
  (∑ i : Axis, z.1 i) / 2

theorem two_mul_halfTotal (z : D6Control) :
    2 * halfTotal z = ∑ i : Axis, z.1 i := by
  rcases z.2 with ⟨k, hk⟩
  simp only [EvenAxisTotal] at hk
  unfold halfTotal
  omega

/-- Explicit six-seam section of the current map.  Only seams
`5, 8, 12, 13, 22, 24` are used. -/
def d6Section (z : D6Control) : SeamCurrent :=
  ![0, 0, 0, 0, 0,
    -z.1 5,
    0, 0,
    halfTotal z - z.1 0,
    0, 0, 0,
    halfTotal z - z.1 2 - z.1 3,
    halfTotal z - z.1 1 - z.1 4 - z.1 5,
    0, 0, 0, 0, 0, 0, 0, 0,
    halfTotal z - z.1 2 - z.1 3 - z.1 4,
    0,
    z.1 3,
    0, 0, 0, 0, 0]

/-- The section is exact over the integers. -/
theorem seamAxisCurrent_d6Section (z : D6Control) :
    seamAxisCurrent (d6Section z) = z.1 := by
  have hhalf := two_mul_halfTotal z
  rw [seamAxisCurrent_eq_table]
  funext i
  fin_cases i <;>
    simp [seamAxisTable, d6Section, Fin.sum_univ_succ] at hhalf ⊢ <;>
    omega

/-- Exact image theorem: thirty-seam currents reach precisely `D6`. -/
theorem exists_seamCurrent_iff_even (z : Axis → ℤ) :
    (∃ c : SeamCurrent, seamAxisCurrent c = z) ↔ EvenAxisTotal z := by
  constructor
  · rintro ⟨c, rfl⟩
    exact seamAxisCurrent_even c
  · intro hz
    let q : D6Control := ⟨z, hz⟩
    exact ⟨d6Section q, seamAxisCurrent_d6Section q⟩

/-- Every doubled six-axis record is in the seam-current image. -/
theorem two_mul_control_in_image (z : Axis → ℤ) :
    ∃ c : SeamCurrent, seamAxisCurrent c = fun i ↦ 2 * z i := by
  apply (exists_seamCurrent_iff_even _).2
  refine ⟨∑ i : Axis, z i, ?_⟩
  rw [← Finset.mul_sum]
  ring

/-! ## The parity cokernel, without quotient overstatement -/

/-- The residual parity coordinate. -/
def axisParity (z : Axis → ℤ) : ℤ :=
  (∑ i : Axis, z i) % 2

theorem axisParity_zero_or_one (z : Axis → ℤ) :
    axisParity z = 0 ∨ axisParity z = 1 := by
  exact Int.emod_two_eq_zero_or_one (∑ i : Axis, z i)

/-- Vanishing parity is exactly membership in the seam-current image. -/
theorem axisParity_eq_zero_iff_in_image (z : Axis → ℤ) :
    axisParity z = 0 ↔ ∃ c : SeamCurrent, seamAxisCurrent c = z := by
  rw [exists_seamCurrent_iff_even]
  exact Int.even_iff.symm

/-- One odd coordinate represents the nontrivial parity coset. -/
def oddUnit : Axis → ℤ := ![1, 0, 0, 0, 0, 0]

theorem oddUnit_parity : axisParity oddUnit = 1 := by
  decide

/-- There are exactly two image cosets: an arbitrary control is in `D6`, or
subtracting the fixed odd unit puts it in `D6`. -/
theorem image_or_oddUnit_coset (z : Axis → ℤ) :
    (∃ c : SeamCurrent, seamAxisCurrent c = z) ∨
      (∃ c : SeamCurrent, seamAxisCurrent c = z - oddUnit) := by
  rcases axisParity_zero_or_one z with hzero | hone
  · exact Or.inl ((axisParity_eq_zero_iff_in_image z).1 hzero)
  · right
    apply (axisParity_eq_zero_iff_in_image (z - oddUnit)).1
    simp [axisParity, oddUnit, Fin.sum_univ_succ] at hone ⊢
    have hdiv := Int.emod_add_mul_ediv (∑ i : Axis, z i) 2
    omega

/-! ## Dense carrier image and completion -/

/-- The standard dense controls from the Gram construction lie in `D6`. -/
theorem denseIntegerControl_even (p : Fin 3 → ℤ × ℤ) :
    EvenAxisTotal (denseIntegerControl p) := by
  refine ⟨(p 0).2 + (p 1).2 + (p 2).2, ?_⟩
  simp [denseIntegerControl, Fin.sum_univ_succ]
  ring

/-- Exact `D6` wrapper equipped below with the response-selected Gram metric,
not the inherited discrete metric of integer tuples. -/
@[ext]
structure D6Point where
  control : Axis → ℤ
  even_total : EvenAxisTotal control

/-- Inclusion into the previously proved integer-control carrier. -/
def d6IntegerPoint (p : D6Point) : IntegerFramePoint := ⟨p.control⟩

/-- `D6` remains injective in the exact Gram-frame image. -/
theorem d6IntegerPoint_injective : Function.Injective d6IntegerPoint := by
  intro p q h
  ext i
  exact congrFun (congrArg IntegerFramePoint.control h) i

/-- Position of an even-sum seam quotient record in the selected carrier. -/
noncomputable def d6Position (p : D6Point) : EuclideanVec3 :=
  pointEuclideanFrame (d6IntegerPoint p)

theorem d6Position_injective : Function.Injective d6Position :=
  pointEuclideanFrame_injective.comp d6IntegerPoint_injective

/-- The even-sum seam-current quotient is still dense in the same
three-dimensional carrier. -/
theorem d6Position_denseRange : DenseRange d6Position := by
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
  have hframe : DenseRange
      (fun p : D6Point ↦ pointFrame (d6IntegerPoint p)) := by
    apply hproduct.mono
    rintro _ ⟨p, rfl⟩
    let q : D6Point := ⟨denseIntegerControl p, denseIntegerControl_even p⟩
    exact ⟨q, integerFrame_denseIntegerControl p⟩
  have hbase : DenseRange
      (fun p : D6Point ↦
        euclideanEquivVec3.symm (pointFrame (d6IntegerPoint p))) :=
    euclideanEquivVec3.symm.surjective.denseRange.comp hframe
      euclideanEquivVec3.symm.continuous
  have hscale : Function.Surjective
      (fun x : EuclideanVec3 ↦ rawRadius⁻¹ • x) := by
    intro y
    refine ⟨rawRadius • y, ?_⟩
    change rawRadius⁻¹ • (rawRadius • y) = y
    rw [← mul_smul]
    simp [rawRadius_ne_zero]
  simpa only [d6Position, pointEuclideanFrame, Function.comp_apply] using
    hscale.denseRange.comp hbase (by fun_prop)

/-- Pull back the response Gram metric to the exact seam quotient. -/
noncomputable instance d6PointMetric : MetricSpace D6Point :=
  MetricSpace.induced d6Position d6Position_injective inferInstance

theorem d6Position_isometry : Isometry d6Position :=
  MetricSpace.isometry_induced d6Position d6Position_injective

/-- The same Euclidean three-carrier is an abstract completion of the exact
even-sum seam-current quotient. -/
noncomputable def d6Completion : AbstractCompletion D6Point where
  space := EuclideanVec3
  coe := d6Position
  uniformStruct := inferInstance
  complete := inferInstance
  separation := inferInstance
  isUniformInducing := d6Position_isometry.isUniformInducing
  dense := d6Position_denseRange

/-- Canonical comparison equivalence from the metric completion of the seam
quotient to the response-selected Euclidean carrier. -/
noncomputable def d6CompletionEquivEuclidean3 :
    UniformSpace.Completion D6Point ≃ᵤ EuclideanVec3 :=
  AbstractCompletion.compareEquiv UniformSpace.Completion.cPkg d6Completion

end OPH.SeamCurrentCarrierQuotient

/- Axiom audit: finite incidence checks, exact integer arithmetic, and the
existing response-selected completion theorems only. -/

#print axioms OPH.SeamCurrentCarrierQuotient.seam_table_complete
#print axioms OPH.SeamCurrentCarrierQuotient.rationalSeamBoundary_range_finrank
#print axioms OPH.SeamCurrentCarrierQuotient.rationalSeamBoundary_kernel_finrank
#print axioms OPH.SeamCurrentCarrierQuotient.exists_seamCurrent_iff_even
#print axioms OPH.SeamCurrentCarrierQuotient.axisParity_eq_zero_iff_in_image
#print axioms OPH.SeamCurrentCarrierQuotient.d6Position_denseRange
#print axioms OPH.SeamCurrentCarrierQuotient.d6CompletionEquivEuclidean3
