import EventAlgebra.PartitionPinching

/-!
# Block averaging onto the commutative partition span

A projective partition induces a second canonical expectation besides the
block pinching of `EventAlgebra.PartitionPinching`. The **partition span**

  `D = span{Pᵢ}`

is the linear span of the projectors. It is closed under multiplication and
conjugate transposition, it is commutative, it contains the identity, and it
sits inside the commutant `N = {X : X Pᵢ = Pᵢ X}`; in fact every element of
`D` commutes with every element of `N`, so `D` lies in the center of `N`.
The **partition average**

  `partitionAverage part X = ∑ i, (Tr Pᵢ)⁻¹ · Tr(X Pᵢ) • Pᵢ`

is the trace-preserving projection onto `D`. Zero projectors are harmless
throughout: Lean's zero-totalized field inverse makes their terms vanish on
both sides of every identity, so no nonzero-block hypothesis appears.

The module proves:

* **membership, fixed points, exact range**: the average lands in `D`, fixes
  `D` pointwise, and its bundled linear map has range exactly `D`;
* **map structure**: linear, unital, positive, trace preserving, idempotent;
* **tower laws**: `𝔸 ∘ 𝔼 = 𝔸` and `𝔼 ∘ 𝔸 = 𝔸` for the partition pinching
  `𝔼` and the partition average `𝔸`;
* **statistics preservation**: the average leaves every Born weight of a
  partition member unchanged;
* **trace duality and uniqueness**: against `D` the average is invisible to
  the trace pairing, and it is the unique `D`-valued linear map with that
  property;
* **conditioning compatibility**: for a partition member of nonzero weight,
  the L\"uders update and the partition average commute, and both composites
  equal the normalized projector `(Tr Pⱼ)⁻¹ • Pⱼ`.

No claim is made that `D` equals the center of the commutant (only the
inclusion is proved), and no completely-positive channel bundle is claimed.

## Tagging convention

As in `EventAlgebra.Basic`: each lemma is tagged **algebra-only** or
**trace-dependent** according to whether its content passes through the
trace pairing.
-/

namespace EventAlgebra

open Matrix
open scoped ComplexOrder

variable {n k : ℕ}

namespace ProjectivePartition

variable (part : ProjectivePartition n k)

/-- The **partition span** `D = span{Pᵢ}`: the linear span of the
projectors of the partition. The closure lemmas below show that it is a
unital, star-closed, commutative subalgebra of the matrix algebra. -/
noncomputable def span : Submodule ℂ (Matrix (Fin n) (Fin n) ℂ) :=
  Submodule.span ℂ (Set.range part.proj)

/-- **Algebra-only.** Every projector of the partition belongs to the
partition span. -/
theorem proj_mem_span (i : Fin k) : part.proj i ∈ part.span :=
  Submodule.subset_span ⟨i, rfl⟩

/-- **Algebra-only.** The identity belongs to the partition span, because
the partition is complete. -/
theorem one_mem_span : (1 : Matrix (Fin n) (Fin n) ℂ) ∈ part.span := by
  rw [← part.complete]
  exact Submodule.sum_mem _ fun i _ => part.proj_mem_span i

/-- **Algebra-only.** The partition span is closed under conjugate
transposition. -/
theorem conjTranspose_mem_span {X : Matrix (Fin n) (Fin n) ℂ}
    (hX : X ∈ part.span) : Xᴴ ∈ part.span := by
  induction hX using Submodule.span_induction with
  | mem Y hY =>
      obtain ⟨i, rfl⟩ := hY
      rw [(part.isEvent i).1.eq]
      exact part.proj_mem_span i
  | zero =>
      rw [conjTranspose_zero]
      exact Submodule.zero_mem _
  | add x y hx hy ihx ihy =>
      rw [conjTranspose_add]
      exact Submodule.add_mem _ ihx ihy
  | smul c x hx ih =>
      rw [conjTranspose_smul]
      exact Submodule.smul_mem _ _ ih

/-- **Algebra-only.** The partition span is closed under multiplication. -/
theorem mul_mem_span {X Y : Matrix (Fin n) (Fin n) ℂ}
    (hX : X ∈ part.span) (hY : Y ∈ part.span) : X * Y ∈ part.span := by
  refine Submodule.span_induction₂
    (p := fun x y _ _ => x * y ∈ part.span) ?_ ?_ ?_ ?_ ?_ ?_ ?_ hX hY
  · rintro x y ⟨i, rfl⟩ ⟨j, rfl⟩
    rw [part.proj_mul_proj]
    split_ifs
    · exact part.proj_mem_span i
    · exact Submodule.zero_mem _
  · intro y _
    rw [zero_mul]
    exact Submodule.zero_mem _
  · intro x _
    rw [mul_zero]
    exact Submodule.zero_mem _
  · intro x y z _ _ _ h1 h2
    rw [add_mul]
    exact Submodule.add_mem _ h1 h2
  · intro x y z _ _ _ h1 h2
    rw [mul_add]
    exact Submodule.add_mem _ h1 h2
  · intro r x y _ _ h
    rw [smul_mul_assoc]
    exact Submodule.smul_mem _ _ h
  · intro r x y _ _ h
    rw [mul_smul_comm]
    exact Submodule.smul_mem _ _ h

/-- **Algebra-only.** The partition span is commutative. -/
theorem span_mul_comm {X Y : Matrix (Fin n) (Fin n) ℂ}
    (hX : X ∈ part.span) (hY : Y ∈ part.span) : X * Y = Y * X := by
  refine Submodule.span_induction₂
    (p := fun x y _ _ => x * y = y * x) ?_ ?_ ?_ ?_ ?_ ?_ ?_ hX hY
  · rintro x y ⟨i, rfl⟩ ⟨j, rfl⟩
    exact part.proj_commute i j
  · intro y _
    rw [zero_mul, mul_zero]
  · intro x _
    rw [mul_zero, zero_mul]
  · intro x y z _ _ _ h1 h2
    rw [add_mul, mul_add, h1, h2]
  · intro x y z _ _ _ h1 h2
    rw [mul_add, add_mul, h1, h2]
  · intro r x y _ _ h
    rw [smul_mul_assoc, mul_smul_comm, h]
  · intro r x y _ _ h
    rw [mul_smul_comm, smul_mul_assoc, h]

/-- **Algebra-only.** Every element of the partition span commutes with
every projector of the partition. -/
theorem commute_proj_of_mem_span {X : Matrix (Fin n) (Fin n) ℂ}
    (hX : X ∈ part.span) (i : Fin k) :
    X * part.proj i = part.proj i * X := by
  induction hX using Submodule.span_induction with
  | mem Y hY =>
      obtain ⟨j, rfl⟩ := hY
      exact part.proj_commute j i
  | zero => rw [zero_mul, mul_zero]
  | add x y hx hy ihx ihy => rw [add_mul, mul_add, ihx, ihy]
  | smul c x hx ih => rw [smul_mul_assoc, mul_smul_comm, ih]

/-- **Algebra-only.** The partition span is contained in the partition
commutant: `D ⊆ N`. -/
theorem span_le_commutant {X : Matrix (Fin n) (Fin n) ℂ}
    (hX : X ∈ part.span) : X ∈ part.commutant :=
  (ProjectivePartition.mem_commutant_iff part).mpr
    (part.commute_proj_of_mem_span hX)

/-- **Algebra-only.** Every element of the partition span commutes with
every element of the commutant: `D` lies in the center of `N`. -/
theorem mul_comm_of_mem_span_of_mem_commutant
    {X Y : Matrix (Fin n) (Fin n) ℂ}
    (hX : X ∈ part.span) (hY : Y ∈ part.commutant) : X * Y = Y * X := by
  have hY' := (ProjectivePartition.mem_commutant_iff part).mp hY
  induction hX using Submodule.span_induction with
  | mem Z hZ =>
      obtain ⟨i, rfl⟩ := hZ
      exact (hY' i).symm
  | zero => rw [zero_mul, mul_zero]
  | add x y hx hy ihx ihy => rw [add_mul, mul_add, ihx, ihy]
  | smul c x hx ih => rw [smul_mul_assoc, mul_smul_comm, ih]

end ProjectivePartition

/-- The **partition average**: the trace-preserving projection onto the
partition span, `X ↦ ∑ i, (Tr Pᵢ)⁻¹ · Tr(X Pᵢ) • Pᵢ`. Lean's
zero-totalized inverse makes the terms of zero projectors vanish. -/
noncomputable def partitionAverage (part : ProjectivePartition n k)
    (X : Matrix (Fin n) (Fin n) ℂ) : Matrix (Fin n) (Fin n) ℂ :=
  ∑ i, (((part.proj i).trace)⁻¹ * bornWeight X (part.proj i)) • part.proj i

/-- The partition average bundled as a complex-linear map. -/
noncomputable def partitionAverageLinearMap (part : ProjectivePartition n k) :
    Matrix (Fin n) (Fin n) ℂ →ₗ[ℂ] Matrix (Fin n) (Fin n) ℂ where
  toFun := partitionAverage part
  map_add' X Y := by
    simp only [partitionAverage, bornWeight, add_mul, trace_add, mul_add,
      add_smul, Finset.sum_add_distrib]
  map_smul' c X := by
    simp only [partitionAverage, bornWeight_smul, Finset.smul_sum, smul_smul,
      RingHom.id_apply]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [mul_left_comm]

@[simp]
theorem partitionAverageLinearMap_apply (part : ProjectivePartition n k)
    (X : Matrix (Fin n) (Fin n) ℂ) :
    partitionAverageLinearMap part X = partitionAverage part X :=
  rfl

variable (part : ProjectivePartition n k)

/-- **Trace-dependent.** Additivity of the partition average. -/
theorem partitionAverage_add (X Y : Matrix (Fin n) (Fin n) ℂ) :
    partitionAverage part (X + Y) =
      partitionAverage part X + partitionAverage part Y :=
  (partitionAverageLinearMap part).map_add X Y

/-- **Trace-dependent.** Homogeneity of the partition average. -/
theorem partitionAverage_smul (c : ℂ) (X : Matrix (Fin n) (Fin n) ℂ) :
    partitionAverage part (c • X) = c • partitionAverage part X :=
  (partitionAverageLinearMap part).map_smul c X

/-- **Trace-dependent.** The partition average lands in the partition
span. -/
theorem partitionAverage_mem_span (X : Matrix (Fin n) (Fin n) ℂ) :
    partitionAverage part X ∈ part.span :=
  Submodule.sum_mem _ fun i _ =>
    Submodule.smul_mem _ _ (part.proj_mem_span i)

/-- **Trace-dependent.** The partition average fixes every projector of
the partition. The zero-projector case degenerates on both sides. -/
theorem partitionAverage_proj (j : Fin k) :
    partitionAverage part (part.proj j) = part.proj j := by
  rw [partitionAverage,
    Finset.sum_eq_single j
      (fun i _ hij => by
        have hbw : bornWeight (part.proj j) (part.proj i) = 0 := by
          rw [bornWeight, part.proj_mul_proj]
          simp [Ne.symm hij]
        rw [hbw, mul_zero, zero_smul])
      (fun h => absurd (Finset.mem_univ j) h)]
  have hbw : bornWeight (part.proj j) (part.proj j) = (part.proj j).trace := by
    rw [bornWeight, (part.isEvent j).2]
  rw [hbw]
  rcases eq_or_ne ((part.proj j).trace) 0 with h0 | h0
  · have hzero : part.proj j = 0 :=
      ((part.isEvent j).posSemidef.trace_eq_zero_iff).mp h0
    rw [hzero, smul_zero]
  · rw [inv_mul_cancel₀ h0, one_smul]

/-- **Trace-dependent.** The partition average fixes the partition span
pointwise. -/
theorem partitionAverage_fixes {X : Matrix (Fin n) (Fin n) ℂ}
    (hX : X ∈ part.span) : partitionAverage part X = X := by
  induction hX using Submodule.span_induction with
  | mem Y hY =>
      obtain ⟨j, rfl⟩ := hY
      exact partitionAverage_proj part j
  | zero => simp [partitionAverage, bornWeight]
  | add x y hx hy ihx ihy => rw [partitionAverage_add, ihx, ihy]
  | smul c x hx ih => rw [partitionAverage_smul, ih]

/-- **Trace-dependent.** Exact fixed-point characterization of the
partition average in terms of the partition span. -/
theorem partitionAverage_eq_self_iff_mem_span
    (X : Matrix (Fin n) (Fin n) ℂ) :
    partitionAverage part X = X ↔ X ∈ part.span := by
  constructor
  · intro h
    rw [← h]
    exact partitionAverage_mem_span part X
  · exact partitionAverage_fixes part

/-- **Trace-dependent.** The partition average is idempotent. -/
theorem partitionAverage_idem (X : Matrix (Fin n) (Fin n) ℂ) :
    partitionAverage part (partitionAverage part X) =
      partitionAverage part X :=
  partitionAverage_fixes part (partitionAverage_mem_span part X)

/-- **Trace-dependent.** The partition average is unital. -/
theorem partitionAverage_unital : partitionAverage part 1 = 1 :=
  partitionAverage_fixes part part.one_mem_span

/-- **Trace-dependent.** The linear range of the bundled partition average
is exactly the partition span. -/
theorem range_partitionAverageLinearMap :
    LinearMap.range (partitionAverageLinearMap part) = part.span := by
  ext X
  constructor
  · rintro ⟨Y, rfl⟩
    exact partitionAverage_mem_span part Y
  · intro hX
    exact ⟨X, partitionAverage_fixes part hX⟩

/-- **Trace-dependent.** The partition average preserves the trace of
every matrix. -/
theorem trace_partitionAverage (X : Matrix (Fin n) (Fin n) ℂ) :
    (partitionAverage part X).trace = X.trace := by
  rw [partitionAverage, trace_sum]
  have hterm : ∀ i : Fin k,
      ((((part.proj i).trace)⁻¹ * bornWeight X (part.proj i)) •
        part.proj i).trace = bornWeight X (part.proj i) := by
    intro i
    rw [trace_smul, smul_eq_mul]
    rcases eq_or_ne ((part.proj i).trace) 0 with h0 | h0
    · have hzero : part.proj i = 0 :=
        ((part.isEvent i).posSemidef.trace_eq_zero_iff).mp h0
      rw [h0, mul_zero, hzero, bornWeight, mul_zero, trace_zero]
    · rw [mul_comm (((part.proj i).trace)⁻¹) _, mul_assoc,
        inv_mul_cancel₀ h0, mul_one]
  rw [Finset.sum_congr rfl fun i _ => hterm i, ← bornWeight_sum,
    part.complete, bornWeight, mul_one]

/-- **Trace-dependent.** The partition average preserves positive
semidefiniteness: the coefficients are nonnegative and the projectors are
positive. -/
theorem partitionAverage_posSemidef {X : Matrix (Fin n) (Fin n) ℂ}
    (hX : X.PosSemidef) : (partitionAverage part X).PosSemidef := by
  refine posSemidef_sum _ fun i _ => ?_
  have htr : (0 : ℂ) ≤ (part.proj i).trace :=
    (part.isEvent i).posSemidef.trace_nonneg
  have hinv : (0 : ℂ) ≤ ((part.proj i).trace)⁻¹ := by
    rcases eq_or_ne ((part.proj i).trace) 0 with h0 | h0
    · rw [h0, _root_.inv_zero]
    · exact (RCLike.inv_pos.mpr (lt_of_le_of_ne htr (Ne.symm h0))).le
  have hbw : (0 : ℂ) ≤ bornWeight X (part.proj i) :=
    bornWeight_nonneg hX (part.isEvent i)
  exact ((part.isEvent i).posSemidef).smul (mul_nonneg hinv hbw)

/-- **Trace-dependent.** The partition average maps states to states. -/
theorem partitionAverage_isState {ρ : Matrix (Fin n) (Fin n) ℂ}
    (hρ : IsState ρ) : IsState (partitionAverage part ρ) :=
  ⟨partitionAverage_posSemidef part hρ.1,
    (trace_partitionAverage part ρ).trans hρ.2⟩

/-- **Trace-dependent.** First tower law: averaging after pinching is
averaging, `𝔸 ∘ 𝔼 = 𝔸`. -/
theorem partitionAverage_partitionPinching (X : Matrix (Fin n) (Fin n) ℂ) :
    partitionAverage part (partitionPinching part X) =
      partitionAverage part X := by
  rw [partitionAverage, partitionAverage]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [bornWeight_partitionPinching part X (part.proj i)
    fun j => part.proj_commute i j]

/-- **Trace-dependent.** Second tower law: pinching after averaging is
averaging, `𝔼 ∘ 𝔸 = 𝔸`, because the average lands in the span, which the
pinching fixes. -/
theorem partitionPinching_partitionAverage (X : Matrix (Fin n) (Fin n) ℂ) :
    partitionPinching part (partitionAverage part X) =
      partitionAverage part X :=
  partitionPinching_fixes part
    (part.commute_proj_of_mem_span (partitionAverage_mem_span part X))

/-- **Trace-dependent.** Against any element of the partition span, the
partition average is invisible to the trace pairing. This is the defining
conditional-expectation property of the average relative to `D`. -/
theorem trace_partitionAverage_mul_of_mem
    {C : Matrix (Fin n) (Fin n) ℂ} (X : Matrix (Fin n) (Fin n) ℂ)
    (hC : C ∈ part.span) :
    (partitionAverage part X * C).trace = (X * C).trace := by
  induction hC using Submodule.span_induction with
  | mem Y hY =>
      obtain ⟨j, rfl⟩ := hY
      rw [partitionAverage, Finset.sum_mul, trace_sum,
        Finset.sum_eq_single j
          (fun i _ hij => by
            rw [smul_mul_assoc, part.proj_mul_proj]
            simp [hij])
          (fun h => absurd (Finset.mem_univ j) h)]
      rw [smul_mul_assoc, (part.isEvent j).2, trace_smul, smul_eq_mul]
      rcases eq_or_ne ((part.proj j).trace) 0 with h0 | h0
      · have hzero : part.proj j = 0 :=
          ((part.isEvent j).posSemidef.trace_eq_zero_iff).mp h0
        rw [h0, _root_.inv_zero, zero_mul, zero_mul, hzero, mul_zero,
          trace_zero]
      · rw [mul_comm (((part.proj j).trace)⁻¹) _, mul_assoc,
          inv_mul_cancel₀ h0, mul_one, bornWeight]
  | zero => rw [mul_zero, mul_zero]
  | add x y hx hy ihx ihy =>
      rw [mul_add, mul_add, trace_add, trace_add, ihx, ihy]
  | smul c x hx ih =>
      rw [mul_smul_comm, mul_smul_comm, trace_smul, trace_smul, ih]

/-- **Trace-dependent.** The partition average preserves the Born weight
of every partition member: the classical shadow carries the exact Born
statistics of the state. -/
theorem bornWeight_partitionAverage (ρ : Matrix (Fin n) (Fin n) ℂ)
    (j : Fin k) :
    bornWeight (partitionAverage part ρ) (part.proj j) =
      bornWeight ρ (part.proj j) :=
  trace_partitionAverage_mul_of_mem part ρ (part.proj_mem_span j)

/-- **Trace-dependent.** Pointwise trace-duality characterization: any
span-valued matrix that is trace-compatible with `X` against the whole
partition span equals the partition average of `X`. -/
theorem partitionAverage_unique {X Y : Matrix (Fin n) (Fin n) ℂ}
    (hYmem : Y ∈ part.span)
    (hYtr : ∀ C, C ∈ part.span → (Y * C).trace = (X * C).trace) :
    Y = partitionAverage part X := by
  set G := Y - partitionAverage part X with hG
  have hGmem : G ∈ part.span :=
    Submodule.sub_mem _ hYmem (partitionAverage_mem_span part X)
  have hGH : Gᴴ ∈ part.span := part.conjTranspose_mem_span hGmem
  have htr0 : (G * Gᴴ).trace = 0 := by
    have h1 : (Y * Gᴴ).trace = (X * Gᴴ).trace := hYtr Gᴴ hGH
    have h2 : (partitionAverage part X * Gᴴ).trace = (X * Gᴴ).trace :=
      trace_partitionAverage_mul_of_mem part X hGH
    rw [hG, sub_mul, trace_sub, h1, h2, sub_self]
  have hG0 : G = 0 := trace_mul_conjTranspose_self_eq_zero_iff.mp htr0
  rw [← sub_eq_zero]
  exact hG0

/-- **Trace-dependent.** Map-level uniqueness: a complex-linear map whose
values lie in the partition span and which has the same trace pairing
against every span element is the bundled partition average. -/
theorem partitionAverageLinearMap_unique
    (F : Matrix (Fin n) (Fin n) ℂ →ₗ[ℂ] Matrix (Fin n) (Fin n) ℂ)
    (hRange : ∀ X, F X ∈ part.span)
    (hTrace : ∀ X C, C ∈ part.span → (F X * C).trace = (X * C).trace) :
    F = partitionAverageLinearMap part := by
  apply LinearMap.ext
  intro X
  exact partitionAverage_unique part (hRange X) (fun C hC => hTrace X C hC)

/-- **Trace-dependent.** Averaging a L\"uders update by a partition member
of nonzero weight collapses it to the normalized projector: the classical
record of the conditioned state is the indicator of the observed outcome. -/
theorem partitionAverage_luedersUpdate_proj (ρ : Matrix (Fin n) (Fin n) ℂ)
    (j : Fin k) (hw : bornWeight ρ (part.proj j) ≠ 0) :
    partitionAverage part (luedersUpdate ρ (part.proj j)) =
      ((part.proj j).trace)⁻¹ • part.proj j := by
  rw [partitionAverage,
    Finset.sum_eq_single j
      (fun i _ hij => by
        have hzero : bornWeight (luedersUpdate ρ (part.proj j))
            (part.proj i) = 0 := by
          rw [luedersUpdate, bornWeight, smul_mul_assoc, trace_smul,
            smul_eq_mul]
          have hmat : part.proj j * ρ * part.proj j * part.proj i = 0 := by
            rw [mul_assoc (part.proj j * ρ),
              part.orthogonal j i (Ne.symm hij), mul_zero]
          rw [hmat, trace_zero, mul_zero]
        rw [hzero, mul_zero, zero_smul])
      (fun h => absurd (Finset.mem_univ j) h)]
  rw [bornWeight_luedersUpdate_self (part.isEvent j) hw, mul_one]

/-- **Trace-dependent.** Conditioning the partition average of a state on
a partition member of nonzero weight also produces the normalized
projector. -/
theorem luedersUpdate_partitionAverage_proj (ρ : Matrix (Fin n) (Fin n) ℂ)
    (j : Fin k) (hw : bornWeight ρ (part.proj j) ≠ 0) :
    luedersUpdate (partitionAverage part ρ) (part.proj j) =
      ((part.proj j).trace)⁻¹ • part.proj j := by
  have hsand : part.proj j * partitionAverage part ρ * part.proj j =
      (((part.proj j).trace)⁻¹ * bornWeight ρ (part.proj j)) •
        part.proj j := by
    rw [partitionAverage, Finset.mul_sum, Finset.sum_mul,
      Finset.sum_eq_single j
        (fun i _ hij => by
          rw [mul_smul_comm, smul_mul_assoc]
          have hmat : part.proj j * part.proj i * part.proj j = 0 := by
            rw [part.orthogonal j i (Ne.symm hij), zero_mul]
          rw [hmat, smul_zero])
        (fun h => absurd (Finset.mem_univ j) h)]
    rw [mul_smul_comm, smul_mul_assoc, (part.isEvent j).2,
      (part.isEvent j).2]
  rw [luedersUpdate, bornWeight_partitionAverage part ρ j, hsand, smul_smul]
  congr 1
  rw [mul_comm (((part.proj j).trace)⁻¹) _, ← mul_assoc,
    inv_mul_cancel₀ hw, one_mul]

/-- **Trace-dependent.** **Classical conditioning**: for a partition member
of nonzero weight, the partition average commutes with the L\"uders
update. Both composites are the normalized projector of the observed
outcome. -/
theorem partitionAverage_luedersUpdate (ρ : Matrix (Fin n) (Fin n) ℂ)
    (j : Fin k) (hw : bornWeight ρ (part.proj j) ≠ 0) :
    partitionAverage part (luedersUpdate ρ (part.proj j)) =
      luedersUpdate (partitionAverage part ρ) (part.proj j) :=
  (partitionAverage_luedersUpdate_proj part ρ j hw).trans
    (luedersUpdate_partitionAverage_proj part ρ j hw).symm

-- Axiom audit: each must report only `[propext, Classical.choice, Quot.sound]`.
#print axioms ProjectivePartition.proj_mem_span
#print axioms ProjectivePartition.one_mem_span
#print axioms ProjectivePartition.conjTranspose_mem_span
#print axioms ProjectivePartition.mul_mem_span
#print axioms ProjectivePartition.span_mul_comm
#print axioms ProjectivePartition.commute_proj_of_mem_span
#print axioms ProjectivePartition.span_le_commutant
#print axioms ProjectivePartition.mul_comm_of_mem_span_of_mem_commutant
#print axioms partitionAverageLinearMap_apply
#print axioms partitionAverage_add
#print axioms partitionAverage_smul
#print axioms partitionAverage_mem_span
#print axioms partitionAverage_proj
#print axioms partitionAverage_fixes
#print axioms partitionAverage_eq_self_iff_mem_span
#print axioms partitionAverage_idem
#print axioms partitionAverage_unital
#print axioms range_partitionAverageLinearMap
#print axioms trace_partitionAverage
#print axioms partitionAverage_posSemidef
#print axioms partitionAverage_isState
#print axioms partitionAverage_partitionPinching
#print axioms partitionPinching_partitionAverage
#print axioms trace_partitionAverage_mul_of_mem
#print axioms bornWeight_partitionAverage
#print axioms partitionAverage_unique
#print axioms partitionAverageLinearMap_unique
#print axioms partitionAverage_luedersUpdate_proj
#print axioms luedersUpdate_partitionAverage_proj
#print axioms partitionAverage_luedersUpdate

end EventAlgebra
