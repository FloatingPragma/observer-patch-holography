import Mathlib
import A5PortModule
import A5PrimitivePortPrediction

open scoped BigOperators

namespace OPH.PrimitivePortTranslationBridge

/-!
# Conditional primitive-port translation bridge

This module places a declared twelve-port finite-range translation law on an
already continuous three-dimensional scalar field.  It proves inverse pairing,
plane-wave diagonalization, nonnegativity, port-relabeling covariance, and
passive orthogonal-frame covariance.  An explicit moment-packet premise then
gives the frozen fourth- and sixth-order coefficient polynomial.

The theorem does not select the translation law from repair dynamics, construct
the spatial frame from the source, fix the positive scale `a`, identify a photon
or another physical sector, supply a physical time equation, or prove a frame,
boost, nuisance, or exclusivity statement.  The abstract shifts are not a site
orbit or a locally finite lattice.
-/

/-- Cartesian three-space represented by three real coordinates. -/
abbrev Vec3 := Fin 3 → ℝ

/-- Exact Cartesian dot product. -/
noncomputable def dot (x y : Vec3) : ℝ :=
  ∑ d : Fin 3, x d * y d

theorem dot_add_right (x y z : Vec3) :
    dot x (y + z) = dot x y + dot x z := by
  unfold dot
  rw [← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro d _
  simp
  ring

theorem dot_smul_right (a : ℝ) (x y : Vec3) :
    dot x (a • y) = a * dot x y := by
  unfold dot
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro d _
  simp
  ring

theorem dot_neg_right (x y : Vec3) : dot x (-y) = -dot x y := by
  simp [dot]

theorem dot_neg_left (x y : Vec3) : dot (-x) y = -dot x y := by
  simp [dot]

/-! ## Paired translations and their reciprocal symbol -/

/-- Translation of a complex scalar field by a declared displacement. -/
noncomputable def shift (v : Vec3) (f : Vec3 → ℂ) : Vec3 → ℂ :=
  fun x => f (x + v)

/-- Symmetrized finite difference for a displacement and its inverse. -/
noncomputable def pairedDifference
    (v : Vec3) (f : Vec3 → ℂ) : Vec3 → ℂ :=
  fun x => f x - (f (x + v) + f (x - v)) / 2

/-- Opposite translations are exact inverses. -/
theorem shift_inverse (v : Vec3) (f : Vec3 → ℂ) :
    shift (-v) (shift v f) = f := by
  funext x
  simp [shift]

/-- Reversing the signed displacement leaves the paired difference fixed. -/
theorem pairedDifference_neg (v : Vec3) (f : Vec3 → ℂ) :
    pairedDifference (-v) f = pairedDifference v f := by
  funext x
  simp [pairedDifference, sub_eq_add_neg, add_comm]

/-- Scalar plane wave on the conditional translation space. -/
noncomputable def planeWave (k : Vec3) : Vec3 → ℂ :=
  fun x => Complex.exp ((dot k x : ℂ) * Complex.I)

/-- Translation acts diagonally on a plane wave. -/
theorem planeWave_shift (k x v : Vec3) :
    planeWave k (x + v) =
      Complex.exp ((dot k v : ℂ) * Complex.I) * planeWave k x := by
  unfold planeWave
  rw [dot_add_right, Complex.ofReal_add, add_mul, Complex.exp_add]
  ring

/-- The symmetrized inverse pair has cosine eigenvalue `1 - cos θ`. -/
theorem pairedDifference_planeWave (k x v : Vec3) :
    pairedDifference v (planeWave k) x =
      ((1 - Real.cos (dot k v) : ℝ) : ℂ) * planeWave k x := by
  unfold pairedDifference
  rw [show planeWave k (x + v) =
      Complex.exp ((dot k v : ℂ) * Complex.I) * planeWave k x from
    planeWave_shift k x v]
  rw [show planeWave k (x - v) =
      Complex.exp ((dot k (-v) : ℂ) * Complex.I) * planeWave k x by
    rw [sub_eq_add_neg]
    exact planeWave_shift k x (-v)]
  rw [dot_neg_right]
  rw [Complex.exp_ofReal_mul_I, Complex.exp_ofReal_mul_I]
  simp only [Real.cos_neg, Real.sin_neg, Complex.ofReal_neg]
  push_cast
  ring

/-- Equal-weight twelve-port finite-range translation operator. -/
noncomputable def primitivePortOperator
    (a : ℝ) (u : Fin 12 → Vec3) (f : Vec3 → ℂ) : Vec3 → ℂ :=
  fun x => (1 / (2 * a ^ 2) : ℂ) *
    ∑ i : Fin 12, pairedDifference (a • u i) f x

/-- Relabeling the complete port orbit leaves the operator unchanged. -/
theorem primitivePortOperator_reindex
    (a : ℝ) (u : Fin 12 → Vec3) (f : Vec3 → ℂ)
    (g : Equiv.Perm (Fin 12)) :
    primitivePortOperator a (fun i => u (g i)) f =
      primitivePortOperator a u f := by
  funext x
  unfold primitivePortOperator
  congr 1
  exact g.bijective.sum_comp
    (fun i : Fin 12 => pairedDifference (a • u i) f x)

/-- Phase sampled by a port direction. -/
noncomputable def portPhase
    (a : ℝ) (k : Vec3) (u : Fin 12 → Vec3) (i : Fin 12) : ℝ :=
  a * dot k (u i)

/-- Exact twelve-port cosine symbol. -/
noncomputable def cosineSymbol
    (a : ℝ) (k : Vec3) (u : Fin 12 → Vec3) : ℝ :=
  (1 / (2 * a ^ 2)) *
    ∑ i : Fin 12, (1 - Real.cos (portPhase a k u i))

/-- Every nonzero-scale cosine symbol is nonnegative. -/
theorem cosineSymbol_nonnegative
    (a : ℝ) (ha : a ≠ 0) (k : Vec3) (u : Fin 12 → Vec3) :
    0 ≤ cosineSymbol a k u := by
  have hfactor : 0 ≤ (1 / (2 * a ^ 2) : ℝ) := by positivity
  have hsum : 0 ≤ ∑ i : Fin 12,
      (1 - Real.cos (portPhase a k u i)) := by
    exact Finset.sum_nonneg fun i _ => sub_nonneg.mpr (Real.cos_le_one _)
  exact mul_nonneg hfactor hsum

/-- Reversing the wavevector leaves the symbol fixed. -/
theorem cosineSymbol_wavevector_neg
    (a : ℝ) (k : Vec3) (u : Fin 12 → Vec3) :
    cosineSymbol a (-k) u = cosineSymbol a k u := by
  simp [cosineSymbol, portPhase, dot_neg_left, Real.cos_neg]

/-- Reversing every port direction leaves the symbol fixed. -/
theorem cosineSymbol_frame_neg
    (a : ℝ) (k : Vec3) (u : Fin 12 → Vec3) :
    cosineSymbol a k (fun i => -(u i)) = cosineSymbol a k u := by
  simp [cosineSymbol, portPhase, dot_neg_right, Real.cos_neg]

/-- The symbol depends on the complete port multiset, not its labels. -/
theorem cosineSymbol_reindex
    (a : ℝ) (k : Vec3) (u : Fin 12 → Vec3)
    (g : Equiv.Perm (Fin 12)) :
    cosineSymbol a k (fun i => u (g i)) = cosineSymbol a k u := by
  unfold cosineSymbol portPhase
  congr 1
  exact g.bijective.sum_comp
    (fun i : Fin 12 => 1 - Real.cos (a * dot k (u i)))

/-- A registered proper-carrier permutation preserves the complete sum by
relabeling.  This theorem alone is not an active spatial rotation. -/
theorem cosineSymbol_P60_covariant
    (a : ℝ) (k : Vec3) (u : Fin 12 → Vec3)
    (g : Equiv.Perm (Fin 12)) (_hg : g ∈ OPH.A5PortModule.P60) :
    cosineSymbol a k (fun i => u (g i)) = cosineSymbol a k u :=
  cosineSymbol_reindex a k u g

/-- A simultaneous orthogonal frame change and port relabeling preserves the
symbol. -/
theorem cosineSymbol_orthogonal_covariant
    (a : ℝ) (k : Vec3) (u : Fin 12 → Vec3)
    (R : Vec3 ≃ₗ[ℝ] Vec3) (g : Equiv.Perm (Fin 12))
    (hOrth : ∀ x y, dot (R x) (R y) = dot x y) :
    cosineSymbol a (R k) (fun i => R (u (g i))) = cosineSymbol a k u := by
  unfold cosineSymbol portPhase
  simp_rw [hOrth]
  congr 1
  exact g.bijective.sum_comp
    (fun i : Fin 12 => 1 - Real.cos (a * dot k (u i)))

/-- Conditional active carrier covariance when an orthogonal map realizes the
declared port permutation on the same frame. -/
theorem cosineSymbol_active_carrier_covariant
    (a : ℝ) (k : Vec3) (u : Fin 12 → Vec3)
    (R : Vec3 ≃ₗ[ℝ] Vec3) (g : Equiv.Perm (Fin 12))
    (hOrth : ∀ x y, dot (R x) (R y) = dot x y)
    (hR : ∀ i, R (u i) = u (g i)) :
    cosineSymbol a (R k) u = cosineSymbol a k u := by
  have hframe : (fun i => R (u (g.symm i))) = u := by
    funext i
    simpa using hR (g.symm i)
  have hcov := cosineSymbol_orthogonal_covariant a k u R g.symm hOrth
  rw [hframe] at hcov
  exact hcov

/-- The finite-range operator is diagonal on plane waves with the real cosine
symbol as eigenvalue. -/
theorem primitivePortOperator_planeWave
    (a : ℝ) (k x : Vec3) (u : Fin 12 → Vec3) :
    primitivePortOperator a u (planeWave k) x =
      (cosineSymbol a k u : ℂ) * planeWave k x := by
  unfold primitivePortOperator
  simp_rw [pairedDifference_planeWave]
  simp_rw [dot_smul_right]
  rw [← Finset.sum_mul]
  unfold cosineSymbol portPhase
  push_cast
  ring

/-! ## Conditional frozen moment packet -/

/-- Even directional port moment. -/
noncomputable def evenMoment
    (order : ℕ) (n : Vec3) (u : Fin 12 → Vec3) : ℝ :=
  ∑ i : Fin 12, dot n (u i) ^ order

/-- Squared radius in the declared Cartesian frame. -/
noncomputable def radialSquare (n : Vec3) : ℝ :=
  ∑ d : Fin 3, (n d) ^ 2

/-- Sixth-order moment expression for the normalized cosine symbol. -/
noncomputable def sixthOrderSymbol
    (a k m2 m4 m6 : ℝ) : ℝ :=
  (1 / (2 * a ^ 2)) *
    (((a * k) ^ 2 / 2) * m2 -
      ((a * k) ^ 4 / 24) * m4 +
      ((a * k) ^ 6 / 720) * m6)

/-- Degree-six directional moment truncation.  This definition contains no
claim about the analytic remainder of the exact cosine symbol. -/
noncomputable def portSymbolTaylorSix
    (a k : ℝ) (n : Vec3) (u : Fin 12 → Vec3) : ℝ :=
  sixthOrderSymbol a k
    (evenMoment 2 n u) (evenMoment 4 n u) (evenMoment 6 n u)

/-- Exact premises needed to specialize the moment truncation to the frozen
icosahedral coefficient packet along one declared unit direction. -/
structure FrozenMomentPacket
    (n : Vec3) (u : Fin 12 → Vec3) (i6 : ℝ) : Prop where
  unitDirection : radialSquare n = 1
  secondMoment : evenMoment 2 n u = 4
  fourthMoment : evenMoment 4 n u = 12 / 5
  sixthMoment : evenMoment 6 n u = 12 / 7 + (64 / 175) * i6

/-- The frozen moment values imply the exact coefficient polynomial. -/
theorem sixthOrderSymbol_frozen_expansion
    (a k i6 : ℝ) (ha : a ≠ 0) :
    sixthOrderSymbol a k 4 (12 / 5) (12 / 7 + (64 / 175) * i6) =
      k ^ 2 - (a ^ 2 / 20) * k ^ 4 +
        (a ^ 4 / 840) * k ^ 6 +
        ((2 * a ^ 4 / 7875) * k ^ 6) * i6 := by
  unfold sixthOrderSymbol
  field_simp [ha]
  ring

/-- A supplied exact moment packet gives the frozen sixth-order expansion. -/
theorem frozen_taylor_expansion_of_packet
    (a k i6 : ℝ) (n : Vec3) (u : Fin 12 → Vec3)
    (ha : a ≠ 0) (h : FrozenMomentPacket n u i6) :
    portSymbolTaylorSix a k n u =
      k ^ 2 - (a ^ 2 / 20) * k ^ 4 +
        (a ^ 4 / 840) * k ^ 6 +
        ((2 * a ^ 4 / 7875) * k ^ 6) * i6 := by
  unfold portSymbolTaylorSix
  rw [h.secondMoment, h.fourthMoment, h.sixthMoment]
  exact sixthOrderSymbol_frozen_expansion a k i6 ha

/-- The coefficient tuple agrees with the frozen arithmetic source. -/
theorem frozen_coefficient_tuple :
    OPH.A5PrimitivePortPrediction.C4 = -1 / 20 ∧
      OPH.A5PrimitivePortPrediction.B0 = 1 / 840 ∧
      OPH.A5PrimitivePortPrediction.B6 = 2 / 7875 := by
  norm_num [OPH.A5PrimitivePortPrediction.C4,
    OPH.A5PrimitivePortPrediction.B0,
    OPH.A5PrimitivePortPrediction.B6]

/-- Scale-free coefficient relations inherited from the frozen source. -/
theorem frozen_scale_free_relations :
    OPH.A5PrimitivePortPrediction.B0 /
        OPH.A5PrimitivePortPrediction.C4 ^ 2 = 10 / 21 ∧
      OPH.A5PrimitivePortPrediction.B6 /
        OPH.A5PrimitivePortPrediction.C4 ^ 2 = 32 / 315 ∧
      OPH.A5PrimitivePortPrediction.B6 /
        OPH.A5PrimitivePortPrediction.B0 = 16 / 75 := by
  exact ⟨OPH.A5PrimitivePortPrediction.b0_over_c4_squared,
    OPH.A5PrimitivePortPrediction.b6_over_c4_squared,
    OPH.A5PrimitivePortPrediction.b6_over_b0⟩

#print axioms shift_inverse
#print axioms pairedDifference_neg
#print axioms planeWave_shift
#print axioms pairedDifference_planeWave
#print axioms primitivePortOperator_reindex
#print axioms cosineSymbol_nonnegative
#print axioms cosineSymbol_wavevector_neg
#print axioms cosineSymbol_frame_neg
#print axioms cosineSymbol_reindex
#print axioms cosineSymbol_P60_covariant
#print axioms cosineSymbol_orthogonal_covariant
#print axioms cosineSymbol_active_carrier_covariant
#print axioms primitivePortOperator_planeWave
#print axioms sixthOrderSymbol_frozen_expansion
#print axioms frozen_taylor_expansion_of_packet
#print axioms frozen_coefficient_tuple
#print axioms frozen_scale_free_relations

end OPH.PrimitivePortTranslationBridge
