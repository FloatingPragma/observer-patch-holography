import SeamCurrentAuxiliaryOscillatorLift
import Mathlib.Tactic

/-! A signed incidence lift realizes the completion oscillator on each Fourier fibre.
The second phase coordinate is the divergence of the link amplitudes. -/
namespace OPH.TransverseIncidence

open OPH.SeamCurrentDirichletGenerator OPH.SeamCurrentHomogeneousAction
open scoped BigOperators

noncomputable def squareWeight (x : ℝ) := ((1-Real.cos x)^2 + (Real.sin x)^2)/60

theorem weight_eq (x : ℝ) : squareWeight x = (1-Real.cos x)/30 := by
  unfold squareWeight
  nlinarith [Real.sin_sq_add_cos_sq x]

noncomputable def edgePhase (k : OPH.SeamCurrentDirichletGenerator.Vec3) (e : Fin 30) :=
  completionStepScale * OPH.PrimitivePortTranslationBridge.dot k (unitCarrierSeamDirection e)

theorem incidence_symbol (k : OPH.SeamCurrentDirichletGenerator.Vec3) :
    ∑ e, squareWeight (edgePhase k e) = completionFourierSymbol k := by
  simp_rw [weight_eq]
  rw [← Finset.sum_div, Finset.sum_sub_distrib]
  simp [completionFourierSymbol, edgePhase]
  ring

/-- Link amplitude in the phase convention j = i times the stored amplitude. -/
noncomputable def generator (b : Fin 30 → ℂ) (s : ℂ × (Fin 30 → ℂ)) :
    ℂ × (Fin 30 → ℂ) := (-(∑ e, (star (b e)) * s.2 e), fun e => b e * s.1)

theorem second_phase (b : Fin 30 → ℂ) (s : ℂ × (Fin 30 → ℂ)) :
    (generator b (generator b s)).1 = -(∑ e, star (b e) * b e) * s.1 := by
  simp [generator, Finset.sum_mul, mul_assoc]

/-- The real and imaginary parts of a signed difference with equal contact
weight 1/sqrt(60) give this squared magnitude. -/
noncomputable def edgeAmplitude (x : ℝ) : ℂ :=
  ⟨(1-Real.cos x)/Real.sqrt 60, -Real.sin x/Real.sqrt 60⟩

theorem edge_norm (x : ℝ) :
    star (edgeAmplitude x) * edgeAmplitude x = (squareWeight x : ℂ) := by
  change (starRingEnd ℂ) _ * _ = _
  rw [mul_comm, Complex.mul_conj]
  apply congrArg Complex.ofReal
  rw [edgeAmplitude, Complex.normSq_mk, div_mul_div_comm, div_mul_div_comm,
    Real.mul_self_sqrt (by norm_num : (0:ℝ) ≤ 60)]
  unfold squareWeight
  ring

theorem transverse_oscillator (k : OPH.SeamCurrentDirichletGenerator.Vec3) (s : ℂ × (Fin 30 → ℂ)) :
    (generator (fun e => edgeAmplitude (edgePhase k e))
      (generator (fun e => edgeAmplitude (edgePhase k e)) s)).1 =
      -(completionFourierSymbol k : ℂ) * s.1 := by
  rw [second_phase]
  simp_rw [edge_norm]
  rw [← Complex.ofReal_sum, incidence_symbol]

/-- Restore the spatial character to the vertex reading. The second temporal
derivative cancels the *existing completion Dirichlet operator* at every
carrier point, for any initial vertex and link amplitudes. -/
theorem vertex_wave_equation (k x : OPH.SeamCurrentDirichletGenerator.Vec3)
    (s : ℂ × (Fin 30 → ℂ)) :
    (generator (fun e => edgeAmplitude (edgePhase k e))
      (generator (fun e => edgeAmplitude (edgePhase k e)) s)).1 *
        OPH.PrimitivePortTranslationBridge.planeWave k x +
      s.1 * completionComplexDirichletGenerator
        (OPH.PrimitivePortTranslationBridge.planeWave k) x = 0 := by
  rw [transverse_oscillator, completionComplexDirichletGenerator_planeWave]
  ring

theorem two_independent_layers (a b : ℝ) :
    (1-a)*(1-b) + a*(1-b) + b = 1 := by ring

/-- Unconditional detection after two independent opportunities with a
direction-independent survival s for each encountered layer. -/
theorem two_layer_survival (a b s : ℝ) :
    (1-a)*(1-b) + a*(1-b)*s + (1-a)*b*s + a*b*s^2 =
      (1-a+a*s)*(1-b+b*s) := by ring

theorem measured_visibility :
    let a : ℝ := 1/4
    let b : ℝ := 3/4
    let s : ℝ := 2/3
    let clear := (1-a)*(1-b)
    let first := a*(1-b)*s
    let last := (1-a)*b*s + a*b*s^2
    clear = 3/16 ∧ first = 1/24 ∧ last = 11/24 ∧
      clear/(clear+first+last) = 3/11 := by norm_num

#print axioms measured_visibility
#print axioms transverse_oscillator
#print axioms vertex_wave_equation
end OPH.TransverseIncidence
