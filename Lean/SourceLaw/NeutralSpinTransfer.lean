import Mathlib

/-! The common donor supplies the same phase on every converted charged
address. These identities apply to a monitored neutral conversion whose
ancillary records are identical for all addresses. The independently tested
pulse and readout determine whether a physical instrument obeys the condition. -/
namespace SourceLaw.NeutralSpinTransfer

/-- A common unit-modulus conversion amplitude preserves every coherence,
including when an address index carries an unread reference coordinate. -/
theorem preserve_matrix_element (c x y : ℂ) (hc : c * star c = 1) :
    (c*x) * star (c*y) = x * star y := by
  rw [star_mul]
  calc
    (c*x)*(star y*star c) = (c*star c)*(x*star y) := by ring
    _ = x*star y := by rw [hc, one_mul]

/-- The compensated magnetic turn can be applied after charge conversion.
Its sign relative to an oppositely ordered turn transfers to the neutral
code whenever the two turns act on its surviving vacancy coordinate. -/
theorem transfer_order {V : Type*} [AddCommGroup V] [Module ℂ V]
    (r s : V →ₗ[ℂ] V) (v : V) (phase : ℂ)
    (h : r (s v) = -(s (r v))) :
    r (s (phase • v)) = -(s (r (phase • v))) := by
  simp only [map_smul, h, smul_neg]

/-- With a uniformly spread rank-two code over twelve ports, the
port-dephasing entanglement-fidelity expression equals one twelfth. -/
theorem uniform_dephasing_read (p : Fin 12 → ℝ)
    (h : ∀ j, p j = 1/6) : (∑ j, (p j)^2)/4 = 1/12 := by
  simp_rw [h]
  norm_num

#print axioms preserve_matrix_element
#print axioms transfer_order
#print axioms uniform_dephasing_read

end SourceLaw.NeutralSpinTransfer
