import ObserverPatchHolography.EinsteinBranch.QuantitativeSmallBall
import ObserverPatchHolography.EinsteinBranch.FixedCapacityWLaw

/-! Exact finite-readout and conditional scalar transfer identities.
The imported modules own the regional tensor bound and capacity w-law.
These results do not construct the geometric premises or a cosmic capacity. -/

namespace OPH.GravityCapacityTransfer

theorem lapse_area_response (p q : ℝ) :
    ((1-q)*(9/10) + q*(11/10)) - ((1-p)*(9/10) + p*(11/10)) =
      (-99/200) * (((1-q)*(20/9) + q*(20/11)) -
                   ((1-p)*(20/9) + p*(20/11))) := by
  ring

theorem poisson_source_error_bound {e0 e1 e2 e3 epsilon : ℝ}
    (h0 : |e0| ≤ epsilon) (h1 : |e1| ≤ epsilon)
    (h2 : |e2| ≤ epsilon) (h3 : |e3| ≤ epsilon) :
    |(e0+e1+e2+e3)/2| ≤ 2*epsilon := by
  obtain ⟨l0, u0⟩ := abs_le.mp h0
  obtain ⟨l1, u1⟩ := abs_le.mp h1
  obtain ⟨l2, u2⟩ := abs_le.mp h2
  obtain ⟨l3, u3⟩ := abs_le.mp h3
  apply abs_le.mpr
  constructor <;> linarith

theorem capacity_epoch_invariant (N kappa : ℝ) (hN : N ≠ 0) :
    N*(kappa/N)=kappa := by
  field_simp

theorem density_scale_invariant (N kappa c : ℝ) (hN : N ≠ 0) (hc : c ≠ 0) :
    (c*kappa)/(c*N)=kappa/N := by
  field_simp

#print axioms lapse_area_response
#print axioms poisson_source_error_bound
#print axioms capacity_epoch_invariant
#print axioms density_scale_invariant

end OPH.GravityCapacityTransfer
