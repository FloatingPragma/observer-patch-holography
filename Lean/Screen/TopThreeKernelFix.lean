import Mathlib

/-!
Kernel-checked S3 ambivalence. Plain `decide` proves this six-element fact
without a native-code axiom in the trust report.
-/

namespace OPH.TopThreeKernelFix

abbrev S3 := Equiv.Perm (Fin 3)

theorem s3_ambivalent : ∀ g : S3, ∃ x : S3, x * g * x⁻¹ = g⁻¹ := by
  decide

#print axioms s3_ambivalent

end OPH.TopThreeKernelFix
