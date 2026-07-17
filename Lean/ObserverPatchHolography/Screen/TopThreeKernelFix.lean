import Mathlib

/-!
Kernel-evaluated replacement for the supplied `TopThree.lean` use of
`native_decide`.  The original theorem was true, but `native_decide` adds a
native-code axiom to the trust report.  Plain `decide` proves this six-element
fact in Lean's kernel.
-/

namespace OPH.TopThreeKernelFix

abbrev S3 := Equiv.Perm (Fin 3)

theorem s3_ambivalent : ∀ g : S3, ∃ x : S3, x * g * x⁻¹ = g⁻¹ := by
  decide

#print axioms s3_ambivalent

end OPH.TopThreeKernelFix
