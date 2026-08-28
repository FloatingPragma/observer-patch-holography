import OphGap.OphData

/-!
# The OPH frozen visible seam complex satisfies the certificate

`GapScout.Oph.edges` is the 11,816-edge signed seam graph on 8,662 carriers
(freeze sha256 `a0be6fc6…`, BFS-relabelled), and `GapScout.Oph.cycle` is the receipt's
11-edge negative-cycle witness. `check_true` is verified by the Lean KERNEL
(`decide +kernel`), so `oph_gram_posDef` is a fully kernel-checked
theorem about this concrete graph.
-/

open Matrix

namespace GapScout.Oph

theorem n_eq : n = 8662 := rfl
theorem edges_length : edges.length = 11816 := by decide +kernel

/-- Every seam sign is `-1` (the repository's uniform reversing convention). -/
theorem all_negative : edges.all (fun e => e.2.2) = true := by decide +kernel

set_option maxRecDepth 100000 in
/-- THE CERTIFICATE PASSES, kernel-checked. -/
theorem check_true : check n edges cycle = true := by decide +kernel

/-- Main result: the signed Laplacian `Dᵀ D` of the frozen OPH graph is positive definite
(`unfold GramPosDef` for the literal statement). -/
theorem oph_gram_posDef : GramPosDef n edges cycle check_true :=
  gramPosDef_of_check check_true

/-- Every eigenvalue is strictly positive: the spectral gap (least eigenvalue) is positive. -/
theorem oph_eigenvalues_pos : EigenvaluesPos n edges cycle check_true :=
  eigenvaluesPos_of_check check_true

/-- Kernel of the signed derivative is trivial. -/
theorem oph_kernel_trivial : KernelTrivial n edges cycle check_true :=
  kernelTrivial_of_check' check_true

/-! ## NEGATIVE CONTROL: flip ONE sign on a witness edge. The same certificate must be refused. -/

def tampered : List Edge := flipAt tamperIndex edges

set_option maxRecDepth 100000 in
/-- The tampered certificate is REFUSED, kernel-checked. -/
theorem tampered_refused : check n tampered cycle = false := by decide +kernel

end GapScout.Oph
