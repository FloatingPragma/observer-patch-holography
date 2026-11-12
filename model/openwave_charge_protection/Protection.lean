import Mathlib.Data.Int.Basic
import Lean.Elab.Tactic.Omega

/-! Integer charge identities conditional on an integer-valued shell degree.
The actual shell map and its homotopy are supplied by the Python certificate
and mathematical proof, not formalized by these arithmetic lemmas. -/

namespace OPH.WindingChargeProtection

theorem unit_degree_has_no_equal_halves (q : ℤ) : q + q ≠ 1 := by omega

theorem negative_unit_degree_has_no_equal_halves (q : ℤ) : q + q ≠ -1 := by omega

theorem integer_curl_step_preserves_sum (a b j : ℤ) :
    (a+j)+(b-j)=a+b := by omega

theorem same_degree_different_charge_forbids_dictionary
    (f : ℤ → ℤ) (w : ℤ) (h0 : f w = 0) (h1 : f w = 1) : False := by omega

theorem integer_winding_permits_unequal_signed_split : (2 : ℤ) + (-1) = 1 := by omega

#print axioms unit_degree_has_no_equal_halves
#print axioms negative_unit_degree_has_no_equal_halves
#print axioms integer_curl_step_preserves_sum
#print axioms same_degree_different_charge_forbids_dictionary

end OPH.WindingChargeProtection
