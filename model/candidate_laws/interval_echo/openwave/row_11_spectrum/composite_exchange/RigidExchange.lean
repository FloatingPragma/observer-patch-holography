import Mathlib.Algebra.BigOperators.Group.List.Basic
import Mathlib.Algebra.Group.Defs
import Mathlib.Tactic.Ring

namespace OPH.RigidCompositeExchange

/-- At half-step i the two copies traverse complementary halves of the path.
Their product accumulates exactly the transport phase of one full winding.
Each entry may itself include all port, species and orientation factors of a
rigid internal packet. -/
theorem matched_winding {K : Type*} [CommMonoid K] (steps : List (K × K)) :
    (steps.map (fun p => p.1 * p.2)).prod =
    (steps.map Prod.fst).prod * (steps.map Prod.snd).prod := by
  induction steps with
  | nil => simp
  | cons p tail ih =>
    simp only [List.map_cons, List.prod_cons, ih]
    ac_rfl

/-- An ungraded block flip is identity on two copies of any common
internal pure state with amplitudes u, including reflection-odd states. -/
theorem identical_packet_flip {A K : Type*} [CommMonoid K]
    (u : A → K) (left right : A) :
    u right * u left = u left * u right := mul_comm _ _

/-- A reflection of just one odd packet produces a minus in the same
interferometer, although swapping the whole copies preserves the product. -/
theorem separate_odd_reflection {A K : Type*} [CommRing K]
    (r : A → A) (u : A → K) (odd : ∀ a, u (r a) = -(u a))
    (left right : A) :
    u (r right) * u left = -(u left * u right) := by
  rw [odd, neg_mul, mul_comm (u right)]

end OPH.RigidCompositeExchange
