import EventAlgebra.NativeBinarySourceEndpoints
import EventAlgebra.NativeL2CountMoment

set_option maxRecDepth 4096
set_option maxHeartbeats 2000000

/-!
Activation and assumption witnesses for the limited source contract.
These witnesses separate the generic mathematical lemmas from the executable
66-state source/consumer census; they do not claim a kernel proof of that census.
-/

namespace EventAlgebra.NativeBinarySource.Witnesses

/-- A serialization witness, accepted only after kernel evaluation of the
actual source parser. It is not a replacement input for source decoding. -/
private def parsedProfile : Program where
  difference := .bin "-" .loadA .loadB
  current := .bin "//" (.bin "-" (.var "d")
    (.bin "*" (.bin "-" (.bin "*" (.lit 2) (.var "coin")) (.lit 1))
      (.bin "%" (.abs (.var "d")) (.lit 2)))) (.lit 2)
  subtractA := true
  updateA := .var "current"
  subtractB := false
  updateB := .var "current"
  recordCoin := .var "coin"
  recordCurrent := .var "current"

theorem original_parse_shape : parseProgram Endpoints.sourceText = some parsedProfile := by decide

/-- The coin guard is a binary-carrier/protocol input condition. It is not
needed for the bare endpoint placement equation on binary endpoint loads. -/
theorem source_binary_without_coin_guard (x y coin : Int) (tail : List Int)
    (hx : x=0 ∨ x=1) (hy : y=0 ∨ y=1) :
    Endpoints.runSource (x::y::tail) 0 1 coin =
      some (Endpoints.endpointA x y coin :: Endpoints.endpointB x y coin :: tail,
            (coin,x-Endpoints.endpointA x y coin)) := by
  have hm : (-1-(2*coin-1))/2 = -coin := by omega
  have hp : (1-(2*coin-1))/2 = 1-coin := by omega
  rcases hx with rfl | rfl <;> rcases hy with rfl | rfl <;>
    simp [Endpoints.runSource,original_parse_shape,parsedProfile,evalProgram,evalExpr,
      Endpoints.endpointA,Endpoints.endpointB,hm,hp] <;> simp only [sub_eq_add_neg]

theorem source_preconditions_inhabited :
    ∃ x y coin : Int, (x=0 ∨ x=1) ∧ (y=0 ∨ y=1) ∧ (coin=0 ∨ coin=1) := by
  exact ⟨0,1,1,Or.inl rfl,Or.inr rfl,Or.inr rfl⟩

theorem signed_source_activation :
    Endpoints.runSource [0,1,7,-11] 0 1 1 = some ([1,0,7,-11],(1,-1)) := by
  simpa [Endpoints.endpointA,Endpoints.endpointB] using
    Endpoints.source_binary_endpoint 0 1 1 [7,-11]
      (Or.inl rfl) (Or.inr rfl) (Or.inr rfl)

theorem binary_carrier_activation :
    Endpoints.runSource [0,1,0,0,0,0,0,0,0,0,0,1] 0 1 1 =
      some ([1,0,0,0,0,0,0,0,0,0,0,1],(1,-1)) := by
  simpa [Endpoints.endpointA,Endpoints.endpointB] using
    Endpoints.source_binary_endpoint 0 1 1 [0,0,0,0,0,0,0,0,0,1]
      (Or.inl rfl) (Or.inr rfl) (Or.inr rfl)

theorem collision_is_not_the_binary_placement_law :
    Endpoints.runSource [2,0] 0 1 0 ≠ some ([0,1],(0,2)) := by decide

theorem collision_at_second_endpoint :
    Endpoints.runSource [0,2] 0 1 0 ≠ some ([0,1],(0,0)) := by decide

theorem outside_coin_leaves_binary_carrier :
    Endpoints.runSource [1,0] 0 1 2 = some ([2,-1],(2,-1)) := by decide

def gramWitness : Matrix (Fin 1) (Fin 1) ℝ := fun _ _ => 2

theorem gram_positive_activation :
    ((4 : ℝ)⁻¹ • (gramWitness*gramWitness)).PosSemidef := by
  exact L2.normalized_square_positive gramWitness rfl 4 (by norm_num)

theorem gram_trace_activation :
    ((4 : ℝ)⁻¹ • (gramWitness*gramWitness)).trace = 1 := by
  apply L2.normalized_trace gramWitness 4 (by norm_num)
  norm_num [gramWitness,Matrix.trace,Matrix.mul_apply]

theorem ratio_activation :
    ((1/60 : ℚ)*104/(66*60^2)) / ((1/60)*(171720/(66*60^2))) = 13/21465 := by
  have h := L2.count_first_ratio 104 171720 66 60 (1/60)
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)
  calc
    ((1/60 : ℚ)*104/(66*60^2)) / ((1/60)*(171720/(66*60^2))) = 104/171720 := h
    _ = 13/21465 := by norm_num

/-- Algebraic equality also holds at total=0. The nonzero total guard is kept
in the public helper as a stated normalization-domain condition, not advertised
as necessary for this bare field identity. The actual source return has Z>0. -/
theorem ratio_without_total_guard (tail total population denominator coefficient : ℚ)
    (hPopulation : population ≠ 0) (hDenominator : denominator ≠ 0)
    (hCoefficient : coefficient ≠ 0) :
    (coefficient*tail/(population*denominator^2)) /
      (coefficient*(total/(population*denominator^2))) = tail/total := by
  by_cases hTotal : total=0
  · simp [hTotal]
  · exact L2.count_first_ratio tail total population denominator coefficient
      hTotal hPopulation hDenominator hCoefficient

theorem zero_coefficient_failure :
    ((0 : ℚ)*104/(66*60^2)) / (0*(171720/(66*60^2))) ≠ 13/21465 := by norm_num

theorem zero_population_failure :
    ((1/60 : ℚ)*104/(0*60^2)) / ((1/60)*(171720/(0*60^2))) ≠ 13/21465 := by norm_num

theorem zero_attempt_denominator_failure :
    ((1/60 : ℚ)*104/(66*0^2)) / ((1/60)*(171720/(66*0^2))) ≠ 13/21465 := by norm_num

theorem zero_normalization_failure :
    ((0 : ℝ) • (gramWitness*gramWitness)).trace ≠ 1 := by simp

def nonsymmetricWitness : Matrix (Fin 2) (Fin 2) ℝ := !![1,1;0,1]

theorem nonsymmetric_square_failure :
    ¬ (((2 : ℝ)⁻¹ • (nonsymmetricWitness*nonsymmetricWitness)).PosSemidef) := by
  intro h
  have equality := congrFun (congrFun h.1 0) 1
  norm_num [nonsymmetricWitness,Matrix.mul_apply,Matrix.conjTranspose_apply,
    Fin.sum_univ_two] at equality

#print axioms original_parse_shape
#print axioms source_binary_without_coin_guard
#print axioms source_preconditions_inhabited
#print axioms signed_source_activation
#print axioms binary_carrier_activation
#print axioms collision_is_not_the_binary_placement_law
#print axioms collision_at_second_endpoint
#print axioms outside_coin_leaves_binary_carrier
#print axioms gram_positive_activation
#print axioms gram_trace_activation
#print axioms ratio_activation
#print axioms ratio_without_total_guard
#print axioms zero_coefficient_failure
#print axioms zero_population_failure
#print axioms zero_attempt_denominator_failure
#print axioms zero_normalization_failure
#print axioms nonsymmetric_square_failure

end EventAlgebra.NativeBinarySource.Witnesses
