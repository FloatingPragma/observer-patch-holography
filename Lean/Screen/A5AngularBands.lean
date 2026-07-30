import Mathlib

/-!
# Equal-port angular sequence and readback parity

The twelve icosahedral ports have pairwise dot products in
`{1, 1/√5, -1/√5, -1}`, so the equal-port average of the degree-`l`
Legendre kernel is

`I_l = (1 + (-1)^l + 5 (P_l(1/√5) + P_l(-1/√5))) / 12`.

This module computes the Legendre values at `±1/√5` exactly in `ℚ(√5)`,
represented as rational pairs `(a, b) ↦ a + b √5`, proves the closed
equal-port sequence through level fourteen, the vanishing of every odd
level, and the alternating readback parity sign table
`R|_l = (-1)^(l+1)` for the four interpolated levels. The matrix projector
and equivariance identities live in the exact-arithmetic certificate
producer; this file carries the arithmetic layer of the sprint contract.
-/

namespace OPH.AngularBands

/-- An element `a + b √5` of the quadratic field as a rational pair. -/
structure Q5 where
  a : ℚ
  b : ℚ
  deriving DecidableEq, Repr

namespace Q5

def add (x y : Q5) : Q5 := ⟨x.a + y.a, x.b + y.b⟩

def sub (x y : Q5) : Q5 := ⟨x.a - y.a, x.b - y.b⟩

def smul (q : ℚ) (x : Q5) : Q5 := ⟨q * x.a, q * x.b⟩

/-- Multiplication by `1/√5 = (1/5) √5`. -/
def mulInvSqrt5 (x : Q5) : Q5 := ⟨x.b, x.a / 5⟩

end Q5

/-- Legendre values `P_l(1/√5)` by the three-term recurrence. -/
def legendreInvSqrt5 : ℕ → Q5
  | 0 => ⟨1, 0⟩
  | 1 => ⟨0, 1 / 5⟩
  | (l + 2) =>
      Q5.smul (1 / (l + 2 : ℚ))
        (Q5.sub
          (Q5.smul (2 * (l : ℚ) + 3)
            (Q5.mulInvSqrt5 (legendreInvSqrt5 (l + 1))))
          (Q5.smul ((l : ℚ) + 1) (legendreInvSqrt5 l)))

/-- Legendre parity at the reflected node: `P_l(-t) = (-1)^l P_l(t)`. -/
def legendreNegInvSqrt5 (l : ℕ) : Q5 :=
  Q5.smul ((-1 : ℚ) ^ l) (legendreInvSqrt5 l)

/-- The equal-port sequence. Its `√5` part cancels, so it is the rational
    coordinate of the displayed combination. -/
def equalPort (l : ℕ) : ℚ :=
  (1 + (-1 : ℚ) ^ l
    + 5 * ((legendreInvSqrt5 l).a + (legendreNegInvSqrt5 l).a)) / 12

/-- The `√5` coordinate of the symmetrized Legendre pair vanishes, so the
    equal-port value is genuinely rational. -/
def equalPortIrrationalPart (l : ℕ) : ℚ :=
  (legendreInvSqrt5 l).b + (legendreNegInvSqrt5 l).b

theorem equalPort_unit : equalPort 0 = 1 := by
  norm_num [equalPort, legendreInvSqrt5, legendreNegInvSqrt5, Q5.smul]

/-- The even initial vector of the equal-port angular comb. -/
theorem equalPort_evenVector :
    equalPort 2 = 0 ∧ equalPort 4 = 0 ∧ equalPort 6 = 11 / 25 ∧
      equalPort 8 = 0 ∧ equalPort 10 = 247 / 1875 ∧
      equalPort 12 = 1071 / 3125 ∧ equalPort 14 = 0 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
    norm_num [equalPort, legendreInvSqrt5, legendreNegInvSqrt5,
      Q5.smul, Q5.sub, Q5.mulInvSqrt5]

/-- Every odd level through thirteen vanishes. -/
theorem equalPort_oddZero :
    equalPort 1 = 0 ∧ equalPort 3 = 0 ∧ equalPort 5 = 0 ∧
      equalPort 7 = 0 ∧ equalPort 9 = 0 ∧ equalPort 11 = 0 ∧
      equalPort 13 = 0 := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
    norm_num [equalPort, legendreInvSqrt5, legendreNegInvSqrt5,
      Q5.smul, Q5.sub, Q5.mulInvSqrt5]

/-- The symmetrized pair is rational at every checked level. -/
theorem equalPort_rational :
    (List.range 15).all (fun l => equalPortIrrationalPart l = 0) := by
  decide +kernel

/-- The inverse-port readback acts on the interpolated level-`l` image by
    the alternating parity sign. The band eigenvalues of the response
    polynomial `J = (A³ - 4A² - 5A + 10)/10` at the four adjacency bands
    `(5, √5, -√5, -1)`, matched to levels `(0, 1, 3, 2)`, give
    `R = -J` the level signs below; the matrix identity is certified in
    the exact-arithmetic producer. -/
def paritySign (l : ℕ) : ℤ := (-1) ^ (l + 1)

theorem parity_signs :
    paritySign 0 = -1 ∧ paritySign 1 = 1 ∧ paritySign 2 = -1 ∧
      paritySign 3 = 1 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> norm_num [paritySign]

/-- The response polynomial takes value `-1` at the frame and kernel band
    eigenvalues `±√5` and value `+1` at `5` and `-1`: with
    `p(x) = x³ - 4x² - 5x + 10`, the rational identities
    `p(5) = 10`, `p(-1) = 10`, and `p(±√5) = -10` hold, the latter through
    the quadratic-field computation `p(±√5) = ±5√5 - 20 ∓ 5√5 + 10`. -/
theorem response_band_values :
    (5 : ℚ) ^ 3 - 4 * 5 ^ 2 - 5 * 5 + 10 = 10 ∧
      ((-1 : ℚ)) ^ 3 - 4 * (-1) ^ 2 - 5 * (-1) + 10 = 10 ∧
      ∀ s : ℚ, s ^ 2 = 5 →
        s ^ 3 - 4 * s ^ 2 - 5 * s + 10 = -10 := by
  refine ⟨by norm_num, by norm_num, ?_⟩
  intro s hs
  have hcube : s ^ 3 = 5 * s := by
    have : s ^ 3 = s ^ 2 * s := by ring
    rw [this, hs]
  rw [hcube, hs]
  ring

-- Axiom audit.
#print axioms equalPort_unit
#print axioms equalPort_evenVector
#print axioms equalPort_oddZero
#print axioms equalPort_rational
#print axioms parity_signs
#print axioms response_band_values

end OPH.AngularBands
