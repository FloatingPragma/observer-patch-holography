import Lean.Elab.Tactic.Omega

/-! An order-three rotation of a twelve-port cell with no fixed port has four
three-port orbits. Each internal seam contributes one head. For a single
compensated mobile charge the conserved internal reference vector sums to
30 minus the charge. This lemma checks the integer obstruction after the
orbit and sum hypotheses are supplied by the geometry and contact checks. -/
namespace SourceLaw.WrittenChargeSymmetry

structure Orbit where
  first : Int
  second : Int
  third : Int

def Orbit.total (x : Orbit) : Int := x.first + x.second + x.third
def Orbit.fixedByTurn (x : Orbit) : Prop := x.first = x.second ∧ x.second = x.third

theorem fixedOrbitTotal (x : Orbit) (h : x.fixedByTurn) : x.total = 3 * x.first := by
  unfold Orbit.fixedByTurn at h
  unfold Orbit.total
  omega

theorem no_fixed_reference_for_single_charge (a b c d : Orbit) (q : Int)
    (charge : q = -1 ∨ q = 1)
    (gauss : a.total + b.total + c.total + d.total = 30 - q) :
    ¬ (a.fixedByTurn ∧ b.fixedByTurn ∧ c.fixedByTurn ∧ d.fixedByTurn) := by
  rintro ⟨ha, hb, hc, hd⟩
  rw [fixedOrbitTotal a ha, fixedOrbitTotal b hb, fixedOrbitTotal c hc,
      fixedOrbitTotal d hd] at gauss
  omega

end SourceLaw.WrittenChargeSymmetry
