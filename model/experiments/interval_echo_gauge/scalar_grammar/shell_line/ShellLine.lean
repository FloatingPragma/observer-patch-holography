import ElectroweakBreakingComposition

namespace OPH.ScalarShellLine

open OPH.ElectroweakBreakingComposition

/-- The chart uses a nonzero reference for every finite doublet. -/
def chartWeight (h : Doublet) : ℝ := 1 + normSqSum h

theorem chartWeight_pos (h : Doublet) : 0 < chartWeight h := by
  unfold chartWeight
  linarith [normSqSum_nonneg h]

/-- A pair of independent signed reads, divided by the two-reference frequency. -/
theorem product_reads_potential (lam r : ℝ) (h : Doublet) :
    lam * ((normSqSum h - r^2) / chartWeight h)^2 /
        (1 / chartWeight h)^2 = potential lam r h := by
  have hn := ne_of_gt (chartWeight_pos h)
  unfold potential
  field_simp

/-- Each complex component mean is its amplitude times its reference
frequency. The three nonzero frequencies cancel in the oriented multilinear
read, including when a component itself vanishes. -/
theorem quadratures_read_yukawa (y : ℝ) (h ψ : Doublet) (χ p q s : ℂ)
    (hp : p ≠ 0) (hq : q ≠ 0) (hs : s ≠ 0) :
    ((y : ℂ) * ((ψ 0 * p) * (h 1 * q) - (ψ 1 * p) * (h 0 * q)) * (χ * s)) /
        (p * q * s) = yukawaLine y h ψ χ := by
  unfold yukawaLine epsPair
  field_simp

/-- Applying the same change of weak frame to both factors multiplies the alternating
pairing by its determinant, which is one for the special-unitary block. -/
theorem epsilon_frame (a b c d : ℂ) (ψ h : Doublet) :
    (a * ψ 0 + b * ψ 1) * (c * h 0 + d * h 1) -
      (c * ψ 0 + d * ψ 1) * (a * h 0 + b * h 1) =
      (a*d - b*c) * epsPair ψ h := by
  unfold epsPair
  ring

theorem lower_vacuum_read (y r : ℝ) (ψ : Doublet) (χ p q s : ℂ)
    (hp : p ≠ 0) (hq : q ≠ 0) (hs : s ≠ 0) :
    ((y : ℂ) * ((ψ 0 * p) * (vacuum r 1 * q) - (ψ 1 * p) * (vacuum r 0 * q)) * (χ * s)) /
        (p * q * s) = ((y * r : ℝ) : ℂ) * (ψ 0 * χ) := by
  rw [quadratures_read_yukawa y (vacuum r) ψ χ p q s hp hq hs, yukawa_mass_at_vacuum]

end OPH.ScalarShellLine
