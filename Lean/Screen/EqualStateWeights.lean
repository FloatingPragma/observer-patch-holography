import Mathlib
import A5PortAction

namespace OPH.EqualStateWeights

open OPH.A5PortAction (perms app)

/-! # Equal state weights from invariance, uniqueness, and transitivity

Machine-checked spine for issue #610. The realized A3 state assigns weight
`1/12` to every port projection whenever

1. the deck action preserves the feasible family and the A3 objective
   (the A3 exactness clauses: reference family, observer cover, and weights
   are constructed from quotient-visible A1 data by a rule natural under
   admissible presentation equivalence),
2. the information projection is unique on the feasible family, and
3. the action is transitive on ports.

The theorems here are exactly those three steps. Strict convexity of the
weighted divergence, which supplies hypothesis 2, and the concrete
invariance measurements, which supply hypothesis 1, live in the executable
certificate `code/a5_closure/equal_state_weights_certificate.py`. The
pairwise transitivity of the listed sixty rotations is checked by kernel
`decide` against `A5PortAction.perms`.

BOUNDARY. Nothing here selects the reference family or proves equal
central-block algebra dimensions; the block-trace value measured by the
carrier certificates is a distinct realization property. Dropping any
hypothesis admits the countermodels recorded in the certificate. -/

/-- Pairwise transitivity of the listed sixty rotations: every ordered port
pair is connected by a listed deck rotation. -/
theorem pairwise_transitive_on_ports :
    ((List.range 12).all fun s => (List.range 12).all fun t =>
      perms.any fun p => app p s == t) = true := by decide

/-- Invariance plus uniqueness forces a symmetric minimizer: if the feasible
set and objective are invariant under a set of port maps whose action is
transitive, then the unique minimizer is constant across ports. The maps
need not be listed as bijections; invariance and transitivity carry the
whole argument. -/
theorem unique_invariant_minimizer_constant
    (K : Set ((Fin 12) → ℝ)) (F : ((Fin 12) → ℝ) → ℝ)
    (G : Set ((Fin 12) → (Fin 12)))
    (hK : ∀ g ∈ G, ∀ x ∈ K, (x ∘ g) ∈ K)
    (hF : ∀ g ∈ G, ∀ x ∈ K, F (x ∘ g) = F x)
    (xstar : (Fin 12) → ℝ)
    (hmem : xstar ∈ K)
    (hmin : ∀ y ∈ K, F xstar ≤ F y)
    (huniq : ∀ y ∈ K, (∀ z ∈ K, F y ≤ F z) → y = xstar)
    (htrans : ∀ p q : Fin 12, ∃ g ∈ G, g p = q) :
    ∀ p q : Fin 12, xstar p = xstar q := by
  intro p q
  obtain ⟨g, hgG, hgpq⟩ := htrans p q
  have hKg : (xstar ∘ g) ∈ K := hK g hgG xstar hmem
  have hFg : F (xstar ∘ g) = F xstar := hF g hgG xstar hmem
  have hming : ∀ z ∈ K, F (xstar ∘ g) ≤ F z := by
    intro z hz
    calc F (xstar ∘ g) = F xstar := hFg
      _ ≤ F z := hmin z hz
  have heq : (xstar ∘ g) = xstar := huniq (xstar ∘ g) hKg hming
  have hp := congrFun heq p
  simp only [Function.comp_apply, hgpq] at hp
  exact hp.symm ▸ rfl

/-- A constant normalized port weight vector assigns exactly `1/12`
everywhere. -/
theorem constant_normalized_is_one_twelfth
    (x : (Fin 12) → ℝ)
    (hconst : ∀ p q : Fin 12, x p = x q)
    (hsum : (Finset.univ.sum x) = 1) :
    ∀ p : Fin 12, x p = 1 / 12 := by
  intro p
  have hall : ∀ q ∈ (Finset.univ : Finset (Fin 12)), x q = x p := by
    intro q _
    exact hconst q p
  have hsum' : (Finset.univ.sum x) = 12 * x p := by
    rw [Finset.sum_congr rfl hall]
    simp [Finset.sum_const, Finset.card_univ]
  have : (12 : ℝ) * x p = 1 := by rw [← hsum'] ; exact hsum
  linarith

/-- The composed statement for issue #610: invariance of the feasible family
and objective under a transitive port action, plus uniqueness of the
information projection and normalization of the port weights, force the
realized state weight `1/12` at every port. -/
theorem equal_state_weights
    (K : Set ((Fin 12) → ℝ)) (F : ((Fin 12) → ℝ) → ℝ)
    (G : Set ((Fin 12) → (Fin 12)))
    (hK : ∀ g ∈ G, ∀ x ∈ K, (x ∘ g) ∈ K)
    (hF : ∀ g ∈ G, ∀ x ∈ K, F (x ∘ g) = F x)
    (xstar : (Fin 12) → ℝ)
    (hmem : xstar ∈ K)
    (hmin : ∀ y ∈ K, F xstar ≤ F y)
    (huniq : ∀ y ∈ K, (∀ z ∈ K, F y ≤ F z) → y = xstar)
    (htrans : ∀ p q : Fin 12, ∃ g ∈ G, g p = q)
    (hsum : (Finset.univ.sum xstar) = 1) :
    ∀ p : Fin 12, xstar p = 1 / 12 :=
  constant_normalized_is_one_twelfth xstar
    (unique_invariant_minimizer_constant K F G hK hF xstar hmem hmin huniq
      htrans)
    hsum

end OPH.EqualStateWeights
