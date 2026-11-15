import Geometry.SourceConstrainedSelection
import Geometry.SourceReadSelection

/-!
# Constrained A3 schedule to native read bridge

The feasible family may contain correlated schedules and the reference need
not be uniform or feasible. Its scalar weights describe a finite history law
when the source grammar supplies normalization. The theorem below needs only
nonnegativity, convexity, faithfulness and an attained KL minimum. The full
history atoms are scored here; a coarse cover needs the separate positive
failure-readout representation theorem.
-/

set_option autoImplicit false

namespace OPH.SourceConstrainedRead
noncomputable section
open OPH.SourceConstrainedSelection OPH.SourceReadSelection
open ObserverPatchHolography.RepairWordSchedule

variable {Move ι R Value : Type*} [Fintype Move] [DecidableEq ι]

/-- A single feasible positive-weight cut-avoiding word survives constrained
information projection and prevents a guaranteed correct two-input read.
No IID, uniform-reference or full-simplex restriction on the finite word is
assumed. -/
theorem constrained_remote_read_obstruction (k : ℕ)
    (K : Set (Word Move k → ℝ)) (hK : Convex ℝ K)
    (hKnonneg : ∀ p ∈ K, ∀ w, 0 ≤ p w)
    (reference selected : Word Move k → ℝ) (hr : ∀ w, 0 < reference w)
    (hp : selected ∈ K)
    (hmin : ∀ q ∈ K, weightedKL (fun _ => 1) reference selected ≤
      weightedKL (fun _ => 1) reference q)
    (q : Word Move k → ℝ) (hq : q ∈ K) (w : Word Move k) (hqw : 0 < q w)
    (edge : Move → ι × ι) (source : ι → Prop)
    (hcut : ∀ i, NoCross source (edge (w i)))
    (readPort : R → ι) (hread : ∀ r, ¬ source (readPort r))
    (x y : ι → ℝ) (hxy : ∀ i, ¬ source i → x i=y i)
    (decode : Word Move k → (ℕ → R → ℝ) → Value)
    (a b : Value) (hab : a ≠ b) :
    0 < selected w ∧
      ¬ (decode w (transcript (List.ofFn (edge ∘ w)) x readPort)=a ∧
         decode w (transcript (List.ofFn (edge ∘ w)) y readPort)=b) := by
  refine ⟨weightedKL_minimizer_support K hK hKnonneg _ reference selected
    (by intro v; norm_num) hr hp hmin q hq w hqw, ?_⟩
  apply cut_transcript_not_both_correct source readPort hread _ _ x y hxy (decode w) a b hab
  intro e he
  obtain ⟨i,rfl⟩ := List.mem_ofFn.mp he
  exact hcut i

/-- Guaranteed remote correctness forces the *entire* convex feasible
family to exclude every cut-avoiding word. A favorable optimizer cannot
replace that grammar-level obligation. This is only a necessary condition. -/
theorem guaranteed_read_excludes_feasible_cut_words (k : ℕ)
    (K : Set (Word Move k → ℝ)) (hK : Convex ℝ K)
    (hKnonneg : ∀ p ∈ K, ∀ w, 0 ≤ p w)
    (reference selected : Word Move k → ℝ) (hr : ∀ w, 0 < reference w)
    (hp : selected ∈ K)
    (hmin : ∀ q ∈ K, weightedKL (fun _ => 1) reference selected ≤
      weightedKL (fun _ => 1) reference q)
    (edge : Move → ι × ι) (source : ι → Prop)
    (readPort : R → ι) (hread : ∀ r, ¬ source (readPort r))
    (x y : ι → ℝ) (hxy : ∀ i, ¬ source i → x i=y i)
    (decode : Word Move k → (ℕ → R → ℝ) → Value)
    (a b : Value) (hab : a ≠ b)
    (correct : ∀ w, 0 < selected w →
      decode w (transcript (List.ofFn (edge ∘ w)) x readPort)=a ∧
      decode w (transcript (List.ofFn (edge ∘ w)) y readPort)=b)
    (q : Word Move k → ℝ) (hq : q ∈ K) (w : Word Move k)
    (hcut : ∀ i, NoCross source (edge (w i))) : q w=0 := by
  apply le_antisymm _ (hKnonneg q hq w)
  apply le_of_not_gt
  intro hqw
  have h := constrained_remote_read_obstruction k K hK hKnonneg reference selected hr hp hmin
    q hq w hqw edge source hcut readPort hread x y hxy decode a b hab
  exact h.2 (correct w h.1)

end
end OPH.SourceConstrainedRead
