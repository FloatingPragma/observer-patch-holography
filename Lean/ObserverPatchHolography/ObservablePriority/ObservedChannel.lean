import ObserverPatchHolography.ObservablePriority.NativeColumns

/-! Trace, coherence and pushforward to a patch's addressed observation. -/
namespace StrictQuantum.Observed
open StrictQuantum
open scoped BigOperators
noncomputable section

variable {E V : Type} [DecidableEq V]

abbrev density (V : Type) := Loads V → Loads V → ℂ

theorem projection_idempotent (q : Screen E V) (es : List E) (M : E → Bool) (ψ : StateVector V) :
    project q es M (project q es M ψ) = project q es M ψ := by
  funext x
  by_cases h : flags (scan q es x).2 = es.map M <;> simp [project,h]

theorem kraus_complex_linear (q : Screen E V) (es : List E) (M : E → Bool) (z : ℂ)
    (ψ φ : StateVector V) (x : Loads V) (c : Coins (scan q es x).2) :
    kraus q es M (fun t => z * ψ t + φ t) x c = z * kraus q es M ψ x c + kraus q es M φ x c := by
  simp only [kraus]
  ring

/-- A positive Gram input on the physical screen and *any* spectator system
has a positive Gram factor after this branch. The spectator indices are
not measured or modified. Zero rows off the physical image extend the factor. -/
theorem spectator_gram {I T : Type} [Fintype T] (q : Screen E V) (es : List E)
    (M : E → Bool) (F : Loads V → I → T → ℂ) (x y : Loads V) (a b : I) :
    (Native.matrixCoefficient q M es x : ℂ) *
      (∑ t, F x a t * star (F y b t)) * star (Native.matrixCoefficient q M es y : ℂ) =
    ∑ t, ((Native.matrixCoefficient q M es x : ℂ) * F x a t) *
      star ((Native.matrixCoefficient q M es y : ℂ) * F y b t) := by
  rw [Finset.mul_sum, Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro t _
  rw [star_mul]
  ring

/-- Matrix unit rule for K_M ρ K_M†, at two *physical* output coordinates
(see physical_injective). The coin coefficients are positive real numbers. -/
def entry (q : Screen E V) (es : List E) (M : E → Bool) (ρ : density V)
    (x y : Loads V) (_ : Coins (scan q es x).2) (_ : Coins (scan q es y).2) : ℂ :=
  (Native.matrixCoefficient q M es x : ℂ) * ρ x y *
    star (Native.matrixCoefficient q M es y : ℂ)

theorem same_outcome_coherence (q : Screen E V) (es : List E) (M : E → Bool)
    (ρ : density V) (x y : Loads V) (c : Coins (scan q es x).2) (d : Coins (scan q es y).2)
    (hx : flags (scan q es x).2 = es.map M) (hy : flags (scan q es y).2 = es.map M) :
    entry q es M ρ x y c d =
      (weight (scan q es x).2 : ℂ) * ρ x y * (weight (scan q es y).2 : ℂ) := by
  simp [entry, Native.matrix_outcome_projection, hx, hy]

theorem offdiagonal_survives (q : Screen E V) (es : List E) (M : E → Bool)
    (ρ : density V) (x y : Loads V) (c : Coins (scan q es x).2) (d : Coins (scan q es y).2)
    (hx : flags (scan q es x).2 = es.map M) (hy : flags (scan q es y).2 = es.map M)
    (hne : ρ x y ≠ 0) : entry q es M ρ x y c d ≠ 0 := by
  rw [same_outcome_coherence q es M ρ x y c d hx hy]
  have ax : (weight (scan q es x).2 : ℂ) ≠ 0 := by exact_mod_cast (weight_pos _).ne'
  have ay : (weight (scan q es y).2 : ℂ) ≠ 0 := by exact_mod_cast (weight_pos _).ne'
  exact mul_ne_zero (mul_ne_zero ax hne) ay

theorem wrong_outcome_zero (q : Screen E V) (es : List E) (M : E → Bool)
    (ρ : density V) (x y : Loads V) (c : Coins (scan q es x).2) (d : Coins (scan q es y).2)
    (h : flags (scan q es x).2 ≠ es.map M ∨ flags (scan q es y).2 ≠ es.map M) :
    entry q es M ρ x y c d = 0 := by
  rcases h with h | h <;> simp [entry, Native.matrix_outcome_projection, h]

theorem trace_one_column (q : Screen E V) (es : List E) (M : E → Bool)
    (ρ : density V) (x : Loads V) :
    (∑ c : Coins (scan q es x).2, entry q es M ρ x x c c) =
      if flags (scan q es x).2 = es.map M then ρ x x else 0 := by
  by_cases h : flags (scan q es x).2 = es.map M
  · simp only [h, if_true, same_outcome_coherence q es M ρ x x _ _ h h]
    calc
      _ = ∑ _ : Coins (scan q es x).2,
          ρ x x * ((weight (scan q es x).2 * weight (scan q es x).2 : ℝ) : ℂ) := by
        apply Finset.sum_congr rfl
        intro c _
        push_cast
        ring
      _ = _ := by
        rw [← Finset.mul_sum, ← Complex.ofReal_sum, normalized]
        simp
  · simp [entry, Native.matrix_outcome_projection, h]

/-- Only words that occur on the blank input space are needed. This also makes
the normalization independent of unused Boolean coordinates outside the screen. -/
def outcomes [Fintype V] (q : Screen E V) (es : List E) : Finset (List Bool) :=
  Finset.univ.image (fun x : Loads V => flags (scan q es x).2)

def mass [Fintype V] (q : Screen E V) (es : List E) (w : List Bool)
    (d : Loads V → ℂ) : ℂ :=
  ∑ x, if flags (scan q es x).2 = w then d x else 0

theorem branch_trace [Fintype V] (q : Screen E V) (es : List E) (M : E → Bool)
    (ρ : density V) :
    (∑ x, ∑ c : Coins (scan q es x).2, entry q es M ρ x x c c) =
      mass q es (es.map M) (fun x => ρ x x) := by
  simp [mass, trace_one_column]

theorem marginal [Fintype V] {B : Type} [DecidableEq B] (q : Screen E V) (es : List E)
    (observe : List Bool → B) (b : B) (d : Loads V → ℂ) :
    (∑ w ∈ (outcomes q es).filter (fun w => observe w = b), mass q es w d) =
      ∑ x, if observe (flags (scan q es x).2) = b then d x else 0 := by
  unfold mass
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro x _
  have hw : flags (scan q es x).2 ∈ outcomes q es := Finset.mem_image.mpr ⟨x, Finset.mem_univ _, rfl⟩
  simp [hw]

theorem total_trace [Fintype V] (q : Screen E V) (es : List E) (ρ : density V) :
    (∑ w ∈ outcomes q es, mass q es w (fun x => ρ x x)) = ∑ x, ρ x x := by
  simpa using marginal q es (fun _ => ()) () (fun x => ρ x x)

#print axioms same_outcome_coherence
#print axioms projection_idempotent
#print axioms kraus_complex_linear
#print axioms spectator_gram
#print axioms offdiagonal_survives
#print axioms wrong_outcome_zero
#print axioms branch_trace
#print axioms marginal
#print axioms total_trace

end
end StrictQuantum.Observed
