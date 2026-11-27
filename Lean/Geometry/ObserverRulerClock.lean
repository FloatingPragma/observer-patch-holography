import Mathlib.Tactic

/-!
# Pulse periods and a ruler on a receiver's clock

`flight t d` means that a labelled pulse emitted at tick `t` can arrive
after `d` ticks. If material contacts only close and a transfer takes one
tick, replaying any successful itinerary earlier gives the same duration.
All period statements below concern received crests.
-/

namespace OPH.ObserverRulerClock

def Earliest (flight : ℕ → ℕ → Prop) (t d : ℕ) : Prop :=
  flight t d ∧ ∀ j, flight t j → d ≤ j

/-- The received period is the emitted period plus the increase in minimum flight time. -/
theorem received_period {flight : ℕ → ℕ → Prop}
    (retime : ∀ s t d, s ≤ t → flight t d → flight s d)
    {s p a b : ℕ} (first : Earliest flight s a) (second : Earliest flight (s+p) b) :
    a ≤ b ∧ (s+p+b) - (s+a) = p + (b-a) := by
  have h := first.2 b (retime s (s+p) b (by omega) second.1)
  omega

/-- Translation-invariant flights give the emitted period at reception, at every cadence. -/
theorem stationary_period {flight : ℕ → ℕ → Prop} {s p a b : ℕ}
    (stationary : ∀ d, flight s d ↔ flight (s+p) d)
    (first : Earliest flight s a) (second : Earliest flight (s+p) b) :
    (s+p+b) - (s+a) = p := by
  have h₁ := first.2 b ((stationary b).mpr second.1)
  have h₂ := second.2 a ((stationary a).mp first.1)
  omega

/-- Two unit-neighbour sources received together issued their crests within one tick.
The last hop of each bilateral length certificate allows entry into the other
source one tick before its emission. The first-arrival clock is strictly ordered
on successful emissions. A missing arrival is represented by `none`. -/
theorem unit_ruler_emission_gap {A B : ℕ → Option ℕ} {a b t : ℕ}
    (ha : A a = some t) (hb : B b = some t)
    (fifoA : ∀ i j x y, i < j → A i = some x → A j = some y → x < y)
    (fifoB : ∀ i j x y, i < j → B i = some x → B j = some y → x < y)
    (viaB : 0 < b → ∃ r, A (b-1) = some r ∧ r ≤ t)
    (viaA : 0 < a → ∃ r, B (a-1) = some r ∧ r ≤ t) :
    a ≤ b+1 ∧ b ≤ a+1 := by
  constructor
  · by_contra! h
    obtain ⟨r, hr, hle⟩ := viaA (by omega)
    have := fifoB b (a-1) t r (by omega) hb hr
    omega
  · by_contra! h
    obtain ⟨r, hr, hle⟩ := viaB (by omega)
    have := fifoA a (b-1) t r (by omega) ha hr
    omega

/-- A complete array with an unknown emission offset distinguishes exactly the
arrival words that do not differ by a common translation. -/
theorem same_relative_arrivals_iff {ι : Type*} (reference : ι) (A B : ι → ℤ) :
    (∀ j, A j - A reference = B j - B reference) ↔
    ∃ k, ∀ j, B j = A j + k := by
  constructor
  · intro h
    refine ⟨B reference - A reference, ?_⟩
    intro j
    have := h j
    omega
  · rintro ⟨k, h⟩ j
    rw [h j, h reference]
    omega

end OPH.ObserverRulerClock
