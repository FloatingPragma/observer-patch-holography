import ObserverPatchHolography.EinsteinBranch.ReceiverGravity

/-! Operational identities for an open material crossing and for the ordinary
population transported by an overlap source. All rates use native units. -/
namespace OPH.EinsteinBranch.ContactCensus

noncomputable def first (d x y : ℝ) := (17*d+34*x+65*y)/4
noncomputable def second (d x y : ℝ) := (65*d+65*x+130*y)/4
noncomputable def beat (d x y : ℝ) := -24*d-(31*x+65*y)/2

/-- Both release accelerations and the clock can coincide for a passive
gradient and an ordinary transfer preserving the total of each load type. -/
theorem signed_transfer (d : ℝ) :
    first 0 (-31*d/3) (17*d/3) = first d 0 0 ∧
    second 0 (-31*d/3) (17*d/3) = second d 0 0 ∧
    beat 0 (-31*d/3) (17*d/3) = beat d 0 0 := by
  unfold first second beat
  constructor <;> [ring; constructor <;> ring]

/-- A third, finite swap read of vacancy identifies the written-degree
gradient even for unknown signed changes in both ordinary neighborhoods. -/
theorem recover_written (d x y : ℝ) :
    -4/17 * first d x y + 124/1105 * second d x y + 3/17*(d+x+y) = d := by
  unfold first second
  ring

theorem recover_ordinary (d x y : ℝ) :
    4/17 * first d x y -192/1105 * second d x y +31/17*(d+x+y) = x ∧
    4/65 * second d x y - (d+x+y) = y := by
  unfold first second
  constructor <;> ring

theorem clock_and_release (d x y : ℝ) :
    beat d x y = -2*(second d x y - first d x y) := by
  unfold first second beat
  ring

/-- At a uniform, zero-current launch, H times the launch measures the onsite
potential. The second Born probability derivative measures its local graph
Laplacian. The kinetic coefficient on each active bond is one half. -/
theorem uniform_population {n : ℕ} (site : ℝ) (neighbor : Fin n → ℝ)
    (amplitude : ℝ) :
    2*((site*amplitude)^2 - amplitude *
      (site^2*amplitude + amplitude*(∑ i, (site-neighbor i))/2)) =
      amplitude^2*(∑ i, (neighbor i-site)) := by
  have h : (∑ i, (neighbor i-site)) = -(∑ i, (site-neighbor i)) := by
    rw [← Finset.sum_neg_distrib]
    apply Finset.sum_congr rfl
    intro i _
    ring
  rw [h]
  ring

/-- If the material census changes by b and every written contact has been
restored, the three independent readers all give the same ordinary dose. -/
theorem restored_coupling (b : ℝ) :
    first 0 b 0 / (17/2) = b ∧ second 0 b 0 / (65/4) = b ∧
    beat 0 b 0 / (-31/2) = b := by
  unfold first second beat
  constructor <;> [ring; constructor <;> ring]

/-- Free evolution preserves an energy-block averaged density at every
delay. The two hypotheses are commutation and a unitary right inverse. -/
theorem stationary_under_delay {R : Type*} [Monoid R] (u ρ inverse : R)
    (comm : u*ρ = ρ*u) (unitary : u*inverse=1) :
    u*ρ*inverse = ρ := by rw [comm, mul_assoc, unitary, mul_one]

/-- A public reversal removes a fixed instrumental offset. -/
theorem reversal (forward reverse : ℝ) :
    forward-(forward+reverse)/2 = (forward-reverse)/2 := by ring

/-- The off-source screened equation of a transitive carrier makes the
early ordinary density negative for a positive odd-visit response. -/
theorem exterior_depletion (α p lap : ℝ) (ha : 0 < α) (hα : α < 1)
    (hp : 0 < p) (bal : (1-α)*p+α*lap/12=0) : lap < 0 := by
  by_contra! h
  have hm : 0 < (1-α)*p := mul_pos (by linarith) hp
  have hn : 0 ≤ α*lap := mul_nonneg (le_of_lt ha) h
  linarith

/-- Terminal positions exhaust both ordinary preparations. Any excess mass
in a subset is compensated by the remainder of the finite board. -/
theorem total_redistribution {n : ℕ} (shared fresh : Fin n → ℝ)
    (hn : ∑ i, shared i = ∑ i, fresh i) : (∑ i, (shared i - fresh i)) = 0 := by
  rw [Finset.sum_sub_distrib, hn, sub_self]

/-- A cell projection commuting with restored propagation gives the same
observable after any free delay, for every input in that sector. -/
theorem conserved_cell_read {R : Type*} [Monoid R] (u inverse proj : R)
    (comm : proj*u = u*proj) (left : inverse*u=1) :
    inverse*(proj*u) = proj := by rw [comm, ← mul_assoc, left, one_mul]

/-- A common area coefficient for two fixed-volume interventions requires
the determinant of their measured area and entropy tangents to vanish. -/
theorem common_coefficient (η a s ao so : ℝ) (hr : η*a+s=0) (ho : η*ao+so=0) :
    a*so-s*ao=0 := by
  linear_combination a*ho - ao*hr

/-- The two preparations differ only in admission by the serial record.
A restored one-load kinetic bandwidth of five and rest energy forty bound
the stock-normalized homogeneous response for every admitted layout. -/
theorem restored_energy_ceiling (admitted quiet : ℝ) (hi : admitted ≤ 5)
    (lo : 0 ≤ quiet) : ((101/144 - 25/36)*(admitted-quiet))/40 ≤ 1/1152 := by
  linarith

/-- Any additional regional entropy on the same two measured geometries
must supply this contrast, even if the ordinary intervention changes it. -/
theorem missing_entropy (η a s ao so x xo : ℝ) (hm : ao ≠ 0)
    (hr : η*a + s + x = 0) (ho : η*ao + so + xo = 0) :
    x - (a/ao)*xo = -s + (a/ao)*so := by
  field_simp
  linear_combination ao*hr - a*ho

end OPH.EinsteinBranch.ContactCensus
