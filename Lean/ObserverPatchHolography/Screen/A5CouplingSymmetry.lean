import A5PortAction

namespace OPH.A5CouplingSymmetry

open OPH.A5PortAction

/-! # Universal cap coupling from icosahedral symmetry

The Einstein-branch coupling clause asks that per-cap coupling ratios be one
constant across the axis-cap family.  Measured on the present heuristic
source dynamics this fails (issue #576: 68 percent spread).  This file proves
that the clause is a theorem, not a measurement, for any source law with the
icosahedral symmetry of the carrier: for every port-data assignment, the
rotation-group average of the cap sum is the same for every port, and for
any invariant assignment the cap sums themselves are equal.  Consequently
the per-cap ratios of any two group-averaged (or invariant) readouts, such
as entropy over stress flux, are equal across the family with zero spread.

Everything is stated over `ℚ` at the `Nat` port level of `A5PortAction`,
with the finite facts (cap images, composition law, right-translation
bijection, transitivity) checked by kernel `decide` and the quantified
statements over arbitrary data proved symbolically from them.

Boundary: this reduces the universality clause to one hypothesis, that the
source law is `A5`-equivariant.  Whether a given implemented law satisfies
that hypothesis is a finite check on the law, owned by the simulator; no
physical coupling value, vacuum, or promotion is implied. -/

/-- The cap of a port: the port together with its neighbor ring. -/
def capPorts (a : Nat) : List Nat := a :: natNeighbors a

/-- The cap sum of a port-data assignment. -/
def capSum (d : Nat → ℚ) (a : Nat) : ℚ := ((capPorts a).map d).sum

/-- The rotation-group average of the cap sum. -/
def groupAverage (d : Nat → ℚ) (a : Nat) : ℚ :=
  (perms.map (fun p => capSum (fun k => d (app p k)) a)).sum

/-! ## Finite facts, by kernel decide -/

/-- Every cap member is a valid port. -/
theorem capPorts_lt : ∀ a, a < 12 → ∀ x ∈ capPorts a, x < 12 := by decide

/-- Rotations map caps to caps: the image of the cap of `a` is a permutation
of the cap of the image of `a`. -/
theorem cap_image_perm :
    ∀ p ∈ perms, ∀ a, a < 12 →
      ((capPorts a).map (app p)).Perm (capPorts (app p a)) := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 4096 in
/-- The listed composition is the pointwise composition on ports. -/
theorem app_comp :
    ∀ p ∈ perms, ∀ q ∈ perms, ∀ k, k < 12 →
      app (comp p q) k = app p (app q k) := by decide

set_option maxHeartbeats 8000000 in
set_option maxRecDepth 8192 in
/-- Right translation by a listed rotation permutes the listed group. -/
theorem right_translate_perm :
    ∀ p ∈ perms, (perms.map (fun q => comp q p)).Perm perms := by decide

/-- The action reaches every port from every port. -/
theorem transitive_pairs :
    ∀ a, a < 12 → ∀ b, b < 12 → ∃ p, p ∈ perms ∧ app p a = b := by decide

/-! ## Symbolic consequences for arbitrary data -/

/-- Equivariance of the cap sum: evaluating at a rotated port equals
evaluating the rotated data at the original port. -/
theorem capSum_equivariant
    (d : Nat → ℚ) {p : List Nat} (hp : p ∈ perms) {a : Nat} (ha : a < 12) :
    capSum d (app p a) = capSum (fun k => d (app p k)) a := by
  unfold capSum
  have himage := cap_image_perm p hp a ha
  have hmap : ((capPorts a).map (app p)).map d
      = (capPorts a).map (fun k => d (app p k)) := by
    simp [List.map_map, Function.comp]
  calc ((capPorts (app p a)).map d).sum
      = (((capPorts a).map (app p)).map d).sum :=
        ((himage.map d).sum_eq).symm
    _ = ((capPorts a).map (fun k => d (app p k))).sum := by rw [hmap]

/-- Invariant data have equal cap sums on every port: universality with zero
spread for `A5`-invariant readouts. -/
theorem capSum_const_of_invariant
    (d : Nat → ℚ)
    (hinv : ∀ p ∈ perms, ∀ k, k < 12 → d (app p k) = d k)
    {a b : Nat} (ha : a < 12) (hb : b < 12) :
    capSum d a = capSum d b := by
  obtain ⟨p, hp, hab⟩ := transitive_pairs a ha b hb
  have hrot : capSum d (app p a) = capSum (fun k => d (app p k)) a :=
    capSum_equivariant d hp ha
  have hpoint : (capPorts a).map (fun k => d (app p k))
      = (capPorts a).map d := by
    refine List.map_congr_left ?_
    intro x hx
    exact hinv p hp x (capPorts_lt a ha x hx)
  have : capSum (fun k => d (app p k)) a = capSum d a := by
    unfold capSum
    rw [hpoint]
  rw [hab] at hrot
  rw [this] at hrot
  exact hrot.symm


/-- The group-averaged cap sum is the same for every port, for every data
assignment: universality of averaged readouts is unconditional. -/
theorem groupAverage_port_independent
    (d : Nat → ℚ) {a b : Nat} (ha : a < 12) (hb : b < 12) :
    groupAverage d a = groupAverage d b := by
  obtain ⟨p₀, hp₀, hab⟩ := transitive_pairs a ha b hb
  unfold groupAverage
  have hstep : ∀ q ∈ perms,
      capSum (fun k => d (app q k)) b
        = capSum (fun k => d (app (comp q p₀) k)) a := by
    intro q hq
    have hrot :
        capSum (fun k => d (app q k)) (app p₀ a)
          = capSum (fun k => d (app q (app p₀ k))) a :=
      capSum_equivariant (fun k => d (app q k)) hp₀ ha
    have hpoint : (capPorts a).map (fun k => d (app q (app p₀ k)))
        = (capPorts a).map (fun k => d (app (comp q p₀) k)) := by
      refine List.map_congr_left ?_
      intro x hx
      rw [app_comp q hq p₀ hp₀ x (capPorts_lt a ha x hx)]
    have hswap :
        capSum (fun k => d (app q (app p₀ k))) a
          = capSum (fun k => d (app (comp q p₀) k)) a := by
      unfold capSum
      rw [hpoint]
    rw [hab] at hrot
    exact hrot.trans hswap
  have hcongr :
      (perms.map (fun q => capSum (fun k => d (app q k)) b))
        = perms.map (fun q => capSum (fun k => d (app (comp q p₀) k)) a) := by
    refine List.map_congr_left ?_
    intro q hq
    exact hstep q hq
  have htrans :
      perms.map (fun q => capSum (fun k => d (app (comp q p₀) k)) a)
        = (perms.map (fun q => comp q p₀)).map
            (fun r => capSum (fun k => d (app r k)) a) := by
    simp [List.map_map, Function.comp]
  have hperm := (right_translate_perm p₀ hp₀).map
    (fun r => capSum (fun k => d (app r k)) a)
  calc (perms.map (fun q => capSum (fun k => d (app q k)) a)).sum
      = ((perms.map (fun q => comp q p₀)).map
          (fun r => capSum (fun k => d (app r k)) a)).sum :=
        (hperm.sum_eq).symm
    _ = (perms.map
          (fun q => capSum (fun k => d (app (comp q p₀) k)) a)).sum := by
        rw [htrans]
    _ = (perms.map (fun q => capSum (fun k => d (app q k)) b)).sum := by
        rw [hcongr]

/-- Universal coupling: for any two group-averaged readouts, for example an
entropy readout over a stress-flux readout, the per-cap ratio is the same
at every port.  The Einstein-branch universality clause holds by symmetry
for every `A5`-equivariant source law. -/
theorem coupling_ratio_universal
    (entropy flux : Nat → ℚ) {a b : Nat} (ha : a < 12) (hb : b < 12) :
    groupAverage entropy a * groupAverage flux b
      = groupAverage entropy b * groupAverage flux a := by
  rw [groupAverage_port_independent entropy ha hb,
    groupAverage_port_independent flux ha hb]

end OPH.A5CouplingSymmetry
