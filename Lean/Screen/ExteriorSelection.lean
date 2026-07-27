import Mathlib

namespace OPH.ExteriorSelection

/-! # Exterior matter selection and statistics forcing

Machine checks for the #314/#567 source-bound closures.

1. The exhaustive selection scan: over the derived primitive block charges
   `(a, b) = (-2, 3)` in `q = 6Y` units, the ten nontrivial isotypic
   components of the exterior module admit exactly two nonempty chiral
   anomaly-free subsets, namely the two fermionic-parity sectors, and they
   are exchanged by charge conjugation.  The parity grading of the matter
   object is therefore an output of the scan.
2. The section obstruction: a group whose nonidentity square roots of the
   identity are unique admits no injective homomorphism from the Klein four
   group, so the measured transport double cover (unique involution `-1`)
   is non-split over every Klein four-subgroup of the deck.
3. The polarization commutant: on the line-class group `ZMod 6 x ZMod 6`
   with the Dirac pairing, any class mutually local with every electric
   class is electric, so the electric lattice is the unique maximal
   mutually-local lattice containing the realized electric lines.
4. The sector-lift menu: the members of the order-`d` subgroups of `ZMod 6`
   are exactly the divisor menu used by the four-form discrimination.

BOUNDARY.  These are finite statements about the declared component table,
abstract groups with unique involutions, and `ZMod 6` arithmetic.  They do
not construct continuum fields, four-dimensional instanton sectors, theta
periodicity, laboratory observables, or quantum field theory. -/

/-- Complex dimension of the color factor per component. -/
def colorDim : Fin 10 → ℤ := ![3, 1, 3, 3, 1, 3, 3, 1, 3, 1]

/-- Complex dimension of the weak factor per component. -/
def weakDim : Fin 10 → ℤ := ![1, 2, 1, 2, 1, 1, 2, 1, 1, 2]

/-- Integer charge `q = 6Y` per component over the derived pair `(-2, 3)`. -/
def charge : Fin 10 → ℤ := ![-2, 3, -4, 1, 6, 4, -1, -6, 2, -3]

/-- Whether the component carries the color triplet or antitriplet. -/
def isColored : Fin 10 → Bool :=
  ![true, false, true, true, false, true, true, false, true, false]

/-- Whether the component carries the weak doublet. -/
def isDoublet : Fin 10 → Bool :=
  ![false, true, false, true, false, false, true, false, false, true]

/-- Fermionic parity: `true` on the odd exterior degrees. -/
def odd : Fin 10 → Bool :=
  ![true, true, false, false, false, true, true, true, false, false]

/-- Charge conjugation on the component index set. -/
def conj : Fin 10 → Fin 10 := ![8, 9, 5, 6, 7, 2, 3, 4, 0, 1]

/-- Bitmask membership for a subset of the ten components. -/
def mem (m : Nat) (i : Fin 10) : Bool := Nat.testBit m i.val

/-- Masked sum of an integer weight over the selected components. -/
def maskSum (f : Fin 10 → ℤ) (m : Nat) : ℤ :=
  (List.finRange 10).foldr (fun i acc => if mem m i then f i + acc else acc) 0

/-- Gravitational anomaly of a selection. -/
def grav (m : Nat) : ℤ := maskSum (fun i => colorDim i * weakDim i * charge i) m

/-- `SU(3)^2 U(1)` anomaly of a selection. -/
def su3 (m : Nat) : ℤ :=
  maskSum (fun i => if isColored i then weakDim i * charge i else 0) m

/-- `SU(2)^2 U(1)` anomaly of a selection. -/
def su2 (m : Nat) : ℤ :=
  maskSum (fun i => if isDoublet i then colorDim i * charge i else 0) m

/-- Cubic `U(1)^3` anomaly of a selection. -/
def cube (m : Nat) : ℤ :=
  maskSum (fun i => colorDim i * weakDim i * charge i ^ 3) m

/-- Chirality: no component appears together with its conjugate. -/
def chiral (m : Nat) : Bool :=
  (List.finRange 10).all fun i => !(mem m i && mem m (conj i))

/-- All four anomaly constraints vanish. -/
def anomalyFree (m : Nat) : Bool :=
  grav m == 0 && su3 m == 0 && su2 m == 0 && cube m == 0

/-- The even-parity sector `{2, 3, 4, 8, 9}` as a mask. -/
def evenMask : Nat := 796

/-- The odd-parity sector `{0, 1, 5, 6, 7}` as a mask. -/
def oddMask : Nat := 227

theorem evenMask_is_even_sector :
    ∀ i : Fin 10, mem evenMask i = !odd i := by decide

theorem oddMask_is_odd_sector :
    ∀ i : Fin 10, mem oddMask i = odd i := by decide

set_option maxRecDepth 8192 in
/-- **The selection theorem.**  Among all 1024 subsets of the ten
components, the nonempty chiral anomaly-free selections are exactly the two
fermionic-parity sectors. -/
theorem selection_unique :
    ∀ m : Fin 1024,
      m.val ≠ 0 → chiral m.val = true → anomalyFree m.val = true →
        m.val = evenMask ∨ m.val = oddMask := by decide +kernel

/-- Both parity sectors do survive the scan. -/
theorem parity_sectors_survive :
    chiral evenMask = true ∧ anomalyFree evenMask = true ∧
      chiral oddMask = true ∧ anomalyFree oddMask = true := by decide

/-- Charge conjugation exchanges the two survivors. -/
theorem conj_exchanges_survivors :
    ∀ i : Fin 10, mem evenMask i = mem oddMask (conj i) := by decide

/-- Witten parity is automatic on both survivors: each carries four weak
doublet slots. -/
theorem witten_automatic :
    maskSum (fun i => if isDoublet i then colorDim i else 0) evenMask = 4 ∧
      maskSum (fun i => if isDoublet i then colorDim i else 0) oddMask = 4 := by
  decide

/-- **Section obstruction.**  A group in which nonidentity square roots of
the identity are unique admits no injective homomorphism from the Klein
four group.  Applied to the measured transport double cover, whose unique
involution is `-1`, this is the non-splitness of the deck lift over every
Klein four-subgroup. -/
theorem no_klein_section_of_unique_involution {G : Type*} [Group G]
    (huniq : ∀ x y : G, x ≠ 1 → x ^ 2 = 1 → y ≠ 1 → y ^ 2 = 1 → x = y)
    (f : Multiplicative (ZMod 2 × ZMod 2) →* G) :
    ¬ Function.Injective f := by
  intro hinj
  set a : Multiplicative (ZMod 2 × ZMod 2) :=
    Multiplicative.ofAdd ((1, 0) : ZMod 2 × ZMod 2) with ha
  set b : Multiplicative (ZMod 2 × ZMod 2) :=
    Multiplicative.ofAdd ((0, 1) : ZMod 2 × ZMod 2) with hb
  have hane : a ≠ b := by decide
  have hsq : ∀ x : Multiplicative (ZMod 2 × ZMod 2), x ^ 2 = 1 := by decide
  have hfa2 : f a ^ 2 = 1 := by rw [← map_pow, hsq, map_one]
  have hfb2 : f b ^ 2 = 1 := by rw [← map_pow, hsq, map_one]
  have hfa1 : f a ≠ 1 := by
    intro h
    exact (by decide : a ≠ (1 : Multiplicative (ZMod 2 × ZMod 2)))
      (hinj (h.trans (map_one f).symm))
  have hfb1 : f b ≠ 1 := by
    intro h
    exact (by decide : b ≠ (1 : Multiplicative (ZMod 2 × ZMod 2)))
      (hinj (h.trans (map_one f).symm))
  exact hane (hinj (huniq (f a) (f b) hfa1 hfa2 hfb1 hfb2))

/-- Dirac pairing on the line-class group. -/
def dirac (x y : ZMod 6 × ZMod 6) : ZMod 6 := x.1 * y.2 - x.2 * y.1

/-- **Polarization commutant.**  Any line class mutually local with every
electric class is electric, so the electric lattice is the unique maximal
mutually-local lattice containing the realized electric lines. -/
theorem electric_commutant :
    ∀ x : ZMod 6 × ZMod 6, (∀ e : ZMod 6, dirac x (e, 0) = 0) → x.2 = 0 := by
  decide

/-- The electric lattice is itself mutually local. -/
theorem electric_isotropic :
    ∀ e e' : ZMod 6, dirac (e, 0) (e', 0) = 0 := by decide

/-- **Sector-lift menu.**  The order-`d` subgroup members of the measured
class group are exactly the divisor menu the four-form discrimination
uses. -/
theorem lift_menu :
    (∀ c : ZMod 6, (1 : ZMod 6) * c = 0 ↔ c = 0) ∧
      (∀ c : ZMod 6, 2 * c = 0 ↔ c = 0 ∨ c = 3) ∧
      (∀ c : ZMod 6, 3 * c = 0 ↔ c = 0 ∨ c = 2 ∨ c = 4) ∧
      (∀ c : ZMod 6, 6 * c = 0) := by decide

#print axioms selection_unique
#print axioms parity_sectors_survive
#print axioms conj_exchanges_survivors
#print axioms witten_automatic
#print axioms no_klein_section_of_unique_involution
#print axioms electric_commutant
#print axioms lift_menu

end OPH.ExteriorSelection
