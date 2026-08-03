import CapFirstLaw
import ObserverPatchHolography.EinsteinBranch.Entropy

namespace OPH.Thermodynamics

/-!
# Discharging the Einstein-branch first-law premises from thermodynamics

The Einstein branch consumes `FiniteFirstLawPremises`: the entropy
differential equals the modular pairing, the modular generator splits
as `2π B + Z`, and the central pairing equals the edge differential.
This module instantiates that premise package on the thermodynamic
model itself. The variation space is the tangent space of the
probability simplex, the kernel of the total-mass functional. On it
the differential of Shannon entropy at the faithful reference `τ` is
the pairing with `-(1 + log τ)`, the modular pairing is the pairing
with `K = -log τ`, and the two agree exactly because mass-preserving
variations kill the constant `1`. The central split distributes
pointwise from `K = 2π B + Z`, and the model realizes the central
charge as the edge differential. All three premise equalities are
theorems here, so `finite_bulk_edge_central_first_law` applies with
no assumed input beyond the split of the reference. One repair step
supplies a canonical variation: its displacement is mass preserving
by stochasticity, and for a central charge measurable through the
repaired visible datum the central pairing vanishes on it, so the
first-order entropy change of repair is purely bulk-modular.

The exact, non-differential counterpart is `cap_firstLaw_exact`,
restated here as `shannon_diff_eq_pairing_sub_kl`: the entropy
difference is the modular pairing of the displacement minus the
relative entropy, so the linear functionals of this module are the
first-order coefficients of an exact identity rather than a formal
linearization.
-/

variable (Ω : Type*) [Fintype Ω]

/-- The total-mass functional on finite signed variations. -/
noncomputable def massFunctional : (Ω → ℝ) →ₗ[ℝ] ℝ where
  toFun δ := ∑ x, δ x
  map_add' a b := Finset.sum_add_distrib
  map_smul' c a := by
    simp [Finset.mul_sum]

/-- The tangent space of the probability simplex: variations with zero
total mass. -/
noncomputable abbrev MassZero : Submodule ℝ (Ω → ℝ) :=
  LinearMap.ker (massFunctional Ω)

variable {Ω}

@[simp] theorem massFunctional_apply (δ : Ω → ℝ) :
    massFunctional Ω δ = ∑ x, δ x := rfl

variable (Ω) in
/-- Pairing a variation with a fixed observable. -/
noncomputable def pairing (w : Ω → ℝ) : (Ω → ℝ) →ₗ[ℝ] ℝ where
  toFun δ := ∑ x, δ x * w x
  map_add' a b := by
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun x _ => by
      simp only [Pi.add_apply]
      ring
  map_smul' c a := by
    simp only [RingHom.id_apply, smul_eq_mul, Finset.mul_sum]
    exact Finset.sum_congr rfl fun x _ => by
      simp [Pi.smul_apply, smul_eq_mul]
      ring

@[simp] theorem pairing_apply (w δ : Ω → ℝ) :
    pairing Ω w δ = ∑ x, δ x * w x := rfl

/-- The exact identity behind the linearization: the entropy
difference to the reference is the modular pairing of the displacement
minus the relative entropy. -/
theorem shannon_diff_eq_pairing_sub_kl [DecidableEq Ω]
    (p tau : Ω → ℝ) (hτ : ∀ x, 0 < tau x) :
    shannon p - shannon tau
      = pairing Ω (fun x => -Real.log (tau x)) (fun x => p x - tau x)
        - kl p tau := by
  have h := cap_firstLaw_exact p tau hτ
  have hpair : pairing Ω (fun x => -Real.log (tau x))
      (fun x => p x - tau x)
      = (∑ x, p x * (-Real.log (tau x)))
        - ∑ x, tau x * (-Real.log (tau x)) := by
    rw [pairing_apply, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun x _ => by ring
  rw [hpair]
  linarith

variable (Ω) in
/-- The thermodynamic instantiation of the Einstein-branch first-law
data on the simplex tangent space: entropy differential, modular
pairing, bulk and central pairings, and the edge differential realized
by the central charge. -/
noncomputable def thermoFirstLawData (tau Bc Z : Ω → ℝ) :
    OPH.EinsteinBranch.FiniteFirstLawData (MassZero Ω) where
  entropyDifferential :=
    (pairing Ω (fun x => -(1 + Real.log (tau x)))).comp
      (MassZero Ω).subtype
  modularPairing :=
    (pairing Ω (fun x => -Real.log (tau x))).comp
      (MassZero Ω).subtype
  bulkPairing := (pairing Ω Bc).comp (MassZero Ω).subtype
  centralPairing := (pairing Ω Z).comp (MassZero Ω).subtype
  edgeDifferential := (pairing Ω Z).comp (MassZero Ω).subtype

/-- **The Einstein-branch first-law premises hold in the thermodynamic
model.** The entropy differential equals the modular pairing on
mass-preserving variations, the modular pairing splits through
`K = 2π B + Z`, and the central pairing is the edge differential by
construction. -/
theorem thermoFirstLawData_passes (tau Bc Z : Ω → ℝ)
    (hsplit : ∀ x, -Real.log (tau x) = 2 * Real.pi * Bc x + Z x) :
    OPH.EinsteinBranch.FiniteFirstLawPremises
      (thermoFirstLawData Ω tau Bc Z) := by
  refine ⟨?_, ?_, rfl⟩
  · apply LinearMap.ext
    rintro ⟨δ, hδ⟩
    rw [LinearMap.mem_ker, massFunctional_apply] at hδ
    show pairing Ω (fun x => -(1 + Real.log (tau x))) δ
      = pairing Ω (fun x => -Real.log (tau x)) δ
    rw [pairing_apply, pairing_apply]
    have hterm : ∀ x ∈ Finset.univ,
        δ x * -(1 + Real.log (tau x))
          = -δ x + δ x * -Real.log (tau x) := by
      intro x _
      ring
    rw [Finset.sum_congr rfl hterm, Finset.sum_add_distrib,
      Finset.sum_neg_distrib, hδ, neg_zero, zero_add]
  · apply LinearMap.ext
    rintro ⟨δ, hδ⟩
    show pairing Ω (fun x => -Real.log (tau x)) δ
      = ((2 * Real.pi) • (pairing Ω Bc).comp (MassZero Ω).subtype
          + (pairing Ω Z).comp (MassZero Ω).subtype) ⟨δ, hδ⟩
    rw [LinearMap.add_apply, LinearMap.smul_apply]
    show pairing Ω (fun x => -Real.log (tau x)) δ
      = (2 * Real.pi) • pairing Ω Bc δ + pairing Ω Z δ
    rw [smul_eq_mul, pairing_apply, pairing_apply, pairing_apply,
      Finset.mul_sum, ← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro x _
    rw [hsplit x]
    ring

/-- The composed Einstein-branch first law, instantiated: on every
mass-preserving variation the entropy differential is
`2π` times the bulk pairing plus the edge differential. -/
theorem thermo_first_law_on_simplex_tangent (tau Bc Z : Ω → ℝ)
    (hsplit : ∀ x, -Real.log (tau x) = 2 * Real.pi * Bc x + Z x)
    (δρ : MassZero Ω) :
    (thermoFirstLawData Ω tau Bc Z).entropyDifferential δρ
      = 2 * Real.pi * (thermoFirstLawData Ω tau Bc Z).bulkPairing δρ
        + (thermoFirstLawData Ω tau Bc Z).edgeDifferential δρ :=
  OPH.EinsteinBranch.finite_bulk_edge_central_first_law
    (thermoFirstLawData Ω tau Bc Z)
    (thermoFirstLawData_passes tau Bc Z hsplit) δρ

section RepairVariation

variable {B : Type*} [DecidableEq B] [DecidableEq Ω]

/-- One repair step supplies a canonical mass-preserving variation:
the displacement of any state under the repair kernel lies in the
simplex tangent space. -/
theorem repair_variation_mem_massZero (π : Ω → ℝ) (b : Ω → B)
    (hπ : ∀ x, 0 < π x) (p : Ω → ℝ) :
    (fun y => push p (heatBath π b) y - p y) ∈ MassZero Ω := by
  rw [LinearMap.mem_ker, massFunctional_apply,
    Finset.sum_sub_distrib,
    push_total p (heatBath π b) (heatBath_row_sum hπ), sub_self]

/-- For a central charge measurable through the repaired visible
datum, the central pairing vanishes on the repair displacement: the
first-order entropy change of one repair step is purely bulk-modular.
-/
theorem repair_variation_central_pairing_zero (π : Ω → ℝ) (b : Ω → B)
    (hπ : ∀ x, 0 < π x) (Z : Ω → ℝ)
    (hZ : ∀ x y, b x = b y → Z x = Z y) (p : Ω → ℝ) :
    pairing Ω Z (fun y => push p (heatBath π b) y - p y) = 0 := by
  rw [pairing_apply]
  have hsplit : ∀ y ∈ Finset.univ,
      (push p (heatBath π b) y - p y) * Z y
        = push p (heatBath π b) y * Z y - p y * Z y := by
    intro y _
    ring
  rw [Finset.sum_congr rfl hsplit, Finset.sum_sub_distrib,
    push_heatBath_fixes_mean hπ Z hZ p, sub_self]

end RepairVariation

end OPH.Thermodynamics

#print axioms OPH.Thermodynamics.shannon_diff_eq_pairing_sub_kl
#print axioms OPH.Thermodynamics.thermoFirstLawData_passes
#print axioms OPH.Thermodynamics.thermo_first_law_on_simplex_tangent
#print axioms OPH.Thermodynamics.repair_variation_mem_massZero
#print axioms OPH.Thermodynamics.repair_variation_central_pairing_zero
