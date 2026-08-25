import Geometry.ScreenCarrierMapCandidate
import CertifiedScaledStepInstrument

set_option autoImplicit false

open scoped BigOperators

namespace OPH.CarrierDynamicsCompatibility

open OPH.ScreenCarrierMapCandidate
open OPH.SeamCurrentCarrierQuotient
open OPH.LocalFaceMaxwellAction
open OPH.DiscreteCoulombGreen
open OPH.TemporalMaxwellEvolution
open OPH.ScaledMaxwellStability
open OPH.CertifiedScaledStepInstrument
open OPH.DiscreteRefinement (Barycentric)
open OPH.C1Lorentz (Spatial Herm2 spatialDot)

/-!
# Dynamics compatibility of the screen carrier-map candidate

The carrier-map candidate of `Geometry/ScreenCarrierMapCandidate.lean`
assigns exact `ℤ[φ]` rays to the twelve committed ports and realizes the
two committed generator permutations `genA` and `genB` by the exact
rotation matrices `rotAZ` and `rotBZ`.  This file proves the dynamics
half of completion item 2 of the 2026-08-24 deep audit for issue #740:
the incidence symmetries realized geometrically through the candidate
embedding transport the committed screen evolution data exactly.  Every
statement is compatibility of the CANDIDATE; no theorem identifies the
candidate as the physical carrier.

WHAT IS PROVED.

* Transport structure.  `DynamicsTransport` packages a port bijection, a
  seam bijection with a seam sign, and a face bijection, constrained by
  two exact clauses on the committed tables: the endpoint clause ties the
  seam bijection to the port bijection through the committed seam table
  with the declared sign, and the incidence clause ties the face
  bijection to the seam bijection through the committed signed face
  incidence.  Both clauses are proved by kernel computation for the two
  generator transports (`genATransport`, `genBTransport`, tables searched
  offline and re-verified here row by row), hold for the identity, and
  are closed under composition (`DynamicsTransport.comp`).

* Operator intertwining.  For every inhabitant of the structure, the
  committed port coboundary, face curvature, face codifferential, and
  incidence boundary intertwine exactly with the induced pullbacks
  (`coboundary_intertwine`, `curvature_intertwine`,
  `codifferential_intertwine`, `boundary_intertwine`), and the seam,
  face, and port pairings are invariant (`seamInner_map`,
  `seamEnergy_map`, `faceInner_map`, `portInner_map`).

* Evolution transport, at every step `h` uniformly.  The scaled electric
  field and the magnetic field transport exactly
  (`electric_intertwine`, `magnetic_intertwine`); the scaled Ampere
  update holds for a transported history exactly when it holds for the
  original (`ampere_intertwine`, both directions); the Gauss constraint
  transports with the pulled-back charge (`gauss_intertwine`, both
  directions); the scaled staggered form is invariant
  (`energy_intertwine`); and the transport commutes with the committed
  time-dependent gauge transformation (`gauge_transport_commute`).

* Sixty-element coverage.  `wordTransport` assigns a transport to every
  word in the two generators, with port component equal to the committed
  word action (`wordTransport_port`).  Through the committed word
  certificate `perms_generated_by_two_rotations`, every row of the
  committed sixty-permutation list is the port component of a transport
  (`listed_perm_transport`), and every such row therefore transports the
  scaled evolution and its staggered form exactly
  (`committed_perm_dynamics_compatibility`).

* Geometric realization for the generators through the candidate
  embedding.  At the port level the transports act by the exact rotation
  matrices with the declared factor two (`rotA_ray_vsc`,
  `rotB_ray_vsc`); at the seam level the rotation carries the unordered
  endpoint-ray pair of a committed seam to the pair of its image seam
  (`genA_seam_geometric`, `genB_seam_geometric`); at the face level the
  rotation carries the candidate face-image multiset to the image-face
  multiset (`genA_face_geometric`, `genB_face_geometric`).  The
  geometric realization is proved for the generating set; words act
  geometrically through composition of the two matrices, and this file
  keeps the word level combinatorial.

* Certified-step receipts transport.  At the committed certified step
  `4/5` of `Screen/CertifiedScaledStepInstrument.lean`, every transport
  carries zero-current solutions to zero-current solutions with the same
  staggered form at every step count and the same rational electric
  energy bound (`certified_step_transport`).

* Refinement compatibility at the committed barycentric level.  The
  matrix transport of the carrier extension commutes with the committed
  refinement by exact scaling (`carrier_transport_refine`, algebraic
  `SameRay` corollary `carrier_transport_refine_sameRay`, including its
  explicitly degenerate zero cases); the aligned face-level
  transport identity holds for both generators (`bary_transport_A`,
  `bary_transport_B`), with the coordinate alignment given by exact
  tables re-verified against the committed face table (`alignA_sound`,
  `alignB_sound`); and the two routes, refine then transport and
  transport then refine, agree exactly (`bary_transport_A_routes`,
  `bary_transport_B_routes`).  The committed corpus carries no
  frequency-`n` incidence complex, so the barycentric level is the level
  at which refinement compatibility is provable.

* Interaction candidate, declared.  `interactionCandidate` couples the
  worldline increment of the canonical Lorentz module to the seam
  electric field through the embedded seam direction vectors.  Gauge
  invariance is exact at every declared relative normalization
  (`interaction_gauge_invariant`); the zero-normalization value vanishes
  identically, recovering the decoupled direct-sum regime of the
  committed joint action as the regression limit
  (`interaction_zero_coupling`); and the normalization is a declared
  free constant, proved non-forced: the invariance clause holds at every
  value while the functional separates values on an explicit
  configuration (`interaction_normalization_not_forced`).

WHAT IS NOT SUPPLIED (audit item (c)).  No metric calibration, unit,
clock rate, or physical scale is attached: the intertwining theorems
hold at every step `h` uniformly, so they select no step and no
calibration.  No interaction is derived: the interaction candidate is a
declared coupling whose shape is chosen, whose relative normalization is
a free constant proved non-forced, and for which no equation of motion,
no force law, no back-reaction, and no derivation from the three axioms
is stated; equivariance of the interaction candidate under the
transports is not stated.  No observer readout is constructed: no
theorem maps a field history or a worldline to an observation.  The
physical propagation, source-attachment, and calibration premise rows of
the committed register stay open and are not consumed.  Compatibility
does not select the embedding or a frame: two distinct transports
satisfy every clause (`transport_selection_not_forced`), and the parent
module exhibits two distinct embeddings satisfying its clauses, so the
clauses of this file are a selection surface for candidates and force no
identification of ports with physical directions.

FALSIFIER.  The module fails if some listed permutation row has no
transport, if an intertwining identity misses a sign or a term, if the
staggered form changes under some transport, if a generator transport
disagrees with the geometric action of its rotation matrix on ports,
seams, or faces, if an alignment table row disagrees with the committed
face table, if the two barycentric routes disagree, or if the gauge
variation of the interaction candidate is nonzero.

Axiom audit.  Every proof composes committed receipts with exact
arithmetic and kernel `decide` checks on the committed integer tables;
no project axiom and no native decision procedure is used.  The audit
lines at the end of the file show at most `propext`, `Classical.choice`,
and `Quot.sound`.
-/

/-! ## The transport structure

A `DynamicsTransport` is a triple of bijections on ports, seams, and
faces together with a seam sign, constrained by two exact clauses on the
committed tables.  The endpoint clause states that the seam bijection is
the port bijection read through the committed seam table, with the sign
recording whether the committed smaller-to-larger orientation is
preserved.  The incidence clause states that the face bijection is the
seam bijection read through the committed signed face incidence.  The
structure is constrained, without stipulation: an inhabitant must
reproduce the committed tables exactly, and the generator inhabitants
below are verified row by row by kernel computation. -/

/-- A transport of the committed screen index data: bijections on ports,
seams, and faces with a seam orientation sign, tied to the committed
seam table and the committed signed face incidence by exact clauses. -/
structure DynamicsTransport where
  /-- The port bijection. -/
  port : Fin 12 ≃ Fin 12
  /-- The seam bijection. -/
  seam : Fin 30 ≃ Fin 30
  /-- The face bijection. -/
  face : Fin 20 ≃ Fin 20
  /-- The seam orientation sign. -/
  sign : Fin 30 → ℤ
  /-- Endpoint clause: the seam bijection is the port bijection read
  through the committed seam table, with the sign recording whether the
  committed orientation is preserved. -/
  endpoint_compat : ∀ e : Fin 30,
    (sign e = 1 ∧ seamLeft (seam e) = port (seamLeft e) ∧
        seamRight (seam e) = port (seamRight e)) ∨
      (sign e = -1 ∧ seamLeft (seam e) = port (seamRight e) ∧
        seamRight (seam e) = port (seamLeft e))
  /-- Incidence clause: the face bijection is the seam bijection read
  through the committed signed face incidence. -/
  incidence_compat : ∀ (f : Fin 20) (e : Fin 30),
    faceIncidenceZ f e * sign e = faceIncidenceZ (face f) (seam e)

/-- The seam sign of every inhabitant is a unit. -/
theorem DynamicsTransport.sign_unit (T : DynamicsTransport) (e : Fin 30) :
    T.sign e = 1 ∨ T.sign e = -1 := by
  rcases T.endpoint_compat e with ⟨hs, -, -⟩ | ⟨hs, -, -⟩
  · exact Or.inl hs
  · exact Or.inr hs

/-- The real seam sign squares to one. -/
theorem DynamicsTransport.sign_mul_self (T : DynamicsTransport) (e : Fin 30) :
    (T.sign e : ℝ) * (T.sign e : ℝ) = 1 := by
  rcases T.sign_unit e with hs | hs <;> rw [hs] <;> norm_num

noncomputable section

/-- Signed seam pullback of a seam configuration. -/
def DynamicsTransport.mapSeam (T : DynamicsTransport) (x : Fin 30 → ℝ) :
    Fin 30 → ℝ :=
  fun e ↦ (T.sign e : ℝ) * x (T.seam e)

/-- Port pullback of a port configuration. -/
def DynamicsTransport.mapPort (T : DynamicsTransport) (x : Fin 12 → ℝ) :
    Fin 12 → ℝ :=
  fun p ↦ x (T.port p)

/-- Face pullback of a face configuration. -/
def DynamicsTransport.mapFace (T : DynamicsTransport) (F : Fin 20 → ℝ) :
    Fin 20 → ℝ :=
  fun f ↦ F (T.face f)

/-- Stepwise seam pullback of a history. -/
def DynamicsTransport.mapSeamH (T : DynamicsTransport) (A : ℕ → Fin 30 → ℝ) :
    ℕ → Fin 30 → ℝ :=
  fun n ↦ T.mapSeam (A n)

/-- Stepwise port pullback of a history. -/
def DynamicsTransport.mapPortH (T : DynamicsTransport) (φ : ℕ → Fin 12 → ℝ) :
    ℕ → Fin 12 → ℝ :=
  fun n ↦ T.mapPort (φ n)

theorem DynamicsTransport.mapSeamH_apply (T : DynamicsTransport)
    (A : ℕ → Fin 30 → ℝ) (n : ℕ) : T.mapSeamH A n = T.mapSeam (A n) := rfl

theorem DynamicsTransport.mapPortH_apply (T : DynamicsTransport)
    (φ : ℕ → Fin 12 → ℝ) (n : ℕ) : T.mapPortH φ n = T.mapPort (φ n) := rfl

theorem DynamicsTransport.mapSeam_add (T : DynamicsTransport)
    (x y : Fin 30 → ℝ) : T.mapSeam (x + y) = T.mapSeam x + T.mapSeam y := by
  funext e
  simp only [DynamicsTransport.mapSeam, Pi.add_apply]
  ring

theorem DynamicsTransport.mapSeam_sub (T : DynamicsTransport)
    (x y : Fin 30 → ℝ) : T.mapSeam (x - y) = T.mapSeam x - T.mapSeam y := by
  funext e
  simp only [DynamicsTransport.mapSeam, Pi.sub_apply]
  ring

theorem DynamicsTransport.mapSeam_smul (T : DynamicsTransport) (c : ℝ)
    (x : Fin 30 → ℝ) : T.mapSeam (c • x) = c • T.mapSeam x := by
  funext e
  simp only [DynamicsTransport.mapSeam, Pi.smul_apply, smul_eq_mul]
  ring

theorem DynamicsTransport.mapSeam_zero (T : DynamicsTransport) :
    T.mapSeam 0 = 0 := by
  funext e
  simp [DynamicsTransport.mapSeam]

/-- The signed seam pullback is injective. -/
theorem DynamicsTransport.mapSeam_injective (T : DynamicsTransport) :
    Function.Injective T.mapSeam := by
  intro x y hxy
  funext e
  have h := congrFun hxy (T.seam.symm e)
  simp only [DynamicsTransport.mapSeam, Equiv.apply_symm_apply] at h
  rcases T.sign_unit (T.seam.symm e) with hs | hs <;> rw [hs] at h <;>
    push_cast at h <;> linarith

/-- The port pullback is injective. -/
theorem DynamicsTransport.mapPort_injective (T : DynamicsTransport) :
    Function.Injective T.mapPort := by
  intro x y hxy
  funext p
  have h := congrFun hxy (T.port.symm p)
  simpa only [DynamicsTransport.mapPort, Equiv.apply_symm_apply] using h

/-! ## Operator intertwining for every inhabitant -/

/-- The committed port coboundary intertwines the port and seam
pullbacks. -/
theorem coboundary_intertwine (T : DynamicsTransport) (φ : Fin 12 → ℝ) :
    realCoboundary (T.mapPort φ) = T.mapSeam (realCoboundary φ) := by
  funext e
  simp only [DynamicsTransport.mapSeam]
  rw [realCoboundary_apply, realCoboundary_apply]
  simp only [DynamicsTransport.mapPort]
  rcases T.endpoint_compat e with ⟨hs, hl, hr⟩ | ⟨hs, hl, hr⟩ <;>
    rw [hs, hl, hr] <;> push_cast <;> ring

/-- The committed face curvature intertwines the seam and face
pullbacks. -/
theorem curvature_intertwine (T : DynamicsTransport) (x : Fin 30 → ℝ) :
    faceCurvature (T.mapSeam x) = T.mapFace (faceCurvature x) := by
  funext f
  simp only [DynamicsTransport.mapFace]
  rw [faceCurvature_apply, faceCurvature_apply,
    ← Equiv.sum_comp T.seam (fun e ↦ faceIncidenceR (T.face f) e * x e)]
  refine Finset.sum_congr rfl fun e _ ↦ ?_
  simp only [DynamicsTransport.mapSeam]
  have hr : faceIncidenceR f e * (T.sign e : ℝ) =
      faceIncidenceR (T.face f) (T.seam e) := by
    unfold faceIncidenceR
    exact_mod_cast T.incidence_compat f e
  rw [← hr]
  ring

/-- The committed face codifferential intertwines the face and seam
pullbacks. -/
theorem codifferential_intertwine (T : DynamicsTransport) (F : Fin 20 → ℝ) :
    faceCodifferential (T.mapFace F) = T.mapSeam (faceCodifferential F) := by
  funext e
  simp only [DynamicsTransport.mapSeam]
  rw [faceCodifferential_apply, faceCodifferential_apply, Finset.mul_sum,
    ← Equiv.sum_comp T.face
      (fun f ↦ (T.sign e : ℝ) * (faceIncidenceR f (T.seam e) * F f))]
  refine Finset.sum_congr rfl fun f _ ↦ ?_
  simp only [DynamicsTransport.mapFace]
  have hr : faceIncidenceR f e * (T.sign e : ℝ) =
      faceIncidenceR (T.face f) (T.seam e) := by
    unfold faceIncidenceR
    exact_mod_cast T.incidence_compat f e
  have hs := T.sign_mul_self e
  rw [← hr]
  linear_combination (-(faceIncidenceR f e * F (T.face f))) * hs

/-- The committed incidence boundary intertwines the seam and port
pullbacks. -/
theorem boundary_intertwine (T : DynamicsTransport) (x : Fin 30 → ℝ) :
    realBoundary (T.mapSeam x) = T.mapPort (realBoundary x) := by
  funext p
  simp only [DynamicsTransport.mapPort]
  rw [realBoundary_apply, realBoundary_apply,
    ← Equiv.sum_comp T.seam (fun e ↦
      (if T.port p = seamRight e then x e else 0) -
        (if T.port p = seamLeft e then x e else 0))]
  refine Finset.sum_congr rfl fun e _ ↦ ?_
  simp only [DynamicsTransport.mapSeam]
  rcases T.endpoint_compat e with ⟨hs, hl, hr⟩ | ⟨hs, hl, hr⟩ <;>
    rw [hs, hl, hr] <;>
    simp only [Equiv.apply_eq_iff_eq] <;>
    push_cast <;>
    split_ifs <;>
    ring

/-- The committed seam pairing is invariant under the seam pullback. -/
theorem seamInner_map (T : DynamicsTransport) (x y : Fin 30 → ℝ) :
    realSeamInner (T.mapSeam x) (T.mapSeam y) = realSeamInner x y := by
  unfold realSeamInner
  rw [← Equiv.sum_comp T.seam (fun e ↦ x e * y e)]
  refine Finset.sum_congr rfl fun e _ ↦ ?_
  simp only [DynamicsTransport.mapSeam]
  linear_combination (x (T.seam e) * y (T.seam e)) * T.sign_mul_self e

/-- The committed seam energy is invariant under the seam pullback. -/
theorem seamEnergy_map (T : DynamicsTransport) (x : Fin 30 → ℝ) :
    realSeamEnergy (T.mapSeam x) = realSeamEnergy x := by
  unfold realSeamEnergy
  rw [← Equiv.sum_comp T.seam (fun e ↦ x e ^ 2)]
  refine Finset.sum_congr rfl fun e _ ↦ ?_
  simp only [DynamicsTransport.mapSeam]
  linear_combination (x (T.seam e) ^ 2) * T.sign_mul_self e

/-- The committed face pairing is invariant under the face pullback. -/
theorem faceInner_map (T : DynamicsTransport) (F H : Fin 20 → ℝ) :
    faceInner (T.mapFace F) (T.mapFace H) = faceInner F H := by
  unfold faceInner
  rw [← Equiv.sum_comp T.face (fun f ↦ F f * H f)]
  exact Finset.sum_congr rfl fun f _ ↦ rfl

/-- The committed port pairing is invariant under the port pullback. -/
theorem portInner_map (T : DynamicsTransport) (x y : Fin 12 → ℝ) :
    realPortInner (T.mapPort x) (T.mapPort y) = realPortInner x y := by
  unfold realPortInner
  rw [← Equiv.sum_comp T.port (fun p ↦ x p * y p)]
  exact Finset.sum_congr rfl fun p _ ↦ rfl

/-! ## Evolution transport for every inhabitant, at every step -/

/-- The scaled electric field of a transported history is the seam
pullback of the electric field. -/
theorem electric_intertwine (T : DynamicsTransport) (h : ℝ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    electricFieldScaled h (T.mapSeamH A) (T.mapPortH φ) n =
      T.mapSeam (electricFieldScaled h A φ n) := by
  unfold electricFieldScaled
  rw [T.mapSeamH_apply, T.mapSeamH_apply, T.mapPortH_apply,
    coboundary_intertwine]
  funext e
  simp only [DynamicsTransport.mapSeam, Pi.sub_apply, Pi.neg_apply,
    Pi.smul_apply, smul_eq_mul]
  ring

/-- The magnetic field of a transported history is the face pullback of
the magnetic field. -/
theorem magnetic_intertwine (T : DynamicsTransport) (A : ℕ → Fin 30 → ℝ)
    (n : ℕ) :
    magneticField (T.mapSeamH A) n = T.mapFace (magneticField A n) := by
  unfold magneticField
  rw [T.mapSeamH_apply]
  exact curvature_intertwine T (A n)

/-- The scaled Ampere update holds for a transported history exactly when
it holds for the original history, at every step. -/
theorem ampere_intertwine (T : DynamicsTransport) (h : ℝ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ) :
    AmpereEvolutionScaled h (T.mapSeamH A) (T.mapPortH φ) (T.mapSeamH J) ↔
      AmpereEvolutionScaled h A φ J := by
  unfold AmpereEvolutionScaled
  constructor
  · intro H n
    have hn := H n
    rw [electric_intertwine, electric_intertwine, magnetic_intertwine,
      codifferential_intertwine, T.mapSeamH_apply, ← T.mapSeam_sub,
      ← T.mapSeam_sub, ← T.mapSeam_smul] at hn
    exact T.mapSeam_injective hn
  · intro H n
    rw [electric_intertwine, electric_intertwine, magnetic_intertwine,
      codifferential_intertwine, T.mapSeamH_apply, ← T.mapSeam_sub,
      ← T.mapSeam_sub, ← T.mapSeam_smul]
    exact congrArg T.mapSeam (H n)

/-- The Gauss constraint holds for a transported history with the
pulled-back charge exactly when it holds for the original. -/
theorem gauss_intertwine (T : DynamicsTransport) (h : ℝ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (ρ : Fin 12 → ℝ) (n : ℕ) :
    realBoundary (electricFieldScaled h (T.mapSeamH A) (T.mapPortH φ) n) =
        T.mapPort ρ ↔
      realBoundary (electricFieldScaled h A φ n) = ρ := by
  rw [electric_intertwine, boundary_intertwine]
  exact ⟨fun hH ↦ T.mapPort_injective hH, fun hH ↦ congrArg T.mapPort hH⟩

/-- The scaled staggered form is invariant under every transport. -/
theorem energy_intertwine (T : DynamicsTransport) (h : ℝ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    fieldEnergyScaled h (T.mapSeamH A) (T.mapPortH φ) n =
      fieldEnergyScaled h A φ n := by
  unfold fieldEnergyScaled
  rw [electric_intertwine, seamEnergy_map, magnetic_intertwine,
    magnetic_intertwine, faceInner_map]

/-- The transport commutes with the committed time-dependent gauge
transformation. -/
theorem gauge_transport_commute (T : DynamicsTransport)
    (A : ℕ → Fin 30 → ℝ) (χ : ℕ → Fin 12 → ℝ) :
    T.mapSeamH (gaugeTransformA A χ) =
      gaugeTransformA (T.mapSeamH A) (T.mapPortH χ) := by
  funext n
  show T.mapSeam (A n + realCoboundary (χ n)) =
    T.mapSeam (A n) + realCoboundary (T.mapPort (χ n))
  rw [T.mapSeam_add, coboundary_intertwine]

end

/-! ## The two generator transports

The committed generator permutations `genA` and `genB` of the candidate
module induce seam and face permutations through the committed tables.
The tables below were searched offline; every fidelity clause is
re-verified here by kernel computation, so the tables carry no trusted
content. -/

/-- Builds a `Fin` equivalence from a table and its inverse table. -/
def finEquivOfTables {n : ℕ} (f g : Fin n → Fin n)
    (hgf : ∀ x, g (f x) = x) (hfg : ∀ x, f (g x) = x) : Fin n ≃ Fin n :=
  ⟨f, g, hgf, hfg⟩

/-- Inverse table of the generator permutation `genA`. -/
def genAInv : Fin 12 → Fin 12 := ![1, 2, 0, 5, 3, 4, 7, 8, 6, 11, 9, 10]

/-- Inverse table of the generator permutation `genB`. -/
def genBInv : Fin 12 → Fin 12 := ![0, 3, 1, 6, 2, 7, 4, 9, 5, 10, 8, 11]

/-- Seam permutation induced by `genA` through the committed seam
table. -/
def seamPermA : Fin 30 → Fin 30 :=
  ![1, 5, 9, 10, 11, 0, 3, 2, 4, 7, 6, 8, 16, 15, 17, 19, 18, 20, 12, 13,
    14, 25, 26, 22, 21, 24, 23, 29, 27, 28]

/-- Inverse table of `seamPermA`. -/
def seamPermAInv : Fin 30 → Fin 30 :=
  ![5, 0, 7, 6, 8, 1, 10, 9, 11, 2, 3, 4, 18, 19, 20, 13, 12, 14, 16, 15,
    17, 24, 23, 26, 25, 21, 22, 28, 29, 27]

/-- Orientation sign of `seamPermA` against the committed
smaller-to-larger seam orientation. -/
def seamSignA : Fin 30 → ℤ :=
  ![-1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, -1, -1]

/-- Face permutation induced by `genA` through the committed face
table. -/
def facePermA : Fin 20 → Fin 20 :=
  ![0, 2, 5, 8, 9, 1, 4, 3, 7, 6, 13, 12, 15, 14, 10, 11, 18, 16, 17, 19]

/-- Inverse table of `facePermA`. -/
def facePermAInv : Fin 20 → Fin 20 :=
  ![0, 5, 1, 7, 6, 2, 9, 8, 3, 4, 14, 15, 11, 10, 13, 12, 17, 18, 16, 19]

/-- Seam permutation induced by `genB` through the committed seam
table. -/
def seamPermB : Fin 30 → Fin 30 :=
  ![1, 3, 0, 4, 2, 9, 5, 11, 10, 15, 16, 17, 6, 7, 8, 12, 22, 21, 19, 25,
    26, 13, 14, 18, 20, 27, 29, 23, 24, 28]

/-- Inverse table of `seamPermB`. -/
def seamPermBInv : Fin 30 → Fin 30 :=
  ![2, 0, 4, 1, 3, 6, 12, 13, 14, 5, 8, 7, 15, 21, 22, 9, 10, 11, 23, 18,
    24, 17, 16, 27, 28, 19, 20, 25, 29, 26]

/-- Orientation sign of `seamPermB`. -/
def seamSignB : Fin 30 → ℤ :=
  ![1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, -1, 1, 1, 1,
    1, 1, 1, -1, 1, 1, 1, 1]

/-- Face permutation induced by `genB` through the committed face
table. -/
def facePermB : Fin 20 → Fin 20 :=
  ![2, 0, 4, 1, 3, 8, 5, 9, 12, 13, 6, 7, 10, 16, 15, 18, 11, 14, 19, 17]

/-- Inverse table of `facePermB`. -/
def facePermBInv : Fin 20 → Fin 20 :=
  ![1, 3, 0, 4, 2, 6, 10, 11, 5, 7, 12, 16, 8, 9, 17, 14, 13, 19, 15, 18]

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- The transport induced by the first generator `genA`, with every
fidelity clause verified by kernel computation against the committed
tables. -/
def genATransport : DynamicsTransport where
  port := finEquivOfTables genA genAInv (by decide) (by decide)
  seam := finEquivOfTables seamPermA seamPermAInv (by decide) (by decide)
  face := finEquivOfTables facePermA facePermAInv (by decide) (by decide)
  sign := seamSignA
  endpoint_compat := by decide
  incidence_compat := by decide

set_option maxHeartbeats 4000000 in
set_option maxRecDepth 16384 in
/-- The transport induced by the second generator `genB`, with every
fidelity clause verified by kernel computation against the committed
tables. -/
def genBTransport : DynamicsTransport where
  port := finEquivOfTables genB genBInv (by decide) (by decide)
  seam := finEquivOfTables seamPermB seamPermBInv (by decide) (by decide)
  face := finEquivOfTables facePermB facePermBInv (by decide) (by decide)
  sign := seamSignB
  endpoint_compat := by decide
  incidence_compat := by decide

theorem genATransport_port (k : Fin 12) : genATransport.port k = genA k := rfl

theorem genBTransport_port (k : Fin 12) : genBTransport.port k = genB k := rfl

theorem genATransport_seam (e : Fin 30) :
    genATransport.seam e = seamPermA e := rfl

theorem genBTransport_seam (e : Fin 30) :
    genBTransport.seam e = seamPermB e := rfl

theorem genATransport_face (f : Fin 20) :
    genATransport.face f = facePermA f := rfl

theorem genBTransport_face (f : Fin 20) :
    genBTransport.face f = facePermB f := rfl

/-! ## Identity, composition, and word transports -/

/-- The identity transport. -/
def idTransport : DynamicsTransport where
  port := Equiv.refl _
  seam := Equiv.refl _
  face := Equiv.refl _
  sign := fun _ ↦ 1
  endpoint_compat := fun _ ↦ Or.inl ⟨rfl, rfl, rfl⟩
  incidence_compat := fun _ _ ↦ mul_one _

/-- Composition of transports: the first argument acts first.  The
compatibility clauses are closed under composition. -/
def DynamicsTransport.comp (T U : DynamicsTransport) : DynamicsTransport where
  port := T.port.trans U.port
  seam := T.seam.trans U.seam
  face := T.face.trans U.face
  sign := fun e ↦ T.sign e * U.sign (T.seam e)
  endpoint_compat := by
    intro e
    simp only [Equiv.trans_apply]
    rcases T.endpoint_compat e with ⟨hs, hl, hr⟩ | ⟨hs, hl, hr⟩ <;>
      rcases U.endpoint_compat (T.seam e) with ⟨hs', hl', hr'⟩ |
        ⟨hs', hl', hr'⟩
    · exact Or.inl ⟨by rw [hs, hs']; ring, by rw [hl', hl],
        by rw [hr', hr]⟩
    · exact Or.inr ⟨by rw [hs, hs']; ring, by rw [hl', hr],
        by rw [hr', hl]⟩
    · exact Or.inr ⟨by rw [hs, hs']; ring, by rw [hl', hl],
        by rw [hr', hr]⟩
    · exact Or.inl ⟨by rw [hs, hs']; ring, by rw [hl', hr],
        by rw [hr', hl]⟩
  incidence_compat := by
    intro f e
    simp only [Equiv.trans_apply]
    rw [← U.incidence_compat (T.face f) (T.seam e), ← T.incidence_compat f e]
    ring

/-- The transport of one generator letter. -/
def letterTransport (b : Bool) : DynamicsTransport :=
  if b then genATransport else genBTransport

/-- The transport of a word in the two generators, letters acting left to
right as in the committed `wordApply`. -/
def wordTransport : List Bool → DynamicsTransport
  | [] => idTransport
  | b :: w => (letterTransport b).comp (wordTransport w)

/-- The port component of a word transport is the committed word
action. -/
theorem wordTransport_port (w : List Bool) :
    ∀ k : Fin 12, (wordTransport w).port k = wordApply w k := by
  induction w with
  | nil => intro k; rfl
  | cons b w ih =>
    intro k
    cases b
    · show (wordTransport w).port (genBTransport.port k) = wordApply w (genB k)
      rw [ih, genBTransport_port]
    · show (wordTransport w).port (genATransport.port k) = wordApply w (genA k)
      rw [ih, genATransport_port]

/-- Every row of the committed sixty-permutation list is the port
component of a transport.  The words come from the committed word
certificate of the candidate module. -/
theorem listed_perm_transport (row : List ℕ)
    (hrow : row ∈ OPH.A5PortAction.perms) :
    ∃ T : DynamicsTransport,
      (List.finRange 12).map (fun k ↦ ((T.port k : Fin 12) : ℕ)) = row := by
  rw [perms_generated_by_two_rotations] at hrow
  obtain ⟨w, -, rfl⟩ := List.mem_map.mp hrow
  refine ⟨wordTransport w, ?_⟩
  refine List.map_congr_left fun k _ ↦ ?_
  exact congrArg Fin.val (wordTransport_port w k)

/-- Every row of the committed sixty-permutation list transports the
scaled evolution and its staggered form exactly. -/
theorem committed_perm_dynamics_compatibility (row : List ℕ)
    (hrow : row ∈ OPH.A5PortAction.perms) :
    ∃ T : DynamicsTransport,
      ((List.finRange 12).map (fun k ↦ ((T.port k : Fin 12) : ℕ)) = row) ∧
        (∀ (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
            (J : ℕ → Fin 30 → ℝ),
          AmpereEvolutionScaled h (T.mapSeamH A) (T.mapPortH φ)
              (T.mapSeamH J) ↔
            AmpereEvolutionScaled h A φ J) ∧
        ∀ (h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ),
          fieldEnergyScaled h (T.mapSeamH A) (T.mapPortH φ) n =
            fieldEnergyScaled h A φ n := by
  obtain ⟨T, hT⟩ := listed_perm_transport row hrow
  exact ⟨T, hT, fun h A φ J ↦ ampere_intertwine T h A φ J,
    fun h A φ n ↦ energy_intertwine T h A φ n⟩

/-- Two distinct transports inhabit the structure, so the compatibility
clauses select a family and force no identification of ports with
directions. -/
theorem transport_selection_not_forced :
    ∃ T U : DynamicsTransport, T.port 0 ≠ U.port 0 :=
  ⟨idTransport, genATransport, by decide⟩

/-! ## Geometric realization of the generator transports

The generator transports act on ports exactly as the committed rotation
matrices act on the candidate rays, with the declared factor two.  The
statements below transport that realization to seams (unordered
endpoint-ray pairs) and faces (ray multisets). -/

/-- Doubling in `ℤ[φ]`: multiplication by the constant `2` is the scaling
`zsc 2`. -/
theorem zmul_two_eq_zsc (x : Zphi) : zmul (2, 0) x = zsc 2 x := by
  obtain ⟨a, b⟩ := x
  apply Prod.ext
  · show 2 * a + 0 * b = ((2 : ℕ) : ℤ) * a
    push_cast
    ring
  · show 2 * b + 0 * a + 0 * b = ((2 : ℕ) : ℤ) * b
    push_cast
    ring

/-- Function-level form of the committed generator equivariance: `rotAZ`
carries the candidate ray of a port to twice the ray of its `genA`
image. -/
theorem rotA_ray_vsc (i : Fin 12) :
    matVecZ rotAZ (candidateRayZ i) = vsc 2 (candidateRayZ (genA i)) := by
  funext k
  rw [candidate_equivariant_rotA i k]
  exact zmul_two_eq_zsc _

/-- Function-level form for the second generator. -/
theorem rotB_ray_vsc (i : Fin 12) :
    matVecZ rotBZ (candidateRayZ i) = vsc 2 (candidateRayZ (genB i)) := by
  funext k
  rw [candidate_equivariant_rotB i k]
  exact zmul_two_eq_zsc _

/-- Seam-level geometric realization of the first generator transport:
the rotation carries the unordered endpoint-ray pair of a committed seam
to the pair of its image seam, with the declared factor two. -/
theorem genA_seam_geometric (e : Fin 30) :
    s(matVecZ rotAZ (candidateRayZ (seamLeft e)),
        matVecZ rotAZ (candidateRayZ (seamRight e))) =
      s(vsc 2 (candidateRayZ (seamLeft (genATransport.seam e))),
        vsc 2 (candidateRayZ (seamRight (genATransport.seam e)))) := by
  simp only [rotA_ray_vsc]
  rcases genATransport.endpoint_compat e with ⟨-, hl, hr⟩ | ⟨-, hl, hr⟩
  · rw [hl, hr, genATransport_port, genATransport_port]
  · rw [hl, hr]
    exact Sym2.eq_swap

/-- Seam-level geometric realization of the second generator
transport. -/
theorem genB_seam_geometric (e : Fin 30) :
    s(matVecZ rotBZ (candidateRayZ (seamLeft e)),
        matVecZ rotBZ (candidateRayZ (seamRight e))) =
      s(vsc 2 (candidateRayZ (seamLeft (genBTransport.seam e))),
        vsc 2 (candidateRayZ (seamRight (genBTransport.seam e)))) := by
  simp only [rotB_ray_vsc]
  rcases genBTransport.endpoint_compat e with ⟨-, hl, hr⟩ | ⟨-, hl, hr⟩
  · rw [hl, hr, genBTransport_port, genBTransport_port]
  · rw [hl, hr]
    exact Sym2.eq_swap

/-- The committed face-port triples transport along `genA` to the triples
of the image faces. -/
theorem facePortTriple_genA :
    ∀ f : Fin 20, (facePortTriple f).map genA = facePortTriple (facePermA f) := by
  decide

/-- The committed face-port triples transport along `genB`. -/
theorem facePortTriple_genB :
    ∀ f : Fin 20, (facePortTriple f).map genB = facePortTriple (facePermB f) := by
  decide

/-- Face-level geometric realization of the first generator transport:
the rotation carries the candidate face-image multiset to the image-face
multiset, with the declared factor two. -/
theorem genA_face_geometric (f : Fin 20) :
    (canonicalCandidate.faceImage f).map (matVecZ rotAZ) =
      (canonicalCandidate.faceImage (facePermA f)).map (vsc 2) := by
  show ((facePortTriple f).map candidateRayZ).map (matVecZ rotAZ) =
    ((facePortTriple (facePermA f)).map candidateRayZ).map (vsc 2)
  rw [← facePortTriple_genA f, Multiset.map_map, Multiset.map_map,
    Multiset.map_map]
  refine Multiset.map_congr rfl fun i _ ↦ ?_
  exact rotA_ray_vsc i

/-- Face-level geometric realization of the second generator
transport. -/
theorem genB_face_geometric (f : Fin 20) :
    (canonicalCandidate.faceImage f).map (matVecZ rotBZ) =
      (canonicalCandidate.faceImage (facePermB f)).map (vsc 2) := by
  show ((facePortTriple f).map candidateRayZ).map (matVecZ rotBZ) =
    ((facePortTriple (facePermB f)).map candidateRayZ).map (vsc 2)
  rw [← facePortTriple_genB f, Multiset.map_map, Multiset.map_map,
    Multiset.map_map]
  refine Multiset.map_congr rfl fun i _ ↦ ?_
  exact rotB_ray_vsc i

/-! ## Refinement compatibility at the committed barycentric level

The committed refinement acts on same-parent barycentric coordinates by
denominator multiplication, and the candidate module extends the carrier
to barycentric points of committed faces.  The theorems below transport
that extension: the matrix transport commutes with refinement by exact
scaling, and the generator transports carry the carrier extension of a
face to the carrier extension of its image face through an exact
coordinate alignment, so the two routes, refine then transport and
transport then refine, agree exactly. -/

theorem vsc_add_vec (n : ℕ) (u v : VecZ) : vsc n (u + v) = vsc n u + vsc n v := by
  funext k
  exact zsc_add n (u k) (v k)

theorem vsc_vsc_vec (m n : ℕ) (v : VecZ) : vsc m (vsc n v) = vsc (m * n) v := by
  funext k
  exact zsc_zsc m n (v k)

theorem vsc_swap (m n : ℕ) (v : VecZ) : vsc m (vsc n v) = vsc n (vsc m v) := by
  rw [vsc_vsc_vec, vsc_vsc_vec, Nat.mul_comm]

theorem matVec_add (M : Fin 3 → Fin 3 → Zphi) (u v : VecZ) :
    matVecZ M (u + v) = matVecZ M u + matVecZ M v := by
  funext r
  show zmul (M r 0) (u 0 + v 0) + zmul (M r 1) (u 1 + v 1) +
      zmul (M r 2) (u 2 + v 2) =
    (zmul (M r 0) (u 0) + zmul (M r 1) (u 1) + zmul (M r 2) (u 2)) +
      (zmul (M r 0) (v 0) + zmul (M r 1) (v 1) + zmul (M r 2) (v 2))
  rw [zmul_add, zmul_add, zmul_add]
  abel

theorem matVec_vsc (M : Fin 3 → Fin 3 → Zphi) (n : ℕ) (v : VecZ) :
    matVecZ M (vsc n v) = vsc n (matVecZ M v) := by
  funext r
  show zmul (M r 0) (zsc n (v 0)) + zmul (M r 1) (zsc n (v 1)) +
      zmul (M r 2) (zsc n (v 2)) =
    zsc n (zmul (M r 0) (v 0) + zmul (M r 1) (v 1) + zmul (M r 2) (v 2))
  rw [zmul_zsc, zmul_zsc, zmul_zsc, zsc_add, zsc_add]

/-- The matrix transport of the carrier extension commutes with the
committed refinement by exact scaling, for every matrix and every port
assignment. -/
theorem carrier_transport_refine (M : Fin 3 → Fin 3 → Zphi)
    (pm : Fin 12 → VecZ) (m : ℕ) (f : Fin 20) (x : Barycentric) :
    matVecZ M (baryCarrierZ pm f (OPH.DiscreteRefinement.refine m x)) =
      vsc m (matVecZ M (baryCarrierZ pm f x)) := by
  rw [baryCarrier_refine, matVec_vsc]

noncomputable section

/-- Algebraic `SameRay` corollary under nonnegative scaling.  It also holds
degenerately when `m = 0` or when the transported vector is zero; this theorem
alone therefore does not assert that either vector determines a geometric
ray. -/
theorem carrier_transport_refine_sameRay (M : Fin 3 → Fin 3 → Zphi)
    (pm : Fin 12 → VecZ) (m : ℕ) (f : Fin 20) (x : Barycentric) :
    SameRay ℝ
      (evalVec (matVecZ M (baryCarrierZ pm f x)))
      (evalVec (matVecZ M
        (baryCarrierZ pm f (OPH.DiscreteRefinement.refine m x)))) := by
  rw [carrier_transport_refine, evalVec_vsc]
  exact SameRay.sameRay_nonneg_smul_right _ (Nat.cast_nonneg m)

end

/-- The barycentric coordinates of a point, indexed. -/
def coordAt (x : Barycentric) : Fin 3 → ℕ := ![x.i, x.j, x.k]

/-- The ordered vertex ports of a committed face, indexed. -/
def vertAt (f : Fin 20) : Fin 3 → Fin 12 :=
  ![(faceVertices f).1, (faceVertices f).2.1, (faceVertices f).2.2]

/-- Indexed form of the carrier extension. -/
theorem baryCarrier_vert (pm : Fin 12 → VecZ) (f : Fin 20) (x : Barycentric) :
    baryCarrierZ pm f x =
      vsc (coordAt x 0) (pm (vertAt f 0)) + vsc (coordAt x 1) (pm (vertAt f 1)) +
        vsc (coordAt x 2) (pm (vertAt f 2)) := rfl

/-- Coordinate alignment of the first generator transport: position `j`
of the image face carries the source coordinate at position
`(alignA f).j`. -/
def alignA : Fin 20 → Fin 3 × Fin 3 × Fin 3 :=
  ![(1, 2, 0), (2, 0, 1), (1, 2, 0), (0, 1, 2), (0, 1, 2), (0, 1, 2),
    (0, 1, 2), (0, 1, 2), (0, 1, 2), (0, 1, 2), (0, 1, 2), (0, 1, 2),
    (0, 1, 2), (0, 1, 2), (0, 1, 2), (0, 1, 2), (0, 1, 2), (0, 1, 2),
    (0, 1, 2), (2, 0, 1)]

/-- Coordinate alignment of the second generator transport. -/
def alignB : Fin 20 → Fin 3 × Fin 3 × Fin 3 :=
  ![(0, 1, 2), (0, 1, 2), (0, 1, 2), (0, 1, 2), (0, 1, 2), (0, 1, 2),
    (1, 2, 0), (0, 1, 2), (0, 1, 2), (0, 1, 2), (0, 1, 2), (0, 1, 2),
    (2, 0, 1), (0, 1, 2), (1, 2, 0), (0, 1, 2), (0, 1, 2), (0, 1, 2),
    (2, 0, 1), (0, 1, 2)]

/-- The six coordinate permutations. -/
def permTriples : List (Fin 3 × Fin 3 × Fin 3) :=
  [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]

/-- The alignment table of the first generator matches the committed face
table exactly: every image-face vertex is the `genA` image of the aligned
source vertex. -/
theorem alignA_sound :
    ∀ f : Fin 20,
      vertAt (facePermA f) 0 = genA (vertAt f (alignA f).1) ∧
        vertAt (facePermA f) 1 = genA (vertAt f (alignA f).2.1) ∧
          vertAt (facePermA f) 2 = genA (vertAt f (alignA f).2.2) := by
  decide

/-- The alignment table of the second generator matches the committed
face table exactly. -/
theorem alignB_sound :
    ∀ f : Fin 20,
      vertAt (facePermB f) 0 = genB (vertAt f (alignB f).1) ∧
        vertAt (facePermB f) 1 = genB (vertAt f (alignB f).2.1) ∧
          vertAt (facePermB f) 2 = genB (vertAt f (alignB f).2.2) := by
  decide

/-- Every alignment row of the first generator is a coordinate
permutation. -/
theorem alignA_perm : ∀ f : Fin 20, alignA f ∈ permTriples := by decide

/-- Every alignment row of the second generator is a coordinate
permutation. -/
theorem alignB_perm : ∀ f : Fin 20, alignB f ∈ permTriples := by decide

/-- Barycentric coordinate transport of the first generator. -/
def baryAlignA (f : Fin 20) (x : Barycentric) : Barycentric :=
  ⟨coordAt x (alignA f).1, coordAt x (alignA f).2.1, coordAt x (alignA f).2.2⟩

/-- Barycentric coordinate transport of the second generator. -/
def baryAlignB (f : Fin 20) (x : Barycentric) : Barycentric :=
  ⟨coordAt x (alignB f).1, coordAt x (alignB f).2.1, coordAt x (alignB f).2.2⟩

theorem baryAlignA_c0 (f : Fin 20) (x : Barycentric) :
    coordAt (baryAlignA f x) 0 = coordAt x (alignA f).1 := rfl

theorem baryAlignA_c1 (f : Fin 20) (x : Barycentric) :
    coordAt (baryAlignA f x) 1 = coordAt x (alignA f).2.1 := rfl

theorem baryAlignA_c2 (f : Fin 20) (x : Barycentric) :
    coordAt (baryAlignA f x) 2 = coordAt x (alignA f).2.2 := rfl

theorem baryAlignB_c0 (f : Fin 20) (x : Barycentric) :
    coordAt (baryAlignB f x) 0 = coordAt x (alignB f).1 := rfl

theorem baryAlignB_c1 (f : Fin 20) (x : Barycentric) :
    coordAt (baryAlignB f x) 1 = coordAt x (alignB f).2.1 := rfl

theorem baryAlignB_c2 (f : Fin 20) (x : Barycentric) :
    coordAt (baryAlignB f x) 2 = coordAt x (alignB f).2.2 := rfl

theorem coordAt_refine (m : ℕ) (x : Barycentric) (j : Fin 3) :
    coordAt (OPH.DiscreteRefinement.refine m x) j = m * coordAt x j := by
  fin_cases j <;> rfl

/-- The coordinate transport commutes with the committed refinement. -/
theorem baryAlignA_refine (f : Fin 20) (m : ℕ) (x : Barycentric) :
    baryAlignA f (OPH.DiscreteRefinement.refine m x) =
      OPH.DiscreteRefinement.refine m (baryAlignA f x) := by
  show (⟨coordAt (OPH.DiscreteRefinement.refine m x) (alignA f).1,
      coordAt (OPH.DiscreteRefinement.refine m x) (alignA f).2.1,
      coordAt (OPH.DiscreteRefinement.refine m x) (alignA f).2.2⟩ :
        Barycentric) = _
  rw [coordAt_refine, coordAt_refine, coordAt_refine]
  rfl

/-- The coordinate transport of the second generator commutes with the
committed refinement. -/
theorem baryAlignB_refine (f : Fin 20) (m : ℕ) (x : Barycentric) :
    baryAlignB f (OPH.DiscreteRefinement.refine m x) =
      OPH.DiscreteRefinement.refine m (baryAlignB f x) := by
  show (⟨coordAt (OPH.DiscreteRefinement.refine m x) (alignB f).1,
      coordAt (OPH.DiscreteRefinement.refine m x) (alignB f).2.1,
      coordAt (OPH.DiscreteRefinement.refine m x) (alignB f).2.2⟩ :
        Barycentric) = _
  rw [coordAt_refine, coordAt_refine, coordAt_refine]
  rfl

/-- Aligned face-level transport of the carrier extension for the first
generator: the rotation carries the carrier extension of a face to twice
the carrier extension of the image face at the aligned coordinates. -/
theorem bary_transport_A (f : Fin 20) (x : Barycentric) :
    matVecZ rotAZ (baryCarrierZ candidateRayZ f x) =
      vsc 2 (baryCarrierZ candidateRayZ (facePermA f) (baryAlignA f x)) := by
  obtain ⟨h0, h1, h2⟩ := alignA_sound f
  have hmem := alignA_perm f
  rw [baryCarrier_vert candidateRayZ f x,
    baryCarrier_vert candidateRayZ (facePermA f) (baryAlignA f x)]
  simp only [matVec_add, matVec_vsc, rotA_ray_vsc]
  rw [vsc_swap (coordAt x 0) 2, vsc_swap (coordAt x 1) 2,
    vsc_swap (coordAt x 2) 2, ← vsc_add_vec, ← vsc_add_vec]
  refine congrArg (vsc 2) ?_
  rw [h0, h1, h2, baryAlignA_c0, baryAlignA_c1, baryAlignA_c2]
  simp [permTriples] at hmem
  rcases hmem with h | h | h | h | h | h <;> simp only [h] <;> abel

/-- Aligned face-level transport of the carrier extension for the second
generator. -/
theorem bary_transport_B (f : Fin 20) (x : Barycentric) :
    matVecZ rotBZ (baryCarrierZ candidateRayZ f x) =
      vsc 2 (baryCarrierZ candidateRayZ (facePermB f) (baryAlignB f x)) := by
  obtain ⟨h0, h1, h2⟩ := alignB_sound f
  have hmem := alignB_perm f
  rw [baryCarrier_vert candidateRayZ f x,
    baryCarrier_vert candidateRayZ (facePermB f) (baryAlignB f x)]
  simp only [matVec_add, matVec_vsc, rotB_ray_vsc]
  rw [vsc_swap (coordAt x 0) 2, vsc_swap (coordAt x 1) 2,
    vsc_swap (coordAt x 2) 2, ← vsc_add_vec, ← vsc_add_vec]
  refine congrArg (vsc 2) ?_
  rw [h0, h1, h2, baryAlignB_c0, baryAlignB_c1, baryAlignB_c2]
  simp [permTriples] at hmem
  rcases hmem with h | h | h | h | h | h <;> simp only [h] <;> abel

/-- The two routes agree exactly for the first generator: refining and
then transporting equals transporting and then refining. -/
theorem bary_transport_A_routes (f : Fin 20) (m : ℕ) (x : Barycentric) :
    matVecZ rotAZ
        (baryCarrierZ candidateRayZ f (OPH.DiscreteRefinement.refine m x)) =
      vsc 2 (baryCarrierZ candidateRayZ (facePermA f)
        (OPH.DiscreteRefinement.refine m (baryAlignA f x))) := by
  rw [bary_transport_A f (OPH.DiscreteRefinement.refine m x), baryAlignA_refine]

/-- The two routes agree exactly for the second generator. -/
theorem bary_transport_B_routes (f : Fin 20) (m : ℕ) (x : Barycentric) :
    matVecZ rotBZ
        (baryCarrierZ candidateRayZ f (OPH.DiscreteRefinement.refine m x)) =
      vsc 2 (baryCarrierZ candidateRayZ (facePermB f)
        (OPH.DiscreteRefinement.refine m (baryAlignB f x))) := by
  rw [bary_transport_B f (OPH.DiscreteRefinement.refine m x), baryAlignB_refine]

/-! ## Certified-step receipts transport

The strongest committed evolution object connected here is the certified
scaled-step instrument at the committed step `4/5`: every transport
carries its zero-current solutions to zero-current solutions with the
same staggered form and the same rational electric energy bound. -/

noncomputable section

/-- At the committed certified step `4/5`, every transport carries a
zero-current solution to a zero-current solution, preserves the staggered
form at every step count, and inherits the committed rational electric
energy bound against the original initial form. -/
theorem certified_step_transport (T : DynamicsTransport)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ)
    (hAmp : AmpereEvolutionScaled (4 / 5) A φ (fun _ ↦ 0)) :
    AmpereEvolutionScaled (4 / 5) (T.mapSeamH A) (T.mapPortH φ)
        (fun _ ↦ 0) ∧
      (∀ n : ℕ, fieldEnergyScaled (4 / 5) (T.mapSeamH A) (T.mapPortH φ) n =
        fieldEnergyScaled (4 / 5) A φ n) ∧
      ∀ n : ℕ,
        realSeamEnergy
            (electricFieldScaled (4 / 5) (T.mapSeamH A) (T.mapPortH φ) n) ≤
          (25 / 2) * fieldEnergyScaled (4 / 5) A φ 0 := by
  have hzero : T.mapSeamH (fun _ ↦ (0 : Fin 30 → ℝ)) =
      fun _ ↦ (0 : Fin 30 → ℝ) := by
    funext n
    exact T.mapSeam_zero
  have hAmp' : AmpereEvolutionScaled (4 / 5) (T.mapSeamH A) (T.mapPortH φ)
      (fun _ ↦ 0) := by
    rw [← hzero]
    exact (ampere_intertwine T (4 / 5) A φ (fun _ ↦ 0)).mpr hAmp
  refine ⟨hAmp', fun n ↦ energy_intertwine T (4 / 5) A φ n, fun n ↦ ?_⟩
  have hb := (certifiedStep_energy_bounds_rational (T.mapSeamH A)
    (T.mapPortH φ) hAmp' n).1
  rw [energy_intertwine] at hb
  exact hb

end

/-! ## A gauge-invariant interaction candidate

The coupling below pairs the worldline increment of the canonical
Lorentz module with the seam electric field through the embedded seam
direction vectors of the candidate.  The coupling shape is declared, and
the relative normalization is a declared free constant: gauge invariance
holds at every value, the zero value recovers the decoupled direct-sum
regime of the committed joint action, and an explicit configuration
separates distinct values, so nothing proved here forces the
normalization. -/

/-- The exact seam direction vector of the candidate embedding: the ray
of the larger endpoint minus the ray of the smaller endpoint, aligned
with the committed smaller-to-larger seam orientation. -/
def seamVectorZ (e : Fin 30) : VecZ :=
  fun k ↦ zsub (candidateRayZ (seamRight e) k) (candidateRayZ (seamLeft e) k)

noncomputable section

/-- The real seam direction vector. -/
def seamVector (e : Fin 30) : Spatial := evalVec (seamVectorZ e)

/-- The interaction candidate: the declared relative normalization `κ`
times the window sum over steps and seams of the seam electric field
paired with the worldline increment through the embedded seam
direction. -/
def interactionCandidate (κ h : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (x : ℕ → Herm2) (N : ℕ) : ℝ :=
  κ * ∑ n ∈ Finset.range N, ∑ e : Fin 30,
      electricFieldScaled h A φ n e *
        spatialDot (seamVector e) ((x (n + 1)).2 - (x n).2)

/-- Exact gauge invariance of the interaction candidate, at every
declared relative normalization. -/
theorem interaction_gauge_invariant (κ h : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ χ : ℕ → Fin 12 → ℝ) (x : ℕ → Herm2) (N : ℕ) :
    interactionCandidate κ h (gaugeTransformA A χ)
        (gaugeTransformPhiScaled h φ χ) x N =
      interactionCandidate κ h A φ x N := by
  unfold interactionCandidate
  congr 1
  refine Finset.sum_congr rfl fun n _ ↦ Finset.sum_congr rfl fun e _ ↦ ?_
  rw [electricFieldScaled_gauge_invariant]

/-- Zero-coupling regression limit: at zero relative normalization the
interaction candidate vanishes identically, recovering the decoupled
direct-sum regime of the committed joint action. -/
theorem interaction_zero_coupling (h : ℝ) (A : ℕ → Fin 30 → ℝ)
    (φ : ℕ → Fin 12 → ℝ) (x : ℕ → Herm2) (N : ℕ) :
    interactionCandidate 0 h A φ x N = 0 := by
  unfold interactionCandidate
  rw [zero_mul]

/-- A pulse history: one seam potential switched on after the first
step. -/
def pulseHistory : ℕ → Fin 30 → ℝ :=
  fun n e ↦ if n = 0 then 0 else if e = 0 then -1 else 0

/-- A step worldline: the spatial part moves by the embedded direction of
seam zero across the first step. -/
def stepWorldline : ℕ → Herm2 :=
  fun n ↦ ((0 : ℝ), if n = 0 then (0 : Spatial) else seamVector 0)

/-- The electric field of the pulse history at the first step is the
indicator of seam zero. -/
theorem pulse_electric :
    electricFieldScaled 1 pulseHistory (fun _ ↦ 0) 0 =
      fun e ↦ if e = 0 then (1 : ℝ) else 0 := by
  funext e
  show -(1⁻¹ • (pulseHistory 1 - pulseHistory 0)) e -
      realCoboundary (0 : Fin 12 → ℝ) e = if e = 0 then (1 : ℝ) else 0
  rw [map_zero]
  simp only [Pi.zero_apply, Pi.smul_apply, Pi.sub_apply,
    smul_eq_mul, inv_one, one_mul, sub_zero]
  unfold pulseHistory
  norm_num
  split_ifs <;> norm_num

/-- The exact self-pairing of the embedded direction of seam zero. -/
theorem seamVec0_selfdot :
    dotZ (seamVectorZ 0) (seamVectorZ 0) = ((4 : ℤ), (0 : ℤ)) := by decide

/-- The interaction candidate at normalization one on the pulse
configuration evaluates to four. -/
theorem interaction_witness_value :
    interactionCandidate 1 1 pulseHistory (fun _ ↦ 0) stepWorldline 1 = 4 := by
  unfold interactionCandidate
  rw [Finset.sum_range_one, one_mul]
  have hx : (stepWorldline (0 + 1)).2 - (stepWorldline 0).2 = seamVector 0 := by
    show (if (1 : ℕ) = 0 then (0 : Spatial) else seamVector 0) -
      (if (0 : ℕ) = 0 then (0 : Spatial) else seamVector 0) = seamVector 0
    norm_num
  rw [hx, pulse_electric]
  simp only [ite_mul, one_mul, zero_mul]
  rw [Finset.sum_ite_eq' Finset.univ (0 : Fin 30)
    (fun e ↦ spatialDot (seamVector e) (seamVector 0)),
    if_pos (Finset.mem_univ _)]
  unfold seamVector
  rw [← evalPhi_dotZ, seamVec0_selfdot]
  show ((4 : ℤ) : ℝ) + ((0 : ℤ) : ℝ) * Real.goldenRatio = 4
  push_cast
  ring

/-- The relative normalization is declared, not forced: gauge invariance
holds at every value, while the pulse configuration separates the values
one and zero. -/
theorem interaction_normalization_not_forced :
    (∀ (κ h : ℝ) (A : ℕ → Fin 30 → ℝ) (φ χ : ℕ → Fin 12 → ℝ)
        (x : ℕ → Herm2) (N : ℕ),
      interactionCandidate κ h (gaugeTransformA A χ)
          (gaugeTransformPhiScaled h φ χ) x N =
        interactionCandidate κ h A φ x N) ∧
      interactionCandidate 1 1 pulseHistory (fun _ ↦ 0) stepWorldline 1 ≠
        interactionCandidate 0 1 pulseHistory (fun _ ↦ 0) stepWorldline 1 := by
  refine ⟨interaction_gauge_invariant, ?_⟩
  rw [interaction_witness_value, interaction_zero_coupling]
  norm_num

end

end OPH.CarrierDynamicsCompatibility

/- Axiom audit: standard axioms only (`propext`, `Classical.choice`,
`Quot.sound`); no `sorry`, no `native_decide`, no project axiom. -/

#print axioms OPH.CarrierDynamicsCompatibility.coboundary_intertwine
#print axioms OPH.CarrierDynamicsCompatibility.curvature_intertwine
#print axioms OPH.CarrierDynamicsCompatibility.codifferential_intertwine
#print axioms OPH.CarrierDynamicsCompatibility.boundary_intertwine
#print axioms OPH.CarrierDynamicsCompatibility.seamInner_map
#print axioms OPH.CarrierDynamicsCompatibility.seamEnergy_map
#print axioms OPH.CarrierDynamicsCompatibility.faceInner_map
#print axioms OPH.CarrierDynamicsCompatibility.portInner_map
#print axioms OPH.CarrierDynamicsCompatibility.electric_intertwine
#print axioms OPH.CarrierDynamicsCompatibility.magnetic_intertwine
#print axioms OPH.CarrierDynamicsCompatibility.ampere_intertwine
#print axioms OPH.CarrierDynamicsCompatibility.gauss_intertwine
#print axioms OPH.CarrierDynamicsCompatibility.energy_intertwine
#print axioms OPH.CarrierDynamicsCompatibility.gauge_transport_commute
#print axioms OPH.CarrierDynamicsCompatibility.genATransport
#print axioms OPH.CarrierDynamicsCompatibility.genBTransport
#print axioms OPH.CarrierDynamicsCompatibility.wordTransport_port
#print axioms OPH.CarrierDynamicsCompatibility.listed_perm_transport
#print axioms OPH.CarrierDynamicsCompatibility.committed_perm_dynamics_compatibility
#print axioms OPH.CarrierDynamicsCompatibility.transport_selection_not_forced
#print axioms OPH.CarrierDynamicsCompatibility.rotA_ray_vsc
#print axioms OPH.CarrierDynamicsCompatibility.rotB_ray_vsc
#print axioms OPH.CarrierDynamicsCompatibility.genA_seam_geometric
#print axioms OPH.CarrierDynamicsCompatibility.genB_seam_geometric
#print axioms OPH.CarrierDynamicsCompatibility.genA_face_geometric
#print axioms OPH.CarrierDynamicsCompatibility.genB_face_geometric
#print axioms OPH.CarrierDynamicsCompatibility.carrier_transport_refine
#print axioms OPH.CarrierDynamicsCompatibility.carrier_transport_refine_sameRay
#print axioms OPH.CarrierDynamicsCompatibility.alignA_sound
#print axioms OPH.CarrierDynamicsCompatibility.alignB_sound
#print axioms OPH.CarrierDynamicsCompatibility.bary_transport_A
#print axioms OPH.CarrierDynamicsCompatibility.bary_transport_B
#print axioms OPH.CarrierDynamicsCompatibility.bary_transport_A_routes
#print axioms OPH.CarrierDynamicsCompatibility.bary_transport_B_routes
#print axioms OPH.CarrierDynamicsCompatibility.certified_step_transport
#print axioms OPH.CarrierDynamicsCompatibility.interaction_gauge_invariant
#print axioms OPH.CarrierDynamicsCompatibility.interaction_zero_coupling
#print axioms OPH.CarrierDynamicsCompatibility.interaction_witness_value
#print axioms OPH.CarrierDynamicsCompatibility.interaction_normalization_not_forced
