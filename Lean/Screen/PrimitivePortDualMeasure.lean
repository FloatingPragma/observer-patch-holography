import PortFrameGram

open scoped BigOperators

namespace OPH.PrimitivePortDualMeasure

open OPH.PortFrameGram (adj)

/-!
# Finite primitive-port dual measure

This module checks the exact finite part of the issue #664 dual-measure
packet directly on the source-pinned twelve-port incidence table.  An
unordered base face is represented once by its increasing vertex triple.
The incidence table has twenty such triangular faces and each port belongs to
five of them.  Equal normalized face weights and the declared barycentric
one-third share per incident vertex then give port-dual weight `1/12`.

The result concerns a finite normalized support measure.  It does not identify
a dual sector with a physical pixel, interpret the weights as continuum area,
prove a refinement limit, select a physical length, or identify any support
radius with a translation hop.
-/

/-- Canonically ordered candidate for one triangular base face. -/
structure FaceTriple where
  first : Fin 12
  second : Fin 12
  third : Fin 12
deriving DecidableEq, Fintype

/-- A canonical base face is an increasing triple whose three port pairs are
    adjacent in the source-pinned incidence graph. -/
def isBaseFace (f : FaceTriple) : Bool :=
  f.first < f.second &&
    f.second < f.third &&
    adj f.first f.second &&
    adj f.first f.third &&
    adj f.second f.third

/-- The finite set of triangular base faces reconstructed from port
    incidence. -/
def baseFaces : Finset FaceTriple :=
  Finset.univ.filter fun f => isBaseFace f = true

/-- Membership of a port in a canonical triangular face. -/
def incident (p : Fin 12) (f : FaceTriple) : Bool :=
  p = f.first || p = f.second || p = f.third

/-- The source-pinned incidence graph has exactly twenty triangular faces. -/
theorem base_face_count : baseFaces.card = 20 := by
  decide

/-- Every source port is incident on exactly five base faces. -/
theorem incident_face_count_five (p : Fin 12) :
    (baseFaces.filter fun f => incident p f = true).card = 5 := by
  fin_cases p <;> decide

/-- Equal normalized weights on the twenty reconstructed base faces are
    exactly `1/20` on every face. -/
theorem equal_normalized_face_weight_is_one_twentieth
    (faceWeight : FaceTriple → ℚ)
    (hconstant : ∀ f ∈ baseFaces, ∀ g ∈ baseFaces,
      faceWeight f = faceWeight g)
    (hnormalized : ∑ f ∈ baseFaces, faceWeight f = 1) :
    ∀ f ∈ baseFaces, faceWeight f = 1 / 20 := by
  intro f hf
  have hsum : ∑ g ∈ baseFaces, faceWeight g =
      baseFaces.card * faceWeight f := by
    calc
      ∑ g ∈ baseFaces, faceWeight g =
          ∑ _g ∈ baseFaces, faceWeight f := by
        apply Finset.sum_congr rfl
        intro g hg
        exact hconstant g hg f hf
      _ = baseFaces.card * faceWeight f := by simp
  rw [hsum, base_face_count] at hnormalized
  norm_num at hnormalized ⊢
  linarith

/-- Port-dual weight at a declared per-vertex share: each incident triangular
    face contributes that share of its normalized face weight. -/
def portDualWeight (vertexShare : ℚ) (faceWeight : FaceTriple → ℚ)
    (p : Fin 12) : ℚ :=
  ∑ f ∈ baseFaces.filter fun f => incident p f = true,
    vertexShare * faceWeight f

/-- Equal normalized face weights, twenty base faces, five incident faces per
    port, and the barycentric one-third share imply exact normalized port-dual
    weight `1/12`. -/
theorem equal_face_weights_imply_port_dual_one_twelfth
    (vertexShare : ℚ) (faceWeight : FaceTriple → ℚ)
    (hbarycentric : vertexShare = 1 / 3)
    (hconstant : ∀ f ∈ baseFaces, ∀ g ∈ baseFaces,
      faceWeight f = faceWeight g)
    (hnormalized : ∑ f ∈ baseFaces, faceWeight f = 1) :
    ∀ p : Fin 12, portDualWeight vertexShare faceWeight p = 1 / 12 := by
  intro p
  have hface := equal_normalized_face_weight_is_one_twentieth
    faceWeight hconstant hnormalized
  have hsum : portDualWeight vertexShare faceWeight p =
      (baseFaces.filter fun f => incident p f = true).card *
        ((1 / 3 : ℚ) * (1 / 20 : ℚ)) := by
    unfold portDualWeight
    rw [hbarycentric]
    calc
      ∑ f ∈ baseFaces.filter (fun f => incident p f = true),
          (1 / 3) * faceWeight f =
          ∑ _f ∈ baseFaces.filter (fun f => incident p f = true),
            (1 / 3 : ℚ) * (1 / 20 : ℚ) := by
        apply Finset.sum_congr rfl
        intro f hf
        rw [hface f (Finset.mem_filter.mp hf).1]
      _ = (baseFaces.filter fun f => incident p f = true).card *
          ((1 / 3 : ℚ) * (1 / 20 : ℚ)) := by simp
  rw [hsum, incident_face_count_five]
  norm_num

#print axioms base_face_count
#print axioms incident_face_count_five
#print axioms equal_normalized_face_weight_is_one_twentieth
#print axioms equal_face_weights_imply_port_dual_one_twelfth

end OPH.PrimitivePortDualMeasure
