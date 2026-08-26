import CarrierModeOscillators
import A5PortAction

open scoped BigOperators Matrix

namespace OPH.CarrierModeEquivariance

open OPH.A5PortAction
open OPH.SeamCurrentCarrierQuotient
open OPH.DiscreteCoulombGreen
open OPH.LocalFaceMaxwellAction
open OPH.TemporalMaxwellEvolution
open OPH.ScaledMaxwellStability
open OPH.CarrierModeOscillators

/-!
# Icosahedral equivariance of the carrier's local operator and the irrep
reading of its mode spectrum

STATUS.  Kernel `decide` checks on the committed finite tables (the sixty
listed port permutations of `A5PortAction.perms`, the thirty seams of
`SeamCurrentCarrierQuotient`, the twenty oriented faces of
`LocalFaceMaxwellAction`, the five integer projector multiples of
`ScaledMaxwellStability`), transported to exact real linear algebra.  The
carrier, its seam orientation, its face orientation, the step `h`, and the
kinetic term are declared in the imported files.  No register row is
discharged.

WHAT IS PROVED.
1. Induced seam action.  For each listed row `p`, `seamAct p e` is the seam
   whose endpoints are the port images of the endpoints of `e`, with sign
   `+1` if the stored smaller-to-larger orientation of the image agrees
   with the image of the stored orientation of `e` and `-1` if reversed.
   The sign is `±1` on every seam of every row (`seamSign_pm_one`), the
   seam map is injective, hence bijective on `Fin 30`
   (`seamPerm_injective`, `seamEquiv`), and the endpoints of the image seam
   are the port images in the order fixed by the sign
   (`seamPerm_endpoints`).
2. Induced face action.  `faceAct p f` is the face whose vertex triple is
   the image triple of `f` up to cyclic order, with sign `+1` for the same
   cyclic order and `-1` for the reversed one.  On every face of every
   listed row the sign is `+1` (`faceSign_eq_one`): the listed rotations
   preserve the committed face orientation.  The face map is injective
   (`facePerm_injective`, `faceEquiv`).
3. Equivariance.  `faceIncidenceZ (p f) (p e) = sign_f sign_e
   faceIncidenceZ f e` over all sixty rows, twenty faces, thirty seams
   (`incidence_equivariant`, checked in full, no generating pair).  With
   the pullback actions `pullSeam p A e = sign_e A (p e)` and
   `pullFace p F f = sign_f F (p f)`, the curvature `C`, the codifferential
   `Cᵀ`, the local operator `CᵀC`, and the face normal operator `C Cᵀ`
   commute with the actions (`faceCurvature_pull`,
   `faceCodifferential_pull`, `localMaxwellOperator_pull`,
   `faceNormal_pull`); the port coboundary commutes with the port pullback
   (`realCoboundary_pull`).
4. Consequences.  Eigenvectors of `CᵀC` and of `C Cᵀ` map to eigenvectors
   of the same eigenvalue (`eigen_pull`, `faceEigen_pull`); the five
   committed projectors commute with the face action
   (`projZeroR_pull` to `projGoldenR_pull`, from the integer identity
   `projectors_equivariant`), so each projector image is invariant
   (`projector_images_invariant`).  Mode histories `cos (n θ) • v` and
   `sin (n θ) • v` transport to mode histories of the same eigenvalue
   (`cosHistory_pull`, `sinHistory_pull`, `cosHistory_pull_ampere`,
   `modeOscillator_pull`).  The scaled Ampere evolution is equivariant for
   every step, potential pair, and current (`ampereEvolutionScaled_pull`),
   in particular the zero-current temporal-gauge solutions map to such
   solutions (`ampereEvolutionScaled_pull_free`), and the staggered energy
   is invariant (`fieldEnergyScaled_pull`, from `realSeamEnergy_pull`,
   `faceInner_pull`).
5. Irrep reading at exact scope.  The face pullbacks obey the right action
   laws `pullFace (comp p q) = pullFace q ∘ pullFace p` and
   `pullFace id = id` (`pullFace_comp`, `pullFace_id`, from the table
   identities `facePerm_comp`, `facePerm_id` checked over all `60 × 60`
   pairs), and through the listed inverse `invPerm` the left action laws
   `leftFace_comp`, `leftFace_id`; the seam pullbacks obey the same laws
   with the sign cocycle (`pullSeam_comp`, `pullSeam_id`).  The image of
   `projZeroR` is the constant face vectors (`projZeroR_const`, trace `1`
   by `projector_traces`), and the constant vector is fixed pointwise by
   every row (`pullFace_const`).  Element orders are `1, 2, 3, 5` with
   class sizes `1, 15, 20, 24` (`elemOrder_spec`, `order_counts`).  The
   character theorem `projector_characters` (integer) and
   `projector_characters_R` (real trace of the action matrix against each
   projector) give, as functions of the element order, `1` on `projZero`,
   `5, 1, -1, 0` on `projTwo` (the five-dimensional irrep), `4, 0, 1, -1`
   on `projThree` and on `projFive` (the four-dimensional irrep), and
   `6, -2, 0, 1` on `projGolden` (the sum of the two three-dimensional
   irreps, whose values `φ` and `1 - φ` on the two order-five classes add
   to `1`).  The character norms `∑_g χ(g)²` are `60` on the `2`, `3`, `5`
   sectors and `120` on the golden sector (`character_norms`).
   Observation, outside the theorems: by the orthogonality relations these
   norms read as irreducibility of the `5`, `4`, `4` sectors and as two
   irreducible constituents of the golden sector, so the traces
   `1, 5, 4, 4, 6` of `CarrierModeOscillators.projector_traces` decompose as
   `1 ⊕ 5 ⊕ 4 ⊕ 4 ⊕ (3 ⊕ 3')` in the `A5` irreps `1, 3, 3', 4, 5`; the
   two four-dimensional eigenspaces (`λ = 3` and `λ = 5`) carry the same
   irrep, and the golden sector splits into the two three-dimensional
   irreps, one per eigenvalue `3 ± √5`, an inference from the character
   values of the two order-five classes that is stated here as observation.
   The checked golden values `6, -2, 0, 1` with norm `120` exclude, by the
   same orthogonality relations, the readings `1 ⊕ 5` (norm `120` but order
   two value `2`, not `-2`) and `3 ⊕ 3` (norm `240`), so `3 ⊕ 3'` is the
   only reading consistent with the checked values; which eigenvalue
   carries `3` and which `3'` is the part stated as observation.

PRIOR WORK.  `Lean/Geometry/CarrierDynamicsCompatibility.lean` proves,
for every inhabitant of its `DynamicsTransport` structure, the coboundary,
curvature, codifferential, and boundary intertwining
(`coboundary_intertwine`, `curvature_intertwine`,
`codifferential_intertwine`, `boundary_intertwine`), invariance of the
seam, face, and port pairings and of the seam energy (`seamInner_map`,
`faceInner_map`, `portInner_map`, `seamEnergy_map`), transport of the
electric and magnetic fields, of the scaled Ampere update in both
directions for every step and current, of the Gauss constraint, and of
the staggered energy (`electric_intertwine`, `magnetic_intertwine`,
`ampere_intertwine`, `gauss_intertwine`, `energy_intertwine`), with
sixty-row coverage through the word certificate
(`listed_perm_transport`, `committed_perm_dynamics_compatibility`).  The
statements of sections 3 and 4 below (`realCoboundary_pull`,
`faceCurvature_pull`, `faceCodifferential_pull`, `realSeamEnergy_pull`,
`faceInner_pull`, `electricFieldScaled_pull`, `magneticField_pull`,
`ampereEvolutionScaled_pull`, `fieldEnergyScaled_pull`) restate that
content for the explicit total action `seamAct`, `faceAct` of each listed
row, obtained here by kernel `decide` on the row itself with no word
certificate, and are cited as such.  The content of this file beyond the
prior module is: the explicit seam and face tables for all sixty rows,
the face sign `+1` on every face (`faceSign_eq_one`), the signed-table
form of the incidence equivariance (`incidence_equivariant`), the
composition and inverse laws of the induced tables (`facePerm_comp`,
`pullFace_comp`, `leftFace_comp`, `pullSeam_comp`), projector
equivariance and eigenspace invariance (`projectors_equivariant`,
`projector_images_invariant`, `eigen_pull`, `faceEigen_pull`), mode
history transport (`cosHistory_pull`, `modeOscillator_pull`), element
orders and class counts (`elemOrder_spec`, `order_counts`), and the
characters and norms (`projector_characters_R`, `character_norms`).

ROWS TOUCHED (none discharged).  Source clock and duration row: the step
`h` is declared with no unit.  Physical spacetime attachment row: the
listed group is an incidence automorphism group of the declared carrier;
its identification with a physical rotation group is open.  Light-signal
row: the identification of a mode with a physical oscillation is open; the
theorems here sort the modes by symmetry type.  Coupled-action row: the
kinetic term is declared in `ScaledMaxwellStability`.  Laboratory clock
and energy calibration import: untouched.  Gravitation-route energy
identification: untouched.

NEGATIVES CITED.  Legendre non-identifiability
(`Lean/Variational/RealizedHistoryLegendreNoGo.lean`): realized histories
select no velocity curvature or Legendre map, so every Lagrangian shape is
a declared enrichment; the equivariance here is of the declared evolution.

CONVENTIONS.  Rows of `perms` act on ports by `app p k = p.getD k 0`;
`comp p q` is `p` after `q`.  Seams carry the committed smaller-to-larger
orientation; faces carry the committed cyclic orientation of
`faceVertices`.  Pullback actions are right actions: `pullFace (comp p q)
= pullFace q ∘ pullFace p`; the left action is `leftFace p = pullFace
(invPerm p)`.  The action matrix `faceActZ p` has entry `sign_i` at
`(i, p i)`.

FALSIFIER.  A row with a seam or face sign `0`, a face sign `-1`, an
incidence entry `faceIncidenceZ (p f) (p e)` differing from
`sign_f sign_e faceIncidenceZ f e`, a projector entry moved by a row, an
element order outside `{1, 2, 3, 5}`, or a trace off the stated character
value would make the corresponding kernel check fail.

Axiom audit.  The `#print axioms` lines at the end of the file show at most
`propext`, `Classical.choice`, and `Quot.sound`.
-/

/-! ## 1. Induced seam and face actions of the listed port permutations -/

/-- Port image of a port under a listed permutation, as a `Fin 12` (the
`% 12` is the identity on every listed row, `app_lt_twelve`). -/
def portPerm (p : List Nat) (k : Fin 12) : Fin 12 :=
  ⟨app p k.val % 12, Nat.mod_lt _ (by omega)⟩

/-- The seam whose endpoints are the images of the endpoints of seam `e`,
with the orientation sign: `+1` if the image seam's stored smaller-to-larger
orientation agrees with the image of the stored orientation of `e`, `-1`
if it is reversed.  The fallback `(e, 0)` is never taken on the listed
rows (`seamSign_pm_one`). -/
def seamAct (p : List Nat) (e : Fin 30) : Fin 30 × ℤ :=
  ((List.finRange 30).findSome? fun d =>
    if (seamLeft d).val = app p (seamLeft e).val ∧
        (seamRight d).val = app p (seamRight e).val then some (d, 1)
    else if (seamLeft d).val = app p (seamRight e).val ∧
        (seamRight d).val = app p (seamLeft e).val then some (d, -1)
    else none).getD (e, 0)

def seamPerm (p : List Nat) (e : Fin 30) : Fin 30 := (seamAct p e).1
def seamSign (p : List Nat) (e : Fin 30) : ℤ := (seamAct p e).2

/-- Cyclic-order test on vertex triples. -/
def cyclicEq (a b : Nat × Nat × Nat) : Bool :=
  a == b || a == (b.2.1, b.2.2, b.1) || a == (b.2.2, b.1, b.2.1)

/-- The image vertex triple of face `f` under `p`. -/
def faceImage (p : List Nat) (f : Fin 20) : Nat × Nat × Nat :=
  let t := faceVertices f
  (app p t.1.val, app p t.2.1.val, app p t.2.2.val)

/-- Nat-level vertex triple of a face. -/
def faceNat (g : Fin 20) : Nat × Nat × Nat :=
  let t := faceVertices g
  (t.1.val, t.2.1.val, t.2.2.val)

/-- The face whose vertex triple is the image triple of `f` up to cyclic
order, with sign `+1` if the cyclic orders agree and `-1` if the image
triple is the reversed cyclic order.  The fallback is never taken on the
listed rows (`faceSign_eq_one`). -/
def faceAct (p : List Nat) (f : Fin 20) : Fin 20 × ℤ :=
  ((List.finRange 20).findSome? fun g =>
    let t := faceNat g
    if cyclicEq (faceImage p f) t then some (g, 1)
    else if cyclicEq (faceImage p f) (t.1, t.2.2, t.2.1) then some (g, -1)
    else none).getD (f, 0)

def facePerm (p : List Nat) (f : Fin 20) : Fin 20 := (faceAct p f).1
def faceSign (p : List Nat) (f : Fin 20) : ℤ := (faceAct p f).2

/-- Every listed row sends every port index below twelve. -/
theorem app_lt_twelve : ∀ p ∈ perms, ∀ k : Fin 12, app p k.val < 12 := by
  decide +kernel

/-- **(1)** Every listed permutation maps seams to seams with a definite
orientation sign: the fallback is never taken. -/
theorem seamSign_pm_one : ∀ p ∈ perms, ∀ e : Fin 30, seamSign p e = 1 ∨ seamSign p e = -1 := by
  decide +kernel

/-- **(1)** The induced seam map of every listed permutation is injective,
hence bijective on `Fin 30`. -/
theorem seamPerm_injective : ∀ p ∈ perms, Function.Injective (seamPerm p) := by
  decide +kernel

/-- **(2)** Every listed permutation maps faces to faces with the fallback
never taken; the sign is `+1` on every face of every listed rotation: the
listed group preserves the committed face orientation. -/
theorem faceSign_eq_one : ∀ p ∈ perms, ∀ f : Fin 20, faceSign p f = 1 := by
  decide +kernel

/-- **(2)** The induced face map of every listed permutation is injective. -/
theorem facePerm_injective : ∀ p ∈ perms, Function.Injective (facePerm p) := by
  decide +kernel

/-! ## 2. Kernel-checked equivariance, composition law, orders, characters -/

/-- Lookup in a table built over `List.finRange`. -/
theorem getD_map_finRange {α : Type*} {n : ℕ} (f : Fin n → α) (d : α) (i : Fin n) :
    ((List.finRange n).map f).getD i.val d = f i := by
  have hl : i.val < (List.finRange n).length := by simp [i.isLt]
  rw [List.getD_eq_getElem?_getD, List.getElem?_map, List.getElem?_eq_getElem hl,
    List.getElem_finRange]
  simp

/-- **(3)** The signed face-seam incidence is invariant under the
simultaneous induced action: for every listed `p`, face `f`, seam `e`,
`faceIncidenceZ (p f) (p e) = sign_f(p) sign_e(p) faceIncidenceZ f e`.
Checked over all sixty listed rows, twenty faces, thirty seams. -/
theorem incidence_equivariant :
    ∀ p ∈ perms, ∀ f : Fin 20, ∀ e : Fin 30,
      faceIncidenceZ (facePerm p f) (seamPerm p e) =
        faceSign p f * seamSign p e * faceIncidenceZ f e := by
  decide +kernel

/-- The image seam's endpoints are the port images of the endpoints, in the
order fixed by the sign. -/
theorem seamPerm_endpoints :
    ∀ p ∈ perms, ∀ e : Fin 30,
      (seamSign p e = 1 → seamLeft (seamPerm p e) = portPerm p (seamLeft e) ∧
        seamRight (seamPerm p e) = portPerm p (seamRight e)) ∧
      (seamSign p e = -1 → seamLeft (seamPerm p e) = portPerm p (seamRight e) ∧
        seamRight (seamPerm p e) = portPerm p (seamLeft e)) := by
  decide +kernel

/-- The identity row induces the identity on seams with sign `+1`. -/
theorem seamAct_id : ∀ e : Fin 30, seamPerm (List.range 12) e = e ∧
    seamSign (List.range 12) e = 1 := by
  decide +kernel

/-- The identity row induces the identity on faces. -/
theorem facePerm_id : ∀ f : Fin 20, facePerm (List.range 12) f = f := by
  decide +kernel

/-- Closure of the listed set, in the form used below. -/
theorem comp_mem : ∀ p ∈ perms, ∀ q ∈ perms, comp p q ∈ perms := by
  decide +kernel

/-- Nat-level table of the induced face action (computed once per row
inside the kernel checks below). -/
def faceTab (p : List Nat) : List Nat :=
  (List.finRange 20).map fun f ↦ (facePerm p f).val

theorem faceTab_getD (p : List Nat) (f : Fin 20) :
    (faceTab p).getD f.val 0 = (facePerm p f).val :=
  getD_map_finRange _ _ f

/-- `app` of a composite is the composite of `app`, below twelve. -/
theorem app_comp (p q : List Nat) (k : Nat) (hk : k < 12) :
    app (comp p q) k = app p (app q k) := by
  have hl : k < ((List.range 12).map fun k ↦ app p (app q k)).length := by simp [hk]
  show ((List.range 12).map fun k ↦ app p (app q k)).getD k 0 = _
  rw [List.getD_eq_getElem?_getD, List.getElem?_eq_getElem hl, List.getElem_map,
    List.getElem_range]
  rfl

/-- The port action is functorial on listed rows. -/
theorem portPerm_comp (p : List Nat) (_hp : p ∈ perms) (q : List Nat) (hq : q ∈ perms)
    (k : Fin 12) : portPerm (comp p q) k = portPerm p (portPerm q k) := by
  apply Fin.ext
  simp only [portPerm, app_comp p q k.val k.isLt, Nat.mod_eq_of_lt (app_lt_twelve q hq k)]

/-- Endpoint form of the seam action, as a conditional on the sign. -/
theorem seamPerm_left : ∀ p ∈ perms, ∀ e : Fin 30, seamLeft (seamPerm p e) =
    if seamSign p e = 1 then portPerm p (seamLeft e) else portPerm p (seamRight e) := by
  decide +kernel

theorem seamPerm_right : ∀ p ∈ perms, ∀ e : Fin 30, seamRight (seamPerm p e) =
    if seamSign p e = 1 then portPerm p (seamRight e) else portPerm p (seamLeft e) := by
  decide +kernel

theorem seam_eq_of_endpoints (a b : Fin 30) (h1 : seamLeft a = seamLeft b)
    (h2 : seamRight a = seamRight b) : a = b :=
  seam_table_injective (Prod.ext h1 h2)

theorem seam_not_swapped (a b : Fin 30) (h1 : seamLeft a = seamRight b)
    (h2 : seamRight a = seamLeft b) : False := by
  have ha := (seam_table_sound a).1
  have hb := (seam_table_sound b).1
  rw [h1, h2] at ha
  exact lt_asymm ha hb

set_option linter.unusedSimpArgs false in
/-- **(5)** Composition law on seams: the seam map and sign of `comp p q`
(`p` after `q`) are the composite map and the product of signs.  Proved
from the endpoint form of the sixty rows, the functoriality of the port
action, and the uniqueness of a seam with given oriented endpoints, with
no `60 × 60` table. -/
theorem seamAct_comp (p : List Nat) (hp : p ∈ perms) (q : List Nat) (hq : q ∈ perms)
    (e : Fin 30) :
    seamPerm (comp p q) e = seamPerm p (seamPerm q e) ∧
      seamSign (comp p q) e = seamSign p (seamPerm q e) * seamSign q e := by
  have hr := comp_mem p hp q hq
  have L1 := seamPerm_left q hq e
  have R1 := seamPerm_right q hq e
  have L2 := seamPerm_left p hp (seamPerm q e)
  have R2 := seamPerm_right p hp (seamPerm q e)
  have L3 := seamPerm_left (comp p q) hr e
  have R3 := seamPerm_right (comp p q) hr e
  rw [portPerm_comp p hp q hq, portPerm_comp p hp q hq] at L3 R3
  rcases seamSign_pm_one q hq e with s1 | s1 <;>
  rcases seamSign_pm_one p hp (seamPerm q e) with s2 | s2 <;>
  rcases seamSign_pm_one (comp p q) hr e with s3 | s3 <;>
  norm_num [s1, s2, s3] at L1 R1 L2 R2 L3 R3 ⊢ <;>
  first
    | exact seam_eq_of_endpoints (seamPerm (comp p q) e) (seamPerm p (seamPerm q e))
        (by simp only [L1, R1, L2, R2, L3, R3]) (by simp only [L1, R1, L2, R2, L3, R3])
    | exact seam_not_swapped (seamPerm (comp p q) e) (seamPerm p (seamPerm q e))
        (by simp only [L1, R1, L2, R2, L3, R3]) (by simp only [L1, R1, L2, R2, L3, R3])

/-- A face is determined by its signed incidence row. -/
theorem face_row_injective : ∀ g g' : Fin 20,
    (∀ e : Fin 30, faceIncidenceZ g e = faceIncidenceZ g' e) → g = g' := by
  decide +kernel

/-- **(5)** Composition law on faces, from the seam law and equivariance:
the incidence rows of `(p q) f` and `p (q f)` agree on every seam. -/
theorem facePerm_comp (p : List Nat) (hp : p ∈ perms) (q : List Nat) (hq : q ∈ perms)
    (f : Fin 20) : facePerm (comp p q) f = facePerm p (facePerm q f) := by
  have hr := comp_mem p hp q hq
  apply face_row_injective
  intro e'
  obtain ⟨e1, rfl⟩ := (Finite.injective_iff_bijective.mp (seamPerm_injective p hp)).2 e'
  obtain ⟨e, rfl⟩ := (Finite.injective_iff_bijective.mp (seamPerm_injective q hq)).2 e1
  obtain ⟨hc1, hc2⟩ := seamAct_comp p hp q hq e
  have hL : faceIncidenceZ (facePerm (comp p q) f) (seamPerm p (seamPerm q e)) =
      seamSign (comp p q) e * faceIncidenceZ f e := by
    rw [← hc1, incidence_equivariant _ hr f e, faceSign_eq_one _ hr]; ring
  have hR : faceIncidenceZ (facePerm p (facePerm q f)) (seamPerm p (seamPerm q e)) =
      seamSign p (seamPerm q e) * seamSign q e * faceIncidenceZ f e := by
    rw [incidence_equivariant p hp, incidence_equivariant q hq, faceSign_eq_one p hp,
      faceSign_eq_one q hq]; ring
  rw [hL, hR, hc2]

/-- A listed inverse, by search; `inv_spec` checks it. -/
def invPerm (p : List Nat) : List Nat :=
  (perms.find? fun q => comp p q == List.range 12).getD (List.range 12)

theorem inv_spec : ∀ p ∈ perms, invPerm p ∈ perms ∧
    comp p (invPerm p) = List.range 12 ∧ comp (invPerm p) p = List.range 12 := by
  decide +kernel

theorem invPerm_id : invPerm (List.range 12) = List.range 12 := by
  decide +kernel

/-- Iterated composition. -/
def permPow (p : List Nat) : Nat → List Nat
  | 0 => List.range 12
  | n + 1 => comp p (permPow p n)

/-- The order of a listed element, by trial among `1, 2, 3, 5`. -/
def elemOrder (p : List Nat) : Nat :=
  if p = List.range 12 then 1
  else if permPow p 2 = List.range 12 then 2
  else if permPow p 3 = List.range 12 then 3
  else 5

theorem elemOrder_spec : ∀ p ∈ perms, permPow p (elemOrder p) = List.range 12 := by
  decide +kernel

/-- Class sizes by order: `1, 15, 20, 24` (the two order-five classes of
size twelve are not separated here; none of the characters below
distinguishes them). -/
theorem order_counts :
    (perms.filter fun p => elemOrder p = 1).length = 1 ∧
    (perms.filter fun p => elemOrder p = 2).length = 15 ∧
    (perms.filter fun p => elemOrder p = 3).length = 20 ∧
    (perms.filter fun p => elemOrder p = 5).length = 24 := by
  decide +kernel

/-- Integer trace of the induced signed face permutation against an integer
face matrix `Q`: `∑ f, sign_f · Q (p f) f`. -/
def traceZ (Q : Matrix (Fin 20) (Fin 20) ℤ) (p : List Nat) : ℤ :=
  ∑ f : Fin 20, faceSign p f * Q (facePerm p f) f

/-- The `A5` characters as functions of the element order (both order-five
classes share these values): trivial `1`; the four-dimensional irrep
`4, 0, 1, -1`; the five-dimensional irrep `5, 1, -1, 0`; the sum of the
two three-dimensional irreps `6, -2, 0, 1`. -/
def chiOne (_ : Nat) : ℤ := 1
def chiFour : Nat → ℤ
  | 1 => 4 | 2 => 0 | 3 => 1 | _ => -1
def chiFive : Nat → ℤ
  | 1 => 5 | 2 => 1 | 3 => -1 | _ => 0
def chiThreePair : Nat → ℤ
  | 1 => 6 | 2 => -2 | 3 => 0 | _ => 1

/-- **(5)** Character values on the five projector images: the trace of the
induced face action composed with each integer projector multiple equals
the normalisation times the `A5` character of the matching irrep, over all
sixty rows. -/
theorem projector_characters :
    ∀ p ∈ perms,
      traceZ projZeroZ p = 20 * chiOne (elemOrder p) ∧
      traceZ projTwoZ p = 12 * chiFive (elemOrder p) ∧
      traceZ projThreeZ p = 10 * chiFour (elemOrder p) ∧
      traceZ projFiveZ p = 30 * chiFour (elemOrder p) ∧
      traceZ projGoldenZ p = 10 * chiThreePair (elemOrder p) := by
  decide +kernel

/-- Character norms: `∑_g χ(g)²` equals `60` times the squared
normalisation for the three irreducible sectors and `120` times it for the
golden sector, the value `⟨χ, χ⟩ = 1` (respectively `2`) that reads as
irreducibility (respectively two irreducible constituents) by the
orthogonality relations, which are cited, not proved here. -/
theorem character_norms :
    (perms.map fun p => traceZ projTwoZ p ^ 2).sum = 60 * 12 ^ 2 ∧
    (perms.map fun p => traceZ projThreeZ p ^ 2).sum = 60 * 10 ^ 2 ∧
    (perms.map fun p => traceZ projFiveZ p ^ 2).sum = 60 * 30 ^ 2 ∧
    (perms.map fun p => traceZ projGoldenZ p ^ 2).sum = 120 * 10 ^ 2 := by
  decide +kernel

/-- List table of an integer face matrix and its lookup. -/
def matTab (Q : Matrix (Fin 20) (Fin 20) ℤ) : List (List ℤ) :=
  (List.finRange 20).map fun i ↦ (List.finRange 20).map fun j ↦ Q i j

def matLookup (T : List (List ℤ)) (a b : Nat) : ℤ := (T.getD a []).getD b 0

theorem matTab_lookup (Q : Matrix (Fin 20) (Fin 20) ℤ) (i j : Fin 20) :
    matLookup (matTab Q) i.val j.val = Q i j := by
  unfold matLookup matTab
  rw [getD_map_finRange, getD_map_finRange]

/-- Table form of the invariance `Q (p a) (p b) = Q a b`. -/
def projCheck (Q : Matrix (Fin 20) (Fin 20) ℤ) (p : List Nat) : Bool :=
  (List.range 20).all fun a ↦ (List.range 20).all fun b ↦
    matLookup (matTab Q) ((faceTab p).getD a 0) ((faceTab p).getD b 0) == matLookup (matTab Q) a b

theorem proj_checks : ∀ p ∈ perms,
    (projCheck projZeroZ p && projCheck projTwoZ p && projCheck projThreeZ p &&
      projCheck projFiveZ p && projCheck projGoldenZ p) = true := by
  decide +kernel

theorem projCheck_entry (Q : Matrix (Fin 20) (Fin 20) ℤ) (p : List Nat)
    (h : projCheck Q p = true) (i j : Fin 20) : Q (facePerm p i) (facePerm p j) = Q i j := by
  unfold projCheck at h
  rw [List.all_eq_true] at h
  have h1 := h i.val (List.mem_range.mpr i.isLt)
  rw [List.all_eq_true] at h1
  have h2 := h1 j.val (List.mem_range.mpr j.isLt)
  simp only [faceTab_getD, matTab_lookup, beq_iff_eq] at h2
  exact h2

/-- **(4)/(5)** Each integer projector multiple commutes with the induced
signed face permutation: `Q (p i) (p j) = sign_i sign_j Q i j`. -/
theorem projectors_equivariant (p : List Nat) (hp : p ∈ perms) (i j : Fin 20) :
      projZeroZ (facePerm p i) (facePerm p j) = faceSign p i * faceSign p j * projZeroZ i j ∧
      projTwoZ (facePerm p i) (facePerm p j) = faceSign p i * faceSign p j * projTwoZ i j ∧
      projThreeZ (facePerm p i) (facePerm p j) = faceSign p i * faceSign p j * projThreeZ i j ∧
      projFiveZ (facePerm p i) (facePerm p j) = faceSign p i * faceSign p j * projFiveZ i j ∧
      projGoldenZ (facePerm p i) (facePerm p j) =
        faceSign p i * faceSign p j * projGoldenZ i j := by
  have h := proj_checks p hp
  simp only [Bool.and_eq_true] at h
  obtain ⟨⟨⟨⟨h0, h2⟩, h3⟩, h5⟩, hg⟩ := h
  rw [faceSign_eq_one p hp i, faceSign_eq_one p hp j]
  simp only [one_mul]
  exact ⟨projCheck_entry _ p h0 i j, projCheck_entry _ p h2 i j, projCheck_entry _ p h3 i j,
    projCheck_entry _ p h5 i j, projCheck_entry _ p hg i j⟩

/-! ## 3. Real transport: signed permutation actions on seam and face fields -/

/-- Pullback action on seam fields: `(p · A) e = sign_e(p) · A (p e)`. -/
def pullSeam (p : List Nat) (A : Fin 30 → ℝ) : Fin 30 → ℝ :=
  fun e ↦ (seamSign p e : ℝ) * A (seamPerm p e)

/-- Pullback action on face fields: `(p · F) f = sign_f(p) · F (p f)`. -/
def pullFace (p : List Nat) (F : Fin 20 → ℝ) : Fin 20 → ℝ :=
  fun f ↦ (faceSign p f : ℝ) * F (facePerm p f)

/-- Pullback action on port fields. -/
def pullPort (p : List Nat) (φ : Fin 12 → ℝ) : Fin 12 → ℝ :=
  fun k ↦ φ (portPerm p k)

/-- The seam map of a listed row as an equivalence of `Fin 30`. -/
noncomputable def seamEquiv (p : List Nat) (hp : p ∈ perms) : Fin 30 ≃ Fin 30 :=
  Equiv.ofBijective (seamPerm p) (Finite.injective_iff_bijective.mp (seamPerm_injective p hp))

/-- The face map of a listed row as an equivalence of `Fin 20`. -/
noncomputable def faceEquiv (p : List Nat) (hp : p ∈ perms) : Fin 20 ≃ Fin 20 :=
  Equiv.ofBijective (facePerm p) (Finite.injective_iff_bijective.mp (facePerm_injective p hp))

theorem seamEquiv_apply (p : List Nat) (hp : p ∈ perms) (e : Fin 30) :
    seamEquiv p hp e = seamPerm p e := rfl

theorem faceEquiv_apply (p : List Nat) (hp : p ∈ perms) (f : Fin 20) :
    faceEquiv p hp f = facePerm p f := rfl

theorem seamSign_sq (p : List Nat) (hp : p ∈ perms) (e : Fin 30) :
    (seamSign p e : ℝ) * (seamSign p e : ℝ) = 1 := by
  rcases seamSign_pm_one p hp e with h | h <;> rw [h] <;> norm_num

theorem faceSign_sq (p : List Nat) (hp : p ∈ perms) (f : Fin 20) :
    (faceSign p f : ℝ) * (faceSign p f : ℝ) = 1 := by
  rw [faceSign_eq_one p hp f]; norm_num

theorem faceIncidenceR_equivariant (p : List Nat) (hp : p ∈ perms) (f : Fin 20) (e : Fin 30) :
    faceIncidenceR (facePerm p f) (seamPerm p e) =
      (faceSign p f : ℝ) * (seamSign p e : ℝ) * faceIncidenceR f e := by
  unfold faceIncidenceR
  rw [incidence_equivariant p hp f e]
  push_cast
  ring

/-- **(3)** The face curvature `C` commutes with the induced actions. -/
theorem faceCurvature_pull (p : List Nat) (hp : p ∈ perms) (A : Fin 30 → ℝ) :
    faceCurvature (pullSeam p A) = pullFace p (faceCurvature A) := by
  funext f
  simp only [pullFace, faceCurvature_apply, pullSeam]
  have hr : (∑ e' : Fin 30, faceIncidenceR (facePerm p f) e' * A e') =
      ∑ e : Fin 30, faceIncidenceR (facePerm p f) (seamPerm p e) * A (seamPerm p e) :=
    (Equiv.sum_comp (seamEquiv p hp) (fun e' ↦ faceIncidenceR (facePerm p f) e' * A e')).symm
  rw [hr, Finset.mul_sum]
  refine Finset.sum_congr rfl fun e _ ↦ ?_
  rw [faceIncidenceR_equivariant p hp f e]
  have h := faceSign_sq p hp f
  linear_combination (-(A (seamPerm p e) * (seamSign p e : ℝ) * faceIncidenceR f e)) * h

/-- **(3)** The face codifferential `Cᵀ` commutes with the induced actions. -/
theorem faceCodifferential_pull (p : List Nat) (hp : p ∈ perms) (F : Fin 20 → ℝ) :
    faceCodifferential (pullFace p F) = pullSeam p (faceCodifferential F) := by
  funext e
  simp only [pullSeam, faceCodifferential_apply, pullFace]
  have hr : (∑ f' : Fin 20, faceIncidenceR f' (seamPerm p e) * F f') =
      ∑ f : Fin 20, faceIncidenceR (facePerm p f) (seamPerm p e) * F (facePerm p f) :=
    (Equiv.sum_comp (faceEquiv p hp) (fun f' ↦ faceIncidenceR f' (seamPerm p e) * F f')).symm
  rw [hr, Finset.mul_sum]
  refine Finset.sum_congr rfl fun f _ ↦ ?_
  rw [faceIncidenceR_equivariant p hp f e]
  have h := seamSign_sq p hp e
  linear_combination (-(F (facePerm p f) * (faceSign p f : ℝ) * faceIncidenceR f e)) * h

/-- **(3)** The local operator `CᵀC` commutes with the induced seam action. -/
theorem localMaxwellOperator_pull (p : List Nat) (hp : p ∈ perms) (A : Fin 30 → ℝ) :
    localMaxwellOperator (pullSeam p A) = pullSeam p (localMaxwellOperator A) := by
  show faceCodifferential (faceCurvature (pullSeam p A)) =
    pullSeam p (faceCodifferential (faceCurvature A))
  rw [faceCurvature_pull p hp, faceCodifferential_pull p hp]

/-- **(3)** The face normal operator `C Cᵀ` commutes with the induced face action. -/
theorem faceNormal_pull (p : List Nat) (hp : p ∈ perms) (F : Fin 20 → ℝ) :
    faceNormalR.mulVec (pullFace p F) = pullFace p (faceNormalR.mulVec F) := by
  rw [← faceNormal_mulVec, ← faceNormal_mulVec, faceCodifferential_pull p hp,
    faceCurvature_pull p hp]

/-- The seam pullback as a linear map. -/
def pullSeamL (p : List Nat) : (Fin 30 → ℝ) →ₗ[ℝ] (Fin 30 → ℝ) where
  toFun := pullSeam p
  map_add' A B := by funext e; simp only [pullSeam, Pi.add_apply]; ring
  map_smul' a A := by
    funext e; simp only [pullSeam, Pi.smul_apply, smul_eq_mul, RingHom.id_apply]; ring

theorem pullSeamL_apply (p : List Nat) (A : Fin 30 → ℝ) : pullSeamL p A = pullSeam p A := rfl

theorem pullSeam_smul (p : List Nat) (a : ℝ) (A : Fin 30 → ℝ) :
    pullSeam p (a • A) = a • pullSeam p A := (pullSeamL p).map_smul a A

theorem pullSeam_zero (p : List Nat) : pullSeam p (0 : Fin 30 → ℝ) = 0 := (pullSeamL p).map_zero

/-- **(4)** Every eigenvector of `CᵀC` maps to an eigenvector of the same eigenvalue. -/
theorem eigen_pull (p : List Nat) (hp : p ∈ perms) (lam : ℝ) (v : Fin 30 → ℝ)
    (hv : localMaxwellOperator v = lam • v) :
    localMaxwellOperator (pullSeam p v) = lam • pullSeam p v := by
  rw [localMaxwellOperator_pull p hp, hv, pullSeam_smul]

/-- **(4)** Every eigenvector of `C Cᵀ` maps to an eigenvector of the same eigenvalue. -/
theorem faceEigen_pull (p : List Nat) (hp : p ∈ perms) (lam : ℝ) (w : Fin 20 → ℝ)
    (hw : faceNormalR.mulVec w = lam • w) :
    faceNormalR.mulVec (pullFace p w) = lam • pullFace p w := by
  rw [faceNormal_pull p hp, hw]
  funext f; simp only [pullFace, Pi.smul_apply, smul_eq_mul]; ring

/-- Mode histories transport: `p · (cos (n θ) • v) = cos (n θ) • (p · v)`. -/
theorem cosHistory_pull (p : List Nat) (θ : ℝ) (v : Fin 30 → ℝ) (n : ℕ) :
    pullSeam p (cosHistory θ v n) = cosHistory θ (pullSeam p v) n := by
  unfold cosHistory scalarHistory; rw [pullSeam_smul]

theorem sinHistory_pull (p : List Nat) (θ : ℝ) (v : Fin 30 → ℝ) (n : ℕ) :
    pullSeam p (sinHistory θ v n) = sinHistory θ (pullSeam p v) n := by
  unfold sinHistory scalarHistory; rw [pullSeam_smul]

/-- **(4)** The image of a mode history under the action is a mode history of
the same eigenvalue: it solves the zero-current scaled evolution. -/
theorem cosHistory_pull_ampere (p : List Nat) (hp : p ∈ perms) (h : ℝ) (hh : h ≠ 0)
    (lam : ℝ) (h0 : 0 ≤ h ^ 2 * lam) (h4 : h ^ 2 * lam ≤ 4) (v : Fin 30 → ℝ)
    (hv : localMaxwellOperator v = lam • v) :
    AmpereEvolutionScaled h (fun n ↦ pullSeam p (cosHistory (modeAngle h lam) v n))
      (fun _ ↦ 0) (fun _ ↦ 0) := by
  have : (fun n ↦ pullSeam p (cosHistory (modeAngle h lam) v n)) =
      cosHistory (modeAngle h lam) (pullSeam p v) := by
    funext n; exact cosHistory_pull p _ v n
  rw [this]
  exact cosHistory_ampere h hh lam h0 h4 _ (eigen_pull p hp lam v hv)

/-- **(4)** The full oscillator packet transports. -/
theorem modeOscillator_pull (p : List Nat) (hp : p ∈ perms) (h : ℝ) (hh : h ≠ 0)
    (lam : ℝ) (h0 : 0 ≤ h ^ 2 * lam) (h4 : h ^ 2 * lam ≤ 4) (v : Fin 30 → ℝ)
    (hv : localMaxwellOperator v = lam • v) : ModeOscillator h lam (pullSeam p v) :=
  modeOscillator h hh lam h0 h4 _ (eigen_pull p hp lam v hv)

/-! ## 4. Energy invariance, gauge and current transport, evolution equivariance -/

theorem realSeamEnergy_pull (p : List Nat) (hp : p ∈ perms) (A : Fin 30 → ℝ) :
    realSeamEnergy (pullSeam p A) = realSeamEnergy A := by
  unfold realSeamEnergy
  rw [← Equiv.sum_comp (seamEquiv p hp) (fun e ↦ A e ^ 2)]
  refine Finset.sum_congr rfl fun e _ ↦ ?_
  simp only [pullSeam, seamEquiv_apply]
  have h := seamSign_sq p hp e
  linear_combination (A (seamPerm p e) ^ 2) * h

theorem realSeamInner_pull (p : List Nat) (hp : p ∈ perms) (A B : Fin 30 → ℝ) :
    realSeamInner (pullSeam p A) (pullSeam p B) = realSeamInner A B := by
  unfold realSeamInner
  rw [← Equiv.sum_comp (seamEquiv p hp) (fun e ↦ A e * B e)]
  refine Finset.sum_congr rfl fun e _ ↦ ?_
  simp only [pullSeam, seamEquiv_apply]
  have h := seamSign_sq p hp e
  linear_combination (A (seamPerm p e) * B (seamPerm p e)) * h

theorem faceInner_pull (p : List Nat) (hp : p ∈ perms) (F G : Fin 20 → ℝ) :
    faceInner (pullFace p F) (pullFace p G) = faceInner F G := by
  unfold faceInner
  rw [← Equiv.sum_comp (faceEquiv p hp) (fun f ↦ F f * G f)]
  refine Finset.sum_congr rfl fun f _ ↦ ?_
  simp only [pullFace, faceEquiv_apply]
  have h := faceSign_sq p hp f
  linear_combination (F (facePerm p f) * G (facePerm p f)) * h

/-- The port coboundary `d` commutes with the induced actions. -/
theorem realCoboundary_pull (p : List Nat) (hp : p ∈ perms) (φ : Fin 12 → ℝ) :
    realCoboundary (pullPort p φ) = pullSeam p (realCoboundary φ) := by
  funext e
  simp only [realCoboundary_apply, pullPort, pullSeam]
  rcases seamSign_pm_one p hp e with h | h
  · obtain ⟨hl, hr⟩ := (seamPerm_endpoints p hp e).1 h
    rw [h, hl, hr]; push_cast; ring
  · obtain ⟨hl, hr⟩ := (seamPerm_endpoints p hp e).2 h
    rw [h, hl, hr]; push_cast; ring

theorem pullPort_zero (p : List Nat) : pullPort p (0 : Fin 12 → ℝ) = 0 := by
  funext k; rfl

theorem electricFieldScaled_pull (p : List Nat) (hp : p ∈ perms) (h : ℝ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    electricFieldScaled h (fun m ↦ pullSeam p (A m)) (fun m ↦ pullPort p (φ m)) n =
      pullSeam p (electricFieldScaled h A φ n) := by
  unfold electricFieldScaled
  show -(h⁻¹ • (pullSeam p (A (n + 1)) - pullSeam p (A n))) - realCoboundary (pullPort p (φ n)) = _
  rw [realCoboundary_pull p hp, ← pullSeamL_apply, ← pullSeamL_apply, ← pullSeamL_apply,
    ← pullSeamL_apply, map_sub, map_neg, map_smul, map_sub]

theorem magneticField_pull (p : List Nat) (hp : p ∈ perms) (A : ℕ → Fin 30 → ℝ) (n : ℕ) :
    magneticField (fun m ↦ pullSeam p (A m)) n = pullFace p (magneticField A n) := by
  unfold magneticField
  exact faceCurvature_pull p hp (A n)

/-- **(4)** The scaled Ampere evolution is equivariant: the induced action
carries a solution with potentials `(A, φ)` and current `J` to a solution
with the transported potentials and current. -/
theorem ampereEvolutionScaled_pull (p : List Nat) (hp : p ∈ perms) (h : ℝ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (J : ℕ → Fin 30 → ℝ)
    (hAmp : AmpereEvolutionScaled h A φ J) :
    AmpereEvolutionScaled h (fun m ↦ pullSeam p (A m)) (fun m ↦ pullPort p (φ m))
      (fun m ↦ pullSeam p (J m)) := by
  intro n
  rw [electricFieldScaled_pull p hp, electricFieldScaled_pull p hp, magneticField_pull p hp,
    faceCodifferential_pull p hp]
  simp only [← pullSeamL_apply, ← map_sub, ← map_smul]
  exact congrArg (pullSeamL p) (hAmp n)

/-- **(4)** Zero-current, temporal-gauge solutions map to zero-current,
temporal-gauge solutions. -/
theorem ampereEvolutionScaled_pull_free (p : List Nat) (hp : p ∈ perms) (h : ℝ)
    (A : ℕ → Fin 30 → ℝ) (hAmp : AmpereEvolutionScaled h A (fun _ ↦ 0) (fun _ ↦ 0)) :
    AmpereEvolutionScaled h (fun m ↦ pullSeam p (A m)) (fun _ ↦ 0) (fun _ ↦ 0) := by
  have h1 := ampereEvolutionScaled_pull p hp h A (fun _ ↦ 0) (fun _ ↦ 0) hAmp
  simpa only [pullPort_zero, pullSeam_zero] using h1

/-- **(4)** The staggered field energy is invariant under the induced action. -/
theorem fieldEnergyScaled_pull (p : List Nat) (hp : p ∈ perms) (h : ℝ)
    (A : ℕ → Fin 30 → ℝ) (φ : ℕ → Fin 12 → ℝ) (n : ℕ) :
    fieldEnergyScaled h (fun m ↦ pullSeam p (A m)) (fun m ↦ pullPort p (φ m)) n =
      fieldEnergyScaled h A φ n := by
  unfold fieldEnergyScaled
  rw [electricFieldScaled_pull p hp, magneticField_pull p hp, magneticField_pull p hp,
    realSeamEnergy_pull p hp, faceInner_pull p hp]

/-! ## 5. Projector images, the fixed vector, the group laws, the characters -/

/-- An integer face matrix satisfying the signed equivariance identity
commutes, after real scaling, with the induced face action. -/
theorem scaledProj_pull (p : List Nat) (hp : p ∈ perms) (Q : Matrix (Fin 20) (Fin 20) ℤ)
    (d : ℤ) (hQ : ∀ i j : Fin 20, Q (facePerm p i) (facePerm p j) =
      faceSign p i * faceSign p j * Q i j) (F : Fin 20 → ℝ) :
    (scaledProj Q d).mulVec (pullFace p F) = pullFace p ((scaledProj Q d).mulVec F) := by
  unfold scaledProj
  rw [Matrix.smul_mulVec, Matrix.smul_mulVec]
  have hc : (castZ Q).mulVec (pullFace p F) = pullFace p ((castZ Q).mulVec F) := by
    funext i
    simp only [Matrix.mulVec, dotProduct, castZ_apply, pullFace]
    have hr : (∑ j : Fin 20, (Q (facePerm p i) j : ℝ) * F j) =
        ∑ j : Fin 20, (Q (facePerm p i) (facePerm p j) : ℝ) * F (facePerm p j) :=
      (Equiv.sum_comp (faceEquiv p hp) (fun j ↦ (Q (facePerm p i) j : ℝ) * F j)).symm
    rw [hr, Finset.mul_sum]
    refine Finset.sum_congr rfl fun j _ ↦ ?_
    rw [hQ i j]; push_cast
    have h := faceSign_sq p hp i
    linear_combination (-((Q i j : ℝ) * (faceSign p j : ℝ) * F (facePerm p j))) * h
  rw [hc]
  funext f
  simp only [pullFace, Pi.smul_apply, smul_eq_mul]
  ring

/-- **(4)** The five committed projectors commute with the induced face action. -/
theorem projZeroR_pull (p : List Nat) (hp : p ∈ perms) (F : Fin 20 → ℝ) :
    projZeroR.mulVec (pullFace p F) = pullFace p (projZeroR.mulVec F) :=
  scaledProj_pull p hp projZeroZ 20 (fun i j ↦ (projectors_equivariant p hp i j).1) F

theorem projTwoR_pull (p : List Nat) (hp : p ∈ perms) (F : Fin 20 → ℝ) :
    projTwoR.mulVec (pullFace p F) = pullFace p (projTwoR.mulVec F) :=
  scaledProj_pull p hp projTwoZ 12 (fun i j ↦ (projectors_equivariant p hp i j).2.1) F

theorem projThreeR_pull (p : List Nat) (hp : p ∈ perms) (F : Fin 20 → ℝ) :
    projThreeR.mulVec (pullFace p F) = pullFace p (projThreeR.mulVec F) :=
  scaledProj_pull p hp projThreeZ 10 (fun i j ↦ (projectors_equivariant p hp i j).2.2.1) F

theorem projFiveR_pull (p : List Nat) (hp : p ∈ perms) (F : Fin 20 → ℝ) :
    projFiveR.mulVec (pullFace p F) = pullFace p (projFiveR.mulVec F) :=
  scaledProj_pull p hp projFiveZ 30 (fun i j ↦ (projectors_equivariant p hp i j).2.2.2.1) F

theorem projGoldenR_pull (p : List Nat) (hp : p ∈ perms) (F : Fin 20 → ℝ) :
    projGoldenR.mulVec (pullFace p F) = pullFace p (projGoldenR.mulVec F) :=
  scaledProj_pull p hp projGoldenZ 10 (fun i j ↦ (projectors_equivariant p hp i j).2.2.2.2) F

/-- **(4)** A projector image (its fixed-point set) is invariant under any
operator commuting with the projector. -/
theorem fixed_of_commute (P : Matrix (Fin 20) (Fin 20) ℝ) (T : (Fin 20 → ℝ) → (Fin 20 → ℝ))
    (hT : ∀ F, P.mulVec (T F) = T (P.mulVec F)) (F : Fin 20 → ℝ) (hF : P.mulVec F = F) :
    P.mulVec (T F) = T F := by
  rw [hT, hF]

/-- **(4)** The five projector images are invariant under every listed row. -/
theorem projector_images_invariant (p : List Nat) (hp : p ∈ perms) (F : Fin 20 → ℝ) :
    (projZeroR.mulVec F = F → projZeroR.mulVec (pullFace p F) = pullFace p F) ∧
    (projTwoR.mulVec F = F → projTwoR.mulVec (pullFace p F) = pullFace p F) ∧
    (projThreeR.mulVec F = F → projThreeR.mulVec (pullFace p F) = pullFace p F) ∧
    (projFiveR.mulVec F = F → projFiveR.mulVec (pullFace p F) = pullFace p F) ∧
    (projGoldenR.mulVec F = F → projGoldenR.mulVec (pullFace p F) = pullFace p F) :=
  ⟨fixed_of_commute _ _ (projZeroR_pull p hp) F, fixed_of_commute _ _ (projTwoR_pull p hp) F,
    fixed_of_commute _ _ (projThreeR_pull p hp) F, fixed_of_commute _ _ (projFiveR_pull p hp) F,
    fixed_of_commute _ _ (projGoldenR_pull p hp) F⟩

/-- Row sums of `projZeroZ` are `20`: it is the all-ones matrix. -/
theorem projZeroZ_row_sum : ∀ i : Fin 20, (∑ j : Fin 20, projZeroZ i j) = 20 := by
  decide +kernel

/-- **(5)** The constant face vector is the image of `projZeroR`. -/
theorem projZeroR_const (c : ℝ) : projZeroR.mulVec (fun _ ↦ c) = fun _ ↦ c := by
  funext i
  unfold projZeroR scaledProj
  rw [Matrix.smul_mulVec]
  simp only [Pi.smul_apply, smul_eq_mul, Matrix.mulVec, dotProduct, castZ_apply]
  rw [← Finset.sum_mul, ← Int.cast_sum, projZeroZ_row_sum i]
  push_cast; ring

/-- **(5)** The constant face vector is fixed pointwise by every listed row:
with the face signs all `+1`, `p · (fun _ ↦ c) = fun _ ↦ c`. -/
theorem pullFace_const (p : List Nat) (hp : p ∈ perms) (c : ℝ) :
    pullFace p (fun _ ↦ c) = fun _ ↦ c := by
  funext f
  simp only [pullFace, faceSign_eq_one p hp f]
  push_cast; ring

/-- **(5)** Right action law on face fields: `(p q) · F = q · (p · F)`. -/
theorem pullFace_comp (p : List Nat) (hp : p ∈ perms) (q : List Nat) (hq : q ∈ perms)
    (F : Fin 20 → ℝ) : pullFace (comp p q) F = pullFace q (pullFace p F) := by
  funext f
  simp only [pullFace]
  rw [facePerm_comp p hp q hq f, faceSign_eq_one _ (comp_mem p hp q hq) f,
    faceSign_eq_one p hp, faceSign_eq_one q hq]
  push_cast; ring

/-- **(5)** Right action law on seam fields, with the sign cocycle. -/
theorem pullSeam_comp (p : List Nat) (hp : p ∈ perms) (q : List Nat) (hq : q ∈ perms)
    (A : Fin 30 → ℝ) : pullSeam (comp p q) A = pullSeam q (pullSeam p A) := by
  funext e
  simp only [pullSeam]
  obtain ⟨h1, h2⟩ := seamAct_comp p hp q hq e
  rw [h1, h2]
  push_cast; ring

theorem pullFace_id (F : Fin 20 → ℝ) : pullFace (List.range 12) F = F := by
  funext f
  simp only [pullFace, facePerm_id f, faceSign_eq_one _ id_mem f]
  push_cast; ring

theorem pullSeam_id (A : Fin 30 → ℝ) : pullSeam (List.range 12) A = A := by
  funext e
  obtain ⟨h1, h2⟩ := seamAct_id e
  simp only [pullSeam, h1, h2]
  push_cast; ring

/-- **(5)** Left action on face fields through the listed inverse:
`leftFace p := pullFace (invPerm p)` obeys `leftFace (p q) = leftFace p ∘ leftFace q`
and `leftFace id = id`; together with `pullFace_comp` this is the group
action of the listed order-sixty group on `Fin 20 → ℝ`. -/
def leftFace (p : List Nat) (F : Fin 20 → ℝ) : Fin 20 → ℝ := pullFace (invPerm p) F

theorem pullFace_inv_left (p : List Nat) (hp : p ∈ perms) (F : Fin 20 → ℝ) :
    pullFace (invPerm p) (pullFace p F) = F := by
  rw [← pullFace_comp p hp (invPerm p) (inv_spec p hp).1, (inv_spec p hp).2.1, pullFace_id]

theorem pullFace_inv_right (p : List Nat) (hp : p ∈ perms) (F : Fin 20 → ℝ) :
    pullFace p (pullFace (invPerm p) F) = F := by
  rw [← pullFace_comp (invPerm p) (inv_spec p hp).1 p hp, (inv_spec p hp).2.2, pullFace_id]

theorem leftFace_comp (p : List Nat) (hp : p ∈ perms) (q : List Nat) (hq : q ∈ perms)
    (F : Fin 20 → ℝ) : leftFace (comp p q) F = leftFace p (leftFace q F) := by
  unfold leftFace
  have hr := comp_mem p hp q hq
  conv_rhs => rw [← pullFace_inv_right (comp p q) hr F]
  rw [pullFace_comp p hp q hq, pullFace_inv_left q hq, pullFace_inv_left p hp]

theorem leftFace_id (F : Fin 20 → ℝ) : leftFace (List.range 12) F = F := by
  unfold leftFace
  rw [invPerm_id, pullFace_id]

/-- The induced signed face permutation as an integer matrix. -/
def faceActZ (p : List Nat) : Matrix (Fin 20) (Fin 20) ℤ :=
  Matrix.of fun i j ↦ if facePerm p i = j then faceSign p i else 0

theorem faceActZ_mulVec (p : List Nat) (F : Fin 20 → ℝ) :
    (castZ (faceActZ p)).mulVec F = pullFace p F := by
  funext i
  simp only [Matrix.mulVec, dotProduct, castZ_apply, faceActZ, Matrix.of_apply, pullFace]
  simp [Finset.sum_ite_eq, ite_mul]

/-- **(5)** The real trace of the induced action against a scaled projector
is the integer trace over the normalisation. -/
theorem trace_faceAct_scaledProj (p : List Nat) (Q : Matrix (Fin 20) (Fin 20) ℤ) (d : ℤ) :
    Matrix.trace (castZ (faceActZ p) * scaledProj Q d) = ((d : ℝ))⁻¹ * (traceZ Q p : ℝ) := by
  unfold scaledProj traceZ
  rw [Matrix.mul_smul, Matrix.trace_smul, smul_eq_mul]
  congr 1
  simp only [Matrix.trace, Matrix.diag, Matrix.mul_apply, castZ_apply, faceActZ, Matrix.of_apply]
  push_cast
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  simp [ite_mul, Finset.sum_ite_eq]

/-- **(5)** Character theorem: on each projector image, the trace of the
induced action equals the `A5` character value, as a function of the
element order, for every listed row. -/
theorem projector_characters_R (p : List Nat) (hp : p ∈ perms) :
    Matrix.trace (castZ (faceActZ p) * projZeroR) = chiOne (elemOrder p) ∧
    Matrix.trace (castZ (faceActZ p) * projTwoR) = chiFive (elemOrder p) ∧
    Matrix.trace (castZ (faceActZ p) * projThreeR) = chiFour (elemOrder p) ∧
    Matrix.trace (castZ (faceActZ p) * projFiveR) = chiFour (elemOrder p) ∧
    Matrix.trace (castZ (faceActZ p) * projGoldenR) = chiThreePair (elemOrder p) := by
  obtain ⟨h0, h2, h3, h5, hg⟩ := projector_characters p hp
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · unfold projZeroR; rw [trace_faceAct_scaledProj, h0]; push_cast; ring
  · unfold projTwoR; rw [trace_faceAct_scaledProj, h2]; push_cast; ring
  · unfold projThreeR; rw [trace_faceAct_scaledProj, h3]; push_cast; ring
  · unfold projFiveR; rw [trace_faceAct_scaledProj, h5]; push_cast; ring
  · unfold projGoldenR; rw [trace_faceAct_scaledProj, hg]; push_cast; ring

end OPH.CarrierModeEquivariance

#print axioms OPH.CarrierModeEquivariance.seamPerm_injective
#print axioms OPH.CarrierModeEquivariance.faceSign_eq_one
#print axioms OPH.CarrierModeEquivariance.incidence_equivariant
#print axioms OPH.CarrierModeEquivariance.localMaxwellOperator_pull
#print axioms OPH.CarrierModeEquivariance.faceNormal_pull
#print axioms OPH.CarrierModeEquivariance.projector_images_invariant
#print axioms OPH.CarrierModeEquivariance.cosHistory_pull_ampere
#print axioms OPH.CarrierModeEquivariance.ampereEvolutionScaled_pull
#print axioms OPH.CarrierModeEquivariance.fieldEnergyScaled_pull
#print axioms OPH.CarrierModeEquivariance.pullFace_comp
#print axioms OPH.CarrierModeEquivariance.leftFace_comp
#print axioms OPH.CarrierModeEquivariance.pullFace_const
#print axioms OPH.CarrierModeEquivariance.projector_characters_R
#print axioms OPH.CarrierModeEquivariance.character_norms
