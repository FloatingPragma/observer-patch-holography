import GoldenSectorIrreducibility

open scoped BigOperators Matrix

namespace OPH.GoldenSectorComplexIrreducibility

open OPH.A5PortAction
open OPH.LocalFaceMaxwellAction
open OPH.ScaledMaxwellStability
open OPH.CarrierModeEquivariance
open OPH.ScreenCarrierMapCandidate
open OPH.GoldenSectorCharacters
open OPH.GoldenSectorIrreducibility

/-!
# Complex irreducibility, inequivalence, and real-type commutant of the
two golden pieces

STATUS.  Base change of the committed real Burnside certificate of
`Screen/GoldenSectorIrreducibility.lean` to `ℂ`, plus one further kernel
check (a listed row of order five with the two listed character values),
transported to exact complex linear algebra through `Matrix.map
Complex.ofReal` and to the language of submodules through
`Matrix.toLin'`.  The carrier, its face orientation, the operator
`N = C Cᵀ`, and the two golden projector tables are declared in the
imported files.  No register row is discharged.

WHAT IS PROVED.  For each of the two certificates `plusCert`
(`goldenPlusZ`) and `minusCert` (`goldenMinusZ`), with complexified
objects `projC = projR.map ofReal`, `ARC`, `BRC`, `gC p = (gR p).map
ofReal`, `WC = range (toLin' projC)`, and transported tables
`gtC b = BRC * gC (elems b) * ARC`:
1. Complex Burnside span.  The complex span of the nine transported
   tables of the same nine listed rows is all of `Matrix (Fin 3)
   (Fin 3) ℂ` (`span_gtC_eq_top`): every complex `3 × 3` matrix is the
   explicit combination `∑ b, ((∑ q, Y q · N q b) / e) • gtC b` with the
   committed inverse certificate `N`, `e` of the real file read through
   `ℝ → ℂ` (`eq_sum_gtC`, `NR_mul_M_C`); no complex certificate is
   recomputed.  The complex span `SC` of the sixty matrices `gC p *
   projC` has dimension nine (`finrank_SC`) and maps onto the full
   transported endomorphism algebra (`map_PhiC_SC_eq_top`).
2. Complex irreducibility.  Every `ℂ`-submodule `U ≤ WC` with
   `gC p *ᵥ u ∈ U` for all listed `p` and all `u ∈ U` is `⊥` or `WC`
   (`complex_irreducible`; instances `goldenPlus_complex_irreducible`,
   `goldenMinus_complex_irreducible`).  Burnside argument: the transported
   image of `U` is stable under the nine tables, hence under every complex
   `3 × 3` matrix (`mulVec_mem_of_gtC`), hence trivial
   (`trivial_of_forall_mulVec`, by a rank-one matrix `v ⊗ eᵢ / uᵢ` at a
   nonzero coordinate `uᵢ`, which replaces the real proof's `u ⬝ u ≠ 0`
   since complex vectors can be isotropic).  `finrank ℂ WC = 3`
   (`finrank_WC`).
3. Inequivalence.  The listed row `p5 = [0, 2, 4, 1, 6, 8, 3, 5, 10, 7,
   9, 11]` (the second of the nine certificate rows, `p5_eq_elems`) has
   order five (`p5_order`) and character values `20 (1 - φ)` on
   `goldenPlusZ` and `20 φ` on `goldenMinusZ` (`chi_plus_p5`,
   `chi_minus_p5`); the transported real traces are `1 - φ` and `φ`
   (`trace_plus_p5`,
   `trace_minus_p5`, from `trace_transported`: the trace of the
   transported action of any listed row is `evalPhi (traceZphi Q p) /
   20`).  Hence: no invertible complex `3 × 3` matrix intertwines the
   transported complex actions of the two pieces on all listed rows
   (`no_unit_intertwiner`; an intertwiner conjugates one action to the
   other and preserves the trace, `trace_ne_of_order_five` with
   `chi_differ_iff`); every intertwiner between them is zero
   (`intertwiner_eq_zero`, Schur form: the kernel of an intertwiner is
   invariant under the `λ₊` action, hence `⊥` or `⊤` by item 2, and `⊥`
   forces invertibility); and every complex `20 × 20` matrix `T` with
   `projC₋ * T * projC₊ = T` commuting with every `gC p` is zero
   (`intertwiner_eq_zero_ambient`).  Instance:
   `complex_pieces_inequivalent`.
4. Commutant.  A real `3 × 3` matrix commuting with the nine transported
   tables of a piece is a real scalar (`commutant_scalar_real`,
   `commutant_scalar_real_perms`), and a real `20 × 20` matrix `T` with
   `T * P = T = P * T` commuting with every listed `gR p` is a real
   multiple of `P` (`commutant_on_piece_real`): the endomorphism algebra
   of each golden piece over `ℝ` is `ℝ · 1`, the real-type case of the
   Frobenius–Schur trichotomy.  The complex analogue is
   `commutant_scalar_complex`.  Route: a matrix commuting with a spanning
   set commutes with every matrix, in particular with every matrix unit,
   and `Matrix.mem_range_scalar_of_commute_single` gives the scalar.

NOT PROVED HERE (inferences, stated for scope).  The identification of
the listed group with the abstract alternating group `A5` and of the two
complexified pieces with the two Galois-conjugate three-dimensional
irreducible characters of the abstract character table remain outside
this file; the theorems concern the sixty listed rows and the two
declared projector tables.  No physical interpretation of a golden mode
is asserted.  Absolute irreducibility in the sense of every field
extension is not stated in Lean; the statement here is for `ℂ`.  `ℂ` is
algebraically closed, so irreducibility over `ℂ` is absolute
irreducibility of the real representation; the statement for an
arbitrary field extension is simply not formalised.

PRIOR WORK.  `Screen/GoldenSectorIrreducibility.lean` proves the real
statements (`finrank_S`, `span_gt_eq_top`, `irreducible`,
`goldenPlus_irreducible`, `goldenMinus_irreducible`) and states complex
irreducibility and the character identification as inferences; this file
replaces the first inference by theorems.  `Screen/GoldenSectorCharacters.lean`
supplies the character values (`chi_plus_values`, `chi_minus_values`,
`chi_differ_iff`) and `traceZphi_eq`.  `QFT/GaugeIrreducibleBorn.lean`
proves a `2 × 2` complex commutant statement (`commutant_scalar`) for a
different pair of generators by direct entry computation.  The content
of this file beyond the prior modules is: the complex base change of the
span certificate, the complex irreducibility theorem for both pieces, the
zero-intertwiner theorem between the two complexified pieces, and the
real-scalar commutant.

ROWS TOUCHED (none discharged).  Premise row PR-53 (physical
finite-carrier propagation, frame, and comparison attachment): cited at
scope, untouched; the theorems classify the declared incidence action on
the two golden eigenspaces over `ℂ` and identify their endomorphism
algebras, and attach no propagation, frame, clock or length scale, or
readout.  Premise row PR-52 (observer-to-physical-spacetime and causal
attachment), the physical spacetime attachment row: cited at scope,
untouched; the identification of the listed group with a physical
rotation group is open.  Both rows are the open rows of the parent claim
OPH-EM-GOLDEN-SECTOR-IRREDUCIBILITY and stay open here.  Source clock
and duration row: untouched (no step appears).  Light-signal row:
untouched; the identification of a golden mode with a physical
oscillation is open.  Coupled-action row, laboratory clock and energy
calibration import, gravitation-route energy identification: untouched.

NEGATIVES CITED.  Legendre non-identifiability
(`Lean/Variational/RealizedHistoryLegendreNoGo.lean`), at scope only:
the operator `N` whose eigenspaces are classified here is part of a
declared evolution.

CONVENTIONS.  As in the real file: `Zphi = ℤ × ℤ` with `(a, b) = a +
bφ`, `φ² = φ + 1`; `gR p = castZ (faceActZ p)`; matrix positions
`(c, a)` listed at `3c + a`; `gt C b = BR C * gR (elems b) * AR C`.
Complexification is entrywise `Complex.ofReal` (`Matrix.map`); every
complex identity below is the image of a real identity under the ring
homomorphism `Complex.ofRealHom` (`map_mul_ofReal`, `map_one_ofReal`),
so the complex objects carry no new table.  The transported complex
action of a listed row `p` is `PhiC C (gPC C p) = BRC C * gC p * ARC C`
(`PhiC_gPC`), which for the nine certificate rows is `gtC C b`
(`gtC_eq_PhiC`).  The face action inherits the right-action composition
law of `Screen/CarrierModeEquivariance.lean` (`pullFace (comp p q) =
pullFace q ∘ pullFace p`, so `gR (comp p q) = gR q * gR p`); the
families `p ↦ PhiC C (gPC C p)` are therefore anti-representations of
`comp`.  No theorem in this file depends on the composition law: every
statement is about the set of sixty listed matrices.

FALSIFIER.  A row `p5` off the listed membership, off order five, or with
character values off `(20, -20)` and `(0, 20)` fails the kernel checks
`p5_mem`, `p5_order`, `chi_plus_p5`, `chi_minus_p5`; the remaining
statements are transports of the committed certificates and fail with
them.

Axiom audit.  The `#print axioms` lines at the end of the file show at
most `propext`, `Classical.choice`, and `Quot.sound`.  No `native_decide`.
-/

noncomputable section

/-! ## 1. Complex transport of the committed real certificate -/

/-- Entrywise complexification commutes with matrix products. -/
theorem map_mul_ofReal {m n o : Type*} [Fintype n] (L : Matrix m n ℝ) (M : Matrix n o ℝ) :
    (L * M).map Complex.ofReal = L.map Complex.ofReal * M.map Complex.ofReal :=
  Matrix.map_mul (f := Complex.ofRealHom)

theorem map_one_ofReal {n : Type*} [DecidableEq n] :
    (1 : Matrix n n ℝ).map Complex.ofReal = 1 :=
  Matrix.map_one Complex.ofReal Complex.ofReal_zero Complex.ofReal_one

/-- The complexified projector of a certificate. -/
def projC (C : GoldenCert) : Matrix (Fin 20) (Fin 20) ℂ := (projR C).map Complex.ofReal
/-- The complexified transport of the first three projector columns. -/
def ARC (C : GoldenCert) : Matrix (Fin 20) (Fin 3) ℂ := (AR C).map Complex.ofReal
/-- The complexified left inverse table, at scale `1/20`. -/
def BRC (C : GoldenCert) : Matrix (Fin 3) (Fin 20) ℂ := (BR C).map Complex.ofReal
/-- The complexified action matrix of a listed row. -/
def gC (p : List Nat) : Matrix (Fin 20) (Fin 20) ℂ := (gR p).map Complex.ofReal

theorem projC_plus : projC plusCert = goldenPlusR.map Complex.ofReal := rfl
theorem projC_minus : projC minusCert = goldenMinusR.map Complex.ofReal := rfl

theorem BRC_mul_ARC (C : GoldenCert) : BRC C * ARC C = 1 := by
  rw [BRC, ARC, ← map_mul_ofReal, BR_mul_AR, map_one_ofReal]

theorem ARC_mul_BRC (C : GoldenCert) : ARC C * BRC C = projC C := by
  rw [ARC, BRC, ← map_mul_ofReal, AR_mul_BR]; rfl

theorem projC_idem (C : GoldenCert) : projC C * projC C = projC C := by
  rw [← ARC_mul_BRC, Matrix.mul_assoc, ← Matrix.mul_assoc (BRC C), BRC_mul_ARC, Matrix.one_mul]

theorem projC_mul_ARC (C : GoldenCert) : projC C * ARC C = ARC C := by
  rw [← ARC_mul_BRC, Matrix.mul_assoc, BRC_mul_ARC, Matrix.mul_one]

theorem BRC_mul_projC (C : GoldenCert) : BRC C * projC C = BRC C := by
  rw [← ARC_mul_BRC, ← Matrix.mul_assoc, BRC_mul_ARC, Matrix.one_mul]

/-- Complex commutation of every listed row with the projector. -/
theorem gC_comm_projC (C : GoldenCert) (p : List Nat) (hp : p ∈ perms) :
    gC p * projC C = projC C * gC p := by
  rw [gC, projC, ← map_mul_ofReal, ← map_mul_ofReal, gR_comm_projR C p hp]

/-- The complexified transported action table of a certificate row. -/
def gtC (C : GoldenCert) (b : Fin 3 × Fin 3) : Matrix (Fin 3) (Fin 3) ℂ :=
  BRC C * gC (elems b) * ARC C

theorem gtC_eq_map (C : GoldenCert) (b : Fin 3 × Fin 3) :
    gtC C b = (gt C b).map Complex.ofReal := by
  rw [gtC, gt, map_mul_ofReal, map_mul_ofReal]; rfl

theorem gtC_eq (C : GoldenCert) (b : Fin 3 × Fin 3) :
    gtC C b = Matrix.of fun c a ↦ ((evalPhi (C.M b (c, a)) : ℝ) : ℂ) := by
  rw [gtC_eq_map, gt_eq]
  ext c a
  simp [Matrix.map_apply]

/-! ## 2. The complex golden piece -/

/-- The complexified golden piece: the range of the complexified projector. -/
def WC (C : GoldenCert) : Submodule ℂ (Fin 20 → ℂ) := LinearMap.range (Matrix.toLin' (projC C))

theorem mem_WC_iff (C : GoldenCert) (w : Fin 20 → ℂ) : w ∈ WC C ↔ (projC C).mulVec w = w := by
  constructor
  · rintro ⟨v, rfl⟩
    rw [Matrix.toLin'_apply, Matrix.mulVec_mulVec, projC_idem]
  · intro h
    exact ⟨w, by rw [Matrix.toLin'_apply, h]⟩

theorem projC_mulVec_of_mem_WC (C : GoldenCert) {w : Fin 20 → ℂ} (hw : w ∈ WC C) :
    (projC C).mulVec w = w := (mem_WC_iff C w).1 hw

theorem WC_eq_range_ARC (C : GoldenCert) : WC C = LinearMap.range (Matrix.toLin' (ARC C)) := by
  apply le_antisymm
  · rintro _ ⟨v, rfl⟩
    refine ⟨(BRC C).mulVec v, ?_⟩
    rw [Matrix.toLin'_apply, Matrix.toLin'_apply, Matrix.mulVec_mulVec, ARC_mul_BRC]
  · rintro _ ⟨v, rfl⟩
    refine ⟨(ARC C).mulVec v, ?_⟩
    rw [Matrix.toLin'_apply, Matrix.toLin'_apply, Matrix.mulVec_mulVec, projC_mul_ARC]

theorem toLin'_ARC_injective (C : GoldenCert) : Function.Injective (Matrix.toLin' (ARC C)) := by
  intro v v' h
  have h' := congrArg (fun x ↦ (BRC C).mulVec x) h
  simp only [Matrix.toLin'_apply, Matrix.mulVec_mulVec, BRC_mul_ARC, Matrix.one_mulVec] at h'
  exact h'

/-- `finrank ℂ WC = 3`. -/
theorem finrank_WC (C : GoldenCert) : Module.finrank ℂ (WC C) = 3 := by
  rw [WC_eq_range_ARC, LinearMap.finrank_range_of_inj (toLin'_ARC_injective C),
    Module.finrank_fin_fun]

/-- `WC` is invariant under every listed row. -/
theorem WC_invariant (C : GoldenCert) (p : List Nat) (hp : p ∈ perms) {w : Fin 20 → ℂ}
    (hw : w ∈ WC C) : (gC p).mulVec w ∈ WC C := by
  rw [mem_WC_iff, Matrix.mulVec_mulVec, ← gC_comm_projC C p hp, ← Matrix.mulVec_mulVec,
    projC_mulVec_of_mem_WC C hw]

/-! ## 3. Base change of the span certificate -/

theorem eC_ne (C : GoldenCert) : ((evalPhi C.e : ℝ) : ℂ) ≠ 0 :=
  Complex.ofReal_ne_zero.mpr (eR_ne C)

/-- Complex form of the committed certificate `N · M = e · 1`: the same
coefficients, read through `ℝ → ℂ`. -/
theorem NR_mul_M_C (C : GoldenCert) (q q' : Fin 3 × Fin 3) :
    (∑ b : Fin 3 × Fin 3, ((NR C q b : ℝ) : ℂ) * ((evalPhi (C.M b q') : ℝ) : ℂ)) =
      if q = q' then ((evalPhi C.e : ℝ) : ℂ) else 0 := by
  have h := NR_mul_M C q q'
  by_cases hq : q = q'
  · rw [if_pos hq] at h ⊢
    have h' := congrArg (fun x : ℝ ↦ (x : ℂ)) h
    push_cast at h'
    exact h'
  · rw [if_neg hq] at h ⊢
    have h' := congrArg (fun x : ℝ ↦ (x : ℂ)) h
    push_cast at h'
    exact h'

/-- **(1)** Every complex `3 × 3` matrix is the explicit combination of
the nine complexified transported tables with the committed coefficients. -/
theorem eq_sum_gtC (C : GoldenCert) (Y : Matrix (Fin 3) (Fin 3) ℂ) :
    Y = ∑ b : Fin 3 × Fin 3,
      ((∑ q : Fin 3 × Fin 3, Y q.1 q.2 * ((NR C q b : ℝ) : ℂ)) / ((evalPhi C.e : ℝ) : ℂ)) •
        gtC C b := by
  ext c a
  simp only [Matrix.sum_apply, Matrix.smul_apply, gtC_eq, Matrix.of_apply, smul_eq_mul]
  have h : (∑ b : Fin 3 × Fin 3,
      (∑ q : Fin 3 × Fin 3, Y q.1 q.2 * ((NR C q b : ℝ) : ℂ)) / ((evalPhi C.e : ℝ) : ℂ) *
        ((evalPhi (C.M b (c, a)) : ℝ) : ℂ)) =
      (∑ q : Fin 3 × Fin 3, Y q.1 q.2 *
        ∑ b : Fin 3 × Fin 3, ((NR C q b : ℝ) : ℂ) * ((evalPhi (C.M b (c, a)) : ℝ) : ℂ)) /
        ((evalPhi C.e : ℝ) : ℂ) := by
    simp only [Finset.sum_div, div_mul_eq_mul_div, Finset.sum_mul, Finset.mul_sum]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun q _ ↦ Finset.sum_congr rfl fun b _ ↦ ?_
    ring
  rw [h]
  simp only [NR_mul_M_C, mul_ite, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, if_true]
  rw [mul_div_assoc, div_self (eC_ne C), mul_one]

/-- **(1)** The complex span of the nine complexified transported tables is
all of `Matrix (Fin 3) (Fin 3) ℂ`. -/
theorem span_gtC_eq_top (C : GoldenCert) : Submodule.span ℂ (Set.range (gtC C)) = ⊤ := by
  rw [eq_top_iff]
  intro Y _
  rw [eq_sum_gtC C Y]
  exact Submodule.sum_mem _ fun b _ ↦ Submodule.smul_mem _ _ (Submodule.subset_span ⟨b, rfl⟩)

/-! ## 4. The complex Burnside span at the carrier level -/

/-- The generator `gC p * projC` of the complex Burnside span. -/
def gPC (C : GoldenCert) (p : List Nat) : Matrix (Fin 20) (Fin 20) ℂ := gC p * projC C

/-- The complex Burnside span of the sixty matrices `gC p * projC`. -/
def SC (C : GoldenCert) : Submodule ℂ (Matrix (Fin 20) (Fin 20) ℂ) :=
  Submodule.span ℂ {X | ∃ p ∈ perms, X = gPC C p}

/-- The transport of a complex `20 × 20` matrix to the basis `ARC` of `WC`. -/
def PhiC (C : GoldenCert) : Matrix (Fin 20) (Fin 20) ℂ →ₗ[ℂ] Matrix (Fin 3) (Fin 3) ℂ where
  toFun X := BRC C * X * ARC C
  map_add' X Y := by simp only [Matrix.mul_add, Matrix.add_mul]
  map_smul' c X := by simp only [Matrix.mul_smul, Matrix.smul_mul, RingHom.id_apply]

theorem PhiC_apply (C : GoldenCert) (X : Matrix (Fin 20) (Fin 20) ℂ) :
    PhiC C X = BRC C * X * ARC C := rfl

theorem sandwich_of_mem_SC (C : GoldenCert) {X : Matrix (Fin 20) (Fin 20) ℂ} (hX : X ∈ SC C) :
    projC C * X * projC C = X := by
  induction hX using Submodule.span_induction with
  | mem X hX =>
    obtain ⟨p, hp, rfl⟩ := hX
    unfold gPC
    rw [← Matrix.mul_assoc, ← gC_comm_projC C p hp, Matrix.mul_assoc, projC_idem,
      Matrix.mul_assoc, projC_idem]
  | zero => simp
  | add X Y _ _ hX hY => rw [Matrix.mul_add, Matrix.add_mul, hX, hY]
  | smul c X _ hX => rw [Matrix.mul_smul, Matrix.smul_mul, hX]

theorem PhiC_injOn_SC (C : GoldenCert) {X : Matrix (Fin 20) (Fin 20) ℂ} (hX : X ∈ SC C)
    (h : PhiC C X = 0) : X = 0 := by
  have h2 : ARC C * PhiC C X * BRC C = 0 := by rw [h]; simp
  rw [PhiC_apply, ← Matrix.mul_assoc, ← Matrix.mul_assoc, ARC_mul_BRC, Matrix.mul_assoc,
    Matrix.mul_assoc, ARC_mul_BRC, ← Matrix.mul_assoc, sandwich_of_mem_SC C hX] at h2
  exact h2

theorem finrank_matrix_three_C : Module.finrank ℂ (Matrix (Fin 3) (Fin 3) ℂ) = 9 := by
  rw [Module.finrank_matrix]
  simp

theorem finrank_SC_le (C : GoldenCert) : Module.finrank ℂ (SC C) ≤ 9 := by
  have hinj : Function.Injective ((PhiC C).domRestrict (SC C)) := by
    rw [← LinearMap.ker_eq_bot, LinearMap.ker_eq_bot']
    intro x hx
    apply Subtype.ext
    exact PhiC_injOn_SC C x.2 hx
  rw [← finrank_matrix_three_C]
  exact LinearMap.finrank_le_finrank_of_injective hinj

theorem PhiC_gPC (C : GoldenCert) (p : List Nat) (_hp : p ∈ perms) :
    PhiC C (gPC C p) = BRC C * gC p * ARC C := by
  rw [PhiC_apply, gPC, Matrix.mul_assoc, Matrix.mul_assoc, projC_mul_ARC, ← Matrix.mul_assoc]

theorem gtC_eq_PhiC (C : GoldenCert) (b : Fin 3 × Fin 3) :
    gtC C b = PhiC C (gPC C (elems b)) := by
  rw [PhiC_gPC C _ (elems_mem b)]; rfl

theorem gPC_mem_SC (C : GoldenCert) (p : List Nat) (hp : p ∈ perms) : gPC C p ∈ SC C :=
  Submodule.subset_span ⟨p, hp, rfl⟩

/-- **(1)** `PhiC` maps the complex Burnside span onto all complex `3 × 3`
matrices. -/
theorem map_PhiC_SC_eq_top (C : GoldenCert) : (SC C).map (PhiC C) = ⊤ := by
  rw [eq_top_iff, ← span_gtC_eq_top, Submodule.span_le]
  rintro _ ⟨b, rfl⟩
  exact ⟨gPC C (elems b), gPC_mem_SC C _ (elems_mem b), (gtC_eq_PhiC C b).symm⟩

theorem finrank_SC_ge (C : GoldenCert) : 9 ≤ Module.finrank ℂ (SC C) := by
  have h := Submodule.finrank_map_le (PhiC C) (SC C)
  rw [map_PhiC_SC_eq_top, finrank_top, finrank_matrix_three_C] at h
  exact h

/-- **(1)** `finrank ℂ SC = 9`: the complex Burnside span of a golden piece
(a subspace of the complex `20 × 20` matrices) maps isomorphically under
`PhiC` onto the full complex endomorphism algebra of the three-dimensional
piece. -/
theorem finrank_SC (C : GoldenCert) : Module.finrank ℂ (SC C) = 9 :=
  le_antisymm (finrank_SC_le C) (finrank_SC_ge C)

/-! ## 5. Complex irreducibility -/

/-- A complex submodule of `Fin 3 → ℂ` stable under the nine complexified
tables is stable under every complex `3 × 3` matrix. -/
theorem mulVec_mem_of_gtC (C : GoldenCert) (U' : Submodule ℂ (Fin 3 → ℂ))
    (h : ∀ b, ∀ u ∈ U', (gtC C b).mulVec u ∈ U') (Y : Matrix (Fin 3) (Fin 3) ℂ) {u : Fin 3 → ℂ}
    (hu : u ∈ U') : Y.mulVec u ∈ U' := by
  have hY : Y ∈ Submodule.span ℂ (Set.range (gtC C)) := by rw [span_gtC_eq_top]; trivial
  induction hY using Submodule.span_induction with
  | mem Z hZ =>
    obtain ⟨b, rfl⟩ := hZ
    exact h b u hu
  | zero => rw [Matrix.zero_mulVec]; exact U'.zero_mem
  | add Z Z' _ _ hZ hZ' => rw [Matrix.add_mulVec]; exact U'.add_mem hZ hZ'
  | smul c Z _ hZ => rw [Matrix.smul_mulVec]; exact U'.smul_mem c hZ

/-- A complex submodule of `Fin 3 → ℂ` stable under every matrix is trivial.
The rank-one matrix `v ⊗ eᵢ / uᵢ` at a nonzero coordinate `uᵢ` carries
`u` to `v`. -/
theorem trivial_of_forall_mulVec (U' : Submodule ℂ (Fin 3 → ℂ))
    (h : ∀ Y : Matrix (Fin 3) (Fin 3) ℂ, ∀ u ∈ U', Y.mulVec u ∈ U') : U' = ⊥ ∨ U' = ⊤ := by
  by_cases hbot : U' = ⊥
  · exact Or.inl hbot
  right
  obtain ⟨u, hu, hne⟩ := Submodule.exists_mem_ne_zero_of_ne_bot hbot
  obtain ⟨i, hi⟩ : ∃ i, u i ≠ 0 := Function.ne_iff.mp hne
  rw [eq_top_iff]
  intro v _
  let Y : Matrix (Fin 3) (Fin 3) ℂ := Matrix.of fun j k ↦ if k = i then v j / u i else 0
  have hYu : Y.mulVec u = v := by
    funext j
    simp only [Matrix.mulVec, dotProduct, Y, Matrix.of_apply, ite_mul, zero_mul,
      Finset.sum_ite_eq', Finset.mem_univ, if_true]
    exact div_mul_cancel₀ (v j) hi
  rw [← hYu]
  exact h Y u hu

/-- Transported form: a complex submodule of `Fin 3 → ℂ` stable under the
nine complexified tables is trivial. -/
theorem transported_irreducible (C : GoldenCert) (U' : Submodule ℂ (Fin 3 → ℂ))
    (h : ∀ b, ∀ u ∈ U', (gtC C b).mulVec u ∈ U') : U' = ⊥ ∨ U' = ⊤ :=
  trivial_of_forall_mulVec U' fun Y _ hu ↦ mulVec_mem_of_gtC C U' h Y hu

/-- **(2)** Complex irreducibility: a complex submodule of the complexified
golden piece invariant under every listed row is `⊥` or the whole piece. -/
theorem complex_irreducible (C : GoldenCert) (U : Submodule ℂ (Fin 20 → ℂ)) (hU : U ≤ WC C)
    (hinv : ∀ p ∈ perms, ∀ u ∈ U, (gC p).mulVec u ∈ U) : U = ⊥ ∨ U = WC C := by
  have hstab : ∀ b, ∀ u' ∈ U.map (Matrix.toLin' (BRC C)),
      (gtC C b).mulVec u' ∈ U.map (Matrix.toLin' (BRC C)) := by
    rintro b _ ⟨u, hu, rfl⟩
    refine ⟨(gC (elems b)).mulVec u, hinv _ (elems_mem b) u hu, ?_⟩
    simp only [Matrix.toLin'_apply, gtC, Matrix.mulVec_mulVec]
    rw [Matrix.mul_assoc, ARC_mul_BRC, ← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec,
      projC_mulVec_of_mem_WC C (hU hu), Matrix.mulVec_mulVec]
  rcases transported_irreducible C _ hstab with hb | ht
  · left
    rw [eq_bot_iff]
    intro u hu
    have hBu : (BRC C).mulVec u = 0 := by
      have : Matrix.toLin' (BRC C) u ∈ U.map (Matrix.toLin' (BRC C)) := ⟨u, hu, rfl⟩
      rw [hb, Submodule.mem_bot] at this
      exact this
    have : u = 0 := by
      rw [← projC_mulVec_of_mem_WC C (hU hu), ← ARC_mul_BRC, ← Matrix.mulVec_mulVec, hBu,
        Matrix.mulVec_zero]
    rw [this]; exact Submodule.zero_mem _
  · right
    refine le_antisymm hU ?_
    intro w hw
    have : Matrix.toLin' (BRC C) w ∈ U.map (Matrix.toLin' (BRC C)) := by rw [ht]; trivial
    obtain ⟨u, hu, huw⟩ := this
    have hw' : w = u := by
      rw [← projC_mulVec_of_mem_WC C hw, ← projC_mulVec_of_mem_WC C (hU hu), ← ARC_mul_BRC,
        ← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec]
      simp only [Matrix.toLin'_apply] at huw
      rw [huw]
    rw [hw']; exact hu

/-! ## 6. The two complexified golden pieces -/

/-- The complexified `λ₊` piece. -/
def WCPlus : Submodule ℂ (Fin 20 → ℂ) :=
  LinearMap.range (Matrix.toLin' (goldenPlusR.map Complex.ofReal))
/-- The complexified `λ₋` piece. -/
def WCMinus : Submodule ℂ (Fin 20 → ℂ) :=
  LinearMap.range (Matrix.toLin' (goldenMinusR.map Complex.ofReal))

theorem WCPlus_eq : WCPlus = WC plusCert := rfl
theorem WCMinus_eq : WCMinus = WC minusCert := rfl

theorem finrank_WCPlus : Module.finrank ℂ WCPlus = 3 := finrank_WC plusCert
theorem finrank_WCMinus : Module.finrank ℂ WCMinus = 3 := finrank_WC minusCert

/-- Non-vacuity of the invariance hypothesis: the whole piece satisfies it. -/
theorem WCPlus_invariant (p : List Nat) (hp : p ∈ perms) {w : Fin 20 → ℂ} (hw : w ∈ WCPlus) :
    ((castZ (faceActZ p)).map Complex.ofReal).mulVec w ∈ WCPlus := WC_invariant plusCert p hp hw
theorem WCMinus_invariant (p : List Nat) (hp : p ∈ perms) {w : Fin 20 → ℂ} (hw : w ∈ WCMinus) :
    ((castZ (faceActZ p)).map Complex.ofReal).mulVec w ∈ WCMinus := WC_invariant minusCert p hp hw

/-- **(1)** The complex Burnside span of the `λ₊` piece has dimension nine. -/
theorem finrank_SCPlus :
    Module.finrank ℂ (Submodule.span ℂ
      {X : Matrix (Fin 20) (Fin 20) ℂ | ∃ p ∈ perms,
        X = (castZ (faceActZ p)).map Complex.ofReal * goldenPlusR.map Complex.ofReal}) = 9 :=
  finrank_SC plusCert
/-- **(1)** The complex Burnside span of the `λ₋` piece has dimension nine. -/
theorem finrank_SCMinus :
    Module.finrank ℂ (Submodule.span ℂ
      {X : Matrix (Fin 20) (Fin 20) ℂ | ∃ p ∈ perms,
        X = (castZ (faceActZ p)).map Complex.ofReal * goldenMinusR.map Complex.ofReal}) = 9 :=
  finrank_SC minusCert

/-- **(2)** The complexified `λ₊` golden piece is irreducible under the
listed group. -/
theorem goldenPlus_complex_irreducible (U : Submodule ℂ (Fin 20 → ℂ)) (hU : U ≤ WCPlus)
    (hinv : ∀ p ∈ perms, ∀ u ∈ U, ((castZ (faceActZ p)).map Complex.ofReal).mulVec u ∈ U) :
    U = ⊥ ∨ U = WCPlus :=
  complex_irreducible plusCert U hU hinv

/-- **(2)** The complexified `λ₋` golden piece is irreducible under the
listed group. -/
theorem goldenMinus_complex_irreducible (U : Submodule ℂ (Fin 20 → ℂ)) (hU : U ≤ WCMinus)
    (hinv : ∀ p ∈ perms, ∀ u ∈ U, ((castZ (faceActZ p)).map Complex.ofReal).mulVec u ∈ U) :
    U = ⊥ ∨ U = WCMinus :=
  complex_irreducible minusCert U hU hinv

/-! ## 7. Traces of the transported actions and an order-five row -/

/-- The trace of the transported real action of a listed row is the
character value `traceZphi Q p` read through `evalPhi`, at scale `1/20`. -/
theorem trace_transported (C : GoldenCert) (p : List Nat) (hp : p ∈ perms) :
    Matrix.trace (Phi C (gP C p)) = evalPhi (traceZphi C.Q p) / 20 := by
  rw [Phi_gP C p hp, Matrix.trace_mul_comm, ← Matrix.mul_assoc, AR_mul_BR, Matrix.trace_mul_comm]
  simp only [Matrix.trace, Matrix.diag_apply, Matrix.mul_apply, gR, projR, castZ_apply, faceActZ,
    Matrix.of_apply, Int.cast_ite, Int.cast_zero, ite_mul, zero_mul, Finset.sum_ite_eq,
    Finset.mem_univ, if_true, faceSign_eq_one p hp, Int.cast_one, one_mul]
  rw [traceZphi_eq C.Q p hp, evalPhi_sum, Finset.sum_div]

/-- The transported traces of the two pieces differ on every listed row of
order five. -/
theorem trace_ne_of_order_five (p : List Nat) (hp : p ∈ perms) (h5 : elemOrder p = 5) :
    Matrix.trace (Phi plusCert (gP plusCert p)) ≠ Matrix.trace (Phi minusCert (gP minusCert p)) := by
  rw [trace_transported _ _ hp, trace_transported _ _ hp]
  intro h
  have h' : evalPhi (traceZphi goldenPlusZ p) = evalPhi (traceZphi goldenMinusZ p) :=
    (div_left_inj' (by norm_num : (20 : ℝ) ≠ 0)).1 h
  exact (chi_differ_iff p hp).1 (evalPhi_injective h') h5

/-- A listed row of order five: the second of the nine certificate rows. -/
def p5 : List Nat := [0, 2, 4, 1, 6, 8, 3, 5, 10, 7, 9, 11]

theorem p5_mem : p5 ∈ perms := by decide
theorem p5_eq_elems : p5 = elems (0, 1) := by decide
theorem p5_order : elemOrder p5 = 5 := by decide +kernel
theorem chi_plus_p5 : traceZphi goldenPlusZ p5 = (20, -20) := by decide +kernel
theorem chi_minus_p5 : traceZphi goldenMinusZ p5 = (0, 20) := by decide +kernel

/-- Numeric check: the transported `λ₊` action of `p5` has trace `1 - φ`. -/
theorem trace_plus_p5 : Matrix.trace (Phi plusCert (gP plusCert p5)) = 1 - Real.goldenRatio := by
  rw [trace_transported _ _ p5_mem]
  show evalPhi (traceZphi goldenPlusZ p5) / 20 = _
  rw [chi_plus_p5]
  simp only [evalPhi]
  push_cast
  ring

/-- Numeric check: the transported `λ₋` action of `p5` has trace `φ`. -/
theorem trace_minus_p5 :
    Matrix.trace (Phi minusCert (gP minusCert p5)) = Real.goldenRatio := by
  rw [trace_transported _ _ p5_mem]
  show evalPhi (traceZphi goldenMinusZ p5) / 20 = _
  rw [chi_minus_p5]
  simp only [evalPhi]
  push_cast
  ring

/-! ## 8. Inequivalence of the two complexified pieces -/

theorem PhiC_gPC_eq_map (C : GoldenCert) (p : List Nat) (hp : p ∈ perms) :
    PhiC C (gPC C p) = (Phi C (gP C p)).map Complex.ofReal := by
  rw [PhiC_gPC C p hp, Phi_gP C p hp, map_mul_ofReal, map_mul_ofReal]; rfl

theorem trace_map_ofReal (X : Matrix (Fin 3) (Fin 3) ℝ) :
    Matrix.trace (X.map Complex.ofReal) = ((Matrix.trace X : ℝ) : ℂ) := by
  simp only [Matrix.trace, Matrix.diag_apply, Matrix.map_apply, Complex.ofReal_sum]

/-- **(3)** No invertible complex `3 × 3` matrix intertwines the nine
complexified transported tables of the two pieces: an intertwiner would
conjugate one table to the other and preserve the trace, but the traces
on `p5` are `1 - φ` and `φ`. -/
theorem no_unit_intertwiner_nine :
    ¬ ∃ T : Matrix (Fin 3) (Fin 3) ℂ, IsUnit T ∧
      ∀ b, T * gtC plusCert b = gtC minusCert b * T := by
  rintro ⟨T, hT, h⟩
  obtain ⟨u, rfl⟩ := hT
  have hb := h (0, 1)
  have e1 : gtC plusCert (0, 1) = (Phi plusCert (gP plusCert p5)).map Complex.ofReal := by
    rw [gtC_eq_map, gt, Phi_gP _ _ p5_mem, p5_eq_elems]
  have e2 : gtC minusCert (0, 1) = (Phi minusCert (gP minusCert p5)).map Complex.ofReal := by
    rw [gtC_eq_map, gt, Phi_gP _ _ p5_mem, p5_eq_elems]
  have htr : Matrix.trace (gtC plusCert (0, 1)) = Matrix.trace (gtC minusCert (0, 1)) := by
    have h1 : gtC plusCert (0, 1) = (↑u⁻¹ : Matrix (Fin 3) (Fin 3) ℂ) *
        (gtC minusCert (0, 1) * (↑u : Matrix (Fin 3) (Fin 3) ℂ)) := by
      rw [← hb, ← Matrix.mul_assoc, Units.inv_mul, Matrix.one_mul]
    rw [h1, Matrix.trace_mul_comm, Matrix.mul_assoc, Units.mul_inv, Matrix.mul_one]
  rw [e1, e2, trace_map_ofReal, trace_map_ofReal] at htr
  exact trace_ne_of_order_five p5 p5_mem p5_order (Complex.ofReal_injective htr)

/-- **(3)** No invertible complex `3 × 3` matrix intertwines the
complexified transported actions of the two pieces on all listed rows. -/
theorem no_unit_intertwiner :
    ¬ ∃ T : Matrix (Fin 3) (Fin 3) ℂ, IsUnit T ∧
      ∀ p ∈ perms, T * PhiC plusCert (gPC plusCert p) = PhiC minusCert (gPC minusCert p) * T := by
  rintro ⟨T, hT, h⟩
  refine no_unit_intertwiner_nine ⟨T, hT, fun b ↦ ?_⟩
  rw [gtC_eq_PhiC, gtC_eq_PhiC]
  exact h _ (elems_mem b)

/-- **(3)** Schur form: every complex `3 × 3` matrix intertwining the nine
complexified transported tables of the two pieces is zero. -/
theorem intertwiner_eq_zero_nine (T : Matrix (Fin 3) (Fin 3) ℂ)
    (h : ∀ b, T * gtC plusCert b = gtC minusCert b * T) : T = 0 := by
  by_contra hT0
  have hker : ∀ b, ∀ u ∈ LinearMap.ker (Matrix.toLin' T),
      (gtC plusCert b).mulVec u ∈ LinearMap.ker (Matrix.toLin' T) := by
    intro b u hu
    rw [LinearMap.mem_ker, Matrix.toLin'_apply] at hu ⊢
    rw [Matrix.mulVec_mulVec, h b, ← Matrix.mulVec_mulVec, hu, Matrix.mulVec_zero]
  rcases transported_irreducible plusCert _ hker with hb | ht
  · have hinj : Function.Injective T.mulVec := by
      have hk := (LinearMap.ker_eq_bot (f := Matrix.toLin' T)).1 hb
      intro x y hxy
      apply hk
      simpa [Matrix.toLin'_apply] using hxy
    exact no_unit_intertwiner_nine ⟨T, Matrix.mulVec_injective_iff_isUnit.1 hinj, h⟩
  · exact hT0 (Matrix.toLin'.map_eq_zero_iff.1 ((LinearMap.ker_eq_top (f := Matrix.toLin' T)).1 ht))

/-- **(3)** Schur form on all listed rows: every complex `3 × 3` matrix
intertwining the complexified transported actions of the two pieces is
zero. -/
theorem intertwiner_eq_zero (T : Matrix (Fin 3) (Fin 3) ℂ)
    (h : ∀ p ∈ perms, T * PhiC plusCert (gPC plusCert p) = PhiC minusCert (gPC minusCert p) * T) :
    T = 0 :=
  intertwiner_eq_zero_nine T fun b ↦ by
    rw [gtC_eq_PhiC, gtC_eq_PhiC]
    exact h _ (elems_mem b)

/-- **(3)** Ambient Schur form: a complex `20 × 20` matrix carrying the
`λ₊` piece into the `λ₋` piece and vanishing off it (`projC₋ * T * projC₊ =
T`) that commutes with every listed row is zero. -/
theorem intertwiner_eq_zero_ambient (T : Matrix (Fin 20) (Fin 20) ℂ)
    (hT : projC minusCert * T * projC plusCert = T)
    (h : ∀ p ∈ perms, T * gC p = gC p * T) : T = 0 := by
  have hTP : T * projC plusCert = T := by
    conv_lhs => rw [← hT]
    rw [Matrix.mul_assoc, projC_idem]
    exact hT
  have hPT : projC minusCert * T = T := by
    conv_lhs => rw [← hT]
    rw [← Matrix.mul_assoc, ← Matrix.mul_assoc, projC_idem]
    exact hT
  have hY : ∀ b, (BRC minusCert * T * ARC plusCert) * gtC plusCert b =
      gtC minusCert b * (BRC minusCert * T * ARC plusCert) := by
    intro b
    have e1 : BRC minusCert * T * ARC plusCert * gtC plusCert b =
        BRC minusCert * T * gC (elems b) * ARC plusCert := by
      rw [show BRC minusCert * T * ARC plusCert * gtC plusCert b =
          BRC minusCert * (T * (ARC plusCert * BRC plusCert)) * gC (elems b) * ARC plusCert by
        unfold gtC; simp only [Matrix.mul_assoc]]
      rw [ARC_mul_BRC, hTP]
    have e2 : gtC minusCert b * (BRC minusCert * T * ARC plusCert) =
        BRC minusCert * T * gC (elems b) * ARC plusCert := by
      rw [show gtC minusCert b * (BRC minusCert * T * ARC plusCert) =
          BRC minusCert * (gC (elems b) * ((ARC minusCert * BRC minusCert) * T)) * ARC plusCert by
        unfold gtC; simp only [Matrix.mul_assoc]]
      rw [ARC_mul_BRC, hPT, ← h _ (elems_mem b), ← Matrix.mul_assoc]
    rw [e1, e2]
  have hY0 := intertwiner_eq_zero_nine _ hY
  rw [← hT, ← ARC_mul_BRC minusCert, ← ARC_mul_BRC plusCert,
    show ARC minusCert * BRC minusCert * T * (ARC plusCert * BRC plusCert) =
      ARC minusCert * (BRC minusCert * T * ARC plusCert) * BRC plusCert by
        simp only [Matrix.mul_assoc],
    hY0, Matrix.mul_zero, Matrix.zero_mul]

/-- **(3)** The two complexified golden pieces are inequivalent: the
transported actions on all listed rows admit no invertible intertwiner,
and every intertwiner is zero. -/
theorem complex_pieces_inequivalent :
    (¬ ∃ T : Matrix (Fin 3) (Fin 3) ℂ, IsUnit T ∧
      ∀ p ∈ perms, T * PhiC plusCert (gPC plusCert p) = PhiC minusCert (gPC minusCert p) * T) ∧
    (∀ T : Matrix (Fin 3) (Fin 3) ℂ,
      (∀ p ∈ perms, T * PhiC plusCert (gPC plusCert p) = PhiC minusCert (gPC minusCert p) * T) →
        T = 0) :=
  ⟨no_unit_intertwiner, intertwiner_eq_zero⟩

/-- Non-vacuity of the intertwiner hypothesis shape: the identity
intertwines each piece with itself. -/
theorem self_intertwiner (C : GoldenCert) :
    ∀ p ∈ perms, (1 : Matrix (Fin 3) (Fin 3) ℂ) * PhiC C (gPC C p) = PhiC C (gPC C p) * 1 := by
  intro p _
  rw [Matrix.one_mul, Matrix.mul_one]

/-! ## 9. Commutant: real type -/

/-- **(4)** A real `3 × 3` matrix commuting with the nine transported
tables of a piece is a real scalar. -/
theorem commutant_scalar_real (C : GoldenCert) (Y : Matrix (Fin 3) (Fin 3) ℝ)
    (h : ∀ b, Y * gt C b = gt C b * Y) : ∃ c : ℝ, Y = c • 1 := by
  have hall : ∀ Z : Matrix (Fin 3) (Fin 3) ℝ, Y * Z = Z * Y := by
    intro Z
    have hZ : Z ∈ Submodule.span ℝ (Set.range (gt C)) := by rw [span_gt_eq_top]; trivial
    induction hZ using Submodule.span_induction with
    | mem Z hZ =>
      obtain ⟨b, rfl⟩ := hZ
      exact h b
    | zero => simp
    | add Z Z' _ _ hZ hZ' => rw [Matrix.mul_add, Matrix.add_mul, hZ, hZ']
    | smul c Z _ hZ => rw [Matrix.mul_smul, Matrix.smul_mul, hZ]
  obtain ⟨c, hc⟩ := Matrix.mem_range_scalar_of_commute_single (M := Y)
    fun i j _ ↦ (hall (Matrix.single i j 1)).symm
  exact ⟨c, by rw [← hc, Matrix.scalar_apply, Matrix.smul_one_eq_diagonal]⟩

/-- **(4)** A real `3 × 3` matrix commuting with the transported actions of
all listed rows of a piece is a real scalar. -/
theorem commutant_scalar_real_perms (C : GoldenCert) (Y : Matrix (Fin 3) (Fin 3) ℝ)
    (h : ∀ p ∈ perms, Y * Phi C (gP C p) = Phi C (gP C p) * Y) : ∃ c : ℝ, Y = c • 1 :=
  commutant_scalar_real C Y fun b ↦ by
    have := h _ (elems_mem b)
    rwa [Phi_gP _ _ (elems_mem b)] at this

/-- **(4)** A real `20 × 20` matrix supported on a golden piece
(`T * P = T = P * T`) that commutes with every listed row is a real
multiple of the projector `P`: the real endomorphism algebra of the piece
is `ℝ · 1`. -/
theorem commutant_on_piece_real (C : GoldenCert) (T : Matrix (Fin 20) (Fin 20) ℝ)
    (hTP : T * projR C = T) (hPT : projR C * T = T)
    (h : ∀ p ∈ perms, T * gR p = gR p * T) : ∃ c : ℝ, T = c • projR C := by
  have hY : ∀ b, (BR C * T * AR C) * gt C b = gt C b * (BR C * T * AR C) := by
    intro b
    have e1 : BR C * T * AR C * gt C b = BR C * T * gR (elems b) * AR C := by
      rw [show BR C * T * AR C * gt C b =
          BR C * (T * (AR C * BR C)) * gR (elems b) * AR C by
        unfold gt; simp only [Matrix.mul_assoc]]
      rw [AR_mul_BR, hTP]
    have e2 : gt C b * (BR C * T * AR C) = BR C * T * gR (elems b) * AR C := by
      rw [show gt C b * (BR C * T * AR C) =
          BR C * (gR (elems b) * ((AR C * BR C) * T)) * AR C by
        unfold gt; simp only [Matrix.mul_assoc]]
      rw [AR_mul_BR, hPT, ← h _ (elems_mem b), ← Matrix.mul_assoc]
    rw [e1, e2]
  obtain ⟨c, hc⟩ := commutant_scalar_real C _ hY
  refine ⟨c, ?_⟩
  rw [← hPT, ← hTP, ← AR_mul_BR C,
    show AR C * BR C * (T * (AR C * BR C)) = AR C * (BR C * T * AR C) * BR C by
      simp only [Matrix.mul_assoc],
    hc, Matrix.mul_smul, Matrix.mul_one, Matrix.smul_mul]

/-- **(4)** Complex analogue: a complex `3 × 3` matrix commuting with the
nine complexified transported tables of a piece is a complex scalar. -/
theorem commutant_scalar_complex (C : GoldenCert) (Y : Matrix (Fin 3) (Fin 3) ℂ)
    (h : ∀ b, Y * gtC C b = gtC C b * Y) : ∃ c : ℂ, Y = c • 1 := by
  have hall : ∀ Z : Matrix (Fin 3) (Fin 3) ℂ, Y * Z = Z * Y := by
    intro Z
    have hZ : Z ∈ Submodule.span ℂ (Set.range (gtC C)) := by rw [span_gtC_eq_top]; trivial
    induction hZ using Submodule.span_induction with
    | mem Z hZ =>
      obtain ⟨b, rfl⟩ := hZ
      exact h b
    | zero => simp
    | add Z Z' _ _ hZ hZ' => rw [Matrix.mul_add, Matrix.add_mul, hZ, hZ']
    | smul c Z _ hZ => rw [Matrix.mul_smul, Matrix.smul_mul, hZ]
  obtain ⟨c, hc⟩ := Matrix.mem_range_scalar_of_commute_single (M := Y)
    fun i j _ ↦ (hall (Matrix.single i j 1)).symm
  exact ⟨c, by rw [← hc, Matrix.scalar_apply, Matrix.smul_one_eq_diagonal]⟩

/-- Instances of the real commutant statement for the two pieces. -/
theorem goldenPlus_commutant_real (Y : Matrix (Fin 3) (Fin 3) ℝ)
    (h : ∀ p ∈ perms, Y * Phi plusCert (gP plusCert p) = Phi plusCert (gP plusCert p) * Y) :
    ∃ c : ℝ, Y = c • 1 :=
  commutant_scalar_real_perms plusCert Y h

theorem goldenMinus_commutant_real (Y : Matrix (Fin 3) (Fin 3) ℝ)
    (h : ∀ p ∈ perms, Y * Phi minusCert (gP minusCert p) = Phi minusCert (gP minusCert p) * Y) :
    ∃ c : ℝ, Y = c • 1 :=
  commutant_scalar_real_perms minusCert Y h

end

end OPH.GoldenSectorComplexIrreducibility

#print axioms OPH.GoldenSectorComplexIrreducibility.span_gtC_eq_top
#print axioms OPH.GoldenSectorComplexIrreducibility.map_PhiC_SC_eq_top
#print axioms OPH.GoldenSectorComplexIrreducibility.finrank_SC
#print axioms OPH.GoldenSectorComplexIrreducibility.finrank_SCPlus
#print axioms OPH.GoldenSectorComplexIrreducibility.finrank_SCMinus
#print axioms OPH.GoldenSectorComplexIrreducibility.finrank_WCPlus
#print axioms OPH.GoldenSectorComplexIrreducibility.finrank_WCMinus
#print axioms OPH.GoldenSectorComplexIrreducibility.complex_irreducible
#print axioms OPH.GoldenSectorComplexIrreducibility.goldenPlus_complex_irreducible
#print axioms OPH.GoldenSectorComplexIrreducibility.goldenMinus_complex_irreducible
#print axioms OPH.GoldenSectorComplexIrreducibility.trace_transported
#print axioms OPH.GoldenSectorComplexIrreducibility.p5_order
#print axioms OPH.GoldenSectorComplexIrreducibility.trace_plus_p5
#print axioms OPH.GoldenSectorComplexIrreducibility.trace_minus_p5
#print axioms OPH.GoldenSectorComplexIrreducibility.no_unit_intertwiner
#print axioms OPH.GoldenSectorComplexIrreducibility.intertwiner_eq_zero
#print axioms OPH.GoldenSectorComplexIrreducibility.intertwiner_eq_zero_ambient
#print axioms OPH.GoldenSectorComplexIrreducibility.complex_pieces_inequivalent
#print axioms OPH.GoldenSectorComplexIrreducibility.commutant_scalar_real
#print axioms OPH.GoldenSectorComplexIrreducibility.commutant_on_piece_real
#print axioms OPH.GoldenSectorComplexIrreducibility.commutant_scalar_complex
#print axioms OPH.GoldenSectorComplexIrreducibility.goldenPlus_commutant_real
#print axioms OPH.GoldenSectorComplexIrreducibility.goldenMinus_commutant_real
