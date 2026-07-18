import ObserverPatchHolography.CollarStates

/-!
# State-side collar states — T1: the reversal exhibit (#544)

T0 (`CollarStates.lean`) showed the paper's current admissible-channel
requirements cannot force the central-interface collar clause: the
identity channel keeps the non-central family `{uuC, XXC}` at closure
defect zero.  T1 sharpens the boundary from the other side: over ℂ the
ℤ-obstruction of `no_integral_flux_retraction` *reverses* — the
Hilbert–Schmidt-orthogonal conditional expectation

  `Eflux : m ↦ (tr m)/4 • 1 + (tr (uuC·m))/4 • uuC`

onto the flux subalgebra `span{1, uuC}` exists, is an admissible channel
(linear, positive, trace-preserving), fixes the flux sector pointwise,
and genuinely annihilates the coupling `XXC`.  And *even so*, coarse-
graining the realized MaxEnt branch of `{uuC, XXC}` by `Eflux` lands
inside the same exponential family with closure defect exactly zero:
the transport identity

  `Eflux (gibbsM ![uuC, XXC] lam) = gibbsM ![uuC, XXC] ![lam 0, 0]`

(coordinate erasure of the coupling multiplier) exhibits the coarse-
grained state as the central member at `lam' = ![lam 0, 0]`.  Possession
of the excluding map still does not *force* the clause — it merely
*deselects* the coupling on the coarse side.  This is the force-vs-select
gap: T1 re-characterizes #544, it does not close it; strengthening the
admissible class at the axiom level (T2) stays out of scope.

Everything here rides T0's banked spectral machinery (`specProjP`,
`gibbsM_spectral`, the family Klein inequality); no new analytic
machinery, no `Matrix.exp`, no Kubo–Mori.
-/

namespace OPH

open Matrix
open scoped ComplexOrder

/-! ## T1-S1 — the flux conditional expectation and the transport identity -/

/-- The flux conditional expectation as a raw linear map: the
    Hilbert–Schmidt-orthogonal projection of the collar algebra onto the
    flux subalgebra `span{1, uuC}` (the map whose ℤ-integral avatar was
    impossible by `no_integral_flux_retraction`). -/
noncomputable def EfluxL : CollarC →ₗ[ℂ] CollarC where
  toFun m := (m.trace / 4) • (1 : CollarC) + ((uuC * m).trace / 4) • uuC
  map_add' m n := by
    rw [Matrix.trace_add, mul_add, Matrix.trace_add, add_div, add_div,
      add_smul, add_smul]
    abel
  map_smul' c m := by
    simp only [Matrix.mul_smul, Matrix.trace_smul, smul_eq_mul,
      RingHom.id_apply, smul_add, smul_smul, mul_div_assoc]

@[simp] theorem EfluxL_apply (m : CollarC) :
    EfluxL m = (m.trace / 4) • (1 : CollarC) + ((uuC * m).trace / 4) • uuC :=
  rfl

/-- `uuC` acts by its sign on each joint spectral projection. -/
theorem uuC_mul_specProjP (p : Fin 2 × Fin 2) :
    uuC * specProjP p = ((sgnR p.1 : ℝ) : ℂ) • specProjP p := by
  rw [specProjP, ← mul_assoc, uuC_mul_quP, Matrix.smul_mul]

theorem uuC_mul_specProjP_trace (p : Fin 2 × Fin 2) :
    (uuC * specProjP p).trace = ((sgnR p.1 : ℝ) : ℂ) := by
  rw [uuC_mul_specProjP, Matrix.trace_smul, specProjP_trace, smul_eq_mul,
    mul_one]

/-- The action of the flux expectation on the joint spectral basis:
    each sector projection is sent to half the corresponding `uuC`
    half-projection — the coupling direction is forgotten, the charge
    direction is kept. -/
theorem Eflux_specProjP (p : Fin 2 × Fin 2) :
    EfluxL (specProjP p) = (2 : ℂ)⁻¹ • quP p.1 := by
  rw [EfluxL_apply, specProjP_trace, uuC_mul_specProjP_trace, quP,
    smul_smul, smul_add, smul_smul]
  congr 1
  · congr 1
    norm_num
  · congr 1
    rw [div_eq_mul_inv, mul_comm]
    norm_num [mul_comm]

/-- The `uuC` half-projections are the fibre sums of the joint spectral
    projections. -/
theorem sum_specProjP_fibre (s : Fin 2) :
    ∑ t : Fin 2, specProjP (s, t) = quP s := by
  simp only [specProjP]
  rw [← Finset.mul_sum, sum_qxP, mul_one]

/-- Coordinate erasure at the eigenvalue level: killing the second
    multiplier removes the coupling contribution to every sector
    eigenvalue. -/
theorem eigval_erase (lam : Fin 2 → ℝ) (p : Fin 2 × Fin 2) :
    eigval ![lam 0, 0] p = lam 0 * sgnR p.1 := by
  simp [eigval]

/-- The Gibbs weights of the erased family are constant along the
    coupling fibre. -/
theorem gibbsProb_erase (lam : Fin 2 → ℝ) (p : Fin 2 × Fin 2) :
    gibbsProb ![lam 0, 0] p
      = Real.exp (-(lam 0 * sgnR p.1)) / partitionZ ![lam 0, 0] := by
  rw [gibbsProb, eigval_erase]

/-- The partition function of the commuting family factorizes over the
    two spectral directions. -/
theorem partitionZ_factor (lam : Fin 2 → ℝ) :
    partitionZ lam
      = (∑ s : Fin 2, Real.exp (-(lam 0 * sgnR s)))
        * (∑ t : Fin 2, Real.exp (-(lam 1 * sgnR t))) := by
  rw [partitionZ, Fintype.sum_prod_type, Finset.sum_mul_sum]
  refine Finset.sum_congr rfl fun s _ => Finset.sum_congr rfl fun t _ => ?_
  rw [eigval, neg_add, Real.exp_add]

/-- **The marginal identity** — the scalar heart of the transport: the
    halved coupling-fibre marginal of the joint Gibbs weights equals the
    erased-family Gibbs weight.  Pure `Real.exp` arithmetic on the
    factorized partition functions. -/
theorem gibbsProb_marginal (lam : Fin 2 → ℝ) (s : Fin 2) :
    (∑ t : Fin 2, gibbsProb lam (s, t)) / 2
      = Real.exp (-(lam 0 * sgnR s)) / partitionZ ![lam 0, 0] := by
  have hZU : (0 : ℝ) < ∑ u : Fin 2, Real.exp (-(lam 0 * sgnR u)) :=
    Finset.sum_pos (fun _ _ => Real.exp_pos _) Finset.univ_nonempty
  have hZX : (0 : ℝ) < ∑ t : Fin 2, Real.exp (-(lam 1 * sgnR t)) :=
    Finset.sum_pos (fun _ _ => Real.exp_pos _) Finset.univ_nonempty
  have hterm : ∀ t : Fin 2, gibbsProb lam (s, t)
      = Real.exp (-(lam 0 * sgnR s)) * Real.exp (-(lam 1 * sgnR t))
        / partitionZ lam := by
    intro t
    rw [gibbsProb, eigval, neg_add, Real.exp_add]
  have hnum : ∑ t : Fin 2, gibbsProb lam (s, t)
      = Real.exp (-(lam 0 * sgnR s))
        * (∑ t : Fin 2, Real.exp (-(lam 1 * sgnR t))) / partitionZ lam := by
    rw [Finset.sum_congr rfl fun t _ => hterm t, ← Finset.sum_div,
      ← Finset.mul_sum]
  have hZ' : partitionZ ![lam 0, 0]
      = (∑ u : Fin 2, Real.exp (-(lam 0 * sgnR u))) * 2 := by
    rw [partitionZ_factor]
    norm_num [Fin.sum_univ_two, Matrix.cons_val_zero, Matrix.cons_val_one,
      Matrix.head_cons]
  rw [hnum, partitionZ_factor lam, hZ']
  field_simp [hZU.ne', hZX.ne']

/-- **T1 transport identity (the crux)**: coarse-graining the realized
    Gibbs state of the non-central family by the flux expectation lands
    exactly on the family member with the coupling multiplier erased,
    `lam' = ![lam 0, 0]`.  The I-projection infimum is therefore hit at
    an explicit central point — no optimization, no Kubo–Mori. -/
theorem Eflux_transport (lam : Fin 2 → ℝ) :
    EfluxL (gibbsM SFam lam) = gibbsM SFam ![lam 0, 0] := by
  rw [gibbsM_spectral lam, map_sum, gibbsM_spectral ![lam 0, 0],
    Fintype.sum_prod_type, Fintype.sum_prod_type]
  refine Finset.sum_congr rfl fun s _ => ?_
  have hL : ∑ t : Fin 2,
        EfluxL (((gibbsProb lam (s, t) : ℝ) : ℂ) • specProjP (s, t))
      = (((∑ t : Fin 2, gibbsProb lam (s, t)) / 2 : ℝ) : ℂ) • quP s := by
    have step : ∀ t : Fin 2,
        EfluxL (((gibbsProb lam (s, t) : ℝ) : ℂ) • specProjP (s, t))
          = ((gibbsProb lam (s, t) / 2 : ℝ) : ℂ) • quP s := by
      intro t
      rw [map_smul, Eflux_specProjP, smul_smul]
      congr 1
      rw [Complex.ofReal_div, div_eq_mul_inv]
      norm_num
    rw [Finset.sum_congr rfl fun t _ => step t, ← Finset.sum_smul]
    congr 1
    rw [← Complex.ofReal_sum]
    congr 1
    exact (Finset.sum_div _ _ _).symm
  have hR : ∑ t : Fin 2,
        ((gibbsProb ![lam 0, 0] (s, t) : ℝ) : ℂ) • specProjP (s, t)
      = ((Real.exp (-(lam 0 * sgnR s)) / partitionZ ![lam 0, 0] : ℝ) : ℂ)
          • quP s := by
    have step : ∀ t : Fin 2,
        ((gibbsProb ![lam 0, 0] (s, t) : ℝ) : ℂ) • specProjP (s, t)
          = ((Real.exp (-(lam 0 * sgnR s)) / partitionZ ![lam 0, 0] : ℝ) : ℂ)
              • specProjP (s, t) := by
      intro t
      rw [gibbsProb_erase]
    rw [Finset.sum_congr rfl fun t _ => step t, ← Finset.smul_sum,
      sum_specProjP_fibre]
  rw [hL, hR, gibbsProb_marginal]

/-! ## T1-S2 — the action table of the flux expectation

The five lemmas of the scoping report plus the sixth (the `uuC·XXC`
cross-term, which the Gibbs exponential carries): the expectation is
unital, fixes the charge, maps the ℤ-impossible sector projector to
itself, and annihilates both the coupling and the charge–coupling
cross-term. -/

theorem Eflux_unital : EfluxL 1 = 1 := by
  rw [EfluxL_apply, mul_one, trace_one_collarC, uuC_trace]
  norm_num

theorem Eflux_fixes_uu : EfluxL uuC = uuC := by
  rw [EfluxL_apply, uuC_trace, uuC_sq, trace_one_collarC]
  norm_num

theorem Eflux_kills_XX : EfluxL XXC = 0 := by
  rw [EfluxL_apply, XXC_trace, uuC_mul_XXC_trace]
  norm_num

/-- The sixth lemma (not in the report's five): the charge–coupling
    cross-term is also annihilated — `uuC·(uuC·XXC) = XXC` is traceless,
    so the Gibbs exponential's cross-term cannot strand any payload. -/
theorem Eflux_kills_uu_mul_XX : EfluxL (uuC * XXC) = 0 := by
  rw [EfluxL_apply, uuC_mul_XXC_trace, ← mul_assoc, uuC_sq, one_mul,
    XXC_trace]
  norm_num

/-- The ℤ-impossible value, now legal: the flux expectation fixes the
    sector projector — a linear consequence of unitality and charge
    fixing, not an entrywise computation. -/
theorem Eflux_maps_pPlus : EfluxL pPlusC = pPlusC := by
  rw [pPlusC, map_smul, map_add, Eflux_unital, Eflux_fixes_uu]

/-! ## T1-S3 — positivity, trace preservation, and the channel

Positivity via the spectral form `E(m) = ∑ₛ (tr(quPₛ·m)/2) • quPₛ`:
each coefficient is the trace of a projector-compression of `m`, hence
nonnegative.  Positivity + trace preservation suffice for admissibility;
complete positivity is true but not required (no Choi theatre). -/

theorem quP_idem (s : Fin 2) : quP s * quP s = quP s := by
  rw [quP_mul, if_pos rfl]

theorem quP_posSemidef (s : Fin 2) : (quP s).PosSemidef := by
  have h : quP s = (quP s)ᴴ * quP s := by
    rw [(quP_isHermitian s).eq, quP_idem]
  rw [h]
  exact Matrix.posSemidef_conjTranspose_mul_self _

/-- The spectral form of the flux expectation on the `uuC`
    half-projections. -/
theorem EfluxL_spectral (m : CollarC) :
    EfluxL m = ∑ s : Fin 2, ((quP s * m).trace / 2) • quP s := by
  have htr : ∀ s : Fin 2, (quP s * m).trace
      = (2 : ℂ)⁻¹ * m.trace
        + ((2 : ℂ)⁻¹ * ((sgnR s : ℝ) : ℂ)) * (uuC * m).trace := by
    intro s
    rw [quP, Matrix.smul_mul, add_mul, one_mul, Matrix.smul_mul,
      Matrix.trace_smul, Matrix.trace_add, Matrix.trace_smul]
    simp only [smul_eq_mul]
    ring
  rw [Fin.sum_univ_two, htr 0, htr 1, EfluxL_apply, quP, quP]
  simp only [sgnR_zero, sgnR_one, Complex.ofReal_one, Complex.ofReal_neg]
  match_scalars <;> ring

/-- Compression of a state by a `uuC` half-projection has nonnegative
    trace. -/
theorem trace_quP_mul_nonneg {m : CollarC} (hm : m.PosSemidef) (s : Fin 2) :
    0 ≤ (quP s * m).trace := by
  have key : (quP s * m).trace = (quP s * m * quP s).trace := by
    rw [Matrix.trace_mul_cycle (quP s) m (quP s), quP_idem]
  have hconj : (quP s * m * quP s).PosSemidef := by
    have h := hm.conjTranspose_mul_mul_same (B := quP s)
    rwa [(quP_isHermitian s).eq] at h
  rw [key]
  exact hconj.trace_nonneg

/-- Nonnegative complex scalars preserve positive semidefiniteness. -/
theorem posSemidef_smul_of_nonneg {c : ℂ} (hc : 0 ≤ c) {M : CollarC}
    (hM : M.PosSemidef) : (c • M).PosSemidef := by
  rw [Complex.nonneg_iff] at hc
  have hcr : c = ((c.re : ℝ) : ℂ) := Complex.ext (by simp) (by simp [hc.2])
  rw [hcr]
  exact posSemidef_ofReal_smul hc.1 hM

theorem EfluxL_pos {m : CollarC} (hm : m.PosSemidef) :
    (EfluxL m).PosSemidef := by
  have h2 : (0 : ℂ) ≤ (2 : ℂ)⁻¹ := by
    rw [Complex.nonneg_iff]
    norm_num
  rw [EfluxL_spectral, Fin.sum_univ_two]
  refine Matrix.PosSemidef.add ?_ ?_ <;>
  · refine posSemidef_smul_of_nonneg ?_ (quP_posSemidef _)
    rw [div_eq_mul_inv]
    exact complex_mul_nonneg (trace_quP_mul_nonneg hm _) h2

theorem EfluxL_trace (m : CollarC) : (EfluxL m).trace = m.trace := by
  rw [EfluxL_apply, Matrix.trace_add, Matrix.trace_smul, Matrix.trace_smul,
    trace_one_collarC, uuC_trace, smul_eq_mul, smul_eq_mul]
  ring

/-- **The flux expectation is an admissible channel** — linear, positive,
    trace-preserving; exactly the paper's named requirements, nothing
    smuggled. -/
noncomputable def EfluxChannel : AdmissibleChannel where
  Φ := EfluxL
  pos := fun _ hm => EfluxL_pos hm
  tracePreserving := EfluxL_trace

@[simp] theorem EfluxChannel_apply (m : CollarC) :
    EfluxChannel.Φ m = EfluxL m := rfl

theorem XXC_ne_zero : XXC ≠ 0 := by
  intro h
  have h1 : XXC (0, 1) (1, 0) = (1 : ℂ) := by
    have hx : XX (0, 1) (1, 0) = 1 := by decide
    simp [XXC, hx]
  rw [h] at h1
  simp at h1

/-- **Deselection, not exclusion** (false-receipt guard): the flux
    expectation genuinely kills the coupling that the identity channel
    keeps fixed — the T1 payload below cannot be read as T0 in a wig. -/
theorem EfluxChannel_deselects_XXC :
    EfluxChannel.Φ XXC = 0 ∧ idChannel.Φ XXC = XXC ∧ XXC ≠ 0 :=
  ⟨Eflux_kills_XX, rfl, XXC_ne_zero⟩

/-! ## T1-S4 — the reversal receipt and the payload -/

/-- The ℂ-span of `{1, uuC}` is a subring — `uuC² = 1` closes
    multiplication.  (Over ℤ this same closure was the wall the
    retraction could not cross; over ℂ it is the home of `Eflux`.) -/
def fluxSpanC : Subring CollarC where
  carrier := {m | ∃ a b : ℂ, m = a • 1 + b • uuC}
  zero_mem' := ⟨0, 0, by simp⟩
  one_mem' := ⟨1, 0, by simp⟩
  add_mem' := by
    rintro x y ⟨a1, b1, rfl⟩ ⟨a2, b2, rfl⟩
    exact ⟨a1 + a2, b1 + b2, by module⟩
  neg_mem' := by
    rintro x ⟨a1, b1, rfl⟩
    exact ⟨-a1, -b1, by module⟩
  mul_mem' := by
    rintro x y ⟨a1, b1, rfl⟩ ⟨a2, b2, rfl⟩
    refine ⟨a1 * a2 + b1 * b2, a1 * b2 + b1 * a2, ?_⟩
    rw [add_mul, mul_add, mul_add, smul_mul_smul_comm, smul_mul_smul_comm,
      smul_mul_smul_comm, smul_mul_smul_comm, uuC_sq]
    simp only [one_mul, mul_one]
    module

/-- The flux expectation fixes the flux sector pointwise: every element
    of `Z(K)` is an actual fixed point, not merely a member of the
    range. -/
theorem EfluxL_fixes_fluxC : ∀ m ∈ FluxC, EfluxL m = m := by
  rintro m ⟨hK, -⟩
  have hle : Subring.closure {uuC} ≤ fluxSpanC := by
    rw [Subring.closure_le]
    intro g hg
    rw [Set.mem_singleton_iff] at hg
    subst hg
    exact ⟨0, 1, by simp⟩
  obtain ⟨a, b, rfl⟩ := hle hK
  rw [map_add, map_smul, map_smul, Eflux_unital, Eflux_fixes_uu]

/-- **The reversal receipt**: over ℂ the excluding map exists — an
    admissible channel that fixes the flux sector pointwise, is unital,
    and annihilates the coupling.  The witness is the concrete
    `EfluxChannel`, not an unrelated existential.  This is the ℂ-negation
    of the ℤ-obstruction `no_integral_flux_retraction`. -/
theorem flux_expectation_exists_over_C :
    ∃ E : AdmissibleChannel,
      (∀ m ∈ FluxC, E.Φ m = m) ∧ E.Φ 1 = 1 ∧ E.Φ XXC = 0 :=
  ⟨EfluxChannel, EfluxL_fixes_fluxC, Eflux_unital, Eflux_kills_XX⟩

/-- **T1 payload — possession of the excluding map still does not
    force.**  Coarse-graining the realized MaxEnt branch of the
    non-central family `{uuC, XXC}` by the genuine flux expectation
    leaves the closure defect exactly zero at every point of the branch:
    the coarse-grained state IS the central family member at
    `lam' = ![lam 0, 0]` (the transport identity), so the I-projection
    infimum is attained at zero.  Together with
    `EfluxChannel_deselects_XXC` this is deselection, not exclusion —
    the excluding map re-characterizes the boundary of #544; it does not
    close it. -/
theorem Eflux_does_not_force (lam : Fin 2 → ℝ) :
    closureDefect ![uuC, XXC] EfluxChannel lam = 0 := by
  have hSF : (![uuC, XXC] : Fin 2 → CollarC) = SFam := rfl
  rw [hSF]
  have hbdd : BddBelow (Set.range fun lam' =>
      relEntropyM (EfluxChannel.Φ (gibbsM SFam lam)) (gibbsM SFam lam')) := by
    refine ⟨0, ?_⟩
    rintro x ⟨lam', rfl⟩
    simp only [EfluxChannel_apply]
    rw [Eflux_transport]
    exact relEntropyM_gibbs_nonneg _ lam'
  refine le_antisymm ?_ ?_
  · have h := closureDefect_le_relEntropy_apply hbdd ![lam 0, 0]
    simp only [EfluxChannel_apply] at h
    rwa [Eflux_transport, relEntropyM_self] at h
  · refine le_ciInf fun lam' => ?_
    simp only [EfluxChannel_apply]
    rw [Eflux_transport]
    exact relEntropyM_gibbs_nonneg _ lam'

/-! ## Axiom audit

Expected footprint for every entry: `[propext, Classical.choice,
Quot.sound]`.  No `sorry`, no `native_decide`, no extra axioms. -/

#print axioms Eflux_specProjP
#print axioms gibbsProb_marginal
#print axioms Eflux_transport
#print axioms Eflux_unital
#print axioms Eflux_fixes_uu
#print axioms Eflux_kills_XX
#print axioms Eflux_kills_uu_mul_XX
#print axioms Eflux_maps_pPlus
#print axioms EfluxL_pos
#print axioms EfluxL_trace
#print axioms EfluxChannel
#print axioms EfluxChannel_deselects_XXC
#print axioms EfluxL_fixes_fluxC
#print axioms flux_expectation_exists_over_C
#print axioms Eflux_does_not_force

end OPH
