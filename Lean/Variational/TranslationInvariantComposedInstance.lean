import Variational.SourceToHamiltonianComposed

namespace OPH.Variational

/-!
# A translation-invariant fully-discharged composed instance
(V3, issue #731)

A declared synthetic two-state kernel, symmetric with stay weight `2/3`
and move weight `1/3`, discharges every field of the symmetry-extended
composed bundle `ComposedSymmetryData` with exact log-of-rational
literals.  The kernel is declared, and is not extracted from the
retained locally hash-pinned run; the instance is a satisfiability
receipt for the symmetry-extended conditional, and no source receipt is
claimed for it.

The enrichment is the translation-invariant member
`quadTwoPoint (2 log 2) (fun _ => -log (3/2))`, the quadratic class at
stiffness `2 log 2` with constant potential.  A translation-invariant
two-point enrichment depends on the record difference alone, so the
corner action it reproduces takes one value on the two stay pairs and
one value on the two move pairs; a strictly positive row-stochastic
kernel with such a corner action has equal stay weights and equal move
weights, which is the symmetric-kernel shape.  The declared kernel is
that shape, and the corner identity holds exactly:
`log (3/2) = -log (2/3)` on the diagonal and
`log 2 + log (3/2) = log 3 = -log (1/3)` off it
(`symLagrangian_corner`).  The committed mixing chain admits no
translation-invariant member of its curvature family: every member
carries the record-dependent base term `chainLogBase`, which the shift
does not preserve, so this instance is stated on the declared kernel and
on nothing else.

The embedded constant-0 history is simultaneously the global mode of
the path law (`constZero_global_mode`), a real fixed-endpoint
single-site minimizer at every interior junction
(`symConstZero_realMin`), a discrete Hamilton trajectory, and a carrier
of a constant Noether translation current along the whole chain
(`symComposed_receipt`).  The conserved value is `J = 0` exactly, which
is degenerate: the committed nonzero-current witness stays the affine
free path of `free_translation_nonzero_noether_witness`, which is not
an embedded chain history, and no nonzero current on a realized history
is claimed.

Boundary.  Register rows PR-05 (step-uniform reference, in the tilt
clause), PR-06 (embedding and enrichment fields), and PR-45 (derivative
and symmetry fields) are consumed as declared data; the kernel itself
is a declared witness under the same rows.  Finite exact statements; no
physical units, clocks, amplitudes, or continuum limit.
-/

/-- The constant state-0 history at every length. -/
def constZeroHistory (M : ℕ) : Fin (M + 1) → Fin 2 := fun _ => 0

/-- The declared symmetric two-state kernel: stay weight `2/3`, move
weight `1/3`.  A synthetic satisfiability witness, not a source
extraction. -/
noncomputable def symKernel : Fin 2 → Fin 2 → ℝ :=
  fun i j => if i = j then 2 / 3 else 1 / 3

/-- The translation-invariant enrichment: the quadratic class at
stiffness `2 log 2` with constant potential `-log (3/2)`. -/
noncomputable def symLagrangian : ℝ → ℝ → ℝ :=
  quadTwoPoint (2 * Real.log 2) fun _ => -Real.log (3 / 2)

/-- The first-slot derivative of the translation-invariant
enrichment. -/
noncomputable def symD1 (x y : ℝ) : ℝ := 2 * Real.log 2 * (x - y)

/-- The second-slot (momentum) derivative of the translation-invariant
enrichment. -/
noncomputable def symD2 (x y : ℝ) : ℝ := 2 * Real.log 2 * (y - x)

/-- The exact velocity solver of the translation-invariant
enrichment. -/
noncomputable def symVel (x p : ℝ) : ℝ := x + p / (2 * Real.log 2)

theorem symKernel_apply_00 : symKernel 0 0 = 2 / 3 := by
  unfold symKernel
  rw [if_pos rfl]

theorem symKernel_apply_01 : symKernel 0 1 = 1 / 3 := by
  unfold symKernel
  rw [if_neg (by decide)]

theorem symKernel_apply_10 : symKernel 1 0 = 1 / 3 := by
  unfold symKernel
  rw [if_neg (by decide)]

theorem symKernel_apply_11 : symKernel 1 1 = 2 / 3 := by
  unfold symKernel
  rw [if_pos rfl]

/-- The declared kernel is strictly positive. -/
theorem symKernel_pos (i j : Fin 2) : 0 < symKernel i j := by
  unfold symKernel
  split <;> norm_num

/-- Every entry of the declared kernel is at most the stay weight
`2/3`. -/
theorem symKernel_le (i j : Fin 2) : symKernel i j ≤ 2 / 3 := by
  unfold symKernel
  split <;> norm_num

/-- The embedded constant-0 history sits at the real value `0` on every
record. -/
theorem chainEmb_constZero (M : ℕ) (n : Fin (M + 1)) :
    chainEmb (constZeroHistory M) n = 0 := by
  simp [chainEmb, constZeroHistory]

/-- Diagonal value of the enrichment: `log (3/2)` at every equal
pair. -/
theorem symLagrangian_diag (x : ℝ) :
    symLagrangian x x = Real.log (3 / 2) := by
  simp only [symLagrangian, quadTwoPoint, quadLagrangian]
  ring

/-- Off-diagonal value at the embedded move `0 → 1`: `log 3`. -/
theorem symLagrangian_zero_one : symLagrangian 0 1 = Real.log 3 := by
  simp only [symLagrangian, quadTwoPoint, quadLagrangian]
  have hlog3 : Real.log 2 + Real.log ((3 : ℝ) / 2) = Real.log 3 := by
    rw [← Real.log_mul (by norm_num) (by norm_num)]
    norm_num
  have hval : (2 : ℝ) * Real.log 2 * ((1 : ℝ) - 0) ^ 2 / 2
        - -Real.log ((3 : ℝ) / 2)
      = Real.log 2 + Real.log ((3 : ℝ) / 2) := by ring
  linarith [hlog3, hval]

/-- Off-diagonal value at the embedded move `1 → 0`: `log 3`. -/
theorem symLagrangian_one_zero : symLagrangian 1 0 = Real.log 3 := by
  simp only [symLagrangian, quadTwoPoint, quadLagrangian]
  have hlog3 : Real.log 2 + Real.log ((3 : ℝ) / 2) = Real.log 3 := by
    rw [← Real.log_mul (by norm_num) (by norm_num)]
    norm_num
  have hval : (2 : ℝ) * Real.log 2 * ((0 : ℝ) - 1) ^ 2 / 2
        - -Real.log ((3 : ℝ) / 2)
      = Real.log 2 + Real.log ((3 : ℝ) / 2) := by ring
  linarith [hlog3, hval]

/-- **The PR-06 corner identity for the declared kernel.**  The
translation-invariant enrichment reproduces the derived corner action
`-log symKernel` at every embedded record pair, with exact log
literals. -/
theorem symLagrangian_corner (i j : Fin 2) :
    symLagrangian ((i : ℕ) : ℝ) ((j : ℕ) : ℝ)
      = -Real.log (symKernel i j) := by
  have h23 : -Real.log ((2 : ℝ) / 3) = Real.log ((3 : ℝ) / 2) := by
    rw [show ((2 : ℝ) / 3) = ((3 : ℝ) / 2)⁻¹ by norm_num, Real.log_inv,
      neg_neg]
  have h13 : -Real.log ((1 : ℝ) / 3) = Real.log (3 : ℝ) := by
    rw [show ((1 : ℝ) / 3) = ((3 : ℝ))⁻¹ by norm_num, Real.log_inv,
      neg_neg]
  fin_cases i <;> fin_cases j
  · norm_num [symKernel]
    rw [symLagrangian_diag]
    exact h23.symm
  · norm_num [symKernel]
    rw [symLagrangian_zero_one]
    exact h13.symm
  · norm_num [symKernel]
    rw [symLagrangian_one_zero]
    exact h13.symm
  · norm_num [symKernel]
    rw [symLagrangian_diag]
    exact h23.symm

/-- The uncurried derivative packet of the quadratic class at constant
potential, from the general polynomial packet. -/
theorem quadTwoPointConst_hasFDerivAt (a v : ℝ) (p : ℝ × ℝ) :
    HasFDerivAt (Function.uncurry (quadTwoPoint a fun _ => v))
      ((a * (p.1 - p.2)) • ContinuousLinearMap.fst ℝ ℝ ℝ
        + (a * (p.2 - p.1)) • ContinuousLinearMap.snd ℝ ℝ ℝ) p := by
  have hfun : Function.uncurry (quadTwoPoint a fun _ => v)
      = fun q : ℝ × ℝ =>
          -v + 0 * q.1 + 0 * q.2 + -a * (q.1 * q.2)
            + (a / 2) * (q.1 * q.1) + (a / 2) * (q.2 * q.2) := by
    funext q
    simp only [Function.uncurry, quadTwoPoint, quadLagrangian]
    ring
  have h := polyTwoPoint_hasFDerivAt (-v) 0 0 (-a) (a / 2) (a / 2) p
  have h1 : a * (p.1 - p.2) = 0 + -a * p.2 + 2 * (a / 2) * p.1 := by
    ring
  have h2 : a * (p.2 - p.1) = 0 + -a * p.1 + 2 * (a / 2) * p.2 := by
    ring
  rw [hfun, h1, h2]
  exact h

/-- **The embedded constant-0 history is a real fixed-endpoint
minimizer** of the translation-invariant enrichment at every length,
interior junction, and real replacement value: the variation gap is the
exact nonnegative quadratic `2 log 2 * x ^ 2`. -/
theorem symConstZero_realMin (M : ℕ) {k m : Fin M}
    (hkm : k.succ = m.castSucc) (x : ℝ) :
    localAction symLagrangian (chainEmb (constZeroHistory M))
      ≤ localAction symLagrangian
          (Function.update (chainEmb (constZeroHistory M)) k.succ x) := by
  have hdiff := localAction_update_diff symLagrangian
    (chainEmb (constZeroHistory M)) hkm x
  simp only [chainEmb_constZero] at hdiff
  have hgap : (symLagrangian 0 x - symLagrangian 0 0)
        + (symLagrangian x 0 - symLagrangian 0 0)
      = 2 * Real.log 2 * x ^ 2 := by
    simp only [symLagrangian, quadTwoPoint, quadLagrangian]
    ring
  have hnn : 0 ≤ 2 * Real.log 2 * x ^ 2 :=
    mul_nonneg
      (le_of_lt (mul_pos two_pos (Real.log_pos (by norm_num))))
      (sq_nonneg x)
  linarith [hdiff, hgap, hnn]

/-- **The constant-0 history is the global mode at every length**: the
entrywise bound by the stay weight `2/3` dominates every path of the
same length, with no variation restriction. -/
theorem constZero_global_mode (M : ℕ) (t : Fin (M + 1) → Fin 2) :
    pathTransitionWeight symKernel t
      ≤ pathTransitionWeight symKernel (constZeroHistory M) := by
  have heq : pathTransitionWeight symKernel (constZeroHistory M)
      = ∏ _n : Fin M, (2 / 3 : ℝ) := by
    unfold pathTransitionWeight
    refine Finset.prod_congr rfl fun n _ => ?_
    show symKernel 0 0 = 2 / 3
    exact symKernel_apply_00
  rw [heq]
  unfold pathTransitionWeight
  exact Finset.prod_le_prod (fun i _ => le_of_lt (symKernel_pos _ _))
    (fun i _ => symKernel_le _ _)

/-- The translation Noether momentum of the embedded constant-0 history
is `0` exactly, on every segment and at every length. -/
theorem symConstZero_noether_zero (M : ℕ) (seg : Fin M) :
    segmentMomentum symD2 (fun _ => 1)
      (chainEmb (constZeroHistory M)) seg = 0 := by
  simp [segmentMomentum, symD2, chainEmb, constZeroHistory]

/-- Exact Hamilton-step literals at the constant-0 state: the solved
velocity and the momentum balance both read `0`. -/
theorem symHamilton_literals :
    symVel 0 (symD2 0 0) = 0
      ∧ symD2 0 0 + symD1 0 (symVel 0 (symD2 0 0)) = 0 := by
  constructor <;> simp [symVel, symD2, symD1]

/-- **The symmetry-extended composed instance.**  The declared
symmetric kernel with the uniform initial law, the alphabet embedding,
the translation-invariant enrichment with its corner identity (PR-06
fields), the derivative packet, momentum monotonicity, exact solver,
and the unit-translation symmetry block (PR-45 fields), at every path
length. -/
noncomputable def symComposedData (M : ℕ) : ComposedSymmetryData (Fin 2)
    where
  pi := fun _ => 1 / 2
  P := symKernel
  n := M
  pi_pos := fun _ => by norm_num
  P_pos := symKernel_pos
  pi_sum := by
    rw [Fin.sum_univ_two]
    norm_num
  P_rows := fun i => by
    fin_cases i <;> norm_num [Fin.sum_univ_two, symKernel]
  emb := fun i => ((i : ℕ) : ℝ)
  Lenr := symLagrangian
  corner := symLagrangian_corner
  d1 := symD1
  d2 := symD2
  vel := symVel
  deriv_packet := quadTwoPointConst_hasFDerivAt (2 * Real.log 2)
    (-Real.log (3 / 2))
  mono := fun x u v huv => by
    simp only [symD2]
    exact mul_lt_mul_of_pos_left (by linarith)
      (mul_pos two_pos (Real.log_pos (by norm_num)))
  sec := fun x p => by
    have h2 : (2 : ℝ) * Real.log 2 ≠ 0 :=
      ne_of_gt (mul_pos two_pos (Real.log_pos (by norm_num)))
    simp only [symD2, symVel]
    field_simp
    ring
  T := fun s x => x + s
  xi := fun _ => 1
  T0 := fun x => add_zero x
  Td := fun x => by simpa using (hasDerivAt_id (0 : ℝ)).const_add x
  inv := fun s x y => by
    simp only [symLagrangian, quadTwoPoint, quadLagrangian]
    ring

/-- **The fully inhabited conditional composed receipt (issue #731).**  At the
declared symmetric kernel and every length, with the embedded
constant-0 history: the enriched local action equals the derived
log-transition action of the same kernel; the path law is the
exponential tilt of the registered step-uniform reference (row PR-05)
by that action at multiplier one; the history is a most probable
history among interior single-site alphabet variations; the embedded
path satisfies the discrete Hamilton equations of the
translation-invariant enrichment at every interior junction; and the
translation Noether momentum is one constant along the whole chain,
with exact value `0` on every segment.  One history under one declared bundle
carries the derived action, the mode, the real extremizer, the Hamilton
flow, and the conserved current.  This does not discharge source selection of
the reference, real enrichment, clock, or physical action. -/
theorem symComposed_receipt (M : ℕ) :
    localAction symLagrangian (chainEmb (constZeroHistory M))
        = InformationProjection.logTransitionAction symKernel M
            (constZeroHistory M)
      ∧ InformationProjection.markovPathLaw (fun _ => (1 / 2 : ℝ))
            symKernel M
          = InformationProjection.tilt
              (InformationProjection.stepUniformRef
                (fun _ => (1 / 2 : ℝ)) M)
              (InformationProjection.logTransitionAction symKernel M) 1
      ∧ (∀ (k m : Fin M), k.succ = m.castSucc → ∀ x : Fin 2,
          InformationProjection.markovPathLaw (fun _ => (1 / 2 : ℝ))
              symKernel M
              (Function.update (constZeroHistory M) k.succ x)
            ≤ InformationProjection.markovPathLaw (fun _ => (1 / 2 : ℝ))
                symKernel M (constZeroHistory M))
      ∧ (∀ (k m : Fin M), k.succ = m.castSucc →
          chainEmb (constZeroHistory M) m.castSucc
              = symVel (chainEmb (constZeroHistory M) k.castSucc)
                  (symD2 (chainEmb (constZeroHistory M) k.castSucc)
                    (chainEmb (constZeroHistory M) k.succ))
            ∧ symD2 (chainEmb (constZeroHistory M) k.castSucc)
                  (chainEmb (constZeroHistory M) k.succ)
                + symD1 (chainEmb (constZeroHistory M) m.castSucc)
                    (symVel (chainEmb (constZeroHistory M) m.castSucc)
                      (symD2 (chainEmb (constZeroHistory M) m.castSucc)
                        (chainEmb (constZeroHistory M) m.succ))) = 0)
      ∧ (∃ J : ℝ, ∀ seg : Fin M,
          segmentMomentum symD2 (fun _ => 1)
            (chainEmb (constZeroHistory M)) seg = J)
      ∧ ∀ seg : Fin M,
          segmentMomentum symD2 (fun _ => 1)
            (chainEmb (constZeroHistory M)) seg = 0 := by
  have hmin : ∀ (k m : Fin M), k.succ = m.castSucc → ∀ x : ℝ,
      localAction symLagrangian (chainEmb (constZeroHistory M))
        ≤ localAction symLagrangian
            (Function.update (chainEmb (constZeroHistory M)) k.succ x) :=
    fun k m hkm x => symConstZero_realMin M hkm x
  have h := source_to_hamiltonian_composed
    (symComposedData M).toComposedMechanicsData (constZeroHistory M) hmin
  have hN := source_to_hamiltonian_composed_noether (symComposedData M)
    (constZeroHistory M) hmin
  exact ⟨h.1, h.2.1, h.2.2.1, h.2.2.2.1, hN,
    symConstZero_noether_zero M⟩

end OPH.Variational

#print axioms OPH.Variational.symKernel_apply_00
#print axioms OPH.Variational.symKernel_apply_01
#print axioms OPH.Variational.symKernel_apply_10
#print axioms OPH.Variational.symKernel_apply_11
#print axioms OPH.Variational.symKernel_pos
#print axioms OPH.Variational.symKernel_le
#print axioms OPH.Variational.chainEmb_constZero
#print axioms OPH.Variational.symLagrangian_diag
#print axioms OPH.Variational.symLagrangian_zero_one
#print axioms OPH.Variational.symLagrangian_one_zero
#print axioms OPH.Variational.symLagrangian_corner
#print axioms OPH.Variational.quadTwoPointConst_hasFDerivAt
#print axioms OPH.Variational.symConstZero_realMin
#print axioms OPH.Variational.constZero_global_mode
#print axioms OPH.Variational.symConstZero_noether_zero
#print axioms OPH.Variational.symHamilton_literals
#print axioms OPH.Variational.symComposed_receipt
