# Specification of the Capacity Readback Map F

Formal acceptance specification for any candidate construction of the OPH capacity
readback map F, the object behind consistency requirement C8
([CONSISTENCY_STACK.md](../../docs/CONSISTENCY_STACK.md)), uniqueness lemma L2, and generator
G2 of the closure-ledger dependency table. G2 is the single object whose construction and
certified execution moves ledger rows CL-3, CL-4, and CL-7
([CLOSURE_LEDGER.md](../../docs/CLOSURE_LEDGER.md)). This file states what a candidate F must
be and what it must pass. It constructs nothing.

**Status: specification only. No candidate F is constructed. CL-7 is open.**

## 1. Object

The closure equation of C8 is

```
N_CRC = F(N_CRC)
```

with F the capacity readback map

```
F(N) = Cap_read( Obs( nf( U_N ) ) ).
```

Sources of the definition:

- `paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.tex`,
  Definition `def:self-closure-density` (Cosmic record-closure capacity) and Remark
  `rem:self-closure-counting-target` (Capacity fixed-point reading); the D6 lane of
  Theorem `thm:lambda-stack` consumes the fixed point.
- `paper/tex_fragments/OBSERVERS_SYNTHESIS_SECTIONS.tex`, global self-closure section
  (the finite readback form `F_r(N) = Cap_read(Obs(nf_{r,N}(U_{r,N})))`, the contraction
  condition `0 <= F'(N) <= kappa < 1`, and the count representation).

### 1.1 Domain and codomain

- **Domain.** An admissible interval `I ⊂ (0, ∞)` of active horizon capacities, measured
  in nats: `N = log dim H_∂,N`, where `H_∂,N` is the observer-facing horizon record
  Hilbert space at capacity `N`. At finite repair cutoff `r` the finite map `F_r` is
  defined on the finite normal-form grammar; the specified `F` is the cofinal refinement
  limit of `F_r` on `I`, and existence of that limit on `I` is part of the construction
  obligation (synthesis section, refinement-limit paragraph).
- **Codomain.** `(0, ∞)` in the same nat units. The self-map requirement `F(I) ⊆ I` is
  a certified property (Section 3), not a definition-level assumption.

### 1.2 The three factors

`F = Cap_read ∘ Obs ∘ nf`, applied to the universe candidate `U_N` (the OPH universe
candidate supplied with active boundary capacity `N`: horizon size, accessible records,
dilution rate, checkpoint continuation, observer-supporting capacity, late-time de Sitter
boundary; synthesis section, trial-capacity paragraph).

1. **nf. Quotient normal form.** The quotient normal-form map on the declared
   terminating and confluent repair surface, supplied by the consensus branch
   (`def:self-closure-density`; consensus paper carries the termination/confluence
   package). `nf(U_N)` is the public quotient of the candidate: the record structure
   that survives overlap agreement and repair. Well-definedness of `nf` is exactly
   termination plus confluence of the repair rewriting at the declared cutoff.
2. **Obs. Stable observer-sector selection.** `Obs` selects the stable self-reading
   observer sector of the quotient: the subfederation of normal forms that support at
   least one stable observer/checkpoint chain and carry the recovered local package
   (`def:self-closure-density`; the `Ω^sc_N` membership predicate lists the same four
   clauses). `Obs` is a selection, not an optimization: its output is determined by the
   declared stability predicate, with no tunable threshold.
3. **Cap_read. Reconstructed capacity readout.** `Cap_read` returns the capacity the
   selected sector reconstructs from inside: from internal geometry, causal
   accessibility, stable records, and self-read closure (synthesis section). Observers
   in `U_N` never see the supplied parameter `N`; `Cap_read` is the inside estimate of
   the boundary capacity, in the same nat normalization `log dim`.

The physical content of the map: supplied raw Hilbert-space room becomes readable public
capacity only insofar as it is reachable, recordable, self-readable, stable under repair,
coupled to observer continuation, and visible on the quotient (synthesis section,
contraction-meaning paragraph). `F(N) < N` slack and `F(N) > N` deficit are both
admissible off the fixed point; the closure equation demands the point with neither.

### 1.3 Count representation and the density log|Ω^sc_N| − N

The screen-normalized count representation of the same target
(`def:self-closure-density`; synthesis section, finite-count paragraphs): let `Ω^sc_N` be
the terminal normal forms at capacity `N` that are repair-closed, support at least one
stable observer/checkpoint subfederation, carry the recovered local package, and whose
own horizon record surface reads back capacity `N`. Define

```
Π(N) = |Ω^sc_N| / dim H_∂,N = |Ω^sc_N| · e^(−N),
ℓ(N) = log Π(N) = log|Ω^sc_N| − N.
```

The subtraction of `N` is not a tunable penalty; it is division by the full screen
Hilbert-space dimension `e^N`. The input-free selector target is

```
N_star = MAR argmax_N [ log|Ω^sc_N| − N ],
```

with stationarity map `T_η(N) = N + η·ℓ'(N)` and fixed-point condition `ℓ'(N) = 0`,
carried by the derivative-sign certificate `H_N(N) = ℓ'(N)`, `H_N(N₀) = 0`, `H_N' < 0`
on the admissible interval, or the stronger concavity bound `−M ≤ ℓ'' ≤ −m < 0`
(synthesis section, pressure-certificate paragraph).

A candidate F therefore carries a **two-representation coherence obligation**: the
readback fixed point `N_CRC = F(N_CRC)` and the count-density stationary point `N_star`
are two displays of one target, and the construction must either certify
`N_CRC = N_star` on the admissible interval or register the discrepancy as a new closure
row. Silent coexistence of two values is excluded by SL-4 (one N; CL-3 reading rule).

## 2. Required properties

A candidate construction is admissible for certification only if all of P1–P5 hold.

- **P1. Well-definedness at finite N.** At every finite cutoff `r` and admissible `N`:
  the repair surface terminates and is confluent (so `nf` is a function), the stability
  predicate of `Obs` is decidable on the quotient, and `Cap_read` is total on the
  selected sector. `F_r` is then a total function of `N`, and the refinement limit `F`
  exists on `I` (cofinal limit, synthesis section). A construction that is only defined
  on a measure-zero or unstated subset of `I` is not a candidate.
- **P2. Monotonicity regime.** On the admissible interval, `F` is differentiable (or
  interval-Lipschitz) with `0 ≤ F'(N)`; supplied capacity never decreases readable
  capacity. The regime statement is the synthesis-section condition
  `0 ≤ F'(N) ≤ κ < 1` on `I`.
- **P3. Boundedness and growth bounds.** Bounds sufficient for a contraction interval
  to exist: a bracketing pair `a < b` in `I` with `F(a) ≥ a` and `F(b) ≤ b` together
  with P2 gives `F([a,b]) ⊆ [a,b]`; the derivative ceiling `κ < 1` on `[a,b]` then makes
  `F` a Banach contraction there (L2). A sufficient growth form: if `Cap_read` reads a
  declared normalization times the log of the sector count, and the sector count grows
  at most exponentially with rate `ρ` in `N`, then `F' ≤ (normalization)·ρ`, and the
  ceiling is a computable product. Whatever the form, the bound must be emitted with the
  construction, not asserted.
- **P4. Count-density coherence.** The obligation of Section 1.3: certified agreement
  of the readback fixed point with the `MAR argmax` of `ℓ(N)`, or a registered
  discrepancy row. The derivative-sign certificate for `H_N` is the accepted carrier.
  For the coupled candidate, the reduction of this obligation to the single count-side
  premise CP-4 (stationarity at load balance), and its degenerate closure under the
  aggregate reading of membership clause 4, are recorded in
  [P4_COUPLED_COHERENCE_REDUCTION_2026-07-17.md](P4_COUPLED_COHERENCE_REDUCTION_2026-07-17.md).
- **P5. Non-triviality.** `F` is not the identity and not a constant pinned to any
  measured value. A map with `F' ≡ 1` on an interval certifies nothing (every point is
  fixed); a constant map wired to the SL-4 estimate is a blindness violation (Section
  3), not a construction.

## 3. Blindness requirement (dependency-cone audit V-08)

The evaluation cone of a candidate F must not contain:

- measured Λ, measured G, or any quantity computed from them;
- the SL-4 estimate `N = 3.31e122` ([STRANGE_LOOP_PRINCIPLES.md](../../docs/STRANGE_LOOP_PRINCIPLES.md), Estimated values),
  including every code alias of it (e.g. `DEFAULT_N_SCR`);
- the CL-3 electroweak-bridge value `N_EW = π·exp(6π/(P·α_U)) = 3.53e122`;
- any ledger residual, basin location, or downstream display derived from the above.

The Λ-estimated capacity is the basin the executed output is compared against, stage 1
of the basin-then-contract protocol (STRANGE_LOOP_PRINCIPLES.md rule 3; CLOSURE_LEDGER.md protocol
table, N row). It is never an input to the evaluation. Direction of inference is
permanent (STRANGE_LOOP_PRINCIPLES.md rule 1): a map revised after seeing the basin converts the basin to
an input and the landing to a retrodiction (rule 4). Enforcement is mechanical: the V-08
dependency-cone audit over the candidate's import and data cone, executed before the
first comparison. The declared detuning constant `P` (SL-3) may appear in the cone only
where the construction derives it from the source side; the CL-3 bridge *value* may not.

## 4. Acceptance tests

CL-7 closes only when a candidate implementation passes all of A1–A7. The certificate
schema mirrors the P-side machinery in
[`code/P_derivation/`](../P_derivation/) (`fixed_point_certificate.py`,
`runtime/fine_structure_interval_certificate_current.json`).

- **A1. Finite executability.** `F_r(N)` executes at declared finite cutoffs and
  returns identical output on repeated runs (deterministic; no timestamps or run
  identifiers in artifact content).
- **A2. Refinement stability.** Emitted evidence that `F_r → F` on `I`: monotone or
  Cauchy behavior of `F_r(N)` in `r` at declared probe capacities, with stated bounds.
- **A3. Self-map enclosure.** Interval-evaluated `F(I') ⊆ I'` for a stated certificate
  interval `I' ⊆ I`, computed with outward rounding.
- **A4. Contraction.** Interval-evaluated derivative (or Lipschitz) enclosure on `I'`
  with certified `sup |F'| ≤ L < 1`. A3 + A4 constitute the stage-2 Banach certificate:
  existence and uniqueness of `N_CRC` in `I'` (L2).
- **A5. Fixed-point enclosure.** A certified enclosure `[N_lo, N_hi] ∋ N_CRC` with
  stated width, produced by interval iteration inside `I'`.
- **A6. Blindness audit.** V-08 dependency-cone audit record over the evaluation cone
  (Section 3), attached to the certificate, produced before the basin comparison.
- **A7. Landing verdict.** The blind output compared once against the SL-4 basin
  (stage 3). Inside: CL-7 closes, and CL-3/CL-4 become evaluable at the one certified N.
  Outside: the verdict is recorded permanently (STRANGE_LOOP_PRINCIPLES.md rule 7); the row does not
  close by relabeling.

Certificate JSON schema (minimum keys, mirroring the P-side artifact):

```
artifact                      identifier string
status                        pass / rejected / conditional
interval_backend              { library, precision, rounding }
contraction_certificate       { interval {lo,hi}, image {lo,hi}, self_map_pass,
                                derivative_enclosure {lo,hi}, lipschitz_bound_L,
                                lipschitz_pass, banach_pass }
fixed_point                   { enclosure {lo,hi}, width }
count_density                 P4 coherence record or registered discrepancy
blindness                     { inputs, reads_measured_lambda, reads_sl4_estimate,
                                dependency_cone }
promotion_allowed             boolean; theorem-grade promotion requires a
                              directed-rounding backend (Arb/MPFI class), matching the
                              P-side promotion gate
```

## 5. What this specification does not do

It does not construct `U_N`, `nf`, `Obs`, or `Cap_read` for the physical grammar; it
does not evaluate F; it does not move any ledger row. The executable schema
demonstration in this directory ([`toy_readback.py`](toy_readback.py)) runs the full
certificate pipeline on a declared toy model with no physical content, so that the
acceptance machinery exists and is tested before any candidate does.

**Status: specification only. No candidate F is constructed. CL-7 open; CL-3 and CL-4
wait on the same object (generator G2).**
