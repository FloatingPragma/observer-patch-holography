# OPH Closure Ledger

The equations the strange-loop principle set ([STRANGE_LOOP_PRINCIPLES.md](STRANGE_LOOP_PRINCIPLES.md)) requires to hold exactly, with
their evaluation status. The proof spine that consumes these rows is [PROOF_SPINE.md](PROOF_SPINE.md). A residual bears on the hypothesis only after the map, physical observable, convention, and uncertainty packet needed for the comparison are complete. Rows are permanent; verdicts append. Last audited: 2026-07-16 (manual; CI regeneration is fix item SLP-02).

| # | Closure | Required | Current status | Residual (relative) | Pull / evaluation status | Resolution path |
|---|---|---|---|---|---|---|
| CL-1 | P-loop, source chain: P = φ + √π/A_T(P) with source-only transport | exact after the declared map is completed | incomplete declared-map output: α⁻¹ = 136.994835177413 (interval-certified unique fixed point, width 7.2×10⁻²⁴); the Ward-projected transport and same-scheme remainder are absent | 3.0×10⁻⁴ diagnostic gap | not a physical pull; map incomplete | register a detached successor whose complete method is blind to withheld data or clean-room audited, emit a target-blind transport function over P, and certify the completed CL-1 map (SLP-03) |
| CL-2 | P-loop with gauge-width term, self-consistent map | exact after the declared map is completed | incomplete declared-map output: certified unique fixed point α⁻¹ = 137.035660136946577 vs 137.035999177(21). The previously displayed 137.0359595 is a mixed-provenance packet, not a fixed point of one declared map | 2.5×10⁻⁶ diagnostic gap | not a physical pull; map incomplete | use the same target-blind transport function but solve and certify CL-2 independently; its point target differs from CL-1 by α_U(P) |
| CL-3 | One capacity (SL-4): N_EW = π·exp(6π/(P·α_U)) vs N from Λ | exact only after F and premises CP-1 to CP-3 are discharged | conditional bridge value 3.53×10¹²² vs nominal Λ-readout value 3.31×10¹²²; F and CP-1 to CP-3 remain open | nominal 6.6×10⁻²; posterior propagation pending | not evaluated; the joint Planck posterior has not been propagated | construct F unconditionally, discharge CP-1 to CP-3, then propagate the Planck posterior through the Λ readout (SLP-04) |
| CL-4 | Hierarchy bridge at the one N: α_U·log(N/π) = 6π/P | conditional on the same capacity premises | nominal 11.5546 vs 11.5573 at the Λ-located N | nominal 2.3×10⁻⁴ (log scale); posterior propagation pending | not evaluated | same object as CL-3 |
| CL-5 | Forward electroweak emission at SL-3 pixel: (M_W, M_Z) | a common physical observable and convention | not evaluable. The emitted (80.330, 91.119) GeV pair is a running/tree chart coordinate. The quoted (80.3692(133), 91.1880(20)) GeV references are stale PDG 2025 mass-dependent-width Breit-Wigner parameters, not pole masses. No derived map places emission and target in one observable convention | not evaluable | not evaluable | derive the readout observable and a complete scheme map, choose a current target in that convention, freeze the comparison, then evaluate it (SLP-06) |
| CL-6 | Printed pair identity: α_root = (P_fwd − φ)/√π | exact | holds to ≥ 35 digits after converged precision-100 reruns of `../code/P_derivation/runtime/full_p_alpha_report_current.json` and `../code/P_derivation/runtime/p_closure_trunk_current.json`; the converged α⁻¹ = 136.994835177413 supersedes the earlier printed tail 136.994835164622 beyond digit 9 | 3.1×10⁻³⁶ (report), 1.2×10⁻³⁸ (trunk) | n/a | closed by converged rerun; CI test `../code/P_derivation/test_printed_pair_identity.py` (A-03) |
| CL-7 | Capacity readback map: F(N) = N | exact | conditionally constructed: G2-GAP-1 proves the readback fixed point equals the bridge capacity modulo premises CP-1 (balance), CP-2 (inversion form), CP-3 (averaging carrier); certified conditional fixed point 3.5321315434×10¹²² at relative width 1.6×10⁻²⁵; unconditional construction open | n/a | n/a | prove the port-load/D10 coupling (G2-GAP-1); declared-structure candidates constructed and excluded 2026-07-14 (`../code/capacity_readback/F_CONSTRUCTION_2026-07-14.md`) |

## Reading rules

- CL-1 and CL-2 consume the same open transport function but are distinct maps.
  Historical v2,
  `../falsification/frozen_targets/hadronic_closure_target_2026-07-14_v2.json`,
  is immutable provenance only: its pass rule compares the emitted total
  `Delta_source` with a residual and cannot govern promotion. It also gives one
  scalar target to both maps although CL-2 adds alpha_U(P). Correct point
  diagnostics are total 8.7280337037… for CL-1 and 8.6869093675… for CL-2;
  their difference is alpha_U(P) = 0.0411243362…. The residual diagnostics are
  0.0414658610… and 0.0003415248…, respectively.
- The v3 corrective contract is
  `../falsification/frozen_targets/hadronic_closure_target_2026-07-16_v3.json`.
  It separates map and coordinate types and requires a target-blind function or
  interval over P. It is a post-target-access erratum scaffold and can never be
  activated in place. A detached successor must bind a complete method frozen
  before withheld data or produced by an audited clean-room operator; v3
  cannot score V1 or any future payload.
- The 2026-07-16 V1 grid is exploratory, not blind and not promoted. Its source embeds
  canon target values in a compare-only block, its directing session had target access,
  and its min/max grid envelope is not a rigorous bracket. Its containment verdict used
  the correct total endpoint and S coordinates despite the residual/total naming bug.
  Its width is 1.17×10⁸ times v2's historical, invalid tolerance. V3
  registers no active tolerance or pass/fail rule.
- CL-3 supersedes all "two capacities" and "resonance at logarithmic depth" language:
  SL-4 admits one N. The displayed 6.6% is a conditional nominal mismatch, not an
  unconditional contradiction, until F and CP-1 to CP-3 close. A final comparison must
  carry the joint Planck posterior through the Λ-to-N map.
- CL-5's exact-match companion lane is calibration and does not appear in this ledger.
  The 96-cell one-loop sweep is a chart-coordinate diagnostic only. The 2026-07-16
  two-loop packet combines an MSSM one-loop baseline with SM two-loop increments and
  excludes only that hybrid prescription. The pole packet is an incomplete,
  scale-dependent partial PRTS/Feynman-gauge prescription with open tadpole/FJ and vev
  schemes. Its JKV check audits a slope term, not the full finite packet. Neither packet
  supplies a physical CL-5 comparison.
- A row closes only after its map and comparison object are complete. A frozen,
  target-blind computation and its contraction certificate are required where the
  protocol calls for them. A row never closes by relabeling.

## Protocol status per coordinate (basin-then-contract, SLP-10)

| Coordinate | Stage 1: basin | Stage 2: contraction certificate | Stage 3: landing | Stage 4: sharpened digits |
|---|---|---|---|---|
| P | located (CODATA α) | certificate for the two incomplete declared maps: `../code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json` (direct Banach on an alpha interval, mpmath.iv outward rounding, centered form, L ≤ 0.0724; both readout modes; declared map at declared cutoffs, edge-sum tails bounded); global at-most-one on the declared domain (α ∈ [0.005, 0.01], both modes, sup \|g′\| < 1 on all 256 subdivision pieces, worst L ≤ 0.3041, exceptional set empty): `../code/P_derivation/runtime/p_global_uniqueness_certificate_2026-07-14.json` | not evaluable until Ward-projected transport completes the map (CL-1/CL-2) | not armed |
| N | located (Planck Λ + bridge) | missing (F not constructed) | not evaluable | not armed |
