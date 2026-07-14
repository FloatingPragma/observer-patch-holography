# OPH Closure Ledger

The equations the strange-loop principle set ([STRANGE_LOOP_PRINCIPLES.md](STRANGE_LOOP_PRINCIPLES.md)) requires to hold exactly, with
current residuals. This table is the primary evidential surface of the program: under
SL-0, a closure residual is a measurement of the hypothesis. Rows are permanent; verdicts
append. Last recomputed: 2026-07-14 (manual; CI regeneration is fix item SLP-02).

| # | Closure | Required | Current status | Residual (relative) | Residual (measurement σ) | Resolution path |
|---|---|---|---|---|---|---|
| CL-1 | P-loop, source chain: P = φ + √π/A_T(P) with source-only transport | exact | forward map contracts to α⁻¹ = 136.994835177413 (interval-certified unique fixed point, width 7.2×10⁻²⁴), outside the SL-3 basin | 3.0×10⁻⁴ | ≈2.0×10⁶ | Ward-projected hadronic transport, blind, vs frozen target (SLP-03) |
| CL-2 | P-loop with gauge-width term, self-consistent map | exact | certified unique fixed point α⁻¹ = 137.035660136946577 vs 137.035999177(21). Verdict appended 2026-07-14: the previously displayed 137.0359595 is a mixed-provenance packet (inner value at P_fwd plus α_U at the SL-3 pixel), not a fixed point of any single declared map; it stays on record as a display packet only | 2.5×10⁻⁶ | ≈1.6×10⁴ | same missing hadronic term |
| CL-3 | One capacity (SL-4): N_EW = π·exp(6π/(P·α_U)) vs N from Λ | exact | 3.53×10¹²² vs 3.31×10¹²² | 6.8×10⁻² | — (dominated by model, not Λ error) | derive the EW-bridge/capacity connection or carry as open contradiction (SLP-04) |
| CL-4 | Hierarchy bridge at the one N: α_U·log(N/π) = 6π/P | exact | 11.5546 vs 11.5573 at the Λ-estimated N | 2.3×10⁻⁴ (log scale) | — | same object as CL-3 |
| CL-5 | Forward electroweak emission at SL-3 pixel: (M_W, M_Z) | exact | (80.330, 91.119) vs (80.3692(133), 91.1876(21)) | 4.9×10⁻⁴ / 7.5×10⁻⁴ | 2.9σ / ≈33σ | repair forward chain from axioms with preregistered revisions (SLP-06); suspects: frozen MSSM coefficients, QT1–QT5, β_EW |
| CL-6 | Printed pair identity: α_root = (P_fwd − φ)/√π | exact | holds to ≥ 35 digits after converged precision-100 reruns of `runtime/full_p_alpha_report_current.json` and `runtime/p_closure_trunk_current.json`; the converged α⁻¹ = 136.994835177413 supersedes the earlier printed tail 136.994835164622 beyond digit 9 | 3.1×10⁻³⁶ (report), 1.2×10⁻³⁸ (trunk) | — | closed by converged rerun; CI test `code/P_derivation/test_printed_pair_identity.py` (A-03) |
| CL-7 | Capacity readback map: F(N) = N | exact | F not yet constructed | — | — | construct and execute F; stage-2 certificate (C-01, SLP-10) |

## Reading rules

- CL-1 and CL-2 are the same open term seen from two points on the chain. The completed
  hadronic transport must supply Δ_source = 0.0414658… α⁻¹ units at the endpoint
  (equivalently S_required = 0.8954001…); the target is frozen in
  `../falsification/frozen_targets/hadronic_closure_target_2026-07-14.json` and the
  computation must not read it.
- CL-3 supersedes all "two capacities" and "resonance at logarithmic depth" language:
  SL-4 admits one N.
- CL-5's exact-match companion lane (frozen measured W/Z pair) is calibration and does not
  appear in this ledger.
- A row closes only by a blind computation landing inside the stated basin, with its
  stage-2 contraction certificate. A row never closes by relabeling.

## Protocol status per coordinate (basin-then-contract, SLP-10)

| Coordinate | Stage 1: basin | Stage 2: contraction certificate | Stage 3: landing | Stage 4: sharpened digits |
|---|---|---|---|---|
| P | located (CODATA α) | certificate: `code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json` (direct Banach on an alpha interval, mpmath.iv outward rounding, centered form, L ≤ 0.0724; both readout modes; declared map at declared cutoffs, edge-sum tails bounded) | outside basin (CL-1) pending hadronic term | not armed |
| N | located (Planck Λ + bridge) | missing (F not constructed) | not evaluable | not armed |
