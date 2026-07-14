# Addendum 01 to hadronic_closure_target_2026-07-14.json

Date: 2026-07-14 (same day, after the interval contraction certificate landed).
The frozen target file is not modified by this addendum; its sha256 remains
`7cedad0a7281c74ca0fb1105120c991aeab2f3c45bf86adbbfd560c6324fb985` and its
OpenTimestamps proofs remain valid.

## What this addendum records

1. The interval contraction certificate
   `reverse-engineering-reality/code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json`
   established that the previously printed source root tail
   `136.994835164621649…` came from an unconverged solver run and is
   superseded beyond digit 9 by the certified converged root
   `136.994835177413…` (ledger row CL-6, closed). Values inside the frozen
   target that were transcribed from the pre-certificate chain inherit that
   stale tail at and beyond the corresponding digit position.
2. Magnitude assessment: the supersession is a shift of ≈1.3×10⁻⁸ in
   inverse-alpha units. The quantities the target tests are of order
   4×10⁻² inverse-alpha units, and the pass criterion propagates the CODATA
   alpha uncertainty (2.1×10⁻⁸ on alpha_inv). The stale-tail shift is
   therefore below the leading test scale by six orders of magnitude and
   comparable to the measurement-uncertainty floor.
3. Evaluation rule: the one-shot public comparison declared in the target
   proceeds against a v2 target recomputed from the certified converged
   chain, frozen and externally timestamped BEFORE any payload work begins.
   The v1 target stays on record untouched; v2 supersedes it for evaluation.
   If payload work has begun before a v2 exists, v1 governs and this
   addendum stands as the recorded precision caveat.
4. The blindness requirement is unchanged: the payload computation must not
   read v1, v2, this addendum's numeric content, CODATA alpha, the SL-3
   pixel estimate, or anything downstream (dependency-cone audit V-08).

## Status

- v2 target: PENDING — recompute the target quantities from the certified
  chain at certificate-grade precision, freeze, timestamp.
- This addendum is itself timestamped alongside its commit.
