# Frozen stage2 task

Inputs: original stage1 `proposal.json`, `source_columns.npz`, all source pins,
the local code and `STAGE2_DESIGN.md` bound in `PREREG_FLOW.json`. Deterministic
complete finite calculations; no seeds, data fitting or measured constants.

1. Reconstruct both full240 source kernels from the pinned literal endpoint
   table. Preserve the held seam(0,5), original1260-context denominator and
   all columns. This recheck uses8 literal endpoint calls and enumerates
   302400 source columns per law; it is not302400 new literal calls.
2. Exact rational Gram and residual infinity/one-norm bounds must each be
   at most10^-9. Every proposed eigenvalue is between the derived104/105
   floor and1. No dropped modes, repaired proposal or eigensolver trust.
3. Compile all240 rank-one phase words for native times1,4 and held-seam
   time4. Every emitted word has at most65536 pulses. Full M6 primitive
   Frobenius error is at most10^-11, pulse allowance2e-11, split component
   allowance2e-12, scalar log-phase enclosure error at most10^-12.
   Every complete flow half-diamond bound, including padding aborts, must be
   below1/5000. Stored-map gain is retained rather than normalized away.
4. Execute the source swap-parity reference instrument under the design's
   exact all-column/selected-full-scratch policy. Exact plus/minus masses are
   239/240 and1/240 after successful selectors; valid plus referenceT/239.
   Reference complete half-diamond bound must be below10^-4. Selected complete
   parity laws must agree within10^-8; dirty scratch probability below10^-16.
   All selector and reset timeouts remain explicit.
5. Native coordinate and original energy-ray pointer records for both C/E
   terminal readers must have primary outside-cell mass below10^-7 and
   f-ray click above1-10^-7 at times1,4. The same original detector must show
   held-seam outside-cell mass above10^-3 at time4. A physical source-X escape
   plant must register above1-10^-7. These numerical arrays are diagnostics;
   certified physical-instrument comparisons use the full error budgets.
6. Two-setting source/reference interference at native time4 must have
   contrast above1-10^-7. Retaining a source bit4 which-cell read must give
   absolute contrast below10^-7. Source-read/reset-to-f with the original
   mixed reference must also give absolute contrast below10^-7, despite its
   preserved f-ray click. These are complete ideal-control laws with the
   inherited native reset/pointer error added, not discarded records.
7. Save every emitted word's deterministic hash, all primitive matrices and
   exact enclosures, source/processor/padding/reference arrays, all512 final
   outcomes per interference setting, raw pointer laws and elapsed runtime.
   Failure creates a retained attempt with its exception; do not overwrite it.
8. Independent headline recomputation must reproduce the source/eigenray,
   held-law change and reported native detector quantities in under60seconds.
   It does not assert full million-pulse reproduction. A same-group reviewer
   must challenge the certificate and scientific meaning before sealing.

Command from this directory:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 run_flow.py --record --attempt stage2_001
```

All identities are separately labelled: positivity, normalization, eigenray
stationarity under its declared spectral calculus and reference algebra do
not establish emergent physical time, dynamic selection or rest mass.
