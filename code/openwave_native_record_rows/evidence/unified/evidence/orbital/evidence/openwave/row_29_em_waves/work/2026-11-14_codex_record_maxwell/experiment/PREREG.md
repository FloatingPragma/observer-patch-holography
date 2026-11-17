# Preregistration before the first execution

The machine timestamp and hashes are in `PREREG.json`. The root coordinator
must acknowledge the shared ledger entry before `probe.py` runs. No output has
been inspected when these thresholds are written.

Use exactly the candidate, source, compression and resources in `CANDIDATE.md`.
Run the frozen `probe.py` once, retain any failure and every raw branch/tape.
Then repeat once and require byte-identical receipt and compressed records.
No fitted parameter, discarded seed or selected successful branch is allowed.

Fixed checks:

1. Graph size is 216 cells, 2,592 vertices and 16,848 unique seams.
2. All 101,088 three-slot native branches conserve the unit and agree with
   current-only replay; 384 full-vector checks agree with that compression.
3. All twelve port-conditioned macro transition counts agree exactly, are
   symmetric, and have row count 2W. Their Laplacian eigensystem residual is
   below 1e-10. No target Green matrix is used to construct the transfer.
4. Exactly 65,536 independently scheduled native events are retained across
   four disclosed preparations. The reader rebuilds all occupation counts
   solely from the current tapes and declared initial positions. Its six
   antisymmetric readouts are retained without a physical accuracy claim.
5. Evaluate the finite-block population covariance at n=10^7 and 10^9 by the
   exact geometric-sum identity, numerically diagonalizing only the
   branch-derived 216-state macro transfer. Record both. For n=10^9 and each
   spatial index (1,0,0), (1,1,0), (1,1,1), require: separated-slice covariance
   at lags1,2,3 agrees with the periodic native Green target within1%; the
   slice recurrence residual is below1%; the fitted energy agrees within1%;
   the magnetic covariance has one eigenvalue of magnitude below1e-10 and
   two positive eigenvalues; reflected electric eigenvalues agree with them
   within1% of the larger magnetic eigenvalue.
6. The unchanged decoder on IID-position resampling has separated covariance
   below1e-12, while native lag2 covariance is above1e-5 at each momentum.

These gates concern a finite conditional population calculation and record
replay. They do not claim a measured quantum phase, finite-block Gaussianity,
an operational physical clock, a source intervention experiment, a full
OpenWave criterion or review by another group. No row-level gate is passed.

Resources: local Python3 plus installed NumPy; no new package; one process;
target CPU below30 seconds, stop/reconsider if above120 seconds. Fixed seed
2026111429. Floating population values are reported as numerical evaluations
of the exact finite transfer formula, not exact rational probabilities.
