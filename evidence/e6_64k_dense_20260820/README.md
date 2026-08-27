# OPH-FPE 65,536-patch finite-run archive

This directory is a curated, hash-bound archive of the OPH-FPE run
`e6_64k_dense_20260820`. The simulator instantiates a bounded data path meant
to exercise observer-like self-reading: each finite patch row has local state
and twelve ports, neighboring ports exchange gauge-transported values, local
readback drives repair moves, and settled rows commit records. It then
constructs finite observer views from those records and evaluates diagnostic
and promotion gates. The dedicated source-to-observer contract marks the
observer-like self-reading receipt false because the required read-after-write
and feedback-loop clauses do not close. A historical aggregate in
`AUTO_THEOREM_UNIVERSE_SUMMARY.json` nevertheless reports the same-named
receipt as true after reducing it to the presence of observer rows. That
semantic collision is a false positive; the dedicated contract controls, and
the archive verifier requires the disagreement to remain visible.

The source code is the public
[`oph-physics-sim`](https://github.com/muellerberndt/oph-physics-sim)
repository at commit
[`b52196b296435d704b14d005d1f69caaaa662f97`](https://github.com/muellerberndt/oph-physics-sim/tree/b52196b296435d704b14d005d1f69caaaa662f97).
`config.yml`, the seed material, source arrays, reports, traces, observer rows,
and compact numerical arrays retained here are listed byte-for-byte in
`archive_manifest.json`.

## What this run reports

The configured finite regulator contains 65,536 patches, twelve local ports per
patch, 390,924 routed seams, and 128 cycles. It materializes 2,048 patch-observer
rows with 96-patch neighborhoods and four cap-observer rows. The separate
theorem-replay source state contains 326,047 covariant mismatches. In the
driven trace, the cycle-0 update begins from 326,053, leaves 324,739, reaches
178,280 after cycle 50, and first reaches zero at cycle 96. The
committed-record count is 28,888 at cycle 100, first
reaches 65,536 at cycle 107, and remains there at cycle 110. The driven trace
has 69 cross-cycle mismatch injections, so its global Lyapunov-descent receipt
is false.

The separate replay report records 326,047 accepted strict-descent moves and no
strict-descent violation. It checks 512 edge-slot diamonds, including 256
shared-node pairs, 256 disjoint pairs, and four frame relabelings without a
reported violation. Those checks do not establish finite consensus. The exact
endpoint-branch argument provides a structural nonconfluence witness and
certifies at least two distinct terminal quotient orbits. Its exact-check hash
list is empty; only one sampled terminal hash was materialized before the
schedule replay stopped fail-closed. `FINITE_CONSENSUS_THEOREM_RECEIPT` is
therefore false.

The global point set and adjacency are explicitly labeled
`legacy_fibonacci_knn_control`. They are a numerical control surface rather
than an earned icosahedral refinement tower or a carrier-to-support
realization. The twelve-port patch state is instantiated and summarized, but
the compact output profile did not write the full patch-state artifact:
`echosahedral_patch_state_report.json` records `artifact.written = false` and
`artifact.reason = compact_output_profile`. The finite replay source state is
present separately as `finite_consensus_source_state.npz`.

The run does not pass the finite-consensus, endogenous modular-clock, strict 3D
bulk, 3+1-dimensional event-manifold, physical CMB, dynamic dark-transport, or
cosmological-perturbation gates. Visualization-only assumptions in
`simulation_assumption_manifest.json` do not change those failed gates. The
generated `README_OPH_UNIVERSE_PACK.md` is preserved as source output; final
claim status must be read from the fail-closed replay, emergence, geometry, and
physical-promotion reports. Its positive observer-like line is affected by the
same aggregate naming defect described above.

## Custody boundary

The original local run tree has 281 files and 1,300,918,262 bytes. An audit
computed the sorted relative-path inventory digest
`8c71ec07c6f549ffa5c9cd0a7de732c0097f9e0cb0bfdfe044cec9fd8d9541f1`.
That digest is informational because the full tree and exact inventory-line
serialization are not archived here. This 45-file selection contains
66,309,382 bytes and has its own verifier-enforced inventory digest.

The 591,815,800-byte `observer_views.jsonl` is retained losslessly as
`observer_views.jsonl.zst`. Its uncompressed SHA-256 is
`3a38b6dfbd480898c58b66b46dba34af500701c047f0052feb679a319a3ca7e0`;
it contains 2,052 rows (2,048 patch observers and four cap observers).

Two provenance limitations remain explicit. The simulator checkout was clean,
while the run manifest records the contemporaneous
`reverse-engineering-reality` checkout as dirty at commit
`cc3abe7b7c3abea294735cd14328a9529b7ad76c`. The dirty-worktree digest is
retained, but no changed-file inventory accompanies it. No stdout/stderr run
log survives. The reports, configuration, seeds, arrays, and row data are
therefore the available primary custody record.

## Verify

From this directory, with Python 3 and the `zstd` executable installed:

```bash
python3 verify_archive.py
```

The verifier checks the exact inventory, every archived byte count and SHA-256,
the compressed observer file and its uncompressed stream, the main finite-run
counts, the trace checkpoints, provenance, internal consistency of the
reported structural witness and terminal-orbit fields, and the false
physical-promotion gates. It is a custody and report-consistency verifier; it
does not independently rerun the repair kernel or derive the incidence
witness from the archived arrays. It exits nonzero if a check fails or if
`zstd` is unavailable.

## Reproduce the simulator run

In a clean clone of the simulator at the pinned commit, install the declared
package and run the archived configuration:

```bash
git checkout b52196b296435d704b14d005d1f69caaaa662f97
python3 -m pip install -e '.[dev]'
oph-fpe run-oph-universe \
  --config configs/e6_axiom_manifest_64k_dense_observers.yml \
  --out-dir runs \
  --run-id e6_64k_dense_20260820 \
  --max-screen-points 65536 \
  --max-observers 2048
```

The seed (`20260805`) is pinned inside the configuration. This is the
reconstructed source-level invocation supported by the retained configuration
and CLI, since the original shell command and stdout log were not retained.
Reproduction should be judged by the archived hashes and report fields, with
platform and dependency differences disclosed rather than silently accepted.
