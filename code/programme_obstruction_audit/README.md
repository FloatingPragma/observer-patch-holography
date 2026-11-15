# Finite source and attachment controls

This package reproduces finite mathematical controls for the fixed source
formulation identified in `provenance/fixed_sources.json`. It checks scoped
attachment and non-entailment cases; the supplied analytic scalar-seam source
is an abstract witness under its declared typed interpretation. These
controls do not determine a physical model, establish an empirical exclusion,
or replace the analytic proof.

## Run the seven offline controls

From this directory, with an ordinary Python interpreter and assertions enabled:

```sh
python3 run_all.py
```

Only the Python standard library and packaged files are required. No network or external scientific input files are read. The entry point uses `sys.executable` for all seven subprocesses, validates package SHA-256 hashes, checks the immutable response fixture against its hardcoded SHA-256 and Git blob SHA-1, validates the fixed target/source provenance, and stops at the first failed assertion, process failure or 120-second per-script timeout. It prints a concise verdict and the path to a full `summary.json`, saved alongside all individual JSON receipts in a fresh temporary directory. Use `--output-dir /path/to/empty-directory` to choose a destination outside this package; package destinations and existing nonempty destinations are rejected before any control executes. A failed run therefore cannot retain a previous run's success summary. Do not use `python -O`; these scientific controls rely on assertions. The publication run used **Python 3.14.7**, exact rational arithmetic and explicit Q(sqrt5) pairs; other supported Python versions can differ only in runtime metadata and their reported script/environment fields, not the exact mathematical results.

| Script | Prescribed inputs and controls | Meaning and stopping boundary |
| --- | --- | --- |
| `channel_author.py` | Rational 8x8 binary seam channel, four Kraus operators, private-controlled and nonnormal-exchange rotations; all 64 matrix units; completion n=2,3; two private inputs | Exact no-feedback control, CPTP/completion and feedback difference 144/625. Toy channel class only; no OPH source membership. Stop after these cases. |
| `channel_independent.py` | Independent private-first tensor order, weighted Kraus reconstruction, all 64 matrix units, completion projector identities and 64x4 private commutator equations | Independent exact feedback/restricted-class validation; scalar-only private commutant excludes every single-product exponential without logarithm-branch assumptions. |
| `joint_author.py` | Fixed twelve-vector Q(sqrt5) response fixture, at most60 stored paths,66 binary configurations, four oscillator cutoffs and three clock settings; s=2 units and idle controls | Corrected24 finite gates for carrier/triplet/CCR/clock attachments. Does not certify full source membership. |
| `joint_independent.py` | Independent D formula/rank, all 20 Cayley factors and60 path words,3960 event branches,431 A3 omission-coordinate controls, posterior and three current-unit controls | 907 prescribed gates; overlaps author/event checks and is not an extra907 scientific findings. Separate reconstruction imports no author implementation. |
| `kinetic_control.py` | A=I+path3 stiffness; e=0,+1/16,-1/16,-1/4; gamma=1,3/2 | Eight rational cases. Singular -1/4 is rejected. Same q-Gibbs marginal is an analytic integration result; full reference/vacuum changes and same-capacity inertia premise fails. |
| `attempt_independent.py` | Incidence reconstructed separately;66 states x30 seams x2 coins; s=2, idle half and unequal edge-weight controls | Exact named-process components B,nu,kappa and mean law. The kappa-weighted ratio needs unit atomic current; 3B/nu is current-unit invariant. No full-core independence assertion. |
| `compact_classifier.py` | Exact partitions of11,12 using imported simple dimensions3,8,10 | Unmodified upstream tiny classifier; finite arithmetic only, not reproof of compact-simple classification or physical gauge realization. |

The matrix-unit controls cover arbitrary complex operators by complex linearity despite real rational basis entries. Conjugated completion at n=2,3 is a regression check; the all-n identity follows analytically by unitary conjugation of the native decomposition. The unit/scheduler controls deliberately change their declared convention/law and do not falsify the frozen relation. No stochastic search, fitting or absence-of-witness inference occurs.

Individual commands are also portable:

```sh
python3 channel_author.py > author-result.json
python3 channel_independent.py > independent-result.json
python3 joint_author.py data/response.json > joint-result.json
python3 joint_independent.py > independent-joint-result.json
python3 kinetic_control.py > kinetic-result.json
python3 attempt_independent.py > attempt-result.json
python3 compact_classifier.py > classifier-result.json
```

`receipts/` contains the final control execution and two retained upstream
replays. `provenance/package_sha256.json` covers every package file except
itself. `provenance/packaging_adaptations.json` preserves original script
hashes and describes the adaptations: the author joint receipt uses a stable
fixture label plus actual-input SHA256; the independent joint fixture path is
package relative and its retained-evidence sidecar write is removed; and the
kinetic control writes to stdout. Arithmetic, cases and scientific scope are
unchanged. Frozen receipts retain their original producer hashes and
environment metadata; current replay receipts describe the current producer.

Run `python3 -B -W error test_run_all.py` for the package-preservation regressions. They copy the package to a temporary directory, run all seven controls twice with explicit and default fresh output destinations, verify that package file names, bytes and modification times remain unchanged, and check that unsafe output destinations are refused without changing existing files.

## Reproduce the separately scoped upstream components

The retained `upstream_components.json` and `upstream_galois.json` were executed on fixed source bytes with **Python 3.14.7, SymPy 1.14.0, mpmath 1.3.0**. They are not counts to add to the offline reconstructions. The seven components equal the upstream retained receipt: response, algebra, constraint grammar, central channels, four-level support, records and refinement. Both structural Galois laws pass response/observable rank12, restricted rank11, even/odd rank6,66 brackets,20 factors,60 actions,720 covariance pairs and3600 compositions. Whole-law conjugation includes the quadrupole and fixes i; it does not map arbitrary PSD preparations/effects.

To rerun those components, provide a **separate exact fixed-source fixture or checkout**, and an environment with those pinned dependencies. This optional step has explicit external-source/dependency requirements; the offline command above does not. The helper validates the source paths/hashes listed in `fixed_sources.json` and `additional_source_pins.json`, then uses the same upstream component calls and canonical JSON comparison as the original wrapper. It does not fetch anything or alter the source fixture. From a repository root, where `FIXED_SOURCE` names that separately prepared fixture:

```sh
python3 code/programme_obstruction_audit/replay_upstream_components.py --source-root "$FIXED_SOURCE" > upstream-component-replay.json
python3 "$FIXED_SOURCE/code/a5_closure/galois_source_response_control.py" > upstream-galois-replay.json
```

The exact Galois script has Git blob SHA-1 `ecb5f9af4b9ab9d3ccde131e7228e174d4385c11`, SHA-256 `71cfa71dce4080febc61e09b3675f8fc8822294acdbacfbd1f87012e9e0e618d`. It validates its four source pins internally. The packaged response snapshot has Git blob SHA-1 `e56b96a9429ac618aea9de7cd5ffa53a13354ff6`, SHA-256 `226ab001894f357d3732d8a898fbae80e28b9e16ca41104a6f293190b441e053`. Replaying a newer checkout by assumption is invalid; use the frozen hashes. The packaged optional component helper was also executed successfully on that fixed fixture; its receipt is `receipts/upstream_components_portable.json`. No full upstream package or installed dependency tree is vendored here.

These replays **did not run** top-level recursive source custody, producer regeneration, the full hostile-input suite, N4's separate 187-check numerical packet, scalar closure rebuilding or Lean compilation. The optional helper verifies the finite fixture paths it actually consumes; that is not a claim that the omitted top-level upstream verifier ran. All-level tower/operational model membership rests on the separately challenged analytic proof, not four finite support levels or a passing receipt. Historical source receipt bytes remain untouched. Runner integration controls additionally reject a modified fixture even when its manifest is edited consistently, and reject optimized Python before executing the scientific controls; `receipts/runner_refusal_controls.json` records those justified refusals.

## Independent review and evidence boundary

`provenance/review_status.json` records methods, reviewers, note hashes and
accepted corrections. Author replays establish reproducibility. Independent
proofs use a separate implementation, tensor order or algebraic
reconstruction, and source-clause review challenges admissibility and
quantifiers. The scalar-seam source is a claimed full **abstract** model under
its typed interpretation. Finite controls do not establish physical source
membership or empirical selection.

The fixed source revision and verified Git objects are recorded in
`provenance/fixed_sources.json` and `provenance/local_object_verification.json`.
The generated active-surface inventory is maintained by the repository's
canonical scanner. These controls do not change the adopted basis, claim
statuses, premise map or frozen evidence. They contain no measured
observations, laboratory preparations or resource lower bounds.
