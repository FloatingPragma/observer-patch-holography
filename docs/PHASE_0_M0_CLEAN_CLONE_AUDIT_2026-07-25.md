# Phase 0 M0 Clean-Clone Audit

Date: 2026-07-25

## Acceptance boundary

Milestone M0 covers the registered publication artifacts, mandatory scientific
receipt suite, fail-closed claim and issue projections, external-data pins, and
the Phase 0 simulator producer inventory. A clean checkout can validate and
rebuild those surfaces without untracked scientific inputs.

M0 does not establish every quantitative or physical claim in the corpus.
Conditional branches, comparison rows, non-identifiability results, absent
physical producers, and later source-binding gates retain their declared
status. Issue #544 and all Phase 1 work are outside this audit.

No release was published. The visible release identifier remains `r1577`, and
no release or version metadata was changed.

## Pushed source states

- `reverse-engineering-reality`: implementation and exact clean-clone CI state
  `6166cc304fcd1aca86ea79262e8de2b24e2952ed`
- `oph-physics-sim`: simulator inventory and fail-closed producer state
  `a29a747cb13503cf930d8b1363d72faa2e00984e`
- The commit containing this audit carries the live issue-ledger snapshot and
  its generated claims scoreboard, proof-spine projection, and compression
  scorecard after all Phase 0 issue closures.

## Issue outcomes

| Issue | Result |
|---|---|
| #514 | Closed. Manifest membership and exact paper-ID-to-path binding are source-derived. Missing, extra, stray, cross-swapped, malformed, and stale publication artifacts fail closed. Credential rotation is waived by owner decision. |
| #512 | Closed. Registry, novelty matrix, falsification matrix, dependency graph, issue ledger, proof spine, compression scorecard, and claims scoreboard enforce exact coverage and generated synchronization. |
| #507 | Closed. One documented mandatory command is the CI command. Optional cloud and hardware lanes cannot enter mandatory science through broad discovery. |
| #518 | Closed with a no-go boundary. Exact same-antecedent countermodels establish non-identifiability for the two audited targets. Back-solved identities are non-promoting. |
| #517 | Closed. The synchronized proof packet covers selector, Byzantine safety, orphan-lock transitions, confluence, refinement moduli, product-pseudometric finiteness, and layer separation. |
| #542 | Closed. Every registered paper root and the book pass the warning, tracked-drift, and two-pass byte-identity gates. |
| #553 | Closed. Clean-checkout source and simulator campaigns pass with frozen provenance and no hidden scientific inputs. |
| #325 | Closed. The class-H evidence packet is synchronized across schema, verifier, fixtures, and documentation. Sufficiency remains relative to the declared threat model. |

## Source verification

The mandatory command is:

```bash
python3 tools/run_mandatory_suite.py
```

The local run on `6166cc30` produced:

- 25 source-bound canonical book diagram renderings validated
- 14 release-manifest tests passed
- 12 deterministic publication-gate tests passed
- 33 Phase 0 receipt and projection tests passed
- 1,281 mandatory tests collected
- 126 audit-fixture tests passed
- 9 A5 closure-ledger tests passed

The external-data registry validates nine artifacts, nine deterministic
loaders, ten upstream pins, and nine explicit `NOASSERTION` license records.
The Planck approximation artifact has SHA-256
`70357bb5100974d14b8c3291ccfcc2385f165a5181143c5c62acd27f812eb00f`.
Its loader writes explicit UTF-8 bytes, so Windows newline translation cannot
alter the artifact.

Clean-checkout continuous integration on `6166cc30`:

- [Mandatory Suite, Ubuntu and Windows](https://github.com/FloatingPragma/observer-patch-holography/actions/runs/30145476385)
- [Publication Build Gate, Ubuntu](https://github.com/FloatingPragma/observer-patch-holography/actions/runs/30145576220)

The publication gate validates source-bound book assets, builds every
registered paper and the book, enforces the warning budget, rejects tracked
artifact drift, and compares a second complete build byte for byte. The book
SHA-256 is
`5d1b8a05cebfa351ebf3333f881159b6852813bf4902d75fb944e927ff1c631a`.

The issue-517 proof-obligations receipt recomputes at SHA-256
`bd94e1b9dee69c9b61386357429996b0c523266d64ba55dfbd85f439939354aa`.
Ten focused tests pass, including platform-independent repository paths. The
independent adversarial review exercised 114 mandatory-field deletion probes
plus scenario and certificate mutations.

The class-H evidence packet has 56 passing focused tests. Its integrity-valid
reference fixture is explicitly nonphysical and returns `INSUFFICIENT`.

## Simulator verification

The generated producer inventory at `a29a747` classifies every audited Phase 0
row as produced, produced but non-promoting, declared input, retired
compatibility input, or unavailable producer. K1 uses a gauge-covariant
production mismatch and perturb-resettle replay. H2 hash-token features have
zero claim-bearing weight. Legacy Einstein sidecars are declaration-only and
cannot promote branch entry.

The isolated simulator clone, with no positive-geometry sibling checkout,
passes the exact 19-file Phase 0 producer and certificate suite:

```bash
python3 tools/build_producer_inventory.py --check
python3 -m pytest -q \
  tests/test_producer_inventory.py \
  tests/test_wzh_source_closure_backend.py \
  tests/test_paper_side_realized_branch.py \
  tests/test_theorem_contract.py \
  tests/test_bulk_proof_certificate.py \
  tests/test_einstein_tower_producer.py \
  tests/test_modular_normalization_producer.py \
  tests/test_gns_tower_producer.py \
  tests/test_event_manifold_producer.py \
  tests/test_stress_coupling_producer.py \
  tests/test_h3_worldline_stitch_producer.py \
  tests/test_subjective_observers.py \
  tests/test_neutral_bulk.py \
  tests/test_positive_geometry_kernel.py \
  tests/test_kernel_dispatcher.py \
  tests/test_run_bundle.py \
  tests/test_covariant_overlap.py \
  tests/test_modular_probe_sector_replay.py \
  tests/test_bw_array.py
```

```text
187 passed, 4 skipped
```

The four skips are declared external-plugin integration tests. Executed
fallback and dispatch tests emit `DECLARED_PLUGIN_UNAVAILABLE` and
`EXACT_GENERIC_FALLBACK`; they confer no acceleration, readout, or physical
prediction.

## Final ledger invariant

The live GitHub issue ledger is generated separately from the offline
scientific build. `tools/build_open_problem_ledger.py --check-live` matches the
final Phase 0 issue state. `tools/build_scoreboard.py --check` binds the same
snapshot into the claims scoreboard, proof spine, and compression scorecard.
The 33 focused ledger and projection tests pass.
