# Reproducing the mandatory scientific suite

This is the clean-clone path into the OPH scientific receipt suite. The first
commands verify that the claim registry is internally connected and that the
public test collection imports without error. Individual evidence families
then have their own theorem, certificate, or experimental acceptance rule.

## Runnable observer and repair-law model

The live [model implementation](model/README.md) has a separate clean-clone
reproduction gate. From the repository root:

```sh
python -m pip install -r model/requirements.txt
python -m model.verify
python -m pytest -q -W error model
```

The verifier checks the artifact manifests and replays every receipt declared
in `model/verify.py`. Its replay inventory is narrower than the set of all
retained results: native-link selection and holonomy reconstruction, for
example, are replayed by their model regression tests. Run both commands;
neither artifact hashing nor one replay command replaces the other tests.
The standalone source-field and particle checks used by CI are listed below.

The finite observer contributes the exact native observer, its
quantum mixing and record capacity, a 20,000-trial sequential run, complete
record recovery, the P-dependent binary detector with 200,000 trials and
the observer-local incident-tick calculation. The field detectors
contribute the scalar detector with native feedback, its retained two-read
continuation, its local evolution and the source-compiled two-site
apparatus. The prediction package and the source-executed C/E experiment
contribute the instrument comparison. The common interacting source
contributes its finite execution, the infinite-occupation history, the
complete sharp-output comparison, the observer readout, smooth-detector and
thermal-clock bridges, the causal cones and the native field relations.
Its corrected native-memory history adds a source encoder and recovery,
complete logical/reference comparisons, and the matched passive failure control.
The [mixed-state observer](model/common_source/OPEN_HISTORY.md) additionally
executes the finite source environment through the four-read continuation,
checks complete quantum outputs and fresh-draw controls, and attaches the
same corrected/passive native memory with separate error and resource budgets.
The [smooth-history certificate](model/common_source/SMOOTH_HISTORY.md)
compares all 16 retained work outputs with the infinite history for each of
the three native-current preparations, including later feedback and the
gauge clock probe. Its state-specific bound is separate from the broader
uniform occupation theorem and from numerical or physical service errors.
The native coherent-continuation and source-selection receipts retain their
distinct phase-sensitive observation contracts.
The spatial Bose and canonical scalar branches add complete interacting
histories and approximation certificates. The shared scalar observations
add spatial and occupation comparisons, the retrospective published-data
readout check and the separate source-executed calibrated readout.
The [observer-record reconstruction](model/observer_reconstruction/README.md)
adds executed free-field and clock records, record-only spatial/spectral
identification, withheld detector predictions and an exact native response
geometry. Its normal-mode occupation bound and noncommuting-phase control
retain the supplied-action and finite-identification scope.
The [binary geometry](model/relative_observer_geometry/README.md) checks a
current-only calibration, quantum occupation response and finite propagation
on the clocks' native sector. The [host-clock comparison](model/native_clock/README.md)
replays actual recorded native calls and recalculates their held-out timing
tests; reproduction checks the retained measurements rather than requiring
another machine to reproduce their timestamps.
The [joint clock bridge](model/native_clock/CLOCK_BRIDGE.md) retains zero-hit
and movement-count correlations, propagates exact finite count uncertainty,
and preserves the recorded workload-control rejection. The
[field refinement](model/observer_time_fields/REFINEMENT.md) also evaluates
canonical fixed-number preparations with their conservation correlations.
[Temporal records](model/relative_observer_geometry/README.md#temporal-hilbert-space-reconstructed-from-field-records)
reconstruct an exact elementary cyclic transfer and a separate, bounded
two-vector compression on the actual screens. Their preparation and clock
domains remain distinct; these checks do not establish physical time.
The [Bell-record replay](model/predictions/BELL_RECORDS.md) reproduces a public
fixed-count comparison against local factorized classical responses. It is
a retrospective boundary test, not a prospective OPH validation.
The tests independently check endpoint arithmetic,
all retained histories, classical-record preservation, quantum coherence,
controls, capacity and
the stopping-time null. The exact native observer and P-response runtime
use the Python standard library; the complete test and local-tick route
uses the pinned model requirements. Individual commands and scope are in
[observer/](model/observer/README.md),
[observer/p_response/](model/observer/p_response/README.md),
[observer/scalar_response/](model/observer/scalar_response/README.md), and
[local_tick/](model/local_tick/README.md).

The complete smooth-history comparison can also be replayed directly:

```sh
python -m model.common_source.smooth_history --check model/common_source/smooth_history_results.json
```

Focused field, temporal, clock-conversion and public-data replays are:

```sh
python -m model.observer_time_fields.canonical_sector --check model/observer_time_fields/canonical_sector_results.json
python -m model.relative_observer_geometry.temporal_records --check model/relative_observer_geometry/temporal_results.json
python -m model.native_clock.clock_bridge --check model/native_clock/clock_bridge_results.json
python -m model.predictions.bell_records --check model/predictions/bell_results.json
```

The scalar adapter executes the supplied field instrument and its conditional
branches, with the original-action vacuum and native record-only continuation.
Its retained [prediction protocol](model/observer/scalar_response/PREDICTION.md)
compares the feedback statistic with null and reversal controls and transfers
the existing scalar refinement error to the complete history law. These are
conditional software tests; the scalar action, Born readout and model-time
convention are supplied, and independent physical comparison is separate.

The [Runnable observer model workflow](.github/workflows/observer-model.yml)
runs the verifier and model regressions on Linux for changes to the model and
its source dependencies, and on Windows nightly, on release tags and on manual
dispatch. It also executes the following
standalone source-field and particle checks in the model environment:

```sh
python code/source_scalar_time_refinement/verify_source_scalar_time_refinement.py
python -m model.common_source.repair_field_attachment --check model/common_source/repair_field_attachment_results.json
python -m model.common_source.native_staggered_field --check model/common_source/native_staggered_field_results.json
python -m model.common_source.native_field_history --check model/common_source/native_field_history_results.json
python -m model.common_source.native_charged_spectrum --check model/common_source/native_charged_spectrum_results.json
python -m model.common_source.native_code_transport --check model/common_source/native_code_transport_results.json
python -m model.common_source.native_charged_spectroscopy --check model/common_source/native_charged_spectroscopy_results.json
python -m model.native_equilibration.wedge_modular --check model/native_equilibration/wedge_modular_results.json
python -B model/observer_particles/run.py --verify-only
python -B model/observer_particles/run.py
```

Together these commands reproduce the workflow's separate gates; running
only `model.verify` does not. They supplement the source-frozen mandatory runner below.
Neither that runner nor a bare `python -m pytest` includes the new `model/`
tests: invoke the model commands explicitly when changing this implementation.
The finite observer construction is executable with declared apparatus and
accepted closure inputs; no hardware fabrication, calibrated physical clock
or whole-universe implementation is inferred from a passing run.

## Additional component calculations

The bounded modal Maxwell-shaped factorization is also an exact, laptop-scale
algebra check rather than a simulation campaign:

```bash
python -m pytest -q \
  code/electromagnetism/test_modal_maxwell_factorization.py \
  tools/test_modal_maxwell_factorization_surfaces.py
```

The coupled-action controls can also be reproduced together:

```bash
python -m pytest -q \
  code/electromagnetism/test_whitney_charged_enclosure.py \
  code/electromagnetism/test_whitney_magnetic_continuum.py \
  code/electromagnetism/test_neutral_packet_observables.py \
  code/electromagnetism/test_neutral_packet_projection.py \
  code/electromagnetism/test_whitney_quantum_packet.py
```

These check a rigorous fixed-mesh time enclosure, numerical controls for the
analytic prescribed-magnetic-field continuum theorem, full-dimensional
neutral-state preparation, and the closed-form neutral projection. Their separate acceptance rules do not certify
physical calibration or useful interacting quantum propagation.

The neutral-packet controls preserve the projection probability and the
scalar nodal squared radius through exact-rational invariants and a budgeted
high-precision evaluation.
The pointwise projection evaluates the closed form of the circle integral,
and the legacy `nodes` argument is validated and has no effect. The
[Whitney numerical methods](code/electromagnetism/WHITNEY_NUMERICS.md)
gives the evaluator, its numerical contract, the independent replay and the
receipt regeneration commands. The dedicated neutral-packet workflow
executes these controls on Linux and Windows.

## Native observer dynamics and retained large-run statistics

The compact [observer-dynamics evidence](evidence/observer_dynamics_20260925/README.md)
contains the pinned simulator source needed for standalone reproduction:

```sh
python3 code/observer_dynamics/verify.py
python3 code/observer_dynamics/verify.py --replay-refinement
python3 evidence/observer_dynamics_20260925/lean/verify.py --check
```

The first command checks immutable inputs, recomputes bounded source and
clock certificates, and runs the numerical and mutation controls. The second
also reconstructs all 61 refinement fresh-noise streams. The third compiles
the corresponding Lean modules and checks their axiom receipt. The mandatory
suite includes the first command. Dependencies and exact interpretation
boundaries are in [the executable interface](code/observer_dynamics/README.md).
The large terminal arrays are represented by their retained statistics and
cross-pins; this command does not replay their microscopic repair histories.
The conditional CMB comparison retains its supplied source and background.

## Finite source, preparation and common-history interfaces

These four packages check bounded observer-like systems with local state,
read ports, retained records and declared control operations:

```sh
python3 code/native_geometric_source/verify.py
python3 code/common_history_packet/verify.py
python3 code/source_scalar_preparation/verify_preparation.py
python3 code/source_density_geometry/verify.py
```

They respectively test fixed-geometry source identifiability, a classical
readout/restoration history with a useful complete error bound, finite local
quantum-state preparation on the supplied scalar action, and relative action
measure with calibrated density controls. Their package contracts distinguish
the mathematical operations from source selection and physical calibration.
The common-history verifier also checks the retained summary and closed
source/evidence inventory. The preparation certificate specifies a quantum
control channel; it is not a sampled quantum execution.

The finite algebraic implications can be compiled independently:

```sh
cd Lean
lake build Geometry.NativeGeometricSourceIdentification CommonHistoryFeedback Geometry.SourceActionMeasure
```

The general spectral error estimate and oscillator pulse argument are
analytical; these finite Lean modules do not formalize those entire proofs.

## Passive-memory requirements of the earlier source lane

The exact nonlinear-record classification, finite native-zero reset bound,
fresh-zero approximate construction and ancillary quadratic budget have a
separate retained-evidence verifier:

```sh
python code/source_passive_memory/verify.py
python -m pytest -q code/source_passive_memory
cd Lean
lake build Geometry.SourcePassiveMemoryAxiomAudit
```

The [contract and scope](code/source_passive_memory/README.md) distinguish
native pair-mean executions from equal-image proof witnesses and supplied
preparations. These controls do not derive that lane's first milestone
(labelled M1 in its package) or a physical energy or clock.
The dedicated workflow runs on Windows and Linux; the frozen mandatory
runner is unchanged.

## Source-lane routing and live storage

The conditional routing and live-storage refinement of the same lane has a complete small
replay and an explicit kernel audit:

```bash
python code/source_routing_refinement/verify.py
python -m pytest -q code/source_routing_refinement/test_refinement.py
cd Lean
lake build Geometry.SourceRoutingRefinementAxiomAudit
```

The Python controls run in the dedicated Linux/Windows
`source-routing-refinement.yml` workflow; the source-frozen mandatory runner
is preserved byte for byte. See
[the derivation and scope](code/source_routing_refinement/README.md): the
23-scalar-slot bound excludes precision, controller and retained audit history,
and does not derive the supplied feedback law from canonical dynamics.

## Environment

- CPython 3.12 or newer for the full suite (continuous integration runs 3.12;
  the observer model alone runs on 3.11 or newer).
- A clean virtual environment for each pin set: the mandatory suite installs
  `requirements.txt` and the observer model installs `model/requirements.txt`.
  Their numpy and SciPy pins differ, and receipt checks are version
  sensitive, so the two gates run in separate environments.
- Tectonic 0.15.0 and Pandoc 3.8.3 for the publication artifacts.
- Ghostscript, `rsvg-convert`, and `pdftotext` on `PATH` for the book and
  publication validation. On Ubuntu these are supplied by `ghostscript`,
  `librsvg2-bin`, and `poppler-utils`.
- `xz` only for the optional NuFIT 6.1 profile replay described below. The
  profile files are external inputs and are not required by the mandatory
  clean-clone suite.

## Mandatory suite

```bash
python -m venv .venv
. .venv/bin/activate            # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python tools/run_mandatory_suite.py
```

`requirements.txt` pins the core dependencies. The default command runs the
complete standard suite. Its FZ-11 Lean replay test runs only when `lake` is
on `PATH` and Mathlib's oleans exist under `Lean/.lake`; on a clone without a
built Mathlib it is skipped, and the dedicated Lean CI lane replays it. On every push and PR, CI
(`.github/workflows/mandatory-suite.yml`) runs the same ordered steps in ten
isolated partitions on each operating system. Each partition has its own
clean checkout and the existing 30-minute job limit. All ten partitions must
succeed; a failed, cancelled or skipped partition cannot produce a passing
aggregate check. The runner also supports other partition counts; a complete
two-part local replay uses:

```bash
python tools/run_mandatory_suite.py --shard-index 0 --shard-count 2
python tools/run_mandatory_suite.py --shard-index 1 --shard-count 2
```

The indices are zero-based. Their ordered union is exactly the default suite;
running only one partition does not verify the complete suite. Partition flags
apply only to standard mode; full and certificate runs use their unpartitioned
commands. Five designated
replay/mutation-scan steps (listed in `HEAVY_STEP_TITLES` inside the runner)
are deferred to

```bash
python tools/run_mandatory_suite.py --full
```

which CI enforces nightly and on demand, and which the release checklist
runs before publication. The standard run prints exactly which steps it
skipped. The suite
both collects and executes: claim-registry validation, generated scientific
register validation, external-data provenance/hash/license-boundary validation,
release-manifest validation and its regression tests, a clean `--collect-only`
pass over `code/`, the scientific validation fixtures in
`code/audit/` (which includes the scope guard proving no cloud or hardware
lane is silently collected), the A5 closure ledger checks, and the earliest proof and
non-identifiability receipts (the suite's step titles call them Phase-0). One fixture is excluded from that step:
`code/audit/test_e4_absence_guards.py` needs the pinned Mathlib sources and
runs in the Lean CI workflow instead, where the toolchain is provisioned.

The exact certificate suites (port-current, matter-lift, axis-center-descent
and matter-menu spectral-ledger, ~26 minutes) run through the same runner in their own CI workflow
(`.github/workflows/certificate-suites.yml`) whenever `code/a5_closure/`
changes and nightly on both platforms:

```bash
python tools/run_mandatory_suite.py --certificates    # mandatory + certificates
python tools/run_mandatory_suite.py --certificates-only
```

## Optional lanes (opt-in extras)

Each optional lane keeps its own requirements file and stays out of the
mandatory collection unless explicitly enabled:

- IBM / Qiskit hardware lane:
  `pip install -r code/ibm_quantum_cloud/requirements-ibm.txt`, then
  `OPH_RUN_IBM=1 python -m pytest code/ibm_quantum_cloud`. A direct
  invocation without the opt-in, or with an incomplete extras installation,
  exits with the missing requirement instead of reporting an empty test run.
- Legacy particle helpers: set `OPH_RUN_LEGACY_D10=1` and
  `OPH_LEGACY_PARTICLE_DIR` (see `code/particles/conftest.py`).

## Scope

The mandatory suite is **collectable and executable** from a clean clone. The
acceptance bar for this path is a green `python tools/run_mandatory_suite.py`:
claim registry, release manifest, scientific-register sync, a clean `--collect-only`
run with zero import errors, and the executed validation fixtures. The generated
registers live in `docs/registers/`; the generation line under each page title
names its JSON source and the `tools/build_*.py` generator, and the suite fails when a committed
page differs from its regeneration. Each research note in `docs/research/` is
pinned by path and SHA-256 in the receipt of the code lane that verifies it, so a
change to a note rebuilds that receipt in the same commit.

Full test execution (`python -m pytest code`) is **not** expected to be green
from a clean clone, so it is not the documented gate here. A bare
`python -m pytest` collects `code/` only (`testpaths` in `pytest.ini`); the
regression tests under `tools/` run through the mandatory suite. The live
observer implementation requires the separate `python -m model.verify` and
`python -m pytest -q -W error model` gate above. Individual scientific
test outcomes are tracked as their own issues, and some are not reproducible
from the public checkout alone. In particular:

- Some byte- and value-level receipt checks are sensitive to platform line
  endings and to `numpy`/`scipy` versions.

Run the full suite for extended scientific validation, not as a clean-clone
pass/fail gate.

## Finite Core Checks

The compact exact evidence route is:

```bash
python3 -m pytest -q \
  code/a5_closure/test_audit.py \
  code/particles/calibration/test_wz_experimental_convention.py \
  code/particles/calibration/test_wz_survival_boundaries.py \
  code/capacity_readback/test_correctable_public_record_capacity.py \
  code/capacity_readback/test_reversible_public_checkpoint_packet.py \
  code/consensus/test_reference_architecture_benchmark_suite.py \
  code/consensus/test_verified_tree_packet_net.py
```

Run the independent strict-one-loop W/Z receipt package separately:

```bash
python3 code/particles/calibration/strict_one_loop_pole_map/run_all.py
```

This regenerates the conditional fixture receipt, runs the adversarial suite,
validates both JSON Schemas, checks the receipt without importing the producer,
and verifies the package manifest.

These checks cover the twelve-port algebra validation, exact physical-boundary
controls for the A5/SM and W/Z lanes, exact public-record capacity, the
reversible reference packet, and the finite consensus packets. They do not
claim a physical three-family attachment, an OPH-native W/Z pole, the missing
physical $N$ packet, or the continuum Einstein tower.

## Lean proofs

The Lean 4 / Mathlib workspace under `Lean/` holds ten libraries (`ObserverPatchHolography`, `EventAlgebra`, `OPHScreen` in `Screen/`, `OPHTime`, `OPHComputation`, `OPHThermodynamics`, `OPHConstruction`, `SourceLaw`, the opt-in `OphGap` certificate, and the standalone `ObservableNormalForms` package), each sorry-free. Thirteen finite proofs use `native_decide`; their native-evaluation axioms are tracked by `tools/check_lean_native_decide_inventory.py`. Rebuild everything with:

```bash
cd Lean
lake exe cache get
lake build
cd ObservableNormalForms
lake exe cache get
lake build
```

CI (`.github/workflows/lean-ci.yml`) runs both builds with a resumable cache,
rejects any `sorry`/`admit`/global-axiom regression, and replays the
Einstein-branch axiom check. `Lean/README.md` documents the layout;
`Lean/docs/PROOF_INDEX.md` maps theorems to paper statements.

## Paper review and publication builds

A new paper release follows the owner's approval of a recorded model
milestone under the
[model research policy](model/CONSTRUCTION.md#research-programme-and-paper-release-hold);
the [paper set](docs/PAPER_SET.md) records the milestone and the role of each
paper. Correctness fixes, errata, review builds and their validation continue
between releases. The publication commands below require that maintainer
decision or an authorized scoped corrective release.

With the pinned publication tools above installed, rebuild every registered
paper, including every extra and cosmology source, the warnings gate, the
local review manifest, and the reader-facing book from the repository root:

```bash
python3 tools/refresh_paper_release.py --preview
python3 tools/build_book_pdf.py
```

This chains `tools/build_tex_papers.py`, the build-warnings gate, manifest
regeneration, and manifest validation, so a rebuilt PDF can never be committed
with stale manifest hashes. This review pass may retain the visible release
identifier from an existing Git tag. It does not publish or replace that tag.
The CI build uses this mode so draft PDFs can be committed and inspected
without a version bump.

After that decision and review, prepare a publication candidate with a
distinct release identifier:

```bash
python3 tools/bump_paper_release.py
python3 tools/refresh_paper_release.py --publication
```

`tools/bump_paper_release.py` also rewrites repository URLs pinned to the
current release tag (`blob/<id>/`, `tree/<id>/`) in the `paper/`, `flagship/`,
`extra/`, `cosmology/` and `book/` TeX sources to the new tag
(disable with `--keep-release-urls`). Sources that cite an earlier tag
are caught up without a bump by
`python3 tools/bump_paper_release.py --retag-urls-only --from-release-id <old tag>`
(`--dry-run` lists the files).

The publication mode fails when the selected identifier exists as a local or
remote Git tag. It also rebuilds the canonical book and stamps its hash, size,
path, and selected release ID into the shared manifest before final
validation. Publication remains a separate maintainer action after the
candidate is committed, pushed, and inspected.

The manually dispatched `Release Channel Integrity` workflow is a
post-publication integrity check. Do not use it to validate a same-release preview. A
preview manifest describes the local bytes under review, while the tagged
GitHub Release remains fixed until the maintainer publishes a new release.

Both builders derive `SOURCE_DATE_EPOCH` from the visible date in
`paper/release_info.tex`, force UTC, avoid host-font selection, and retain the
logs consumed by the warning gate. The preview CI performs the sequence
twice on a clean Ubuntu runner and rejects any paper or book PDF whose SHA-256
changes on the second pass. A local check can make the same comparison with
`sha256sum` (or `shasum -a 256` on macOS).

The book's SVGs have source-bound canonical PDF renderings under
`assets/book_pdf_renderings/`. A normal build validates the exact SVG and PDF
inventory plus both SHA-256 digests, then stages those committed bytes. This
keeps librsvg, Pango, Cairo, and host-font differences outside the clean-clone
build boundary. After editing a book SVG, regenerate and validate the
renderings explicitly:

```bash
python3 tools/generate_book_pdf_assets.py
python3 tools/book_pdf_assets.py
```

The regeneration command requires `rsvg-convert` and Ghostscript. The normal
book build requires neither tool.

A clean clone must also retain no tracked publication drift after the first
rebuild:

```bash
git diff --exit-code -- paper flagship extra cosmology book/reverse-engineering-reality-book.pdf
```

The `Paper Preview Build` workflow enforces this check, compares the committed
book with a scratch rebuild, and rejects generated paper PDFs missing from
Git, including ignored files. The cosmology companions obey the same
reproducibility checks. A committed PDF and
manifest pair therefore cannot substitute another paper's bytes at the
expected path; the source rebuild restores the correct artifact and makes the
job fail.

## External comparison data

The mandatory suite validates
`code/audit/external_data_provenance_registry.json` with:

```bash
python3 tools/check_external_data_provenance.py
```

The registry pins each retained local artifact and loader by repository path,
byte count, and SHA-256; records its publisher, version, HTTPS source, and
license status; and distinguishes three boundaries:

- deterministic generators from hand-transcribed published constants
  (KNT19/PDG, the Planck Table 2 plus CODATA-G Gaussian approximation,
  the CODATA-2022 inverse-alpha comparison fixture, and the PDG-2026 W/Z
  running-width fixture);
- live PDG API snapshots whose normalized local artifacts are frozen but whose
  raw response bodies were not archived; and
- hash-pinned external NuFIT tables and the five-source Bouchard--Donagi
  literature packet, which are deliberately not vendored because their
  redistribution licenses are `NOASSERTION`.

The validator requires the complete twelve-artifact inventory; deleting a
dataset entry is itself a gate failure. Those declared gaps are data lineage,
not hidden build inputs. The paper and book build does not consume them. To
replay the optional NuFIT score, obtain the profile files listed in
`code/particles/neutrino/nufit61_sources.json`, keep their `.xz` bytes
unchanged, and pass the normal-ordering files to:

```bash
python3 code/particles/neutrino/score_neutrino_nufit61.py \
  --tb-off-no /path/to/v61.release-TBoff-NO.txt.xz \
  --tb-yes-no /path/to/v61.release-TByes-NO.txt.xz
```

The scorer checks the registered byte counts and SHA-256 values before parsing.
No credential, cloud cache, or untracked fixture is accepted as a mandatory
scientific input.
