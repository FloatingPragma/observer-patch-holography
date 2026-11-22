# Reproducing the mandatory scientific suite

This is the clean-clone path into the OPH scientific receipt suite. The first
commands verify that the claim registry is internally connected and that the
public test collection imports without error. Individual evidence families
then have their own theorem, certificate, or experimental acceptance rule.

## Executable model

The [executable model](model/README.md) has a fixed carrier, the spherical echosahedral federation in
`model/executable_model`, and candidate laws in `model/candidate_laws`. The main candidate is
[Interval Echo](model/candidate_laws/interval_echo/README.md). Its [specification](model/candidate_laws/interval_echo/SPEC.md) defines the material clock and the prepared service protocol. From the repository root run

```sh
python -m pip install -r model/requirements.txt
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -m model                              # 20-cell material clock and paired reunion cohorts
python -m model --full                       # all three sizes, stationary records, repair forks and port currents
python -m model run --law contrast_memory_spherical_v1   # reference binding and reader replay
python -m model run --law scalar_channel
python -m model evaluate --check
python -m pytest -q -W error model
```

The main replay compares fresh clock decisions and reunion words with the deposited results. Full verification includes every stationary record in the archive, the two-seam repair outcomes and the guarded port-current response. A mismatch returns a nonzero exit status. Use `--full --seed 42 --output /tmp/echo-42` to generate an independent realization. These commands run in a standalone RER checkout.

The station experiment tests whether ranges inferred from simulated observer records predict held clock and pulse measurements. Its three-dimensional fit uses the recorded times and station labels. [Source-free chart and held flights](model/experiments/interval_echo_space/README.md#can-a-source-free-chart-predict-clocks-and-pulses-at-other-stations)

## Theory certificates

These packages check conditional theorems on their own declared carriers. Their replay establishes no transfer
to a candidate law without an explicit source map:

```sh
PYTHONPATH=code python -m pytest -q -W error \
  code/carrier_selection_reduction code/defect_particle_classes \
  code/finite_instrument_calibration code/matched_einstein_limit \
  code/maxwell_local_selection code/maxwell_weak_field \
  code/support_refinement
```

## Additional component calculations

The bounded modal Maxwell-shaped factorization is also an exact, laptop-scale
algebra check rather than a simulation campaign:

```bash
python -m pytest -q code/electromagnetism/test_modal_maxwell_factorization.py
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
receipt regeneration commands. The `neutral_packet` lane of the
[CI workflow](tools/CI.md) executes these controls on Linux.

## Environment

- CPython 3.12 or newer for the full suite (continuous integration runs 3.12;
  the executable model alone runs on 3.11 or newer).
- A clean virtual environment for each pin set: the mandatory suite installs
  `requirements.txt` and the executable model installs `model/requirements.txt`.
  Their NumPy and SciPy pins differ, so the two gates run in separate environments.
- Tectonic 0.15.0 and Pandoc 3.8.3 for the publication artifacts.
- Ghostscript, `rsvg-convert`, and `pdftotext` on `PATH` for the book and
  publication validation. On Ubuntu these are supplied by `ghostscript`,
  `librsvg2-bin`, and `poppler-utils`.

The full model gate runs on native x86-64 Linux in GitHub Actions. Exact seeded draws and long-time numerical receipts can depend on the host math libraries even with the same Python package versions. On Apple Silicon, the following Linux container command reproduces the station-chart tests and the archived electron-wall and loaded-clock counts. The complete native-runner command is in the executable model section above.

```sh
docker run --rm --platform linux/amd64 -v "$PWD:/src:ro" -w /src \
  -e PYTHONDONTWRITEBYTECODE=1 -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 \
  python:3.12-slim sh -ec 'python -m pip install -r model/requirements.txt &&
    python -m pytest -p no:cacheprovider -q -W error \
      model/experiments/interval_echo_space/test_station_chart.py \
      model/experiments/seam_program/test_blind_chart.py \
      model/benchmarks/eight_slot_torus/electron_wall/test_electron_wall.py::test_published_counts_all_settings \
      model/benchmarks/eight_slot_torus/formed_clocks/test_loaded_clock_archive.py::test_published_source_and_signed_counts_replay'
```

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
(`.github/workflows/ci.yml`) runs the same ordered steps in ten
isolated partitions. Each partition has its own clean checkout and a hard
job timeout. All ten partitions must
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
commands. CI runs the whole standard list in six shards. The suite
both collects and executes: claim-registry validation, generated scientific
register validation, external-data provenance/hash/license-boundary validation,
release-manifest validation and its regression tests, a clean `--collect-only`
pass over `code/`, the scientific validation fixtures in
`code/audit/` (which includes the scope guard proving no cloud or hardware
lane is silently collected), the A5 closure ledger checks, and the earliest proof and
non-identifiability receipts (the suite's step titles call them Phase-0). One fixture is excluded from that step:
`code/audit/test_e4_absence_guards.py` needs the pinned Mathlib sources and
runs in the Lean CI workflow instead, where the toolchain is provisioned.

The exact certificate suites (port-current and matter-lift, about 26 minutes) run on demand through the same
runner; the canonical smoke receipts run on every push in `.github/workflows/ci.yml`:

```bash
python tools/run_mandatory_suite.py --certificates    # mandatory + certificates
python tools/run_mandatory_suite.py --certificates-only
```

## Scope

The mandatory suite is **collectable and executable** from a clean clone. The
acceptance bar for this path is a green `python tools/run_mandatory_suite.py`:
claim registry, release manifest, scientific-register sync, a clean `--collect-only`
run with zero import errors, and the executed validation fixtures. The generated
registers live in `docs/registers/`; the generation line under each page title
names its JSON source and the `tools/build_*.py` generator, and the suite fails when a committed
page differs from its regeneration. Scientific receipts bind executable sources, formal proofs and machine inputs.
Research notes explain those results; they are not scientific inputs. Explicit
retained documentary contracts identify their authenticated input bytes and do
not assert that those bytes validate the current explanatory prose.

Full test execution (`python -m pytest code`) is **not** expected to be green
from a clean clone, so it is not the documented gate here. A bare
`python -m pytest` collects `code/` only (`testpaths` in `pytest.ini`); the
regression tests under `tools/` run through the mandatory suite. The
executable model uses the full model replay, candidate evaluation and
`python -m pytest -q -W error model` as specified above. Individual scientific
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
  code/capacity_readback/test_correctable_public_record_capacity.py \
  code/capacity_readback/test_reversible_public_checkpoint_packet.py \
  code/consensus/test_reference_architecture_benchmark_suite.py \
  code/consensus/test_verified_tree_packet_net.py
```

These checks cover the twelve-port algebra validation, exact public-record
capacity, the reversible reference packet, and the finite consensus packets.
They do not claim a physical three-family attachment, the physical $N$
packet, or the continuum Einstein tower.

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

CI (`.github/workflows/lean.yml`) runs both builds with a resumable cache,
rejects any `sorry`/`admit`/global-axiom regression, and replays the
Einstein-branch axiom check. `Lean/README.md` documents the layout;
`Lean/docs/PROOF_INDEX.md` maps theorems to paper statements.

## Paper review and publication builds

A new paper release follows the owner's approval of a recorded model
milestone of the [executable model](model/README.md); the
[paper set](docs/PAPER_SET.md) records the role of each paper. Correctness fixes, errata, review builds and their validation continue
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
  (KNT19/PDG hadronic data, the Planck Table 2 plus CODATA-G Gaussian
  approximation, and the CODATA-2022 inverse-alpha comparison fixture);
- live PDG API snapshots whose normalized local artifacts are frozen but whose
  raw response bodies were not archived; and
- hash-pinned external literature packets, which are deliberately not
  vendored because their redistribution licenses are `NOASSERTION`.

The validator requires the complete artifact inventory; deleting a dataset
entry is itself a gate failure. Those declared gaps are data lineage, not
hidden build inputs. The paper and book build does not consume them.
No credential, cloud cache, or untracked fixture is accepted as a mandatory
scientific input.
