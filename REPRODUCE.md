# Reproducing the mandatory scientific suite

This is the clean-clone path into the OPH scientific receipt suite. The first
commands verify that the claim registry is internally connected and that the
public test collection imports without error. Individual evidence families
then have their own theorem, certificate, or experimental acceptance rule.

The bounded modal Maxwell-shaped factorization is also an exact, laptop-scale
algebra check rather than a simulation campaign:

```bash
python -m pytest -q \
  code/electromagnetism/test_modal_maxwell_factorization.py \
  tools/test_modal_maxwell_factorization_surfaces.py
```

## Environment

- CPython 3.12 or newer (verified on 3.12 and 3.13).
- A clean virtual environment.
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

`requirements.txt` pins the core dependencies. The runner is the single
documented mandatory command, and it is the exact command CI
(`.github/workflows/mandatory-suite.yml`) enforces on every push and PR. The
default run is the standard suite; five long-running replay/mutation-scan
steps (listed in `HEAVY_STEP_TITLES` inside the runner, together most of the
suite's runtime) are deferred to

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
lane is silently collected), the A5 closure ledger checks, and the Phase-0
proof/non-identifiability receipts. One fixture is excluded from that step:
`code/audit/test_e4_absence_guards.py` needs the pinned Mathlib sources and
runs in the Lean CI workflow instead, where the toolchain is provisioned.

The exact certificate suites (#566 port-current, #314 matter-lift, ~26
minutes) run through the same runner in their own CI workflow
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
run with zero import errors, and the executed validation fixtures.

Full test execution (`python -m pytest code`) is **not** expected to be green
from a clean clone, so it is not the documented gate here. Individual scientific
test outcomes are tracked as their own issues, and some are not reproducible
from the public checkout alone. In particular:

- Two runtime-surface tests in
  `code/particles/test_compute_current_output_table_runtime_surface.py` require
  the untracked sibling tree `../arXiv/RC1/ancillary/code/particles`, which a
  clean clone does not provide.
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

The Lean 4 / Mathlib workspace under `Lean/` holds four libraries
(`ObserverPatchHolography`, `EventAlgebra`, `OPHScreen` in `Screen/`, and the
standalone `ObservableNormalForms` package), each sorry-free with standard
axioms only. Rebuild everything with:

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

With the pinned publication tools above installed, rebuild every registered
paper, the warnings gate, the local review manifest, and the reader-facing
book from the repository root:

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

After review, prepare a publication candidate with a new release identifier:

```bash
python3 tools/bump_paper_release.py
python3 tools/refresh_paper_release.py --publication
```

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
git diff --exit-code -- paper flagship extra book/reverse-engineering-reality-book.pdf
```

The `Paper Preview Build` workflow enforces this check. A committed PDF and
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

The validator requires the complete nine-artifact inventory; deleting a
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
