# Reproducing the mandatory scientific suite

This is the clean-clone path into the OPH scientific receipt suite. The first
commands verify that the claim registry is internally connected and that the
public test collection imports without error. Individual evidence families
then have their own theorem, certificate, or experimental acceptance rule.

## Environment

- CPython 3.12 or newer (verified on 3.12 and 3.13).
- A clean virtual environment.

## Mandatory suite

```bash
python -m venv .venv
. .venv/bin/activate            # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python tools/run_mandatory_suite.py
```

`requirements.txt` pins the core dependencies. The runner is the single
documented mandatory command, and it is the exact command CI
(`.github/workflows/mandatory-suite.yml`) enforces on every push and PR. It
both collects and executes: claim-registry validation (including live-gate
and wording checks), release-manifest validation and its regression tests,
the generated claims-scoreboard drift check, a clean `--collect-only` pass
over `code/`, the audit fixture suite in `code/audit/` (which includes the
scope guard proving no cloud or hardware lane is silently collected), and
the A5 closure ledger audit.

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
  `OPH_RUN_IBM=1 python -m pytest code/ibm_quantum_cloud`.
- Legacy particle helpers: set `OPH_RUN_LEGACY_D10=1` and
  `OPH_LEGACY_PARTICLE_DIR` (see `code/particles/conftest.py`).

## Scope

The mandatory suite is **collectable and executable** from a clean clone. The
acceptance bar for this path is a green `python tools/run_mandatory_suite.py`:
claim registry, release manifest, scoreboard sync, a clean `--collect-only`
run with zero import errors, and the executed audit fixtures.

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

Run the full suite for scientific auditing, not as a clean-clone pass/fail gate.

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

These checks cover the twelve-port algebra audit, exact physical-boundary
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
Einstein-branch axiom audit. `Lean/README.md` documents the layout;
`Lean/docs/PROOF_INDEX.md` maps theorems to paper statements.

## Paper Build

With [Tectonic](https://tectonic-typesetting.github.io/) installed, rebuild the
complete paper stack, warnings gate, and release manifest in one pass from the
repository root:

```bash
python3 tools/refresh_paper_release.py
```

This chains `tools/build_tex_papers.py`, the build-warnings gate, manifest
regeneration, and manifest validation, so a rebuilt PDF can never be committed
with stale manifest hashes. Every generated paper displays the shared release
identifier from `paper/release_info.tex`; bumping that identifier
(`tools/bump_paper_release.py`) is a separate, deliberate release step.
