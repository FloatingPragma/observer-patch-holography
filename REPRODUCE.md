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
python tools/check_claim_registry.py
python -m pytest --collect-only code
```

`requirements.txt` pins the core dependencies. Two of them, `mpmath` and
`sympy`, are imported across `code/` and require explicit declarations. With
them declared, the mandatory command above collects the repository test set
and exits 0.

## Optional lanes (opt-in extras)

Each optional lane keeps its own requirements file and stays out of the
mandatory collection unless explicitly enabled:

- IBM / Qiskit hardware lane:
  `pip install -r code/ibm_quantum_cloud/requirements-ibm.txt`, then
  `OPH_RUN_IBM=1 python -m pytest code/ibm_quantum_cloud`.
- CAMB / Boltzmann cosmology:
  `pip install -r code/dark_matter/requirements-boltzmann.txt`.
- Legacy particle helpers: set `OPH_RUN_LEGACY_D10=1` and
  `OPH_LEGACY_PARTICLE_DIR` (see `code/particles/conftest.py`).

## Scope

The mandatory suite is **collectable** from a clean clone. The acceptance bar
for this path is a green claim registry plus a clean `--collect-only` run with
zero import errors.

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

## Paper Build

With [Tectonic](https://tectonic-typesetting.github.io/) installed, rebuild the
complete paper stack from the repository root:

```bash
python3 tools/build_tex_papers.py
```

Every generated paper displays the shared release identifier from
`paper/release_info.tex`. The manifest records the corresponding PDF hashes.
