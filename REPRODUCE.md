# Reproducing the mandatory scientific suite

This is the single clean-clone path for the mandatory OPH scientific receipt
suite. It does not prove any physics claim. It makes the existing machine
receipts runnable and auditable from a fresh checkout, so a green claim
registry is not mistaken for a green scientific reproduction.

## Environment

- CPython 3.11 or newer (verified on 3.12).
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
`sympy`, are imported across `code/` but were previously undeclared, so a fresh
clone failed at pytest collection before any test ran. With them declared, the
mandatory command above collects 918 tests and exits 0.

## Optional lanes (opt-in extras)

Each optional lane keeps its own requirements file and stays out of the
mandatory collection unless explicitly enabled:

- IBM / Qiskit hardware lane:
  `pip install -r code/ibm_quantum_cloud/requirements-ibm.txt`, then
  `OPH_RUN_IBM=1 python -m pytest code/ibm_quantum_cloud`.
- CAMB / Boltzmann cosmology:
  `pip install -r code/dark_matter/requirements-boltzmann.txt`.
- Legacy arXiv D10 helpers: set `OPH_RUN_LEGACY_D10=1` and
  `OPH_LEGACY_PARTICLE_DIR` (see `code/particles/conftest.py`).

## Scope

This PR makes the mandatory suite **collectable** from a clean clone. The
acceptance bar for this path is a green claim registry plus a clean
`--collect-only` (zero import errors, 918 tests).

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
