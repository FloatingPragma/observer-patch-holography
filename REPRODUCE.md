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
python -m pytest code
```

`requirements.txt` pins the core dependencies. Two of them, `mpmath` and
`sympy`, are imported across `code/` but were previously undeclared, so a fresh
clone failed at pytest collection before any test ran.

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

The command above collects with zero import errors from a clean clone.
Individual scientific test outcomes are separate from collection health and are
tracked as their own issues.
