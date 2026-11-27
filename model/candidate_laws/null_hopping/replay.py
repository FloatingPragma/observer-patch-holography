"""First seam reading after a free release on one cell, with matched preparation and dwell."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector, pair_outcomes
from model.executable_model import cell
from . import candidate


def pair_readout(sector, kappa):
    initial = np.zeros(sector.P)
    initial[sector.index(0, 1)] = 1
    state = expm_multiply(-1j * sector.hamiltonian(kappa), initial)
    words = pair_outcomes(sector, state)
    return {"no_pair": float(words[0]), "by_seam": words[1:].tolist()}


def experiment():
    board = cell()
    sector = PairSector(board)
    return {
        "start": [0, 1],
        "exchange_time": 1,
        "seams": [list(edge) for edge in board.seams],
        "kappa_0": pair_readout(sector, candidate.KAPPA),
        "kappa_8": pair_readout(sector, 8),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="compare every outcome with the stored probabilities")
    args = parser.parse_args(argv)
    actual = experiment()
    if not args.check:
        print(json.dumps(actual, indent=2, allow_nan=False))
        return 0
    stored = json.loads(Path(__file__).with_name("results.json").read_text(encoding="utf-8"))
    for key in ("start", "exchange_time", "seams"):
        if actual[key] != stored[key]:
            raise ValueError(f"pair readout {key} differs")
    for key in ("kappa_0", "kappa_8"):
        result, expected = actual[key], stored[key]
        if not math.isclose(result["no_pair"], expected["no_pair"], rel_tol=0, abs_tol=1e-12):
            raise ValueError(f"{key}: no-pair probability differs")
        np.testing.assert_allclose(result["by_seam"], expected["by_seam"], rtol=0, atol=1e-12)
    print("Both complete 31-outcome distributions agree with the stored probabilities.")
    return 0
