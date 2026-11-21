"""Compare the independent replay with the registered law in this RER checkout."""

import sys
from pathlib import Path

import numpy as np
from scipy import sparse

from row_04_check import screen


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.executable_model.carrier import federation


def main():
    _, _, exchange, adjacency = screen(HERE / "FEDERATION_LEVEL0.npz")
    sector = PairSector(federation(0))
    for kappa in (candidate.KAPPA, 0):
        independent = exchange - 2 * kappa * sparse.diags(adjacency)
        if (independent != sector.hamiltonian(kappa)).nnz:
            raise ValueError("The independent and native sector matrices differ")
    law = candidate.native_law
    echo = law.echo()
    outgoing = np.arange(law.DIM) % len(law.SRC) == law.SRC.index((1, 1))
    for loads in ((0, 0), (0, 1), (1, 0), (1, 1)):
        state = np.eye(law.DIM)[:, law.SRC.index(loads)]
        pair = loads == (1, 1)
        written = echo @ state
        for flag in (True, False):
            recovered = echo.T @ (written * (outgoing == flag))
            np.testing.assert_allclose(recovered, state * (pair == flag), atol=1e-12)
    print("Both 28,680-state generators agree entrywise. All eight terminal branches agree.")


if __name__ == "__main__":
    main()
