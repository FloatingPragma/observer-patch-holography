"""Finite signed angular panels on the same 80-cell equilibrium preparation."""

import json

import numpy as np
from scipy.sparse.linalg import expm_multiply

from matched_einstein_limit.rate_geometry.static_directions import design
from ..measure import pairing, ranging
from .geometry import HERE


def transfer():
    directions, _, _, _ = design()
    board = pairing.formation.federation(1)
    # Local angular readout enrolls only the first cell's 12 native ports.
    adjacency = board.adjacency()[:12, :12].toarray()
    identity = np.eye(12)
    weighted = 2 * identity + adjacency
    quadrupole = (
        directions[:, :, None] * directions[:, None, :] - np.eye(3) / 3
    ).reshape(12, -1)
    np.testing.assert_allclose(
        adjacency @ directions, np.sqrt(5) * directions, atol=1e-12
    )
    np.testing.assert_allclose(adjacency @ quadrupole, -quadrupole, atol=1e-12)
    basis = np.column_stack([np.ones(12), directions, quadrupole])
    ranks = [
        int(np.linalg.matrix_rank(matrix @ basis, tol=1e-10))
        for matrix in (identity + adjacency, weighted)
    ]
    return directions, weighted, ranks


def acquire():
    directions, weighted, ranks = transfer()
    receipt = json.loads(
        (HERE.parents[1] / "full_algebra" / "results.json").read_text()
    )
    state = np.asarray(receipt["x_real"]) + 1j * np.asarray(receipt["x_imag"])
    board = pairing.formation.federation(1)
    source = pairing.formation.region(board)
    emitter = source[0]
    dwell = 1.5
    starts = 2**26
    occupations = np.zeros(board.n)
    for odd, block in enumerate(state):
        for charge in (-1, 1):
            h = ranging.generator(board, source, charge, 8, bool(odd))
            initial = np.zeros((h.shape[0], 4), complex)
            initial[4 * emitter : 4 * emitter + 4] = np.eye(4)
            waves = expm_multiply(-1j * dwell * h, initial, traceA=0).reshape(
                board.n, 4, 4
            )
            occupations += (
                np.einsum("xai,ij,xaj->x", waves, block, waves.conj()).real / 2
            )
    np.testing.assert_allclose(sum(occupations), 1.0, atol=5e-13, rtol=0.0)
    probabilities = np.append(occupations[:12], sum(occupations[12:]))
    counts = np.random.default_rng(1370519).multinomial(starts, probabilities)
    # Bounded per-start score: one for the pole, one half for each of its
    # neighbors, zero at all other ports (including the outside bin).
    panels = (weighted @ counts[:12]) / (2 * starts)
    expected = weighted @ probabilities[:12] / 2
    error = np.sqrt(np.log(2 * 12 / 0.01) / (2 * starts)) + 1.81e-6
    assert np.max(abs(panels - expected)) < error
    opp = (directions @ directions.T).argmin(axis=1)
    even = (panels + panels[opp]) / 2
    odd = (panels - panels[opp]) / 2
    return dict(
        cells=board.cells,
        launch_port=int(emitter),
        dwell=dwell,
        starts=starts,
        counts_including_outside=counts.tolist(),
        weighted_panel_score=panels.tolist(),
        simultaneous_error=float(error),
        maximum_odd_amplitude=float(np.max(abs(odd))),
        even_panel_spread=float(np.ptp(even)),
        ordinary_rank=ranks[0],
        pole_weighted_rank=ranks[1],
    )


if __name__ == "__main__":
    print(json.dumps(acquire(), indent=2))
