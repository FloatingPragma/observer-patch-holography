"""Enumerate every prospective setting as development evidence on the fixed board."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from .. import gateway_independent as peer
from ..experiment import gateway_table
from ..gateway_acquisition import prepare_layout, terminal_port_effects
from ..gateway_bound import SCORE, arbitrary_sweep_floor, board_charts
from ..gateway_dephased_bound import dephased_home_ceiling
from ..gateway_release import ARMS
from ..gateway_transport_bound import swap_home_ceiling
from .experiment import terms


def survey():
    spec = terms()
    gates = spec["gates"]
    board = federation(0)
    h = candidate.generator(board)
    heat = -board.laplacian() / 2
    # Keep full-board transport. Removing a common phase preserves every endpoint effect.
    h = h - sparse.eye(board.n) * h.diagonal().mean()
    times = np.arange(15, 39)
    initial = np.eye(board.n)
    unitary = expm_multiply(-1j * h, initial, start=15, stop=38, num=24, traceA=0)
    stochastic = expm_multiply(heat, initial, start=15, stop=38, num=24)
    np.testing.assert_allclose(terminal_port_effects(), np.eye(12), atol=1e-13)
    output = np.empty((20, 12, 2, 24, 4, 14))
    cross_error = 0.0
    for home, chart in enumerate(gateway_table(board)):
        ids = slice(12 * home, 12 * (home + 1))
        local = unitary[:, ids, ids]
        sharp = abs(local) ** 2
        prep = prepare_layout(tuple(chart))
        # Unnormalized clear columns are common to all start orders.
        formed = abs(local @ prep.accepted_vectors) ** 2
        columns = abs(prep.accepted_vectors) ** 2
        for start in range(12):
            for d, direction in enumerate((-1, 1)):
                weights = np.zeros(12)
                reach = 1.0
                for attempt in range(48):
                    port = (start + direction * attempt) % 12
                    weights[port] += reach
                    reach *= prep.first_click[:, :, port].sum()
                diagonal = columns @ weights
                inside = (
                    formed @ weights,
                    sharp @ diagonal,
                    sharp.mean(axis=2),
                    stochastic[:, ids, ids] @ diagonal,
                )
                for arm, values in enumerate(inside):
                    failure = 0.0 if ARMS[arm] == "unfiltered_uniform" else reach
                    row = output[home, start, d, :, arm]
                    row[:, :12] = values
                    row[:, 12] = 1 - failure - values.sum(axis=1)
                    row[:, 13] = failure
        # Separate implementation uses dense raw-edge spectral flight, literal gateway
        # queries, Poisson summation and a sum over every exterior destination.
        start, d, t = home % 12, home % 2, 15 + home
        separate = peer.probabilities(peer.BOARD, home, start, (-1, 1)[d], t)
        for j, arm in enumerate(ARMS):
            delta = max(abs(output[home, start, d, t - 15, j] - separate[arm]))
            cross_error = max(cross_error, float(delta))
    if (
        not np.isfinite(output).all()
        or output.min() < -2e-10
        or output.max() > 1 + 2e-10
        or abs(output.sum(axis=-1) - 1).max() > gates["normalization_tolerance"]
        or cross_error > gates["independent_bin_tolerance"]
    ):
        raise ArithmeticError("Invalid grid or independent check")
    # Small negative complements can arise only at roundoff. Return the full computed grid.
    homes = output[..., :12].sum(axis=-1)
    coherent, dephased, _, swap = np.moveaxis(homes, -1, 0)
    charts = board_charts()
    per_chart = {c: arbitrary_sweep_floor(c) for c in set().union(*charts.values())}
    floors = {
        str(n): min(per_chart[c] for c in chartset) for n, chartset in charts.items()
    }
    dephased_bound = dephased_home_ceiling()[-1]
    swap_bound = swap_home_ceiling()[-1]
    if (
        min(floors.values()) / SCORE <= gates["exact_joint_home_strictly_above"]
        or coherent.min() <= gates["coherent_terminal_home_strictly_above"]
        or swap.max() >= gates["stochastic_swap_home_strictly_below"]
        or (coherent - dephased).min() < gates["minimum_coherent_minus_dephased"]
        or min(floors.values()) / SCORE - float(dephased_bound)
        < gates["minimum_coherent_minus_dephased"]
        or float(swap_bound) >= gates["stochastic_swap_home_strictly_below"]
    ):
        raise ArithmeticError(
            "A response or certificate failed a prospective threshold"
        )
    summary = dict(
        kind="exhaustive development evidence",
        official_confirmation=False,
        settings=int(np.prod(output.shape[:4])),
        bins=int(output.size),
        shape=list(output.shape),
        directions=[-1, 1],
        times=times.tolist(),
        arms=ARMS,
        all_time_floor_ppt_by_board=floors,
        grid_minimum_coherent=float(coherent.min()),
        grid_maximum_stochastic=float(swap.max()),
        grid_minimum_coherent_advantage=float((coherent - dephased).min()),
        dephased_ceiling=float(dephased_bound),
        stochastic_ceiling=float(swap_bound),
        raw_edge_checked_settings=20,
        raw_edge_checked_bins=1120,
        max_raw_edge_bin_error=cross_error,
    )
    return summary, output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bins", type=Path, help="optional compressed full-precision grid"
    )
    args = parser.parse_args()
    result, grid = survey()
    if args.bins:
        with args.bins.open("xb") as stream:
            np.savez_compressed(stream, probabilities=grid)
    print(json.dumps(result, indent=2, allow_nan=False))
