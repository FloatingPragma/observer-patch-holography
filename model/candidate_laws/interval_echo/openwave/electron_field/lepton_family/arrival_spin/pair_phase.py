"""Compare opposite physical half-turn orders on the blank-memory remote pair-phase preparation."""

import argparse
from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import connected_components

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[6]))
# Standalone archive replay starts outside the repository.
# ruff: noqa: E402
from model.candidate_laws.interval_echo.openwave.row_30_branch import (
    pair_phase as flight,
)
from model.executable_model import cell
from model.executable_model.symmetry import rotations
from model.experiments.interval_echo_angular.charged_orientation.measure import (
    compile_turn,
)


def addresses(sector, station):
    u, v, w, x = flight.rails(station)
    return [sector.index(u, w), sector.index(v, x)]


def geometric_orders():
    # In a protected blank contact an addressed pi pulse swaps the two loads, including 00 and 11.
    law = flight.candidate.native_law
    exchange, contrast = law.local_term_parts()
    pulse = expm(-1j * np.pi * (exchange + flight.candidate.KAPPA * contrast))
    for word in ((0, 0), (0, 1), (1, 0), (1, 1)):
        expected = np.zeros(law.DIM)
        expected[law.SRC.index(word[::-1])] = 1
        np.testing.assert_allclose(
            pulse[:, law.SRC.index(word)], expected, atol=1e-12, rtol=0
        )
    sector = flight.PairSector(cell())
    n = sector.P
    # Off-diagonal exchanges connect all resolved pair addresses in the protected template.
    graph = sector.exchange().copy()
    graph.setdiag(0)
    graph.eliminate_zeros()
    components = connected_components(graph, directed=False, return_labels=False)
    if components != 1:
        raise ArithmeticError("Disconnected pair sector")
    relabelings = []
    for rotation in sorted(rotations()):
        if rotation == tuple(range(12)) or any(
            rotation[rotation[i]] != i for i in range(12)
        ):
            continue
        port_registers = np.arange(12)
        for a, b in compile_turn(rotation):
            port_registers[[a, b]] = port_registers[[b, a]]
        image = np.argsort(port_registers)
        np.testing.assert_array_equal(image, rotation)
        perm = np.array(
            [sector.index(image[a], image[b]) for a, b in zip(sector.aa, sector.bb)]
        )
        np.testing.assert_array_equal(np.sort(perm), np.arange(n))
        relabelings.append((rotation, perm))
    orders = []
    for (a, first), (b, second) in combinations(relabelings, 2):
        if all(a[b[i]] == b[a[i]] for i in range(12)):
            xy, yx = second[first], first[second]
            np.testing.assert_array_equal(xy, yx)
            orders.append((np.argsort(xy), np.argsort(yx)))
    if len(orders) != 15:
        raise ArithmeticError("Incomplete perpendicular-axis census")
    return sector, orders, int(components)


def run():
    local, orders, components = geometric_orders()
    board = flight.federation(0)
    pairs = flight.PairSector(board)
    launch, analyzer = flight.apparatus()
    archive = json.loads(flight.PREDICTIONS.read_text())
    observed = {
        (source, dest): bins
        for source, dest, law, t, *bins in archive["receiver"]
        if law == "native" and t == 256
    }
    sectors = []
    for source in (2, 3):
        initial = np.zeros((pairs.P, 2))
        initial[addresses(pairs, source), np.arange(2)] = 1
        (states,) = flight.evolve(
            pairs.hamiltonian(flight.candidate.KAPPA),
            flight.candidate.KAPPA,
            initial,
            [256],
        )
        for dest in flight.receivers(board, source)[1:]:
            k = states[addresses(pairs, dest)]
            received = k @ launch
            records = flight.observe(
                pairs, dest, states @ launch, np.zeros(2), analyzer
            )
            np.testing.assert_allclose(
                records, observed[source, dest], atol=2e-9, rtol=0
            )
            full = np.zeros((local.P, 2), complex)
            full[addresses(local, 0)] = received
            admitted = np.sum(abs(full) ** 2, axis=0)
            np.testing.assert_allclose(
                admitted, np.array(records)[:, -2:].sum(axis=1), atol=2e-13
            )
            dark = []
            for xy, yx in orders:
                left, right = full[xy], full[yx]
                np.testing.assert_allclose(
                    np.sum(abs((left + right) / 2) ** 2, axis=0), admitted, atol=2e-13
                )
                dark.append(np.sum(abs((left - right) / 2) ** 2, axis=0))
            # An ideal input reference maximally entangled with the two sharp source words.
            # K is the unnormalized terminal filter given a prepared input; exterior occupancy is vacuum.
            sv = np.linalg.svd(k, compute_uv=False)
            sectors.append(
                {
                    "source": source,
                    "receiver": dest,
                    "exposure": 256,
                    "native_plus_minus_admission_per_enrollment": admitted.tolist(),
                    "maximum_destructive_weight_over_15_axis_pairs": np.max(
                        dark, axis=0
                    ).tolist(),
                    "source_word_to_receiver_word_singular_values": sv.tolist(),
                    "computed_reference_partial_transpose_min_given_prepared_input": float(
                        -sv.prod() / 2
                    ),
                }
            )
    return {
        "two_load_local_addresses": local.P,
        "exchange_graph_components": components,
        "perpendicular_axis_pairs": len(orders),
        "routes": sectors,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = run()
    with args.output.open("x") as out:
        json.dump(data, out, indent=2, allow_nan=False)
        out.write("\n")
