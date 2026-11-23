"""Measure a single mobile defect in the cell's guarded colour band."""

import argparse
import hashlib
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.linalg import block_diag
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import candidate as contrast
from model.candidate_laws.interval_echo import candidate, currents
from model.executable_model import cell, federation
from model.experiments.interval_echo_encounter.conjugation import defect_generator
from model.experiments.interval_echo_gauge.__main__ import coordinates, ordered_word

HERE = Path(__file__).resolve().parent
TIMES = (1 / 16, 1 / 8, 1 / 4, 1 / 2, 1.0, 2.0, 4.0)
SHOTS = 2**28
SEED = 861221


def colour_meter():
    car = cell()
    _, source, h, poles, bands, _, _, _, addresses, _ = currents._service(
        car, False, 1.0
    )
    matrices = np.asarray(addresses)
    if [
        np.linalg.matrix_rank(coordinates(block))
        for block in (matrices[:, :3, :3], matrices[:, 3:, 3:], matrices)
    ] != [9, 3, 12]:
        raise ValueError("Wrong one-load current representation")
    basis = np.column_stack(bands)
    colour = bands[0] @ bands[0].T
    antipodes = currents._geometry(car)[0]
    exact = (
        (np.eye(12) - np.eye(12)[antipodes])
        @ (np.eye(12) + car.adjacency().toarray() / np.sqrt(5))
        / 4
    )
    np.testing.assert_allclose(colour, exact, atol=1e-13)
    np.testing.assert_allclose(h @ colour, poles[0]["energy"] * colour, atol=1e-13)
    axis = bands[0].T @ source[:, 0]
    axis /= np.linalg.norm(axis)
    target = block_diag(np.eye(3) / 3 - np.outer(axis, axis), np.zeros((3, 3)))
    coeff = np.linalg.lstsq(
        coordinates(matrices).T, coordinates([target])[0], rcond=1e-12
    )[0]
    np.testing.assert_allclose(
        np.einsum("i,ijk->jk", coeff, matrices), target, atol=1e-12
    )
    word = (
        np.eye(12)
        + basis @ (ordered_word(matrices, coeff, 2 * np.pi, 512) - np.eye(6)) @ basis.T
    )
    phase = np.exp(2j * np.pi / 3)
    error = np.linalg.norm(word - np.eye(12) - (phase - 1) * exact, 2)
    return (
        exact,
        word,
        antipodes,
        dict(
            triplet_energy=poles[0]["energy"],
            sharp_port_colour_fraction=0.25,
            compiled_centre_error=float(error),
            pulse_subdivisions=512,
        ),
    )


def gateway_constraints(board, antipodes):
    """Exact port-incidence elimination of the tangent conditions P_i H psi=0."""
    endpoints = {}
    cell_ports = [[] for _ in range(board.cells)]
    for edge, (a, b) in enumerate(board.external):
        for port in (a, b):
            cell_ports[port // 12].append(port % 12)
            endpoints[port] = edge
    linked = nx.Graph()
    linked.add_nodes_from(range(len(board.external)))
    forced_ports = set()
    dependent = 0
    for owner, local in enumerate(cell_ports):
        if len(local) != 3:
            raise ValueError("Expected three distinct native gateways per cell")
        pairs = [(p, q) for p in local for q in local if p < q and antipodes[p] == q]
        base = 12 * owner
        if not pairs:
            forced_ports.update(base + p for p in local)
        else:
            ((p, q),) = pairs
            dependent += 1
            (third,) = set(local) - {p, q}
            forced_ports.add(base + third)
            linked.add_edge(endpoints[base + p], endpoints[base + q])
    forced_edges = {endpoints[port] for port in forced_ports}
    # At alpha=1 this endpoint difference is the negative at the far end.
    # Equality of the antipodal coefficients transfers a zero throughout this graph.
    undecided_empty = [
        sorted(c) for c in nx.connected_components(linked) if c.isdisjoint(forced_edges)
    ]
    open_graph = nx.Graph()
    open_graph.add_nodes_from(range(board.cells))
    open_graph.add_edges_from(
        (a // 12, b // 12)
        for e, (a, b) in enumerate(board.external)
        if e not in forced_edges
    )
    components = list(nx.connected_components(open_graph))
    return dict(
        bridges=len(board.external),
        cells_with_antipodal_gateways=dependent,
        locally_forced_bridges=len(forced_edges),
        empty_background_undecided_components=undecided_empty,
        largest_possible_flow_component_for_signed_sectors=max(map(len, components)),
        largest_diameter_in_unconstrained_edges=max(
            nx.diameter(open_graph.subgraph(c)) for c in components
        ),
    )


def operator(board, arm, kappa):
    if arm == "one_load":
        h = contrast.generator(board, kappa)
        sign = 0
    else:
        sign = -1 if arm == "vacancy" else 1
        h = defect_generator(board, sign, kappa, 0)
    # Subtract a spatially constant phase. Compare to a separate direct reduction.
    constant = 6 * kappa + 2 * kappa * sign
    reduced = h - constant * sparse.eye(board.n, format="csr")
    check = board.laplacian() / 2 - sparse.diags(2 * kappa * sign * (board.degree - 5))
    np.testing.assert_allclose((reduced - check).data, 0, atol=1e-12)
    return reduced


def flight(board, level, arm, arm_index, projector, word, kappa, exchange_control):
    gate = min(a for a, b in board.external if a < 12)
    localized = np.eye(12)[:, gate] - projector[:, gate]
    localized /= np.linalg.norm(localized)
    launches = np.zeros((board.n, 2), complex)
    launches[:12] = np.column_stack((localized, np.ones(12) / np.sqrt(12)))
    gen = -1j * operator(board, arm, kappa)
    state = launches.copy()
    prev = 0
    rows = []
    for appointment, t in enumerate(TIMES):
        state = expm_multiply((t - prev) * gen, state)
        prev = t
        for source in range(2):
            outside = state[12:, source].reshape(-1, 12)
            arrived = np.linalg.norm(outside) ** 2
            colour = np.linalg.norm(outside @ projector.T) ** 2
            dark = np.linalg.norm(outside @ (word - np.eye(12)).T) ** 2 / 4
            p = np.array([colour, arrived - colour, 1 - arrived])
            if min(p) < -2e-13:
                raise ArithmeticError("Invalid observer probabilities")
            rng = np.random.default_rng(
                [SEED, level, arm_index, appointment, source, exchange_control]
            )
            counts = rng.multinomial(SHOTS, np.maximum(0, p) / np.maximum(0, p).sum())
            centre_counts = rng.multinomial(SHOTS, [dark, arrived - dark, 1 - arrived])
            rows.append(
                dict(
                    exposure=t,
                    source=["projected_gateway", "uniform_cell"][source],
                    outside_probability=float(arrived),
                    outside_colour_probability=float(colour),
                    colour_given_arrival=float(colour / arrived),
                    centre_dark_probability=float(dark),
                    energy_read=counts.tolist(),
                    centre_read=centre_counts.tolist(),
                )
            )
    # Cut the bridges after taking the same diagonal entries; intercell transfer then vanishes.
    interior = gen[:12, :12]
    halted = expm_multiply(4 * interior, launches[:12])
    np.testing.assert_allclose(np.linalg.norm(halted, axis=0), 1, atol=1e-12)
    return rows


def run():
    projector, word, antipodes, calibration = colour_meter()
    rows = []
    for level in range(3):
        board = federation(level)
        local = gateway_constraints(board, antipodes)
        arms = {}
        for i, arm in enumerate(("one_load", "vacancy", "excess")):
            arms[arm] = flight(
                board, level, arm, i, projector, word, candidate.KAPPA, 0
            )
        arms["exchange_only"] = flight(
            board, level, "one_load", 3, projector, word, 0.0, 1
        )
        rows.append(dict(cells=board.cells, constraints=local, arms=arms))
    return dict(
        protocol_sha256=hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        shots=SHOTS,
        seed=SEED,
        outcomes=["outside_in_band", "outside_in_complement", "home"],
        centre_outcomes=["outside_dark", "outside_bright", "home"],
        simultaneous_99_percent_probability_radius=float(
            np.sqrt(np.log(2 * 1008 / 0.01) / (2 * SHOTS))
        ),
        calibration=calibration,
        boards=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as file:
        json.dump(result, file, indent=2, allow_nan=False)
        file.write("\n")
