"""Send observed source symbols to the comparison station with native material swaps."""

import networkx as nx
import numpy as np
from scipy.linalg import expm

from . import reference_arrival as source

local = source.local
SEED = 28122925
TRIALS = 2**15


def transfer(path, bit, unitary, n):
    if not bit:
        return 0
    amplitude = np.zeros(n, complex)
    amplitude[path[0]] = 1
    for i, j in zip(path, path[1:]):
        amplitude[[i, j]] = unitary @ amplitude[[i, j]]
    result = abs(amplitude[path[-1]]) ** 2
    assert abs(result - round(result)) < 1e-12
    return round(result)


def run():
    rng = np.random.default_rng(SEED)
    board = local.survey.bulk.federation(0)
    sources = local.inputs("shared_contact")["serial_route"]["cells"]
    codes, common = local.survey.routes.codebook()
    graph = nx.Graph()
    for a, b in sorted(board.seams):
        e = a % 12, b % 12
        if a // 12 != b // 12 or (
            common[e] if a // 12 in sources else e not in codes[0]
        ):
            graph.add_edge(a, b)
    station = 0
    paths = []
    for v in sources:
        alternatives = [
            nx.shortest_path(graph, g, station)
            for g in sorted(graph)
            if g // 12 == v and nx.has_path(graph, g, station)
        ]
        paths.append(min(alternatives, key=lambda t: (len(t), t)))
    law = local.survey.routes.volume.candidate.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + 8 * contrast
    blank = np.eye(11)[:, 0]
    columns = np.kron(
        blank[:, None],
        np.eye(9)[
            :, [law.SRC.index((1, 0)), law.SRC.index((0, 1)), law.SRC.index((0, 0))]
        ],
    )
    u = expm(-1j * np.pi * h)
    projected = columns.conj().T @ u @ columns
    np.testing.assert_allclose(u @ columns, columns @ projected, atol=3e-13)
    vacuum = projected[2, 2]
    assert abs(abs(vacuum) - 1) < 1e-12
    op = projected[:2, :2] / vacuum
    np.testing.assert_allclose(abs(op), [[0, 1], [1, 0]], atol=3e-13)
    letters = local.inputs("pair_clock_information")["local_reader"]["calibrations"][
        "fixed_contact"
    ]["sharp_records_in_calibration_order"]
    alphabet = []
    for path in paths:
        received = []
        for energy, tag in letters:
            # Seven public integer bits report summed energy. One reports
            # arrival. A new blank carrier is used in each time slot.
            bits = [(energy >> j) & 1 for j in range(7)] + [tag]
            word = [transfer(path, bit, op, board.n) for bit in bits]
            received.append([sum((b << j) for j, b in enumerate(word[:7])), word[7]])
        assert received == letters
        alphabet.append(received)
    symbols = np.array([u + v for u in alphabet[0] for v in alphabet[1]])
    ref = local.survey.standards(board)[0][0].ravel() / 2592
    camera = local.inputs("reference_arrival")
    # Independently execute finite-length readouts and actual stopped courier
    # itineraries. The high-resolution receiver comparisons have separate,
    # batched cohorts and are not inferred from these path counts.
    steps, _, _ = source.coherent.steady.exp.flight.routes.layout(board)
    windows = {}
    for key, arm in camera["arms"].items():
        target = rng.choice(
            81, TRIALS, p=np.array(arm["source_program"]).ravel() / 2592
        )
        budget = source.panel_lengths(TRIALS, rng)
        live = np.arange(TRIALS)
        success = np.zeros(TRIALS, bool)
        count = 0
        while len(live):
            other = rng.choice(81, len(live), p=ref)
            count += len(live)
            live = live[np.any(symbols[target[live]] != symbols[other], axis=1)]
            budget[live] -= 1
            done = budget[live] == 0
            success[live[done]] = True
            live = live[~done]
        cohort = np.flatnonzero(success)
        current = np.zeros(TRIALS, dtype=int)
        history = np.zeros(TRIALS, dtype=int)
        while len(cohort):
            # The stopped service writes at each recorded station, including
            # the dispatching and terminal visits, under the original parity law.
            history[cohort] ^= (current[cohort] == 5).astype(int) * 2 + (
                current[cohort] == 19
            ).astype(int)
            cohort = cohort[rng.integers(32, size=len(cohort)) != 0]
            current[cohort] = steps[current[cohort], rng.integers(12, size=len(cohort))]
        windows[key] = {
            "comparisons": count,
            "accepted": int(success.sum()),
            "flag_counts_all_starts": np.bincount(history, minlength=4).tolist(),
        }
    return {
        "seed": SEED,
        "starts_per_arm": TRIALS,
        "collector": station,
        "paths": paths,
        "received_energy_and_arrival": alphabet,
        "complete_symbol_pairs": len(np.unique(symbols, axis=0)),
        "forward_paths": windows,
        "inputs": source.coherent.steady.exp.save_inputs(
            "ARCHIVED_OVERLAP.md",
            "reference_arrival.json",
            "shared_contact.json",
            "pair_clock_information.json",
        ),
    }


if __name__ == "__main__":
    z = run()
    source.coherent.steady.exp.emit("reference_bus", z)
    print(z["paths"], z["forward_paths"])
