"""Positive formation spectrum of a charged vacancy with native zero-contrast marks."""

import argparse
import base64
from itertools import combinations
import json
from math import pi, sqrt
from pathlib import Path
import zlib

import numpy as np
from scipy.optimize import minimize_scalar
import networkx as nx
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_01_charge.measure import (
    sharp_contrast,
)
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
SHOTS = 2**15
TIMES = np.arange(1024) / 4


def seam_certificate():
    law = candidate.native_law
    blank, zero = (
        list(range(9)),
        [
            9 + 9 * i + j
            for i, (_, d) in enumerate(law.LABELS)
            if d == 0
            for j in range(9)
        ],
    )
    # Exact integer write pattern, 2 sqrt(2) factors removed from projector products.
    a = law._write_pattern().astype(np.int64)
    b = 2 * np.eye(90, dtype=np.int64) - a @ a.T
    refl, d2 = law._label_ops()
    fields = (np.eye(10, dtype=int) - refl, 2 * d2)
    for twice_g, actual in zip(fields, law.local_term_parts(), strict=True):
        twice_g = np.kron(twice_g, np.eye(9, dtype=int))
        assert not (a.T @ twice_g @ b).any()
        h8 = np.zeros((99, 99), dtype=np.int64)
        h8[:9, :9] = 2 * a.T @ twice_g @ a
        h8[9:, 9:] = b @ twice_g @ b
        assert not h8[:, zero].any()
        np.testing.assert_allclose(h8 / 8, actual, atol=2e-14, rtol=0)
    pulse = law.echo()[:, law.SRC.index((1, 1))]
    assert sum(abs(pulse[zero]) ** 2) > 1 - 2e-15
    charge_read = [sharp_contrast(x, 0) for x in range(3)]
    assert charge_read == [0, 1, 2]
    return dict(
        complete_native_columns=law.DIM,
        blank_columns=len(blank),
        silent_columns=len(zero),
        filled_echo_support=np.flatnonzero(abs(pulse) > 1e-12).tolist(),
        signed_load_against_zero=charge_read,
    )


def path_fixture(board):
    """Deposit just two live edges inside a degree-five triple of the first cell."""
    for hub in range(12):
        if board.degree[hub] != 5:
            continue
        leaves = [q for q in board.neighbours[hub] if q < 12 and board.degree[q] == 5]
        for left, right in combinations(leaves, 2):
            if right not in board.neighbours[left]:
                return [left, hub, right]
    raise ValueError("No degree-five induced path")


def equal_sign_components(board):
    gateways = [p for p in range(board.n) if board.degree[p] == 6]
    graph = nx.Graph()
    graph.add_nodes_from(gateways)
    graph.add_edges_from(
        (a, b) for a, b in board.seams if board.degree[a] == board.degree[b] == 6
    )
    return dict(
        maximum_gateway_degree=max(dict(graph.degree()).values()),
        maximum_gateway_component=max(map(len, nx.connected_components(graph))),
        gateway_graph_is_forest=nx.is_forest(graph),
    )


def local_matrix(board, path, sign, kappa):
    """Assemble actual unmarked seam blocks on the vacancy/excess sector above the host."""
    law = candidate.native_law
    ex, contrast = law.local_term_parts()
    term = ex + kappa * contrast
    i11 = law.SRC.index((1, 1))
    ordered = [law.SRC.index((1 + sign, 1)), law.SRC.index((1, 1 + sign))]
    seam = term[np.ix_(ordered, ordered)] - np.eye(2) * term[i11, i11]
    h = np.diag(kappa * ((1 + sign) ** 2 - 1) * (6 - board.degree[path])).astype(float)
    for a, b in [(0, 1), (1, 2)]:
        if path[b] not in board.neighbours[path[a]]:
            raise ArithmeticError("Source needs live material contacts")
        h[np.ix_([a, b], [a, b])] += seam
    return h


def read_frequencies(signal):
    spacing = TIMES[1] - TIMES[0]
    grid = 2 * pi * np.fft.rfftfreq(len(signal), spacing)
    peak = 1 + np.argmax(abs(np.fft.rfft(signal - np.mean(signal)))[1:])

    def residual(omega):
        design = np.column_stack(
            (np.ones(len(signal)), np.cos(omega * TIMES), np.sin(omega * TIMES))
        )
        weights = np.linalg.lstsq(design, signal, rcond=None)[0]
        return np.linalg.norm(signal - design @ weights) ** 2

    fit = minimize_scalar(
        residual,
        bounds=(grid[max(0, peak - 1)], grid[min(len(grid) - 1, peak + 1)]),
        method="bounded",
        options={"xatol": 1e-11},
    )
    return float(fit.x)


def bank(h, rng):
    energies, basis = np.linalg.eigh(h)
    state = (np.exp(-1j * TIMES[:, None] * energies) * basis[0]) @ basis.T
    probs = abs(state) ** 2
    counts = []
    for prob in probs:
        # These sources have exact support on the path. Simulate the full
        # terminal read keeping the outside outcome and deterministic dark ports.
        p = np.round(prob, 12)
        row = rng.multinomial(SHOTS, p / p.sum())
        counts.append([*row, 0])
    counts = np.array(counts, dtype="<u4")
    fast = read_frequencies(counts[:, 1] / SHOTS)
    slow = read_frequencies((counts[:, 0].astype(float) - counts[:, 2]) / SHOTS)
    return dict(
        host_relative_spectrum=energies.tolist(),
        predicted_fast_gap=float(energies[2] - energies[0]),
        predicted_slow_gap=float(energies[1] - energies[0]),
        fitted_hub_gap=fast,
        fitted_leaf_imbalance_gap=slow,
        counts_shape=list(counts.shape),
        counts_u32_le_zlib_b64=base64.b64encode(
            zlib.compress(counts.tobytes(), 9)
        ).decode(),
    )


def measure():
    k = candidate.KAPPA
    if k != 8:
        raise ValueError("Sharp numeric bounds assume the registered coupling")
    alpha = (19 - sqrt(297)) / 4
    beta = alpha / (1 - 2 * alpha)
    assert abs(8 - 2 * beta - alpha) < 2e-14
    assert all(8 * (d - 1) - d * beta >= alpha - 2e-14 for d in range(2, 7))
    boards = [federation(i) for i in range(3)]
    paths = [path_fixture(b) for b in boards]
    banks = []
    for label, strength in [("native", k), ("half_contrast", k / 2)]:
        for sign in (-1, 1):
            h = local_matrix(boards[0], paths[0], sign, strength)
            for b, p in zip(boards[1:], paths[1:], strict=True):
                np.testing.assert_allclose(
                    local_matrix(b, p, sign, strength), h, atol=1e-13
                )
            banks.append(
                dict(
                    law=label,
                    sign=sign,
                    **bank(
                        h, np.random.default_rng(11500 + 100 * int(strength) + sign)
                    ),
                )
            )
    return dict(
        schema=1,
        source="one signed defect on a unit-filled host with clean private references",
        local=seam_certificate(),
        boards=[
            dict(
                cells=b.cells,
                live_path=p,
                marked_seams=len(b.seams) - 2,
                **equal_sign_components(b),
            )
            for b, p in zip(boards, paths, strict=True)
        ],
        vacancy_positive_floor=alpha,
        vacancy_upper=54,
        positive_vacancy_ratio_ceiling=54 / alpha,
        excess_positive_floor=8,
        excess_upper=69,
        sharp_leaf_payment=beta,
        starts_per_wait=SHOTS,
        waits=dict(count=len(TIMES), step=0.25),
        banks=banks,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    document = json.dumps(measure(), indent=2) + "\n"
    if parser.parse_args().write:
        (HERE / "results.json").write_text(document)
    else:
        print(document, end="")


if __name__ == "__main__":
    main()
