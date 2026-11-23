"""Predict the full local instrument for a freely released joint phase."""

import argparse
import json
from itertools import permutations
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_04.__main__ import evolve
from model.executable_model import federation
from model.executable_model.carrier import CELL_SEAMS, seam_digest
from tools.strict_json import loads

HERE = Path(__file__).resolve().parent
PREDICTIONS = HERE / "pair_phase_predictions.json"
SOURCES = (0, 2, 3)
TIMES = (0, 16, 64, 256)
TOL = 2e-9
OUTCOMES = [
    "source refused",
    "first rail refused",
    "second rail refused",
    "unequal source-side ports",
    "equal",
    "unequal",
]
LAWS = ("native", "exchange", "independent")


def rails(cell):
    contacts = set(CELL_SEAMS)
    ports = next(
        (u, v, w, x)
        for u, v, w, x in permutations(range(12), 4)
        if all(tuple(sorted(e)) in contacts for e in ((u, v), (w, x), (u, w), (v, x)))
    )
    return tuple(12 * cell + p for p in ports)


def receivers(board, source):
    others = {p // 12 for edge in board.external if source in {q // 12 for q in edge} for p in edge}
    return [source, *sorted(others - {source})]


def apparatus():
    """Obtain source amplitudes and the two analyzer rows from the native blank block."""
    law = candidate.native_law
    echo = law.echo()
    blank = np.eye(law.DIM)[:, : law.N_SRC]
    for magnitude in (0, 1):
        flag = np.repeat([False, *(abs(d) == magnitude for _, d in law.LABELS)], law.N_SRC)
        observed = echo.T @ (flag[:, None] * (echo @ blank))
        expected = blank * [abs(a - b) == magnitude for a, b in law.SRC]
        np.testing.assert_allclose(observed, expected, rtol=0, atol=1e-13)
    # In each admitted rail the first basis state occupies its first port.
    columns = [law.SRC.index(word) for word in ((1, 0), (0, 1))]
    exchange, contrast = law.local_term_parts()

    def splitter(duration):
        local = expm(-1j * duration * (exchange + candidate.KAPPA * contrast)) @ blank[:, columns]
        matrix = local[columns]
        np.testing.assert_allclose(local, blank[:, columns] @ matrix, rtol=0, atol=1e-12)
        return matrix

    half = splitter(np.pi / 2)
    inputs = []
    for second, sign in ((splitter(3 * np.pi / 2), 1), (half, -1)):
        prepared = np.kron(half[:, 0], second[:, 0])
        source = prepared[[0, 3]]
        np.testing.assert_allclose(abs(prepared) ** 2, 0.25, rtol=0, atol=1e-13)
        np.testing.assert_allclose(source[1], sign * source[0], rtol=0, atol=1e-13)
        inputs.append(source)
    analyzer = np.kron(half, half)[:, [0, 3]]
    return np.array(inputs).T, analyzer


def independent_flights(sector, ports, times):
    """Symmetric product evolution, with the full doubled-port sector retained."""
    values, modes = eigh(sector.car.laplacian().toarray() / 2)
    for time in times:
        one = (modes * np.exp(-1j * time * values)) @ modes[list(ports)].T
        pairs, doubles = [], []
        for i, j in ((0, 2), (1, 3)):
            pairs.append(
                one[sector.aa, i] * one[sector.bb, j] + one[sector.bb, i] * one[sector.aa, j]
            )
            doubles.append(np.sqrt(2) * one[:, i] * one[:, j])
        yield np.array(pairs).T, np.array(doubles).T


def probabilities(bins):
    result = np.asarray(bins, dtype=float)
    if (
        not np.isfinite(result).all()
        or result.min() < -1e-13
        or np.max(abs(result.sum(axis=0) - 1)) > 2e-10
    ):
        raise ArithmeticError("Incomplete probability distribution")
    return np.maximum(result, 0).T.round(14).tolist()


def observe(sector, receiver, states, double_weights, analyzer):
    a, b, c, d = rails(receiver)
    # One unit on each rail is necessary. A double occupation fails the first query.
    in_first = np.isin(sector.aa, (a, b)).astype(int) + np.isin(sector.bb, (a, b)) == 1
    in_second = np.isin(sector.aa, (c, d)).astype(int) + np.isin(sector.bb, (c, d)) == 1
    equal_ids = [sector.index(a, c), sector.index(b, d)]
    crossed_ids = [sector.index(a, d), sector.index(b, c)]
    weight = abs(states) ** 2
    mixed = abs(analyzer @ states[equal_ids]) ** 2
    return probabilities(
        [
            np.full(2, 0.5),
            weight[~in_first].sum(axis=0) + double_weights,
            weight[in_first & ~in_second].sum(axis=0),
            weight[crossed_ids].sum(axis=0),
            mixed[[0, 3]].sum(axis=0),
            mixed[[1, 2]].sum(axis=0),
        ]
    )


def forecast():
    board = federation(0)
    pairs = PairSector(board)
    inputs, analyzer = apparatus()
    result = dict(
        board=seam_digest(board),
        ports=list(rails(0)),
        signs=["plus", "minus"],
        outcomes=OUTCOMES,
        association_outcomes=["source refused", "adjacent", "separated", "double"],
        association=[],
        receiver=[],
        all_time_adjacent_given_launch={},
    )
    for source in SOURCES:
        u, v, w, x = ports = rails(source)
        ids = [pairs.index(u, w), pairs.index(v, x)]
        basis = np.zeros((pairs.P, 2))
        basis[ids, np.arange(2)] = 1
        # For the accepted state, E = <K> - 2*kappa initially and K is positive.
        kinetic = basis.T @ (pairs.exchange() @ basis)
        mean = np.diag(inputs.conj().T @ kinetic @ inputs).real / 0.5
        result["all_time_adjacent_given_launch"][str(source)] = (
            (1 - mean / (2 * candidate.KAPPA)).round(14).tolist()
        )
        for law in LAWS:
            if law == "independent":
                flights = independent_flights(pairs, ports, TIMES)
            else:
                strength = candidate.KAPPA if law == "native" else 0
                states = evolve(pairs.hamiltonian(strength), strength, basis, TIMES)
                flights = ((state, np.zeros((0, 2))) for state in states)
            for t, (columns, doubles) in zip(TIMES, flights, strict=True):
                gram = columns.conj().T @ columns + doubles.conj().T @ doubles
                np.testing.assert_allclose(gram, np.eye(2), rtol=0, atol=2e-10)
                states = columns @ inputs
                double_weights = (abs(doubles @ inputs) ** 2).sum(axis=0)
                weight = abs(states) ** 2
                global_bins = probabilities(
                    [
                        np.full(2, 0.5),
                        weight[pairs.adjacent].sum(axis=0),
                        weight[~pairs.adjacent].sum(axis=0),
                        double_weights,
                    ]
                )
                result["association"].append([source, law, t, *global_bins])
                for dest in receivers(board, source):
                    recorded = observe(pairs, dest, states, double_weights, analyzer)
                    result["receiver"].append([source, dest, law, t, *recorded])
    return result


def verify_forecast(actual, expected):
    if set(actual) != set(expected):
        raise ValueError("Prediction schema differs")
    for field in actual.keys() - {"association", "receiver"}:
        if actual[field] != expected[field]:
            raise ValueError(f"Wrong {field}")
    for field in ("association", "receiver"):
        for row, stored in zip(actual[field], expected[field], strict=True):
            if json.dumps(row[:-2]) != json.dumps(stored[:-2]):
                raise ValueError("Wrong appointment")
            a, b = np.asarray(row[-2:]), np.asarray(stored[-2:])
            if a.shape != b.shape or b.dtype.kind != "f":
                raise ValueError("Wrong bin layout")
            np.testing.assert_allclose(a, b, rtol=0, atol=TOL, equal_nan=False)


def store_result(result, path):
    # Keep one setting per line so each probability and its address can be reviewed together.
    fields = []
    for key, value in result.items():
        if key in ("association", "receiver"):
            body = "[\n" + ",\n".join(json.dumps(v, allow_nan=False) for v in value) + "\n]"
        else:
            body = json.dumps(value, allow_nan=False)
        fields.append(json.dumps(key) + ": " + body)
    with path.open("x", encoding="utf-8") as stream:
        stream.write("{\n" + ",\n".join(fields) + "\n}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write predictions to a new file")
    args = parser.parse_args()
    result = forecast()
    if args.output:
        store_result(result, args.output)
    else:
        verify_forecast(result, loads(PREDICTIONS.read_text(encoding="utf-8")))
    print(len(result["receiver"]), "receiver appointments")
