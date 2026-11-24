"""Terminal scattering of opposite signed probes on the freely arriving bound pair."""

import argparse
import json
from itertools import combinations, product
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.experiments.interval_echo_encounter.conjugation import free_seam
from . import measure as flight

HERE = Path(__file__).resolve().parent
N = 2**40
CAL = 2**35
SIGNS = (-1, 1)
BINS = np.arange(-2, 3)
K = candidate.KAPPA
WORDS = list(product(range(3), repeat=3))


def select_rails(board):
    layouts = []
    for p in range(board.n):
        edges = [
            (q, r)
            for q in board.neighbours[p]
            for r in board.neighbours[q]
            if r != p and board.degree[r] == 5
        ]
        lays = next(c for c in combinations(edges, 3) if len(set(sum(c, ()))) == 6)
        layouts.append(lays)
    return layouts


def instrument(dwell, contact_strength):
    """The effective ternary effects of the actual three-port circuit and signed echo."""
    law = candidate.native_law
    u = law.echo()
    blank = np.eye(law.DIM, law.N_SRC)
    split = np.kron(np.eye(3), expm(-0.5j * np.pi * free_seam(K)))
    recombine = np.kron(np.eye(3), expm(-1.5j * np.pi * free_seam(K)))
    ref = np.diag([np.exp(-1j * np.pi * r**2 / 2) for _, _, r in WORDS])
    contact = np.kron(expm(-1j * dwell * free_seam(contact_strength)), np.eye(3))
    unitary = recombine @ ref @ contact @ split
    native_reads = []
    for d in BINS:
        written = np.repeat([False, *(v == d for _, v in law.LABELS)], law.N_SRC)
        restored = u.T @ (written[:, None] * (u @ blank))
        expected = np.diag([a - b == d for a, b in law.SRC])
        np.testing.assert_allclose(restored, blank @ expected, atol=2e-14)
        native_reads.append(np.kron(np.eye(3), blank.T @ restored))
    curves = []
    for sign in SIGNS:
        columns = unitary[:, [WORDS.index((s, 1 + sign, 1)) for s in range(3)]]
        probs = []
        for effect in native_reads:
            # Diagonality permits the full address sum after the nondemolition availability tests.
            read = columns.conj().T @ effect @ columns
            values = np.diag(read).real
            np.testing.assert_allclose(read, np.diag(values), atol=1e-13)
            probs.append(values)
        curves.append(np.array(probs).T)
    curves = np.array(curves)
    np.testing.assert_allclose(curves.sum(axis=2), 1, atol=2e-13)
    return curves


def pathways(board, layouts, x, y, weights):
    """For each port class combine all geometry-selected detector addresses and sign histories."""
    panels = {}
    for degree in (5, 6):
        pool = np.flatnonzero(board.degree == degree)
        amounts = np.zeros((250, 3))
        for p in pool:
            state = 1 + (x == p).astype(int) - (y == p)
            waiting = np.ones(len(weights), bool)
            history = np.ones(len(weights), int)
            for q, r in layouts[p]:
                sign = (x == q).astype(int) - (y == q) - (x == r) + (y == r)
                history[waiting] = 5 * history[waiting] + sign[waiting] + 2
                valid = waiting & (sign == 0)
                keys = 3 * history[valid] + state[valid]
                amounts += np.bincount(
                    keys, weights[valid], minlength=amounts.size
                ).reshape(amounts.shape)
                waiting &= ~valid
            if waiting.any():
                raise ValueError("An enrolled query has no host-only probe rails")
        amounts /= len(pool)
        np.testing.assert_allclose(amounts.sum(), 1, atol=5e-13)
        panels[degree] = amounts
    return panels


def sample(shots, probs, rng):
    probs = np.maximum(np.asarray(probs), 0)
    np.testing.assert_allclose(probs.sum(), 1, rtol=0, atol=2e-12)
    probs = probs.round(14)
    return rng.multinomial(shots, (probs / probs.sum()).ravel()).reshape(probs.shape)


def shots_panel(amounts, curves, admission, arrival, rng):
    rows = []
    used = np.flatnonzero(amounts.sum(axis=1))
    for curve in curves:
        probs = amounts[used] @ curve
        prep, transit, accepted = sample(
            N, [1 - admission, admission * (1 - arrival), admission * arrival], rng
        )
        counts = sample(accepted, probs, rng)
        rows.append(
            {
                "filter_rejected": int(prep),
                "no_arrival": int(transit),
                "availability_codes": used.tolist(),
                "terminal_counts": counts.tolist(),
            }
        )
    return rows


def observables(counts):
    means, variances = [], []
    for sign, bins in zip(SIGNS, counts, strict=True):
        bins = np.array(bins)
        prob = bins / bins.sum()
        score = (BINS == sign).astype(int) - (BINS == -sign)
        mean = prob @ score
        means.append(mean)
        variances.append((prob @ score**2 - mean**2) / bins.sum())
    minus, plus = means
    return np.array([(plus - minus) / 2, (plus + minus) / 2]), sum(variances) / 4


def decode(data, calibrations, scale_calibrations):
    cal = [observables(c) for c in calibrations]
    blank, var_blank = cal[1]
    unit = [observables(c) for c in scale_calibrations]
    scale = (unit[2][0][0] - unit[0][0][0]) / 2
    var_scale = (unit[2][1] + unit[0][1]) / 4
    panel = {}
    acc, var = np.zeros(2), 0.0
    total = 0
    for deg, settings in data.items():
        size = settings["size"]
        mean, v = observables(
            [np.sum(r["terminal_counts"], axis=0) for r in settings["probes"]]
        )
        result = size * (mean - blank) / scale
        se = np.sqrt(size**2 * (v + var_blank) + result**2 * var_scale) / abs(scale)
        panel[str(deg)] = {
            "odd_even_in_probe_units": result.tolist(),
            "sampling_se": se.tolist(),
        }
        acc += size * mean
        var += size**2 * v
        total += size
    net = (acc - total * blank) / scale
    se = np.sqrt(var + total**2 * var_blank + net**2 * var_scale) / abs(scale)
    return {
        "classes": panel,
        "whole_odd_even": net.tolist(),
        "sampling_se": se.tolist(),
        "single_excess_odd_scale": float(scale),
    }


def sources(board, layouts):
    x, y, _, _, h, cost = flight.defect.generators(board)
    reduced = flight.PairReduction(h, cost)
    energies, basis, residual, _ = reduced.lowest(h, 24)
    prior = json.loads((HERE / "results.json").read_text())
    np.testing.assert_allclose(energies, prior["lines"], rtol=0, atol=3e-12)
    graph = nx.Graph((a // 12, b // 12) for a, b in board.external)
    prepared = []
    for row in prior["flights"]:
        if row["admission"] < 0.001:
            continue
        sharp = flight.defect.ordered_index(*row["ports"], board.n)
        coeff = basis[sharp]
        admission = coeff @ coeff
        state = basis @ (
            np.exp(-1j * (energies - energies[0]) * 1e6) * coeff / np.sqrt(admission)
        )
        dist = nx.single_source_shortest_path_length(graph, row["owner"])
        away = np.array([dist[i] >= 2 for i in range(board.cells)])
        flag = (x // 12 == y // 12) & away[x // 12]
        mass = abs(state[flag]) ** 2
        arrival = mass.sum()
        np.testing.assert_allclose(
            arrival,
            row["appointments"][-1]["shared_owner_two_or_more"],
            rtol=0,
            atol=2e-8,
        )
        prepared.append(
            (
                f"pair_from_{row['owner']}",
                admission,
                arrival,
                pathways(board, layouts, x[flag], y[flag], mass / arrival),
            )
        )
    for name, xp, yp in (
        ("single_excess", 0, -1),
        ("single_vacancy", -1, 1),
        ("host", -1, -1),
    ):
        prepared.append(
            (
                name,
                1.0,
                1.0,
                pathways(
                    board, layouts, np.array([xp]), np.array([yp]), np.array([1.0])
                ),
            )
        )
    return prepared, float(residual)


def run():
    board = flight.federation(0)
    layouts = select_rails(board)
    prepared, residual = sources(board, layouts)
    rng = np.random.default_rng(381040)
    panels = []
    for strength, dwell in ((K, 1 / 32), (K, 1 / 64), (0.0, 1 / 32)):
        curves = instrument(dwell, strength)
        cal = [
            [sample(CAL, curves[i, s], rng).tolist() for i in range(2)]
            for s in range(3)
        ]
        rows = []
        for name, admission, arrival, amounts in prepared:
            records = {
                deg: {
                    "size": int(sum(board.degree == deg)),
                    "probes": shots_panel(a, curves, admission, arrival, rng),
                }
                for deg, a in amounts.items()
            }
            # Hold the original single-charge unit fixed in the changed-contact panel.
            unit = panels[0]["calibration_counts"] if strength == 0 else cal
            scale = (observables(cal[2])[0][0] - observables(cal[0])[0][0]) / 2
            rows.append(
                {
                    "source": name,
                    "filter_probability": admission,
                    "arrival_given_filter": arrival,
                    "records": records,
                    "read": decode(records, cal, unit),
                }
            )
        panels.append(
            {
                "strength_at_measured_contact": strength,
                "dwell": dwell,
                "calibration_counts": cal,
                "calibrated_odd_scale": float(scale),
                "sharp_diagonal_energies": np.diag(free_seam(strength))
                .reshape(3, 3)
                .tolist(),
                "panels": rows,
            }
        )
    return {
        "source_file": "results.json",
        "seed": 381040,
        "full_generator_residual": residual,
        "enrollments_per_probe_sign_and_degree": N,
        "calibration_enrollments": CAL,
        "availability_history": "start at 1, append signed contrast d by multiplying by 5 and adding d+2; zero contrast opens probe",
        "signed_bins": BINS.tolist(),
        "probe_signs": list(SIGNS),
        "studies": panels,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(
            json.loads(json.dumps(result)),
            json.loads((HERE / "probe.json").read_text()),
            "pair_probe",
        )
        print("Reproduced both signed probes on the freely arriving bound pair.")
