"""Two written clocks with one freely propagating excess on the registered host."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_03_clock.written_coin import measure as clock
from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare
from model.executable_model import federation
from model.experiments.material_inertia.dark_collar import label

HERE = Path(__file__).resolve().parent
ENDPOINTS = (0, 5)
SHOTS, TIMER_SHOTS, TIMES = clock.SHOTS, clock.TIMER_SHOTS, clock.TIMES
# Covers each of the four joint bins, marginals and branch events on every bank.
EPS = float(np.sqrt(np.log(2 * 150000 / 0.01) / (2 * SHOTS)))


def native_blocks(k):
    """Certify the full incoming columns, including every off-sector row."""
    law = candidate.native_law
    ex, co = law.local_term_parts()
    h = ex + k * co
    np.testing.assert_allclose(h[:, :9], np.vstack([law.seam_term(k), np.zeros((90, 9))]), atol=1e-13)
    for parity in (1, -1):
        energy = 4 * k + (1 - parity) / 2
        for coin in (1, -1):
            for loads in ((1, 1), (2, 1), (1, 2)):
                vector = np.kron(label(2, parity, coin), np.eye(9)[law.SRC.index(loads)])
                value = 0 if loads == (1, 1) and coin == 1 else energy
                np.testing.assert_allclose(h @ vector, value * vector, atol=1e-13)
    timer = [law.SRC.index(pair) for pair in ((1, 0), (0, 1))]
    np.testing.assert_allclose(
        h[np.ix_(timer, timer)], (k + .5) * np.eye(2) - .5 * np.array([[0, 1], [1, 0]]),
        atol=1e-13,
    )


def board_generator(board, k):
    # The written seam acts separately; the degree reference retains the
    # original host degrees. Subtract the common fully filled host energy.
    adjacency = board.adjacency().toarray()
    a, b = ENDPOINTS
    assert adjacency[a, b] == 1
    adjacency[a, b] = adjacency[b, a] = 0
    degree = adjacency.sum(axis=1)
    return np.diag((.5 + k) * degree + 3 * k * (6 - board.degree)) - .5 * adjacency


def kernels(h, starts):
    energies, vectors = eigh(h, driver="evd")
    phases = np.exp(-1j * TIMES[:, None] * energies)
    return {site: (phases * vectors[site]) @ vectors.T for site in starts}


def complete_probabilities(sym, anti, near):
    same = abs(sym + anti) ** 2 / 4
    flip = abs(sym - anti) ** 2 / 4
    joint = np.stack([same[:, ~near].sum(1), same[:, near].sum(1),
                      flip[:, ~near].sum(1), flip[:, near].sum(1)], axis=1)
    np.testing.assert_allclose(joint.sum(axis=1), 1, atol=3e-13)
    return joint, [np.sum(abs(v[:, near]) ** 2, axis=1) for v in (sym, anti)]


def enrolled(rng, p):
    return np.array([clock.enroll(rng, row) for row in p])


def fit(timer, flip):
    return clock.fit(np.column_stack([timer, SHOTS - flip, flip, np.zeros(len(flip), dtype=int)]))


def forecast(timer, flips, frequency):
    t = 16 * np.arcsin(np.sqrt(timer / TIMER_SHOTS))
    design = np.column_stack([np.ones(len(t)), np.cos(frequency * t), np.sin(frequency * t)])
    p = flips / SHOTS
    coeff = np.linalg.lstsq(design[:97], p[:97], rcond=None)[0]
    return float(np.sqrt(np.mean((design[97:] @ coeff - p[97:]) ** 2)))


def source_meter(board, launches):
    """Meter the ordinary excess on the filled host before clock installation."""
    law = candidate.native_law
    exchange, _ = law.local_term_parts()
    echo = law.echo()
    rng = np.random.default_rng(27112507)
    records = {}
    probes = {}
    for distance, site in launches.items():
        odd_count = 0
        contrast_sum = 0
        for neighbour in board.neighbours[site]:
            loads = (2, 1) if site < neighbour else (1, 2)
            start = np.eye(law.DIM)[:, law.SRC.index(loads)]
            # On the blank code the exchange term is an orthogonal projector.
            p_odd = np.linalg.norm(exchange @ start) ** 2
            parity_clicks = int(rng.binomial(SHOTS, round(p_odd, 13)))
            written = (echo @ start)[law.N_SRC:].reshape(law.N_LAB, law.N_SRC)
            p_magnitude = [np.sum(abs(written[[abs(d) == m for _, d in law.LABELS]]) ** 2)
                           for m in (0, 1, 2)]
            magnitudes = clock.enroll(rng, p_magnitude)
            records[f"energy_{site}_{neighbour}"] = np.array([parity_clicks, *magnitudes])
            odd_count += parity_clicks
            contrast_sum += magnitudes[1] + 4 * magnitudes[2]
        # Local sharp port reads give 2 in the excess arm, 1 in the matched
        # filled reference. Their squared-load difference is three.
        refs = 6 - int(board.degree[site])
        energy = odd_count / SHOTS + candidate.KAPPA * (contrast_sum / SHOTS + 3 * refs)
        probes[str(distance)] = {"site": site, "degree": int(board.degree[site]),
                                  "excess_site_load": 2, "reference_site_load": 1,
                                  "formation_energy": float(energy)}
    return probes, records


def run():
    board = federation(0)
    distance = shortest_path(board.adjacency(), directed=False, indices=ENDPOINTS, unweighted=True).min(0)
    launches = {int(d): int(np.flatnonzero(distance == d)[0]) for d in (0, 1, 2, 3, distance.max())}
    mask = np.zeros(board.n, dtype=bool)
    mask[list(ENDPOINTS)] = True
    rng = np.random.default_rng(27112506)
    raw = {}
    panels = {}
    for k in (candidate.KAPPA, 0):
        native_blocks(k)
        h = board_generator(board, k)
        common = kernels(h, list(launches.values()))
        for parity in (1, -1):
            gap = 4 * k + (1 - parity) / 2
            moving = kernels(h + np.diag(mask * gap), list(launches.values()))
            phase = np.exp(-1j * gap * TIMES[:, None])
            arms = {"blank": (np.ones_like(phase), phase, np.array([False]))}
            for dist, site in launches.items():
                arms[f"free{dist}"] = (moving[site], common[site] * phase, mask)
                # With the source parked all other port columns are empty.
                # The parked degree reference and open filled contacts give
                # the same global phase on both arms and drop out.
                near_site = mask[site]
                arms[f"park{dist}"] = (phase if near_site else np.ones_like(phase),
                                       phase, np.array([near_site]))
            for name, (sym, anti, near) in arms.items():
                key = f"k{int(k)}_{'even' if parity == 1 else 'odd'}_{name}"
                joint, arrivals = complete_probabilities(sym, anti, near)
                raw[key] = enrolled(rng, joint)
                timer = rng.binomial(TIMER_SHOTS, np.sin(TIMES / 16) ** 2)
                raw[key + "_timer"] = timer
                flips = raw[key][:, 2:].sum(axis=1)
                panels[key] = {"fit": fit(timer, flips)}
                if name.startswith("free"):
                    selected = []
                    for branch, p in zip(("sym", "anti"), arrivals):
                        bank = enrolled(rng, np.column_stack([np.full(len(TIMES), .5),
                                                              (1 - p) / 2, p / 2]))
                        raw[key + "_" + branch] = bank
                        selected.append(bank[:, 2] / bank[:, 1:].sum(axis=1))
                    diff = abs(selected[0] - selected[1])
                    loc = int(diff.argmax())
                    panels[key]["backaction"] = {
                        "appointment": loc, "symmetric_arrival": float(selected[0][loc]),
                        "antisymmetric_arrival": float(selected[1][loc]),
                        "absolute_difference": float(diff[loc]),
                    }
    for name in ["blank"] + [f"{mode}{d}" for mode in ("free", "park") for d in launches]:
        ek, ok = f"k8_even_{name}", f"k8_odd_{name}"
        src = panels[ek]["fit"]
        valid = src["frequency"] is not None and src["held_rms"] <= .025
        predicted = None
        if valid:
            ratio = panels["k8_odd_blank"]["fit"]["frequency"] / panels["k8_even_blank"]["fit"]["frequency"]
            w = src["frequency"] * ratio
            predicted = {"frequency": w, "held_rms": forecast(raw[ok + "_timer"], raw[ok][:, 2:].sum(1), w)}
        panels[ok]["even_clock_forecast"] = predicted
    for parity in ("even", "odd"):
        for k in (8, 0):
            blank = raw[f"k{k}_{parity}_blank"][:, 2:].sum(1) / SHOTS
            for d in launches:
                key = f"k{k}_{parity}_free{d}"
                flip = raw[key][:, 2:].sum(1) / SHOTS
                panels[key]["max_flip_change_from_blank"] = float(max(abs(flip - blank)))
    energies, energy_counts = source_meter(board, launches)
    raw.update(energy_counts)
    result = {
        "source_only_energy_meter": energies, "host_ports": board.n, "endpoints": list(ENDPOINTS),
        "launch_by_native_separation": {str(d): s for d, s in launches.items()},
        "shots": SHOTS, "timer_shots": TIMER_SHOTS,
        "joint_columns": ["same_away", "same_at_clock", "flip_away", "flip_at_clock"],
        "branch_columns": ["branch_rejected", "selected_away", "selected_at_clock"],
        "simultaneous_99_single_probability_radius": EPS, "panels": panels,
    }
    return result, raw


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result, raw = run()
    archive = HERE / "counts.npz"
    summary = HERE / "results.json"
    if args.write:
        np.savez_compressed(archive, **raw)
        summary.write_text(json.dumps(result, indent=2) + "\n")
    else:
        with np.load(archive) as expected:
            assert set(raw) == set(expected.files)
            for key, array in raw.items():
                np.testing.assert_array_equal(array, expected[key])
        compare(result, json.loads(summary.read_text()), "free source clock")
    for key, row in result["panels"].items():
        if key.startswith("k8_") and ("_free" in key or "_blank" in key):
            print(key, row)


if __name__ == "__main__":
    main()
