"""Test the constant-density round radial product on clocked, freely evolving stations."""

from hashlib import sha256
from itertools import combinations, product
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.optimize import least_squares
from scipy.sparse.linalg import expm_multiply

from . import station_chart as station
from .shell_source import Board, check, local_source
from .shells import hands, protocol

HERE = Path(__file__).resolve().parent
BOOK = HERE / "product_protocol.json"


def classify(counts):
    symmetric = counts + counts.T
    values = np.sort(symmetric[np.triu_indices(len(symmetric), 1)])
    gaps = np.diff(values)
    index = int(np.argmax(gaps))
    check(
        gaps[index] > 0 and np.sum(gaps == gaps[index]) == 1,
        "angular classes unresolved",
    )
    incidence = symmetric > values[index]
    np.fill_diagonal(incidence, False)
    return incidence


def angular_read(board, book, rng):
    n = board.car.cells
    start = np.repeat(np.eye(n), 12, axis=0) / np.sqrt(12.0)
    h, _, offset = board.generator("native")
    h = h - offset * sparse.eye(h.shape[0], format="csr")
    counts = []
    for factor in (1, 2):
        action = -1j * factor * book["angular_step"] * h
        wave = expm_multiply(action, start, traceA=action.diagonal().sum())
        probabilities = (abs(wave.reshape(n, 12, n)) ** 2).sum(axis=1).T
        check(
            np.max(abs(probabilities.sum(axis=1) - 1)) < 1e-10, "lost angular records"
        )
        counts.append(
            np.array(
                [
                    rng.multinomial(book["angular_shots"], p / p.sum())
                    for p in probabilities
                ]
            )
        )
    incidence = classify(counts[0])
    check(
        np.array_equal(incidence, classify(counts[1])), "held angular forecast failed"
    )
    laplacian = np.diag(incidence.sum(axis=1).astype(float)) - incidence
    eigen, basis = np.linalg.eigh(laplacian)
    check(eigen[1] > 1e-8, "disconnected angular record")
    band = basis[:, np.isclose(eigen, eigen[1], atol=1e-9, rtol=0)]
    lengths = np.linalg.norm(band, axis=1)
    check(band.shape[1] == 3 and min(lengths) > 1e-8, "angular triplet unresolved")
    directions = band / lengths[:, None]
    # A disconnected generator is block diagonal, so intercell arrivals vanish at every exposure.
    parked, _, _ = board.generator("external_exchange_off")
    nonzero = parked.tocoo()
    check(
        np.all(nonzero.row // 12 == nonzero.col // 12),
        "disconnected control couples cells",
    )
    return incidence, directions, counts, float(eigen[1])


def menu(incidence, home):
    cells = [home, *np.flatnonzero(incidence[home]).tolist()]
    check(len(cells) == 4, "angular anchor count")
    for c in cells[1:4]:
        cells.append(min(set(np.flatnonzero(incidence[c]).tolist()).difference(cells)))
    return cells


def triangle_certificate(angles, paths, readings, tol, bins):
    """Relax each unknown line coordinate difference over an entire screen-radius interval."""
    readings = np.asarray(readings, dtype=float)
    idx = dict(zip(paths, range(len(paths)), strict=True))
    theta = np.array([angles[i, j] for i, j in paths])
    check(np.all(theta > 0) and np.all(readings > 0), "unresolved station pair")
    # Numerical padding applied outwards. Conclusions are conditional on the reconstructed angles.
    pad = 1e-9 * max(1.0, readings.max())
    low, high = (1 - tol) * readings - pad, (1 + tol) * readings + pad
    limit = float(np.min(high / theta)) + pad
    boundaries = np.linspace(0.0, limit, bins + 1)
    low_angular = np.maximum(0.0, boundaries[:-1, None] * theta - pad)
    high_angular = boundaries[1:, None] * theta + pad
    excluded = np.any(low_angular > high, axis=1)
    lower = np.maximum(0.0, np.sqrt(np.maximum(0.0, low**2 - high_angular**2)) - pad)
    upper = np.sqrt(np.maximum(0.0, high**2 - low_angular**2)) + pad
    witnesses = []
    for tri in combinations(range(len(angles)), 3):
        keys = list(combinations(tri, 2))
        if not all(k in idx for k in keys):
            continue
        slots = [idx[k] for k in keys]
        possible = np.zeros(bins, bool)
        # On a real line, one absolute difference is the sum of the other two.
        for k in range(3):
            a, b, c = [slots[(k + j) % 3] for j in range(3)]
            possible |= (lower[:, a] <= upper[:, b] + upper[:, c]) & (
                upper[:, a] >= lower[:, b] + lower[:, c]
            )
        ruled = ~possible & ~excluded
        if np.any(ruled):
            witnesses.append(dict(stations=tri, intervals=int(ruled.sum())))
        excluded |= ~possible
    return dict(
        excludes_every_radius=bool(excluded.all()),
        tested_intervals=bins,
        surviving_intervals=int((~excluded).sum()),
        radius_upper=float(limit),
        padding=float(pad),
        triangle_witnesses=witnesses,
    )


def fit(angles, training):
    """Fit only the supplied training distances. Reflection and radial origin are free."""
    left, right = np.array(station.TRAIN).T
    theta = angles[left, right]
    d = np.asarray(training)
    check(np.all(d > 0) and np.all(theta > 0), "invalid training ranges")

    def residual(x):
        radial = np.r_[0.0, x[1:]]
        return np.hypot(x[0] * theta, radial[left] - radial[right]) / d - 1

    def jacobian(x):
        radial = np.r_[0.0, x[1:]]
        delta = radial[left] - radial[right]
        norm = np.maximum(np.hypot(x[0] * theta, delta), 1e-14) * d
        out = np.zeros((len(d), 8))
        slots = np.arange(len(d))
        out[:, 0] = x[0] * theta**2 / norm
        out[slots, 1 + left] = delta / norm
        out[slots, 1 + right] = -delta / norm
        return out[:, [0, 2, 3, 4, 5, 6, 7]]

    scale = 0.6 * min(d / theta)
    best_cost = np.inf
    for signs in product((-1.0, 1.0), repeat=5):
        initial = np.r_[scale, 0.5 * d.mean() * np.array([1.0, *signs])]
        candidate = least_squares(
            residual,
            initial,
            jac=jacobian,
            bounds=(np.r_[0.0, np.full(6, -np.inf)], np.inf),
        )
        if candidate.cost < best_cost:
            best_cost, best = candidate.cost, candidate.x
    radial = np.r_[0.0, best[1:]]
    predictions = np.hypot(best[0] * angles, radial[:, None] - radial)
    return predictions, float(best[0]), radial.tolist(), float(max(abs(residual(best))))


def score(angles, observation, settings, book):
    d = observation["training_ranges"]
    prediction, scale, radial, placement = fit(angles, d)
    anchors = [d[station.TRAIN.index(pair)] for pair in station.ANCHOR_PAIRS]
    intercept, slope = np.linalg.lstsq(
        np.column_stack((np.ones(6), anchors)),
        observation["training_pulse"],
        rcond=None,
    )[0]
    held = []
    for (i, j), echo, pulse in zip(
        station.HELD_PAIRS,
        observation["held_echo"],
        observation["held_pulse"],
        strict=True,
    ):
        pred = prediction[i, j]
        held.append(
            dict(
                stations=[i, j],
                echo=echo,
                pulse=pulse,
                echo_prediction=float(pred),
                pulse_prediction=float(intercept + slope * pred),
                echo_relative_error=float(abs(pred / echo - 1)),
                pulse_relative_error=float(abs((intercept + slope * pred) / pulse - 1)),
            )
        )
    tol = book["relative_tolerance"]
    phase = observation["stationary_clock_max_error"]
    ok = (
        scale > 0
        and placement <= tol
        and slope > 0
        and phase <= settings["clock_tolerance_radians"]
        and all(
            max(h["echo_relative_error"], h["pulse_relative_error"]) <= tol
            for h in held
        )
    )
    return dict(
        passed=bool(ok),
        fitted_screen_radius=scale,
        radial=radial,
        maximum_training_error=placement,
        stationary_clock_error=phase,
        held=held,
        training_certificate=triangle_certificate(
            angles, station.TRAIN, d, tol, book["radius_bins"]
        ),
        complete_certificate=triangle_certificate(
            angles,
            station.TRAIN + station.HELD_PAIRS,
            d + observation["held_echo"],
            tol,
            book["radius_bins"],
        ),
    )


def control(settings, book):
    vectors = np.array(
        [
            [5, 0, 0],
            [1, 4, 0],
            [2, 0, 3],
            [3, -1, -2],
            [0, 1, 4],
            [-2, 3, 2],
            [3, 2, -1],
        ],
        float,
    )
    vectors /= np.linalg.norm(vectors, axis=1)[:, None]
    angles = np.arccos(np.clip(vectors @ vectors.T, -1.0, 1.0))
    r = np.array([0, 1, 2, -0.8, -1.3, 0.3, 2.8])
    distance = np.hypot(3 * angles, r[:, None] - r)

    def take(pairs):
        return [float(distance[i, j]) for i, j in pairs]

    result = score(
        angles,
        dict(
            training_ranges=take(station.TRAIN),
            held_echo=take(station.HELD_PAIRS),
            training_pulse=[1 + 0.8 * d for d in take(station.ANCHOR_PAIRS)],
            held_pulse=[1 + 0.8 * d for d in take(station.HELD_PAIRS)],
            stationary_clock_max_error=0.0,
        ),
        settings,
        book,
    )
    check(
        result["passed"]
        and not result["complete_certificate"]["excludes_every_radius"],
        "round-product control rejected",
    )
    return result


def run(output):
    book = json.loads(BOOK.read_text())
    settings = json.loads(station.PROTOCOL.read_text())
    matching = protocol("shell_protocol.json")["matching"]
    archived = json.loads((HERE / "station_chart.json").read_text())
    archived_by_cell = {q["cells"]: q for q in archived["boards"]}
    hand_term, _ = local_source(1.0, 8.0)
    signs = hands(
        hand_term, settings["step"] * np.arange(2 * settings["appointments"] + 1)
    )
    result = dict(
        protocol_sha256=sha256(BOOK.read_bytes()).hexdigest(),
        control=control(settings, book),
        boards=[],
    )
    result["source_sha256"] = {
        name: sha256((HERE / name).read_bytes()).hexdigest()
        for name in (
            "radial_product.py",
            "station_chart.py",
            "station_chart_protocol.json",
            "station_chart.json",
            "shell_sources.json",
            "shell_source.py",
            "shells.py",
        )
    }
    records = {}
    for level in settings["levels"]:
        board = Board(level, matching)
        n = board.car.cells
        incidence, directions, counts, eigen = angular_read(
            board, book, np.random.default_rng([book["seed"], level, 0])
        )
        for k, c in enumerate(counts):
            records[f"angular_{level}_{k}"] = c
        menus = []
        for group, home in enumerate((0, int(book["far_marker_fraction"] * n))):
            cells = menu(incidence, home)
            angular = directions[cells]
            angles = np.arccos(np.clip(angular @ angular.T, -1.0, 1.0))
            if group == 0:
                prior = archived_by_cell[n]
                check(
                    cells == prior["station_cells"], "original station reconstruction"
                )
                population, real = prior["population"], prior["acquisitions"]
            else:
                leg, intensity, _ = station.flights(
                    board, [12 * c + 4 for c in cells], settings, "native"
                )
                population = station.acquisition(
                    leg, intensity, signs, settings, None, {}
                )
                real = []
                for replica in range(settings["replicas"]):
                    data = {}
                    real.append(
                        station.acquisition(
                            leg,
                            intensity,
                            signs,
                            settings,
                            [book["seed"], level, group, replica],
                            data,
                        )
                    )
                    records.update(
                        {f"fresh_{level}_{replica}_{k}": v for k, v in data.items()}
                    )
            menus.append(
                dict(
                    stations=cells,
                    angles=angles.tolist(),
                    population=score(angles, population, settings, book),
                    acquisitions=[score(angles, o, settings, book) for o in real],
                )
            )
        result["boards"].append(
            dict(
                cells=n,
                first_angular_eigenvalue=eigen,
                angular_rank=directions.shape[1],
                held_angular_classes_match=True,
                menus=menus,
            )
        )
    output.mkdir(parents=True, exist_ok=True)
    path = output / "radial_product_counts.npz"
    np.savez_compressed(path, **records)
    result["counts_sha256"] = sha256(path.read_bytes()).hexdigest()
    (output / "radial_product.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args().output)
