"""Calibrate a source-free three-dimensional chart and predict independent flights."""

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy import sparse
from scipy.linalg import eigh
from scipy.sparse.linalg import expm_multiply

from .shell_source import Board, check, local_source
from .shells import hands, protocol, read_hand

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "station_chart_protocol.json"
PIN = "35ce62d3cd940c53071d79fe4338d921ede647261743a01213c47140564942e2"
ANCHORS = tuple(range(4))
HELD = (4, 5, 6)
ANCHOR_PAIRS = list(combinations(ANCHORS, 2))
HELD_PAIRS = list(combinations(HELD, 2))
ALL_PAIRS = list(combinations(range(7), 2))
TRAIN = [p for p in ALL_PAIRS if p not in HELD_PAIRS]


def station_cells(board):
    neighbours = [set() for _ in range(board.car.cells)]
    for a, b in board.car.external:
        left, right = a // 12, b // 12
        neighbours[left].add(right)
        neighbours[right].add(left)
    cells = [0, *sorted(neighbours[0])]
    for anchor in cells[1:4]:
        cells.append(min(neighbours[anchor].difference(cells)))
    check(len(cells) == len(set(cells)) == 7, "station placements must be distinct")
    return cells


def flights(board, ports, spec, member):
    h, _, offset = board.generator(member)
    h = h - offset * sparse.eye(h.shape[0], format="csr")
    initial = np.zeros((h.shape[0], len(ports)))
    initial[ports, np.arange(len(ports))] = 1
    n, dt = spec["appointments"], spec["step"]
    u = expm_multiply(
        -1j * h,
        initial,
        start=0,
        stop=n * dt,
        num=n + 1,
        traceA=-1j * h.diagonal().sum(),
    )[:, ports, :]
    intensity = abs(u[1:].transpose(2, 1, 0)) ** 2
    leg = np.zeros((7, 7, n + 1))
    for target in range(7):
        # After an unsuccessful read, subtract the amplitude that reached the target.
        first = np.zeros((n + 1, 7), complex)
        for k in range(1, n + 1):
            first[k] = u[k, target] - u[1:k, target, target] @ first[k - 1 : 0 : -1]
        leg[:, target, :-1] = abs(first[1:].T) ** 2
    leg[:, :, -1] = 1 - leg[:, :, :-1].sum(axis=2)
    check(np.min(leg) >= 0 and np.max(intensity) <= 1 + 1e-12, "invalid flight law")
    error = None
    if board.level == 0 and member == "native":
        values, modes = eigh(h.toarray())
        unitary = (modes * np.exp(-1j * dt * values)) @ modes.T
        error = 0.0
        for j, port in enumerate(ports):
            survivor, free = initial.copy(), initial.copy()
            for k in range(n):
                survivor, free = unitary @ survivor, unitary @ free
                error = max(
                    error,
                    np.max(abs(abs(survivor[port]) ** 2 - leg[:, j, k])),
                    np.max(abs(abs(free[port]) ** 2 - intensity[:, j, k])),
                )
                survivor[port] = 0
            error = max(
                error, np.max(abs((abs(survivor) ** 2).sum(axis=0) - leg[:, j, -1]))
            )
        check(error < 1e-10, "independent native propagators disagree")
    return leg, intensity, error


def normalized_mean(weights, times):
    total = np.sum(weights)
    return float(weights @ times / total) if total else float("nan")


def chart(ranges, anchor_pulse, held_echo, held_pulse, phase_error, spec):
    tol = spec["relative_tolerance"]
    measurements = np.array([*ranges, *anchor_pulse, *held_echo, *held_pulse])
    if (
        not np.isfinite(measurements).all()
        or min(measurements) <= 0
        or not np.isfinite(phase_error)
    ):
        return dict(passed=False, status="unresolved flight or clock")
    d = np.full((7, 7), np.nan)
    np.fill_diagonal(d, 0)
    for (i, j), value in zip(TRAIN, ranges, strict=True):
        d[i, j] = d[j, i] = value
    radii = d[0, 1:4] ** 2
    gram = (radii[:, None] + radii - d[1:4, 1:4] ** 2) / 2
    eigen = np.linalg.eigvalsh(gram)
    out = dict(
        anchor_gram_eigenvalues=eigen.tolist(), stationary_clock_max_error=phase_error
    )
    if eigen[0] <= spec["minimum_anchor_eigenvalue_ratio"] * eigen[-1]:
        return dict(
            **out,
            passed=False,
            status="invalid anchor tetrahedron",
            observed_held_echo=held_echo,
            observed_held_pulse=held_pulse,
        )
    basis = np.linalg.cholesky(gram)
    coords = np.vstack((np.zeros(3), basis, np.zeros((3, 3))))
    residuals = []
    for i in HELD:
        rhs = (d[i, 0] ** 2 + radii - d[i, 1:4] ** 2) / 2
        coords[i] = np.linalg.solve(basis, rhs)
        inferred = np.linalg.norm(coords[i] - coords[:4], axis=1)
        residuals.append((inferred / d[i, :4] - 1).tolist())
    design = np.array([[1, d[i, j]] for i, j in ANCHOR_PAIRS])
    intercept, slope = np.linalg.lstsq(design, anchor_pulse, rcond=None)[0]
    rows = []
    for (i, j), echo, pulse in zip(HELD_PAIRS, held_echo, held_pulse, strict=True):
        prediction = np.linalg.norm(coords[i] - coords[j])
        pulse_prediction = intercept + slope * prediction
        rows.append(
            dict(
                stations=[i, j],
                echo=float(echo),
                echo_prediction=float(prediction),
                echo_relative_error=float(abs(prediction / echo - 1)),
                pulse=float(pulse),
                pulse_prediction=float(pulse_prediction),
                pulse_relative_error=float(abs(pulse_prediction / pulse - 1)),
            )
        )
    radius_error = float(np.max(np.abs(residuals)))
    ok = (
        radius_error <= tol
        and slope > 0
        and phase_error <= spec["clock_tolerance_radians"]
        and all(
            max(row["echo_relative_error"], row["pulse_relative_error"]) <= tol
            for row in rows
        )
    )
    return dict(
        **out,
        coordinates=coords.tolist(),
        held_anchor_fractional_residuals=residuals,
        max_placement_error=radius_error,
        pulse_offset=float(intercept),
        pulse_slope=float(slope),
        held=rows,
        passed=bool(ok),
        status="evaluated",
    )


def rank_bound(ranges, echo):
    """Lower bound on the relative distortion of any Euclidean three-space fit."""
    d = np.zeros((7, 7))
    for (i, j), value in zip(TRAIN + HELD_PAIRS, [*ranges, *echo], strict=True):
        d[i, j] = d[j, i] = value
    largest = 0.0
    for subset in combinations(range(7), 5):
        squared = d[np.ix_(subset, subset)] ** 2
        center = np.eye(5) - np.ones((5, 5)) / 5
        eig = np.linalg.eigvalsh(-center @ squared @ center / 2)
        # A 3D Gram is PSD of rank <=3. Entrywise squared-distance error is
        # bounded by (2e+e^2)D^2; centering is a contraction in operator norm.
        obstruction = max(0.0, eig[-4], -eig[0])
        largest = max(
            largest, np.sqrt(1 + 2 * obstruction / np.linalg.norm(squared, 2)) - 1
        )
    return float(largest)


def corner_certificates(population, tolerance, length):
    """Check all five-station subsets after the held experiment is scored."""
    from ..interval_echo_einstein_station.refinement import angular_comparison_bound
    from ..interval_echo_einstein_station.volume import five_site_rank_bound

    distances = np.zeros((7, 7))
    readings = population["training_ranges"] + population["held_echo"]
    for (i, j), reading in zip(TRAIN + HELD_PAIRS, readings, strict=True):
        distances[i, j] = distances[j, i] = reading
    witnesses = []
    for subset in combinations(range(7), 5):
        # The archived box verifier reads positions 0,1,2,4,5 of a six-site menu.
        other = min(set(range(7)).difference(subset))
        order = [*subset[:3], other, *subset[3:]]
        panel = distances[np.ix_(order, order)]
        receipt = five_site_rank_bound(panel, panel, tolerance)
        normalized = distances[np.ix_(subset, subset)] / length
        curved_exact = float(
            angular_comparison_bound(normalized, np.zeros_like(normalized))
        )
        curved_tolerance = float(
            angular_comparison_bound(normalized, tolerance * normalized)
        )
        witnesses.append(
            dict(
                stations=list(subset),
                box=receipt,
                curved_exact_gram_floor=curved_exact,
                curved_tolerance_gram_floor=curved_tolerance,
            )
        )
    return dict(
        tolerance=tolerance,
        five_station_subsets=witnesses,
        excluding_subsets=sum(w["box"]["excludes_three_dimensions"] for w in witnesses),
        curvature_comparison_length=length,
        curved_exact_excluding_subsets=sum(
            w["curved_exact_gram_floor"] > 1e-9 for w in witnesses
        ),
        curved_tolerance_excluding_subsets=sum(
            w["curved_tolerance_gram_floor"] > 1e-9 for w in witnesses
        ),
    )


def control(spec):
    p = np.array(
        [[0, 0, 0], [3, 0, 0], [0, 4, 0], [0, 0, 5], [1, 2, 1], [-2, 1, 3], [2, -1, 1]]
    )
    d = np.linalg.norm(p[:, None] - p, axis=2)
    ranges = [d[i, j] for i, j in TRAIN]
    echo = [d[i, j] for i, j in HELD_PAIRS]
    result = chart(
        ranges,
        [1 + 0.8 * d[i, j] for i, j in ANCHOR_PAIRS],
        echo,
        [1 + 0.8 * x for x in echo],
        0.0,
        spec,
    )
    check(
        result["passed"] and rank_bound(ranges, echo) < 1e-12,
        "flat chart control failed",
    )
    return result


def acquisition(leg, intensity, signs, spec, seed, archive):
    """Population readouts for seed=None, otherwise disjoint finite ensembles."""
    n, shots = spec["appointments"], spec["shots"]
    rng = None if seed is None else np.random.default_rng(seed)

    def multi(total, probabilities):
        return rng.multinomial(total, probabilities / probabilities.sum())

    def standards():
        if rng is None:
            return np.broadcast_to(signs, (7, *signs.shape)).copy()
        return np.array([[multi(shots, prob) for prob in signs] for _ in range(7)])

    cal, fresh = standards(), standards()
    archive.update(calibration_hands=cal, independent_hands=fresh)
    clock = np.array([read_hand(row) for row in cal])
    anchor_phase = clock[:4].mean(axis=0)
    fresh_phase = np.array([read_hand(row) for row in fresh])
    phase_error = float(np.max(abs(fresh_phase[4:] - anchor_phase)))
    ranges = []
    for i, j in TRAIN:
        readings = []
        for a, b in ((i, j), (j, i)):
            counts = leg[a, b] if rng is None else multi(shots, leg[a, b])
            if rng is not None:
                check(counts.sum() == shots, "lost first-hit or timeout records")
            archive[f"first_{a}_{b}"] = counts
            readings.append(
                normalized_mean(counts[:-1], clock[b, 1 : n + 1] - clock[a, 0])
            )
        ranges.append(float(np.mean(readings)))

    def pulse(i, j):
        means = []
        for a, b in ((i, j), (j, i)):
            hit = intensity[a, b]
            if rng is None:
                yes, total_hands = hit, signs[1 : n + 1]
            else:
                panels = np.array(
                    [
                        multi(shots, np.outer([p, 1 - p], hand).ravel()).reshape(2, 4)
                        for p, hand in zip(hit, signs[1 : n + 1], strict=True)
                    ]
                )
                check(
                    np.all(panels.sum(axis=(1, 2)) == shots),
                    "lost pulse or quiet records",
                )
                archive[f"free_{a}_{b}"] = panels
                yes = panels[:, 0].sum(axis=1)
                total_hands = panels.sum(axis=1)
            phase = read_hand(total_hands)
            means.append(normalized_mean(yes, phase - fresh_phase[a, 0]))
        return float(np.mean(means))

    anchor_pulse = [pulse(i, j) for i, j in ANCHOR_PAIRS]
    held_pulse = [pulse(i, j) for i, j in HELD_PAIRS]
    echoes = []
    for i, j in HELD_PAIRS:
        timings = []
        for a, b in ((i, j), (j, i)):
            outward, back = leg[a, b], leg[b, a]
            duration = np.convolve(outward[:-1], back[:-1])
            success = duration.sum()
            fail_out = outward[-1]
            fail_back = outward[:-1].sum() * back[-1]
            np.testing.assert_allclose(success + fail_out + fail_back, 1, atol=1e-13)
            if rng is None:
                counts, terminal = duration, fresh[a, 2:]
            else:
                quota = spec["return_quota"]
                failures = int(rng.negative_binomial(quota, success))
                missed = int(rng.binomial(failures, fail_out / (fail_out + fail_back)))
                counts = multi(quota, duration)
                panels = np.array(
                    [multi(int(k), p) for k, p in zip(counts, signs[2:], strict=True)]
                )
                check(panels.sum() == quota, "incomplete returned-hand enrollment")
                archive[f"echo_hands_{a}_{b}"] = panels
                archive[f"echo_timeouts_{a}_{b}"] = np.array(
                    [missed, failures - missed]
                )
                # Matched reference hands cover rare duration bins. Factorization
                # of the position and hand channels is checked by local_source.
                terminal = fresh[a, 2:] + panels
            timings.append(
                normalized_mean(counts, (read_hand(terminal) - fresh_phase[a, 0]) / 2)
            )
        echoes.append(float(np.mean(timings)))
    result = chart(ranges, anchor_pulse, echoes, held_pulse, phase_error, spec)
    result.update(
        training_ranges=ranges,
        training_pulse=anchor_pulse,
        held_echo=echoes,
        held_pulse=held_pulse,
        euclidean_distortion_lower_bound=rank_bound(ranges, echoes),
    )
    return result


def run(output):
    check(sha256(PROTOCOL.read_bytes()).hexdigest() == PIN, "station protocol changed")
    spec = json.loads(PROTOCOL.read_text())
    matching = protocol("shell_protocol.json")["matching"]
    sign_term, _ = local_source(1.0, 8.0)
    signs = hands(sign_term, spec["step"] * np.arange(2 * spec["appointments"] + 1))
    result = dict(
        protocol_sha256=PIN,
        flat_control=control(spec),
        boards=[],
        python=platform.python_version(),
        numpy=np.__version__,
        scipy=scipy.__version__,
        source_sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                Path(__file__),
                HERE / "shells.py",
                HERE / "shell_source.py",
                HERE / "shell_protocol.json",
                HERE.parent / "interval_echo_einstein_station" / "volume.py",
                HERE.parent / "interval_echo_einstein_station" / "refinement.py",
                HERE.parents[1] / "candidate_laws" / "interval_echo" / "candidate.py",
                HERE.parents[1]
                / "candidate_laws"
                / "contrast_memory_spherical_v1"
                / "law.py",
                HERE.parents[1] / "executable_model" / "carrier.py",
            )
        },
    )
    records = {}
    for level in spec["levels"]:
        board = Board(level, matching)
        cells = station_cells(board)
        ports = [12 * c + 4 for c in cells]
        leg, intensity, error = flights(board, ports, spec, "native")
        parked, quiet, _ = flights(board, ports, spec, "external_exchange_off")
        cross = ~np.eye(7, dtype=bool)
        check(
            np.max(parked[cross, :-1]) == np.max(quiet[cross]) == 0,
            "disconnected control reached another station",
        )
        exact = acquisition(leg, intensity, signs, spec, None, {})
        rows = []
        for rep in range(spec["replicas"]):
            archive = {}
            seed = [spec["seed"], level, rep]
            rows.append(acquisition(leg, intensity, signs, spec, seed, archive))
            for key, value in archive.items():
                records[f"l{level}_r{rep}_{key}"] = value
        result["boards"].append(
            dict(
                cells=board.car.cells,
                station_cells=cells,
                local_port=4,
                independent_propagator_error=error,
                exchange_off_max_interstation_probability=0.0,
                population=exact,
                acquisitions=rows,
                post_run_corner_check=corner_certificates(
                    exact,
                    spec["relative_tolerance"],
                    spec["step"] * spec["appointments"],
                ),
            )
        )
    result["source_free_pass"] = all(
        row["passed"]
        for b in result["boards"]
        for row in [b["population"], *b["acquisitions"]]
    )
    output.mkdir(parents=True, exist_ok=True)
    count_path = output / "station_chart_counts.npz"
    np.savez_compressed(count_path, **records)
    result["counts_sha256"] = sha256(count_path.read_bytes()).hexdigest()
    (output / "station_chart.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args().output)
