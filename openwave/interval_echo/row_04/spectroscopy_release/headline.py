"""Replay free signed-history survival from three public binary inputs."""

import argparse
import hashlib
import json
from itertools import combinations
from math import asin, expm1, log, log1p, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import expm


BOARD = "openwave/interval_echo/row_04/FEDERATION_LEVEL0.npz"
SOURCE = "model/candidate_laws/interval_echo/openwave"
WORDS = f"{SOURCE}/row_30_selective/analyzers.json"
CLOCK = f"{SOURCE}/row_04_spectral/clock.json"
DIGESTS = {
    BOARD: "99870398e3d516f899c600e9bf25ba53c7f513142ba8facbe1cedd94b23f915b",
    WORDS: "bca73380e3f4960b8b4aa1268023573845e9f40463c8545af35ed550e7f17a69",
    CLOCK: "0259cd71489f34564d6ef4e2edb91b8f420de947911ebfe35ec733f7336a398d",
}


def physical_histories(board):
    edges = [tuple(sorted(map(int, row))) for row in board["edges"]]
    adjacent = [set() for _ in board["degree"]]
    for left, right in edges:
        adjacent[left].add(right)
        adjacent[right].add(left)
    degree = np.array([len(ports) for ports in adjacent])
    np.testing.assert_array_equal(degree, board["degree"])
    paths = [(0,)]
    joins = []
    for i, path in enumerate(paths):
        for port in sorted(adjacent[path[-1]]):
            if port > path[-1]:
                joins.append((i, len(paths)))
                paths.append((*path, port))
    a = np.zeros((len(paths), len(paths)), dtype=int)
    for parent, child in joins:
        a[parent, child] = a[child, parent] = 1
    ends = np.array([path[-1] for path in paths])
    graph_degree = np.diag(a.sum(axis=1))
    return paths, ends, graph_degree - a, 2 * np.diag(6 - degree[ends]) - graph_degree - a


def eigenspaces(eigenvalues):
    return np.split(np.arange(len(eigenvalues)), np.flatnonzero(np.diff(eigenvalues) > 1e-8) + 1)


def select(paths, laplacian, flight):
    frequency, v = np.linalg.eigh(4 * flight)
    wavelength, w = np.linalg.eigh(laplacian)
    matches = []
    for energy in eigenspaces(frequency):
        for wave in eigenspaces(wavelength):
            overlap = np.linalg.norm(v[:, energy].T @ w[:, wave], ord=2)
            if overlap**2 >= 0.9:
                matches.append((energy, wave))
    [(energy, wave)] = matches
    np.testing.assert_allclose(frequency[energy], [4, 4], atol=1e-12, rtol=0)
    np.testing.assert_allclose(wavelength[wave], [1, 1], atol=1e-12, rtol=0)
    # The donor chart orthogonalizes successive projected coordinate vectors.
    # For disjoint signed sibling pairs this orders the rays by first support.
    rays = []
    for i, j in combinations(range(len(paths)), 2):
        if paths[i][:-1] == paths[j][:-1]:
            ray = np.zeros(len(paths), dtype=int)
            ray[[i, j]] = [1, -1]
            if np.array_equal(flight @ ray, ray) and np.array_equal(laplacian @ ray, ray):
                rays.append(ray)
    normalized = np.array(rays) / sqrt(2)
    np.testing.assert_allclose(normalized.T @ normalized, w[:, wave] @ w[:, wave].T,
                               atol=2e-13, rtol=0)
    return list(zip(map(int, wave), normalized, strict=True))


def read_clock(counts, times, shots, settings):
    if len(times) != len(counts) or any(type(k) is not int or not 0 <= k <= shots for k in counts):
        raise ValueError("Malformed clock record")
    error = sqrt(log(2 * len(times) * settings / 0.01) / (2 * shots)) + 5.1e-14
    lower, upper = 0.0, 128.0
    for t, k in zip(times, counts, strict=True):
        if t * upper > pi:
            return None
        p = k / shots
        lower = max(lower, 2 / t * asin(sqrt(max(0.0, p - error))))
        upper = min(upper, 2 / t * asin(sqrt(min(1.0, p + error))))
        if lower > upper:
            return None
        if p > 1 / 8:
            return min(upper, max(lower, 2 / t * asin(sqrt(p)))), lower, upper
    return None


def replay(root):
    for name, digest in DIGESTS.items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Input digest mismatch for {name}")
    with np.load(root / BOARD, allow_pickle=False) as data:
        paths, ends, lap, flight = physical_histories(data)
        stations = len(data["edges"]) + 1
    rays = select(paths, lap, flight)
    clock = json.loads((root / CLOCK).read_bytes())
    words = json.loads((root / WORDS).read_bytes())
    exhaustion = -expm1(stations * log1p(-2.0**-96))
    lowest, highest_swap, max_drift = 1.0, 0.0, 0.0
    unresolved, releases = [], 0
    for bank in clock["runs"]:
        strength = bank["strength_in_simulator"]
        read = read_clock(bank["counts"], clock["appointments"],
                          clock["shots_per_appointment"], len(clock["runs"]))
        if read is None:
            unresolved.append(strength)
            continue
        estimate, lower, upper = read
        if not lower <= abs(strength) <= upper:
            raise ArithmeticError("The exposure calibration missed this simulated strength")
        ratio = strength / estimate
        pulse = expm(-4j * ratio * words["step"] * flight)
        for column, ray in rays:
            psi = np.eye(len(paths), dtype=complex)[0]
            for angles in words["pulses"][column]["angles"]:
                psi = pulse @ psi
                for port, angle in zip(words["ports"], angles, strict=True):
                    psi[ends == port] *= np.exp(-1j * ratio * angle)
            launch = (1 - exhaustion) * abs(ray @ psi)**2
            lowest = min(lowest, launch)
            for wait in (1, 7, 37, 150):
                terminal = (1 - exhaustion) * abs(ray @ (expm(-0.5j * strength * wait * flight) @ psi))**2
                swap = expm(-0.5 * abs(strength) * wait * lap) @ (abs(psi)**2)
                highest_swap = max(highest_swap, (1 - exhaustion) * (ray**2 @ swap))
                max_drift = max(max_drift, abs(launch - terminal))
            releases += 1

    # A commuting control on the same graph passes the donor's spectroscopic
    # selection exactly. Its first nonconstant wave can still spread over ports.
    k2, guides = np.linalg.eigh(lap)
    first = eigenspaces(k2)[1]
    if len(first) != 1:
        raise ArithmeticError("The commuting-control wave needs a degenerate-space read")
    mass = np.bincount(ends, weights=guides[:, first[0]]**2)
    aperture = sum(sorted(mass)[-2:])
    # The allowance covers propagation roundoff and the sub-ulp difference
    # between unity and the original-start commissioning success probability.
    lower_weight = lowest - 2e-9
    if (unresolved != [0] or releases != 32 or lower_weight <= 0.9999994989
            or highest_swap > 0.5 + 2e-12 or max_drift > 2e-10 or aperture >= 0.98):
        raise ArithmeticError("A signed-return, clock or spatial-control gate failed")
    return {"source_columns": [col for col, _ in rays], "resolved_releases": releases,
            "unresolved_strengths": unresolved, "commissioning_exhaustion": exhaustion,
            "minimum_calculated_signed_weight": round(float(lowest), 12),
            "minimum_calculated_click_probability": round(float(lowest / 8), 12),
            "minimum_signed_weight_after_allowance": round(float(lower_weight), 12),
            "maximum_observed_swap_weight": round(float(highest_swap), 12),
            "largest_signed_weight_drift_below": 2e-10,
            "commuting_control_maximum_two_port_probability": round(float(aperture), 12)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("oph_root", type=Path, help="Directory containing the three pinned inputs")
    print(json.dumps(replay(parser.parse_args().oph_root), indent=2, allow_nan=False))
