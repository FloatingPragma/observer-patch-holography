"""Compare observer records of geometry, formation and recoil without mass inputs."""

import argparse
import json
from itertools import permutations
from math import acos, asin, atan2, ceil, floor, log, pi, sqrt
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from ..phase_feedback import measure as phase_read

SHOTS = phase_read.SHOTS
SLACK = 5e-7
RECOIL_SETTINGS = ((1 / 64, 0, 0), (0, 1 / 32, 0))


def quadrature(response, rng, acceptance=1.0):
    """Read two ports and keep every failed return."""
    probabilities = np.array(
        [
            abs(1 + response) ** 2 / 4,
            abs(1 - response) ** 2 / 4,
            max(0.0, (1 - abs(response) ** 2) / 2),
        ]
    )
    probabilities *= acceptance / probabilities.sum()
    probabilities[-1] += 1 - acceptance
    return rng.multinomial(SHOTS, probabilities).tolist()


def ramsey_records(symbol, states):
    """Read the formation phase relative to the same compensated idle arm."""
    eigenvalues, basis = np.linalg.eigh(symbol)
    weights = abs(basis.conj().T @ states) ** 2
    rng = np.random.default_rng(726041)
    preparations = []
    for one_branch in weights.T:
        rows = []
        for exponent in phase_read.POWERS:
            ticks = 2**exponent
            z = one_branch @ np.exp(-1j * ticks * eigenvalues)
            rows.append(
                dict(ticks=ticks, counts=[quadrature(z, rng), quadrature(1j * z, rng)])
            )
        preparations.append(
            dict(
                energy_interval_99=phase_read.infer(rows, 1.0, SLACK), observations=rows
            )
        )
    return dict(
        shots_per_quadrature=SHOTS,
        seed=726041,
        numerical_slack=SLACK,
        preparations=preparations,
    )


def geometry_record(unitary, detectors, preparation):
    """One shot bank with three independent detector names and an inactive outcome."""
    response = abs(detectors.conj().T @ unitary @ preparation) ** 2
    rng = np.random.default_rng(726042)
    counts = rng.multinomial(SHOTS, np.r_[response, 1 - response.sum()]).tolist()
    return dict(shots=SHOTS, seed=726042, counts=counts)


def terminal_record(probabilities):
    rng = np.random.default_rng(726043)
    returns = [int(rng.binomial(SHOTS, value)) for value in probabilities]
    return dict(
        wait=1,
        shots=SHOTS,
        seed=726043,
        returned_and_absent=[[n, SHOTS - n] for n in returns],
    )


def native_phase_banks(generator, modes):
    """Retain both quadratures of the equal three-input mixture on native flight."""
    rng = np.random.default_rng(726045)
    prepared = np.zeros((generator.shape[0], 3), dtype=complex)
    prepared[: len(modes)] = modes
    banks = []
    for stride in (1 / 8, 1 / 4, 1.0):
        flights = expm_multiply(
            -1j * generator, prepared, start=0, stop=6 * stride, num=7
        )
        overlap = np.einsum("ij,tij->t", prepared.conj(), flights) / 3
        counts = [[quadrature(val, rng), quadrature(1j * val, rng)] for val in overlap]
        banks.append(dict(stride=stride, counts=counts))
    return dict(shots_per_quadrature=SHOTS, seed=726045, banks=banks)


def native_frequency_check(records):
    """A four-by-four minor excludes an unfiltered sum of just three phases."""
    banks = records["banks"]
    coordinate_radius = sqrt(2 * log(2 * 7 * 2 * len(banks) / 0.01) / SHOTS) + SLACK
    matrix_error = 4 * sqrt(2) * coordinate_radius
    moments = np.add.outer(np.arange(4), np.arange(4))
    output = []
    for bank in banks:
        coordinates = np.array(bank["counts"])
        phases = (coordinates[..., 0] - coordinates[..., 1]) / records[
            "shots_per_quadrature"
        ]
        response = phases[:, 0] - 1j * phases[:, 1]
        fourth = np.linalg.svd(response[moments], compute_uv=False)[-1]
        output.append(
            dict(
                stride=bank["stride"],
                fourth_singular_value=float(fourth),
                lower_bound=max(0.0, float(fourth) - matrix_error),
            )
        )
    return output


def recoil_record(symbol, states, momentum, rng):
    """Read one temporal cosine bank per tag and independent spatial phase banks."""
    eigenvalues, basis = np.linalg.eigh(symbol)
    # A creation attempt above the filled negative band supplies this preparation.
    positive = eigenvalues > 0
    weights = abs(basis[:, positive].conj().T @ states) ** 2
    acceptance = np.minimum(weights.sum(axis=0), 1)
    injection = [int(rng.binomial(SHOTS, x)) for x in acceptance]
    weights = weights / weights.sum(axis=0)
    temporal = []
    for exponent in phase_read.POWERS:
        ticks = 2**exponent
        z = np.exp(-1j * ticks * eigenvalues[positive]) @ weights
        temporal.append(
            dict(
                ticks=ticks,
                counts=[
                    quadrature(response, rng, accepted)
                    for response, accepted in zip(z, acceptance, strict=True)
                ],
            )
        )
    spatial = []
    for baseline in (8, 16):
        phase = np.exp(-1j * baseline * np.linalg.norm(momentum))
        spatial.append(
            dict(
                baseline=baseline,
                counts=[quadrature(phase, rng), quadrature(1j * phase, rng)],
            )
        )
    return dict(
        shots_per_quadrature=SHOTS,
        injection_counts=[[n, SHOTS - n] for n in injection],
        temporal=temporal,
        spatial=spatial,
    )


def frequency_from_cosines(record, branch, radius):
    """Fold all cosine aliases into the previous interval before each new dwell."""
    lo, hi = 0.0, 2.0
    acceptance = record["injection_counts"][branch][0] / SHOTS
    assert acceptance > radius
    for row in record["temporal"]:
        left, right, _ = row["counts"][branch]
        x = (left - right) / SHOTS
        endpoints = [
            z / p
            for z in (x - radius, x + radius)
            for p in (acceptance - radius, acceptance + radius)
        ]
        a = acos(min(1.0, max(endpoints)))
        b = acos(max(-1.0, min(endpoints)))
        ticks = row["ticks"]
        intervals = []
        for start, stop in ((a, b), (-b, -a)):
            first = ceil((ticks * lo - stop) / (2 * pi))
            last = floor((ticks * hi - start) / (2 * pi))
            for winding in range(first, last + 1):
                intervals.append(
                    (
                        max(lo, (start + 2 * pi * winding) / ticks),
                        min(hi, (stop + 2 * pi * winding) / ticks),
                    )
                )
        if not intervals:
            raise ArithmeticError("Inconsistent recoil phase records.")
        lo = min(interval[0] for interval in intervals)
        hi = max(interval[1] for interval in intervals)
    return [lo, hi]


def wavelength_momentum(rows, radius):
    bounds = []
    for row in rows:
        x, y = ((left - right) / SHOTS for left, right, _ in row["counts"])
        center = atan2(y, x)
        half = asin(sqrt(2) * radius / np.hypot(x, y))
        baseline = row["baseline"]
        bounds.append(((center - half) / baseline, (center + half) / baseline))
    return intersect(bounds)


def intersect(bounds):
    return [max(lo for lo, _ in bounds), min(hi for _, hi in bounds)]


def quotient(bounds):
    lo, hi = np.array(bounds).T
    return [
        float(lo.sum() / np.sqrt(hi).sum() ** 2),
        float(hi.sum() / np.sqrt(lo).sum() ** 2),
    ]


def geometry_intervals(geometry):
    # Invert |frequency - p| <= sqrt(2*p*r) + 2*r/3 from Bernstein's inequality.
    # Three marginal binary variables are covered together by a union bound.
    r = log(6 / 0.01) / geometry["shots"]
    slack = sqrt(r / 2)
    intervals = []
    for count in geometry["counts"][:3]:
        f = count / geometry["shots"]
        lo = max(0.0, sqrt(max(0.0, f - r / 6)) - slack) ** 2
        hi = min(1.0, (sqrt(f + 7 * r / 6) + slack) ** 2)
        intervals.append([lo, hi])
    return intervals


def reconstruct(geometry, formation, recoil, registered_return):
    """Identify a common-scale triple using only named detector outcomes."""
    probabilities = geometry_intervals(geometry)
    frequencies = [
        phase_read.infer(row["observations"], 1.0, formation["numerical_slack"])
        for row in formation["preparations"]
    ]
    geometric_q, formation_q = quotient(probabilities), quotient(frequencies)
    assert geometric_q[0] < 2 / 3 < geometric_q[1]
    assert formation_q[0] < 2 / 3 < formation_q[1]
    candidates = []
    for names in permutations(range(3)):
        scale = intersect(
            [
                (f[0] / probabilities[i][1], f[1] / probabilities[i][0])
                for f, i in zip(frequencies, names, strict=True)
            ]
        )
        if scale[0] < scale[1]:
            candidates.append(
                dict(geometry_channel_for_each_formation_tag=names, scale=scale)
            )
    assert len(candidates) == 1

    banks = len(recoil) * (3 * len(phase_read.POWERS) + 7)
    radius = sqrt(2 * log(2 * banks / 0.01) / SHOTS) + SLACK
    recoil_reads = []
    for record in recoil:
        k0, k1 = wavelength_momentum(record["spatial"], radius)
        assert 0 < k0 < k1
        fits = [frequency_from_cosines(record, tag, radius) for tag in range(3)]
        coefficients = []
        for (e0, e1), (m0, m1) in zip(fits, frequencies, strict=True):
            coefficients.append([(e0**2 - m1**2) / k1**2, (e1**2 - m0**2) / k0**2])
        recoil_reads.append(
            dict(
                momentum_interval=[k0, k1],
                energy_intervals=fits,
                speed_squared=coefficients,
            )
        )
    common = intersect(
        [interval for read in recoil_reads for interval in read["speed_squared"]]
    )
    assert common[0] < 1 < common[1]
    sorted_tags = sorted(range(3), key=lambda tag: frequencies[tag][0])
    measured = [
        [
            frequencies[i][0] / frequencies[sorted_tags[0]][1],
            frequencies[i][1] / frequencies[sorted_tags[0]][0],
        ]
        for i in sorted_tags[1:]
    ]
    assert all(
        rest[1] < later[0]
        for rest, later in zip(
            [frequencies[i] for i in sorted_tags[:-1]],
            [frequencies[i] for i in sorted_tags[1:]],
            strict=True,
        )
    )
    return dict(
        geometry_probability_intervals=probabilities,
        mapping=candidates[0],
        geometry_quotient_interval=geometric_q,
        formation_quotient_interval=formation_q,
        formation_tags_light_to_heavy=sorted_tags,
        measured_mass_ratios=measured,
        common_speed_squared_interval=common,
        recoil=recoil_reads,
        registered_return_fractions=[
            n / registered_return["shots"]
            for n, _ in registered_return["returned_and_absent"]
        ],
    )


def decode_receipt(path):
    with path.open() as handle:
        receipt = json.load(handle)
    observers = receipt["observer_records"]
    field = receipt["field"]
    return dict(
        reconstruction=reconstruct(
            observers["geometry"],
            field["rest_interferometer"],
            field["recoil_interferometers"],
            observers["native_free_return"],
        ),
        native_frequency_check=native_frequency_check(observers["native_phase_counts"]),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt", type=Path)
    print(json.dumps(decode_receipt(parser.parse_args().receipt), indent=2))
