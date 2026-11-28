"""Calibrate and read the compensated rest pole on a sharp negative arrival."""

import argparse
from hashlib import sha256
import json
from math import pi, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal

from .....replay import compare
from ..measure import (
    FanPath,
    archived_fixture,
    arrival,
    conjugacy,
    counts,
    divide,
    generator,
    interval,
    mass,
    phase,
    read_phase,
    trace_route,
)
from .exact import certificate

HERE = Path(__file__).resolve().parent
WAITS = (0.5, *(2**j for j in range(12)))
SIZE, STAGES, STEP = 2**15, 6, 1 / 8
SEED = 104510521753788


def clock(rng, energies, weights, gate, pad):
    energies = np.asarray(energies, dtype=float)
    previous = [0, 4]
    banks = []
    for t in WAITS:
        signal = np.dot(weights, np.exp(-1j * t * energies))
        bank = read_phase(rng, signal, gate)
        angle = phase(*bank, pad)
        n = round((sum(previous) * t - sum(angle)) / (4 * pi))
        band = [(a + n * 2 * pi) / t for a in angle]
        assert (previous[1] - previous[0]) * t < 2 * pi
        previous = [max(previous[0], band[0]), min(previous[1], band[1])]
        assert previous[0] < previous[1]
        banks.append({"wait": t, "counts": bank, "frequency_interval": previous})
    return banks


def spectral(h):
    e, v = eigh_tridiagonal(h.diagonal(), h.diagonal(1))
    orthogonal_error = np.linalg.norm(v.T @ v - np.eye(len(e)), ord="fro")
    residual = np.linalg.norm(h @ v - v * e, ord="fro")
    assert orthogonal_error < 5e-13 and residual < 2e-12
    # The polar factor changes the frame by at most twice its orthogonality
    # defect. With ||H|| and max|e|<=24 the backward bound is below 1e-10.
    radius = 1e-10
    assert np.max(np.asarray(abs(h).sum(axis=1))) <= 24
    assert e[-1] < 24
    return e, v, radius


def accepted(energies, center):
    x = (energies - center) * STEP / (2 * pi)
    return (np.sinc(SIZE * x) / np.sinc(x)) ** (2 * STAGES)


def pole_bounds(rng, p, gate, pad, unit_loss, tail):
    bank = counts(rng, [p, 1 - p], gate)
    read = interval(bank[0], sum(bank[:-1]), pad)
    weighted = [max(0, read[0] - tail), read[1]]
    weight = [weighted[0], min(1, weighted[1] / (1 - unit_loss))]
    return {
        "accept_miss_and_gate_refusal": bank,
        "acceptance_given_gate": read,
        "pole_weight_given_gate": weight,
    }


def run():
    rng = np.random.default_rng(SEED)
    with archived_fixture() as archive:
        route = FanPath(archive, 4, -1, True)
        positive = FanPath(archive, 4, 1, True)
    captured, words, source_error = arrival(route, conjugacy(positive, route)[0])
    gate = mass(captured)
    pad = 4 * source_error / (gate - source_error) + 1e-9
    ports, degrees, _ = trace_route(route, captured, words)
    h = generator(degrees)
    np.testing.assert_array_equal(h.diagonal(1), -4)
    algebra = certificate(h.diagonal().astype(int).tolist(), 4)
    e, v, radius = spectral(h)
    anchor = 10 - 2 * sqrt(17)
    trial = np.tile([1, 1 - anchor / 4, 1], len(ports) // 3)
    trial /= np.linalg.norm(trial)
    assert np.linalg.norm(h @ trial - anchor * trial) < 1e-13
    ring = route.field.ring_laplacian(route.field.words[route.heads])
    np.testing.assert_allclose(ring @ route.vacuum, 0, atol=1e-14)
    np.testing.assert_allclose(ring @ route.code[:, 0], route.code[:, 0], atol=1e-14)
    # The interferometric X,Y signal on the neutral reference is the
    # relative phase of its prepared even and odd branches.
    neutral = clock(rng, [1], [1], 1, 1e-10)
    reads = clock(rng, e, abs(v.T @ trial) ** 2, gate, pad + WAITS[-1] * radius)
    frequency = reads[-1]["frequency_interval"]
    base = neutral[-1]["frequency_interval"]
    center = sum(frequency) / 2
    delta = (frequency[1] - frequency[0]) / 2
    assert frequency[0] < anchor < frequency[1]
    unit_loss = STAGES * ((SIZE * STEP * delta) ** 2) / 12
    # All undesired detunings lie between 0.0016-delta and 24. The
    # argument of the denominator sine is below pi/2, so checking its
    # lower endpoint suffices. sin(x)>=x-x^3/6 gives the rational guard.
    x = (algebra["rest_gap_lower"] - delta) * STEP / 2
    assert x - x**3 / 6 > 0.000099
    tail = (1 / (SIZE * 0.000099)) ** (2 * STAGES)
    probability = accepted(e, center)
    numerical = 2 * STAGES * SIZE * STEP * radius
    results = {}
    for name, p in (
        ("sharp_1045", abs(v[0]) ** 2 @ probability),
        ("native_swap_1052", abs(v[1]) ** 2 @ probability),
        ("filled_compensated", accepted(np.array([0]), center)[0]),
    ):
        results[name] = pole_bounds(
            rng, float(p), gate, pad + numerical, unit_loss, tail
        )
    a = results["sharp_1045"]["pole_weight_given_gate"]
    b = results["native_swap_1052"]["pole_weight_given_gate"]
    blank = results["filled_compensated"]["pole_weight_given_gate"]
    assert a[0] > b[1] > b[0] > blank[1] >= blank[0] == 0
    assert a[0] < abs(trial[0]) ** 2 < a[1]
    assert b[0] < abs(trial[1]) ** 2 < b[1]
    differences = e[:, None] - e[None, :]
    i, j = np.unravel_index(np.argmin(abs(differences - anchor)), differences.shape)
    nearest = float(differences[i, j])
    error = 2 * radius + delta
    certified_separation = float(np.min(abs(differences - center))) - error
    assert certified_separation > 0.00014
    witness = HERE.parents[3] / "fan_path" / "local_launch" / "witness.json"
    return {
        "source_arrival_witness_sha256": sha256(witness.read_bytes()).hexdigest(),
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "signed_two_word_gate": words,
        "gate_probability_interval": [gate - source_error, gate + source_error],
        "clock_bypass_and_write_refusals": [0, 0],
        "neutral_ring_clock": neutral,
        "compensated_formation_clock": reads,
        "rest_frequency_in_neutral_clock_units": divide(frequency, base),
        "filter": {
            "center_measured_in_appointment_units": center,
            "coherent_wait_appointments_per_stage": SIZE,
            "stages": STAGES,
            "appointment_spacing": STEP,
            "pole_attenuation_upper": unit_loss,
            "off_pole_response_upper": tail,
            "numerical_probability_allowance": numerical,
        },
        "source_banks": results,
        "finite_autonomous_clock": {
            "modular_and_gap_certificate": algebra,
            "nearest_positive_bohr_pair_indices": [int(i), int(j)],
            "nearest_bohr_frequency_numeric": nearest,
            "measured_rest_to_any_bohr_line_minimum_separation": certified_separation,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = HERE / "records.json"
    result = run()
    if args.write:
        receipt.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    else:
        compare(result, json.loads(receipt.read_text()), "pole_clock")
        print(
            "Verified sharp-arrival pole weights and finite autonomous-clock exclusion."
        )


if __name__ == "__main__":
    main()
