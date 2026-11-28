"""Independently read the approximate long-wave shell on the finite open guide."""

import argparse
from hashlib import sha256
import json
from math import acos, pi
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal

from .....replay import compare
from ..measure import (
    CELLS,
    FanPath,
    archived_fixture,
    arrival,
    calibration,
    conjugacy,
    counts,
    difference_of_pointer,
    divide,
    generator,
    interval,
    mass,
    phase,
    prepare,
    read_phase,
    symbol,
    trace_route,
    translation,
)
from .measure import clock, spectral, WAITS

HERE = Path(__file__).resolve().parent
SEED = 1045105217548
ORDERS = (2, 3, 4)


def standing(order):
    k = pi * order / CELLS
    _, u, _ = symbol(k, 4, 8)
    j = np.arange(CELLS)
    direct = np.exp(1j * k * j[:, None]) * u
    mirrored = np.exp(-1j * k * (j[:, None] + 1)) * u[::-1]
    packet = (direct + mirrored).ravel()
    return prepare(packet / np.linalg.norm(packet))


def wavelength(rng, packet, gate, pad):
    across = translation(rng, packet, gate, pad)
    assert (
        across["y"]["quadrature_interval"][0]
        < 0
        < across["y"]["quadrature_interval"][1]
    )
    ends = packet.reshape(CELLS, 3)[[0, -1]]
    plus = float((abs(ends[:, 1]) ** 2 + abs(ends[:, 0] + ends[:, 2]) ** 2 / 2).sum())
    minus = float((abs(ends[:, 0] - ends[:, 2]) ** 2 / 2).sum())
    bank = counts(rng, [plus, minus, 1 - plus - minus], gate)
    reflected = difference_of_pointer(bank, pad)
    c = [
        a + b / 2
        for a, b in zip(across["x"]["quadrature_interval"], reflected, strict=True)
    ]
    assert -1 < c[0] < c[1] < 1
    return {
        "translation_counts": across["x"]["counts"],
        "orthogonal_quadrature_counts": across["y"]["counts"],
        "end_reflection_counts": bank,
        "direct_cosine_interval": c,
        "wave_number_interval": [acos(c[1]), acos(c[0])],
    }


def owner_counts(rng, packet, gate):
    return counts(rng, (abs(packet.reshape(CELLS, 3)) ** 2).sum(axis=1), gate)


def run():
    rng = np.random.default_rng(SEED)
    raw = (HERE / "records.json").read_bytes()
    prior = json.loads(raw)
    rest = prior["rest_frequency_in_neutral_clock_units"]
    clock_band = prior["neutral_ring_clock"][-1]["frequency_interval"]
    with archived_fixture() as archive:
        path = FanPath(archive, 4, -1, True)
        positive = FanPath(archive, 4, 1, True)
    capture, words, abs_pad = arrival(path, conjugacy(positive, path)[0])
    gate = mass(capture)
    pad = 4 * abs_pad / (gate - abs_pad) + 1e-9
    _, degrees, _ = trace_route(path, capture, words)
    h = generator(degrees)
    e, v, radius = spectral(h)
    cut = generator(degrees, removed=True)
    ce, cv = eigh_tridiagonal(cut.diagonal(), cut.diagonal(1))
    cut_coo = cut.tocoo()
    assert np.all(cut_coo.row // 3 == cut_coo.col // 3)
    calib = calibration(path, rng)
    rlo, rhi = calib["transfer_rate_interval"]
    vlo, vhi = calib["degree_cost_interval"]
    elo, ehi = rest
    d = [
        (3 * rlo - ehi) * (3 * rlo + vlo - 2 * ehi),
        (3 * rhi - elo) * (3 * rhi + vhi - 2 * elo),
    ]
    c2 = divide([2 * elo * rlo**3, 2 * ehi * rhi**3], d)
    rows = []
    for order in ORDERS:
        packet = standing(order)
        # The exact reflected eigenvector has unit phase visibility.
        eigenvalue = float(np.vdot(packet, h @ packet).real)
        assert np.linalg.norm(h @ packet - eigenvalue * packet) < 2e-12
        w = wavelength(rng, packet, gate, pad)
        phase_reads = clock(
            rng, e, abs(v.T @ packet) ** 2, gate, pad + WAITS[-1] * radius
        )
        formation = divide(phase_reads[-1]["frequency_interval"], clock_band)
        a, b = w["wave_number_interval"]
        lo, hi = formation
        residual = [lo**2 - ehi**2 - c2[1] * b**2, hi**2 - elo**2 - c2[0] * a**2]
        tolerance = 0.5 * a**4
        assert 0 < a < b < 0.27
        assert max(abs(q) for q in residual) < tolerance
        holdout = {
            "order_label": order,
            "wavelength": w,
            "compensated_frequency": phase_reads,
            "frequency_in_neutral_clock": formation,
            "squared_frequency_residual": residual,
            "quartic_tolerance": tolerance,
        }
        if order in (2, 3):
            state = cv @ (np.exp(-1j * 4 * ce) * (cv.T @ packet))
            # With the cells cut, addressed material stays on its respective
            # owner even in the presence of all grade-five and within-cell terms.
            np.testing.assert_allclose(
                (abs(state.reshape(CELLS, 3)) ** 2).sum(axis=1),
                (abs(packet.reshape(CELLS, 3)) ** 2).sum(axis=1),
                atol=1e-13,
                rtol=0,
            )
            t = 0.5
            evolved = np.dot(abs(cv.T @ packet) ** 2, np.exp(-1j * t * ce))
            bank = read_phase(rng, evolved, gate)
            cut_frequency = divide([p / t for p in phase(*bank, pad)], clock_band)
            assert cut_frequency[1] < formation[0]
            before, after = (
                owner_counts(rng, packet, gate),
                owner_counts(rng, state, gate),
            )
            differences = []
            for start, end in zip(before[:-1], after[:-1], strict=True):
                a_owner = interval(start, sum(before[:-1]), pad)
                b_owner = interval(end, sum(after[:-1]), pad)
                difference = [b_owner[0] - a_owner[1], b_owner[1] - a_owner[0]]
                assert difference[0] <= 0 <= difference[1]
                differences.append(difference)
            holdout["intercell_cut"] = {
                "owner_counts_at_0_and_4": [before, after],
                "owner_probability_difference_intervals": differences,
                "half_wait_phase_counts": bank,
                "frequency_pointer_interval": cut_frequency,
            }
        rows.append(holdout)
    return {
        "rest_receipt_sha256": sha256(raw).hexdigest(),
        "protocol_sha256": sha256((HERE / "LONG_WAVE.md").read_bytes()).hexdigest(),
        "original_gate_words": words,
        "local_calibration": calib,
        "rest_frequency": rest,
        "near_rest_squared_speed_from_contacts_and_anchor": c2,
        "holdouts": rows,
    }


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    data = run()
    out = HERE / "long_wave.json"
    if args.write:
        out.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    else:
        compare(data, json.loads(out.read_text()), "long_wave")
        print("Verified three separately clocked long-wave shell comparisons.")


if __name__ == "__main__":
    main()
