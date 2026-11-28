"""Physical vacancy counts on the independent winding-spectroscopy material."""

import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import least_squares

from model.experiments.interval_echo_angular.charged_orientation.winding_memory import (
    mixed_masses as source,
)

HERE = Path(__file__).resolve().parent
N, SHOTS, FFT = 144, 262144, 2304
HOLDS = np.arange(1, 113) / 4
TRAIN = HOLDS <= 4
LIMIT = 0.015
PROTOCOL = "1f45da04de263ebe1d3d165d09f934cb2db1d40cb539139888d64a1596b1096f"


def generator(motif):
    pattern = motif * (N // len(motif))
    word, offset = source.lifted(pattern)
    d, h = source.coefficients(word, "even")
    energy, basis = eigh_tridiagonal(np.array(d) / 8, np.array(h) / 8)
    return (np.arange(len(energy)) + offset) % N, energy, basis


def release(ports, energy, basis, initial, times):
    amplitude = (np.exp(-1j * np.outer(times, energy)) * (basis.T @ initial)) @ basis.T
    density = np.zeros((len(times), N))
    for i, port in enumerate(ports):
        density[:, port] += abs(amplitude[:, i]) ** 2
    np.testing.assert_allclose(density.sum(axis=1), 1, atol=1e-11)
    return density


def acquire(density, rng):
    return np.array([rng.multinomial(SHOTS, p / p.sum()) for p in density])


def symbol(k):
    # A site and the B to its right share a blank, with a written link to B on the left.
    z = -0.5 - 4 * np.exp(-1j * k)
    return z, abs(z)


def bulk_prediction(splitting):
    k = 2 * np.pi * np.fft.fftfreq(FFT)
    z, exact = symbol(k)
    phase = np.outer(HOLDS, splitting)
    a = np.fft.ifft(np.cos(phase), axis=1)
    b = np.fft.ifft(-1j * np.sin(phase) * z.conj() / exact, axis=1)
    physical = np.zeros((len(HOLDS), N))
    for component, parity in ((a, 0), (b, 1)):
        by_cell = (abs(component) ** 2).reshape(len(HOLDS), -1, N // 2).sum(axis=1)
        physical[:, parity::2] = np.roll(by_cell, N // 4, axis=1)
    return physical


def fit_rival(training):
    # q = pi-k, with its rest at k=pi.
    x = 1 + np.cos(2 * np.pi * np.fft.fftfreq(512))

    def residual(p):
        phases = np.outer(HOLDS[TRAIN], p[0] + p[1] * x)
        return (np.cos(phases) ** 2).mean(axis=1) - training

    mesh = [(a, b) for a in np.linspace(3, 4.5, 31) for b in np.linspace(0, 1.5, 31)]
    first = sorted(mesh, key=lambda p: np.linalg.norm(residual(p)))[:5]
    fits = [least_squares(residual, p, bounds=([3, 0], [4.5, 1.5])) for p in first]
    return min(fits, key=lambda fit: fit.cost).x


def errors(counts, prediction):
    diff = counts[~TRAIN] / SHOTS - prediction[~TRAIN]
    even = np.max(abs(diff[:, ::2].sum(axis=1)))
    ports = np.max(abs(diff))
    return dict(
        max_even_error=even, max_port_error=ports, passed=bool(max(even, ports) < LIMIT)
    )


def impulse(raw):
    rng = np.random.default_rng(300144)
    stocks = {}
    for name, motif in (("alternating", "01"), ("triple_blank_16", "01" * 7 + "00")):
        ports, energy, basis = generator(motif)
        initial = np.zeros(len(energy))
        initial[np.flatnonzero(ports == 72).item()] = 1
        counts = acquire(release(ports, energy, basis, initial, HOLDS), rng)
        raw[name + "_impulse_counts"] = counts
        raw[name + "_matched_host_levels"] = energy
        stocks[name] = counts
    observed = stocks["alternating"][:, ::2].sum(axis=1) / SHOTS
    params = fit_rival(observed[TRAIN])
    k = 2 * np.pi * np.fft.fftfreq(FFT)
    native = bulk_prediction(symbol(k)[1])
    rival = bulk_prediction(params[0] + params[1] * (1 + np.cos(k)))
    measured = {}
    for name, counts in stocks.items():
        measured[name] = dict(
            lattice=errors(counts, native), fitted_hopping=errors(counts, rival)
        )
    return dict(
        shots_per_hold=SHOTS,
        calibration_holds=int(TRAIN.sum()),
        held_holds=int((~TRAIN).sum()),
        tolerance=LIMIT,
        fitted_hopping=params.tolist(),
        stocks=measured,
    )


def packet(ports, q, branch):
    # Interior support excludes both multiply visited ports at the winding ends.
    k = np.pi - q
    phase, gap = symbol(k)
    cells = ports // 2
    envelope = np.exp(-((cells - 36) ** 2) / 100) * np.exp(1j * k * cells)
    state = envelope * np.where(ports % 2 == 0, 1, branch * phase.conj() / gap)
    state[(ports < 4) | (ports >= 140)] = 0
    return state / np.linalg.norm(state)


def group_velocity(ports, state):
    waves = []
    for p in (0, 1):
        src = np.zeros(FFT, dtype=complex)
        mask = (ports % 2 == p) & (abs(state) > 0)
        src[ports[mask] // 2] = state[mask]
        waves.append(np.fft.fft(src) / np.sqrt(FFT))
    k = 2 * np.pi * np.fft.fftfreq(FFT)
    z, gap = symbol(k)
    plus, minus = [
        abs((waves[0] + s * z / gap * waves[1]) / np.sqrt(2)) ** 2 for s in (1, -1)
    ]
    np.testing.assert_allclose(np.sum(plus + minus), 1, atol=2e-13)
    return np.sum((plus - minus) * (-2 * np.sin(k) / gap)), plus.sum(), minus.sum()


def fringe(state, ports, pair, rng):
    u, v = [state[np.flatnonzero(ports == 2 * cell).item()] for cell in pair]
    trials = []
    for phase in (1, -1j):
        p = np.array([abs(u + phase * v) ** 2 / 2, abs(u - phase * v) ** 2 / 2])
        trials.append(rng.multinomial(SHOTS, [*p, 1 - p.sum()]))
    x, y = [(row[0] - row[1]) / SHOTS for row in trials]
    q = np.angle(-np.exp(-1j * np.arctan2(y, x)))
    return q, trials


def packets(raw):
    ports, energy, basis = generator("01")
    rng = np.random.default_rng(301144)
    rows, hist, fringes = [], [], []
    for q in (0.2, -0.2, 0.4, -0.4, 0.6, -0.6, 1, -1, 1.6, -1.6):
        for branch in (-1, 1):
            state = packet(ports, q, branch)
            counts = acquire(release(ports, energy, basis, state, [0, 12]), rng)
            speed = ((counts[1] - counts[0]) @ (np.arange(N) // 2)) / (SHOTS * 12)
            observed_q, local = [], []
            for pair in ((33, 34), (37, 38)):
                reading, trials = fringe(state, ports, pair, rng)
                observed_q.append(reading)
                local.append(trials)
            group, plus, minus = group_velocity(ports, state)
            point = -branch * 2 * np.sin(q) / np.sqrt(3.5**2 + 8 * np.sin(q / 2) ** 2)
            continuum = -branch * 2 * q / np.sqrt(3.5**2 + 2 * q**2)
            rows.append(
                dict(
                    q=q,
                    branch=branch,
                    source_fringe_q=observed_q,
                    spectral_weight=plus if branch == 1 else minus,
                    measured_velocity=speed,
                    predicted_velocity=group,
                    point_velocity=point,
                    continuum_point_velocity=continuum,
                    passed=bool(
                        np.sign(speed) == np.sign(group) and abs(speed - group) < 0.025
                    ),
                )
            )
            hist.append(counts)
            fringes.append(local)
    raw["packet_counts"] = np.array(hist)
    raw["fringe_counts"] = np.array(fringes)
    return rows


def discriminant(motif):
    # Trace of a cyclic Jacobi transfer product with integer coefficient polynomials.
    diag, hop = source.coefficients(motif, "even", cyclic=True)
    ident, zero = [1], [0]

    def mul(a, b):
        out = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] += x * y
        return out

    def add(a, b):
        return [
            (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
            for i in range(max(len(a), len(b)))
        ]

    a, b, c, d = ident, zero, zero, ident
    den = 1
    for i, h in enumerate(hop):
        linear, previous = [-diag[i], 8], [-hop[i - 1]]
        a, b, c, d = (
            add(mul(linear, a), mul(previous, c)),
            add(mul(linear, b), mul(previous, d)),
            mul([h], a),
            mul([h], b),
        )
        den *= h
    trace = add(a, d)
    assert len(trace) == len(motif) + 1 and trace[-1] == 8 ** len(motif)
    if len(motif) == 2:
        assert [2 * x for x in trace] == [den * x for x in (4, -9, 1)]
    return dict(
        primitive_period=len(motif),
        trace_degree=len(trace) - 1,
        leading_numerator=str(trace[-1]),
        denominator=str(den),
    )


def run():
    assert hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest() == PROTOCOL
    raw = {}
    terms = source.contacts()
    route, _ = source.material.geometry()
    stock = source.stock_routes(route, ["01", "01" * 7 + "00"])
    response = impulse(raw)
    movement = packets(raw)
    archive = io.BytesIO()
    np.savez_compressed(archive, **raw)
    data = archive.getvalue()
    out = dict(
        protocol_sha256=PROTOCOL,
        counts_sha256=hashlib.sha256(data).hexdigest(),
        stock=stock,
        source_receipt="5301ed528344e89a6115cda2bdcb99b3feccba28",
        contact_terms_in_eighths=terms,
        physical_ports=route,
        minimum_formation=float(raw["alternating_matched_host_levels"].min()),
        maximum_formation=float(raw["alternating_matched_host_levels"].max()),
        winding_states=len(raw["alternating_matched_host_levels"]),
        impulse=response,
        packets=movement,
        transfer_traces=[
            discriminant(s)
            for s in [
                "01",
                *["01" * ((p - 2) // 2) + "00" for p in (12, 16, 18, 24, 36, 48, 72)],
            ]
        ],
    )
    return json.dumps(out, indent=2, sort_keys=True, default=int) + "\n", data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    summary, raw = run()
    for name, expected in (("results.json", summary.encode()), ("counts.npz", raw)):
        path = HERE / name
        if args.write:
            path.write_bytes(expected)
        else:
            assert expected == path.read_bytes(), f"Replay mismatch for {name}"
    print("Written" if args.write else "Verified", "winding charge records")


if __name__ == "__main__":
    main()
