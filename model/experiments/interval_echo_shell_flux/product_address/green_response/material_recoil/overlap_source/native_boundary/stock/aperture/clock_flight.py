"""Transfer independently measured paired clock frequencies to free material records."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import clock_covariance as clock
from . import occupied_release as flight

MEAN_AGES = (4, 16, 64)
TAGS = (*map(str, flight.AGES), "uniform_0_64")


def native_spectrum(h, psi, masks, receivers, omega):
    e, v = np.linalg.eigh(h)
    weights = v * (v.T @ psi)
    ages = np.array(flight.AGES)
    waves = weights @ np.exp(-1j * np.outer(e, ages))
    average = ((weights @ np.sinc((e[:, None] - e) * ages[-1] / np.pi)) * weights).sum(
        1
    )
    p = np.column_stack((abs(waves) ** 2, average))
    assert p.min() >= -1e-13
    np.testing.assert_allclose(p.sum(0), 1, atol=5e-13, rtol=0)
    projections = receivers @ weights
    spectrum = np.array(
        [projections @ (1j / (t * (omega - e) + 1j)) for t in MEAN_AGES]
    )
    orth = np.linalg.norm(v.T @ v - np.eye(len(e)), "fro")
    residual = np.linalg.norm(h @ v - v * e, "fro")
    return masks @ p, spectrum, float(4 * orth + 2 * ages[-1] * residual)


def boundary_response(h0, marks, psi, receivers, omega, scale):
    """Use only the ordinary propagator and the values on the four marked sites."""
    sites = np.flatnonzero(marks.any(0))
    ports = np.eye(len(h0))[:, sites]
    rhs = np.column_stack((psi, ports))
    response = np.empty((3, 3, len(MEAN_AGES), len(receivers)), complex)
    error = 0.0
    for k, age in enumerate(MEAN_AGES):
        z = omega + 1j / age
        blank = z * np.eye(len(h0)) - h0
        green = np.linalg.solve(blank, rhs)
        free, legs = green[:, 0], green[:, 1:]
        for i, a in enumerate(clock.LINE):
            for j, b in enumerate(clock.LINE):
                strength = scale * (a * marks[0, sites] + b * marks[1, sites])
                inside = np.linalg.solve(
                    np.eye(len(sites)) - legs[sites] * strength, free[sites]
                )
                routed = free + legs @ (strength * inside)
                residual = blank @ routed - psi
                residual[sites] -= strength * routed[sites]
                error = max(error, float(np.linalg.norm(residual)))
                response[i, j, k] = 1j / age * (receivers @ routed)
    assert error < 1e-10
    return response, error


def calculate(board, edges, rad, clamp, load, scale):
    h, marks = flight.apparatus(board, edges, load, clamp, 1)
    other, _ = flight.apparatus(board, edges, load, clamp, 2)
    active = np.arange(board.n) != clamp
    adj = board.adjacency().toarray()[clamp, active]
    expected = 15.5 if load == 0 else -17
    np.testing.assert_array_equal((other - h) / 2, np.diag(expected * adj))
    h0 = h / 2
    psi = np.full(len(h), 1 / np.sqrt(len(h)))
    masks = np.array([rad[active] <= r for r in flight.CUTS])
    receivers = np.concatenate((masks, ~masks)).astype(float)
    receivers /= np.linalg.norm(receivers, axis=1, keepdims=True)
    omega = float(psi @ h0 @ psi)
    response, residual = boundary_response(h0, marks, psi, receivers, omega, scale)
    p, full = [], []
    spectral_error = 0.0
    for a, b in flight.SECTORS:
        mat = h0 + np.diag(scale * (a * marks[0] + b * marks[1]))
        rates, spectrum, bound = native_spectrum(mat, psi, masks, receivers, omega)
        p.append(rates)
        full.append(spectrum)
        spectral_error = max(spectral_error, bound)
    p = np.reshape(p, (3, 3, len(masks), len(TAGS)))
    full = np.reshape(full, response.shape)
    np.testing.assert_allclose(full, response, atol=1e-11, rtol=0)
    assert spectral_error < 2e-9
    return (
        p,
        full,
        response,
        dict(
            tuning=omega,
            resolvent_residual=residual,
            spectral_error_check=spectral_error,
        ),
    )


def forecast(clock_bank, inverse, source_delta, reading):
    from_clocks = clock.estimate(clock_bank, inverse, inverse, reading)
    sharp = flight.CONDITIONAL @ reading.ravel()
    from_source = np.sort(source_delta * (flight.MIX @ sharp))
    # Include a numerical margin beyond the checked native spectral residual.
    pad = np.array([-1, 1]) * 5e-9
    return from_clocks + pad, from_source + pad


def measurement(law, native, model, field, clock_bank, inverse, delta, rng, signed):
    weights = law.reshape(2, -1)
    results = {}
    scale = 2 if signed else 1
    for key in np.ndindex(native.shape[2:]):

        def column(array):
            return array[(slice(None), slice(None), *key)]

        truth = weights @ column(native).ravel()
        prob = (1 + truth) / 2 if signed else truth
        hits = rng.binomial(clock.N, prob)
        diff = scale * (hits[1] - hits[0]) / clock.N
        held = diff + np.array([-2, 2]) * scale * clock.EPS
        direct, source = forecast(clock_bank, inverse, delta, column(model))
        gravity, gravity_source = forecast(clock_bank, inverse, delta, column(field))
        for predicted in (direct, source):
            assert max(predicted[0], held[0]) <= min(predicted[1], held[1])
        excluded = max(gravity[0], held[0]) > min(gravity[1], held[1])
        results[",".join(map(str, key))] = dict(
            hits=hits.tolist(),
            held=held.tolist(),
            local_clocks=direct.tolist(),
            source=source.tolist(),
            ordinary_gain=gravity.tolist(),
            source_at_ordinary_gain=gravity_source.tolist(),
            excluded=bool(excluded),
        )
    return results


def run():
    board, old, edges, rad, _ = flight.chart()
    rng = np.random.default_rng(191156)
    clamp, law, delta, bank, inverse, ratio, report = clock.experiment(
        board, old, edges, rng
    )
    out = {}
    for label, load in (("vacancy", 0), ("excess", 2)):
        native, timing, compressed, checked = calculate(
            board, edges, rad, clamp, load, 1
        )
        if load == 0:
            scaled, field_response = native, compressed
            ordinary_check = checked
        else:
            scaled, _, field_response, ordinary_check = calculate(
                board, edges, rad, clamp, load, float(ratio)
            )
        p = measurement(law, native, native, scaled, bank, inverse, delta, rng, False)
        real = measurement(
            law,
            timing.real,
            compressed.real,
            field_response.real,
            bank,
            inverse,
            delta,
            rng,
            True,
        )
        imag = measurement(
            law,
            timing.imag,
            compressed.imag,
            field_response.imag,
            bank,
            inverse,
            delta,
            rng,
            True,
        )
        out[label] = dict(
            population=p,
            ramsey_real=real,
            ramsey_imag=imag,
            native_check=checked,
            ordinary_check=ordinary_check,
        )
    return dict(
        sha256={
            f: sha256((flight.HERE / f).read_bytes()).hexdigest()
            for f in ("CORRELATED_CLOCKS.md", "pair_contact.json")
        },
        trials=clock.N,
        clock=report,
        delta=delta.tolist(),
        cuts=list(flight.CUTS),
        flight_ages=list(TAGS),
        mean_ages=list(MEAN_AGES),
        receiver_order="inside each cut, then outside each cut",
        flight=out,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    path = flight.HERE / "clock_flight.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        flight.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), result
        )
    for label, out in result["flight"].items():
        print(
            label,
            {
                k: sum(r["excluded"] for r in out[k].values())
                for k in ("population", "ramsey_real", "ramsey_imag")
            },
            out["population"]["1,3"],
        )
