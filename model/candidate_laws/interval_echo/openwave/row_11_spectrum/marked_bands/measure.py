"""Prepare three charged orbitals and compare their formation energies and recoil records."""

import argparse
import json

import numpy as np
from scipy.linalg import eigh
from .material import (
    HERE,
    enclosure,
    formation_certificate,
    geometry,
    route_matrix,
    synthesize,
)
from ..conjugate_flight.measure import binomial_interval
from ..conjugate_flight.pulses import native_word

SHOTS = 2**26
CELLS = 48
TIMES = [0, *[2.0**j for j in range(-1, 20)]]
DELTA = 2 * np.pi / CELLS
# Bonferroni allowance for clock, energy and all Fourier quadratures in both charges.
RISK = 0.01 / 1050


def basis(q):
    phase = np.exp(1j * q * np.arange(CELLS)) / np.sqrt(CELLS)
    return np.kron(phase[:, None], np.eye(3))


def orbital(h, q, band):
    b = basis(q)
    matrix = b.conj().T @ h @ b
    np.testing.assert_allclose(h @ b, b @ matrix, atol=2e-13)
    _, u = eigh(matrix)
    v = u[:, band]
    return b @ (v * np.exp(-1j * np.angle(v[np.argmax(abs(v))])))


def counts(rng, state):
    p = abs(state) ** 2
    return rng.multinomial(SHOTS, p / p.sum())


def read_energy(state, rng, bond_word, sign):
    loc = counts(rng, state)
    # Choose a random adjacent bond and interrogate it with the established
    # three-rail native inverse decoder. The middle-port hit tags odd exchange.
    n = len(state)
    analyzer = np.column_stack(
        [
            native_word(v.astype(complex), bond_word, 8, sign, inverse=True)
            for v in np.eye(3)
        ]
    )[1]
    prob = np.mean(
        abs(analyzer @ np.stack([state, np.roll(state, -1), np.roll(state, -2)])) ** 2
    )
    odd = int(rng.binomial(SHOTS, prob))
    mid = int(loc[1::3].sum())
    ci = (
        16
        + 32 * binomial_interval(mid, SHOTS, RISK)
        + n * binomial_interval(odd, SHOTS, RISK)
    )
    return dict(
        ports=loc.tolist(),
        odd_bond=odd,
        energy=16 + 32 * mid / SHOTS + n * odd / SHOTS,
        energy_interval=ci.tolist(),
    )


def read_recoil(eigen, state, rng, clock):
    ev, vec = eigen
    initial = vec.conj().T @ state
    harmonic = np.exp(1j * DELTA * np.arange(CELLS))
    records, estimates = [], []
    eps = np.sqrt(4 * np.log(2 / RISK) / SHOTS)
    for t in TIMES:
        pop = (
            counts(rng, vec @ (np.exp(-1j * t * ev) * initial))
            .reshape(CELLS, 3)
            .sum(axis=1)
        )
        z = harmonic @ pop / SHOTS
        records.append(pop.tolist())
        estimates.append(z)
    angles = np.angle(np.asarray(estimates) / estimates[0])
    phase_error = np.arcsin(eps / np.abs(estimates))
    assert max(phase_error) < 0.03
    # Adjacent exposures double. Select the wrap from the preceding read,
    # beginning at half a tick where the one-step change is below pi.
    phase = angles[1]
    for raw in angles[2:]:
        phase = raw + 2 * np.pi * round((2 * phase - raw) / (2 * np.pi))
    err = phase_error[0] + phase_error[-1]
    ci = [
        (phase - err) / (TIMES[-1] * clock[1]),
        (phase + err) / (TIMES[-1] * clock[0]),
    ]
    return dict(
        by_cell=records,
        gap=phase / TIMES[-1],
        gap_interval=ci,
        weakest_visibility=float(min(np.abs(estimates))),
    )


def flight(h, eigen, q, band, rng):
    # Finite source envelope with a fixed, independently computed three-rail port profile.
    x = np.arange(CELLS) - 20
    v = orbital(h, q, band) * np.repeat(np.exp(-x * x / 36), 3)
    state = synthesize(v / np.linalg.norm(v))
    ev, vec = eigen
    t = 12288 if band == 2 else 256
    z0, z1 = [
        counts(rng, s).reshape(CELLS, 3).sum(axis=1)
        for s in (state, vec @ (np.exp(-1j * t * ev) * (vec.T @ state)))
    ]
    harmonic = np.exp(1j * DELTA * x)
    drift = np.angle((harmonic @ z1) / (harmonic @ z0)) / DELTA
    return dict(time=t, before=z0.tolist(), after=z1.tolist(), displacement=drift)


def run():
    ports, marks = geometry()
    word = json.loads((HERE.parent / "conjugate_flight/pulses.json").read_text())[
        "bond_minus"
    ]
    rng = np.random.default_rng(204961)
    # Independent isolated blank dimer: at a scheduled tick t it crosses with sin(t/2)^2.
    timer = int(rng.binomial(SHOTS, np.sin(0.5) ** 2))
    clock = 2 * np.arcsin(np.sqrt(binomial_interval(timer, SHOTS, RISK)))
    cases = []
    for sign in (-1, 1):
        h = route_matrix(8, sign)
        eig = eigh(h)
        for band in range(3):
            rest_q = np.pi if band == 1 else 0.0
            v = orbital(h, rest_q, band)
            w = orbital(h, rest_q + DELTA, band)
            rested = synthesize(v)
            source = synthesize((v + w) / np.sqrt(2))
            assert np.linalg.norm(rested - v) < 2e-12
            assert np.linalg.norm(source - (v + w) / np.sqrt(2)) < 3e-12
            cases.append(
                dict(
                    sign=sign,
                    band=band,
                    rest_q=rest_q,
                    energy=read_energy(rested, rng, word, sign),
                    recoil=read_recoil(eig, source, rng, clock),
                    flight=[
                        flight(h, eig, rest_q + s * np.pi / 2, band, rng)
                        for s in [-1, 1]
                    ],
                )
            )
    # Same topology, mark and pulse *times*, with a different contrast coefficient.
    # k=40 gives an extra 2*pi on the fixed phase hold in the orienting pair source.
    changed = route_matrix(40, -1)
    change_control = []
    for band in range(3):
        q = np.pi if band == 1 else 0
        nominal = orbital(h, q, band)
        emitted = synthesize(nominal, actual_k=40, source_k=8)
        pop = counts(rng, emitted)
        change_control.append(
            dict(
                band=band,
                ports=pop.tolist(),
                weight_on_nominal=float(abs(np.vdot(nominal, emitted)) ** 2),
                changed_mean_native=float(np.vdot(emitted, changed @ emitted).real),
            )
        )
    return dict(
        shots=SHOTS,
        times=TIMES,
        route=ports,
        side_seams=marks,
        enclosures=[enclosure(s) for s in (-1, 1)],
        formation_bound=formation_certificate(),
        timer_crossed=timer,
        timer_interval=clock.tolist(),
        cases=cases,
        same_pulses_kappa40=change_control,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = json.dumps(run(), indent=2) + "\n"
    if args.write:
        (HERE / "results.json").write_text(data)
    else:
        print(data, end="")
