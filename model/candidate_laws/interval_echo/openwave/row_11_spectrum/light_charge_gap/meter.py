"""Prepare each finite line with native pulses and read energy from local counts."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.stats import beta
from .measure import bench, running_coefficients


def contact(state, bond, time, local):
    segment = slice(bond, bond + 2)
    state[segment] = expm(-1j * time * local) @ state[segment]


def phase(state, site, angle, local, private):
    # Route the chosen coordinate to the existing private contact and back.
    for bond in reversed(range(site)):
        contact(state, bond, np.pi, local)
    hold = ((-angle * np.sign(private)) % (2 * np.pi)) / abs(private)
    state[0] *= np.exp(-1j * private * hold)
    for bond in range(site):
        contact(state, bond, np.pi, local)


def emit(target, local, private):
    state = np.array([1, 0, 0, 0], complex)
    for j in range(len(target) - 1):
        tail = np.linalg.norm(target[j + 1 :])
        contact(state, j, 2 * np.arctan2(tail, abs(target[j])), local)
        phase(state, j, np.angle(target[j]) - np.angle(state[j]), local, private)
    phase(state, 3, np.angle(target[3]) - np.angle(state[3]), local, private)
    return state


def interval(k, n, tail):
    lo = 0.0 if k == 0 else beta.ppf(tail, k, n - k + 1)
    hi = 1.0 if k == n else beta.isf(tail, k + 1, n - k)
    return np.array([lo, hi])


def read():
    rng = np.random.default_rng(871217)
    shots = 2**30
    tail = 0.01 / (2 * 8 * 4)
    rows = []
    for sign in [-1, 1]:
        local, private = running_coefficients(sign)
        h = bench(sign)[3][True]
        energies, waves = np.linalg.eigh(h)
        tip_cost = 8 + private
        gateway_cost = 16
        for index, (energy, wave) in enumerate(zip(energies, waves.T, strict=True)):
            state = emit(wave, local, private)
            np.testing.assert_allclose(state, wave, atol=2e-12)
            state = expm(-7j * h) @ state
            ports = rng.multinomial(shots, abs(state) ** 2 / np.linalg.norm(state) ** 2)
            gateways = int(sum(ports[1:3]))
            odd = []
            for j in range(3):
                probe = state.copy()
                phase(probe, j + 1, np.pi / 2, local, private)
                contact(probe, j, np.pi / 2, local)
                p = abs(probe[j]) ** 2
                np.testing.assert_allclose(
                    p, abs(state[j] - state[j + 1]) ** 2 / 2, atol=2e-12
                )
                odd.append(int(rng.binomial(shots, p)))
            base = tip_cost + (gateway_cost - tip_cost) * gateways / shots
            estimate = base + sum(odd) / shots
            band = np.sort(
                tip_cost + (gateway_cost - tip_cost) * interval(gateways, shots, tail)
            )
            for count in odd:
                band += interval(count, shots, tail)
            rows.append(
                dict(
                    sign=sign,
                    level=index,
                    exact=energy,
                    port_counts=ports.tolist(),
                    odd_bond_counts=odd,
                    measured=estimate,
                    simultaneous_99=band.tolist(),
                )
            )
    return dict(seed=871217, starts_per_bank=shots, free_hold=7, rows=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = json.dumps(read(), indent=2) + "\n"
    if parser.parse_args().write:
        Path(__file__).with_name("meter.json").write_text(data)
    else:
        print(data, end="")
