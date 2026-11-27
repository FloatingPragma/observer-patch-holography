"""Separated Ramsey, position and signed-impulse readers for a compensated material pair."""

import argparse
import json
from pathlib import Path

import numpy as np

from .guide import Flight, physical_rings, signed_pulse_audit
from .release import release_report

SHOTS = 2**20
TIMES = (0, 64, 256, 768)


def probabilities(values):
    p = np.maximum(np.real_if_close(values), 0)
    np.testing.assert_allclose(sum(p), 1, atol=2e-10, rtol=0)
    return p / sum(p)


def draws(rng, p):
    return rng.multinomial(SHOTS, probabilities(p)).tolist()


def ramsey(amplitude, rng):
    results = []
    for value in (amplitude.real, amplitude.imag):
        results.append(draws(rng, [(1 - value) / 2, (1 + value) / 2]))
    return results


def read_quadratures(counts):
    x, y = [(plus - minus) / (plus + minus) for minus, plus in counts]
    return complex(x, y)


def calibration(sign, rng):
    rest, shifted = [Flight(20, character, 8, sign) for character in (10, 11)]
    host_tick = 1 / 32
    rest_clock = ramsey(np.vdot(rest.initial, rest.at(host_tick)), rng)
    tick = 256
    relative = np.vdot(rest.initial, rest.at(tick)).conjugate() * np.vdot(shifted.initial, shifted.at(tick))
    clock = ramsey(relative, rng)
    inverse_mass = -np.angle(read_quadratures(clock)) / (tick * (1 - np.cos(2 * np.pi / 20)))
    return dict(host_tick=host_tick, rest_quadratures=rest_clock,
                measured_rest_cost=-np.angle(read_quadratures(rest_clock)) / host_tick,
                displacement_tick=tick, displacement_quadratures=clock, inverse_mass=float(inverse_mass),
                mass=1 / inverse_mass, exact_return_visibility=float(abs(relative)))


def position_panel(size, j, sign, kappa, inverse_mass, rng):
    k0 = size // 2
    q = 2 * np.pi * j / size
    rest, moving = [Flight(size, k, kappa, sign) for k in (k0, k0 + j)]
    phase = np.exp(1j * q * np.arange(size))
    steps = []
    for t in TIMES:
        a, b = rest.at(t), moving.at(t)
        full = (a[None, :] + phase[:, None] * b[None, :]) / np.sqrt(2 * size)
        p = abs(full)**2
        counts = draws(rng, p.sum(axis=1))
        harmonic = 2 * np.vdot(phase, counts) / SHOTS
        if t == 0:
            initial_harmonic = harmonic
        ratio = harmonic / initial_harmonic
        forecast = np.exp(-1j * t * (1 - np.cos(q)) * inverse_mass)
        steps.append(dict(time=t, vacancy_position_counts=counts,
                          adjacency_probability=float(p[:, [0, -1]].sum()),
                          harmonic_ratio=[float(ratio.real), float(ratio.imag)],
                          mass_discrepancy=float(abs(ratio - forecast)), frozen_discrepancy=float(abs(ratio - 1)),
                          reversed_discrepancy=float(abs(ratio - forecast.conjugate()))))
    return dict(size=size, offset=j, orientation=sign, contrast=kappa, steps=steps)


def impulse_panel(size, stock, sign, rng):
    character = size // 2 + 1
    x = np.arange(size)[:, None]
    wave = np.exp(2j * np.pi * character * x / size) / np.sqrt(size)
    if stock == 'pair':
        wave = wave * Flight(size, character, 8, sign).at(11)[None, :]
        y = (x + np.arange(1, size)[None, :]) % size
        signed, unsigned = y - x, y + x
    else:
        charge = 1 if stock == 'excess' else -1
        wave *= np.exp(-11j * (17 - np.cos(2 * np.pi * character / size)))
        signed, unsigned = charge * x, x
    results = {}
    for label, generator in [('inactive', 0), ('signed', signed), ('unsigned', unsigned)]:
        after = wave * np.exp(2j * np.pi / size * generator)
        weights = (abs(np.fft.fft(after, axis=0))**2).sum(axis=1) / size
        results[label] = draws(rng, weights)
    return dict(size=size, stock=stock, orientation=sign, incoming=character, terminal_character_counts=results)


def report():
    rng = np.random.default_rng(190027)
    owners, rings, station = physical_rings()
    calibration_data = {str(s): calibration(s, rng) for s in (1, -1)}
    held = [position_panel(n, j, s, kappa, calibration_data[str(s)]['inverse_mass'], rng)
            for n in (30, 40) for j in (-3, -2, -1, 1, 2, 3) for s in (1, -1) for kappa in (8, 0)]
    impulse = [impulse_panel(n, stock, s, rng) for n in rings
               for stock, s in [('pair', 1), ('pair', -1), ('vacancy', 0), ('excess', 0)]]
    native = [step for row in held if row['contrast'] == 8 for step in row['steps']]
    changed = [step for row in held if row['contrast'] == 0 for step in row['steps']]
    summary = dict(max_native_mass_error=max(s['mass_discrepancy'] for s in native),
                   max_native_frozen_error=max(s['frozen_discrepancy'] for s in native),
                   max_native_reversed_error=max(s['reversed_discrepancy'] for s in native),
                   max_changed_mass_error=max(s['mass_discrepancy'] for s in changed),
                   min_native_adjacency=min(s['adjacency_probability'] for s in native),
                   min_changed_adjacency=min(s['adjacency_probability'] for s in changed))
    return dict(seed=190027, shots_per_setting=SHOTS, owner_cycle=owners, rings=rings, station=station,
                material_impulse_pulse_checks={n: signed_pulse_audit(n) for n in rings},
                calibration=calibration_data, positions=held, impulses=impulse, summary=summary,
                full_release=release_report(rings[20]['ports'], rng))


if __name__ == '__main__':
    cli = argparse.ArgumentParser()
    cli.add_argument('--output', type=Path, required=True)
    args = cli.parse_args()
    args.output.write_text(json.dumps(report(), indent=1, allow_nan=False) + '\n')
