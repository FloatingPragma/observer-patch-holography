"""Read the reference component invisible to the two terminal clock settings."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import clock
from model.executable_model import federation
from .measure import EXPOSURE, HERE, PERIOD

SHOTS = 4096
PORTS = (18, 19, 20)
VACANCY = 21


def contact(car, first, second, duration):
    """A protected contact with one double, blank memories and a spectator rail."""
    if PORTS[second] not in car.neighbours[PORTS[first]]:
        raise ValueError('unavailable material contact')
    law = clock.law
    h0, h1 = law.local_term_parts()
    h = h0 + clock.KAPPA * h1
    # The double outside the two addressed rails presents (0,0) to this seam.
    states = [(2 * int(i == first), 2 * int(i == second)) for i in range(3)]
    basis = np.eye(law.DIM)[:, [law.SRC.index(pair) for pair in states]]
    reduced = basis.T @ h @ basis
    np.testing.assert_allclose(h @ basis, basis @ reduced, atol=1e-12, rtol=0)
    return expm(-1j * duration * reduced)


def pointer(car, rail, vacancy):
    """Compute the local yes effect for a sharp read against a vacant neighbour."""
    if vacancy in PORTS or vacancy not in car.neighbours[PORTS[rail]]:
        raise ValueError('the pointer needs a neighbouring vacant rail')
    law = clock.law
    echo = law.echo()
    balanced = np.arange(law.DIM) % len(law.SRC) == law.SRC.index((1, 1))
    columns = [law.SRC.index((load, 0)) for load in (0, 2)]
    yes = np.sum(abs(echo[balanced][:, columns])**2, axis=0)
    np.testing.assert_allclose(yes, [0, 1], atol=1e-12, rtol=0)
    # The spectator double positions are orthogonal and survive the read.
    return np.diag([yes[int(i == rail)] for i in range(3)])


def instruments(car):
    half_a = contact(car, 0, 2, np.pi / 2)
    half_b = contact(car, 1, 2, np.pi / 2)
    # Launch from the third rail, visit the second and then the first rail.
    raw = (half_a @ half_b)[:, 2]
    third_occupied = pointer(car, 2, VACANCY)
    ready = float(1 - np.vdot(raw, third_occupied @ raw).real)
    np.testing.assert_allclose(ready, 3 / 4, atol=1e-12, rtol=0)
    launches = [(1., np.array([1., 0., 0.])), (ready, np.array([*raw[:2], 0.]) / ready**.5)]
    first_occupied = pointer(car, 0, 22)
    settings = [(first_occupied, np.eye(3)), (first_occupied, contact(car, 0, 1, np.pi / 2)),
                (third_occupied, half_a @ half_b)]
    pulse = contact(car, 0, 1, EXPOSURE)
    return launches, settings, pulse


def run():
    car = federation(0)
    launches, settings, pulse = instruments(car)
    rows = []
    for prep, (ready, initial) in enumerate(launches):
        phases = [np.linalg.matrix_power(pulse, n) @ initial for n in range(PERIOD)]
        rho = sum(np.outer(state[:2], state[:2].conj()) for state in phases) / PERIOD
        target = np.eye(2) / 2 + prep * np.array([[0, 1/3], [1/3, 0]])
        np.testing.assert_allclose(rho, target, atol=1e-12, rtol=0)
        if prep:
            weights, axes = np.linalg.eigh(rho)
            np.testing.assert_allclose(weights, [1/6, 5/6], atol=1e-12, rtol=0)
            # Check equality up to the central phase on the whole rail algebra.
            modular = (axes * np.exp(1j * np.log(weights) * EXPOSURE / np.log(5))) @ axes.conj().T
            relative = modular.conj().T @ pulse[:2, :2]
            np.testing.assert_allclose(relative, np.eye(2) * relative[0, 0], atol=1e-12, rtol=0)
        readings = []
        for analyzer, (effect, before) in enumerate(settings):
            output = [before @ state for state in phases]
            p = [float(np.vdot(state, effect @ state).real) for state in output]
            theory = Fraction(1, 2) if analyzer < 2 else [Fraction(3, 8), Fraction(13, 24)][prep]
            np.testing.assert_allclose(np.mean(p), float(theory), atol=1e-12, rtol=0)
            registers = clock.Registers(car, [0x524546, prep, analyzer])
            registers.prepare_front(0)
            failed = 0
            # Retain every supplied phase and outcome, including failed heralds.
            phase_counts = np.zeros((PERIOD, 2), dtype=int)
            for _ in range(SHOTS):
                if registers.rng.random() >= ready:
                    failed += 1
                    continue
                n = sum(registers.fair(0, {0}) << i for i in range(5))
                # Local sharp double-vacancy read after coherent contact(s).
                hit = int(registers.rng.random() < p[n])
                phase_counts[n, hit] += 1
            measured = int(phase_counts[:, 1].sum()) / (SHOTS - failed)
            # Six conditional readouts and six readiness frequencies.
            def limit(n):
                return np.sqrt(np.log(2 * 12 / .01) / (2 * n))
            if abs(measured - float(theory)) > limit(SHOTS - failed) or abs(1 - failed / SHOTS - ready) > limit(SHOTS):
                raise ArithmeticError('herald or local detector frequency missed its prediction')
            readings.append({'analyzer': analyzer, 'expected_first_or_third': str(theory),
                             'failed_heralds': failed, 'phase_by_miss_hit': phase_counts.tolist()})
        rows.append({'preparation': ['sharp first rail', 'heralded from third rail'][prep],
                     'success_probability': ['1', '3/4'][prep], 'settings': readings})
    return {'board_cells': car.cells, 'rails': list(PORTS), 'third_rail_read_vacancy': VACANCY,
            'shots_per_setting': SHOTS, 'records': rows}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='write a fresh record at this path')
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open('x') as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write('\n')
    elif result != json.loads((HERE / 'reference_records.json').read_text()):
        raise SystemExit('reference readouts differ from retained counts')
    for prep in result['records']:
        print(prep['preparation'], [(s['failed_heralds'], sum(row[1] for row in s['phase_by_miss_hit']))
                                    for s in prep['settings']])
