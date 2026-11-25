"""Compare the stationary Interval Echo hand with its phase-averaged reference."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import clock
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
EXPOSURE = np.pi / 16
PERIOD = 32
SHOTS = 1024
SEED = 0x50484153


def rail_certificate():
    law = clock.law
    h0, h1 = law.local_term_parts()
    h_full = h0 + clock.KAPPA * h1
    basis = np.eye(law.DIM)[:, [law.SRC.index(pair) for pair in ((2, 0), (0, 2))]]
    generator = basis.T @ h_full @ basis
    np.testing.assert_allclose(h_full @ basis, basis @ generator, atol=1e-12, rtol=0)
    np.testing.assert_allclose(generator, [[32.5, -.5], [-.5, 32.5]], atol=1e-12, rtol=0)
    energy, eigvec = np.linalg.eigh(generator)
    start = np.diag([1., 0.])
    pulse = clock._dimer_pulse(EXPOSURE)
    orbit = []
    state = start
    for _ in range(PERIOD):
        orbit.append(state)
        state = pulse @ state @ pulse.conj().T
    np.testing.assert_allclose(state, start, atol=1e-12, rtol=0)
    reference = np.mean(orbit, axis=0)
    np.testing.assert_allclose(reference, np.eye(2) / 2, atol=1e-12, rtol=0)
    # Shift the upper energy by one sampling frequency. This compares
    # generators of the hand, without changing the native sampling apparatus.
    alias = (eigvec * np.exp(-1j * EXPOSURE * (energy + [0, PERIOD]))) @ eigvec.T
    np.testing.assert_allclose(alias, pulse, atol=1e-12, rtol=0)
    contrast = []
    for upper_shift in (0, PERIOD):
        half = (eigvec * np.exp(-.5j * EXPOSURE * (energy + [0, upper_shift]))) @ eigvec.T
        contrast.append(float(abs(half[0, 0])**2))
    return {
        'rail_generator': generator.real.round(12).tolist(),
        'uniform_phase_reference': reference.real.round(12).tolist(),
        'native_and_alias_gaps': [1, 1 + PERIOD],
        'half_tick_first_rail_probabilities_native_alias': np.round(contrast, 12).tolist(),
    }


def acquire(car, services, analyzer):
    # Each separated bank receives an independent stream. A trial starts
    # with a fresh sharp hand, the same public digit and the native fair pair.
    counts = Counter()
    for launch in range(SHOTS):
        registers = clock.Registers(car, [SEED, services, analyzer, launch])
        owner = {0, 1}
        registers.prepare_front(0)
        hand = registers.prepare_oscillator(1, owner)
        digit = clock._digits(1, clock.EXTRA_DIGITS[1:])[0]
        registers.initialize_digits([digit])
        accepted = 0
        for _ in range(services):
            registers.publish_tick(digit, registers.fair(0, owner), owner)
            accepted += registers.read_tick(digit, owner)
            hand.on_tick(registers, digit, owner)
        detector = hand.terminal(registers, owner, analyzer)['bit']
        counts[accepted, detector] += 1
    return [[n, counts[n, 0], counts[n, 1]] for n in range(services + 1)]


def summarize(bank, services, analyzer):
    shift = analyzer * np.pi / 2
    shots = sum(no + yes for _, no, yes in bank)
    measured = sum(yes for _, _, yes in bank) / shots
    predicted = sum((no + yes) * (1 + np.cos(EXPOSURE * n + shift)) / 2
                    for n, no, yes in bank) / shots
    # On each opportunity the fair source produces either identity or one pulse.
    visibility = np.cos(EXPOSURE / 2)**services
    unread = (1 + visibility * np.cos(services * EXPOSURE / 2 + shift)) / 2
    # Four banks x two comparisons, family-wise false rejection <= 1%.
    limit = np.sqrt(np.log(2 * 8 / .01) / (2 * shots))
    if shots != SHOTS or max(abs(measured - predicted), abs(measured - unread)) > limit:
        raise ValueError('registered observer bank differs from the fair-source prediction')
    return {'opportunities': services, 'analyzer': analyzer,
            'by_accepted_count_misses_hits': bank, 'observed_first_rail': measured,
            'conditional_prediction_averaged_over_records': round(float(predicted), 12),
            'unread_opportunity_prediction': round(float(unread), 12),
            'unread_visibility': round(float(visibility), 12)}


def measure():
    car = federation(0)
    cert = rail_certificate()
    banks = [summarize(acquire(car, m, a), m, a) for m in (32, 64) for a in (0, 1)]
    if banks[0]['observed_first_rail'] >= .5 or banks[2]['observed_first_rail'] <= .5:
        raise ValueError('first-rail measurements failed to resolve the hand reversal')
    return {'candidate': 'interval_echo', 'board_cells': car.cells, 'seed': SEED,
            'shots_per_bank': SHOTS, 'rail_certificate': cert, 'banks': banks}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='write a fresh record at this path')
    args = parser.parse_args()
    receipt = measure()
    if args.output:
        with args.output.open('x') as output:
            json.dump(receipt, output, indent=2, allow_nan=False)
            output.write('\n')
    elif receipt != json.loads((HERE / 'records.json').read_text()):
        raise SystemExit('observer bank does not agree with the retained record')
    print([(b['opportunities'], b['analyzer'], b['observed_first_rail']) for b in receipt['banks']])
