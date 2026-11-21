"""Compare an independent coherent seam hand to Interval Echo at guided reunion."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import clock, cohorts
from model.executable_model import federation
from model.executable_model.carrier import seam_digest

HERE = Path(__file__).resolve().parent
SEED = 2026112200
PATTERNS = ((2, 6), (6, 2), (4, 4))
RULES = ("positive_primary", "accepted_tick")
BLOCKS, BLOCK_SIZE = 32, 32768
PULSE = np.pi / 16
TOLERANCE = 0.006


def native_probe():
    """Read the first rail against an empty neighbour with the native agreement pointer."""
    law = clock.law
    h0, h1 = law.local_term_parts()
    h = h0 + clock.KAPPA * h1
    rails = [law.SRC.index((2, 0)), law.SRC.index((0, 2))]
    basis = np.eye(law.DIM)[:, rails]
    block = basis.T @ h @ basis
    if np.linalg.norm(h @ basis - basis @ block) > 1e-11:
        raise ArithmeticError("the dimer leaves the two coherent rails")
    pointer = np.arange(law.DIM) % len(law.SRC) == law.SRC.index((1, 1))
    # Read the first rail together with its own sharp vacant neighbour.
    occupancy = [float(np.sum(np.abs(law.echo()[pointer, law.SRC.index((load, 0))]) ** 2)) for load in (2, 0)]
    if not np.allclose(occupancy, (1, 0), rtol=0, atol=1e-12):
        raise ArithmeticError("the signed read does not distinguish the rails")
    transfer = expm(-1j * np.pi * h)
    for a in (0, 2):
        for b in (0, 2):
            column = transfer[:, law.SRC.index((a, b))]
            expected = np.eye(law.DIM)[:, law.SRC.index((b, a))]
            if np.linalg.norm(column - expected) > 1e-11:
                raise ArithmeticError("the transfer changes a relative phase")

    tables = []
    for analyzer in (0, np.pi / 2):
        probabilities = []
        for marks in range(1 + max(map(sum, PATTERNS))):
            state = expm(-1j * (marks * PULSE + analyzer) * block)[:, 0]
            probabilities.append(np.clip(occupancy @ (np.abs(state) ** 2), 0, 1))
        tables.append(probabilities)
    tables = np.array(tables)
    zero_phase = np.angle((2 * tables[0, 0] - 1) + 1j * (1 - 2 * tables[1, 0]))
    if abs(zero_phase) > 1e-12 or abs(tables[0, 0] - 1) > 1e-12:
        raise ArithmeticError("coupling-off or rest preparation failed")
    return tables


class RailLabels(clock.Registers):
    """Apply the registered physical handoff sequence to distinct transport labels."""

    def __init__(self, car):
        super().__init__(car, seed=0)
        self.loads = np.arange(car.n)

    def position(self, labels):
        return [int(np.flatnonzero(self.loads == label).item()) for label in labels]


def copy_at_reunion(reg, shared, moving, home, donor):
    """Exercise both physically possible terminal bits at the observers' overlap."""
    owner = moving if donor == "moving" else home
    seam = (12 * shared + 8, 12 * shared + 9, 12 * shared + 6)
    reg.initialize_digits([seam])
    result = []
    for value in (0, 1):
        outcomes = reg.signed_record(seam, value, owner, (moving, home))
        if outcomes != [value, value]:
            raise ArithmeticError("reunion owners disagree on the dimer pointer")
        result.append(outcomes)
    return result


def witness(car, pattern, seed, probe):
    record = clock.material_replay(car, pattern, age_bits=3, seed=seed)
    rail = record['guide']
    home = set(record['home_private']) | set(rail[-2:])
    moving = set(rail[:-1])
    bank, parked = rail[-3], record['home_private'][-1]
    tips = (12 * bank + 6, 12 * bank + 7, 12 * bank + 10)
    routed = RailLabels(car)
    for c, owner in ((bank, moving), (parked, home)):
        routed.native(12 * c + 6, 12 * c + 7, owner)
        routed.native(12 * c + 6, 12 * c + 10, owner)
    position, age = 0, 0
    chronological = []
    for event in record['local_events']:
        age += 1
        if event['handoff']:
            routed.shift(rail, right=(position == 0))
            position = 1 - position
            expected = [12 * rail[-3 + position] + k for k in (6, 7, 10)]
            if routed.position(tips) != expected:
                raise ArithmeticError("coherent rails lost or interchanged during transfer")
            age = 0
        if event['age'] != age:
            raise ArithmeticError("inconsistent stored interval")
        chronological.append([int(event[k]) for k in ('handoff', 'home_primary', 'primary', 'tick')])
    if position or routed.position(tips) != list(tips):
        raise ArithmeticError("dimer failed to return to its original patch")
    if sum(row[3] for row in chronological) != record['reunion']['traveller']:
        raise ArithmeticError("the material hand and the local accepted history disagree")
    positive_home = sum(row[1] for row in chronological)
    positive_moving = sum(row[2] for row in chronological)
    posted_moving = sum(row[3] for row in chronological)
    # Transfers have carried both columns with the same phase; the local seams
    # between transfers therefore compose at their actual accumulated exposures.
    fields = {'home': positive_home, 'moving_primary': positive_moving, 'moving_accepted': posted_moving}
    terminal = {key: [round(float(p), 12) for p in probe[:, marks]] for key, marks in fields.items()}
    receipts = {donor: copy_at_reunion(clock.Registers(car, seed), record['shared_cell'], moving, home, donor)
                for donor in ('home', 'moving')}
    return {'seed': seed, 'home_probe_cell': parked, 'moving_probe_birth_cell': bank,
            'moving_probe_away_cell': rail[-2], 'reunion_cell': record['shared_cell'],
            'native_transfer_contacts': routed.gates['pi_contacts'],
            'two_owner_readings_by_terminal_pointer': receipts, 'local_events': chronological,
            'terminal_first_rail_probabilities': terminal, 'published_words': record['reunion']}


def joint_marks(handoffs):
    """Probability of retained raw and accepted words after the measured handoffs."""
    distribution = {(0, 0): 1.}
    age = 0
    for moved in handoffs:
        age += 1
        suppressed = cohorts.suppression(age, 3) if moved else 0.
        if moved:
            other = clock.interval_response(age, age_bits=3)['closing_candidate_suppression']
            if abs(suppressed - other) > 1e-13:
                raise ArithmeticError("finite first-return sources disagree")
            age = 0
        following = Counter()
        for (raw, accepted), probability in distribution.items():
            following[raw, accepted] += probability / 2
            following[raw + 1, accepted] += probability * suppressed / 2
            following[raw + 1, accepted + 1] += probability * (1 - suppressed) / 2
        distribution = following
    if age or abs(sum(distribution.values()) - 1) > 1e-13:
        raise ArithmeticError("incomplete itinerary or probability")
    support = sorted(distribution)
    return np.array(support), np.array([distribution[k] for k in support])


def read_ratios(mean):
    """Input comprises terminal observer words and occupancy bits for both quadratures."""
    hands = mean[:, 1].sum() / mean[:, 0].sum()
    waves = (2 * mean[0, 2:] - 1) + 1j * (1 - 2 * mean[1, 2:])
    home, moving = np.angle(waves)
    if home <= 0 or min(abs(waves)) < .9:
        raise ArithmeticError("stationary phase or visibility control failed")
    dimer = moving / home
    return np.array([hands, dimer, dimer - hands, home, moving, *np.abs(waves)])


def cohort(pattern_index, rule_index, handoffs, probe, *, seed=SEED, block_size=BLOCK_SIZE):
    support, weights = joint_marks(handoffs)
    sums = np.empty((2, BLOCKS, 4), dtype=np.int64)
    for setting in range(2):
        rng = np.random.default_rng(np.random.SeedSequence([seed, 777, pattern_index, rule_index, setting]))
        for b in range(BLOCKS):
            raw, accepted = support[rng.choice(len(support), size=block_size, p=weights)].T
            home = rng.binomial(len(handoffs), .5, size=block_size)
            moving = (raw, accepted)[rule_index]
            home_bit = rng.random(block_size) < probe[setting, home]
            moving_bit = rng.random(block_size) < probe[setting, moving]
            sums[setting, b] = [home.sum(), accepted.sum(), home_bit.sum(), moving_bit.sum()]
    means = sums.mean(axis=1) / block_size
    full = read_ratios(means)
    omitted = (sums.sum(axis=1)[:, None] - sums) / ((BLOCKS - 1) * block_size)
    leave_one = np.array([read_ratios(omitted[:, b]) for b in range(BLOCKS)])
    errors = np.sqrt((BLOCKS - 1) * ((leave_one - leave_one.mean(0))**2).mean(0))
    delta, uncertainty = full[2], 4 * errors[2]
    if abs(delta) + uncertainty <= TOLERANCE:
        verdict = "agrees_within_tolerance"
    elif abs(delta) - uncertainty > TOLERANCE:
        verdict = "disagrees_beyond_tolerance"
    else:
        verdict = "unresolved"
    keys = ('registered_ratio', 'dimer_ratio', 'dimer_minus_registered', 'home_phase', 'moving_phase',
            'home_visibility', 'moving_visibility')
    return {'readings': {k: {'value': round(float(v), 10), 'se': round(float(e), 10)}
                         for k, v, e in zip(keys, full, errors, strict=True)},
            'verdict': verdict, 'block_sums': sums.tolist()}


def material_ensemble(pattern, seed, probe):
    boards = []
    sequence = None
    for level in range(3):
        car = federation(level)
        rows = [witness(car, pattern, seed + 100 * level + j, probe) for j in range(5)]
        for row in rows:
            observed = tuple(e[0] for e in row['local_events'])
            if sequence is not None and sequence != observed:
                raise ArithmeticError("cannot pool different material traces")
            sequence = observed
        boards.append({'cells': car.cells, 'seams_sha256': seam_digest(car), 'witnesses': rows})
    return sequence, boards


def run():
    protocol_hash = sha256((HERE / 'PROTOCOL.md').read_bytes()).hexdigest()
    probe = native_probe()
    studies = []
    for i, pattern in enumerate(PATTERNS):
        sequence, boards = material_ensemble(pattern, SEED + 10 * i, probe)
        groups = {rule: cohort(i, ri, sequence, probe) for ri, rule in enumerate(RULES)}
        studies.append({'prepared_services': list(pattern), 'boards': boards, 'readouts': groups})
    return {'schema': 'oph.interval_echo.coincident_seam_hands.v1', 'protocol_sha256': protocol_hash,
            'preparations_per_setting': BLOCKS * BLOCK_SIZE, 'dimer_pulse': 'pi/16',
            'block_fields': ['home_elapsed_word', 'traveller_elapsed_word', 'home_first_rail', 'moving_first_rail'],
            'local_event_fields': ['handoff', 'home_primary', 'traveller_primary', 'traveller_posted'],
            'analyzer_settings': ['0', 'pi/2'], 'studies': studies}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='write a fresh receipt')
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open('x') as f:
            f.write(json.dumps(result, separators=(',', ':'), allow_nan=False) + '\n')
    elif json.loads((HERE / 'results.json').read_text()) != result:
        raise SystemExit('receipt differs from fresh run')
    for row in result['studies']:
        print(row['prepared_services'], [(rule, arm['verdict'], arm['readings']['dimer_minus_registered'])
                                         for rule, arm in row['readouts'].items()])


if __name__ == '__main__':
    main()
