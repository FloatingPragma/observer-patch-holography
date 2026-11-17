"""Infer a nonuniform IID seam scheduler from visible pair-clock laws.

The estimator reads mean count ratios and terminal coincidences, never a
hidden attempt counter, source seam table, load preparation or coordinates.
Its output predicts marked field laws under a separately declared stationary
binary preparation. Exact ensemble identification does not certify IID
scheduling or identify a positive edge from finitely many absent events.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations
from pathlib import Path

import sympy as sp

from ..laboratory_clock import event_geometry as event
from ..observer import native
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*event.SOURCE_PATHS,
    'model/relative_observer_geometry/clock_reconstruction.py',
    'model/relative_observer_geometry/__init__.py'))))
PORTS = 12


def _exact(token, *, positive=False):
    if type(token) is not str:
        raise ValueError('a canonical exact probability or mean string is required')
    try:
        value = F(token)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError('a canonical rational is required') from error
    if str(value) != token or value < 0 or (positive and not value):
        raise ValueError('a nonnegative canonical rational is required')
    return value


def _weights(weights):
    if type(weights) is not dict or not weights:
        raise ValueError('a normalized seam probability dictionary is required')
    for seam, probability in weights.items():
        if (type(seam) is not tuple or len(seam) != 2
                or any(type(x) is not int or x not in range(PORTS) for x in seam)
                or seam[0] >= seam[1] or type(probability) is not F
                or _exact(str(probability), positive=True) != probability):
            raise ValueError('increasing native pairs and positive exact weights required')
    if sum(weights.values(), F(0)) != 1:
        raise ValueError('seam probabilities must sum to one')
    return dict(sorted(weights.items()))


def clock_records(weights):
    """Source-side generator; only its pair-law records reach recover()."""
    weights = _weights(weights)
    exposure = [sum((p for seam, p in weights.items() if i in seam), F(0))
                for i in range(PORTS)]
    if not all(exposure):
        raise ValueError('every declared clock must tick with positive probability')
    return [{'ports': [a, b], 'mean_B_per_A': str(exposure[b]/exposure[a]),
             'terminal_coincidence': str(weights.get(tuple(sorted((a, b))), F(0))/exposure[a])}
            for a in range(PORTS) for b in range(PORTS) if a != b]


def recover(records):
    """Unique IID single-seam law consistent with all oriented clock pairs.

    Each event ticks exactly its two endpoints, hence sum_i exposure_i=2.
    This fixes probability normalization; it introduces no physical time unit.
    """
    if type(records) not in (list, tuple) or len(records) != PORTS*(PORTS-1):
        raise ValueError('all 132 oriented pair-clock laws are required')
    pairs = {}
    for row in records:
        if type(row) is not dict or set(row) != {'ports', 'mean_B_per_A', 'terminal_coincidence'}:
            raise ValueError('only pair labels, mean count and coincidence are admitted')
        ports = row['ports']
        if (type(ports) not in (list, tuple) or len(ports) != 2
                or any(type(x) is not int or x not in range(PORTS) for x in ports)
                or ports[0] == ports[1] or tuple(ports) in pairs):
            raise ValueError('unique oriented pairs of distinct native ports required')
        mean = _exact(row['mean_B_per_A'], positive=True)
        coincidence = _exact(row['terminal_coincidence'])
        if coincidence > 1:
            raise ValueError('coincidence must be a probability')
        pairs[tuple(ports)] = mean, coincidence
    relative = [F(1)]+[pairs[0, i][0] for i in range(1, PORTS)]
    exposure = tuple(2*x/sum(relative, F(0)) for x in relative)
    weights = {}
    for a, b in combinations(range(PORTS), 2):
        forward = exposure[a]*pairs[a, b][1]
        backward = exposure[b]*pairs[b, a][1]
        if forward != backward:
            raise ValueError('opposite coincidence laws disagree on the shared event')
        if forward:
            weights[a, b] = forward
    for (a, b), (mean, _) in pairs.items():
        if mean != exposure[b]/exposure[a]:
            raise ValueError('mean clock ratios have no common exposure vector')
    if any(sum((p for seam, p in weights.items() if i in seam), F(0)) != exposure[i]
           for i in range(PORTS)):
        raise ValueError('coincidence rows do not account for every local tick')
    weights = _weights(weights)
    reached = {0}
    while True:
        enlarged = reached | {v for edge in weights if any(u in reached for u in edge) for v in edge}
        if enlarged == reached:
            break
        reached = enlarged
    if len(reached) != PORTS:
        raise ValueError('the recovered positive support is disconnected')
    return weights, exposure


def structural_gram(records):
    """Recover unweighted structural geometry separately from weighted kinetics."""
    weights, _ = recover(records)
    laplacian = [[F(0) for _ in range(PORTS)] for _ in range(PORTS)]
    for a, b in weights:
        laplacian[a][a] += 1; laplacian[b][b] += 1
        laplacian[a][b] -= 1; laplacian[b][a] -= 1
    checks, projector = event.geometry._band_certificate(laplacian)
    if not all(checks.values()) or any(laplacian[i][i] != 5 for i in range(PORTS)):
        raise ValueError('the inferred graph lacks the declared native rank-three certificate')
    return event.geometry._encode_q5([
        [event.geometry._qscale(x, 4) for x in row] for row in projector])


def _arguments(a, b, z, ticks):
    event._ports(a, b)
    event._integer(ticks, 1, 8)
    return event._argument(z)


def renewal_pgf(records, a, b, z=F(1), ticks=1):
    z = _arguments(a, b, z, ticks)
    weights, exposure = recover(records)
    c = weights.get(tuple(sorted((a, b))), F(0))
    x, y = exposure[a], exposure[b]
    return ((x-c+c*z)/(x+y-c-(y-c)*z))**ticks


def predict_field(records, a, b, z=F(1), ticks=2):
    """Predict the complete initial-A/final-B bit PGF using inferred weights."""
    z = _arguments(a, b, z, ticks)
    weights, _ = recover(records)
    remote, incident = sp.zeros(PORTS), sp.zeros(PORTS)
    for (x, y), probability in weights.items():
        v = sp.zeros(PORTS, 1); v[x] = 1; v[y] = -1
        step = sp.Rational(probability)*sp.Rational(z)**int(b in (x, y))*(sp.eye(PORTS)-v*v.T/2)
        if a in (x, y):
            incident += step
        else:
            remote += step
    stopped = (incident*(sp.eye(PORTS)-remote).inv(method='DM'))**ticks
    mass = renewal_pgf(records, a, b, z, ticks)
    if any(F(x) != mass for x in sp.ones(1, PORTS)*stopped):
        raise ArithmeticError('count generating mass disagrees with marked density transfer')
    joint = mass/36+F(event.CHI*(stopped*event.P0)[b, a])
    values = {(1, 1): joint, (1, 0): mass/6-joint,
              (0, 1): mass/6-joint, (0, 0): 2*mass/3+joint}
    if any(p < 0 for p in values.values()) or sum(values.values()) != mass:
        raise ArithmeticError('field generating weights lost positivity or total mass')
    return values


@lru_cache(maxsize=8)
def _native_field(items, a, b, z_pair, ticks):
    """Independent 66-state native endpoint oracle, with all waiting tails."""
    z = F(*z_pair)
    states = tuple(native.states(binary=True))
    index = {state: i for i, state in enumerate(states)}
    remote, incident = sp.zeros(66), sp.zeros(66)
    for column, state in enumerate(states):
        for seam, numerator, denominator in items:
            probability = F(numerator, denominator)
            target = incident if a in seam else remote
            for coin in (0, 1):
                after, (_, current) = native.native(state, seam, coin)
                for port in (a, b):
                    if after[port] != state[port]+current*(int(port == seam[1])-int(port == seam[0])):
                        raise ArithmeticError('local current readback failed')
                target[index[after], column] += sp.Rational(probability/2)*sp.Rational(z)**int(b in seam)
    population = sp.Matrix(66, 2, lambda i, bit: sp.Rational(int(states[i][a] == bit), 66))
    inverse = (sp.eye(66)-remote).inv(method='DM')
    for _ in range(ticks):
        population = incident*inverse*population
    result = []
    for first in (0, 1):
        for last in (0, 1):
            value = F(sum(population[i, first] for i in range(66) if states[i][b] == last))
            result.append(((first, last), value.numerator, value.denominator))
    return tuple(result)


def native_field(weights, a, b, z=F(1), ticks=2):
    z = _arguments(a, b, z, ticks)
    items = tuple((edge, p.numerator, p.denominator) for edge, p in _weights(weights).items())
    return {key: F(numerator, denominator)
            for key, numerator, denominator in _native_field(items, a, b, (z.numerator, z.denominator), ticks)}


def periodic_control(seams, a=0, b=1):
    """Stationary cyclic schedule: identical pair means, different interval law.

    A uniform initial phase is conditioned on an A tick (its Palm law).
    Each seam appears once per period. This is a dependent scheduling control,
    not an alternative IID generator admitted by recover().
    """
    event._ports(a, b)
    if len(set(seams)) != len(seams) or not all(any(p in e for e in seams) for p in range(PORTS)):
        raise ValueError('a complete positive, unrepeated seam period is required')
    intervals = []
    for start, seam in enumerate(seams):
        if a not in seam:
            continue
        count = 0
        for offset in range(1, len(seams)+1):
            edge = seams[(start+offset) % len(seams)]
            count += int(b in edge)
            if a in edge:
                intervals.append((count, int(b in edge)))
                break
    return {'mean_B_per_A': F(sum(x for x, _ in intervals), len(intervals)),
            'terminal_coincidence': F(sum(c for _, c in intervals), len(intervals)),
            'zero_B_count_probability': F(sum(x == 0 for x, _ in intervals), len(intervals))}


def run():
    seams = tuple(sorted(tuple(sorted(e)) for e in event._source().seams))
    uniform = {edge: F(1, len(seams)) for edge in seams}
    biased = {edge: F(i+1, 465) for i, edge in enumerate(seams)}
    rows = []
    for label, weights in (('uniform', uniform), ('biased', biased)):
        records = clock_records(weights)
        recovered, exposure = recover(records)
        if recovered != weights or structural_gram(records) != structural_gram(clock_records(uniform)):
            raise ArithmeticError('clock inference failed scheduler recovery or structural invariance')
        predictions = []
        for z, ticks in ((F(1, 2), 1), (F(1), 2)):
            observed = native_field(weights, 0, 1, z, ticks)
            forecast = predict_field(records, 0, 1, z, ticks)
            if observed != forecast:
                raise ArithmeticError('weighted field prediction fails independent native execution')
            predictions.append({'A_ticks': ticks, 'B_count_PGF_argument': str(z),
                'joint_bit_generating_weights': [
                    {'initial_A': i, 'final_B': j, 'weight': str(p)}
                    for (i, j), p in sorted(forecast.items())]})
        rows.append({'schedule': label, 'seam_weights': [
            {'ports': list(e), 'probability': str(p)} for e, p in recovered.items()],
            'incident_exposures': list(map(str, exposure)), 'oriented_pair_laws': len(records),
            'predictions_at_ports_0_1': predictions})
    changed = predict_field(clock_records(biased), 0, 1) != predict_field(clock_records(uniform), 0, 1)
    if not changed:
        raise ArithmeticError('the biased scheduler failed to change the field prediction')
    cyclic = periodic_control(seams)
    uniform_pair = clock_records(uniform)[0]
    if any(str(cyclic[key]) != uniform_pair[key] for key in ('mean_B_per_A', 'terminal_coincidence')):
        raise ArithmeticError('the dependent control changed its matched pair statistics')
    iid_zero = renewal_pgf(clock_records(uniform), 0, 1, F(0))
    if cyclic['zero_B_count_probability'] == iid_zero:
        raise ArithmeticError('the temporal-dependence control failed to distinguish interval laws')
    return {'schema': 'oph.clock_scheduler_reconstruction.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'normalization': 'each single-seam event ticks exactly two ports; sum of exposures is two',
        'reconstruction': 'a_i=2*mu_0i/sum_j(mu_0j), mu_00=1; p_ij=a_i*coincidence_ij',
        'same_structural_gram_different_kinetics': changed,
        'independent_full_state_comparisons': 4, 'cases': rows,
        'dependent_schedule_control': {
            'cyclic_Palm_statistics': {key: str(value) for key, value in cyclic.items()},
            'IID_zero_B_count_probability': str(iid_zero),
            'pair_means_and_terminal_coincidences_do_not_establish_IID': True},
        'scope': {'exact_ensemble_laws': True, 'uniform_scheduler_required': False,
            'IID_exogenous_scheduler_required': True, 'positive_connected_support_required': True,
            'hidden_attempt_counter_read': False, 'clock_geometry_requires_load_preparation': False,
            'field_prediction_requires_stationary_binary_total_two': True,
            'finite_sample_support_identification_claimed': False,
            'temporal_independence_inferred_from_pair_means': False,
            'physical_spacetime_identification_claimed': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Weighted clock reconstruction and native field replay verified.')
    else:
        print(canonical(result))


if __name__ == '__main__':
    main()
