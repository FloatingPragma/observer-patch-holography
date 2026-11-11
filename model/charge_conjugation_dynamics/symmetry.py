"""Exact local symmetry, paired native histories and finite bias bounds.

The conjugation is a symmetry of classical load/class records. It does not
identify those records with a baryon number or the cosmological baryon density.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path

from model.defect_integer_lift import model as native

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def local_step(state, coin):
    """Return both loads, both classes, current and accepted-twist flag."""
    xa, xb, k1, k2 = state
    if any(type(x) is not int or x not in range(6) for x in state) or type(coin) is not int or coin not in (0, 1):
        raise ValueError('loads/classes in 0..5 and a binary coin required')
    a, b, current = native.attempt(xa, xb, coin)
    applied = bool(current and native.accepts(k1, k2, current, xa - xb))
    return (a, b, (k1 + applied * current) % 6,
            (k2 - applied * current) % 6, current, applied)


def conjugate(state):
    xa, xb, k1, k2 = state
    return 5 - xa, 5 - xb, -k1 % 6, -k2 % 6


def local_certificate():
    cases = applied = 0
    for state in product(range(6), repeat=4):
        for coin in (0, 1):
            out = local_step(state, coin)
            mirror = local_step(conjugate(state), 1 - coin)
            expected = (*conjugate(out[:4]), -out[4], out[5])
            if mirror != expected:
                raise ArithmeticError((state, coin, out, mirror))
            cases += 1
            applied += out[5]
    return {'complete_local_branches': cases, 'accepted_nonzero_twists': applied,
            'every_branch_conjugates': True, 'no_operations_and_rejections_retained': True}


def bias_bound(initial_asymmetry, defects):
    """Bound TV(mu_n,C*mu_n) by maximal coupling at each complete step.

    defects[t] bounds sup_x TV(K_t(Cx), C*K_t(x)). All arguments are exact
    rationals. A Bernoulli(p) repair coin has defect at most abs(2*p-1).
    """
    defects = tuple(defects)
    values = (initial_asymmetry, *defects)
    if any(type(x) not in (int, F) or not 0 <= x <= 1 for x in values):
        raise ValueError('exact rational probabilities required')
    overlap = 1 - F(initial_asymmetry)
    for defect in defects:
        overlap *= 1 - defect
    return 1 - overlap


def selection_bound(asymmetry, acceptance):
    """Odd conditional mean bound for a conjugation-even selected event."""
    if type(acceptance) not in (int, F) or not 0 < acceptance <= 1:
        raise ValueError('positive exact acceptance probability required')
    return min(F(1), bias_bound(asymmetry, ()) / acceptance)


def _state(fed, loads):
    return {'loads': list(loads), 'classes': [0] * fed.faces,
            'values': [0] * fed.faces, 'link': [0] * fed.seams}


def _attempt(fed, state, seam, coin):
    a, b = fed.end_a[seam], fed.end_b[seam]
    p, m = fed.face_plus[seam], fed.face_minus[seam]
    before = (state['loads'][a], state['loads'][b], state['classes'][p], state['classes'][m])
    na, nb, kp, km, current, applied = local_step(before, coin)
    state['loads'][a], state['loads'][b] = na, nb
    state['classes'][p], state['classes'][m] = kp, km
    if applied:
        state['values'][p] += current
        state['values'][m] -= current
        state['link'][seam] += current
    return (seam, coin, current, applied)


def paired_history(level=0, sweeps=12):
    """Replay complete native schedules with conjugate preparation and coins.

    This is a pathwise certificate, with all attempted records retained in a
    digest. Each individual deterministic seed is a sample, not an ensemble.
    """
    fed = native.Federation(level)
    loads = native.prepare(fed)
    direct, mirror = _state(fed, loads), _state(fed, [5 - x for x in loads])
    digests = hashlib.sha256(), hashlib.sha256()
    counts, curve, max_abs = 0, [], 0
    for sweep, (seams, coins) in enumerate(native.schedule(fed, sweeps)):
        for seam, coin in zip(seams, coins):
            a = _attempt(fed, direct, seam, coin)
            b = _attempt(fed, mirror, seam, 1 - coin)
            if b != (a[0], 1 - a[1], -a[2], a[3]):
                raise ArithmeticError('record conjugation failed')
            for digest, record in zip(digests, (a, b)):
                digest.update((json.dumps(record, separators=(',', ':')) + '\n').encode())
            counts += 1
        for key in ('values', 'link'):
            if mirror[key] != [-x for x in direct[key]]:
                raise ArithmeticError(f'{key} symmetry failed')
        if (mirror['loads'] != [5 - x for x in direct['loads']] or
                mirror['classes'] != [-x % 6 for x in direct['classes']]):
            raise ArithmeticError('state conjugation failed')
        classes = direct['classes']
        odd = [classes.count(k) - classes.count(6 - k) for k in (1, 2)]
        other = [mirror['classes'].count(k) - mirror['classes'].count(6 - k) for k in (1, 2)]
        if other != [-x for x in odd]:
            raise ArithmeticError('odd observable failed')
        max_abs = max(max_abs, *(abs(x) for x in odd))
        curve.append({'attempts': (sweep + 1) * fed.seams,
                      'single_history_odd_counts': odd, 'paired_mean_odd_counts': [0, 0]})
    if native.curl(fed, direct['link']) != direct['values']:
        raise ArithmeticError('integer lift differs from native curl')
    return {'level': level, 'faces': fed.faces, 'attempts_per_history': counts,
            'complete_record_sha256': [x.hexdigest() for x in digests],
            'curve': curve, 'maximum_single_history_odd_count': max_abs,
            'paired_mean_zero': True, 'total_integer_lift': sum(direct['values'])}


def abundance_kernel_certificate():
    """Supply the same conjugation premise for the exact population source."""
    from model.native_abundance_epochs import model as epochs

    def flip(state):
        return (*[5 - x for x in state[:3]], -state[3] % 6)

    branches = 0
    for state in epochs.STATES:
        for seam, coin in product(range(3), range(2)):
            target, record = epochs.event(state, seam, coin)
            other, reflected = epochs.event(flip(state), seam, 1 - coin)
            expected = dict(record, coin=1 - coin, current=-record['current'])
            if other != flip(target) or reflected != expected:
                raise ArithmeticError('population-source record conjugation failed')
            branches += 1
        reflected_row = {flip(s): weight for s, weight in epochs.kernel()[state].items()}
        if reflected_row != epochs.kernel()[flip(state)]:
            raise ArithmeticError('population kernel is not equivariant')
    p = epochs.initial_distribution()
    if any(p.get(flip(s), F(0)) != p.get(s, F(0)) for s in epochs.STATES):
        raise ArithmeticError('population preparation is not even')
    return {'states': len(epochs.STATES), 'complete_record_branches': branches,
            'equivariant_kernel': True, 'even_preparation': True,
            'all_finite_times_odd_mean_zero': True}


def report():
    source_paths = ['model/defect_integer_lift/model.py', 'model/defect_integer_lift/federation.json',
                    'Lean/Dynamics/ChargeConjugation.lean', 'Lean/Dynamics/ProtectedCharge.lean',
                    'Lean/Dynamics/PublicMarkov.lean']
    source_paths += ['model/native_abundance_epochs/model.py']
    source_paths += [p.relative_to(ROOT).as_posix() for p in sorted(HERE.iterdir())
                     if p.suffix in ('.py', '.md')]
    return {'schema': 'oph.charge_conjugation_dynamics.v1',
            'sources': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in source_paths},
            'local': local_certificate(),
            'native_paired_history': paired_history(),
            'abundance_kernel': abundance_kernel_certificate(),
            'bias_examples': [{'attempts': n, 'coin_probability': '501/1000',
                               'odd_mean_upper_bound': str(bias_bound(0, [F(1, 500)] * n))}
                              for n in (1, 10, 100)],
            'scope': {'observable': 'Odd signed class counts divided by the number of faces.',
                      'symmetry': 'Load complement, class/lift sign reversal and complementary coins.',
                      'ensemble': 'Conjugation-symmetric preparation and independent unbiased coins; schedule equivariant.',
                      'clock': 'Complete native attempts, including no-operation and rejected twists.',
                      'physical': 'No baryon identification, baryon-to-photon ratio, temperature or physical duration.',
                      'formal': 'Generic finite-kernel symmetry theorem in Lean; native premises checked by exact enumeration.'}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    options = parser.add_mutually_exclusive_group()
    options.add_argument('--write', action='store_true')
    options.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = report()
    path = HERE / 'results.json'
    if args.write:
        path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    if args.check and json.loads(path.read_text()) != value:
        raise ValueError('charge-conjugation receipt differs from replay')
    print(json.dumps({'local': value['local'], 'paired_history': value['native_paired_history']}, indent=2))
