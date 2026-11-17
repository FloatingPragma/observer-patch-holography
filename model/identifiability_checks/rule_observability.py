"""Exact C/E observability under native-record-controlled quantum work.

Let D_N retain the diagonal native blocks of a joint operator on N x W x R,
without measuring or discarding W or the untouched reference R. For every
seam and common record r=(coin,current), the actual partial injections give

    D_N I_C,r = D_N I_E,r = I_C,r D_N = I_E,r D_N.

Here I_E,r is a sum of CP maps, not a coherent sum of contrast Kraus maps.
Partial injectivity is essential: an off-diagonal native matrix unit cannot
become a diagonal output inside a recorded branch. This file checks every
matrix unit on all thirty seams of the full 78-state carrier.

An arbitrary work-only CP instrument commutes with D_N. Induction therefore
gives equal D_N blocks at every common history, even for correlated quantum
inputs, work measurements, record-dependent seams and bounded stopping.
Tracing N gives identical complete records AND conditional W,R operators.
Randomized policies are convex mixtures of the same identity. The statement
is dimension-independent in W,R, does not require IID trials or a diagonal
native preparation, and concerns any finite common-access protocol.

The access contract excludes additional coherent native controls or recovery,
coherent N-W gates, native coherent readout, and feedback on E's extra
contrast. The complete native quantum outputs can differ. An existing
interference witness demonstrates exactly that boundary.

This specializes quantum-instrument composition and no-signalling to the
actual C/E records; ordinary trace preservation alone would not suffice in
the presence of record-conditioned work. It removes an independence premise
from this observational equivalence, not a physical apparatus-identification
premise. The empirical oscillator apparatus is separate from these native
records. No laboratory action, noise channel or probability map is added.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from model.observer import native, quantum_mixing
from model.predictions import observer_instrument as witness
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
STATES = native.states()
INDEX = {state: i for i, state in enumerate(STATES)}
RECORDS = tuple((coin, current) for coin in (0, 1) for current in (-1, 0, 1))
ZERO = (F(0), F(0))
ONE = (F(1), F(0))
I = (F(0), F(1))
SOURCE_PATHS = ('model/observer/native.py', 'model/observer/quantum_mixing.py',
                'model/predictions/observer_instrument.py',
                'model/identifiability_checks/rule_observability.py')


def _add(target, key, value):
    old = target.get(key, ZERO)
    new = (old[0]+value[0], old[1]+value[1])
    if new == ZERO:
        target.pop(key, None)
    else:
        target[key] = new


def _multiply(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def _conjugate(a):
    return (a[0], -a[1])


def _scale(a, value):
    return (a[0]*value, a[1]*value)


def _common(law, raw):
    coin, value = raw
    current = value if law == 'C' else (value-(2*coin-1)*(abs(value) % 2))//2
    return coin, current


def columns(law, seam):
    """Actual full-carrier native columns; each coin has squared weight 1/2."""
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('candidate C or E required')
    if (type(seam) is not tuple or len(seam) != 2
            or any(type(v) is not int for v in seam)
            or seam not in quantum_mixing.seams()):
        raise ValueError('an oriented seam from the declared carrier is required')
    result = []
    for state in STATES:
        row = []
        for coin in (0, 1):
            output, current = native.native(state, seam, coin)
            raw = current if law == 'C' else (coin, state[seam[0]]-state[seam[1]])
            if _common(law, raw) != current:
                raise ArithmeticError('contrast does not recover the same native current')
            row.append((raw, INDEX[output]))
        result.append(tuple(row))
    return tuple(result)


def _diagonal_image(law, first, second):
    result = defaultdict(F)
    for (raw, out), (other_raw, other_out) in zip(first, second, strict=True):
        if raw == other_raw and out == other_out:
            result[_common(law, raw), out] += F(1, 2)
    return dict(result)


def seam_certificate(seam):
    """Check the diagonal-block intertwiner, including every off-diagonal unit."""
    images = {law: columns(law, seam) for law in ('C', 'E')}
    counts = {}
    for law, rows in images.items():
        supports = defaultdict(list)
        for source, row in enumerate(rows):
            for raw, target in row:
                supports[raw].append((source, target))
        for entries in supports.values():
            if len({i for i, _ in entries}) != len(entries) or len({j for _, j in entries}) != len(entries):
                raise ArithmeticError('a recorded native branch lost partial injectivity')
        counts[law] = len(supports)
    for a in range(len(STATES)):
        for b in range(len(STATES)):
            c = _diagonal_image('C', images['C'][a], images['C'][b])
            e = _diagonal_image('E', images['E'][a], images['E'][b])
            expected = {} if a != b else {
                _common('C', raw): target for raw, target in images['C'][a]}
            expected = {(record, target): F(1, 2) for record, target in expected.items()}
            if c != e or c != expected:
                raise ArithmeticError('the common-record diagonal-block identity failed')
    return {'seam': list(seam), 'recorded_partial_injections': counts,
            'matrix_units_checked': len(STATES)**2}


def intertwining_certificate():
    rows = [seam_certificate(seam) for seam in quantum_mixing.seams()]
    return {'native_dimension': len(STATES), 'seams': rows,
            'matrix_units_checked': sum(row['matrix_units_checked'] for row in rows),
            'arbitrary_initial_native_work_reference_correlations': True,
            'arbitrary_work_reference_dimensions': True,
            'work_instruments_act_only_on_work': True,
            'reference_is_untouched': True,
            'native_quantum_output_must_be_traced_or_diagonally_read': True,
            'extra_contrast_feedback_admitted': False,
            'additional_coherent_native_controls_admitted': False,
            'coherent_native_work_coupling_admitted': False}


def native_step(density, law, seam):
    """Recorded native CP maps on sparse exact (native,work,reference) entries."""
    images = columns(law, seam)
    result = {record: {} for record in RECORDS}
    for (a, b), value in density.items():
        for (raw, out), (other_raw, other_out) in zip(images[a[0]], images[b[0]], strict=True):
            if raw == other_raw:
                _add(result[_common(law, raw)],
                     ((out, *a[1:]), (other_out, *b[1:])), _scale(value, F(1, 2)))
    return result


def work_map(density, matrix, weight=F(1)):
    """Apply sqrt(weight)*matrix to W only; R and N are unchanged."""
    result = {}
    for (a, b), value in density.items():
        for i in (0, 1):
            for j in (0, 1):
                coefficient = _multiply(matrix[i][a[1]], _conjugate(matrix[j][b[1]]))
                _add(result, ((a[0], i, a[2]), (b[0], j, b[2])),
                     _scale(_multiply(coefficient, value), weight))
    return result


IDENTITY = ((ONE, ZERO), (ZERO, ONE))
HADAMARD_NUMERATOR = ((ONE, ONE), (ONE, (F(-1), F(0))))
PHASE = ((ONE, ZERO), (ZERO, I))


def work_read(density, axis):
    if axis not in ('X', 'Z'):
        raise ValueError('X or Z work read required')
    if axis == 'Z':
        matrices = (((ONE, ZERO), (ZERO, ZERO)), ((ZERO, ZERO), (ZERO, ONE)))
        return tuple(work_map(density, matrix) for matrix in matrices)
    matrices = tuple(((ONE, (F(sign), F(0))), ((F(sign), F(0)), ONE)) for sign in (1, -1))
    return tuple(work_map(density, matrix, F(1, 4)) for matrix in matrices)


def trace(density):
    value = ZERO
    for (a, b), entry in density.items():
        if a == b:
            value = (value[0]+entry[0], value[1]+entry[1])
    if value[1]:
        raise ArithmeticError('an operator has nonreal trace')
    return value[0]


def reduced_work_reference(density):
    result = {}
    for (a, b), value in density.items():
        if a[0] == b[0]:
            _add(result, (a[1:], b[1:]), value)
    return result


def initial_state():
    """A normalized pure N-W-R state with complex cross-native coherences."""
    amplitudes = {(witness.X, 0, 0): ONE, (witness.X, 1, 1): I,
                  (witness.Z, 0, 1): ONE, (witness.Z, 1, 0): (F(-1), F(0))}
    return {(a, b): _scale(_multiply(x, _conjugate(y)), F(1, 4))
            for a, x in amplitudes.items() for b, y in amplitudes.items()}


def execute(law, initial=None):
    """Adaptive exact witness, retaining 504 record words including zeros.

    Native current controls an H/S work pulse. Its X/Z outcome chooses the
    second seam. The second current controls a further work pulse. A Z work
    read stops on zero or invokes a third repair on one. No outcome is dropped.
    This is an access-class witness, not an oscillator or laboratory model.
    """
    initial = initial_state() if initial is None else initial
    result = {}
    for (c1, j1), first in native_step(initial, law, (0, 1)).items():
        first_controlled = work_map(first, HADAMARD_NUMERATOR, F(1, 2)) if j1 else work_map(first, PHASE)
        for a, measured in enumerate(work_read(first_controlled, 'X' if c1 == 0 else 'Z')):
            seam = (0, 1) if a == 0 else (0, 5)
            for (c2, j2), second in native_step(measured, law, seam).items():
                second_controlled = work_map(second, PHASE) if j2 else work_map(second, HADAMARD_NUMERATOR, F(1, 2))
                for b, final in enumerate(work_read(second_controlled, 'Z')):
                    prefix = (c1, j1, a, c2, j2, b)
                    if b == 0:
                        result[(*prefix, 'halt')] = final
                    else:
                        for (c3, j3), terminal in native_step(final, law, (1, 5)).items():
                            result[(*prefix, c3, j3)] = terminal
    if len(result) != 504 or sum(map(trace, result.values()), F()) != trace(initial):
        raise ArithmeticError('the complete adaptive witness lost a branch or trace')
    return result


def coherent_access_control():
    """The existing coherent native detector lies outside the no-effect class."""
    laws = {}
    for law in ('C', 'E'):
        law_rows = {}
        for setting in witness.SETTINGS:
            joint = {((a, 0, 0), (b, 0, 0)): value
                     for (a, b), value in witness.prepared_density(setting).items()}
            actual = native_step(joint, law, (0, 1))
            plus = F()
            for block in actual.values():
                native_density = {(a[0], b[0]): value for (a, b), value in block.items()}
                plus += witness.trace(witness.readout_density(native_density, 1))
            if plus != witness.prediction(law, setting):
                raise ArithmeticError('the coherent-access control changed the existing instrument')
            law_rows[setting] = str(plus)
        laws[law] = law_rows
    return {'plus_probabilities': laws, 'native_coherent_readout_added': True,
            'existing_preparation_and_detector_unchanged_between_candidates': True,
            'physical_apparatus_identified': False}


def run():
    certificate = intertwining_certificate()
    outputs = {law: execute(law) for law in ('C', 'E')}
    reduced = {law: {record: reduced_work_reference(rho) for record, rho in rows.items()}
               for law, rows in outputs.items()}
    if reduced['C'] != reduced['E']:
        raise ArithmeticError('record-controlled work distinguishes the native rules')
    squared_difference = F()
    for record in outputs['C']:
        c, e = outputs['C'][record], outputs['E'][record]
        for key in c.keys() | e.keys():
            x, y = c.get(key, ZERO), e.get(key, ZERO)
            squared_difference += sum((a-b)**2 for a, b in zip(x, y, strict=True))
    if not squared_difference:
        raise ArithmeticError('the witness no longer retains a hidden native quantum difference')
    rows = [{'record': list(record), 'probability': str(trace(rho)),
             'work_reference_density': [list(a)+list(b)+list(map(str, value))
                 for (a, b), value in sorted(rho.items())]}
            for record, rho in reduced['C'].items()]
    return {'schema': 'oph.rule_observability.v1', 'intertwining': certificate,
            'adaptive_execution': {'complete_record_count': len(rows),
                'nonzero_record_count': sum(bool(rho) for rho in reduced['C'].values()),
                'common_law_and_every_conditional_work_reference_operator_equal': True,
                'complete_probability': str(sum((trace(rho) for rho in reduced['C'].values()), F())),
                'retained_native_quantum_block_difference_HS_squared': str(squared_difference),
                'initial_state_includes_native_work_reference_correlations': True,
                'records': rows},
            'coherent_access_control': coherent_access_control(),
            'scope': {'empirical_likelihoods_constructed': False,
                'oscillator_noise_or_action_inserted': False,
                'native_rule_empirically_selected': False,
                'identity_applies_to_common_record_controlled_work': True,
                'complete_native_quantum_instruments_identical': False},
            'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                              for name in SOURCE_PATHS}}


def compare(expected, actual):
    """Exact JSON comparison, distinguishing Boolean/count/float aliases."""
    def canonical(value):
        return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)
    if canonical(expected) != canonical(actual):
        raise ValueError('the exact observability receipt does not replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--check', type=Path)
    mode.add_argument('--write', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare(strict_json.load(args.check), result)
        print('Exact C/E reduced-observer equivalence and coherent boundary replayed.')
    else:
        payload = json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n'
        if args.write:
            with args.write.open('x', encoding='utf-8') as stream:
                stream.write(payload)
        else:
            print(payload, end='')


if __name__ == '__main__':
    main()
