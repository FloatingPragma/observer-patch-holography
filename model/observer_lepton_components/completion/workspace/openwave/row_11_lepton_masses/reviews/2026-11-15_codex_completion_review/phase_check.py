#!/usr/bin/env python3
"""Independent Fraction-only full-kernel pole check; imports no producer code."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib
import json
import sys
import time

HERE = Path(__file__).resolve().parent
PACKET = HERE / 'paused_originals/20261115_codex_intrinsic_phase_process'
STATES = tuple(product(range(4), range(8)))
PI = (F(5, 12), F(5, 12), F(1, 12), F(1, 12))
RHO = tuple(PI[c] / 8 for c, b in STATES)


def kernel(e, d, p):
    """Sum the two independently refreshed code outcomes and original flips."""
    out = []
    for c, b in STATES:
        row = []
        for target, z in STATES:
            total = F(0)
            for mid in range(4):
                before = (1-d) * (c == mid) + d * PI[mid]
                after = (1-d) * (mid == target) + d * PI[target]
                # Code order00,01,10,11 gives J=0,0,+1,-1.
                probabilities = (e/2, e*p/2 if mid != 2 else 0,
                                 e*p/2 if mid == 3 else 0)
                branch = before * after
                for bit, probability in enumerate(probabilities):
                    branch *= probability if (b ^ z) & (1 << bit) else 1-probability
                total += branch
            row.append(total)
        out.append(row)
    return out


def action(matrix, vector):
    return [sum((a*b for a, b in zip(row, vector)), F(0)) for row in matrix]


def inner(a, b):
    return sum((w*x*y for w, x, y in zip(RHO, a, b)), F(0))


def check_case(e, d, p, expected=None):
    matrix = kernel(e, d, p)
    if expected is not None:
        assert matrix == [[F(x) for x in row] for row in expected]
    assert all(sum(row) == 1 for row in matrix)
    assert all(sum(RHO[i]*matrix[i][j] for i in range(32)) == RHO[j] for j in range(32))
    assert all(RHO[i]*matrix[i][j] == RHO[j]*matrix[j][i] for i in range(32) for j in range(32))
    f = [F((c == 0) - (c == 1)) for c, b in STATES]
    f1 = [v * (-1 if b & 2 else 1) for v, (c, b) in zip(f, STATES)]
    ref = [F(-1 if b & 1 else 1) for c, b in STATES]
    poles = ((f, (1-d)**2, F(5, 6)),
             (f1, (1-d)**2*(1-e*p), F(5, 6)),
             (ref, 1-e, F(1)))
    correlations = []
    for observable, eigenvalue, norm in poles:
        assert inner(observable, observable) == norm
        assert action(matrix, observable) == [eigenvalue*x for x in observable]
        vector = observable
        values = []
        for lag in range(9):
            value = inner(observable, vector)
            assert value == norm*eigenvalue**lag
            values.append(str(value))
            vector = action(matrix, vector)
        correlations.append(values)
    # This rejects the earlier tempting standalone phase-gap interpretation.
    assert action(matrix, f1) != [(1-e*p)*x for x in f1]
    relative_ratio = poles[1][1] / poles[0][1]
    assert relative_ratio == 1-e*p
    if p == 0:
        assert poles[0][1] == poles[1][1]
    if p == 1:
        assert relative_ratio == poles[2][1]
    return {'epsilon': str(e), 'delta': str(d), 'p': str(p),
            'lambda_refresh': str(poles[0][1]), 'lambda_phase': str(poles[1][1]),
            'lambda_reference': str(poles[2][1]),
            'relative_pole_ratio': str(relative_ratio),
            'norms': ['5/6', '5/6', '1'], 'correlations_lag0_to8': correlations,
            'standalone_phase_eigenvalue_rejected': True,
            'full_receipt_matrix_matched': expected is not None}


def main():
    assert not sys.flags.optimize, 'assertions must be enabled'
    start = time.perf_counter()
    receipt = json.loads((PACKET/'RECEIPT.json').read_text())
    assert receipt['original_count_probability'] == '2/9'
    assert list(map(F, receipt['refresh_law'])) == list(PI)
    # Derive count and refresh weights from complete stored native records.
    for slot in ('1', '2'):
        table = receipt['count_pair_tables'][slot]
        assert len(table) == 3600 and len({tuple(row[:2]) for row in table}) == 3600
        assert all(row[4] == row[2] ^ row[3] and row[5] == row[2]*row[3] for row in table)
        assert F(sum(row[4] for row in table), 3600) == F(2, 9)
        assert F(sum(row[5] for row in table), 3600) == F(1, 72)
    for role in ('refresh_before', 'refresh_after'):
        records = receipt['source_roles'][role]
        assert len(records) == 120
        codes = ((0, 0), (0, 1), (1, 0), (1, 1))
        assert tuple(F(sum(tuple(r['record_code']) == c for r in records), 120) for c in codes) == PI
    cases = {}
    for k in (1, 2, 3, 6):
        old = receipt['refinements'][str(k)]
        e = F(1, 2**k)
        d = e**3
        cases[f'k{k}_native'] = check_case(e, d, F(2, 9), old['kernel32'])
        for name, p in (('count_HOLD', F(0)), ('count_AND', F(1, 72)), ('count_one_plant', F(1))):
            cases[f'k{k}_{name}'] = check_case(e, d, p, old['controls'][name]['full_kernel'])
    for e, d, p in ((F(1,3), F(1,7), F(2,9)), (F(2,5), F(1,11), F(1,72)),
                    (F(3,7), F(2,13), F(0)), (F(4,9), F(3,17), F(1))):
        cases[f'algebraic_{e}_{d}_{p}'] = check_case(e, d, p)
    result = {'schema': 'row11.phase_pole_completion.v1', 'all_pass': True,
              'cases': cases, 'full_matrix_receipt_cases': 16,
              'additional_algebraic_cases': 4, 'observable_lag_checks': 20*3*9,
              'producer_code_imported': False, 'new_native_calls': 0,
              'old_receipt_sha256': hashlib.sha256((PACKET/'RECEIPT.json').read_bytes()).hexdigest(),
              'scope': 'Exact relative-pole lemma and recorded stationary correlations, not coherent phase execution or physical mass.',
              'seconds': time.perf_counter()-start}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
