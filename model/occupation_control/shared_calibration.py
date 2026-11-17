"""Exact calibration-only obstruction for all fixed transfer intervals.

For every deposited central value y_i and certified apparatus interval [l_i,u_i],
allow one common gain g >= 0, one common offset b, and remaining allowance e >= 0:

    g*l_i + b - e <= y_i <= g*u_i + b + e.

Minimizing e gives the most permissive calibration-only explanation within this
affine family. The optimizing g,b are diagnostic witnesses, never adopted
calibrations. Each endpoint is allowed to vary independently inside its own
interval, an additional relaxation of the actual common apparatus dynamics.

The producer uses a floating LP only to propose active constraints. Exact
rational primal feasibility and a nonnegative dual witness prove the answer.
The verifier checks all constraints and strong dual equality without calling
an optimizer, including the g=0 and e=0 boundaries. This is retrospective
arithmetic on published central values, not a confidence limit or OPH test.

The same certificate also proves the optimum when each potential has its own
gain and offset: its dual is supported within one case, while its common
primal is feasible for every case. The endpoint-only dual additionally gives
a necessary weighted inequality for heterogeneous remaining allowances.
"""
from __future__ import annotations

from copy import deepcopy
import argparse
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
from pathlib import Path

from tools import strict_json
from . import transfer
try:
    from tools.pin_scope import repo_sha256
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import repo_sha256

ROOT = Path(__file__).resolve().parents[2]
SOURCE = 'model/occupation_control/transfer_results.json.gz'
SOURCE_SHA256 = repo_sha256('model/occupation_control/transfer_results.json.gz')
SCHEMA = 'oph.occupation_control.shared_calibration.v2'
CASES = tuple('abcde')
CONTRACT = {
    'model': 'observed in gain*[certified_lower,certified_upper]+offset+[-allowance,allowance]',
    'calibration_families': ['one common gain and offset', 'independent gain and offset per case'],
    'uniform_allowance': 'one remaining nonnegative allowance across all 125 endpoints',
    'heterogeneous_allowances': 'nonnegative endpoint allowances; the retained weighted sum is a necessary lower bound',
    'gain_domain': 'nonnegative; includes zero as a conservative boundary',
    'offset_domain': 'all real numbers',
    'units': 'deposited dimensionless two-point position-estimator units',
    'classification': 'retrospective calibration-only diagnostic',
    'optimized_calibration_adopted': False,
    'predictions_changed': False,
    'experimental_uncertainty_applied': None,
    'statistical_rejection_claimed': False,
    'full_uncertainty_agreement_established': False,
    'OPH_rejection_claimed': False,
    'boundary': ('Only affine output gains and offsets, common or independent per case, are varied. Probe-coordinate, '
                 'preparation, timing, Hamiltonian and physical-service changes are outside '
                 'this family. Missing experimental and calibration uncertainties remain missing.'),
}


def _digest(value):
    try:
        raw = json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
    except (TypeError, ValueError) as error:
        raise ValueError('Finite JSON required') from error
    return hashlib.sha256(raw).hexdigest()


def _keys(value, keys, label):
    if type(value) is not dict or set(value) != set(keys):
        raise ValueError(f'Exact {label} fields required')


def _rational(value):
    if type(value) is not str:
        raise ValueError('Canonical rational string required')
    try:
        result = F(value)
    except (ValueError, ZeroDivisionError, OverflowError) as error:
        raise ValueError('Canonical rational string required') from error
    if str(result) != value:
        raise ValueError('Canonical rational string required')
    return result


def source_rows(*, root=ROOT):
    """Authenticate the retained own-path enclosures; never rerun their ODEs."""
    path = Path(root)/SOURCE
    if hashlib.sha256(path.read_bytes()).hexdigest() != SOURCE_SHA256:
        raise ValueError('Pinned transfer certificate changed')
    retained = transfer.read_receipt(path)
    rows = retained['comparison']['runs'][0]['held_out_rows']
    order = [(case, str(F(step, 5))) for case in 'abcde' for step in range(2, 51, 2)]
    if [(row['case'], row['time_ms']) for row in rows] != order:
        raise ValueError('All 125 endpoints in their fixed order are required')
    if any(row['status'] != 'certified_interval' for row in rows):
        raise ValueError('Every endpoint requires its own certified interval')
    return [{'id': f"{row['case']}@{row['time_ms']}",
             'lower': row['certified_interval'][0], 'upper': row['certified_interval'][1],
             'observed': str(F(row['observed_decimal']))} for row in rows]


def constraints(rows):
    """Produce A*x <= rhs for x=(gain,offset,allowance), preserving row order."""
    if type(rows) is not list or not rows:
        raise ValueError('Nonempty ordered endpoint list required')
    result, seen = [], set()
    for row in rows:
        _keys(row, ('id', 'lower', 'upper', 'observed'), 'endpoint')
        name = row['id']
        if type(name) is not str or not name or name in seen or ':' in name:
            raise ValueError('Distinct endpoint identifiers required')
        seen.add(name)
        lo, hi, value = map(_rational, (row['lower'], row['upper'], row['observed']))
        if lo > hi:
            raise ValueError('Reversed certified interval')
        result.extend([(name+':lo', (lo, F(1), F(-1)), value),
                       (name+':hi', (-hi, F(-1), F(-1)), -value)])
    result.extend([('gain_nonnegative', (F(-1), F(0), F(0)), F(0)),
                   ('allowance_nonnegative', (F(0), F(0), F(-1)), F(0))])
    return result


def _binding(rows):
    return {'endpoint_count': len(rows), 'endpoints_sha256': _digest(rows),
            'constraint_order': [name for name, _, _ in constraints(rows)]}


def verify_exact(rows, certificate):
    """Independently prove a supplied primal/dual optimum with exact arithmetic.

    A*x<=rhs, lambda>=0, sum(lambda*A)=-c imply
    c*x>=-sum(lambda*rhs). Equality certifies optimality. Boundary constraints
    are explicit members of A, so no strictly-positive-gain assumption enters.
    """
    _keys(certificate, ('input_binding', 'primal', 'dual'), 'certificate')
    if _digest(certificate['input_binding']) != _digest(_binding(rows)):
        raise ValueError('Complete ordered constraint binding differs')
    primal = certificate['primal']
    _keys(primal, ('gain', 'offset', 'allowance'), 'primal')
    point = tuple(_rational(primal[key]) for key in ('gain', 'offset', 'allowance'))
    all_constraints = constraints(rows)
    if any(sum(a*x for a, x in zip(coeff, point, strict=True)) > rhs
           for _, coeff, rhs in all_constraints):
        raise ValueError('Primal violates a complete endpoint or boundary constraint')
    dual = certificate['dual']
    if type(dual) is not list or not dual:
        raise ValueError('Nonempty sparse dual witness required')
    positions = {name: index for index, (name, _, _) in enumerate(all_constraints)}
    previous, vector, bound = -1, [F(0), F(0), F(0)], F(0)
    for term in dual:
        _keys(term, ('constraint', 'weight'), 'dual term')
        name = term['constraint']
        if type(name) is not str or name not in positions or positions[name] <= previous:
            raise ValueError('Dual constraints must be unique and in source order')
        previous = positions[name]
        weight = _rational(term['weight'])
        if weight <= 0:
            raise ValueError('Sparse dual weights must be strictly positive')
        _, coefficients, rhs = all_constraints[previous]
        vector = [v+weight*a for v, a in zip(vector, coefficients, strict=True)]
        bound -= weight*rhs
    if vector != [F(0), F(0), F(-1)] or bound != point[2]:
        raise ValueError('Dual cancellation or exact objective equality failed')
    return point[2]


def _case_groups(rows):
    """Bind the stronger claim to every deposited case and exact time label."""
    constraints(rows)
    expected = [f'{case}@{F(step, 5)}' for case in CASES for step in range(2, 51, 2)]
    if [row['id'] for row in rows] != expected:
        raise ValueError('All five fixed cases and their 25 exact endpoint times are required')
    return {case: rows[25*index:25*(index+1)] for index, case in enumerate(CASES)}


def _case_bindings(groups):
    return [{'case': case, 'endpoint_count': len(rows), 'endpoints_sha256': _digest(rows)}
            for case, rows in groups.items()]


def _heterogeneous_terms(rows, dual):
    """Aggregate endpoint coefficients only after exact dual verification.

    Each endpoint inequality has -e_i in place of -e. Cancellation of gain
    and offset therefore gives sum_i weight_i*e_i >= the same dual bound.
    Boundary multipliers cannot be silently converted into endpoint errors.
    """
    weights = {}
    endpoints = {row['id'] for row in rows}
    for term in dual:
        name, separator, side = term['constraint'].rpartition(':')
        if not separator or side not in ('lo', 'hi') or name not in endpoints:
            raise ValueError('Heterogeneous allowance witness must use only selected-case endpoints')
        weights[name] = weights.get(name, F())+_rational(term['weight'])
    if sum(weights.values(), F()) != 1:
        raise ValueError('Endpoint allowance weights must sum exactly to one')
    return [{'endpoint': row['id'], 'weight': str(weights[row['id']])}
            for row in rows if row['id'] in weights]


def independent_case_certificate(rows, common, *, witness_case='a'):
    """Reuse the existing exact optimum; do not solve a second optimization."""
    groups = _case_groups(rows)
    optimum = verify_exact(rows, common)
    if type(witness_case) is not str or witness_case not in groups:
        raise ValueError('A deposited witness case is required')
    selected = groups[witness_case]
    restricted = dict(common, input_binding=_binding(selected))
    if verify_exact(selected, restricted) != optimum:
        raise ValueError('Selected case does not certify the common optimum')
    return {
        'case_bindings': _case_bindings(groups),
        'witness_case': witness_case,
        'uniform_allowance': str(optimum),
        'heterogeneous_allowance': {
            'terms': _heterogeneous_terms(selected, common['dual']),
            'lower_bound': str(optimum),
        },
    }


def verify_independent_cases(rows, common, certificate):
    """Prove the enlarged-family optimum and endpoint-error inequality.

    The common primal satisfies every case, so assigning that same feasible
    pair to all five independent pairs is an upper witness. Restricting the
    exact dual to one complete case gives the matching lower witness even
    when every other gain/offset is unrestricted. No per-case minimizer is
    inferred or fitted, and no experimental error values are assumed.
    """
    _keys(certificate, ('case_bindings', 'witness_case', 'uniform_allowance',
                        'heterogeneous_allowance'), 'independent-case certificate')
    groups = _case_groups(rows)
    if _digest(certificate['case_bindings']) != _digest(_case_bindings(groups)):
        raise ValueError('Complete case binding differs')
    case = certificate['witness_case']
    if type(case) is not str or case not in groups:
        raise ValueError('A deposited witness case is required')
    optimum = verify_exact(rows, common)
    selected = groups[case]
    restricted = dict(common, input_binding=_binding(selected))
    if (verify_exact(selected, restricted) != optimum
            or _rational(certificate['uniform_allowance']) != optimum):
        raise ValueError('Independent-case objective equality failed')
    allowance = certificate['heterogeneous_allowance']
    _keys(allowance, ('terms', 'lower_bound'), 'heterogeneous allowance')
    if (_digest(allowance['terms']) != _digest(_heterogeneous_terms(selected, common['dual']))
            or _rational(allowance['lower_bound']) != optimum):
        raise ValueError('Heterogeneous allowance inequality differs from exact dual')
    return optimum


def _solve_square(matrix, rhs):
    """Small exact Gaussian elimination, used only by the producer."""
    count = len(rhs)
    augmented = [list(row)+[value] for row, value in zip(matrix, rhs, strict=True)]
    for column in range(count):
        pivot = next((i for i in range(column, count) if augmented[i][column]), None)
        if pivot is None:
            raise ValueError('Singular proposed active set')
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value/divisor for value in augmented[column]]
        for row in range(count):
            if row != column:
                factor = augmented[row][column]
                augmented[row] = [a-factor*b for a, b in
                                  zip(augmented[row], augmented[column], strict=True)]
    return tuple(row[-1] for row in augmented)


def solve_exact(rows):
    """Propose a vertex numerically, then retain only an exact verified optimum."""
    from scipy.optimize import linprog
    conditions = constraints(rows)
    proposal = linprog([0., 0., 1.],
                       A_ub=[[float(v) for v in coeff] for _, coeff, _ in conditions],
                       b_ub=[float(rhs) for _, _, rhs in conditions],
                       bounds=[(None, None)]*3, method='highs')
    if not proposal.success:
        raise ValueError('No LP proposal found; no exact result claimed')
    active = [i for i, slack in enumerate(proposal.ineqlin.residual) if abs(slack) < 1e-7]
    for indices in combinations(active, 3):
        chosen = [conditions[i] for i in indices]
        matrix = [item[1] for item in chosen]
        try:
            point = _solve_square(matrix, [item[2] for item in chosen])
            weights = _solve_square(list(zip(*matrix, strict=True)), [F(0), F(0), F(-1)])
            if any(weight < 0 for weight in weights):
                continue
            certificate = {
                'input_binding': _binding(rows),
                'primal': dict(zip(('gain', 'offset', 'allowance'), map(str, point), strict=True)),
                'dual': [{'constraint': item[0], 'weight': str(weight)}
                         for item, weight in zip(chosen, weights, strict=True) if weight],
            }
            verify_exact(rows, certificate)
            return certificate
        except ValueError:
            continue
    raise ValueError('Proposed active set has no exact certificate; no optimum claimed')


def produce():
    rows = source_rows()
    common = solve_exact(rows)
    return {'schema': SCHEMA, 'source': {'path': SOURCE, 'sha256': SOURCE_SHA256},
            'contract': deepcopy(CONTRACT), 'certificate': common,
            'independent_cases': independent_case_certificate(rows, common)}


def verify(receipt):
    _keys(receipt, ('schema', 'source', 'contract', 'certificate', 'independent_cases'), 'receipt')
    if (receipt['schema'] != SCHEMA
            or _digest(receipt['source']) != _digest({'path': SOURCE, 'sha256': SOURCE_SHA256})
            or _digest(receipt['contract']) != _digest(CONTRACT)):
        raise ValueError('Source or diagnostic claim contract changed')
    return verify_independent_cases(source_rows(), receipt['certificate'], receipt['independent_cases'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--output', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.output:
        result = produce()
        verify(result)
        with args.output.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, sort_keys=True, indent=2, allow_nan=False)
            stream.write('\n')
    else:
        optimum = verify(strict_json.load(args.check))
        print(f'Exact all-125 common or independent-case affine allowance: {float(optimum):.12f}')
        print('Calibration-only diagnostic; no fitted calibration adopted or statistical rejection.')


if __name__ == '__main__':
    main()
