"""Exact passive-repair versus modular-flow tests on native observer algebras.

Modular unitaries preserve distance from the tracial state. The admitted
unread native laws contract that distance instead. A separate Palm-event
experiment keeps every waiting/timeout record and exposes the distinction
even inside the modular centralizer. Active coherent controllers are outside
this passive comparison.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as sp

from ..observer import native, quantum_mixing as mixing
from ..record_laws import burn_in
from ..repair_identification import core as exact
from tools import strict_json
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
STATES = tuple(native.states(binary=True))
INDEX = {state: i for i, state in enumerate(STATES)}
D = len(STATES)
SOURCE_PATHS = tuple(sorted(set((*burn_in.SOURCE_PATHS, *exact.SOURCE_PATHS,
    'model/relative_modular_dynamics/__init__.py',
    'model/relative_modular_dynamics/passive.py',
    'model/repair_identification/__init__.py'))))


def _integer(value, low, high):
    if type(value) is not int or not low <= value <= high:
        raise ValueError('an original integer in the declared finite range is required')
    return value


def _fraction(value, low=F(0), high=F(1)):
    value = exact._fraction(value)
    if not low <= value <= high:
        raise ValueError('an exact fraction in the declared range is required')
    return value


def operator(value, dimension=D):
    """Canonical sparse rational-complex operator; no shared Fraction aliases."""
    _integer(dimension, 1, D)
    if type(value) is not dict:
        raise ValueError('a sparse rational-complex operator is required')
    result = {}
    for key, coefficient in value.items():
        if (type(key) is not tuple or len(key) != 2
                or any(type(i) is not int or not 0 <= i < dimension for i in key)
                or type(coefficient) is not tuple or len(coefficient) != 2):
            raise ValueError('matrix-unit indices and exact real/imaginary pairs required')
        parts = tuple(exact._fraction(x) for x in coefficient)
        if any(parts):
            result[key] = parts
    return result


def identity(dimension=D):
    _integer(dimension, 1, D)
    return {(i, i): (F(1), F(0)) for i in range(dimension)}


def _matrix(value, dimension):
    return sp.ImmutableSparseMatrix(dimension, dimension,
        {key: sp.Rational(a)+sp.I*sp.Rational(b) for key, (a, b) in operator(value, dimension).items()})


def _sparse(matrix):
    return {key: tuple(F(x) for x in value.as_real_imag()) for key, value in matrix.todok().items()}


def inspect_algebra(basis, reference, heisenberg_map):
    """Exact algebra and necessary modular-isometry checks, without logarithms.

    The caller supplies a linear map; this function does not certify its
    complete positivity. It checks its values on the entire declared algebra
    by a basis, including possible leakage outside that algebra. Passing the
    norm checks is necessary, not sufficient, for modular dynamics.
    """
    if type(reference) is not tuple or type(basis) is not tuple or not basis or not callable(heisenberg_map):
        raise ValueError('a finite operator basis, faithful diagonal reference and linear map required')
    d = _integer(len(reference), 1, D)
    reference = tuple(_fraction(x) for x in reference)
    if not all(reference) or sum(reference) != 1:
        raise ValueError('the diagonal modular reference must be faithful and normalized')
    matrices = tuple(_matrix(x, d) for x in basis)
    columns = sp.ImmutableSparseMatrix(d*d, len(matrices),
        {(i*d+j, k): x for k, matrix in enumerate(matrices) for (i, j), x in matrix.todok().items()})
    gram = columns.H*columns
    if gram.det() == 0:
        raise ValueError('the declared operator basis is linearly dependent')
    inverse = gram.inv(method='DM')

    def project(matrix):
        vector = sp.ImmutableSparseMatrix(d*d, 1, {(i*d+j, 0): x for (i, j), x in matrix.todok().items()})
        return sp.ImmutableSparseMatrix((columns*(inverse*(columns.H*vector))).reshape(d, d))

    def norm_squared(matrix):
        return F(sum(sp.conjugate(x)*x for x in matrix.todok().values()))

    def member(matrix):
        return norm_squared(matrix-project(matrix)) == 0

    unit = sp.eye(d)
    if not member(unit):
        raise ValueError('the observer algebra must contain the ambient identity')
    if any(not member(x.H) for x in matrices):
        raise ValueError('the observer space is not adjoint closed')
    if any(not member(x*y) for x in matrices for y in matrices):
        raise ValueError('the observer space is not multiplicatively closed')
    if _matrix(heisenberg_map(identity(d)), d) != unit:
        raise ValueError('the supplied Heisenberg map is not unital')
    frequencies, invariant = set(), True
    rows = []
    for number, matrix in enumerate(matrices):
        pieces = defaultdict(dict)
        for (i, j), value in matrix.todok().items():
            ratio = reference[i]/reference[j]
            pieces[ratio][i, j] = value
        for ratio, entries in pieces.items():
            frequencies.add(ratio)
            invariant &= member(sp.ImmutableSparseMatrix(d, d, entries))
        centered = matrix-sp.trace(matrix)*unit/d
        image = _matrix(heisenberg_map(_sparse(centered)), d)
        rows.append({'basis_index': number, 'centered_HS_norm_squared': str(norm_squared(centered)),
            'image_HS_norm_squared': str(norm_squared(image)),
            'compressed_image_HS_norm_squared': str(norm_squared(project(image))),
            'leakage_HS_norm_squared': str(norm_squared(image-project(image))),
            'strict_HS_loss': norm_squared(image) < norm_squared(centered)})
    loss = any(row['strict_HS_loss'] for row in rows)
    return {'dimension': len(basis), 'ambient_dimension': d,
        'unital_star_algebra_checked': True, 'modular_invariant_for_every_parameter': invariant,
        'modular_eigenvalue_ratios': list(map(str, sorted(frequencies))),
        'native_map_preserves_algebra': all(F(row['leakage_HS_norm_squared']) == 0 for row in rows),
        'strict_loss_obstructs_every_unitary_modular_parameter': loss,
        'passing_isometry_check_is_a_sufficiency_claim': False, 'CPTP_admissibility_checked_here': False,
        'basis_checks': rows}


@lru_cache(maxsize=2)
def _branches(law):
    """Independent actual primitive branches, cached as integer tuples only."""
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('native record law C or E required')
    groups = defaultdict(dict)
    for seam_number, seam in enumerate(mixing.seams()):
        for coin in (0, 1):
            for i, state in enumerate(STATES):
                after, (record_coin, current) = native.native(state, seam, coin)
                if record_coin != coin or current != state[seam[0]]-after[seam[0]]:
                    raise ArithmeticError('native current convention changed')
                label = current if law == 'C' else state[seam[0]]-state[seam[1]]
                groups[seam_number, coin, label, abs(current)][i] = INDEX[after]
    if any(len(set(mapping.values())) != len(mapping) for mapping in groups.values()):
        raise ArithmeticError('a committed branch is not a partial injection')
    return tuple((key[-1], tuple(sorted(mapping.items()))) for key, mapping in sorted(groups.items()))


def native_map(law, value, *, adjoint=False, current=None):
    if type(adjoint) is not bool or (current is not None and (type(current) is not int or current not in (0, 1))):
        raise ValueError('an exact adjoint switch and optional absolute-current bit are required')
    value = operator(value)
    result = defaultdict(lambda: [F(0), F(0)])
    for bit, entries in _branches(law):
        if current is not None and bit != current:
            continue
        mapping = {b: a for a, b in entries} if adjoint else dict(entries)
        for (i, j), parts in value.items():
            if i in mapping and j in mapping:
                target = result[mapping[i], mapping[j]]
                for k in (0, 1):
                    target[k] += parts[k]/60
    return operator({key: tuple(value) for key, value in result.items()})


def orbit_gap(initial_tracial_distance, attempts, coverage=F(1)):
    """Lower bound against every source modular unitary, at any chosen parameter.

    With a reference, the supplied initial distance must be to tau tensor the
    untouched reference marginal. Record-conditioned phases or retaining the
    warm-up tape are not covered by this single-orbit comparison.
    """
    distance = _fraction(initial_tracial_distance)
    _integer(attempts, 0, 10000)
    coverage = _fraction(coverage)
    if coverage == 0:
        raise ValueError('strictly positive seam coverage is required')
    return max(F(0), distance-burn_in.binary_bound(attempts, coverage))


def observer_HS_gap(attempts, coverage=F(1)):
    """Relative HS gap for a centered observer, including algebra compression."""
    _integer(attempts, 0, 10000)
    coverage = _fraction(coverage)
    if coverage == 0:
        raise ValueError('strictly positive seam coverage is required')
    return 1-(1-coverage*(1-burn_in.R_UPPER))**attempts


def palm_data():
    moving = native_map('C', identity(), current=1)
    degree = tuple(int(60*moving[i, i][0]) for i in range(D))
    if Counter(degree) != {8: 30, 10: 36} or sum(degree) != 600:
        raise ArithmeticError('the native movement degree census changed')
    return degree, tuple(F(x, 600) for x in degree)


def event_history(cap, law='C', *, modular_control=False):
    """Complete herald/degree/wait/degree records, with matched timeout control.

    The control keeps each initial degree's native waiting law but preserves
    that degree through any modular phase. It is an explicit comparison law,
    not a claim that native waiting implements that phase. Every record is
    unconditional; unsuccessful heralds and timeout are never normalized away.
    """
    _integer(cap, 0, 32)
    if type(modular_control) is not bool:
        raise ValueError('a Boolean modular comparison switch is required')
    degree, _ = palm_data()
    initial = {(i, i): (F(1, D), F(0)) for i in range(D)}
    herald = native_map(law, initial, current=1)
    mass = sum(parts[0] for (i, j), parts in herald.items() if i == j)
    result = {('failed_herald',): 1-mass}
    for first in (8, 10):
        state = {key: value for key, value in herald.items() if degree[key[0]] == first}
        for number in range(1, cap+1):
            event = native_map(law, state, current=1)
            total = sum(parts[0] for (i, j), parts in event.items() if i == j)
            for last in (8, 10):
                weight = (total if first == last else F(0)) if modular_control else sum(
                    parts[0] for (i, j), parts in event.items() if i == j and degree[i] == last)
                result['complete', first, number, last] = weight
            state = native_map(law, state, current=0)
        result['timeout', first, cap] = sum(parts[0] for (i, j), parts in state.items() if i == j)
    if sum(result.values()) != 1:
        raise ArithmeticError('the complete degree history lost probability')
    return result


def run():
    degree, palm = palm_data()
    unit = identity()
    projector = {(i, i): (F(1), F(0)) for i, value in enumerate(degree) if value == 8}
    i, j = next((i, j) for i in range(D) for j in range(D) if degree[i] == 8 and degree[j] == 10)
    block = (unit, {(i, i): (F(1), F(0))}, {(j, j): (F(1), F(0))},
             {(i, j): (F(1), F(0)), (j, i): (F(1), F(0))},
             {(i, j): (F(0), F(1)), (j, i): (F(0), F(-1))})
    algebras, histories = {}, []
    for law in ('C', 'E'):
        mapping = lambda value, law=law: native_map(law, value, adjoint=True)
        algebras[law] = {'scalar': inspect_algebra((unit,), palm, mapping),
            'degree_centralizer': inspect_algebra((unit, projector), palm, mapping),
            'coherent_block': inspect_algebra(block, palm, mapping)}
        if (algebras[law]['scalar']['strict_loss_obstructs_every_unitary_modular_parameter']
                or not all(algebras[law][name]['strict_loss_obstructs_every_unitary_modular_parameter']
                           for name in ('degree_centralizer', 'coherent_block'))):
            raise ArithmeticError('the source algebra obstruction controls failed')
        for cap in (0, 1, 2, 8):
            observed = event_history(cap, law)
            target = event_history(cap, law, modular_control=True)
            distance = sum(abs(observed[key]-target[key]) for key in observed)/2
            expected = (2-F(13, 15)**cap-F(5, 6)**cap)/33
            if distance != expected or any(observed[key] != target[key] for key in observed if key[0] != 'complete'):
                raise ArithmeticError('the native degree-change law or matched timeouts changed')
            histories.append({'law': law, 'cap': cap, 'records_including_zero_probability': len(observed),
                'unconditional_record_TV_from_every_matched_modular_phase': str(distance),
                'herald_fail_and_timeout_probabilities_match': True,
                'native_history': [{'record': list(key), 'probability': str(value)} for key, value in sorted(observed.items())]})
    tolerance = F(1, 1000)
    attempts = burn_in.sufficient_attempts(tolerance)
    if attempts != 221 or burn_in.binary_bound(attempts) > tolerance:
        raise ArithmeticError('the upstream law-uniform warm-up budget changed')
    radius = sum(abs(p-F(1, D)) for p in palm)/2
    if radius != F(3, 55):
        raise ArithmeticError('Palm and attempt references were conflated')
    return {'schema': 'oph.passive_modular_observer.v1',
        'source_sha256': bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}),
        'algebra_audits': algebras, 'coherent_block_indices': [i, j],
        'all_law_orbit_bound': {'attempts': attempts, 'trace_error_upper': str(tolerance),
            'pure_input_minimum_distance_from_any_modular_orbit': str(F(65, 66)-tolerance),
            'Palm_reference_minimum_distance_from_its_stationary_modular_orbit': str(radius-tolerance),
            'relative_observer_HS_gap_after_one_attempt_lower': str(observer_HS_gap(1)),
            'arbitrary_untouched_reference_version': 'replace initial distance by dist(rho_SR,tau tensor rho_R)',
            'parameter_choice': 'any deterministic choice at each attempt count; no record-conditioned orbit mixture'},
        'Palm_event_histories': histories,
        'completed_event_unconditional_TV_limit': '2/33',
        'completed_event_conditional_TV_limit_given_herald': '2/5',
        'scope': {'actual_native_C_E_maps_independently_enumerated': True,
            'passive_attempt_obstruction_all_admitted_current_tables': True,
            'full_reference_tracial_modular_flow_is_trivial': True,
            'Palm_stationarity_implies_modular_dynamics': False,
            'scalar_algebra_passes_but_contains_no_nontrivial_clock': True,
            'fixed_stationary_state_alone_identifies_dynamics': False,
            'active_coherent_controllers_or_feedback_gates_excluded': True,
            'attempt_bound_substituted_for_stopping_time_bound': False,
            'physical_modular_clock_identified': False}, 'all_pass': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        if strict_json.load(args.check) != result:
            raise ValueError('the passive modular compatibility receipt does not replay')
        print('Passive native modular obstruction and complete event records replayed.')
    else:
        print(json.dumps(result, sort_keys=True, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
