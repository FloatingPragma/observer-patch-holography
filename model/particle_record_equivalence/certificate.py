"""Source-bound record recovery, recurrent memory and duty certificates.

The finite algebra is exact. Eigenvalue/readout numbers are floating diagnostics.
No empirical reference values enter the producer.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as s

from ..common_source import native_code_transport as code
from ..common_source import native_family_record_memory as memory
from ..common_source import native_gauge_object as geometry
from ..native_flavor_orbit import intervals as iv, spectrum as certified_spectrum

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SUPPORT = (2, 3, 9, 10)


def clean(matrix):
    return s.Matrix(matrix).applyfunc(lambda x: s.simplify(s.expand(x)))


def entries(matrix):
    return [[str(x) for x in row] for row in clean(matrix).tolist()]


def validate_images(images, *, dimension=6):
    """Admit the actual uniform recorded error on all logical columns."""
    if len(images) != dimension:
        raise ValueError('every logical column is required')
    lengths = {len(row) for row in images}
    if len(lengths) != 1:
        raise ValueError('branch probability must be independent of logical input')
    flat = [address for row in images for address in row]
    if len(flat) != len(set(flat)):
        raise ValueError('recorded logical columns must have orthogonal images')
    if any(type(i) is not int or not 0 <= i < 512 for i in flat):
        raise ValueError('native output addresses must lie in the complete carrier')
    count = next(iter(lengths))
    if count > 6:
        raise ValueError('a branch cannot exceed the six source terms')
    return Fraction(count, 12)


@lru_cache(maxsize=1)
def recovery_certificate():
    """Check recorded isometry identities for both codes and all native seams.

    Uniform code terms have amplitude 1/sqrt(6); a fair placement coin
    contributes 1/sqrt(2). Disjoint equal-cardinality branch images imply
    (K_y E)^* K_y E = (count/12) I6 exactly. On a positive branch, V_y is
    the normalized image, and R_y=E V_y^* on that image recovers the code.
    """
    branches = positive = 0
    weights = set()
    for carrier in code.PAIR_CODES.values():
        for law, nraw in (('C', 4), ('E', 8)):
            for seam in geometry.geometry()['seams']:
                total = Fraction(0)
                for coin in (0, 1):
                    for raw in range(nraw):
                        images = carrier.branch_images(law, seam, coin, raw)
                        weight = validate_images(images)
                        assert weight == carrier.branch_weight(law, seam, coin, raw)
                        total += weight
                        weights.add(str(weight))
                        branches += 1
                        positive += bool(weight)
                assert total == 1
    return {
        'codes': 2, 'seams_per_code': 30, 'laws': ['C', 'E'],
        'logical_dimension': 6, 'carrier_dimension': 512,
        'branches_checked_including_zero': branches,
        'positive_branches': positive, 'exact_branch_weights': sorted(weights),
        'exact_identity': 'R_y K_y E = sqrt(w_y) E',
        'record_independence_of_logical_state': True,
        'arbitrary_reference_identity': True,
        'compiled_service_error_set_to_zero': False,
        'scope': 'Ideal recorded recovery on the specified matching code; source compilation and physical service errors retain their separate budgets.'}


def branch_destinations():
    data = memory.architecture()
    return [(tuple(row[:4]), row[-1], tuple(memory._destination(
        old, row[4], row[2] == 0 and row[3] in (4, 5),
        data['reset'], data['gC'], data['gE']) for old in range(12)))
        for row in data['records']]


def compression_certificate(*, destinations=None):
    """Check the actual all-outcome recurrent quantum-instrument embedding.

    A raw/old branch is sqrt(w)|destination><old| tensor U_destination(h).
    Matching its old/destination labels proves the quantum identity for all h.
    All raw labels and every original old-memory outcome remain explicit;
    outcomes outside SUPPORT have zero probability on the admitted code.
    """
    data = memory.architecture()
    E = s.zeros(12, 4)
    for j, old in enumerate(SUPPORT):
        E[old, j] = 1
    source_rows = branch_destinations()
    rows = source_rows if destinations is None else destinations
    if destinations is not None and rows != source_rows:
        raise ValueError('supplied source branch table differs from the native weights or destinations')
    if len(rows) != 128 or tuple(row[0] for row in rows) != memory.RAW_LABELS:
        raise ValueError('all128 source raw labels are required in source order')
    if sum(row[1] for row in rows) != 1 or any(w < 0 for _, w, _ in rows):
        raise ValueError('complete nonnegative branch weights are required')
    small = s.zeros(4)
    identities = zero_slots = 0
    for label, w, dest in rows:
        if len(dest) != 12:
            raise ValueError('all12 source memory inputs are required')
        for old in range(12):
            if old not in SUPPORT or w == 0:
                zero_slots += 1
                continue
            new = dest[old]
            if new not in SUPPORT:
                raise ValueError('a positive branch leaves the recurrent encoding')
            i, j = SUPPORT.index(new), SUPPORT.index(old)
            A = s.zeros(12)
            A[new, old] = w
            B = s.zeros(4)
            B[i, j] = w
            if A*E != E*B:
                raise ArithmeticError('outcome-specific memory intertwining failed')
            small += B
            identities += 1
    if any(sum(small[:, j]) != 1 for j in range(4)):
        raise ArithmeticError('compressed quantum instrument is incomplete')
    pi = E.T*data['stationary']
    if E*pi != data['stationary'] or small*pi != pi:
        raise ArithmeticError('stationary initialization is incompatible')
    if data['T']*E != E*small:
        raise ArithmeticError('unread memory embedding failed')
    assert E.T*E == s.eye(4)
    startup_domain = tuple(old for old in range(12)
                           if all(w == 0 or dest[old] in SUPPORT for _, w, dest in rows))
    zero_start = s.zeros(4, 1)
    for _, w, dest in rows:
        if w:
            zero_start[SUPPORT.index(dest[0])] += w
    assert startup_domain == (0, 1, 2, 3, 9, 10)
    assert zero_start == s.Matrix([s.Rational(71, 72), 0, s.Rational(1, 72), 0])
    return {
        'full_memory_addresses': 12, 'compressed_memory_addresses': 4,
        'embedding_ports': list(SUPPORT),
        'complete_raw_labels': 128, 'complete_raw_old_outcomes': 1536,
        'positive_branch_intertwinings': identities,
        'zero_outcomes_on_admitted_code': zero_slots,
        'transition_exact': entries(small), 'stationary_exact': entries(pi),
        'admitted_initial_states': 'Every state supported in the four-address memory code, including arbitrary work/reference correlations; reset2 and stationary memory are included.',
        'full_quantum_identity': 'K_y E = E B_y for every raw/old outcome, with identical U_destination(h); equality holds for every real h.',
        'complete_history_identity': True,
        'one_step_startup_domain': list(startup_domain),
        'zero_start_memory_after_first_event': entries(zero_start),
        'zero_start_scope': 'Initial address0 is admitted after retaining its full first-event instrument and original old0 record; its continuation then uses the same four-address encoding.',
        'physical_capacity_universality': False,
        'scope': 'An encoding of the same declared source law. Transient initial memory outside the code, altered schedules, clocks or physical capacities are not identified.'}


@lru_cache(maxsize=1)
def duty_certificate():
    data = memory.architecture()
    policies = ('always', 'C', 'E', 'CE')
    Z = {policy: s.zeros(3) for policy in policies}
    for record in data['records']:
        cc, cr, ec, er, c, w = record
        e = ec == 0 and er in (4, 5)
        for old in range(12):
            new = memory._destination(old, c, e, data['reset'], data['gC'], data['gE'])
            term = w*data['stationary'][old]*data['projectors'][new]
            enabled = {'always': True, 'C': c, 'E': e, 'CE': c and e}
            for policy in policies:
                if enabled[policy]:
                    Z[policy] += term
    Z = {policy: clean(value) for policy, value in Z.items()}
    assert clean(Z['always'] - data['Y']) == s.zeros(3)
    assert clean(Z['C'] - (data['Y']-(1-data['p'])*data['projectors'][2])) == s.zeros(3)
    expected_traces = {'always': s.Integer(1), 'C': data['p'], 'E': data['q'], 'CE': data['p']*data['q']}
    result = {}
    for policy, z in Z.items():
        trace = s.trace(z)
        assert trace == expected_traces[policy]
        eigenvalues = np.linalg.eigvalsh(np.asarray(z.evalf(30), float))
        rank = int(z.rank())
        result[policy] = {'matrix_exact': entries(z), 'trace_exact': str(trace),
                          'rank_exact': rank,
                          'eigenvalues_diagnostic': eigenvalues.tolist(),
                          'ratios_diagnostic': (eigenvalues/eigenvalues[0]).tolist() if rank == 3 else None,
                          'normalized_trace_square_exact': str(s.simplify(s.trace(z*z)/trace**2))}
    # Different normalized quadratic traces exclude similarity even after
    # the one common scalar needed to normalize the duty.
    assert result['always']['normalized_trace_square_exact'] != result['C']['normalized_trace_square_exact']
    return {'policies': result,
            'C_identity_exact': 'Z_C=Y-(1-p)P2',
            'C_is_p_times_Y': False,
            'C_and_always_differ_beyond_common_clock_scale': True,
            'conditioning_on_C_divides_Z_C_by_p': True,
            'conditioning_changes_eigenvalue_ratios': False,
            'duty_changes_source_generator': True,
            'physical_mass_map_derived': False}


def certified_family_roots(matrix):
    """Reuse the finite-stack rational LDL inertia certificate for trace-one Y."""
    boxes = certified_spectrum.eigenvalue_boxes(iv.matrix(matrix), rank=3)
    return [(s.Rational(v.lo.numerator, v.lo.denominator),
             s.Rational(v.hi.numerator, v.hi.denominator)) for v in boxes]


def readout_certificate():
    data = memory.architecture()
    # Independent exact finite mixing, including the final observed transition.
    p0 = s.zeros(12, 1); p0[2] = 1
    p9 = data['T']**(memory.BURNIN+1)*p0
    difference = p9-data['stationary']
    l1 = sum(abs(x) for x in difference)
    upper = 2*data['p']**(memory.BURNIN+1)
    assert l1 <= upper
    Yn = clean(sum((p9[k]*data['projectors'][k] for k in range(12)), s.zeros(3)))
    values = np.linalg.eigvalsh(np.asarray(data['Y'].evalf(30), float))
    finite_values = np.linalg.eigvalsh(np.asarray(Yn.evalf(30), float))
    delta = float(upper)
    root_intervals = certified_family_roots(data['Y'])
    assert root_intervals[0][0] > upper
    exact_ratios = [((lo-upper)/(root_intervals[0][1]+upper),
                     (hi+upper)/(root_intervals[0][0]-upper))
                    for lo, hi in root_intervals[1:]]
    intervals = [[float(lo), float(hi)] for lo, hi in exact_ratios]
    records = memory.forward_records()
    budget = memory.readout_error_budget()
    decoded = memory.decode_records(records['probabilities'], steps=records['steps'],
                                    operator_error=budget['operator_error_upper'])
    inverse_error = float(np.linalg.norm(decoded['estimate']-np.asarray(Yn, float), 2))
    if inverse_error > budget['numerical_operator_allowance']:
        raise ArithmeticError('two-clock complete-record inverse failed')
    return {'burnin': memory.BURNIN, 'total_memory_transitions': memory.BURNIN+1,
            'exact_memory_l1_error': str(l1), 'exact_memory_l1_upper': str(upper),
            'operator_bias_upper': delta,
            'stationary_ratios_diagnostic': (values/values[0]).tolist(),
            'finite_ratios_diagnostic': (finite_values/finite_values[0]).tolist(),
            'family_enclosure_method': 'rational interval LDL inertia from native_flavor_orbit.spectrum.eigenvalue_boxes',
            'family_root_enclosures_exact': [[str(lo), str(hi)] for lo, hi in root_intervals],
            'finite_preparation_ratio_enclosures_exact': [[str(lo), str(hi)] for lo, hi in exact_ratios],
            'mixing_only_stationary_ratio_diagnostic_intervals': intervals,
            'ratio_endpoint_floats_are_outward_rounded': False,
            'complete_pointer_records_retained': int(records['probabilities'].size),
            'two_clock_inverse_error_diagnostic': inverse_error,
            'two_clock_numerical_allowance': budget['numerical_operator_allowance'],
            'inverse_has_zero_exact_finite_pulse_bias': True,
            'physical_service_errors_or_shot_noise_certified': False}


def pole_certificate():
    z = s.Symbol('z')
    # An unobserved carrier complement may have extra eigenvalues. Its
    # compressed logical response must not falsely inherit those poles.
    B = s.diag(s.Rational(1, 2), s.Rational(1, 3))
    T = s.diag(s.Rational(1, 2), s.Rational(1, 3), s.Rational(1, 5))
    E = s.Matrix([[1, 0], [0, 1], [0, 0]])
    L = E.T
    assert L*E == s.eye(2) and T*E == E*B
    recovered = clean(L*(s.eye(3)-z*T).inv()*E)
    assert recovered == (s.eye(2)-z*B).inv()
    visible = (s.ones(1, 2)*recovered*s.ones(2, 1))[0]
    masked = (s.Matrix([[1, 0]])*recovered*s.Matrix([1, 0]))[0]
    return {'exact_logical_resolvent': entries(recovered),
            'full_visible_response': str(s.factor(visible)),
            'masked_response': str(s.factor(masked)),
            'masked_poles': ['3'],
            'unreachable_carrier_pole': '5',
            'faithful_recovery_moves_logical_poles': False,
            'lossy_observation_can_remove_residues': True,
            'changed_feedback_can_change_dynamics': True}


def build():
    paths = ['model/common_source/native_code_transport.py',
             'model/common_source/record_code.py',
             'model/common_source/native_repair_mass.py',
             'model/common_source/native_family_record_memory.py',
             'model/common_source/native_gauge_object.py',
             'model/particle_record_equivalence/certificate.py',
             'model/native_flavor_orbit/spectrum.py',
             'model/native_flavor_orbit/intervals.py']
    return {'schema': 'oph.particle_record_equivalence.v1',
            'scope': 'Exact ideal record and encoding identities for a declared finite source; no empirical correction is fitted.',
            'recovery': recovery_certificate(),
            'memory_compression': compression_certificate(),
            'duty': duty_certificate(), 'readout': readout_certificate(),
            'poles': pole_certificate(),
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = build()
    output = HERE/'results.json'
    if args.check:
        old = json.loads(output.read_text())
        # Only explicitly labelled numerical diagnostics use tolerance.
        def compare(a, b, path=''):
            if type(a) is float and type(b) in (int, float):
                if not np.isclose(a, b, atol=2e-12, rtol=1e-11):
                    raise AssertionError(path)
            elif isinstance(a, dict):
                assert set(a) == set(b), path
                for k in a: compare(a[k], b[k], path+'.'+k)
            elif isinstance(a, list):
                assert len(a) == len(b), path
                for i, (x, y) in enumerate(zip(a, b)): compare(x, y, path+f'[{i}]')
            else: assert a == b, path
        compare(result, old)
        print('particle record equivalence: receipt verified')
    else:
        output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
        print(output)


if __name__ == '__main__':
    main()
