"""A local completion invisible on the protected seven-momentum field.

The comparison fixes the internal action, mass kernel, modular clock and
low-momentum metric speed. A two-hop lattice symbol still agrees exactly
on every retained momentum. Equality concerns the restricted observable
algebra: the larger lattice has additional modes, states and resource costs.
An additional momentum pair and neutral pair-quadrature read separate the
two declared kernels, with every selector and clock outcome retained.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import sympy as s
from scipy.linalg import expm

from tools import strict_json
from . import action, gaussian_source as g, kinetic, protected_momentum as protected
from ..halfline import pi_interval
from ...calibrated_measurement.bounds import log_interval, root_interval
from ...characteristic_constraints.bounds import trig_interval
from ...scalar_observations.experimental_readout.stats import concentration_certificate


P0 = s.Rational(1, 3)
SECOND_SHELL = s.Rational(2, 3)
PAIR_SITES = (1, 6)
INTERVAL = F(1, 4)
ACCEPTED = F(15, 16)
ROOT = Path(__file__).resolve().parents[3]


def completion_coefficients():
    c = action.source_coefficients()['Omega']/2
    total = 2*c/(3*s.sqrt(3))
    return c, 2*total-c/s.pi, c/s.pi-total


def lattice_component(momentum):
    """Finite-range-two symbol on a six-cell torus of coordinate side6*pi."""
    if (not isinstance(momentum, s.Basic) or momentum.is_real is not True
            or momentum.atoms(s.Float)):
        raise ValueError('an exact real symbolic coordinate momentum required')
    _, first, second = completion_coefficients()
    return first*s.sin(s.pi*momentum)+second*s.sin(2*s.pi*momentum)


def effective_momentum(momentum, law):
    if law not in ('continuum', 'two_hop_lattice') or type(law) is not str:
        raise ValueError('one of the two fixed kinetic laws required')
    if type(momentum) is not tuple or len(momentum) != 3:
        raise ValueError('three exact coordinate momenta required')
    c = completion_coefficients()[0]
    values = tuple(s.sympify(x) for x in momentum)
    if any(x.is_real is not True or x.atoms(s.Float) for x in values):
        raise ValueError('exact real momenta required before source conversion')
    return tuple(s.simplify(lattice_component(x)/c) if law == 'two_hop_lattice' else x for x in values)


def exact_certificate():
    c, first, second = completion_coefficients()
    p = s.Symbol('p', real=True)
    derivative = s.simplify(s.diff(lattice_component(p), p).subs(p, 0))
    retained = [(s.Integer(0),)*3]
    for axis in range(3):
        for sign in (-1, 1):
            vector = [s.Integer(0)]*3
            vector[axis] = sign*P0
            retained.append(tuple(vector))
    if derivative != c or any(effective_momentum(p, 'two_hop_lattice') != p for p in retained):
        raise ArithmeticError('local completion failed the principal speed or retained field')
    # B<0 follows from pi>3 and sqrt3<2; A+2B=c/pi>0. Therefore
    # sin(pi*p)*(A+2B*cos(pi*p)) has no additional scalar zeros.
    if s.simplify(first+2*second-c/s.pi) != 0:
        raise ArithmeticError('local completion lost its positive sine factor')
    negative = expression_interval(s.simplify(second/c))
    minimum = expression_interval(1/s.pi)
    corner = expression_interval(3-8*s.pi/(3*s.sqrt(3)))
    if negative[1] >= 0 or minimum[0] <= 0 or corner[1] >= 0:
        raise ArithmeticError('exact local-completion sign proof failed')
    second_value = s.simplify(lattice_component(SECOND_SHELL)/c)
    if s.simplify(second_value-(1-s.sqrt(3)/s.pi)) != 0:
        raise ArithmeticError('second shell did not give the separating symbol')
    nodes = [{'coordinate_momentum': list(bits),
              'kinetic_Weyl_chirality': -1 if sum(bits) % 2 == 0 else 1}
             for bits in __import__('itertools').product((0, 1), repeat=3)]
    return {'local_symbol': 'A*sin(pi*p)+B*sin(2*pi*p)',
        'coefficients': {'c': str(c), 'A': str(first), 'B': str(second)},
        'periodic_coordinate_side': '6*pi', 'cells_per_axis': 6, 'coordinate_spacing': 'pi',
        'same_low_momentum_derivative': str(derivative),
        'retained_momenta': [[str(x) for x in p] for p in retained],
        'exact_retained_kinetic_and_Nambu_matrices_equal': True,
        'same_internal_pairing_current_and_clock_operators': True,
        'same_complete_restricted_quantum_histories_for_all_references': True,
        'restriction_contract': 'trace exterior modes and their records; compare the same retained/reference input marginal; retained controls and stopping use no exterior records or failures',
        'retained_continuation_or_stop_uses_exterior_records_or_failures': False,
        'whole_lattice_quantum_outputs_or_source_resources_equal': False,
        'retained_CAR_modes': 630, 'local_lattice_CAR_modes': 6**3*90,
        'unretained_local_lattice_CAR_modes': 6**3*90-630,
        'kinetic_Weyl_nodes': nodes,
        'left_and_right_kinetic_nodes': [4, 4],
        'node_proof': 'B<0; A+2*B=c/pi>0; sin(pi*p) is the only zero factor; corner derivative signs alternate',
        'node_sign_intervals_exact': {'B_over_c': list(map(str, negative)),
            'minimum_positive_sine_factor_over_c': list(map(str, minimum)),
            'corner_derivative_over_c': list(map(str, corner))},
        'established_lattice_doubling_reference': 'https://doi.org/10.1016/0550-3213(81)90361-8',
        'second_shell': {'momentum': '2/3', 'continuum_kinetic_over_c': '2/3',
                         'two_hop_kinetic_over_c': str(second_value)},
        'physical_chiral_gauge_regulator_or_measured_universe_selected': False}


def _add(a, b):
    return a[0]+b[0], a[1]+b[1]


def _neg(a):
    return -a[1], -a[0]


def _mul(a, b):
    products = [x*y for x in a for y in b]
    return min(products), max(products)


def _inverse(a):
    if a[0] <= 0 <= a[1]:
        raise ValueError('interval inverse crosses zero')
    return 1/a[1], 1/a[0]


def _point(value):
    return F(value), F(value)


def _compact(interval):
    grid = 10**28
    return F((interval[0]*grid).__floor__(), grid), F((interval[1]*grid).__ceil__(), grid)


@lru_cache(maxsize=256)
def expression_interval(value):
    """Directed rational evaluation of this fixed source's log/radical data."""
    value = s.sympify(value)
    if value.is_Rational:
        return _point(F(int(value.p), int(value.q)))
    if value == s.pi:
        return _compact(pi_interval())
    if value.is_Add:
        result = _point(0)
        for item in value.args:
            result = _add(result, expression_interval(item))
        return _compact(result)
    if value.is_Mul:
        result = _point(1)
        for item in value.args:
            result = _mul(result, expression_interval(item))
        return _compact(result)
    if value.is_Pow and value.exp.is_Integer:
        exponent = int(value.exp)
        base = expression_interval(value.base)
        if exponent < 0:
            base, exponent = _inverse(base), -exponent
        result = _point(1)
        for _ in range(exponent):
            result = _mul(result, base)
        return _compact(result)
    if value.is_Pow and value.exp == s.Rational(1, 2):
        lo, hi = expression_interval(value.base)
        return root_interval(lo)[0], root_interval(hi)[1]
    if value.func == s.log:
        lo, hi = expression_interval(value.args[0])
        def logarithm(x):
            return log_interval(x) if x >= 1 else _neg(log_interval(1/x))
        return _compact((logarithm(lo)[0], logarithm(hi)[1]))
    raise ValueError('unsupported exact source expression: '+str(value))


def _trig(interval, *, cosine=False):
    midpoint, radius = sum(interval)/2, (interval[1]-interval[0])/2
    lo, hi = trig_interval(midpoint, cosine=cosine)
    return _compact((max(F(-1), lo-radius), min(F(1), hi+radius)))


def _cadd(a, b):
    return _add(a[0], b[0]), _add(a[1], b[1])


def _cmul(a, b):
    return _add(_mul(a[0], b[0]), _neg(_mul(a[1], b[1]))), _add(_mul(a[0], b[1]), _mul(a[1], b[0]))


def _conjugate(a):
    return a[0], _neg(a[1])


def _matrix_vector(matrix, vector):
    return tuple(_cadd(_cmul(row[0], vector[0]), _cmul(row[1], vector[1])) for row in matrix)


def _selected_mass():
    labels = action.mode_labels()
    chosen = next(row['index'] for row in labels if row['spin'] == 0
                  and row['internal_type'] == 'L' and row['family'] == 0
                  and s.Rational(row['electric_Q']) == -1)
    kernel = action.mass_kernel()
    partners = [i for i in range(45) if kernel[chosen, i] != 0]
    if len(partners) != 1:
        raise ArithmeticError('source charged mode no longer has one trial pair')
    clean = action.source_coefficients()['Omega']*action.source_coefficients()['J']*action.family_feature()[0, 0]
    if s.simplify(kernel[chosen, partners[0]]-clean) != 0:
        raise ArithmeticError('declared source mass and retained charged pair differ')
    return clean


def ideal_event_interval(law):
    """Exact continuum score: completed, clock_X=1, neutral pair read=0.

    The first pair number read supplies the SAME conditional feedback as
    the full Fock experiment. The second number read is summed, not dropped
    by conditioning. Two-mode even pair Hamiltonians commute on this z axis.
    """
    omega = expression_interval(action.source_coefficients()['Omega'])
    mass = expression_interval(_selected_mass())
    p = effective_momentum((0, 0, SECOND_SHELL), law)[2]
    frequency = expression_interval(completion_coefficients()[0]*p)
    energy2 = _add(_mul(mass, mass), _mul(frequency, frequency))
    energy = root_interval(energy2[0])[0], root_interval(energy2[1])[1]
    argument = _mul(energy, _point(INTERVAL))
    cosine, sine = _trig(argument, cosine=True), _trig(argument)
    scale = _mul(sine, _inverse(energy))
    z, x = _mul(frequency, scale), _mul(mass, scale)
    unitary = (((cosine, _neg(z)), (_point(0), _neg(x))),
               ((_point(0), _neg(x)), (cosine, z)))
    first = (unitary[0][0], unitary[1][0])
    moving = _mul(x, x)
    probabilities = (_add(_point(1), _neg(moving)), moving)
    pair_zero = _point(0)
    for bit in (0, 1):
        feedback = first if not bit else (first[0], (_neg(first[1][1]), first[1][0]))
        after = _matrix_vector(unitary, feedback)
        cross = _cmul(_conjugate(after[0]), after[1])[0]
        prob = _mul(_add(_point(1), _mul(_point(2), cross)), _point(F(1, 2)))
        pair_zero = _add(pair_zero, _mul(probabilities[bit], prob))
    calibration = protected.model.ideal_clock_calibration()
    contrast = F(calibration['populations'][0])-F(calibration['populations'][1])
    clock_cosine = _trig(_mul(omega, _point(2*INTERVAL)), cosine=True)
    clock_one = _mul(_add(_point(1), _neg(_mul(_point(contrast), clock_cosine))), _point(F(1, 2)))
    return _compact(_mul(_point(ACCEPTED**2), _mul(pair_zero, clock_one)))


@lru_cache(maxsize=1)
def _pair_detector():
    """Actual source neutral-pair X read, including fermion routing signs."""
    modes, a, b = 8, *PAIR_SITES
    state = np.eye(1 << modes, dtype=complex)
    words = []
    resource_keys = g._ADDITIVE_RESOURCES
    resources = {key: 0 for key in resource_keys}
    swap = g.local_gate('fswap')
    route = [(j-1, j) for j in range(b, a+1, -1)]
    calls = [*((swap, sites) for sites in route),
             (g.local_gate('pairing', -1j, math.pi/4), (a, a+1)),
             *((swap, sites) for sites in reversed(route))]
    for gate, sites in calls:
        global_sites = tuple(630+i for i in sites)
        literal = protected.decoded_source_word(gate.word, global_sites)
        lowered = protected._lower_decoded_word(literal, global_sites)
        if lowered != list(gate.word):
            raise ArithmeticError('new pair detector touches a wrong global source slot')
        matrix = protected.compiler.execute(lowered, 2)
        state = protected.model.apply_subsystem(state, matrix, sites, modes)
        cost = protected.compiler.resources(literal)
        for key in resource_keys:
            resources[key] += cost[key]
        words.extend(literal)
    h, delta = np.zeros((8, 8), complex), np.zeros((8, 8), complex)
    delta[a, b], delta[b, a] = -1j, 1j
    target = expm(-1j*math.pi/4*protected._fock_hamiltonian(h, delta))
    residual = float(np.linalg.norm(state-target, 2))
    if not math.isfinite(residual) or residual > 2e-10:
        raise ArithmeticError('source neutral pair detector differs from independent CAR target')
    return g._freeze(state), g._freeze(target), tuple(words), dict(resources), residual


def pair_detector():
    unitary, target, words, resources, residual = _pair_detector()
    return unitary, target, tuple({**event, 'targets': list(event['targets'])} for event in words), dict(resources), residual


def additional_read(packet, *, continuous=False):
    detector = pair_detector()[1 if continuous else 0]
    rows, arrays = [], []
    for old in packet['terminal']:
        changed = detector@old['probe_vector']
        outputs = protected._number_vectors(changed, PAIR_SITES[0], continuous=continuous)
        spectator = math.prod(float(np.vdot(vector, vector).real)
                              for key, vector in packet['spectator_states'].items() if key[2] == old['elapsed_cycles'])
        for bit, vector in enumerate(outputs):
            probability = (float(F(old['selector_weight_exact']))*spectator
                           *float(np.vdot(vector, vector).real)
                           *float(np.vdot(old['clock_reference_columns'], old['clock_reference_columns']).real))
            rows.append({key: value for key, value in old.items()
                         if key not in ('probe_vector', 'clock_reference_columns', 'probability')}
                        | {'pair_X': bit, 'probability': probability})
            arrays.append(g._freeze(vector))
    mass = math.fsum(row['probability'] for row in rows)
    if len(rows) != 80 or not math.isfinite(mass) or abs(mass-1) > 2e-7:
        raise ArithmeticError('additional complete neutral pair read lost an outcome or mass')
    score = math.fsum(row['probability'] for row in rows
                     if row['flag'] == 'completed' and row['clock_X'] == 1 and row['pair_X'] == 0)
    return {'terminal': rows, 'total_probability': mass, 'scored_joint_event_probability': score}, np.stack(arrays)


def execute(*, steps=16):
    protected._steps(steps)
    proof = exact_certificate()
    ideals = {law: ideal_event_interval(law) for law in ('continuum', 'two_hop_lattice')}
    separation = ideals['continuum'][0]-ideals['two_hop_lattice'][1]
    if separation <= F(1, 250):
        raise ArithmeticError('exact joint event separation missed its declared design floor')
    reports, arrays = {}, {}
    for law in ideals:
        physical = (s.Integer(0), s.Integer(0), SECOND_SHELL)
        source_p = tuple(float(value.evalf(40)) for value in effective_momentum(physical, law))
        block = protected.field_block(source_p, INTERVAL, steps=steps)
        actual = protected.protocol((block,), 0, INTERVAL)
        target = protected.protocol((block,), 0, INTERVAL, continuous=True)
        scored, probe_outputs = additional_read(actual)
        compared, reference_outputs = additional_read(target, continuous=True)
        score = compared['scored_joint_event_probability']
        if not float(ideals[law][0])-2e-11 <= score <= float(ideals[law][1])+2e-11:
            raise ArithmeticError('independent full Fock/clock law disagrees with exact two-pair interval')
        arrays[law+'_probe_outputs'] = probe_outputs
        arrays[law+'_clock_reference_columns'] = np.stack([row['clock_reference_columns'] for row in actual['terminal'] for _ in (0, 1)])
        for factor_index, factor in enumerate(block.factors):
            blank = np.eye(len(factor.unitary), dtype=complex)[:, 0]
            arrays[f'{law}_factor{factor_index}'] = np.stack((blank, factor.unitary@blank, factor.unitary@factor.unitary@blank))
        output_error = float(np.max(np.abs(probe_outputs-reference_outputs)))
        score_error = abs(scored['scored_joint_event_probability']-score)
        if not math.isfinite(output_error) or output_error > .01:
            raise ArithmeticError('complete source pair outputs exceed numerical admission')
        reports[law] = {'physical_coordinate_momentum': ['0', '0', '2/3'],
            'source_kinetic_vector': [str(x) for x in effective_momentum(physical, law)],
            'CAR_modes_in_executed_momentum_pair': 180,
            'all_internal_family_spin_and_opposite_momenta_retained': True,
            'complete_law': scored, 'ideal_joint_event_interval_exact': list(map(str, ideals[law])),
            'source_joint_event_difference_diagnostic': score_error,
            'complete_probe_output_max_entry_difference_diagnostic': output_error,
            'source_field_all_input_operator_error_diagnostic': block.tensor_operator_error,
            'source_field_resources_per_cycle': block.resources,
            'same_reference_clock_calibration': actual['clock_calibration'],
            'same_first_read_feedback': 'first=1 applies the same neutral i phase to mode1'}
    margin = F(1, 250)/2-F(1, 2000)
    width = float(2*margin)
    if F(width) > 2*margin:
        width = math.nextafter(width, -math.inf)
    sampling = concentration_certificate(1, 1_200_000, width)
    if F(sampling['familywise_error_upper']) >= F(1, 100):
        raise ArithmeticError('fixed all-attempt decision design missed one percent')
    return {'schema': 'oph.chiral_candidate.dispersion_identifiability.v1',
        'source_split_count': steps,
        'exact_restricted_field_certificate': proof, 'second_shell_source_comparison': reports,
        'ideal_joint_event_separation_lower_exact': str(separation),
        'fixed_decision_design': {'attempted_independent_trials': 1_200_000,
            'score': 'completed AND clock_X=1 AND pair_X=0; all other outcomes count as0',
            'threshold_exact': str((ideals['continuum'][0]+ideals['two_hop_lattice'][1])/2),
            'minimum_ideal_probability_gap': '1/250',
            'assumed_complete_event_error_per_model': '1/2000',
            'Bernoulli_probability_decision_margin_exact': str(margin),
            'sampling_bounded_observable': 'Y=2*event-1 in [-1,1]; certificate radius is at most twice the Bernoulli probability margin',
            'sampling_certificate': sampling,
            'source_numeric_or_physical_service_error_certified_below_assumed_error': False},
        'neutral_pair_detector': {'sites': list(PAIR_SITES), 'global_native_CAR_indices': [631, 636],
            'coherent_rotation_source_resources': pair_detector()[3],
            'complete_number_pointer_source_resources': protected.compiler.resources(g._occupation_pointer()[0]),
            'source_CAR_operator_error_diagnostic': pair_detector()[4]},
        'scope': {'first_shell_comparison_is_exact_on_restricted_observable_algebra': True,
            'larger_lattice_exterior_outputs_or_resource_costs_reproduced': False,
            'new180_mode_pair_executes_all_Fock_columns_and80_clock_detector_outputs': True,
            'new_pair_can_be_appended_to_the_fixed630_mode_field_as_disjoint180_mode_factor': True,
            'full810_mode_native_source_error_budget_reexecuted': False,
            'universal_local_chiral_regulator_or_physical_dictionary_derived': False,
            'physical_clock_momentum_or_detector_service_calibrated': False,
            'laboratory_observations_or_empirical_fit_used': False,
            'finite_numeric_diagnostics_are_outward_certificates': False}}, arrays


def source_bindings():
    names = set(protected.source_bindings())|{
        'model/common_source/chiral_candidate/dispersion_identifiability.py',
        'model/common_source/chiral_candidate/test_dispersion_identifiability.py',
        'model/common_source/halfline.py', 'model/calibrated_measurement/bounds.py',
        'model/characteristic_constraints/bounds.py',
        'model/scalar_observations/experimental_readout/stats.py'}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def canonical_report(report):
    result = deepcopy(report)
    for row in result['second_shell_source_comparison'].values():
        for key, gate in (('source_joint_event_difference_diagnostic', 1e-4),
                          ('complete_probe_output_max_entry_difference_diagnostic', .01),
                          ('source_field_all_input_operator_error_diagnostic', .01)):
            value = row.pop(key)
            if not math.isfinite(value) or not 0 <= value <= gate:
                raise ArithmeticError('fresh second-shell source diagnostic exceeds admission: '+key)
            row[key.replace('_diagnostic', '_admission')] = gate
    detector = result['neutral_pair_detector']
    residual = detector.pop('source_CAR_operator_error_diagnostic')
    if not math.isfinite(residual) or not 0 <= residual <= 2e-10:
        raise ArithmeticError('fresh neutral pair detector exceeds admission')
    detector['source_CAR_operator_error_admission'] = 2e-10
    result['backend_contract'] = {'raw_backend_residuals_authenticated': False,
        'complete_complex_output_comparison_tolerance': 2e-9,
        'fixed_admissions_checked_on_owned_fresh_execution': True,
        'admissions_are_outward_or_physical_certificates': False}
    return json.loads(json.dumps(result, allow_nan=False))


def check(path):
    """Own source replay and compare every retained complex factor and record."""
    path = Path(path)
    document = strict_json.load(path)
    json.dumps(document, allow_nan=False)
    if document.get('schema') != 'oph.chiral_candidate.dispersion_identifiability.v1':
        raise ValueError('wrong dispersion identifiability schema')
    if document.get('source_bindings') != source_bindings():
        raise ValueError('dispersion identifiability source custody changed')
    artifact = document['artifact']
    if set(artifact) != {'filename', 'sha256'} or Path(artifact['filename']).name != artifact['filename']:
        raise ValueError('one adjacent quantum artifact and exact hash required')
    array_path = path.with_name(artifact['filename'])
    if hashlib.sha256(array_path.read_bytes()).hexdigest() != artifact['sha256']:
        raise ValueError('dispersion quantum artifact custody changed')
    with np.load(array_path, allow_pickle=False) as archive:
        if len(archive.files) != len(set(archive.files)):
            raise ValueError('duplicate dispersion quantum outputs')
        stored = dict(archive)
    current, arrays = execute(steps=document['source_split_count'])
    fresh = canonical_report(current)
    fresh['source_bindings'], fresh['artifact'] = source_bindings(), artifact
    def compare(left, right, key='root'):
        if type(left) is not type(right):
            raise ValueError('dispersion receipt type mismatch: '+key)
        if type(right) is dict:
            if left.keys() != right.keys():
                raise ValueError('dispersion receipt menu differs: '+key)
            for name in right:
                compare(left[name], right[name], key+'.'+name)
        elif type(right) is list:
            if len(left) != len(right):
                raise ValueError('dispersion receipt length differs: '+key)
            for index, (a, b) in enumerate(zip(left, right, strict=True)):
                compare(a, b, key+f'[{index}]')
        elif type(right) is float:
            if not math.isfinite(left) or not math.isfinite(right):
                raise ValueError('nonfinite dispersion evidence')
            if 'admission' in key or 'tolerance' in key:
                if left != right:
                    raise ValueError('dispersion admission changed: '+key)
            elif not math.isclose(left, right, rel_tol=2e-9, abs_tol=2e-9):
                raise ValueError('dispersion numerical evidence differs: '+key)
        elif left != right:
            raise ValueError('dispersion exact evidence differs: '+key)
    compare(document, fresh)
    if stored.keys() != arrays.keys():
        raise ValueError('dispersion quantum output menu differs')
    for name, value in arrays.items():
        retained = stored[name]
        if (retained.shape != value.shape or retained.dtype.kind != 'c'
                or not np.isfinite(retained).all()
                or not np.allclose(retained, value, atol=2e-9, rtol=2e-9)):
            raise ValueError('complete complex dispersion output differs: '+name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--output', type=Path)
    choice.add_argument('--check', type=Path)
    parser.add_argument('--steps', type=int, default=16)
    args = parser.parse_args()
    if args.check:
        check(args.check)
        print('Verified restricted-field ambiguity and complete second-shell source outputs')
        return
    companion = args.output.with_suffix('.npz')
    if args.output.exists() or companion.exists():
        parser.error('choose fresh comparison and companion paths')
    report, arrays = execute(steps=args.steps)
    report = canonical_report(report)
    report['source_bindings'] = source_bindings()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with companion.open('xb') as stream:
        np.savez_compressed(stream, **arrays)
    report['artifact'] = {'filename': companion.name, 'sha256': hashlib.sha256(companion.read_bytes()).hexdigest()}
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({law: row['complete_law']['scored_joint_event_probability']
                      for law, row in report['second_shell_source_comparison'].items()}))


if __name__ == '__main__':
    main()
