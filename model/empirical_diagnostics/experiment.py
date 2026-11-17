"""Execute fixed preparation probes without receiving empirical observations."""
from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..anharmonic_observations import source as primitives
from ..anharmonic_transfer import calibration
from ..occupation_control import intervals as iv
from ..occupation_control.experiment import thermal_state, _add, _total
from ..occupation_control.solver import evolve
from ..observer.source_scalar import compiler
from ..primitive_observer import anharmonic, runtime
from ..scalar_observations.experimental_readout import executable_source as native


FOLDER = Path(__file__).resolve().parent
ROOT = FOLDER.parents[1]


def _real(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f'finite real {name} required')
    return float(value)


def gaussian_characteristic(alpha, thermal, u, v):
    alpha, thermal, u, v = (_real(x, n) for x, n in
                            ((alpha, 'alpha'), (thermal, 'thermal'), (u, 'u'), (v, 'v')))
    if thermal < 0:
        raise ValueError('nonnegative thermal occupation required')
    return np.exp(-(thermal+.5)*(u*u+v*v)+2j*alpha*v)


def identify_gaussian(value, h):
    """Unique inverse on -pi/(2h) < alpha <= pi/(2h), for trusted beta."""
    h = _real(h, 'probe coordinate')
    if h <= 0 or isinstance(value, (bool, np.bool_)):
        raise ValueError('positive coordinate and nonzero finite complex read required')
    value = complex(value)
    if not math.isfinite(value.real) or not math.isfinite(value.imag) or not 0 < abs(value) <= 1:
        raise ValueError('nonzero finite characteristic value in the unit disk required')
    angle = math.atan2(value.imag, value.real)
    if angle == -math.pi:
        angle = math.pi
    alpha = angle/(2*h)
    thermal = -math.log(abs(value))/(h*h)-.5
    if thermal < -1e-12:
        raise ValueError('read lies outside the displaced-thermal family')
    return {'alpha_real': alpha, 'thermal_occupation': thermal,
            'alpha_alias_spacing': math.pi/h,
            'phase_branch': '-pi/(2h) < alpha <= pi/(2h)'}


def coherent_mixture_weight(alpha, thermal, h, second_alpha=-1.5):
    target = gaussian_characteristic(alpha, thermal, 0., h).imag
    denominator = gaussian_characteristic(second_alpha, 0., 0., h).imag
    if denominator == 0:
        raise ValueError('the specified coherent component cannot match this read')
    weight = target/denominator
    if not 0 < weight < 1:
        raise ValueError('the specified pair cannot give a nontrivial positive mixture')
    return float(weight)


def _law(history):
    law = np.zeros(4)
    for branch in history.branches:
        law[branch.record[-1]] += float(np.trace(branch.density).real)
    if (not np.isfinite(law).all() or np.min(law) < -2e-10
            or abs(float(law.sum())-1) > 2e-9):
        raise ArithmeticError('complete source law is not normalized')
    return law


def native_component(dimension, alpha, thermal, probes):
    """Native preparation and four-outcome instrument in the shared runtime."""
    states, weights, preparation = native.preparation(dimension, alpha, thermal)
    density = np.zeros((2*dimension, 2*dimension), complex)
    density[:dimension, :dimension] = (states*weights) @ states.conj().T
    spin, spin_metadata = anharmonic._spin_record(dimension)
    records = []
    for probe in probes:
        if probe['quadrature'] not in ('real', 'imag'):
            raise ValueError('real or imaginary quadrature required')
        theta = 0. if probe['quadrature'] == 'real' else -math.pi/2
        detector, detector_metadata = primitives.detector(dimension, probe['u'], probe['v'], theta)
        groups = tuple((np.kron(np.eye(2), detector[bit]) @ spin[choice],)
                       for choice in (0, 1) for bit in (0, 1))
        instrument = runtime.Instrument(kraus=groups, visible=True,
            binding={'spin': spin_metadata, 'detector': detector_metadata, 'probe': probe})
        program = runtime.Program('fixed_preparation_readout', 2*dimension, density,
            (runtime.Event('spin_and_detector_read', instruments=(instrument,)),),
            binding={'preparation': preparation, 'no_evolution': True})
        history = runtime.execute(program)
        law = _law(history)
        records.append({'name': probe['name'], 'joint_probabilities': law.reshape(2, 2).tolist(),
                        'mean': float(law[0]+law[2]-law[1]-law[3]),
                        'mass': float(law.sum()),
                        'source_detector': detector_metadata,
                        'complete_quantum_branches': len(history.branches)})
    return {'preparation': preparation, 'spin_pointer': spin_metadata, 'records': records}


def semantic_component(dimension, alpha, thermal, probes):
    """Independent infinite-occupation certificate, not native source execution.

    Only axis probes occur in the locked menu. The pointer is initially |+>.
    Its off-diagonal trace after opposite half-displacements is chi(beta)/2.
    Measuring X or -Y gives its real or imaginary part respectively.
    """
    state, initial = thermal_state(dimension, thermal)
    state, preparation = evolve(state, cells=8, order=12, force=alpha, force_imaginary=True)
    pointer_cache = {}
    records = []
    for probe in probes:
        if probe['quadrature'] not in ('real', 'imag'):
            raise ValueError('real or imaginary quadrature required')
        u, v = probe['u'], probe['v']
        if u and v:
            raise ValueError('the certified diagnostic supports axis probes only')
        if (u, v) not in pointer_cache:
            joint = np.block([[state*.5, state*.5], [state*.5, state*.5]])
            target_half = iv.complex_scale(iv.complex_exact(state), iv.exact(.5))
            embedding = iv.complex_sub(iv.complex_exact(state*.5), target_half)
            embedding_error = float(iv.up(2*iv.trace_upper(*embedding)))
            joint, detector_budget = evolve(joint, cells=8, order=12, sectors=4,
                force=u/2 if u else -v/2, force_imaginary=bool(u), pointer_sign=True)
            detector_budget['pointer_embedding_rounding'] = embedding_error
            budget = _add([initial, preparation, detector_budget])
            pointer_cache[u, v] = joint, budget
        joint, budget = pointer_cache[u, v]
        masses, signals = [], []
        for spin in (0, 1):
            pairs = [(spin*dimension+i, (2+spin)*dimension+i) for i in range(dimension)]
            masses.append(sum((Fraction(joint[i, i].real)+Fraction(joint[j, j].real)
                               for i, j in pairs), Fraction()))
            signals.append(2*sum((Fraction(joint[i, j].real if probe['quadrature'] == 'real'
                                         else joint[i, j].imag) for i, j in pairs), Fraction()))
        exact_law = [(mass+sign*signal)/2 for mass, signal in zip(masses, signals, strict=True)
                     for sign in (1, -1)]
        law = [float(x) for x in exact_law]
        rounding = sum(abs(Fraction(x)-y) for x, y in zip(law, exact_law, strict=True))/2
        epsilon = _total(budget)
        tv = iv.sum_up([epsilon, iv.rational(rounding)[1]])
        finite = _total({key: value for key, value in budget.items()
                         if key not in ('thermal_tail', 'occupation_defect')})
        records.append({'name': probe['name'], 'joint_probabilities': np.asarray(law).reshape(2, 2).tolist(),
                        'total_variation_error_upper': tv,
                        'unrestricted_vs_exact_finite_half_trace_upper': iv.sum_up([epsilon, finite]),
                        'budget': budget})
    return records


def _source_bindings():
    paths = [FOLDER/'protocol.json', FOLDER/'experiment.py', FOLDER/'__main__.py',
             ROOT/'model/anharmonic_transfer/protocol.json',
             ROOT/'model/anharmonic_transfer/calibration.py',
             ROOT/'model/anharmonic_transfer/data/calibration_receipt.json.gz',
             ROOT/'model/anharmonic_observations/source.py',
             ROOT/'model/primitive_observer/runtime.py', ROOT/'model/primitive_observer/anharmonic.py',
             ROOT/'model/common_source/open_protocol.py', ROOT/'model/common_source/open_environment.py',
             ROOT/'model/common_source/protocol.py',
             ROOT/'model/occupation_control/experiment.py', ROOT/'model/occupation_control/solver.py',
             ROOT/'model/occupation_control/intervals.py',
             ROOT/'model/occupation_control/__main__.py',
             ROOT/'model/scalar_observations/experimental_readout/executable_source.py',
             ROOT/'model/scalar_observations/experimental_readout/canonical_source.py',
             ROOT/'model/scalar_observations/source.py',
             ROOT/'model/observer/source_scalar/compiler.py', ROOT/'tools/strict_json.py',
             ROOT/'model/requirements.txt', ROOT/'code/source_selection_model/response.json',
             ROOT/'code/m1_source_realization/receipt.json']
    return {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def produce(*, dimension=None, certificate=True):
    protocol = strict_json.load(FOLDER/'protocol.json')
    transfer = strict_json.load(ROOT/'model/anharmonic_transfer/protocol.json')
    supplied = calibration.transferred_parameters(transfer)
    alpha, thermal = supplied['alpha_real'], supplied['thermal_occupation']
    if alpha != protocol['gaussian_preparation']['alpha_real'] or thermal != protocol['gaussian_preparation']['thermal_occupation']:
        raise ValueError('diagnostic detached from the frozen preparation')
    dimension = protocol['dimension'] if dimension is None else native._dimension(dimension)
    probes, h = protocol['probe_menu'], protocol['probe_h']
    alternative_alpha = protocol['non_gaussian_control']['coherent_displacements'][1]
    weight = coherent_mixture_weight(alpha, thermal, h, alternative_alpha)
    components = [('gaussian', alpha, thermal), ('vacuum', 0., 0.), ('displaced_vacuum', alternative_alpha, 0.)]
    executed, certified = {}, {}
    for name, a, n in components:
        executed[name] = native_component(dimension, a, n, probes)
        if certificate:
            certified[name] = semantic_component(dimension, a, n, probes)
    records = []
    for index, probe in enumerate(probes):
        gaussian = np.asarray(executed['gaussian']['records'][index]['joint_probabilities'])
        mixture = ((1-weight)*np.asarray(executed['vacuum']['records'][index]['joint_probabilities'])
                   +weight*np.asarray(executed['displaced_vacuum']['records'][index]['joint_probabilities']))
        expected_g = gaussian_characteristic(alpha, thermal, probe['u'], probe['v'])
        expected_m = ((1-weight)*gaussian_characteristic(0., 0., probe['u'], probe['v'])
                      +weight*gaussian_characteristic(alternative_alpha, 0., probe['u'], probe['v']))
        means = [float(law[0, 0]+law[1, 0]-law[0, 1]-law[1, 1]) for law in (gaussian, mixture)]
        ideal_means = [float(value.real if probe['quadrature'] == 'real' else value.imag)
                       for value in (expected_g, expected_m)]
        row = {'name': probe['name'], 'gaussian_joint_probabilities': gaussian.tolist(),
               'non_gaussian_joint_probabilities': mixture.tolist(),
               'gaussian_mean': means[0], 'non_gaussian_mean': means[1],
               'ideal_unrestricted_means': ideal_means,
               'unrestricted_formula_residuals': [abs(a-b) for a, b in zip(means, ideal_means, strict=True)],
               'native_total_variation_contrast': float(np.sum(abs(gaussian-mixture))/2)}
        if certificate:
            semantic_g = certified['gaussian'][index]
            semantic_m = ((1-weight)*np.asarray(certified['vacuum'][index]['joint_probabilities'])
                          +weight*np.asarray(certified['displaced_vacuum'][index]['joint_probabilities']))
            mixture_bound = iv.sum_up([
                iv.mul(iv.rational(1-Fraction(weight)), iv.exact(certified['vacuum'][index]['total_variation_error_upper']))[1],
                iv.mul(iv.exact(weight), iv.exact(certified['displaced_vacuum'][index]['total_variation_error_upper']))[1],
                # Exact arithmetic probabilities and binary mixture weights are
                # rounded in forming the displayed combined law. Charge each
                # operation separately with intervals.
                _mixture_rounding(certified['vacuum'][index]['joint_probabilities'],
                                  certified['displaced_vacuum'][index]['joint_probabilities'], weight, semantic_m)])
            contrast = sum(abs(Fraction(float(a))-Fraction(float(b))) for a, b in
                           zip(np.asarray(semantic_g['joint_probabilities']).flat, semantic_m.flat, strict=True))/2
            lower = max(Fraction(0), contrast-Fraction(semantic_g['total_variation_error_upper'])-Fraction(mixture_bound))
            row['independent_semantic_check'] = {
                'gaussian_tv_error_upper': semantic_g['total_variation_error_upper'],
                'non_gaussian_tv_error_upper': mixture_bound,
                'native_vs_semantic_gaussian_max_probability': float(np.max(abs(gaussian-np.asarray(semantic_g['joint_probabilities'])))),
                'native_vs_semantic_non_gaussian_max_probability': float(np.max(abs(mixture-semantic_m))),
                'unrestricted_contrast_lower': iv.rational(lower)[0],
                'bound_applies_to_semantic_path_not_floating_source_execution': True}
        records.append(row)
    lookup = {row['name']: row for row in records}
    z = complex(lookup['position_h_real']['gaussian_mean'], lookup['position_h_imag']['gaussian_mean'])
    z2 = complex(lookup['position_2h_real']['gaussian_mean'], lookup['position_2h_imag']['gaussian_mean'])
    zm = complex(lookup['momentum_h_real']['gaussian_mean'], lookup['momentum_h_imag']['gaussian_mean'])
    checks = {'all_complete_laws_normalized': all(abs(sum(sum(row[key], []))-1) < 2e-9
                 for row in records for key in ('gaussian_joint_probabilities', 'non_gaussian_joint_probabilities')),
              'old_imaginary_probe_blind': lookup['position_h_imag']['native_total_variation_contrast'] < 2e-9,
              'extra_real_probe_separates': lookup['position_h_real']['native_total_variation_contrast'] > .02,
              'gaussian_inverse_recovers_preparation': abs(identify_gaussian(z, h)['alpha_real']-alpha) < 2e-9
                   and abs(identify_gaussian(z, h)['thermal_occupation']-thermal) < 2e-9,
              'gaussian_crosschecks_hold': abs(z2-abs(z)**2*z*z) < 2e-9 and abs(zm-abs(z)) < 2e-9}
    if certificate:
        checks['native_matches_independent_semantic_path'] = all(
            math.isfinite(row['independent_semantic_check'][key])
            and row['independent_semantic_check'][key] < 2e-10
            for row in records for key in ('native_vs_semantic_gaussian_max_probability',
                                           'native_vs_semantic_non_gaussian_max_probability'))
        checks['discriminator_survives_certified_occupation_error'] = (
            lookup['position_h_real']['independent_semantic_check']['unrestricted_contrast_lower'] > .07)
    if not all(checks.values()):
        raise ArithmeticError(f'preparation diagnostic failed: {checks}')
    return {'schema': 'oph.empirical_diagnostics.results.v1', 'dimension': dimension,
            'source_controls': compiler.source_summary(),
            'source_bindings': _source_bindings(), 'protocol': protocol,
            'frozen_preparation': supplied, 'non_gaussian_second_component_weight': weight,
            'source_executions': executed, 'semantic_certificates': certified,
            'records': records, 'recovered_gaussian': identify_gaussian(z, h),
            'checks': checks, 'all_pass': all(checks.values()),
            'scope': {'new_empirical_observations_loaded': False, 'old_score_changed': False,
                      'preparation_and_probe_coordinates_source_derived': False,
                      'semantic_occupation_error_certified': certificate,
                      'native_floating_error_certified': False,
                      'external_probe_scale_still_required': True}}


def _mixture_rounding(first, second, weight, actual):
    exact = [(1-Fraction(weight))*Fraction(a)+Fraction(weight)*Fraction(b)
             for a, b in zip(np.asarray(first).flat, np.asarray(second).flat, strict=True)]
    error = sum(abs(Fraction(float(x))-y) for x, y in zip(actual.flat, exact, strict=True))/2
    return iv.rational(error)[1]
