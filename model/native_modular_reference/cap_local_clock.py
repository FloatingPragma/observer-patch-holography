"""Bare-cap clock interventions without selecting an exterior configuration.

On the binary total-two shell, J_ab=sum_k(|a,k><b,k|-|b,k><a,k|)
with k outside {a,b} belongs to the bare two-port algebra. Its exponential
R_ab=exp(pi J_ab/4) rotates all ten spectator pairs equally. If a=1 and
b=7 or 8, both R_ab and the common reference's modular unitary belong to the
corresponding core-containing cap algebra. Complete native reads on (a,b)
also belong to that algebra. Their unread operation preserves every outside
observable, including arbitrary outside/reference coherences. Conditioning
on a local outcome may steer an initially correlated outside system.

The common state is the SAME locally prepared density, not a selected pure
pair. Its three core-charge weights are 2/225,1/36,2/45. Ramsey rotation,
modular action and inverse rotation followed by the native current give

 E[j]=17/200 cos(s log(25/8))+1/120 cos(s log(8/5)), P(j!=0)=181/900.

Both caps have this isolated read law. Their continuing joint history has
its own disturbed conditional states; no reset or equality of second
marginals is assumed. The second preparation rotation is fixed; only its
analyzer sign is selected by the first retained signed current. A closed-Gray uniformly controlled Y decomposition
compiles each two-address factor through the existing source services.
Complete original 512 arrays certify the full cap rotation. Encoded source
pulses may address the full register: locality belongs to the completed
operation on the native shell, not to its intermediate implementation.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from . import local_history as history, common_caps as common, local_preparation as local
from ..native_interface import construction as interface
from ..observer.source_scalar import compiler

full = common.full
ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('cap_local_clock_results.json')
SEAMS = history.SEAMS
PARAMETER = F(1)
ROTATION_GATE = F(1, 10**8)
WORD_GATE = F(1, 10**9)


def _probe(index, sign=1):
    if type(index) is not int or index not in (0, 1) or type(sign) is not int or sign not in (-1, 1):
        raise ValueError('one cap index and an integer rotation sign required')
    return index, sign


def pairs(index):
    """Every oriented occupation transfer, retaining all ten spectators."""
    _probe(index)
    states = common.regional.source()[0].states
    codes = full.support()[0]
    a, b = SEAMS[index]
    result = []
    for i, x in enumerate(states):
        if x[a] == 1 and x[b] == 0:
            y = list(x); y[a], y[b] = 0, 1
            j = states.index(tuple(y))
            result.append((i, j, codes[i], codes[j]))
    if len(result) != 10 or len({v for row in result for v in row[:2]}) != 20:
        raise ArithmeticError('native local transfer lost a spectator pair')
    return tuple(result)


def pair_word(a, b, sign):
    """Exact-target Gray conditional-Y word on the complete nine-bit chart."""
    if (type(a) is not int or type(b) is not int or not 0 <= a < 512
            or not 0 <= b < 512 or a == b or type(sign) is not int or sign not in (-1, 1)):
        raise ValueError('two distinct original code addresses and signed rotation required')
    different = a ^ b
    pivot = (different & -different).bit_length()-1
    controls = tuple(j for j in range(9) if j != pivot)
    fan = [pulse for bit in controls if different >> bit & 1
           for pulse in compiler.cnot(pivot, bit)]
    transformed_a = a ^ ((different ^ (1 << pivot)) if a >> pivot & 1 else 0)
    orientation = 1 if not (transformed_a >> pivot & 1) else -1
    control_word = sum(((transformed_a >> bit) & 1) << j for j, bit in enumerate(controls))
    tape = list(fan)
    for k in range(256):
        gray = k ^ (k >> 1)
        angle = sign*orientation*math.pi/1024*(-1)**((control_word & gray).bit_count())
        tape.append(compiler.pulse('Y', angle, [pivot]))
        successor = (k+1) % 256
        changed = gray ^ (successor ^ (successor >> 1))
        tape.extend(compiler.cnot(controls[changed.bit_length()-1], pivot))
    tape.extend(compiler.inverse(fan))
    return tape


def rotation_word(index, sign):
    _probe(index, sign)
    return [pulse for _, _, a, b in pairs(index) for pulse in pair_word(a, b, sign)]


@lru_cache(maxsize=4)
def _rotation(index, sign):
    # Execute every declared factor on ALL original columns, then compose
    # those unmodified source arrays. This evaluates the same pulse product
    # in blocks and avoids propagating floating subnormals through a long
    # monolithic tape. No polar correction, projection or normalization is
    # used; the complete final array is certified against its exact target.
    a = np.eye(512, dtype=complex)
    for _, _, left, right in pairs(index):
        a = compiler.execute(pair_word(left, right, sign), 9)@a
    return a.dtype.str, a.tobytes()


def rotation(index, sign):
    _probe(index, sign)
    dtype, data = _rotation(index, sign)
    return np.frombuffer(data, dtype=np.dtype(dtype)).reshape(512, 512)


def ideal_rotation(index, sign):
    _probe(index, sign)
    matrix = np.eye(512, dtype=complex)
    for _, _, a, b in pairs(index):
        matrix[np.ix_((a, b), (a, b))] = np.array([[1, sign], [-sign, 1]])/math.sqrt(2)
    return matrix


def certify_rotation(index, sign, array=None):
    """Rational error bound for every original entry, including padding."""
    _probe(index, sign)
    array = rotation(index, sign) if array is None else array
    if type(array) is not np.ndarray:
        raise ValueError('an original base array required')
    array = full.exact._array(array, (512, 512))
    lo, hi = full.exact._sqrt_interval(F(1, 2))
    middle, radius = (lo+hi)/2, (hi-lo)/2
    special = {}
    for _, _, a, b in pairs(index):
        special.update({(a, a): middle, (b, b): middle,
                        (a, b): sign*middle, (b, a): -sign*middle})
    squared = F()
    for i, row in enumerate(array):
        for j, value in enumerate(row):
            re, im = full.exact._component(value)
            wanted = special.get((i, j), F(int(i == j)))
            squared += (re-wanted)**2+im*im
    return full.exact._sqrt_interval(squared)[1]+full.exact._sqrt_interval(F(40))[1]*radius


def _in_algebra(entries, cap):
    supports = common.regional.basis_supports(cap, completed=False)
    covered = set()
    for support in supports:
        value = entries.get(support[0], 0)
        if any(entries.get(address, 0) != value for address in support):
            return False
        covered.update(support)
    return all(address in covered for address in entries)


def exact_algebra():
    """Exact generators and Kraus supports, not a floating locality check."""
    sigma = local.reference()
    rows = []
    for index, seam in enumerate(SEAMS):
        generator = {address: sign for i, j, _, _ in pairs(index)
                     for address, sign in (((i, j), 1), ((j, i), -1))}
        if not _in_algebra(generator, seam):
            raise ArithmeticError('the uniform generator is not a bare two-port operation')
        # The skew generator squares to minus its active projection, hence
        # exp(theta J)=I+sin(theta)J+(1-cos(theta))J^2 is in the same algebra.
        cap = common.CAPS[index]
        if not _in_algebra({(i, i): w for i, w in enumerate(sigma)}, cap):
            raise ArithmeticError('the same modular reference is not in the cap algebra')
        laws = {}
        for law in ('C', 'E'):
            branches = interface.fixed_maps(law, seam)
            for branch in branches:
                # Every nonzero Kraus entry is the same scalar 1/sqrt2.
                entries = {tuple(map(int, address)): 1 for address in np.argwhere(branch.matrix != 0)}
                if not _in_algebra(entries, seam):
                    raise ArithmeticError('native current instrument is not a bare local instrument')
            laws[law] = len(branches)
        i, j, _, _ = pairs(index)[0]
        selected = {(i, j): 1, (j, i): -1}
        if _in_algebra(selected, cap):
            raise ArithmeticError('selected-exterior countercontrol unexpectedly became cap local')
        rows.append({'cap': list(cap), 'native_read_seam': list(seam),
                     'spectator_pairs': 10, 'native_kraus_branches': laws,
                     'selected_exterior_generator_fails_bare_cap_membership': True})
    return {'rows': rows, 'completed_net_operations_in_bare_cap_algebras': True,
            'outside_marginal_preserved_by_unread_instrument': True,
            'conditioned_outside_state_may_be_steered': True,
            'intermediate_encoded_pulses_proved_physically_local': False}


def current_law(s):
    """The ideal first current distribution on the common prepared reference."""
    s = full.parameter(s)
    mean = F(17, 200)*math.cos(float(s)*math.log(25/8))+F(1, 120)*math.cos(float(s)*math.log(8/5))
    moving = F(181, 900)
    return {-1: (float(moving)-mean)/2, 0: 1-float(moving), 1: (float(moving)+mean)/2}


def calibration_certificate():
    """A finite-window inverse for the complete-trial expected current.

    Preparation failures and a pre-read export abort count as current zero.
    The ideal mean is q*m(s), q=5/33. On 1/2<=s<=1 both phases are in (0,pi/2).
    Using sin(x)>=2*x/pi, pi<22/7 and exact lower logarithm enclosures gives
    -d(q*m)/ds >1/200. Thus the inverse mean-readout map is 200-Lipschitz.
    This is a mean calibration, not monotonicity of individual random records.
    Outside the declared window the cosine response has orientation aliases.
    """
    lower = (F(113, 100), F(47, 100))
    for ratio, bound in zip((F(25, 8), F(8, 5)), lower, strict=True):
        lo, hi = common.logarithm_interval(ratio)
        if not bound < lo < hi < F(3, 2):
            raise ArithmeticError('the exact clock calibration phase bounds failed')
    slope = F(5, 33)*F(7, 22)*(F(17, 200)*lower[0]**2+F(1, 120)*lower[1]**2)
    if slope <= F(1, 200):
        raise ArithmeticError('the common clock inverse lost its finite slope')
    return {'modular_parameter_window': ['1/2', '1'],
            'ideal_herald_probability': '5/33',
            'complete_trial_expected_current': '(5/33)*(17/200*cos(s*log(25/8))+1/120*cos(s*log(8/5)))',
            'absolute_mean_derivative_lower': str(slope),
            'inverse_mean_readout_lipschitz_upper': 200,
            'individual_random_records_strictly_monotone': False,
            'global_orientation_or_phase_aliases_removed': False,
            'empirical_sampling_uncertainty_or_physical_service_fidelity_supplied': False}


@lru_cache(maxsize=1)
def _source_budget():
    for index in (0, 1):
        for sign in (-1, 1):
            full.corner.report._gate(certify_rotation(index, sign), ROTATION_GATE, 'uniform bare-cap rotation')
    full.corner.report._gate(common.certify_diagonal(PARAMETER, history.ratios(), history.clock()),
                            WORD_GATE, 'shared local modular word')
    error, gain = F(), F(1)
    for eps in (ROTATION_GATE, WORD_GATE, ROTATION_GATE)*2:
        error = (1+eps)**2*error+eps*(2+eps)/2
        gain *= (1+eps)**2
    pulses = sum(2*len(rotation_word(i, 1))+len(history.word()) for i in (0, 1))
    return tuple({'complete_history_half_diamond_upper': error, 'channel_norm_upper': gain,
                  'maximum_source_pulses': pulses, 'maximum_proper_code_service_events': 2*pulses,
                  'native_current_reads': 2, 'native_export_checks': 2,
                  'exterior_configuration_selection': False,
                  'initial_pure_probe_or_projective_work_bit_read': False,
                  'unmodified_full_source_factor_arrays_composed': True}.items())


def source_budget():
    return dict(_source_budget())


@lru_cache(maxsize=16)
def _stage(index, sign, raw, s):
    r = rotation if raw else ideal_rotation
    u = history.clock(s) if raw else np.diag([complex(float(w))**(-1j*float(s)) for w in history.ratios()])
    a = r(index, sign)@u@r(index, -1)
    return a.dtype.str, a.tobytes()


def _read(density, index, sign, law, raw, s):
    dtype, data = _stage(index, sign, raw, s)
    stage = np.frombuffer(data, dtype=np.dtype(dtype)).reshape(512, 512)
    return _read_stage(density, index, law, stage)


def _read_stage(density, index, law, stage):
    """Complete native/export instrument for an already certified source stage.

    Internal composition helper: callers own stage admission and its error
    certificate. All padding output is retained as a reference-bearing abort.
    """
    ref = density.shape[0]//66
    codes = np.array(full.support()[0])
    outside = np.array([i for i in range(512) if i not in set(codes)])
    action = np.kron(stage[np.ix_(codes, codes)], np.eye(ref))
    after = action@density@action.conj().T
    escape = np.kron(stage[np.ix_(outside, codes)], np.eye(ref))
    aborted = (escape@density@escape.conj().T).reshape(len(outside), ref, len(outside), ref)
    abort = np.einsum('aras->rs', aborted)
    outputs = {}
    for branch in interface.fixed_maps(law, SEAMS[index]):
        k = np.kron(branch.matrix, np.eye(ref))
        outputs[(branch.coin, branch.label, branch.current)] = k@after@k.conj().T
    return outputs, abort


def execute(density, law, *, raw=True, s=PARAMETER, feedback=True):
    density = history._density(density)
    if (type(law) is not str or law not in ('C', 'E') or type(raw) is not bool
            or type(feedback) is not bool):
        raise ValueError('one native law and original Boolean switches required')
    s = full.parameter(s)
    if raw and s != PARAMETER:
        raise ValueError('raw histories are certified at the retained parameter only')
    def stage(index, sign):
        dtype, data = _stage(index, sign, raw, s)
        return np.frombuffer(data, dtype=np.dtype(dtype)).reshape(512, 512)
    return _execute_stages(density, law, stage, feedback=feedback)


def _execute_stages(density, law, stage, *, feedback=True):
    """Compose admitted stages without resetting the carrier between reads."""
    first, abort = _read_stage(density, 0, law, stage(0, 1))
    result = {('abort', 0): abort}
    for record, rho in first.items():
        sign = -1 if feedback and record[2] < 0 else 1
        second, abort = _read_stage(rho, 1, law, stage(1, sign))
        result['abort', 1, *record] = abort
        for later, final in second.items():
            result['read', *record, *later] = final
    if any(not np.isfinite(rho).all() for rho in result.values()):
        raise ArithmeticError('complete cap-local history is nonfinite')
    return result


def prepared_history(law, *, raw=True):
    maps = local.population_maps(law=law)
    result = {('preparation', flag): np.diag(matrix[:, 0]).astype(complex)
              for flag, matrix in maps.items() if flag != 'success'}
    result.update(execute(np.diag(maps['success'][:, 0]).astype(complex), law, raw=raw))
    return result


def target_history(law):
    sigma = np.diag(list(map(float, local.reference()))).astype(complex)
    raw_reference = local.incident.reference(local.CORE)[0]
    failure = np.diag([float(F(1, 66)-F(5, 33)*w) for w in raw_reference]).astype(complex)
    result = {('preparation', 'herald_failure'): failure,
              ('preparation', 'warmup_timeout'): np.zeros((66, 66), complex),
              ('preparation', 'completion_timeout'): np.zeros((66, 66), complex)}
    result.update(execute(float(F(5, 33))*sigma, law, raw=False))
    return result


def run():
    budget = source_budget()
    prep = local.budget()
    total_bound = budget['channel_norm_upper']*prep['complete_flagged_reference_halftrace_upper']+F(5, 33)*budget['complete_history_half_diamond_upper']
    sigma = np.diag(list(map(float, local.reference()))).astype(complex)
    rows = []
    for law in ('C', 'E'):
        raw, ideal = execute(sigma, law), execute(sigma, law, raw=False)
        p_feedback = history.probabilities(ideal)
        p_fixed = history.probabilities(execute(sigma, law, raw=False, feedback=False))
        feedback_effect = sum(abs(p_feedback[k]-p_fixed[k]) for k in p_feedback)/2
        if feedback_effect <= 1e-3:
            raise ArithmeticError('retained native records lost analyzer feedback')
        prepared = prepared_history(law)
        discrepancy = history.distance(prepared, target_history(law))
        if discrepancy > float(total_bound) or history.distance(raw, ideal) > float(budget['complete_history_half_diamond_upper']):
            raise ArithmeticError('complete local operation exceeds its certificate')
        first = []
        for index in (0, 1):
            read, abort = _read(sigma, index, 1, law, False, PARAMETER)
            law_p = {j: sum(float(np.trace(rho).real) for key, rho in read.items() if key[2] == j) for j in (-1, 0, 1)}
            if max(abs(law_p[j]-current_law(PARAMETER)[j]) for j in law_p) > 2e-14 or np.linalg.norm(abort)>1e-14:
                raise ArithmeticError('local current does not have the counted two-gap law')
            first.append({str(j): p for j, p in law_p.items()})
        rows.append({'law': law, 'isolated_cap_current_laws': first,
                     'complete_prepared_history_halftrace_diagnostic': discrepancy,
                     'retained_current_feedback_total_variation_diagnostic': feedback_effect,
                     'complete_prepared_history_mass_diagnostic': sum(history.probabilities(prepared).values()),
                     'zero_clock_quantum_distance_diagnostic': history.distance(ideal, execute(sigma, law, raw=False, s=F()))})
    paths = sorted(set((*common.SOURCE_PATHS, *interface.SOURCE_PATHS, *local.SOURCE_PATHS,
                        'model/native_modular_reference/local_history.py',
                        'model/native_modular_reference/cap_local_clock.py')))
    return {'schema': 'oph.native_modular_reference.cap_local_clock.v1', 'all_pass': True,
            'parameter': str(PARAMETER), 'source_budget': full.preparation._compact(budget),
            'preparation': local._json(prep), 'exact_algebra': exact_algebra(),
            'mean_clock_calibration': calibration_certificate(),
            'complete_preparation_and_history_halftrace_upper': full.compiled._outward(total_bound),
            'source_clock_vs_actual_prepared_modular_halfdiamond_upper': full.compiled._outward(
                WORD_GATE*(2+WORD_GATE)/2+prep['prepared_modular_generator_operator_error_upper']),
            'ideal_mean_current': '17/200*cos(s*log(25/8))+1/120*cos(s*log(8/5))',
            'ideal_movement_probability': '181/900', 'histories': rows,
            'source_words': {'clock': hashlib.sha256(full.corner.canonical(history.word()).encode()).hexdigest(),
                'rotations': [hashlib.sha256(full.corner.canonical(rotation_word(i, sign)).encode()).hexdigest()
                              for i in (0, 1) for sign in (-1, 1)]},
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
            'scope': {'all_exterior_states_and_untouched_references': True,
                      'ideal_cap_algebra_action_at_completed_protocol_boundaries': True,
                      'coherent_encoded_control_services_remain_admitted': True,
                      'intermediate_pulse_locality_or_physical_signal_cone': False,
                      'conditional_remote_state_independence': False,
                      'physical_clock_temperature_or_spacetime_identification': False}}


def compare_receipt(supplied, actual):
    full.corner.compare_receipt(supplied, actual)
    for receipt in (supplied, actual):
        for row in receipt['histories']:
            for values in row['isolated_cap_current_laws']:
                ps = list(values.values())
                if (set(values) != {'-1', '0', '1'} or any(type(p) is not float or not math.isfinite(p) or not 0<=p<=1 for p in ps)
                        or abs(sum(ps)-1)>2e-13):
                    raise ValueError('a complete nonnegative local current law required')
            for key, value in row.items():
                if key.endswith('_diagnostic') and (type(value) is not float or not math.isfinite(value) or value<0):
                    raise ValueError('nonnegative finite complete-history diagnostics required')
            if abs(row['complete_prepared_history_mass_diagnostic']-1)>2e-10:
                raise ValueError('complete prepared history mass changed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path); group.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = strict_json.loads(full.preparation.canonical(run()))
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Bare-cap clock instruments and complete local history verified.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(full.corner.canonical(result).rstrip('\n')+'\n')
    else:
        print(full.corner.canonical(result))


if __name__ == '__main__':
    main()
