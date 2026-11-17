"""Native open-field ground state and explicitly closed-history Gibbs states.

The transfer acts on reconstructed configuration amplitudes, not on the
native quantum working algebra. Endpoint closure is an additional full-state
classical read service. Every failure/timeout remains in the native output.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations
import json
from pathlib import Path
from typing import NamedTuple

import numpy as np
import sympy as sp

from tools import strict_json
from model.observer import native
from model.record_laws import burn_in
from model.native_clock import clock as native_clock
from . import core

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*core.SOURCE_PATHS, *burn_in.SOURCE_PATHS, *native_clock.SOURCE_PATHS,
    'model/anharmonic_transfer/__init__.py',
    'model/anharmonic_transfer/environment.py',
    'model/common_source/smooth_protocol.py',
    'model/primitive_observer/__init__.py',
    'model/relational_event_clock/__init__.py',
    'model/relativistic_scalar/__init__.py',
    'model/relativistic_scalar/lattice.py',
    'model/scalar_observations/__init__.py',
    'model/scalar_observations/experimental_readout/__init__.py',
    'model/scalar_observations/experimental_readout/canonical_source.py',
    'model/scalar_observations/experimental_readout/executable_source.py',
    'model/scalar_observations/lattice.py',
    'model/scalar_observations/source.py',
    'model/observer_time_fields/thermal.py',
    'Lean/QFT/StructuralInheritance.lean'))))
DIMENSION = 66
TRANSFER_RATE_UPPER = F(35237, 38000)


def _integer(value, minimum, maximum, name):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(f'{name} must be an original integer in [{minimum},{maximum}]')
    return value


def _fraction(value):
    return F(int(value.p), int(value.q)) if isinstance(value, sp.Rational) else F(value)


def _array(rows):
    return np.array(rows, dtype=object)


def _tuples(matrix):
    return tuple(tuple(int(x) for x in row) for row in matrix)


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def interval(cap=None):
    """Unnormalized success/timeout kernels at one true incident interval."""
    data = core.source()
    if cap is None:
        kernel = core.completed()
        return kernel.counts, None, kernel.denominator
    _integer(cap, 1, 32, 'attempt cap')
    a, b = _array(data.remote_counts), _array(data.incident_counts)
    power = np.eye(DIMENSION, dtype=object)
    success = np.zeros((DIMENSION, DIMENSION), dtype=object)
    for j in range(cap):
        success += (b@power)*data.denominator**(cap-j-1)
        power = a@power
    denominator = data.denominator**cap
    mass = denominator-(data.denominator-len(data.incident_edges)*2)**cap
    assert all(sum(row) == mass for row in success)
    assert all(sum(row) == denominator-mass for row in power)
    return _tuples(success), _tuples(power), denominator


class LoopLaw(NamedTuple):
    ports: tuple
    cap: int | None
    feedback: bool
    records: tuple
    counts: tuple
    denominator: int

    def probabilities(self):
        """Complete public closure/field outcomes; endpoint tables are retained."""
        out = {}
        for (kind, word), matrix in zip(self.records, self.counts, strict=True):
            total = sum(map(sum, matrix))
            if kind == 'timeout':
                out[kind, word] = F(total, self.denominator)
            else:
                returned = sum(matrix[i][i] for i in range(DIMENSION))
                out['return', word] = F(returned, self.denominator)
                out['different', word] = F(total-returned, self.denominator)
        return out


def _protocol(ports, cap, feedback):
    if (type(ports) is not tuple or not 1 <= len(ports) <= 4
            or any(type(p) is not int or p not in core.source().region for p in ports)):
        raise ValueError('one through four read ports in the declared region required')
    if cap is not None:
        _integer(cap, 1, 32, 'attempt cap')
    if type(feedback) is not bool:
        raise ValueError('an original Boolean feedback policy required')


def loop_law(ports=(0, 1), *, cap=None, feedback=False):
    """All words and original/final configurations, before closure conditioning.

    A full initial configuration is saved solely for the final equality test.
    Local field bits follow the usual disclosed-bit/current dictionary. The
    optional existing internal repair occurs after a first read of one and
    itself costs an additional native incident tick.
    """
    _protocol(ports, cap, feedback)
    return _loop_law(ports, cap, feedback)


@lru_cache(maxsize=12, typed=True)
def _loop_law(ports, cap, feedback):
    rows, waits, step_denominator = interval(cap)
    success = _array(rows)
    timeout = None if waits is None else _array(waits)
    internal = _array(core.source().internal_counts)
    active = {(): np.eye(DIMENSION, dtype=object)}
    outputs = {}
    denominator = DIMENSION
    for cursor, port in enumerate(ports):
        factor = 2 if feedback and cursor == 0 else 1
        full_step = step_denominator*factor
        outputs = {key: matrix*full_step for key, matrix in outputs.items()}
        following = {}
        mask = np.array([state[port] for state in core.source().states], dtype=object)
        for word, matrix in active.items():
            if timeout is not None:
                outputs['timeout', word] = (timeout@matrix)*factor
            evolved = success@matrix
            for bit in (0, 1):
                chosen = evolved*(mask if bit else 1-mask)[:, None]
                if factor == 2:
                    chosen = internal@chosen if bit else 2*chosen
                following[(*word, bit)] = chosen
        active = following
        denominator *= full_step
    outputs.update({('completed', word): matrix for word, matrix in active.items()})
    records = tuple(sorted(outputs))
    counts = tuple(_tuples(outputs[key]) for key in records)
    if sum(sum(map(sum, matrix)) for matrix in counts) != denominator:
        raise ArithmeticError('closure experiment lost a completion or timeout')
    return LoopLaw(ports, cap, feedback, records, counts, denominator)


@lru_cache(maxsize=1)
def _native_branches():
    """Independent endpoint/coin expansion, without producer transfer counts."""
    states = tuple(state for state in native.states() if max(state) == 1)
    index = {state: i for i, state in enumerate(states)}
    seams = tuple(core.quantum.mixing.seams())
    branches = []
    for state in states:
        row = []
        for seam in seams:
            for coin in (0, 1):
                final, (_, current) = native.native(state, seam, coin)
                # The field updates use actual signed currents, including zero.
                for port in (0, 1):
                    expected = state[port]+current*(int(port == seam[1])-int(port == seam[0]))
                    if expected != final[port]:
                        raise ArithmeticError('the independent current oracle lost local readback')
                row.append((bool(set(seam)&{0, 1}), index[final]))
        branches.append(tuple(row))
    if states != core.source().states or len(seams) != 30:
        raise ArithmeticError('oracle and producer do not name the same native source')
    feedback = tuple(tuple(index[native.native(state, (0, 1), coin)[0]] for coin in (0, 1))
                     for state in states)
    return states, tuple(branches), feedback


def native_oracle(ports=(0, 1), *, cap=2, feedback=False):
    """Direct bounded native histories, including every endpoint and timeout.

    This routine never calls interval(), completed(), or a matrix recurrence.
    Integer padding only puts differently stopped words over one denominator.
    """
    _protocol(ports, cap, feedback)
    if cap is None or cap > 4 or len(ports) > 3:
        raise ValueError('the independent finite oracle accepts caps 1..4 and at most three reads')
    states, transitions, internal = _native_branches()
    active = {(i, i, ()): 1 for i in range(DIMENSION)}
    outputs = defaultdict(int)
    denominator = DIMENSION
    for cursor, port in enumerate(ports):
        factor = 2 if feedback and cursor == 0 else 1
        outputs = defaultdict(int, {key: n*60**cap*factor for key, n in outputs.items()})
        waiting = active
        completed = defaultdict(int)
        for attempt in range(1, cap+1):
            following = defaultdict(int)
            for (origin, state, word), count in waiting.items():
                for incident, final in transitions[state]:
                    if incident:
                        bit = states[final][port]
                        coefficient = count*60**(cap-attempt)
                        if factor == 2 and bit:
                            for target in internal[final]:
                                completed[origin, target, (*word, bit)] += coefficient
                        else:
                            completed[origin, final, (*word, bit)] += coefficient*factor
                    else:
                        following[origin, final, word] += count
            waiting = following
        for (origin, final, word), count in waiting.items():
            outputs['timeout', word, origin, final] += count*factor
        active = completed
        denominator *= 60**cap*factor
    for (origin, final, word), count in active.items():
        outputs['completed', word, origin, final] += count
    # Include zero-probability field words as well as every reachable endpoint.
    expected_records = []
    from itertools import product
    for depth in range(len(ports)):
        expected_records.extend(('timeout', tuple(word)) for word in product((0, 1), repeat=depth))
    expected_records.extend(('completed', tuple(word)) for word in product((0, 1), repeat=len(ports)))
    records = tuple(sorted(expected_records))
    counts = tuple(tuple(tuple(outputs[kind, word, origin, final] for origin in range(DIMENSION))
                         for final in range(DIMENSION)) for kind, word in records)
    if sum(sum(map(sum, matrix)) for matrix in counts) != denominator:
        raise ArithmeticError('independent native history lost mass')
    return LoopLaw(ports, cap, feedback, records, counts, denominator)


@lru_cache(maxsize=5)
def _powers(length):
    kernel = core.completed()
    power = np.eye(DIMENSION, dtype=object)
    values = [power]
    for _ in range(length):
        power = _array(kernel.counts)@power
        values.append(power)
    return tuple(values)


def open_moments(maximum=4):
    _integer(maximum, 1, 8, 'maximum lag')
    kernel = core.completed()
    field = [6*state[0]-1 for state in core.source().states]
    moments = []
    for n, matrix in enumerate(_powers(maximum)):
        numerator = sum(field[i]*int(matrix[i, j])*field[j]
                        for i in range(DIMENSION) for j in range(DIMENSION))
        moments.append(F(numerator, 36*DIMENSION*kernel.denominator**n))
    if moments[0] != F(5, 36) or moments[1] != F(115, 1188):
        raise ArithmeticError('the native local-field KMS obstruction changed')
    return tuple(moments)


def kms_defect_lower(beta):
    if type(beta) not in (int, F) or beta <= 0:
        raise ValueError('a positive exact inverse-temperature parameter required')
    return min(F(beta), F(1))*(open_moments(1)[0]-open_moments(1)[1])


def thermal_cycle(length=2):
    """Exact unforced Euclidean circle correlators, with beta=loop length."""
    _integer(length, 1, 8, 'loop length')
    powers = _powers(length)
    z_count = sum(int(powers[length][i, i]) for i in range(DIMENSION))
    partition = F(z_count, core.completed().denominator**length)
    occupied = [i for i, state in enumerate(core.source().states) if state[0]]
    mean = F(sum(int(powers[length][i, i]) for i in occupied), z_count)
    correlation = tuple(F(sum(int(powers[length-n][i, j])*int(powers[n][j, i])
                              for i in occupied for j in occupied), z_count)
                        for n in range(length+1))
    if correlation != tuple(reversed(correlation)) or correlation[0] != mean:
        raise ArithmeticError('native closed cylinders lost the Euclidean KMS wrap')
    return {'length': length, 'partition': partition, 'return_probability': partition/DIMENSION,
            'local_occupation_mean': mean, 'two_point_cycle': correlation,
            'open_mean': F(1, 6), 'mean_differs_from_open': mean != F(1, 6)}


@lru_cache(maxsize=1)
def _inverse_counts():
    data, kernel = core.source(), core.completed()
    inverse = (60*sp.eye(DIMENSION)-sp.Matrix(data.remote_counts))*sp.Matrix(data.incident_counts).inv(method='DM')
    denominator = int(sp.ilcm(*(entry.q for entry in inverse)))
    counts = np.array([[int(entry*denominator) for entry in row] for row in inverse.tolist()], dtype=object)
    original = _array(kernel.counts)
    identity = kernel.denominator*denominator*np.eye(DIMENSION, dtype=object)
    if not np.array_equal(original@counts, identity) or not np.array_equal(counts@original, identity):
        raise ArithmeticError('the exact native transfer inverse failed')
    return _tuples(counts), denominator


def boundary_certificate(length=2):
    """Finite KMS trace boundary checked in the similar native coordinates.

    K is self-adjoint in the D^-1 metric. Its density notation is not an
    assertion of ordinary-matrix positivity of the nonsymmetric K^length.
    The positive state is T^length/Tr(T^length) in reconstructed coordinates.
    """
    _integer(length, 1, 4, 'KMS loop length')
    inverse, inverse_denominator = _inverse_counts()
    power = _powers(length)[length]
    reverse = np.eye(DIMENSION, dtype=object)
    for _ in range(length):
        reverse = _array(inverse)@reverse
    f0 = np.array([state[0] for state in core.source().states], dtype=object)
    f1 = np.array([state[1] for state in core.source().states], dtype=object)
    checked = []
    for name, a, b in (('F0,F0', np.diag(f0), np.diag(f0)),
                       ('F0,F1', np.diag(f0), np.diag(f1)),
                       ('2*internal_repair,F0', _array(core.source().internal_counts), np.diag(f0))):
        left = sum(int(x) for x in np.diag(power@a@power@b@reverse))
        right = sum(int(x) for x in np.diag(power@b@a))
        if left != (core.completed().denominator*inverse_denominator)**length*right:
            raise ArithmeticError('a source-derived Gibbs KMS boundary identity failed')
        checked.append(name)
    return {'beta':length, 'tested_insertions':checked,
            'exact_two_sided_transfer_inverse_checked':True,
            'trace_cyclicity_supplies_the_full_matrix_algebra_boundary':True,
            'ordinary_Gibbs_density':'T^beta/Tr(T^beta)',
            'native_coordinate_star':'X -> D X^* D^-1',
            'Heisenberg_convention':'alpha_t(X)=exp(itH) X exp(-itH), H=-log(T)'}


def ground_limit_bound(length):
    _integer(length, 0, 10000, 'loop length')
    # Reuse the original native exclusion gap; D<=19/30 on the same graph.
    if TRANSFER_RATE_UPPER != 1-(1-burn_in.R_UPPER)/F(19, 30):
        raise ArithmeticError('the incident-transfer gap comparison changed')
    excess = 65*TRANSFER_RATE_UPPER**length
    return excess/(1+excess)


def cap_bound(length, cap):
    _integer(length, 1, 8, 'loop length')
    _integer(cap, 1, 1024, 'attempt cap')
    timeout = 1-(1-F(7, 10)**cap)**length
    acceptance = thermal_cycle(length)['return_probability']
    return {'complete_output_error': timeout, 'return_probability': acceptance,
            'capped_return_lower': max(F(), acceptance-timeout),
            'conditional_return_error': min(F(1), timeout/acceptance)}


def reflection_certificate():
    """Exact weighted Gram checks for four native half-history vectors."""
    data, kernel = core.source(), core.completed()
    d = sp.eye(DIMENSION)-sp.Matrix(data.remote_counts)/data.denominator
    inverse = d.inv(method='DM')
    k = sp.Matrix(kernel.counts)/kernel.denominator
    f0 = sp.diag(*(state[0] for state in data.states))
    f1 = sp.diag(*(state[1] for state in data.states))
    one = sp.ones(DIMENSION, 1)
    vectors = sp.Matrix.hstack(one, f0*one, k*f0*one, f1*k*f0*one)
    if inverse*k != k.T*inverse or d*f0 != f0*d or d*f1 != f1*d:
        raise ArithmeticError('the positive weighted native field dictionary changed')
    grams = []
    for offset in (0, 1):
        matrix = data.tick_probability/DIMENSION*vectors.T*inverse*(k**offset)*vectors
        if matrix != matrix.T:
            raise ArithmeticError('reflection Gram is not symmetric')
        minors = []
        for size in range(1, 5):
            for subset in combinations(range(4), size):
                minor = matrix.extract(subset, subset).det(method='domain-ge')
                if minor < 0:
                    raise ArithmeticError('a native reflected half-history is negative')
                minors.append(_fraction(minor))
        grams.append({'offset': offset, 'matrix': tuple(tuple(_fraction(x) for x in row) for row in matrix.tolist()),
                      'principal_minor_count': len(minors), 'strictly_positive_minor_count': sum(x > 0 for x in minors)})
    return {'half_history_vectors': ['1', 'F0 1', 'K F0 1', 'F1 K F0 1'],
            'metric': 'p D^-1 /66', 'grams': grams,
            'ground_vector_fixed': k*one == one,
            'all_local_words_have_positive_reflection_grams': True,
            'native_full_configuration_kernel_is_classically_reversible': k == k.T}


def local_closure_control():
    """Returning local loads is weaker than returning the full source state."""
    law = loop_law()
    states = core.source().states
    full, local = {}, {}
    for (_, word), matrix in zip(law.records, law.counts, strict=True):
        full[word] = F(sum(matrix[i][i] for i in range(DIMENSION)), law.denominator)
        local[word] = F(sum(matrix[i][j] for i in range(DIMENSION) for j in range(DIMENSION)
                            if states[i][:2] == states[j][:2]), law.denominator)
    p, q = sum(full.values(), F()), sum(local.values(), F())
    distance = sum((abs(full[word]/p-local[word]/q) for word in full), F())/2
    if not q > p or distance <= 0:
        raise ArithmeticError('local-return versus full-return control disappeared')
    return {'full_return_probability': p, 'local_return_probability': q,
            'conditional_field_distance': distance,
            'primitive_local_return_does_not_supply_the_Gibbs_closure': True}


def _upper(value, denominator=10**12):
    return F(-(-value.numerator*denominator//value.denominator), denominator)


@lru_cache(maxsize=1)
def class_source():
    kernel = native_clock.transfer.native_counts()
    reduction = native_clock.lumping(kernel)
    expected = ((F(14,15), F(1,15), F()), (F(1,15), F(9,10), F(1,30)),
                (F(), F(1,6), F(5,6)))
    if reduction['matrix'] != expected or reduction['sizes'] != (30,30,6):
        raise ArithmeticError('the active modular probe and class closure source differ')
    return kernel, reduction


def class_loop_law(length=2, *, herald=True):
    """Complete coarse class records, including initial rejection and nonreturn."""
    _integer(length, 1, 8, 'class loop length')
    if type(herald) is not bool:
        raise ValueError('an original Boolean initial filter switch required')
    _, reduction = class_source()
    prior = (F(1,11),)*3 if herald else reduction['stationary']
    active = {(i,i,()): prior[i] for i in range(3)}
    for _ in range(length):
        following = defaultdict(F)
        for (origin, state, word), weight in active.items():
            for final, probability in enumerate(reduction['matrix'][state]):
                following[origin,final,(*word,final)] += weight*probability
        active = following
    output = defaultdict(F)
    output['initial_reject',()] = F(8,11) if herald else F()
    for (origin, final, word), weight in active.items():
        output['return' if origin==final else 'different',word] += weight
    # Every complete word is retained, even when its probability is zero.
    from itertools import product
    for word in product(range(3),repeat=length):
        for status in ('return','different'):output[status,word] += F()
    if sum(output.values(),F())!=1:
        raise ArithmeticError('class closure omitted a rejection or failed return')
    return dict(sorted(output.items()))


def class_native_oracle(length=2, *, selection='first'):
    """Native endpoints after a deterministic six-states-per-class herald."""
    _integer(length, 1, 4, 'native class oracle length')
    if selection not in ('first','last'):
        raise ValueError('first or last six states per class required')
    kernel, reduction = class_source()
    classes = reduction['classes']
    chosen = []
    for label in range(3):
        members = sorted((i for i,c in enumerate(classes) if c==label),key=lambda i:kernel.states[i])
        chosen.extend(members[:6] if selection=='first' else members[-6:])
    index = {state:i for i,state in enumerate(kernel.states)}
    direct = []
    for state in kernel.states:
        direct.append(tuple(index[native.native(state,seam,coin)[0]]
                            for seam in kernel.seams for coin in (0,1)))
    active = {(classes[i],i,()):1 for i in chosen}
    for _ in range(length):
        following = defaultdict(int)
        for (origin,state,word),count in active.items():
            for final in direct[state]:
                following[origin,final,(*word,classes[final])] += count
        active = following
    denominator = 66*60**length
    output = defaultdict(F)
    output['initial_reject',()] = F(8,11)
    for (origin,final,word),count in active.items():
        output['return' if origin==classes[final] else 'different',word] += F(count,denominator)
    from itertools import product
    for word in product(range(3),repeat=length):
        for status in ('return','different'):output[status,word] += F()
    return dict(sorted(output.items()))


def class_cycle(length=2):
    _integer(length, 1, 4, 'class cycle length')
    _, reduction = class_source()
    p = sp.Matrix(reduction['matrix'])
    power = p**length
    partition = _fraction(sp.trace(power))
    law = class_loop_law(length)
    acceptance = sum((v for (status,_),v in law.items() if status=='return'),F())
    if acceptance != partition/11:
        raise ArithmeticError('the heralded native class loop lost its Gibbs trace')
    for selection in ('first','last'):
        if law != class_native_oracle(length,selection=selection):
            raise ArithmeticError('class strong lumping did not survive the initial filter')
    missing = class_loop_law(length,herald=False)
    other = sum((v for (status,_),v in missing.items() if status=='return'),F())
    difference = sum((abs(v/acceptance-missing[key]/other)
                      for key,v in law.items() if key[0]=='return'),F())/2
    if difference <= 0:
        raise ArithmeticError('dropping the explicit uniform-class herald no longer changes the Gibbs words')
    return {'native_attempt_length':length, 'inverse_temperature_beta':length,
            'initial_class_sizes':reduction['sizes'],'selected_native_states_per_class':6,
            'initial_filter_success':F(3,11),'initial_filter_failure':F(8,11),
            'partition':partition,'combined_return_probability':acceptance,
            'combined_different_probability':F(3,11)-acceptance,
            'conditional_Gibbs_class_diagonal':tuple(_fraction(power[i,i])/partition for i in range(3)),
            'complete_class_word_law':[{'status':s,'word':w,'probability':v} for (s,w),v in law.items()],
            'first_and_last_six_per_class_native_oracles_agree':True,
            'missing_herald_conditional_distance':difference,
            'reference':'sigma_L=S3^L/Tr(S3^L), S3=diag(pi)^1/2 P diag(pi)^-1/2',
            'active_probe_state_modular_parameter_for_U_k':'s=k/L with U_k=exp(+ikH3)',
            'physical_coherent_code_Gibbs_preparation_claimed':False}


def _summary(law):
    probabilities = law.probabilities()
    totals = {status: sum((p for (s, _), p in probabilities.items() if s == status), F())
              for status in ('return', 'different', 'timeout')}
    return {'ports': law.ports, 'attempt_cap': law.cap, 'feedback': law.feedback,
            'full_endpoint_block_count': len(law.records), 'endpoint_table_entries': len(law.records)*DIMENSION**2,
            'counts_sha256': _digest(law.counts), 'denominator': str(law.denominator),
            'status_probabilities': totals,
            'complete_word_probabilities': [{'status': status, 'word': word, 'probability': probability}
                                            for (status, word), probability in sorted(probabilities.items())]}


def _value(value):
    if type(value) is F:
        return str(value)
    if type(value) in (list, tuple):
        return [_value(x) for x in value]
    if type(value) is dict:
        return {k: _value(v) for k, v in value.items()}
    if value is None or type(value) in (str, int, bool):
        return value
    raise TypeError(f'unsupported receipt value {type(value).__name__}')


def run():
    verified = []
    for cap in (1, 2):
        for feedback in (False, True):
            law = loop_law(cap=cap, feedback=feedback)
            oracle = native_oracle(cap=cap, feedback=feedback)
            if law != oracle:
                raise ArithmeticError('independent native closure oracle disagrees')
            verified.append(_summary(law))
    exact = loop_law()
    changed = loop_law(feedback=True)
    p, q = exact.probabilities(), changed.probabilities()
    feedback_distance = sum((abs(p[key]-q[key]) for key in p), F())/2
    if feedback_distance <= 0:
        raise ArithmeticError('the equilibrium/intervention control disappeared')
    cap = loop_law(cap=2)
    capped, complete = cap.probabilities(), exact.probabilities()
    cap_tv = sum((abs(capped.get(key,F())-complete.get(key,F())) for key in set(capped)|set(complete)),F())/2
    if cap_tv != cap_bound(2,2)['complete_output_error']:
        raise ArithmeticError('complete capped history error did not retain timeout')
    budget = next(n for n in range(1000) if ground_limit_bound(n) <= F(1,1000))
    return _value({'schema':'oph.observer_time_fields.thermal.v1',
        'source_sha256':{path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
        'source':{'region':core.source().region,'dimension':66,'tick_probability':F(3,10),
                  'effective_transfer_lower':core.spectral_certificate()['effective_transfer_lower']},
        'open_ground_state':{'centered_field':'n0-1/6','moments':open_moments(),
             'finite_beta_KMS_defect_at_beta_one':kms_defect_lower(1),
             'all_positive_beta_defect_lower':'min(beta,1)*25/594',
             'finite_beta_KMS_on_the_dynamical_field_algebra':False},
        'reflection':reflection_certificate(),
        'closed_cycles':[thermal_cycle(n) for n in (1,2,3,4)],
        'finite_KMS_boundary':boundary_certificate(),
        'same_three_class_modular_reference':[class_cycle(length) for length in (1,2,3)],
        'uncapped_complete_output':_summary(exact), 'native_oracle_checks':verified,
        'cap_bounds':[{'cap':h,**cap_bound(2,h)} for h in (2,16,32,64)],
        'feedback_complete_output':_summary(changed), 'feedback_complete_distance':feedback_distance,
        'local_closure_control':local_closure_control(),
        'ground_limit':{'nonconstant_transfer_upper':TRANSFER_RATE_UPPER,
            'length_for_error_at_most_1_over_1000':budget,'error_upper':_upper(ground_limit_bound(budget))},
        'scope':{'full_initial_final_configuration_closure_service_required':True,
            'local_bit_and_current_readback_required':True, 'uniform_sector_preparation_declared':True,
            'closure_failure_and_wait_timeouts_retained':True,
            'conditional_closed_field_words_equal_source_derived_Gibbs_correlators':True,
            'faithful_Gibbs_reference_reconstructs_the_declared_closed_cylinders':True,
            'local_records_uniquely_identify_the_full_Gibbs_density':False,
            'native_quantum_M66_density_identified_with_transfer_Gibbs_state':False,
            'unconditioned_stationary_native_observer_is_finite_temperature_KMS':False,
            'feedback_histories_assumed_unforced_equilibrium':False,
            'native_incident_waiting_identified_with_modular_unitary_flow':False,
            'physical_temperature_seconds_or_Lorentzian_geometry_identified':False}})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--write',type=Path);mode.add_argument('--check',type=Path)
    args=parser.parse_args();result=run()
    if args.check:
        from model.observer.scalar_response.experiment import compare_receipt
        compare_receipt(strict_json.load(args.check),result)
        print('Native open-state obstruction, closed Gibbs words, closure failures and cap errors replayed.')
    elif args.write:
        with args.write.open('x',encoding='utf-8',newline='\n') as stream:
            json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
    else:
        print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':main()
