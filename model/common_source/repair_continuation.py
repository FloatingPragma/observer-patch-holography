"""Finite coherent continuation of the fixed recorded native repair.

The two-hop diagonal phase is an explicitly supplied coherent service. Its
effective logical generator is calculated from the existing code and repair;
neither the phase service nor its scheduling is selected by C/E repair alone.

This program explicitly CONTINUES after each one-proposal timeout, retaining
it as an idle record. The previous globally stopped history is a separate
control. Every interval applies exp(-i delta H2), then a complete proposal,
then its existing record-selected recovery. All raw records remain available.
We bound the work/reference marginal; no full classical-record transcript is
silently identified with a product environment or with a unitary channel.

The actual512 completion, not just its code columns or an alternative66
Householder completion, preserves the36-dimensional matching-support space.
For each matching its six-row Kraus blocks G_s are obtained from the COMPLETE
existing Gray mixers U6*U_count.T and ordered native input subsets. All rows
outside that space vanish algebraically. The cross transfer
B=sum_h sqrt(p_h) A_h obeys B*V=V and B.H*V=V. Indeed the complete
Kraus stack W=sum_h |h> tensor A_h is an isometry and the code's record
vector e=sum_h sqrt(p_h)|h> is normalized. Thus B=(<e| tensor I)W has norm
at most1; norm saturation on every code vector gives B.H*V=V as well.
On the reachable complement,
M=Q*B*Q is a contraction, and directed exact radical arithmetic establishes
max_j ||M_j^4||F<3/4 for C and max_j ||M_j^2||F<2/3 for E. Consequently
sum_(n>=1)||M^n|| <=15 (C), <=5 (E). No syndrome reset or postselection occurs.

Let J(rho)=V*rho*V.H, h=V.H*H2*V, K=Q*H2*V and
Y=sum_(n>=1)M^n*K. Set D(rho)=-i(Y*rho*V.H-V*rho*Y.H).
Then (I-Phi)D=Phi*(-i[H2,J])-J*(-i[h,.]). For F=Phi*Ad(exp(-i delta H2))
and G=Ad(exp(-i delta h)), the dressed embedding E=J+delta*D cancels the
first-order defect. CPTP contraction and integral exponential remainders give
||F^N J-J G^N||_diamond <= (2*t*b+t^2*c)/N,
c=(d^2+d0^2)/2+b*(d+d0), b=2*sum_bound*kappa,
where d=diam(H2 restricted to the reachable36 rows)=17/512,
d0=diam(h)=15/3072 and ||K||<kappa=11/1024. The full512 padding extension
has a different diameter and is not used in this invariant-space estimate.
This trace-norm bound includes arbitrary untouched reference systems. E is
only a proof device; it is neither prepared nor assumed positive.

Nearest theorem: Burgarth et al., Quantum4,289(2020), Theorem1, m=1,
https://arxiv.org/html/1809.09570v3#S4. The mapping is a finite CPTP kick and
the bounded Hamiltonian GKLS generator -i[H2,.]. Our added result is its
actual-source reachable algebra and explicit finite-N code/reference bound;
we do not infer its hypotheses from an observed endpoint identity.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import sympy as s

from . import parameter_selection as selection, record_code as code
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler
from ..calibrated_measurement.bounds import root_interval
from tools import strict_json


ROOT = Path(__file__).resolve().parents[2]
NUMERIC_ADMISSION = 1e-8
SOURCE_REPLAY = 2e-9
MAX_STEPS = 1 << 24


def _law(law):
    if type(law) is not str or law not in ("C", "E"):
        raise ValueError("the exact C or E instrument name required")
    return law


@lru_cache(maxsize=6)
def _mixer_exact(count):
    """The complete existing Gray state-preparation unitary, not its first column."""
    if type(count) is not int or not 1 <= count <= 6:
        raise ValueError("one through six surviving matching entries required")
    gray = tuple(j ^ (j >> 1) for j in range(8))
    residual = [1 if gray[j] < count else 0 for j in range(8)]
    rotations = []
    for row in range(7, 0, -1):
        x, y = residual[row-1], residual[row]
        norm = s.sqrt(x*x+y*y)
        cosine, sine = (x/norm, y/norm) if norm else (s.S.One, s.S.Zero)
        residual[row-1], residual[row] = norm, s.S.Zero
        rotations.append((gray[row-1], gray[row], cosine, sine))
    result = s.eye(8)
    for a, b, cosine, sine in reversed(rotations):
        old_a, old_b = result[a, :], result[b, :]
        result[a, :], result[b, :] = cosine*old_a-sine*old_b, sine*old_a+cosine*old_b
    result = result.applyfunc(s.simplify)
    if result.T*result != s.eye(8) or result[:, 0] != s.Matrix([1/s.sqrt(count)]*count+[0]*(8-count)):
        raise ArithmeticError("complete existing mixer identities failed")
    return s.ImmutableMatrix(result)


def _mix_columns(columns, matrix):
    return np.einsum("ab,hblc->halc", matrix,
                     columns.reshape(8, 8, 8, -1)).reshape(512, -1)


def _source_raw_maps(law, seam, coin):
    """Exact integer action of the existing source dilation on all512 addresses."""
    _law(law)
    route = np.array(code._route_permutation(seam))
    inverse = np.argsort(route)
    result = {raw: [] for raw in range(1 << circuit.RECORD_BITS[law])}
    for before in range(512):
        local = int(route[before] % 8)
        for a, b in circuit.REPAIR_SWAPS["C", coin]:
            if local == a:
                local = b
            elif local == b:
                local = a
        if law == "E" and (local >> coin) & 1:
            local ^= 32
        raw = local >> 3
        after = int(inverse[(local % 8)+8*(route[before]//8)])
        result[raw].append((before, after))
    return result


def _source_recover(columns, law, seam, coin, raw):
    count = len(code.branch_images(law, seam, coin, raw)[0])
    if not count or raw not in circuit.record_values(law, coin):
        return columns.copy()
    branch = np.array(code._branch_permutation(law, seam, coin, raw))
    encoded = np.array(code._encoder_permutation())
    value = columns[branch]
    value = _mix_columns(value, np.asarray(_mixer_exact(count), float).T)
    value = _mix_columns(value, np.asarray(_mixer_exact(6), float))
    output = np.empty_like(value)
    output[encoded] = value
    return output


def source_cross_transfer(law):
    """Overlap transfer from the ACTUAL existing512 unitary completion."""
    _law(law)
    weights, timeout = selection.completion_weights()
    result = float(timeout)*np.eye(512)
    for seam, weight in weights.items():
        for coin in (0, 1):
            for raw, mapping in _source_raw_maps(law, seam, coin).items():
                probability = code.branch_weight(law, seam, coin, raw)
                if not probability:
                    continue
                incoming = np.zeros((512, 512))
                for before, after in mapping:
                    incoming[after, before] = 1/math.sqrt(2)
                result += float(weight)*math.sqrt(float(probability))*_source_recover(
                    incoming, law, seam, coin, raw)
    return result


@lru_cache(maxsize=32)
def _source_block_exact(mask):
    if type(mask) is not int or not 0 < mask < 64:
        raise ValueError("a nonempty six-entry native subset mask required")
    selected = tuple(j for j in range(6) if mask >> j & 1)
    count = len(selected)
    change = _mixer_exact(6)*_mixer_exact(count).T
    if change[6:8, :count] != s.zeros(2, count):
        raise ArithmeticError("actual512 completion escapes the matching-support space")
    result = s.zeros(6)
    for rank, column in enumerate(selected):
        result[:, column] = change[:6, rank]
    effect = s.diag(*(int(j in selected) for j in range(6)))
    target = s.ones(6, 1)/s.sqrt(6)
    if ((result.T*result-effect).applyfunc(s.simplify) != s.zeros(6)
            or (result*target-s.sqrt(s.Rational(count, 6))*target).applyfunc(s.simplify) != s.zeros(6, 1)):
        raise ArithmeticError("actual complete-mixer Kraus identity failed")
    return s.ImmutableMatrix(result)


@lru_cache(maxsize=2)
def _model(law):
    _law(law)
    geometry = selection.code_geometry()
    weights, timeout = selection.completion_weights()
    binary_addresses = {circuit.encode(state) for state in geometry['states']}
    rows, maps, effects, cross = [], [], [s.zeros(6) for _ in range(6)], [s.Rational(timeout)*s.eye(6) for _ in range(6)]
    full_integer_maps_checked = 0
    for seam, weight in weights.items():
        native_maps = selection._raw_maps(law, seam, geometry)
        for coin in (0, 1):
            source_maps = _source_raw_maps(law, seam, coin)
            for raw, mapping in source_maps.items():
                expected = native_maps[coin, raw]
                actual = {geometry['index'][circuit.decode(before)]: geometry['index'][circuit.decode(after)]
                          for before, after in mapping if before in binary_addresses}
                if actual != expected:
                    raise ArithmeticError("512-source native dilation disagrees on binary66 inputs")
                full_integer_maps_checked += len(mapping)
                subsets = tuple(tuple(j for j, state in enumerate(support) if state in expected)
                                for support in geometry['supports'])
                if len(set(map(len, subsets))) != 1:
                    raise ArithmeticError("native record reveals a logical label")
                count = len(subsets[0])
                probability = weight*F(count, 12)
                rows.append({'seam': list(seam), 'coin': coin, 'raw': raw,
                             'code_probability_exact': str(probability), 'zero': not count,
                             'raw_abort': raw not in circuit.record_values(law, coin)})
                if not count:
                    continue
                blocks = []
                for logical, subset in enumerate(subsets):
                    mask = sum(1 << j for j in subset)
                    block = _source_block_exact(mask)
                    # Bind the ordered subset chart to the actual completed
                    # branch permutation, not just its uniform code vector.
                    inverse_chart = np.argsort(np.array(code._branch_permutation(law, seam, coin, raw)))
                    for rank, position in enumerate(subset):
                        before = geometry['supports'][logical][position]
                        address = circuit.encode(geometry['states'][expected[before]])
                        if inverse_chart[address] != logical+8*rank:
                            raise ArithmeticError("actual ordered syndrome chart differs")
                    effects[logical] += s.Rational(weight, 2)*s.diag(*(int(j in subset) for j in range(6)))
                    cross[logical] += s.Rational(weight)*s.sqrt(s.Rational(count, 24))*block
                    blocks.append(np.asarray(block, float)*math.sqrt(float(weight)/2))
                maps.append(np.array(blocks))
    for effect in effects:
        if effect+s.Rational(timeout)*s.eye(6) != s.eye(6):
            raise ArithmeticError("complete corrected channel is not trace preserving")
    maps.append(np.array([math.sqrt(float(timeout))*np.eye(6)]*6))
    data = np.array(maps)
    data.flags.writeable = False
    q = s.eye(6)-s.ones(6)/6
    cross = tuple(s.ImmutableMatrix((q*block*q).applyfunc(s.simplify)) for block in cross)
    return {'maps': data, 'cross': cross, 'records': rows,
            'source_integer_input_checks': full_integer_maps_checked,
            'timeout': timeout, 'weights': weights}


def _plus(a, b):
    return a[0]+b[0], a[1]+b[1]


def _times(a, b):
    corners = tuple(x*y for x in a for y in b)
    return min(corners), max(corners)


@lru_cache(maxsize=1024)
def _radical(expression):
    value = s.sympify(expression)
    if value.is_Rational:
        point = F(int(value.p), int(value.q))
        return point, point
    if value.is_Add:
        result = F(), F()
        for term in value.args:
            result = _plus(result, _radical(term))
        return result
    if value.is_Mul:
        result = F(1), F(1)
        for term in value.args:
            result = _times(result, _radical(term))
        return result
    if value.is_Pow and value.exp == s.Rational(1, 2):
        lower, upper = _radical(value.base)
        if lower < 0:
            raise ValueError("negative radical base")
        return root_interval(lower)[0], root_interval(upper)[1]
    raise ValueError("unsupported exact mixer radical")


def _interval_product(left, right):
    return [[sum_intervals(_times(left[i][k], right[k][j]) for k in range(6))
             for j in range(6)] for i in range(6)]


def sum_intervals(values):
    result = F(), F()
    for value in values:
        result = _plus(result, value)
    return result


@lru_cache(maxsize=2)
def _contraction_data(law):
    """Cache immutable exact admissions and diagnostic values only."""
    model = _model(law)
    power, upper, geometric_sum = (4, F(3, 4), 15) if law == 'C' else (2, F(2, 3), 5)
    rows = []
    for block in model['cross']:
        interval = [[_radical(block[i, j]) for j in range(6)] for i in range(6)]
        product = [[(F(int(i == j)), F(int(i == j))) for j in range(6)] for i in range(6)]
        for _ in range(power):
            product = _interval_product(product, interval)
        frobenius_squared_upper = sum((max(abs(a), abs(b))**2 for row in product for a, b in row), F())
        if frobenius_squared_upper >= upper**2:
            raise ArithmeticError("actual512 reachable cross-transfer contraction is unproved")
        rows.append(float(root_interval(frobenius_squared_upper)[1]))
    return power, str(upper), tuple(rows), str(geometric_sum)


def contraction_certificate(law):
    """Fresh public receipt; caller mutation cannot alter a certified bound."""
    _law(law)
    power, upper, rows, geometric_sum = _contraction_data(law)
    return {'power': power, 'strict_Frobenius_upper_exact': upper,
            'directed_Frobenius_upper_diagnostics': list(rows),
            'sum_norm_powers_from_one_upper_exact': geometric_sum,
            'proof': 'CPTP/code identity gives ||M||<=1; blockwise ||M^m||F<r gives sum(n>=1)||M^n||<=m/(1-r)-1',
            'outward_rational_arithmetic': True}


def logical_hamiltonian():
    return selection.diagonal_cost(2)['compression']


def _time_steps(time, steps):
    if (type(time) is not F or type(time.numerator) is not int or type(time.denominator) is not int
            or time.denominator <= 0 or not 0 <= time <= 32768):
        raise ValueError("an original Fraction coherent parameter in[0,32768] required")
    if type(steps) is not int or not 1 <= steps <= MAX_STEPS:
        raise ValueError("one through2^24 declared intervals required")
    return time, steps


def error_bound(law, time, steps):
    """Completely bounded TRACE norm; half trace distance is at most half this."""
    _law(law)
    time, steps = _time_steps(time, steps)
    certificate = contraction_certificate(law)
    energies = selection.diagonal_cost(2)
    geometry = selection.code_geometry()
    reachable = [energies['values'][i] for row in geometry['supports'] for i in row]
    d = F(max(reachable)-min(reachable))
    means = tuple(energies['compression'].diagonal())
    d0, kappa = F(max(means)-min(means)), F(11, 1024)
    if max(energies['variance']) >= s.Rational(kappa*kappa):
        raise ArithmeticError("exact code-to-syndrome coupling upper bound failed")
    b = 2*F(certificate['sum_norm_powers_from_one_upper_exact'])*kappa
    constant = (d*d+d0*d0)/2+b*(d+d0)
    return {'trace_norm_upper_exact': str(min(F(2), (2*time*b+time*time*constant)/steps)),
            'untruncated_trace_norm_upper_exact': str((2*time*b+time*time*constant)/steps),
            'd_reachable_exact': str(d), 'd0_exact': str(d0), 'kappa_strict_upper_exact': str(kappa),
            'dressing_map_cb_norm_upper_exact': str(b), 'quadratic_constant_exact': str(constant),
            'arbitrary_untouched_reference_included': True,
            'scope': 'ideal-source reachable carrier/reference marginal; full native records remain retained'}


@lru_cache(maxsize=72)
def _pair_transfer(law, left, right):
    result = sum((np.kron(row[left], row[right].conj()) for row in _model(law)['maps']), np.zeros((36, 36)))
    result.flags.writeable = False
    return result


def continuation_blocks(law, time, steps, *, absorbing=False):
    """All36 input operator images on the closed36 matching-support rows.

    Axes are output matching,left local index,output matching,right local
    index, stored as(block_left,block_right,local_left,local_right). Input
    |j><k| occupies only output block(j,k). No code projection is performed.
    """
    _law(law)
    time, steps = _time_steps(time, steps)
    if type(absorbing) is not bool:
        raise ValueError("explicit Boolean terminal policy required")
    geometry = selection.code_geometry()
    energy = selection.diagonal_cost(2)['values']
    values = np.array([[float(energy[i]) for i in row] for row in geometry['supports']])
    delta = float(time/steps)
    result = np.empty((6, 6, 6, 6), complex)
    timeout = float(_model(law)['timeout'])
    for left in range(6):
        for right in range(left, 6):
            phase = np.exp(-1j*delta*(values[left, :, None]-values[right, None, :])).ravel()
            base = _pair_transfer(law, left, right)
            if absorbing:
                live = (base-timeout*np.eye(36))*phase
                full = np.block([[live, np.zeros((36, 36))],
                                 [timeout*np.diag(phase), np.eye(36)]])
                output = np.linalg.matrix_power(full, steps)@np.r_[np.ones(36)/6, np.zeros(36)]
                vector = output[:36]+output[36:]
            else:
                vector = np.linalg.matrix_power(base*phase, steps)@(np.ones(36)/6)
            result[left, right] = vector.reshape(6, 6)
            result[right, left] = result[left, right].conj().T
    return result


def _choi(blocks):
    # An isometric chart of the full216-dimensional Choi state: output's
    # matching label always equals the untouched reference's input label.
    return blocks.transpose(0, 2, 1, 3).reshape(36, 36)/6


def _trace_norm(matrix):
    return float(np.sum(abs(np.linalg.eigvalsh(matrix))))


def _ideal_choi(time):
    means = np.asarray(list(logical_hamiltonian().diagonal()), float).ravel()
    vector = np.repeat(np.exp(-1j*float(time)*means), 6)/6
    return np.outer(vector, vector.conj())


def ramsey(blocks):
    rows = []
    for spectator in range(3):
        a, b = 2*spectator, 2*spectator+1
        visibility = np.sum(blocks[a, b])/6
        population = float((np.sum(blocks[a, a])+np.sum(blocks[b, b])).real/12)
        rows.append({'spectator': spectator, 'X_mean': float(visibility.real),
                     'Y_mean': float(-visibility.imag), 'code_pair_probability': population,
                     'outside_code_pair_probability': 1-population})
    return rows


def two_time_control(law, time=F(1024), steps=65536):
    """Two complete three-outcome X instruments; outside-code outcomes retained."""
    blocks = continuation_blocks(law, time, steps)
    transfers = {}
    geometry = selection.code_geometry()
    energies = selection.diagonal_cost(2)['values']
    values = np.array([[float(energies[i]) for i in row] for row in geometry['supports']])
    for a in (0, 1):
        for b in (0, 1):
            phase = np.exp(-1j*float(time/steps)*(values[a, :, None]-values[b, None, :])).ravel()
            transfers[a, b] = np.linalg.matrix_power(_pair_transfer(law, a, b)*phase, steps)
    def evolve(rho):
        output = np.zeros_like(rho)
        for a in (0, 1):
            for b in (0, 1):
                output[6*a:6*a+6, 6*b:6*b+6] = (transfers[a, b]@rho[6*a:6*a+6, 6*b:6*b+6].ravel()).reshape(6, 6)
        return output
    vectors = [np.r_[np.ones(6), sign*np.ones(6)]/math.sqrt(12) for sign in (1, -1)]
    projectors = [np.outer(v, v) for v in vectors]
    projectors.append(np.eye(12)-sum(projectors))
    before = np.block([[blocks[a, b]/2 for b in (0, 1)] for a in (0, 1)])
    probabilities = np.array([[float(np.trace(Q@evolve(P@before@P)).real) for Q in projectors]
                              for P in projectors])
    gap = float(logical_hamiltonian()[1, 1]-logical_hamiltonian()[0, 0])
    cosine = math.cos(gap*float(time))
    ideal = np.zeros((3, 3))
    for i, first in enumerate((1, -1)):
        for j, second in enumerate((1, -1)):
            ideal[i, j] = (1+first*cosine)*(1+first*second*cosine)/4
    if abs(probabilities.sum()-1) > NUMERIC_ADMISSION or probabilities.min() < -NUMERIC_ADMISSION:
        raise ArithmeticError("complete two-time probabilities lost positivity/normalization")
    return {'coherent_parameter_per_interval': str(time), 'steps_per_interval': steps,
            'outcomes': ['plus', 'minus', 'outside_code_pair'], 'joint_probabilities': probabilities.tolist(),
            'ideal_joint_probabilities': ideal.tolist(), 'L1_difference': float(np.sum(abs(probabilities-ideal))),
            'universal_L1_upper_exact': str(min(F(2), 2*F(error_bound(law, time, steps)['trace_norm_upper_exact']))),
            'no_outcome_postselection_or_probability_repair': True}


def implementation_bound(law, time, steps, epsilon_step):
    """Conditional service-error accumulation; no physical precision inferred."""
    if (type(epsilon_step) is not F or type(epsilon_step.numerator) is not int
            or type(epsilon_step.denominator) is not int or epsilon_step.denominator <= 0
            or not 0 <= epsilon_step <= 2):
        raise ValueError("an exact per-step diamond trace-norm allowance in[0,2] required")
    bound = error_bound(law, time, steps)
    ideal = F(bound['untruncated_trace_norm_upper_exact'])
    return {'ideal_target_bound_exact': str(ideal), 'assumed_per_step_diamond_error_exact': str(epsilon_step),
            'accumulated_service_error_exact': str(steps*epsilon_step),
            'combined_trace_norm_upper_exact': str(min(F(2), ideal+steps*epsilon_step)),
            'hypothesis': 'actual whole interval is CPTP on full512 including abort/service flags; full-output diamond error is uniformly bounded on reachable inputs and arbitrary references, including all padding leakage; initial preparation/readout errors are additional',
            'literal_source_diagnostics_or_physical_errors_certify_this_allowance': False}


def finite_rate_modes(law):
    """Wrapped response poles of this finite program, not physical masses."""
    _law(law)
    geometry = selection.code_geometry()
    energy = selection.diagonal_cost(2)['values']
    values = np.array([[float(energy[i]) for i in row] for row in geometry['supports']])
    rows = []
    for delta in (F(1, 4), F(1), F(4)):
        channels = []
        for spectator in range(3):
            a, b = 2*spectator, 2*spectator+1
            phase = np.exp(-1j*float(delta)*(values[a, :, None]-values[b, None, :])).ravel()
            transfer = _pair_transfer(law, a, b)*phase
            eigenvalues, vectors = np.linalg.eig(transfer)
            order = np.argsort(abs(eigenvalues-1))
            index = order[0]
            value = eigenvalues[index]
            separation = float(abs(eigenvalues[order[1]]-1)-abs(value-1))
            residual = float(np.linalg.norm(transfer@vectors[:, index]-value*vectors[:, index]))
            if not 0 < abs(value) < 1 or separation < 0.01 or residual > 1e-10:
                raise ArithmeticError("finite-rate pole is unresolved, noncontractive or aliased to another mode")
            channels.append({'spectator': spectator, 'pole_real': float(value.real), 'pole_imag': float(value.imag),
                'wrapped_phase_rate': math.atan2(value.imag, value.real)/float(delta),
                'attenuation_rate': -math.log(abs(value))/float(delta),
                'nearest_to_one_mode_separation': separation, 'eigenvector_residual': residual})
        rows.append({'coherent_parameter_per_interval': str(delta), 'channels': channels})
    return {'time_parameter': 'declared coherent phase per interval', 'controls': rows,
            'finite_rate_corrections_include_damping': True,
            'wrapped_phase_rates_are_unwrapped_physical_masses': False,
            'three_spectator_channels_are_three_independently_occupiable_species': False,
            'comparison': 'limit signed gaps(1/256,-1/1024,-1/3072); finite rates are dissipative corrections of this supplied H2 program'}


def absorbing_control(law, time=F(1024), steps=256):
    """The old terminal policy is not restarted or conditioned on survival."""
    _time_steps(time, steps)
    timeout = selection.completion_weights()[1]
    blocks = continuation_blocks(law, time, steps, absorbing=True)
    state = _choi(blocks)
    expected_intervals = (1-(1-timeout)**steps)/timeout
    return {'steps': steps, 'total_programmed_coherent_parameter': str(time),
            'survival_probability_exact': str((1-timeout)**steps),
            'mean_executed_intervals_exact': str(expected_intervals),
            'mean_accumulated_coherent_parameter_exact': str(time*expected_intervals/steps),
            'clock_parameter_upper_exact': str(time/(steps*timeout)),
            'trace_norm_to_unchanged_input_upper_exact': str(min(F(2), F(15, 256)*time/(steps*timeout))),
            'Choi_trace_norm_to_unchanged_input_diagnostic': _trace_norm(state-_ideal_choi(F())),
            'limit_at_fixed_total_coherent_parameter': 'unchanged input; mean executed coherent parameter<=t/(N*timeout)->0',
            'timeout_and_stopping_index_retained': True, 'timeout_restarted': False}


def _phase_word(time):
    values = tuple(F(value) for value in selection._probe_energies())
    coefficients = list(values)
    size = 1
    while size < 512:
        for first in range(0, 512, 2*size):
            for offset in range(size):
                a, b = coefficients[first+offset], coefficients[first+offset+size]
                coefficients[first+offset], coefficients[first+offset+size] = a+b, a-b
        size *= 2
    word = []
    for mask, coefficient in enumerate(coefficients):
        label = ['Z' if mask >> bit & 1 else 'I' for bit in range(9)]
        word.extend(compiler.pauli_rotation(label, float(-time*coefficient/512)))
    return word


def _source_targets():
    ideal = code.ideal_code_matrix()
    reachable = np.zeros((512, 36), complex)
    for j, state in enumerate(state for row in code.SUPPORTS for state in row):
        reachable[circuit.encode(state), j] = 1
    phase = np.exp(-1j*float(F(1, 64))*np.array(list(map(float, selection._probe_energies()))))
    targets = {'encoded': ideal, 'phase_outputs': phase[:, None]*np.column_stack((ideal, reachable))}
    for count in (1, 4, 5, 6):
        targets[f'mixer_{count}'] = np.asarray(_mixer_exact(count), complex)
    for law, seam in (('C', (2, 3)), ('E', (2, 10))):
        output = []
        for coin in (0, 1):
            for raw, mapping in _source_raw_maps(law, seam, coin).items():
                incoming = np.zeros((512, 42), complex)
                for before, after in mapping:
                    incoming[after] = targets['phase_outputs'][before]/math.sqrt(2)
                output.append(_source_recover(incoming, law, seam, coin, raw))
        targets['native_'+law+'_outputs'] = np.array(output)
    return targets


def verify_source_arrays(arrays):
    """Retained literal arrays independently meet the complete ideal targets."""
    targets = _source_targets()
    if type(arrays) is not dict or arrays.keys() != targets.keys():
        raise ValueError("complete source output menu changed")
    errors = {}
    for name, target in targets.items():
        value = arrays[name]
        if (np.ma.isMaskedArray(value) or not isinstance(value, np.ndarray)
                or value.dtype != np.dtype(complex) or value.shape != target.shape
                or not np.isfinite(value).all()):
            raise ValueError("source output schema changed: "+name)
        difference = (value-target).reshape(-1, target.shape[-1])
        error = float(np.linalg.norm(difference, 2))
        if error > NUMERIC_ADMISSION:
            raise ValueError("retained source operator exceeds admission: "+name)
        errors[name] = error
    return errors


def source_execution():
    arrays = {'encoded': np.array(code.code_matrix())}
    for count in (1, 4, 5, 6):
        target = np.array([1/math.sqrt(count)]*count+[0.]*(8-count))
        arrays[f'mixer_{count}'] = compiler.execute(compiler.compile_real_state(target), 3)
    word = _phase_word(F(1, 64))
    reachable = np.zeros((512, 36), complex)
    for j, state in enumerate(state for row in code.SUPPORTS for state in row):
        reachable[circuit.encode(state), j] = 1
    arrays['phase_outputs'] = compiler.execute(word, 9, np.column_stack((arrays['encoded'], reachable)))
    resources = {}
    differences = {}
    geometry = selection.code_geometry()
    binary_addresses = [circuit.encode(state) for state in geometry['states']]
    targets = _source_targets()
    for law, seam in (('C', (2, 3)), ('E', (2, 10))):
        branches = code.repair_branches(law, seam, arrays['phase_outputs'])
        arrays['native_'+law+'_outputs'] = np.stack([branch['amplitude'] for branch in branches])
        difference = 0.
        for branch, actual in zip(branches, targets['native_'+law+'_outputs'], strict=True):
            raw, coin = branch['record_code'], branch['coin']
            if not code.branch_weight(law, seam, coin, raw):
                continue
            native_map = selection._raw_maps(law, seam, geometry)[coin, raw]
            incoming = np.zeros((66, 42), complex)
            for before, after in native_map.items():
                incoming[after] = targets['phase_outputs'][binary_addresses[before]]/math.sqrt(2)
            alternative = code.declared_recovery_matrix66(law, seam, coin, raw)@incoming
            lifted = np.zeros((512, 42), complex)
            lifted[binary_addresses] = alternative
            difference = max(difference, float(np.linalg.norm(actual-lifted, 2)))
        differences[law] = difference
        resources[law] = {'seam': list(seam), 'all_raw_outputs': len(branches),
            'native_repair': code.repair_resources(law, seam),
            'maximum_recovery_pulses': max(branch['recovery_resources']['native_pulses'] for branch in branches)}
    if min(differences.values()) <= 1e-7:
        raise ArithmeticError("off-code completion comparison lost its nontrivial control")
    errors = verify_source_arrays(arrays)
    phase_cost = compiler.resources(word)
    maximum_repair = max(code.repair_resources(law, seam)['native_pulses']+
                         max(code.recovery_resources(law, seam, coin, raw)['native_pulses']
                             for coin in (0, 1) for raw in range(1 << circuit.RECORD_BITS[law]))
                         for law in ('C', 'E') for seam in _model(law)['weights'])
    # One existing selector proposal draws its five-bit address, of which
    # twenty accepted seams have multiplicity1 or2 and five addresses timeout.
    selector = selection.preparation.selectors(1)
    return {'probe_coherent_parameter': '1/64', 'all512_rows_and_six_reference_columns_retained': True,
            'additional_reachable_syndrome_input_columns_retained': 36,
            'total_literal_input_columns_per_native_branch': 42,
            'reachable_basis_columns_are_operator_queries_not_free_state_preparations': True,
            'all_eight_input_columns_of_complete_existing_mixers_executed': [1, 4, 5, 6],
            'retained_source_operator_errors': errors, 'native_executions': resources,
            'actual512_versus_alternative66_output_operator_difference': differences,
            'finite_time_theorem_uses_actual512_completion': True,
            'phase_word_resources': phase_cost, 'maximum_native_repair_plus_recovery_pulses': maximum_repair,
            'one_proposal_selector': {key: value for key, value in selector.items()
                                      if key in ('completion_timeout_probability', 'completion_accepted_seam_weights')},
            'per_interval_phase_pulses': phase_cost['native_pulses'],
            'per_interval_selector_auxiliary_native_attempts_upper': selector['completion_auxiliary_calls_per_slot_upper'],
            'fresh_selector_coin_and_raw_record_bits_per_interval_upper': {'C': 8, 'E': 9},
            'selector_source_history_executed_here': False,
            'N_interval_pulse_requirement': 'N*(phase_word_pulses+worst_native_repair_and_recovery_pulses+selector_preparation_and_readout_pulses)',
            'record_register_requirement': 'fresh selector address, coin and rawC/E record per interval; no archive overwrite; total storage grows linearly withN',
            'proper_load_export_service_requirement': 'two services per emitted pulse; physical initialization/storage/readout costs additional',
            'full_N_interval_pulse_or_protection_history_executed': False,
            'numerical_source_admission': NUMERIC_ADMISSION,
            'fresh_source_output_operator_replay_tolerance': SOURCE_REPLAY,
            'physical_service_error_certified': False}, arrays


def clock_winding_control():
    plus = np.array([1., 1.])/math.sqrt(2)
    values = []
    for time in (F(1), F(1, 2)):
        word = [compiler.pulse('I', -math.pi*float(time), [0]),
                compiler.pulse('Z', math.pi*float(time), [0])]
        source = compiler.execute(word, 1, plus)
        values.append({'coherent_parameter': str(time),
                       'X_mean': float(2*(source[0].conjugate()*source[1]).real),
                       'zero_generator_X_mean': 1.})
    return {'source_winding_Hamiltonian': '2*pi*|1><1|',
            'endpoint_unitary_equals_identity_at_coherent_parameter_one': True,
            'intermediate_quadrature_records_distinguish_it': values,
            'clock_parameter': 'supplied coherent program parameter, not native proposal index or laboratory seconds',
            'an_identity_endpoint_selects_a_unique_unwrapped_generator': False}


def _json(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json(child) for key, child in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(child) for child in value]
    return value


def source_bindings():
    paths = set(selection.source_bindings())|{
        'model/common_source/repair_continuation.py', 'model/common_source/test_repair_continuation.py',
        'model/calibrated_measurement/bounds.py'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def run():
    source_report, arrays = source_execution()
    laws = {}
    for law in ('C', 'E'):
        model = _model(law)
        rows = []
        for steps in (1, 16, 256, 4096, 65536):
            blocks = continuation_blocks(law, F(1024), steps)
            state = _choi(blocks)
            if (abs(np.trace(state)-1) > NUMERIC_ADMISSION
                    or np.min(np.linalg.eigvalsh(state)) < -NUMERIC_ADMISSION):
                raise ArithmeticError("finite continuing CPTP marginal failed its numerical admission")
            rows.append({'steps': steps, 'coherent_parameter': '1024',
                'Choi_trace_defect_diagnostic': float(abs(np.trace(state)-1)),
                'Choi_trace_norm_error_diagnostic': _trace_norm(state-_ideal_choi(F(1024))),
                'ideal_all_reference_bound': error_bound(law, F(1024), steps),
                'Ramsey_three_spectator_controls': ramsey(blocks)})
        laws[law] = {'accepted_native_seams': len(model['weights']),
            'complete_raw_record_menu': model['records'], 'timeout_idle_probability_exact': str(model['timeout']),
            'all512_integer_source_inputs_checked': model['source_integer_input_checks'],
            'actual_complete_source_mixer_reachable_rows': 36,
            'contraction': contraction_certificate(law), 'finite_time_controls': rows,
            'finite_rate_response_modes': finite_rate_modes(law),
            'two_time_record_control': two_time_control(law),
            'ideal_target_bound_at_2pow24_intervals': error_bound(law, F(1024), MAX_STEPS),
            'conditional_service_tradeoff_example': implementation_bound(law, F(1024), MAX_STEPS, F(1, 10**10)),
            'old_absorbing_terminal_control': absorbing_control(law)}
    report = {'schema': 'oph.repair_continuation.v1',
        'classification': 'controlled finite-source Zeno continuation of the existing corrected native cycle',
        'logical_hamiltonian_exact': selection._serialize(logical_hamiltonian()),
        'signed_one_CAR_three_spectator_gaps_exact': ['1/256', '-1/1024', '-1/3072'],
        'schedule': {'interval_order': ['two-hop coherent phase', 'one completion proposal', 'committed-record recovery'],
            'one_proposal_timeout': 'idle outcome retained; next interval explicitly allowed by this NEW continuing program',
            'earlier_global_terminal_program_is_not_reinterpreted': True,
            'N': 'declared number of intervals', 'time_parameter': 'supplied dimensionless coherent phase parametert; delta=t/N',
            'arbitrarily_fast_correction_or_phase_service_is_a_limit_assumption': True,
            'all_raw_records_and_timeout_addresses_retained': True},
        'laws': laws, 'source_execution': source_report, 'clock_winding_control': clock_winding_control(),
        'nearest_theorem': {'source': 'https://arxiv.org/html/1809.09570v3#S4', 'theorem': 'Theorem1, m=1',
            'mapping': 'finite-dimensional CPTP kickPhi; bounded Hamiltonian GKLS generator-i[H2,.]',
            'addition': 'actual512 reachable closure, directed contraction and explicit arbitrary-reference finite-N bound'},
        'scope': {'two_hop_phase_and_continuing_schedule_supplied': True,
            'logical_generator_derived_from_actual_complete_repair': True,
            'syndrome_postselection_or_reset_added': False,
            'empirical_masses_or_calibrated_packet_used': False,
            'one_CAR_mode_with_three_spectator_values_is_three_independent_species': False,
            '12_to_3_to_1_splittings_are_a_measured_hierarchy': False,
            'full_classical_record_transcript_converges_to_a_product_unitary_environment': False,
            'source_replay_numerics_or_finite_N_grid_are_outward_certificates': False,
            'directed_ideal_contraction_and_ideal_error_bound_are_certificates': True,
            'physical_clock_or_service_accuracy_or_interacting_mass_law_derived': False},
        'source_bindings': source_bindings()}
    return _json(report), arrays


def _compare(old, new, path='root'):
    if type(old) is not type(new):
        raise ValueError('receipt type changed: '+path)
    if isinstance(old, dict):
        if old.keys() != new.keys():
            raise ValueError('receipt keys changed: '+path)
        for key in old:
            _compare(old[key], new[key], path+'.'+key)
    elif isinstance(old, list):
        if len(old) != len(new):
            raise ValueError('receipt list changed: '+path)
        for index, (a, b) in enumerate(zip(old, new, strict=True)):
            _compare(a, b, path+f'[{index}]')
    elif type(old) is float:
        exact = path.endswith(('numerical_source_admission', 'fresh_source_output_operator_replay_tolerance'))
        if (not math.isfinite(old) or not math.isfinite(new)
                or (old != new if exact else not math.isclose(old, new, abs_tol=2e-8, rel_tol=2e-8))):
            raise ValueError('receipt numerical replay changed: '+path)
    elif old != new:
        raise ValueError('receipt exact replay changed: '+path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--output', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    report, arrays = run()
    path = args.output or args.check
    artifact = path.with_suffix('.npz')
    if args.output:
        if path.exists() or artifact.exists():
            parser.error('fresh receipt and array paths required')
        with artifact.open('xb') as stream:
            np.savez_compressed(stream, **arrays)
        report['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
        path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n')
    else:
        retained = strict_json.load(path)
        if retained.pop('artifact', None) != {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}:
            raise ValueError('literal source artifact custody changed')
        with np.load(artifact, allow_pickle=False) as saved:
            if len(saved.files) != len(set(saved.files)) or set(saved.files) != set(arrays):
                raise ValueError('literal source array menu changed')
            original = {key: saved[key] for key in saved.files}
            verify_source_arrays(original)
            for key, fresh in arrays.items():
                if np.linalg.norm((original[key]-fresh).reshape(-1, fresh.shape[-1]), 2) > SOURCE_REPLAY:
                    raise ValueError('literal source numerical replay changed: '+key)
        _compare(retained, report)
        print('actual-source coherent repair continuation, all-record maps and finite-N certificates: PASS')


if __name__ == '__main__':
    main()
