"""Native execution and matched comparison for one interacting band extension."""
from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import ceil, factorial, prod, sqrt
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.optimize import least_squares
from scipy.sparse.linalg import expm_multiply

from . import algebra as a
from ..observer.source_scalar import compiler
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..scalar_observations import source as inherited_source
from ..scalar_observations.protocol import native_records
from ..observer.native import native as native_repair
from ..observer.quantum_mixing import seams as native_seams
from ..relative_observer_geometry.temporal_records import incident_clock_bound


ROOT = Path(__file__).resolve().parents[2]
INTERVAL = F(1, 16)
STEPS = 64
SOURCE_PATHS = (
    'model/interacting_band_extension/algebra.py',
    'model/interacting_band_extension/construction.py',
    'model/observer/source_scalar/compiler.py',
    'model/scalar_observations/source.py',
    'model/scalar_observations/protocol.py',
    'model/scalar_observations/bounds.py',
    'model/relativistic_scalar/lattice.py',
    'model/observer/native.py',
    'model/observer/quantum_mixing.py',
    'model/common_source/slow_geometry.py',
    'model/common_source/records.py',
    'model/relative_observer_geometry/temporal_records.py',
    'evidence/native_geometric_source_20260925/spec.json',
    compiler.RESPONSE, compiler.CONTROL_RECEIPT)


def normalized(matrix):
    return np.array([[float(value)*sqrt(factorial(i)/factorial(j))
                      for j, value in enumerate(row)] for i, row in enumerate(matrix)])


@dataclass(frozen=True)
class Term:
    sites: tuple
    factors: tuple
    coefficient: F
    label: str

    def matrix(self):
        result = np.array([[1.]])
        for factor in reversed(self.factors):
            result = np.kron(result, normalized(factor))
        return float(self.coefficient)*result

    def norm_bound(self):
        return abs(self.coefficient)*prod(a.row_norm(matrix) for matrix in self.factors)


def terms(extended=True, matched=False, coupling=1):
    if type(extended) is not bool or type(matched) is not bool or (extended and matched):
        raise ValueError('admitted original or matched band selector required')
    g = a.rational(coupling)
    dims = a.DIMENSIONS if extended else a.LOW_DIMENSIONS
    result = []
    for mode, d in enumerate(dims):
        frequency_squared = F(1 + sum(x*x for x in a.WAVES[mode]))+(g/4 if matched else 0)
        quartic = g/(24 if mode == 0 else 16)
        p2, q2, q4 = (a.polynomial(d, 2, True), a.polynomial(d, 2), a.polynomial(d, 4))
        matrix = tuple(tuple(p2[i][j]/2+frequency_squared*q2[i][j]/2+quartic*q4[i][j]
                             for j in range(d)) for i in range(d))
        result.append(Term((mode,), (matrix,), F(1), f'onsite_{mode}'))
    coefficients = a.declared_quartic() if extended else a.quartic_coefficients(False)
    for powers, coefficient in sorted(coefficients.items()):
        selected = tuple(i for i, power in enumerate(powers) if power)
        if len(selected) == 1:
            continue  # Already included in its complete onsite generator.
        matrices = tuple(a.polynomial(dims[i], powers[i]) for i in selected)
        result.append(Term(selected, matrices, g*coefficient/24,
                           'quartic_'+'_'.join(f'{i}^{powers[i]}' for i in selected)))
    return tuple(result)


class Lattice:
    """Finite tensor code, with sparse independent full-Hamiltonian reference."""

    def __init__(self, extended=True, matched=False, coupling=1):
        self.extended, self.matched = extended, matched
        self.coupling = a.rational(coupling)
        self.dimensions = a.DIMENSIONS if extended else a.LOW_DIMENSIONS
        self.dimension = prod(self.dimensions)
        self.terms = terms(extended, matched, self.coupling)
        self.hamiltonian = sum((self.embed(term) for term in self.terms),
                               sparse.csr_matrix((self.dimension, self.dimension)))
        self.hamiltonian.eliminate_zeros()
        self.trace = float(self.hamiltonian.diagonal().sum())

    def embed(self, term):
        factors = dict(zip(term.sites, term.factors, strict=True))
        result = sparse.csr_matrix([[float(term.coefficient)]])
        for mode in reversed(range(len(self.dimensions))):
            factor = (sparse.csr_matrix(normalized(factors[mode])) if mode in factors
                      else sparse.eye(self.dimensions[mode], format='csr'))
            result = sparse.kron(result, factor, format='csr')
        return result

    def apply(self, state, sites, matrix):
        value = np.asarray(state, dtype=complex)
        if value.ndim not in (1, 2) or value.shape[0] != self.dimension:
            raise ValueError('one amplitude row per field configuration required')
        sites = tuple(sites)
        if not sites or len(set(sites)) != len(sites) or any(type(i) is not int or not 0 <= i < len(self.dimensions) for i in sites):
            raise ValueError('distinct admitted local tensor indices required')
        local_size = prod(self.dimensions[i] for i in sites)
        matrix = np.asarray(matrix)
        if matrix.shape != (local_size, local_size) or not np.isfinite(matrix).all():
            raise ValueError('a finite matrix matching the local code is required')
        tensor = value.reshape(self.dimensions+value.shape[1:], order='F')
        moved = np.moveaxis(tensor, sites, range(len(sites)))
        changed = (matrix@moved.reshape((local_size, -1), order='F')).reshape(moved.shape, order='F')
        return np.moveaxis(changed, range(len(sites)), sites).reshape(value.shape, order='F')

    def blank(self):
        result = np.zeros(self.dimension, complex)
        result[0] = 1
        return result

    def evolve(self, state, duration):
        t = float(a.rational(duration))
        return expm_multiply(-1j*t*self.hamiltonian, state, traceA=-1j*t*self.trace)

    def source_evolve(self, state, duration, steps):
        a.integer(steps)
        dt = a.rational(duration)/(2*steps)
        operations = [(term.sites, execute_term(term, dt)[0]) for term in self.terms]
        for _ in range(steps):
            for sites, matrix in operations+operations[::-1]:
                state = self.apply(state, sites, matrix)
        return state


@lru_cache(maxsize=256)
def execute_term(term, duration):
    """Contract the emitted local source pulses, never substituting a target."""
    duration = a.rational(duration)
    changes, eigenvalues, offset = [], [], 0
    for matrix in term.factors:
        values, basis = inherited_source._basis(normalized(matrix))
        bits = len(matrix).bit_length()-1
        changes.append(compiler.remap(compiler.compile_real_orthogonal(basis),
                                      list(range(offset, offset+bits))))
        eigenvalues.append(values)
        offset += bits
    values = np.array([float(term.coefficient)])
    for vector in reversed(eigenvalues):
        values = np.kron(values, vector)
    stages = [(f'inverse_factor_{i}', compiler.inverse(change)) for i, change in enumerate(changes)]
    stages.append(('product_phase', compiler.compile_diagonal(-float(duration)*values)))
    stages.extend((f'factor_{i}', change) for i, change in reversed(list(enumerate(changes))))
    target = expm(-1j*float(duration)*term.matrix())
    return inherited_source._execute(stages, offset, target)


def splitting_bound(lattice, duration, steps):
    """Uniform all-input operator bound for the local palindrome."""
    time = a.rational(duration)
    a.integer(steps)
    # Center each local term; the identical total scalar phase cancels in
    # the exact/product comparison. Norm bounds remain exact rationals.
    budget = F(0)
    for term in lattice.terms:
        if len(term.sites) == 1:
            matrix = term.factors[0]
            scalar = matrix[0][0]
            centered = tuple(tuple(value-(scalar if i == j else 0)
                                   for j, value in enumerate(row)) for i, row in enumerate(matrix))
            budget += abs(term.coefficient)*a.row_norm(centered)
        else:
            budget += term.norm_bound()
    error = budget**3*time**3/(3*steps**2)
    return {'centered_generator_norm_sum_upper': budget,
            'one_interval_operator_error_upper': min(F(2), error),
            'two_interval_complete_output_trace_distance_upper': min(F(1), 2*error)}


def history(lattice, sign=1, *, source=True, feedback=True, interval=INTERVAL,
            steps=STEPS, reset_high=False, return_states=False, initial=None):
    """Complete low two-bit law with the same disturbed field throughout."""
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError('preparation sign must be plus or minus one')
    for flag in (source, feedback, reset_high, return_states):
        if type(flag) is not bool:
            raise ValueError('original boolean protocol flags required')
    d = lattice.dimensions[1]
    q = normalized(a.polynomial(d, 1))
    if source:
        kick = inherited_source.feedback(d, sign/4, 2)[0]
        kraus = inherited_source.instrument(d, 1, np.pi/2, 2)[0]
        control = inherited_source.feedback(d, 1/4, 2)[0]
        evolve = lambda v: lattice.source_evolve(v, interval, steps)
    else:
        kick, control = expm(1j*sign*q/4), expm(1j*q/4)
        values, basis = np.linalg.eigh(q)
        kraus = tuple((basis*function((values+np.pi/2)/2))@basis.T
                      for function in (np.cos, np.sin))
        evolve = lambda v: lattice.evolve(v, interval)
    prepared = lattice.apply(lattice.blank() if initial is None else initial, (1,), kick)
    first_state = evolve(prepared)
    branches = [lattice.apply(first_state, (1,), matrix) for matrix in kraus]
    # The reset control is a genuine trace-and-reprepare channel: retain
    # every high occupation as a Kraus label, rather than selecting vacuum.
    if reset_high:
        if not lattice.extended:
            raise ValueError('high reset needs the extended field')
        low_size = prod(a.LOW_DIMENSIONS)
        reset = []
        for branch in branches:
            columns = branch.reshape((low_size, -1), order='F')
            lifted = np.zeros((lattice.dimension, columns.shape[1]), complex)
            lifted[:low_size, :] = columns
            reset.append(lifted)
        branches = reset
    if feedback:
        branches[1] = lattice.apply(branches[1], (1,), control)
    evolved = [evolve(branch) for branch in branches]
    rows, states = [], []
    for first, branch in enumerate(evolved):
        for last, matrix in enumerate(kraus):
            terminal = lattice.apply(branch, (1,), matrix)
            rows.append({'bits': [first, last], 'probability': float(np.linalg.norm(terminal)**2)})
            states.append(terminal)
    return (rows, states) if return_states else rows


class ProductBand(Lattice):
    """Two radial modes of the same angular-constant product field.

    The interval, coupling, occupation cutoff and smooth-background
    interpretation are declared. This is not a closed nonlinear mode sector.
    """

    def __init__(self, coupling=1, occupation=4):
        self.coupling = g = a.rational(coupling)
        d = a.integer(occupation, minimum=2)
        self.extended, self.matched = False, False
        self.dimensions, self.dimension = (d, d), d*d
        q2, q4, p2 = a.polynomial(d, 2), a.polynomial(d, 4), a.polynomial(d, 2, True)
        terms = []
        for index, frequency in enumerate((1, 2)):
            matrix = tuple(tuple(p2[i][j]/2+frequency**2*q2[i][j]/2+g*q4[i][j]/16
                                 for j in range(d)) for i in range(d))
            terms.append(Term((index,), (matrix,), F(1), f'product_onsite_{index}'))
        terms.append(Term((0, 1), (q2, q2), g/4, 'product_local_quartic'))
        self.terms = tuple(terms)
        self.hamiltonian = sum((self.embed(term) for term in self.terms),
                               sparse.csr_matrix((d*d, d*d)))
        self.trace = float(self.hamiltonian.diagonal().sum())


def product_gibbs_preparation(lattice, beta=F(1, 4)):
    """Source preparation of this interacting Hamiltonian's purification.

    Keep all work/reference amplitudes. Numerical synthesis guards and
    inherited exact-code failure semantics are explicit; no postselection
    or state normalization is used to repair the returned amplitudes.
    Occupations two and four keep the energy-labelled reference, including
    the unchanged original four-occupation word. Occupation eight prepares
    the canonical square-root density amplitude in the work Hermite basis
    on both registers; it needs no six-qubit arbitrary basis compiler.
    """
    beta = a.rational(beta, positive=True)
    if (not isinstance(lattice, ProductBand)
            or lattice.dimensions not in ((2, 2), (4, 4), (8, 8))):
        raise ValueError('a product band with occupation two, four or eight is required')
    dimension = lattice.dimension
    workbits = dimension.bit_length()-1
    values, basis = inherited_source._basis(lattice.hamiltonian.toarray())
    weights = np.exp(-float(beta)*(values-values[0])); weights /= weights.sum()
    if not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ArithmeticError('faithful Gibbs weights exceed the numerical preparation range')
    if lattice.dimensions == (8, 8):
        target = (basis*np.sqrt(weights))@basis.T
        word = compiler.compile_real_state(target.T.reshape(-1))
    else:
        word = compiler.compile_real_state(np.sqrt(weights))
        for bit in range(workbits):
            word += compiler.cnot(bit, bit+workbits)
        word += compiler.compile_real_orthogonal(basis)
        target = basis@np.diag(np.sqrt(weights))
    blank = np.zeros(dimension*dimension, complex); blank[0] = 1
    raw = compiler.execute(word, 2*workbits, blank).reshape(dimension, dimension).T
    error = float(np.linalg.norm(raw-target))
    if not np.isfinite(error) or error > 1e-9:
        raise ArithmeticError('the interacting Gibbs preparation exceeded its source guard')
    metadata = {'inverse_temperature': beta,
        'retained_work_addresses': dimension, 'retained_reference_addresses': dimension,
        'source_resources': compiler.resources(word),
        'raw_norm_squared': float(np.linalg.norm(raw)**2),
        'purification_vector_error_numerical': error,
        'purification_vector_error_guard': F(1, 10**9),
        'ideal_minimum_density_eigenvalue': float(weights.min()),
        'physical_preparation_error_certified': False}
    if lattice.dimensions == (8, 8):
        metadata['reference_basis'] = 'canonical_fixed_frequency_two_Hermite'
        metadata['preparation_source_word'] = 'twelve_qubit_real_amplitude_tree'
    return raw, target, metadata


@lru_cache(maxsize=1)
def product_incident_tapes():
    """All two-attempt native tapes; null incidents still trigger the clock.

    The timing carrier is independent of the protected field register. Its
    native seam/coin schedule pauses during compiled field services. The
    endpoint updates themselves are the unchanged native repair operation.
    """
    initial, port = (0,)*10+(1, 1), 10
    result = []
    for first_edge, first_seam in enumerate(native_seams()):
        for first_coin in (0, 1):
            middle, (_, first_current) = native_repair(initial, first_seam, first_coin)
            first = (first_edge, first_coin, first_current, int(port in first_seam))
            for second_edge, second_seam in enumerate(native_seams()):
                for second_coin in (0, 1):
                    final, (_, second_current) = native_repair(middle, second_seam, second_coin)
                    second = (second_edge, second_coin, second_current, int(port in second_seam))
                    result.append((first, second, final))
    return tuple(result)


def product_incident_clock_bound(level, *, degree=5, occupation=4,
                                 interval=INTERVAL, steps=1, tolerance=None,
                                 preparation_error=None, service_error=None):
    """Same interacting H, native count clock and complete instrument bound.

    Errors refer to normalized complete instruments, including failures and
    an arbitrary reference, not only successful work matrices. The default
    decreasing service errors are requirements, not hardware certificates.
    Both compared instruments retain the SAME native clock tape; closeness
    of timestamps to model time is a separate probability statement.
    """
    a.integer(level, minimum=0)
    if level > 100 or degree not in (5, 6) or type(degree) is not int:
        raise ValueError('an admitted source level and native port degree required')
    a.integer(steps)
    interval = a.rational(interval, positive=True)
    lattice = ProductBand(occupation=occupation)
    rate = 2*degree*4**level
    ticks_per_interval = ceil(rate*interval)
    ticks = 2*ticks_per_interval
    rounded_horizon = F(ticks, rate)
    eta = F(1, 2**(level//2+1)) if tolerance is None else a.rational(tolerance, positive=True)
    clock = incident_clock_bound(level, degree, rounded_horizon, eta)
    local_operations = 2*steps*len(lattice.terms)*ticks+4
    prep = F(1, rate**2) if preparation_error is None else a.rational(preparation_error)
    service = (F(1, rate**2*local_operations) if service_error is None
               else a.rational(service_error))
    split = splitting_bound(lattice, F(1, rate), steps)
    # Round an exact upper bound upward; avoid growing decimal-root fractions
    # in a clock schedule whose purpose is an explicit conservative budget.
    norm = F(ceil(split['centered_generator_norm_sum_upper']))
    total_split = ticks*min(F(2), norm**3/(3*rate**3*steps**2))
    time_error = 2*norm*(F(ticks_per_interval, rate)-interval)
    word_error = prep+local_operations*service
    failure = F(clock['failure_probability_upper'])
    return {'level': level, 'occupation_dimension_per_mode': occupation,
        'native_clock': clock, 'incident_rate_per_normalized_attempt_unit': rate,
        'model_time_per_zero_inclusive_incident': F(1, rate),
        'requested_model_time_per_read_interval': interval,
        'incident_ticks_per_read_interval': ticks_per_interval,
        'maximum_interacting_ticks': ticks,
        'palindrome_steps_per_incident': steps,
        'centered_generator_norm_sum_upper': norm,
        'maximum_guarded_local_word_calls': local_operations,
        'complete_preparation_half_trace_error_required': prep,
        'per_word_complete_instrument_half_diamond_error_required': service,
        'total_source_instrument_error_upper': word_error,
        'same_H_product_formula_error_upper': total_split,
        'nominal_read_time_rounding_error_upper': time_error,
        'same_clock_tape_complete_work_reference_record_error_upper': min(F(1), word_error+total_split),
        'nominal_work_law_with_retained_clock_tape_error_upper': min(F(1), word_error+total_split+time_error+failure),
        'timeout_probability_upper': failure,
        'model_time_and_native_arrival_time_discrepancy_upper': eta+F(1, rate),
        'modular_increment_per_incident_at_beta_one_quarter': -F(4, rate),
        'source_services_and_storage_accuracy_certified': False,
        'autonomous_clock_field_coupling_derived': False,
        'deterministic_timestamp_total_variation_claimed': False,
        'same_occupation_H_at_every_clock_level': True}


def _product_clock_prefixes(lattice, initial, *, source, interval, steps, feedback=True, sign=1):
    """Keep the work/reference state at every possible timeout prefix."""
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError('an original initial kick sign of minus or plus one is required')
    d = lattice.dimensions[1]
    q = normalized(a.polynomial(d, 1))
    if source:
        kick = inherited_source.feedback(d, sign/4, 2)[0]
        kraus = inherited_source.instrument(d, 1, np.pi/2, 2)[0]
        control = inherited_source.feedback(d, .25, 2)[0]
        advance = lambda value: lattice.source_evolve(value, interval, steps)
    else:
        kick, control = expm(1j*sign*q/4), expm(1j*q/4)
        eigenvalues, basis = np.linalg.eigh(q)
        kraus = tuple((basis*function((eigenvalues+np.pi/2)/2))@basis.T
                      for function in (np.cos, np.sin))
        advance = lambda value: lattice.evolve(value, interval)
    branches = {(): lattice.apply(initial, (1,), kick)}
    result = [branches]
    for tick in range(2):
        output = {}
        for prefix, value in branches.items():
            advanced = advance(value)
            for bit, matrix in enumerate(kraus):
                terminal = lattice.apply(advanced, (1,), matrix)
                if tick == 0 and bit == 1 and feedback:
                    terminal = lattice.apply(terminal, (1,), control)
                output[prefix+(bit,)] = terminal
        branches = output
        result.append(branches)
    return tuple(result)


def _pure_output_distance(left, right):
    """Stable half trace norm for two unnormalized rank-one outputs."""
    x, y = left.reshape(-1), right.reshape(-1)
    xx, yy = float(np.vdot(x, x).real), float(np.vdot(y, y).real)
    if xx == 0:
        return yy/2
    perpendicular = y-x*(np.vdot(x, y)/xx)
    return sqrt((xx-yy)**2+4*xx*float(np.vdot(perpendicular, perpendicular).real))/2


def product_incident_clock_control():
    """Execute this interacting Gibbs history on every native timing tape.

    Shared work computations are grouped only when the emitted pulse words
    coincide. No native record, timeout, measurement or reference output is
    postselected. The finite fixture has two attempts, not two forced ticks.
    """
    lattice, interval, steps = ProductBand(), F(1, 16), 64
    raw, target, preparation = product_gibbs_preparation(lattice)
    actual = _product_clock_prefixes(lattice, raw, source=True, interval=interval, steps=steps)
    ideal = _product_clock_prefixes(lattice, target, source=False, interval=interval, steps=steps)
    no_feedback = _product_clock_prefixes(lattice, target, source=False, interval=interval,
                                         steps=steps, feedback=False)
    tapes = product_incident_tapes()
    counts, missed_counts, null_ticks = [0]*3, [0]*3, 0
    for first, second, _ in tapes:
        counts[first[3]+second[3]] += 1
        missed_counts[sum(row[3] and row[2] != 0 for row in (first, second))] += 1
        null_ticks += sum(row[3] and row[2] == 0 for row in (first, second))
    weights = [F(count, len(tapes)) for count in counts]
    if weights != [F(25, 36), F(5, 18), F(1, 36)] or not null_ticks:
        raise ArithmeticError('the full native incident/null clock law changed')
    complete_error = sum(float(weight)*sum(_pure_output_distance(value, ideal[tick][bits])
                            for bits, value in actual[tick].items())
                         for tick, weight in enumerate(weights))
    operation_count = 4*steps*len(lattice.terms)+4
    guard = F(251, 10**12)
    gain = operation_count*guard
    vector_guard = (preparation['purification_vector_error_guard']+gain)/(1-gain)
    numerical_guard = vector_guard*(2+vector_guard)/2
    split = splitting_bound(lattice, interval, steps)['two_interval_complete_output_trace_distance_upper']
    if complete_error > float(split+numerical_guard):
        raise ArithmeticError('the joint incident-clock/interacting-field source law exceeded its bound')
    records = [{'incident_ticks': tick, 'timeout': tick < 2, 'field_bits': list(bits),
                'probability': float(weights[tick])*float(np.linalg.norm(value)**2)}
               for tick, branches in enumerate(actual) for bits, value in branches.items()]
    feedback_error = sum(float(weights[tick])*sum(_pure_output_distance(value, no_feedback[tick][bits])
                           for bits, value in branches.items())
                         for tick, branches in enumerate(ideal))
    resources = {key: preparation['source_resources'][key] for key in inherited_source.COUNT_KEYS}
    for term in lattice.terms:
        meta = execute_term(term, interval/(2*steps))[1]
        for key in resources:
            resources[key] += 4*steps*meta['resources'][key]
    for multiplier, meta in ((2, inherited_source.instrument(4, 1, np.pi/2, 2)[1]),
                              (2, inherited_source.feedback(4, .25, 2)[1])):
        for key in resources:
            resources[key] += multiplier*meta['resources'][key]
    return {'scope': 'Active zero-inclusive native incident clock for the same two-mode interacting Gibbs history.',
        'timing_carrier_ports': 12, 'timing_carrier_initial_loads': (0,)*10+(1, 1),
        'clock_port': 10, 'native_attempts_per_history': 2,
        'native_seams': len(native_seams()), 'incident_seams': 5,
        'native_tape_count': len(tapes),
        'native_tape_sha256': hashlib.sha256(canonical(tapes).encode()).hexdigest(),
        'native_tape_weight': F(1, len(tapes)),
        'tape_schema': '(edge index, placement coin, signed current, incident bit) twice, final timing-carrier loads',
        'tape_counts_by_incident_ticks': counts,
        'complete_native_tape_and_field_prefix_outcomes': sum(count*2**tick for tick, count in enumerate(counts)),
        'null_incident_occurrences_across_all_tapes': null_ticks,
        'wrong_nonzero_current_tick_counts': missed_counts,
        'nonzero_current_tick_marginal_tv': sum(abs(F(x-y, len(tapes))) for x, y in zip(counts, missed_counts))/2,
        'model_time_per_incident': interval, 'palindrome_steps_per_incident': steps,
        'finite_fixture_attempt_normalization_for_rate_sixteen': 96,
        'fixture_is_full_refined_screen_clock': False,
        'modular_parameter_per_incident_at_beta_one_quarter': -F(1, 4),
        'complete_clock_and_field_record_law': records,
        'complete_record_probability_sum': sum(row['probability'] for row in records),
        'retained_work_addresses': 16, 'retained_reference_addresses': 16,
        'timeout_probability': weights[0]+weights[1],
        'complete_native_tape_work_reference_record_trace_error_numerical': complete_error,
        'ideal_splitting_complete_output_bound': split,
        'complete_source_numerical_guard': numerical_guard,
        'feedback_removal_complete_output_distance': feedback_error,
        'worst_branch_source_resources': resources,
        'additional_native_repairs_per_history': 2,
        'source_pauses_during_field_services': True,
        'protected_persistent_work_register_required': True,
        'all_tapes_timeout_prefixes_and_reference_outputs_retained': True,
        'physical_service_or_lab_second_calibration_certified': False,
        'clock_refinement_bounds': [product_incident_clock_bound(level) for level in (2, 5, 8)]}


@dataclass(frozen=True)
class ProductControllerState:
    """Classical registers of the fixed finite interacting protocol.

    The instruction pointer is part of the declared apparatus. Neither a
    host iteration nor a quantum amplitude is an input to its transition.
    Native records include null currents; field outcomes are retained bits.
    """
    instruction: str = 'kick'
    cursor: int = 0
    attempts: int = 0
    bits: tuple = ()
    loads: tuple = (0,)*10+(1, 1)
    tape: tuple = ()
    failure: tuple = ()


def product_controller_transition(state, outcome=None, *, steps=STEPS):
    """One time-independent classical control transition.

    ``idle`` emits a retained identity outcome. A service failure is an
    absorbing output, never a request to normalize or retry its work state.
    Only the ready instruction enables a native timing-carrier repair.
    """
    a.integer(steps)
    if (not isinstance(state, ProductControllerState)
            or state.instruction not in ('kick', 'ready', 'evolve', 'read', 'feedback', 'halt', 'failure')
            or type(state.attempts) is not int or not 0 <= state.attempts <= 2
            or type(state.tape) is not tuple or type(state.bits) is not tuple
            or type(state.loads) is not tuple or type(state.failure) is not tuple
            or len(state.tape) != state.attempts or len(state.bits) > 2
            or any(type(bit) is not int or bit not in (0, 1) for bit in state.bits)
            or len(state.loads) != 12 or any(type(x) is not int or x not in (0, 1) for x in state.loads)
            or sum(state.loads) != 2 or type(state.cursor) is not int
            or not 0 <= state.cursor < 6*steps):
        raise ValueError('valid finite classical controller registers required')
    # The admitted classical domain is a consistent checkpoint, not every
    # bit pattern of a hardware register. Reverse the retained currents and
    # replay the two-event maximum tape to check its actual native meaning.
    loads = list(state.loads)
    for record in reversed(state.tape):
        if (type(record) is not tuple or len(record) != 4
                or any(type(x) is not int for x in record)
                or not 0 <= record[0] < 30 or record[1] not in (0, 1)
                or record[2] not in (-1, 0, 1) or record[3] not in (0, 1)):
            raise ValueError('a complete native tape is required')
        seam = native_seams()[record[0]]
        if record[3] != int(10 in seam):
            raise ValueError('the retained incidence must match its native seam')
        loads[seam[0]] += record[2]; loads[seam[1]] -= record[2]
        if any(value not in (0, 1) for value in loads):
            raise ValueError('native records must reverse to the admitted binary sector')
    for edge, coin, current, _ in state.tape:
        loads, (_, actual) = native_repair(tuple(loads), native_seams()[edge], coin)
        if actual != current:
            raise ValueError('retained currents must be actual native outcomes')
    if tuple(loads) != state.loads:
        raise ValueError('retained native loads and tape disagree')
    phase = state.instruction
    if phase == 'failure':
        if (len(state.failure) != 3 or state.failure[0] not in ('kick', 'evolve', 'read', 'feedback')
                or state.failure[1] != state.cursor or type(state.failure[2]) is not int
                or not 2 <= state.failure[2] < 6):
            raise ValueError('an absorbing failure needs its source instruction and complement label')
        phase = state.failure[0]
    elif state.failure:
        raise ValueError('a live instruction cannot carry a terminal failure')
    incidents = sum(record[3] for record in state.tape)
    pending = int(phase in ('evolve', 'read'))
    if (incidents != len(state.bits)+pending
            or (phase in ('evolve', 'read') and len(state.bits) >= 2)
            or (phase == 'feedback' and state.bits != (1,))
            or (phase == 'kick' and (state.attempts or state.bits))
            or (phase == 'halt' and state.attempts != 2)
            or (phase != 'evolve' and state.cursor)):
        raise ValueError('instruction, incident count and retained field prefix must be consistent')
    if outcome == 'idle':
        return state
    if state.instruction in ('halt', 'failure'):
        if outcome is not None:
            raise ValueError('an absorbing controller has no further active outcome')
        return state
    if isinstance(outcome, tuple) and len(outcome) == 2 and outcome[0] == 'failure':
        if state.instruction not in ('kick', 'evolve', 'read', 'feedback'):
            raise ValueError('only an active source service can emit a service failure')
        if type(outcome[1]) is not int or not 2 <= outcome[1] < 6:
            raise ValueError('the retained six-level source complement label is required')
        return replace(state, instruction='failure', failure=(state.instruction, state.cursor, outcome[1]))
    if state.instruction == 'ready':
        if state.attempts == 2:
            if outcome is not None:
                raise ValueError('the finite attempt budget is exhausted')
            return replace(state, instruction='halt')
        if (not isinstance(outcome, tuple) or len(outcome) != 2
                or any(type(x) is not int for x in outcome)
                or not 0 <= outcome[0] < 30 or outcome[1] not in (0, 1)):
            raise ValueError('a native seam index and placement coin are required')
        edge, coin = outcome
        seam = native_seams()[edge]
        loads, (_, current) = native_repair(state.loads, seam, coin)
        incident = int(10 in seam)
        record = (edge, coin, current, incident)
        return replace(state, instruction='evolve' if incident else 'ready', cursor=0,
                       attempts=state.attempts+1, loads=loads, tape=state.tape+(record,))
    if state.instruction == 'read':
        if type(outcome) is not int or outcome not in (0, 1):
            raise ValueError('a committed binary pointer outcome is required')
        bits = state.bits+(outcome,)
        return replace(state, instruction='feedback' if bits == (1,) else 'ready', bits=bits)
    if outcome is not None:
        raise ValueError('this source instruction has no classical choice')
    if state.instruction in ('kick', 'feedback'):
        return replace(state, instruction='ready')
    following = state.cursor+1
    return replace(state, instruction='read' if following == 6*steps else 'evolve',
                   cursor=0 if following == 6*steps else following)


class ProductControllerInstrument:
    """A single homogeneous instrument on control, timer and work registers.

    Its local service matrices are contractions of the existing emitted
    source words. The ideal channel uses their exact code identities; the
    numerical representation retains the inherited guarded roundoff policy.
    A static read-only instruction table, source services and protected
    registers remain construction inputs, not consequences of passive repair.
    The initial kick has the declared sign; record-selected feedback is
    always +1/4. Native transitions and durable prefixes are unchanged.
    """

    def __init__(self, *, occupation=4, sign=1, interval=INTERVAL, steps=STEPS):
        if type(occupation) is not int or occupation not in (2, 4, 8):
            raise ValueError('an original occupation of two, four or eight is required')
        if type(sign) is not int or sign not in (-1, 1):
            raise ValueError('an original initial kick sign of minus or plus one is required')
        self.occupation, self.sign = occupation, sign
        self.lattice = ProductBand(occupation=occupation)
        self.interval, self.steps = a.rational(interval, positive=True), a.integer(steps)
        self.kick, self.kick_metadata = inherited_source.feedback(occupation, sign/4, 2)
        self.feedback, self.feedback_metadata = inherited_source.feedback(occupation, .25, 2)
        self.kraus, self.read_metadata = inherited_source.instrument(occupation, 1, np.pi/2, 2)
        local = [(term.sites, *execute_term(term, self.interval/(2*self.steps)))
                 for term in self.lattice.terms]
        self.cycle = tuple(local+local[::-1])
        self.prefix_transition_classes = set()

    def _checked_branches(self, state, branches):
        for following, _, _ in branches:
            if (following.tape[:len(state.tape)] != state.tape
                    or following.bits[:len(state.bits)] != state.bits):
                raise ArithmeticError('a controller transition changed an already committed prefix')
            # Busy words copy native data, so their instruction classes are
            # checked once; every actually enumerated ready branch is also
            # checked on its full native tape before this diagnostic grouping.
            self.prefix_transition_classes.add((state.instruction, state.cursor,
                state.attempts, state.bits, following.instruction, following.cursor,
                following.attempts, following.bits))
        return branches

    def branches(self, state):
        """Return complete active outcomes (next control, work map, weight).

        Native weights are 1/60 and multiply amplitudes by their square
        roots. Born probabilities reside in the read Kraus maps. Invalid
        proper-code exports are the inherited complement channel, routed
        by ``product_controller_transition`` to the failure instruction.
        They have exactly zero weight on the declared valid-code input.
        At an absorbing failure, ``None`` denotes identity on its retained
        failure-sector payload, whose dimension need not equal the work code.
        """
        product_controller_transition(state, 'idle', steps=self.steps)
        if state.instruction == 'failure':
            return self._checked_branches(state, ((product_controller_transition(state, steps=self.steps), None, F(1)),))
        identity = np.eye(self.lattice.dimension, dtype=complex)
        if state.instruction == 'ready' and state.attempts < 2:
            return self._checked_branches(state, tuple((product_controller_transition(state, (edge, coin), steps=self.steps),
                          identity, F(1, 60)) for edge in range(30) for coin in (0, 1)))
        if state.instruction == 'read':
            return self._checked_branches(state, tuple((product_controller_transition(state, bit, steps=self.steps),
                          self.lattice.apply(identity, (1,), matrix), F(1))
                         for bit, matrix in enumerate(self.kraus)))
        if state.instruction in ('kick', 'feedback'):
            matrix = self.kick if state.instruction == 'kick' else self.feedback
            operation = self.lattice.apply(identity, (1,), matrix)
        elif state.instruction == 'evolve':
            sites, matrix, _ = self.cycle[state.cursor % len(self.cycle)]
            operation = self.lattice.apply(identity, sites, matrix)
        else:
            operation = identity
        return self._checked_branches(state, ((product_controller_transition(state, steps=self.steps), operation, F(1)),))

    @lru_cache(maxsize=8)
    def service_closure(self, instruction, bits):
        """Exact common-subroutine reuse, ending at the next ready state.

        This contracts the fixed transition table on every work-input
        columns. It has no access to the native state or retained tape.
        Each call records its instruction count, including conditional reads
        and feedback; this is not substitution by an ideal Hamiltonian.
        """
        if instruction == 'kick' and bits == ():
            start = ProductControllerState()
        elif instruction == 'evolve' and len(bits) < 2:
            # A valid representative checkpoint for a service that is
            # independent of its timer contents. This is an actual native
            # tape, not fabricated currents or an inaccessible field value.
            edge = next(i for i, seam in enumerate(native_seams()) if 10 in seam)
            loads, tape = ProductControllerState().loads, ()
            for _ in range(len(bits)+1):
                loads, (_, current) = native_repair(loads, native_seams()[edge], 0)
                tape += ((edge, 0, current, 1),)
            start = ProductControllerState(instruction='evolve', bits=bits,
                attempts=len(tape), loads=loads, tape=tape)
        else:
            raise ValueError('a reachable source subroutine entry is required')
        todo = [(start, np.eye(self.lattice.dimension, dtype=complex), 0)]
        complete = []
        while todo:
            state, operation, count = todo.pop()
            if state.instruction == 'ready':
                complete.append((state.bits, operation, count))
                continue
            for following, local, weight in self.branches(state):
                if weight != 1:
                    raise ArithmeticError('a busy source subroutine attempted native repair')
                todo.append((following, local@operation, count+1))
        return tuple(sorted(complete, key=lambda row: row[0]))

    def completed(self, initial_loads=(0,)*10+(1, 1)):
        """All terminal native/field records and their all-input Kraus maps."""
        start = ProductControllerState(loads=tuple(initial_loads))
        todo = [(start, np.eye(self.lattice.dimension, dtype=complex), F(1), 0)]
        complete = []
        while todo:
            state, operation, weight, count = todo.pop()
            if state.instruction in ('halt', 'failure'):
                complete.append((state, operation, weight, count))
            elif state.instruction in ('kick', 'evolve'):
                for bits, local, length in self.service_closure(state.instruction, state.bits):
                    following = replace(state, instruction='ready', cursor=0, bits=bits)
                    todo.append((following, local@operation, weight, count+length))
            else:
                for following, local, probability in self.branches(state):
                    todo.append((following, local@operation, weight*probability, count+1))
        return tuple(sorted(complete, key=lambda row: (row[0].tape, row[0].bits)))


def product_controller_idle_bound(progress_steps, invocations, idle_probability=F(1, 2)):
    """Finite unfinished-output bound for the SAME homogeneous instrument.

    Laziness means q Id+(1-q) T, independent of control and quantum state.
    Idle records and every unfinished work/reference state are retained.
    The Chernoff bound compares a binomial progress count with the worst
    admitted active length; no physical duration is assigned to invocations.
    """
    a.integer(progress_steps); a.integer(invocations, minimum=0)
    q = a.rational(idle_probability)
    if q >= 1:
        raise ValueError('idle probability must be below one')
    mean = (1-q)*invocations
    if mean < progress_steps:
        bound = F(1)
    else:
        # exp(-x) <= 1/(1+x) keeps a conservative rational certificate.
        exponent = (mean-progress_steps)**2/(2*mean)
        bound = 1/(1+exponent)
    return {'maximum_active_instructions': progress_steps, 'controller_invocations': invocations,
            'idle_probability': q, 'unfinished_probability_upper': bound,
            'expected_invocations_upper': F(progress_steps)/(1-q),
            'unfinished_states_and_idle_records_retained': True,
            'controller_invocations_are_native_clock_ticks': False}


def product_controller_export(state, dimension, *, steps=STEPS):
    """Route the existing full six-level export to success or retained failure.

    Success leaves the surrounding service instruction active. Its internal
    source-word cursor advances within the supplied service implementation;
    the outer controller advances only after that word completes.
    """
    if not isinstance(state, ProductControllerState) or state.instruction not in ('kick', 'evolve', 'read', 'feedback'):
        raise ValueError('an active controller source service is required')
    exports = compiler.transfer_kraus(dimension)
    return ((state, exports[0], 'code'),)+tuple(
        (product_controller_transition(state, ('failure', label), steps=steps), matrix, label)
        for label, matrix in zip(range(dimension, 6), exports[1:], strict=True))


def product_homogeneous_controller_control():
    """Audit the fixed controller against the separately defined history law."""
    program = ProductControllerInstrument()
    completed = program.completed()
    direct = _product_clock_prefixes(program.lattice, np.eye(16), source=True,
                                    interval=INTERVAL, steps=STEPS)
    total, maximum_error, lengths, tape_counts = np.zeros((16, 16), complex), 0., {}, {}
    raw, _, preparation = product_gibbs_preparation(program.lattice)
    complete_law, rows = [], []
    for state, operation, weight, count in completed:
        target = direct[len(state.bits)][state.bits]
        maximum_error = max(maximum_error, float(np.linalg.norm(operation-target, 2)))
        total += float(weight)*(operation.conj().T@operation)
        lengths[count] = lengths.get(count, 0)+1
        tape_counts[state.tape] = tape_counts.get(state.tape, 0)+1
        complete_law.append((state.tape, state.bits, float(weight)*float(np.linalg.norm(operation@raw)**2)))
        rows.append((state.tape, state.bits, state.loads, count))
    completeness = float(np.linalg.norm(total-np.eye(16), 2))
    expected_tapes = {(first, second): final for first, second, final in product_incident_tapes()}
    if (len(completed) != 4900 or set(tape_counts) != set(expected_tapes)
            or any(state.loads != expected_tapes[state.tape] for state, _, _, _ in completed)
            or maximum_error > 2e-12 or completeness > 2e-11):
        raise ArithmeticError('the homogeneous controller changed the complete all-input history instrument')
    failures = {}
    for dimension in (2, 4):
        exports = product_controller_export(ProductControllerState(), dimension)
        failures[str(dimension)] = {
            'complete_export_residual': float(np.linalg.norm(sum(k.conj().T@k for _, k, _ in exports)-np.eye(6), 2)),
            'complement_labels': list(range(dimension, 6)),
            'valid_code_failure_probability': 0,
            'failure_instruction_absorbing': all(product_controller_transition(
                product_controller_transition(ProductControllerState(), ('failure', label)))
                .instruction == 'failure' for label in range(dimension, 6))}
    law_rows = []
    for tick in range(3):
        for bits in direct[tick]:
            law_rows.append({'incident_ticks': tick, 'timeout': tick < 2, 'field_bits': list(bits),
                'probability': sum(probability for _, outcome, probability in complete_law if outcome == bits)})
    maximum_length = max(lengths)
    # Persistence is checked after commitment, not during record creation.
    # Include idle and genuine export failure transitions on reachable busy
    # checkpoints and the absorbing terminal states.
    ready = product_controller_transition(ProductControllerState())
    first_edge = next(i for i, seam in enumerate(native_seams()) if 10 in seam)
    written = product_controller_transition(ready, (first_edge, 0))
    probes = [written, completed[-1][0]]
    for state in probes:
        program._checked_branches(state, ((product_controller_transition(state, 'idle'), None, F(1)),))
        if state.instruction == 'evolve':
            for failed, _, label in product_controller_export(state, 4):
                if label != 'code':
                    program._checked_branches(state, ((failed, None, F(1)),))
                    program.branches(failed)
        else:
            program.branches(state)
    prefix_certificate = {
        'checked_instruction_transition_classes': len(program.prefix_transition_classes),
        'instruction_transition_classes_sha256': hashlib.sha256(canonical(sorted(program.prefix_transition_classes)).encode()).hexdigest(),
        'all_emitted_transitions_preserve_existing_native_and_field_prefixes': True,
        'postcommit_public_projectors_satisfy_dual_fixed_point_identity': True,
        'generator': 'T* - I on the invariant postcommit controller domain',
        'readback': 'sum_p Tr(Q_p X Q_p)/rank(Q_p) Q_p; block ranks may differ',
        'state': 'actual first native write after the source Gibbs preparation and kick',
        'new_record_creation_globally_preserves_that_projector': False,
        'blank_to_first_selected_record_probability': F(1, 60),
        'OperationalObserver_interface_realized_on_this_program': True,
        'Lean_inhabitant_compiled': False,
        'unused_buffer_refinement_architecture_supplied': True,
        'arbitrary_quantum_work_observables_claimed_durable': False}
    resources = {key: preparation['source_resources'][key] for key in inherited_source.COUNT_KEYS}
    for _, _, metadata in program.cycle:
        for key in resources:
            resources[key] += 2*STEPS*metadata['resources'][key]
    for metadata in (program.kick_metadata, program.read_metadata):
        for key in resources:
            resources[key] += 2*metadata['resources'][key]
    return {'scope': 'One finite homogeneous controller for the existing native-clock interacting field; supplied controller architecture and source services.',
        'complete_native_tapes': len(tape_counts), 'complete_tape_field_outputs': len(completed),
        'retained_work_input_columns': 16, 'retained_Gibbs_reference_addresses': 16,
        'controller_registers': ['instruction', 'cursor', 'attempts', 'field bits', 'native loads', 'native tape', 'failure'],
        'external_iteration_or_quantum_amplitudes_read_by_controller': False,
        'maximum_active_instructions': maximum_length,
        'worst_branch_source_resources_including_preparation': resources,
        'additional_native_repairs_per_history': 2,
        'source_word_services_are_atomic_at_this_controller_level': True,
        'terminal_instruction_counts': [{'instructions': count, 'outcomes': amount} for count, amount in sorted(lengths.items())],
        'complete_controller_record_sha256': hashlib.sha256(canonical(rows).encode()).hexdigest(),
        'all_input_Kraus_comparison_residual': maximum_error,
        'complete_all_input_channel_residual': completeness,
        'complete_clock_field_record_law': law_rows,
        'total_output_probability': sum(row['probability'] for row in law_rows),
        'native_repairs_are_enabled_only_at_ready_instruction': True,
        'postcommit_prefix_observer': prefix_certificate,
        'service_complement_outputs': failures,
        'finite_idle_bounds': [product_controller_idle_bound(maximum_length, multiple*maximum_length)
                               for multiple in (2, 4, 8)],
        'Gibbs_preparation_is_a_source_prefix': preparation['source_resources'],
        'all_input_equality_covers_arbitrary_references_and_classical_native_correlations': True,
        'coherent_native_load_superposition_claimed': False,
        'passive_repair_selected_the_controller_or_interacting_H': False,
        'controller_hardware_storage_and_source_service_accuracy_derived': False,
        'microinstruction_count_calibrated_as_physical_time': False}


def _product_history_maps(advance, kick, reads, feedback):
    """Complete ideal work maps, retaining every stopped read prefix."""
    branches = {(): kick}
    prefixes = [branches]
    for tick in range(2):
        following = {}
        for bits, operation in branches.items():
            for bit, read in enumerate(reads):
                value = read@advance@operation
                if tick == 0 and bit == 1:
                    value = feedback@value
                following[bits+(bit,)] = value
        branches = following
        prefixes.append(branches)
    return tuple(prefixes)


def product_modular_history_control():
    """One reference clock for the complete controlled field instrument.

    The fixed Gibbs reference determines the modular evolution. Its ideal
    channel identity includes interventions and arbitrary work references;
    the executed source word agrees only within its splitting/service
    budget. A conditioned work state is not silently used as that reference.
    Translated detector menus and relative modular cocycles are mathematical
    controls, not additional compiled source protocols.
    """
    program = ProductControllerInstrument()
    lattice, beta, interval = program.lattice, F(1, 4), program.interval
    raw, target, preparation = product_gibbs_preparation(lattice, beta)
    rho = target@target.conj().T
    h = lattice.hamiltonian.toarray()
    values, basis = np.linalg.eigh(rho)
    logarithm = (basis*np.log(values))@basis.conj().T
    modular_h = -logarithm/float(beta)
    scalar = float(np.trace(modular_h-h).real/16)
    s = -float(interval/beta)
    # Schrödinger propagation is rho^(-is), adjoint to sigma_s.
    modular_advance = expm(-1j*s*logarithm)
    advance = expm(-1j*float(interval)*h)
    phase = np.exp(1j*float(interval)*scalar)
    q = normalized(a.polynomial(4, 1))
    q_values, q_basis = np.linalg.eigh(q)
    local_reads = tuple((q_basis*f((q_values+np.pi/2)/2))@q_basis.T
                        for f in (np.cos, np.sin))
    identity = np.eye(16)
    kick = lattice.apply(identity, (1,), expm(1j*q/4))
    reads = tuple(lattice.apply(identity, (1,), matrix) for matrix in local_reads)
    physical = _product_history_maps(advance, kick, reads, kick)
    modular = _product_history_maps(modular_advance, kick, reads, kick)
    ideal_residual = max(float(np.linalg.norm(physical[tick][bits]
                            -phase**tick*operation, 2))
                         for tick, row in enumerate(modular) for bits, operation in row.items())
    completed = program.completed()
    difference_gram = np.zeros((16, 16), complex)
    source_gram = np.zeros_like(difference_gram)
    output_error = 0.
    for state, operation, weight, _ in completed:
        tick = len(state.bits)
        expected = phase**tick*modular[tick][state.bits]
        difference = operation-expected
        difference_gram += float(weight)*(difference.conj().T@difference)
        source_gram += float(weight)*(operation.conj().T@operation)
        output_error += float(weight)*_pure_output_distance(operation@raw, expected@target)
    stinespring_error = sqrt(max(0., float(np.linalg.eigvalsh(difference_gram)[-1])))
    # This is a numerical all-reference bound, not a hardware certificate.
    channel_error = stinespring_error*(1+sqrt(float(np.linalg.eigvalsh(source_gram)[-1])))/2
    count = 4*program.steps*len(lattice.terms)+4
    gain = count*F(251, 10**12)
    vector_guard = (preparation['purification_vector_error_guard']+gain)/(1-gain)
    source_guard = vector_guard*(2+vector_guard)/2
    split = splitting_bound(lattice, interval, program.steps)['two_interval_complete_output_trace_distance_upper']

    # A common time change conjugates the whole intervention menu, including
    # the conditional feedback. It is a channel identity on every input.
    translation = expm(-1j*float(interval)*h)
    translate = lambda value: translation@value@translation.conj().T
    shifted = _product_history_maps(advance, translate(kick),
                                   tuple(map(translate, reads)), translate(kick))
    wrong_feedback = _product_history_maps(advance, translate(kick),
                                          tuple(map(translate, reads)), kick)
    covariance = max(float(np.linalg.norm(shifted[tick][bits]-translate(operation), 2))
                     for tick, row in enumerate(physical) for bits, operation in row.items())
    feedback_mismatch = max(float(np.linalg.norm(wrong_feedback[tick][bits]
                               -translate(operation), 2))
                            for tick, row in enumerate(physical) for bits, operation in row.items())

    # The selected posterior is faithful, but disturbance changes its
    # modular generator. The finite relative cocycle gives the exact repair
    # of the representation without identifying the two reference states.
    branch = physical[1][(1,)]@target
    posterior = branch@branch.conj().T/float(np.linalg.norm(branch)**2)
    masses, vectors = np.linalg.eigh(posterior)
    if masses.min() <= 0:
        raise ArithmeticError('the selected posterior is not faithful')
    log_posterior = (vectors*np.log(masses))@vectors.conj().T
    ref_unitary = expm(1j*s*logarithm)
    post_unitary = expm(1j*s*log_posterior)
    observable = lattice.apply(identity, (1,), q)
    ref_action = ref_unitary@observable@ref_unitary.conj().T
    post_action = post_unitary@observable@post_unitary.conj().T
    cocycle = post_unitary@ref_unitary.conj().T
    cocycle_residual = float(np.linalg.norm(cocycle.conj().T@post_action@cocycle-ref_action, 2))
    posterior_mismatch = float(np.linalg.norm(post_action-ref_action, 2))
    second = s/2
    relative = lambda t: expm(1j*t*log_posterior)@expm(-1j*t*logarithm)
    cocycle_law = float(np.linalg.norm(relative(s+second)
        -cocycle@ref_unitary@relative(second)@ref_unitary.conj().T, 2))
    parity = lattice.apply(identity, (0,), np.diag([1, -1, 1, -1]))
    parity_residual = max(float(np.linalg.norm(parity@value-value@parity, 2))
                          for value in (h, rho, kick, *reads, posterior, cocycle))
    if (ideal_residual > 1e-11 or covariance > 1e-11 or cocycle_residual > 1e-11
            or cocycle_law > 1e-11 or parity_residual > 1e-11
            or feedback_mismatch < 1e-5 or posterior_mismatch < 1e-5
            or max(channel_error, output_error) > float(split+source_guard)):
        raise ArithmeticError('the complete fixed-reference modular history contract failed')
    return {'scope': 'Same finite product program and fixed Gibbs reference; compiled errors retained, no posterior-reference substitution.',
        'inverse_temperature': beta, 'model_time_per_incident': interval,
        'fixed_reference_modular_parameter_per_incident': -interval/beta,
        'complete_native_tapes': 3600, 'complete_native_and_field_outputs': len(completed),
        'retained_work_input_columns': 16, 'retained_reference_addresses': 16,
        'ideal_all_input_history_Kraus_residual_after_scalar_phase': ideal_residual,
        'source_to_modular_Stinespring_error_numerical': stinespring_error,
        'source_to_modular_all_reference_half_diamond_bound_numerical': channel_error,
        'source_prepared_complete_output_half_trace_error_numerical': output_error,
        'source_splitting_and_service_complete_output_bound': split+source_guard,
        'full_intervention_menu_covariance_residual': covariance,
        'untranslated_feedback_covariance_error': feedback_mismatch,
        'conditioned_reference_minimum_eigenvalue_numerical': float(masses.min()),
        'conditioned_reference_modular_action_mismatch': posterior_mismatch,
        'relative_modular_cocycle_action_residual': cocycle_residual,
        'relative_modular_cocycle_law_residual': cocycle_law,
        'first_mode_parity_commutator_residual': parity_residual,
        'same_native_clock_tape_retained_in_both_channels': True,
        'fixed_reference_replaced_by_measurement_posterior': False,
        'translated_menu_or_cocycle_compiled_as_extra_protocol': False,
        'Gibbs_reference_and_controller_selected_by_passive_repair': False,
        'physical_time_or_own_algebra_modular_flow_for_classical_records_identified': False}


def product_interaction_control():
    """Conditional same-product interaction, Gibbs clock and full readouts."""
    lattice, free = ProductBand(), ProductBand(0)
    beta, interval, steps = F(1, 4), F(1, 16), 64
    raw, target, preparation = product_gibbs_preparation(lattice, beta)
    actual, reference, actual_states, reference_states = [], [], [], []
    for sign in (-1, 1):
        rows, states = history(lattice, sign, interval=interval, steps=steps,
                               initial=raw, return_states=True)
        ideal_rows, ideal_states = history(lattice, sign, source=False, interval=interval,
                                           initial=target, return_states=True)
        actual.append(rows); actual_states.append(states)
        reference.append(ideal_rows); reference_states.append(ideal_states)
    def distance(left, right):
        return sum(float(np.linalg.norm(np.outer(x.reshape(-1), x.reshape(-1).conj())
                                      -np.outer(y.reshape(-1), y.reshape(-1).conj()), ord='nuc'))/2
                   for x, y in zip(left, right, strict=True))
    complete_error = max(distance(x, y) for x, y in zip(actual_states, reference_states, strict=True))
    split = splitting_bound(lattice, interval, steps)
    certificate = split['two_interval_complete_output_trace_distance_upper']
    # Telescope the inherited fixed raw-target and polar-correction guards,
    # including the complete measurement dilations and both possible kicks.
    # This is a floating-replay allowance, not a physical service certificate.
    operation_count = 4*steps*len(lattice.terms)+4
    local_guard = F(251, 10**12)  # 2e-10 raw residual + 5.1e-11 correction.
    gain = operation_count*local_guard
    vector_guard = (preparation['purification_vector_error_guard']+gain)/(1-gain)
    numerical_guard = vector_guard*(2+vector_guard)/2
    if complete_error > float(certificate+numerical_guard):
        raise ArithmeticError('the full interacting work/reference/record law exceeded its bound')
    # Compare the same input and detectors after removing only the interaction.
    free_rows = [history(free, sign, source=False, interval=interval, initial=target)
                 for sign in (-1, 1)]
    no_feedback = [history(lattice, sign, source=False, feedback=False,
                           interval=interval, initial=target) for sign in (-1, 1)]
    h = lattice.hamiltonian.toarray(); rho = target@target.conj().T
    q = np.kron(normalized(a.polynomial(4, 1)), np.eye(4))
    t = float(interval)
    evolved = expm(1j*t*h)@q@expm(-1j*t*h)
    masses, vectors = np.linalg.eigh(rho)
    modular = (vectors*np.exp(-1j*t*np.log(masses)/float(beta)))@vectors.T
    modular_residual = float(np.linalg.norm(modular@q@modular.conj().T-evolved))
    raw_masses, raw_vectors = np.linalg.eigh(raw@raw.conj().T)
    raw_modular = (raw_vectors*np.exp(-1j*t*np.log(raw_masses)/float(beta)))@raw_vectors.conj().T
    raw_modular_residual = float(np.linalg.norm(raw_modular@q@raw_modular.conj().T-evolved))
    free_values, free_basis = np.linalg.eigh(free.hamiltonian.toarray())
    free_weights = np.exp(-float(beta)*free_values); free_weights /= free_weights.sum()
    free_rho = (free_basis*free_weights)@free_basis.T
    free_state_nonstationarity = float(np.linalg.norm(h@free_rho-free_rho@h))
    if modular_residual > 1e-10 or raw_modular_residual > 1e-9 or free_state_nonstationarity < 1e-5:
        raise ArithmeticError('the interacting clock/state identity or mismatched-state control failed')
    native = native_records(actual)
    resources = dict(preparation['source_resources'])
    counts = {key: resources[key] for key in inherited_source.COUNT_KEYS}
    for term in lattice.terms:
        meta = execute_term(term, interval/(2*steps))[1]
        for key in counts:
            counts[key] += 4*steps*meta['resources'][key]
    for count, meta in ((2, inherited_source.instrument(4, 1, np.pi/2, 2)[1]),
                        (2, inherited_source.feedback(4, .25, 2)[1])):
        for key in counts:
            counts[key] += count*meta['resources'][key]
    return {'scope': 'Two-mode interacting effective field on the declared radial product; finite ultraviolet and occupation cutoffs.',
        'spatial_domain': '(0,pi) times the same native angular screen',
        'field_modes': ['sqrt(2) sin(r) times angular constant', 'sqrt(2) sin(2r) times angular constant'],
        'spatial_measure': 'dr/pi times normalized native angular measure',
        'quartic_polynomial': '3/2*(Q1^4+Q2^4)+6*Q1^2*Q2^2',
        'effective_coupling': lattice.coupling, 'oscillator_basis_frequency': 2,
        'covariant_coupling_map': 'g_eff=lambda/pi integral mu_screen^2 sqrt(det a_screen) dx',
        'occupation_dimensions': lattice.dimensions, 'preparation': preparation,
        'model_time_per_interval': interval, 'steps_per_interval': steps,
        'source_histories_by_native_sign': actual,
        'independent_hamiltonian_histories_by_native_sign': reference,
        'complete_work_reference_record_trace_error_numerical': complete_error,
        'ideal_splitting_complete_output_bound': certificate,
        'complete_source_numerical_guard': numerical_guard,
        'guarded_local_operation_count': operation_count,
        'local_raw_plus_polar_operator_guard': local_guard,
        'interaction_removal_record_tv': max(tv(x, y) for x, y in zip(reference, free_rows)),
        'feedback_removal_record_tv': max(tv(x, y) for x, y in zip(reference, no_feedback)),
        'ideal_Gibbs_modular_identity_residual': modular_residual,
        'raw_prepared_marginal_modular_identity_residual': raw_modular_residual,
        'free_Gibbs_state_interacting_commutator_norm': free_state_nonstationarity,
        'modular_clock': 'sigma_s=alpha_(-beta*s) for this interacting H and its ideal Gibbs marginal',
        'clock_parameter': 'Supplied Hamiltonian model time; no event-counter or seconds calibration in this protocol.',
        'native_history_count': len(native),
        'native_history_probability_sum': sum(row['probability'] for row in native),
        'worst_branch_source_resources': counts,
        'additional_native_repairs_per_history': 2,
        'physical_service_error_certified': False,
        'conditional_smooth_background_requires_smooth_coefficient_density': True,
        'coupling_and_temperature_selected_by_repair': False,
        'finite_band_is_invariant_under_untruncated_local_interaction': False,
        'UV_limit_or_cutoff_microcausality_certified': False,
        'observed_QFT_or_laboratory_clock_identified': False}


def product_parameter_recovery(source_histories, raw_preparation):
    """Synthetic inverse check with an independent padded-oscillator oracle.

    Fit coupling, inverse temperature and radial frequency scale from the
    first interval only. A second interval is withheld. Rounded estimates
    are retained under explicit numerical guards, not statistical intervals.
    """
    annihilation = np.diag(np.sqrt(np.arange(1, 8)), 1)
    q = (annihilation+annihilation.T)/2
    p = 1j*(annihilation.T-annihilation)
    q2, q4, p2 = (matrix[:4, :4] for matrix in (q@q, np.linalg.matrix_power(q, 4), p@p))

    def oracle(parameters, interval=F(1, 16)):
        coupling, beta, scale = parameters
        h1, h2 = (p2/2+(scale*k)**2*q2/2+coupling*q4/16 for k in (1, 2))
        h = np.kron(np.eye(4), h1)+np.kron(h2, np.eye(4))+coupling*np.kron(q2, q2)/4
        lattice = ProductBand()
        lattice.hamiltonian, lattice.trace = sparse.csr_matrix(h), float(np.trace(h).real)
        values, vectors = np.linalg.eigh(h)
        weights = np.exp(-beta*(values-values[0])); weights /= weights.sum()
        initial = vectors*np.sqrt(weights)
        return np.array([[row['probability'] for row in history(
            lattice, sign, source=False, initial=initial, interval=interval)] for sign in (-1, 1)])

    target = np.array([[row['probability'] for row in rows] for rows in source_histories])
    select = lambda probabilities: probabilities[:, :3].ravel()
    starts = ([.6, .4, .8], [1.6, .15, 1.2], [2., .8, .7])
    fits = [least_squares(lambda x: select(oracle(x)-target), start,
                         bounds=([.01, .03, .5], [3., 2., 2.]),
                         xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=150)
            for start in starts]
    expected = np.array([1., .25, 1.])  # Used only for the synthetic acceptance check.
    if any(not fit.success or np.max(abs(fit.x-expected)) >= 1e-4
           or np.max(abs(fit.fun)) >= 1e-8 for fit in fits):
        raise ArithmeticError('the synthetic observer records failed parameter recovery')
    estimate = fits[0].x
    step = 1e-5
    jacobian = np.column_stack([(select(oracle(estimate+np.eye(3)[i]*step))
                                -select(oracle(estimate-np.eye(3)[i]*step)))/(2*step)
                               for i in range(3)])
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    withheld = np.array([[row['probability'] for row in history(
        ProductBand(), sign, initial=raw_preparation, interval=F(1, 8))] for sign in (-1, 1)])
    if singular_values[-1] <= 1e-5 or np.max(abs(oracle(estimate, F(1, 8))-withheld)) >= 1e-8:
        raise ArithmeticError('the local identifiability or withheld-interval check failed')
    return {'data_kind': 'synthetic source-executed complete probabilities, without sampling noise',
            'parameters': ['effective coupling', 'inverse temperature', 'radial frequency scale'],
            'estimates_rounded_to_1e_minus_4': np.round(estimate, 4).tolist(),
            'three_distinct_initial_guesses_agree_with_declared_inputs_within': F(1, 10000),
            'training_max_probability_residual_guard': F(1, 10**8),
            'training_model_interval': F(1, 16), 'withheld_model_interval': F(1, 8),
            'withheld_max_probability_error_guard': F(1, 10**8),
            'local_sensitivity_singular_values_rounded_to_1e_minus_6': np.round(singular_values, 6).tolist(),
            'global_identifiability_or_finite_sample_confidence_proved': False,
            'radial_topology_or_angular_tensor_identified': False,
            'absolute_energy_time_scale_identified': False,
            'physical_parameter_validation': False}


def tv(left, right):
    if [row['bits'] for row in left] != [row['bits'] for row in right]:
        raise ValueError('record alphabets must agree')
    return sum(abs(x['probability']-y['probability']) for x, y in zip(left, right, strict=True))/2


def reduced_output_distance(enlarged_states, low_states):
    """Full low quantum state jointly with each recorded bit pair.

    A thin QR evaluates the trace norm of a difference of rank-at-most17
    positive matrices; no branch is normalized or selected.
    """
    if len(enlarged_states) != len(low_states):
        raise ValueError('complete terminal record alphabets must agree')
    low_size = prod(a.LOW_DIMENSIONS)
    result = 0.
    for larger, smaller in zip(enlarged_states, low_states, strict=True):
        if np.shape(larger) != (prod(a.DIMENSIONS),) or np.shape(smaller) != (low_size,):
            raise ValueError('pure terminal branches in the declared field codes required')
        factors = np.column_stack((larger.reshape((low_size, -1), order='F'), smaller))
        _, triangular = np.linalg.qr(factors, mode='reduced')
        signs = np.ones(factors.shape[1])
        signs[-1] = -1
        matrix = (triangular*signs)@triangular.conj().T
        result += float(np.sum(abs(np.linalg.eigvalsh(matrix)))/2)
    return result


def exact_apply(lattice, vector, centered=False):
    """Sparse rational action on unnormalized Bargmann basis coefficients."""
    result = {}
    for occupations, value in vector.items():
        for term in lattice.terms:
            entries = [(occupations, value*term.coefficient)]
            for mode, matrix in zip(term.sites, term.factors, strict=True):
                following = []
                for state, coefficient in entries:
                    for row in range(len(matrix)):
                        factor = matrix[row][state[mode]]
                        if factor:
                            key = state[:mode]+(row,)+state[mode+1:]
                            following.append((key, coefficient*factor))
                entries = following
            for key, coefficient in entries:
                result[key] = result.get(key, F(0))+coefficient
    if centered:
        scalar = sum((term.coefficient*prod(matrix[0][0] for matrix in term.factors)
                      for term in lattice.terms), F(0))
        for key, value in vector.items():
            result[key] = result.get(key, F(0))-scalar*value
    return {key: value for key, value in result.items() if value}


def exact_norm_squared(vector):
    return sum((value*value*prod(factorial(n) for n in key) for key, value in vector.items()), F(0))


@lru_cache(maxsize=1)
def high_parity_instrument():
    """Source CNOTs write the two high occupations' parity into a blank bit."""
    parity = np.array([(index % 4+index//4) % 2 for index in range(16)])
    yes, no = np.diag(parity), np.diag(1-parity)
    target = np.block([[no, yes], [yes, no]])
    word = compiler.cnot(0, 4)+compiler.cnot(2, 4)
    matrix, metadata = inherited_source._execute([('high_parity_pointer', word)], 5, target)
    return (matrix[:16, :16], matrix[16:, :16]), metadata


def leakage_witness(lattice, steps=32, time=F(1, 1024), execute=True):
    """Positive odd-high occupation from the resonant interaction alone."""
    if not lattice.extended:
        raise ValueError('the leakage witness requires both new modes')
    time = a.rational(time, True)
    blank = {(0,)*9: F(1)}
    first = exact_apply(lattice, blank, centered=True)
    odd = {key: value for key, value in first.items() if (key[7]+key[8]) % 2}
    leading = exact_norm_squared(odd)
    second_squared = exact_norm_squared(exact_apply(lattice, first, centered=True))
    # sqrt(1/768) >=1/28, used as an exactly checked lower enclosure.
    root_lower = F(1, 28)*lattice.coupling
    if root_lower**2 > leading:
        raise ArithmeticError('invalid leading-amplitude lower bound')
    remainder = time*time*a.sqrt_upper(second_squared)/2
    split = splitting_bound(lattice, time, steps)['one_interval_operator_error_upper']
    ideal_margin = time*root_lower-remainder
    source_margin = ideal_margin-split
    result = {'model_time': time, 'steps': steps,
              'odd_high_first_derivative_norm_squared_exact': leading,
              'centered_second_derivative_norm_squared_exact': second_squared,
              'taylor_vector_remainder_upper': remainder,
              'source_splitting_vector_error_upper': split,
              'continuous_odd_high_probability_lower': max(F(0), ideal_margin)**2,
              'ideal_source_odd_high_probability_lower': max(F(0), source_margin)**2,
              'no_interaction_odd_high_probability_exact': '0',
              'high_read_is_separate_from_low_only_history_theorem': True}
    if execute:
        state = lattice.source_evolve(lattice.blank(), time, steps)
        reference = lattice.evolve(lattice.blank(), time)
        kraus, read_metadata = high_parity_instrument()
        probabilities = [float(np.linalg.norm(lattice.apply(state, (7, 8), matrix))**2)
                         for matrix in kraus]
        indices = np.arange(lattice.dimension)
        low = prod(a.LOW_DIMENSIONS)
        odd_mask = ((indices//low) % 4+(indices//(4*low)) % 4) % 2 == 1
        result.update({'source_odd_high_probability': probabilities[1],
                       'source_even_odd_probabilities': probabilities,
                       'source_parity_pointer': read_metadata,
                       'reference_odd_high_probability': float(np.linalg.norm(reference[odd_mask])**2)})
    return result


def checks_and_run():
    from ..common_source.records import (product_access_certificate,
                                        product_central_weight_certificate)
    if a.quartic_coefficients() != a.declared_quartic():
        raise ArithmeticError('extended quartic failed independent exact Fourier integration')
    extended, matched, bare = Lattice(), Lattice(False, True), Lattice(False, False)
    low_size = matched.dimension
    block = extended.hamiltonian[:low_size, :low_size]
    constant = a.residual_bound()['matched_scalar_increment']
    residual = block-matched.hamiltonian-float(constant)*sparse.eye(low_size)
    projected_identity_residual = float(max(abs(residual.data), default=0))
    outputs, terminal_states = {}, {}
    for name, lattice in (('extended', extended), ('matched_low', matched), ('bare_low', bare)):
        reference = [history(lattice, sign, source=False) for sign in (-1, 1)]
        actual = None
        if name != 'bare_low':
            executions = [history(lattice, sign, return_states=True) for sign in (-1, 1)]
            actual = [execution[0] for execution in executions]
            terminal_states[name] = [execution[1] for execution in executions]
        outputs[name] = {'reference_histories_by_native_sign': reference,
                         'source_histories_by_native_sign': actual,
                         'dimension': lattice.dimension,
                         'splitting_certificate': splitting_bound(lattice, INTERVAL, STEPS)}
    high_reset = [history(extended, sign, source=False, reset_high=True) for sign in (-1, 1)]
    no_feedback = [history(extended, sign, source=False, feedback=False) for sign in (-1, 1)]
    reference = outputs['extended']['reference_histories_by_native_sign']
    actual = outputs['extended']['source_histories_by_native_sign']
    low = outputs['matched_low']['source_histories_by_native_sign']
    matching = a.history_bound(2*INTERVAL)
    product_extended = outputs['extended']['splitting_certificate']['two_interval_complete_output_trace_distance_upper']
    product_low = outputs['matched_low']['splitting_certificate']['two_interval_complete_output_trace_distance_upper']
    combined = min(F(1), matching['reduced_complete_output_trace_distance_upper']+product_extended+product_low)
    witness = leakage_witness(extended)
    resource_counts = {key: 0 for key in inherited_source.COUNT_KEYS}
    maximum_local_residual = 0.
    for term in extended.terms:
        metadata = execute_term(term, INTERVAL/(2*STEPS))[1]
        maximum_local_residual = max(maximum_local_residual, metadata['full_phase_target_residual'])
        for key in resource_counts:
            resource_counts[key] += 4*STEPS*metadata['resources'][key]
    for multiplicity, metadata in ((1, inherited_source.feedback(4, .25, 2)[1]),
                                   (1, inherited_source.feedback(4, .25, 2)[1]),
                                   (2, inherited_source.instrument(4, 1, np.pi/2, 2)[1])):
        maximum_local_residual = max(maximum_local_residual, metadata['full_phase_target_residual'])
        for key in resource_counts:
            resource_counts[key] += multiplicity*metadata['resources'][key]
    native = native_records(actual)
    comparisons = {
        'source_to_same_enlarged_reference_tv': max(tv(x, y) for x, y in zip(actual, reference)),
        'source_enlarged_to_matched_low_tv': max(tv(x, y) for x, y in zip(actual, low)),
        'source_enlarged_to_matched_complete_low_quantum_record_trace_distance': max(
            reduced_output_distance(x, y) for x, y in zip(terminal_states['extended'], terminal_states['matched_low'])),
        'reference_enlarged_to_matched_low_tv': max(tv(x, y) for x, y in zip(reference, outputs['matched_low']['reference_histories_by_native_sign'])),
        'reference_enlarged_to_unmatched_bare_low_tv': max(tv(x, y) for x, y in zip(reference, outputs['bare_low']['reference_histories_by_native_sign'])),
        'trace_and_reset_high_history_tv': max(tv(x, y) for x, y in zip(reference, high_reset)),
        'feedback_deletion_history_tv': max(tv(x, y) for x, y in zip(reference, no_feedback)),
        'ideal_source_enlarged_to_matched_source_trace_distance_upper': combined,
    }
    checks = {
        'exact_quartic_integral': a.quartic_coefficients() == a.declared_quartic(),
        'vacuum_block_is_matched_low_action_plus_scalar': projected_identity_residual < 1e-12,
        'all_complete_source_laws_normalized': all(abs(sum(row['probability'] for row in law)-1) < 2e-9 for law in actual+low),
        'source_error_within_analytic_split_budget': comparisons['source_to_same_enlarged_reference_tv'] < float(product_extended)+2e-9,
        'matched_history_within_analytic_combined_budget': comparisons['source_enlarged_to_matched_low_tv'] < float(combined)+2e-9,
        'matched_complete_quantum_output_within_analytic_combined_budget': comparisons['source_enlarged_to_matched_complete_low_quantum_record_trace_distance'] < float(combined)+2e-9,
        'strictly_positive_certified_source_resonance': witness['ideal_source_odd_high_probability_lower'] > 0,
        'source_resonance_above_certified_lower_bound': witness['source_odd_high_probability'] >= float(witness['ideal_source_odd_high_probability_lower'])-1e-15,
        'continuing_high_memory_is_observable': comparisons['trace_and_reset_high_history_tv'] > 1e-10,
        'finite_history_matching_bound_is_nonvacuous': combined < 1,
    }
    if not all(checks.values()):
        raise ArithmeticError(f'band-extension checks failed: {checks}')
    product = product_interaction_control()
    product['synthetic_parameter_recovery'] = product_parameter_recovery(
        product['source_histories_by_native_sign'], product_gibbs_preparation(ProductBand())[0])
    product['native_incident_clock'] = product_incident_clock_control()
    product['homogeneous_native_clock_controller'] = product_homogeneous_controller_control()
    product['joint_fixed_reference_modular_history'] = product_modular_history_control()
    product['operational_work_access'] = product_access_certificate()
    product['central_weight_identification'] = product_central_weight_certificate()
    product['central_weight_held_out_readout'] = a.product_central_holdout_control()
    product['detector_KMS'] = a.product_detector_kms_control(product)
    product['occupation_refinement'] = a.product_occupation_refinement()
    return a.serializable({
        'schema': 'oph.interacting_band_extension.v1', 'all_pass': True,
        'scope': 'One supplied interacting Fourier-band extension, finite occupation, matched low-only histories; no ultraviolet, Lorentzian or empirical completion.',
        'source_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_PATHS},
        'modes': a.MODES, 'occupation_dimensions': a.DIMENSIONS,
        'grid_without_quartic_aliasing_minimum': 13,
        'complete_quartic_monomials': [{'powers': key, 'coefficient': value} for key, value in sorted(a.quartic_coefficients().items())],
        'model_time_per_interval': INTERVAL, 'steps_per_interval': STEPS,
        'matched_projection_matrix_residual': projected_identity_residual,
        'matching_certificate': matching, 'executed_comparisons': comparisons,
        'nine_mode_spatial_grid_certificates': [a.grid_history_bound(n, 2*INTERVAL)
                                               for n in (13, 33, 65, 129)],
        'models': outputs, 'high_reset_reference_histories': high_reset,
        'no_feedback_reference_histories': no_feedback,
        'resonant_leakage_witness': witness,
        'same_product_interaction': product,
        'complete_native_history_count': len(native),
        'complete_native_history_probability_sum': sum(row['probability'] for row in native),
        'complete_native_history_law': native,
        'worst_branch_source_resources': resource_counts,
        'additional_native_repairs_per_history': 2,
        'maximum_full_phase_local_source_residual': maximum_local_residual,
        'checks': checks,
        'acceptance': {'one_controlled_interacting_band_extension': True,
                       'quadratic_bound_for_reduced_low_quantum_and_record_output': True,
                       'high_mode_reset_or_postselection_used_in_principal_protocol': False,
                       'occupation_limit_of_enlarged_model_certified': False,
                       'ultraviolet_band_limit_or_renormalized_qft_constructed': False,
                       'finite_occupation_code_translation_covariance_claimed': False,
                       'lorentz_invariance_or_empirical_clock_identified': False}})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--output', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('--output requires a new path')
    receipt = checks_and_run()
    if args.output:
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(canonical(receipt))
    if args.check:
        compare_receipt(json.loads(args.check.read_text(encoding='utf-8')), receipt)
    display = dict(receipt['executed_comparisons'])
    key = 'ideal_source_enlarged_to_matched_source_trace_distance_upper'
    display[key] = float(F(display[key]))
    print(json.dumps({'all_pass': receipt['all_pass'], 'checks': receipt['checks'],
                      'comparison': display}, indent=2))


if __name__ == '__main__':
    main()
