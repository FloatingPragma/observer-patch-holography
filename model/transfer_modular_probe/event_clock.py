"""A native incident record triggers one fixed source-derived modular step.

The clock carrier is paused while a disjoint probe-code service executes.
Source pulses are not clock-carrier attempts. This is an active observer
coupling, not a claim that passive dissipative waiting becomes modular flow.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path
from types import MappingProxyType

import numpy as np

from . import construction as c
from ..observer import native

ROOT = c.ROOT
_SOURCE_LIBRARY = 'evidence/observer_dynamics_20260925/oph-physics-sim'
SOURCE_PATHS = tuple(sorted(set((*c.SOURCE_PATHS,
    *(f'{_SOURCE_LIBRARY}/{path}' for path in (
        'oph_exact/__init__.py', 'oph_fpe/__init__.py', 'oph_fpe/claims.py',
        'oph_fpe/constants/__init__.py', 'oph_fpe/constants/oph_pixel.py',
        'oph_fpe/core/__init__.py', 'oph_fpe/core/array_geometry.py',
        'oph_fpe/core/echosahedral_dynamics.py', 'oph_fpe/core/echosahedral_federation.py',
        'oph_fpe/core/graph.py', 'oph_fpe/core/icosahedral.py',
        'oph_fpe/core/patch_state.py', 'oph_fpe/core/patchnet.py',
        'oph_fpe/core/pixel_scale.py', 'oph_fpe/core/screen_microphysics.py',
        'oph_fpe/core/screen_ports.py', 'oph_fpe/dynamics/__init__.py',
        'oph_fpe/dynamics/anneal.py', 'oph_fpe/dynamics/consensus_certificate.py',
        'oph_fpe/dynamics/kernel_dispatcher.py', 'oph_fpe/dynamics/modular.py',
        'oph_fpe/dynamics/positive_geometry.py', 'oph_fpe/dynamics/repair.py',
        'oph_fpe/evidence/__init__.py', 'oph_fpe/evidence/bundle.py',
        'oph_fpe/evidence/cross_repo_artifacts.py', 'oph_fpe/evidence/hashes.py',
        'oph_fpe/evidence/particle_input_policy.py', 'oph_fpe/evidence/production_envelope.py',
        'oph_fpe/evidence/verifier.py', 'oph_fpe/finite_groups.py',
        'oph_fpe/gauge/__init__.py', 'oph_fpe/gauge/a5_sm_certificate.py',
        'oph_fpe/gauge/covariant_overlap.py', 'oph_fpe/gauge/higgs_carrier.py',
        'oph_fpe/gauge/mar_sieve.py', 'oph_fpe/gauge/physical_a5_sm_requirements.py',
        'oph_fpe/gauge/repair_projection.py', 'oph_fpe/gauge/source_to_eft.py',
        'oph_fpe/gauge/yang_mills_gap.py', 'oph_fpe/groups/__init__.py',
        'oph_fpe/groups/base.py', 'oph_fpe/groups/clock.py',
        'oph_fpe/groups/registry.py', 'oph_fpe/groups/s3.py', 'oph_fpe/groups/z2.py')),
    'model/transfer_modular_probe/event_clock.py',
    'model/__init__.py', 'model/anharmonic_transfer/__init__.py',
    'model/common_source/__init__.py', 'model/common_source/composition.py',
    'model/common_source/native_join.py', 'model/common_source/operational.py',
    'model/common_source/records.py', 'model/common_source/thermal_clock.py',
    'model/native_clock/__init__.py', 'model/observer/__init__.py',
    'model/observer/p_response/__init__.py', 'model/observer/p_response/experiment.py',
    'model/observer/scalar_response/__init__.py',
    'model/observer/source_discrimination/__init__.py',
    'model/observer/source_scalar/__init__.py', 'model/observer/source_scalar/experiment.py',
    'model/predictions/__init__.py', 'model/predictions/observer_instrument.py',
    'model/predictions/protocol.py', 'model/primitive_observer/__init__.py',
    'model/relativistic_scalar/__init__.py', 'model/relativistic_scalar/lattice.py',
    'model/scalar_observations/__init__.py',
    'model/scalar_observations/experimental_readout/__init__.py',
    'model/scalar_observations/lattice.py', 'model/scalar_observations/source.py'))))
MAX_GATES = 8
COMPLETE_ALLOWANCE = F(1, 10**10)


def _cap(value):
    if type(value) is not int or value not in (1, 2):
        raise ValueError('the native-attempt service cap is one or two')
    return value


@lru_cache(maxsize=1)
def native_tables():
    """Actual branch counts and complete native-C partial-permutation checks."""
    source = c.transfer.native_counts()
    c.transfer.exact_audit(source)
    index = {state: i for i, state in enumerate(source.states)}
    marks = tuple((seam, current) for seam in source.seams if 0 in seam for current in (-1, 0, 1))
    matrices = {mark: np.zeros((66, 66), dtype=np.int64) for mark in marks}
    remote = np.zeros((66, 66), dtype=np.int64)
    for seam in source.seams:
        groups = {}
        for coin in (0, 1):
            for i, state in enumerate(source.states):
                after, (_, current) = native.native(state, seam, coin)
                j = index[after]
                groups.setdefault((coin, current), {})[i] = j
                if 0 in seam:
                    matrices[seam, current][i, j] += 1
                else:
                    remote[i, j] += 1
        # With amplitude 1/sqrt2 each recorded native-C Kraus is a partial
        # permutation, and the complete seam effect is exactly the identity.
        if (any(len(set(group.values())) != len(group) for group in groups.values())
                or any(sum(i in group for group in groups.values()) != 2 for i in range(66))):
            raise ArithmeticError('complete native-C instrument structure failed')
    if (not np.all(remote.sum(axis=1) == 50)
            or not np.all(sum(matrices.values()).sum(axis=1) == 10)):
        raise ArithmeticError('the clock-port exposure is not one sixth')
    for matrix in (remote, *matrices.values()):
        matrix.flags.writeable = False
    return source, remote, MappingProxyType(matrices)


def stopped_tables(cap=2):
    """Exact integer completed local-record maps and the unrenormalized timeout."""
    cap = _cap(cap)
    source, remote, marks = native_tables()
    result = {mark: np.zeros((66, 66), dtype=np.int64) for mark in marks}
    power = np.eye(66, dtype=np.int64)
    for attempt in range(cap):
        for mark, matrix in marks.items():
            result[mark] += 60**(cap-attempt-1)*(power@matrix)
        power = power@remote
    if not np.all((power+sum(result.values())).sum(axis=1) == 60**cap):
        raise ArithmeticError('a stopped interval lost an outcome or timeout')
    return source, result, power, 60**cap


def endpoint_oracle(cap=2):
    """Independent floor/ceiling path recursion, without native() or resolvents."""
    cap = _cap(cap)
    source = c.transfer.native_counts()
    index = {state: i for i, state in enumerate(source.states)}
    result = {(seam, current): np.zeros((66, 66), dtype=np.int64)
              for seam in source.seams if 0 in seam for current in (-1, 0, 1)}
    timeout = np.zeros((66, 66), dtype=np.int64)
    def visit(origin, state, remaining):
        if not remaining:
            timeout[origin, index[state]] += 1
            return
        for seam in source.seams:
            a, b = seam
            total = state[a]+state[b]
            for coin in (0, 1):
                after = list(state)
                after[a], after[b] = ((total//2, (total+1)//2) if coin == 0
                                      else ((total+1)//2, total//2))
                after = tuple(after)
                if 0 in seam:
                    result[seam, state[a]-after[a]][origin, index[after]] += 60**(remaining-1)
                else:
                    visit(origin, after, remaining-1)
    for i, state in enumerate(source.states):
        visit(i, state, cap)
    return result, timeout


@lru_cache(maxsize=2)
def _gates(mode):
    mode = c._mode(mode)
    word = c.balanced_word(mode)
    b, br = c._source(word, 'balanced', j=mode)
    undo, ur = c._source(c.compiler.inverse(word), 'balanced', j=mode)
    advance, ar = c._source(c.diagonal_word(1), 'phase', s=1)
    # Reuse validated shared instruments, including complete four-column maps.
    return tuple(c._unitary(matrix, resources, name) for matrix, resources, name in
        ((b, br, 'balanced'), (undo, ur, 'inverse_balanced'), (advance, ar, 'one_modular_tick')))


def triggered_word(record):
    """The controller reads only seam/current; no elapsed-time/phase argument."""
    _, _, marks = native_tables()
    if (type(record) is not tuple or len(record) != 2 or type(record[0]) is not tuple
            or type(record[1]) is not int or record not in marks):
        raise ValueError('an actual port-zero incident seam/current record is required')
    return c.diagonal_word(1)


def probe_outputs(mode, *, current=0, feedback=True, ideal=False, reset=False):
    """Complete continuing probe factors, including all padding and backaction."""
    mode = c._mode(mode)
    if (type(current) is not int or current not in (-1, 0, 1)
            or any(type(x) is not bool for x in (feedback, ideal, reset))):
        raise ValueError('an exact current and explicit protocol/control flags are required')
    balanced, undo, advance = _gates(mode)
    b, inverse, u = (x.kraus[0][0] for x in (balanced, undo, advance))
    if ideal:
        b = np.eye(4, dtype=complex)
        b[0, 0] = b[0, mode] = b[mode, 0] = 1/math.sqrt(2)
        b[mode, mode] = -1/math.sqrt(2)
        inverse = b
        u = np.diag([complex(math.cos(float(sum(e)/2)), math.sin(float(sum(e)/2)))
                     for e in c.energies()]+[1])
    initial = b@np.array([1, 0, 0, 0], dtype=complex)
    first, second = {}, {}
    analyzed = b@u@initial
    for a in range(4):
        v = np.zeros(4, dtype=complex); v[a] = analyzed[a]
        v = inverse@v
        if feedback and current > 0:
            v = b@v  # Fixed noncommuting intervention, triggered by native current.
        if reset:
            # A countercontrol only: a fresh preparation with the old weight.
            v = initial*np.linalg.norm(v)
        first[a] = v
        measured = b@u@v
        for d in range(4):
            out = np.zeros(4, dtype=complex); out[d] = measured[d]
            second[a, d] = inverse@out
    return initial, first, second


def density_oracle(mode, *, current=0, feedback=True):
    """Independent ideal density evolution from the closed spectral roots."""
    mode = c._mode(mode)
    if type(current) is not int or current not in (-1, 0, 1) or type(feedback) is not bool:
        raise ValueError('a declared native current and explicit feedback flag are required')
    b = np.eye(4, dtype=complex)
    for a, d in ((0, 0), (0, mode), (mode, 0)):
        b[a, d] = math.sqrt(.5)
    b[mode, mode] = -math.sqrt(.5)
    values = (1., (25+math.sqrt(3))/30, (25-math.sqrt(3))/30, 1.)
    u = np.diag(np.exp(-1j*np.log(values)))
    blank = np.zeros((4, 4), complex); blank[0, 0] = 1
    initial = b@blank@b.conj().T
    before = b@u@initial@u.conj().T@b.conj().T
    first, second = {}, {}
    for a in range(4):
        projector = np.diag([int(i == a) for i in range(4)])
        rho = b@projector@before@projector@b.conj().T
        if feedback and current > 0:
            rho = b@rho@b.conj().T
        first[a] = rho
        measured = b@u@rho@u.conj().T@b.conj().T
        for d in range(4):
            pointer = np.diag([int(i == d) for i in range(4)])
            second[a, d] = b@pointer@measured@pointer@b.conj().T
    return initial, first, second


def complete_output(mode, cap=2, *, feedback=True, ideal=False, reset=False, oracle=False):
    """CQ execution witness: exact native terminal populations and raw probe factors.

    The native witness is stationary diagonal. The composition error theorem
    additionally covers coherent native inputs and arbitrary references because
    the native Kraus instruments and all four probe columns are controlled.
    """
    cap = _cap(cap)
    if type(oracle) is not bool:
        raise ValueError('an explicit oracle flag is required')
    source, maps, timeout, denominator = stopped_tables(cap)
    if oracle:
        maps, timeout = endpoint_oracle(cap)
    initial, _, _ = probe_outputs(mode, feedback=feedback, ideal=ideal)
    prior = np.ones(66, dtype=np.int64)
    # Each output is (record, exact integer native terminal weights,
    # denominator, unnormalized original-binary probe vector).
    output = [(('timeout',), prior@timeout, 66*denominator, initial)]
    for mark, matrix in maps.items():
        first_native = prior@matrix
        _, first, second = probe_outputs(mode, current=mark[1], feedback=feedback,
                                         ideal=ideal, reset=reset)
        for read in range(4):
            output.append(((mark, read, 'timeout'), first_native@timeout,
                           66*denominator**2, first[read]))
            for following, next_matrix in maps.items():
                native_state = first_native@next_matrix
                for final in range(4):
                    output.append(((mark, read, following, final), native_state,
                                   66*denominator**2, second[read, final]))
    return output


def _probability(branch):
    _, native_state, denominator, vector = branch
    return float(F(int(native_state.sum()), denominator))*float(np.vdot(vector, vector).real)


def _joint_diagnostic(actual, ideal):
    if len(actual) != len(ideal):
        raise ArithmeticError('a complete output branch was lost')
    distance = 0.
    for a, b in zip(actual, ideal, strict=True):
        if a[0] != b[0] or a[2] != b[2] or not np.array_equal(a[1], b[1]):
            raise ArithmeticError('native endpoint and stopped-transfer joint records differ')
        difference = np.outer(a[3], a[3].conj())-np.outer(b[3], b[3].conj())
        distance += float(F(int(a[1].sum()), a[2]))*np.abs(np.linalg.eigvalsh(difference)).sum()/2
    return float(distance)


def density_oracle_diagnostic(actual, mode):
    oracles = {current: density_oracle(mode, current=current) for current in (-1, 0, 1)}
    total = 0.
    for record, weights, denominator, vector in actual:
        if len(record) == 1:
            rho = oracles[0][0]
        elif len(record) == 3:
            rho = oracles[record[0][1]][1][record[1]]
        else:
            rho = oracles[record[0][1]][2][record[1], record[3]]
        delta = np.outer(vector, vector.conj())-rho
        total += float(F(int(weights.sum()), denominator))*np.abs(np.linalg.eigvalsh(delta)).sum()/2
    return float(total)


def multiplicative_independence():
    """Exact ingredients of the conjugation proof in README.md."""
    c.spectrum_audit()
    norm = F(25*25-3, 30*30)
    if norm != F(311, 450) or not 0 < norm < 1:
        raise ArithmeticError('the conjugate eigenvalue norm changed')
    return {'quadratic_field': 'Q(sqrt(3))', 'conjugate_product': str(norm),
        'conjugate_ratio_strictly_greater_than_one': True,
        'integer_relation_lambda_plus_power_m_times_lambda_minus_power_n_equals_one': 'm=n=0',
        'nonzero_log_frequency_ratio_is_irrational': True,
        'joint_phase_pair_has_no_exact_nonzero_real_parameter_alias': True,
        'single_measurement_determines_elapsed_parameter': False}


def _round_interval(interval, bits=90):
    scale = 2**bits
    lo, hi = interval
    return F((lo*scale).numerator//(lo*scale).denominator, scale), F(-((-hi*scale).numerator//(-hi*scale).denominator), scale)


def _product(a, b):
    values = [x*y for x in a for y in b]
    return min(values), max(values)


def _complex_product(a, b):
    ac, bd = _product(a[0], b[0]), _product(a[1], b[1])
    ad, bc = _product(a[0], b[1]), _product(a[1], b[0])
    return (_round_interval((ac[0]-bd[1], ac[1]-bd[0])),
            _round_interval((ad[0]+bc[0], ad[1]+bc[1])))


def phase_power(ticks, mode):
    """Exact outward dyadic interval powering; long windows need no float log."""
    if type(ticks) is not int or not 0 <= ticks <= 2**20:
        raise ValueError('an integer clock reading from zero through 2^20 is required')
    mode = c._mode(mode)
    base = c.phase_targets(1)[mode]
    result = ((F(1), F(1)), (F(0), F(0)))
    while ticks:
        if ticks & 1:
            result = _complex_product(result, base)
        base = _complex_product(base, base)
        ticks //= 2
    return result


def chord_interval(ticks, mode):
    real, imag = phase_power(ticks, mode)
    real = (real[0]-1, real[1]-1)
    lower = upper = F(0)
    for lo, hi in (real, imag):
        lower += 0 if lo <= 0 <= hi else min(lo*lo, hi*hi)
        upper += max(lo*lo, hi*hi)
    return c.exact._sqrt_interval(lower)[0], c.exact._sqrt_interval(upper)[1]


def resolution_audit():
    # The existing complete Y analysis is certified as well. Repeating the
    # same U(1), rather than compiling a new chosen phase, uses at most n+5
    # source gates including preparation and both analysis undo operations.
    for mode in (1, 2):
        c.program(mode, (1,), quadrature='y')
    if c._sequence_allowance(69) >= F(1, 10**9):
        raise ArithmeticError('the fixed-step window source bound failed')
    separation = min(max(chord_interval(n, j)[0] for j in (1, 2)) for n in range(1, 65))
    if separation <= F(1, 4):
        raise ArithmeticError('the declared 64-tick joint phase separation failed')
    # Find an approximate recurrence numerically, then certify that one integer
    # independently using exact interval arithmetic. It is not an exact alias.
    energies = [float(sum(x)/2) for x in c.energies()[1:]]
    recurrence = min(range(65, 4097), key=lambda n:
        max(abs(complex(math.cos(n*e), math.sin(n*e))-1) for e in energies))
    recurrence_upper = max(chord_interval(recurrence, j)[1] for j in (1, 2))
    outward = F(math.ceil(float(recurrence_upper)*10**6)+1, 10**6)
    if outward >= F(1, 8) or recurrence_upper > outward:
        raise ArithmeticError('the finite-resolution recurrence control failed')
    return {'tick_window': [0, 64], 'joint_complex_phase_separation_lower': '1/4',
        'readout': 'two quadrature ensemble means for each mode; no one-shot exact reading',
        'four_ideal_probability_mean_errors_at_most': '1/32',
        'independent_shots_per_quadrature_sufficient_failure_bound': 'min(1,1024/shots)',
        'finite_shot_bound_requires_repeated_known_preparations': True,
        'fixed_U1_window_maximum_source_gates_per_quadrature': 69,
        'fixed_U1_window_complete_source_half_diamond_bound': '1/1000000000',
        'certified_approximate_recurrence_ticks': recurrence,
        'approximate_recurrence_max_chord_upper': str(outward),
        'pigeonhole_recurrence': 'for integer Q>=1, some 1<=n<=Q^2 has both phase chords<=2pi/Q',
        'unbounded_time_perfect_resolution_claimed': False}


def run():
    if c._sequence_allowance(MAX_GATES) > COMPLETE_ALLOWANCE:
        raise ArithmeticError('complete active-clock source allowance is insufficient')
    rows = []
    for cap in (1, 2):
        _, maps, timeout, denominator = stopped_tables(cap)
        direct, direct_timeout = endpoint_oracle(cap)
        if not np.array_equal(timeout, direct_timeout) or any(
                not np.array_equal(matrix, direct[mark]) for mark, matrix in maps.items()):
            raise ArithmeticError('stopped maps disagree with independent native endpoint paths')
        for mode in (1, 2):
            actual = complete_output(mode, cap)
            ideal = complete_output(mode, cap, ideal=True)
            distance = _joint_diagnostic(actual, ideal)
            independent = density_oracle_diagnostic(actual, mode)
            total = sum(map(_probability, actual))
            if max(distance, independent, abs(total-1)) > 1e-11:
                raise ArithmeticError('the complete controller-driven source output failed its diagnostic')
            digest = hashlib.sha256()
            for record, weights, den, _ in actual:
                digest.update(c.canonical([record, [int(x) for x in weights], den]).encode())
            rows.append({'cap_per_incident_interval': cap, 'mode': mode,
                'complete_joint_records': len(actual), 'native_terminal_record_sha256': digest.hexdigest(),
                'first_interval_timeout_probability_exact': str(F(5, 6)**cap),
                'both_ticks_probability_exact': str((1-F(5, 6)**cap)**2),
                'raw_output_total_weight_diagnostic': total,
                'complete_classical_native_quantum_probe_half_trace_diagnostic': distance,
                'independent_density_oracle_half_trace_diagnostic': independent,
                'complete_all_input_reference_source_half_diamond_bound': str(COMPLETE_ALLOWANCE)})
    ordinary = complete_output(2, feedback=False, ideal=True)
    reset = complete_output(2, feedback=False, ideal=True, reset=True)
    event = lambda branch: len(branch[0]) == 4 and branch[0][-1] == 0
    no_reset_difference = sum(_probability(b) for b in reset if event(b))-sum(_probability(b) for b in ordinary if event(b))
    real = c.phase_targets(1)[2][0]
    qlo, qhi = (1+real[0])/2, (1+real[1])/2
    lower = F(11, 36)**2*(1-qhi)*(2*qlo-1)
    if lower <= F(1, 1000) or no_reset_difference < float(lower)-1e-14:
        raise ArithmeticError('the continuing-clock reset control lost its certified contrast')
    return {'schema': 'oph.transfer_modular_probe.event_clock.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'contract': {'clock_event': 'native port-zero incident attempt, including zero current',
            'phase_advance_per_event': 'one fixed derived U(1), with no variable phase input',
            'carrier_paused_during_disjoint_probe_source_services': True,
            'service_cap_counter_controls_timeout_not_phase': True,
            'source_preparation_and_coherent_probe_services_supplied': True,
            'current_positive_feedback': 'the same compiled balanced intervention between reads',
            'native_execution_witness': 'stationary diagonal 66-state carrier; complete terminal populations',
            'error_theorem': 'complete native C quantum instrument and arbitrary untouched references',
            'native_clock_events_and_probe_source_services_are_distinct': True},
        'complete_history_checks': rows,
        'source_resources': {str(j): [dict(g.resources) for g in _gates(j)] for j in (1, 2)},
        'maximum_source_resources_per_mode': {str(j): dict(c.runtime.add_resources(
            *(_gates(j)[0].resources,)*4, *(_gates(j)[1].resources,)*2,
            *(_gates(j)[2].resources,)*2)) for j in (1, 2)},
        'maximum_native_attempts_and_pointer_reads': [4, 2],
        'no_reset_control': {'complete_success_and_second_zero_event_difference_lower': '1/1000',
            'ideal_difference_diagnostic': no_reset_difference,
            'two_raw_compiled_protocols_event_weight_error_upper': str(4*COMPLETE_ALLOWANCE)},
        'frequency_independence': multiplicative_independence(),
        'finite_resolution': resolution_audit(),
        'acceptance': {'native_records_drive_the_same_derived_modular_unitary': True,
            'complete_sequential_records_interventions_and_timeouts_retained': True,
            'source_gate_error_includes_arbitrary_references_and_raw_norm_gain': True,
            'full_coherent_native_input_trajectory_executed': False,
            'passive_native_repair_waiting_is_modular_flow': False,
            'coupling_or_reference_uniquely_selected_by_OPH': False,
            'physical_service_error_or_observed_time_identification_established': False}}


def main():
    parser = argparse.ArgumentParser(__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--check', type=Path)
    modes.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        c.compare_receipt(result, c.strict_json.load(args.check))
        print('Native-record-driven modular clock, complete histories and resolution bounds replayed.')
    elif args.output:
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(c.canonical(result))
    else:
        print(c.canonical(result), end='')


if __name__ == '__main__':
    main()
