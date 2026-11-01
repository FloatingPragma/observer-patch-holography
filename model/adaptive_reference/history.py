"""Complete record access and literal-source bounds for adaptive continuation."""
from __future__ import annotations

from fractions import Fraction as F
from itertools import product

import numpy as np

from ..continuing_reference import preparation
from . import controller as c, instrument


def length(value):
    if type(value) is not int or not 0 <= value <= 1_000_000:
        raise ValueError('a finite history length from zero to one million is required')
    return value


def branch(memory, records):
    """An exact complete corner history, including every zero/outside slot.

    The physical Kraus operator is matrix/sqrt(2)^len(records). The outside
    outcome is zero on this particular invariant corner; its zero is never
    extrapolated to the full Fock space.
    """
    checkpoint = c.Checkpoint(c.bits(memory))
    if type(records) not in (list, tuple) or len(records) > 128:
        raise ValueError('at most128 explicit complete records required')
    matrix = np.eye(3, dtype=object)
    for record in records:
        record = c.label(record)
        op = instrument.corner(checkpoint.memory[0])[c.LABELS.index(record)]
        matrix = op@matrix
        checkpoint = checkpoint.commit(record)
    return checkpoint, matrix, 2**len(records)


def enumerate_positive(ticks=8):
    """Enumerate every nonzero basis history; zero slots remain addressable."""
    length(ticks)
    if ticks > 12:
        raise ValueError('explicit positive-history enumeration is bounded at12 ticks')
    paths = [(memory, position, c.Checkpoint(memory), position, (), 1)
             for memory in product(range(2), repeat=2) for position in range(3)]
    for _ in range(ticks):
        updated = []
        for initial, initial_position, checkpoint, position, records, sign in paths:
            for r, op in zip(c.LABELS, instrument.corner(checkpoint.memory[0]), strict=True):
                for target in np.flatnonzero(op[:, position]):
                    updated.append((initial, initial_position, checkpoint.commit(r), int(target),
                                    records+(r,), sign*int(op[target, position])))
        paths = updated
    counts = np.zeros(12, dtype=np.int64)
    for initial, initial_position, checkpoint, position, records, sign in paths:
        if checkpoint.recover_initial_memory() != initial:
            raise ArithmeticError('complete record did not recover the initial memory')
        # Recover the field basis input from each reverse native current too.
        loads = [int(j in instrument.PAIRS[position]) for j in range(3)]
        memory = checkpoint.memory
        for record in reversed(records):
            memory = c.retreat(memory, record)
            seam = c.selector(memory)[1]
            left, right = (c.PORTS.index(p) for p in seam)
            raw = record[2] % 4
            current = {0: 0, 1: 1, 2: -1}[raw]
            loads[left] += current
            loads[right] -= current
        if tuple(j for j, n in enumerate(loads) if n) != instrument.PAIRS[initial_position]:
            raise ArithmeticError('complete history lost the initial field configuration')
        counts[3*(2*checkpoint.memory[0]+checkpoint.memory[1])+position] += 1
    if len(paths) != 12*2**ticks or np.any(counts != 2**ticks):
        raise ArithmeticError('exact complete path ensemble is not stationary')
    return {'ticks': ticks, 'positive_basis_histories': len(paths),
            'formal_record_words_per_initial_state': f'17^{ticks}',
            'all_zero_words_accessible_with_branch': True,
            'exact_output_population': [F(int(x), 12*2**ticks) for x in counts],
            'initial_memory_and_field_basis_recovered_from_final_state_and_history': True,
            'full_history_channel_forgets_initial_information': False}


def continuing_theorem():
    return {
        'actual_literal_source_success_reference': 'sigma_M=diag(a,b)^tensor M=a^M exp(gamma N), gamma=log(b/a)',
        'stationary_field_controller_success_block': 'Omega=sigma_M tensor I4/4, trace=(a+b)^M',
        'faithful_algebra': 'full proper finite-Fock field algebra tensor classical four-state controller algebra',
        'premises': ['complete selected repair: sum L*L=sum LL*=I',
                     'every selected L commutes with N',
                     'field Hamiltonian H commutes with N',
                     'two-register update (c,d)->(d XOR f(r),c)',
                     'independent exact native-balanced controller seed',
                     'all source-export abort sectors are frozen with their controller states'],
        'Hamiltonian_extension': 'separate unitary layers selected by the current working memory, using arbitrary real durations and original number-preserving H, before or after the complete repair-plus-memory update; stationary state theorem only',
        'proof': 'at output(u,v), sum_r L_r^(v) sigma_M L_r^(v)*/4=sigma_M/4',
        'arbitrary_finite_complete_histories_preserve_marginal': True,
        'traced_for_stationarity': 'seed records, successful source occupation readouts, and appended native records',
        'source_abort_mass': '1-(a+b)^M, retained outside the all-proper field block',
        'abort_archive': 'each globally aborted branch emits a distinct bypass symbol and freezes field and working memory at every tick',
        'no_success_normalization': True,
        'full_growing_history_is_stationary': False,
        'modular_generator_on_success_algebra': '-gamma N-M log(a)+log(4), controller-central scalar retained',
        'modular_automorphism': 'Ad exp(i t gamma N), identical in every classical controller sector',
        'complete_record_covariance': '(exp(i t gamma N) tensor I_memory tensor I_record) V = V (exp(i t gamma N) tensor I_memory)',
        'new_physical_energy_temperature_clock_identification': False,
        'arbitrary_feedback_policies_are_covered': False,
        'simultaneous_record_conditioned_branch_state_stationarity': False,
    }


def resource_certificate(modes=12, ticks=1000, region=4):
    ticks = length(ticks)
    source = preparation.resource_certificate(modes, region)
    return {'modes': modes, 'ticks': ticks, 'working_controller_bits': 2,
            'native_balanced_seed_events': 2,
            'field_preparation_native_source_pulses': modes,
            'field_preparation_export_records': f'6^{modes}',
            'seed_archive_formal_raw_words': 64,
            'appended_proper_repair_record_words': f'17^{ticks}',
            'global_abort_archive': 'one repeated bypass word; all initial local abort labels retained',
            'maximum_tick_record_alphabet_with_abort_bypass': 18,
            'record_storage_bits_upper': 6+5*ticks,
            'record_storage_is_unbounded_as_ticks_increases': True,
            'conditional_XOR_and_swap_updates': ticks,
            'controller_quantum_pulse_compilation_claimed': False,
            'native_seed_or_memory_rounding_error': 0,
            'literal_field_preparation_complete_trace_error': source['complete_output_trace_norm_error_to_ideal_number_reference'],
            'same_trace_error_after_any_ideal_balanced_adaptive_history': source['complete_output_trace_norm_error_to_ideal_number_reference'],
            'local_reference_log_error': source['initial_local_modular_operator_error'],
            'local_bound_uniform_in_total_modes_and_ideal_adaptive_history_length': True,
            'local_order': 'trace complete outside outputs, then project the local proper code',
            'compiled_per_tick_errors': 'add 2 sum_j epsilon_j in complete-output trace norm for complete isometry errors epsilon_j; no hardware error assumed zero',
            'finite_corner_arbitrary_input_mixing_bound': f'2*(61/64)^{ticks//8}',
            'mixing_bound_applies_to_full_export_history': False}
