"""Source-executed ordered-cover feedback; all native outputs are retained."""
from __future__ import annotations

import math
import numpy as np

from ..common_source import native_routing as routing, record_code
from ..exchange_sector import statistics as exchange, statistics_service as source
from ..observer.source_scalar import compiler
from . import protection


SIZE = 1024
LIMIT = 2e-7


def controlled_parity_quarters(seam, current):
    """Controlled Z is nontrivial only at negative ordered-lift columns."""
    return (0,)*512 + source.phase_quarters(seam, current, -1)


def memory_word(seam, current):
    """Y(-pi/4) controlled-Z Y(pi/4) = controlled-X, exactly."""
    phases = controlled_parity_quarters(seam, current)
    if not any(phases):
        return []
    return ([compiler.pulse('Y', math.pi/4, [9])]
            + source.diagonal_word(phases)
            + [compiler.pulse('Y', -math.pi/4, [9])])


def apply_memory(seam, current, state):
    return compiler.execute(memory_word(seam, current), 10, source._state(state, 10))


def repair_branches(law, seam, state):
    """Execute native C/E, then its record-selected ordering update.

    Source routing treats the memory as an untouched reference.  Feedback
    then acts on the full ten-bit register, with an identity padding and
    collision extension.  No branch normalization or output trimming occurs.
    """
    state = source._state(state, 10)
    vector, width = state.ndim == 1, 1 if state.ndim == 1 else state.shape[1]
    columns = state.reshape(2, 512, width).transpose(1, 0, 2).reshape(512, 2*width)
    rows = source._bare_repair(law, seam, columns)
    for row in rows:
        value = row['amplitude'].reshape(512, 2, width).transpose(1, 0, 2).reshape(1024, width)
        if row['current'] in (-1, 1):
            value = apply_memory(seam, row['current'], value)
        row['amplitude'] = value[:, 0] if vector else value
    source._validate_grid(law, rows, shape=state.shape)
    return rows


def source_basis():
    value = np.zeros((SIZE, protection.COVER_DIMENSION), complex)
    for r in (0, 1):
        for i, pair in enumerate(exchange.PAIRS):
            value[r*512+source._address(pair), r*66+i] = 1
    return value


def exact_target(row):
    value = np.zeros((SIZE, protection.COVER_DIMENSION), complex)
    for a, b in protection.lifted_pairs(row):
        r, i = divmod(b, 66)
        value[r*512+source._address(exchange.PAIRS[i]), a] = 1/math.sqrt(2)
    return value


def _resource_rows(law, seam):
    base = record_code.repair_resources(law, seam)
    result = []
    for row in protection.branches(law):
        if row.seam != seam:
            continue
        labels = source.circuit.record_values(law, row.coin)
        raw_value = labels.get(row.raw)
        current = raw_value if law == 'C' or raw_value is None else (
            raw_value-(2*row.coin-1)*(abs(raw_value) % 2))//2
        word = memory_word(seam, current) if current in (-1, 1) else []
        feedback = source._counts(word)
        result.append({'coin': row.coin, 'raw': row.raw, 'abort': row.abort,
                       **source._failure_account(source._add(base, feedback)),
                       'orientation_feedback_pulses': len(word)})
    return result


def audit_feedback():
    """Join all exact compiled phases to every native transported-order column."""
    # Exact full-address audit of every feedback word, including scalar
    # terms and restored parity wires.  Also compare the resulting controlled
    # phase to each complete native branch, rather than relying on the word
    # producer to be its own parity oracle.
    compiled = {}
    addresses = {source._address(pair) for pair in exchange.PAIRS}
    for seam in exchange.SEAMS:
        for current in (-1, 0, 1):
            phases = controlled_parity_quarters(seam, current)
            if any(phases[:512]) or any(phases[512+i] for i in range(512) if i not in addresses):
                raise ArithmeticError('feedback changed its identity padding or control extension')
            source._audit_diagonal_word(phases, source.diagonal_word(phases))
            compiled[seam, current] = phases
    checked = 0
    for law in ('C', 'E'):
        for row in protection.branches(law):
            for i, j, sign in zip(row.inputs, row.outputs, row.signs, strict=True):
                current = exchange.STATES[i][row.seam[0]]-exchange.STATES[j][row.seam[0]]
                phase = compiled[row.seam, current][512+source._address(exchange.PAIRS[j])]
                if phase != 2*int(sign == -1):
                    raise ArithmeticError('compiled feedback does not transport native token order')
                checked += 1
    return checked


def certificate():
    # Numerical source execution independently checks all132 columns on one
    # witness seam; the exact audit covers every native seam and raw column.
    checked = audit_feedback()
    seam = exchange.WITNESS_SEAM
    expected_basis = source_basis()
    results = {}
    for law in ('C', 'E'):
        rows = repair_branches(law, seam, expected_basis)
        exact = {(row.coin, row.raw): row for row in protection.branches(law) if row.seam == seam}
        squared_error = sum(float(np.linalg.norm(
            row['amplitude']-exact_target(exact[row['coin'], row['record_code']]))**2) for row in rows)
        if math.sqrt(squared_error) > LIMIT:
            raise ArithmeticError('complete source-executed cover instrument failed its132-column check')
        # Frobenius error bounds operator error, hence all unknown source,
        # memory and untouched-reference inputs.  We retain a fixed admission
        # ceiling, not platform-sensitive residual decimal digits.
        results[law] = {'raw_slots': len(rows), 'abort_slots': sum(row['abort'] for row in rows),
                        'terminal_addresses_per_slot': SIZE, 'all_joint_input_columns': 132,
                        'combined_Stinespring_Frobenius_error_below': str(LIMIT),
                        'complete_raw_output_and_padding_checked': True,
                        'native_and_feedback_resource_rows': _resource_rows(law, seam)}
    return {'padded_source_bits': 9, 'persistent_orientation_bits': 1,
            'witness_seam': list(seam), 'all_seam_current_feedback_words_audited': 90,
            'native_source_columns_joined_to_exact_feedback': checked,
            'addresses_per_exact_diagonal_word_audit': SIZE,
            'native_raw_records_kept': True, 'proper_code_failure_complements_kept': True,
            'ideal_proper_code_failure_probability': '0',
            'sharp_memory_preparation_from_zero': 'Y(-pi/4) for +; Y(pi/4) for -',
            'sharp_memory_readout': 'Y(pi/4), then both computational outcomes',
            'preparation_and_readout_pulses_each': 1,
            'fresh_coin_bits_per_attempt': 1, 'fresh_record_bits': {'C': 2, 'E': 3},
            'growing_classical_tape_bits_per_attempt': {'C': 3, 'E': 4},
            'initial_source_state_for_column_certificate': 'arbitrary supplied joint source/memory input',
            'physical_preparation_storage_energy_costs_calibrated': False,
            'floating_execution_is_not_interval_or_hardware_proof': True,
            'laws': results}
