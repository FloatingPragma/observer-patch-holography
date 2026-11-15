"""Actual C10 gauge characters and a source-compiled phase instrument.

This tests an alternative topology/current map. Gauge-orbit winding is not
identified with a spatial defect degree. Proper-code coherent processor
services implement the apparatus; passive repair does not synthesize them.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[3]
WORKSPACE = BASE.parent
RER = WORKSPACE / 'reverse-engineering-reality'
sys.path[:0] = [str(BASE / 'next/C10'), str(RER)]
from joint_source import JointNativeSource
from model.observer.source_scalar import compiler


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def degree(samples, *, minimum_modulus=1e-9):
    """Polygonal winding, admitted only for a closed, nonzero resolved loop."""
    samples = np.asarray(samples, complex)
    require(samples.ndim == 1 and len(samples) > 2, 'loop samples required')
    require(np.max(abs(samples)) < math.inf, 'finite samples required')
    require(abs(samples[0] - samples[-1]) < 1e-10, 'loop is not closed')
    require(np.min(abs(samples)) > minimum_modulus, 'loop meets zero')
    increments = np.angle(samples[1:] / samples[:-1])
    require(np.max(abs(increments)) < math.pi / 2, 'phase sampling is unresolved')
    value = float(sum(increments) / (2 * math.pi))
    require(abs(value - round(value)) < 1e-9, 'winding is not integral')
    return int(round(value))


def phase_word(theta):
    """Controlled exp(i theta(N-3)), with the pointer on bit six.

    Six controlled occupation phases and one reference-offset phase use
    existing proper-code source words. The six record-fibre labels tensor
    identity. The compound operation preserves the C10 encoded state space;
    individual processor pulses are not C10 native repair events.
    """
    word = compiler.remap(compiler.compile_diagonal([0., -3 * theta]), [6])
    for bit in range(6):
        local = compiler.compile_diagonal([0., 0., 0., theta])
        word += compiler.remap(local, [bit, 6])
    return word


def source_meter(theta, axis):
    require(axis in ('X', 'Y'), 'X or Y pointer axis required')
    prepare = [compiler.pulse('H', math.pi / 2, [6])]
    read = [] if axis == 'X' else compiler.remap(
        compiler.compile_diagonal([0., -math.pi / 2]), [6])
    word = prepare + phase_word(theta) + read + compiler.inverse(prepare)
    actual = compiler.execute(word, 7)
    # Both pointer outcomes are retained on every occupation input column.
    operators = [actual[64 * bit:64 * (bit + 1), :64] for bit in (0, 1)]
    offdiagonal = max(float(np.max(abs(k - np.diag(np.diag(k))))) for k in operators)
    completeness = float(np.max(abs(sum(k.conj().T @ k for k in operators) - np.eye(64))))
    return operators, word, offdiagonal, completeness


def evaluate():
    contract_path = BASE / 'next/C10/CANONICAL_SOURCE_CONTRACT.json'
    contract = json.loads(contract_path.read_text())
    for relative, digest in contract['source_sha256'].items():
        require(sha(BASE / relative) == digest, 'canonical source drift: ' + relative)
    joint_receipt = json.loads((BASE / 'next/C10/JOINT_SOURCE.json').read_text())
    for relative, digest in joint_receipt['source_sha256'].items():
        require(sha(BASE / relative) == digest, 'joint source drift: ' + relative)
    model = JointNativeSource()
    archive_path = BASE / 'next/C10/joint_source_states.npz'
    archive = np.load(archive_path)
    require(np.array_equal(archive['hamiltonian'], model.hamiltonian), 'actual H changed')
    masks = np.array([mask for mask, _ in model.states], int)
    numbers = np.array([int(mask).bit_count() for mask in masks])
    charge = numbers - 3
    require(np.array_equal(charge, -archive['boundary_B']), 'boundary charge differs')
    h_error = float(np.max(abs((charge[:, None] - charge[None, :]) * model.hamiltonian)))
    native_error = 0.
    native_operators = 0
    for seams in model.field.operators_by_law.values():
        for operators in seams:
            for _, k in operators:
                native_error = max(native_error, float(np.max(abs((charge[:, None] - charge[None, :]) * k))))
                native_operators += 1
    reservoir_error = 0.
    export_error = 0.
    export_failure_ranks = {}
    for dimension in (2, 4):
        export = compiler.transfer_kraus(dimension)
        export_error = max(export_error, float(np.max(abs(
            sum(k.conj().T @ k for k in export) - np.eye(6)))))
        export_failure_ranks[str(dimension)] = len(export) - 1
    branches = []
    for meta, k in model.source.reservoir_instrument():
        delta = meta['reservoir_initial'] - meta['reservoir_final']
        reservoir_error = max(reservoir_error, float(np.max(abs((charge[:, None] - charge[None, :] - delta) * k))))
        key = f"global_allN_reservoir_{meta['reservoir_initial']}_{meta['reservoir_final']}"
        state = archive[key]
        p = float(np.trace(state).real)
        require(np.max(abs((charge[:, None] - delta) * state)) < 1e-12, 'branch not sharp in actual charge')
        branches.append({**meta, 'charge_relative_to_reference': delta,
                         'archive_key': key, 'unconditional_probability': p})

    # The derivative relation follows from the CAR action on each actual
    # occupation mask. The generator is obtained by summing occupied modes;
    # the reference fixes the otherwise undetectable central offset.
    car_edges = 0
    for mask in range(64):
        for mode in range(6):
            if not (mask >> mode) & 1:
                require((mask | (1 << mode)).bit_count() - mask.bit_count() == 1,
                        'creator gauge weight')
                car_edges += 1

    angles = np.linspace(0., 2 * math.pi, 33)
    expected_degrees = {q: degree(np.exp(1j * q * angles)) for q in range(-3, 4)}
    reconstructed = {q: np.zeros(len(angles), complex) for q in range(-3, 4)}
    complete_records = []
    max_off = max_complete = max_probability_error = 0.
    example_words = {}
    source_q = np.array([mask.bit_count() - 3 for mask in range(64)])
    for angle_index, theta in enumerate(angles):
        for axis in ('X', 'Y'):
            ops, word, off, complete = source_meter(float(theta), axis)
            max_off = max(max_off, off)
            max_complete = max(max_complete, complete)
            diagonals = [np.diag(k) for k in ops]
            probabilities = [abs(d) ** 2 for d in diagonals]
            target = np.cos(source_q * theta) if axis == 'X' else np.sin(source_q * theta)
            max_probability_error = max(max_probability_error,
                float(np.max(abs(probabilities[0] - (1 + target) / 2))),
                float(np.max(abs(probabilities[1] - (1 - target) / 2))))
            for q in range(-3, 4):
                ids = np.flatnonzero(source_q == q)
                contrast = float(np.mean(probabilities[0][ids] - probabilities[1][ids]))
                reconstructed[q][angle_index] += contrast * (1 if axis == 'X' else 1j)
            if angle_index == 3:
                example_words[axis] = {'theta': float(theta), 'word': word,
                                       'resources': compiler.resources(word)}
            for branch in branches:
                state = archive[branch['archive_key']]
                for bit, amplitude in enumerate(diagonals):
                    # The actual C10 record fibre is unchanged. This diagonal
                    # operator gives the full subnormalized output state.
                    a = amplitude[masks]
                    post = a[:, None] * state * a.conj()[None, :]
                    p = float(np.trace(post).real) / 2
                    complete_records.append({'angle_index': angle_index, 'axis': axis,
                        'pointer': bit, 'reservoir_initial': branch['reservoir_initial'],
                        'reservoir_final': branch['reservoir_final'],
                        'joint_probability_including_fair_axis': p,
                        'poststate_definition': 'diag(a_axis,bit(mask)) rho_reservoir diag(a)* / 2'})
    actual_degrees = {q: degree(values) for q, values in reconstructed.items()}
    sums = [sum(row['joint_probability_including_fair_axis'] for row in complete_records
                if row['angle_index'] == j) for j in range(len(angles))]
    rejected = {}
    for name, curve in {'half_charge_nonclosed': np.exp(.5j * angles),
                        'equal_zero_one_mixture_has_zero': (1 + np.exp(1j * angles)) / 2}.items():
        try:
            degree(curve)
        except ValueError as error:
            rejected[name] = str(error)
        else:
            raise AssertionError('inadmissible loop accepted: ' + name)
    # Joint input of charge one and a neutral reference has sharp total one;
    # equal expected daughter charge can be produced without half eigenvalues.
    split_control = {'state': '(|1,0>+|0,1>)/sqrt(2)',
                     'scope': 'Tensor-product character control; no spatial fission dynamics or preparation is asserted.',
                     'local_expectations': [.5, .5],
                     'joint_charge_outcomes': [{'charges': [1, 0], 'probability': .5},
                                               {'charges': [0, 1], 'probability': .5}],
                     'both_local_sharp_half_charges': False}
    checks = {'same_actual_384_source': model.n == 384,
              'charge_equals_retained_minus_B': bool(np.array_equal(charge, -archive['boundary_B'])),
              'same_H_and_every_native_branch_preserve_charge': h_error < 1e-12 and native_error < 1e-12,
              'all_actual_reservoir_branches_intertwine_charge': reservoir_error < 1e-12,
              'compiled_actual_gauge_character_degrees': actual_degrees == expected_degrees == {q: q for q in range(-3, 4)},
              'complete_M6_exports_include_failure_complements': export_error == 0.,
              'complete_source_compiled_pointer_instrument': max_complete < 1e-10 and max_off < 1e-10,
              'all_four_branches_and_both_axes_and_pointer_values': len(complete_records) == 528 and max(abs(x - 1) for x in sums) < 1e-10,
              'compiled_probabilities_match_actual_charge_action': max_probability_error < 1e-10,
              'nonclosed_half_and_zero_crossing_mixture_rejected': len(rejected) == 2}
    require(all(checks.values()), 'gauge-loop certificate failed')
    paths = [Path(__file__), contract_path, archive_path,
             RER / 'model/observer/source_scalar/compiler.py',
             RER / 'code/source_selection_model/response.json',
             RER / 'code/m1_source_realization/receipt.json']
    paths += [BASE / rel for rel in set(contract['source_sha256']) | set(joint_receipt['source_sha256'])]
    return {'schema': 'oph.row1.gauge_loop.v1', 'checks': checks,
        'construction_id': contract['construction_id'],
        'source_charge': 'Q=N-3=-B, the actual relative source number; no electric SI unit assigned',
        'native_operators_checked': native_operators, 'CAR_creator_edges_checked': car_edges,
        'all_charge_character_degrees': actual_degrees, 'actual_reservoir_branches': branches,
        'phase_settings': angles.tolist(),
        'reconstructed_character_loops': {str(q): [[float(z.real), float(z.imag)] for z in values] for q, values in reconstructed.items()},
        'complete_records': complete_records, 'example_executed_source_words': example_words,
        'proper_code_export_failure_ranks': export_failure_ranks,
        'errors': {'H_charge_commutator': h_error, 'native_charge_intertwiner': native_error,
                   'reservoir_charge_intertwiner': reservoir_error, 'compiled_offdiagonal': max_off,
                   'compiled_instrument_completeness': max_complete, 'compiled_probability': max_probability_error,
                   'full_export_completeness': export_error},
        'rejecting_controls': rejected, 'fractional_expectation_control': split_control,
        'exact_statement': 'Each sharp actual relative-number sector carries the gauge character exp(i q theta), whose winding is q. Every internal H/native branch preserves it; reservoir branches shift it by their actual retained integer current. Character tensor products add winding, excluding two equal sharp half characters in the same representation.',
        'source_services': ['Existing twelve-generator M6 exponentials execute all pointer/gauge words.',
                           'Coherent proper-code load/export, blank pointer and analog angles are declared services.',
                           'All four reservoir outcomes, axis choices and pointer outcomes are retained.',
                           'Proper-code export has its existing full complement failure effects; valid code has zero ideal failure. Invalid C10 encodings are outside physical admission, not accepted charge states.',
                           'The compound apparatus preserves encoded Gauss; its intermediate processor pulses are not native C10 repair events.'],
        'remaining_row1_boundary': ['Gauge-orbit topology is not identified with spatial defect degree.',
                                   'The primitive CAR gauge representation and reference preparation are declared architecture.',
                                   'This implements an observable on C10 using declared coherent source services; it does not derive C10 H or CAR content from passive repair.',
                                   'No physical electric dictionary, apparatus implementation from committed repair alone, or maintainer row award follows.'],
        'numerics': 'Binary64 source-pulse execution; exact source identities justify ideal gates, no outward rounding certificate.',
        'source_sha256': {str(p.relative_to(WORKSPACE)): sha(p) for p in paths},
        'full_criterion_met': False, 'official_rows_changed': []}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = evaluate()
    if args.write:
        (HERE / 'GAUGE_LOOP.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'checks': result['checks'], 'errors': result['errors'],
                      'degrees': result['all_charge_character_degrees']}, indent=2))
