"""Frozen-source replay of the finite native scalar construction.

Verification is read-only. Recording requires an explicit unused output path
and an exact production-source freeze. Scientific gates run under Python -O.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

from tools.pin_scope import bound, bound_pins
from .carrier import load_source
from .compiler import require

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]

# Receipt gates match instruments.TOL and the stricter compiler, fusion and
# contrast checks. They are independent of the backend replay allowance.
INSTRUMENT_TOL = 2e-11
INSTRUMENT_ERRORS = frozenset((
    'decode_residual', 'decoded_matrix_error', 'pulse_unitarity_error',
    'parity_reference_error', 'source_projector_error', 'frame_isometry_error',
    'largest_circulation_occupation_error', 'largest_probability_roundoff',
    'work_balance_residual', 'largest_source_work_balance_residual'))
STRICT_ERRORS = frozenset((
    'write_isometry_error', 'full_involution_error', 'blank_memory_leakage',
    'full_99_state_seam_error', 'full_99_reversed_seam_error', 'state_norm_error',
    'selected_code_Gram_error', 'collective_SO3_rotation_error',
    'logical_Pauli_involution_and_anticommutation_error',
    'retained_colored_complement_amplitude_error'))
# Additional fusion diagnostics use the same 1e-12 receipt ceiling.
FUSION_DIAGNOSTICS = frozenset((
    'local_Clebsch_transform_full_unitarity_error', 'code_output_vector_error',
    'passive_transported_chart_amplitude_error', 'physical_pair_transform_row_isometry_error'))


def _finite(value):
    try:
        return type(value) in (int, float) and math.isfinite(value)
    except OverflowError:
        return False


def _relation(value, expected, path, *operands):
    """Check arithmetic within one receipt, allowing only summation roundoff.

    This is not a second backend comparison. The ULP allowance accounts for
    subtraction of source energies and NumPy versus Python reductions, and is
    capped below the instrument gate even for large operands.
    """
    require(all(_finite(x) for x in (value, expected, *operands)),
            'Invalid numerical relation: '+path)
    tolerance = min(INSTRUMENT_TOL / 16,
                    16 * max(math.ulp(float(x)) for x in (1, value, expected, *operands)))
    require(abs(value-expected) <= tolerance, 'Inconsistent saved relation: '+path)


def _receipt_fields(record, path):
    """Check producer gates and linked fields of the fixed receipt schema."""
    for key, value in record.items():
        field = path+'.'+key
        if key in ('probability', 'pointer_probability', 'code_probability', 'dirty_probability',
                   'success_probability_after_feedforward', 'no_adjacent_pair_probability',
                   'largest_pointer_probability_difference'):
            require(_finite(value) and 0 <= value <= 1, 'Invalid scalar probability: '+field)
        if key in INSTRUMENT_ERRORS:
            require(_finite(value) and 0 <= value <= INSTRUMENT_TOL,
                    'Scientific field exceeds its gate: '+field)
        elif key in STRICT_ERRORS | FUSION_DIAGNOSTICS or key == 'complex_hopping_formula_dense_error':
            bound = 1e-13 if key == 'complex_hopping_formula_dense_error' else 1e-12
            require(_finite(value) and 0 <= value < bound,
                    'Scientific field exceeds its gate: '+field)
        elif key == 'source_particle_number_error':
            require(_finite(value) and 0 <= value <= 2*math.sqrt(150)*INSTRUMENT_TOL,
                    'Scientific field exceeds its gate: '+field)
        elif key == 'source_particle_number_tolerance':
            require(_finite(value) and value == 2*math.sqrt(150)*INSTRUMENT_TOL,
                    'Changed particle-number gate: '+field)
        elif key == 'full_78_state_operator_error':
            require(_finite(value) and value == 0, 'Scientific field exceeds its gate: '+field)
        elif key == 'code_errors':
            require(type(value) is list and len(value) == 3 and
                    all(_finite(v) and 0 <= v < 1e-12 for v in value),
                    'Scientific field exceeds its gate: '+field)
        elif key in ('largest_independent_owner_rotation_change',
                     'load_only_orientation_control_error', 'largest_pointer_probability_difference'):
            lower = {'largest_independent_owner_rotation_change': .1,
                     'load_only_orientation_control_error': .5,
                     'largest_pointer_probability_difference': .01}[key]
            require(_finite(value) and value > lower, 'Scientific field exceeds its gate: '+field)
        elif key == 'all_vacuum_codeword_amplitude_real':
            require(_finite(value) and abs(value-1) < 1e-12,
                    'Scientific field exceeds its gate: '+field)
        elif key == 'raw_single_weak_owner_positive_control':
            require(_finite(value) and value == 1, 'Scientific field exceeds its gate: '+field)
    if 'steps_per_star' in record and record['steps_per_star'] == 256:
        require(_finite(record['distance_from_exact_sea']) and
                0 <= record['distance_from_exact_sea'] < .0058,
                'Scientific field exceeds its gate: '+path+'.distance_from_exact_sea')
    if 'source_energy_before' in record and 'source_energy_after' in record:
        before, after = record['source_energy_before'], record['source_energy_after']
        work = record['external_driving_work']
        require(all(_finite(v) for v in (before, after, work)), 'Invalid source energy ledger: '+path)
        change = after-before
        residual = abs(change-work)
        require(residual <= INSTRUMENT_TOL, 'Source energy/work balance exceeds its gate: '+path)
        _relation(record['source_energy_change'], change, path+'.source_energy_change', before, after)
        _relation(record['work_balance_residual'], residual, path+'.work_balance_residual', before, after)
        if record['stage'] in ('pulse_1', 'pulse_2'):
            _relation(work, record['switch_on_work']+record['switch_off_work'],
                      path+'.external_driving_work', record['switch_on_work'], record['switch_off_work'])
        else:
            require(work == 0, 'Non-pulse source work must be zero: '+path)
    if 'net_source_energy_change' in record:
        stages = record['stages']
        require([s['stage'] for s in stages] ==
                ['pulse_1', 'free_interval', 'fuel_phase', 'pulse_2', 'terminal_pointer'],
                'Incomplete energy ledger: '+path)
        for previous, following in zip(stages, stages[1:]):
            _relation(following['source_energy_before'], previous['source_energy_after'],
                      path+'.stage_continuity')
        before, after = stages[0]['source_energy_before'], stages[-1]['source_energy_after']
        change = after-before
        work = sum(s['external_driving_work'] for s in stages)
        residual = abs(change-work)
        require(residual <= INSTRUMENT_TOL, 'Source energy/work balance exceeds its gate: '+path)
        _relation(record['net_source_energy_change'], change, path+'.net_source_energy_change', before, after)
        _relation(record['total_external_driving_work'], work, path+'.total_external_driving_work')
        _relation(record['work_balance_residual'], residual, path+'.work_balance_residual', before, after)
    if 'raw_final_fuel_probability' in record:
        raw = record['raw_final_fuel_probability']
        require(_finite(raw) and -INSTRUMENT_TOL <= raw <= 1+INSTRUMENT_TOL,
                'Raw fuel probability exceeds its gate: '+path)
        law = record['pointer_probabilities']
        p = min(1., max(0., raw))
        require(len(law) == 4 and law[1] == law[2] == 0, 'Forbidden pointer support: '+path)
        for value, expected in zip(law, (1-p, 0, 0, p)):
            _relation(value, expected, path+'.pointer_probabilities')
    if 'checks' in record and 'records' in record and 'target' in record:
        records = record['records']
        roundoff = max(abs(min(1., max(0., r['raw_final_fuel_probability']))-
                           r['raw_final_fuel_probability']) for r in records)
        residual = max(v['work_balance_residual'] for r in records
                       for v in [r['source_energy_ledger'], *r['source_energy_ledger']['stages']])
        _relation(record['checks']['largest_probability_roundoff'], roundoff, path+'.largest_probability_roundoff')
        _relation(record['checks']['largest_source_work_balance_residual'], residual,
                  path+'.largest_source_work_balance_residual')
    if 'unnormalized_weights' in record:
        weights = record['unnormalized_weights']
        law = record['probabilities']
        require(type(weights) is list and len(weights) == len(law) and
                all(_finite(v) and v >= 0 for v in weights), 'Invalid contrast weights: '+path)
        total = sum(weights)
        require(abs(total-1) < 1e-12 and _finite(record['numerical_normalization']) and
                abs(record['numerical_normalization']-1) < 1e-12,
                'Contrast normalization exceeds its gate: '+path)
        _relation(record['numerical_normalization'], total, path+'.numerical_normalization')
        _relation(record['state_norm_error'], abs(record['numerical_normalization']-1), path+'.state_norm_error')
        for weight, probability in zip(weights, law):
            _relation(probability, weight/record['numerical_normalization'], path+'.probabilities')
    if 'outcomes' in record and 'kappa' in record and 'records' in record:
        edges = {tuple(map(int, edge)) for edge in load_source()['edges']}
        for r in record['records']:
            require(len(r['probabilities']) == len(record['outcomes']), 'Incomplete contrast law: '+path)
            p = sum(prob for outcome, prob in zip(record['outcomes'], r['probabilities'])
                    if tuple(outcome) not in edges)
            _relation(r['no_adjacent_pair_probability'], p, path+'.no_adjacent_pair_probability')
        peak = max(r['no_adjacent_pair_probability'] for r in record['records'])
        require(peak <= .4 if record['kappa'] == 8 else peak > .4,
                'Contrast control exceeds its gate: '+path)


def source_pins():
    paths = [p.relative_to(ROOT).as_posix() for p in HERE.glob('*.py')]
    paths += ['model/candidate_laws/contrast_memory_spherical_v1/replay.py', 'model/candidate_laws/scalar_channel/protocol.json',
              'model/executable_model/carrier.py', 'model/candidate_laws/contrast_memory_spherical_v1/law.py',
              'model/candidate_laws/contrast_memory_spherical_v1/sectors.py',
              'model/executable_model/carriers/CELL.npz', 'tools/pin_scope.py', 'code/pin_scope.py']
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
            for p in sorted(bound(paths))}


def load_json(path):
    def finite_float(token):
        value = float(token)
        require(math.isfinite(value), 'Nonfinite JSON number: '+token)
        return value
    def pairs(items):
        out = {}
        for k, v in items:
            require(k not in out, 'Duplicate JSON key: '+k)
            out[k] = v
        return out
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=pairs,
                      parse_float=finite_float,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError('Nonfinite JSON: '+x)))


def compare(saved, actual, path='result', exact_zero=True):
    if type(actual) is dict:
        require(type(saved) is dict and saved.keys() == actual.keys(), 'Different keys: '+path)
        for k in actual:
            compare(saved[k], actual[k], path+'.'+k)
        _receipt_fields(saved, path)
        _receipt_fields(actual, path)
    elif type(actual) is list:
        require(type(saved) is list and len(saved) == len(actual), 'Different sequence: '+path)
        if path.endswith(('probabilities', 'complete_basis_law')):
            for law in (saved, actual):
                require(all(_finite(p) and 0 <= p <= 1
                            for p in law) and abs(sum(law)-1) < 1e-12,
                        'Invalid probability law: '+path)
        for i, (a, b) in enumerate(zip(saved, actual)):
            # For the blank two-pointer instrument only 10 and 01 are
            # forbidden. A valid 00 or 11 weight can round to zero without
            # changing its membership in the allowed outcome support.
            structural = not (path.endswith('pointer_probabilities') and
                              len(actual) == 4 and i in (0, 3))
            compare(a, b, path+f'[{i}]', exact_zero=structural)
    elif type(actual) is float:
        require(_finite(saved) and _finite(actual),
                'Invalid numerical value: '+path)
        # Forbidden outcome support is exact where the support is structural:
        # the blank two-pointer law and the integer phase witness. Calculated
        # differences of energies, and the contrast law that an eigensolver
        # produces, can be zero in one backend and roundoff in another.
        structural_law = any(key+'[' in path for key in
                             ('pointer_probabilities', 'complete_basis_law'))
        require(saved == 0 if structural_law and exact_zero and actual == 0 else
                abs(saved-actual) <= 2e-11*(1+abs(actual)),
                'Numerical replay differs: '+path)
    else:
        require(type(saved) is type(actual) and saved == actual, 'Changed value: '+path)


def _plain(value):
    if isinstance(value, dict):
        return {str(k): _plain(v) for k,v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def validate_protocol():
    """Reject a frozen declaration that differs from this fixed implementation."""
    protocol = load_json(HERE/'protocol.json')
    require(protocol['schema'] == 'oph.native_base_protocol.v1', 'Unknown protocol schema')
    expected = {
        'scalar': {'rest_level': 1, 'half_link_coefficient': '1/sqrt(2)',
                   'penalty': 1024, 'loop_sign_seed': 193019,
                   'default_ramp_substeps_per_star': 16,
                   'extended_ramp_substeps_per_star': 256,
                   'total_ramp_action_per_star': 16},
        'contrast': {'coupling': 8, 'competing_coupling': 0, 'cell_load': 2,
                     'sample_actions': [0, 0.125, 0.5, 1, 2, 4, 8]},
        'apparatus': {'pulse_strength': 64, 'pulse_action': 'pi/4',
                      'queries': [0, 4, 8], 'source_action_per_count': '1/8',
                      'phases': ['0', 'pi/2', 'pi', '3*pi/2'],
                      'targets': ['single_seam', 'circulation', 'five_seam', 'complex_seed']},
        'architecture': {'code_bits': 432, 'source_modes': 150, 'fuel_modes': 12,
                         'parity_reference_modes': 12, 'additional_local_scratch_bits': 12,
                         'additional_local_pointer_bits': 2,
                         'total_coherent_bits_including_scratch_and_pointers': 446,
                         'controller': 'supplied finite local classical program',
                         'record_policy': 'retain all loop, diagnostic and failure records; no archive reuse',
                         'whole_star_code_gates': 'supplied full unitary extensions, identity on complement',
                         'owner_initialization': 'supplied vacuum on each owner'},
    }
    for name, values in expected.items():
        compare(protocol[name], values, 'protocol.'+name)
    return protocol


def scientific_result(extended=False):
    protocol = validate_protocol()
    from . import replay
    from .routed import run as routed
    from .instruments import run as instruments
    from .gauge_capacity import (compiler_capacity, current_algebra_dimensions,
                                 fock_capacity, phase_response_witness)
    from model.candidate_laws.contrast_memory_spherical_v1.replay import run as contrast, primary_join
    d = load_source()
    with np.load(ROOT/'model/executable_model/carriers/CELL.npz', allow_pickle=False) as original:
        require(np.array_equal(d['edges'], original['edges']), 'Imported carrier edges differ')
        require(np.array_equal(d['wedges'], original['wedges']), 'Imported carrier wedges differ')
    result = replay()
    result['routed_preparation'] = routed()
    result['instruments'] = instruments(extended=extended)
    result['capacity'] = compiler_capacity()
    result['seven_species_capacity'] = fock_capacity(7)
    result['raw_currents'] = current_algebra_dimensions()
    phase = phase_response_witness()['records']
    result['phase_records'] = {k: {key: [int(p) for p in v] if isinstance(v, tuple) else int(v)
                                 for key,v in r.items()} for k,r in phase.items()}
    result['contrast_alternative'] = contrast()
    result['contrast_primary_join'] = primary_join()
    return _plain(dict(schema='oph.native_base_replay.v1', extended=extended,
                       protocol=protocol, scientific=result))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended', action='store_true')
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--record', type=Path)
    modes.add_argument('--verify', type=Path)
    parser.add_argument('--freeze', type=Path, default=HERE/'FREEZE.json')
    args = parser.parse_args(argv)
    freeze = load_json(args.freeze)
    require(bound_pins(freeze['source_sha256']) == source_pins(), 'Production source freeze differs')
    require(freeze['schema'] == 'oph.native_base_freeze.v1', 'Unknown freeze schema')
    validate_protocol()
    if args.record:
        require(not args.record.exists(), 'Refusing to replace recorded evidence')
    result = scientific_result(args.extended)
    result['source_sha256'] = source_pins()
    # A new receipt has the same domain gates as a replayed receipt.
    compare(result, result)
    path = args.verify or HERE/('extended.json' if args.extended else 'results.json')
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        with args.record.open('x', encoding='utf-8', newline='\n') as f:
            json.dump(result, f, sort_keys=True, separators=(',', ':'), allow_nan=False)
            f.write('\n')
        print(json.dumps({'recorded': str(args.record), 'extended': args.extended}))
    else:
        saved = load_json(path)
        # Schema, values, types, structural zeros, failed controls and complete
        # outcome arrays are all compared; no unknown fields are ignored.
        compare(saved, result)
        print(json.dumps({'verified': True, 'extended': args.extended,
                          'construction': 'OPH-NSC-1'}))
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, KeyError, OSError) as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1)
