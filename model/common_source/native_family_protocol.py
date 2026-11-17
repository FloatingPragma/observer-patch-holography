"""Bounded source-bound receipt for the declared current/family/field protocol.

The baseline experiment executes complete correlated classical-memory records
and one logical projector-feedback pulse at two clocks. The interaction-gated
channels instead supply stationary coefficients to a declared static mean
feedback service. Their complete detector records are retained separately;
no equality of those two finite processes is asserted. Both services, the
interaction duration, current filter, Higgs background and Weyl symbol remain
declared architecture. This receipt contains no empirical mass packet and
establishes no physical pole masses or interacting continuum field theory.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import sympy as s

from tools import strict_json
from . import native_family_record_memory as memory
from . import native_filtered_current_control as control
from . import native_gauge_family_records as gates
from . import native_record_yukawa_field as field

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = Path(__file__).with_name('native_family_protocol_results.json')
SCHEMA = 'oph.native_family_protocol.v1'
REPLAY_ADMISSION = 3e-11
RATIO_REPLAY_RELATIVE = 1e-8
CHANNELS = ('charged_lepton', 'down', 'up')
DUTIES = ('always', 'C', 'E', 'CE')
ARRAY_SCHEMA = {
    'memory_probabilities': ((2, 9, 1536, 3), 'float64'),
    'memory_labels': ((1536, 6), 'int64'),
    'memory_burnin_memories': ((9, 12), 'float64'),
    'memory_transition': ((12, 12), 'float64'),
    'memory_raw_labels': ((128, 4), 'int64'),
    'memory_raw_weights': ((128, 2), 'int64'),
    'memory_pulses': ((2, 12, 6, 6), 'complex128'),
    'filtered_current_generators': ((12, 6, 6), 'complex128'),
    'filtered_current_projector': ((6, 6), 'float64'),
    'current_kraus': ((13, 15, 15), 'complex128'),
    'current_covariance_errors': ((12,), 'float64'),
    'gate_binary_labels': ((3, 512, 8), 'int64'),
    'gate_binary_weights': ((3, 512), 'float64'),
    'gate_left_support': ((3, 15, 15), 'float64'),
    'gate_right_support': ((3, 15, 15), 'float64'),
    'gate_Higgs_components': ((3, 2, 15, 15), 'complex128'),
    'gate_coefficients': ((3, 3, 3), 'float64'),
    'gate_rates': ((3, 2), 'float64'),
    'catalog_gate_rates': ((3, 3, 2), 'float64'),
    'catalog_stationary': ((3, 3, 12), 'float64'),
    'catalog_coefficients': ((3, 3, 3, 3), 'float64'),
    'duty_coefficients': ((4, 3, 3, 3), 'float64'),
    'duty_projector_weights': ((4, 3, 12), 'float64'),
    'duty_traces': ((4, 3), 'float64'),
    'duty_ranks': ((4, 3), 'int64'),
    'static_probabilities': ((3, 3, 3, 2, 3), 'float64'),
    'static_survival': ((3, 3, 2), 'float64'),
    'static_recovered': ((3, 3, 3), 'complex128'),
    'field_kernel': ((45, 45), 'complex128'),
    'field_bdg': ((180, 180), 'complex128'),
    'field_spectrum': ((180,), 'float64')}


def _json(value):
    if isinstance(value, np.ndarray):
        return _json(value.tolist())
    if isinstance(value, np.generic):
        return _json(value.item())
    if isinstance(value, dict):
        return {key: _json(child) for key, child in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json(child) for child in value]
    return value


def _compare(actual, expected, path='receipt'):
    if type(actual) is not type(expected):
        raise ValueError(path+' original type differs')
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            raise ValueError(path+' complete keys differ')
        for key in expected:
            _compare(actual[key], expected[key], path+'.'+key)
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise ValueError(path+' complete list differs')
        for index, (left, right) in enumerate(zip(actual, expected, strict=True)):
            _compare(left, right, path+f'[{index}]')
    elif type(expected) is float:
        # Small eigenvalue denominators amplify platform LAPACK roundoff in
        # displayed ratios. Complete primary arrays keep absolute admission;
        # this relative allowance applies only to scalar diagnostics. It is
        # neither a physical/theory uncertainty nor an empirical fit tolerance.
        relative = RATIO_REPLAY_RELATIVE if 'ratio' in path or 'span' in path else REPLAY_ADMISSION
        allowance = REPLAY_ADMISSION+relative*abs(expected)
        if not math.isfinite(actual) or abs(actual-expected) > allowance:
            raise ValueError(path+' finite numerical value differs')
    elif actual != expected:
        raise ValueError(path+' exact claim differs')


def source_bindings():
    paths = (set(control.response.source_bindings()) | set(gates.dressing.source_bindings())
             | set(memory.feedback.transport.source_bindings()))
    names = ('native_filtered_current_control', 'native_family_record_memory',
             'native_gauge_family_records', 'native_record_yukawa_field', 'native_family_protocol')
    paths.update('model/common_source/'+name+suffix for name in names
                 for suffix in ('.py',))
    paths.update('model/common_source/test_'+name+'.py' for name in names)
    paths.update(('model/common_source/native_repair_mass_readout.py',
                  'model/common_source/native_family_selection.py',
                  'model/common_source/native_matter_theorems.py',
                  'model/common_source/chiral_candidate/kinetic.py', 'model/requirements.txt',
                  'tools/strict_json.py'))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def build():
    """Run every declared finite producer; keep complete numerical laws in arrays."""
    # Derive the shared exact control/exterior objects first. Their owning
    # functions cache immutable algebra; complete records are still freshly
    # produced and validated on every build/check.
    current = deepcopy(control.run())
    baseline, baseline_records = memory.execute()
    memory.validate_records(baseline_records)
    data = memory.architecture()
    exact = current['exact']
    # These operators and all commutators are recomputed by control.run().
    # Retain the exact invertible remaps in JSON; large operator tables use
    # the binary artifact rather than a second matrix serialization.
    exact['filtered_port_generator_count'] = len(exact.pop('filtered_port_generators'))
    exact['commutator_closure_count'] = len(exact.pop('commutator_closure'))
    current_data = control.exact_architecture()
    interaction = gates.summary()
    instrument = gates.instrument()
    channels = gates.candidate_channels()
    if tuple(sorted(channels)) != CHANNELS:
        raise ArithmeticError('complete native action-channel menu changed')
    records = [field.forward_records(channels[name]['Y']) for name in CHANNELS]
    recovered = [field.validate_records(channels[name]['Y'], row)
                 for name, row in zip(CHANNELS, records, strict=True)]
    binary = [gates.complete_binary_records(tuple(map(s.sympify, channels[name]['charge_origin'])))
              for name in CHANNELS]
    catalog_models = []
    for angle in gates.protocol_catalog().values():
        by_origin = {tuple(map(str, row['gate']['origin'])): row for row in gates.family_models(angle)}
        catalog_models.append([by_origin[channels[name]['charge_origin']] for name in CHANNELS])
    duty_models = []
    duties = gates.duty_coefficients()
    if tuple(duties) != DUTIES:
        raise ArithmeticError('complete unconditional pulse-duty menu changed')
    for policy in DUTIES:
        by_channel = {row['action_channel']: row for row in duties[policy]}
        duty_models.append([by_channel[name] for name in CHANNELS])
    field_report = field.summary()
    K = field.assemble_kernel(channels)
    diagnostic = field_report['free_field_diagnostic']
    H = field.kinetic.bdg_symbol(diagnostic['momentum'], K, speed=diagnostic['coordinate_speed'])
    arrays = {
        'memory_probabilities': baseline_records['probabilities'],
        'memory_labels': baseline_records['labels'],
        'memory_burnin_memories': baseline_records['burnin_memories'],
        'memory_transition': np.asarray(data['T'], float),
        'memory_raw_labels': np.asarray([row[:4] for row in data['records']], np.int64),
        'memory_raw_weights': np.asarray([[int(row[-1].p), int(row[-1].q)] for row in data['records']], np.int64),
        'memory_pulses': np.asarray([memory.feedback_unitaries(h) for h in memory.READOUT_STEPS]),
        'filtered_current_generators': np.asarray(current_data['generators'], complex),
        'filtered_current_projector': np.asarray(current_data['P'], float),
        'current_kraus': instrument['kraus'],
        'current_covariance_errors': gates.covariance_errors(instrument['kraus']),
        'gate_binary_labels': np.asarray([[row[:-1] for row in rows] for rows in binary], np.int64),
        'gate_binary_weights': np.asarray([[float(row[-1]) for row in rows] for rows in binary]),
        'gate_left_support': np.asarray([channels[name]['left_projector'] for name in CHANNELS]),
        'gate_right_support': np.asarray([channels[name]['right_projector'] for name in CHANNELS]),
        'gate_Higgs_components': np.asarray([channels[name]['Higgs_components'] for name in CHANNELS]),
        'gate_coefficients': np.asarray([channels[name]['Y'] for name in CHANNELS]),
        'gate_rates': np.asarray([[channels[name]['p'], channels[name]['q']] for name in CHANNELS]),
        'catalog_gate_rates': np.asarray([[[float(row['gate']['p']), float(row['gate']['q'])]
                                           for row in models] for models in catalog_models]),
        'catalog_stationary': np.asarray([[list(row['family']['stationary']) for row in models]
                                           for models in catalog_models], float),
        'catalog_coefficients': np.asarray([[row['family']['Y'] for row in models]
                                             for models in catalog_models], float),
        'duty_coefficients': np.asarray([[row['coefficient'] for row in models] for models in duty_models], float),
        'duty_projector_weights': np.asarray([[row['projector_weights'] for row in models] for models in duty_models], float),
        'duty_traces': np.asarray([[float(row['trace']) for row in models] for models in duty_models]),
        'duty_ranks': np.asarray([[row['family_rank_exact'] for row in models] for models in duty_models], np.int64),
        'static_probabilities': np.asarray([row['probabilities'] for row in records]),
        'static_survival': np.asarray([row['survival'] for row in records]),
        'static_recovered': np.asarray(recovered, complex),
        'field_kernel': np.asarray(K, complex), 'field_bdg': np.asarray(H, complex),
        'field_spectrum': np.linalg.eigvalsh(H)}
    report = {'schema': SCHEMA, 'filtered_current': current, 'baseline_memory': baseline,
              'interaction_records': interaction, 'field': field_report,
              'array_order': {'action_channels': list(CHANNELS),
                              'geometry_protocols': list(gates.protocol_catalog()),
                              'pulse_duties': list(DUTIES),
                              'memory_clocks': list(memory.READOUT_STEPS),
                              'static_detector_angle': field.PULSE_ANGLE},
              'record_processes': {
                  'baseline': 'finite record-driven memory preparation plus actual one-pulse CQ feedback at each fixed clock',
                  'interaction_gated': 'complete current instruments and stationary family coefficients; correlated gated quantum history not executed here',
                  'field_inverse': 'complete detector records of a declared static mean-feedback service',
                  'static_mean_service_identified_with_finite_CQ_memory': False},
              'scope': {'empirical_mass_packet_consumed': False,
                        'all_geometry_catalog_entries_retained': True,
                        'target_selected_catalog_winner': False,
                        'source_protocol_is_preregistered_physical_prediction': False,
                        'logical_services_are_literal_compiled_source_histories': False,
                        'numerical_acceptance_is_outward_physical_error_certificate': False,
                        'dimensionless_coefficients_are_physical_pole_masses': False,
                        'physical_spectrum_or_interacting_continuum_derived': False}}
    report = _json(report)
    validate_arrays(report, arrays)
    return report, arrays


def _schema(arrays):
    if type(arrays) is not dict or set(arrays) != set(ARRAY_SCHEMA):
        raise ValueError('complete bounded protocol array menu required')
    for key, (shape, dtype) in ARRAY_SCHEMA.items():
        value = arrays[key]
        if (type(value) is not np.ndarray or value.shape != shape or value.dtype != np.dtype(dtype)
                or not np.isfinite(value).all()):
            raise ValueError('invalid complete protocol array '+key)


def validate_arrays(report, arrays):
    """Own-array checks and source binding, in addition to fresh replay at check()."""
    _schema(arrays)
    if report.get('schema') != SCHEMA or report.get('array_order') != {
            'action_channels': list(CHANNELS), 'geometry_protocols': list(gates.protocol_catalog()),
            'pulse_duties': list(DUTIES),
            'memory_clocks': list(memory.READOUT_STEPS), 'static_detector_angle': field.PULSE_ANGLE}:
        raise ValueError('declared complete protocol or settings changed')
    baseline_records = {'probabilities': arrays['memory_probabilities'], 'labels': arrays['memory_labels'],
                        'burnin_memories': arrays['memory_burnin_memories'],
                        'burnin': memory.BURNIN, 'steps': memory.READOUT_STEPS}
    memory.validate_records(baseline_records)
    channels = gates.candidate_channels()
    for index, name in enumerate(CHANNELS):
        records = {'probabilities': arrays['static_probabilities'][index],
                   'survival': arrays['static_survival'][index], 'angle': field.PULSE_ANGLE}
        decoded = field.validate_records(channels[name]['Y'], records)
        if np.max(abs(decoded-arrays['static_recovered'][index])) > REPLAY_ADMISSION:
            raise ValueError('retained static inverse differs from its own complete records')
        L, R, C = (arrays[key][index] for key in ('gate_left_support', 'gate_right_support', 'gate_Higgs_components'))
        if (max(np.linalg.norm(L@L-L), np.linalg.norm(R@R-R), np.linalg.norm(L@R)) > REPLAY_ADMISSION
                or any(np.linalg.norm((L+R)@component-component) > REPLAY_ADMISSION for component in C)):
            raise ValueError('retained full Higgs supports lost their domain')
    gates.validate_kraus(arrays['current_kraus'])
    if np.max(gates.covariance_errors(arrays['current_kraus'])) > REPLAY_ADMISSION:
        raise ValueError('retained complete current dilation lost covariance')
    if (np.max(abs(arrays['gate_binary_weights'].sum(axis=1)-1)) > REPLAY_ADMISSION
            or np.min(arrays['gate_binary_weights']) < -REPLAY_ADMISSION):
        raise ValueError('all binary interaction records are required without postselection')
    diagnostic = field.dispersion_certificate(arrays['field_kernel'])
    if diagnostic['dispersion_spectral_error'] > REPLAY_ADMISSION:
        raise ValueError('complete declared Weyl symbol fails dispersion')
    _compare(diagnostic, report['field']['free_field_diagnostic'], 'receipt.field.free_field_diagnostic')
    H = arrays['field_bdg']
    if (np.linalg.norm(H-H.conj().T, 2) > REPLAY_ADMISSION
            or np.max(abs(np.linalg.eigvalsh(H)-arrays['field_spectrum'])) > REPLAY_ADMISSION):
        raise ValueError('complete retained180-dimensional field spectrum changed')


def check(path):
    path = Path(path)
    pins = source_bindings()
    retained = strict_json.load(path)
    if type(retained) is not dict or retained.get('source_bindings') != pins:
        raise ValueError('source bindings differ')
    artifact = retained.get('artifact')
    if (type(artifact) is not dict or set(artifact) != {'path', 'sha256'}
            or artifact['path'] != path.with_suffix('.npz').name):
        raise ValueError('the paired local array artifact is required')
    target = path.with_suffix('.npz')
    if hashlib.sha256(target.read_bytes()).hexdigest() != artifact['sha256']:
        raise ValueError('array artifact custody differs')
    with np.load(target, allow_pickle=False) as archive:
        if len(archive.files) != len(set(archive.files)):
            raise ValueError('duplicate array keys')
        arrays = {key: archive[key] for key in archive.files}
    report = {key: value for key, value in retained.items() if key not in ('source_bindings', 'artifact')}
    validate_arrays(report, arrays)
    fresh, source_arrays = build()
    _compare(report, fresh)
    for key, (_, dtype) in ARRAY_SCHEMA.items():
        if np.dtype(dtype).kind in 'iu':
            if not np.array_equal(arrays[key], source_arrays[key]):
                raise ValueError('source input labels or exact integers differ: '+key)
        elif np.max(abs(arrays[key]-source_arrays[key])) > REPLAY_ADMISSION:
            raise ValueError('complete retained law differs from fresh source: '+key)
    if source_bindings() != pins:
        raise RuntimeError('source changed during replay')
    return True


def write(path):
    path = Path(path)
    if path.suffix != '.json':
        raise ValueError('a JSON receipt destination is required')
    pins = source_bindings()
    report, arrays = build()
    if source_bindings() != pins:
        raise RuntimeError('source changed during generation')
    artifact = path.with_suffix('.npz')
    np.savez_compressed(artifact, **arrays)
    report.update(source_bindings=pins,
                  artifact={'path': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()})
    path.write_text(json.dumps(report, sort_keys=True, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--output', type=Path)
    action.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.check:
        check(args.check)
        print('native family protocol replay: PASS')
    else:
        write(args.output)
        print('native family protocol receipt generated')


if __name__ == '__main__':
    main()
