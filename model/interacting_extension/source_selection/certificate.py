"""Replay native references and clocks with exact receipts and bounded diagnostics."""
from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import os
import json
from pathlib import Path
import tempfile

import numpy as np

from .._source import ROOT, PIN, verify_sources
from . import native_record_dynamics as native

HERE = Path(__file__).resolve().parent
TOLERANCE = 3e-12
TOLERANCE_EXACT = '3/1000000000000'
OUTPUT = HERE / 'results.json'
ARRAYS = HERE / 'complete_native_instruments.npz'


def _finite_tree(value):
    if isinstance(value, dict):
        for child in value.values():
            _finite_tree(child)
    elif isinstance(value, list):
        for child in value:
            _finite_tree(child)
    elif isinstance(value, float) and not np.isfinite(value):
        raise ArithmeticError('nonfinite native diagnostic')


def _residual(value, name):
    if not np.isfinite(value) or abs(value) > TOLERANCE:
        raise ArithmeticError('native diagnostic exceeded its bound: '+name)


def stable_report(raw, arrays, repo=ROOT):
    """Keep exact data and verified bounds; floating values live in the array packet."""
    _finite_tree(raw)
    # Integer census keys acquire their JSON string spelling at this boundary.
    raw = json.loads(json.dumps(raw, allow_nan=False))
    result = {key: deepcopy(raw[key]) for key in
              ('source_commit', 'exact_native_preparation', 'exact_activity_quotient',
               'capacity_theorem', 'source_status')}
    result['schema'] = 'oph.native_record_dynamics.certificate.v1'
    result['floating_diagnostic_absolute_tolerance'] = TOLERANCE_EXACT
    result['candidates'] = {}
    for name, row in raw['candidates'].items():
        clean = {}
        for key, value in row.items():
            if isinstance(value, float):
                _residual(value, name+'/'+key)
                clean[key+'_within_bound'] = True
            elif key not in ('tilts', 'record_phase_audit', 'record_conditioned_geometric_reference'):
                clean[key] = deepcopy(value)
        clean['tilts'] = []
        for tilted in row['tilts']:
            for key, value in tilted.items():
                if key.endswith('_residual'):
                    _residual(value, name+'/'+key)
            if (tilted['top_eigenvalue_multiplicity_at_2e-12'] != 1
                    or tilted['tilted_transfer_minimum_eigenvalue'] <= 0
                    or tilted['all_operator_eigenvalues_count'] != 4356):
                raise ArithmeticError('incomplete positive native transfer spectrum')
            clean['tilts'].append({
                'z': str(Fraction(tilted['z'])),
                'complete_filter_record_labels': tilted['complete_filter_record_labels'],
                'all_operator_eigenvalues_count': tilted['all_operator_eigenvalues_count'],
                'positive_transfer_diagnostic': True,
                'simple_perron_diagnostic': True,
                'trace_preservation_stationarity_and_herald_residuals_within_bound': True,
                'euclidean_gap_is_a_physical_particle_mass': False})
        phase = row['record_phase_audit']
        for key in ('complete_committed_law_difference', 'activity_clock_law_difference',
                    'commit_first_pointer_law_difference'):
            _residual(phase[key], name+'/'+key)
        if phase['complete_coherent_pointer_total_variation'] <= 1e-3:
            raise ArithmeticError('coherent-pointer phase witness disappeared')
        clean['record_phase_audit'] = {
            'committed_and_count_laws_phase_invariant_within_bound': True,
            'coherent_pointer_total_variation_greater_than': '1/1000',
            'all_pointer_outcomes': phase['all_pointer_outcomes'],
            'seam_and_coin_randomizers_remain_committed': True}
        ref = row['record_conditioned_geometric_reference']
        clean_ref = {}
        for key, value in ref.items():
            if isinstance(value, float):
                _residual(value, name+'/'+key)
                clean_ref[key+'_within_bound'] = True
            elif key != 'two_record_conditional_spectra_groups_at_11_decimals_diagnostic':
                clean_ref[key] = deepcopy(value)
        clean['record_conditioned_geometric_reference'] = clean_ref
        result['candidates'][name] = clean
    certificates = []
    for row in raw['exact_perron_root_certificates']:
        if Fraction(row['floating_perron_value_absolute_error_upper_rational']) > Fraction(1, 10**14):
            raise ArithmeticError('floating Perron approximation exceeded its exact enclosure bound')
        clean = {key: value for key, value in row.items()
                 if not key.startswith('floating_perron_value_absolute_error_upper_')}
        clean['floating_perron_value_absolute_error_upper_rational'] = '1/100000000000000'
        certificates.append(clean)
    result['exact_perron_root_certificates'] = certificates
    names = sorted(raw['candidates']['C']['metadata']['source_hashes'])
    result['source_sha256'] = verify_sources(repo, names)
    implementation = [HERE/'native_record_dynamics.py', HERE/'certificate.py',
                      HERE.parent/'_source.py', HERE.parent/'source_manifest.json']
    # Package-relative keys, the convention model/verify.py uses for this field.
    result['implementation_sha256'] = {
        os.path.relpath(path, HERE).replace(os.sep, '/'): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in implementation}
    result['array_contract'] = {}
    for name, array in sorted(arrays.items()):
        if not np.isfinite(array).all():
            raise ArithmeticError('nonfinite native array: '+name)
        item = {'shape': list(array.shape), 'kind': array.dtype.kind}
        if array.dtype.kind in 'iu':
            item['integer_values_sha256'] = hashlib.sha256(
                array.astype('<i8').tobytes(order='C')).hexdigest()
        else:
            item['comparison_absolute_tolerance'] = TOLERANCE_EXACT
        result['array_contract'][name] = item
    return result


def replay(repo=ROOT):
    with tempfile.TemporaryDirectory(prefix='oph-native-record-') as directory:
        raw = native.audit(repo, directory)
        with np.load(Path(directory)/'complete_native_instruments.npz', allow_pickle=False) as packet:
            arrays = {name: packet[name].copy() for name in packet.files}
    return stable_report(raw, arrays, repo), arrays


def compare_arrays(expected, actual):
    if set(expected) != set(actual):
        raise ValueError('native array keys differ')
    for name in expected:
        left, right = expected[name], actual[name]
        if left.shape != right.shape or left.dtype.kind != right.dtype.kind:
            raise ValueError('native array contract differs: '+name)
        if not np.isfinite(left).all() or not np.isfinite(right).all():
            raise ValueError('nonfinite native array: '+name)
        if left.dtype.kind in 'iu':
            equal = np.array_equal(left, right)
        else:
            equal = np.max(np.abs(left-right), initial=0) <= TOLERANCE
        if not equal:
            raise ValueError('native array values differ: '+name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=ROOT)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    result, arrays = replay(args.repo)
    if args.write is not None:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
        np.savez_compressed(ARRAYS, **arrays)
    else:
        if json.loads(args.check.read_text()) != result:
            raise ValueError('native reference and clock certificate differs')
        with np.load(ARRAYS, allow_pickle=False) as packet:
            compare_arrays({name: packet[name] for name in packet.files}, arrays)
    print('Verified native references, correlated counts and complete heralded dynamics.')


if __name__ == '__main__':
    main()
