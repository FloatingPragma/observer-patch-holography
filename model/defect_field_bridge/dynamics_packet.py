"""Frozen programmable quasiparticle controls and complete-output replay."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from tools import strict_json
from ..defect_field_benchmark import quasiparticles as qp
from ..defect_protection import packet as protected_packet
from . import dynamics


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CONTROL = HERE/'dynamics_controls.json'
OUTPUT = HERE/'dynamics_results.json'
ARRAYS = HERE/'dynamics_results.npz'


def control():
    value = strict_json.load(CONTROL)
    if value.get('schema') != 'oph.colored_quasiparticle_controls.v1':
        raise ValueError('wrong frozen quasiparticle control schema')
    data = qp.quasiparticles()
    energies = np.asarray([float(x) for x in value['energies_GeV']])
    charge6 = np.asarray(value['charge6'])
    if energies.shape != (90,) or charge6.shape != (90,):
        raise ValueError('all90 programmed modes must remain')
    if not np.array_equal(charge6, data['charge6']):
        raise ValueError('programmed electric charge chart differs from the quasiparticle proof')
    if (np.max(abs(energies-data['energies'])) > 1e-8
            or abs(float(value['vacuum_energy_GeV'])-data['vacuum_energy']) > 1e-8):
        raise ValueError('frozen programmable spectrum left its numerical approximation admission')
    if value['ramsey_indices'] != list(data['ramsey_indices']):
        raise ValueError('frozen Ramsey energy/charge blocks changed')
    return value, energies, charge6/6


def source_bindings():
    values = protected_packet.source_bindings()
    files = [HERE/'dynamics.py', HERE/'dynamics_packet.py', HERE/'test_dynamics.py',
             HERE/'bridge.py', HERE/'colored.py', HERE/'REPRESENTATION.md',
             HERE/'DYNAMICS.md', CONTROL, Path(qp.__file__),
             ROOT/'model/common_source/calibrated_matter.py',
             ROOT/'model/common_source/calibrated_matter_parameters.json']
    for path in files:
        values[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return dict(sorted(values.items()))


def build():
    frozen, energies, charges = control()
    pins = source_bindings()
    report, arrays = dynamics.execute(energies, charges,
        ramsey_indices=tuple(frozen['ramsey_indices']), partner=frozen['partner_index'],
        duration=float(frozen['duration_GeV_inverse']),
        electromagnetic=float(frozen['electromagnetic_phase']),
        vacuum_energy=float(frozen['vacuum_energy_GeV']))
    report['programmable_control'] = frozen
    report['program_to_fresh_spectrum_numerical_admission_GeV'] = '1e-8'
    report['per_eigenvalue_program_quantization_GeV'] = '1e-9'
    report['control_status'] = 'rounded and frozen programmed parameters from explicitly calibrated MSbar(MZ) matter inputs; no refit'
    report['native_work_qubits'] = 9
    report['ordering_qubits'] = 1
    report['internal_work_qubits'] = 14
    report['fresh_internal_pointer_bits'] = 14
    report['fresh_neutral_charge_pointer_bits'] = 1
    report['terminal_decoded_ordering_read_bits'] = 1
    report['shared_internal_array_aliases'] = {
        row['array_prefix']+'_internal_factors':
            'shared_internal_'+('abort' if row['abort'] else str(row['current']))
        for rows in report['laws'].values() for row in rows}
    report['source_sha256'] = pins
    if source_bindings() != pins:
        raise RuntimeError('sources changed during colored dynamics replay')
    validate_arrays(report, arrays)
    return report, arrays


def storage_arrays(report, arrays):
    """Alias only tensors constructed by the same exact control formula."""
    output = {}
    aliases = report['shared_internal_array_aliases']
    for key, value in arrays.items():
        stored = aliases.get(key, key)
        if stored in output and not np.array_equal(output[stored], value):
            raise ValueError('shared control produced inconsistent full internal tensors')
        output[stored] = value
    return output


def load_arrays(report, path=ARRAYS):
    with np.load(path, allow_pickle=False) as data:
        stored = {key: data[key] for key in data.files}
    aliases = report['shared_internal_array_aliases']
    shared = set(aliases.values())
    if not shared.issubset(stored):
        raise ValueError('missing shared complete internal tensor')
    result = {key: value for key, value in stored.items() if key not in shared}
    result.update({key: stored[value] for key, value in aliases.items()})
    return result


def validate_arrays(report, arrays):
    expected = {row['array_prefix']+suffix for rows in report['laws'].values() for row in rows
                for suffix in ('_source_factors', '_internal_factors', '_all_pointer_probabilities')}
    if set(arrays) != expected:
        raise ValueError('missing or extra complete colored output array')
    for rows in report['laws'].values():
        total = 0.
        for row in rows:
            key = row['array_prefix']
            vectors, matrices, probabilities = [np.asarray(arrays[key+suffix]) for suffix in
                ('_source_factors', '_internal_factors', '_all_pointer_probabilities')]
            if (vectors.shape != (2, 512, 4) or matrices.shape != (2, 4, 128, 128)
                    or probabilities.shape != (2, 2, 128, 128)
                    or not all(np.isfinite(x).all() for x in (vectors, matrices, probabilities))
                    or np.min(probabilities) < 0):
                raise ValueError('invalid complete colored quantum output')
            for orientation in (0, 1):
                for charge in (0, 1):
                    actual = dynamics.pointer_probabilities(vectors[orientation], matrices[charge])
                    if np.max(abs(actual-probabilities[orientation, charge])) > 2e-9:
                        raise ValueError('saved pointer probabilities disagree with the quantum output')
            total += float(probabilities.sum())
        if abs(total-1) > dynamics.ADMISSION:
            raise ValueError('complete colored record tree is not normalized')


def compare_arrays(actual, expected):
    if set(actual) != set(expected):
        raise ValueError('complete colored array set changed')
    for key in actual:
        a, b = np.asarray(actual[key]), np.asarray(expected[key])
        if a.shape != b.shape or not np.isfinite(a).all() or np.max(abs(a-b)) > 2e-8:
            raise ValueError('complete colored source output changed: '+key)


def compare_json(actual, expected, path='report'):
    if isinstance(expected, float):
        if (type(actual) is not float or not np.isfinite(actual)
                or abs(actual-expected) > 2e-10*max(1., abs(expected))):
            raise ValueError('colored numeric report differs: '+path)
    elif type(actual) is not type(expected):
        raise ValueError('colored report type differs: '+path)
    elif isinstance(expected, dict):
        if actual.keys() != expected.keys():
            raise ValueError('colored report fields differ: '+path)
        for key in expected:
            compare_json(actual[key], expected[key], path+'.'+key)
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise ValueError('colored report list differs: '+path)
        for i, (a, b) in enumerate(zip(actual, expected, strict=True)):
            compare_json(a, b, path+'.'+str(i))
    elif actual != expected:
        raise ValueError('colored exact report differs: '+path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    report, arrays = build()
    if args.write:
        OUTPUT.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n', encoding='utf-8')
        np.savez_compressed(ARRAYS, **storage_arrays(report, arrays))
    elif args.check:
        retained = strict_json.load(OUTPUT)
        saved = load_arrays(retained)
        validate_arrays(retained, saved)
        compare_json(retained, report)
        compare_arrays(saved, arrays)
    print('Verified complete native colored quasiparticle dynamics and pointer outputs.')


if __name__ == '__main__':
    main()
