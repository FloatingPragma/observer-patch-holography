"""Independent numerical spectral audit of the actual native quark blocks.

The producer uses exact ground-field moments and polynomial residues. This
audit orthonormalizes the positive Gram forms, diagonalizes the literal
source/target Hamiltonians and solves their complex resolvents directly.
Floating checks do not replace the producer's rational interval proof.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

PIN = '3bcacf1276ddb3ead35be2b9adcd63b6643e2d35'


def orthonormalize(block):
    gram = np.asarray(block['gram'], dtype=np.complex128)
    h = np.asarray(block['H'], dtype=np.complex128)
    frame = np.linalg.cholesky(gram).conj().T
    converted = frame @ h @ np.linalg.inv(frame)
    error = float(np.max(np.abs(converted-converted.conj().T)))
    if error > 2e-10:
        raise ArithmeticError('physical metric transform is not Hermitian')
    return frame, converted, error


def resolvents(h, image, points):
    return np.array([np.vdot(image, np.linalg.solve(z*np.eye(len(h))-h, image)) for z in points])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rer', required=True, type=Path)
    parser.add_argument('--receipt', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    root = args.rer.resolve()
    if subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip() != PIN:
        raise ValueError('pinned source required')
    subprocess.run(['git', '-C', str(root), 'diff', '--exit-code', PIN, '--', 'model', 'code'], check=True)
    sys.path.insert(0, str(root))
    from model.native_quark_response.local import transfer
    from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import source
    source_block = source(1)
    source_frame, source_h, source_error = orthonormalize(source_block)
    source_energies, source_vectors = np.linalg.eigh(source_h)
    ground_coordinate = np.linalg.solve(source_frame, source_vectors[:, 0])
    channels = {}
    for name in ('uc', 'dc'):
        block = transfer(name)
        frame, h, error = orthonormalize(block)
        image = frame @ np.asarray(block['map'], dtype=np.complex128) @ ground_coordinate
        energies, vectors = np.linalg.eigh(h)
        weights = np.abs(vectors.conj().T @ image)**2
        channels[name] = {'h': h, 'image': image, 'energies': energies,
                          'weights': weights, 'hermitian_error': error}
    rng = np.random.default_rng(917346)
    points = rng.uniform(-8, 8, 17)+1j*rng.uniform(.5, 2, 17)
    up = resolvents(channels['uc']['h'], channels['uc']['image'], points)
    down = resolvents(channels['dc']['h'], channels['dc']['image'], points)
    defect = float(np.max(np.abs(up-down)))
    if defect > 2e-10:
        raise ArithmeticError('direct native quark resolvents disagree')
    changed_h = channels['dc']['h']+.01*np.eye(len(channels['dc']['h']))
    mutation_defect = float(np.max(np.abs(up-resolvents(changed_h, channels['dc']['image'], points))))
    if mutation_defect < 1e-6:
        raise ArithmeticError('true energy-translation mutation did not fire')

    receipt = json.loads(args.receipt.read_text())
    if receipt['source_commit'] != PIN or receipt['full_criterion_met'] is not False:
        raise ValueError('native receipt changed source or promoted a physical criterion')
    rows = receipt['complete_positive_records']
    total_lower = sum((Fraction(row['weight_interval'][0]) for row in rows), Fraction())
    total_upper = sum((Fraction(row['weight_interval'][1]) for row in rows), Fraction())
    if len(rows) != 13 or not total_lower < Fraction(1, 2) < total_upper:
        raise ArithmeticError('saved positive measure does not enclose total weight one half')
    maximum_energy_error = maximum_weight_error = 0.0
    for row in rows:
        e = sum(map(Fraction, row['excitation_interval']))/2
        w = sum(map(Fraction, row['weight_interval']))/2
        energy = float(e)+source_energies[0]
        for channel in channels.values():
            errors = np.abs(channel['energies']-energy)
            maximum_energy_error = max(maximum_energy_error, float(np.min(errors)))
            matching = errors < 1e-8
            maximum_weight_error = max(maximum_weight_error, abs(float(sum(channel['weights'][matching]))-float(w)))
    if maximum_energy_error > 2e-9 or maximum_weight_error > 2e-9:
        raise ArithmeticError('direct spectral decomposition disagrees with exact positive records')
    norms = {name: float(np.vdot(c['image'], c['image']).real) for name, c in channels.items()}
    if any(abs(value-.5) > 2e-10 for value in norms.values()):
        raise ArithmeticError('complete raw success weight is not one half')
    result = {
        'schema': 'oph.openwave.native_quark_adversarial_review.v1',
        'source_commit': PIN, 'computational_pass': True,
        'receipt_sha256': hashlib.sha256(args.receipt.read_bytes()).hexdigest(),
        'auditor_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'method': 'Direct metric-orthonormal Hamiltonian diagonalization and17 independent complex resolvent solves; different from exact polynomial/residue producer.',
        'checks': {'all_saved_records_match_direct_source_spectrum': True,
                   'both_raw_success_weights_one_half': True,
                   'uc_dc_resolvents_equal_within_declared_tolerance': True,
                   'rational_total_weight_encloses_one_half': True},
        'mutations': {'true_down_energy_translation_detected': True},
        'diagnostics': {'source_dimension': len(source_h),
                        'target_dimensions': {k: len(v['h']) for k, v in channels.items()},
                        'max_metric_Hermiticity_error': max(source_error, *(x['hermitian_error'] for x in channels.values())),
                        'max_resolvent_difference': defect, 'energy_shift_control_difference': mutation_defect,
                        'max_energy_difference': maximum_energy_error, 'max_weight_difference': maximum_weight_error,
                        'raw_success_weights': norms, 'comparison_seed': 917346},
        'verdicts': {'equal_native_up_down_measures': 'CONFIRMED',
                     'rational_receipt_spectral_record_consistency': 'CONFIRMED',
                     'distinct_physical_quark_masses_or_confinement': 'NOT_ESTABLISHED'},
        'scope': 'Independent floating audit on the actual pinned finite Gram/Hamiltonian/probe blocks; exact interval validity remains with the producer proof. Raw successful spectral measures have weight1/2, with pump failure1/2 retained in the source instrument. No physical quark mass scale or confinement follows.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print('PASS independent native quark spectral audit')


if __name__ == '__main__':
    main()
