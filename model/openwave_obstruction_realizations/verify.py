#!/usr/bin/env python3
"""Authenticate the frozen sources and execute selected scientific controls."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
REL = Path('openwave-tuning/next/OBSTRUCTION_REPAIR')


def read(path):
    return json.loads(path.read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def custody(root, manifest):
    for name, item in manifest['files'].items():
        p = Path(name)
        require(not p.is_absolute() and '..' not in p.parts, 'unsafe export path')
        data = (root / p).read_bytes()
        require(len(data) == item['bytes'], f'byte count changed: {name}')
        require(hashlib.sha256(data).hexdigest() == item['sha256'],
                f'source changed: {name}')


def scientific_checks(root):
    base = root / REL
    locality = read(base / 'continuum_locality_bound/RESULT.json')
    require(locality['exact_envelope_inequalities_pass'], 'kernel inequalities')
    bound = Fraction(locality['cross_observer_norm_upper_exact'])
    c3 = Fraction(locality['full_three_dimensional_C3_rational_upper'])
    require(bound == 8*c3*9/Fraction(100000)**3, 'complete-mode locality bound')
    cell = read(base / 'cell_entropy_limit/RESULT.json')
    require(cell['all_modes_retained'], 'cell modes omitted')
    require(max(x['entropy_identity_error'] for x in cell['records']) < 1e-10,
            'finite cell entropy identity')
    require(not cell['continuum_error_rate_certified'], 'unproved cell error rate')
    consumers = read(base / 'vacuum_trace_review/CONSUMERS.json')
    require(len(consumers['consumers']) == 4, 'vacuum cases omitted')
    for row in consumers['consumers']:
        require(row['larger_cut_eigen_count_upper'] < 2, 'vacuum eigen count')
        require(Fraction(row['neutral_gap_lower_exact']) == Fraction(1, 6),
                'neutral gap')
    sharp = next(x for x in consumers['consumers'] if x['case'] == 'zero_K4')
    require(sharp['ground_fidelity_lower'] > .9999987, 'vacuum fidelity')
    require(sharp['actual_vacuum_Green_error_upper'] < .000794, 'vacuum response')
    charge = read(base / 'row1_v3/RESULT.json')
    require(charge['registered_order'] == list(range(7, 30))+list(range(7)),
            'charge policy')
    require(charge['computed_joint_normalizations'] == 4788
            and charge['all_computed_normalizations_one'], 'complete charge laws')
    require(sum(charge['pointer_outcomes'].values()) ==
            charge['coverage']['regional_reads'], 'pointer outcomes omitted')
    require(not charge['channel_completeness_authenticated']
            and not charge['physical_charge_identification_solved'], 'charge scope')
    spin = read(base / 'BARYON_SPIN_REVIEW.json')
    require(spin['all_checks_pass'] and spin['exact_matrix_certificate_count'] == 6,
            'full-source spin certificates')
    mobility = read(base / 'baryon_mobility/RESULT.json')
    require(mobility['complete_product_dimension'] == 8464
            and mobility['fermion_parity_retained']
            and not mobility['full625_spectral_transfer'], 'composite-motion scope')
    ledger = read(base / 'ISSUES.json')
    require(len(ledger['failed_constructions']) >= 46, 'failed constructions omitted')
    require(not ledger['official_rows_changed'], 'unearned platform award')
    result = {'kernel_bound': str(bound), 'cell_controls': len(cell['records']),
            'vacuum_cases': len(consumers['consumers']),
            'charge_laws': charge['computed_joint_normalizations'],
            'spin_certificates': spin['exact_matrix_certificate_count'],
            'singlet_product_dimension': mobility['complete_product_dimension'],
            'failed_constructions': len(ledger['failed_constructions'])}
    neutral = base/'species_leptons/neutral_paths/EXACT_POLE_INTERVALS.json'
    if neutral.exists():
        poles = read(neutral)
        lower, upper = map(Fraction, poles['splitting_ratio_interval_exact'])
        require(Fraction(33, 10000) < lower <= upper < Fraction(34, 10000),
                'scale-independent neutral exclusion')
        require(len(poles['all63_actual_pole_intervals_exact']) == 63,
                'neutral poles omitted')
        current = read(neutral.parent/'CURRENT_CERTIFICATE.json')
        require(current['all63_bands_and_18_records_retained']
                and Fraction(current['complete_five_outcome_record_error_exact'])
                    < Fraction(48, 10**24), 'complete neutral current records')
        result['neutral_pole_intervals'] = 63
        result['neutral_ratio_excludes_overall_scale_repair'] = True
    gravity = base/'gravity/joint_continuum_entropy/L40_h0p5/RESULT.json'
    if gravity.exists():
        selection = read(gravity)
        errors = []
        for solution in selection['solutions']:
            require(solution['rank'] == 6, 'gravity fit rank')
            require(max(abs(x['total_residual']) for x in solution['training']) < 1e-12,
                    'gravity training equations')
            errors.append(max(abs(x['relative_residual']) for x in solution['held_out']))
        require(len(errors) == 2 and max(errors) > .31 and min(errors) > .09,
                'held-out geometric failure omitted')
        result['gravity_held_out_residuals'] = errors
    return result


def tamper_control():
    with tempfile.TemporaryDirectory(prefix='oph-custody-control-') as tmp:
        root = Path(tmp)
        (root/'source.py').write_bytes(b'complete source\n')
        data = (root/'source.py').read_bytes()
        manifest = {'files': {'source.py': {'bytes': len(data),
                     'sha256': hashlib.sha256(data).hexdigest()}}}
        custody(root, manifest)
        (root/'source.py').write_bytes(b'altered  source\n')
        try:
            custody(root, manifest)
        except ValueError:
            return True
        raise ValueError('altered source accepted')


def replay(group, manifest):
    work = Path(tempfile.mkdtemp(prefix='oph-obstruction-replay-'))
    shutil.copytree(HERE/'inputs', work/'inputs')
    base = work/'inputs'/REL
    commands = []
    if group in ('quick', 'all'):
        commands += ['continuum_locality_bound/check.py', 'cell_entropy_limit/check.py',
                     'vacuum_trace_review/consumers.py']
    if group in ('trace', 'all'):
        commands += ['vacuum_trace_review/run.py', 'vacuum_trace_review/consumers.py']
    if group in ('charge', 'all'):
        commands += ['row1_v3/verify.py']
    if group in ('additions', 'all'):
        commands += ['species_leptons/neutral_paths/exact_endpoints.py',
                     'species_leptons/neutral_paths/current_certificate.py',
                     'gravity/joint_continuum_entropy/review.py']
    environment = os.environ.copy()
    environment.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
                       MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    receipts = []
    print(f'Replay directory: {work}', flush=True)
    for index, relative in enumerate(commands):
        script = base/relative
        log = work/f'{index}-{script.parent.name}-{script.stem}.log'
        print(f'Executing {relative}', flush=True)
        with log.open('w') as stream:
            result = subprocess.run([sys.executable, str(script)], cwd=script.parent,
                                    env=environment, stdout=stream,
                                    stderr=subprocess.STDOUT, timeout=3600)
        require(result.returncode == 0, f'replay failed; see {log}')
        receipts.append({'script': relative, 'exit_code': result.returncode,
                         'log': str(log)})
    checks = scientific_checks(work/'inputs')
    custody(HERE/'inputs', manifest)
    return {'directory': str(work), 'commands': receipts, 'checks': checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', choices=('quick', 'trace', 'charge', 'additions', 'all'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if (HERE/'MANIFEST.json').exists():
        custody(HERE, read(HERE/'MANIFEST.json'))
    manifest = read(HERE/'EXPORT_MANIFEST.json')
    require(not manifest['missing_pins'] and not manifest['changed_pins'],
            'unresolved source custody')
    custody(HERE/'inputs', manifest)
    result = {'schema': 'oph.obstruction_realizations.verification.v1',
              'export_manifest_sha256': hashlib.sha256(
                  (HERE/'EXPORT_MANIFEST.json').read_bytes()).hexdigest(),
              'files_checked': len(manifest['files']),
              'mathematical_receipts': scientific_checks(HERE/'inputs'),
              'changed_byte_control_rejected': tamper_control(),
              'physical_rows_solved': False, 'official_rows_changed': False}
    if args.replay:
        result['replay'] = replay(args.replay, manifest)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
