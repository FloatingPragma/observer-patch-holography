"""Cold deterministic replay of the composite-search integer certificate."""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from tools import strict_json
from tools.build_completion_program import integration_source_paths

from ..centered_particle_spectrum import packet as inherited
from . import certificate

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = Path(__file__).resolve().parent
OUTPUT = DIRECTORY/'results.json'
TABLE = DIRECTORY/'lines.csv'
WITNESS = DIRECTORY/'witness.npz'


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def load_witness(path=WITNESS):
    with np.load(path, allow_pickle=False) as archive:
        return {name: archive[name] for name in archive.files}


def bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/charged_composite_search/README.md',
        'source_receipts': ['model/centered_particle_spectrum/results.json'],
        'independent_checks': []})
    paths.update({'model/charged_composite_search/witness.npz',
                  'model/centered_particle_spectrum/witness.npz'})
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def validate(value):
    if value.get('schema') != 'oph.charged_composite_search.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid composite-search receipt')
    if value.get('all_pass_scope') != 'exact finite composite spectrum and source sum rules; inherited custody reported separately':
        raise ValueError('incorrect certificate scope')
    s = value['spectrum']
    for name, expected in (('complete_block_lines', 138), ('strictly_positive_lines', 138),
                           ('weight_positivity_unresolved_lines', 0), ('omitted_block_lines', 0),
                           ('eligible_empirical_comparisons_added', 0)):
        if type(s[name]) is not int or s[name] != expected:
            raise ValueError('composite search candidate accounting changed')
    if s['weight_outside_reducing_block'] != '0' or len(s['rows']) != 138:
        raise ValueError('incomplete source response coverage')
    for i, row in enumerate(s['rows']):
        if row['index'] != i or row['strictly_positive_weight_proved'] is not True:
            raise ValueError('missing certified composite line')
        for name in ('energy_interval', 'weight_interval'):
            lo, hi = map(certificate.F, row[name])
            if lo <= 0 or hi < lo:
                raise ValueError('invalid positive line interval')
    for name in ('all_dark_charged_states_exhausted', 'finite_gaps_identified_as_relativistic_rest_masses',
                 'physical_particle_species_identified', 'native_detector_executed',
                 'autonomous_localized_defect_proved', 'continuum_or_gravity_conclusion'):
        if s[name] is not False:
            raise ValueError('composite search promoted beyond proved scope')
    if s['source_selection']['all_linear_CAR_responses_orthogonal'] is not True:
        raise ValueError('independence from linear probes missing')
    if s['independent_sum_rules']['independent_of_composite_eigenvectors'] is not True:
        raise ValueError('independent energy moment check missing')


def produce():
    pins = bindings()
    result = {
        'schema': 'oph.charged_composite_search.v1', 'all_pass': True,
        'all_pass_scope': 'exact finite composite spectrum and source sum rules; inherited custody reported separately',
        'interface': {
            'state': 'unique centered N8 hidden-spin-zero eta-zero ground of the open 2x2x2 Weyl/contact cube',
            'algebra': 'finite 16-mode CAR algebra; complete origin-stabilizer S=3/2 eta=1/2 N7 reducing block',
            'operation': 'A_x=(sum of hidden S_y+ over native neighbors of x) c_xdown, with A_x/3 and complete complement as a mathematical instrument',
            'clock': 'original centered source generator K=H-N/2; hbar_model=1; no SI calibration',
            'refinement': 'one finite eight-site cube at c=1/8 U=1; no volume or continuum limit'},
        'spectrum': certificate.verify(load_witness()),
        'inherited_inputs': inherited.inherited_inputs(),
        'source_bindings': pins}
    validate(result)
    if pins != bindings():
        raise ValueError('source changed during composite replay')
    return result


def report():
    pins = bindings()
    with tempfile.TemporaryDirectory(prefix='oph-composite-') as temp:
        destination = Path(temp)/'worker.json'
        environment = dict(os.environ)
        environment['OPENBLAS_NUM_THREADS'] = environment['OMP_NUM_THREADS'] = '1'
        process = subprocess.run([sys.executable, '-W', 'error', '-X',
                                  'pycache_prefix='+str(Path(temp)/'pycache'),
                                  '-m', 'model.charged_composite_search.packet',
                                  '--worker', str(destination)], cwd=ROOT, env=environment,
                                 capture_output=True, text=True)
        if process.returncode:
            raise RuntimeError('cold composite replay failed:\n'+process.stderr)
        result = strict_json.loads(destination.read_text())
    if pins != bindings() or result['source_bindings'] != pins:
        raise ValueError('source changed across cold composite replay')
    validate(result)
    return result


def table(value):
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator='\n')
    writer.writerow(('index', 'energy_lower', 'energy_upper', 'weight_lower', 'weight_upper'))
    for row in value['spectrum']['rows']:
        writer.writerow((row['index'], *row['energy_interval'], *row['weight_interval']))
    return buffer.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    group.add_argument('--worker', type=Path)
    group.add_argument('--propose', action='store_true')
    args = parser.parse_args()
    if args.propose:
        from .proposal import propose
        witness = propose()
        certificate.verify(witness)
        np.savez_compressed(WITNESS, **witness)
        print('integer proposal admitted; use --write to bind the receipt')
        return
    if args.worker:
        args.worker.write_text(canonical(produce()))
        return
    value = report()
    if args.write:
        OUTPUT.write_text(canonical(value))
        TABLE.write_text(table(value))
    else:
        retained = strict_json.loads(OUTPUT.read_text())
        validate(retained)
        if canonical(value) != canonical(retained) or TABLE.read_text() != table(value):
            raise ValueError('composite receipt or complete line table changed')
    print(json.dumps({'all_pass': True, 'positive_composite_lines': 138,
                      'independent_energy_sum_rules': True,
                      'source_bindings': len(value['source_bindings'])}))


if __name__ == '__main__':
    main()
