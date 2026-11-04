"""Cold, source-bound receipt for the conditional Majorana action join."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from . import native, readout, spectrum, topology

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'
PARENTS = ('model/common_source/native_matter_attachment_results.json',
           'model/common_source/calibrated_matter_results.json',
           readout.RECEIPT)
DIRECT = ('model/conditional_standard_model/spectrum.py',
          'model/conditional_standard_model/defect_constraints.py',
          'model/conditional_standard_model/PROOF.md',
          'model/common_source/test_calibrated_matter.py',
          'model/common_source/native_repair_mass_results.npz',
          'tools/strict_json.py', 'model/requirements.txt')


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    inherited = {}
    for name in PARENTS:
        record = strict_json.load(ROOT/name)
        bindings = record.get('source_bindings', record.get('source_sha256'))
        if type(bindings) is not dict or not bindings:
            raise ValueError('parent lacks a complete source manifest')
        for path, digest in bindings.items():
            target = (ROOT/path).resolve()
            if not target.is_relative_to(ROOT.resolve()):
                raise ValueError('source binding escapes repository')
            if hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                raise ValueError('stale inherited source: '+path)
            if path in inherited and inherited[path] != digest:
                raise ValueError('inconsistent inherited source: '+path)
            inherited[path] = digest
    paths = set(inherited) | set(PARENTS) | set(DIRECT)
    paths.update(str(path.relative_to(ROOT)) for path in PACKAGE.iterdir()
                 if path.is_file() and path.suffix in ('.py', '.md'))
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def scope():
    return {'conditional_native_current_and_spin_exterior_architecture': True,
            'canonical_Weyl_kinetic_and_same_Higgs_action_are_inputs': True,
            'Majorana_coefficient_and_empirical_calibration_are_inputs': True,
            'all_45_internal_and_90_CAR_and_180_Nambu_rows_retained': True,
            'all_neutrino_ranks_zero_through_three_admitted': True,
            'connected_PMNS_required_for_generation_charge_identification': True,
            'massless_exclusion_applied_to_calibrated_zero_lightest_mass': True,
            'source_selects_physical_couplings_or_Majorana_phases': False,
            'interacting_continuum_and_physical_pole_masses_constructed': False,
            'physical_defect_family_assignment_constructed': False,
            'exact_topological_charge_preserved_by_native_repair_proved': False,
            'Einstein_dynamics_derived': False,
            'integrated_gate_accepted': False}


def _produce():
    pins = source_bindings()
    result = {'schema': 'oph.majorana_defect_bridge.v1', 'all_pass': True,
              'issues': [1141, 1142, 1143],
              'native_current_Higgs_dictionary': native.certificate(),
              'concrete_native_coefficient_readout': readout.certificate(),
              'actual_action_spectrum': spectrum.catalogue(),
              'full_momentum_symbol': spectrum.universal_symbol_certificate(),
              'additive_topological_constraints': topology.certificate(),
              'scope': scope(), 'canonical_physical_premises_discharged': [],
              'source_bindings': pins}
    if pins != source_bindings():
        raise ValueError('source changed during the cold calculation')
    return result


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-majorana-defect-') as temporary:
        path = Path(temporary)/'worker.json'
        environment = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        child = subprocess.run([sys.executable, '-W', 'error', '-m', 'model.majorana_defect_bridge.packet',
                                '--worker', str(path)], cwd=ROOT, env=environment, check=False,
                               capture_output=True, text=True, timeout=600)
        if child.returncode:
            raise RuntimeError('cold Majorana worker failed:\n'+child.stderr)
        value = strict_json.load(path)
    if value.get('source_bindings') != pins or pins != source_bindings():
        raise ValueError('source changed across the cold calculation')
    return value


def validate(value, replay=None):
    if canonical(value) != canonical(report() if replay is None else replay):
        raise ValueError('Majorana bridge differs from the complete fresh calculation')
    if value.get('source_bindings') != source_bindings():
        raise ValueError('Majorana bridge source manifest is incomplete or stale')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choices = parser.add_mutually_exclusive_group(required=True)
    choices.add_argument('--write', nargs='?', const=str(OUTPUT))
    choices.add_argument('--check', nargs='?', const=str(OUTPUT))
    choices.add_argument('--worker')
    args = parser.parse_args()
    if args.worker:
        Path(args.worker).write_text(canonical(_produce()), encoding='utf-8')
    elif args.write:
        Path(args.write).write_text(canonical(report()), encoding='utf-8')
        print('wrote '+args.write)
    else:
        validate(strict_json.load(args.check))
        print('Majorana source action and topology replay passed')


if __name__ == '__main__':
    main()
