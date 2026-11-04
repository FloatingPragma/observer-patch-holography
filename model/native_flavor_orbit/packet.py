"""Cold exact replay of the full native flavor orbit and source coefficients."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from . import algebra, correlated, spectrum

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'
PARENTS = ('model/common_source/native_family_protocol_results.json',
           'model/common_source/native_matter_attachment_results.json',
           'model/common_source/native_repair_mass_results.json')
DIRECT = ('tools/strict_json.py', 'model/requirements.txt')


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


def interface():
    return {'source_endpoints': [
                'native_matter_attachment.representation: twelve actual currents on selected M15',
                'native_repair_mass.native_records: all128 corrected C/E records',
                'native_family_record_memory: all12 old memory inputs and reset2',
                'native_gauge_family_records: complete13-outcome current instrument and all512 binary/native records'],
            'declared_readable_inputs': [
                'auxiliary spin/exterior M15 realization and three gauge-trivial family copies',
                'canonical common kinetic I45 and unit neutral Higgs norm',
                'current interaction, independent two-leg preparation domain and angle pi/5',
                'source half-turn controller, reset port2 and coherent projector feedback',
                'static mean-feedback service for the joint Higgs action assembly',
                'unconditional always/C/E/CE pulse-duty policy'],
            'observable_maps': [
                'complete current jump/no-jump effect restricted to both Higgs-component domains',
                'unconditioned stationary first-order feedback coefficient from every retained record',
                'canonical same-Higgs K45 singular spectrum with every zero and multiplicity',
                'left spectral-projector overlaps giving squared charged-current mixing moduli'],
            'empirical_coefficients_or_pole_windows_consumed': [],
            'runtime_empirical_separation_proves_target_independent_design_history': False,
            'physical_identifications_not_supplied': [
                'readable source architecture is the physical source',
                'stationary coefficient is an interacting relativistic pole mass',
                'dimensionless source step is a calibrated physical clock',
                'this auxiliary matter current acts on the localized radial defect states']}


def scope():
    return {'complete_same_current_canonical_kinetic_unitary_embedding_class': True,
            'source_gated_distinct_spectra_and_nontrivial_real_mixing_derived': True,
            'all_four_duty_policies_keep_complete_records_and_zero_directions': True,
            'all_reported_spectrum_and_mixing_intervals_are_rational_certificates': True,
            'charged_current_CP_invariant_zero_for_all_real_projector_coefficients': True,
            'complete_finite_histories_have_declared_antiunitary_grading': True,
            'finite_unread_channel_nonunitarity_on_positive_gate_domain': True,
            'full_twelve_memory_channel_has_matrix_logarithm': False,
            'history_grading_is_a_universal_physical_CP_theorem': False,
            'gauge_covariance_alone_forces_zero_CP': False,
            'all_possible_source_feedback_laws_classified': False,
            'source_selected_physical_couplings_or_particle_map': False,
            'compiled_current_instrument_native_history_executed': False,
            'static_mean_action_equals_finite_correlated_memory_channel': False,
            'interacting_relativistic_continuum_constructed': False,
            'physical_pole_masses_or_absolute_units_derived': False,
            'repair_protection_or_Einstein_dynamics_derived': False,
            'integrated_gate_accepted': False}


def _produce():
    pins = source_bindings()
    result = {'schema': 'oph.native_flavor_orbit.v1', 'all_pass': True,
              'issues': [1141, 1142, 1143], 'interface': interface(),
              'complete_current_commutant': algebra.commutant_certificate(),
              'gauge_invariant_coefficient_freedom': algebra.coefficient_freedom_certificate(),
              'fixed_repair_embedding_invariants': algebra.fixed_repair_certificate(),
              'complete_memory_theorem': algebra.memory_certificate(),
              'current_gated_joint_spectrum': spectrum.catalogue(),
              'finite_correlated_memory': correlated.certificate(),
              'scope': scope(), 'canonical_physical_premises_discharged': [],
              'source_bindings': pins}
    if pins != source_bindings():
        raise ValueError('source changed during the cold calculation')
    return result


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-native-flavor-') as temporary:
        path = Path(temporary)/'worker.json'
        environment = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        child = subprocess.run([sys.executable, '-W', 'error', '-m', 'model.native_flavor_orbit.packet',
                                '--worker', str(path)], cwd=ROOT, env=environment, check=False,
                               capture_output=True, text=True, timeout=600)
        if child.returncode:
            raise RuntimeError('cold native flavor worker failed:\n'+child.stderr)
        value = strict_json.load(path)
    if value.get('source_bindings') != pins or pins != source_bindings():
        raise ValueError('source changed across the cold calculation')
    return value


def validate(value, replay=None):
    if canonical(value) != canonical(report() if replay is None else replay):
        raise ValueError('native flavor receipt differs from the complete fresh calculation')
    if value.get('source_bindings') != source_bindings():
        raise ValueError('native flavor source manifest is incomplete or stale')


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
        validate(strict_json.load(Path(args.check)))
        print('native flavor full-current and complete-memory replay passed')


if __name__ == '__main__':
    main()
