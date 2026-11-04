"""Source-bound replay of the complete finite defect/current/mass embedding."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import algebra, controls

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'
PARENTS = ('model/radial_defect_gap/results.json',
           'model/common_source/native_matter_attachment_results.json',
           'model/mcpr_native_memory_bridge/results.json')
DIRECT = ('tools/strict_json.py', 'model/requirements.txt',
          'code/P_derivation/paper_math.py', 'model/conditional_standard_model/spectrum.py',
          'model/common_source/calibrated_matter.py')


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    inherited = {}
    for name in PARENTS:
        value = strict_json.load(ROOT/name)
        pins = value.get('source_bindings', value.get('source_sha256'))
        if not isinstance(pins, dict) or not pins:
            raise ValueError('parent must provide complete source custody')
        for path, digest in pins.items():
            target = (ROOT/path).resolve()
            if not target.is_relative_to(ROOT.resolve()):
                raise ValueError('source binding escapes repository')
            if hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                raise ValueError('stale inherited source: '+path)
            if path in inherited and inherited[path] != digest:
                raise ValueError('conflicting inherited source: '+path)
            inherited[path] = digest
    paths = set(inherited) | set(PARENTS) | set(DIRECT)
    paths.update(str(path.relative_to(ROOT)) for path in PACKAGE.iterdir()
                 if path.is_file() and path.suffix in ('.py', '.md'))
    paths.update(integration_source_paths(ROOT, {
        'domain_and_proof': 'model/defect_current_bridge/PROOF.md',
        'source_receipts': list(PARENTS), 'independent_checks': []}))
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def scope():
    return {'chosen_internal_carrier_and_pairing_are_explicit_design_inputs': True,
            'same_defect_CAR_states_native_currents_and_memory_mass_operator': True,
            'every_source_site_spin_family_and_Nambu_direction_retained': True,
            'Nambu_holes_are_redundancy_not_extra_physical_mirror_species': True,
            'finite_quadratic_dynamics_and_even_state_map_exact': True,
            'chosen_full_record_pulse_service_has_controlled_mean_Hamiltonian_limit': True,
            'single_quasiparticle_core_band_isolation_proved': True,
            'P_dependent_lepton_formula_realized_on_same_core': True,
            'historical_MCPR_architecture_target_independence_proved': False,
            'three_families_forced_by_topological_index': False,
            'fermion_creation_changes_texture_winding': False,
            'unmodified_site_repair_preserves_core': False,
            'full_interacting_Fock_band_isolation_proved': False,
            'relativistic_dispersion_or_pole_mass_identification_proved': False,
            'physical_particle_species_identified': False,
            'local_Gauss_constraint_or_dynamical_gauge_links_derived': False,
            'physical_clock_or_GeV_calibration_derived': False,
            'massless_interacting_continuum_constructed': False,
            'Einstein_geometry_constructed': False,
            'integrated_gate_accepted': False}


def produce():
    pins = source_bindings()
    result = {'schema': 'oph.defect_current_bridge.v1', 'all_pass': True,
              'issues': [1139, 1140, 1141, 1142, 1143],
              'design': {'physical_carrier': 'CAR(ell2(Lambda) tensor C2 tensor M15 tensor F3)',
                         'three_family_copies': 'declared',
                         'texture': 'supplied radial degree-one background',
                         'energy_unit': 'chosen v(P); radial hopping v/a=1 and texture mass64 in this unit',
                         'lepton_mass_operator': 'chosen coherent pairing service implements the MCPR native C-duty stationary-mean matrix times its source gain',
                         'quark_operator': 'declared Stage5 epsilon=1/6, powers (6,3,0) and (6,4,2)',
                         'neutral_couplings': 'three minimal neutrino directions set to zero',
                         'finite_memory_channel_equals_quadratic_Hamiltonian': False},
              'exact_identities': algebra.certificate(),
              'complete_record_pulse_service': algebra.pulse_service_certificate(),
              'source_cases': [controls.source_case(branch, index)
                               for branch in ('forward', 'measured_endpoint') for index in (0, 1)],
              'independent_countercontrols': controls.independent_control(),
              'scope': scope(), 'canonical_physical_premises_discharged': [],
              'source_bindings': pins}
    if pins != source_bindings():
        raise ValueError('source changed during calculation')
    return result


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-defect-current-') as temporary:
        path = Path(temporary)/'worker.json'
        environment = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        child = subprocess.run([sys.executable, '-W', 'error', '-m', 'model.defect_current_bridge.packet',
                                '--worker', str(path)], cwd=ROOT, env=environment, check=False,
                               capture_output=True, text=True, timeout=600)
        if child.returncode:
            raise RuntimeError('cold defect-current worker failed:\n'+child.stderr)
        value = strict_json.load(path)
    if value.get('source_bindings') != pins or pins != source_bindings():
        raise ValueError('source changed across calculation')
    return value


def validate(value, replay=None):
    if canonical(value) != canonical(report() if replay is None else replay):
        raise ValueError('defect-current receipt differs from complete replay')
    if value.get('source_bindings') != source_bindings():
        raise ValueError('defect-current manifest is incomplete or stale')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choices = parser.add_mutually_exclusive_group(required=True)
    choices.add_argument('--write', nargs='?', const=str(OUTPUT))
    choices.add_argument('--check', nargs='?', const=str(OUTPUT))
    choices.add_argument('--worker')
    args = parser.parse_args()
    if args.worker:
        Path(args.worker).write_text(canonical(produce()), encoding='utf-8')
    elif args.write:
        Path(args.write).write_text(canonical(report()), encoding='utf-8')
        print('wrote '+args.write)
    else:
        validate(strict_json.load(args.check))
        print('complete defect-current-mass replay passed')


if __name__ == '__main__':
    main()
