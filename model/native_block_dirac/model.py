"""Produce exact encoding, clock, complete records and source custody."""
import argparse
import ast
import hashlib
import json
from pathlib import Path

from .source import ROOT, HERE, PARENT, PARENT_SHA, PREPARATION, admit, constituent_bounds
from .algebra import four_mode_code, spinor_checks, complete_two_tests
from .controller import exact_noncommuting_pulse_control
from .reference import spatial_group, complete_lift_control
from .resources import finite_member
from model.particle_common_contract.provenance import import_closure

RESEARCH_REFERENCE = HERE/'research_reference.json'
RESEARCH_SHA = 'df277440b874d997024a9446f4d84c2c086f572dccfbb283c39bc78861a87fec'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_pins():
    parent, _ = admit()
    inherited = dict(parent['input_sha256']); inherited[PARENT] = PARENT_SHA
    inherited[PREPARATION] = digest(ROOT/PREPARATION)
    seeds = set(inherited)
    seeds.update({'model/particle_common_contract/provenance.py', 'model/requirements.txt'})
    seeds.update(str(path.relative_to(ROOT)) for path in HERE.iterdir() if path.suffix in ('.py', '.md'))
    seeds.add(str(RESEARCH_REFERENCE.relative_to(ROOT)))
    for path in HERE.glob('*.py'):
        for node in ast.walk(ast.parse(path.read_text())):
            names = [alias.name for alias in node.names] if isinstance(node, ast.Import) else (
                [node.module] if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module else [])
            for name in names:
                if name.startswith('model.'):
                    stem = name.replace('.', '/')
                    for candidate in (stem+'.py', stem+'/__init__.py'):
                        if (ROOT/candidate).is_file(): seeds.add(candidate)
    pins = {name: digest(ROOT/name) for name in sorted(import_closure(seeds))}
    for name, expected in inherited.items():
        if pins[name] != expected: raise ValueError('source pin changed: '+name)
    return pins


def local_pins():
    return {path.name: digest(path) for path in sorted(HERE.iterdir()) if path.suffix in ('.py', '.md')}


def validate_pins(value):
    for mapping, root in ((value['input_sha256'], ROOT), (value['local_sha256'], HERE)):
        for name, expected in mapping.items():
            relative = Path(name)
            if relative.is_absolute() or '..' in relative.parts or digest(root/name) != expected:
                raise ValueError('source custody mismatch: '+name)


def certificate():
    result = {'schema': 'oph.native_block_dirac.v1',
            'input_sha256': source_pins(), 'local_sha256': local_pins(),
            'parent_native_join_receipt_sha256': PARENT_SHA,
            'complete_native_constituent': constituent_bounds(),
            'canonical_occupation_encoding': four_mode_code(),
            'Dirac_frame_and_boost_algebra': spinor_checks(),
            'finite_spatial_reference': spatial_group(),
            'complete_reference_lift_control': complete_lift_control(),
            'noncommuting_pulse_clock_control': exact_noncommuting_pulse_control(),
            'literal_complete_active_circuit': complete_two_tests(),
            'literal_zero_separation_control': complete_two_tests(False),
            'finite_member': finite_member(),
            'field_comparison': {
                'fixed_window_frequency_error': '0<=omega-Omega_a<=a^2 P^4/(24m)',
                'finite_volume_field_error': 'sqrt(3) h ||grad f||+(3h/2)(2/m+T)||f||+T a^2(P+sqrt(3)h/2)^4||f||/(24m)',
                'complete_record_TV_error': 'epsilon_preparation+P(all attempts fail)+sum_j |theta_j| epsilon_fj/Cj+sum_j tau_j eta_controller,j',
                'same_source_state_reference_and_clock': True,
                'unconditional_finite_refinement_uses_retained_attempt_amplification': True,
                'continuum_field': 'free massive Dirac CAR field with induced spinor action and causal anticommutator',
                'microscopic_square_root_spatial_couplers': 'full range on each connected finite lattice',
                'finite_momentum_cutoff_strict_microcausality': False,
                'interacting_continuum_or_electromagnetic_Gauss_law_claimed': False,
                'SM_weak_chirality_or_observed_mass_identity_claimed': False,
                'mass_parameter_is_empirically_calibrated': False,
                'passive_native_controller_or_efficient_full_source_execution_claimed': False}}
    if digest(RESEARCH_REFERENCE) != RESEARCH_SHA:
        raise ValueError('the immutable independent free Dirac receipt changed')
    ignored = {'schema', 'input_sha256', 'local_sha256'}
    scientific = lambda value: {k: v for k, v in value.items() if k not in ignored}
    if scientific(result) != scientific(json.loads(RESEARCH_REFERENCE.read_text())):
        raise ValueError('public free Dirac science differs from frozen research')
    result['independent_research_comparison'] = {
        'reference_sha256': RESEARCH_SHA, 'all_scientific_fields_equal': True,
        'ignored_only_provenance_fields': sorted(ignored), 'fresh_calculation_was_replayed': True}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true'); mode.add_argument('--check', action='store_true')
    args = parser.parse_args(); path = HERE/'results.json'
    if args.check:
        saved = json.loads(path.read_text()); validate_pins(saved)
    result = certificate()
    if args.write: path.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    elif result != saved: raise ValueError('native free Dirac receipt differs')
    print('PASS: actual native CAR code, driven same-clock records, free field comparison and complete source custody')


if __name__ == '__main__':
    main()
