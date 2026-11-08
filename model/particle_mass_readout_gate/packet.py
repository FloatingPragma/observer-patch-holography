"""Generate a free-branch census while physical particle mass rows stay empty."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from model.continuum_free_limit import packet as free
from model.spin_canonical_energy import packet as spin
from . import algebra, bounds, charge_bridge

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'table.json'
MARKDOWN = HERE/'table.md'
PARENTS = ('model/continuum_free_limit/results.json',
           'model/spin_canonical_energy/results.json')


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/particle_mass_readout_gate/README.md',
        'source_receipts': list(PARENTS),
        'independent_checks': ['model/protected_spin_limit/README.md',
                               'model/strong_coupling_refinement/README.md'],
    })
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def checked_parents():
    free_value = strict_json.load(ROOT/PARENTS[0])
    free.validate(free_value)
    if canonical(free_value) != canonical(free.report()):
        raise ValueError('retained free-limit producer or proof changed')
    spin_value = strict_json.load(ROOT/PARENTS[1])
    spin.validate(spin_value)
    return free_value, spin_value


def build():
    before = source_bindings()
    free_value, spin_value = checked_parents()
    f, s = free_value['theorem'], spin_value['theorem']
    if (f['massless_Weyl_field_count'] != 8 or f['rest_masses_in_this_branch'] != 'all zero'
            or f['physical_species_or_unit_dictionary_identified'] is not False
            or s['physical_species_or_particle_masses_identified'] is not False):
        raise ValueError('parent particle-identification scope changed')
    rows = [{'kinetic_valley': row['valley'],
             'charge_branches': ['particle_addition', 'hole'],
             'scope': 'free model one-particle and hole subspaces only',
             'dispersion': 'E=v|p|',
             'invariant_rest_energy_squared': '0',
             'empirical_species': None,
             'physical_particle_identified': False}
            for row in algebra.casimir_checks()]
    result = {
        'schema': 'oph.particle_mass_readout_gate.v1',
        'title': 'Free-branch census',
        'source_bindings': before,
        'parent_theorems': {PARENTS[0]: f, PARENTS[1]: s},
        'parent_complete_producers_replayed': list(PARENTS),
        'all_ancestral_expensive_producers_replayed': False,
        'free_branch_convention': {
            'state': f['reference'], 'refinement': f['refinement'],
            'clock': 'K=H-(U/2)N in its own centered ground GNS representation',
            'momentum': 'common continuous spatial translations of all eight continuum valleys',
            'units': 'hbar=1 and declared v>0; no physical unit calibration',
            'one_particle_and_hole_scope_only': True,
            'entire_Fock_space_has_zero_invariant_mass': False,
        },
        'free_branch_census': rows,
        'physical_particle_mass_table': [],
        'mass_ratios': [],
        'zero_over_zero_ratio_defined': False,
        'ordered_spin_branch': {
            'clock': s['scaled_clock'],
            'mass_readout_eligible': False,
            'invariant_rest_energy_squared': None,
            'sharp_mass_shell_proved': s['sharp_mass_shell_proved'],
            'Lorentz_covariance_proved': s['Lorentz_covariance_proved'],
            'merged_with_free_branch': False,
        },
        'physical_admission': {
            'source_capacity_selects_physical_coupling_trajectory': False,
            'controlled_native_preparation_and_complete_detector_join': False,
            'empirical_species_dictionary_and_unit_calibration': False,
            'physical_particle_prediction_admitted': False,
            'finite_H_B_selection_gap_used_as_rest_mass': False,
        },
        'protected_trajectory_clock_controls': [bounds.charged_clock_bound(e) for e in (2, 4, 8, 16)],
        'finite_energy_charged_weight_controls': [bounds.charged_weight_bound(e) for e in (2, 4, 8, 16)],
        'finite_charge_operator_bridge': charge_bridge.certificate(),
        'verification_checks_pass': True,
    }
    if source_bindings() != before:
        raise RuntimeError('census inputs changed during replay')
    return result


def validate(value):
    if type(value) is not dict or value.get('schema') != 'oph.particle_mass_readout_gate.v1':
        raise ValueError('invalid free-branch census')
    if canonical(value) != canonical(build()):
        raise ValueError('census, source, convention or physical admission changed')


def render(value):
    lines = ['# Free-branch census', '',
             '**Physical particle mass table: empty.**', '',
             'The rows are kinetic valleys of the declared centered free model, not identified particles of nature.',
             'Both particle-addition and hole branches are retained. The zero parameter applies only to their',
             'one-particle subspaces. It is not a rest-frame measurement or a statement about multiparticle states.', '',
             '| Kinetic valley | Free dispersion | Invariant rest-energy squared (model control) | Physical species |',
             '| --- | --- | --- | --- |']
    for row in value['free_branch_census']:
        valley = '('+','.join(map(str, row['kinetic_valley']))+')'
        lines.append(f"| {valley} | {row['dispersion']} | {row['invariant_rest_energy_squared']} | Unidentified |")
    lines += ['', 'Convention: centered ground, K=H-(U/2)N, hbar=1, declared v>0. No physical unit calibration.',
              'The ordered-spin branch has no admitted relativistic mass readout and is not combined with these rows.',
              'Ratios of the zero control parameters are undefined. Finite selection gaps are not entered as masses.', '',
              'Generated by `python -m model.particle_mass_readout_gate --write`.', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    mode.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    if args.check is not None:
        value = strict_json.load(args.check)
        validate(value)
        if MARKDOWN.read_text() != render(value):
            raise ValueError('generated census markdown differs')
    else:
        value = build()
        args.write.write_text(canonical(value))
        if args.write.resolve() == OUTPUT:
            MARKDOWN.write_text(render(value))
    print('Free-branch census verified; physical particle mass table empty')


if __name__ == '__main__':
    main()
