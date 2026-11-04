"""Cold reproducible mass-shell calculation with unchanged parent custody."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import model
from tools import strict_json
from . import algebra, bounds, exterior, units

ROOT = Path(model.__file__).resolve().parent.parent
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'
PARENTS = (
    'model/native_mass_channels/results.json',
    'model/bounded_defect_index/results.json',
    'model/radial_interacting_core/results.json',
)
OWN_FILES = (
    '__init__.py', '__main__.py', 'algebra.py', 'bounds.py', 'exterior.py',
    'units.py', 'packet.py', 'test_dispersion.py', 'test_exterior.py',
    'test_units.py', 'test_packet.py', 'PROOF.md', 'CALIBRATION_PROOF.md', 'README.md',
)
DIRECT_SOURCES = (
    'model/relativistic_sector_audit/sector.py',
    'model/contact_scattering/onshell.py',
    'model/periodic_momentum/momentum.py',
    'model/bounded_defect_index/finite.py',
    'model/bounded_defect_index/PROOF.md',
    'model/radial_interacting_core/core.py',
    'model/radial_interacting_core/PROOF.md',
)


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def _path(name):
    prefix = 'model/source_mass_shell/'
    return PACKAGE/name[len(prefix):] if name.startswith(prefix) else ROOT/name


def source_bindings():
    inherited = {}
    for name in PARENTS:
        parent = strict_json.load(ROOT/name)
        if parent.get('all_pass') is not True or type(parent.get('source_bindings')) is not dict:
            raise ValueError('source parent receipt has no complete passing manifest')
        for path, digest in parent['source_bindings'].items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
                raise ValueError('source parent has a stale source binding: '+path)
            if path in inherited and inherited[path] != digest:
                raise ValueError('parent source manifests disagree')
            inherited[path] = digest
    paths = set(inherited) | set(PARENTS) | set(DIRECT_SOURCES)
    paths.update('model/source_mass_shell/'+name for name in OWN_FILES)
    return {name: hashlib.sha256(_path(name).read_bytes()).hexdigest() for name in sorted(paths)}


def numerical_admission():
    from fractions import Fraction as F
    data = exterior.diagnostic(8, 10000, (F(1, 3), F(1, 4), F(1, 2)), mass=64)
    for key in ('norm_error', 'fiber_eigenvector_error', 'kinetic_envelope_residual',
                'texture_residual', 'actual_radial_residual'):
        if type(data[key]) is not float or not math.isfinite(data[key]) or data[key] < 0:
            raise ArithmeticError('nonfinite, negative or untyped numerical diagnostic')
    if max(data['norm_error'], data['fiber_eigenvector_error']) > 1e-10:
        raise ArithmeticError('floating packet norm or fiber check failed')
    for measured, bound in [('kinetic_envelope_residual', 'kinetic_envelope_error_upper'),
                             ('texture_residual', 'radial_texture_error_upper'),
                             ('actual_radial_residual', 'lattice_eigenvector_residual_upper')]:
        if data[measured] > float(F(data['bounds'][bound])):
            raise ArithmeticError('literal radial packet exceeded its proved residual bound')
    return {'cells_per_axis': 8, 'displacement': 10000,
            'support_sites': 4096, 'evaluated_rows_including_halo': 11664,
            'norm_and_fiber_admission_upper': '1/10000000000',
            'norm_and_fiber_checks_pass': True,
            'literal_radial_residuals_below_proved_bounds': True,
            'raw_platform_float_diagnostics_promoted_to_exact_evidence': False,
            'all_source_packet_and_halo_rows_retained': True}


def report():
    pins = source_bindings()
    value = {
        'schema': 'oph.source_mass_shell.v1', 'all_pass': True,
        'algebra': algebra.clifford_certificate(),
        'dispersion': bounds.certificate(),
        'source_packet_control': numerical_admission(),
        'physical_units': units.certificate(),
        'scope': {
            'same_bounded_radial_one_body_exterior_packets': True,
            'all_eight_valleys_and_sixteen_fiber_rows': True,
            'supplied_mass_and_profile': True,
            'original_one_body_clock_sign_and_origin': True,
            'asymptotic_bulk_Weyl_sequences': True,
            'light_speed_matching_is_explicit_physical_identification': True,
            'radial_global_translation_or_momentum_symmetry': False,
            'moving_localized_defect_constructed': False,
            'time_dependent_radial_scattering_wave_operators_constructed': False,
            'interacting_radial_charged_poles_identified_with_bulk_shell': False,
            'interacting_relativistic_continuum': False,
            'one_common_massless_bulk_and_gapped_defect_regime': False,
            'finite_strong_mass_core_continued_to_zero_spacing': False,
            'absolute_time_unit_or_GeV_mass_derived': False,
            'physical_species_identified': False,
            'native_repair_protection_proved': False,
        },
        'physical_comparison_count_added': 0,
        'integrated_gate_accepted': False,
        'source_bindings': pins,
    }
    if source_bindings() != pins:
        raise ValueError('source files changed during the calculation')
    return value


def validate(value):
    if canonical(value) != canonical(report()):
        raise ValueError('mass-shell receipt differs from the full fresh calculation')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    operation = parser.add_mutually_exclusive_group(required=True)
    operation.add_argument('--write', nargs='?', const=str(OUTPUT))
    operation.add_argument('--check', nargs='?', const=str(OUTPUT))
    args = parser.parse_args()
    if args.write:
        path = Path(args.write)
        path.write_text(canonical(report()))
        print('wrote '+str(path))
    else:
        validate(strict_json.load(Path(args.check)))
        print('source mass-shell replay passed')


if __name__ == '__main__':
    main()
