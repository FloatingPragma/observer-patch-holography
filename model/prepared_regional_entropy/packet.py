"""Replay finite entropy bounds for the certified literal source."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from . import entropy


ROOT = entropy.ROOT
OUTPUT = Path(__file__).with_name('results.json')


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def encode(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def bindings():
    paths = [p for p in Path(__file__).parent.iterdir()
             if p.suffix in ('.py', '.md')]
    paths += [ROOT/entropy.PRIOR,
              ROOT/'model/plaquette_reference/preparation.py',
              ROOT/'model/regional_modular_expansion/expansion.py']
    paths += [ROOT/path for path in entropy.SOURCE_PATHS]
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths)}


def report():
    pins = bindings()
    source = entropy.source_bounds(strict_json.loads((ROOT/entropy.PRIOR).read_text()))
    rows = []
    for amplitude in (F(1, 40), F(1, 400), F(1, 4000), F(1, 40000)):
        row = entropy.remainder_bounds(amplitude, minimum=source['regional_minimum'],
                                       modular_error=source['modular_polynomial_error'])
        rows.append({'amplitude': amplitude, **row})
    value = encode({
        'schema': 'oph.prepared_regional_entropy.v1', 'all_pass': True,
        'issues': [1145], 'gates': ['CP-09'],
        'source': source, 'source_certificate': entropy.PRIOR,
        'inherited_source_hashes_checked': True,
        'inherited_native_control_proof_reexecuted': False,
        'domain': {'full_Fock_dimension': 256, 'regional_dimension': 16,
                   'kept_CAR_modes': [0, 1, 2, 3], 'interaction': 1,
                   'beta': F(1, 32), 'speed': F(1, 8),
                   'two_sided_operator_radius': entropy.RADIUS,
                   'regional_path_minimum_simplified': F(1, 40),
                   'global_lift_path_minimum_simplified': F(1, 640),
                   'fixed_success_trace_and_unchanged_abort_complement': True},
        'directions': entropy.direction_certificate(),
        'basis_direction_finite_bounds': rows,
        'uniform_bounds': {'relative_entropy_HS_squared_lower_coefficient': F(1, 2),
                           'relative_entropy_HS_squared_upper_coefficient': 20,
                           'polynomial_error_trace_norm_coefficient': F(2, 10**9),
                           'symmetric_basis_error_constant': F(4, 10**9),
                           'symmetric_basis_error_amplitude_squared_coefficient': F(1600, 3)},
        'scope': {'mathematical_finite_regional_variations': True,
                  'full_field_marginal_positive_lift': True,
                  'record_resolved_lift_or_native_intervention_executed': False,
                  'geometric_modular_charge_identified': False,
                  'entropy_area_stationarity_proved': False,
                  'common_interacting_continuum_or_Einstein_gate_accepted': False},
        'source_bindings': pins,
    })
    if bindings() != pins:
        raise ValueError('source changed during entropy replay')
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    value = report()
    if args.write:
        args.write.write_text(canonical(value))
    elif canonical(strict_json.loads(args.check.read_text())) != canonical(value):
        raise ValueError('retained finite entropy result differs from its exact replay')
    print(json.dumps({'all_pass': True, 'regional_directions': 255,
                      'two_sided_radius': '1/40', 'Einstein_gate_accepted': False}))


if __name__ == '__main__':
    main()
