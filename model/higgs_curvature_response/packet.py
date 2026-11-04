"""Replay the same-Higgs metric, Wald and generalized-volume calculation."""
import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from . import balls, source, ward

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'
PARENTS = ('model/common_source/native_matter_attachment_results.json',
           'model/common_source/calibrated_matter_results.json')
DIRECT = ('model/__init__.py', 'model/requirements.txt', 'tools/strict_json.py',
          'model/conditional_standard_model/__init__.py',
          'model/conditional_standard_model/spectrum.py',
          'model/common_source/calibrated_matter.py',
          'model/matched_einstein_limit/PROOF.md')


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    inherited = {}
    for name in PARENTS:
        parent = strict_json.load(ROOT/name)
        bindings = parent.get('source_bindings', parent.get('source_sha256'))
        if type(bindings) is not dict or not bindings:
            raise ValueError('parent receipt lacks complete source custody')
        for path, digest in bindings.items():
            target = (ROOT/path).resolve()
            if not target.is_relative_to(ROOT) or hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                raise ValueError('stale or external parent source: '+path)
            if path in inherited and inherited[path] != digest:
                raise ValueError('parent source manifests disagree')
            inherited[path] = digest
    paths = set(inherited)|set(PARENTS)|set(DIRECT)
    paths.update(str(path.relative_to(ROOT)) for path in PACKAGE.iterdir()
                 if path.is_file() and path.suffix in ('.py', '.md'))
    if any(name.startswith(('tracking/', 'model/completion/')) for name in paths):
        raise ValueError('central registers cannot enter this scientific source closure')
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def report():
    bindings = source_bindings()
    value = {
        'schema': 'oph.higgs_curvature_response.v1',
        'all_pass': True,
        'same_higgs_action': source.certificate(),
        'metric_variation': ward.metric_variation_certificate(),
        'wald_and_generalized_volume': ward.wald_certificate(),
        'curved_force_balance': ward.curved_ward_certificate(),
        'regional_first_variation': balls.certificate(),
        'scope': {
            'same_conditional_SM_Higgs_and_Yukawa_tensors': True,
            'curved_extension_and_nonminimal_coupling_specified': True,
            'metric_stress_and_Wald_coefficient_from_same_action': True,
            'generalized_volume_from_same_curvature_derivative': True,
            'exact_first_variation_around_flat_constant_F_reference': True,
            'explicit_fourth_derivative_remainder_hypotheses': True,
            'conditional_linearized_scalar_tensor_entropy_join': True,
            'ordinary_fixed_volume_premise_instantiated': False,
            'source_regional_entropy_or_stationarity_constructed': False,
            'finite_parity_entropy_promoted_to_Wald_entropy': False,
            'scalar_stress_conserved_without_gauge_or_Yukawa_exchange': False,
            'full_fermionic_quantum_Ward_identity_proved': False,
            'native_curved_metric_and_common_proper_clock_constructed': False,
            'nonlinear_Einstein_equation_derived': False,
            'curvature_coupling_or_Newton_scale_source_selected': False,
        },
        'canonical_physical_premises_discharged': [],
        'physical_comparison_count_added': 0,
        'integrated_gate_accepted': False,
        'source_bindings': bindings,
    }
    if source_bindings() != bindings:
        raise ValueError('source changed during the calculation')
    return value


def validate(value):
    if canonical(value) != canonical(report()):
        raise ValueError('Higgs curvature receipt differs from the complete fresh calculation')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    operations = parser.add_mutually_exclusive_group(required=True)
    operations.add_argument('--write', nargs='?', const=str(OUTPUT))
    operations.add_argument('--check', nargs='?', const=str(OUTPUT))
    args = parser.parse_args()
    if args.write:
        Path(args.write).write_text(canonical(report()))
        print('wrote '+args.write)
    else:
        validate(strict_json.load(Path(args.check)))
        print('same-Higgs metric and generalized-volume replay passed')


if __name__ == '__main__':
    main()
