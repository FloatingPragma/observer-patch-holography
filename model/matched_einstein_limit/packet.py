"""Replay exact algebra and the source-bound conditional Lean theorem."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

from tools import strict_json
from . import frame

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
LEAN = 'Lean/ObserverPatchHolography/EinsteinBranch/QuantitativeSmallBall.lean'
THEOREMS = ('tensor_bound', 'residual_bound', 'small_ball_tensor_bound',
            'normalized_error_tendsto', 'tensor_tendsto_zero', 'lambda_constant')


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def source_bindings():
    paths = {str(path.relative_to(ROOT)) for path in Path(__file__).parent.iterdir()
             if path.suffix in ('.py', '.md')}
    paths.update(('model/__init__.py', 'tools/strict_json.py', 'model/requirements.txt',
                  'Lean/lakefile.lean', 'Lean/lake-manifest.json', 'Lean/lean-toolchain', LEAN))
    pending = [LEAN]
    while pending:
        path = pending.pop()
        for module in re.findall(r'^import\s+(\S+)', (ROOT/path).read_text(), re.MULTILINE):
            local = 'Lean/'+module.replace('.', '/')+'.lean'
            if (ROOT/local).is_file() and local not in paths:
                paths.add(local)
                pending.append(local)
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def scope():
    return {
        'status': 'CONDITIONAL_MATCHED_REGIONAL_LIMIT',
        'same_source_state_and_refinement_required': True,
        'actual_metric_curvature_and_action_stress_identification_required': True,
        'common_proper_clock_and_fixed_volume_regions_required': True,
        'regional_modular_and_entropy_area_estimates_required': True,
        'all_ten_probes_share_one_isotropic_vacuum_reference': True,
        'positive_spacetime_constant_area_coefficient_limit_required': True,
        'actual_finite_deficits_required_for_nonlinear_conclusion': True,
        'single_vacuum_tangents_imply_only_linearized_equation': True,
        'independent_Ward_and_Bianchi_required_for_constant_Lambda': True,
        'separate_null_tensor_balance_hypothesis_required': False,
        'regional_source_premises_instantiated': False,
        'reciprocal_thermal_coefficient_is_Newton_coefficient': False,
        'finite_profiles_prove_refinement': False,
        'Einstein_dynamics_derived_from_executable_source': False,
        'physical_gate_accepted': False,
    }


def parse_lean_output(stdout, stderr=''):
    if 'sorryAx' in stdout+stderr or 'warning:' in stdout+stderr:
        raise ValueError('Lean proof contains an admission or warning')
    approved = {'propext', 'Classical.choice', 'Quot.sound'}
    actual = {}
    for name in THEOREMS:
        pattern = (re.escape("'OPH.EinsteinBranch.Quantitative."+name)
                   + r"' depends on axioms:\s*\[([^\]]*)\]")
        matches = re.findall(pattern, stdout)
        if len(matches) != 1:
            raise ValueError('missing Lean theorem audit: '+name)
        axioms = [entry.strip() for entry in matches[0].split(',') if entry.strip()]
        if len(axioms) != len(set(axioms)) or not set(axioms) <= approved:
            raise ValueError('unapproved Lean theorem axiom: '+name)
        actual[name] = sorted(axioms)
    return actual


def lean_replay():
    pins = source_bindings()
    subprocess.run(['lake', 'build', 'ObserverPatchHolography.EinsteinBranch.QuantitativeSmallBall'],
                   cwd=ROOT/'Lean', check=True, capture_output=True, text=True, timeout=600)
    completed = subprocess.run(
        ['lake', 'env', 'lean', LEAN.removeprefix('Lean/')], cwd=ROOT/'Lean',
        check=True, capture_output=True, text=True, timeout=240)
    actual = parse_lean_output(completed.stdout, completed.stderr)
    if source_bindings() != pins:
        raise ValueError('source changed during Lean build or kernel audit')
    return {'module': LEAN, 'fresh_source_execution': True,
            'audited_theorems': list(THEOREMS),
            'actual_axioms_by_theorem': actual,
            'constructs_source_premise_instance': False}


def lean_source_certificate():
    source = (ROOT/LEAN).read_text()
    for name in THEOREMS:
        if 'theorem '+name+' ' not in source or 'assert_no_sorry '+name not in source:
            raise ValueError('missing Lean theorem or admission audit: '+name)
    if re.search(r'^\s*(axiom|sorry|admit)\b|\bnative_decide\b', source, re.MULTILINE):
        raise ValueError('unsupported Lean admission or native proof shortcut')
    return {'module': LEAN, 'audited_theorems': list(THEOREMS),
            'source_bound': True,
            'Lean_kernel_reexecuted_by_Python_receipt': False,
            'kernel_check_command': 'python -m model.matched_einstein_limit --lean-check',
            'constructs_source_premise_instance': False}


def report():
    pins = source_bindings()
    value = {'schema': 'oph.matched_einstein_limit.v1', 'all_pass': True,
             'issues': [1145, 1146], 'gates': ['CP-09', 'CP-10'],
             'frame_and_bounds': frame.certificate(), 'lean': lean_source_certificate(),
             'scope': scope(), 'source_bindings': pins}
    if source_bindings() != pins:
        raise ValueError('source changed during conditional Einstein replay')
    return value


def validate(value, expected):
    if value.get('scope') != scope():
        raise ValueError('conditional regional interface or scope changed')
    if canonical(value) != canonical(expected):
        raise ValueError('conditional Einstein receipt differs from fresh replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--lean-check', action='store_true')
    args = parser.parse_args()
    if args.lean_check:
        print(canonical(lean_replay()), end='')
        return
    value = report()
    if args.write:
        args.write.write_text(canonical(value))
    else:
        validate(strict_json.load(args.check), value)
    print('matched regional Einstein implication: PASS; source premises uninstantiated')


if __name__ == '__main__':
    main()
