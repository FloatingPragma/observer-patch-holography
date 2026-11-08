"""Source-bound replay of the complete spin response-memory certificate."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
from tools import strict_json
from . import algebra, asymptotic, charge_detector, charge_memory, finite, zero_field

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
PROOFS = ('model/spin_canonical_field/README.md',
          'model/spin_canonical_energy/README.md',
          'model/spin_current_spectrum/README.md',
          'model/spin_shell_diagnostic/README.md',
          'model/spin_shell_readout/PROOF.md',
          'model/spin_relaxed_response/README.md',
          'model/native_spin_observer/README.md',
          'model/spin_operational_position/PROOF.md',
          'model/native_ordered_ground/README.md',
          'model/protected_defect_catalog/PROOF.md')


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def canonical(value):
    return json.dumps(encode(value), sort_keys=True, indent=2, allow_nan=False)+'\n'


def bindings():
    legacy = 'evidence/observer_dynamics_20260925/oph-physics-sim/'
    names = set(PROOFS)|{finite.WITNESS, zero_field.WITNESS,
                        'code/source_selection_model/response.json',
                        'code/m1_source_realization/receipt.json',
                        legacy+'oph_exact/carrier.py', legacy+'oph_exact/__init__.py'}
    names.update(path.relative_to(ROOT).as_posix() for path in HERE.iterdir()
                 if path.suffix in ('.py', '.md'))
    pending = list(names)
    visited = set()

    def add(module):
        parts = module.split('.')
        prefix = legacy if parts[0] in ('oph_exact', 'oph_fpe') else ''
        stem = prefix+'/'.join(parts)
        candidates = [stem+'.py', stem+'/__init__.py']
        candidates += [prefix+'/'.join(parts[:i])+'/__init__.py' for i in range(1, len(parts))]
        for name in candidates:
            if (ROOT/name).is_file() and name not in names:
                names.add(name)
                pending.append(name)
    while pending:
        name = pending.pop()
        if name in visited or not name.endswith('.py'):
            continue
        visited.add(name)
        module = name.removeprefix(legacy)[:-3].replace('/', '.')
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                target = importlib.util.resolve_name(target, package) if node.level else target
                add(target)
                for alias in node.names:
                    add(target+'.'+alias.name)
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def build():
    before = bindings()
    finite.verified_source.cache_clear()
    charge_detector.native_two_site_execution.cache_clear()
    value = {'schema': 'oph.spin_dynamical_memory.v1', 'source_bindings': before,
        'exact_spectral_controls': algebra.exact_controls(),
        'complete_source_memory': finite.memory_certificate(),
        'complete_source_noncentral_history': finite.bounded_history_certificate(),
        'complete_same_box_zero_field_memory': zero_field.build(),
        'complete_source_charge_resolved_memory': charge_memory.build(),
        'complete_coherent_charge_class_readout': charge_detector.build(),
        'complete_memory_score_error_example': algebra.readout_memory_error(
            F(1, 4), F(28), F(1, 10**8), F(1, 10**8)),
        'actual_source_asymptotic_filter_bounds': asymptotic.source_filter_bounds(),
        'two_speed_local_field_implication_control': asymptotic.two_speed_control(),
        'scope': {
            'source_ground_witness_and_complete_Hamiltonian_replayed': True,
            'all_256_finite_spectral_levels_retained': True,
            'exact_response_complement_memory_and_causal_equation': True,
            'strictly_positive_finite_source_memory': True,
            'strictly_positive_finite_noncentral_history': True,
            'arbitrary_reference_complete_control_error_bound': True,
            'same_bounded_history_filtered_continuum_memory_identification': True,
            'declared_coherent_controller_interface_required': True,
            'native_controlled_history_tape_executed': False,
            'infinite_ground_preparation_effective_modulus': False,
            'continuum_memory_numerically_evaluated': False,
            'noncentral_history_survives_field_scaling': False,
            'sharp_mass_shell_Lorentz_covariance_or_scattering_proved': False,
            'finite_gap_identified_as_physical_particle_mass': False}}
    if bindings() != before:
        raise RuntimeError('dynamical source changed during replay')
    return encode(value)


def validate(value):
    if canonical(value) != canonical(build()):
        raise ValueError('dynamical memory, source binding, complete outcomes or scope differ')
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    group.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    args = parser.parse_args()
    if args.write:
        args.write.write_text(canonical(build()))
    else:
        validate(strict_json.load(args.check))
    print('PASS: complete spin dynamical memory and noncentral source history')
