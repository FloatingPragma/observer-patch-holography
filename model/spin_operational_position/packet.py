"""Independent replay of source coordinate synthesis and complete native output."""
import ast
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path

from tools import strict_json
from . import bounds, circuits, experiment, geometry, generic

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'


def bindings():
    names = {p.relative_to(ROOT).as_posix() for p in HERE.iterdir()
             if p.suffix in ('.py', '.md') or p.name == 'reflectors.json'}
    names.update(('model/native_ordered_ground/witness.npz',
                  'model/native_spin_observer/README.md',
                  'model/native_ordered_ground/README.md',
                  'model/spin_covariant_localization/PROOF.md'))
    pending = list(names)
    visited = set()

    def add(module):
        if not module.startswith(('model.', 'tools.')):
            return
        parts = module.split('.')
        stem = '/'.join(parts)
        for name in [stem+'.py', stem+'/__init__.py']+[
                '/'.join(parts[:i])+'/__init__.py' for i in range(1, len(parts))]:
            if (ROOT/name).is_file() and name not in names:
                names.add(name)
                pending.append(name)

    while pending:
        name = pending.pop()
        if name in visited or not name.endswith('.py'):
            continue
        visited.add(name)
        module = name[:-3].replace('/', '.')
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
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(names)}


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def execution_certificate(run):
    probabilities = run['complete_record_probabilities']
    if (run['native_pulses'] != 49228 or run['first_failure_slots'] != 172368
            or len(probabilities) != 16 or min(probabilities) < 0
            or run['complete_output_distance_diagnostic'] >= 2e-8
            or run['all_first_failure_probability'] >= 1e-20
            or run['normalization_residual'] >= 2e-8
            or abs(probabilities[0]-.5) >= 1e-7
            or sum(probabilities[1:8]) >= 1e-20):
        raise ArithmeticError('complete native detector execution failed its admission')
    # Scientific receipt stores stable bounds; binary64 residuals are checked
    # on every replay instead of being promoted to exact arithmetic data.
    return {'decoded_core_executed': True, 'full_twenty_bit_word_executed': False,
            'work_bits': 12, 'full_quantum_output_rows': 4096,
            'retained_reference_columns': 2, 'native_pulses': run['native_pulses'],
            'word_sha256': run['word_sha256'], 'first_failure_slots': run['first_failure_slots'],
            'all_sixteen_binary_labels_retained': True,
            'complete_output_distance_diagnostic_below': '1/50000000',
            'first_failure_probability_diagnostic_below': '1/100000000000000000000',
            'normalization_residual_below': '1/50000000',
            'dark_probability_interval': ['0.4999999', '0.5000001'],
            'unused_binary_label_probability_below': '1/100000000000000000000',
            'discarded_failure_or_accepted_branch_normalization': False}


@lru_cache(maxsize=1)
def _build():
    before = bindings()
    run = experiment.run()
    normal = ((1, 2), (0, 0), (-3, 1), (2, -4))
    _, generic_record = generic.rank_one_phase(normal, F(3, 7), 48)
    result = {'schema': 'oph.spin_operational_position.v1', 'source_bindings': before,
              'exact_source_geometry': geometry.certify(),
              'literal_native_detector': circuits.plan()[2],
              'complete_native_execution': execution_certificate(run),
              'complex_rational_rank_one_example': generic_record,
              'generic_finite_instrument_examples': [bounds.generic_instrument_plan(q, m, F(1, 10**6))
                  for q, m in ((1, 3), (8, 9), (24, 17))],
              'scope': {'complete_finite_cube_coordinate_source_word_constructed': True,
                        'arbitrary_reference_probabilities_and_luders_backaction_bounded': True,
                        'generic_computable_finite_instrument_under_declared_services': True,
                        'refined_dyadic_interface_distinct_from_binary64_execution': True,
                        'continuum_detector_convergence_proved_here': False,
                        'second_moment_tail_control_proved_here': False,
                        'bare_repair_alone_measures_internal_spin': False,
                        'actual_protected_topological_label_attached': False,
                        'reflection_coordinate_is_physical_particle_position': False,
                        'particle_mass_spectrum_or_lorentz_limit_identified': False,
                        'full_encoded_twenty_bit_native_execution_claimed': False}}
    if bindings() != before:
        raise RuntimeError('source changed during operational detector replay')
    return result


def build():
    value = deepcopy(_build())
    if value['source_bindings'] != bindings():
        raise RuntimeError('source changed after the cached detector replay')
    return value


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path))) != canonical(build()):
        raise ValueError('source detector proof, program, output or scope differs')
    return True
