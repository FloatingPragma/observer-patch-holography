"""Replay exact spin-response identities, bounds and complete source intervals."""
import ast
import hashlib
import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path
from tools import strict_json
from . import algebra, bounds, finite

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
PROOFS = ('model/spin_canonical_field/README.md',
          'model/spin_canonical_energy/README.md',
          'model/spin_current_spectrum/README.md',
          'model/spin_packet_readout/README.md',
          'model/native_ordered_ground/README.md')


def canonical(value):
    return json.dumps(bounds.encode(value), sort_keys=True, indent=2, allow_nan=False)+'\n'


def bindings():
    names = set(PROOFS)
    names.add(finite.WITNESS)
    names.update(p.relative_to(ROOT).as_posix() for p in HERE.iterdir() if p.suffix in ('.py', '.md'))
    pending = list(names); visited = set()
    def add(module):
        if not module.startswith(('model.', 'tools.')):
            return
        parts = module.split('.'); stem = '/'.join(parts)
        for name in [stem+'.py', stem+'/__init__.py']+['/'.join(parts[:i])+'/__init__.py' for i in range(1, len(parts))]:
            if (ROOT/name).is_file() and name not in names:
                names.add(name); pending.append(name)
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
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def build():
    before = bindings()
    value = {'schema': 'oph.spin_relaxed_response.v1', 'source_bindings': before,
             'exact_variational_controls': algebra.exact_checks(),
             'actual_continuum_response_bounds': [bounds.continuum_response(c, F(1, 1000), r)
                 for c in ('P', 'Q') for r in (F(1), F(1, 2), F(1, 4))],
             'complete_source_resolvent_response': finite.finite_response(),
             'source_only_filter_example': bounds.filter_budget(F(1, 4), 1024, F(1, 4096), 28, F(1, 10**6)),
             'complete_score_error_example': bounds.finite_readout_error(
                 F(1, 4), F(1, 10**6), F(1, 10**6), far_label_probability=F(1, 1000)),
             'bounded_source_probe_error_example': bounds.bounded_probe_error(
                 F(1, 4), F(1, 8), F(1, 10**8)),
             'scope': {'same_ordered_spin_continuum_Q_and_P_vectors': True,
                       'relaxed_tangent_energy_exact_variational_formula': True,
                       'nonzero_finite_relaxed_response_both_channels': True,
                       'ground_specific_inverse_response_coefficient_numerically_derived_in_continuum': False,
                       'bare_twist_coefficient_identified_with_relaxed_stiffness': False,
                       'first_energy_uniform_integrability_assumed': False,
                       'regularized_K_plus_eta_form_optimum_attained': True,
                       'actual_Stone_susceptibility_from_bounded_source_ordered_limits': True,
                       'admitted_complete_controlled_probe_encoding_required': True,
                       'raw_unbounded_source_vector_limit_identified_without_projection_control': False,
                       'unregularized_optimum_attained_without_domain_hypothesis': False,
                       'sharp_shell_Lorentz_covariance_or_scattering_derived': False,
                       'physical_particle_masses_identified': False,
                       'effective_infinite_ground_preparation_modulus_derived': False,
                       'all_inherited_producers_replayed': False}}
    if bindings() != before:
        raise RuntimeError('response source changed during replay')
    return bounds.encode(value)


def validate(value):
    if canonical(value) != canonical(build()):
        raise ValueError('response arithmetic, source binding or scientific scope differs')
    return True


def verify(path=OUTPUT):
    return validate(strict_json.load(Path(path)))
