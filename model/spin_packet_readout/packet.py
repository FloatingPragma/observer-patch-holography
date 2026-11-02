"""Current narrow proof custody and replay of the new arithmetic/finite controls."""
import ast
import hashlib
import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path
from . import bounds, coefficients, controls, envelope, filters, finite

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
FRAGMENT = 'paper/tex_fragments/SPIN_PACKET_READOUT.tex'
PROOFS = ('model/spin_goldstone_sector/README.md',
          'model/native_ordered_ground/README.md',
          'model/native_spin_observer/README.md',
          'model/native_spin_observer/precision.py',
          'model/native_joint_spin_spectroscopy/README.md',
          'model/spin_canonical_energy/README.md')


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def bindings():
    paths = {p.relative_to(ROOT).as_posix() for p in HERE.glob('*.py')}
    paths |= {str((HERE/'README.md').relative_to(ROOT)), FRAGMENT, finite.WITNESS}
    pending = list(paths)
    def add(module):
        if not module.startswith(('model.', 'tools.')):
            return
        stem = module.replace('.', '/')
        pieces = stem.split('/')
        for candidate in [stem+'.py', stem+'/__init__.py']+['/'.join(pieces[:i])+'/__init__.py' for i in range(1, len(pieces))]:
            if (ROOT/candidate).is_file() and candidate not in paths:
                paths.add(candidate)
                pending.append(candidate)
    while pending:
        name = pending.pop()
        if not name.endswith('.py'):
            continue
        module = name[:-3].replace('/', '.')
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node, ast.Import):
                for entry in node.names:
                    add(entry.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                target = importlib.util.resolve_name(target, package) if node.level else target
                add(target)
                for entry in node.names:
                    add(target+'.'+entry.name)
    paths.update(PROOFS)
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def build():
    before = bindings()
    table = finite.finite_readout()
    quadrature = bounds.quadrature_plan(1, 2, F(1, 4), 1, F(1, 10000), F(1, 4))
    labels = 10*(quadrature['cells']+2)
    native = bounds.native_plan(8, 12, 24, 8, labels, F(quadrature['time_cutoff']),
                                1, F(1, 4), F(1, 100000), F(1, 1000000))
    eta = quadrature['total_time_filter_error_upper']/5+native['maximum_conjugation_product_error']+native['native_complete_isometry_error_upper']
    q = min(F(row['ideal_LCU_success_probability_interval'][0]) for row in table['rows'])
    conditional = bounds.conditioned_errors(q, eta, F(1, 10**8))
    if conditional['conditional_trace_error_upper'] >= F(1, 100):
        raise ArithmeticError('finite demonstration failed its conditional error target')
    result = {
        'schema': 'oph.spin_packet_readout.v1',
        'source_bindings': before,
        'new_checks_pass': True,
        'explicit_Goldstone_filter': filters.encode(filters.goldstone_filter()),
        'explicit_band_limited_envelope': envelope.certificate(),
        'explicit_envelope_success_examples': filters.encode([
            bounds.packet_success(p, F(envelope.certificate()['sampling_LCU_bound_Lambda']))
            for p in (F(1), F(1, 2), F(1, 4))]),
        'finite_source_readout': table,
        'finite_native_resource_plan': filters.encode({'quadrature': quadrature,
            'native': native, 'success_amplitude_error_upper': eta,
            'input_state_trace_error_allowance': F(1, 10**8),
            'conditional': conditional}),
        'computed_complex_coefficient': filters.encode(coefficients.midpoint_coefficient(
            F(1, 3), F(1, 16), 1, 2, F(1, 4), F(1, 10**18))),
        'computed_spatial_coefficients': filters.encode(coefficients.contractive_coefficients(
            [envelope.coefficient_rectangle(site, F(1, 2))
             for site in ((0, 0, 0), (1, 0, 0), (-1, 0, 0))], F(1, 10**16))),
        'complete_instrument_controls': controls.exact_control(),
        'theorem': {
            'source': 'unchanged cubic antiferromagnetic spin source; selected finite H_B variant explicitly separate',
            'explicit_filter_equal_one_on_certified_band': True,
            'all_LCU_rejections_retained': True,
            'complete_positive_energy_POVM_composition_available': True,
            'coarse_energy_and_reject_l1_error': '4*success_operator_error+input_trace_error+detector_error',
            'absolute_ground_energy_not_needed_for_Heisenberg_filter': True,
            'success_lower_scales_as_p_to_four': True,
            'native_coherent_control_overhead_included': True,
            'fixed_binary64_extrapolated_to_arbitrary_precision': False,
            'finite_ground_input_source_computable': True,
            'effective_box_size_for_infinite_ordered_ground_proved': False,
            'variational_packet_identified_with_continuum_Stone_field': False,
            'new_native_execution_performed': False,
            'sharp_mass_shell_or_physical_species_identified': False,
        },
        'custody_scope': 'Current runtime import closure, exact retained witness and named mathematical proof inputs. No full inherited receipt replay or repository-wide CI claim.',
    }
    if before != bindings():
        raise RuntimeError('source changed during packet construction')
    return result


def validate(value, expected=None):
    if canonical(value) != canonical(build() if expected is None else expected):
        raise ValueError('readout arithmetic, source bindings, or scope changed')
