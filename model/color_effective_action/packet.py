"""Reproduce the finite action and its full-source custody."""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path

from tools import strict_json
from model.color_loop_readout.packet import bindings as parent_bindings
from model.color_loop_readout.source import Controls, source_replay
from . import action, link, rotor, nonabelian

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def bindings():
    names = set(parent_bindings())
    names.update(p.relative_to(ROOT).as_posix() for p in HERE.iterdir()
                 if p.suffix in ('.py', '.md'))
    names.update((
        'model/native_gauge_strings/strings.py',
        'model/common_source/native_gauge_continuum.py',
        'model/particle_gauge_band/QUANTUM_LINK_PROOF.md',
        'model/particle_gauge_band/quantum_link_completion.py',
        'model/interacting_extension/chiral_lattice/PROOF.md',
    ))
    pending = list(names)
    visited = set()
    def add(module):
        if not module.startswith(('model.', 'tools.')):
            return
        parts = module.split('.')
        paths = ['/'.join(parts)+'.py', '/'.join(parts)+'/__init__.py']
        paths += ['/'.join(parts[:i])+'/__init__.py' for i in range(1, len(parts))]
        for path in paths:
            if (ROOT/path).is_file() and path not in names:
                names.add(path)
                pending.append(path)
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
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(names)}


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def build():
    before = bindings()
    controls = (Controls(), Controls(action.F(3, 2), action.F(2, 3), action.F(5, 4),
                                    action.F(7, 5), action.F(3, 2)))
    result = {
        'schema':'oph.color_effective_action.v1',
        'source_bindings':before,
        'literal_source_replay':source_replay(),
        'native_link_coordinate_test':link.certificate(),
        'arbitrary_control_component_identities':action.component_certificate(),
        'full_source_moment_certificates':[action.exact_moment_certificate(c) for c in controls],
        'effective_actions':[action.witness(action.F(mu)) for mu in (64, 128, 256)],
        'finite_action_family_criteria':action.uniform_family_remainders(
            action.F(64), action.Elimination().d_bound, action.F(1)),
        'distinct_rotor_regulator':{
            'literal_source_columns':[rotor.source_replay(q, cutoff)
                                      for q in (-4, -3, 1, 2, 6) for cutoff in (4, 12, 48)],
            'cutoff_uniform_comparisons':[{'charge6':q, 'cutoff':cutoff,
                                         **rotor.cutoff_error(q, cutoff)}
                                        for q in (-4, -3, 1, 2, 6) for cutoff in (4, 12, 48)],
            'complete_infinite_flux_actions':[rotor.action_witness(q) for q in (-4, -3, 1, 2, 6)],
            'same_as_native60_Higgs_source':False,
            'nonabelian_QCD_identification':False,
        },
        'conditional_nonabelian_coordinate_regulator':{
            'source_normalization_and_growth':nonabelian.growth_certificate(),
            'same_input_source_identity_test':nonabelian.completion_response_witness(),
            'fixed_graph_bounds':[nonabelian.cutoff_bounds(action.F(radius)) for radius in (8, 16, 32, 64)],
            'same_source_capacity_attachment_derived':False,
            'physical_alpha_s':None,
        },
        'premises':{
            'action':'the fixed two-vertex three-family Higgs-join source with identity family controls',
            'space':'all 438 Gauss states in Qstar0, split into 6 scalar-vacuum and 432 occupied-scalar states',
            'state':'the original empty source input gives the finite reflection-positive boundary state',
            'clock':'the same specified H; mu has the same model-energy units',
            'readout':'literal simulator specification and the declared source operators',
        },
        'scope':{
            'positive_source_Euclidean_transfer_derived':True,
            'reflection_positive_boundary_state_equals_stationary_ground':False,
            'full_complement_integrated_without_postselection':True,
            'exact_six_dimensional_Schur_action_derived':True,
            'operator_norm_and_all_entry_enclosures':True,
            'complete_filter_operator_instrument_proved':True,
            'filter_instrument_native_pulse_execution_claimed':False,
            'source_generated_additional_gauge_invariant_operators':True,
            'electric_magnetic_projection_residual_discarded':False,
            'direct_Wilson_coordinate_substitution_valid':False,
            'collective_Wilson_reconstruction_excluded':False,
            'spatial_refining_source_family_derived':False,
            'physical_scale_or_QCD_scheme_matching_derived':False,
            'physical_alpha_s':None,
            'physical_particle_masses':None,
        },
    }
    if bindings() != before:
        raise RuntimeError('source bytes changed during calculation')
    return result


def verify(path):
    if canonical(strict_json.load(Path(path))) != canonical(build()):
        raise ValueError('finite color action, source, interval or scope differs')
