"""Reproduce the full-source color-loop coefficient certificate."""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path

from tools import strict_json
from model.interacting_extension._source import verify_sources
from model.interacting_extension.higgs_join.native_inputs import SOURCE_NAMES
from . import algebra, readout, source

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def bindings():
    names = {p.relative_to(ROOT).as_posix() for p in HERE.iterdir()
             if p.suffix in ('.py', '.md')}
    names.update(verify_sources(ROOT, SOURCE_NAMES))
    names.add('model/interacting_extension/source_manifest.json')
    names.add('model/requirements.txt')
    for name in ('joined_basis.json', 'joined_exact_operators.json',
                 'joined_operator_receipt.json', 'joined_operators.npz', 'PROOF.md'):
        names.add('model/interacting_extension/higgs_join/'+name)
    pending = list(names)
    visited = set()
    def add(module):
        if not module.startswith(('model.', 'tools.')):
            return
        parts = module.split('.')
        stem = '/'.join(parts)
        paths = [stem+'.py', stem+'/__init__.py']
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
    result = {
        'schema':'oph.color_loop_readout.v1',
        'source_bindings':before,
        'literal_source_replay':source.source_replay(),
        'full_space_moment_identities':algebra.certificate(),
        'finite_time_complete_experiments':[readout.witness(q) for q in (0, 1)],
        'premise_consumption':{
            'action':'the declared two-site three-family charged-lepton Higgs join',
            'family_controls':'identity hopping and Yukawa matrices; zero phases',
            'state':'empty Q0 input and its actual reference-neutral addition image in Q1',
            'clock':'the same finite Hamiltonian evolution time in every setting',
            'instrument':'complete gauge-invariant Higgs-occupation and two-link-block PVM',
            'statistical_resource_bound':'independent repetitions with the stated probability tolerance',
        },
        'scope':{
            'complete_438_and_1596_Gauss_spaces_retained':True,
            'all_Higgs_and_noncolored_record_outcomes_retained':True,
            'full_interacting_generator_used_at_every_time':True,
            'arbitrary_untouched_reference_in_instrument_error_bound':True,
            'exact_rational_probability_and_coefficient_enclosures':True,
            'default_bare_electric_over_abs_magnetic':'1',
            'default_color_Casimir_coefficient_over_abs_magnetic':'2',
            'physical_coupling_or_particle_data_used':False,
            'flux_instrument_compiled_from_bare_repair_records':False,
            'action_controls_and_clock_selected_from_bare_repair':False,
            'colored_chiral_matter_and_refining_3plus1_field_constructed':False,
            'renormalized_QCD_coupling_identified':False,
            'physical_alpha_s':None,
            'physical_particle_masses':None,
        },
    }
    if bindings() != before:
        raise RuntimeError('source bytes changed during certificate generation')
    return result


def verify(path):
    if canonical(strict_json.load(Path(path))) != canonical(build()):
        raise ValueError('color-loop source, complete probabilities, error bounds or scope differ')
