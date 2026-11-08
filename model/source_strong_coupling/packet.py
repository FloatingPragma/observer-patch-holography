"""Bind the exact native coefficient readout and its physical scope."""
from pathlib import Path
import ast
import hashlib
import importlib.util
import json

from tools import strict_json
from . import algebra, readout

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'


def bindings():
    names={p.relative_to(ROOT).as_posix() for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    names.update(('code/source_selection_model/response.json',
                  'model/native_quantum_gauss/PROOF.md',
                  'model/native_gauge_strings/PROOF.md',
                  'model/native_gauge_strings/strings.py',
                  'Lean/Screen/GaugeKineticInvariantForms.lean',
                  'code/particles/qcd/derive_lambda_qcd_source_transmutation.py',
                  'code/particles/qcd/test_lambda_qcd_source_transmutation.py',
                  'code/particles/calibration/legacy_d10/particle_masses_paper_d10_d11.py',
                  'code/particles/runs/qcd/hadron_source_backend/qcd_ensemble/source_action.json',
                  'code/particles/runs/qcd/hadron_source_backend/qcd_ensemble/source_parameter_map.json',
                  'model/requirements.txt'))
    # The native port stack imports these retained simulator modules through
    # its configured evidence root. Bind the whole two source packages,
    # including package initializers, rather than omitting path-based imports.
    retained=ROOT/'evidence/observer_dynamics_20260925/oph-physics-sim'
    for package in ('oph_exact','oph_fpe'):
        names.update(p.relative_to(ROOT).as_posix() for p in (retained/package).rglob('*.py'))
    pending=list(names); visited=set()
    def add(module):
        if not module.startswith(('model.','tools.')):return
        parts=module.split('.'); stem='/'.join(parts)
        paths=[stem+'.py',stem+'/__init__.py']+['/'.join(parts[:i])+'/__init__.py' for i in range(1,len(parts))]
        for path in paths:
            if (ROOT/path).is_file() and path not in names:
                names.add(path);pending.append(path)
    while pending:
        name=pending.pop()
        if name in visited or not name.endswith('.py'):continue
        visited.add(name);module=name[:-3].replace('/','.')
        package=module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node,ast.Import):
                for alias in node.names:add(alias.name)
            elif isinstance(node,ast.ImportFrom):
                target='.'*node.level+(node.module or '')
                target=importlib.util.resolve_name(target,package) if node.level else target
                add(target)
                for alias in node.names:add(target+'.'+alias.name)
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(names)}


def declared_inputs():
    """Report existing inputs separately; they do not enter the native witness."""
    folder=ROOT/'code/particles/runs/qcd/hadron_source_backend/qcd_ensemble'
    action=strict_json.load(folder/'source_action.json')
    mapping=strict_json.load(folder/'source_parameter_map.json')
    filename=ROOT/'code/particles/qcd/derive_lambda_qcd_source_transmutation.py'
    source=filename.read_text();values={}
    for node in ast.parse(source).body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ('ALPHA_3_MZ','MU_Z_OVER_E',
                         'E_STAR_DISPLAY_GEV','MB_MB_GEV','MC_MC_GEV'):
                    values[target.id]=ast.get_source_segment(source,node.value)
    if len(values)!=5:raise ArithmeticError('conditional running inputs changed')
    return {'source_action_status':action['status'],
            'source_parameter_map_status':mapping['status'],
            'source_action_prediction_promoted':action['promotion_allowed'],
            'source_parameter_prediction_promoted':mapping['promotion_allowed'],
            'conditional_running_inputs':values,
            'declared_inputs_enter_native_coefficient_witness':False,
            'alpha_s_prediction_derived_from_native_defaults':False}


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def build():
    before=bindings()
    result={'schema':'oph.source_strong_coupling.v1','source_bindings':before,
            'native_color_and_electric_action':algebra.source_certificate(),
            'complete_two_time_readout':readout.witness(),
            'inherited_input_boundary':declared_inputs(),
            'scope':{'fixed_native_J_and_epsilon_readable':True,
                     'actual_source_Gauss_and_electric_action_replayed':True,
                     'all_position_outcomes_and_reference_backaction_retained':True,
                     'exact_rational_probability_and_coefficient_enclosures':True,
                     'bare_repair_selects_this_auxiliary_link_design':False,
                     'source_defaults_identified_with_renormalized_QCD':False,
                     'QCD_continuum_and_scheme_matching_derived':False,
                     'empirical_strong_coupling_used_to_choose_a_coefficient':False,
                     'Hubbard_U_over_c_identified_with_alpha_s':False,
                     'physical_alpha_s_prediction':None}}
    if bindings()!=before:raise RuntimeError('source changed during verification')
    return result


def verify(path):
    if canonical(strict_json.load(Path(path)))!=canonical(build()):
        raise ValueError('source strong-coupling receipt, source or scope differs')
