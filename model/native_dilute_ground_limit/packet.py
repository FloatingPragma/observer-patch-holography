"""Replay source-ground dilution, complete endpoint laws and density closure."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

from tools import strict_json
from model.native_fair_local_preparation import packet as inherited
from . import dilution as d

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'


def encode(value):
    if type(value)is F:return str(value)
    if type(value)in (tuple,list):return [encode(x)for x in value]
    if type(value)is dict and all(type(k)is str for k in value):return {k:encode(v)for k,v in value.items()}
    if value is None or type(value)in (str,int,bool):return value
    raise TypeError('exact typed ground-dilution receipt required')


def canonical(value):return json.dumps(encode(value),sort_keys=True,indent=2,allow_nan=False)+'\n'


def source_bindings():
    paths=set(inherited.source_bindings())
    paths.update(p.relative_to(ROOT).as_posix()for p in HERE.iterdir()if p.suffix in ('.py','.md'))
    paths.update(('model/federated_observer/SCALING.md','model/native_spatial_scaling/PROOF.md',
                  'model/native_exchange_clock/PROOF.md'))
    pending=list(paths);visited=set()
    def add(module):
        if not module.startswith(('model.','tools.')):return
        parts=module.split('.');stem='/'.join(parts)
        for name in [stem+'.py',stem+'/__init__.py']+['/'.join(parts[:i])+'/__init__.py'for i in range(1,len(parts))]:
            if (ROOT/name).is_file()and name not in paths:paths.add(name);pending.append(name)
    while pending:
        name=pending.pop()
        if name in visited or not name.endswith('.py'):continue
        visited.add(name);module=name[:-3].replace('/','.')
        package=module.removesuffix('.__init__')if module.endswith('.__init__')else module.rpartition('.')[0]
        add(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node,ast.Import):
                for alias in node.names:add(alias.name)
            elif isinstance(node,ast.ImportFrom):
                target='.'*node.level+(node.module or '')
                target=importlib.util.resolve_name(target,package)if node.level else target
                add(target)
                for alias in node.names:add(target+'.'+alias.name)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()for name in sorted(paths)}


def produce():
    before=source_bindings()
    examples=[]
    for ports in (240,960,240<<400):
        examples.append({'ports':ports,'total':2,
            'local_observable_bound':d.observable_bounds(ports,2,1,2),
            'complete_history_bound':d.history_bound(ports,2,(2,2,4))})
    value={'schema':'oph.native_dilute_ground_limit.v1','source_bindings':before,
           'source_ground_dilution':examples,
           'literal_native_endpoint':d.native_endpoint_certificate(),
           'exact_density_ground_closure':d.density_closure_certificate(),
           'necessary_visibility':d.visibility_requirement(240,2,F(1),F(1,10)),
           'finite_density_control':{'ports':240,'total':120,'local_number_variance':'1/4'},
           'growing_norm_control':{'ports':240,'total':2,'observable':'sqrt(V) n0',
                                  'variance':F(2*238,240),'uniform_norm_bound':False},
           'scope':{'pure_source_Dicke_ground_required':True,
                    'each_fine_local_Kraus_preserves_occupation':True,
                    'complete_local_instrument_includes_vacuum_and_all_failure_slots':True,
                    'fixed_finite_complete_source_and_record_protocol':True,
                    'no_history_postselection':True,'uniform_over_controlled_clock_rescaling':True,
                    'sqrt_upper_grid_bits':128,'analytic_asymptotic_bound_has_rounding_floor':False,
                    'native_passive_repair_identified_with_unitary_flow':False,
                    'bounded_fixed_local_dilute_field_has_nonzero_spectral_residue':False,
                    'finite_density_growing_support_norm_or_charged_fields_excluded':False,
                    'spatial_density_scaling_dictionary_changed':False,
                    'relativistic_interacting_theory_ruled_out_in_all_dictionaries':False}}
    if source_bindings()!=before:raise ArithmeticError('source bytes changed during dilution replay')
    return encode(value)


def check_payload(retained,fresh):
    if canonical(retained)!=canonical(fresh):raise ArithmeticError('ground-dilution receipt differs from fresh owner')
    return True


def verify(path=OUTPUT):return check_payload(strict_json.load(Path(path)),produce())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.write:
        value=produce();OUTPUT.write_text(canonical(value),encoding='utf-8')
        check_payload(strict_json.load(OUTPUT),value)
        print('Native dilute-ground receipt written.')
    else:verify();print('Native dilute-ground receipt verified.')


if __name__=='__main__':main()
