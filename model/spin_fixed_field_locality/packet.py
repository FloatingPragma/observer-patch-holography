"""Source-bound fixed positive-field marginal certificates."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

from tools import strict_json
from . import bounds,finite

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).parent


def serialize(value):
    if type(value)is F:return str(value)
    if type(value)is float:return format(value,'.10f')
    if isinstance(value,dict):return {key:serialize(v)for key,v in value.items()}
    if isinstance(value,(list,tuple)):return [serialize(v)for v in value]
    return value


def bindings():
    paths=set(HERE.glob('*.py'));queue=list(paths)
    def include(path):
        if path.is_file()and path not in paths:paths.add(path);queue.append(path)
    while queue:
        path=queue.pop();package='.'.join(path.relative_to(ROOT).with_suffix('').parts[:-1])
        for node in ast.walk(ast.parse(path.read_text())):
            names=[]
            if isinstance(node,ast.Import):names=[a.name for a in node.names]
            elif isinstance(node,ast.ImportFrom):
                name='.'*node.level+(node.module or '')
                base=importlib.util.resolve_name(name,package)if node.level else name
                names=[base]+[base+'.'+a.name for a in node.names]
            for name in names:
                parts=name.split('.');include(ROOT.joinpath(*parts).with_suffix('.py'))
                for j in range(1,len(parts)+1):include(ROOT.joinpath(*parts[:j])/'__init__.py')
    paths.update((HERE/'README.md',HERE/'PROOF.md'))
    paths.update(ROOT/p for p in ('model/native_ordered_ground/README.md',
        'model/spin_canonical_field/README.md','model/native_spin_clock_join/README.md',
        'model/protected_spin_limit/README.md'))
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted(paths)}


def claims():
    return {'actual_source_Lee_Yang_hypotheses_instantiated':True,
        'fixed_positive_field_uniform_boundary_path_gap':True,
        'Gaussian_flow_bound_derived_for_actual_ground':True,
        'exterior_volume_independent_marginal_modulus':True,
        'rational_source_box_search_terminates':True,
        'complete_local_instrument_and_independent_reference_bound':True,
        'imported_Lee_Yang_theorem_reproved_or_Lean_checked_here':False,
        'arbitrary_source_reference_correlations_covered':False,
        'arbitrary_boundary_fields_covered':False,
        'zero_field_ordered_state_selection_modulus_supplied':False,
        'positive_field_gap_identified_as_particle_mass':False,
        'locality_bound_identified_as_physical_light_speed':False,
        'large_planned_source_box_actually_simulated':False}


def build():
    before=bindings()
    value={'schema':'oph.spin_fixed_field_locality.v1','all_pass':True,
        'source_bindings':before,
        'external_theorem':{'url':'https://arxiv.org/html/2607.10765v1',
            'authors':['Chaithanya Rayudu','Jun Takahashi'],
            'title':'Spectral gap of Lee-Yang Hamiltonians','version':'2607.10765v1',
            'consumed':'Definition1 and Theorem11; real nonnegative longitudinal coefficients and strictly positive uniform field'},
        'source_attachment':finite.exact_source_attachment(),
        'finite_source_derivative_controls':serialize(finite.certificate()),
        'box_certificates':serialize([bounds.choose_box(F(1),b,1,27,eps)
            for b in (F(1),F(1,4),F(1,16))for eps in (F(1,100),F(1,10000))]),
        'claims':claims()}
    if bindings()!=before:raise ValueError('source changed during fixed-field marginal checks')
    return value


def validate(value):
    expected=build()
    if json.dumps(value,sort_keys=True,allow_nan=False)!=json.dumps(expected,sort_keys=True,allow_nan=False):
        raise ValueError('source binding, bound, source control or fixed-field scope changed')
    return expected


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',nargs='?',type=Path,const=HERE/'results.json')
    group.add_argument('--check',nargs='?',type=Path,const=HERE/'results.json')
    args=parser.parse_args()
    if args.write:args.write.write_text(json.dumps(build(),indent=2,sort_keys=True,allow_nan=False)+'\n')
    else:validate(strict_json.load(args.check))
    print('PASS: fixed positive-field source marginal modulus')


if __name__=='__main__':main()
