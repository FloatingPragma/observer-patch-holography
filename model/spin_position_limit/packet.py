"""Verify finite source controls and exact Gram probability certificates."""
from pathlib import Path
from fractions import Fraction as F
import ast
import hashlib
import importlib.util
import json
from tools import strict_json
from . import bounds, finite

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'


def bindings():
    names={p.relative_to(ROOT).as_posix() for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    names.update(('model/spin_covariant_localization/CONTINUUM_POSITION.md',
                  'model/spin_canonical_field/README.md',
                  'model/native_spin_observer/README.md',
                  'model/spin_operational_position/PROOF.md',
                  'model/spin_operational_position/generic.py',
                  'model/spin_operational_position/bounds.py',
                  'model/spin_operational_position/results.json',
                  'model/native_ordered_ground/witness.npz',
                  'model/requirements.txt'))
    pending=list(names);visited=set()
    def add(module):
        if not module.startswith(('model.','tools.')):return
        parts=module.split('.');stem='/'.join(parts)
        candidates=[stem+'.py',stem+'/__init__.py']+['/'.join(parts[:i])+'/__init__.py' for i in range(1,len(parts))]
        for path in candidates:
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


def encode(value):
    if type(value)is F:return str(value)
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [encode(v) for v in value]
    if value is None or type(value)in(str,int,bool):return value
    raise TypeError('exact typed receipt values required')


def canonical(value):
    return json.dumps(encode(value),sort_keys=True,indent=2,allow_nan=False)+'\n'


def build():
    before=bindings()
    result={'schema':'oph.spin_position_limit.v1','source_bindings':before,
            'finite_source_controls':finite.controls(),
            'bounded_probe_error_example':bounds.bounded_probe_errors(3,F(1,10),F(1,1000000),F(1,1000000)),
            'exact_error_example':bounds.probability_error(F(1,100),F(1,1000000),F(1),F(1,1000000),F(1,100000)),
            'scope':{'finite_polar_sections_converge_to_same_field_position':True,
                     'same_bounded_history_source_limit_used':True,
                     'full_position_partition_and_dark_outcome_retained':True,
                     'fixed_history_statistical_convergence_proved_in_text':True,
                     'finite_coherent_observer_service_required':True,
                     'bare_repair_implements_spin_position':False,
                     'effective_ground_or_scaling_modulus_supplied':False,
                     'infinite_dimensional_uniform_diamond_convergence':False,
                     'second_moment_tail_control_supplied':False,
                     'finite_checks_prove_infinite_theorem':False,
                     'mass_shell_or_empirical_species_identification':False}}
    if bindings()!=before:raise RuntimeError('source changed during Gram verification')
    return result


def verify(path):
    if canonical(strict_json.load(Path(path)))!=canonical(build()):
        raise ValueError('position limit receipt, source or scope differs')
