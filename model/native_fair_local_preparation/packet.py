"""Replay exact local preparation and bind its native sources."""
import argparse
import ast
import importlib.util
import hashlib
import json
from pathlib import Path
from tools import strict_json
from model.native_sector_readout import packet as inherited
from model.native_effect_synthesis import packet as coherent_sources
from . import preparation as p, coupling, witness

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
SOURCE_FILES=tuple('model/native_fair_local_preparation/'+name for name in (
    '__init__.py','preparation.py','coupling.py','witness.py','packet.py','test_preparation.py','test_packet.py','test_coupling.py','PROOF.md','COUPLING.md','README.md'))+(
    'model/native_phase_experiment/README.md',
    'Lean/ObserverPatchHolography/Locality/NoSignalling.lean',
    'tools/strict_json.py','model/requirements.txt')


def source_bindings():
    paths=set(SOURCE_FILES)|set(inherited.source_bindings())|set(coherent_sources.source_bindings())
    pending=list(paths);visited=set()
    def add(module):
        if not module.startswith(('model.','tools.')):return
        parts=module.split('.');stem='/'.join(parts)
        for name in [stem+'.py',stem+'/__init__.py']+['/'.join(parts[:i])+'/__init__.py' for i in range(1,len(parts))]:
            if (ROOT/name).is_file() and name not in paths:paths.add(name);pending.append(name)
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
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def endpoint_control():
    checked=0
    for edge in p.SEAMS:
        a,b=edge
        for bits in range(4):
            mask=((bits&1)<<a)|((bits>>1)<<b)
            observed=sorted(p.readout.native.image(mask,edge,c,12)[0] for c in (0,1))
            swapped=mask^((1<<a)|(1<<b)) if bits in (1,2) else mask
            if observed!=sorted((mask,swapped)):
                raise ArithmeticError('fair native endpoint law differs from the certified population map')
            checked+=1
    return {'literal_endpoint_inputs':checked,'both_coins_per_input':2,'all30_seams':True}


def attachments():
    rows=[]
    for level in (0,1,2,3,200):
        value=p.program(level=level,total=2,blocks=32)
        rows.append({'level':level,'carrier':value['carrier'],'program_digest':p.source.digest(value),
                     'resources':value['declared_resources'],'literal_internal_word_admitted':True})
    return rows


def produce():
    bindings=source_bindings()
    examples=[p.source_query(total=q,mask=(1<<q)-1) for q in range(13)]
    examples.append(p.source_query(total=2,blocks=1,mask=3))
    for level in (1,2):examples.append(p.source_query(level=level,total=2,cell=1,mask=5))
    stop=p.source.execute(p.program(total=0,blocks=0),{(0,0,0,0):p.F(1)},max_branches=25)
    high=p.source.execute(p.program(level=3,total=0,blocks=32),{(0,0,0,0):p.F(1)})
    if stop.status!='RESOURCE_STOP' or stop.attempted_steps!=1 or high.status!='RESOURCE_STOP_UNEVALUATED_CARRIER' or high.attempted_steps!=0:
        raise ArithmeticError('the original expanded host boundary changed')
    result={'schema':'oph.native_fair_local_preparation.v1','endpoint_attachment':endpoint_control(),
            'exact_all_sector_minorization':p.minorization(),
            'occupation_wire_intertwining':coupling.intertwining_certificate(),
            'declared_coherent_services':coupling.service_certificate(),
            'conditional_complete_phase_witness':witness.certificate(),
            'chosen_budget':p.budget(32),'maximal_budget':p.budget(68),
            'source_admission':attachments(),'authenticated_raw_words':examples,
            'boundaries':{'expanded_control_branch_budget':25,'native_branch_ceiling':p.source.MAX_EVALUATED_BRANCHES,'expanded_tree_status':stop.status,'expanded_tree_attempts':stop.attempted_steps,
                          'higher_level_status':high.status,'higher_level_attempts':high.attempted_steps,
                          'complete_archive_independence_distance_at_q6':str(p.archive_separation(6)),
                          'constant_corner_dimension':2,'phase_only_record_control':p.phase_blindness_control()},
            'scope':{'all13_local_number_sectors':True,'arbitrary_untouched_reference':True,
                     'preparation_tree_archive_retained_in_error_comparison':True,
                     'mixing_archive_traced_only_for_comparison':True,
                     'complete_source_output_retains_all_raw_records':True,
                     'conditional_sector_weights_selected_from_input':True,
                     'all_global_canonical_states_prepared':False,
                     'independent_source_full_archive_reference_claimed':False,
                     'proper_code_occupation_wire_coupling_explicitly_supplied':True,
                     'growing_coherent_load_derived_from_passive_repair':False,
                     'physical_vacuum_derived':False,
                     'physical_clock_QFT_or_particle_masses_derived':False,
                     'native_program_grammar_or_host_limits_changed':False},
            'source_bindings':bindings}
    if source_bindings()!=bindings:raise ArithmeticError('preparation source bytes changed during replay')
    return result


def canonical(value):return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path)))!=canonical(produce()):
        raise ValueError('local preparation receipt differs from the exact producer')
    return True


def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.write:OUTPUT.write_text(canonical(produce()))
    else:verify()
    print('Fair local preparation receipt verified.')


if __name__=='__main__':main()
