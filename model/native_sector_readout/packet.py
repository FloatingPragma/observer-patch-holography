"""Reproduce the full-number native configuration-readout certificate."""
from __future__ import annotations
import argparse
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
from fractions import Fraction as F
from math import comb
from tools import strict_json
from . import readout as r
from model.native_local_response import packet as inherited

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
SOURCE_FILES=tuple('model/native_sector_readout/'+name for name in (
    '__init__.py','readout.py','packet.py','test_readout.py','PROOF.md','README.md'))+(
    'tools/strict_json.py','model/requirements.txt')


def source_bindings():
    paths=set(SOURCE_FILES)|set(inherited.source_bindings());pending=list(paths);visited=set()
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
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths) if not p.startswith('docs/')}


def examples():
    rows=[]
    local=r.program();value,family,total=r.specification(local)
    x=(1<<0)|(1<<11);h=r.trace_word(value,x,tuple(i%2 for i in range(11)))
    r.authenticate_trace(value,h)
    rows.append({'carrier':value['carrier'],'ports':12,'total':2,
                 'program_digest':r.source.digest(value),'attempts':11,
                 'actual_full_raw_word_authenticated':True,'weight':str(h['weight'])})
    for total in (0,1,2,120,239,240):
        value=r.program(level=0,total=total);family=r.cone.source_family(value)
        mask=(1<<total)-1;coins=tuple((i+total)%2 for i in range(239))
        h=r.trace_word(value,mask,coins);r.authenticate_trace(value,h)
        if r.reverse_currents(family,h['terminal_mask'],h['records'])[0]!=mask:
            raise ArithmeticError('terminal readback lost the complete current reversal')
        rows.append({'carrier':value['carrier'],'ports':240,'total':total,
                     'program_digest':r.source.digest(value),'attempts':239,
                     'actual_full_raw_word_authenticated':True,'weight':str(h['weight'])})
    return rows


def attachment():
    examples=[]
    for level in (0,1,2):
        ports=240*4**level;value=r.program(level=level,total=ports//2)
        family=r.cone.source_family(value);edges=r.tree(family)
        if len(edges)!=ports-1 or len({i for e in edges for i in e})!=ports:
            raise ArithmeticError('actual tree does not span every declared binary port')
        examples.append({'level':level,'ports':ports,'tree_commands':len(edges),
                         'tree_sha256':r.source.digest(edges),'program_digest':r.source.digest(value),
                         'record_resources':value['declared_resources'],
                         'available_body_AST_nodes':r.source.MAX_NODES-len(edges),
                         'fixed_intercell_program_admitted':True})
    return examples


def small_full_instrument():
    family=r.native.Family((2,3,5),((0,1),(1,2)),(1,1));rows=[]
    totals=[F() for _ in range(8)];counts=[0]*8;ambiguous=0
    for a in range(25):
        for b in range(25):
            weight,pairs,_=r.branch(family,(a,b))
            for x,y in pairs:totals[x]+=weight;counts[x]+=1
            if len(pairs)==2:
                if pairs!=((0,0),(7,7)):raise ArithmeticError('nonconstant tree ambiguity')
                ambiguous+=1
    if totals!=[F(1)]*8 or counts!=[4]*8 or ambiguous!=4:
        raise ArithmeticError('complete all-number detector lost probability or degeneracy')
    for m in range(4):
        nonzero=sum(bool(r.branch(family,(a,b),m)[1]) for a in range(25) for b in range(25))
        if nonzero!=comb(3,m)*4:raise ArithmeticError('known-number fine-word count changed')
        rows.append({'total':m,'nonzero_words':nonzero})
    return {'all_binary_source_columns':8,'formal_words':625,'coin_words_per_input':4,
            'unknown_total_rank_two_words':4,'known_total_sectors':rows,
            'empty_full_coherence_preserved':True}


def boundaries():
    value=r.program();stop=r.source.execute(value,{(0,0,0,0):F(1)})
    if stop.status!='RESOURCE_STOP' or stop.attempted_steps!=3:
        raise ArithmeticError('expanded host boundary changed')
    control=r.source.program(commands=[
        {'op':'repair','selector':{'kind':'fixed','seam':[0,1]},'law':'E','repeat':1},
        {'op':'repair','selector':{'kind':'last_current','positive':[2,3],'otherwise':[2,10]},'law':'C','repeat':1}])
    return {'expanded_status':stop.status,'expanded_attempts':stop.attempted_steps,
            'expanded_raw_blocks':len(stop.branches),'admission_nodes':r.source.MAX_NODES,
            'expanded_steps':r.source.MAX_EVALUATED_STEPS,'expanded_branches':r.source.MAX_EVALUATED_BRANCHES,
            'weighted_adaptive_reset_lower_bound':r.no_reset_bound(control),
            'fair_reset_lower_bound':r.no_reset_bound(r.program(level=0,total=1)),
            'unknown_total_constant_subspace_dimension':2,
            'known_total_is_read_from_source_specification':True}


def produce():
    bindings=source_bindings()
    result={'schema':'oph.native_sector_readout.v1','actual_tree_programs':attachment(),
            'authenticated_sparse_full_histories':examples(),'all_number_control':small_full_instrument(),
            'source_and_preparation_boundaries':boundaries(),
            'scope':{'all_binary_number_sectors_in_general_tree_proof':True,
                     'actual_federation_geometry_levels':[0,1,2],
                     'original_kernel_word_queries_not_normalized_sampler':True,
                     'arbitrary_untouched_reference_preserved':True,
                     'unknown_local_total_empty_full_coherence_retained':True,
                     'outcome_labelled_terminal_basis_preparation':True,
                     'deterministic_fixed_target_reset_derived':False,
                     'initial_census_supplied_to_terminal_readback':False,
                     'new_coherent_register_or_controller_policy_added':False,
                     'generic_host_limits_changed':False,'physical_QFT_or_particle_masses_derived':False},
            'source_bindings':bindings}
    if source_bindings()!=bindings:raise ArithmeticError('source readout inputs changed during replay')
    return result


def canonical(value):return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path)))!=canonical(produce()):
        raise ValueError('full-number native readout receipt differs from replay')
    return True


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);m=parser.add_mutually_exclusive_group(required=True)
    m.add_argument('--write',action='store_true');m.add_argument('--check',action='store_true');args=parser.parse_args()
    if args.write:OUTPUT.write_text(canonical(produce()));print('Native full-number readout receipt written.')
    else:verify();print('Native full-number readout receipt verified.')
