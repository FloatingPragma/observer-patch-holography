"""Replay source-bound local readout, erasure and admission certificates."""
import argparse
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path
from tools import strict_json
from . import readout as r
from model.native_sector_readout import packet as sector
from model.native_rooted_refinement import packet as rooted

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
SOURCE_FILES=tuple('model/native_fair_local_readout/'+name for name in (
    '__init__.py','readout.py','packet.py','test_readout.py','PROOF.md','README.md'))+(
    'Lean/ObserverPatchHolography/Locality/NoSignalling.lean',
    'tools/strict_json.py','model/requirements.txt')


def source_bindings():
    paths=set(SOURCE_FILES)|set(sector.source_bindings())|set(rooted.source_bindings())
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths) if not p.startswith('docs/')}


def encode(value):
    if type(value) is F:return str(value)
    if type(value) in (tuple,list):return [encode(x) for x in value]
    if type(value) is dict and all(type(k) is str for k in value):return {k:encode(v) for k,v in value.items()}
    if value is None or type(value) in (int,str,bool):return value
    raise ValueError('typed exact JSON output required')


def canonical(value):return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def all_number_control():
    ports=(0,1,2);edges=((0,1),(1,2));p=r.plan(ports,edges);steps=3
    blocks={};totals=[F()]*8;weight=F(1,4**steps)
    for word in product(edges,repeat=steps):
        for coins in product((0,1),repeat=steps):
            for initial in range(8):
                state=initial;records=[]
                for t,(edge,coin) in enumerate(zip(word,coins,strict=True),1):
                    a,b=edge;rec=r.record(edge,coin,((state>>a)&1)-((state>>b)&1))
                    records.append((t,rec));state+=(1<<b)*rec.current-(1<<a)*rec.current
                key=tuple(rec for _,rec in records)
                blocks.setdefault(key,[]).append((initial,state));totals[initial]+=weight
                result=r.decode(p,tuple(records),steps)
                if result['status']=='DECODED' and (result['initial_mask'],result['terminal_mask'])!=(initial,state):
                    raise ArithmeticError('record reconstruction failed')
    failures=[F()]*8;constant=[F()]*8;success=0
    for labels,pairs in blocks.items():
        result=r.decode(p,tuple(enumerate(labels,1)),steps)
        if len({x for x,y in pairs})!=len(pairs) or len({y for x,y in pairs})!=len(pairs):
            raise ArithmeticError('native history lost its partial-injection property')
        if result['status']=='DECODED':
            if {x for x,y in pairs}!={result['initial_mask']}:raise ArithmeticError('K_h=K_h P_x failed')
            success+=1
        else:
            for x,y in pairs:failures[x]+=weight
            if result['status']=='CONSTANT_ERASURE':
                if {x for x,y in pairs}!={0,7}:raise ArithmeticError('constant block changed')
                for x,y in pairs:constant[x]+=weight
    pcov=r.coverage(2,2,steps,exact=True)['exact']
    if totals!=[F(1)]*8 or failures!=[F(1) if x in (0,7) else pcov for x in range(8)]:
        raise ArithmeticError('complete-instrument erasure probabilities differ')
    return {'binary_input_columns':8,'local_numbers':list(range(4)),
        'formal_raw_words':(1+24*len(edges))**steps,'live_fine_words':len(blocks),
        'decoded_fine_words':success,'coverage_error':pcov,
        'erasure_by_initial_mask':failures,'constant_erasure_by_initial_mask':constant,
        'completeness_by_initial_mask':totals,
        'all_live_fine_maps_partial_injections_checked':True,
        'decoded_fine_map_initial_projection_identity_checked':True}


def actual_kernel_controls():
    p=r.source_plan(0,tuple(range(12)));fixed=r.fixed_program(0,3,p)
    initial=3+(1<<20);other=3+(1<<21);state=initial;slots=[]
    for edge in p['tree']:
        a,b=edge;coin=len(slots)%2;rec=r.record(edge,coin,((state>>a)&1)-((state>>b)&1),'fixed')
        rows=r.source._outcomes(fixed,fixed['commands'][len(slots)],{})
        slots.append(next(i for i,out in enumerate(rows) if r.inherited._same_record(out.record,rec)))
        state+=(1<<b)*rec.current-(1<<a)*rec.current
    density={(initial,0,other,1):F(2,3),(other,1,initial,0):F(2,3)}
    out=r.raw_query(fixed,p,tuple(slots),density)
    paired=state-(1<<20)+(1<<21)
    expected={(state,0,paired,1):F(2,3*2**11),(paired,1,state,0):F(2,3*2**11)}
    if out['density']!=expected or out['decoded']['initial_mask']!=3:
        raise ArithmeticError('the original kernel lost exterior/reference coherence')
    zero=r.raw_query(fixed,p,(0,)*11,density)
    if zero['density'] or len(zero['archive'])!=10:raise ArithmeticError('zero raw word was omitted')
    uniform=r.uniform_program(0,1,1);rho={(1,0,1,0):F(1)}
    rows=r.source._outcomes(uniform,uniform['commands'][0],rho)
    if len(rows)!=15121 or sum(sum(x for (a,i,b,j),x in row.density.items() if a==b and i==j) for row in rows)!=1:
        raise ArithmeticError('complete original uniform raw alphabet changed')
    slot=next(i for i,row in enumerate(rows) if row.density and row.record.seam in p['incident_seams'])
    queried=r.raw_query(uniform,p,(slot,),rho)
    if queried['density']!=rows[slot].density or queried['memory']!=rows[slot].memory:
        raise ArithmeticError('raw-word query differs from original continuation')
    stop=r.source.execute(fixed,{(initial,0,initial,0):F(1)})
    if stop.status!='RESOURCE_STOP' or stop.attempted_steps!=3:raise ArithmeticError('expanded stop changed')
    contradictory_slots=[];contradictory_records=[]
    for command in fixed['commands']:
        rec=r.record(tuple(command['selector']['seam']),1,1,'fixed')
        options=r.source._outcomes(fixed,command,{})
        contradictory_slots.append(next(i for i,row in enumerate(options) if r.inherited._same_record(row.record,rec)))
        contradictory_records.append(rec)
    contradictory=[]
    for informed in (False,True):
        zero=r.raw_query(fixed,p,tuple(contradictory_slots),density,use_global_total=informed)
        if zero['density'] or zero['decoded']['status']!='ZERO_RAW_WORD' or not r.inherited._same_records(
                zero['archive']+zero['memory'],tuple(contradictory_records)):
            raise ArithmeticError('a contradictory complete binary word lost its zero block or raw records')
        contradictory.append({'uses_declared_global_total':informed,'raw_slots':contradictory_slots,
            'decoded_status':zero['decoded']['status'],'quantum_block_is_zero':True,
            'archive_records':len(zero['archive']),'FIFO_records':len(zero['memory'])})
    return {'fixed_program_digest':r.source.digest(fixed),'fixed_attempts':11,
        'fixed_raw_slots_per_attempt':25,'fixed_reference_word_weight':F(1,2**11),
        'initial_local_mask':3,'decoded_terminal_local_mask':out['decoded']['terminal_mask'],
        'archive_records':10,'FIFO_records':1,'exact_off_diagonal_reference_control':True,
        'zero_raw_word_retained':True,'uniform_program_digest':r.source.digest(uniform),
        'uniform_full_raw_slots':len(rows),'uniform_selected_slot':slot,
        'expanded_status':stop.status,'expanded_attempts':stop.attempted_steps,
        'contradictory_binary_word_controls':contradictory}


def global_specification_controls():
    p=r.source_plan(0,tuple(range(12)));rows=[]
    for total,mask,label in ((1,1<<20,0),(239,((1<<240)-1)^(1<<20),4095)):
        value=r.fixed_program(0,total,p);slots=[]
        for command in value['commands']:
            expected=r.record(tuple(command['selector']['seam']),len(slots)%2,0,'fixed')
            outcomes=r.source._outcomes(value,command,{})
            slots.append(next(i for i,out in enumerate(outcomes) if r.inherited._same_record(out.record,expected)))
        density={(mask,0,mask,0):F(1)}
        ordinary=r.raw_query(value,p,tuple(slots),density)
        informed=r.raw_query(value,p,tuple(slots),density,use_global_total=True)
        if ordinary['decoded']['status']!='CONSTANT_ERASURE' or informed['decoded']['initial_mask']!=label:
            raise ArithmeticError('global-total feasibility failed to resolve the unique constant mask')
        for field in ('density','memory','archive','observed_events'):
            if ordinary[field]!=informed[field]:raise ArithmeticError('label refinement changed the native instrument')
        rows.append({'global_total':total,'program_digest':r.source.digest(value),
            'feasibility':r.feasibility(240,total,12),'initial_local_mask':label,
            'raw_word_weight':F(1,2**11),'complete_quantum_and_record_output_unchanged':True,
            'canonical_ambiguous_weight':r.canonical_erasure(240,total,12,use_global_total=True)['exact']})
    extremes=[]
    for total,label in ((0,0),(240,4095)):
        result=r.decode(p,(),0,population=240,total=total)
        if result['status']!='DECODED' or result['initial_mask']!=label:
            raise ArithmeticError('globally constant source needs no measurement coverage')
        extremes.append({'global_total':total,'tree_hits_required':0,'initial_local_mask':label,
                         'oracle_error':r.error_bound(F(1),F(),global_constant=True)})
    zero_rows=[]
    for total in (0,240):
        value=r.uniform_program(0,total,1);mask=0 if total==0 else (1<<240)-1
        density={(mask,0,mask,0):F(1)}
        outcomes=r.source._outcomes(value,value['commands'][0],density)
        slot=next(i for i,out in enumerate(outcomes) if out.record.kind=='repair'
            and out.record.seam in p['incident_seams'] and out.record.law=='E' and out.record.label==1)
        result=r.raw_query(value,p,(slot,),density,use_global_total=True)
        if result['density'] or result['decoded']['status']!='ZERO_RAW_WORD' or result['memory']!=outcomes[slot].memory:
            raise ArithmeticError('a valid zero-map global-sector word was rejected or lost its record')
        zero_rows.append({'global_total':total,'raw_slot':slot,'decoded_status':result['decoded']['status'],
                          'quantum_block_is_zero':True,'complete_raw_record_retained':True})
    return {'source_bound_unique_constant_controls':rows,'globally_constant_controls':extremes,
            'infeasible_global_sector_zero_word_controls':zero_rows,
            'ambiguous_both_feasible_control':r.feasibility(240,80,12)}


def routed_coverage():
    level=6;outer=r.routing.routed_window(level,5);inner=r.routing.routed_window(level,4)
    inner_set=set(inner['source_ports'])
    indices=tuple(i for i,x in enumerate(outer['source_ports']) if x in inner_set)
    p=r.rooted_plan(outer,indices);V=240*4**level;M=630*4**level;T=64*M
    if len(p['ports'])!=348 or len(p['tree'])!=347:raise ArithmeticError('the finite readout window changed')
    value=r.uniform_program(level,V//3,T)
    cov=r.coverage(M,len(p['tree']),T);eta=r.canonical_erasure(V,V//3,len(p['ports']))
    eta_bound=F(1,10**60);error=r.error_bound(cov['upper'],eta_bound)
    if eta['product_upper']>=eta_bound or error>=F(1,10**16):raise ArithmeticError('finite readout bound failed')
    try:r.fixed_program(level,V//3,p)
    except ValueError as exception:
        if 'intercell' not in str(exception):raise
    else:raise ArithmeticError('high-level fixed intercell programme unexpectedly admitted')
    high=200;ceiling=r.uniform_ceiling(high,80);highvalue=r.uniform_program(high,80,ceiling)
    stopped=r.source.execute(highvalue,{'sentinel':True})
    if stopped.status!='RESOURCE_STOP_UNEVALUATED_CARRIER' or stopped.attempted_steps:
        raise ArithmeticError('unevaluated source host did not preserve its stop')
    return {'level':level,'global_ports':V,'global_seams':M,'global_total':V//3,
        'detector_cell_radius':4,'incident_window_cell_radius':5,'detector_ports':348,
        'tree_edges':347,'incident_seams':len(p['incident_seams']),
        'plan_sha256':r.source.digest(p),'program_digest':r.source.digest(value),
        'attempts':T,'per_seam_attempt_time':F(T,M),'resources':value['declared_resources'],
        'coverage_upper':cov['upper'],'canonical_constant_weight_upper':eta_bound,
        'full_CQ_oracle_error_upper':error,'displayed_error_gate':F(1,10**16),
        'fixed_intercell_rejection_checked':True,'expanded_geometry_executed':False,
        'large_capacity_admission_control':{'level':high,'attempt_ceiling':ceiling,
            'bits_per_record':highvalue['declared_resources']['record_bits_per_slot_upper'],
            'coverage_lower':r.coverage(630*4**high,11,ceiling)['lower'],
            'host_status':stopped.status,'attempted_steps':stopped.attempted_steps}}


def produce():
    before=source_bindings()
    output={'schema':'oph.native_fair_local_readout.v1','source_bindings':before,
        'complete_small_uniform_instrument':all_number_control(),
        'actual_source_kernel_controls':actual_kernel_controls(),
        'declared_global_total_refinement':global_specification_controls(),
        'finite_routed_uniform_readout':routed_coverage(),
        'scope':{'all_local_number_sectors':True,'arbitrary_untouched_reference':True,
            'optional_source_bound_global_total_feasibility':True,'boundary_currents_required':True,'constant_input_coherence_preserved':True,
            'oracle_is_comparison_not_implemented_service':True,
            'coverage_error_is_unconditional':True,'coherent_service_lift_derived':False,
            'all_capacity_admitted_coverage_guarantee':False,
            'expanded_host_limits_changed':False,'physical_clock_or_particle_identification':False}}
    if before!=source_bindings():raise ArithmeticError('readout sources changed during replay')
    return encode(output)


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path)))!=canonical(produce()):
        raise ValueError('local fair-source readout receipt differs from replay')
    return True


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);action=parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write',action='store_true');action.add_argument('--check',action='store_true');args=parser.parse_args()
    if args.write:OUTPUT.write_text(canonical(produce()));print('Native fair local readout receipt written.')
    else:verify();print('Native fair local readout receipt verified.')
