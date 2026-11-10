"""Typed exact-template admission for unchanged complete native producer outputs."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import copy,hashlib,json,math,re
import numpy as np
from model.protected_particle_readout import packet,experiment
from model.cube_eta_protection import repair
from . import kernel,interval_replay,arithmetic

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent

def require(test,message):
    if not test:raise ValueError(message)

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def specification():
    spec=json.loads((HERE/'contract.json').read_text())
    require(digest(ROOT/spec['retained_receipt'])==spec['retained_receipt_sha256'],'retained JSON changed')
    require(digest(ROOT/spec['retained_array'])==spec['retained_array_sha256'],'retained complete NPZ changed')
    require(digest(ROOT/spec['original_producer'])==spec['original_producer_sha256'],'original producer changed')
    old=json.loads((ROOT/spec['retained_receipt']).read_text())
    expected=dict(old['source_bindings'])
    for path,transition in spec['typed_parent_metadata_transitions'].items():
        require(expected[path]==transition['retained'],'parent baseline changed')
        expected[path]=transition['current']
    require(packet.source_bindings()==expected,'unreviewed executable, evidence or source scope changed')
    canonical=json.loads((ROOT/'model/protected_particle_readout/results.json').read_text())
    verify_legacy_projection(old,canonical,expected)
    return spec,old,expected

def verify_legacy_projection(old,canonical,pins):
    require(canonical['source_bindings'] in (old['source_bindings'],pins),'legacy source projection is unreviewed')
    restored=copy.deepcopy(canonical);restored['source_bindings']=old['source_bindings']
    require(restored==old,'legacy scientific fields or archive metadata changed')


def retained():
    spec,value,_=specification()
    with np.load(ROOT/spec['retained_array'],allow_pickle=False)as archive:
        arrays={key:archive[key]for key in archive.files}
    value=copy.deepcopy(value);value.pop('complete_record_archive')
    return value,arrays


def validate(value,arrays,*,budget=None):
    """Validate each raw law independently; no renormalization or zeroing."""
    packet.validate_accounting(value,arrays)
    budget=arithmetic.admitted_budget()if budget is None else budget
    E=F(budget['complete_state_error_upper']);norm=F(budget['normalization_error_upper'])
    terminal=F(budget['one_terminal_probability_upper']);failure=F(budget['one_export_probability_upper'])
    intervals=interval_replay.compare_arrays(arrays)
    for row in intervals:require(F(row['complete_state_error_upper'])<=E,'full state exceeds exact-template enclosure')
    law=value['executed_native_readout'];rows=law['all130_repair_records']
    for row in rows:
        p=F(row['probability'])
        if row['record'][0]=='outside':require(abs(p-1)<=norm,'outside probability violates protected template')
        else:
            require(0<=p<=terminal,'terminal repair probability exceeds certified enclosure')
            seam=tuple(repair.PORTS.index(j)for j in row['native_seam'])
            matrix=dict(repair.matrices(row['law'],seam))[tuple(row['record'])]
            if not np.any(matrix):require(p==0,'structurally zero native outcome changed')
    require(F(law['terminal_native_repair_probability'])<=120*terminal,'total native mass outside enclosure')
    records=law['source_prefix_records']+[x['analyzer_source_record']for x in law['pointer_settings']]
    for row in records:
        table=arrays[row['failure_weight_array']]
        require(float(np.max(table,initial=0))<=float(failure),'source export entry exceeds primitive enclosure')
        require(0<=F(row['first_failure_probability'])<=row['native_abort_slots']*failure,'source export total exceeds enclosure')
        require(0<=row['failure_sum_error_diagnostic']<=row['failure_sum_admission'],'source sum diagnostic exceeds own admission')
    # Pair-code columns2/3 are padding. The actual producer validates this
    # schema; their exact zero values are checked from its emitted tape here.
    from model.cube_eta_protection import source
    from model.protected_particle_readout import circuits
    tapes={'tree':source.tree_word(),'copy':source.copy_word(),'pointer':circuits.pointer_prepare(),
           'field_first':circuits.flatten(circuits.field_binding(F(1,2),32)),
           'analyzer_0':circuits.pointer_analyzer(F()),'analyzer_1':circuits.pointer_analyzer(F(1,2))}
    tapes['field_second']=tapes['field_first']
    for name,tape in tapes.items():
        pair=[j for j,event in enumerate(tape)if len(event['targets'])==2]
        require(np.count_nonzero(arrays[name+'_failure_weights'][pair,2:])==0,'padding is not a physical source output')
    prefix=sum(x['first_failure_probability']for x in law['source_prefix_records'])+law['terminal_native_repair_probability']
    for j,setting in enumerate(law['pointer_settings']):
        require(all(p==0 for p in setting['pointer_probabilities'][2:]),'impossible final export complement changed')
        require(abs(F(setting['complete_probability'])-1)<=norm,'complete normalization exceeds arithmetic bound')
        lo,hi=map(F,intervals[j]['exact_template_pointer_zero_interval'])
        require(lo-norm<=F(setting['pointer_probabilities'][0])<=hi+norm,'pointer event outside exact-template enclosure')
        target=experiment.pointer_outputs(experiment.ideal_analyzed(F(1),experiment.PHASES[j]))
        actual=experiment.pointer_outputs(arrays[f'analyzer_{j}_full_success_state'])
        recomputed=math.sqrt(sum(float(np.linalg.norm(a-b)**2)for a,b in zip(actual,target,strict=True))
                            +prefix+setting['analyzer_source_record']['first_failure_probability'])
        require(abs(recomputed-setting['complete_output_error_diagnostic'])<=packet.REPLAY_ATOL,
                'output-error diagnostic differs from its own full array')
    for j,row in enumerate(law['executed_zero_field_duration_control']):
        desired=(F(1),F())if j==0 else(F(1,2),F(1,2))
        for actual,target in zip(row['all_six_pointer_probabilities'][:2],desired,strict=True):
            require(abs(F(actual)-target)<=norm,'zero-duration pointer outside exact control')
        require(all(p==0 for p in row['all_six_pointer_probabilities'][2:]),'zero-duration export complement changed')
        require(0<=F(row['pointer_source_first_failures'])<=160*failure,'zero-duration source failures exceed enclosure')
    return intervals


def numerical_paths(value):
    base='receipt.executed_native_readout';paths={base+'.terminal_native_repair_probability'}
    paths.update(f'{base}.all130_repair_records[{i}].probability'for i in range(130))
    fields=('first_failure_probability','failure_sum_error_diagnostic')
    for i in range(5):paths.update(f'{base}.source_prefix_records[{i}].{x}'for x in fields)
    for i in range(2):
        p=f'{base}.pointer_settings[{i}]'
        paths.update(p+'.'+x for x in ('complete_probability','complete_output_error_diagnostic'))
        paths.update(p+'.analyzer_source_record.'+x for x in fields)
        paths.update(p+f'.pointer_probabilities[{j}]'for j in range(6))
        p=f'{base}.executed_zero_field_duration_control[{i}]'
        paths.add(p+'.pointer_source_first_failures')
        paths.update(p+f'.all_six_pointer_probabilities[{j}]'for j in range(6))
    return paths


def compare_science(old,fresh):
    """Only independently validated diagnostic fields use the successor rule."""
    _,_,pins=specification();require(fresh['source_bindings']==pins,'fresh source pins changed')
    baseline=copy.deepcopy(old);baseline['source_bindings']=pins
    allowed=numerical_paths(old)
    def walk(a,b,path='receipt'):
        require(type(a)is type(b),'scientific type changed: '+path)
        if type(a)is dict:
            require(a.keys()==b.keys(),'scientific keys changed: '+path)
            for key in a:walk(a[key],b[key],path+'.'+key)
        elif type(a)is list:
            require(len(a)==len(b),'scientific list changed: '+path)
            for i,(x,y)in enumerate(zip(a,b,strict=True)):walk(x,y,f'{path}[{i}]')
        elif path in allowed:
            require(type(a)is float and math.isfinite(a)and math.isfinite(b),'finite diagnostic required')
        else:packet.compare(a,b,path)
    walk(baseline,fresh)


def backend_certificate():
    budget=arithmetic.admitted_budget();a=arithmetic.certificate()
    cap=F(budget['local_matrix_error_cap'])
    for row in a['local_success_block_bounds']:
        require(F(row['exact_dyadic_Frobenius_error_upper'])<=cap,'field local matrix misses backend accuracy contract')
    require(F(a['pointer_success_error_upper'])<=cap,'pointer local matrix misses accuracy contract')
    require(all(F(x)<=cap for x in a['analyzer_success_error_uppers']),'analyzer local matrix misses accuracy contract')
    require(F(a['copied_initial_state_error_upper'])<=F(budget['initial_state_error_cap']),'copied preparation misses accuracy contract')
    require(F(a['V_embedding_error_upper'])<=F(budget['embedding_error_cap']),'repair embedding misses accuracy contract')
    primitives=arithmetic.primitive_certificate()
    require(all(F(row['complete_M6_error_upper'])<=cap for row in primitives),'primitive misses complete-M6 accuracy contract')
    return a,primitives,budget


def report(fresh_pair=None):
    old,arrays=retained();k=kernel.certificate();backend,primitives,budget=backend_certificate()
    historical=validate(old,arrays,budget=budget)
    value,current_arrays=packet.report()if fresh_pair is None else fresh_pair
    compare_science(old,value);current=validate(value,current_arrays,budget=budget)
    rejected=False
    try:
        baseline=copy.deepcopy(old);baseline['source_bindings']=value['source_bindings'];packet.compare(baseline,value)
    except ValueError:rejected=True
    return {'schema':'oph.protected_readout_numerics.v1','scope':'same finite17-bit program; corrected numerical replay contract',
            'physical_statuses_changed':False,'old_checker_and_admission_preserved':True,
            'exact_template_kernel':k,'fixed_numerical_budget':budget,
            'retained_full_state_interval_validation':historical,
            'backend_replay_diagnostics':{'old_checker_rejected_this_replay':rejected,'fresh_full_state_interval_validation':current,
                'arithmetic_certificate':backend,'primitive_certificates':primitives,
                'fresh_raw_complete_law':value['executed_native_readout']},
            'observed_universe_or_hardware_accuracy_claimed':False}


def live_projection():
    spec,old,pins=specification()
    result=copy.deepcopy(old);result['source_bindings']=pins
    return result
