"""Frozen finite controller replay, complete records and strict receipt checks."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from . import controller as c

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
TOL=2e-10
PROB_TOL=2e-11
SOURCE_PATHS=tuple('model/experiments/contrast_clock/'+name for name in
                   ('__init__.py','__main__.py','controller.py','experiment.py','protocol.json'))+tuple(
    'model/candidate_laws/contrast_memory_spherical_v1/'+name for name in ('__init__.py','common.py','law.py','reader.py','sectors.py'))+(
    'model/executable_model/carrier.py',)


def require(condition,message):
    if not condition:
        raise ValueError(message)


def _pairs(items):
    result={}
    for key,value in items:
        require(key not in result,f'duplicate JSON key: {key}')
        result[key]=value
    return result


def _finite_tree(value):
    if isinstance(value,float):
        require(math.isfinite(value),'nonfinite JSON number')
    elif isinstance(value,dict):
        for child in value.values():
            _finite_tree(child)
    elif isinstance(value,list):
        for child in value:
            _finite_tree(child)


def read_json(path):
    value=json.loads(Path(path).read_text(encoding='utf8'),object_pairs_hook=_pairs,
                    parse_constant=lambda s: (_ for _ in ()).throw(ValueError(f'nonfinite {s}')))
    _finite_tree(value)
    return value


def source_pins():
    frozen=read_json(HERE/'FREEZE.json')
    require(set(frozen)==set(SOURCE_PATHS),'source pin coverage')
    for rel,digest in frozen.items():
        require(hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest,f'source mismatch: {rel}')
    return frozen


def protocol():
    p=read_json(HERE/'protocol.json')
    require(p['construction']=='contrast-memory-clock-v1','construction')
    require(p['angle']==1/8 and p['rounds_per_read']==8,'schedule')
    require(p['read_seams']==[[0,1],[0,5]] and p['laws']==[8.,0.,-8.],'source choices')
    require(p['apparatus_controls']==list(c.CONTROLS),'apparatus')
    require(p['empirical_selection'] is False,'empirical ancestry')
    require(p['gates']['operator_residual']==TOL,'operator gate')
    require(p['gates']['probability_residual']==PROB_TOL,'probability gate')
    require(p['gates']['maximum_source_gate_operator_error']==1e-6,'precision gate')
    return p


def public(words):
    return [[(word>>6)&3,(word>>8)&7] for word in words]


def readiness_support(words,step):
    decisions=[(word>>6)&3 for word in words]
    require((decisions[0]!=0)==(step>=668) and (decisions[1]!=0)==(step>=1336),
            'archive readiness outside program support')


def run_case(tape,kappa,control):
    blocks,prefixes=c.run_history(tape,kappa,control,full_operator=True)
    states=c.source_sector()
    start=int(np.flatnonzero((states[:,0]==1)&(states[:,1]==1))[0])
    rows=[]
    complete=np.zeros((78,78),complex)
    for words,K in sorted(blocks.items()):
        complete+=K.conj().T@K
        vector=K[:,start]
        rows.append({'words':list(words),'public':public(words),
                     'all_input_weights':np.sum(abs(K)**2,axis=0).tolist(),
                     'prepared_output_real':vector.real.tolist(),
                     'prepared_output_imag':vector.imag.tolist()})
    residual=float(np.linalg.norm(complete-np.eye(78),2))
    # The nominal clock consumes every prefix. Controls retain complete final
    # instruments plus period/read-stage checkpoints, without a clock claim.
    marks={0,len(tape)}
    for k,step in enumerate(tape,1):
        if step.kind=='reader' or (step.kind=='seam' and step.ports==(0,1) and step.angle==1/16):
            marks.add(k)
    return {'kappa':kappa,'apparatus':control,'branches':rows,
            'complete_instrument_error':residual,
            'prefix_records':[row for row in prefixes if (kappa==8 and control=='nominal') or row['step'] in marks],
            'final_no_pair':prefixes[-1]['no_pair']},prefixes


def public_distribution(case):
    out={}
    for row in case['branches']:
        key=tuple(tuple(x) for x in row['public'])
        weight=sum(x*x+y*y for x,y in zip(row['prepared_output_real'],row['prepared_output_imag']))
        out[key]=out.get(key,0.)+weight
    return out


def tv(a,b):
    return sum(abs(a.get(k,0)-b.get(k,0)) for k in a.keys()|b.keys())/2


def clock_records(prefixes,length,fraction,kind='engineered'):
    count=len(prefixes)-1
    time=fraction*np.pi*length/2
    weights=(c.clock_law(length,time) if kind=='engineered'
             else abs(c.clock_amplitudes(length,time,kind))**2)
    # All cursor positions are retained, including unfinished and return positions.
    # The joint prefix law is weights[k] times prefixes[min(k,count)].
    public_law={}
    for k,w in enumerate(weights):
        for branch in prefixes[min(k,count)]['archive']:
            key=tuple(tuple(p) for p in public(branch['words']))
            public_law[key]=public_law.get(key,0.)+float(w)*branch['weight']
    return {'kind':kind,'transfer_fraction':fraction,'position_weights':weights.tolist(),
            'unfinished_weight':float(weights[:count].sum()),
            'public_records':[{'readings':[list(x) for x in key],'weight':value}
                              for key,value in sorted(public_law.items())]}


def run():
    pins=source_pins()
    p=protocol()
    tape=c.program()
    cases=[]
    nominal_prefixes=None
    for kappa,control in [(k,'nominal') for k in p['laws']]+[(8.,a) for a in p['apparatus_controls'][1:]]:
        case,prefixes=run_case(tape,kappa,control)
        cases.append(case)
        if kappa==8 and control=='nominal':
            nominal_prefixes=prefixes
    n=len(tape)
    length=2*n
    clocks=[clock_records(nominal_prefixes,length,f) for f in p['clock_fractions_of_transfer_time']]
    clocks.extend(clock_records(nominal_prefixes,length,1.,kind) for kind in ('uniform','broken'))
    comparisons=[{'kappa':row['kappa'],'apparatus':row['apparatus'],
                  'public_total_variation':tv(public_distribution(cases[0]),public_distribution(row))}
                 for row in cases[1:]]
    result={'construction':p['construction'],'source_pins':pins,
            'resources':resources(tape),
            'cases':cases,'competing_controls':comparisons,
            'clock_records':clocks,'literal_clock_check':c.direct_clock_check(),
            'claims':{'finite_autonomous_declared_extension':True,'absorbing_halt':False,
                      'physical_clock_or_particle_identification':False,
                      'continuous_time_binding_bound_transferred':False}}
    return result


def resources(tape):
    n=len(tape)
    length=2*n
    capacities=(c.count_table(tape)[-1]+1).tolist()
    source_gate_count=sum(step.kind!='reader' for step in tape)
    return {'source_qutrits':12,'private_seam_memories_11_states':30,
                         'degree_reference_qutrits':12,'reference_memories_11_states':12,
                         'fresh_reader_bits':22,'clock_bits_one_excitation':length+1,
                         'counter_capacities':capacities,
                         'source_gates':source_gate_count,'reader_stages':8,
                         'program_edges':n,'idle_edges':n,'clock_length':length,
                         'total_two_source_dimension':78,
                         'clock_operator_norm':1.,
                         'largest_clock_edge_coupling':float(c.couplings(length).max()),
                         'transfer_action':np.pi*length/2,
                         'completion_window_action':[np.pi*length/3,2*np.pi*length/3],
                         'window_unfinished_bound':math.exp(-length/8),
                         'coherent_gate_error_tv_bound':source_gate_count*1e-6}


def formula_compare(expected,actual):
    """Analytic resource formulas have no absolute numerical-error floor."""
    if type(expected) is float:
        require(type(actual) is float and math.isfinite(actual),'resource numeric type')
        require(abs(expected-actual)<=4*math.ulp(expected),'resource formula mismatch')
    elif isinstance(expected,dict):
        require(type(actual) is dict and expected.keys()==actual.keys(),'resource keys')
        for key in expected: formula_compare(expected[key],actual[key])
    elif isinstance(expected,list):
        require(type(actual) is list and len(expected)==len(actual),'resource list')
        for a,b in zip(expected,actual): formula_compare(a,b)
    else:
        require(type(expected) is type(actual) and expected==actual,'resource value')


def array(values,shape):
    require(isinstance(values,list),'expected list')
    def numbers(v):
        if isinstance(v,list):
            for child in v: numbers(child)
        else:
            require(type(v) in (float,int),'expected numeric array')
    numbers(values)
    a=np.asarray(values,dtype=float)
    require(a.shape==shape and np.isfinite(a).all(),'invalid numeric array')
    return a


def validate(result):
    """Scientific gates recomputed from records, independently of replay tolerance."""
    _finite_tree(result)
    r=result['resources']
    formula_compare(resources(c.program()),r)
    require(r['program_edges']==1336 and r['idle_edges']==1336 and r['clock_length']==2672,'resource count')
    require(r['counter_capacities']==[193]*11+[177],'counter allocation')
    require(len(result['cases'])==6 and len(result['clock_records'])==8,'complete campaign')
    require(len(result['competing_controls'])==5,'complete controls')
    states=c.source_sector()
    counter_words=c.count_table(c.program()).tolist()
    start=int(np.flatnonzero((states[:,0]==1)&(states[:,1]==1))[0])
    nominal_law=None
    for case in result['cases']:
        require(type(case['complete_instrument_error']) is float and 0<=case['complete_instrument_error']<TOL,
                'instrument completeness gate')
        all_weights=np.zeros(78)
        no_pair=0.
        mask=np.array([any(row[a]==row[b]==1 for a,b in c.cell().seams) for row in states])
        seen=set()
        for branch in case['branches']:
            words=branch['words']
            require(len(words)==2 and all(type(w) is int and 0<=w<c.APP for w in words),'archive words')
            require(tuple(words) not in seen,'duplicate archive branch')
            seen.add(tuple(words))
            require(branch['public']==public(words),'public decoder mismatch')
            weights=array(branch['all_input_weights'],(78,))
            require(np.all(weights>=0) and np.all(weights<=1+PROB_TOL),'outcome weights')
            if np.any(weights>0): readiness_support(words,1336)
            v=array(branch['prepared_output_real'],(78,))+1j*array(branch['prepared_output_imag'],(78,))
            require(abs(np.vdot(v,v).real-weights[start])<PROB_TOL,'prepared branch weight mismatch')
            no_pair+=float(np.vdot(v[~mask],v[~mask]).real)
            all_weights+=weights
        require(np.max(abs(all_weights-1))<PROB_TOL,'missing input outcome')
        require(abs(no_pair-case['final_no_pair'])<PROB_TOL,'no-pair mismatch')
        for row in case['prefix_records']:
            require(type(row['step']) is int and 0<=row['step']<=1336,'prefix position')
            require(type(row['no_pair']) is float and 0<=row['no_pair']<=1+PROB_TOL,'prefix no-pair probability')
            if row['step']==0:
                require(row['no_pair']==0.,'initial pair preparation')
            require(len(row['counts'])==12 and all(type(x) is int for x in row['counts']) and
                    row['counts']==counter_words[row['step']],'local counter words')
            seen_prefix=set()
            for branch in row['archive']:
                require(type(branch['weight']) is float and 0<=branch['weight']<=1+PROB_TOL,'prefix weight')
                words=branch['words']
                require(len(words)==2 and all(type(w) is int and 0<=w<c.APP for w in words),'prefix archive words')
                require(tuple(words) not in seen_prefix,'duplicate prefix archive word')
                seen_prefix.add(tuple(words))
                if branch['weight']>0: readiness_support(words,row['step'])
            require(abs(sum(b['weight'] for b in row['archive'])-1)<PROB_TOL,'prefix mass')
        if nominal_law is None:
            nominal_law=public_distribution(case)
            require([row['step'] for row in case['prefix_records']]==list(range(1337)),
                    'missing clock-conditioned prefix')
    for control,case in zip(result['competing_controls'],result['cases'][1:]):
        actual=tv(nominal_law,public_distribution(case))
        require(abs(actual-control['public_total_variation'])<PROB_TOL,'control recomputation')
        require(actual>0.001,'competing control lacks detectable consequence')
    for row in result['clock_records']:
        weights=array(row['position_weights'],(2673,))
        require(np.all(weights>=0) and abs(weights.sum()-1)<PROB_TOL,'clock probability mass')
        require(type(row['unfinished_weight']) is float and
                float(weights[:1336].sum())==row['unfinished_weight'],'unfinished recomputation')
        require(abs(sum(v['weight'] for v in row['public_records'])-1)<PROB_TOL,'clock public mass')
        require(all(type(v['weight']) is float and 0<=v['weight']<=1+PROB_TOL
                    for v in row['public_records']),'clock public weight')
        recomputed={}
        prefixes=result['cases'][0]['prefix_records']
        for k,w in enumerate(weights):
            for branch in prefixes[min(k,1336)]['archive']:
                key=tuple(tuple(p) for p in public(branch['words']))
                recomputed[key]=recomputed.get(key,0.)+float(w)*branch['weight']
        retained={tuple(tuple(x) for x in r['readings']):r['weight'] for r in row['public_records']}
        require(len(retained)==len(row['public_records']),'duplicate clock public record')
        require(retained.keys()==recomputed.keys() and tv(retained,recomputed)<PROB_TOL,
                'clock public joint-law recomputation')
        require(all(retained[key]==0. for key,value in recomputed.items() if value==0.),
                'impossible clock public support')
        unavailable=sum(value for key,value in retained.items() if key[1][0]==0)
        require(abs(unavailable-row['unfinished_weight'])<=PROB_TOL*row['unfinished_weight'],
                'readiness marginal exceeds clock failure')
        if row['kind']=='engineered' and 2/3<=row['transfer_fraction']<=4/3:
            require(row['unfinished_weight']<=math.exp(-2672/8),'completion window gate')
        if row['kind']=='broken':
            require(row['unfinished_weight']>1-PROB_TOL,'broken clock must retain unfinished outcomes')
    errors=array(result['literal_clock_check']['operator_errors'],(3,))
    require(np.all(errors>=0) and np.max(errors)<TOL,'literal clock operator gate')


def compare(expected,actual,path='result',clock_kind=None):
    if type(expected) is float:
        require(type(actual) is float and math.isfinite(actual),f'{path}: float required')
        source_weight=('.all_input_weights[' in path or
                       ('.archive[' in path and path.endswith('.weight')))
        clock_weight=('.position_weights[' in path or path.endswith('.unfinished_weight') or
                      ('.public_records[' in path and path.endswith('.weight')))
        # Source archive zeros come from exact predicate/permutation support.
        # The engineered binomial law also has no eigensolver cancellation.
        # Numerical competing clocks retain their ordinary replay allowance.
        if expected==0. and (source_weight or (clock_weight and clock_kind=='engineered')):
            require(actual==0.,f'{path}: impossible probability support')
        elif clock_weight and clock_kind=='engineered':
            # Preserve resolved tiny tails; only the last eight subnormal ULPs
            # have an absolute allowance for probability underflow rounding.
            require(abs(expected-actual)<=PROB_TOL*abs(expected)+8*math.ulp(0.),
                    f'{path}: relative probability replay mismatch')
        else:
            require(abs(expected-actual)<=PROB_TOL*(1+abs(expected)),f'{path}: replay mismatch')
    elif isinstance(expected,dict):
        require(type(actual) is dict and expected.keys()==actual.keys(),f'{path}: keys')
        if 'position_weights' in expected:
            clock_kind=expected['kind']
        for key in expected: compare(expected[key],actual[key],f'{path}.{key}',clock_kind)
        if 'prepared_output_real' in expected and 'prepared_output_imag' in expected:
            # Test a complex zero jointly. A zero real or imaginary part alone
            # can acquire ordinary floating-point phase roundoff.
            for re,im,ar,ai in zip(expected['prepared_output_real'],expected['prepared_output_imag'],
                                   actual['prepared_output_real'],actual['prepared_output_imag']):
                if re==im==0.:
                    require(ar==ai==0.,f'{path}: impossible prepared-output support')
    elif isinstance(expected,list):
        require(type(actual) is list and len(expected)==len(actual),f'{path}: length')
        for k,(e,a) in enumerate(zip(expected,actual)): compare(e,a,f'{path}[{k}]',clock_kind)
    else:
        require(type(expected) is type(actual) and expected==actual,f'{path}: value')


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('--record',type=Path)
    parser.add_argument('--verify',type=Path,default=HERE/'results.json')
    args=parser.parse_args()
    result=run()
    if args.record:
        with args.record.open('x',encoding='utf8',newline='\n') as stream:
            json.dump(result,stream,separators=(',',':'),allow_nan=False)
            stream.write('\n')
    else:
        expected=read_json(args.verify)
        validate(expected)
        compare(result,expected)
    # A failed confirmatory gate leaves its complete receipt at --record PATH.
    validate(result)
    print(json.dumps({'verified':True,'construction':result['construction'],
                      'program_edges':1336,'cases':6}))


if __name__=='__main__':
    main()
