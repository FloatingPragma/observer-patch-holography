"""Changed-clock control: literal arc words, all counters, unchanged source frame."""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback
import numpy as np
import source_join
from custody import check

HERE=Path(__file__).resolve().parent
OPENWAVE=HERE.parents[2]
AMP=OPENWAVE/'row_07_double_cover/work/2026-11-16_codex_stationary_readout_repair_v2'
sys.path.insert(0,str(AMP))
import protocol as p
from interval_reader import read
from interval_math import real,packed,endpoints

CASES=[]


def require(ok,message):
    if not ok:raise ValueError(message)


def normal(x):
    if isinstance(x,np.ndarray):return {'real':x.real.tolist(),'imag':x.imag.tolist()} if np.iscomplexobj(x) else x.tolist()
    if isinstance(x,complex):return {'real':x.real,'imag':x.imag}
    if isinstance(x,np.generic):return normal(x.item())
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):normal(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [normal(v) for v in x]
    return x


def run(ks):
    check();started=time.monotonic();provider_before=dict(p.old.CALLS);CASES.clear();join=source_join.run();frames={};prepared={};prep_records={};clock_records={}
    sources=join['representatives']
    for sign in (-1,1):
        native=sources[f'native:{sign}'];changed=sources[f'swap:{sign}']
        require(native['prep_columns']==changed['prep_columns'],'same native preparation columns')
        factor,success,record=p.prepare('native',native['prep_columns'],sign,native['prep_blank'])
        prepared[sign]=(factor,success);prep_records[str(sign)]=record
        frames[sign]=p.catalog(32,'native',columns=native['prep_columns'],sign=sign)
    for k in ks:
        for sign in (-1,1):
            for clock_law in ('native','swap'):
                measured=sources[f'{clock_law}:{sign}'];cl=p.clock(k,measured['clock_columns'],sign)
                R=np.eye(33)[measured['spectrum']['R_permutation']]
                require(np.linalg.norm(cl['R']-R)<1e-14,'actual measured R was not used')
                means=list(map(F,measured['spectrum']['mode_means']))
                mu=1+(np.asarray(list(map(float,means)))-1)/(4**k)
                require(np.max(abs(np.linalg.eigvalsh(cl['P'])-np.sort(np.r_[mu,1.])))<1e-12,'derived A*RA source spectrum')
                factor,decoded,entry=p.counter_input(prepared[sign][0],cl)
                clock_records[f'{clock_law}:{sign}:{k}']={'input':entry,'clock':{name:cl[name] for name in ('R','A','P','W','a','decode','Q','epsilon','columns','edges','charge_sign','blank_predicate')},
                    'measured_mode_means':measured['spectrum']['mode_means'],'word_iteration_count':cl['Q']}
                for turns in (1,2):
                    for direction in (-1,1):
                        result=p.run_case(factor,cl,32,2,turns,direction,'native',frames[sign],prepared[sign][1]*decoded,1)
                        result.update(clock_law=clock_law,frame_law='native',preparation_law='native')
                        CASES.append(result)
                        require(result['visibility']<-.98 if turns==1 else result['visibility']>.96,'same native frame sign/return gate')
                        unconditional=result['selected_clock_probability_given_entry']*float(1-F(1816,2**64))
                        require(unconditional>.03,'complete double-calibration acceptance gate')
                        result['selected_clock_probability_unconditional_lower']=unconditional
                        selected=1 if turns%2 else 0;branch=result['returned']['conditional_clocks'][selected]
                        rho=np.asarray(branch['arc_density']['real'])+1j*np.asarray(branch['arc_density']['imag'])
                        vectors=[];zs=[]
                        for mode,mean in enumerate(means):
                            cosine=1+(float(mean)-1)/(4**k);sine=np.sqrt(1-cosine*cosine)
                            a=cl['A']@p.OPTICAL[:,mode];b=(cl['R']@a-cosine*a)/sine
                            vectors.extend([(a-1j*b)/np.sqrt(2),(a+1j*b)/np.sqrt(2)]);zs.extend([cosine+1j*sine,cosine-1j*sine])
                        C=np.column_stack(vectors);spectral=C.conj().T@rho@C;weights=np.diag(spectral).real
                        require(np.linalg.norm(rho-C@spectral@C.conj().T)<1e-9,'full retained four-mode clock state')
                        off=1-weights[0]**2-weights[1]**2
                        require(min(weights)>0 and off<.01,'actual changed-clock pair impurity gate')
                        records={}
                        for stride in (1,8):
                            observed=branch['ramsey_stride'+str(stride)][1]
                            moment=complex(2*observed['x_plus']-1,2*observed['y_plus']-1)**2
                            records[str(stride)]={'x_plus':(1+moment.real)/2,'x_minus':(1-moment.real)/2,
                                'y_plus':(1+moment.imag)/2,'y_minus':(1-moment.imag)/2}
                        # Reader receives records and its fixed certificates only.
                        bands=read(*[{name:real(str(value)) for name,value in records[str(stride)].items()} for stride in (1,8)])
                        require(len(bands)==1,'source-blind all-branch uniqueness')
                        prediction=2*2**k*np.arccos(1+(float(means[0])-1)/(4**k))
                        low,high=map(float,endpoints(bands[0]))
                        require(low<=prediction<=high,'changed-clock prediction after readout')
                        # The native fixed anchor is unchanged for this control.
                        native_prediction=2*2**k*np.arccos(1-(1-float(F(3,5)))/(4**k))
                        if clock_law=='swap':require(not low<=native_prediction<=high,'changed-clock detector did not reject native prediction')
                        result['independent_pair']={'signed_mode_weights':weights.tolist(),'joint_mode_weights':np.outer(weights,weights).tolist(),
                            'pair_off_band':float(off),'records':records,'frequency_interval':packed(bands[0]),
                            'energy_interval_keV':packed(bands[0]*(511*real(5)**.5/2)),
                            'source_prediction_after_readout':float(prediction),'native_prediction_after_readout':float(native_prediction),
                            'same_511_anchor':True,'reader_given_expected_frequency':False}
        print(json.dumps({'k':k,'cases':len(CASES)}),flush=True)
    require(len(CASES)==16*len(ks),'complete changed-clock case coverage')
    for row in CASES:
        if row['clock_law']!='swap':continue
        baseline=next(x for x in CASES if x['clock_law']=='native' and all(x[key]==row[key] for key in ('charge_sign','k','turns','direction')))
        a=list(map(float,endpoints(real_from_pack(row['independent_pair']['frequency_interval']))));b=list(map(float,endpoints(real_from_pack(baseline['independent_pair']['frequency_interval']))))
        require(a[0]>b[1] or b[0]>a[1],'same-reader native/control intervals overlap')
    return normal({'candidate':'A-native-changed-clock-control-v1','passed':True,'k_values':ks,'case_count':len(CASES),
        'source_join':join,'native_preparation':prep_records,'native_frames':frames,'clocks':clock_records,'cases':CASES,
        'seconds':time.monotonic()-started,'literal_provider_endpoint_calls':{k:v-provider_before.get(k,0) for k,v in p.old.CALLS.items()},'word_updates_are_literal_native_calls':False,
        'readout_scope':'Exact represented source record laws; no finite-shot precision claim','external_G6':False})


def real_from_pack(x):
    from interval_math import iv
    def dyadic(r):return (-1 if r[0] else 1)*real(int(r[1]))*real(2)**r[2]
    return iv.mpf([dyadic(x['lo_dyadic']).a,dyadic(x['hi_dyadic']).b])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--k',nargs='+',type=int,default=[0]);group=parser.add_mutually_exclusive_group();group.add_argument('--record-to',type=Path);group.add_argument('--verify',type=Path);args=parser.parse_args()
    require(len(args.k)==len(set(args.k)) and set(args.k)<={0,6},'filed k menu')
    if args.record_to is not None:require(not args.record_to.exists(),'refusing evidence overwrite')
    try:
        result=run(args.k)
        if args.record_to is not None:
            with args.record_to.open('xb') as f:f.write(gzip.compress((json.dumps(result,separators=(',',':'))+'\n').encode(),mtime=0))
        else:
            old=json.loads(gzip.decompress((args.verify or HERE/'HEADLINE.json.gz').read_bytes()));a=dict(result);b=dict(old);a.pop('seconds');b.pop('seconds');require(a==b,'saved changed-clock evidence differs')
        print(json.dumps({k:v for k,v in result.items() if k not in ('source_join','native_preparation','native_frames','clocks','cases')},indent=2))
    except Exception as error:
        if args.record_to is not None:
            failure=args.record_to.with_name(args.record_to.name+'.FAILURE.json.gz')
            with failure.open('xb') as f:f.write(gzip.compress(json.dumps(normal({'error':repr(error),'traceback':traceback.format_exc(),'source_progress':source_join.PROGRESS,'completed_or_current_cases':CASES})).encode(),mtime=0))
        raise


if __name__=='__main__':main()
