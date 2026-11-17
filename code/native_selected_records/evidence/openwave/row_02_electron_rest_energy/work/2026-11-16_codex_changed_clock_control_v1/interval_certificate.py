"""Outward changed-clock observable, with means derived from measured columns."""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import time
from interval_math import *
from interval_reader import read
from filter_math import preparation,arms,returned,pair_records
from counter_core import counter
from custody import check

HERE=Path(__file__).resolve().parent


def rational(x):
    v=F(x);return real(v.numerator)/v.denominator


def mode_means(columns):
    require(len(columns)==18,'measured conditional clock domain')
    tab={(r['input_position'],tuple(r['edge']),r['coin']):r for r in columns}
    require(len(tab)==18,'unique measured columns')
    tickets=((0,1),(0,1),(0,2),(0,2),(1,2));P=[[F(0) for _ in range(3)] for _ in range(3)]
    R=[]
    for p in range(3):
        for h in range(11):
            if h==10:R.append(11*p+10);continue
            row=tab[p,tickets[h//2],h%2];P[row['decoded_output']][p]+=F(1,10)
            R.append(11*row['decoded_output']+(h^int(row['current']!=0)))
    require(sorted(R)==list(range(33)) and all(R[R[i]]==i for i in range(33)),'actual source reversal')
    means=[]
    for v in ((0,1,-1),(-2,1,1)):
        y=[sum(P[i][j]*v[j] for j in range(3)) for i in range(3)];j=next(i for i,x in enumerate(v) if x);mu=y[j]/v[j]
        require(y==[mu*x for x in v],'measured optical clock eigenmode');means.append(mu)
    for k in (0,6):
        weights=[F(1,10*4**k)]*10+[1-F(1,4**k)];AstarRA=[[F(0) for _ in range(3)] for _ in range(3)]
        for p in range(3):
            for h in range(11):
                y,hh=divmod(R[11*p+h],11);require(weights[h]==weights[hh],'actual coherent choice amplitude');AstarRA[y][p]+=weights[h]
        require(AstarRA==[[P[i][j]/(4**k)+(1-F(1,4**k))*int(i==j) for j in range(3)] for i in range(3)],'native R/A source transfer')
    return means


def run(source_path):
    check();started=time.monotonic();data=json.loads(gzip.decompress(source_path.read_bytes()));require(data['passed'],'primary source failed')
    entries=data['source_join']['representatives'];programs=data['native_preparation'];frames=data['native_frames'];rows=[];counter_cache={}
    epsilon=F(2*2*2**22+1,2**52);conditional=epsilon/F(9,10000)
    require(conditional<F('0.00002'),'same full-service pair budget')
    for k in data['k_values']:
        for sign in (-1,1):
            for clock_law in ('native','swap'):
                means=mode_means(entries[f'{clock_law}:{sign}']['clock_columns']);key=(k,tuple(means))
                if key not in counter_cache:counter_cache[key]=counter(k,means)
                gram,z=counter_cache[key]
                position,prep_probability,prep_weights=preparation(programs[str(sign)],gram)
                for turns in (1,2):
                    for direction in (-1,1):
                        pair=arms(frames[str(sign)],32,2,turns,direction,1,False)
                        visibility,selected,weights,analyzer=returned(gram,position,pair,turns,False)
                        cost=selected*(1-real(1)/2**16)**(96*turns)*(1-(real(3)/8)**8)**2
                        unconditional=cost*(1-rational(F(1816,2**64)))
                        require(unconditional>real('.03'),'complete changed-clock probability floor')
                        require(visibility<real('-.98') if turns==1 else visibility>real('.96'),'unchanged frame sign/return')
                        joint,off,records=pair_records(weights,weights,z,k);require(off<real('.01'),'changed-clock pair impurity')
                        bands=read(records['1'],records['8']);require(len(bands)==1,'source-blind changed-clock inverse branches')
                        prediction=2*2**k*acos_point(1-(1-rational(means[0]))/(4**k))
                        require(endpoints(bands[0])[0]<=endpoints(prediction)[0] and endpoints(bands[0])[1]>=endpoints(prediction)[1],'frequency comparison after reading')
                        native_mu=mode_means(entries[f'native:{sign}']['clock_columns'])[0]
                        native_prediction=2*2**k*acos_point(1-(1-rational(native_mu))/(4**k))
                        if clock_law=='swap':require(endpoints(bands[0])[0]>endpoints(native_prediction)[1] or endpoints(bands[0])[1]<endpoints(native_prediction)[0],'changed-clock readout failed to exclude native frequency')
                        original=next(r for r in data['cases'] if (r['clock_law'],r['charge_sign'],r['k'],r['turns'],r['direction'])==(clock_law,sign,k,turns,direction))
                        def discrepancy(ivalue,number):return max(abs(float(v)-number) for v in endpoints(ivalue))
                        differences=[discrepancy(visibility,original['visibility']),discrepancy(cost,original['selected_clock_probability_given_entry'])]
                        for stride in ('1','8'):
                            for name,value in records[stride].items():differences.append(discrepancy(value,original['independent_pair']['records'][stride][name]))
                        require(max(differences)<1e-8,'full literal-word source and outward observable differ')
                        rows.append({'clock_law':clock_law,'charge_sign':sign,'k':k,'turns':turns,'direction':direction,
                            'mode_means_from_records':[str(x) for x in means],'visibility':packed(visibility),
                            'selected_probability_given_entry':packed(cost),'unconditional_probability_lower':packed(unconditional),
                            'signed_mode_weights':[packed(x) for x in weights],'all16_joint_weights':[[packed(x) for x in r] for r in joint],
                            'pair_off_band':packed(off),'records':{s:{n:packed(x) for n,x in rec.items()} for s,rec in records.items()},
                            'frequency_interval':packed(bands[0]),'energy_interval_keV':packed(bands[0]*(511*iv.sqrt(real(5))/2)),
                            'derived_prediction_after_reader':packed(prediction),'native_prediction_after_reader':packed(native_prediction),
                            'maximum_primary_float_discrepancy':max(differences),'same_anchor':True,
                            'conditional_pair_service_error':str(conditional),'reader_binary_allowance':'0.00002',
                            'scope':'Actual clock-law control; exact record-law interval and declared complete-service precision, not finite-shot mass measurement.'})
    require(len(rows)==16*len(data['k_values']),'complete interval case domain')
    return {'candidate':'A-native-changed-clock-control-v1','passed':True,'rows':rows,'case_count':len(rows),
        'source_sha256':hashlib.sha256(source_path.read_bytes()).hexdigest(),'reader_receives_expected_frequency':False,
        'per_source_entry_failure_upper':str(F(1816,2**64)),'service_cap_each':2**22,'pair_service_delta':str(F(1,2**52)),
        'counter_scope':'Complete finite geometric source word from measured R/A means, after exact signed source binding; full bins are in the primary FFT records.',
        'seconds':time.monotonic()-started,'external_G6':False}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,default=HERE/'EVIDENCE.json.gz');p.add_argument('--record-to',type=Path);args=p.parse_args()
    if args.record_to is not None:require(not args.record_to.exists(),'refusing interval overwrite')
    try:result=run(args.source)
    except Exception as error:
        if args.record_to is not None:
            import traceback
            with args.record_to.with_name(args.record_to.name+'.FAILURE.json').open('x') as f:json.dump({'error':repr(error),'traceback':traceback.format_exc()},f,indent=2)
        raise
    if args.record_to is not None:
        with args.record_to.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    else:
        saved=json.loads((HERE/'INTERVALS.json').read_text());a=dict(result);b=dict(saved);a.pop('seconds');b.pop('seconds');require(a==b,'saved interval differs')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
