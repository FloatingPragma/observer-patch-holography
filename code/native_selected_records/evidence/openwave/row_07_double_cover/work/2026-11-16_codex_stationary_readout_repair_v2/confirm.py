"""Default: read-only replay. First evidence requires exclusive --first-record.

The finite coefficient confirmation reproduces disclosed exploratory leads.
It does not award external review or replace the separate interval certificate.
"""
import argparse
import gzip
import hashlib
import json
import math
import time
from fractions import Fraction
from pathlib import Path
import numpy as np
import native_source as old
import protocol as p
import workspaces
import error_budget
from amplified_native import require

ROOT=Path(__file__).resolve().parent
FREEZE_FILE='PREREG_V5.json'


def encode(value):
    if isinstance(value,np.ndarray):return old.encode(value)
    if isinstance(value,complex):return [value.real,value.imag]
    if isinstance(value,Fraction):return str(value)
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {k:encode(v)for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [encode(v)for v in value]
    return value


def compare(a,b,path='$'):
    if isinstance(a,dict):
        require(a.keys()==b.keys(),path+' keys differ')
        for key in a:compare(a[key],b[key],path+'.'+key)
    elif isinstance(a,list):
        require(len(a)==len(b),path+' length differs')
        for i,(x,y)in enumerate(zip(a,b)):compare(x,y,path+'['+str(i)+']')
    elif isinstance(a,int)and not isinstance(a,bool):
        require(type(b)is int and a==b,path+' exact integer/type mismatch')
    elif isinstance(a,float):
        require(type(b)is float,path+' float/type mismatch')
        require(math.isfinite(a)and math.isfinite(b),path+' nonfinite number')
        if path.endswith('.selected_clock_probability_given_entry')or path.endswith('.selected_clock_probability_including_calibration_given_entry'):
            require((a==0)==(b==0)and(a==b or abs(a-b)<=1e-9*max(abs(a),abs(b))),path+' zero-preserving relative cost mismatch')
        else:require(abs(a-b)<=2e-10+2e-12*max(abs(a),abs(b)),path+' numerical difference')
    else:require(a==b,path+' value differs')


def exact_path_semantics(record):
    control=record.get('exact_scalar_control')
    if control is None:return
    scalars=tuple(Fraction(x)for x in control['scalars'])
    _,_,plus,minus=p.scalar_path(scalars,record['M'],record['h'],record['turns'],record['direction'],record['filter_stride'])
    require(control['minus_exact_zero']==(minus==0),'Exact zero/positive path flag changed')
    for index,exact in enumerate((plus,minus)):
        observed=record['returned']['path_probabilities_plus_minus'][index]
        require((observed>0)==(exact>0),'A positive exact branch was lost or a zero branch invented')
        if exact:
            require(abs(observed/float(exact)-1)<1e-10,'Tiny exact branch relative probability changed')
        require(record['returned']['conditional_clocks'][index]['available']==(exact>0),'Conditional clock availability contradicts exact path law')
    selected=1 if record['turns']%2 else 0
    require(record['selected_clock_available']==((plus,minus)[selected]>0),'Fixed-reader continuation availability mismatch')
    require((record['selected_clock_probability_given_entry']>0)==((plus,minus)[selected]>0),'Rare selected-clock cost was rounded to zero')


def verify_freeze():
    freeze=json.loads((ROOT/FREEZE_FILE).read_text())
    for name,expected in freeze['files'].items():
        require(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,'Frozen input/source changed: '+name)


def verify_custody(folder):
    index=json.loads((folder/'INDEX.json').read_text())
    names=[r['path']for r in index['files']]
    require(len(names)==len(set(names)),'Duplicate saved artifact path')
    for row in index['files']:
        require(Path(row['path']).name==row['path'],'Artifact path escapes confirmation folder')
        require(hashlib.sha256((folder/row['path']).read_bytes()).hexdigest()==row['sha256'],'Saved artifact digest mismatch: '+row['path'])


def emit(folder,name,value,first):
    path=folder/(name+'.json.gz');value=encode(value)
    data=(json.dumps(value,separators=(',',':'),allow_nan=False)+'\n').encode()
    if first:
        with path.open('xb')as raw:
            with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0)as stream:stream.write(data)
    else:
        with gzip.open(path,'rt')as stream:saved=json.load(stream)
        if isinstance(saved,dict)and'exact_scalar_control'in saved:exact_path_semantics(saved)
        compare(value,saved)
    return {'path':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def scientific_gates(preparations,rows):
    for group in preparations.values():
        require(group['native']['optical_weights'][0]>.999999,'Native selected-mode gate')
        require(group['swap']['optical_weights'][0]<.7,'Changed-law preparation control did not fail selection')
    expected={(sign,k,stride,turns,direction,law)for sign in(-1,1)for k in range(7)
        for stride in(1,2)for turns in(1,2)for direction in(-1,1)
        for law in(('native','swap','swap_fixed_native_calibration')if k in(0,6)else('native','swap'))}
    key=lambda r:(r['charge_sign'],r['k'],r['filter_stride'],r['turns'],r['direction'],r['law'])
    require(len(rows)==len(expected)and len({key(r)for r in rows})==len(expected)and {key(r)for r in rows}==expected,
        'Confirmation has missing, duplicate or extra scientific cases')
    lookup={key(r):r for r in rows}
    for row in rows:
        if row['law']=='native':
            threshold=.96 if row['filter_stride']==1 else .95
            require(row['visibility']<-.98 if row['turns']==1 else row['visibility']>threshold,'Native stationary sign/return gate')
            require(row['selected_clock_probability_including_calibration_given_entry']>.03,'Native complete selected-clock acceptance gate')
            if row['filter_stride']==1:
                coarse=lookup[row['charge_sign'],row['k'],2,row['turns'],row['direction'],'native']
                require(abs(row['visibility'])>abs(coarse['visibility'])+.001,'Finer applied frames did not improve the visibility')
        elif row['law']=='swap':
            require(row['visibility']>.999999,'Adaptive swap must lack a one-turn sign')
            if row['turns']==1:
                require(row['selected_clock_probability_given_entry']<1e-15,'Same fixed parity reader failed to reject the swap sign control')
        else:
            require(row['visibility']>=0,'Fixed native-calibration swap falsely produces negative sign')
        require(row['geometry']['reference_native_hops']==0,'Reference moved')
        require(row['geometry']['maximum_increment_degrees']<10,'Frames too coarse')
        require(row['geometry']['maximum_applied_increment_degrees']<15,'Applied frames too coarse')


def run(folder,first):
    artifacts=[];preparations={};factors={};success={}
    artifacts.append(emit(folder,'WORKSPACES',workspaces.audit(),first))
    calibration_path=ROOT/'CALIBRATION_INPUT.json'
    calibrations=json.loads(calibration_path.read_text())
    entries=calibrations['entries']
    require(len(entries)==4 and {(e['sign'],e['law'])for e in entries}=={(s,l)for s in(-1,1)for l in('native','swap')},
        'Operational calibration input coverage')
    tables={(e['sign'],e['law']):e['columns']for e in entries}
    blanks={(e['sign'],e['law']):e['actual_blank_records']for e in entries}
    artifacts.append(emit(folder,'CALIBRATION_INPUT',calibrations,first))
    for sign in(-1,1):
        sk=str(sign);factors[sk]={};success[sk]={};preparations[sk]={}
        for law in('native','swap'):
            factors[sk][law],success[sk][law],preparations[sk][law]=p.prepare(law,tables[sign,law],sign,blanks[sign,law])
    artifacts.append(emit(folder,'PREPARATION',preparations,first))
    catalogs={}
    for sign in(-1,1):
        sk=str(sign);catalogs[sk]={'native':p.catalog(32,'native',columns=tables[sign,'native'],sign=sign)}
        for law in('swap','swap_fixed_native_calibration'):
            catalogs[sk][law]=p.catalog(32,law,catalogs[sk]['native'],tables[sign,'swap'],sign)
    artifacts.append(emit(folder,'CATALOG',catalogs,first))
    summary=[]
    for sign in(-1,1):
        sk=str(sign)
        for k in range(7):
            cl=p.clock(k,tables[sign,'native'],sign)
            factor,decoded,entry=p.counter_input(factors[sk]['native'],cl)
            artifacts.append(emit(folder,'INPUT_s'+sk+'_k'+str(k),{'source':entry,'clock':{
                key:cl[key]for key in ('R','A','P','W','a','decode','Q','epsilon','columns','edges','charge_sign','blank_predicate')}},first))
            for stride,h in((1,2),(2,4)):
                for turns in(1,2):
                    for direction in(-1,1):
                        laws=('native','swap','swap_fixed_native_calibration')if k in(0,6)else('native','swap')
                        for law in laws:
                            row=p.run_case(factor,cl,32,h,turns,direction,law,catalogs[sk][law],success[sk]['native']*decoded,stride)
                            count=1 if law=='native'else 2
                            calibration_lower=1-count*Fraction(calibrations['failure_probability_upper_per_census'])
                            row['source_calibration_success_lower']=str(calibration_lower)
                            row['selected_clock_probability_including_calibration_given_entry']=row['selected_clock_probability_given_entry']*float(calibration_lower)
                            row['calibration_inclusive_cost_is_lower_bound']=True
                            row['source_calibration_scope']='Bounded continuing-unit record census; actual source tapes retained in the bound provider packet.'
                            row['service_budget']=error_budget.budget(k,turns,count)
                            name='CASE_s'+sk+'_k'+str(k)+'_r'+str(stride)+'_t'+str(turns)+'_d'+str(direction)+'_'+law
                            artifacts.append(emit(folder,name,row,first))
                            summary.append({key:row[key]for key in ('charge_sign','k','M','h','filter_stride','turns','direction','law','visibility',
                                'selected_path','selected_clock_probability_given_entry','selected_clock_available',
                                'selected_clock_probability_including_calibration_given_entry',
                                'geometry','source_accesses_rotation_including_padding','source_accesses_counter')})
    scientific_gates(preparations,summary)
    result={'candidate':'A-native-finite-amplified-stationary-source-v5','rows':summary,
        'preparation_optical_weights':{s:{law:rec['optical_weights']for law,rec in group.items()}for s,group in preparations.items()},
        'preparation_success':success,'gates_pass':True,'prior_pilots_disclosed':True,
        'reference':'literal stationary; full path law and fixed parity continuation',
        'finite_coefficients':'40-bit nearest rational tangent; finite inactive-selector amplitudes retained',
        'numerical_scope':'Nominal complete-source replay. Interval certification and the explicit implementation envelope are separately required.',
        'classical_entry_scope':'Probabilities condition on the separately guarded common Q1 entry and declared blank/vacuum registers; no old unconditional prefix is invented.',
        'memory_boundary':'Declared coherent source inverse V†, phase/Hadamard/address/record services; 2forward+1inverse source access perR.',
        'literal_native_calls':dict(old.CALLS)}
    artifacts.append(emit(folder,'ERROR_BUDGET',{'spin':error_budget.envelope(),
        'joint_pair_stricter_example':error_budget.envelope(52,2),'sign_assay':error_budget.sign_assay()},first))
    artifacts.append(emit(folder,'SUMMARY',result,first))
    return result,artifacts


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--first-record',action='store_true')
    args=parser.parse_args();folder=ROOT/'CONFIRMATION'
    verify_freeze()
    if args.first_record:folder.mkdir(exist_ok=False)
    else:
        require(folder.is_dir(),'Missing saved confirmation; explicit exclusive first-record required')
        verify_custody(folder)
    start=time.monotonic();result,artifacts=run(folder,args.first_record)
    index={'candidate':result['candidate'],'files':artifacts}
    path=folder/'INDEX.json'
    if args.first_record:
        with path.open('x')as stream:json.dump(index,stream,indent=2);stream.write('\n')
    else:compare(index,json.loads(path.read_text()))
    print(json.dumps({'verified':not args.first_record,'created':args.first_record,
        'cases':len(result['rows']),'gates_pass':result['gates_pass'],'seconds':time.monotonic()-start,
        'minimum_native_success_including_calibration':min(r['selected_clock_probability_including_calibration_given_entry']for r in result['rows']if r['law']=='native'),
        'index_sha256':hashlib.sha256(path.read_bytes()).hexdigest()},indent=2))


if __name__=='__main__':main()
