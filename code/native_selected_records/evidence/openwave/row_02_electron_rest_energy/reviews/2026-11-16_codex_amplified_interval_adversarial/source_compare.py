"""Compare exact-source certificate with full stored producer densities and tapes.

This independently reads source artifacts, imports no producer code and does
not claim that floating source centers are exact interval endpoints.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time
import numpy as np

HERE=Path(__file__).resolve().parent
OPENWAVE=HERE.parents[2]
SOURCE=OPENWAVE/'row_07_double_cover/work/2026-11-16_codex_stationary_readout_repair_v2/CONFIRMATION'
CERT=HERE.parents[1]/'work/2026-11-16_codex_amplified_source_interval_v1'
TOLERANCE=1e-8


def require(ok,message):
    if not ok:raise ValueError(message)


def matrix(x):return np.array(x['real'])+1j*np.array(x['imag'])
def point(x):return sum(map(float,x['display']))/2
def read(path):return json.loads(gzip.decompress(path.read_bytes()))


def run():
    started=time.monotonic();index=json.loads((SOURCE/'INDEX.json').read_text());listed={r['path']:r['sha256'] for r in index['files']}
    require(len(listed)==len(index['files']),'source index duplicate')
    checked=set()
    def source_file(name):
        p=SOURCE/name;blob=p.read_bytes();require(hashlib.sha256(blob).hexdigest()==listed[name],'source file custody '+name);checked.add(name)
        return json.loads(gzip.decompress(blob))
    certificate=json.loads((CERT/'CERTIFICATE.json').read_text());inputs=json.loads((CERT/'SOURCE_INPUTS.json').read_text())
    prep=source_file('PREPARATION.json.gz');catalog=source_file('CATALOG.json.gz')
    for sign,program in inputs['sources'].items():
        for label,descriptor in program['preparation'].items():
            require(all(prep[sign]['native'][label][k]==v for k,v in descriptor.items()),'actual preparation descriptor binding')
        for law,frames in program['catalogs'].items():
            for i,descriptor in enumerate(frames):require(all(catalog[sign][law][i][k]==v for k,v in descriptor.items()),'actual full frame descriptor binding')
    E=np.array([[0,-2/np.sqrt(6)],[1/np.sqrt(2),1/np.sqrt(6)],[-1/np.sqrt(2),1/np.sqrt(6)]])
    clock={};maxima={k:0. for k in ['visibility','selected_cost','signed_weight','optical_subspace','ramsey_pair_probability','first_failure_mass','preparation_success','counter_decode_success','clock_eigenblock']}
    cases=[];target_stages=marker_hops=ramsey_records=0
    for row in certificate['cases']:
        sign,k=row['charge_sign'],row['k'];key=(sign,k)
        if key not in clock:
            entry=source_file(f'INPUT_s{sign}_k{k}.json.gz');cl=entry['clock'];A,R,W=map(matrix,(cl['A'],cl['R'],cl['W']))
            vectors=[];zs=[]
            for mode,lam in enumerate((4,6)):
                p=1-lam/(10*4**k);t=np.sqrt(1-p*p);a=A@E[:,mode];b=(R@a-p*a)/t
                vectors.extend([(a-1j*b)/np.sqrt(2),(a+1j*b)/np.sqrt(2)]);zs.extend([p+1j*t,p-1j*t])
            C=np.column_stack(vectors);zs=np.array(zs)
            maxima['clock_eigenblock']=max(maxima['clock_eigenblock'],float(np.linalg.norm(W@C-C*zs)))
            require(np.linalg.norm(C.conj().T@C-np.eye(4))<TOLERANCE,'source chart orthonormality')
            clock[key]=(C,zs,entry)
            saved=next(r for r in certificate['preparations'] if (r['charge_sign'],r['k'])==key)
            maxima['preparation_success']=max(maxima['preparation_success'],abs(prep[str(sign)]['native']['success']-point(saved['success'])))
            actual=entry['source'];value=actual['counter_success']*actual['decode_success']
            maxima['counter_decode_success']=max(maxima['counter_decode_success'],abs(value-point(saved['counter_decode_success_given_preparation'])))
            for probabilities in [actual['counter_probabilities'],actual['decode_probabilities']]:require(min(probabilities)>=0 and abs(sum(probabilities)-1)<TOLERANCE,'complete counter/decode probabilities')
        C,zs,entry=clock[key];stride,turns,direction,law=(row[n] for n in ('filter_stride','turns','direction','law'))
        name=f'CASE_s{sign}_k{k}_r{stride}_t{turns}_d{direction}_{law}.json.gz';actual=source_file(name)
        require(all(actual[n]==row[n] for n in ('k','charge_sign','law','M','h','turns','direction','filter_stride','selected_path')),'complete source case setting')
        maxima['visibility']=max(maxima['visibility'],abs(actual['visibility']-point(row['visibility'])))
        maxima['selected_cost']=max(maxima['selected_cost'],abs(actual['selected_clock_probability_given_entry']-point(row['selected_clock_probability_given_classical_entry'])))
        require(len(actual['marker_records'])==96*turns,'actual marker coverage')
        fail=0.
        for stage in actual['marker_records']:
            require(sum(h['probability'] for h in stage['histories'])>1-TOLERANCE,'complete marker histories')
            fail+=stage['first_failure_mass'];marker_hops+=1
        expected_stages=2*row['h']*(96*turns//stride)
        got_stages=0
        for stage in actual['target_first_failures']:
            if stage['kind']=='declared_filter_hold':continue
            probs=stage['probabilities'];require(min(probs)>=0 and abs(sum(probs)-1)<TOLERANCE,'complete target flag law')
            require(abs(stage['first_failure_mass']-stage['reach_before']*sum(probs[1:]))<TOLERANCE,'target first-failure mass')
            fail+=stage['first_failure_mass'];target_stages+=1;got_stages+=1
        require(got_stages==expected_stages,'target stage multiplicity')
        maxima['first_failure_mass']=max(maxima['first_failure_mass'],abs(fail+actual['raw_float_rotation_marker_success']-1))
        selected=1 if turns%2 else 0;branch=actual['returned']['conditional_clocks'][selected];rho=matrix(branch['arc_density'])
        require(branch['available'] and branch['path']==row['selected_path'],'actual parity branch')
        require(np.linalg.norm(rho-rho.conj().T)<TOLERANCE and abs(np.trace(rho)-1)<TOLERANCE,'actual normalized Hermitian density')
        spectral=C.conj().T@rho@C;weights=np.diag(spectral).real
        maxima['signed_weight']=max(maxima['signed_weight'],max(abs(weights[i]-point(row['signed_mode_weights'][i])) for i in range(4)))
        maxima['optical_subspace']=max(maxima['optical_subspace'],float(np.linalg.norm(rho-C@spectral@C.conj().T)))
        for read_stride in (1,8):
            observed=branch['ramsey_stride'+str(read_stride)];require([r['time'] for r in observed]==list(range(17)),'all stored Ramsey times')
            for r in observed:
                moment=sum(weights*zs**(read_stride*2**k*r['time']))
                require(max(abs(r['x_plus']-(1+moment.real)/2),abs(r['y_plus']-(1+moment.imag)/2))<TOLERANCE,'actual source full Ramsey law')
                ramsey_records+=1
            first=observed[1];moment=complex(2*first['x_plus']-1,2*first['y_plus']-1)**2
            result={'x_plus':(1+moment.real)/2,'x_minus':(1-moment.real)/2,'y_plus':(1+moment.imag)/2,'y_minus':(1-moment.imag)/2}
            maxima['ramsey_pair_probability']=max(maxima['ramsey_pair_probability'],max(abs(v-point(row['pair_records'][str(read_stride)][n])) for n,v in result.items()))
        cases.append(name)
    require(len(cases)==len(set(cases))==256,'all actual source case multiplicities')
    require(all(x<TOLERANCE for x in maxima.values()),'full-source residual threshold '+str(maxima))
    return {'passed':True,'producer_imports':False,'source_index_sha256':hashlib.sha256((SOURCE/'INDEX.json').read_bytes()).hexdigest(),
        'certificate_sha256':hashlib.sha256((CERT/'CERTIFICATE.json').read_bytes()).hexdigest(),'case_count':len(cases),
        'checked_source_files':len(checked),'target_flag_stages':target_stages,'marker_hops':marker_hops,'stored_ramsey_records':ramsey_records,
        'max_absolute_residuals':maxima,'predeclared_float_comparison_tolerance':TOLERANCE,'seconds':time.monotonic()-started,
        'scope':'Stored full actual source matrices/flags and positive-counter densities agree with the exact finite-source certificate. Floating source centers are not claimed to be exact dyadic interval points. No laboratory/hardware realization or external G6.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--record-to',type=Path);args=p.parse_args()
    if args.record_to is not None and args.record_to.exists():raise FileExistsError(args.record_to)
    result=run()
    if args.record_to is not None:
        with args.record_to.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    else:
        saved=json.loads((HERE/'SOURCE_COMPARISON.json').read_text());a=dict(result);b=dict(saved);a.pop('seconds');b.pop('seconds');require(a==b,'saved comparison differs')
    print(json.dumps(result,indent=2))
