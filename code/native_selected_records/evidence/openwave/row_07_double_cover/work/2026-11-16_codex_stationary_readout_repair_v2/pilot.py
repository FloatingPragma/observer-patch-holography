"""Exploration only: complete primitive instruments, actual old mixed input.
Output is created only with --first-record; default verifies the saved pilot.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import native_source as source
from amplified_native import compile_instrument,routing,OPTICAL,require

ROOT = Path(__file__).resolve().parent
PERMUTATION = np.eye(3)[[1,2,0]]


def input_state(k):
    cl = source.clock(k)
    rho = np.diag([0.,1.,0.])
    logp = 0.
    stages = [(source.contrast(cl['edges'],'native'),[1])]+[(source.instrument(cl['edges'],'native'),[0])]*8
    for ops,keep in stages+[(cl['qpe'],range(1,cl['Q']//2)),(cl['decode'],[0])]:
        rho,p,probabilities = source.channel(rho,ops,keep)
        source.check_probabilities(probabilities)
        logp += math.log10(p)
    # Exact source constant-mode removal fixes this invariant sector. This is
    # the full mixed optical density, not a chosen ideal ray. Record raw drift.
    optical = OPTICAL.T@rho@OPTICAL
    drift = float(np.trace(optical).real-1)
    require(abs(drift) < 1e-10, 'Input optical trace drift')
    optical /= np.trace(optical).real
    reencoder = np.einsum('kxa,ab->kxb',cl['qpe'][1:cl['Q']//2],OPTICAL)
    metric = np.einsum('kxa,kxb->ab',reencoder.conj(),reencoder)
    return optical,metric,logp,drift


def run():
    rows = []
    calibration = []
    prepared = {k:input_state(k) for k in (0,6)}
    for m in (4,8,16,32):
        native_catalog = [compile_instrument(routing(slot,m),'native') for slot in range(3*m)]
        for law in ('native','swap','swap_fixed_native_calibration'):
            if law == 'native':
                catalog = native_catalog
            else:
                catalog = [compile_instrument(routing(slot,m),'swap',
                    native_catalog[slot] if law=='swap_fixed_native_calibration' else None) for slot in range(3*m)]
            basis = PERMUTATION@OPTICAL
            filters = [basis.T@c['K']@basis for c in catalog]
            calibration += [{'M':m,'law':law,'slot':j,'s':c['s'],'phase':c['phase'],
                'preboost':c['preboost'],'norm_record':c['norm_record'],'native_reversals':c['native_reversals_per_filter']}for j,c in enumerate(catalog)]
            for k in (0,6):
                rho,metric,logprep,drift = prepared[k]
                for h in (1,2,4):
                    for turns in (1,2):
                        a = np.eye(2,dtype=complex)
                        b = np.eye(2,dtype=complex)
                        loga = logb = 0.
                        for step in range(1,3*m*turns+1):
                            slot = (2*m+step)%(3*m)
                            a = np.linalg.matrix_power(filters[slot],h)@a
                            b = np.linalg.matrix_power(filters[2*m],h)@b
                            sa,sb = np.linalg.norm(a),np.linalg.norm(b)
                            require(sa>0 and sb>0,'Null scaled product')
                            a /= sa; b /= sb; loga += math.log(sa); logb += math.log(sb)
                        ra,rb = a@rho@a.conj().T,b@rho@b.conj().T
                        cross = float(np.trace(metric@a@rho@b.conj().T).real)
                        aa,bb = float(np.trace(metric@ra).real),float(np.trace(metric@rb).real)
                        require(aa>0 and bb>0,'Invalid scaled pilot probability')
                        left,right = 2*loga+math.log(aa),2*logb+math.log(bb)
                        high = max(left,right)
                        logp = high+math.log(math.exp(left-high)+math.exp(right-high))-math.log(2)
                        normalized_cross = cross/(math.sqrt(aa)*math.sqrt(bb))
                        visibility = normalized_cross/math.cosh((left-right)/2)
                        la,lb = 2*loga+math.log(float(np.trace(ra).real)),2*logb+math.log(float(np.trace(rb).real))
                        top = max(la,lb)
                        log_rotation = top+math.log(math.exp(la-top)+math.exp(lb-top))-math.log(2)
                        p = math.exp(logp)
                        probability_rotation = math.exp(log_rotation)
                        marker = (1-2**-8)**(3*m*turns)
                        selectors = float(source.selector_success(5))**11
                        sign = -1 if turns==1 and law=='native' else 1
                        selected = (1+sign*visibility)/2
                        require(selected>0,'No selected path')
                        rows.append({'M':m,'law':law,'k':k,'h':h,'turns':turns,
                            'visibility':visibility,'normalized_cross':normalized_cross,
                            'rotation_success':probability_rotation,'rotation_reencode_success':p,
                            'log10_rotation_success':log_rotation/math.log(10),
                            'log10_rotation_reencode_success':logp/math.log(10),
                            'arm_product_natural_logs':[loga,logb],
                            'log10_clock_success_given_classical_prefix':logprep+logp/math.log(10)+math.log10(marker*selectors*selected),
                            'input_optical_trace_drift':drift,'native_reversals':sum(catalog[(2*m+j)%(3*m)]['native_reversals_per_filter'] for j in range(1,3*m*turns+1))*h})
    return {'candidate':'A-native-stationary-amplified-filter-v3','exploratory':True,
            'confirmatory_result':False,'actual_old_mixed_input':True,'reference':'literal stationary',
            'input_coordinate_note':'Invariant optical mixed density; raw float trace drift recorded, normalized only as exploratory numerical representation.',
            'full_composed_failure_ledger':False,'rows':rows,'calibration':calibration}


def compare(a,b,path='$'):
    if isinstance(a,dict):
        require(a.keys()==b.keys(),path+' keys')
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        require(len(a)==len(b),path+' length')
        for i,(x,y)in enumerate(zip(a,b)):compare(x,y,path+str(i))
    elif isinstance(a,(int,float))and not isinstance(a,bool):
        require(abs(a-b)<=1e-10+1e-12*max(abs(a),abs(b)),path+' numerical mismatch')
    else:require(a==b,path+' mismatch')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--first-record',action='store_true')
    args = parser.parse_args()
    path = ROOT/'pilots/AMPLIFIED_PILOT.json'
    require(args.first_record or path.is_file(),'Missing saved pilot; explicit --first-record required')
    require(not args.first_record or not path.exists(),'Refusing to replace saved evidence')
    result = run()
    if args.first_record:
        with path.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    else:compare(result,json.loads(path.read_text()))
    print(json.dumps({'verified':not args.first_record,'record_created':args.first_record,
        'rows':len(result['rows']),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'k6_M8_h4':[r for r in result['rows']if r['k']==6 and r['M']==8 and r['h']==4]},indent=2))


if __name__=='__main__':main()
