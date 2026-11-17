"""Exploratory full-source contrast/cooling repair, not confirmatory evidence."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import native_source as source
from amplified_native import compile_instrument,routing,OPTICAL,require
from selected_source import weighted_source,contrast_dilation,amplify
from pilot import compare

ROOT=Path(__file__).resolve().parent
PERMUTATION=np.eye(3)[[1,2,0]]


def preparation(law,h):
    raw=weighted_source(law)
    contrast=amplify(contrast_dilation(raw))
    cool=amplify(raw)
    rho=np.diag([0.,1.,0.]).astype(complex)
    logp=0.; stages=[]
    for name,instrument in [('contrast',contrast)]+[('cool',cool)]*h:
        for orientation in ('forward','inverse'):
            rho,p,prob=source.channel(rho,instrument[orientation],[0])
            source.check_probabilities(prob)
            logp+=math.log10(p)
            stages.append({'stage':name,'orientation':orientation,'success':p,'probabilities':prob.tolist()})
    weights=(OPTICAL.T@rho@OPTICAL).diagonal().real
    return rho,logp,{'law':law,'h':h,'selection_success':10**logp,
        'optical_weights':weights.tolist(),'constant_weight':float(np.trace(rho).real-sum(weights)),
        'contrast_optical_filter':source.encode(contrast['optical_filter']),
        'cool_optical_filter':source.encode(cool['optical_filter']),
        'contrast_phase':contrast['phase'],'cool_phase':cool['phase'],
        'full_stages':stages,'rho':source.encode(rho),
        'native_reversals':contrast['native_reversals_per_filter']+h*cool['native_reversals_per_filter']}


def decode_input(rho,cl):
    logp=0.; stages=[]
    for name,ops,keep in [('positive_counter',cl['qpe'],range(1,cl['Q']//2)),('choice_decode',cl['decode'],[0])]:
        rho,p,prob=source.channel(rho,ops,keep);source.check_probabilities(prob)
        logp+=math.log10(p)
        stages.append({'stage':name,'success':p,'probabilities':prob.tolist()})
    optical=OPTICAL.T@rho@OPTICAL
    drift=float(np.trace(optical).real-1)
    require(abs(drift)<1e-10,'Optical invariant defect')
    optical/=np.trace(optical).real
    reencoder=np.einsum('kxa,ab->kxb',cl['qpe'][1:cl['Q']//2],OPTICAL)
    metric=np.einsum('kxa,kxb->ab',reencoder.conj(),reencoder)
    return optical,metric,logp,{'full_stages':stages,'rho':source.encode(rho),'trace_drift':drift}


def run():
    preparations=[];inputs=[];rows=[]
    clocks={k:source.clock(k)for k in (0,6)}
    basis=PERMUTATION@OPTICAL
    catalogs={}
    for m in (16,32):
        for law in ('native','swap'):
            catalogs[m,law]=[compile_instrument(routing(slot,m),law)for slot in range(3*m)]
    for prep_law in ('native','swap'):
        for cool_h in (16,24,32):
            prepared,logprep,preprecord=preparation(prep_law,cool_h)
            preparations.append(preprecord)
            for k,cl in clocks.items():
                rho,metric,logdecode,inputrecord=decode_input(prepared,cl)
                inputs.append({'prep_law':prep_law,'cool_h':cool_h,'k':k,**inputrecord})
                # Rotation main experiment keeps native preparation; its control
                # changes only the target rotation source, with the same input.
                if prep_law!='native':continue
                for m,rotation_h in ((16,4),(32,2)):
                    for law in ('native','swap'):
                        catalog=catalogs[m,law]
                        filters=[basis.T@c['K']@basis for c in catalog]
                        for turns in (1,2):
                            a=np.eye(2,dtype=complex);b=a.copy()
                            for step in range(1,3*m*turns+1):
                                a=np.linalg.matrix_power(filters[(2*m+step)%(3*m)],rotation_h)@a
                                b=np.linalg.matrix_power(filters[2*m],rotation_h)@b
                            aa=float(np.trace(metric@a@rho@a.conj().T).real)
                            bb=float(np.trace(metric@b@rho@b.conj().T).real)
                            cross=float(np.trace(metric@a@rho@b.conj().T).real)
                            p=(aa+bb)/2;visibility=2*cross/(aa+bb)
                            raw_visibility=visibility
                            if law=='swap':
                                require(abs(visibility-1)<1e-10,'Adaptive scalar control theorem defect')
                                visibility=1.0  # Exact source-derived Kopt=I, not a clipped probability.
                            path_probabilities=[(1+visibility)/2,(1-visibility)/2]
                            selected=path_probabilities[1 if turns==1 else 0]
                            require(p>0 and selected>=0,'Invalid accepted path')
                            marker=(1-2**-16)**(3*m*turns)
                            selectors=float(source.selector_success(5))**2
                            total=10**(logprep+logdecode)*p*selected*marker*selectors
                            rows.append({'k':k,'cool_h':cool_h,'M':m,'rotation_h':rotation_h,
                                'law':law,'turns':turns,'visibility':visibility,'raw_visibility':raw_visibility,
                                'path_probabilities_plus_minus':path_probabilities,'selected_path':'minus'if turns==1 else'plus',
                                'returned_clock_available':selected>0,'rotation_reencode_success':p,
                                'selection_success':10**logprep,'counter_decode_success':10**logdecode,
                                'selected_clock_success_given_classical_prefix':total,
                                'log10_selected_clock_success_given_classical_prefix':math.log10(total)if total else None,
                                'native_reversals_rotation':rotation_h*sum(catalog[(2*m+j)%(3*m)]['native_reversals_per_filter']for j in range(1,3*m*turns+1))})
    return {'candidate':'A-native-amplified-selected-source-v4','exploratory':True,
        'confirmation':False,'new_actual_source':True,'reference':'literal stationary',
        'compiler':'record-derived positive adjoint-pair amplification',
        'counter_unchanged':True,'marker_cap':16,'preparations':preparations,'inputs':inputs,'rows':rows,
        'scope':'Nominal pilot; no finite implementation certification or full rotation first-failure ledger yet.'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--first-record',action='store_true');args=parser.parse_args()
    path=ROOT/'pilots/COOLING_PILOT.json'
    require(args.first_record or path.is_file(),'Saved record absent; explicit first-record required')
    require(not args.first_record or not path.exists(),'Refusing to overwrite saved record')
    result=run()
    if args.first_record:
        with path.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    else:compare(result,json.loads(path.read_text()))
    print(json.dumps({'verified':not args.first_record,'rows':len(result['rows']),
        'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'preparations':[{k:p[k]for k in ('law','h','selection_success','optical_weights')}for p in result['preparations']],
        'k6_h24':[r for r in result['rows']if r['k']==6 and r['cool_h']==24]},indent=2))


if __name__=='__main__':main()
