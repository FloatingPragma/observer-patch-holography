"""Apply every full stored Kraus catalog to every stored input, independently."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time
import numpy as np

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[2]/'row_07_double_cover/work/2026-11-16_codex_stationary_readout_repair_v2/CONFIRMATION'


def require(ok,message):
    if not ok:raise ValueError(message)


def array(x):return np.asarray(x['real'])+1j*np.asarray(x['imag'])


def run():
    started=time.monotonic();index=json.loads((SOURCE/'INDEX.json').read_text());hashes={r['path']:r['sha256'] for r in index['files']}
    def read(name):
        blob=(SOURCE/name).read_bytes();require(hashlib.sha256(blob).hexdigest()==hashes[name],'bound instrument file '+name)
        return json.loads(gzip.decompress(blob))
    preparation=read('PREPARATION.json.gz');catalogs=read('CATALOG.json.gz');cache={};normmax=probmax=reachmax=polymax=0.;catalog_count=stage_count=branch_count=0
    def instrument(key,record):
        nonlocal normmax,polymax,catalog_count
        ops={o:array(record[o]) for o in ('forward','inverse')}
        for values in ops.values():
            normmax=max(normmax,float(np.linalg.norm(np.einsum('kxi,kxj->ij',values.conj(),values)-np.eye(3))))
        P=np.asarray([[float(__import__('fractions').Fraction(x)) for x in row] for row in record['counts']]);T=4*np.linalg.matrix_power(P,3)-3*P if record['preboost'] else P
        z=complex(float(__import__('fractions').Fraction(record['phase_record']['real'])),float(__import__('fractions').Fraction(record['phase_record']['imag'])))
        forward=(2*z-1)*T+(z-1)**2*np.linalg.matrix_power(T,3)
        polymax=max(polymax,float(np.linalg.norm(ops['forward'][0]-forward)),float(np.linalg.norm(ops['inverse'][0]-forward.conj().T)))
        cache[key]=ops;catalog_count+=1
    for sign,groups in preparation.items():
        for law,record in groups.items():
            for name in ('contrast','cool'):instrument(('prep',sign,law,name),record[name])
            reach=1.;failed=0.
            for stage in record['stages']:
                name,orientation=stage['kind'].split(':');F=array(stage['input_factor']);ops=cache['prep',sign,law,name][orientation]
                out=np.einsum('kxa,ab->kxb',ops,F);prob=np.sum(abs(out)**2,axis=(1,2));saved=np.array(stage['probabilities'])
                probmax=max(probmax,float(np.max(abs(prob-saved))));reachmax=max(reachmax,abs(reach-stage['reach_before']))
                failed+=reach*sum(prob[1:]);reach*=prob[0];stage_count+=1;branch_count+=len(prob)
            reachmax=max(reachmax,abs(reach-record['success']),abs(failed+reach-1))
    for sign,groups in catalogs.items():
        for law,frames in groups.items():
            for slot,record in enumerate(frames):instrument(('frame',sign,law,slot),record)
    cases=0
    for name in hashes:
        if not name.startswith('CASE_'):continue
        case=read(name);sign=str(case['charge_sign']);law=case['law'];reach=1.;failed=0.
        stages={step:[] for step in range(1,96*case['turns']+1)}
        for stage in case['target_first_failures']:stages[stage['step']].append(stage)
        for marker in case['marker_records']:
            reachmax=max(reachmax,abs(reach-marker['reach_before']))
            failures=sum(h['probability'] for h in marker['histories'] if h['outcome']=='timeout')
            failed+=reach*failures;reach*=marker['success']
            for stage in stages[marker['step']]:
                reachmax=max(reachmax,abs(reach-stage['reach_before']))
                if stage['kind']=='declared_filter_hold':continue
                F=array(stage['input_factor']);orientation=stage['orientation']
                rops=cache['frame',sign,law,stage['rot_slot']][orientation];fops=cache['frame',sign,law,stage['ref_slot']][orientation]
                r=np.einsum('kxa,ab->kxb',rops,F[:3]);f=np.einsum('kxa,ab->kxb',fops,F[3:])
                prob=np.sum(abs(r)**2+abs(f)**2,axis=(1,2));saved=np.array(stage['probabilities'])
                probmax=max(probmax,float(np.max(abs(prob-saved))))
                failed+=reach*sum(prob[1:]);reach*=prob[0];stage_count+=1;branch_count+=len(prob)
        reachmax=max(reachmax,abs(reach-case['raw_float_rotation_marker_success']),abs(failed+reach-1));cases+=1
    require(cases==256 and catalog_count==584,'complete instrument catalog and case coverage')
    require(max(normmax,probmax,reachmax,polymax)<1e-8,'full CP/source stage mismatch')
    return dict(passed=True,producer_imports=False,source_index_sha256=hashlib.sha256((SOURCE/'INDEX.json').read_bytes()).hexdigest(),
                complete_kraus_catalogs=catalog_count,cases=cases,actual_stage_applications=stage_count,
                branch_probabilities_recomputed=branch_count,maximum_kraus_completeness_residual=normmax,
                maximum_forward_amplitude_polynomial_residual=polymax,maximum_actual_stage_probability_residual=probmax,
                maximum_reach_or_completeness_residual=reachmax,threshold=1e-8,seconds=time.monotonic()-started,
                scope='Full source Kraus and stage-record consistency; independent floating calculation supplements exact unitarity/polynomial and interval proofs. No producer imports or source regeneration.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--record-to',type=Path);args=p.parse_args()
    if args.record_to is not None and args.record_to.exists():raise FileExistsError(args.record_to)
    result=run()
    if args.record_to is not None:
        with args.record_to.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    else:
        saved=json.loads((HERE/'INSTRUMENT_CHECK.json').read_text());a=dict(result);b=dict(saved);a.pop('seconds');b.pop('seconds');require(a==b,'saved instrument verification differs')
    print(json.dumps(result,indent=2))
