#!/usr/bin/env python3
"""Exact field-spin continuity plus complete native source transfer witness."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RER = ROOT / 'reverse-engineering-reality'
sys.path.insert(0, str(RER))
sys.dont_write_bytecode = True
import numpy as np
from scipy.sparse import csr_matrix, diags, eye
from model.defect_field_bridge import full_domain
from model.native_spin_observer import circuits
from model.observer import native
from model.observer.source_scalar import compiler


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dependencies():
    paths = [HERE/'DESIGN.md', HERE/'TASK.md', Path(__file__)]
    paths += [RER/p for p in (
        'model/defect_field_bridge/full_domain.py',
        'model/native_spin_observer/circuits.py',
        'model/observer/native.py', 'model/observer/source_scalar/compiler.py',
        'model/exchange_sector/statistics_service.py',
        'model/common_source/native_routing.py',
        'model/common_source/record_code.py',
        'code/source_selection_model/response.json',
        'code/m1_source_realization/receipt.json')]
    return {str(p.relative_to(ROOT)): sha(p) for p in paths}


def zero(a, label):
    a = a.tocsr()
    a.eliminate_zeros()
    if a.nnz:
        raise AssertionError(label)


def spin(sites, basis=None):
    if basis is None:
        basis = list(range(1 << (2*sites)))
    where = {v:i for i,v in enumerate(basis)}
    out = []
    for site in range(sites):
        rows, cols, xs, ys, zs = [], [], [], [], []
        for col, mask in enumerate(basis):
            occupied = (mask >> (2*site)) & 3
            zs.append(1 if occupied == 1 else -1 if occupied == 2 else 0)
            if occupied in (1, 2):
                rows.append(where[mask ^ (3 << (2*site))])
                cols.append(col)
                xs.append(1)
                ys.append(1 if occupied == 1 else -1)
        shape = (len(basis), len(basis))
        out.append((csr_matrix((xs,(rows,cols)),shape=shape,dtype=np.int64),
                    csr_matrix((ys,(rows,cols)),shape=shape,dtype=np.int64),
                    diags(zs,format='csr',dtype=np.int64)))
    # Physical doubled generators are X, i*Yreal, Z.
    total = tuple(sum((row[a] for row in out),csr_matrix(shape,dtype=np.int64))
                  for a in range(3))
    return out,total,basis


def branches(sites, law, seam, basis):
    where = {mask:i for i,mask in enumerate(basis)}
    table = {(c,r):[] for c in (0,1) for r in range(16)}
    for coin in (0,1):
        for col,mask in enumerate(basis):
            raw,dest,sign,current = full_domain.branch_action(
                mask,law,coin,sites=sites,flavors=2,seam=seam)
            loads = [((mask>>(2*x))&3).bit_count() for x in range(sites)]
            if all(loads[x] <= 1 for x in seam):
                expected,(_,j) = native.native(loads,seam,coin)
                actual = tuple(((dest>>(2*x))&3).bit_count() for x in range(sites))
                assert current == j and actual == expected
            table[coin,raw].append((where[dest],col,sign,current))
    answer = []
    for (coin,raw),rows in table.items():
        shape = (len(basis),len(basis))
        v = csr_matrix(([r[2] for r in rows],([r[0] for r in rows],
                         [r[1] for r in rows])),shape=shape,dtype=np.int64)
        js = set(r[3] for r in rows)
        assert len(js) <= 1
        answer.append((coin,raw,v,next(iter(js),0)))
    return answer


def exact_audit(sites,law,seam,basis=None):
    local,total,basis = spin(sites,basis)
    n = len(basis)
    identity = eye(n,format='csr',dtype=np.int64)
    z = csr_matrix((n,n),dtype=np.int64)
    number = diags([m.bit_count() for m in basis],format='csr',dtype=np.int64)
    tx,ry,tz = total
    zero(tx@ry-ry@tx-2*tz,'su2 xy')
    zero(ry@tz-tz@ry-2*tx,'su2 yz')
    zero(tz@tx-tx@tz+2*ry,'su2 zx')
    one = [i for i,m in enumerate(basis) if m.bit_count()==1]
    casimir = (tx@tx-ry@ry+tz@tz)[one][:,one]
    zero(casimir-3*eye(len(one),format='csr',dtype=np.int64),'one-particle Casimir')
    assert set(tz.diagonal()[one]) == {-1,1}
    complete=z.copy();unital=z.copy();nonzero=0;negative_signs=0;checks=0
    joint = branches(sites,law,seam,basis)
    # Unnormalized spatial trace with up polarization; trace=sites.
    polarized=diags([int(m.bit_count()==1 and (m.bit_length()-1)%2==0)
                    for m in basis],format='csr',dtype=np.int64)
    next_polarized=z.copy();mean_j=0
    complete_records=[]
    for coin,raw,v,current in joint:
        e=v.T@v
        complete+=e;unital+=v@v.T
        nonzero+=int(v.nnz>0);negative_signs+=int(np.count_nonzero(v.data<0))
        zero(v@number-number@v,'number conservation')
        for axis,t in enumerate(total):
            zero(v@t-t@v,'spin conservation')
            for site in range(sites):
                before=local[site][axis]
                delta=v.T@before@v-before@e
                sign=-1 if site==seam[0] else 1 if site==seam[1] else 0
                flux=current*sign*(e@(local[seam[0]][axis]+local[seam[1]][axis]))
                zero(delta-flux,'local spin continuity')
                checks+=1
        next_polarized+=v@polarized@v.T
        mean_j+=current*int((e@polarized).diagonal().sum())
        complete_records.append({'coin':coin,'raw':raw,'current':current,
                                 'domain_rank':int(e.diagonal().sum())})
    zero(complete-2*identity,'complete effects')
    zero(unital-2*identity,'unital channel')
    zero(next_polarized-2*polarized,'polarized stationary state')
    assert mean_j==0
    return {'sites':sites,'law':law,'seam':list(seam),'dimension':n,
            'all_coin_raw_slots':32,'nonzero_slots':nonzero,
            'negative_CAR_columns':negative_signs,'local_continuity_checks':checks,
            'complete_records':complete_records,'one_particle_dimension':len(one),
            'casimir':'3/4','spin_z_spectrum':['-1/2','1/2'],
            'polarized_stationary':True,'stationary_mean_occupation_current':0}


def exact():
    full=[exact_audit(sites,law,seam) for sites,seam in ((2,(0,1)),(4,(0,3)))
          for law in ('C','E')]
    # All66 unordered pairs contain the30 native geometric seams.
    one=[1<<i for i in range(24)]
    twelve=[exact_audit(12,law,(a,b),one) for a in range(12)
            for b in range(a+1,12) for law in ('C','E')]
    # This exact global flavor swap preserves every native effect but flips Sz.
    _,total,basis=spin(2)
    rows=[];sgns=[]
    for mask in basis:
        out=0;sign=1
        for site in range(2):
            x=(mask>>(2*site))&3
            out|=({0:0,1:2,2:1,3:3}[x]<<(2*site))
            if x==3:sign=-sign
        rows.append(out);sgns.append(sign)
    f=csr_matrix((sgns,(rows,list(range(16)))),shape=(16,16),dtype=np.int64)
    zero(f.T@total[2]@f+total[2],'control spin reversal')
    defects=[]
    for law in ('C','E'):
        for coin,raw,v,j in branches(2,law,(0,1),basis):
            changed=f@v
            zero(changed.T@changed-v.T@v,'same record effects')
            d=changed@total[2]-total[2]@changed
            defects.append(int(np.sum(d.data*d.data)))
    assert max(defects)>0
    return {'full_Fock_audits':full,'twelve_port_one_particle_audits':twelve,
            'changed_law':{'all_effects_identical':True,
                           'max_doubled_Sz_intertwiner_squared_Frobenius':max(defects)},
            'exact_integer_checks_pass':True},f.toarray()


def execute_source(control):
    _,total,_=spin(2)
    doubled_spin=[total[0].toarray(),1j*total[1].toarray(),total[2].toarray()]
    word=list(circuits.emit(circuits.fswap(0)+circuits.fswap(2)))
    actual_f=compiler.execute(word,4)
    ferr=float(np.linalg.norm(actual_f-control))
    assert ferr<2e-7
    copy=[event for bit in range(4) for event in compiler.cnot(bit,9+bit)]
    arrays={'native_control_matrix':actual_f}
    results=[]
    expected_spin=np.array([sum(((m>>(2*x))&1)-((m>>(2*x+1))&1)
                               for x in range(2))/2 for m in range(16)])
    pointer_spin=expected_spin[np.arange(8192)>>9]
    for law in ('C','E'):
        actual,target,resources=full_domain.source_isometry(2,2,law)
        err=float(np.linalg.norm(actual-target))
        assert err<2e-7
        intertwiner=[]
        for t in doubled_spin:
            output=np.kron(np.eye(32),t)@actual
            intertwiner.append(float(np.linalg.norm(output-actual@t)))
        assert max(intertwiner)<2e-6
        changed=compiler.execute(word,9,actual)
        rev=float(np.linalg.norm(np.kron(np.eye(32),doubled_spin[2])@changed
                                 +changed@doubled_spin[2]))
        assert rev<2e-6
        measured={}
        for name,state,sign in (('native',actual,1),('spin_reversal',changed,-1)):
            initial=np.zeros((8192,16),complex);initial[:512]=state
            final=compiler.execute(copy,13,initial)
            laws=abs(final)**2
            normerr=float(np.max(abs(laws.sum(axis=0)-1)))
            wrong=np.array([laws[pointer_spin!=sign*expected_spin[col],col].sum()
                            for col in range(16)])
            assert normerr<2e-7 and float(wrong.max())<2e-12
            arrays[f'{law}_{name}_final_amplitudes']=final
            measured[name]={'normalization_error_max':normerr,
                            'wrong_spin_probability_by_input':wrong.tolist(),
                            'mean_spin_by_input':(pointer_spin@laws).tolist(),
                            'all_joint_addresses':8192,'all_Fock_inputs':16}
        arrays[f'{law}_native_dilation']=actual
        results.append({'law':law,'Frobenius_source_error':err,
                        'doubled_spin_intertwiner_errors':intertwiner,
                        'changed_doubled_Sz_anticommuting_error':rev,
                        'field_source_resources':resources,'pointer':measured})
    np.savez_compressed(HERE/'source_arrays.npz',**arrays)
    exports={str(d):{'branches':len(compiler.transfer_kraus(d)),
                    'complete':bool(np.array_equal(sum(k.conj().T@k for k in
                         compiler.transfer_kraus(d)),np.eye(6)))} for d in (2,4)}
    assert all(v['complete'] for v in exports.values())
    return {'cases':results,'control_source_matrix_error':ferr,
            'control_source_pulses':len(word),'pointer_source_pulses':len(copy),
            'complete_proper_code_exports':exports,
            'source_arrays_sha256':sha(HERE/'source_arrays.npz'),
            'numeric_scope':'Binary64 source execution, not outward-rounded certification.'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--freeze',action='store_true')
    parser.add_argument('--source',action='store_true')
    args=parser.parse_args()
    if args.freeze:
        path=HERE/'PREREG.json'
        if path.exists():raise SystemExit('Refusing to replace frozen preregistration')
        path.write_text(json.dumps({'schema':'oph.native_field_spin_prereg.v1',
            'timestamp_utc':datetime.now(timezone.utc).isoformat(),
            'sha256':dependencies()},indent=2,sort_keys=True)+'\n')
        print('Preregistered; no evidence run.')
        return
    registered=json.loads((HERE/'PREREG.json').read_text())
    assert registered['sha256']==dependencies(),'preregistered input changed'
    started=time.monotonic()
    result,control=exact()
    if args.source:result['source_execution']=execute_source(control)
    result.update({'schema':'oph.native_field_spin_current.v1',
        'status':'EXACT_NATIVE_FIELD_SPIN_CONSERVATION_AND_CURRENT',
        'physical_row5_complete':False,'physical_row6_complete':False,
        'preregistration_sha256':sha(HERE/'PREREG.json'),
        'elapsed_seconds':time.monotonic()-started,'source_sha256':dependencies()})
    name='source_results.json' if args.source else 'exact_results.json'
    (HERE/name).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'output':name,
                      'elapsed_seconds':result['elapsed_seconds'],
                      'source_execution':args.source}))


if __name__=='__main__':main()
