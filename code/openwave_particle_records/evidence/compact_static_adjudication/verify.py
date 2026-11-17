#!/usr/bin/env python3
"""Exact complete-coordinate source audit; write only into this packet."""
import argparse,gzip,hashlib,json,time
from pathlib import Path
import numpy as np
import encoded_source as S

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PROVIDER=HERE.parent/'localized_rest_join'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(x):return json.dumps(x,sort_keys=True,indent=2)+'\n'

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true');args=parser.parse_args()
    start=time.perf_counter()
    prereg=json.loads((HERE/'PREREG.json').read_text())
    for f,h in prereg['files'].items():assert digest(ROOT/f)==h,f
    calls=0;original=S.native.native
    def native(*a):
        nonlocal calls
        calls+=1;return original(*a)
    S.native.native=native
    table=S.endpoint_table()
    assert calls==18
    from model.federated_observer import source
    graph=source.federation(0);seams=tuple(graph.seams)
    retained=json.loads((PROVIDER/'evidence/graph.json').read_text())
    assert [list(e) for e in seams]==retained['seams'] and len(seams)==630
    labels,p,r,valid=S.pair_chart();valid_indices=labels[valid]
    assert valid.sum()==57361 and (~valid).sum()==8175
    assert np.array_equal(labels[valid& (labels!=0)],(256*p+r)[valid& (labels!=0)])
    sums={'positive':np.zeros(256,dtype=np.int64),'negative':np.zeros(256,dtype=np.int64),'pair':np.zeros(65536,dtype=np.int64)}
    f=np.zeros(256,dtype=np.int64);g=f.copy();f[[4,5,6,7]]=[1,-1,-1,1];g[[16,17,18,19]]=[1,-1,-1,1];fg=np.outer(f,g).reshape(-1)
    counts=[];offsets=[0];changed_inputs=[];changed_outputs=[];changed_currents=[]
    action=hashlib.sha256();scalar_checks=0;contacts=0;continuity=0
    prep_columns=0
    held=np.empty_like(table)
    for c in (0,1):
        for x in (-1,0,1):
            for y in (-1,0,1):held[c,x+1,y+1]=(x,y,0)
    for ei,seam in enumerate(seams):
        a,b=seam
        spectator=next(v for v in range(240) if v not in seam)
        examples=[0,256*a+b,256*b+a,256*a+spectator,256*spectator+a,256*b+spectator,256*spectator+b,65535,257]
        for c in (0,1):
            out,j,mask=S.pair_event(table,seam,c)
            assert np.array_equal(mask,valid)
            assert np.all(valid[out[valid]])
            assert np.array_equal(out[~valid],labels[~valid]) and np.all(j[~valid]==0)
            zero=valid&(j==0);assert np.array_equal(out[zero],labels[zero])
            for current in (-1,1):
                dest=out[valid&(j==current)];assert len(np.unique(dest))==len(dest)
            counts.append([int(np.count_nonzero(valid&(j==cur))) for cur in (-1,0,1)])
            assert sum(counts[-1])==57361
            selected=valid&(j!=0)
            changed_inputs.append(labels[selected].astype(np.uint16));changed_outputs.append(out[selected]);changed_currents.append(j[selected])
            offsets.append(offsets[-1]+int(selected.sum()))
            action.update(out.astype('<u2').tobytes());action.update(j.tobytes())
            for label in examples:
                expected,current,flag=S.apply_pair_column(table,seam,c,label)
                assert (int(out[label]),int(j[label]))==(expected,current)
                assert flag==('body' if valid[label] else 'invalid')
                scalar_checks+=1
            assert int(out[256*a+b])==0 and int(j[256*a+b])==1
            assert int(out[256*b+a])==0 and int(j[256*b+a])==-1
            contacts+=2
            po=np.where(out==0,-1,out//256);ro=np.where(out==0,-1,out%256)
            for lo,hi in ((4,8),(16,20)):
                before=((p>=lo)&(p<hi)).astype(np.int8)-((r>=lo)&(r<hi)).astype(np.int8)
                after=((po>=lo)&(po<hi)).astype(np.int8)-((ro>=lo)&(ro<hi)).astype(np.int8)
                assert np.all((after-before)[valid]==(-j*(int(lo<=a<hi)-int(lo<=b<hi)))[valid])
                continuity+=int(valid.sum())
            nz=np.flatnonzero(fg);np.add.at(sums['pair'],out[nz],fg[nz])
            for name,vec,sign in (('positive',f,1),('negative',g,-1)):
                so,sj,sv=S.single_event(table,seam,c,sign)
                assert sv.sum()==240 and np.all(sv[so[sv]])
                assert np.array_equal(so[~sv],np.flatnonzero(~sv))
                for cur in (-1,0,1):
                    dest=so[sv&(sj==cur)];assert len(np.unique(dest))==len(dest)
                nn=np.flatnonzero(vec);np.add.at(sums[name],so[nn],vec[nn])
                prep_columns+=256
    assert np.array_equal(sums['positive'],1254*f)
    assert np.array_equal(sums['negative'],1254*g)
    assert np.array_equal(sums['pair'],1248*fg)
    hold_out,hold_j,_=S.pair_event(held,seams[0],0)
    assert np.array_equal(hold_out,labels) and not np.any(hold_j)
    assert 1260!=1254 and 1260!=1248
    # Wrong chart changes a source column while leaving its reader untouched.
    permutation=np.arange(65536);permutation[256*4+16],permutation[256*8+16]=permutation[256*8+16],permutation[256*4+16]
    chart_mutant=0
    for e in seams:
        out,j,_=S.pair_event(table,e,0)
        chart_mutant+=int(np.count_nonzero(permutation[out[permutation]]!=out))
    assert chart_mutant>0
    arrays=np.load(PROVIDER/'evidence/raw_source_arrays.npz')
    for name,vec in (('f',f),('g',g)):
        U=arrays['U_'+name];inverse=arrays['inverse_'+name]
        assert U.shape==(256,256) and inverse.shape==(256,256)
        assert np.linalg.norm(U[:,0]-1j*vec/2)<1e-10
        assert np.linalg.norm(inverse@U-np.eye(256),'fro')<1e-10
    laws=np.load(PROVIDER/'evidence/complete_outcome_laws.npz');single={}
    for law in ('C','E'):
        values={}
        for mode,wanted in (('coherent',1.),('dephased',.25),('sign_plant',.25)):
            joint=laws[f'{law}_energy_{mode}'];assert joint.shape==(256,256)
            marginal=joint.sum(axis=1)
            assert abs(marginal.sum()-1)<5e-10 and abs(marginal[0]-wanted)<5e-10
            values[mode]=float(marginal[0])
        pos=laws[f'{law}_position_coherent'].sum(axis=1)
        assert np.max(abs(pos-(f/2)**2))<5e-10
        single[law]=values
    body_columns=630*2*57361
    result={'schema':'oph.compact.encoded-interface.v1','candidate':'A-signed-compact-coordinate-v1',
      'all_pass':True,'native_calls_total':calls,'native_endpoint_table_calls':18,
      'scalar_reference_columns':scalar_checks,'all_pair_columns_per_context':65536,
      'pair_valid_columns':57361,'pair_invalid_columns':8175,'body_contexts':1260,
      'valid_pair_body_columns_checked':body_columns,'pair_all_columns_checked':1260*65536,
      'single_all_columns_checked':prep_columns,'full_regional_continuity_checks':continuity,
      'contact_columns':contacts,'nonzero_current_records':offsets[-1],
      'complete_action_sha256':action.hexdigest(),'body_effect_sum':'I on valid sector; invalid identity flag on full complement',
      'amplitude_phase_convention':'positive real 1/sqrt(1260), explicitly declared quantum lift',
      'single_energy_certificate':single,'wrong_chart_failing_columns':chart_mutant,
      'held_law_energy':0,'native_record_counts_are_compressed_exact_instrument_not_trajectories':True,
      'external_review_or_platform_award':False,'elapsed_seconds':time.perf_counter()-start}
    if args.record:
        target=HERE/'evidence';target.mkdir(exist_ok=False)
        np.savez_compressed(target/'complete_pair_actions.npz',
            counts=np.array(counts,dtype=np.int64),offsets=np.array(offsets,dtype=np.int64),
            inputs=np.concatenate(changed_inputs),outputs=np.concatenate(changed_outputs),currents=np.concatenate(changed_currents),
            invalid=labels[~valid].astype(np.uint16),seams=np.array(seams,dtype=np.int16),endpoint_table=table)
        result['artifact_sha256']={'evidence/complete_pair_actions.npz':digest(target/'complete_pair_actions.npz')}
        (HERE/'RECEIPT.json').write_text(canonical(result))
    print(canonical(result))
    return result
if __name__=='__main__':main()
