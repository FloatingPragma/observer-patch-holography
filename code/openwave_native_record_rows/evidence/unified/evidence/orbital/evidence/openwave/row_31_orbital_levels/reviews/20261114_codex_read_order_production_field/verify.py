#!/usr/bin/env python3
"""Independent complete adaptive production replay, no producer imports."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import gzip,hashlib,json,math,random
import numpy as np
HERE=Path(__file__).resolve().parent
P=HERE.parent.parent/'work/2026-11-14_codex_production_field/failed_initial'
FACES=((0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),(1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),
       (3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),(4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1))

def make_graph(inner):
    cells=list(product(range(-2,3),repeat=3)) if inner else list(product(range(-1,2),range(-1,2),range(-1,2),range(3)))
    ci={v:i for i,v in enumerate(cells)};edges=[]
    internal=sorted({tuple(sorted((f[j],f[(j+1)%3]))) for f in FACES for j in range(3)})
    for x in cells:
        c=ci[x];edges += [(12*c+a,12*c+b,-1) for a,b in internal]
        for axis in range(len(x)):
            y=list(x);y[axis]+=1
            if inner:
                if y[axis]>2:continue
            elif axis==3:y[axis]%=3
            elif y[axis]>1:y[axis]=-1
            edges += [(12*c+p,12*ci[tuple(y)]+p,axis) for p in range(12)]
    adj=[[] for _ in range(12*len(cells))];pair={}
    for eid,(a,b,_) in enumerate(edges):adj[a].append(eid);adj[b].append(eid);pair[frozenset((a,b))]=eid
    return cells,ci,edges,adj,pair,internal

def main():
    expected_receipt='7df9e4455e35d23643a11201ad87f755498662c83e3f3b71f3bed734e7696be8'
    expected_records='342b68a0cb4a3bb15b7fec5621d8c960a8fe00f2e91041938302c042c6dde9b0'
    assert hashlib.sha256((P/'RECEIPT.json').read_bytes()).hexdigest()==expected_receipt
    assert hashlib.sha256((P/'RECORDS.json.gz').read_bytes()).hexdigest()==expected_records
    raw=gzip.decompress((P/'RECORDS.json.gz').read_bytes());record=json.loads(raw);receipt=json.loads((P/'RECEIPT.json').read_text())
    assert hashlib.sha256(raw).hexdigest()==receipt['raw_uncompressed_sha256']
    oc,oi,oe,oa,op,internal=make_graph(False);ic,ii,ie,ia,ip,_=make_graph(True);W=len(oe);NV=972;n=32768
    census=0
    for name,cells,edges in [('outer',oc,oe),('inner',ic,ie)]:
        assert record[name+'_graph']=={'cells':[list(x) for x in cells],'edges':[list(e) for e in edges]}
        rows=iter(record['branches'][name])
        for eid in range(len(edges)):
            for slot in range(3):
                for c in (0,1):
                    j=int(slot==0)-c if slot<2 else 0;after=[int(x==slot) for x in range(3)];after[0]-=j;after[1]+=j
                    assert next(rows)==[eid,slot,c,j,after.index(1)];census+=1
        assert next(rows,None) is None
    Q=[[F(16,22),F(4,22),F(0)],[F(2,22),F(16,22),F(2,22)],[F(0),F(3,22),F(16,22)]]
    hit=[F(1,22),F(0),F(0)];escape=[F(1,22),F(2,22),F(3,22)];eh=[F(0)]*3;ee=[F(0)]*3;tail=[F(1)]*3
    for _ in range(256):
        eh=[hit[i]+sum(Q[i][j]*eh[j] for j in range(3)) for i in range(3)]
        ee=[escape[i]+sum(Q[i][j]*ee[j] for j in range(3)) for i in range(3)]
        tail=[sum(Q[i][j]*tail[j] for j in range(3)) for i in range(3)]
    assert all(h+e+t==1 for h,e,t in zip(eh,ee,tail))
    assert max(abs(float(x)-y) for x,y in zip(ee,receipt['capped_escape']))<1e-15
    lint=np.zeros((12,12))
    for a,b in internal:lint[a,a]+=1;lint[b,b]+=1;lint[a,b]-=1;lint[b,a]-=1
    vals,vecs=np.linalg.eigh(lint);bands=[];projs=[];uvecs=[]
    for v in vals:
        if not bands or abs(v-bands[-1][0])>1e-8:bands.append([v,1])
        else:bands[-1][1]+=1
    for v,_ in bands:
        p=vecs[:,abs(vals-v)<1e-8]@vecs[:,abs(vals-v)<1e-8].T;projs.append(p);u=p[:,np.argmax(np.diag(p))];uvecs.append(u/np.linalg.norm(u))
    outer_calls=inner_calls=holds=0;stats=[];all_counts=[];band_errors=[];correct_bands=[]
    for law,run in enumerate(record['runs']):
        disabled=bool(law);assert run['disabled']==disabled;rng=random.Random((2026111433,2026111434)[law]);law_counts=[];law_stats=[]
        for channel in run['channels']:
            reader=rng.randrange(NV);assert reader==channel['start'];ipos=12*ii[(0,0,0)]
            origin=oc[reader//12][3];ticks=0;clock=[];occ=[0]*NV;inner_cursor=0;active=observed=0;expected=0.;censors=[]
            def inner_step(eid,coin):
                nonlocal ipos,inner_cursor,inner_calls
                a,b,_=ie[eid];j=int(ipos==a)-coin if ipos in (a,b) else 0
                assert channel['inner_events'][inner_cursor]==[eid,coin,j];inner_cursor+=1;inner_calls+=1
                if j:ipos=b if j==1 else a
                return j
            assert len(channel['macros'])==n
            for step,saved in enumerate(channel['macros'],1):
                eid=rng.randrange(W);coin=rng.randrange(2);a,b,axis=oe[eid];begin=inner_cursor;pc=qc=0;outcome='skip';hold=False
                risk=axis==-1 and reader in (a,b)
                if risk:
                    active+=1;target=oc[reader//12][:3];failed=False
                    for ax in range(3):
                        while ic[ipos//12][ax]!=target[ax]:
                            z=list(ic[ipos//12]);z[ax]+=1 if target[ax]>z[ax] else -1
                            dest=12*ii[tuple(z)]+ipos%12;pe=ip[frozenset((ipos,dest))]
                            for attempt in range(128):
                                j=inner_step(pe,rng.randrange(2));pc+=1
                                if j:break
                            else:failed=True
                            if failed:break
                        if failed:break
                    if failed:outcome='preparation_timeout'
                    else:
                        for attempt in range(257):
                            z=ic[ipos//12]
                            if z==(0,0,0):outcome='hit';break
                            if max(map(abs,z))==2:outcome='escape';break
                            if attempt==256:outcome='probe_timeout';break
                            pe=rng.choice(ia[ipos]);inner_step(pe,rng.randrange(2));qc+=1
                    if 'timeout' in outcome:censors.append([step,outcome])
                    hold=outcome!='escape' and not disabled;orbit=sum(x*x for x in target)
                    expected+=.5 if disabled else (.5*float(ee[orbit-1]) if orbit else 0.)
                if hold:j=0;holds+=1
                else:j=int(reader==a)-coin if reader in (a,b) else 0
                outer_calls+=1
                if j:reader=b if j==1 else a
                if risk and j:observed+=1
                assert saved==[eid,coin,j,int(hold),begin,inner_cursor,pc,qc,outcome]
                occ[reader]+=1
                if oc[reader//12][3]==(origin+ticks+1)%3:ticks+=1;clock.append([step,ticks,reader])
            assert inner_cursor==len(channel['inner_events'])
            assert occ==channel['occupation_counts'] and clock==channel['clock'] and censors==channel['censors']==[]
            assert active==channel['active_internal_proposals'] and observed==channel['observed_active_moves']
            assert expected==channel['exact_preparation_expected_active_moves']
            radius=math.sqrt(max(active,1)*math.log(128*n/.01)/2);assert radius==channel['martingale_radius']
            assert abs(observed-expected)<=radius+active*12/2**128
            y=(math.sqrt(NV/(4*W*n))*(np.array(occ)-n/NV)).reshape(81,12)
            fixed=[]
            for proj,saved in zip(projs,channel['band_readouts']):
                correct=y@proj;fixed.append(correct.tolist());band_errors.append(float(np.max(abs(correct-np.array(saved)))))
                # Preserve and exactly identify the producer's view-mutation defect.
                corrupted=proj.copy();column=np.argmax(np.diag(corrupted));corrupted[:,column]/=np.linalg.norm(corrupted[:,column])
                assert np.max(abs(y@corrupted-np.array(saved)))<1e-14
            correct_bands.append(fixed)
            constant=(y@np.ones(12)/math.sqrt(12)).reshape(3,3,3,3)
            assert np.max(abs(constant-channel['constant_readout']))<1e-14
            open_interval=n-(clock[-1][0] if clock else 0);assert open_interval==channel['open_clock_interval']
            law_counts.append(occ);law_stats.append(dict(ticks=ticks,open_interval=open_interval,active=active,observed_moves=observed,inner_calls=inner_cursor))
        av=np.array([x['constant_readout'] for x in run['channels']]);dirs=(3,0,1,2);curl=[]
        for mu in range(4):
            for nu in range(mu+1,4):curl.append(np.roll(av[nu],-1,dirs[mu])-av[nu]-np.roll(av[mu],-1,dirs[nu])+av[mu])
        assert np.max(abs(np.array(curl)-run['field_strength_readouts']))<1e-14
        stats.append(law_stats);all_counts.append(law_counts)
    assert (census,outer_calls,inner_calls,holds)==(82008,262144,2613,9)
    assert census+768+outer_calls+inner_calls==receipt['native_calls']==347533
    assert receipt['outer_history_preparation_error_bound']==8*n*12/2**128
    spatial=list(product(range(-1,2),repeat=3));si={x:i for i,x in enumerate(spatial)};ls=np.zeros((27,27))
    for x in spatial:
        for ax in range(3):
            y=list(x);y[ax]+=1
            if y[ax]>1:y[ax]=-1
            i,j=si[x],si[tuple(y)];ls[i,i]+=1;ls[j,j]+=1;ls[i,j]-=1;ls[j,i]-=1
    acceptance=np.array([float(ee[sum(v*v for v in x)-1]) if x!=(0,0,0) else 0 for x in spatial]);max_eig=0.
    for law,population in enumerate(receipt['population']):
        for bi,b in enumerate(population['bands'],1):
            lam=b['native_lambda'];sl=ls+lam*np.diag(np.ones(27) if law else acceptance);e,v=np.linalg.eigh(sl)
            max_eig=max(max_eig,float(np.max(abs(e-b['spatial_squared_levels']))));assert max_eig<1e-12
            assert min(e)>0
            s=e[0];g0=(1/s+2/(s+3))/3;g1=(1/s-1/(s+3))/3
            assert abs(b['ground']['first_hit_transform']-1/(1+s))<1e-12
            assert abs(g1/g0-1/(1+s))<1e-12
            for c,observed in zip(all_counts[law],b['actual_ground_mode_record_readouts']):
                y=math.sqrt(NV/(4*W*n))*(np.array(c)-n/NV)
                mode=np.kron(np.array([v[si[z[:3]],0]/math.sqrt(3) for z in oc]),uvecs[bi])
                assert abs(y@mode-observed)<1e-12
    (HERE/'CORRECT_BAND_READOUTS.json').write_text(json.dumps({'order':'source four channels, then disabled four channels; each internal band in eigenvalue order','bands':correct_bands},separators=(',',':'))+'\n')
    out={'pass':False,'blocking':'original band projector corrupted by in-place column normalization; corrected derived readout required','same_group':True,'G6_pass':False,'census_records_checked':census,'outer_macro_records_checked':outer_calls,
         'inner_records_checked':inner_calls,'auxiliary_holds_checked':holds,'all_adaptive_choices_reproduced':True,
         'producer_full_vector_calls_reviewed_as_code_not_reexecuted':768,'all_source_calls_accounted':347533,
         'same_tape_curl_clocks_and_mode_readouts_replayed':True,'original_band_decoder_bug_exactly_reproduced':True,
         'maximum_original_band_readout_error':max(band_errors),'independent_correct_band_readouts_saved':True,
         'stats':stats,'max_finite_spatial_spectrum_difference':max_eig,
         'fully_capped_marginal_exactly_Markov':False,'population_is_actual_long_tape':False,'hydrogen_shell_limit_reviewed_here':False,
         'receipt_sha256':expected_receipt,'records_sha256':expected_records}
    (HERE/'VERIFICATION.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out))

if __name__=='__main__':main()
