#!/usr/bin/env python3
"""Native global production field with retained inner first-hit feedback."""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,importlib.util,itertools,json,math,random,sys
import numpy as np

HERE=Path(__file__).resolve().parent
SEEDS=(2026111433,2026111434)
STEPS=32768
PROBE_CAP=256
ROUTE_CAP=128
POP_N=100_000_000
FACES=((0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),
       (1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),
       (3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),
       (4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1))
INTERNAL=sorted({tuple(sorted((f[j],f[(j+1)%3]))) for f in FACES for j in range(3)})

def enc(x):return (json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()

def graph(inner):
    cells=list(itertools.product(range(-2,3),repeat=3)) if inner else list(itertools.product(range(-1,2),range(-1,2),range(-1,2),range(3)))
    ci={z:i for i,z in enumerate(cells)};edges=[]
    for z in cells:
        c=ci[z];edges.extend((12*c+a,12*c+b,-1) for a,b in INTERNAL)
        for axis in range(len(z)):
            zz=list(z);zz[axis]+=1
            if inner:
                if zz[axis]>2:continue
            else:
                if axis==3:zz[axis]%=3
                elif zz[axis]>1:zz[axis]=-1
            edges.extend((12*c+k,12*ci[tuple(zz)]+k,axis) for k in range(12))
    adj=[[] for _ in range(12*len(cells))];pair={}
    for eid,(a,b,axis) in enumerate(edges):
        adj[a].append(eid);adj[b].append(eid);pair[frozenset((a,b))]=eid
    return dict(cells=cells,ci=ci,edges=edges,adj=adj,pair=pair)

def main():
    path=HERE/'inputs/native.py'
    assert hashlib.sha256(path.read_bytes()).hexdigest()=='dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
    spec=importlib.util.spec_from_file_location('production_native',path)
    mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod);native=mod.native
    outer,inner=graph(False),graph(True)
    assert (len(outer['edges']),len(inner['edges']))==(6318,7350)
    W=len(outer['edges']);NV=12*len(outer['cells'])
    branches={};moves={};full_checks=0;calls=0
    for name,g in [('outer',outer),('inner',inner)]:
        records=[];M=np.zeros((12*len(g['cells']),12*len(g['cells'])),dtype=np.int64)
        for eid,(a,b,axis) in enumerate(g['edges']):
            spectator=next(i for i in range(3) if i not in (a,b))
            for slot in range(3):
                for coin in (0,1):
                    loads=[0,0,0];loads[slot]=1
                    out,rec=native(loads,(0,1),coin);calls+=1
                    target=out.index(1)
                    assert rec==(coin,loads[0]-out[0])
                    records.append([eid,slot,coin,rec[1],target])
                    if slot<2 and slot!=target:M[(a,b)[slot],(a,b)[target]]+=1
                    if eid<64:
                        mapping=(a,b,spectator);x=[0]*len(M);x[mapping[slot]]=1
                        y,r=native(x,(a,b),coin);calls+=1;full_checks+=1
                        assert r==rec and y.index(1)==mapping[target]
        branches[name]=records;moves[name]=M
    census_calls=sum(map(len,branches.values()))
    assert census_calls==82008 and full_checks==768
    lint=np.zeros((12,12),dtype=np.int64)
    for a,b in INTERNAL:
        assert moves['outer'][a,b]==1 and moves['outer'][b,a]==1
        lint[a,a]+=1;lint[b,b]+=1;lint[a,b]-=1;lint[b,a]-=1
    vals,vecs=np.linalg.eigh(lint.astype(float));bands=[];projectors=[];internal_vectors=[]
    for val in vals:
        if not bands or abs(val-bands[-1][0])>1e-8:bands.append([float(val),1])
        else:bands[-1][1]+=1
    assert [d for _,d in bands]==[1,3,5,3]
    for lam,_ in bands:
        select=abs(vals-lam)<1e-8;proj=vecs[:,select]@vecs[:,select].T
        projectors.append(proj);u=proj[:,np.argmax(np.diag(proj))].copy();u/=np.linalg.norm(u);internal_vectors.append(u)
    assert max(abs(x-y) for (x,_),y in zip(bands,(0,5-math.sqrt(5),6,5+math.sqrt(5))))<1e-10
    for (lam,dimension),proj in zip(bands,projectors):
        assert np.max(abs(proj@proj-proj))<1e-12
        assert np.max(abs(proj-proj.T))<1e-12
        assert abs(np.trace(proj)-dimension)<1e-10
        poly=np.eye(12)
        for other,_ in bands:
            if abs(other-lam)>1e-8:poly=poly@(lint-other*np.eye(12))/(lam-other)
        assert np.max(abs(poly-proj))<1e-12
    assert np.max(abs(sum(projectors)-np.eye(12)))<1e-12

    # Exact orbit kernel from every interior starting cell and internal port.
    orbit_rows={}
    for z in itertools.product(range(-1,2),repeat=3):
        orbit=sum(x*x for x in z)
        if orbit==0:continue
        for port in range(12):
            start=12*inner['ci'][z]+port;den=2*len(inner['adj'][start]);counts=[0,0,0];hit=escape=0
            for eid in inner['adj'][start]:
                a,b,_=inner['edges'][eid]
                for coin in (0,1):
                    row=branches['inner'][6*eid+2*(0 if start==a else 1)+coin]
                    end=(a,b)[row[4]];zz=inner['cells'][end//12]
                    if zz==(0,0,0):hit+=1
                    elif max(map(abs,zz))==2:escape+=1
                    else:counts[sum(x*x for x in zz)-1]+=1
            got=tuple(F(x,den) for x in counts)+(F(hit,den),F(escape,den))
            expected=[(F(16,22),F(4,22),0,F(1,22),F(1,22)),
                      (F(2,22),F(16,22),F(2,22),0,F(2,22)),
                      (0,F(3,22),F(16,22),0,F(3,22))][orbit-1]
            assert got==tuple(expected)
            if orbit in orbit_rows:assert orbit_rows[orbit]==got
            orbit_rows[orbit]=got
    qf=[[orbit_rows[i+1][j] for j in range(3)] for i in range(3)]
    hf=[orbit_rows[i+1][3] for i in range(3)];ef=[orbit_rows[i+1][4] for i in range(3)]
    pf=[[F(int(i==j)) for j in range(3)] for i in range(3)]
    eh=[F(0)]*3;ee=[F(0)]*3
    for _ in range(PROBE_CAP):
        eh=[eh[i]+sum(pf[i][j]*hf[j] for j in range(3)) for i in range(3)]
        ee=[ee[i]+sum(pf[i][j]*ef[j] for j in range(3)) for i in range(3)]
        pf=[[sum(pf[i][k]*qf[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    exact_tail=[sum(row) for row in pf]
    assert all(eh[i]+ee[i]+exact_tail[i]==1 for i in range(3))
    Q=np.array(qf,dtype=float);hit_vector=np.array(hf,dtype=float)
    escape_by=np.array(ee,dtype=float);hit_by=np.array(eh,dtype=float);tail=np.array(exact_tail,dtype=float)
    assert max(tail)<(10/11)**PROBE_CAP
    assert max(abs(np.linalg.solve(np.eye(3)-Q,hit_vector)-np.array([5/22,1/11,1/22])))<1e-14

    def event(g,pos,eid,coin):
        nonlocal calls
        a,b,_=g['edges'][eid]
        out,rec=native((int(pos==a),int(pos==b)),(0,1),coin);calls+=1
        if rec[1]:
            assert pos==(a if rec[1]==1 else b);pos=b if rec[1]==1 else a
        return pos,rec[1]

    runs=[];all_inner_calls=0;all_aux=0
    for disabled,seed in enumerate(SEEDS):
        rng=random.Random(seed);channels=[]
        for channel in range(4):
            start=source=reader=rng.randrange(NV)
            ipos=12*inner['ci'][(0,0,0)];ireader=ipos
            origin=outer['cells'][start//12][3];ticks=0;clock=[];occ=[0]*NV
            macros=[];inner_events=[];censors=[];active=0;expected_moves=0.;observed_moves=0
            for step in range(1,STEPS+1):
                eid=rng.randrange(W);coin=rng.randrange(2);a,b,axis=outer['edges'][eid]
                begin=len(inner_events);prep_count=probe_count=0;outcome='skip';hold=False
                at_risk=axis==-1 and reader in (a,b)
                if at_risk:
                    active+=1;target=outer['cells'][reader//12][:3];prep_failed=False
                    # Route from retained inner position; no external state reset.
                    for ax in range(3):
                        while inner['cells'][ireader//12][ax]!=target[ax]:
                            z=list(inner['cells'][ireader//12]);z[ax]+=1 if target[ax]>z[ax] else -1
                            nxt=12*inner['ci'][tuple(z)]+ireader%12
                            peid=inner['pair'][frozenset((ireader,nxt))]
                            for attempt in range(ROUTE_CAP):
                                c=rng.randrange(2);ipos,j=event(inner,ipos,peid,c);inner_events.append([peid,c,j]);prep_count+=1
                                pa,pb,_=inner['edges'][peid]
                                if j:ireader=pb if j==1 else pa
                                assert ipos==ireader
                                if j:break
                            else:prep_failed=True
                            if prep_failed:break
                        if prep_failed:break
                    if prep_failed:outcome='preparation_timeout'
                    else:
                        for attempt in range(PROBE_CAP+1):
                            z=inner['cells'][ireader//12]
                            if z==(0,0,0):outcome='hit';break
                            if max(map(abs,z))==2:outcome='escape';break
                            if attempt==PROBE_CAP:outcome='probe_timeout';break
                            peid=rng.choice(inner['adj'][ireader]);c=rng.randrange(2)
                            ipos,j=event(inner,ipos,peid,c);inner_events.append([peid,c,j]);probe_count+=1
                            pa,pb,_=inner['edges'][peid]
                            if j:ireader=pb if j==1 else pa
                            assert ipos==ireader
                    if 'timeout' in outcome:censors.append([step,outcome])
                    hold=(outcome!='escape') and not disabled
                    orbit=sum(x*x for x in target)
                    expected_moves+=.5 if disabled else (.5*escape_by[orbit-1] if orbit else 0.)
                if hold:
                    out,rec=native((0,0),(0,1),coin);calls+=1;all_aux+=1
                    assert out==(0,0) and rec==(coin,0);j=0
                else:source,j=event(outer,source,eid,coin)
                if j:reader=b if j==1 else a
                assert source==reader
                if at_risk and j:observed_moves+=1
                occ[reader]+=1
                if outer['cells'][reader//12][3]==(origin+ticks+1)%3:
                    ticks+=1;clock.append([step,ticks,reader])
                macros.append([eid,coin,j,int(hold),begin,len(inner_events),prep_count,probe_count,outcome])
            # Independent current-only record reconstruction, with no source state.
            rpos=start;ii=12*inner['ci'][(0,0,0)];recount=[0]*NV;reclock=[];rticks=0
            for step,record in enumerate(macros,1):
                eid,c,j,hold,begin,end,pc,qc,outcome=record
                assert end-begin==pc+qc
                for peid,cc,jj in inner_events[begin:end]:
                    aa,bb,_=inner['edges'][peid]
                    if jj:assert ii==(aa if jj==1 else bb);ii=bb if jj==1 else aa
                a,b,_=outer['edges'][eid]
                if hold:assert j==0
                elif j:assert rpos==(a if j==1 else b);rpos=b if j==1 else a
                recount[rpos]+=1
                if outer['cells'][rpos//12][3]==(origin+rticks+1)%3:
                    rticks+=1;reclock.append([step,rticks,rpos])
            assert recount==occ and reclock==clock and ii==ireader
            # Time-uniform Hoeffding-Azuma at the randomly sized active menu.
            radius=math.sqrt(max(active,1)*math.log(32*4*STEPS/.01)/2)
            assert abs(observed_moves-expected_moves)<=radius+active*12/2**128
            assert len(censors)<=2
            y=math.sqrt(NV/(4*W*STEPS))*(np.array(occ)-STEPS/NV)
            band=[(y.reshape(-1,12)@proj).tolist() for proj in projectors]
            af=(y.reshape(-1,12)@np.ones(12)/math.sqrt(12)).reshape(3,3,3,3)
            channels.append(dict(start=start,macros=macros,inner_events=inner_events,
                occupation_counts=occ,band_readouts=band,constant_readout=af.tolist(),
                clock=clock,open_clock_interval=STEPS-(clock[-1][0] if clock else 0),
                censors=censors,active_internal_proposals=active,
                observed_active_moves=observed_moves,exact_preparation_expected_active_moves=expected_moves,
                martingale_radius=radius))
            all_inner_calls+=len(inner_events)
        # The local curl is read across four independent channels.
        av=np.array([c['constant_readout'] for c in channels]);directions=(3,0,1,2)
        curl=[(np.roll(av[nu],-1,directions[mu])-av[nu]-np.roll(av[mu],-1,directions[nu])+av[mu]).tolist()
              for mu in range(4) for nu in range(mu+1,4)]
        runs.append(dict(disabled=bool(disabled),seed=seed,channels=channels,field_strength_readouts=curl))

    spatial=list(itertools.product(range(-1,2),repeat=3));si={x:i for i,x in enumerate(spatial)}
    ls=np.zeros((27,27))
    for i,x in enumerate(spatial):
        for ax in range(3):
            y=list(x);y[ax]+=1
            if y[ax]>1:y[ax]=-1
            j=si[tuple(y)];ls[i,i]+=1;ls[j,j]+=1;ls[i,j]-=1;ls[j,i]-=1
    acceptance=np.array([escape_by[sum(v*v for v in x)-1] if x!=(0,0,0) else 0. for x in spatial])
    # Assemble the effective complete outer operator directly from native branch transitions.
    leff=[];population=[];max_res=max_cov=max_clock=0.
    for disabled in (False,True):
        mat=np.zeros((NV,NV))
        for eid,(a,b,axis) in enumerate(outer['edges']):
            factor=1. if disabled or axis!=-1 else acceptance[si[outer['cells'][a//12][:3]]]
            assert moves['outer'][a,b]==1 and moves['outer'][b,a]==1
            mat[a,a]+=factor;mat[b,b]+=factor;mat[a,b]-=factor;mat[b,a]-=factor
        leff.append(mat)
        rows=[]
        for band_index,(lam,multiplicity) in enumerate(bands):
            if band_index==0:continue
            sl=ls+lam*np.diag(np.ones(27) if disabled else acceptance)
            sv,su=np.linalg.eigh(sl);u=internal_vectors[band_index]
            assert min(sv)>0
            if disabled:assert abs(sv[0]-lam)<1e-10
            else:assert sv[0]<85*lam/99+1e-10
            ground_mode=[]
            for mode,sval in enumerate(sv):
                cv=[];green=[]
                for t in range(3):
                    val=sval+4*math.sin(math.pi*t/3)**2
                    v=np.kron(np.array([su[si[z[:3]],mode]*np.exp(2j*np.pi*t*z[3]/3)/math.sqrt(3)
                                       for z in outer['cells']]),u)
                    residual=float(np.linalg.norm(mat@v-val*v));max_res=max(max_res,residual)
                    assert residual<1e-10
                    z=1-val/(2*W)
                    correction=z*(-math.expm1(POP_N*math.log(z)))/(2*W*POP_N*(1-z)**2)
                    cv.append(1/val-1/(4*W)-correction);green.append(1/val)
                covariance=[sum(math.cos(2*math.pi*t*lag/3)*cv[t] for t in range(3))/3 for lag in (0,1)]
                g=[sum(math.cos(2*math.pi*t*lag/3)*green[t] for t in range(3))/3 for lag in (0,1)]
                err=max(abs((covariance[lag]+(1/(4*W) if lag==0 else 0))/g[lag]-1) for lag in (0,1))
                max_cov=max(max_cov,err);assert err<.001
                q=np.array([[1-(sval+2)/(2*W),1/(2*W)],[1/(2*W),1-(sval+2)/(2*W)]])
                hit=float(np.linalg.solve(np.eye(2)-q,np.ones(2)/(2*W))[0])
                ce=max(abs(hit-1/(1+sval)),abs(hit-g[1]/g[0]));max_clock=max(max_clock,ce);assert ce<1e-10
                if mode==0:ground_mode=dict(squared_energy=float(sval),slice_energy=2*math.asinh(math.sqrt(sval)/2),
                    finite_cycle_clock_energy=-math.log(hit),first_hit_transform=hit,covariance=covariance,green=g,
                    profile=su[:,0].tolist())
            observed=[]
            for c in runs[int(disabled)]['channels']:
                y=math.sqrt(NV/(4*W*STEPS))*(np.asarray(c['occupation_counts'])-STEPS/NV)
                sm=np.kron(np.array([su[si[z[:3]],0]/math.sqrt(3) for z in outer['cells']]),u)
                observed.append(float(y@sm))
            rows.append(dict(native_lambda=lam,internal_multiplicity=multiplicity,
                spatial_squared_levels=sv.tolist(),ground=ground_mode,actual_ground_mode_record_readouts=observed))
        population.append(dict(disabled=disabled,bands=rows))
    # Internal constant sector is exactly unchanged by the feedback.
    u0=np.ones(12)/math.sqrt(12);basis0=np.kron(np.eye(81),u0[:,None])
    invariant=float(np.max(abs((leff[0]-leff[1])@basis0)));assert invariant<1e-12
    raw=enc(dict(outer_graph={k:outer[k] for k in ['cells','edges']},inner_graph={k:inner[k] for k in ['cells','edges']},
        branches=branches,runs=runs))
    (HERE/'RECORDS.json.gz').write_bytes(gzip.compress(raw,mtime=0))
    receipt=dict(candidate='MFR-native-probe-binding',kind='actual global production source, finite field and clock join',
        native_calls=calls,census_calls=census_calls,full_vector_checks=full_checks,
        actual_outer_macro_calls=8*STEPS,actual_inner_calls=all_inner_calls,auxiliary_holds=all_aux,
        member=dict(outer_cells=[3,3,3,3],inner_boundary_radius=2,a=1,channels_per_law=4,macro_attempts_per_channel=STEPS),
        capped_escape=escape_by.tolist(),capped_hit=hit_by.tolist(),probe_tail=tail.tolist(),
        preparation_error_per_active_attempt=12/2**128,outer_history_preparation_error_bound=8*STEPS*12/2**128,
        exact_population_scope='exact preparation with capped probe; fully capped conditional history error bounded separately',
        population_block_length=POP_N,population_is_actual_long_tape=False,
        population=population,max_source_eigen_residual=max_res,max_covariance_relative_error=max_cov,
        max_clock_error=max_clock,massless_sector_change=invariant,
        actual_clock_ticks=[[len(c['clock']) for c in r['channels']] for r in runs],
        actual_final_open_intervals=[[c['open_clock_interval'] for c in r['channels']] for r in runs],
        actual_active_proposals=[[c['active_internal_proposals'] for c in r['channels']] for r in runs],
        actual_censors=[[c['censors'] for c in r['channels']] for r in runs],
        first_tick_cap_tail=(1-1/(2*W))**STEPS,
        hydrogen_ladder_numerically_tested=False,all_registered_checks='pass',official='untested',
        raw_uncompressed_sha256=hashlib.sha256(raw).hexdigest())
    assert calls==census_calls+full_checks+8*STEPS+all_inner_calls
    (HERE/'RECEIPT.json').write_bytes(enc(receipt))
    print(json.dumps({k:v for k,v in receipt.items() if k!='population'},indent=2,sort_keys=True))

if __name__=='__main__':main()
