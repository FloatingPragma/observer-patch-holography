"""Two mobile signed native copies under one complete feedback clock."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

import field_feedback as fixed

HERE=Path(__file__).resolve().parent


def pair_table(neighbours,rounds=512):
    n=len(neighbours);rho=np.eye(n);rho[0,:]-=1
    phi=fixed.local_field(rho,neighbours,rounds)
    diagonal=np.diag(phi)
    table=(diagonal[:,None]+diagonal[None,:]-phi-phi.T)/2
    assert np.max(abs(table-table.T))<1e-14 and np.max(abs(np.diag(table)))<1e-14
    # Enforce the proved exact symmetry/diagonal identities after roundoff.
    table=(table+table.T)/2;np.fill_diagonal(table,0.)
    return table


def pair_rates(counters,edges,table,kappa=1.,beta=4.,nu=1.):
    if any(sum(c)!=1 or any(v not in (0,1) for v in c) for c in counters):
        raise ValueError('outside two one-load-copy domain')
    x,y=(int(np.argmax(c)) for c in counters)
    rate=np.full((2,len(edges)),nu)
    for side,pos in enumerate((x,y)):
        for e,(a,b) in enumerate(edges):
            if pos in (a,b):
                dest=b if pos==a else a
                value=table[dest,y] if side==0 else table[x,dest]
                rate[side,e]=nu*np.exp(-beta*kappa*(value-table[x,y])/2)
    return rate


def execute(initial,edges,table,kappa,seed,horizon=.5,cap=4096,omit_first_move=False):
    n=len(table);native=fixed.literal_native();rng=np.random.default_rng(seed)
    sources=[tuple(int(i==p) for i in range(n)) for p in initial]
    counters=[list(s) for s in sources];flows=[[0]*len(edges) for _ in range(2)]
    tape=[];t=0.;logp=0.;omitted=False
    for index in range(cap):
        rates=pair_rates(counters,edges,table,kappa);total=float(rates.sum())
        wait=float(rng.exponential(1/total))
        if t+wait>horizon:
            logp-=total*(horizon-t)
            terminal={'kind':'no_further_event','duration':horizon-t,'total_rate':total,
                      'survival_probability':float(np.exp(-total*(horizon-t)))}
            break
        mark=int(rng.choice(rates.size,p=rates.ravel()/total));side,e=divmod(mark,len(edges))
        coin=int(rng.integers(2));a,b=edges[e]
        sources[side],(_,j)=native(sources[side],edges[e],coin)
        old_positions=[int(np.argmax(c)) for c in counters]
        reported=j
        if omit_first_move and j and not omitted:reported=0;omitted=True
        counters[side][a]-=reported;counters[side][b]+=reported;flows[side][e]+=reported
        new_positions=[int(np.argmax(c)) for c in counters]
        delta=float(kappa*(table[tuple(new_positions)]-table[tuple(old_positions)]))
        t+=wait;logp+=-total*wait+np.log(rates[side,e]/2)
        tape.append({'index':index,'time':t,'wait':wait,'copy':side,'seam':e,
                     'coin':coin,'current':j,'observer_current':reported,
                     'seam_rate':float(rates[side,e]),'total_rate':total,
                     'record_energy_increment':delta})
        if reported!=j:
            terminal={'kind':'stopped_at_injected_record_omission','time':t}
            break
    else:terminal={'kind':'cap_exhausted','time':t}
    return {'initial':initial,'kappa':kappa,'seed':seed,'horizon':horizon,
            'records':tape,'terminal':terminal,'log_path_density':logp,
            'final_counters':counters,'signed_flow':[a-b for a,b in zip(*flows)],
            'source_matches_record_counters':all(tuple(c)==s for c,s in zip(counters,sources)),
            'omission_control':omit_first_move,'omission_executed':omitted}


def replay(path,edges,table):
    n=len(table);positions=list(path['initial']);counters=[[int(i==p) for i in range(n)] for p in positions]
    t=0.;logp=0.;work=[0.,0.];signed_flow=[0]*len(edges);omission=False
    initial_energy=path['kappa']*table[tuple(positions)]
    for r in path['records']:
        rates=pair_rates(counters,edges,table,path['kappa']);side,e=r['copy'],r['seam'];a,b=edges[e]
        old=positions[side];new=(a if r['coin'] else b) if old in (a,b) else old
        j=int(old==a)-int(new==a)
        assert j==r['current'] and abs(r['seam_rate']-rates[side,e])<1e-14
        if r['observer_current']!=j:omission=True
        jobs=r['observer_current'];before=[int(np.argmax(c)) for c in counters]
        counters[side][a]-=jobs;counters[side][b]+=jobs
        after=[int(np.argmax(c)) for c in counters]
        du=path['kappa']*(table[tuple(after)]-table[tuple(before)])
        assert abs(du-r['record_energy_increment'])<1e-14
        work[side]+=du;signed_flow[e]+=(1 if side==0 else -1)*jobs
        positions[side]=new;t+=r['wait'];assert abs(t-r['time'])<1e-12
        logp+=-rates.sum()*r['wait']+np.log(rates[side,e]/2)
    if path['terminal']['kind']=='no_further_event':
        rate=pair_rates(counters,edges,table,path['kappa']).sum()
        logp-=rate*(path['horizon']-t)
    assert abs(logp-path['log_path_density'])<1e-10
    final=[int(np.argmax(c)) for c in counters]
    assert abs(sum(work)-(path['kappa']*table[tuple(final)]-initial_energy))<1e-12
    rho0=np.array([int(i==path['initial'][0])-int(i==path['initial'][1]) for i in range(n)])
    divergence=np.zeros(n,dtype=int)
    for f,(a,b) in zip(signed_flow,edges):divergence[a]+=f;divergence[b]-=f
    assert np.array_equal(rho0-divergence,np.array(counters[0])-np.array(counters[1]))
    return {'works_by_copy':work,'energy_telescoping':True,'signed_current_Gauss':True,
            'corrupted_record_detected':omission,'counter_matches_hidden_source':final==positions}


def generator(table,neighbours,kappa,beta=4.):
    n=len(table);rr=[];cc=[];vv=[];work=np.zeros((2,n*n))
    for x in range(n):
        for y in range(n):
            src=x*n+y;out=0.
            for side,p in enumerate((x,y)):
                for dest in neighbours[p]:
                    xx,yy=(int(dest),y) if side==0 else (x,int(dest))
                    du=table[xx,yy]-table[x,y]
                    rate=float(.5*np.exp(-beta*kappa*du/2))
                    rr.append(xx*n+yy);cc.append(src);vv.append(rate);out+=rate
                    work[side,src]+=rate*du
            rr.append(src);cc.append(src);vv.append(-out)
    return coo_matrix((vv,(rr,cc)),shape=(n*n,n*n)).tocsc(),work


def run():
    paths_to_pin=(*fixed.PINS,Path(fixed.__file__))
    before={str(p.relative_to(fixed.ROOT)):fixed.sha(p) for p in paths_to_pin}
    edges,neighbours=fixed.graph();table=pair_table(neighbours);n=len(table)
    # Reciprocal difference: changing either position changes the other's rates.
    maximum_response=0.
    for x in (1,35,95):
        c0=[int(i==x) for i in range(n)]
        r0=pair_rates([c0,[int(i==0) for i in range(n)]],edges,table)
        r1=pair_rates([c0,[int(i==12) for i in range(n)]],edges,table)
        maximum_response=max(maximum_response,float(max(abs(r0[0]-r1[0]))))
    assert maximum_response>0
    p=np.ones(n*n)/(n*n);laws=[]
    for kappa in (0.,1.,2.):
        q,work=generator(table,neighbours,kappa)
        pi=np.exp(-4*kappa*table.ravel());pi/=sum(pi)
        db=q.multiply(pi[None,:]);difference=db-db.T
        error=0. if not difference.nnz else float(max(abs(difference.data)))
        assert error<1e-16 and max(abs(q@pi))<1e-16
        evolved=expm_multiply(.5*q,p)
        assert min(evolved)>=0 and abs(sum(evolved)-1)<2e-13
        split=work@p
        assert abs(split[0]-split[1])<1e-14
        assert max(abs(work[0].reshape(n,n)-work[1].reshape(n,n).T).flat)<1e-14
        laws.append({'kappa':kappa,'complete_joint_position_count':n*n,
                     'complete_joint_probability_sum':float(sum(evolved)),
                     'detailed_balance_residual':error,
                     'unit_energy_current_rate_by_copy':split.tolist(),
                     'initial_mean_unit_energy':float(table.ravel()@p),
                     'mean_unit_energy_at_half_time':float(table.ravel()@evolved),
                     'complete_joint_probability_sha256':hashlib.sha256(evolved.tobytes()).hexdigest(),
                     'each_copy_marginal':evolved.reshape(n,n).sum(axis=1).tolist()})
    assert max(abs(v) for v in laws[0]['unit_energy_current_rate_by_copy'])<1e-14
    assert all(v<0 for v in laws[1]['unit_energy_current_rate_by_copy'])
    paths=[];event_count=0;hold_count=0;copy_moves=[0,0]
    for initial in ((1,35),(95,12),(0,0),(5,0)):
        for kappa in (0.,1.,2.):
            for seed in (173,804):
                path=execute(initial,edges,table,kappa,seed)
                audit=replay(path,edges,table);assert audit['counter_matches_hidden_source']
                path['independent_replay']=audit;paths.append(path)
                event_count+=len(path['records'])
                for r in path['records']:
                    if r['current']:copy_moves[r['copy']]+=1
                    else:hold_count+=1
    assert all(v>0 for v in copy_moves)
    omission=execute((1,35),edges,table,1.,173,horizon=4.,omit_first_move=True)
    omission_audit=replay(omission,edges,table)
    assert omission['omission_executed'] and omission_audit['corrupted_record_detected']
    # Even if a later chance reunion occurred, the first discarded transfer
    # is a conclusive intermediate continuity mismatch; retain complete path.
    first=next(r for r in omission['records'] if r['current']!=r['observer_current'])
    omission['first_corrupt_event']=first['index'];omission['independent_replay']=omission_audit
    # Exact intervention at the first corrupted record: both are valid
    # one-load counters, but the true versus stale settings differ.
    proper=[[int(i==x) for i in range(n)] for x in omission['initial']]
    stale=[c[:] for c in proper]
    for r in omission['records'][:first['index']+1]:
        side,e=r['copy'],r['seam'];a,b=edges[e]
        proper[side][a]-=r['current'];proper[side][b]+=r['current']
        stale[side][a]-=r['observer_current'];stale[side][b]+=r['observer_current']
    omitted_rate_error=float(max(abs(pair_rates(proper,edges,table)-pair_rates(stale,edges,table)).flat))
    assert omitted_rate_error>0
    cap=execute((1,35),edges,table,1.,71,horizon=100.,cap=1)
    assert cap['terminal']['kind']=='cap_exhausted';cap['independent_replay']=replay(cap,edges,table)
    after={str(p.relative_to(fixed.ROOT)):fixed.sha(p) for p in paths_to_pin};assert before==after
    eps=2*(23/24)**512
    result={'schema':'oph.native_reciprocal_feedback.development.v1',
      'checks':{'both_mobile_native_copies_executed':True,'current_only_counter_controller':True,
                'same_clock_and_complete_hold_wait_law':True,'reciprocal_cross_response':True,
                'exchange_symmetric_full_generators':True,'both_copy_energy_work_telescopes':True,
                'all_joint_output_probability_retained':True,'kappa_zero_and_record_omission_controls':True,
                'all_sources_unchanged':True},
      'domain':{'ports_per_copy':n,'seams_per_copy':len(edges),'joint_source_states':n*n,
                'executed_paths':len(paths),'native_events':event_count,'zero_current_records':hold_count,
                'moves_by_copy':copy_moves,'official_rows_changed':[]},
      'architecture':{'charge':'rho=delta_x-delta_y, signed readouts of two declared native one-load copies',
                      'energy':'kappa*R_K(x,y)/2, computed by the same finite neighbour circuit',
                      'reference':'both signed copies move; no pinned compensator',
                      'coincidence':'two carrier loads survive at coincidence; no annihilation rule added',
                      'scheduler':'one joint Poisson clock; fair native coins; beta4, nu1 for each seam',
                      'initial_state':'declared pair positions; evolving source positions never passed to controller'},
      'joint_laws':laws,'max_cross_source_rate_response':maximum_response,
      'record_omission_rate_error':omitted_rate_error,
      'exact_arithmetic_full_path_TV_bound_time_half_kappa1':float(.5*12*np.exp(4)*np.expm1(4*eps)),
      'finite_precision_scope':'Exact truncation bound excludes floating rounding, as in fixed-reference proof.',
      'complete_paths':paths,'current_omission_control':omission,'cap_failure':cap,
      'source_pins':before,'auditor_sha256':fixed.sha(Path(__file__))}
    (HERE/'RECIPROCAL_FEEDBACK.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'domain':result['domain'],'work_rates':[r['unit_energy_current_rate_by_copy'] for r in laws],
                      'cross_response':maximum_response,'omission_rate_error':omitted_rate_error}))


if __name__=='__main__':run()
