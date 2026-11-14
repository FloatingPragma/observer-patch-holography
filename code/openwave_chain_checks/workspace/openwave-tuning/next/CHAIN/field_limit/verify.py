"""Small adversarial checks for the controlled native field adapter."""
from copy import deepcopy
from fractions import Fraction as F
import json
import math

import numpy as np
import scipy.linalg as la

import controlled_field as c


def refuses(call):
    try:
        call()
    except (AssertionError, ValueError):
        return True
    return False


def short_actual_path(fixed, edges):
    initial=[1,35]
    counters=[[int(i==x)for i in range(96)]for x in initial]
    states=[tuple(v)for v in counters]
    records=[]
    native=fixed.literal_native()
    for side, moving in ((0,True),(0,False),(1,True),(1,False)):
        position=states[side].index(1)
        seam=next(e for e,ends in enumerate(edges)if position in ends)
        a,b=edges[seam]
        coin=int(position==b) if moving else int(position==a)
        states[side],(_,current)=native(states[side],edges[seam],coin)
        counters[side][a]-=current;counters[side][b]+=current
        records.append({'copy':side,'seam':seam,'coin':coin,
                        'current':current,'observer_current':current})
    assert [abs(v['current'])for v in records]==[1,0,1,0]
    return {'initial':initial,'records':records,'final_counters':counters,
            'terminal':{'kind':'cap_exhausted','time':.01}}


def run():
    pins=c.dependencies();fixed,n20=c.sources();edges,nb=fixed.graph()
    lap=6*np.eye(len(nb))
    for a,b in edges:lap[a,b]-=1;lap[b,a]-=1
    green=la.inv(lap+np.ones_like(lap)/len(nb))-np.ones_like(lap)/len(nb)
    path=short_actual_path(fixed,edges)
    baseline=c.track_path(path,fixed,edges,nb,green,lap)
    assert baseline['native_records']==4 and baseline['native_holds']==2
    dropped=deepcopy(path);dropped['records'][0]['observer_current']=0
    wrong_sign=deepcopy(path);wrong_sign['records'][0]['current']*=-1
    omitted_cap=deepcopy(path);omitted_cap['terminal']['kind']='success_conditioned'
    mutations={
        'omitted_native_current':refuses(lambda:c.track_path(dropped,fixed,edges,nb,green,lap)),
        'wrong_native_current_sign':refuses(lambda:c.track_path(wrong_sign,fixed,edges,nb,green,lap)),
        'omitted_cap_flag':refuses(lambda:c.track_path(omitted_cap,fixed,edges,nb,green,lap)),
    }
    original=c.local_fixed_point
    def bad_arithmetic(*args,**kwargs):
        out=original(*args,**kwargs);out[0]+=1<<40;return out
    c.local_fixed_point=bad_arithmetic
    mutations['corrupted_local_arithmetic']=refuses(lambda:c.track_path(path,fixed,edges,nb,green,lap))
    c.local_fixed_point=original
    # Exact rounding proof covers negative inputs and ties in both directions.
    rounding=max(abs(F(c.round_ratio(v))-F(v,12))for v in range(-1000,1001))
    assert rounding==F(1,2)
    # Independent full-Laplacian eigenvalues check the reused exact gap proof.
    eigen=la.eigvalsh(lap)
    assert abs(eigen[0])<1e-13 and eigen[1]>.5 and eigen[-1]<=12+1e-13
    result=json.loads((c.HERE/'RESULT.json').read_text())
    assert result['source_pins']==pins and all(result['checks'].values())
    assert [r['cells_per_axis']for r in result['static_refinement']]==[16,32,64]
    assert all(row['wrong_dimension_rejected'] and row['wrong_optical_mode_rejected']
               for row in result['static_refinement'])
    # A complete frozen cap law must be compared before any renormalization.
    for row in result['complete_cap_controls']:
        mu=576*row['time'];C=row['cap']
        probability=math.exp(-mu)*sum(mu**k/math.factorial(k)for k in range(C))
        assert abs(probability-row['conclusive'])<1e-15
        assert abs(1-probability-row['cap_probability'])<1e-15
        assert probability<.99
    mutations['success_only_renormalization']=all(
        abs(r['conclusive']-1)>.01 for r in result['complete_cap_controls'])
    # A constant-potential shift is unobservable and must not be falsely rejected.
    rho=np.array(path['final_counters'][0])-np.array(path['final_counters'][1])
    state=c.local_fixed_point([0]*96,rho,nb,32)
    assert c.residual_squared_certificate(state,rho,nb)==c.residual_squared_certificate(
        [z+(1<<40)for z in state],rho,nb)
    assert all(mutations.values()) and c.dependencies()==pins
    out={'schema':'oph.controlled_native_field_limit.adversarial_review.v1',
         'checks':{'actual_move_hold_readout':True,'exact_negative_rounding':True,
                   'source_gap_numeric_cross_check':True,'all_source_pins_verified':True,
                   'frozen_cap_law_independent_Poisson_check':True,
                   'constant_gauge_shift_preserved':True},
         'mutation_rejections':mutations,
         'rounding_max_exact':str(rounding),
         'result_sha256':c.sha(c.HERE/'RESULT.json'),
         'reviewer_independence':'Same author adversarial review; no external review claimed.'}
    (c.HERE/'REVIEW.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':run()
