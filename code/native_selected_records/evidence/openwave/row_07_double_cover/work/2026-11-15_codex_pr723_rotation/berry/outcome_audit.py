#!/usr/bin/env python3
"""Balanced native paths and all first-rejection outcomes from retained tables."""
from pathlib import Path
import copy
import json
import math
from fractions import Fraction as F
import numpy as np
import run
import independent_check as independent

HERE = Path(__file__).resolve().parent

def logsum(a, b):
    top = max(a, b)
    return top + math.log(math.exp(a - top) + math.exp(b - top))

def maps_from_records(table):
    fine = np.zeros((384, 3, 3))
    counts = np.zeros((3, 3), dtype=int)
    for z, ticket, coin, current, out in table:
        flag = (ticket * 2 + coin) * 3 + current + 1
        fine[flag, out, z] = 1 / np.sqrt(128)
        counts[out, z] += 1
    w = np.ones(384) / np.sqrt(384); u = w.copy(); u[0] -= 1; u /= np.linalg.norm(u)
    complete = fine - 2 * u[:, None, None] * np.einsum('r,rij->ij', u, fine)
    assert np.max(abs(np.einsum('rai,raj->ij', complete, complete) - np.eye(3))) < 1e-12
    return counts / 128, complete

def path_compare(pa, pb):
    assert len(pa) == len(pb)
    a, la = run.transport_optical(pa); b, lb = run.transport_optical(pb)
    ratio = math.exp(-abs(la - lb))
    visibility = 2 * ratio * float(a @ b) / (1 + ratio * ratio)
    assert -1 - 1e-14 <= visibility <= 1 + 1e-14
    visibility = min(1., max(-1., visibility))
    logp = logsum(2 * la, 2 * lb) - math.log(2) - len(pa) * math.log(3)
    return {'events_per_arm': len(pa), 'relative_visibility': run.num(visibility),
            'conditional_minus_probability': run.num((1 - visibility) / 2),
            'log10_original_success_probability': run.num(logp / math.log(10)),
            'log10_attempts_per_success': run.num(-logp / math.log(10)),
            'loop_overlap_initial': run.num(a @ run.INITIAL),
            'reference_overlap_initial': run.num(b @ run.INITIAL),
            'log10_initial_constant_amplitude_for_one_percent_final_leakage': run.num(-2 + min(la, lb) / math.log(10))}

def exact_compare(wa, wb):
    aa = [independent.rational_matrix(w, 'native')[0] for w in wa]
    bb = [independent.rational_matrix(w, 'native')[0] for w in wb]
    a, b = independent.product(aa); c, d = independent.product(bb)
    norm = a*a/3 + b*b + c*c/3 + d*d
    visibility = 2 * (a*c/3 + b*d) / norm
    success = norm / (2 * 3 ** len(wa))
    logp = (math.log(success.numerator) - math.log(success.denominator)) / math.log(10)
    return {'visibility': float(visibility), 'exact_visibility_numerator': str(visibility.numerator),
            'exact_visibility_denominator': str(visibility.denominator),
            'log10_success': logp, 'negative_cross_amplitude_exact': visibility < 0}

def outcomes(wa, wb, maps):
    a = run.INITIAL.copy(); b = run.INITIAL.copy(); la = lb = 0.
    stages = []; failure_mass = 0.
    # Compute source vectors in exact invariant chart at each success step.
    basis = np.array([[1.,1.],[-1.,1.],[0.,-2.]]) / np.array([np.sqrt(2),np.sqrt(6)])
    va = vb = np.array([0., 1.])
    for step, (x, y) in enumerate(zip(wa, wb), 1):
        p, ka = maps[tuple(x)]; q, kb = maps[tuple(y)]
        prior = logsum(2*la, 2*lb) - math.log(2)
        total = logsum(2*la, 2*lb)
        weight_a = math.exp(2*la-total); weight_b = math.exp(2*lb-total)
        ax = np.einsum('rij,j->ri',ka,a); bx = np.einsum('rij,j->ri',kb,b)
        probs = weight_a * np.sum(ax*ax,axis=1) + weight_b * np.sum(bx*bx,axis=1)
        assert abs(sum(probs)-1) < 1e-11
        failure_mass += math.exp(prior) * sum(probs[1:])
        stages.append({'stage':step,'log10_probability_of_reaching_stage':run.num(prior/math.log(10)),
                       'conditional_probabilities_by_eraser_flag0_through383':list(map(run.num,probs))})
        va=(basis.T@p@basis)@va/np.sqrt(3); vb=(basis.T@q@basis)@vb/np.sqrt(3)
        na=np.linalg.norm(va);nb=np.linalg.norm(vb);va/=na;vb/=nb;la+=math.log(na);lb+=math.log(nb)
        a=basis@va;b=basis@vb
    final=logsum(2*la,2*lb)-math.log(2)
    assert abs(failure_mass+math.exp(final)-1)<1e-10
    return {'stages':stages,'all_first_rejections_mass':run.num(failure_mass),
            'log10_final_success_probability':run.num(final/math.log(10)),
            'probability_conservation_error':run.num(abs(failure_mass+math.exp(final)-1))}

def execute():
    receipt=json.loads((HERE/'RECEIPT.json').read_text());records=json.loads((HERE/'RECORDS.json').read_text())
    loop=[tuple(w) for w in receipt['weight_loop']]
    maps={w:maps_from_records(table) for w,table in zip(loop,records['source']['native']['loop'])}
    a=(32,16,16);b=(16,32,16);c=(16,16,32)
    vertex_path=[a,b,a,c,b,c,a]
    zero=[tuple(x[k]+step*(y[k]-x[k])//16 for k in range(3))
          for x,y in zip(vertex_path,vertex_path[1:]) for step in range(1,17)]
    one=[w for w in loop for _ in range(2)];two=loop+loop
    assert sorted(one)==sorted(two)==sorted(zero)
    summaries={};exact={}
    for name,reference in [('two_windings',two),('zero_winding',zero)]:
        summaries[name]=path_compare([maps[w][0] for w in one],[maps[w][0] for w in reference])
        exact[name]=exact_compare(one,reference)
        assert abs(exact[name]['visibility']-summaries[name]['relative_visibility'])<2e-10
        assert abs(exact[name]['log10_success']-summaries[name]['log10_original_success_probability'])<2e-8
    all_outcomes={
        'original48':outcomes(loop,[a]*48,maps),
        'original192':outcomes([w for w in loop for _ in range(4)],[a]*192,maps),
        'balanced96_zero_winding':outcomes(one,zero,maps)}
    corrected=copy.deepcopy(receipt);corrections=[]
    for law, cases in corrected['products'].items():
        for dwell, tests in cases.items():
            for name, test in tests.items():
                key='conditional_minus_probability';p=test[key]
                if p<0 or p>1:
                    assert -1e-14<=p<=1+1e-14
                    test[key]=min(1.,max(0.,p))
                    corrections.append({'law':law,'dwell':dwell,'case':name,'original':p,'corrected':test[key]})
    corrected['probability_roundoff_corrections']=corrections
    result={'no_new_native_calls':True,'same_histogram_and_event_count':True,
            'summaries':summaries,'exact_checks':exact,
            'useful_balanced_fringe_threshold':-.9,
            'balanced_fringe_results':{k:v['relative_visibility']<-.9 for k,v in summaries.items()},
            'probability_roundoff_corrections':corrections,
            'scope':'Parameterized native-source comparison; preparation, erasure, path control and physical rotation mapping supplied.'}
    return result,all_outcomes,corrected

if __name__=='__main__':
    result,outcomes,corrected=execute()
    for name,obj in [('BALANCED_RECEIPT.json',result),('ALL_OUTCOMES.json',outcomes),('CORRECTED_RECEIPT.json',corrected)]:
        with (HERE/name).open('x') as f:f.write(run.canonical(obj))
    print(run.canonical({'summaries':result['summaries'],'balanced_fringe_results':result['balanced_fringe_results'],
                        'roundoff_corrections':len(result['probability_roundoff_corrections'])}),end='')
