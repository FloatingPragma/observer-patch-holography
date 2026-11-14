"""Complete attempt-count instrument for both field-feedback executors."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib
import json

import numpy as np
from scipy.sparse import bmat,diags,csc_matrix
from scipy.sparse.linalg import expm_multiply

import field_feedback as fixed
import reciprocal_feedback as reciprocal

HERE=Path(__file__).resolve().parent


def count_generator(q, scheduled_rates, cap):
    """Each native event, including every holding mark, increments count.

    Levels 0..cap-1 retain position and count; level cap retains position
    and an absorbing cap flag.  This is exactly execute's early-cap rule.
    """
    n=q.shape[0];q=csc_matrix(q);a=q+diags(scheduled_rates)
    assert min(a.data)>=0
    blocks=[[None]*(cap+1) for _ in range(cap+1)]
    for count in range(cap):
        blocks[count][count]=-diags(scheduled_rates)
        blocks[count+1][count]=a
    blocks[cap][cap]=csc_matrix((n,n))
    return bmat(blocks,format='csc'),a


def ordered_wait_reference(a,rates,initial,cap,time):
    """Independent analytic one- and two-event waiting-time integrals."""
    n=len(rates);lx=rates[initial]
    p0=np.zeros(n);p0[initial]=np.exp(-lx*time)
    ax=a[:,initial].toarray().ravel()
    if cap==1:
        absorbed=ax*(-np.expm1(-lx*time))/lx
        return np.r_[p0,absorbed]
    delta=rates-lx
    h=np.empty(n)
    same=abs(delta)<1e-12
    h[same]=time*np.exp(-lx*time)
    h[~same]=np.exp(-lx*time)*(-np.expm1(-delta[~same]*time))/delta[~same]
    p1=ax*h
    integral0=-np.expm1(-lx*time)/lx
    absorbed=a@(ax*(integral0-h)/rates)
    return np.r_[p0,p1,absorbed]


def run():
    dependencies=[Path(fixed.__file__),Path(reciprocal.__file__),fixed.NATIVE]
    before={str(p.relative_to(fixed.ROOT)):fixed.sha(p) for p in dependencies}
    edges,neighbours=fixed.graph();controller=fixed.compile_controller(neighbours)
    table=reciprocal.pair_table(neighbours);n=len(neighbours)
    rows=[]
    for copies in (1,2):
        for kappa in (0.,1.,2.):
            if copies==1:q=csc_matrix(fixed.generator(edges,controller['energies'],kappa))
            else:q=reciprocal.generator(table,neighbours,kappa)[0]
            # Q's diagonal contains move loss only. Twice that loss is the
            # incident scheduled rate; all inactive seams have rate one.
            rates=copies*(len(edges)-6)-2*q.diagonal()
            assert max(rates)<(330 if copies==1 else 660)
            for cap in (1,2):
                augmented,a=count_generator(q,rates,cap)
                assert max(abs(np.asarray(augmented.sum(axis=0))).flat)<1e-12
                initial=35 if copies==1 else 1*n+35
                p=np.zeros(augmented.shape[0]);p[initial]=1.
                time=.001
                out=expm_multiply(time*augmented,p)
                reference=ordered_wait_reference(a,rates,initial,cap,time)
                error=float(max(abs(out-reference)))
                assert error<2e-13 and abs(sum(out)-1)<2e-13 and min(out)>-2e-15
                success=float(sum(out[:cap*q.shape[0]]));failure=float(sum(out[cap*q.shape[0]:]))
                if kappa==0:
                    mu=copies*len(edges)*time
                    expected=np.exp(-mu)*sum(mu**k/factorial(k) for k in range(cap))
                    assert abs(success-expected)<2e-13
                rows.append({'copies':copies,'kappa':kappa,'cap':cap,'horizon':time,
                             'admitted_position_states':q.shape[0],'success_probability':success,
                             'cap_probability':failure,'total_probability':success+failure,
                             'ordered_wait_integral_max_error':error,
                             'full_count_position_law_sha256':hashlib.sha256(out.tobytes()).hexdigest()})
    eps=2*F(23,24)**512;small=8*eps
    eupper=sum((F(1,factorial(k)) for k in range(9)),F())+F(1,factorial(9))/(1-F(1,10))
    assert eupper<F(11,4) and small<1
    assert F(11,4)**2/(1-small)<8
    # The exact effective-resistance maximum-principle bound supplies
    # |Delta U|<=kappa/2; approximation adds 2*kappa*epsilon.
    # beta<=4 and kappa<=2 then imply every incident rate<8.
    bounds=[]
    for copies in (1,2):
        envelope=copies*(len(edges)-6+6*8)
        mu=F(envelope,2);cap=4096
        log_bound=float(cap-float(mu)+cap*np.log(float(mu)/cap))
        bounds.append({'copies':copies,'max_scheduled_rate_exact_integer':envelope,
                       'horizon':'1/2','poisson_mean_upper_exact':str(mu),'event_cap':cap,
                       'cap_tail_upper_exact_expression':f'exp({cap}-{mu})*({mu}/{cap})^{cap}',
                       'natural_log_upper_bound_approx':log_bound,
                       'base10_log_upper_bound_approx':log_bound/np.log(10),
                       'rounding_scope':'Exact-arithmetic scheduler. Executed finite rates also satisfy this envelope with a wide numerical margin; no general library rounding theorem is claimed.'})
    actual=[]
    for cap in (1,2):
        one=fixed.execute(35,edges,controller,1.,89,horizon=.5,cap=cap)
        fixed.independent_replay(one,edges,controller)
        two=reciprocal.execute((1,35),edges,table,1.,89,horizon=.5,cap=cap)
        reciprocal.replay(two,edges,table)
        assert one['terminal']['kind']==two['terminal']['kind']=='cap_exhausted'
        assert len(one['records'])==len(two['records'])==cap
        actual.append({'cap':cap,'single_copy':one,'two_copies':two})
    after={str(p.relative_to(fixed.ROOT)):fixed.sha(p) for p in dependencies};assert before==after
    result={'schema':'oph.native_feedback.complete_capped_law.v1',
      'checks':{'count_augmented_generator_increments_all_holds':True,'complete_capped_laws_sum_one':True,
                'independent_ordered_wait_integrals':True,'kappa_zero_exact_Poisson_control':True,
                'small_cap_actual_executor_replays':True,'rigorous_ideal_scheduler_tail_expression':True,
                'source_bytes_unchanged':True},
      'semigroup_scope':'The original 96- and9216-position expm laws are uncapped source laws. The instrument executed with cap4096 has an additional nonzero cap failure output.',
      'complete_instrument':'Levels0..C-1 have diagonal loss-Lambda and next-level marked native kernel A=Q+diagLambda; levelC is absorbing and retains source position plus cap flag. Fine records retain all seam/coin marks and times. Sum surviving count levels for conclusive position outputs.',
      'count_controls':rows,'tail_bounds':bounds,'actual_small_cap_paths':actual,
      'source_pins':before,'auditor_sha256':fixed.sha(Path(__file__)),'official_rows_changed':[]}
    (HERE/'CAPPED_LAW.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks':result['checks'],'controls':len(rows),'tails':bounds}))


if __name__=='__main__':run()
