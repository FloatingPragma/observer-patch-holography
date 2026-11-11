"""Exact finite attempted-event population law and certified persistent-defect interval."""
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import product
from pathlib import Path

from .model import (STATES, SEAMS, kernel, initial_distribution, evolve, observations, energy,
                    event, balanced, extinct_or_frozen, balance_tail_bound, persistent_probability, survival_probabilities)

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def encode(value):
    if isinstance(value,F):return {'exact':str(value),'value':float(value)}
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,list):return [encode(v) for v in value]
    return value


def certificate():
    matrix=kernel(); checks={}
    checks['all_1296_rows_normalized']=len(matrix)==1296 and all(sum(row.values())==1 for row in matrix.values())
    checks['all_edges_preserve_integer_load']=all(sum(s[:3])==sum(t[:3]) for s,row in matrix.items() for t in row)
    checks['all_edges_preserve_Z6_face_charge']=all((t[3]+(-t[3])%6)%6==0 for row in matrix.values() for t in row)
    checks['all_edges_nonincreasing_combined_energy']=all(energy(t)<=energy(s) for s,row in matrix.items() for t in row)
    checks['terminal_set_is_closed']=all(extinct_or_frozen(t) for s,row in matrix.items()
                                       if extinct_or_frozen(s) for t in row)
    checks['frozen_defects_have_constant_loads']=all(t==s for s,row in matrix.items()
                                                  if s[3] and len(set(s[:3]))==1 for t in row)
    q_gap=[]; lower=[]
    for s,row in matrix.items():
        total=sum(s[:3]); quotient,remainder=divmod(total,3)
        qmin=(3-remainder)*quotient**2+remainder*(quotient+1)**2
        q_gap.append(sum(x*x for x in s[:3])-qmin)
        if not balanced(s):
            lower.append(sum(p for t,p in row.items() if sum(x*x for x in t[:3])<=sum(x*x for x in s[:3])-2))
    checks['load_square_excess_at_most_16']=max(q_gap)==16
    checks['unbalanced_descent_probability_at_least_one_third']=min(lower)>=F(1,3)
    # All six attempt records have weight 1/6, even if their states coincide.
    records=[event(s,a,c)[1] for s in STATES for a,c in product(range(3),range(2))]
    record_counts=Counter('no_operation' if r['no_operation'] else 'applied_twist' if r['applied'] else 'load_only'
                          for r in records)
    checks['all_attempt_records_retained']=sum(record_counts.values())==7776 and all(record_counts.values())
    survival3=survival_probabilities(3);survival32=survival_probabilities(32)
    checks['all_states_terminal_probability_by_three_at_least_1_over_108']=max(survival3.values())<=F(107,108)
    checks['all_states_survival_by_32_at_most_1_over_25']=max(survival32.values())<=F(1,25)

    initial=initial_distribution(); p=initial
    trajectory=[]; total_law=Counter()
    for s,weight in initial.items():total_law[sum(s[:3])]+=weight
    charge_conservation=True; monotone=True; last=None
    for n in range(129):
        obs=observations(p)
        if last is not None:monotone &= obs['mean_energy']<=last
        last=obs['mean_energy']
        now=Counter()
        for s,weight in p.items():now[sum(s[:3])]+=weight
        charge_conservation &= now==total_law
        checks[f'epoch_{n}.normalized']=obs['total_probability']==1 and sum(obs['class_shares'])==1
        checks[f'epoch_{n}.load_balance_tail_bound']=obs['unbalanced_probability']<=balance_tail_bound(n)
        checks[f'epoch_{n}.uniform_terminal_tail_bound']=obs['unresolved_probability']<=F(1,25)**(n//32)
        trajectory.append(dict(attempts=n, **encode(obs), balance_tail_bound=encode(balance_tail_bound(n))))
        if n<128:p=evolve(p)
    checks['total_load_distribution_conserved']=charge_conservation
    checks['mean_combined_energy_monotone']=monotone
    obs=observations(p);lo=obs['frozen_defect_probability'];hi=lo+obs['unresolved_probability']
    eventual=persistent_probability()
    exact=sum(weight*eventual[state] for state,weight in initial.items())
    checks['eventual_abundance_harmonic_on_every_state']=all(
        eventual[state]==sum(weight*eventual[target] for target,weight in row.items())
        for state,row in matrix.items())
    checks['eventual_probabilities_in_unit_interval']=all(0<=value<=1 for value in eventual.values())
    checks['exact_limit_inside_finite_epoch_interval']=lo<=exact<=hi
    checks['nonzero_persistent_defect_lower_bound']=lo>0
    checks['persistent_interval_width_below_1e_8']=hi-lo<F(1,10**8)
    source_paths=[str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')]
    source_paths+=['model/defect_integer_lift/model.py','model/defect_integer_lift/federation.json']
    return dict(schema='oph.native_abundance_epochs.v1',
        scope='Exact native repair and register dynamics on a declared closed two-face triangular cell; no transfer to the full federation or physical thermal abundances.',
        construction=dict(ports=3,oriented_seams=[list(s) for s in SEAMS],faces=2,states=len(STATES),
            initial_loads='independent uniform integers 0..5; initial face flux zero',
            scheduler='each attempt selects one of three seams and one of two coins uniformly',
            clock='number of attempted local events, including no-operation events',
            temperature=None,physical_seconds=None,particle_dictionary=None),
        complete_attempt_records=dict(record_counts),
        uniform_survival_bound=encode(dict(maximum_after_3=max(survival3.values()),
            maximum_after_32=max(survival32.values()), block_attempts=32, contraction=F(1,25),
            formula='Pr(outside terminal set after n attempts) <= (1/25)^floor(n/32) for every initial distribution')),

        persistent_occupied_face_share=encode(dict(exact_limit=exact,lower=lo,upper=hi,width=hi-lo)),
        trajectory=trajectory, checks=checks,all_checks_pass=all(checks.values()),
        input_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(source_paths)})
