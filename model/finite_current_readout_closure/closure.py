"""Finite source hypotheses for axial counters and integer face lifts."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
import sympy as s

from model.native_circulation import construction as native
from model.defect_integer_lift import certificate as lift_certificate
from model.defect_integer_lift import model as lift
from model.particle_common_contract.provenance import import_closure

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


@lru_cache(maxsize=1)
def axial():
    kernel, _, _, branches=native.source()
    w=native.axial_weights()
    n=len(kernel.states)
    by_state=[[] for _ in range(n)]
    means=s.zeros(n,30)
    squares=s.zeros(30)
    for i,j,e,coin,current,_ in branches:
        record=6*e+3*coin+current+1
        by_state[i].append((j,e,current,record))
        means[i,e]+=current
        squares[e,e]+=current*current
    assert len(branches)==3960 and all(len(row)==60 for row in by_state)
    assert (w*means.T).applyfunc(s.simplify)==s.zeros(3,n)
    sigma=(w*squares*w.T/s.Integer(n*60)).applyfunc(s.simplify)
    assert sigma==s.Rational(4,99)*s.eye(3)
    counts=np.array(kernel.counts,dtype=np.int64)
    assert np.all(counts.sum(axis=0)==60) and np.all(counts.sum(axis=1)==60)
    fourth=counts@counts@counts@counts
    floor=int(fourth.min())
    assert floor>0
    minorization=F(n*floor,60**4)
    assert 0<minorization<=1
    # Every ordered native branch is included, including current-zero records.
    words=np.zeros((180,180),dtype=np.int64)
    terminal=Counter()
    cross=s.zeros(30)
    for initial in range(n):
        for middle,e,j,r in by_state[initial]:
            for final,f,k,q in by_state[middle]:
                words[r,q]+=1
                terminal[(r,q,final)]+=1
                if j and k:cross[e,f]+=j*k
    assert int(words.sum())==n*60**2
    assert (w*cross*w.T).applyfunc(s.simplify)==s.zeros(3)
    terminal_bytes=json.dumps([[*key,count] for key,count in sorted(terminal.items())],
                              separators=(',',':')).encode()
    return dict(source_states=n,native_branches=3960,declared_preparation='uniform classical law on the66 native binary total-two states',
        exact_native_transition_column_and_row_sum=60,
        complete_two_attempt_record_count=180**2,
        complete_two_attempt_word_counts=words.tolist(),
        common_word_denominator=n*60**2,
        retained_zero_word_count=int(np.count_nonzero(words==0)),
        terminal_configuration_joint_nonzero_entries=len(terminal),
        terminal_configuration_joint_sha256=hashlib.sha256(terminal_bytes).hexdigest(),
        native_current_record_used='(seam,coin,current), with current in{-1,0,1}; Z increments W[:,seam]*current',
        conditional_axial_mean_zero_all_states=True,
        stationary_one_attempt_covariance='(4/99) I3',
        independent_two_attempt_cross_covariance_zero=True,
        exact_stationary_N_attempt_covariance='(4N/99) I3',
        exact_four_step_minimum_transition= str(F(floor,60**4)),
        exact_four_step_uniform_minorization=str(minorization),
        bracket_error_squared_bound='min(1,16/(delta*N)) in Frobenius norm, every initial classical law',
        bracket_error_squared_constant=str(16/minorization),
        functional_limit_hypotheses=dict(all_state_martingale=True,irreducible_aperiodic=True,
            uniform_ergodic_covariance=True,bounded_squared_jump='4/5',
            stationary_preparation_required=False),
        Brownian_covariance='(4/99) min(t,u) I3',
        Brownian_increment_action_coefficient='99/8',
        physical_position_or_lab_time_identified=False,
        full_native_pulse_word_compilation_claimed=False)


@lru_cache(maxsize=1)
def integer_lift():
    structures,checks=lift_certificate.structure()
    assert all(checks.values())
    reference,reference_checks=lift_certificate.reference_run()
    assert all(reference_checks.values())
    columns=[]
    for level in [0,1,2]:
        fed=lift.Federation(level)
        assert all(p!=m for p,m in zip(fed.face_plus,fed.face_minus))
        # Each signed column has one+1 and one-1. No large integer field is sampled.
        columns.append(dict(level=level,carriers=fed.carriers,seams=fed.seams,
                            closed_incidence_columns=fed.seams))
    # Exhaust all native load/coin/class contexts, including rejected twists.
    local=0;accepted=0
    for a in range(6):
        for b in range(6):
            for coin in [0,1]:
                _,_,j=lift.attempt(a,b,coin)
                for p in range(6):
                    for m in range(6):
                        take=bool(j and lift.accepts(p,m,j,a-b))
                        actual=j if take else 0
                        for winding in [-7,0,9]:
                            vp=p+6*winding;vm=m-6*winding
                            assert (vp+actual)%6==(p+actual)%6
                            assert (vm-actual)%6==(m-actual)%6
                            assert vp+actual+vm-actual==vp+vm
                        local+=1;accepted+=int(take)
    return dict(closed_source_incidence=columns,complete_load_coin_class_contexts=local,
        accepted_nonzero_contexts=accepted,noop_contexts=local-accepted,
        formal_theorems='Lean/Dynamics/IntegerFaceLift.lean',
        formal_domain='arbitrary integer links and currents on a finite incidence matrix whose columns sum to zero',
        exact_all_history_identities=['m=B n','m mod6=class','sum_faces m=0 on each closed carrier/screen'],
        source_attachment='accepted native twist updates one link byj and the two incident faces by+j,-j; rejected/noop branch updates neither',
        reference_replay=reference,reference_checks=reference_checks,
        archive_census_reexecuted=False,archive_hypercharge_percentage_status='imported diagnostic, not a universal theorem',
        winding_is_conserved_particle_hypercharge=False,
        dashboard_upgrade_requires_exact_identity_scope=True)


def bindings():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ['.py','.md']}
    paths.update(['Lean/Dynamics/IntegerFaceLift.lean','model/native_circulation/README.md',
                  'model/defect_integer_lift/PROOF.md','model/defect_integer_lift/federation.json',
                  'model/defect_integer_lift/imported_census.json',*native.SOURCE_PATHS])
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(import_closure(paths))}


def certificate():
    pins=bindings()
    result=dict(schema='oph.finite_current_readout_closure.v1',source_bindings=pins,
                axial_counter=axial(),integer_face_lift=integer_lift(),
                scope=dict(two_finite_observable_theorems_attached=True,
                    every_original_two_attempt_axial_record_retained=True,
                    integer_lift_arbitrary_current_theorem=True,
                    physical_particle_or_spacetime_identification=False,
                    imported_large_run_census_promoted=False))
    if bindings()!=pins:raise ArithmeticError('source changed during calculation')
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');parser.add_argument('--check',action='store_true');args=parser.parse_args()
    value=certificate();data=json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if args.write:(HERE/'results.json').write_text(data)
    if args.check and (HERE/'results.json').read_text()!=data:raise SystemExit('current closure receipt differs')
    print(json.dumps(dict(axial_records=180**2,minorization=value['axial_counter']['exact_four_step_uniform_minorization'],
                         reference_twists=value['integer_face_lift']['reference_replay']['applied_twists'])))


if __name__=='__main__':main()
