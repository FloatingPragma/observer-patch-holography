"""All-class native absorption without a physical freezeout-time identification."""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from model.native_abundance_epochs import model as native
from model.native_thermal_admission.thermal import classification, encode
from model.particle_common_contract.provenance import import_closure

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
CUTS=(0,1,2,4,8,16,32,64,128)


def bindings():
    paths={str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md')}
    paths.update({'model/native_thermal_admission/thermal.py','model/native_thermal_admission/PROOF.md',
        'model/native_thermal_admission/results.json','model/native_abundance_epochs/model.py',
        'model/native_abundance_epochs/PROOF.md','model/native_abundance_epochs/results.json',
        'Lean/Dynamics/ProtectedCharge.lean','Lean/Dynamics/ChargeConjugation.lean'})
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(import_closure(paths))}


def solve_many(matrix, rhs):
    """Rational elimination on native permutation blocks, at most six states."""
    n=len(matrix);width=len(rhs[0]);rows=[list(a)+list(b) for a,b in zip(matrix,rhs)]
    for col in range(n):
        pivot=next((i for i in range(col,n) if rows[i][col]),None)
        if pivot is None:raise ArithmeticError('transient source block is singular')
        rows[col],rows[pivot]=rows[pivot],rows[col]
        scale=rows[col][col];rows[col]=[v/scale for v in rows[col]]
        for i in range(n):
            if i==col:continue
            scale=rows[i][col]
            if scale:rows[i]=[v-scale*w for v,w in zip(rows[i],rows[col])]
    return [row[n:n+width] for row in rows]


@lru_cache(maxsize=1)
def hitting_data():
    source=classification();kernel=source['matrix'];classes=source['closed'];r=len(classes)
    groups=defaultdict(list)
    for state in native.STATES:
        groups[(sum(x*x for x in state[:3]),2*native.COST[state[3]],sum(state[:3]),state[3])].append(state)
    values={};blocks=[]
    for _,group in sorted(groups.items()):
        if all(s in source['class_of'] for s in group):
            for s in group:
                row=[F(0)]*(r+1);row[source['class_of'][s]]=F(1);values[s]=row
            continue
        index={s:i for i,s in enumerate(group)};n=len(group)
        a=[[F(i==j) for j in range(n)] for i in range(n)]
        b=[[F(0)]*r+[F(1)] for _ in group]
        for i,s in enumerate(group):
            for t,p in kernel[s].items():
                if t in index:a[i][index[t]]-=p
                else:
                    for j in range(r+1):b[i][j]+=p*values[t][j]
        values.update(zip(group,solve_many(a,b)));blocks.append((group,a))
    second={s:F(0) for s in source['class_of']}
    joint_time={s:[F(0)]*r for s in source['class_of']}
    for group,a in blocks:
        indices=set(group);b=[]
        for s in group:
            value=2*values[s][-1]-1
            value+=sum(p*second[t] for t,p in kernel[s].items() if t not in indices)
            times=[values[s][j]+sum(p*joint_time[t][j] for t,p in kernel[s].items() if t not in indices)
                   for j in range(r)]
            b.append([value]+times)
        for s,v in zip(group,solve_many(a,b)):
            second[s]=v[0];joint_time[s]=v[1:]
    for s,row in kernel.items():
        if sum(values[s][:-1])!=1 or any(v<0 for v in values[s]):
            raise ArithmeticError('invalid exact absorption law')
        if s not in source['class_of']:
            for j in range(r):
                if values[s][j]!=sum(p*values[t][j] for t,p in row.items()):
                    raise ArithmeticError('absorption harmonic equation fails')
                if joint_time[s][j]!=values[s][j]+sum(p*joint_time[t][j] for t,p in row.items()):
                    raise ArithmeticError('class-conditioned lifetime equation fails')
            if values[s][-1]!=1+sum(p*values[t][-1] for t,p in row.items()):
                raise ArithmeticError('first lifetime equation fails')
            if second[s]!=1+sum(p*(2*values[t][-1]+second[t]) for t,p in row.items()):
                raise ArithmeticError('second lifetime equation fails')
        if second[s]<values[s][-1]**2:raise ArithmeticError('negative lifetime variance')
        if sum(joint_time[s])!=values[s][-1]:raise ArithmeticError('class lifetimes do not sum to total lifetime')
    return dict(classes=classes,absorption={s:values[s][:-1] for s in values},
                mean={s:values[s][-1] for s in values},second=second,joint_time=joint_time,
                maximum_elimination_dimension=max(len(group) for group,_ in blocks))


def validate_distribution(distribution):
    if any(not native.valid_state(s) or type(v) not in (int,F) or v<0 for s,v in distribution.items()) or sum(distribution.values())!=1:
        raise ValueError('an exact normalized native preparation is required')


def terminal_law(distribution):
    validate_distribution(distribution);data=hitting_data()
    return [sum(p*data['absorption'][s][i] for s,p in distribution.items()) for i in range(len(data['classes']))]


def limiting_state(distribution):
    result={s:F(0) for s in native.STATES}
    for part,mass in zip(hitting_data()['classes'],terminal_law(distribution)):
        for s in part:result[s]=mass/len(part)
    return result


def statistics(distribution):
    validate_distribution(distribution)
    oriented=[sum(p for s,p in distribution.items() if s[3]==k) for k in range(6)]
    faces=[(oriented[k]+oriented[(-k)%6])/2 for k in range(6)]
    return dict(oriented_face_class_shares=oriented,two_face_class_shares=faces,
        oriented_1_minus_5=oriented[1]-oriented[5],
        unresolved_probability=sum(p for s,p in distribution.items() if not native.extinct_or_frozen(s)),
        mean_cost=sum(p*native.energy(s) for s,p in distribution.items()),
        total_probability=sum(distribution.values()))


def tv(left,right):
    return sum(abs(left.get(s,0)-right.get(s,0)) for s in set(left)|set(right))/2


def reset(distribution,target,probability):
    validate_distribution(distribution)
    if not native.valid_state(target) or type(probability) not in (int,F) or not 0<=probability<=1:
        raise ValueError('exact reset probability and native target are required')
    result={s:(1-probability)*p for s,p in distribution.items()};result[target]=result.get(target,0)+probability
    return result


def boundary_payoff(values):
    data=hitting_data()
    if len(values)!=len(data['classes']) or any(type(v) not in (int,F) for v in values):
        raise ValueError('one exact payoff per recurrent class required')
    return {s:sum(a*b for a,b in zip(row,values)) for s,row in data['absorption'].items()}


def trajectory(initial, *, boundary_reset=None):
    validate_distribution(initial);p=dict(initial);rows=[];branches=None
    if boundary_reset is None:limit=limiting_state(initial)
    else:
        cut,target,epsilon=boundary_reset
        at=native.evolve(initial,cut);after=reset(at,target,epsilon);limit=limiting_state(after)
    final_classes=[sum(limit[s] for s in part) for part in hitting_data()['classes']]
    for n in range(max(CUTS)+1):
        if boundary_reset is not None and n==boundary_reset[0]:
            epsilon=boundary_reset[2]
            branches={'no_reset':{s:(1-epsilon)*mass for s,mass in p.items()},
                      'reset':{boundary_reset[1]:epsilon}}
            p=reset(p,boundary_reset[1],boundary_reset[2])
        if n in CUTS:
            obs=statistics(p)
            branch_rows=[]
            if branches is not None:
                if any(sum(branch.get(s,0) for branch in branches.values())!=p.get(s,0) for s in native.STATES):
                    raise ArithmeticError('retained boundary records do not reconstruct the complete marginal')
                if sum(branches['no_reset'].values())!=1-boundary_reset[2] or sum(branches['reset'].values())!=boundary_reset[2]:
                    raise ArithmeticError('an original reset-record probability changed')
                branch_rows=[dict(record=label,original_probability=sum(branch.values()),
                    unnormalized_state_law=[branch.get(s,F(0)) for s in native.STATES])
                    for label,branch in branches.items()]
            rows.append(dict(attempts=n,**obs,terminal_class_probability=terminal_law(p),
                complete_state_law=[p.get(s,F(0)) for s in native.STATES],
                complete_boundary_record_state_laws=branch_rows,
                state_TV_to_limit=tv(p,limit)))
        if n<max(CUTS):
            p=native.evolve(p)
            if branches is not None:
                for label,branch in list(branches.items()):
                    advanced=defaultdict(F)
                    for s,mass in branch.items():
                        for t,weight in native.kernel()[s].items():advanced[t]+=mass*weight
                    branches[label]=dict(advanced)
    return dict(initial_state_law=[initial.get(s,F(0)) for s in native.STATES],
        boundary_reset=None if boundary_reset is None else dict(native_attempt=boundary_reset[0],
            target=list(boundary_reset[1]),reset_probability=boundary_reset[2],
            retained_public_records=['reset','no_reset'],additional_native_event=False,
            service='explicit supplied preparation/reset at the input boundary; not native dynamics'),
        terminal_class_probabilities=final_classes,limiting_class_shares=statistics(limit),
        limiting_state_law=[limit[s] for s in native.STATES],cuts=rows)


def universal_count_bound(epsilon):
    if type(epsilon) not in (int,F) or not 0<epsilon<1:raise ValueError('exact tolerance in(0,1) required')
    blocks=0;bound=F(1)
    while bound>epsilon:blocks+=1;bound/=25
    return dict(tolerance=epsilon,native_attempts=32*blocks,unresolved_bound=bound,
                physical_time=None,temperature=None)


def certificate():
    pins=bindings();data=hitting_data();base=native.initial_distribution()
    base_limit=limiting_state(base);base_statistics=statistics(base_limit)
    source_state=(0,1,5,0);asymmetric_state=(0,1,5,1);epsilon=F(1,10)
    base_classes=terminal_law(base);asymmetric_classes=terminal_law({asymmetric_state:F(1)})
    mixed=reset(base,asymmetric_state,epsilon);mixed_classes=terminal_law(mixed)
    if mixed_classes!=[(1-epsilon)*a+epsilon*b for a,b in zip(base_classes,asymmetric_classes)]:
        raise ArithmeticError('exact initial response is not affine')
    response=[b-a for a,b in zip(base_classes,asymmetric_classes)]
    cases={'native_uniform_zero_flux':trajectory(base),
           'fixed_load_zero_flux':trajectory({source_state:F(1)}),
           'fixed_load_with_declared_asymmetric_reset':trajectory({source_state:F(1)},
                boundary_reset=(8,asymmetric_state,epsilon))}
    lifetime=dict(mean=sum(p*data['mean'][s] for s,p in base.items()),
        second_moment=sum(p*data['second'][s] for s,p in base.items()),
        maximum_state_mean=max(data['mean'].values()),
        maximum_state_mean_attainers=[list(s) for s,t in data['mean'].items() if t==max(data['mean'].values())])
    lifetime['variance']=lifetime['second_moment']-lifetime['mean']**2
    lifetime['by_terminal_class']=[dict(id=i,probability=mass,
        joint_first_moment=sum(p*data['joint_time'][s][i] for s,p in base.items()),
        conditional_mean=(sum(p*data['joint_time'][s][i] for s,p in base.items())/mass if mass else None))
        for i,mass in enumerate(base_classes)]
    harmonic_payoff=boundary_payoff([F(int(part[0][3]==1)-int(part[0][3]==5)) for part in data['classes']])
    bound_cases=[universal_count_bound(F(1,10**k)) for k in (3,6,9,12)]
    checks=dict(all_absorption_rows_normalized=True,all_46_harmonic_equations_on_all_states=True,
        first_and_second_lifetime_backward_equations_exact=True,
        all_class_conditioned_lifetime_equations_exact=True,
        all_lifetime_variances_nonnegative=True,
        exact_persistent_share_matches_native_source=1-base_statistics['two_face_class_shares'][0]==F(1007,13824),
        no_oriented_asymmetry_for_even_native_preparation=base_statistics['oriented_1_minus_5']==0,
        asymmetric_boundary_has_nonzero_native_response=harmonic_payoff[asymmetric_state]!=0,
        affine_initial_distribution_response_exact=True,
        reset_record_joint_law_independently_propagated_and_reconstructs_marginal=True,
        boundary_readout_is_harmonic=all(harmonic_payoff[s]==sum(p*harmonic_payoff[t] for t,p in row.items())
            for s,row in native.kernel().items()),
        all_saved_full_state_laws_normalized=all(sum(c['complete_state_law'])==1 for case in cases.values() for c in case['cuts']))
    if not all(checks.values()):raise ArithmeticError('native freezeout response check failed')
    if pins!=bindings():raise ArithmeticError('source changed during absorption calculation')
    return dict(schema='oph.native_freezeout_response.v1',all_checks_pass=True,source_bindings=pins,
        state_order=[list(s) for s in native.STATES],
        recurrent_classes=[dict(id=i,states=[list(s) for s in part]) for i,part in enumerate(data['classes'])],
        all_state_hitting_data=[dict(state=list(s),class_probabilities=data['absorption'][s],
            expected_attempts=data['mean'][s],second_moment_attempts=data['second'][s],
            variance_attempts=data['second'][s]-data['mean'][s]**2,
            joint_class_first_lifetime_moment=data['joint_time'][s]) for s in native.STATES],
        complete_native_records_receipt='model/native_thermal_admission/results.json',
        maximum_exact_elimination_dimension=data['maximum_elimination_dimension'],
        native_uniform_initial_lifetime=lifetime,cases=cases,checks=checks,
        initial_response=dict(asymmetric_input_state=list(asymmetric_state),epsilon=epsilon,
            terminal_class_derivative=response,mixed_terminal_class_probabilities=mixed_classes,
            terminal_TV_response_per_epsilon=sum(abs(x) for x in response)/2,
            all_input_lipschitz_bound=1,sharp_bound_attained_by_distinct_absorbing_states=True),
        freezeout_criterion=dict(observable='future native face-class changes, after the final supplied reset',
            statement='TV of current and limiting face-class laws is at most probability outside the terminal set',
            origin_of_count='initial preparation, or the final supplied boundary reset; not total native count before that reset',
            uniform_bound='25^(-floor(n/32))',source_count_examples=bound_cases,
            full_state_bound='TV(p_n,p_infinity) <= 25^(-floor(floor(n/2)/32)) + (2/3)*2^(-(n-floor(n/2)))',
            calibrated_expansion_clock=False,cosmological_temperature=False),
        acceptance_flags=dict(exact_native_all_class_absorption_and_lifetimes=True,
            boundary_asymmetry_is_supplied=True,spontaneous_baryogenesis_claimed=False,
            native_cosmological_freezeout_or_relic_abundances=False))


def main():
    parser=argparse.ArgumentParser();m=parser.add_mutually_exclusive_group(required=True)
    m.add_argument('--write',action='store_true');m.add_argument('--check',action='store_true');args=parser.parse_args()
    result=json.loads(json.dumps(certificate(),default=encode,allow_nan=False));path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    elif json.loads(path.read_text())!=result:raise ArithmeticError('native freezeout receipt differs')
    print('PASS:46-class exact absorption, complete class histories, lifetimes and boundary response')


if __name__=='__main__':main()
