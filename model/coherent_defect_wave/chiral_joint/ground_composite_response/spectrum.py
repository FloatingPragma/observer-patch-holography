"""Exact projector moments and every ground-composite energy probability."""
from functools import lru_cache
import sympy as s

from .model import moments
from ..family_cell.quotient import add,mul,divide,evaluate,outward
from ..local_reference.spectrum import certificate as reference_spectrum


def projector_coefficients(polynomial,factor):
    """p(H)/(H-E)/p'(E) for a semisimple finite matrix, modulo its root factor."""
    x=polynomial.gen
    return tuple(s.Poly(sum(polynomial.nth(j)*x**(j-i-1)
        for j in range(i+1,polynomial.degree()+1)),x).rem(factor)
        for i in range(polynomial.degree()))


def numerator(polynomial,target_polynomial,ground_factor,target_factor,cross_moments):
    source=projector_coefficients(polynomial,ground_factor)
    target=projector_coefficients(target_polynomial,target_factor)
    values=[[s.Rational(0) for _ in range(target_factor.degree())]
            for _ in range(ground_factor.degree())]
    for j,t in enumerate(target):
        for i,a in enumerate(source):
            moment=cross_moments[j][i]
            if not moment:continue
            for u in range(ground_factor.degree()):
                if not a.nth(u):continue
                for v in range(target_factor.degree()):values[u][v]+=moment*a.nth(u)*t.nth(v)
    return tuple(tuple(row) for row in values)


def evaluate_bivariate(coefficients,source_interval,target_interval):
    result=(s.Rational(0),s.Rational(0))
    for row in reversed(coefficients):
        evaluated=(s.Rational(0),s.Rational(0))
        for coefficient in reversed(row):evaluated=add(mul(evaluated,target_interval),(coefficient,coefficient))
        result=add(mul(result,source_interval),evaluated)
    return result


@lru_cache(maxsize=None)
def response(kind,spin):
    data,p,q,cross=moments(kind,spin);x=p.gen
    ground_factor=s.Poly((x*(x-1))**3-313*(x*(x-1))**2+31708*x*(x-1)-1041468,x)
    assert p.rem(ground_factor).is_zero
    ground_roots=s.polys.polytools.intervals(ground_factor,eps=s.Rational(1,10**110))
    assert len(ground_roots)==6 and all(m==1 for _,m in ground_roots)
    ground_interval=tuple(map(s.Rational,ground_roots[0][0]))
    declared=reference_spectrum()['global_physical_ground_energy_interval']
    assert s.Rational(declared[0])<=ground_interval[0]<=ground_interval[1]<=s.Rational(declared[1])
    source_denominator=p.diff().rem(ground_factor)
    source_denominator_interval=evaluate(source_denominator,ground_interval)
    source_coefficients=projector_coefficients(p,ground_factor)
    source_power=s.eye(32);source_trace=0
    for coefficient in source_coefficients:
        source_trace+=coefficient.as_expr()*s.trace(source_power)
        source_power=data['source_H']*source_power
    assert s.Poly(source_trace,x).rem(ground_factor)==2*source_denominator
    total_numerator=s.Poly(sum(c.as_expr()*moment for c,moment in zip(source_coefficients,cross[0])),x).rem(ground_factor)
    direct_success=mul(divide(evaluate(total_numerator,ground_interval),source_denominator_interval),
                       (s.Rational(1,32),)*2)
    full_target=s.Poly(data['target_H'].charpoly(x).as_expr(),x)
    target_ranks={str(factor.monic().as_expr()):rank for factor,rank in s.factor_list(full_target)[1]}
    target_traces=[];target_power=s.eye(26)
    for _ in range(q.degree()):
        target_traces.append(s.trace(target_power));target_power=data['target_H']*target_power
    records=[];success=(s.Rational(0),s.Rational(0));algebraic=[]
    for factor,_ in s.factor_list(q)[1]:
        factor=factor.monic()
        coefficients=numerator(p,q,ground_factor,factor,cross)
        target_denominator=q.diff().rem(factor)
        rank=target_ranks[str(factor.as_expr())]
        target_trace=s.Poly(sum(c.as_expr()*t for c,t in zip(projector_coefficients(q,factor),target_traces)),x).rem(factor)
        assert target_trace==rank*target_denominator
        algebraic.append({'target_energy_factor':str(factor.as_expr()),
            'numerator_coefficients_source_then_target':[[str(c) for c in row] for row in coefficients],
            'source_denominator':str(source_denominator.as_expr()),
            'target_denominator':str(target_denominator.as_expr()),
            'joint_probability_scale':'1/32'})
        zero=not any(c for row in coefficients for c in row)
        for root_index,(interval,multiplicity) in enumerate(s.polys.polytools.intervals(factor,eps=s.Rational(1,10**110))):
            assert multiplicity==1
            interval=tuple(map(s.Rational,interval))
            if zero:weight=(s.Rational(0),s.Rational(0))
            else:
                den=mul(source_denominator_interval,evaluate(target_denominator,interval))
                weight=mul(divide(evaluate_bivariate(coefficients,ground_interval,interval),den),
                           (s.Rational(1,32),)*2)
                assert 0<weight[0]<=weight[1]<s.Rational(1,4),(kind,spin,factor,root_index,weight)
            success=add(success,weight)
            excitation=(interval[0]-ground_interval[1],interval[1]-ground_interval[0])
            assert excitation[0]>0
            records.append({'target_energy_factor':str(factor.as_expr()),'target_root_index':root_index,
                'target_energy_interval':outward(interval,48),'excitation_energy_interval':outward(excitation,48),
                'joint_probability_interval':outward(weight,48),'exactly_zero_weight':zero,
                'complete_target_spectral_rank':rank})
    records.sort(key=lambda row:s.Rational(row['target_energy_interval'][0]))
    assert len(records)==9 and success[1]<s.Rational(1,8)
    assert sum(row['complete_target_spectral_rank'] for row in records)==26
    assert max(success[0],direct_success[0])<=min(success[1],direct_success[1])
    failure=(s.Rational(1,4)-direct_success[1],s.Rational(1,4)-direct_success[0])
    assert failure[0]>0
    return {'kind':kind,'probe_spin':spin,'setting_probability':'1/4',
        'reference_density':'P_ground/2','success_operator':'J_kind,spin/2',
        'target_raw_L':data['target_raw_L'],'complete_physical_target_dimension':26,
        'complete_target_energy_records':records,'joint_success_probability_interval':outward(direct_success,48),
        'joint_failure_probability_interval':outward(failure,48),
        'outside_complete_target_record':{'joint_probability_interval':['0','0'],'retained':True},
        'exact_source_minimal_polynomial':str(p.as_expr()),'exact_target_minimal_polynomial':str(q.as_expr()),
        'exact_cross_moments_target_rows_source_columns':[[str(v) for v in row] for row in cross],
        'exact_projector_weight_functions':algebraic,
        'source_projector_trace':2,'source_and_target_full_metric_intertwiners':True,
        'failure_operator_applied_on_full_CAR_carrier':True,
        'conditioned_success_renormalization_applied':False,'numerical_eigensolver_used':False}


def certificate():
    settings=[response(kind,spin) for kind in ('remove','add') for spin in (0,1)]
    baseline=moments('remove',0)[1:]
    for setting in settings:
        assert moments(setting['kind'],setting['probe_spin'])[1:]==baseline
        assert setting['complete_target_energy_records']==settings[0]['complete_target_energy_records']
    all_intervals=[row['joint_probability_interval'] for setting in settings
                   for row in setting['complete_target_energy_records']]
    all_intervals.extend(setting['joint_failure_probability_interval'] for setting in settings)
    assert sum(s.Rational(row[0]) for row in all_intervals)<=1<=sum(s.Rational(row[1]) for row in all_intervals)
    return {'recorded_settings':settings,'complete_joint_record_count':44,
            'positive_probability_records':sum(1 for row in all_intervals if s.Rational(row[0])>0),
            'exact_setting_sum_and_failure_completeness':True,
            'probe_spin_independence_on_rotation_invariant_reference':True,
            'exact_charge_conjugate_addition_removal_weights_equal':True,
            'all_success_excitation_energies_strictly_positive':True}
