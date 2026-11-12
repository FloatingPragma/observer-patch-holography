"""Explicit native atomic spectral replay of the full polynomial filter."""
import mpmath as mp
import sympy as s
from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import target
from .isolation import data
from .filter import design


def number(value):
    value=s.Rational(value)
    return mp.mpf(int(value.p))/int(value.q)


def value(energy,parameters=None):
    f=design() if parameters is None else parameters
    z=(2*energy-number(f['upper'])-number(f['c']))/(number(f['upper'])-number(f['c']))
    m=f['degree'];denominator=mp.cosh(m*mp.acosh(number(f['za'])))
    if z<=-1:numerator=(-1)**m*mp.cosh(m*mp.acosh(-z))
    elif z>=1:numerator=mp.cosh(m*mp.acosh(z))
    else:numerator=mp.cos(m*mp.acos(z))
    base=(-1)**m*numerator/denominator
    return 1-(1-base)**f['flattening_power']


def certificate():
    with mp.workdps(130):
        d=data();f=design();x=s.Symbol('x');h=target(2,0)['H']
        polynomial=s.Poly(h.charpoly(x).as_expr(),x)
        e0=sum(map(number,d['e0']))/2
        roots=polynomial.intervals(eps=s.Rational(1,10**110))
        rows=[];weights=[];maximum=mp.mpf(0)
        for index,(interval,rank) in enumerate(roots):
            local=sum(map(number,interval))/2;energy=7*e0+local
            amplitude=value(energy,f);expected=int(index==0)
            error=abs(amplitude-expected);maximum=max(maximum,error)
            assert error<number(f['epsilon'])
            weight=rank*amplitude**2/h.rows;weights.append(weight)
            rows.append({'local_energy_decimal':mp.nstr(local,40),'native_spectral_rank':rank,
                'filter_amplitude_decimal':mp.nstr(amplitude,45),
                'joint_accept_and_energy_probability_decimal':mp.nstr(weight,45),
                'is_selected_neutral_eigenline':bool(expected)})
        success=sum(weights);failure=1-success
        assert sum(rank for _,rank in roots)==h.rows==16 and failure>0
        assert abs(rows[0]['native_spectral_rank']-1)==0
        return {'source':'actual 16-row local physical neutral block, seven exact spectator grounds and blank links, lambda=0 control',
            'nonzero_lambda_full_carrier_filter_numerically_materialized':False,
            'full_interacting_validity_supplied_by_operator_isolation_and_filter_bounds':True,
            'all_native_block_spectral_ranks_retained':16,
            'energy_records':rows,'complete_failure_probability_decimal':mp.nstr(failure,45),
            'complete_acceptance_probability_decimal':mp.nstr(success,45),
            'maximum_native_block_filter_error_decimal':mp.nstr(maximum,25),
            'comparison_precision_decimal_digits':130,
            'floating_replay_is_not_the_rational_filter_error_certificate':True,
            'complete_probability_sum':'1'}
