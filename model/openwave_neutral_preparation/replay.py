"""Actual local sixteen-mode polynomial replay, retaining every complement."""
import mpmath as mp
import sympy as s
from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import target
from model.openwave_neutral_stationary.isolation import data
from .filtering import parameters


def number(x):
    x=s.Rational(x)
    return mp.mpf(int(x.p))/int(x.q)


def value(energy,p):
    u=number(p['upper']);c=number(p['c'])
    x=(2*energy-u-c)/(u-c)
    f=(-1)**p['degree']*mp.chebyt(p['degree'],x)/mp.chebyt(p['degree'],number(p['za']))
    return 1-(1-f)**p['flattening_power']


def certificate():
    p=parameters();d=data();h=target(2,0)['H'];x=s.Symbol('x')
    # Exact factors and multiplicities come from the full native block.
    factors=s.factor_list(h.charpoly(x).as_expr(),x)[1]
    rows=[]
    with mp.workdps(150):
        e0=sum(number(z) for z in d['e0'])/2
        en=sum(number(z) for z in d['en'])/2
        for factor,multiplicity in factors:
            for interval,count in s.Poly(factor,x).intervals(eps=s.Rational(1,10**120)):
                local=sum(number(z) for z in interval)/2
                energy=7*e0+local;amplitude=value(energy,p)
                selected=abs(local-en)<mp.mpf('1e-50')
                if selected:assert 1-number(p['epsilon'])<=amplitude<=1
                else:assert abs(amplitude)<=number(p['epsilon'])
                success=amplitude**2;failure=1-success
                assert failure>=0 and abs(success+failure-1)<mp.mpf('1e-140')
                rows.append({'local_root_interval':[str(z) for z in interval],
                    'multiplicity':multiplicity*count,'selected_neutral_line':bool(selected),
                    'filter_amplitude_decimal':mp.nstr(amplitude,50),
                    'complete_success_probability_decimal':mp.nstr(success,50),
                    'complete_failure_probability_decimal':mp.nstr(failure,50)})
    assert sum(row['multiplicity'] for row in rows)==16
    assert sum(row['multiplicity'] for row in rows if row['selected_neutral_line'])==1
    return {'full_native_block_dimension':16,'all_modes':rows,
        'complete_accept_reject_law_per_mode':'M(E)^2+(1-M(E)^2)=1',
        'numeric_precision_digits':150,'full_many_cell_H_materialized':False,
        'all_source_spectrum_filter_bound':'rational Chebyshev certificate; local replay is an additional native finite block control'}
