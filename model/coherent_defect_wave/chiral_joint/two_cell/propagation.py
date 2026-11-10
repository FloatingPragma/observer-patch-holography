"""Time-controlled transport of a gauge-singlet quantum pair preparation."""
from math import factorial,isqrt
from collections import defaultdict
import sympy as s
from . import hamiltonian as opened
from .basis import local_current_column
from .symmetries import (apply_column,periodic_column,periodic_current_column,
                         periodic_translation,norm2)
from ..model import clean


def sqrt_upper(value,bits=64):
    value=s.Rational(value)
    if value<0:raise ValueError('nonnegative squared norm')
    scale=1<<bits;n=isqrt(int(value.p)*scale*scale//int(value.q))
    if n*n*int(value.q)<int(value.p)*scale*scale:n+=1
    return s.Rational(n,scale)


def pair_preparation(site=1,periodic=False):
    """Actual local quantum Yukawa action on the empty state, normalized."""
    if site not in (0,1):raise ValueError('two source sites')
    local=opened.local_yukawa(0,0)
    assert sum(s.conjugate(v)*v for v in local.values())==4
    state={}
    for (mask,scalar),value in local.items():
        fullmask=mask<<(6*site)
        hl,hr=(scalar,0) if site==0 else (0,scalar)
        key=(fullmask,0,0,hl,hr) if periodic else (fullmask,0,hl,hr)
        state[key]=value/2
    assert norm2(state)==1
    # A complete local preparation unitary is exp[-i pi H_Y,site/4].
    # Its vacuum cyclic subspace has H_Y^2|empty>=4|empty>.
    square=defaultdict(lambda:s.Integer(0))
    for (mask,scalar),value in local.items():
        for key,coefficient in opened.local_yukawa(mask,scalar).items():
            square[key]+=value*coefficient
    assert clean(square)=={(0,0):s.Integer(4)}
    return state


def arrival(state,site):
    """Complete odd-occupation outcome at the initially empty receiver."""
    if site not in (0,1):raise ValueError('receiver site')
    return {key:value for key,value in state.items()
            if ((key[0]>>(6*site))&63).bit_count()%2}


def local_only_column(key):
    """Open source control with the hopping term removed, all other terms kept."""
    mask,link,hl,hr=opened.validate_key(key)
    result=dict(opened.column(key))
    for (target,newlink),value in opened.hopping(mask,link).items():
        out=(target,newlink,hl,hr);result[out]=result.get(out,0)-value
    return clean(result)


def taylor_state(state,action,time,order):
    time=s.sympify(time)
    if time.has(s.Float) or time.is_real is not True:raise ValueError('exact real source time')
    if type(order) is not int or order<0:raise ValueError('nonnegative order')
    output=defaultdict(lambda:s.Integer(0));power=dict(state)
    for n in range(order+1):
        coefficient=(-s.I*time)**n/s.factorial(n)
        for key,value in power.items():output[key]+=coefficient*value
        if n<order:power=apply_column(power,action)
    return clean(output)


def probability_interval(polynomial,receiver,remainder,bits=36):
    p=norm2(arrival(polynomial,receiver));upper_norm=sqrt_upper(p)
    low=max(s.Integer(0),p-2*remainder*upper_norm)
    high=min(s.Integer(1),p+2*remainder*upper_norm+remainder**2)
    scale=1<<bits
    outward_low=s.Rational(s.floor(low*scale),scale)
    outward_high=s.Rational(s.ceiling(high*scale),scale)
    assert 0<=outward_low<=low<=high<=outward_high<=1
    return p,outward_low,outward_high


def transport_control(periodic=False):
    state=pair_preparation(1,periodic)
    action=periodic_column if periodic else opened.column
    gauss=periodic_current_column if periodic else local_current_column
    for site in range(2):
        for current in range(12):
            assert not apply_column(state,lambda k:gauss(k,site,current))
    assert not arrival(state,0)
    first=apply_column(state,action)
    c=norm2(arrival(first,0))
    assert c==(s.Rational(3,16) if periodic else s.Rational(3,4))
    # Twenty half-unit open hopping terms, or forty quarter-unit periodic
    # terms, contribute ten. Two local Yukawa operators contribute sixteen.
    # Scalar/link projectors contribute three open, four periodic.
    bound=s.Integer(30 if periodic else 29)
    time=s.Rational(1,20);order=8
    polynomial=taylor_state(state,action,time,order)
    remainder=(time*bound)**(order+1)/s.factorial(order+1)
    p,lo,hi=probability_interval(polynomial,0,remainder)
    assert lo>0
    for site in range(2):
        for current in range(12):
            assert not apply_column(polynomial,lambda k:gauss(k,site,current))
    result={'source':'periodic-two-oriented-links' if periodic else 'uniform-open-two-cell',
            'full_source_dimension':3686400 if periodic else 614400,
            'preparation':'exp(-i pi H_Y,right/4) on empty physical source, up to global phase',
            'preparation_driven':True,'preparation_work_calculated':False,
            'prepared_state_is_physical_ground':False,
            'prepared_state_L_twoSz_parity':[0,0,'even'],
            'all_24_local_Gauss_actions_on_preparation_and_polynomial_zero':True,
            'detector':'odd matter occupation at initially empty left cell',
            'same_idle_H_preparation_arrival_initial_probability':'0',
            'exact_first_derivative_arrival_norm_squared':str(c),
            'source_time':str(time),'Taylor_order':order,'full_H_norm_bound':str(bound),
            'unitary_operator_norm_remainder_bound':str(remainder),
            'polynomial_source_basis_support':len(polynomial),
            'polynomial_arrival_probability':str(p),
            'complete_arrival_yes_probability_interval':[str(lo),str(hi)],
            'complete_arrival_no_probability_interval':[str(1-hi),str(1-lo)],
            'all_outcomes_retained':True,'finite_quantum_transport_observed':True,
            'physical_particle_rest_mass_or_scattering_claim':False}
    if periodic:
        left=pair_preparation(0,True)
        assert apply_column(state,periodic_translation)==left
        translated=apply_column(polynomial,periodic_translation)
        assert norm2(arrival(translated,1))==p
        # The exact H-T intertwiner gives equality for the whole unitary,
        # not only the displayed Taylor polynomial.
        result['opposite_direction_complete_transport_laws_equal_by_exact_translation']=True
        result['periodic_reference_is_open_source_ground']=False
    else:
        # Removing propagation leaves left occupation even forever: each
        # remaining local term changes it by zero or two, for all columns.
        assert not arrival(apply_column(state,local_only_column),0)
        for mask in range(64):
            for scalar in range(5):
                for (target,_),value in opened.local_yukawa(mask,scalar).items():
                    assert value and (target.bit_count()-mask.bit_count())%2==0
        assert not arrival(taylor_state(state,local_only_column,time,order),0)
        result['no_hopping_control_arrival_probability_at_all_times']='0'
    return result


def certificate():
    return {'open_transport':transport_control(False),'periodic_transport':transport_control(True)}


def certificate_data():
    return certificate()
