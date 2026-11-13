"""Outward Taylor integration of the radiation-plus-event-horizon branch."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
import sympy as s
from model.finite_area_response.response import Interval,exponential
from model.openwave_universal_energy_geometry.source import root_interval
from model.openwave_horizon_capacity.source import log_interval

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
INPUTS=json.loads((HERE/'inputs.json').read_text())


def round_interval(x,bits=64):
    D=2**bits
    return Interval(F((x.lo*D).__floor__(),D),F((x.hi*D).__ceil__(),D))


@lru_cache(maxsize=20000)
def exp_point(x):
    x=F(x);scale=0
    while abs(x)>2:x/=2;scale+=1
    value=exponential(x)
    for _ in range(scale):value=value*value
    return value


def sigmoid(w):
    e=exp_point(F(w))
    return e/(1+e)


@lru_cache(maxsize=20000)
def radiation_ratio(t,om,orr):
    z1=exp_point(t)
    return orr*z1/(om+orr*z1)


def flow(t,w,om,orr):
    x=sigmoid(w);r=radiation_ratio(t,om,orr);root=root_interval(x)
    f=Interval.point(-1)-r-2*root
    fw=-root*(Interval.point(1)-x)
    ft=-r*(Interval.point(1)-r)
    return f,ft+fw*f


def taylor_step(t,bracket,h,om,orr):
    """Monotone second-order Taylor map plus a proved |w'''|<8 bound."""
    h=F(h)
    if not 0<h<=F(1,2):raise ValueError('positive step at most1/2 required')
    endpoints=[]
    for w in (bracket.lo,bracket.hi):
        f,second=flow(t,w,om,orr)
        endpoints.append(Interval.point(w)+h*f+h*h*second/2)
    remainder=F(4,3)*h**3
    return round_interval(Interval(endpoints[0].lo-remainder,endpoints[1].hi+remainder))


def derivative_certificate():
    # r in[0,1], x=sqrt(Omega) in[0,1]. All bounds apply globally.
    # x(1-x^2)<=2/5: maximum2/(3sqrt3)<2/5, since4/27<4/25.
    assert F(4,27)<F(4,25)
    f=F(4);fw=F(2,5);ft=F(1,4);fww=F(2,5);ftt=F(1,4)
    second=ft+fw*f
    third=ftt+fww*f*f+fw*second
    assert third==F(739,100) and third<8
    # d/dw Taylor2 >=1-(2/5)h-(4/5)h^2 >0 for h<=1/2.
    minimum=1-fw*F(1,2)-F(4,5)*F(1,2)**2
    assert minimum==F(3,5)
    return {'logit_flow':'w_t=-(1+r(t)+2 sqrt(sigmoid(w)))',
            'r':'Omega_r0*exp(t)/(Omega_m0+Omega_r0*exp(t))',
            'absolute_third_derivative_upper_exact':str(third),
            'used_third_derivative_budget_exact':'8',
            'taylor_map_derivative_lower_for_step_at_most_half':str(minimum),
            'local_remainder':'(4/3)*h^3',
            'interval_rounding':'outward64binaryplaces plus original30decimalinterval arithmetic'}


def integrate(om,orr,redshifts,step=F(1,128)):
    om,orr,step=map(F,(om,orr,step));x0=1-om-orr
    if not 0<orr<om<1 or not 0<x0<1:raise ValueError('positive flat density fractions required')
    if not redshifts or list(redshifts)!=sorted(set(redshifts)) or redshifts[0]!=0:raise ValueError('ordered distinct redshifts starting at zero required')
    if not 0<step<=F(1,2):raise ValueError('supported positive integration step required')
    initial=log_interval(x0/(1-x0));w=Interval(*initial);t=F(0);rows=[];steps=0
    for z in redshifts:
        tlo,thi=log_interval(F(1)+z)
        while t<tlo:
            h=min(step,tlo-t);w=taylor_step(t,w,h,om,orr);t+=h;steps+=1
        # f lies in[-4,-1]. This also encloses the exact irrational endpoint.
        endpoint=Interval(w.lo-4*(thi-tlo),w.hi)
        a,b=sigmoid(endpoint.lo),sigmoid(endpoint.hi)
        omega=Interval(a.lo,b.hi)
        H_ratio=root_interval(Interval.point(1)/(Interval.point(1)-omega))
        rows.append({'redshift':z,'Omega_dark_interval':omega.encoded(),
            'H_over_same_matter_radiation_background_interval':H_ratio.encoded(),
            'logit_enclosure_width_exact':str(endpoint.hi-endpoint.lo)})
    return rows,steps


def run():
    parent=json.loads((ROOT/'model/openwave_horizon_capacity/results.json').read_text())
    om=F(parent['anchors']['Omega_m']);h0=F(parent['anchors']['H0_km_s_Mpc'])
    orr=F(INPUTS['radiation_fraction_today']);x0=1-om-orr
    cert=derivative_certificate()
    rows,steps=integrate(om,orr,INPUTS['redshifts'],F(INPUTS['maximum_log_redshift_step']))
    z=s.symbols('z',nonnegative=True);m,r=s.symbols('m r',positive=True);x=s.symbols('x',positive=True)
    wdark=-s.Rational(1,3)-s.Rational(2,3)*s.sqrt(x)
    ratio=r*(1+z)/(m+r*(1+z));omega_r=(1-x)*ratio
    flow_a=s.expand(3*x*(omega_r/3+wdark*x-wdark))
    expected=x*(1-x)*(1+ratio+2*s.sqrt(x))
    assert s.simplify(flow_a-expected)==0
    fixed=[]
    for redshift in INPUTS['redshifts']:
        a1=F(1)+redshift
        fixed.append({'redshift':redshift,'Omega_dark_exact':str(x0/(om*a1**3+orr*a1**4+x0))})
    factor=F(INPUTS['radiation_neutrino_energy_factor']);neff=F(INPUTS['nominal_relativistic_neutrino_number'])
    delta_neff=x0/(1-x0)*(1+factor*neff)/factor
    today_boost=root_interval(Interval.point(1)/(1-x0))
    return {'schema':'oph.radiation_capacity_history.v1','verdict':'PASS_VALIDATED_RADIATION_CAPACITY_BRANCHES',
       'inputs':{'H0_km_s_Mpc':str(h0),'Omega_m0':str(om),'Omega_r0':str(orr),'Omega_dark0':str(x0),
          'ancestry':INPUTS['ancestry'],'new_empirical_fits':0},
       'conservation_identity':'dOmega_dark/dlog(a)=Omega_dark*(1-Omega_dark)*(1+r(a)+2sqrt(Omega_dark)), gamma=1',
       'equation_of_state':'w_dark=-1/3-2sqrt(Omega_dark)/3',
       'integration_certificate':cert,'accepted_taylor_steps':steps,
       'future_event_branch':rows,'fixed_capacity_branch':fixed,
       'early_scaling':{'future_event_radiation':'Omega_dark proportional to a^2 asymptotically',
          'future_event_matter':'Omega_dark proportional to a asymptotically',
          'fixed_capacity_radiation':'Omega_dark proportional to a^4 asymptotically',
          'hubble_area_conserved':'Omega_dark remains the supplied current constant; w_dark follows the mixture of conserved matter and radiation'},
       'hubble_capacity_control':{'constant_dark_fraction':str(x0),
          'H_over_same_matter_radiation_background_interval':today_boost.encoded(),
          'radiation_expansion_equivalent_extra_Neff_exact':str(delta_neff),
          'interpretation':'Expansion-equivalent radiation density only, not a CMB perturbation or nucleosynthesis likelihood.'},
       'scope':{'radiation_and_matter_conservation_derived':True,'horizon_selection_derived':False,
          'cosmological_amplitude_predicted':False,'standard_model_thermal_history_derived':False,
          'early_universe_abundances_calculated':False,'independent_observational_validation':False,
          'maximum_redshift':max(INPUTS['redshifts']),
          'radiation_domain':INPUTS['radiation_scope']}}
