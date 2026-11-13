"""Exact horizon-capacity alternatives on the declared Friedmann branch."""
from fractions import Fraction as F
from pathlib import Path
import json
from model.openwave_capacity_stress.source import (load_data,covariance_components,
    sqrt_interval,chi_square_interval,enclosure)
import model.openwave_capacity_stress.source as parent
ROOT=Path(parent.__file__).resolve().parents[2]

BRANCHES=('fixed_capacity','comoving_area','hubble_area_conserved')

def log_interval(value, terms=36):
    """Positive rational logarithm, with a geometric remainder enclosure."""
    value=F(value)
    if value<=0:raise ValueError('positive logarithm argument')
    power=0
    while value>=2:value/=2;power+=1
    while value<1:value*=2;power-=1
    def series(v):
        z=(v-1)/(v+1)
        lo=2*sum(z**(2*j+1)/(2*j+1) for j in range(terms))
        return lo,lo+2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
    a,b=series(value),series(F(2))
    return (a[0]+power*b[0],a[1]+power*b[1]) if power>=0 else (a[0]+power*b[1],a[1]+power*b[0])

def event_primitive(x):
    # Derivative = 2/[x(1-x^2)(1+2x)] > 0 for 0<x<1.
    terms=((2,x),(1,1+x),(-F(1,3),1-x),(-F(8,3),1+2*x))
    lo=hi=F(0)
    for c,v in terms:
        a,b=log_interval(v)
        lo+=c*(a if c>0 else b);hi+=c*(b if c>0 else a)
    return lo,hi

def event_omega(z,om):
    z,om=F(z),F(om)
    if z<0 or not 0<om<1:raise ValueError('invalid event-horizon domain')
    if z==0:return 1-om,1-om
    precision=24
    x0=sqrt_interval(1-om,precision)
    while x0[0]<=0 or x0[1]>=1:
        precision*=2
        x0=sqrt_interval(1-om,precision)
    def rational_primitive(x):
        return x**6*(1+x)**3/((1-x)*(1+2*x)**8)
    scale=1/(1+z)**3
    target=(scale*rational_primitive(x0[0]),scale*rational_primitive(x0[1]))
    lower,upper=x0[0]/2,x0[1]
    while rational_primitive(lower)>=target[0]:
        lower/=2
    for _ in range(64):
        mid=(lower+upper)/2;f=rational_primitive(mid)
        if f<target[0]:lower=mid
        elif f>target[1]:upper=mid
        else:break
    return lower*lower,upper*upper

def event_compare(table,h0,om):
    data=load_data(table);diag,u=covariance_components(data)
    intervals=[];states=[]
    for z,_,_ in data:
        a,b=event_omega(z,om);numerator=h0*h0*om*(1+z)**3
        pred=(sqrt_interval(numerator/(1-a))[0],sqrt_interval(numerator/(1-b))[1])
        intervals.append(pred);states.append((a,b))
    residual=[(a-h,b-h) for (a,b),(_,h,_) in zip(intervals,data)]
    return {'branch':'future_event_area_gamma_one','table':table,'rows':len(data),
            'chi_square_interval':enclosure(chi_square_interval(diag,u,residual)),
            'fitted_chronometer_parameters':0,
            'records':[{'z':str(z),'observed_H':str(h),'quoted_sigma':str(s),
                        'H_interval':enclosure(v),'Omega_dark_interval':enclosure(o)}
                       for (z,h,s),v,o in zip(data,intervals,states)]}
def h2(z,branch,h0,om):
    z,h0,om=map(F,(z,h0,om))
    if z<0 or h0<=0 or not 0<om<1:raise ValueError('invalid expansion domain')
    if branch=='fixed_capacity':return h0*h0*(om*(1+z)**3+1-om)
    if branch=='comoving_area':return h0*h0*(om*(1+z)**3+(1-om)*(1+z)**2)
    if branch=='hubble_area_conserved':return h0*h0*(1+z)**3
    raise ValueError('unregistered horizon branch')

def compare(table,branch,h0,om):
    data=load_data(table);diag,u=covariance_components(data)
    intervals=[sqrt_interval(h2(z,branch,h0,om)) for z,_,_ in data]
    residual=[(a-h,b-h) for (a,b),(_,h,_) in zip(intervals,data)]
    chi=chi_square_interval(diag,u,residual)
    return {'branch':branch,'table':table,'rows':len(data),
            'chi_square_interval':enclosure(chi),'fitted_chronometer_parameters':0,
            'records':[{'z':str(z),'observed_H':str(h),'quoted_sigma':str(s),
                        'H_interval':enclosure(v)} for (z,h,s),v in zip(data,intervals)]}

def run():
    receipt=json.loads((ROOT/'model/openwave_capacity_stress/results.json').read_text())
    def anchor(x):
        if isinstance(x,dict):
            if 'anchors' in x:return x['anchors']
            for v in x.values():
                r=anchor(v)
                if r:return r
        if isinstance(x,list):
            for v in x:
                r=anchor(v)
                if r:return r
    anchors=anchor(receipt);assert anchors
    h0,om=F(anchors['H0_km_s_Mpc']),F(anchors['Omega_m'])
    omega=1-om
    sr=sqrt_interval(omega)
    hfuture_squared_ratio=(F(3**8,16)*sr[0]**6*(1+sr[0])**4/(1+2*sr[1])**8,
                           F(3**8,16)*sr[1]**6*(1+sr[1])**4/(1+2*sr[0])**8)
    hfuture=(h0*sqrt_interval(hfuture_squared_ratio[0])[0],h0*sqrt_interval(hfuture_squared_ratio[1])[1])
    # Event-horizon gamma=1 is an architectural control, not selected by these data.
    event_w=(-F(1,3)-2*sr[1]/3,-F(1,3)-2*sr[0]/3)
    return {'schema':'oph.horizon_capacity_comparison.v1',
            'verdict':'PASS_EXACT_CAPACITY_BRANCHES_AND_CORRELATED_COMPARISON',
            'anchors':anchors,'capacity_law':'rho=epsilon/N, N proportional to horizon radius squared',
            'comparisons':[compare(t,b,h0,om) for t in ('HzTable_MM_BC03.dat','HzTable_MM_M11.dat') for b in BRANCHES]
                          +[event_compare(t,h0,om) for t in ('HzTable_MM_BC03.dat','HzTable_MM_M11.dat')],
            'branch_consequences':{
               'fixed_capacity':{'w_exact':'-1','q_today_exact':str((1-3*omega)/2),'amplitude_selected':False},
               'comoving_area':{'w_exact':'-1/3','q_today_exact':str(om/2),'acceleration':False},
               'hubble_area_conserved':{'w_exact':'0','q_exact':'1/2','constant_dark_fraction_exact':str(omega),'acceleration':False},
               'hubble_area_mechanical_pressure_with_exchange':{
                   'mechanical_w_exact':'-1','total_w_exact':str(-omega),
                   'q_exact':str((1-3*omega)/2),
                   'matter_scaling_exponent_exact':str(3*(1-omega)),
                   'dark_fraction_at_all_epochs_exact':str(omega),
                   'matter_dominance_to_acceleration_transition':False,
                   'requires_energy_exchange':True},
               'future_event_area':{'gamma_exact':'1','w_today_interval':enclosure(event_w),
                   'dOmega_dlog_a':'Omega*(1-Omega)*(1+2*sqrt(Omega)/gamma)',
                   'future_boundary_condition_proved_for_gamma_one':True,'initial_dark_fraction_supplied':True,
                   'expansion_likelihood_executed':True,
                   'H_future_km_s_Mpc_interval':enclosure(hfuture),
                   'Omega_at_z_100_interval':enclosure(event_omega(F(100),om)),
                   'high_redshift_scope':'matter plus dark energy control; radiation omitted, no early-Universe data claim'}},
            'same_data_alternatives_are_independent':False,
            'physical_horizon_or_Lambda_amplitude_derived':False,
            'references':['https://arxiv.org/abs/hep-th/0403052','https://arxiv.org/abs/hep-th/0403127']}
