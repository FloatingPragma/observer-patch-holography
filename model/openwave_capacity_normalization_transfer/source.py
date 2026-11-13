"""One rigorous quadrature replaces seven ODEs in the regular zero-k limit."""
from fractions import Fraction as F
from functools import lru_cache
import sympy as sy
from model.openwave_radiation_capacity.source import Interval,exp_point,sigmoid,taylor_step,root_interval,log_interval
from model.openwave_capacity_transfer.source import background_inputs,integrate
from model.openwave_capacity_stress.source import enclosure

STEP=F(1,128)
TAIL_T=F(12)
REDSHIFTS=(0,1,10,1000)
ONE=Interval.point(1)


def exact_mode():
    m,r,s,b,h,A,M,T,R,hN=sy.symbols('m r s b h A M T R hN')
    Phi=R-T;Tn=Phi+h*T;P=-Tn
    y=s*T;z=b*T+s*Tn;dm=-3*T;dr=-4*T
    bN=-hN*s-(3+h)*b-M*s+2*h*A
    zN=bN*T+2*b*Tn+s*(hN*T+(h-1)*Tn)
    residuals=[P-(-Phi+(3*m*T+4*r*T+s*y)/2),
      -3*Tn-3*P,Tn-(h*T+Phi),-4*Tn-4*P,
      Tn-((1+h)*T+dr/4+Phi),b*T+s*Tn-z,
      zN-(-(3+h)*z-M*y+4*s*P-2*A*Phi),
      3*(P+Phi)+(3*m*dm+3*r*dr+s*z-s*s*Phi+A*y)/2]
    out=[sy.expand(v.subs(b,-(3+h)*s-A).subs(h,-(3*m+4*r+s*s)/2))for v in residuals]
    assert out==[0]*8
    x,u=sy.symbols('x u',nonnegative=True)
    hh=-(3+u)/2+(1+u)*x*x/2+x**3
    xp=x*(1-x*x)*(1+u+2*x)/2;up=-u*(1-u)
    hn=sy.expand(sy.diff(hh,x)*xp+sy.diff(hh,u)*up)
    assert sy.expand(hn-((1+u)*x+3*x*x)*xp-(1-x*x)*u*(1-u)/2)==0
    assert F(4,27)<F(4,25)
    upper=5*F(4,5)+F(1,2)*F(1,4)
    assert upper==F(33,8) and max(F(9),upper-1)==9
    return {'regular_adiabatic_mode':'T=R(H/a) integral_0^a da/H; Phi=R-T, U_m=U_r=T, delta_m=-3T, delta_r=-4T, y=sT, z=s_N T+s(Phi+hT).',
        'normalization':'R=3/2 gives unit potential in the radiation asymptotic regime.',
        'all_seven_evolution_and_Hamiltonian_constraint_residuals':list(map(str,out)),
        'momentum_constraint_also_checked':True,'comoving_dust_density_exact':'0 in the regular k=0 limit',
        'comoving_radiation_density_exact':'0 in the regular k=0 limit','comoving_scalar_displacement_exact':'0',
        'U_definition':'Regular limit of aH theta/k^2; no division by zero is performed.',
        'h_formula':'h=-(3+r_mix)/2+(1+r_mix)x^2/2+x^3, x=sqrt(Omega_dark).',
        'h_N_polynomial_exact':str(hn),'h_interval_exact':['-2','0'],
        'h_N_interval_exact':['0','33/8'],
        'lambda_control_h_N':'One half the variance of dilution exponents 3,4,0 under weights Omega_m,Omega_r,Omega_Lambda; hence 0<=h_N<=2.',
        'F_second_derivative':'For t=-log(a), F=e^(-t)/(H/H0), F_t=(h-1)F and F_tt=[(h-1)^2-h_N]F; hence abs(F_tt)<=9F.'}


def normalization():
    H0,l,Mp,Om,lam,k,a,En,R=sy.symbols('H0 l Mp Om lam k a En R',positive=True)
    rho=3*Mp**2*H0**2*Om;radius=l/H0;Q=k*k/(a*a*H0*H0*En*En)
    assert sy.simplify(rho.subs(H0,lam*H0)-lam**2*rho)==0
    assert sy.simplify(radius.subs(H0,lam*H0)-radius/lam)==0
    assert sy.simplify(Q.subs({H0:lam*H0,k:lam*k},simultaneous=True)-Q)==0
    data=background_inputs();h=F(data['H0_km_s_Mpc'])
    return {'exact_scale_symmetry':'H0->lambda H0, all present densities->lambda^2 rho, horizon length->L/lambda, k->lambda k; dimensionless fractions, scalar field/M_Pl, k/H0 and normalized transfers are unchanged.',
        'scale_residuals_exact':['0','0','0'],
        'amplitude_symmetry':'The complete homogeneous linear perturbation system admits v->A v, including A=0. The regular mode has Phi+T=R; R is a supplied normalization.',
        'examples':[{'H0_km_s_Mpc_exact':str(h*q),'density_multiplier_exact':str(q*q),'horizon_length_multiplier_exact':str(1/q)}for q in(F(1,2),F(1),F(2))],
        'absolute_present_dark_density':'3 M_Pl^2 H0^2 Omega_dark0; H0 and present fractions are inherited observational inputs.',
        'architecture_uniqueness_needed':False,'native_horizon_capacity_or_physical_ruler_relation_supplied':False,
        'absolute_dark_energy_or_primordial_amplitude_fixed_by_dimensionless_flow':False}


def integrand(t,w,kind,om,orr,ox):
    a=exp_point(-t);a2=a*a;a3=a2*a;a4=a2*a2
    if kind=='capacity':
        omega=Interval(sigmoid(w.lo).lo,sigmoid(w.hi).hi)
        return a3*root_interval((ONE-omega)/(om*a+orr))
    if kind=='lambda':return a3/root_interval(om*a+orr+ox*a4)
    raise ValueError('unknown background')


@lru_cache(None)
def validated(step=STEP,tail_t=TAIL_T,redshifts=REDSHIFTS):
    step,tail_t=F(step),F(tail_t)
    if not 0<step<=F(1,2) or tail_t<=0:raise ValueError('positive supported quadrature parameters required')
    if tuple(sorted(set(redshifts)))!=redshifts or redshifts[0]!=0:raise ValueError('increasing redshifts starting at zero required')
    data=background_inputs();om=F(data['Omega_m0']);orr=F(data['Omega_r0']);ox=1-om-orr
    targets={F(log_interval(1+z)[0]):z for z in redshifts}
    if max(targets)>=tail_t:raise ValueError('tail endpoint must precede every requested scale factor')
    count=(tail_t/step).__floor__()
    points=sorted(set([n*step for n in range(count+1)]+[tail_t]+list(targets)))
    w=Interval(*log_interval(ox/(1-ox)));states=[];previous=F(0)
    for j,t in enumerate(points):
        if j:w=taylor_step(previous,w,t-previous,om,orr)
        states.append({'t':t,'w':w,**{kind:integrand(t,w,kind,om,orr,ox)for kind in('capacity','lambda')}})
        previous=t
    amin=exp_point(-tail_t)
    tail=amin.hi**3/(3*root_interval(Interval.point(orr)).lo)
    integrals={kind:Interval(0,tail)for kind in('capacity','lambda')};errors={kind:F(0)for kind in integrals};rows=[]
    for j in range(len(states)-1,-1,-1):
        row=states[j];t=row['t']
        if j<len(states)-1:
            nxt=states[j+1];dt=nxt['t']-t
            for kind in integrals:
                error=Interval.point(F(3,4)*dt**3*row[kind].hi).hi
                integrals[kind]=integrals[kind]+dt*(row[kind]+nxt[kind])/2+Interval(-error,error)
                errors[kind]+=error
        if t not in targets:continue
        z=targets[t];a_inv=F(1+z);uncertainty=log_interval(1+z)[1]-t
        endpoint=Interval(row['w'].lo-4*uncertainty,row['w'].hi)
        omega=Interval(sigmoid(endpoint.lo).lo,sigmoid(endpoint.hi).hi)
        Hs={'capacity':root_interval(Interval.point(om*a_inv**3+orr*a_inv**4)/(ONE-omega)),
            'lambda':root_interval(Interval.point(om*a_inv**3+orr*a_inv**4+ox))}
        result={}
        for kind,H in Hs.items():
            exact_integral=Interval(integrals[kind].lo-uncertainty*row[kind].hi,integrals[kind].hi)
            assert exact_integral.lo>0
            T=F(3,2)*a_inv*H*exact_integral;phi=Interval.point(F(3,2))-T
            assert phi.lo>0 and phi.hi<1
            result[kind]={'metric_transfer_interval':phi.encoded(),'metric_transfer_display':enclosure((phi.lo,phi.hi)),
                          'T_interval':T.encoded(),'integral_interval':exact_integral.encoded(),
                          'H_over_H0_interval':H.encoded(),'quadrature_remainder_upper_exact':str(errors[kind])}
        cp=Interval(*map(F,result['capacity']['metric_transfer_interval']));lp=Interval(*map(F,result['lambda']['metric_transfer_interval']));ratio=cp/lp
        rows.append({'redshift':z,**result,'capacity_to_lambda_metric_ratio_interval':ratio.encoded(),
                     'capacity_to_lambda_metric_ratio_display':enclosure((ratio.lo,ratio.hi))})
    return {'step_exact':str(step),'tail_log_redshift_exact':str(tail_t),'nodes':len(points),
            'initial_integral_tail_upper_exact':str(tail),'background_interval_method':'Original monotone Taylor2 logit enclosure with complete 8*h^3/6 remainder and outward rational arithmetic.',
            'quadrature_method':'Composite trapezoids plus local 3*h^3*F_left_upper/4 from abs(F_tt)<=9F; all node and endpoint intervals retained.',
            'records':sorted(rows,key=lambda row:row['redshift'])}


def diagnostics(cert):
    rows=[]
    for kind in('capacity','lambda'):
        interval=list(map(F,cert['records'][0][kind]['metric_transfer_interval']))
        for k in(1e-4,2e-4):
            result,_=integrate(k,kind=kind,rtol=2e-11)
            last=result['records'][-1];phi=last['metric_potential_transfer']
            assert float(interval[0])<phi<float(interval[1])
            rows.append({'background':kind,'k_over_H0':k,'present_metric_transfer':phi,
                'present_comoving_dust_transfer':last['comoving_matter_density_transfer'],
                'sampled_constraint_residual':result['maximum_sampled_normalized_constraint_residual'],
                'inside_validated_zero_k_metric_interval':True})
    return rows


def run():
    mode=exact_mode();norm=normalization();cert=validated();diag=diagnostics(cert)
    return {'schema':'oph.capacity_normalization_and_validated_metric_transfer.v1',
      'verdict':'PASS_EXACT_NORMALIZATION_FREEDOM_AND_VALIDATED_ADIABATIC_METRIC_TRANSFER',
      'inherited_inputs':background_inputs(),'exact_normalization_theorem':norm,'exact_regular_mode':mode,
      'validated_metric_quadrature':cert,'independent_full_seven_ODE_small_k_diagnostics':diag,
      'power_law_controls':[{'background':'pure radiation','h_exact':'-2','T_exact':'1/2','Phi_exact':'1'},
                            {'background':'pure dust','h_exact':'-3/2','T_exact':'3/5','Phi_exact':'9/10'}],
      'mode_and_observation_scope':{'R_primordial_normalization_exact':'3/2, supplied for the displayed unit-radiation-potential transfer.',
        'strict_zero_k_mode_is_local_density_growth':False,'nonzero_k_error_estimate_supplied':False,
        'full_seven_ODE_diagnostics_are_interval_certificates':False,
        'metric_ratio_is_parent_finite_k_matter_growth_ratio':False,
        'strict_zero_k_adiabatic_mode':'A residual adiabatic coordinate mode; the computed coefficient specifies its regular long-wavelength limit, not a measured spatial perturbation.',
        'native_capacity_to_scalar_stress_horizon_and_clock_map_proved':False,
        'primordial_statistics_recombination_neutrino_free_streaming_or_Boltzmann_transfer':False,
        'absolute_H0_or_dark_energy_or_primordial_amplitude_predicted':False,
        'new_empirical_fit':False,'new_Nature_match':False,'OpenWave_row28_never_fitted_amplitude_obligation_met':False},
      'finite_theorem_shortcut':'Exact adiabatic substitution and constraint propagation reduce seven coupled equations in the regular zero-k limit to one positive integral; the original validated background and an elementary derivative bound certify it.',
      'scope':'A rigorously enclosed normalized metric response of the supplied Einstein/dust/radiation/canonical-scalar comparison and an exact normalization underdetermination theorem. It supplies no physical capacity map or independent cosmological amplitude.'}
