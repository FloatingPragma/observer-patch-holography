"""Einstein/fluid/Klein--Gordon transfer with exact constraint propagation.

The equations are exact for the declared classical ideal-fluid/scalar model.
ODE outputs are floating numerical diagnostics, not rigorous interval solutions
or a fit to CMB data. No native capacity-to-scalar transfer is assumed proved.
"""
from functools import lru_cache
from fractions import Fraction
import json
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, quad

ROOT = Path(__file__).resolve().parents[2]


@lru_cache(None)
def constraint_identity():
    m,r,s,b,h,A,M,Q,Phi,dm,um,dr,ur,y,z = sp.symbols('m r s b h A M Q Phi dm um dr ur y z')
    P = -Phi+(3*m*um+4*r*ur+s*y)/2
    D = 3*m*dm+3*r*dr+s*z-s*s*Phi+A*y
    pressure = r*dr+s*z-s*s*Phi-A*y
    C = Q*Phi+3*(P+Phi)+D/2
    derivatives = {m:(-3-2*h)*m, r:(-4-2*h)*r, s:b, A:M*s-2*h*A,
        Q:-2*(1+h)*Q, Phi:P, dm:-Q*um+3*P, dr:-sp.Rational(4,3)*Q*ur+4*P,
        um:h*um+Phi, ur:(1+h)*ur+dr/4+Phi, y:z,
        z:-(3+h)*z-(Q+M)*y+4*s*P-2*A*Phi}
    def prime(value):
        return sum(sp.diff(value,v)*dv for v,dv in derivatives.items())
    replacements = {b:-(3+h)*s-A}
    gravity = {h:-(3*m+4*r+s*s)/2}
    conservation = sp.expand((prime(C)+(3+2*h)*C).subs(replacements).subs(gravity))
    trace = sp.expand((prime(P)+(4+h)*P+(3+2*h)*Phi-pressure/2).subs(replacements).subs(gravity))
    assert conservation == trace == 0
    return {'Hamiltonian_constraint_propagation_exact': 'C_N=-(3+2h)*C',
            'spatial_trace_equation_exact': 'Phi_NN+(4+h)*Phi_N+(3+2h)*Phi=delta_p/(2*H^2*M_Pl^2)',
            'momentum_constraint': 'Phi_N+Phi=(3*Omega_m*U_m+4*Omega_r*U_r+s*y)/2',
            'symbolic_residuals_exact': ['0','0']}


@lru_cache(None)
def background_inputs():
    return json.loads((ROOT/'model/openwave_radiation_capacity/results.json').read_text())['inputs']


@lru_cache(None)
def mass_function():
    parent=json.loads((ROOT/'model/openwave_capacity_perturbations/results.json').read_text())
    x,r=sp.symbols('x r')
    expression=sp.sympify(parent['canonical_scalar_completion']['V_second_derivative_over_H_squared_polynomial'],locals={'x':x,'r':r})
    return sp.lambdify((x,r),expression,'math')


@lru_cache(None)
def background(kind='capacity'):
    if kind not in ('capacity','lambda'):
        raise ValueError('unknown declared background')
    data=background_inputs()
    om,orad,ox=(float(Fraction(data[name])) for name in ('Omega_m0','Omega_r0','Omega_dark0'))
    # The capacity background is integrated independently from its supplied
    # present fraction; the radiation ratio is an exact analytic function.
    def x_rhs(N,value):
        x=value[0]; r=orad/(om*np.exp(N)+orad)
        return [x*(1-x*x)*(1+r+2*x)/2]
    sol=solve_ivp(x_rhs,(0.,np.log(1e-14)),[np.sqrt(ox)],method='DOP853',
                  rtol=2e-13,atol=1e-15,max_step=.08,dense_output=True)
    assert sol.success
    def evaluate(N):
        a=np.exp(N)
        if kind=='lambda':
            H2=om/a**3+orad/a**4+ox; m=om/a**3/H2; r=orad/a**4/H2
            return m,r,0.,0.,0.,0.,-1.5*m-2*r,H2
        x=float(sol.sol(N)[0]); ratio=orad/(om*a+orad)
        xp=x*(1-x*x)*(1+ratio+2*x)/2
        H2=(om/a**3+orad/a**4)/(1-x*x)
        m=(1-ratio)*(1-x*x); r=ratio*(1-x*x)
        s=np.sqrt(2*x*x*(1-x)); b=s*xp*(1/x-1/(2*(1-x)))
        h=-1.5*m-2*r-s*s/2; A=-b-(3+h)*s; mass=float(mass_function()(x,ratio))
        assert 0<x<1 and 2-1e-12<=mass<=6+1e-12
        return m,r,s,b,A,mass,h,H2
    return evaluate


def components(N,state,k,kind):
    Dm,um,Dr,relative_velocity,y,z,Phi=state
    ur=um+relative_velocity
    m,r,s,b,A,mass,h,H2=background(kind)(N)
    Q=k*k/(np.exp(2*N)*H2)
    P=-Phi+(3*m*um+4*r*ur+s*y)/2
    # Algebraically cancel the order-one velocity terms before evaluation.
    # This is the same 00 constraint, in comoving-density coordinates.
    terms=[Q*Phi,1.5*m*Dm,1.5*r*Dr,s*z/2,-s*s*Phi/2,(A+3*s)*y/2]
    C=sum(terms)
    scale=1+sum(abs(t)for t in terms)
    return m,r,s,b,A,mass,h,H2,Q,P,C,scale


def rhs(N,state,k,kind):
    Dm,um,Dr,relative_velocity,y,z,Phi=state
    ur=um+relative_velocity
    m,r,s,b,A,mass,h,H2,Q,P,C,scale=components(N,state,k,kind)
    # Dm=delta_m+3Um, Dr=delta_r+4Ur, and relative_velocity=Ur-Um.
    # Evolving these small differences directly avoids superhorizon loss of
    # precision; the exact symbolic theorem concerns the equivalent equations.
    return np.array([-Q*um+6*r*relative_velocity+1.5*s*(y-s*um),h*um+Phi,
                     Dr-4*Q*ur/3-6*m*relative_velocity+2*s*(y-s*ur),
                     h*relative_velocity+Dr/4,z,
                     -(3+h)*z-(Q+mass)*y+4*s*P-2*A*Phi,P])


def initial(k,kind,a0):
    N=np.log(a0);m,r,s,b,A,mass,h,H2=background(kind)(N)
    Q=k*k/(a0*a0*H2)
    # The regular k=0 growing mode obeys Phi+T=3/2 and
    # T=(3/2)*(H/a)*integral_0^a da'/H(a'). This fixes a common unit
    # primordial potential, rather than independently normalizing each start.
    lower=1e-14/a0
    integral=quad(lambda t:np.sqrt(H2/background(kind)(N+np.log(t))[7]),
                  lower,1.,epsabs=2e-14,epsrel=2e-13)[0]
    orad=float(Fraction(background_inputs()['Omega_r0']))
    omitted_tail=np.sqrt(H2)*a0*a0*lower**3/(3*np.sqrt(orad))
    assert omitted_tail<1e-14
    T=1.5*integral; Phi=1.5-T
    # The explicitly supplied finite-k correction satisfies 00 exactly.
    state=np.array([0.,T,-2*Q*Phi/(3*r),0.,s*T,b*T+s*(Phi+h*T),Phi])
    assert abs(components(N,state,k,kind)[10])<2e-14
    return state


def integrate(k,kind='capacity',a0=1e-7,rtol=2e-11,defect=0.):
    if k<=0 or not 1e-9<=a0<=1e-4:
        raise ValueError('positive k/H0 and declared early start required')
    N0=np.log(a0);state=initial(k,kind,a0);state[2]+=defect
    result=solve_ivp(lambda N,v:rhs(N,v,k,kind),(N0,0.),state,method='DOP853',
                     rtol=rtol,atol=rtol*1e-3,max_step=.06,dense_output=True)
    assert result.success
    rows=[]
    for redshift in (1000,10,1,0):
        N=-np.log1p(redshift);v=result.sol(N);derivative=rhs(N,v,k,kind)
        c=components(N,v,k,kind);Delta=v[0]
        rows.append({'redshift':redshift,'metric_potential_transfer':float(v[6]),
            'metric_potential_N_transfer':float(c[9]),'comoving_matter_density_transfer':float(Delta),
            'matter_logarithmic_growth':float(derivative[0]/Delta),
            'scalar_field_transfer_over_M_Pl':float(v[4]),
            'Hamiltonian_constraint_normalized_residual':float(c[10]/c[11])})
    samples=np.linspace(N0,0.,601)
    residuals=[abs(components(N,result.sol(N),k,kind)[10]/components(N,result.sol(N),k,kind)[11])for N in samples]
    return {'background':kind,'k_over_H0':k,'initial_scale_factor':a0,
            'rtol':rtol,'atol':rtol*1e-3,'rhs_evaluations':result.nfev,
            'maximum_sampled_normalized_constraint_residual':float(max(residuals)),
            'records':rows}, result


def run():
    identities=constraint_identity();pairs=[]
    for k in (.1,1.,10.,100.):
        coarse,_=integrate(k,rtol=2e-10)
        fine,_=integrate(k,rtol=2e-12)
        earlier,_=integrate(k,a0=1e-8,rtol=2e-12)
        control,_=integrate(k,kind='lambda',rtol=2e-12)
        fields=('metric_potential_transfer','metric_potential_N_transfer','comoving_matter_density_transfer',
                'matter_logarithmic_growth','scalar_field_transfer_over_M_Pl')
        error=max(abs(a[key]-b[key])/max(1.,abs(b[key])) for a,b in zip(coarse['records'],fine['records'],strict=True) for key in fields)
        start_error=max(abs(a[key]-b[key])/max(1.,abs(b[key])) for a,b in zip(earlier['records'],fine['records'],strict=True) for key in fields)
        assert error<2e-6 and start_error<2e-5
        assert fine['maximum_sampled_normalized_constraint_residual']<2e-7
        pairs.append({'capacity':fine,'lambda_control':control,
            'solver_refinement_max_scaled_difference':float(error),
            'earlier_start_max_scaled_difference':float(start_error),
            'present_comoving_matter_transfer_ratio_to_lambda':fine['records'][-1]['comoving_matter_density_transfer']/control['records'][-1]['comoving_matter_density_transfer']})
    # Deliberately violate the 00 constraint. Its nonzero residual must remain
    # visible and obey the exact homogeneous propagation law; it is not reset.
    faulty,solution=integrate(1.,defect=1e-5,rtol=2e-12)
    N0=np.log(1e-7);H20=background()(N0)[7]
    C0=components(N0,solution.sol(N0),1.,'capacity')[10]
    checks=[]
    for N in (N0/2,0.):
        actual=components(N,solution.sol(N),1.,'capacity')[10]
        expected=C0*np.exp(-3*(N-N0))*H20/background()(N)[7]
        checks.append({'log_scale_factor':float(N),'actual_constraint_defect':float(actual),
                       'propagated_initial_defect':float(expected)})
        assert abs(actual-expected)<2e-7*max(1.,abs(expected))
    return {'schema':'oph.capacity_constrained_linear_transfer.v1',
        'verdict':'PASS_EXACT_CONSTRAINT_IDENTITY_AND_NUMERICAL_TRANSFER_CONTROLS',
        'inherited_background_inputs':background_inputs(),
        'declared_effective_theory':'Einstein gravity with conserved ideal dust, ideal radiation without anisotropic stress, and the previously reconstructed canonical scalar',
        'variables':'N=log(a), h=H_N/H, s=phi_N/M_Pl, A=V_phi/(H^2 M_Pl), M=V_phiphi/H^2, Q=k^2/(aH)^2, U_i=aH theta_i/k^2, y=delta_phi/M_Pl, z=y_N',
        'exact_constraint_theorems':identities,
        'initial_data':'regular k=0 adiabatic mode: T=(3/2)*(H/a)*integral_0^a da/H, Phi=3/2-T, U_m=U_r=T, y=sT, z=s_N*T+s*(Phi+h*T); delta_m=-3T, delta_r=-4T-2Q*Phi/(3Omega_r). The finite-k density correction is supplied and tested by moving the start earlier.',
        'initial_integral_tail':'H(a)>=sqrt(Omega_r0)/a^2 bounds the omitted interval below a=1e-14; its dimensionless integral contribution is below1e-14',
        'normalization':'linear transfer per unit asymptotic radiation-era potential, not a physical perturbation of unit amplitude or a primordial-amplitude fit',
        'numerical_coordinates':'integrate Delta_m=delta_m+3U_m, Delta_r=delta_r+4U_r and U_r-U_m directly to avoid cancellation of superhorizon order-one variables; no constraint projection is applied',
        'numerical_transfer_comparisons':pairs,
        'nonzero_constraint_control':checks,
        'numerical_status':'floating ODE diagnostics with independent tolerance and initial-time controls; no rigorous interval integration error certificate',
        'references':['https://arxiv.org/abs/astro-ph/9506072','https://arxiv.org/abs/astro-ph/0304212'],
        'scope':{'full_linear_gravitational_constraints_in_declared_model':True,
            'metric_perturbations_frozen_or_omitted':False,'radiation_and_matter_evolved':True,
            'native_capacity_scalar_or_primordial_state_derived':False,
            'free_streaming_neutrinos_recombination_or_full_Boltzmann_hierarchy':False,
            'CMB_large_scale_structure_likelihood':False,'new_empirical_fit':False,
            'new_Nature_validation':False,'nonlinear_Einstein_or_quantum_continuum_completion':False}}
