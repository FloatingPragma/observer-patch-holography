"""Finite source spectral calculations with a uniform analytic mass limit."""
from fractions import Fraction as F
import numpy as np
from scipy.linalg import eigh
from model.openwave_angular_torque.source import oscillator,operators,BETA


def stable_kernel(data,mass,b,external):
    h,d,e,trap=operators(data,mass=mass,omega=1.)
    A=h-np.kron(BETA,mass*np.eye(data['spatial_dimension']))
    A+=sum((v*x for v,x in zip(b,d)),np.zeros_like(A))+sum((v*x for v,x in zip(external,e)),np.zeros_like(A))
    T=A@A;w,u=eigh(T);w=np.maximum(w,0);roots=np.sqrt(mass*mass+w)
    shift=w/(roots+mass);K=(u*shift)@u.conj().T+trap
    def derivative(E):
        Y=u.conj().T@(A@E+E@A)@u
        return u@(Y/(roots[:,None]+roots[None,:]))@u.conj().T
    return K,tuple(derivative(x)for x in d),tuple(derivative(x)for x in e),A


def axial_subspace(data,sign):
    n=data['spatial_dimension'];w,u=eigh(data['J'][2][:2*n,:2*n]);small=u[:,np.abs(w-sign/2)<1e-10]
    return np.vstack((small,np.zeros_like(small)))


def eigenline(data,mass,b,external,sign):
    K,D,E,A=stable_kernel(data,mass,b,external);U=axial_subspace(data,sign)
    values,vectors=eigh(U.conj().T@K@U);psi=U@vectors[:,0]
    residual=np.linalg.norm(K@psi-values[0]*psi)
    assert residual<1e-10
    return float(values[0]),psi,D,E,residual


def bounds(mass,field=F(1,10)):
    m=F(mass);B=F(field);a0=F(8,5);r=F(3,5);c1=F(5);aB=a0+F(4,5)*abs(B);kappa=1/m**6
    # Full finite operator certificates justify these rational upper bounds.
    eps0=F(5,4)/m;eta=2*eps0/(1-eps0)
    static=12*eta+F(15,2)/m**2
    eps=aB*aB/(2*m);remainder=aB**4/(8*m**3);schur=eps*eps/(1-eps)
    bb=kappa*r*aB/(m*c1-kappa*r*r)
    delta=((aB+r*bb)**2-aB*aB)/(2*m)
    hessian=r*r/m+r*r*(aB+r*bb)**2/m**3
    contraction=kappa*hessian/c1
    rest_shift=delta+7*bb*bb/(2*kappa)
    gap=1-eps0-2*delta
    self_state_distance=2*delta/gap
    probe_norm=F(4,5)*(a0+r*bb)/m
    cross_hessian=F(4,5)*r/m+F(4,5)*r*(a0+r*bb)**2/m**3
    moment_change=2*probe_norm*self_state_distance+cross_hessian*bb
    gstatic=static+4*m*moment_change+(2+static+4*m*moment_change)*(eps0+rest_shift)/m
    static_contraction=kappa*(hessian+2*(r*(a0+r*bb)/m)**2/gap)/c1
    etaB=2*eps/(1-eps);radius=F(1,2)-5*etaB
    duration=64*m;state_error=duration*delta
    assert radius>0 and 5*(2+state_error)*state_error<radius/2
    phase_error=10*(2+state_error)*state_error/radius
    gpre=2*(eps0+rest_shift)/m+4*(m+eps0+rest_shift)*(remainder+schur)/abs(B)
    total_pre=gpre+2*(m+eps0+rest_shift)*phase_error/(abs(B)*duration)
    return {'mass':m,'field':B,'coupling':kappa,'ground_energy_excess_upper':eps0,'static_g_error_from_two':gstatic,
        'frozen_precession_g_error':gpre,'quasistatic_precession_g_error_from_two':total_pre,
        'self_field_norm_upper':bb,'self_kernel_perturbation':delta,'rest_energy_self_shift_bound':rest_shift,
        'self_field_contraction':contraction,'axial_stationary_branch_contraction':static_contraction,
        'vacuum_state_distance_bound':eta,'precession_duration':duration,'uniform_precession_radius_lower':radius,
        'quasistatic_state_error':state_error,'uniform_phase_error':phase_error}


def row(mass):
    data=oscillator();B=.1;kappa=mass**-6;b=np.zeros(3)
    for _ in range(8):
        shift,psi,D,E,res=eigenline(data,mass,b,(0,0,0),1)
        current=-np.array([np.vdot(psi,x@psi).real for x in D]);target=kappa*current/data['magnetic_metric']
        if np.linalg.norm(target-b)<kappa*1e-12:break
        b=target
    field_res=np.linalg.norm(target-b)/kappa
    rest=mass+shift+data['magnetic_metric']*np.dot(b,b)/(2*kappa)
    moment=-np.vdot(psi,E[2]@psi).real;j=np.vdot(psi,data['J'][2]@psi).real
    J2=sum((x@x for x in data['J']),np.zeros_like(data['J'][0]));j2_values,j2_vectors=eigh(J2)
    j2_probabilities=abs(j2_vectors.conj().T@psi)**2
    zero=eigenline(data,mass,np.zeros(3),(0,0,0),1)
    zero_casimir_residual=float(np.linalg.norm(J2@zero[1]-.75*zero[1]));assert zero_casimir_residual<1e-10
    static=2*rest*abs(moment)/abs(j)
    up=eigenline(data,mass,np.zeros(3),(0,0,B),1)
    down=eigenline(data,mass,np.zeros(3),(0,0,B),-1)
    gap=up[0]-down[0];assert gap>0
    overlap=np.vdot(up[1],(data['J'][0]+1j*data['J'][1])@down[1])
    downpsi=down[1]*np.exp(-1j*np.angle(overlap));amplitude=np.vdot(up[1],(data['J'][0]+1j*data['J'][1])@downpsi)/2
    assert abs(amplitude.imag)<1e-12 and amplitude.real>.49
    analytic=bounds(mass);time=float(analytic['precession_duration']);delta=float(analytic['self_kernel_perturbation'])
    state_error=time*delta;radius=float(amplitude.real)
    # J+ norm <= 2(N+1/2); expectation and phase transfer retain every mode.
    expectation_error=5*(2+state_error)*state_error
    assert expectation_error<radius/2
    phase_error=2*expectation_error/radius
    dynamical_g_error=2*(mass+float(analytic['ground_energy_excess_upper']))*phase_error/(B*time)
    central=2*rest*gap/B
    samples=[]
    w,u=eigh(data['J'][0])
    for t in np.linspace(0,time,9):
        v=(np.exp(-1j*up[0]*t)*up[1]+np.exp(-1j*down[0]*t)*downpsi)/np.sqrt(2)
        prob=abs(u.conj().T@v)**2;assert abs(prob.sum()-1)<1e-12
        samples.append({'time':float(t),'all_Jx_values':w.tolist(),'all_Jx_probabilities':prob.tolist(),
            'angular_xy':[float(np.vdot(v,data['J'][a]@v).real)for a in (0,1)]})
    assert abs(central-2)<float(analytic['frozen_precession_g_error'])+dynamical_g_error
    certified_error=float(analytic['quasistatic_precession_g_error_from_two'])
    return {'mass':mass,'coupling':kappa,'external_field':B,'static_stationary_total_energy':rest,
        'static_total_angular_momentum':float(j),'static_magnetic_moment':float(moment),
        'zero_self_field_J_squared_eigenvalue':'3/4','zero_self_field_Casimir_residual':zero_casimir_residual,
        'static_axial_J_squared_complete_values':j2_values.tolist(),'static_axial_J_squared_complete_probabilities':j2_probabilities.tolist(),
        'finite_self_field_exact_Casimir_claimed':False,
        'static_g_supplied_mass':4*mass*abs(moment),'static_g_bound_energy':float(static),
        'stationary_field_amplitude':b.tolist(),'scaled_field_equation_residual':float(field_res),'stationary_eigen_residual':float(res),
        'actual_full_kernel_energy_difference':float(gap),'finite_precession_duration':float(time),'precession_radius':radius,
        'precession_g_bound_energy':float(central),'quasistatic_state_error_bound':state_error,
        'analytic_g_distance_from_two_upper':certified_error,
        'analytic_rational_bounds':{k:str(v)for k,v in analytic.items()},
        'static_self_field_additional_error_scope':'The axial stationary-branch contraction and spectral gap give a self-consistent state and an analytic static g bound. The displayed stationary vectors and moments are numerical evaluations of that branch.',
        'complete_separate_copy_records':samples}


def run():
    from .exact import certificate as exact_certificate
    data=oscillator();K,D,E,A=stable_kernel(data,1.,np.zeros(3),np.zeros(3));T=A@A
    polynomial=np.linalg.norm(T@T-2.5*T);assert polynomial<1e-12
    eig=np.linalg.eigvalsh(T);assert eig[-1]<2.5+1e-12
    probe=operators(data,1.,1.)[2][2];W=A@probe+probe@A
    assert np.linalg.norm(A,2)<1.6 and np.linalg.norm(probe,2)<.8 and np.linalg.norm(W,2)<3
    vacuum=np.zeros(len(A));vacuum[0]=1
    assert abs(np.vdot(vacuum,W@vacuum)-1)<1e-12
    zero_spin_mutation=np.vdot(vacuum,(np.kron(np.eye(4),data['L'][2]))@vacuum).real
    assert zero_spin_mutation==0
    rows=[row(m)for m in (64,256,1024,4096,16384)]
    return {'schema':'oph.rotational_native_g_limit.v1','verdict':'PASS','source':'same complete-shell rotational positive kernel and Maxwell modes as angular_torque; Omega=1, radius=6/5, native q=-1, coupling=m^-6',
        'exact_full_matrix_certificate':exact_certificate(),
        'full_matrix_shortcuts':{'A_squared_minimal_polynomial_residual':float(polynomial),'A_squared_eigenvalues':[0,'5/2'],
            'A_norm_upper':'8/5','uniform_probe_norm_upper':'4/5','squared_kernel_probe_norm_upper':3,
            'vacuum_squared_kernel_probe_expectation':1,'spin_deleted_mutation_expectation':float(zero_spin_mutation)},
        'rows':rows,'theorem':'At fixed complete shell cutoff the self-consistent axial static g and full-kernel quasistatic precession g converge to2 as m grows, with kappa=m^-6. Both use actual J, the same positive kernel, native charge and bound-energy reference.',
        'scope':'The low-momentum parameter route is a declared regulator construction, not a physical mass or alpha fit. The finite-cutoff theorem does not establish the cubic-source continuum limit. Numerical roots and stationary residuals are separate from analytic uniform bounds; no narrow floating-point interval is called a spectral enclosure.'}
