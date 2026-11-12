"""Independent matrix derivatives, source reversals and dipole diagnostics."""
import numpy as np
from scipy import linalg as la
from scipy import sparse as sp
from model.openwave_field_observables.lattice import gradient,curl
from .solver import dirac,source_current,geometry,solve,SPIN,BETA,GAMMA,moment


def absolute(h):
    values,vectors=la.eigh(h)
    return (vectors*abs(values))@vectors.conj().T


def derivative_controls():
    n=3;N=9;a=.8;field=np.zeros((N,)*3+(3,))
    h,phases=dirac(n,a,N,field)
    d=h.toarray();energy,v=la.eigh(d)
    assert min(abs(energy))>=8
    i,j,axis,grid,phase=phases[7];delta=2e-5
    plus=field.copy();minus=field.copy()
    plus[grid+(axis,)]+=delta;minus[grid+(axis,)]-=delta
    hp=dirac(n,a,N,plus)[0].toarray();hm=dirac(n,a,N,minus)[0].toarray()
    numerical_d=(hp-hm)/(2*delta)
    nn=n**3
    derivative=np.zeros_like(d)
    for row in range(4):
        for col in range(4):
            value=-GAMMA[axis][row,col]/2-.5j*.5*BETA[row,col]
            derivative[row*nn+i,col*nn+j]=value
            derivative[col*nn+j,row*nn+i]=value.conjugate()
    assert la.norm(derivative-numerical_d)<1e-9
    transformed=v.conj().T@derivative@v
    ratio=(energy[:,None]+energy[None,:])/(abs(energy[:,None])+abs(energy[None,:]))
    frechet=v@(ratio*transformed)@v.conj().T
    numerical_k=(absolute(hp)-absolute(hm))/(2*delta)
    error=la.norm(frechet-numerical_k)
    assert error<2e-7
    same_sign_errors=[]
    for sign in (-1,1):
        basis=v[:,energy*sign>0]
        same_sign_errors.append(la.norm(basis.conj().T@(frechet-sign*derivative)@basis))
    assert max(same_sign_errors)<1e-12
    # Arbitrary basis mixing inside a degenerate positive eigenspace obeys
    # the same operator identity, rather than only diagonal expectations.
    pos=np.flatnonzero(energy>0);pair=pos[:2]
    assert abs(energy[pair[1]]-energy[pair[0]])<1e-10
    rotation=np.array([[3,4j],[4j,3]],complex)/5
    basis=v[:,pair]@rotation
    assert la.norm(basis.conj().T@(frechet-derivative)@basis)<1e-12
    wrong=derivative.copy()
    for row in range(4):
        value=-GAMMA[axis][row,row]/2
        wrong[row*nn+i,row*nn+j]=value
        wrong[row*nn+j,row*nn+i]=value.conjugate()
    assert la.norm(wrong-numerical_d)>.1
    return {'full_Dirac_rows':len(d),'actual_link_D_derivative_residual':float(la.norm(derivative-numerical_d)),
        'absolute_kernel_Frechet_finite_difference_residual':float(error),
        'both_signed_spectral_subspace_current_identities':True,
        'degenerate_subspace_operator_identity':True,'omitted_Wilson_current_rejected':True}


def gauge_control(result):
    n=result['n'];N=result['field_n'];a=result['spacing']
    rng=np.random.default_rng(724)
    chi=.03*rng.normal(size=(N,)*3)
    changed=result['field']+gradient(chi,a)
    h,phases=dirac(n,a,N,changed,result['mass'],result['trap'],result['wilson'],result['external_B'])
    coords,_,_=geometry(n,a,N);off=(N-n)//2
    phase=np.exp(-1j*np.array([chi[tuple(x+off)] for x in coords]))
    psi=(result['state'].reshape(4,-1)*phase).ravel()
    residual=np.linalg.norm(h@psi-result['eigenvalue']*psi)
    current=source_current(psi,phases,N,result['wilson'])
    error=np.linalg.norm(current-result['current'])
    wrong=source_current(result['state'],phases,N,result['wilson'])
    assert residual<1e-10 and error<1e-12 and np.linalg.norm(wrong-current)>1e-3
    assert np.linalg.norm(curl(changed,a)-curl(result['field'],a))<1e-12
    return {'same_state_and_field_gauge_covariance_residual':float(residual),
        'bond_current_gauge_covariance_residual':float(error),
        'unchanged_state_in_changed_link_gauge_rejected':True}


def dipole_profile(result):
    field=result['field'];N=result['field_n'];a=result['spacing'];center=N//2
    b=curl(field,a);mu=result['moment'][2];kappa=result['coupling']
    rows=[]
    for distance in (4,5,6):
        offset=int(round(distance/a))
        if center+offset>=N:continue
        radius=offset*a
        # xy faces carry B_z at their plaquette centers. Average four adjacent
        # faces to the z-axis vertex before comparing with the dipole profile.
        value=np.mean([b[center+dx,center+dy,center+offset,0] for dx in (0,-1) for dy in (0,-1)])
        expected=kappa*mu/(2*np.pi*radius**3)
        rows.append({'radius':radius,'axis_Bz':float(value),'dipole_from_independent_current_moment':float(expected),
                     'ratio_to_infinite_space_dipole':float(value/expected)})
    return {'moment_from_actual_bonds':result['moment'].tolist(),'periodic_boundary_profile':rows,
        'dipole_coefficient_fit_to_field':False,'infinite_space_or_continuum_profile_asserted':False}


def external_probe(result,step=2e-4):
    common={k:result[k] for k in ('n','spacing','field_n','coupling','mass','trap','wilson','spin_sign')}
    plus=solve(**common,external_B=step,initial_field=result['field'])
    minus=solve(**common,external_B=-step,initial_field=result['field'])
    derivative=(plus['vacuum_relative_rest_energy']-minus['vacuum_relative_rest_energy'])/(2*step)
    predicted=-result['moment'][2]
    assert abs(derivative-predicted)<2e-6
    return {'uniform_external_probe_step':step,'self_consistent_total_energy_derivative':float(derivative),
        'minus_stationary_current_moment':float(predicted),'difference':float(derivative-predicted),
        'g_from_actual_bound_energy_and_probe':float(2*result['vacuum_relative_rest_energy']*derivative/result['spin']),
        'field_and_state_reoptimized_for_each_probe':True,
        'coupled_time_dependent_Larmor_dynamics_computed':False}
