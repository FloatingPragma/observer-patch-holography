"""Precession, complete records and independent finite-time controls."""
from fractions import Fraction as Q
from functools import lru_cache
import numpy as np
from scipy.linalg import eigh,expm
from model.openwave_intrinsic_magnetic.solver import solve,native_reduction,dirac,moment
from model.openwave_intrinsic_magnetic.checks import absolute
from .dynamics import (trajectory,stationary,diagonalize,evolve,current,kernel_apply,
    spin_operators,quasistatic,observe,field_minimizer,transverse)

PARAMS={'n':3,'a':1.,'field_n':9,'mass':8.,'trap':.01,'wilson':.5,'kappa':.1}


def contraction_certificate():
    mass=Q(8);r=Q(1,2);a=Q(1);M=mass+Q(3,100)+6*r/a+3/a
    first_squared=2*(1+r*r);second=a*(1+r)/2
    curvature=((1+(M/mass)**2)*first_squared+M*second)/mass
    sine_lower=Q(3,9)-Q(22,63)**3/6
    lap_lower=4*sine_lower**2/a**2
    ratio=Q(1,10)*curvature/(a**3*lap_lower)
    assert 0<ratio<Q(61,100)
    return {'D_norm_upper':str(M),'first_derivative_norm_squared_upper':str(first_squared),
        'second_derivative_norm_upper':str(second),'absolute_kernel_Hessian_upper':str(curvature),
        'transverse_Laplacian_gap_lower':str(lap_lower),'fixed_point_contraction_upper':str(ratio),
        'contraction_decimal':float(ratio),
        'exact_residual_to_field_error_factor':str(1/(1-ratio)),
        'unique_quasistatic_minimizer_for_every_normalized_state_and_uniform_probe':True,
        'proof_inputs':'3<pi<22/7; sin(x)>=x-x^3/6 for0<=x<=pi/2; all sixteen internal rows retain four identical blocks.'}


def current_controls():
    p=PARAMS;N=p['field_n'];field=np.zeros((N,)*3+(3,));w,v,ph=diagonalize(field,p,.037)
    rng=np.random.default_rng(120321);state=rng.normal(size=len(w))+1j*rng.normal(size=len(w));state/=np.linalg.norm(state)
    j=current(state,w,v,ph,p);left,right,axis,grid,phase=ph[19];step=2e-5
    fp=field.copy();fm=field.copy();fp[grid+(axis,)]+=step;fm[grid+(axis,)]-=step
    hp=dirac(p['n'],p['a'],N,fp,external_B=.037)[0].toarray()
    hm=dirac(p['n'],p['a'],N,fm,external_B=.037)[0].toarray()
    finite_difference=-np.vdot(state,(absolute(hp)-absolute(hm))@state).real/(2*step)
    error=abs(j[grid+(axis,)]-finite_difference)
    naive=-np.vdot(state,(hp-hm)@state).real/(2*step)
    assert error<2e-8 and abs(naive-finite_difference)>1e-4
    # Both signs occur in the prepared test state; deleting one is detected.
    signed=[float(sum(abs(v[:,w*sign>0].conj().T@state)**2)) for sign in (-1,1)]
    initial,A,_=stationary(p);obs,_=observe(initial,A,p,.05,True)
    moment_error=max(abs(row['expectation']-obs['moment'][row['axis']]) for row in obs['magnetic_moment_records'])
    assert moment_error<1e-12
    # Direct field inversion must project nonstationary current transversely.
    noisy=rng.normal(size=(N,)*3+(3,));projected=transverse(noisy,p['a'])
    from model.openwave_field_observables.lattice import divergence_adjoint
    assert np.linalg.norm(divergence_adjoint(projected,p['a']))<1e-12
    assert np.linalg.norm(divergence_adjoint(noisy,p['a']))>1
    return {'arbitrary_state_signed_branch_populations':signed,'current_Frechet_FD_error':float(error),
        'positive_branch_only_current_error':float(abs(naive-finite_difference)),
        'magnetic_projective_records_match_current_moment_error':float(moment_error),
        'longitudinal_current_without_projection_rejected':True}


def frozen_control(duration=2200.):
    p=PARAMS;state,field,_=stationary(p);w,v,ph=diagonalize(field,p,.05)
    h=dirac(p['n'],p['a'],p['field_n'],field,external_B=.05)[0].toarray()
    residual=np.linalg.norm(h@v-v*w,2);unitarity=np.linalg.norm(v.conj().T@v-np.eye(len(w)),2)
    exact=evolve(state,w,v,17.3,p['mass'])
    k=(v*abs(w))@v.conj().T
    direct=expm(-1j*17.3*(k-p['mass']*np.eye(len(w))))@state
    assert np.linalg.norm(exact-direct)<2e-12
    times=np.linspace(0,duration,45);spins=spin_operators(p['n']);rows=[]
    for t in times:
        psi=evolve(state,w,v,t,p['mass']);spin=[float(np.vdot(psi,s@psi).real) for s in spins]
        j=current(psi,w,v,ph,p)
        rows.append({'time':float(t),'spin':spin,'moment':moment(j,p['a']).tolist(),
            'self_field_defect':float(np.linalg.norm(field_minimizer(j,p['a'],p['kappa'])-field))})
    angle=np.unwrap([np.angle(x['spin'][0]+1j*x['spin'][1]) for x in rows])
    return {'scope':'exact spectral evolution of the frozen initial self-field; intentionally not the quasistatic coupled solution',
        'spin_frequency':float(np.polyfit(times,angle,1)[0]),'rows':rows,
        'full_signed_spectral_residual':float(residual),'eigenbasis_unitarity_residual':float(unitarity),
        'linear_unitary_error_bound_from_exact_residual_formula':float(duration*residual+2*unitarity),
        'direct_matrix_exponential_error':float(np.linalg.norm(exact-direct)),
        'floating_residuals_outward_rounded':False}


def public(row):return {k:v for k,v in row.items() if k not in ('final_state','final_field')}


@lru_cache(maxsize=1)
def run():
    algebra=native_reduction();proof=contraction_certificate();controls=current_controls()
    refinements=[trajectory(PARAMS,.05,h,2200.,int(100/h)) for h in (20.,10.,5.)]
    fine=refinements[-1]
    smaller=trajectory(PARAMS,.025,10.,4400.,10)
    reversed_B=trajectory(PARAMS,-.05,10.,2200.,10)
    zero=trajectory(PARAMS,0.,20.,200.,5)
    comparisons=[]
    for coarse in refinements[:-1]:
        fidelity=min(1.,abs(np.vdot(coarse['final_state'],fine['final_state']))**2)
        comparisons.append({'step':coarse['step'],'final_pure_state_trace_distance_diagnostic':float(2*np.sqrt(max(0,1-fidelity))),
            'spin_frequency_difference':float(coarse['spin_angular_frequency']-fine['spin_angular_frequency']),
            'final_field_difference':float(np.linalg.norm(coarse['final_field']-fine['final_field']))})
    assert abs(reversed_B['spin_angular_frequency']+fine['spin_angular_frequency'])<2e-8
    assert abs(zero['spin_angular_frequency'])<1e-10
    assert abs(fine['spin_angular_frequency']-fine['moment_angular_frequency'])<1e-7
    stationary_branches=[solve(n=3,spacing=1.,field_n=9,external_B=.05,spin_sign=sign) for sign in (1,-1)]
    splitting=stationary_branches[0]['vacuum_relative_rest_energy']-stationary_branches[1]['vacuum_relative_rest_energy']
    assert abs(splitting-fine['spin_angular_frequency'])<4e-8
    ground=refinements[0]['rows'][0]
    energy0=refinements[0]['initial_stationary']['energy']
    # Zero-probe field energy belongs in the physical reference convention.
    state,A,_=stationary(PARAMS);zero_observation,_=observe(state,A,PARAMS,0.)
    bound_energy=zero_observation['total_energy'];spin0=refinements[0]['initial_stationary']['spin'][0]
    mu0=refinements[0]['initial_stationary']['moment'][0]
    for row in refinements+[smaller,reversed_B,zero]:
        assert max(x['continuity_residual'] for x in row['rows'])<1e-11
        assert max(x['transverse_Maxwell_residual'] for x in row['rows'])<1e-10
        assert max(x['norm_error'] for x in row['rows'])<2e-10
        assert row['total_energy_drift']<1e-8
    return {'schema':'oph.quasistatic_larmor_dynamics.v1','native_algebra':algebra,
        'source':'same |D_W(A)| matter kernel and transverse field energy; a declared instantaneous quasistatic self-field dynamics, with fixed external uniform B',
        'complete_carrier':'all sixteen native internal rows and both D signs; the prepared particle occupies one complete Dirac block; three spectator blocks and the empty Fock vacuum have zero preparation probability',
        'theorem':proof,'independent_current_controls':controls,
        'time_step_refinement':[public(x) for x in refinements],'refinement_comparisons':comparisons,
        'half_external_field':public(smaller),'reversed_external_field':public(reversed_B),
        'zero_field_stationary_control':public(zero),'frozen_field_control':frozen_control(),
        'energy_split_comparison':{'stationary_total_energies':[float(x['vacuum_relative_rest_energy']) for x in stationary_branches],
            'stationary_total_energy_difference':float(splitting),'actual_spin_frequency':fine['spin_angular_frequency'],
            'actual_moment_frequency':fine['moment_angular_frequency'],'difference':float(fine['spin_angular_frequency']-splitting)},
        'normalizations':{'supplied_mass':8.,'zero_probe_actual_rest_energy':bound_energy,'measured_canonical_spin':spin0,
            'zero_probe_intrinsic_moment':mu0,'static_g_using_mass_and_measured_spin':-2*8*mu0/spin0,
            'static_g_using_bound_energy_and_measured_spin':-2*bound_energy*mu0/spin0,
            'precession_g_using_mass':2*8*fine['spin_angular_frequency']/.05,
            'precession_g_using_bound_energy':2*bound_energy*fine['spin_angular_frequency']/.05,
            'doublet_normalization':'A split two-level doublet precesses at DeltaE approximately2|mu|B. The measured projected canonical spin has length s<1/2; the static convention dividing by s differs from the doublet-precession convention by2s in the weak-field limit. No continuum spatial rotation generator is identified here.'},
        'limits':{'quasistatic_self_field_globally_unique_for_tested_parameters':True,
            'semiclassical_state_dependent_mean_field_law':True,
            'linear_CP_channel_on_arbitrary_entangled_inputs':False,
            'full_propagating_time_dependent_Maxwell_dynamics':False,'longitudinal_electric_backreaction':False,
            'interval_validated_long_time_nonlinear_trajectory_error':False,
            'numerical_step_refinement_is_a_continuum_error_certificate':False,
            'physical_electron_or_precision_g_prediction':False,'primitive_PW119_implementation_transfer':False},
        'verdict':'PASS_ACTUAL_QUASISTATIC_SPIN_AND_DIPOLE_PRECESSION'}
