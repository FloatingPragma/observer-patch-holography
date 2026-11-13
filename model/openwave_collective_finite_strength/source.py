"""Finite-energy collective equilibrium and complete autonomous quench records."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigh
from model.openwave_energy_lapse.source import native
from model.openwave_universal_energy_geometry.source import QSIZE, OMEGA, profile
from .bounds import exp_bounds, root_bounds, cos_bounds, moments, complete_gap, down, up

ROOT = Path(__file__).resolve().parents[2]
Q, W = QSIZE, OMEGA


def uniform_certificate(n, total_energy, minimum_source):
    """All profiles and all positive source/probe energies in the stated simplex."""
    if type(n) is not int or n < 1 or total_energy <= 0 or minimum_source <= 0:
        raise ValueError('positive count and energy domain required')
    M, emin = F(total_energy), F(minimum_source)
    eq = exp_bounds(Q)[1]
    x = 2*M*Q/(n*W)
    mu2, mu4 = moments(n,Q,x*x/(1+x*x))
    gap = W-M*eq*mu2/2
    if gap <= 0:
        raise ValueError('this uniform product comparison does not certify the domain')
    epsilon = M*eq*root_bounds(mu4)[1]/(2*gap)
    if epsilon >= 1:
        raise ValueError('the stated angle bound is not useful on this domain')
    inv_root_n = root_bounds(F(1,n))[1]
    # A-I, rather than A, makes the quadratic observable error proportional to F.
    state_error_per_profile = 2*eq*Q*(epsilon*inv_root_n+epsilon*epsilon)
    k0 = 2*Q*Q/(n*W)
    product_error_per_profile = k0*M*x*x/2+F(n*(n-1),2)*x*x*(Q/n)**2
    C = state_error_per_profile+product_error_per_profile
    vacuum_upper = 1/(1-Q*Q/(2*n))
    relative_mutual = (vacuum_upper-1)+vacuum_upper*Q*Q/(3*n*n)+2*vacuum_upper*C/(k0*emin)
    return {'total_energy_upper_exact':str(M),'positive_source_energy_lower_exact':str(emin),
        'full_geometry_sector_gap_lower_exact':str(down(gap)),
        'ground_angle_upper_exact':str(up(epsilon)),
        'normalized_lapse_linear_remainder_per_profile_upper_exact':str(up(C)),
        'mutual_energy_relative_error_upper_exact':str(up(relative_mutual)),
        'leading_mutual_energy':'-2*q^2*F*Es*Ep/(N*omega)',
        'domain':'Es>=minimum_source, Ep>=0, Es+Ep<=total_energy; every F in[0,1]',
        'complete_ground_source_and_probe_included':True,
        'diagnostic':{'gap_lower':float(gap),'angle_upper':float(epsilon),
                     'lapse_remainder_per_profile':float(C),'mutual_relative_error':float(relative_mutual)}}


def quench_certificate(n, total_energy, time):
    """Duhamel comparison with exact product dynamics, without ground preparation."""
    M,t = F(total_energy),F(time)
    if type(n) is not int or n<1 or M<0 or t<0:
        raise ValueError('positive count and nonnegative energy/time required')
    eq=exp_bounds(Q)[1]
    max_m2=min(F(1),(4*M*Q/(n*W))**2)
    _,mu4=moments(n,Q,max_m2)
    epsilon=min(F(1),t*M*eq*root_bounds(mu4)[1]/2)
    observable=2*eq*Q*(epsilon*root_bounds(F(1,n))[1]+epsilon*epsilon)
    return {'time_exact':str(t),'total_energy_upper_exact':str(M),
        'state_angle_sine_upper_exact':str(up(epsilon)),
        'normalized_lapse_error_per_profile_upper_exact':str(up(observable)),
        'complete_Q_outcome_total_variation_upper_exact':str(up(epsilon)),
        'product_spin_mean':'-2*h*omega/(omega^2+4*h^2)*(1-cos(sqrt(omega^2+4*h^2)*t)); h=q*(Es+F*Ep)/N',
        'coordinate_probabilities':'binomial(N,k)*((1+m(t))/2)^k*((1-m(t))/2)^(N-k), every k=0,...,N',
        'arbitrary_initial_geometry_state':False,
        'initial_geometry':'complete Hg vacuum |+>^N; source/probe energy labels retained',
        'diagnostic':{'state_angle_sine_upper':float(epsilon),'lapse_error_per_profile_upper':float(observable)}}


def coarse_record_certificate(n, total_energy, minimum_source, time=F(355,113)):
    """Two independent complete preparations, with equal sampling weights."""
    M,emin,t=map(F,(total_energy,minimum_source,time))
    c=cos_bounds(W*t)
    timing_bias=(1+c[1])/2
    x=2*M*Q/(n*W)
    linear_factor_error=timing_bias+x*x*(1+W*t/4)
    k0=2*Q*Q/(n*W)
    qcert=quench_certificate(n,M,t)
    dynamical=F(qcert['normalized_lapse_error_per_profile_upper_exact'])/2
    product=F(n*(n-1))*x*x*(Q/n)**2
    total=dynamical+product+k0*M*linear_factor_error
    relative=Q*Q/(3*n*n)+2*total/(k0*emin)
    vacuum_upper=1/(1-Q*Q/(2*n))
    return {'two_observation_delays_exact':['0',str(t)],'weights_exact':['1/2','1/2'],
        'two_independent_vacuum_preparations':True,'all_outcomes_and_both_contexts_retained':True,
        'positive_source_energy_lower_exact':str(emin),'total_source_probe_energy_upper_exact':str(M),
        'source_on_minus_off_normalized_lapse':'-2*q^2*F*Es/(N*omega) times (1+bounded relative error)',
        'relative_error_upper_exact':str(up(relative)),
        'timing_factor_error_upper_exact':str(up(timing_bias)),
        'per_preparation_readout_work_interval_exact':[str(down(F(n)*W/2-M*vacuum_upper)),str(F(n)*W/2)],
        'readout_lapse_normalization':'the complete zero-matter geometry-vacuum expectation cosh(Fq/N)^N',
        'supplied_timing_and_blank_pointer_resources':True,
        'exact_source_coupling_evolution_no_equilibrium_preparation':True,
        'diagnostic_relative_error_upper':float(relative)}


def complete_matrix(n,es,ep=0.,f=.5):
    dimension=2**n
    coordinate=np.array([float(Q)*(n-2*k.bit_count())/n for k in range(dimension)])
    H=np.diag(float(W)*n/2+es*np.exp(coordinate)+ep*np.exp(f*coordinate))
    for k in range(dimension):
        for bit in range(n): H[k,k^(1<<bit)]-=float(W)/2
    return H,coordinate


def finite_controls():
    rows=[]
    for n in (1,2,4,6):
        for es,ep,f in ((F(0),F(0),F(1,2)),(F(1),F(2),F(1,2)),(F(20),F(1),F(1))):
            terms=[(F(1),(es,es)),(f,(ep,ep))]
            cert=complete_gap(n,terms)
            H,q=complete_matrix(n,float(es),float(ep),float(f))
            energies,V=eigh(H)
            assert float(cert['ground'][0])-2e-12<=energies[0]<=float(cert['ground'][1])+2e-12
            actual_gap=energies[1]-energies[0]
            assert float(cert['gap'][0])-2e-12<=actual_gap<=float(cert['gap'][1])+2e-12
            assert cert['gap'][0]>0
            probabilities=[]
            initial=np.ones(2**n)/np.sqrt(2**n)
            for time in (F(0),F(1),F(3)):
                psi=V@(np.exp(-1j*float(time)*energies)*(V.T@initial))
                p=abs(psi)**2
                outcomes=[float(sum(p[k] for k in range(2**n) if n-k.bit_count()==j)) for j in range(n+1)]
                assert abs(sum(outcomes)-1)<2e-12
                energy=float(np.vdot(psi,H@psi).real)
                assert abs(energy-float(initial@H@initial))<2e-11
                work=float(p@np.diag(H))-energy
                probabilities.append({'time':float(time),'all_N_plus_one_coordinate_probabilities':outcomes,
                    'coordinate_values':[float(Q)*F(2*j-n,n) for j in range(n+1)],
                    'energy_conservation_residual':energy-float(initial@H@initial),
                    'coordinate_measurement_work':work,
                    'mean_lapse':float(p@np.exp(float(f)*q)),
                    'certificate':quench_certificate(n,es+ep,time)})
            rows.append({'N':n,'Es_exact':str(es),'Ep_exact':str(ep),'profile_exact':str(f),
                'full_geometry_dimension':2**n,'full_gap_interval_exact':list(map(str,cert['gap'])),
                'full_ground_interval_exact':list(map(str,cert['ground'])),
                'independent_full_bit_matrix_gap':float(actual_gap),'complete_quench_records':probabilities})
    return rows


def run():
    parent=json.loads((ROOT/'model/openwave_collective_gravity_strength/results.json').read_text())
    n=parent['calibrated_member']['registers']
    gaps=native()[2]
    emin=gaps[1].lo
    assert 2*gaps[-1].hi<7
    cert=uniform_certificate(n,F(7),emin)
    broad=uniform_certificate(n,F(10**6),emin)
    assert F(cert['mutual_energy_relative_error_upper_exact'])<F(1,10**18)
    shape=profile(); central=[]
    for distance in (1,2,3):
        idx=np.ravel_multi_index((3+distance,3,3),(7,7,7))
        central.append(F(int(shape['F_integer'][idx]),2**40))
    force_factor=(central[0]+central[2])/(central[0]-central[2])
    force_error=force_factor*F(cert['mutual_energy_relative_error_upper_exact'])
    assert force_error<1
    pi_time=F(355,113)
    sampled=coarse_record_certificate(n,F(7),emin,pi_time)
    assert F(sampled['relative_error_upper_exact'])<F(1,10**11)
    return {'schema':'oph.collective_finite_strength.v1',
      'verdict':'PASS_COMPLETE_FINITE_ENERGY_MUTUAL_FORCE_AND_AUTONOMOUS_QUENCH_BOUNDS',
      'source':{'H':'(omega/2)sum(I-X)+Es exp(Q)+Ep exp(F Q)',
        'Q':'q/N sum Z','q_exact':str(Q),'omega_exact':str(W),
        'full_dimension':'64*2^N for all eight source and eight probe energies',
        'native_energy_intervals':[x.encoded() for x in gaps],
        'vacuum_and_all_49_positive_native_pairs_retained':True},
      'two_Jacobi_theorem':{'spin_multiplicity':'binomial(N,N/2-j)-binomial(N,N/2-j-1)',
        'ground_block':'j=N/2 by connected-hypercube Perron positivity',
        'full_gap_for_N_at_least_two':'min(E1(j=N/2),E0(j=N/2-1))-E0(j=N/2)',
        'N_equals_one':'only j=1/2 exists','no_complementary_spin_sector_dropped':True,
        'diagonal':'omega*N/2+Es*exp(2q*m/N)+Ep*exp(2q*F*m/N)',
        'offdiagonal_squared':'(omega^2/4)*(j-m)*(j+m+1)',
        'certificate':'outward rational exponential bounds, diagonal ordering and exact Sturm inertia'},
      'inherited_calibrated_registers':n,
      'native_finite_energy_certificate':cert,
      'broader_declared_energy_domain_certificate':broad,
      'complete_native_force':{'positive_native_source_probe_pairs':49,
        'three_profile_values_exact':list(map(str,central)),
        'central_force_relative_error_upper_exact':str(up(force_error)),
        'force_is_attractive_for_every_positive_native_pair':True,
        'source_energy_factorization_holds_with_quantified_remainder':True,
        'geometry_ground_preparation_required_for_equilibrium_force':True},
      'autonomous_source_quench':[quench_certificate(n,F(7),t) for t in (F(1),F(3),pi_time)],
      'two_time_coarse_record':sampled,
      'readout_resources':{'initial_energy':'Es*cosh(q/N)^N+Ep*cosh(F*q/N)^N',
        'complete_coordinate_measurement_work':'omega*N/2 - <Hg> before measurement',
        'vacuum_readout_work_source_units_exact':str(F(n)*W/2),
        'huge_carrier_enumerated':False},
      'complete_small_carrier_controls':finite_controls(),
      'scope':{'finite_native_source_and_probe_energies_not_only_zero_derivatives':True,
        'static_force_and_autonomous_quench_distinguished':True,
        'geometry_sector_gap_is_not_the_full_matter_spectrum_gap':True,'independent_G_prediction':False,
        'N_inherited_from_observed_G_calibration':True,'common_strong_and_geometry_source':False,
        'physical_Einstein_or_quantum_continuum':False,'new_Nature_validation':False}}
