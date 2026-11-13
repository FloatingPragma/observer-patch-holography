"""Exact complete-source identities and rational all-N information bounds."""
from fractions import Fraction as F
from math import comb, ceil
from pathlib import Path
import json
import numpy as np
from scipy.linalg import eigh
from model.openwave_collective_gravity_strength.source import coefficient_bounds, inherited_units
from model.openwave_universal_energy_geometry.source import QSIZE as Q, OMEGA as W, profile
from model.openwave_collective_finite_strength.bounds import root_bounds, down, up
from model.openwave_collective_finite_strength.source import complete_matrix

ROOT=Path(__file__).resolve().parents[2]


def domain(n,f=F(1),g=F(1)):
    if type(n) is not int or n<1 or not 0<=f<=1 or not 0<=g<=1:
        raise ValueError('positive integer count and profiles in [0,1] required')


def noise_bounds(n,f,g):
    """Cov(exp(fQ),exp(gQ))/(<exp(fQ)><exp(gQ)>), all outcomes."""
    f,g=F(f),F(g); domain(n,f,g)
    a=Q/n; v=Q*Q/n; x=f*g*v
    lower=x*(1-f*f*a*a/3)*(1-g*g*a*a/3)
    upper=x/(1-x)
    return lower,upper


def information_bounds(n):
    """Ground-state QFI for the continuous test source Es at Es=Ep=0."""
    domain(n);a=Q/n;v=Q*Q/n
    return 4*v*(1-a*a/3)**2/(W*W),4*v/(W*W*(1-v)**2)


def small_control(n):
    """Independent complete bit eigensystem and every coordinate outcome."""
    domain(n);H,q=complete_matrix(n,0,0);energy,V=eigh(H)
    vacuum=np.ones(2**n)/np.sqrt(2**n)
    A=np.exp(q); M=float(vacuum@(A*vacuum))
    covariance=float(vacuum@(q*A*vacuum))
    expected=float(Q)*np.cosh(float(Q)/n)**n*np.tanh(float(Q)/n)
    assert abs(covariance-expected)<2e-15
    amplitude=V.T@(A*vacuum)
    derivative=V[:,1:]@(-amplitude[1:]/energy[1:])
    assert abs(vacuum@derivative)<2e-13
    qfi=float(4*derivative@derivative)
    a=float(Q)/n
    exact_sum=4*np.cosh(a)**(2*n)*sum(comb(n,r)*np.tanh(a)**(2*r)/(r*float(W))**2 for r in range(1,n+1))
    assert abs(qfi-exact_sum)<2e-13
    lower,upper=information_bounds(n)
    assert float(lower)-2e-14<=qfi<=float(upper)+2e-14
    p=[];dp=[]
    for k in range(n+1):
        mask=np.array([n-j.bit_count()==k for j in range(2**n)])
        p.append(float(np.sum(vacuum[mask]**2)))
        dp.append(float(np.sum(2*vacuum[mask]*derivative[mask])))
    classical=sum(d*d/x for d,x in zip(dp,p,strict=True))
    assert abs(classical-qfi)<2e-13
    h=1e-5
    minus=eigh(complete_matrix(n,-h,0)[0])[1][:,0]
    plus=eigh(complete_matrix(n,h,0)[0])[1][:,0]
    # Fix the Perron phases before a finite-difference control.
    minus*=np.sign(minus@vacuum);plus*=np.sign(plus@vacuum)
    residual=float(np.linalg.norm((plus-minus)/(2*h)-derivative))
    assert residual<2e-8
    noise=[]
    for f,g in ((F(1),F(1)),(F(1,2),F(1,3)),(F(0),F(1))):
        Af=np.exp(float(f)*q);Ag=np.exp(float(g)*q)
        value=float(np.mean(Af*Ag)/(np.mean(Af)*np.mean(Ag))-1)
        lo,hi=noise_bounds(n,f,g)
        assert float(lo)-3e-15<=value<=float(hi)+3e-15
        noise.append({'f_exact':str(f),'g_exact':str(g),'interval_exact':[str(lo),str(hi)],'full_bit_control':value})
    quench=[]
    for time in (0.,1.,3.):
        value=float(16*np.cosh(a)**(2*n)*sum(comb(n,r)*np.tanh(a)**(2*r)*np.sin(r*float(W)*time/2)**2/(r*float(W))**2 for r in range(1,n+1)))
        assert value<=4*qfi+2e-14
        quench.append({'time':time,'QFI_at_zero_source':value})
    return {'N':n,'full_bit_dimension':2**n,'all_coordinate_probabilities':p,
        'all_coordinate_probability_derivatives':dp,'cov_Q_expQ':covariance,
        'exact_spectral_QFI':qfi,'complete_coordinate_classical_FI':classical,
        'finite_difference_ground_derivative_error':residual,'normalized_lapse_covariances':noise,
        'quench_information':quench,'coordinate_readout_work_source_units_exact':str(F(n)*W/2)}


def run():
    prior=json.loads((ROOT/'model/openwave_collective_gravity_strength/results.json').read_text())
    n=prior['calibrated_member']['registers'];u=inherited_units()[0]
    v=Q*Q/n;a=Q/n;lo,hi=information_bounds(n)
    work=F(n)*W*u/2;stiffness=profile()['stiffness_exact']
    count=ceil(1/hi)
    rows=[]
    for multiplier in (1,2,4):
        N=multiplier*n;glo,ghi=coefficient_bounds(N)
        rows.append({'N':N,'Q_variance_exact':str(Q*Q/N),
            'conditional_G_GeV_minus2_interval_exact':[str(glo/u**2),str(ghi/u**2)],
            'coordinate_work_GeV_exact':str(F(N)*W*u/2),
            'complete_free_gap_exact':str(W),'native_matter_source_changed':False})
    assert F(rows[1]['conditional_G_GeV_minus2_interval_exact'][1])<F(rows[0]['conditional_G_GeV_minus2_interval_exact'][0])
    rms=root_bounds(v)
    work_noise=work*v
    assert work_noise==Q*Q*W*u/2
    return {'schema':'oph.gravity_noise_normalization.v1',
      'verdict':'PASS_EXACT_RESPONSE_NOISE_WORK_AND_INFORMATION_IDENTITIES',
      'source':{'H_geometry':'omega/2 sum(I-X)','Q':'q/N sum Z','q_exact':str(Q),'omega_exact':str(W),
        'parameter_family':'H(Es)=H_geometry+Es exp(Q), Ep=0; information evaluated at Es=0',
        'initial_quench_state':'the complete zero-source geometry vacuum; no Es-dependent reference',
        'full_carrier':'all 2^N geometry states; matter energy labels are retained external controls'},
      'exact_identities':{'mean_Q':'0','variance_Q':'q^2/N',
        'cov_Q_expQ':'q cosh(q/N)^N tanh(q/N)',
        'dimensionless_Newton':'S Cov(Q,exp(Q))/(2 pi omega)',
        'physical_unit_identity':'G/(c_model/c)=g_N/u^2; speed matching is required to write G=g_N/u^2',
        'normalized_lapse_covariance':'[1+tanh(fq/N)tanh(gq/N)]^N-1',
        'work_times_Q_variance_GeV':'q^2 omega u/2',
        'ground_QFI':'4 cosh(q/N)^(2N) sum_{r=1}^N binomial(N,r)tanh(q/N)^(2r)/(r omega)^2',
        'complete_Q_classical_FI_equals_ground_QFI':True,
        'quench_QFI':'16 cosh(q/N)^(2N) sum_{r=1}^N binomial(N,r)tanh(q/N)^(2r)sin^2(r omega t/2)/(r omega)^2',
        'quench_QFI_upper':'4 times ground_QFI, every real time'},
      'conditional_calibrated_member':{'N':n,'variance_Q_exact':str(v),
        'rms_Q_interval_exact':list(map(str,rms)),
        'normalized_lapse_variance_interval_exact':list(map(str,noise_bounds(n,1,1))),
        'ground_QFI_interval_exact':[str(down(lo)),str(up(hi))],
        'single_copy_locally_unbiased_variance_lower_exact':str(down(1/hi)),
        'independent_ground_preparations_required_for_variance_at_most_one':count,
        'coordinate_work_GeV_exact':str(work),'work_times_Q_variance_GeV_exact':str(work_noise),
        'G_times_work_relative_to_S_q_squared_over_4pi_u_interval_exact':[str(1-a*a/3),str(1/(1-v/2))],
        'diagnostic':{'variance_Q':float(v),'rms_Q':float(rms[1]),'ground_QFI_upper':float(hi),
            'minimum_independent_ground_preparations':float(count),'coordinate_work_GeV':float(work)}},
      'common_spatial_mode':{'covariance_leading_term':'f*g*q^2/N',
        'all_profiles_share_one_random_Q':True,'independent_local_geometry_modes_added':False,
        'leading_covariance_matrix_rank_at_most_one':True},
      'normalization_counterfamily':rows,
      'complete_small_source_controls':[small_control(k) for k in (1,2,4,6,8)],
      'scope':{'observed_G_selects_inherited_N':True,'independent_Newton_prediction':False,
        'independent_covariance_measurement_would_fix_N':True,
        'hardware_uniqueness_required':False,
        'physical_metric_noise_identified':False,
        'information_bound_applies_to_specified_geometry_state_family_only':True,
        'arbitrary_interacting_observer_or_feedback_bound':False,
        'quantum_Cramer_Rao_requires_local_unbiasedness_and_independent_preparations':True,
        'common_strong_geometry_carrier_or_Einstein_continuum':False,'new_Nature_validation':False}}
