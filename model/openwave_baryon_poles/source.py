"""Exact residue bounds for native interpolating fields and complete records."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
import sympy as sy
from model.openwave_baryon_dispersion import source as bands
from model.openwave_hadron_refinement_flow.analytic import constants
from model.openwave_hadron_refinement_flow.parameters import degree_coefficients
from model.openwave_hadron_refinement_flow.certified_rest import sqrt_upper
from model.openwave_capacity_stress.source import enclosure

HERE=Path(__file__).resolve().parent
PHASES=((F(1),F(0)),(F(63,65),F(16,65)),(F(35,37),F(12,37)),(F(3,5),F(4,5)),(F(0),F(1)),(F(-1),F(0)))


def key(proton,cosine,sine):
    return ('proton' if proton else 'neutron')+':'+str(cosine)+':'+str(sine)


@lru_cache(None)
def exact_symbol(cosine,sine,proton=False):
    A,G=bands.matrix(cosine,sine,proton)
    if sine:
        A=A[:5,:5]+sy.I*A[5:,:5]
        G=G[:5,:5]
    J,_,_=bands.repair_matrix(proton)
    A+=2*(1-sy.Rational(cosine))*J
    assert A==A.conjugate().T and G==sy.diag(6,1,1,1,1)
    return A,G


def vector(record):
    D=record['dyadic_denominator']
    assert type(D)is int and D==2**40
    parts=[]
    for name in ('real_integer_coefficients','imaginary_integer_coefficients'):
        row=record[name]
        assert type(row)is list and len(row)==5
        assert all(type(x)is str and str(int(x))==x and abs(int(x))<=D for x in row)
        parts.append([sy.Rational(int(x),D) for x in row])
    return sy.Matrix([a+sy.I*b for a,b in zip(*parts)])


def real_fraction(value):
    value=sy.expand(value)
    assert sy.im(value)==0 and sy.re(value).is_Rational
    return F(sy.re(value))


def trial_moments(A,G,v):
    n=real_fraction((v.conjugate().T*G*v)[0]);assert n>0
    mean=real_fraction((v.conjugate().T*A*v)[0])/n
    residual=A*v-sy.Rational(mean)*G*v
    variance=real_fraction((residual.conjugate().T*G.inv()*residual)[0])/n
    assert variance>=0
    return n,mean,variance


def residue_bounds(compact,complement_floor,coupling_square,schur_error,mean,variance):
    e0,e1=compact[0],compact[1]
    C,b2,eps=map(F,(complement_floor,coupling_square,schur_error))
    assert b2>=0 and eps>=0 and C>e0[1]
    assert e0[1]<e1[0]-eps
    assert e0[0]<=mean<e1[0] and variance>=0
    q=b2/(C-e0[1])**2
    s=eps/(e1[0]-e0[1]-eps)
    tau=sqrt_upper(variance)/(e1[0]-mean)
    assert q<1 and s+tau<1
    lower=(1-q)*(1-(s+tau)**2)
    assert 0<lower<=1
    return {'complete_ground_gap_lower_source':str(e1[0]-eps-e0[1]),
        'full_Q_probability_upper':str(q),'P_ground_sine_upper':str(s),
        'prepared_trial_sine_upper':str(tau),
        'prepared_trial_pole_weight_lower':str(lower),
        'prepared_trial_other_spectral_weight_upper':str(1-lower),
        'five_field_pole_residue_trace_interval':[str(1-q),'1'],
        'original_input_pole_event_probability_lower':str(lower/5)}


@lru_cache(None)
def phase(proton,cosine,sine):
    row=bands.phase_record(cosine,sine,True,proton)
    A,G=exact_symbol(cosine,sine,proton)
    witnesses=json.loads((HERE/'trial_witnesses.json').read_text())
    v=vector(witnesses[key(proton,cosine,sine)])
    n,mean,variance=trial_moments(A,G,v)
    intervals=[tuple(map(F,r))for r in row['all_five_compact_source_intervals']]
    _,_,_,p,_,unit=constants(*degree_coefficients())
    certificate=residue_bounds(intervals,F(row['complete_high_complement_source_floor']),5*p['t']**2,F(row['uniform_Schur_error_source']),mean,variance)
    return {'assignment':'proton'if proton else'neutron','cosine':str(cosine),'sine':str(sine),
        'source_pole_GeV_interval':row['all_five_full_source_GeV_intervals'][0],
        'source_pole_GeV_display':enclosure(tuple(map(F,row['all_five_full_source_GeV_intervals'][0]))),
        'all_five_low_source_GeV_intervals':row['all_five_full_source_GeV_intervals'],
        'entire_high_complement_source_floor':row['complete_high_complement_source_floor'],
        'trial_source_coordinates':witnesses[key(proton,cosine,sine)],
        'trial_Gram_norm_squared':str(n),'trial_source_energy':str(mean),'trial_source_variance':str(variance),
        'trial_GeV_energy':str(unit*mean),'pole':certificate}


def run():
    rows=[phase(proton,c,s)for proton in (False,True)for c,s in PHASES]
    zeros=[r for r in rows if r['cosine']=='1']
    return {'schema':'oph.native_baryon_source_poles.v1',
        'verdict':'PASS_COMPLETE_NATIVE_SOURCE_POLES_AND_PREPARED_RESIDUES',
        'momentum_fibre_records':rows,'rest_source_poles':zeros,
        'matrix_spectral_measure':{'interpolators':'Complete five-mode native CAR triquark frame in the fixed spin/isospin/charge fibre',
            'vacuum_energy_exact':'0','zeroth_moment_exact':'I_5','isolated_ground_residue_rank':1,
            'resolvent':'G(z)=R_0/(z-E_0)+integral_complement dM(E)/(z-E)',
            'spectral_complement_discarded':False,'number_sector_dynamically_reducing':True},
        'complete_observer_instrument':{'input':'P/5 on the complete compact five-mode source, normalized momentum fibre or wave packet',
            'actual_prepared_state':'Fixed rational complex native trial, with its exact Gram normalization and certified compact-ground angle',
            'preparation_success_exact':'1/5','preparation_failure_exact':'4/5',
            'filter':'K_success=|trial><trial|; K_failure=I-K_success on the entire source',
            'nonselective_preparation_source_work_exact':'0 on the declared P/5 input',
            'branch_conditioned_preparation_work_claimed_zero':False,
            'terminal_outcomes':'isolated ground spectral event and its full orthogonal complement, on each preparation branch',
            'terminal_spectral_source_work_exact':'0','arbitrary_idle_external_quantum_reference':True},
        'input_ancestry':'Inherited electron anchor, seven-input strong source and neutron/axial bulk matching; supplied strong length and physical speed in the positive correlated kinetic rule.',
        'mass_refit_performed':False,'source_dimension':'Infinite one-dimensional native chain at fixed spacing, with all relative-coordinate complements.',
        'physical_3plus1_Lorentz_particle_pole_identification':False,
        'scope':'Isolated poles of the declared interacting strong-source momentum-resolved two-point functions. These are not a proof of a relativistic continuum particle dictionary or observed nucleon scattering poles.'}
