"""Exact source-gap clock instruments and their onset/drive controls."""
from fractions import Fraction as F
from pathlib import Path
import json
import numpy as np
from scipy.linalg import expm,sqrtm
from model.openwave_common_mass_source.source import lepton_source,ix,INPUTS
from model.openwave_common_electromagnetism.source import exact_ward
ROOT=Path(__file__).resolve().parents[2]
CONSTANTS=json.loads((ROOT/'model/openwave_lepton_widths/inputs.json').read_text())
UNIT=F(INPUTS['electron_unit_GeV'])
H_SI=F(CONSTANTS['h_J_s_exact'])
GEV_J=10**9*F(CONSTANTS['elementary_charge_C_exact'])


def pair_operators():
    # Native source addition times a neutral fermionic apparatus creation:
    # only |00> and |11> are connected; all four Fock states remain.
    a=np.array([[0,1],[0,0]],complex);z=np.diag([1.,-1.])
    c=np.kron(a,np.eye(2));r=np.kron(z,a)
    N=c.conj().T@c;Nr=r.conj().T@r
    A=c.conj().T@r.conj().T;X=A+A.conj().T;Y=1j*(A-A.conj().T)
    parity=np.diag([1.,-1.,-1.,1.]);assert np.array_equal(parity@A,A@parity)
    assert np.array_equal(A.conj().T@A,np.diag([1.,0.,0.,0.]))
    assert np.array_equal(A@A.conj().T,np.diag([0.,0.,0.,1.]))
    # Common longitudinal W shifts its original quantum background charge.
    source_charge=-N;gauss_reference_charge=N
    assert np.array_equal((source_charge+gauss_reference_charge)@A,A@(source_charge+gauss_reference_charge))
    return N,Nr,A,X,Y


def records(phase,m=1.,reference_gap=0.,differential_drive=0.,dephase=False):
    N,Nr,A,X,Y=pair_operators();H=m*N+reference_gap*Nr+differential_drive*N
    rho=np.zeros((4,4),complex);rho[0,0]=1;time=phase/m
    first=[np.asarray(sqrtm((np.eye(4)+s*X)/2),complex)for s in (-1,1)]
    U=expm(-1j*H*time);out=[];first_energy=0.;last_energy=0.
    for s,K in zip((-1,1),first):
        state=K@rho@K.conj().T;first_energy+=np.trace(state@H).real
        if dephase:state=np.diag(np.diag(state))
        state=U@state@U.conj().T
        for quadrature,O in [('X',X),('Y',Y)]:
            for sign in (-1,1):
                effect=(np.eye(4)+sign*O)/2
                probability=.5*np.trace(effect@state).real
                L=np.asarray(sqrtm(effect),complex);post=.5*L@state@L.conj().T
                last_energy+=np.trace(post@H).real
                out.append({'first_sign':s,'quadrature':quadrature,'final_sign':sign,'probability':float(probability)})
    assert abs(sum(x['probability']for x in out)-1)<2e-14
    assert abs(first_energy-last_energy)<2e-14*max(1,m)
    effective=(m+reference_gap+differential_drive)*time
    for row in out:
        wave=0 if dephase else np.cos(effective) if row['quadrature']=='X' else np.sin(effective)
        # CAR order makes A|00> = |11>; this Y convention has negative sine.
        if row['quadrature']=='Y':wave=-wave
        expected=(1+row['first_sign']*row['final_sign']*wave)/8
        assert abs(row['probability']-expected)<2e-14
    return {'phase_in_mass_units':phase,'source_time_diagnostic':time,'complete_eight_records':out,
        'first_measurement_mean_work_electron_units_diagnostic':float(first_energy),
        'terminal_nonselective_mean_work_electron_units_diagnostic':float(last_energy-first_energy)}


def native_attachment():
    l=lepton_source();T=np.vstack((l['A'],l['B'].conj()));rows=[]
    for f,name in enumerate(('electron','muon','tau')):
        j=ix(0,0,11,f);mass=l['calibration']['charged_exact'][f]/l['calibration']['charged_exact'][0]
        frame=T[:,j];assert abs(np.vdot(frame,frame)-1)<2e-14
        # Exact charged Gaussian branch and literal complete native probe isometry.
        assert l['charge'][j]==-1
        gev=mass*UNIT;cycles=gev*GEV_J/H_SI
        rows.append({'species':name,'mode':j,'source_rest_energy_electron_units_exact':str(mass),
            'rest_energy_GeV_exact':str(gev),'clock_cycles_per_second_exact':str(cycles),
            'clock_cycles_per_second_diagnostic':float(cycles),
            'clock_angular_frequency_per_second_diagnostic':float(cycles)*2*np.pi,
            'native_probe_weights_96':(abs(frame)**2).tolist(),
            'all_rest_spin_and_charge_copies_retained':True,
            'unit_role':'empirical electron anchor and exact SI unit definitions; analyzer phase has retrospective ancestry'})
    return rows


def certificate():
    N,Nr,A,X,Y=pair_operators();m=F(1);H=N
    assert np.array_equal(H@A-A@H,A)
    rho=np.diag([1.,0,0,0]);stationary=[]
    for t in (0.,.3,1.7,11.):
        Ut=expm(-1j*H*t)
        correlation=np.trace(rho@Ut.conj().T@X@Ut@X)
        assert abs(correlation-np.exp(-1j*t))<1e-14
        stationary.append({'time':t,'one_point_X':float(np.trace(rho@Ut.conj().T@X@Ut).real),
            'two_time_correlation_real':float(correlation.real),'two_time_correlation_imag':float(correlation.imag)})
    normal=[records(t)for t in (0.,.3,1.,2.,3.,4.,7.)]
    dephased=records(1.3,dephase=True);drive=records(1.3,differential_drive=.2);reference=records(1.3,reference_gap=.4)
    common=.7*np.eye(4);assert np.array_equal(common@A-A@common,np.zeros((4,4)))
    # Native eigenstate is stationary: no fixed-observable one-point onset.
    assert all(r['one_point_X']==0 for r in stationary)
    assert all(abs(r['probability']-.125)<2e-14 for r in dephased['complete_eight_records'])
    return {'schema':'oph.native_two_time_rest_clock.v1','verdict':'PASS_NATIVE_GAP_CORRELATION_COMPLETE_CLOCK_AND_ONSET_CONTROLS',
        'source':'same192-CAR longitudinal matter member, colourful vacuum, complete one-positive-charged reducing sector, retained Gauss background and a neutral fermionic apparatus memory',
        'exact_native_charge_Ward':exact_ward(),
        'native_rest_clock_rows':native_attachment(),
        'operator_identities':{'A':'W a_j* P_vac W* tensor r*; bounded native addition and full neutral fermionic apparatus memory',
            'reference_Hamiltonian':'0 on both retained apparatus occupations; original Gauss compensator retained',
            'reference_energy_is_declared':True,'full_Fock_generator_reduced_only_by_proof':True,
            'commutator':'[H,A]=m_j A','stationary_two_time_correlation':'<X(t)X(0)>=exp(-i m_j t)',
            'one_point_no_onset':'[rho,H]=0 implies every fixed-observable expectation is constant',
            'source_mass_not_doubled':'frequency m_j; the negative-to-positive 2m_j beat is a different experiment'},
        'stationary_controls':stationary,'complete_sequential_records':normal,
        'dephased_record_control':dephased,'differential_drive_control':drive,'reference_energy_control':reference,
        'rigidity_domain':{'initial_spin_population_and_clock_first_outcome':'same frequency on the entire equal-mass spin doublet',
            'common_scalar_or_total_joint_charge_drive':'commutes with A, hence cannot entrain this correlation',
            'local_differential_scalar_drive_changes_frequency':True,
            'nonzero_reference_energy_changes_frequency':True,
            'arbitrary_drive_frequency_rigidity_proved':False},
        'energy_and_records':{'first_measurement_mean_source_work':'m_j/2, including its retained quantum reference',
            'subsequent_uncontrolled_evolution':'same autonomous source H; no continuing clock drive',
            'terminal_nonselective_mean_work':'0 for the stated equatorial inputs',
            'complete_kinematic_carrier':'full four-mode-occupation coordinates in the displayed reduction; every original native Fock and constraint complement retained',
            'failure_or_complement_outputs':'effects (I plus or minus X)/2 and (I plus or minus Y)/2 assign their explicit one-half effects on the complement'},
        'passive_ground_state_time_crystal_or_self_starting_one_point_clock':False,
        'scope':'The first native measurement creates the conditional coherence and costs m/2. This is an operational two-time rest correlation with declared apparatus reference, not passive one-point onset, arbitrary-drive rigidity, a source-independent laboratory clock, or a proof of the literal self-starting time-crystal criterion.'}
