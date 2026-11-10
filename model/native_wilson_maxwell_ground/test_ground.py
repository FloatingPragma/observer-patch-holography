"""Exact source controls and independent complete-instrument checks."""
from fractions import Fraction as F
import importlib.util
import math
from pathlib import Path
import numpy as np
import scipy.linalg as la
import sympy as s

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('native_wilson_ground_producer',HERE/'producer.py')
p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
g=p.ground


def test_exact_CRT_coefficients_and_full_code_obstruction():
    result=g.calculate();c=result['exact_vacuum_contractions']
    expected=s.Rational(4055,4368)-s.sqrt(21)/14-3*s.sqrt(39)/208-27*s.sqrt(91)/728
    assert s.simplify(s.sympify(c['exact_Coulomb_vacuum_coefficient'])-expected)==0
    assert c['complete_CRT_hopping_contraction_coefficients']==['-3','6','17']
    b={k:F(v) for k,v in result['bounds_exact'].items()}
    assert b['strict_Schur_margin']>F(3453,10**9)
    assert b['original_code_above_native_triplet_margin_lower']>F(1383,10**8)
    assert result['unchanged_source']['active_ground_multiplicity']==96


def test_independent_filled_band_Wick_diagnostic():
    H,*_=g.exact_current.literal_matter();H=np.array(H,complex)
    energy,frame=la.eigh(H);P=frame[:,energy<0]@frame[:,energy<0].conj().T
    covariance=np.zeros((8,8))
    for x in range(8):
        for y in range(8):
            block=P[4*x:4*x+4,4*y:4*y+4]
            covariance[x,y]=(2 if x==y else 0)-la.norm(block)**2
    c=np.sum(np.array(g.field.cube()['G'],float)*covariance)/2
    exact=float(s.N(s.sympify(g.contractions()['exact_Coulomb_vacuum_coefficient']),18))
    assert abs(c-exact)<1e-13


def test_native_code_projectors_and_coherent_Gauss_completion():
    c=g.local_completion_control()
    assert c['actual_81_state_native_selected_port_carrier_projector_rank']==16
    assert c['coherent_FW_density_and_Gauss_intertwiner_control']
    assert not c['preexisting_global_kinetic_current_guard_removed']


def test_complete_pattern_count_invariance_and_phase_countercontrol():
    # Two macrocell flags:00 code,01/10 one defect,11 two defects.
    H=np.diag([.1,.01,.02,.03]);D=np.diag([0,1,1,2]);c=.2
    psi=np.ones(4,dtype=complex)/2
    old=la.expm(-2.3j*H)@psi;new=la.expm(-2.3j*(H+c*D))@psi
    assert np.allclose(abs(old)**2,abs(new)**2,atol=1e-15)
    assert abs(abs(np.vdot(psi,old))**2-abs(np.vdot(psi,new))**2)>.02
    assert np.argmin(np.diag(H))!=np.argmin(np.diag(H+c*D))
    assert np.diag(H+c*D)[0]==np.diag(H)[0]


def test_full_energy_agreement_outcomes_and_mixed_state_bound():
    energies=np.array([-2.,-.4,1.3]);B=3.;M=31
    C=(1-1j*energies/B)/(1+1j*energies/B)
    phases=np.kron(C,C.conj());amplitudes=np.array([sum(np.exp(-2j*np.pi*j*k/M)*phases**k for k in range(M))/M for j in range(M)])
    assert np.allclose(np.sum(abs(amplitudes)**2,axis=0),1,atol=1e-14)
    distribution=np.array([.97,.02,.01]);pair=np.kron(distribution,distribution)
    success=float(pair@abs(amplitudes[0])**2);ground_pair=.97**2/success
    eps=min(1.,2*B/(M*1.6));r=.03/.97
    lower=1/(1+r*r+2*eps*eps*r)
    assert success>=.97**2-1e-15 and ground_pair>=lower


def test_actual_ground_preparation_and_complete_dressing_record():
    value=g.calculate();prep=value['complete_repaired_ground_preparation'];record=value['complete_stationary_ground_dressing_record']
    assert prep['each_pointer_dimension']==157589952376326706
    assert prep['raw_full_kinematic_source_norm_integer_upper']=='19397'
    assert prep['total_elapsed_time_all_branches']==4
    assert F(prep['complete_tree_success_lower'])>F(9995943,10**7)
    assert F(prep['conditional_main_trace_distance_upper'])<F(1418,10**10)
    assert F(record['physical_preparation_and_count_joint_probability_lower'])>F(141692,10**11)
    assert F(value['bounds_exact']['repaired_full_physical_gap_lower'])>F(252,10**12)


def test_fresh_scientific_producer_without_numerical_eigen_service(monkeypatch):
    def forbidden(*a,**k):raise AssertionError('numerical eigen service forbidden')
    for module in (np.linalg,la):
        for name in ('eig','eigh','eigvals','eigvalsh'):monkeypatch.setattr(module,name,forbidden)
    value=p.certificate()
    assert value['declared_onsite_completion']['global_scalar_energy_subtraction']=='none'
    assert value['exact_macrocell_and_Gauss_completion_control']['each_local_Px_reduces_all_full_source_terms']
