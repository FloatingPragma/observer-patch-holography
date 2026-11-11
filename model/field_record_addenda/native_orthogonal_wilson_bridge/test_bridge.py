"""Exact identities, scope countercontrols, and independent numerical diagnostics."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import numpy as np
import pytest
import sympy as s
from . import algebra as a
from . import bridge
from . import vacuum


def test_actual_cube_characteristic_polynomial_and_fixed_charge_covariance():
    H,Js,minimum,powers,edges=a.cube()
    assert H.shape==(32,32) and minimum.degree()==12
    assert all(J**3==J for J in Js)
    result=a.covariance()
    assert result['exact_edge_charge_commutators']==96
    assert not result['canonical_field_and_onsite_charge_operators_changed']
    assert not result['outer_Gauss_projector_or_isometry_Y_changed']
    assert a.fw_identity()['exact_unitarity']
    assert a.rotation_algebra()['exact_free_projection_polynomial_identities_checked']


def test_all_256_charge_corners_and_neutral_countercontrol():
    result=bridge.coulomb_control();rows=result['exact_all_256_code_box_corners']
    assert len(rows)==256
    values=[F(row['Coulomb_over_m_g_squared'])for row in rows]
    assert max(values)==F(116,3)
    assert max(values)>F(59,4)>8
    assert result['all_native_physical_charge_Coulomb_norm_over_m_g_squared']=='464/3'
    # Exactly neutral box corners saturate the separate neutral bound8.
    assert max(F(row['Coulomb_over_m_g_squared'])for row in rows if sum(row['charges'])==0)==8
    assert bridge.edge_phase_control()['complete_four_mode_Fock_dimension']==16


def test_raw_sea_shift_and_complement_completion_are_not_optional():
    value=a.energy_alignment();lo,hi=map(F,value['signed_code_vacuum_shift_over_m'])
    nu=F(value['all_code_charge_perturbation_norm_over_m_upper'])
    assert 9<lo<hi<10 and nu==F(18,64)+F(116,3*64**2)
    assert lo-nu>8
    # No completion, or the beta source's tiny completion, leaves the native
    # complement below the entire perpendicular code.
    assert F(1,1000)<lo-nu
    c=F(value['robust_declared_positive_complement_coefficient'])
    assert c==F(21,2) and c-hi-nu>1
    assert F(value['full_physical_ground_gap_over_m_lower'])>F(1,2)
    assert value['old_code_scalar_minus_Esea_old_retained']
    assert not value['additional_code_scalar_offset_inserted']


def test_independent_filled_band_contractions_and_direct_rotation_diagnostic():
    H,Js,*_=a.cube();H=np.array(H,dtype=complex)
    energy,states=np.linalg.eigh(H);P=states[:,energy<0]@states[:,energy<0].conj().T
    C=np.array(a.field.cube()['Q_vacuum_covariance'],dtype=float)
    J=[np.array(x,dtype=complex)for x in Js]
    sigma=sum(C[e,f]*np.trace(P@J[e]@(np.eye(32)-P)@J[f])for e in range(12)for f in range(12))
    result=vacuum.calculate();lo,hi=map(F,result['new_pair_one_photon_current_square_over_m_squared'])
    assert abs(sigma.imag)<1e-12 and float(lo)<sigma.real<float(hi)
    for x in range(8):assert abs(np.trace(P[4*x:4*x+4,4*x:4*x+4])-2)<1e-12
    oldH,*_=a.matter.cube_matter();ew,ev=np.linalg.eigh(np.array(oldH,dtype=complex))
    Q=ev[:,ew<0]@ev[:,ew<0].conj().T
    M=np.eye(32)-(P-Q)@(P-Q);me,mv=np.linalg.eigh(M)
    V=(P@Q+(np.eye(32)-P)@(np.eye(32)-Q))@((mv/np.sqrt(me))@mv.conj().T)
    assert np.linalg.norm(V@Q@V.conj().T-P)<1e-12
    assert np.linalg.norm(V.conj().T@V-np.eye(32))<1e-12
    assert np.linalg.norm(P-Q,2)<float(F(result['exact_projector_distance_upper']))
    # Rotating the physical state is necessary; fixed old sea is distinct.
    assert np.linalg.norm(P-Q,2)>.5


def test_actual_native_pattern_and_full_raw_physical_ground():
    parent=bridge.parents();alignment=a.energy_alignment()
    controls=bridge.ground.local_completion_control()
    assert controls['actual_81_state_native_selected_port_carrier_projector_rank']==16
    result=bridge.source_gap(parent,alignment)
    mlo=F(result['source_selected_mass_interval'][0])
    assert F(result['full_physical_ground_gap_source_units_lower'])>mlo/2
    assert F(result['full_physical_ground_raw_energy_interval'][1])<0
    assert not result['native_Gauss_violating_kinematic_states_assigned_physical_ground_interpretation']
    assert not result['coupled_ground_vector_or_ground_preparation_numerically_constructed']


def test_complete_imperfect_preparation_record_and_clock():
    parent=bridge.parents();alignment=a.energy_alignment()
    result=bridge.interaction(vacuum.calculate(),parent);e=result['exact']
    assert 0<F(e['actual_code_probability_lower'])<1
    assert F(e['actual_original_attempt_target_probability_lower'])>F(1,10**8)
    assert F(e['actual_original_attempt_connected_covariance_lower'])>0
    assert not result['source_state_is_actual_coupled_ground']
    assert result['new_sea_count_projectors_differ_from_parent_count_projectors']
    r=bridge.resources(parent,alignment)
    assert r['count_pointer_dimension']==(17**2+1)*4+1
    assert r['complete_elapsed_source_time']=='2+1/(64m)'
    assert F(r['gross_two_pulse_switching_work_upper'])==4*F(r['each_control_norm_upper'])
    assert not r['controller_coefficient_workspace_in_capacity']


def test_cold_certificate_forbids_every_eigen_service():
    import os
    import subprocess
    code="""
import sys,runpy,numpy as np,scipy.linalg,sympy as s
from pathlib import Path
def forbidden(*args,**kwargs):raise AssertionError('eigen service invoked')
for module in (np.linalg,scipy.linalg):
 for name in ('eig','eigh','eigvals','eigvalsh'):setattr(module,name,forbidden)
for name in ('eigenvals','eigenvects','diagonalize'):setattr(s.MatrixBase,name,forbidden)
sys.argv=['producer','--check'];runpy.run_module('model.field_record_addenda.native_orthogonal_wilson_bridge.producer',run_name='__main__')
"""
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    result=subprocess.run([sys.executable,'-B','-c',code],cwd=Path(__file__).resolve().parents[3],
                          env=env,capture_output=True,text=True,check=False)
    assert result.returncode==0,result.stdout+result.stderr
    value=json.loads((Path(__file__).parent/'results.json').read_text())
    assert not value['scope']['interacting_relativistic_or_spatial_continuum_limit']
    assert not value['scope']['measured_SM_or_GeV_mass']
