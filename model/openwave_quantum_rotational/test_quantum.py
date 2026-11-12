import json
import numpy as np
from scipy.linalg import eigh
from .source import source,photon,cross,radial_effect
from .experiment import full_fock_order_control,complete_records,quantum_dynamics_control
from .enclosure import certify,rational_bounds


def test_exact_dressed_Casimir_includes_nonzero_quantum_field_share():
    d=source(2);v=d['psi']
    assert np.linalg.norm(d['J2']@v-.75*v)<1e-12
    assert np.linalg.norm(d['J'][2]@v-.5*v)<1e-12
    assert np.vdot(v,d['Jf'][2]@v).real>1e-6
    assert np.linalg.norm(d['H']@d['Jm'][0]-d['Jm'][0]@d['H'])>1e-3
    assert np.linalg.norm(d['H']@d['J'][0]-d['J'][0]@d['H'])<1e-11


def test_operator_valued_full_Fock_lift_rejects_wrong_block_order():
    d=full_fock_order_control()
    assert d['correct_order_minimum']>=-1e-14
    assert d['transposed_field_block_mutation_minimum']<-.49


def test_Maxwell_highest_shell_completion_is_required():
    for cap in(1,2,3):
        f=photon(cap);dim=len(f['basis']);eta=1/8
        raw=sum((eta**2*p@p+q@q/(4*eta**2)for p,q in zip(f['pi'],f['q'])),np.zeros((dim,dim),complex))-1.5*np.eye(dim)
        assert np.linalg.norm(raw+(cap+3)*f['top']/2-f['N'])<1e-13
        assert np.linalg.norm(raw-f['N'])>1


def test_cartesian_individual_caps_break_rotation_algebra():
    basis=[(i,j,k)for i in range(2)for j in range(2)for k in range(2)];where={p:i for i,p in enumerate(basis)};a=[]
    for axis in range(3):
        v=np.zeros((8,8),complex)
        for j,p in enumerate(basis):
            if p[axis]:q=list(p);q[axis]-=1;v[where[tuple(q)],j]=1
        a.append(v)
    J=[sum((-1j*cross(i,j,k)*a[j].conj().T@a[k]for j in range(3)for k in range(3)),np.zeros((8,8),complex))for i in range(3)]
    assert np.linalg.norm(J[0]@J[1]-J[1]@J[0]-1j*J[2])>1


def test_complete_filter_position_field_records_and_readout_work():
    r=complete_records(2)
    assert len(r['complete_preparation_position_photon_records'])==7
    assert abs(sum(x['original_input_probability']for x in r['complete_preparation_position_photon_records'])-1)<1e-12
    assert r['energy_parts']['photon_number_measurement_work']>1e-6
    ball=radial_effect(5);assert np.linalg.eigvalsh(ball)[0]>0 and np.linalg.eigvalsh(ball)[-1]<1
    assert float(rational_bounds(2)['all_state_outside_radius_five_upper'])<1.1e-7


def test_reference_uniform_unitary_and_opposite_torques():
    r=quantum_dynamics_control()
    assert r['unitary_norm_error']<1e-12
    assert r['total_angular_drift']<1e-12
    assert r['frozen_field_torque_mutation']>1e-5


def test_interval_energy_current_and_full_band_gap():
    r=certify(1);d=source(1)
    assert json.loads(json.dumps(r))==r
    assert r['actual_ground_energy_interval'][0]<=d['E']<=r['actual_ground_energy_interval'][1]
    g=4*d['E']*abs(np.vdot(d['psi'],d['current_derivative']@d['psi']).real)
    assert r['rest_energy_normalized_static_g_interval'][0]<=g<=r['rest_energy_normalized_static_g_interval'][1]
    assert r['exact_dressed_line_trace_distance_upper']<1e-10
    assert d['all_energies'][2]>=17
