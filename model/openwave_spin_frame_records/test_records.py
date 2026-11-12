from fractions import Fraction as F
import numpy as np
from scipy.linalg import expm
from .records import SIG,SIN,cube_certificate,occupied_rotation_phase,packet,frame_field_records,pointer_control,full_spin_records,tv
from model.openwave_lepton_mass_shell.spectrum import dirac_hamiltonian
from model.openwave_kinetic_phase_cycle.compiler import rounds_for_error,trace_bound


def test_literal_cube_and_source_ground_rotation():
    assert cube_certificate()['quarter_turns_to_source_parity']==4
    for m in (None,.0000753757,.01558546,.26211695):
        result=occupied_rotation_phase(m)
        assert abs(result['rotation_phase_real']-1)<1e-12


def test_native_packet_autonomous_evolution_and_spin():
    for m in (.0000753757,.01558546,.26211695):
        for axis in ('x','z'):
            initial,w,_=packet(m,0,axis)
            momentum=np.fft.fft(initial,axis=0)/2
            target=[]
            for k,v in zip(SIN,momentum):target.append(expm(-3j*dirac_hamiltonian((0,0,k),m))@v)
            direct=np.fft.ifft(target,axis=0)*2
            actual,weight,_=packet(m,3,axis)
            assert abs(weight-w)<1e-14 and np.linalg.norm(actual-direct)<3e-14
            assert np.linalg.norm(abs(initial)**2-abs(actual)**2)>1e-3
            if axis=='z':assert np.linalg.norm(actual@np.kron(np.eye(2),SIG[2]/2).T-actual/2)<1e-14


def test_complete_two_copy_effects_without_fine_sign_readout():
    m=.23;t=1.7;state,w,_=packet(m,t);v=state[0];p=np.vdot(v,v).real
    for a in range(3):
        for b in range(3):
            probs=[]
            for A in (SIG[a],SIG[b]):
                probs.append({u:(1-p if u==0 else np.vdot(v,np.kron(np.eye(2),(np.eye(2)+u*A)/2)@v).real)for u in (-1,0,1)})
            direct={z:sum(probs[0][u]*probs[1][v]for u in (-1,0,1)for v in (-1,0,1)if u*v==z)for z in (-1,0,1)}
            record=frame_field_records(m,t,0)
            chosen=[r for r in record['complete_two_copy_tensor_records']if r.get('axes')==[a,b]]
            for row in chosen:assert abs(row['probability']-w*w/9*direct[row['product_outcome']])<1e-15
    assert abs(sum(r['probability']for r in record['complete_frame_records'])-1)<1e-14


def test_periods_and_pointer_mutations():
    rows=[frame_field_records(.2,1.3,j)for j in range(9)]
    field=lambda j:rows[j]['complete_two_copy_tensor_records']
    frame=lambda j:rows[j]['complete_frame_records']
    w=rows[0]['preparation_success'];p=rows[0]['plane_success_conditional']
    assert tv(field(0),field(2))<1e-14 and tv(frame(0),frame(4))<1e-14
    assert abs(tv(field(0),field(1))-w*w*p*p/9)<1e-14
    assert abs(tv(frame(0),frame(2))-w*p/3)<1e-14
    # A vector substituted for the apolar tensor fails the claimed pi period.
    assert np.linalg.norm(np.array(rows[0]['native_plane_spin_current'])-rows[2]['native_plane_spin_current'])>.1
    controls=pointer_control(w);p0,pref,p4,_=controls['source_and_idle_reference_records']
    assert p0['complete_probabilities'][2]>.5 and pref['complete_probabilities'][1]>.5
    assert p4['complete_probabilities'][1]>.5
    assert abs(controls['pointer_dephasing_control_probabilities'][1]-w/2)<1e-15


def test_ground_spin_law_by_slater_characteristic_function():
    m=.3;p=np.array([1.,1.,-1.]);H=dirac_hamiltonian(p,m)
    _,B=np.linalg.eigh(H);B=B[:,:2];Sz=np.kron(np.eye(2),SIG[2]/2)
    a=(p[0]**2+p[1]**2)/(2*(p@p+m*m))
    for theta in (.1,.7,2.2):
        exact=np.linalg.det(B.conj().T@expm(1j*theta*Sz)@B)
        assert abs(exact-(1-a+a*np.cos(theta)))<2e-15
    full=full_spin_records([.0001,.01,.2]);assert full['ground_variance']>400
    assert not full['unfiltered_source_is_a_spin_eigenstate']


def test_complete_compiler_error_smaller_than_period_separation():
    n=rounds_for_error(F(1,10**12),F(3),sites=64)
    eps=trace_bound(n,F(3),sites=64)
    assert 2*eps<(1-F(263,1000))**4/4096
