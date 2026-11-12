"""Native intertwiners, complete free-carrier records and failure controls."""
from fractions import Fraction as F
import numpy as np
import sympy as s
from scipy.linalg import expm
from .spectrum import (native_frames,dirac_hamiltonian,exact_massless_charge_control,
                        all_lattice_shells,spatial_records,complete_rest_records)


def test_native_exact_frames_and_global_mass_counterexample():
    value=native_frames();assert len(value['frames'])==3
    control=exact_massless_charge_control()
    assert control['exact_native_colour_generator_trace_checks']==8
    assert not control['total_charge_one_sector_has_positive_mass_floor']


def test_all_valleys_full_spectrum_and_complete_rest_records():
    value=all_lattice_shells();assert len(value['rows'])==64
    assert sum(row['kinetic_norm_squared']==0 for row in value['rows'])==8
    records=complete_rest_records();assert records['record_count']==864
    assert records['zero_probability_records']==720


def test_independent_dense_spatial_evolution_and_full_position_outputs():
    n=4;mass=.21;time=.7;momenta=np.arange(n)*2*np.pi/n
    points=list(__import__('itertools').product(range(n),repeat=3));volume=n**3
    transform=np.array([[np.exp(2j*np.pi*np.dot(x,p)/n)/np.sqrt(volume)for p in points]for x in points])
    from scipy.linalg import block_diag
    blocks=[];plus=[]
    for p in points:
        sinp=np.sin(np.array(p)*2*np.pi/n);sinp[abs(sinp)<1e-14]=0
        h=dirac_hamiltonian(sinp,mass);E=np.sqrt(sinp@sinp+mass*mass)
        blocks.append(h);plus.append((np.eye(4)+h/E)/2)
    change=np.kron(transform,np.eye(4));H=change@block_diag(*blocks)@change.conj().T
    P=change@block_diag(*plus)@change.conj().T
    seed=np.eye(4*volume)[:,0];actual=expm(-1j*time*H)@P@seed
    result=spatial_records(mass,time)
    assert np.max(abs(np.abs(actual.reshape(volume,4))**2-result['unconditional_all_site_spinor_arrival_probabilities']))<2e-14
    assert abs(np.vdot(actual,actual).real+result['reservoir_failure_probability']-1)<2e-14


def test_massive_shell_is_single_particle_not_pair_gap():
    from model.lepton_mass_compiler.compiler import _annihilator
    c=[_annihilator(2,j)for j in range(2)];m=.23;p=.71
    pair=m*(c[0].conj().T@c[1].conj().T);H=p*sum(a.conj().T@a for a in c)+pair+pair.conj().T
    e=np.linalg.eigvalsh(H);gap=np.sqrt(p*p+m*m)
    assert np.allclose(e-e[0],[0,gap,gap,2*gap],atol=1e-14)
    assert not np.isclose(e[-1]-e[0],gap)


def test_initial_projection_tail_is_not_a_propagation_claim():
    t0=spatial_records(.1,0);t1=spatial_records(.1,1)
    a=np.array(t0['unconditional_all_site_spinor_arrival_probabilities']);b=np.array(t1['unconditional_all_site_spinor_arrival_probabilities'])
    assert a[1:].sum()>0 and np.abs(a-b).sum()/2>1e-4
