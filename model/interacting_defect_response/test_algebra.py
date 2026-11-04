from fractions import Fraction as F

import numpy as np
import pytest

from . import algebra,finite


def test_exact_native_symmetry_characteristic_polynomial_and_nonlinearity():
    result = algebra.exact_certificate()
    assert result['native_currents_commuting_with_core_number_interaction'] == 12
    assert result['core_number_conserved_by_pairing'] is False


@pytest.mark.parametrize('mass,coupling',[(1,0),(1,.5),(.1,2),(2,.1)])
def test_independent_full_fock_spectra_and_transition_weights(mass,coupling):
    finite.spectrum_control(mass,coupling)


def test_four_mode_block_is_exact_reducing_sector_with_singly_occupied_spectator():
    c = finite.numpy_car(6)
    number = sum((a.T@a for a in c),np.zeros((64,64)))
    pairs = [c[i].T@c[i+1].T+c[i+1]@c[i] for i in (0,2,4)]
    h = pairs[0]+pairs[1]+.37*pairs[2]+.5*(number-3*np.eye(64))@(number-3*np.eye(64))
    inclusion = np.eye(64)[:,[4*i+2 for i in range(16)]]
    expected,_,_ = finite.numerical_block(1,.5)
    assert np.linalg.norm(h@inclusion-inclusion@expected) < 1e-12


def test_nonzero_connected_density_is_not_a_gaussian_cross_contraction():
    h,_,c = finite.numerical_block(1,.5)
    _,u = np.linalg.eigh(h)
    g = u[:,0]
    assert abs(g@c[0].T@c[2]@g) < 1e-12
    assert abs(g@c[0]@c[2]@g) < 1e-12
    n0,n2 = c[0].T@c[0],c[2].T@c[2]
    cumulant = g@n0@n2@g-(g@n0@g)*(g@n2@g)
    assert abs(cumulant+.5/(4*np.hypot(1,.5))) < 1e-12
    assert abs(cumulant) > .1


def test_onsite_interaction_fails_same_core_invariance():
    result = finite.onsite_control()
    assert float(result['unit_onsite_leakage']) > .02


def test_weighted_onsite_leakage_identity():
    _,w,_,row = finite.construction.radial()
    e = np.kron(w.toarray(),w.toarray())
    sites = np.repeat(np.arange(len(row['density'])),2)
    weights = np.linspace(-1,1,len(row['density']))
    v = np.where(sites[:,None] == sites[None,:],weights[sites][:,None],0).ravel()
    applied = v[:,None]*e
    compressed = e.conj().T@applied
    leak = applied-e@compressed
    rho2 = row['density']**2
    expected = weights**2@rho2-(weights@rho2)**2
    assert np.linalg.norm(leak.conj().T@leak-expected*np.eye(4)) < 1e-12


def test_source_mass_blocks_are_retained_and_interacting_gaps_change():
    rows = finite.source_controls()['lepton_blocks']
    assert len(rows) == 3
    assert all(float(row['single_fermion_gaps'][0]) != float(row['mass']) for row in rows)


def test_radius_uniform_bound_and_exponential_localization():
    a = algebra.radius_uniform_bound(1,F(1,128),F(1,12),1,100)
    b = algebra.radius_uniform_bound(1,F(1,128),F(1,12),1,200)
    assert F(a['arbitrary_initial_CQ_error_upper']) == 2*F(b['arbitrary_initial_CQ_error_upper'])
    assert F(algebra.localization_bound(3)['one_particle_tail_probability_upper']) == F(algebra.localization_bound(2)['one_particle_tail_probability_upper'])/1024


@pytest.mark.parametrize('args',[(1.,1,F(1,12),1,10),(1,-1,F(1,12),1,10),(1,1,1,1,10),(1,1,0,1,0)])
def test_invalid_bound_domains_fail(args):
    with pytest.raises(ValueError):
        algebra.radius_uniform_bound(*args)


def test_complete_nonlinear_memory_control():
    result = finite.memory_control()
    assert result['retained_fine_records_per_input'] == 512


def test_full_history_square_root_bound_is_sharp_for_symmetric_phase_records():
    angle = .31
    psi = np.array([1,1])/np.sqrt(2)
    rho = np.outer(psi,psi)
    outputs=[]
    for sign in (-1,1):
        u = np.diag(np.exp(sign*1j*angle*np.array([1,-1])))
        outputs.append(u@rho@u.conj().T)
    unread = sum(outputs)/2
    delta = np.sum(np.abs(np.linalg.eigvalsh(unread-rho)))
    retained = sum(np.sum(np.abs(np.linalg.eigvalsh(out-rho))) for out in outputs)/2
    assert abs(retained**2-2*delta)<1e-12
    assert retained>delta
