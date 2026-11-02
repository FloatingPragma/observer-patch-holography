import numpy as np
import pytest
from scipy.linalg import eigh
from . import algebra as a, experiment as e
from model.strong_coupling_refinement.algebra import exchange_four


def test_actual_full_domain_repair_has_protected_internal_spin():
    r = a.protected_columns()
    assert r['protected_cube_columns'] == 12*2*2*256
    assert r['small_full_Fock_columns'] == 2*3*2*64


@pytest.mark.parametrize('length,area', [(4, 1), (6, 1), (2, 3), (4, 3)])
def test_twist_translation_minus_sign_on_complete_sz_zero_space(length, area):
    assert a.twist_translation(length, area)['twist_translation_factor'] == '-1'


@pytest.mark.parametrize('law,adaptive', [('C', False), ('E', False), ('C', True)])
def test_complete_choi_history_including_rare_collision_records(law, adaptive):
    r = e.recorded_history(law=law, adaptive=adaptive)
    assert abs(r['complete_probability']-1) < 2e-10
    assert 0 < r['trace_norm_difference'] < r['analytic_trace_bound']
    assert r['nonzero_raw_probability'] > 0
    assert r['zero_children'] > 0


@pytest.mark.parametrize('sites', [4, 6, 8])
def test_twist_gap_control_against_unprojected_spin_diagonalization(sites):
    edges = tuple((x, x+1, 1) for x in range(sites-1))+((0, sites-1, 1),)
    H = exchange_four(sites, edges).toarray()/4
    energy, vector = eigh(H)
    ground = vector[:, 0]
    phase = np.exp(2j*np.pi*np.array([sum(x*(.5-((s >> x) & 1)) for x in range(sites)) for s in range(1 << sites)])/sites)
    plus, minus = phase*ground, phase.conj()*ground
    assert abs(np.vdot(ground, plus)) < 1e-10
    shift = ((np.vdot(plus, H@plus)+np.vdot(minus, H@minus))/2).real-energy[0]
    assert energy[1]-energy[0] <= shift+1e-12 <= np.pi**2/sites+1e-12
