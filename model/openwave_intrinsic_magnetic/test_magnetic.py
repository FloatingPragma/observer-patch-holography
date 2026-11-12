"""Complete-kernel derivatives and stationary source/field controls."""
import numpy as np
import pytest
from scipy import linalg as la
from .solver import native_reduction,dirac,solve,maxwell,source_current
from .checks import derivative_controls,gauge_control,external_probe,dipole_profile


@pytest.fixture(scope='module')
def stationary():
    return solve()


def test_all_internal_rows_and_every_gauged_mass_gap():
    result=native_reduction()
    assert result['native_internal_rows']==16 and result['retained_folded_valley_labels']==8
    rng=np.random.default_rng(844)
    for amplitude in (0,.2,2.):
        field=amplitude*rng.normal(size=(9,9,9,3))
        h,_=dirac(3,.8,9,field)
        spectrum=la.eigvalsh(h.toarray())
        assert min(abs(spectrum))>=8
        assert sum(spectrum<0)==sum(spectrum>0)==54


def test_actual_absolute_kernel_derivative_and_degenerate_spaces():
    result=derivative_controls()
    assert result['both_signed_spectral_subspace_current_identities']
    assert result['degenerate_subspace_operator_identity']
    assert result['omitted_Wilson_current_rejected']


def test_stationary_current_maxwell_and_gauge_law(stationary):
    assert stationary['eigen_residual']<1e-10
    assert stationary['fixed_point_residual']<2e-11
    assert stationary['maxwell_residual']<1e-9
    assert stationary['current_divergence_norm']<1e-10
    assert stationary['current_norm']>.005
    assert stationary['moment'][2]<-.04
    assert 1.4<stationary['g_actual_rest_energy']<1.6
    assert abs(stationary['g_actual_rest_energy']-stationary['g_supplied_mass'])>.04
    result=gauge_control(stationary)
    assert result['unchanged_state_in_changed_link_gauge_rejected']
    assert np.linalg.norm(maxwell(np.zeros_like(stationary['current']),1.,.1))==0


def test_independent_perturbation_and_zero_intrinsic_control(stationary):
    repeated=solve(initial_field=1.4*stationary['field'],relaxation=.7)
    assert np.linalg.norm(repeated['field']-stationary['field'])<2e-9
    naive=solve(wilson=0)
    assert naive['current_norm']<1e-10
    assert abs(naive['moment'][2])<1e-10


def test_actual_external_probe_matches_current_moment(stationary):
    result=external_probe(stationary)
    assert abs(result['difference'])<2e-6
    profile=dipole_profile(stationary)
    assert not profile['dipole_coefficient_fit_to_field']
    assert all(.85<r['ratio_to_infinite_space_dipole']<1.1 for r in profile['periodic_boundary_profile'])


def test_custody_and_invalid_construction(monkeypatch):
    from . import certificate as custody
    with pytest.raises(ValueError):solve(coupling=0)
    with pytest.raises(ValueError):solve(wilson=-1)
    with pytest.raises(ValueError):custody.validate_pins({'input_sha256':{'model/openwave_intrinsic_magnetic/solver.py':'0'*64}})
    monkeypatch.setattr(custody,'DEPENDENCIES',(*custody.DEPENDENCIES,'missing-stationary-magnetic.json'))
    with pytest.raises(FileNotFoundError):custody.inputs()
