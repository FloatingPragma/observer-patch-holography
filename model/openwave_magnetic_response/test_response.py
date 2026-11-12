"""Independent direct spectra, source changes and operator-bound controls."""
import numpy as np
import pytest
from scipy import sparse
from scipy.linalg import sqrtm
from .response import (covariant_momenta,square_blocks,positive_energy_minus_mass,
                       wavepacket,native_algebra,current_and_gauge_controls)


def test_complete_native_internal_multiplicity():
    result=native_algebra()
    assert result['dimension']==16 and result['native_valley_labels']==8
    assert sum(result['joint_sector_ranks'].values())==16


def test_direct_positive_matrix_function_and_analytic_remainder():
    px,py,_,_=covariant_momenta(5,.5,.05)
    _,_,plus,minus=square_blocks(px,py)
    for k in (plus,minus):
        approximation,bound=positive_energy_minus_mass(k,8.,.5)
        direct=sqrtm(64*np.eye(25)+k.toarray())-8*np.eye(25)
        assert np.linalg.norm(direct-approximation.toarray(),2)<bound+1e-12
        assert np.linalg.eigvalsh(k.toarray())[0]>-1e-12


def test_current_and_gauge_response():
    result=current_and_gauge_controls()
    assert result['link_current_is_derivative_of_actual_Peierls_H']
    assert result['untransformed_preparation_in_changed_gauge_rejected']


def test_magnetic_sign_reversal_and_zero_field():
    zero=square_blocks(*covariant_momenta(9,1.,0)[:2])[1]
    assert zero.nnz==0
    state=wavepacket(9,1.,1.)
    responses=[]
    for field in (.03,-.03):
        curvature=square_blocks(*covariant_momenta(9,1.,field)[:2])[1]
        responses.append(np.vdot(state,curvature@state).real)
    assert responses[0]>0 and abs(responses[0]+responses[1])<1e-13


def test_probe_has_vanishing_boundary():
    state=wavepacket(33,1.).reshape(33,33)
    assert np.max(np.abs(state[[0,-1],:]))<1e-16
    assert np.max(np.abs(state[:,[0,-1]]))<1e-16
    assert abs(np.linalg.norm(state)-1)<1e-14


def test_inadmissible_series_and_grid_are_rejected():
    with pytest.raises(ValueError):covariant_momenta(4,1.,.1)
    with pytest.raises(ValueError):positive_energy_minus_mass(sparse.eye(5),1.,.1)
