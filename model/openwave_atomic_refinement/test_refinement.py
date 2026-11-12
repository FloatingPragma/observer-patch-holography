import numpy as np
from .refinement import spectrum,volume_control
from model.openwave_field_observables.orbitals import reduced_hamiltonian


def test_nested_green_comparison():
    assert volume_control()['strictly_more_attractive_in_larger_box']


def test_parity_spectral_residuals_and_volume_order():
    small=spectrum(7,8.);large=spectrum(11,12.)
    for k in range(2):
        assert large['sectors'][k]['energies'][0]<small['sectors'][k]['energies'][0]
        assert max(large['sectors'][k]['residual_norms'])<1e-9


def test_positive_vector_bound_and_wrong_charge_control():
    h,_=reduced_hamiltonian(7,8.)
    v=np.linspace(.1,1.,h.shape[0]);ratios=(h@v)/v
    exact=np.linalg.eigvalsh(h.toarray())[0]
    assert min(ratios)<=exact<=max(ratios)
    wrong,_=reduced_hamiltonian(7,8.,potential_sign=-1.)
    assert np.linalg.eigvalsh(wrong.toarray())[0]>0
