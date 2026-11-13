import numpy as np
from scipy.linalg import eigvalsh
from .cochains import relative_cube, transverse_frame


def test_exact_complete_curl_gradient_and_counts():
    for n in (2, 3, 4, 6):
        d = relative_cube(n)
        assert (d['curl']@d['gradient']).nnz == 0
        assert len(d['edges'])-len(d['vertices']) == d['transverse_dimension']
        assert d['harmonic_one_forms'] == 0


def test_relative_Maxwell_stiffness_has_checked_positive_transverse_domain():
    for n in (3, 4):
        d = transverse_frame(n); t = d['frame']
        k = (d['curl']@t).T@(d['curl']@t)
        e = eigvalsh(k)
        assert e[0] > 8/n**2
        assert abs(e[0]-8*np.sin(np.pi/(2*n))**2) < 1e-12
        assert len(e) == d['transverse_dimension']


def test_all_boundary_Gauss_constraints_and_no_discarded_harmonics():
    from .cochains import complete_transverse_frame
    for n in (2, 3):
        d = complete_transverse_frame(n); t = d['frame']
        assert np.max(abs(d['gradient'].T@t)) < 1e-13
        assert d['gradient'].shape[1] == (n+1)**3
        e = eigvalsh((d['curl']@t).T@(d['curl']@t))
        assert abs(e[0]-8*np.sin(np.pi/(2*(n+1)))**2) < 1e-12
