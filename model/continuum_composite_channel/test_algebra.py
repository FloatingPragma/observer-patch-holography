import sympy as s
from . import algebra


def test_exact_complex_spinor_algebra():
    result=algebra.certificate()
    assert result['composite_CAR_exact']
    assert result['mixed_normal_and_anomalous_formulas_exact']
    assert result['microscopic_mode_added'] is False


def test_full_fock_covariance_against_trace_contraction():
    result=algebra.finite_covariance_control()
    assert result['exact_variance']==result['one_particle_trace_variance']=='1/16'
    assert result['Hilbert_Schmidt_upper_bound']=='1/2'


def test_unscaled_contact_is_not_canonical():
    c=[algebra.annihilator(j) for j in range(2)]
    b=(c[1].T*c[1])*c[0]
    assert algebra.anti(b,b.T)!=s.eye(4)
    assert algebra.anti(b,b.T)==c[1].T*c[1]


def test_local_mixed_defect_is_not_zero_operator():
    c=[algebra.annihilator(j) for j in range(2)]
    d=(2*c[1].T*c[1]-s.eye(4))*c[0]
    mixed=algebra.anti(d,c[0].T)
    assert mixed!=s.zeros(4)
    assert mixed*mixed==s.eye(4)


def test_exact_averaged_variance_in_full_tracial_fock_space():
    # This finite reference is an independent scaling control, not the OPH ground.
    for sites in (1,2,3):
        modes=2*sites;dim=1<<modes
        c=[algebra.annihilator(2*x+1,modes) for x in range(sites)]
        mixed=2*sum((a.T*a for a in c),s.zeros(dim))/sites-s.eye(dim)
        assert s.trace(mixed)==0
        assert s.trace(mixed*mixed)/dim==s.Rational(1,sites)
        assert mixed!=s.zeros(dim)
