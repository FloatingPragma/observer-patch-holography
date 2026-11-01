from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from ..field_symmetry import symmetry as s
from . import connected as c


def test_exact_contact_derivative_and_full_fermionic_Wick_cancellation():
    value=c.local_derivative_certificate()
    assert value['spins_checked']==2
    assert value['anomalous_pairings_retained'] is True
    assert c.wick_pairing_residual()=={}
    # A bosonic all-plus replacement must fail even before choosing a state.
    assert c.wick_pairing_residual((1,1,1))
    assert c.wick_pairing_residual((1,-1,-1))


def test_interacting_atomic_control_has_the_predicted_positive_variance():
    # Independent exact-time diagnostic. The product density is quasifree at
    # t=0, but its contact-evolved field is not Gaussian at unequal times.
    a=[s._annihilator(2,j).matrix().toarray() for j in range(2)]
    n=[x.conj().T@x for x in a];h=n[0]@n[1]
    r=F(1,4)
    density=np.diag([float((r if mask&1 else 1-r)*(r if mask&2 else 1-r))
                     for mask in range(4)])
    for t in (0.001,-0.001,0.02):
        u=expm(1j*t*h);at=u@a[0]@u.conj().T
        anti=at@a[0].conj().T+a[0].conj().T@at
        centered=anti-np.trace(density@anti)*np.eye(4)
        actual=np.trace(density@centered.conj().T@centered).real
        exact=float(r*(1-r))*abs(np.exp(-1j*t)-1)**2
        assert abs(actual-exact)<1e-15
        assert actual>float(c.variance_lower(F(str(t)),0,1))
        d={'a':at,'ad':at.conj().T,'b':a[0],'bd':a[0].conj().T}
        expectation=lambda *words:np.trace(density@np.linalg.multi_dot(words)) if len(words)>1 else np.trace(density@words[0])
        cumulants=[]
        for word in c.WORDS:
            w,x,y,z=(d[key] for key in word)
            cumulants.append(expectation(w,x,y,z)-expectation(w,x)*expectation(y,z)
                             +expectation(w,y)*expectation(x,z)-expectation(w,z)*expectation(x,y))
        assert abs(sum(cumulants)-actual)<1e-15


def test_uniform_second_derivative_bound_and_nonzero_interval_are_rational():
    value=c.constants()
    assert value['second_derivative_norm']==F(169,64)
    assert value['time_radius']==F(4,169)
    assert c.variance_lower(F(4,169))==F(1,114244)
    for speed in (F(0),F(1,8),F(2)):
        for interaction in (F(1,4),F(1),F(3)):
            b=c.constants(speed,interaction);t=b['time_radius']
            # Replace sqrt(7)/8 by the smaller rational 1/4.
            amplitude=t*interaction/4-2*b['second_derivative_norm']*t*t
            assert amplitude==t*interaction/8
            assert amplitude*amplitude==c.variance_lower(t,speed,interaction)
            assert c.variance_lower(-t,speed,interaction)==amplitude*amplitude


def test_paired_auxiliary_quadratic_dynamics_is_a_distinct_zero_control():
    value=c.quadratic_negative_control()
    assert value['paired_and_auxiliary_modes_present'] is True
    assert value['scalar_anticommutators'] is True
    assert value['state_needed_for_scalar_operator_identity'] is False


@pytest.mark.parametrize('call',[
    lambda:c.constants(0.125,1),lambda:c.constants(True,1),
    lambda:c.constants(-1,1),lambda:c.constants(1,0),
    lambda:c.variance_lower(0),lambda:c.variance_lower(F(1)),
    lambda:c.variance_lower(0.001),
])
def test_nonexact_or_outside_proved_domain_rejected(call):
    with pytest.raises(ValueError):call()


def test_connected_field_claim_keeps_continuum_and_encoding_boundaries():
    value=c.certificate()
    assert value['non_Gaussian_multitime_ground_field_proved'] is True
    assert value['linear_CAR_quadratic_evolution_excluded'] is True
    assert value['quadratic_exclusion_requires_Gaussian_reference'] is False
    assert value['nonlinear_observable_encodings_excluded'] is False
    assert value['local_state_faithfulness_assumed'] is False
    assert value['continuum_limit_or_scattering_proved'] is False
    assert value['physical_particle_species_identified'] is False
