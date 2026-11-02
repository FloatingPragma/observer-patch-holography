from fractions import Fraction as F
import pytest
from . import bounds as b


def test_joint_trajectory_closes_all_extensive_errors():
    for E in (2, 54, 300, 10000):
        r = b.trajectory(E)
        assert r['isometry_error'] == F(2, E**3)
        assert r['generator_error'] == F(3, E)
        assert r['complete_record_trace_norm_error'] == min(F(2), F(8*(1+E), E**3)+F(6, E))
        assert r['charged_gap_lower'] > 0
    assert b.cube(8)['complete_record_trace_norm_error'] < F(1, 100)


def test_full_trace_history_bound_and_poisson_are_distinct():
    assert b.history(512, 512**2, 2, F(1, 2), 3) == F(9, 64)
    assert b.trajectory(100, duration=0)['complete_record_trace_norm_error'] == F(8, 100**3)


def test_gap_bound_decays_without_claiming_a_massless_particle():
    for w in (3, 5, 9):
        r = b.torus_gap(w)
        assert r['source_gap_upper'] == F(5, w*w)+F(1, w**6)
        assert r['gap_counts_ground_multiplicity'] is True
        assert r['infinite_GNS_massless_particle_inferred'] is False
    assert b.torus_gap(101)['source_gap_upper'] < F(1, 1000)


def test_infinite_spin_locality_is_used_at_fixed_exchange():
    r = b.local_comparison(160, 0, 1, F(1), F(1), F(1))
    assert r['spin_locality_error'] < F(1, 1000)
    assert r['local_expectation_error'] < F(1, 100)


@pytest.mark.parametrize('kwargs', [{'edges': True}, {'edges': 1}, {'edges': 2, 'exchange': 0},
                                  {'edges': 2, 'duration': -1}, {'edges': 2, 'repair_rate_per_edge': 1.0}])
def test_invalid_scaling_inputs_rejected(kwargs):
    with pytest.raises(ValueError):
        b.trajectory(**kwargs)


@pytest.mark.parametrize('width', [True, 1, 2, 4, 5.0])
def test_twist_family_requires_odd_transverse_width(width):
    with pytest.raises(ValueError):
        b.torus_gap(width)
