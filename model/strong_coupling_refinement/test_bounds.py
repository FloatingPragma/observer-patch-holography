from fractions import Fraction as F
import pytest
from . import bounds


def test_source_cube_at_old_parameters_is_outside_global_expansion_domain():
    with pytest.raises(ValueError):bounds.finite_graph(F(1,8),1,12)


def test_exact_controlled_fixed_graph_limits():
    row=bounds.finite_graph(65536,65536**2,12,1)
    assert row['spin_exchange']==1
    assert row['complete_dynamics_error_upper']==F(327,4096)<F(1,10)
    assert bounds.finite_graph(128,128**2,1,2)['complete_dynamics_error_upper']==F(5,64)


def test_physical_box_refinement_does_not_inherit_fixed_graph_control():
    fixedJ=bounds.refinement_exponents(2)
    assert fixedJ['nonzero_fixed_exchange'] and not fixedJ['global_error_proved_to_vanish']
    controlled=bounds.refinement_exponents(7)
    assert controlled['global_error_proved_to_vanish'] and controlled['exchange_power']==-5


@pytest.mark.parametrize('k',[1,2,64,256,1024,10**9])
def test_constructive_growing_u_schedule_has_exact_vanishing_error(k):
    row=bounds.growing_schedule(k)
    assert row['U']==k*k and row['log2_c']==256*k*k
    assert row['U_squared_mass_upper'] < row['U_sqrt_mass_upper']**2
    assert row['U_over_R']==F(1,k)
    assert max(row['exponential_remainder_powers']) <= -1


def test_gaussian_envelope_matches_direct_exact_evaluator_at_small_k():
    # k=1,2 are small enough to evaluate the full rational Gaussian majorant.
    # Larger k use the proved exponent inequalities, never float underflow.
    from ..continuum_free_limit.bounds import gaussian_low_energy_upper
    for k in (1,2):
        row=bounds.growing_schedule(k)
        mass=gaussian_low_energy_upper(2**row['log2_c'],row['U'],row['R'],F(1,k))
        assert row['U']**2*mass <= row['U_squared_mass_upper']


@pytest.mark.parametrize('bad',[True,0,-1,1.0,'2'])
def test_schedule_rejects_wrong_domains(bad):
    with pytest.raises(ValueError):bounds.growing_schedule(bad)


def test_finite_bound_requires_verified_separation_and_exact_inputs():
    for args in [(1,3,1),(1.,8,1),(True,8,1),(1,8,0)]:
        with pytest.raises(ValueError):bounds.finite_graph(*args)
