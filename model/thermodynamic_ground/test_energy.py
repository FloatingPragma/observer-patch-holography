from fractions import Fraction as F
from itertools import product

import pytest

from model.thermodynamic_ground import energy as e


def test_crossing_bond_formula_against_direct_cubic_enumeration():
    for large,tile in ((2,2),(4,2),(6,2),(8,4)):
        crossing=0
        for x in product(range(large),repeat=3):
            for axis in range(3):
                if x[axis]+1<large and x[axis]//tile!=(x[axis]+1)//tile:crossing+=1
        assert e.intertile_bonds(large,tile)==crossing


def test_finite_energy_uncertainty_and_surface_width_are_both_retained():
    low,high=e.finite_energy_bracket(4,-2,-1,speed=F(1,8))
    assert low==-F(1,8) and high==-F(1,64)
    assert high-low==F(1,64)+F(3,32)
    assert e.finite_energy_bracket(2,-1,-1,speed=0)==(-F(1,8),-F(1,8))


def test_untiled_terms_are_counted_and_vanish_per_volume():
    assert e.remainder_norm_bound(8,4)==0
    assert e.remainder_norm_bound(6,4)==F(7,4)*(6**3-4**3)
    ratios=[e.remainder_norm_bound(m,4)/m**3 for m in (6,10,18,34)]
    assert all(x>y for x,y in zip(ratios,ratios[1:],strict=False))


def test_existing_cube_global_certificate_is_consumed_without_a_new_eigensolver():
    result=e.certificate()
    assert result['finite_example']['all_source_hashes_verified'] is True
    assert result['finite_example']['full_ground_energy']=='-sqrt3/4'
    assert result['finite_ground_solver_run_for_growing_boxes'] is False
    assert result['state_convergence_rate_or_particle_rest_mass'] is False


@pytest.mark.parametrize('call',[
    lambda:e.intertile_bonds(6,4),lambda:e.intertile_bonds(4,3),
    lambda:e.finite_energy_bracket(2,1,0),lambda:e.finite_energy_bracket(2,0,0,speed=.125),
    lambda:e.remainder_norm_bound(2,4),lambda:e.remainder_norm_bound(4,2,interaction=-1),
])
def test_uncertified_domains_rejected(call):
    with pytest.raises(ValueError):call()
