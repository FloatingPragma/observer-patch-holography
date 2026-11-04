from fractions import Fraction as F

import pytest
import sympy as s

from . import algebra, bulk
from ..native_mass_channels import algebra as source


def test_uniform_bulk_dirac_blocks_include_all_massive_and_massless_branches():
    gamma, masses = source.folded_generators()
    m0, magnitude = s.Rational(3), s.Rational(2)
    h = sum((x*g for x, g in zip((1, 2, 3), gamma, strict=True)), s.zeros(16))+m0*masses[0]+magnitude*masses[3]
    all_columns = []
    for sign, _, basis in bulk.dirac_blocks():
        small = basis.H*h*basis
        expected_squared = 14+(magnitude+sign*m0)**2
        assert (small*small-expected_squared*s.eye(4)).applyfunc(s.expand) == s.zeros(4)
        assert s.trace(small) == 0
        all_columns.extend(basis[:, j] for j in range(4))
    assert s.Matrix.hstack(*all_columns).rank() == 16
    critical = 2*masses[0]+2*masses[3]
    assert critical.rank() == 8
    result = bulk.certificate()['bulk']
    assert result['at_m0_equal_abs_m'] == {'massless_Dirac_blocks': 2, 'mass_2abs_m_Dirac_blocks': 2}


@pytest.mark.parametrize('direction', [(2, 3, 6), (0, 0, -1), (1, 0, 0)])
def test_actual_source_rotation_transports_every_kinetic_and_mass_block(direction):
    gamma, masses = source.folded_generators()
    rotation = bulk.spin_rotation(direction)
    length = s.sqrt(sum(x*x for x in direction))
    target = sum((s.Rational(x)*g/length for x, g in zip(direction, masses[1:], strict=True)), s.zeros(16))
    assert (rotation*rotation.H-s.eye(16)).applyfunc(s.simplify) == s.zeros(16)
    assert (rotation*masses[3]*rotation.H-target).applyfunc(s.simplify) == s.zeros(16)
    assert all((rotation*g-g*rotation).applyfunc(s.simplify) == s.zeros(16) for g in gamma+(masses[0],))


def test_coulomb_exterior_form_signs_and_complete_orbital_degeneracies():
    result = bulk.angular_potential_certificate()
    assert [x['rank'] for x in result['form_sectors']] == [2, 6, 6, 2]
    assert result['two_and_three_form_subthreshold_states'] == 0
    radius = s.symbols('r', positive=True)
    for principal in range(1, 5):
        assert sum(2*ell+1 for ell in range(principal)) == principal**2
        for ell in range(principal):
            radial = s.exp(-radius/principal)*radius**ell*s.assoc_laguerre(principal-ell-1, 2*ell+1, 2*radius/principal)
            actual = (-s.diff(radial, radius, 2)-2*s.diff(radial, radius)/radius+
                      ell*(ell+1)*radial/radius**2+(1-2/radius)*radial)
            expected = (1-s.Rational(1, principal**2))*radial
            assert s.simplify(actual-expected) == 0
        level = bulk.coulomb_level(principal)
        assert level['each_signed_energy_multiplicity'] == (0 if principal == 1 else 2*principal**2)


def test_bulk_domain_and_smooth_core_boundary():
    with pytest.raises(ValueError):
        bulk.spin_rotation((0, 0, 0))
    for value in (True, 0, -1, F(2), 2.):
        with pytest.raises(ValueError):
            bulk.coulomb_level(value)
    assert bulk.certificate()['generic_smooth_core_all_excited_levels_solved'] is False
