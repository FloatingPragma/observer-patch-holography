from fractions import Fraction as Q
from unittest.mock import patch

import pytest

from . import certificate as c


def varied():
    return tuple(Q(i+1, 2211) for i in range(66))


def test_every_actual_native_seam_is_its_uniform_fibre_heat_bath():
    certificate = c.heat_bath_certificate()
    assert certificate['source_coin_branches'] == 3960
    assert certificate['connected_token_vertices'] == 66
    assert all(row['singleton_fibres'] == 46 and row['two_point_fibres'] == 10
               for row in certificate['local_heat_baths'])
    assert not certificate['random_scan_idempotent']
    assert certificate['nonidempotence_witness']['P'] != certificate['nonidempotence_witness']['P_squared']


def test_actual_source_mutation_cannot_be_replaced_by_desired_heat_bath():
    original = c.cut.seam_population

    def changed(edge):
        rows = [list(row) for row in original(edge)]
        i = next(i for i, row in enumerate(rows) if row[i] == Q(1, 2))
        j = next(j for j, w in enumerate(rows[i]) if j != i and w)
        rows[i][i] += Q(1, 10)
        rows[i][j] -= Q(1, 10)
        return tuple(map(tuple, rows))

    with patch.object(c.cut, 'seam_population', side_effect=changed):
        with pytest.raises(ValueError, match='not its uniform conditional resampling'):
            c.heat_bath_certificate()


def test_exact_log_arithmetic_and_zero_support_conventions():
    assert c.combine((1, c.log_rational(Q(12, 35))), (1, c.log_rational(Q(35, 12)))) == ()
    assert c.entropy((Q(1), Q())) == ()
    assert c.relative_entropy((Q(1), Q()), (Q(1, 2), Q(1, 2))) == ((2, Q(1)),)
    with pytest.raises(ValueError, match='infinite'):
        c.relative_entropy((Q(1), Q()), (Q(), Q(1)))
    c.log_rational(Q(1))
    with pytest.raises(ValueError, match='no floats or Booleans'):
        c.log_rational(True)


def test_state_and_conditional_row_optimization_share_the_same_reference():
    p = varied()
    u = tuple(Q(1, 66) for _ in range(66))
    assert c.entropy(p) == c.combine((1, c.log_rational(66)), (-1, c.relative_entropy(p, u)))
    for row in ((Q(1, 7), Q(6, 7)), (Q(), Q(1)), (Q(1, 2), Q(1, 2))):
        ref = (Q(1, 2), Q(1, 2))
        assert c.entropy(row) == c.combine((1, c.log_rational(2)), (-1, c.relative_entropy(row, ref)))


@pytest.mark.parametrize('zero_support', [False, True])
def test_local_entropy_identity_and_fixed_schedule_telescoping(zero_support):
    p = tuple(Q(i == 0) for i in range(66)) if zero_support else varied()
    initial, terms = p, []
    for edge in c.born.carrier_seams()[::7]:
        terms.append(c.projection_entropy(p, edge))
        p = c.push(p, c.kernel(edge))
    assert c.combine(*[(1, x) for x in terms]) == c.combine((1, c.entropy(p)), (-1, c.entropy(initial)))


def test_random_scan_requires_the_nonzero_mixing_information_term():
    p = varied()
    weights = tuple(Q(i+1, 465) for i in range(30))
    result = c.random_scan_entropy(p, weights)
    assert result['mixing_information']
    assert result['entropy_change'] != result['mean_local_dissipation']
    assert result['entropy_change'] == c.combine((1, result['mean_local_dissipation']),
                                               (1, result['mixing_information']))
    # An exogenous one-seam scan is the projection case; the mixing term vanishes.
    only = c.random_scan_entropy(p, (Q(1),)+tuple(Q() for _ in range(29)))
    assert only['mixing_information'] == ()


def test_integral_fluctuation_equality_and_necessary_support_boundary():
    P = c.kernel(c.born.carrier_seams()[0])
    assert c.fluctuation_average(varied(), P) == 1
    assert c.fluctuation_average(varied(), c.random_scan()) == 1
    moving = next(i for i in range(66) if P[i][i] == Q(1, 2))
    point = tuple(Q(i == moving) for i in range(66))
    with pytest.raises(ValueError, match='positive'):
        c.fluctuation_average(point, P)
    assert c.fluctuation_average(point, P, require_positive=False) == Q(1, 2)


def test_complete_tape_preserves_input_information_and_adds_the_coin_entropy():
    tree = c.born.source_tree()
    p = varied()
    full = tuple(p[c.INDEX[x]] if x in c.INDEX else Q() for x in c.born.BASIS)
    joint = tuple(c.born.execute(tree, full).values())
    assert c.combine((1, c.entropy(joint)), (-1, c.entropy(p))) == ((2, Q(2)),)
    # Reversibility requires terminal configuration AND tape; tape alone conflates inputs.
    assert len(tree.records) < len(c.born.BASIS)


def test_full_78_state_startup_has_a_strict_entropy_drop():
    result = c.startup_control()
    assert not result['uniform_78_is_stationary']
    assert result['one_seam_entropy_change'] == '-log(3)/26'
    assert result['one_covering_scan_collision_probability'] == '0'


@pytest.mark.parametrize('bad', [True, 1.0, Q(-1)])
def test_exact_probability_inputs_reject_coercions_and_negative_mass(bad):
    with pytest.raises(ValueError):
        c.law((bad,)+tuple(Q() for _ in range(65)))


def test_fluctuation_requires_the_same_preserved_reference():
    collapsed = tuple(tuple(Q(j == 0) for j in range(66)) for _ in range(66))
    with pytest.raises(ValueError, match='column sums one'):
        c.fluctuation_average(varied(), collapsed)


def test_zero_mass_boolean_kernel_entry_is_not_silently_coerced():
    with pytest.raises(ValueError, match='no floats or Booleans'):
        c.sparse_rows(((Q(1), False), (Q(), Q(1))))
