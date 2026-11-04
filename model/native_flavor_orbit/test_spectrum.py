"""Independent high-precision source-spectrum and interval countercontrols."""
from fractions import Fraction as F

import mpmath as mp
import pytest
import sympy as s

from ..common_source import native_gauge_family_records as source
from . import algebra, intervals as iv, spectrum


def inside(value, bounds):
    low, high = map(F, bounds)
    return mp.mpf(low.numerator)/low.denominator <= value <= mp.mpf(high.numerator)/high.denominator


def test_direct_parent_duties_eigenvectors_and_all_45_rows_match_certified_intervals():
    with mp.workdps(70):
        parent = source.duty_coefficients()
        catalogue = spectrum.catalogue()
        for row in catalogue['policies']:
            matrices, eigs, vectors = {}, {}, {}
            for native in parent[row['policy']]:
                name = native['action_channel']
                matrices[name] = mp.matrix([[mp.mpf(str(x.evalf(75))) for x in native['coefficient'].row(i)]
                                            for i in range(3)])
                eigs[name], vectors[name] = mp.eigsy(matrices[name])
            for channel in row['channels']:
                name = channel['channel']
                for i, bounds in enumerate(channel['dimensionless_singular_values']):
                    if i == 0 and channel['rank'] == 2:
                        assert abs(eigs[name][i]) < mp.mpf('1e-65')
                    else:
                        assert inside(eigs[name][i], bounds)
            mixing = vectors['up'].T*vectors['down']
            for i in range(3):
                for j in range(3):
                    probability = min(mp.mpf(1), max(mp.mpf(0), mixing[i, j]**2))
                    assert inside(probability, row['quark_mixing_modulus_squared'][i][j])
            assert len(row['full_K45_squared_spectrum']) == 45
            assert row['full_K45_nullity'] == (3 if row['policy'] in ('always', 'C') else 17)
            assert row['charged_current_CP_commutator_trace_cube_exact'] == '0'


def test_exact_default_current_gates_and_nontrivial_real_mixing():
    result = spectrum.catalogue()
    gates = result['source_gates']
    assert gates['up']['p'] == gates['down']['p']
    assert F(gates['up']['q_interval'][0]) > F(gates['down']['q_interval'][1])
    mixing = result['policies'][0]['quark_mixing_modulus_squared']
    assert any(F(mixing[i][j][0]) > 0 for i in range(3) for j in range(3) if i != j)
    assert algebra.memory_certificate()['negative_factor_uniform_upper_bound'] == '-95/144'


def test_Machin_sine_and_square_root_enclosures_have_independent_precision_checks():
    with mp.workdps(80):
        assert inside(mp.pi, iv.pi_interval().strings())
        for radicand in (F(2), F(3), F(5), F(25, 8), F(9, 8)):
            assert inside(mp.sqrt(mp.mpf(radicand.numerator)/radicand.denominator),
                          iv.sqrt_rational(radicand).strings())
        for denominator in (3, 5, 7):
            x = iv.pi_interval()/denominator
            assert inside(mp.sin(mp.pi/denominator), iv.sine(x).strings())


def test_interval_LDL_checks_inertia_and_rejects_uncertain_pivots():
    matrix = [[iv.I(x) for x in row] for row in ((2, 1, 0), (1, 2, 0), (0, 0, 4))]
    assert [iv.inertia_below(matrix, F(x, 2)) for x in (1, 3, 7, 9)] == [0, 1, 2, 3]
    with pytest.raises(ValueError, match='pivot'):
        iv.inertia_below(matrix, 2)
    with pytest.raises(ValueError):
        iv.I(-1, 1).reciprocal()
    with pytest.raises(ValueError):
        iv.sine(iv.I(3))
    with pytest.raises(ValueError):
        iv.expression(s.Symbol('unbound'))


def test_spectral_projector_overlap_fails_if_probabilities_are_invented():
    first = [[iv.I(x) for x in row] for row in ((1, 0, 0), (0, 2, 0), (0, 0, 3))]
    values = [iv.I(x) for x in (1, 2, 3)]
    overlaps = spectrum.mixing_boxes(first, first, values, values)
    assert [[x.strings() for x in row] for row in overlaps] == [
        [['1', '1'], ['0', '0'], ['0', '0']],
        [['0', '0'], ['1', '1'], ['0', '0']],
        [['0', '0'], ['0', '0'], ['1', '1']]]
