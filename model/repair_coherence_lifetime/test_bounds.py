from fractions import Fraction as F
from itertools import combinations

import pytest

from model.exchange_sector import many_body as native
from . import bounds as b
from . import certificate as c


def retained_coherence(first, second, law):
    """Independent complete-record transport of one native matrix unit."""
    mass = F()
    for seam in native.SEAMS:
        for coin in (0, 1):
            r, x, s, _ = native.column(first, law, seam, coin, -1)
            t, y, u, _ = native.column(second, law, seam, coin, -1)
            if r == t:
                assert x != y
                assert abs(s*u) == 1
                mass += F(1, 2*len(native.SEAMS))
    return mass


@pytest.mark.parametrize('law', ['C', 'E'])
@pytest.mark.parametrize('first,second', [(0, 4095), (0, 1), (1, 2), (31, 47), (1365, 2730)])
def test_formula_matches_actual_full_record_mass(law, first, second):
    assert retained_coherence(first, second, law) == b.survival(first, second, 12, native.SEAMS, law=law)


def test_extremal_cat_is_not_falsely_damped():
    assert b.survival(0, 4095, 12, native.SEAMS) == 1
    assert b.survival(1, 2, 12, native.SEAMS) == F(5, 6)
    assert b.survival(0, 1, 12, native.SEAMS) == F(11, 12)


@pytest.mark.parametrize('vertices,edges', [(3, ((0, 1), (1, 2))), (4, ((0, 1), (0, 3), (1, 2), (2, 3)))])
def test_coarea_all_matrix_units_on_other_graphs(vertices, edges):
    weights = tuple(F(i+1, sum(range(1, len(edges)+1))) for i in range(len(edges)))
    cuts = b.cut_table(vertices, edges, weights)
    floor = min(cuts[1:-1])
    full = (1 << vertices)-1
    for first, second in combinations(range(1 << vertices), 2):
        plus, minus = first & ~second, second & ~first
        loss = (cuts[plus]+cuts[minus])/2
        assert b.survival(first, second, vertices, edges, weights) == 1-loss
        if (first, second) != (0, full):
            assert loss >= floor/2
        if first.bit_count() == second.bit_count():
            assert loss >= floor


def test_exact_native_spectral_and_cut_witnesses():
    assert c.diagonal_gap()['strict_rational_gap_floor'] == '11/240'
    assert c.native_cuts()['ternary_differences_checked'] == 3**12
    assert c.local_overlap_table()['C'][1][2] == '0'


def test_complete_source_inventory():
    assert c.source_columns()['complete_native_columns_checked'] == 2*30*2*4096*2


def test_actual_single_particle_coherence_is_an_eigenoperator():
    first, second = 1, 2
    outputs = {}
    for seam in native.SEAMS:
        for coin in (0, 1):
            r, x, s, _ = native.column(first, 'C', seam, coin, -1)
            t, y, u, _ = native.column(second, 'C', seam, coin, -1)
            if r == t:
                outputs[x, y] = outputs.get((x, y), F())+F(s*u, 60)
    assert outputs == {(first, second): F(5, 6)}


def test_positive_rate_does_not_assume_physical_clock_or_repair():
    cuts = b.cut_table(3, ((0, 1), (1, 2)))
    assert b.coherence_upper(1, cuts, 100, repair_probability=0) == 1
    assert b.coherence_upper(1, cuts, 2, repair_probability=F(1, 2)) == F(9, 16)
    assert b.diffusion_coherence_exponent(1, F(1, 10), 1) == 1200
    assert b.diffusion_coherence_exponent(1, F(1, 20), 1) == 4800


@pytest.mark.parametrize('vertices,edges,weights', [
    (True, ((0, 1),), None), (3, ((0, 1),), None),
    (2, ((0, 1), (0, 1)), None), (3, ((0, 1), (1, 2)), (F(1), F(0))),
    (2, ((0, 1),), (1.0,)), (2, ((1, 0),), None)])
def test_bad_or_disconnected_selector_refused(vertices, edges, weights):
    with pytest.raises(ValueError):
        b.graph(vertices, edges, weights)


def test_missing_edge_or_wrong_record_formula_is_detected(monkeypatch):
    monkeypatch.setattr(native, 'SEAMS', native.SEAMS[:-1])
    with pytest.raises(ArithmeticError):
        c.diagonal_gap()


def test_bound_targets_trace_norm_not_a_particle_mass():
    assert b.mixing_upper(12, 6, 0, F(11, 240), F(1, 6)) == 2
    assert b.mixing_upper(12, 6, 250, F(11, 240), F(1, 6)) < F(1, 100)


@pytest.mark.parametrize('cuts', [(0, 2, 2, 0), (0, F(1, 2), 0), (False, 1, 1, 0)])
def test_malformed_cut_input_cannot_produce_negative_norm_bound(cuts):
    with pytest.raises(ValueError):
        b.coherence_upper(1, cuts, 1)
