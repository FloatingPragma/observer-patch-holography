"""Independent matrix checks and rejected false joins for the cone criterion."""
from fractions import Fraction as F
import json
import sys

import pytest
import sympy as sp

from . import cone as c
from . import certificate


def symbolic(rows):
    return sp.Matrix([[sp.Rational(x[0].numerator, x[0].denominator)
                       +sp.sqrt(5)*sp.Rational(x[1].numerator, x[1].denominator)
                       for x in row] for row in rows])


def test_exact_certificate_and_source_decoder():
    result = c.certificate()
    assert result['native_null_design_rank'] == 9
    assert result['constructive_match']['all_twelve_charges_zero']
    assert result['radial_control']['one_clock_cannot_match_both_radii']
    assert not result['scope']['source_selected_common_causal_field_clock']
    assert not result['scope']['source_radial_action_changed']


def test_inverse_pivots_and_principal_covector_tangent_duality():
    rows = c.m3([[0, 1, 0], [1, 0, 0], [0, 0, 2]])
    assert symbolic(c.inverse(rows)) == symbolic(rows).inv()
    a = c.m3([[3, 1, F(1, 3)], [1, 2, 0], [F(1, 3), 0, 4]])
    j = c.m3([[1, 2, 0], [0, 1, F(1, 2)], [0, 0, 2]])
    k = symbolic(c.field_form(a, j, F(7, 3)))
    expected = sp.diag(-sp.Rational(7, 3), *(sp.zeros(3).diagonal()))
    expected[1:, 1:] = symbolic(j).T*symbolic(a).inv()*symbolic(j)
    assert k == expected
    assert k.det() < 0


def test_nondiagonal_family_matches_at_every_supplied_positive_clock():
    j = c.m3([[2, 1, 0], [0, 1, 1], [0, 0, 3]])
    r2 = c.t.source().q5_sub(c.t.q(0), c.t.metric()[0][0])
    gram = c._multiply(j, c._transpose(j))
    for s in (1, F(3, 2), (F(5), F(1))):
        factor = c.t.source().q5_mul(r2, c.t.source().q5_inv(c.t.q(s)))
        a = tuple(tuple(c.t.source().q5_mul(factor, x) for x in row) for row in gram)
        result = c.match(a, j, s)
        assert set(result['null_charges']) == {c.t.q(0)}
        assert sp.simplify(symbolic(result['form'])-symbolic(c.t.metric())
                           *symbolic([[result['factor']]])[0, 0]) == sp.zeros(4)


def test_same_frame_covariance_keeps_cross_terms_and_distinguishes_wrong_observer_frame():
    r2 = c.t.source().q5_sub(c.t.q(0), c.t.metric()[0][0])
    a = [[1, 0, 0], [0, 4, 0], [0, 0, 9]]
    j = [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
    shear = c.t.matrix([[1, F(1, 3), 0, 0], [0, 1, 0, 0],
                        [0, 0, 1, 0], [0, 0, 0, 1]])
    changed = c.frame_covariance(a, j, r2, shear)
    assert changed['field_form'] == changed['observer_form']
    assert changed['field_form'][0][1] != c.t.q(0)
    with pytest.raises(ValueError):
        c.t.conformal_form(changed['field_form'])
    b = c.t.boost(F(1, 5))
    assert c.frame_covariance(a, j, r2, b)['field_form'] == c.t.metric()


def test_wrong_clock_and_spatial_dictionary_are_rejected():
    r2 = c.t.source().q5_sub(c.t.q(0), c.t.metric()[0][0])
    with pytest.raises(ValueError):
        c.match([[1, 0, 0], [0, 4, 0], [0, 0, 9]],
                [[1, 0, 0], [0, 1, 0], [0, 0, 1]], r2)
    with pytest.raises(ValueError):
        c.match([[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 1)


def test_radial_control_uses_fixed_map_and_cannot_be_fixed_by_one_clock():
    r2 = c.t.source().q5_sub(c.t.q(0), c.t.metric()[0][0])
    assert c.t.conformal_form(c.radial_product_form(1, 1, r2)) == c.t.q(1)
    for s in (r2, c.t.source().q5_mul(r2, c.t.q(F(1, 4)))):
        with pytest.raises(ValueError):
            c.t.conformal_form(c.radial_product_form(2, 1, s))
    # Independent spherical differential: dx=n dr+r dn, n.dn=0.
    direction = (F(3, 5), F(4, 5), 0)
    tangent = (F(-4, 5), F(3, 5), 0)
    matrix = symbolic(c.radial_product_form(2, 1, 1, direction))[1:, 1:]
    n, v = sp.Matrix(direction), sp.Matrix(tangent)
    assert (n.T*matrix*n)[0] == 1
    assert (v.T*matrix*v)[0] == F(1, 4)
    assert (n.T*matrix*v)[0] == 0


@pytest.mark.parametrize('tensor', [
    [[1, 0, 0], [0, 1, 0], [0, 0, 0]],
    [[1, 0, 0], [0, 1, 0], [0, 0, -1]],
    [[1, 1, 0], [0, 1, 0], [0, 0, 1]],
    [[1.0, 0, 0], [0, 1, 0], [0, 0, 1]],
])
def test_nonpositive_nonsymmetric_or_inexact_tensors_rejected(tensor):
    with pytest.raises(ValueError):
        c.field_form(tensor, [[1, 0, 0], [0, 1, 0], [0, 0, 1]], 1)


def test_singular_or_nonpositive_maps_rejected():
    identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    with pytest.raises(ValueError):
        c.field_form(identity, [[1, 0, 0], [0, 1, 0], [0, 0, 0]], 1)
    for s in (0, -1, 1.0):
        with pytest.raises(ValueError):
            c.field_form(identity, identity, s)
    for radius in (0, -1):
        with pytest.raises(ValueError):
            c.radial_product_form(radius, 1, 1)
    with pytest.raises(ValueError):
        c.radial_product_form(1, 1, 1, (1, 1, 0))
    with pytest.raises(ValueError):
        c.inverse([[1, 0, 0], [0, 0, 0], [0, 0, 1]])


def test_receipt_replays_source_and_rejects_an_altered_claim(tmp_path, monkeypatch):
    packet = certificate.build()
    packet['scope']['source_selected_common_causal_field_clock'] = True
    path = tmp_path/'wrong.json'
    path.write_text(json.dumps(packet))
    monkeypatch.setattr(sys, 'argv', ['certificate', '--check', str(path)])
    with pytest.raises(SystemExit, match='certificate differs'):
        certificate.main()
