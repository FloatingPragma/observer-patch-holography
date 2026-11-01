"""Independent source, signed-trace, all-state and rational-tail controls."""
from fractions import Fraction as F
import math

import mpmath as mp
import numpy as np
import pytest
import sympy as s

from model.regional_modular_expansion import expansion as e
from model.regional_modular_expansion import packet
from model.field_thermal.thermal import car_restriction, restrict
from model.full_fock_field.field import FullFock


def test_literal_full_field_and_independent_signed_partial_trace():
    data = e.plaquette()
    h = data['H_numerator'].matrix(data['scale']).toarray()
    actual = FullFock(8).hamiltonian(e.GEOMETRY).toarray()-np.eye(256)
    np.testing.assert_array_equal(h, actual)
    for kept in ((0, 1, 2, 3), (1, 4, 6, 7)):
        exact = np.asarray(e.exact_trace(data['H_numerator'], kept=kept), complex)/data['scale']
        np.testing.assert_array_equal(exact, restrict(actual, car_restriction(8, kept))/16)
    # Off-diagonal odd words check the sign; H alone has even parity.
    rows, cols, signs = [], [], []
    for mask in range(256):
        if mask >> 4 & 1:
            rows.append(mask ^ 16); cols.append(mask)
            signs.append((-1)**(mask & 15).bit_count())
    from scipy import sparse
    operator = e.GI.from_exact(sparse.coo_matrix((signs, (rows, cols)), shape=(256, 256)))
    expected = restrict(operator.matrix().toarray(), car_restriction(8, (1, 4, 6, 7)))/16
    np.testing.assert_array_equal(np.asarray(e.exact_trace(operator, kept=(1, 4, 6, 7)), complex), expected)


def test_exact_noncommuting_cubic_formula_and_interaction_deletion():
    data = e.plaquette()
    A, B, C = data['moments'][:3]
    c1, c2, c3 = data['coefficients'][:3]
    assert c1 == A
    assert c2 == (A*A-B)/2
    assert c3 == C/6-(A*B+B*A)/4+A*A*A/3
    assert c2 == -s.Rational(99, 512)*s.eye(16)
    assert c1*c3-c3*c1 != s.zeros(16)
    terms = e.pauli_terms(e.centered(c3))
    assert terms['IIZZ'] == terms['ZZII'] == '-1/3072'
    assert terms['IIIZ'] == terms['ZIII'] == '1/4096'
    assert terms['IXXI'] == '1/98304'
    free = e.pauli_terms(e.centered(e.plaquette(0)['coefficients'][2]))
    assert set(free) == {'IXXI', 'IYYI', 'XZZX', 'YZZY'}
    assert all(e.centered(c) == s.zeros(16)
               for c in e.plaquette(boundary=False)['coefficients'][1:])


@pytest.mark.parametrize('beta', [F(1, 64), F(1, 32), F(1, 16), F(1, 8)])
def test_independent_complete_gibbs_log_satisfies_operator_bound(beta):
    data = e.plaquette()
    result = packet.numerical_comparison(data, beta)
    assert result['full_operator_polynomial_error'] <= float(e.logarithm_tail(beta, data['bound']))+1e-13
    if beta == e.BETA:
        proof = e.witness(data, beta)
        assert float(proof['exact_difference_lower']) <= result['full_minus_empty_bare_difference'] <= float(proof['exact_difference_upper'])
        assert proof['bare_generator_distance_lower'] > F(27, 10**9)


@pytest.mark.parametrize('x,order', [(F(1, 100), 3), (F(3, 32), 8), (F(1, 2), 8)])
def test_outward_rational_tail_against_independent_high_precision(x, order):
    with mp.workdps(100):
        point = mp.mpf(x.numerator)/x.denominator
        reference = -mp.log(2-mp.exp(point))
        for n, a in enumerate(e.majorant_coefficients(order), 1):
            reference -= mp.mpf(a.numerator)/a.denominator*point**n
        bound = e.logarithm_tail(x, 1, order)
        assert 0 < reference <= mp.mpf(bound.numerator)/bound.denominator
        assert mp.mpf(bound.numerator)/bound.denominator-reference < mp.mpf('1e-19')


def test_noncommuting_entropy_tangent_and_untouched_reference():
    data = e.plaquette()
    h = data['H_numerator'].matrix(data['scale']).toarray()
    values, vectors = np.linalg.eigh(h)
    density = (vectors*np.exp(-float(e.BETA)*values))@vectors.conj().T
    regional = restrict(density, car_restriction(8, e.KEPT))
    regional /= np.trace(regional)
    weights, basis = np.linalg.eigh(regional)
    log_reference = (basis*np.log(weights))@basis.conj().T
    vector = np.zeros(16, complex); vector[0] = 1/math.sqrt(2); vector[3] = 1j/math.sqrt(2)
    direction = np.outer(vector, vector.conj())-regional
    assert np.linalg.norm(direction@regional-regional@direction) > 1e-4
    polynomial = np.asarray(e.polynomial(data), complex)
    error = abs(np.trace(direction@(-log_reference-polynomial)))
    norm = np.linalg.norm(direction, ord='nuc')
    assert error <= float(e.logarithm_tail(e.BETA, data['bound']))*norm+1e-14
    # Tensoring the operator error with identity preserves the bound, rather
    # than dividing it by an arbitrary spectator dimension.
    scalar = np.trace(-log_reference-polynomial)/16
    residue = -log_reference-polynomial-scalar*np.eye(16)
    assert np.linalg.norm(np.kron(residue, np.eye(3)), 2) == pytest.approx(np.linalg.norm(residue, 2), abs=1e-15)


@pytest.mark.parametrize('beta,bound', [(True, 1), (.1, 1), (F(-1, 32), 3), (F(1, 4), 3), (1, -1)])
def test_rational_and_high_temperature_domain_fails_closed(beta, bound):
    with pytest.raises(ValueError):
        e.logarithm_tail(beta, bound)


def test_partial_trace_rejects_sector_truncation_and_malformed_modes():
    from scipy import sparse
    with pytest.raises(ValueError, match='complete Fock'):
        e.exact_trace(e.GI.from_exact(sparse.eye(28)))
    for kept in ((True, 1), (2, 1), (0, 0), (0, 8), ()):
        with pytest.raises(ValueError):
            e.exact_trace(e.plaquette()['H_numerator'], kept=kept)


def test_retained_packet_is_exactly_source_bound():
    from tools import strict_json
    from model.observer.scalar_response.experiment import compare_receipt
    actual = packet.run()
    compare_receipt(strict_json.loads(packet.OUTPUT.read_text()), actual)
    assert sum(actual['domain']['number_sector_dimensions']) == 256
    assert actual['scope']['no_Einstein_or_physical_particle_gate_accepted']


def test_cli_rejects_malformed_supplied_receipt_path(tmp_path, monkeypatch):
    import sys
    broken = tmp_path/'broken.json'
    broken.write_text('{"schema": "wrong"}')
    monkeypatch.setattr(sys, 'argv', ['packet', '--check', str(broken)])
    with pytest.raises((ValueError, AssertionError, KeyError)):
        packet.main()
