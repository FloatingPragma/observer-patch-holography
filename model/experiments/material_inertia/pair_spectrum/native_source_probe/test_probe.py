"""Independent source algebra, native routing, calibration and complete records."""

from fractions import Fraction
import json
from itertools import product
import subprocess
import sys

import numpy as np
import sympy as sp

from model.experiments.material_inertia.pair_spectrum.native_source_probe import circuit, physical, measure as m


def test_source_even_odd_cost_and_full_state_stroboscopic_cancellation():
    t = 8 * sp.pi / 31
    for n in (0, 1, 2):
        q = n - 1
        def cost(source):
            return 8 * (n - source)**2 + (sp.Rational(1, 2) if n != source else 0)
        assert cost(0) - sp.Rational(17, 2) == sp.Rational(31, 4) * q*q + sp.Rational(65, 4) * q
        assert cost(2) - sp.Rational(17, 2) == sp.Rational(31, 4) * q*q - sp.Rational(65, 4) * q
        assert sp.simplify(sp.exp(-sp.I * t * sp.Rational(31, 4) * q*q)) == 1
    meters, audits = circuit.audit()
    assert len(audits) == 2
    assert max(x['full_operator_error'] for x in audits) < 3e-12
    assert max(audits[0]['extreme_refocus_omission_leakage']) > 1e-4
    assert abs(meters[1][-1]['X'][2] - meters[1][1]['X'][2]) < 1e-12
    assert abs(meters[1][-1]['Y'][2] - meters[1][1]['Y'][2]) > 0.8


def test_actual_y_reference_and_inverse_on_every_ternary_word():
    from model.experiments.interval_echo_consensus.signed_effects import fixed_pointer as p
    u = circuit.apply([p.pulse('reference', 2, 2, Fraction(1, 64))], np.eye(27))
    expected = np.diag(np.exp(-1j * np.pi * circuit.STATES[:, 2]**2 / 8))
    np.testing.assert_allclose(u, expected, atol=2e-13)
    H = circuit.apply(p.neutral_basis(1, 2), np.eye(27))
    inverse = circuit.apply(p.inverse(p.neutral_basis(1, 2)), np.eye(27))
    np.testing.assert_allclose(H @ inverse, np.eye(27), atol=3e-13)


def test_all_physical_routes_and_compensated_multitarget_standards():
    meters, _ = circuit.audit()
    board = m.flight.federation(0)
    profiles = m.regions(board)
    assert len(profiles) == 60
    result = physical.audit(board, profiles, meters)
    m.compare(result, json.loads(json.dumps(result)))
    assert result['all_routes_checked'] == 3600
    assert result['maximum_physical_standard_error'] < 2e-11
    ports = result['physical_ports']
    assert len(set(ports)) == 5
    for values, q in zip(result['physical_inputs'], result['region_charges'], strict=True):
        assert sum(values) == 5
        assert sum(values[:2]) - 2 == q


def test_complete_source_and_eligibility_counts_and_nontrivial_full_strength_controls():
    data = json.loads((m.HERE / 'results.json').read_text())
    for launch in data['launches']:
        for name, panels in launch['panels'].items():
            assert len(panels) == 60
            for panel in panels:
                np.testing.assert_allclose(sum(panel['admitted_region_charge']),
                                           panel['eligible_probability'], atol=2e-12)
                for record in panel['complete_records']:
                    assert len(record['counts']) == 14
                    assert sum(record['counts']) == m.SHOTS
                    assert all(p >= 0 for p in record['probabilities'])
        native = max(p['accepted_charged_fraction'] for p in launch['panels']['native_arrival'])
        model = max(p['accepted_charged_fraction'] for p in launch['panels']['matched_marginal_model'])
        changed = max(p['accepted_charged_fraction'] for p in launch['panels']['exchange_only_four_more_units'])
        assert native < 0.0011 and model > 0.25 and changed > 0.25


def test_invalid_probabilities_and_optimized_mismatch_reject():
    import pytest
    with pytest.raises(ArithmeticError):
        m.sample([-1e-13, 1 + 1e-13], np.random.default_rng(1))
    code = f"from {m.__package__} import measure; measure.compare({{'count':0}},{{'count':1}})"
    result = subprocess.run([sys.executable, '-O', '-B', '-c', code], capture_output=True, text=True)
    assert result.returncode != 0 and 'ArithmeticError' in result.stderr
