import json
from pathlib import Path

import numpy as np
import pytest

from . import certificate, enclosure, station


def test_native_swap_loading_and_unread_sign():
    gate = station.phase.marked_gate(1,8).reshape(4,4)
    z = np.diag([1,-1])
    expected = -np.kron(z,z) @ np.eye(4)[[0,2,1,3]]
    np.testing.assert_allclose(gate,expected,atol=1e-14,rtol=0)
    np.testing.assert_allclose(gate @ gate,np.eye(4),atol=1e-14,rtol=0)
    station.iso.ranging.ordinary_write()


def test_full_source_station_and_held_instruments():
    report = station.run()
    assert report['target_volume'] == pytest.approx(2.7634718367887,abs=2e-12)
    assert report['reports'][6]['max_residual'] < 1e-10
    assert report['constrained_eigenvalues'][-1] < -3
    assert report['hessian_check'] < 1e-8
    assert report['preparation']['min_kept'] > .007
    assert report['preparation']['min_dephased'] > .003
    assert report['preparation']['weight_sum'] == pytest.approx(1)
    assert report['preparation']['max_error'] < 1e-12
    assert abs(report['held']['hidden_entropy']) < 1e-12
    assert abs(report['held']['shape_entropy']) > 0.6
    a,b = report['held']['finite']
    assert b['max_volume_error'] < 1e-12 and abs(b['residual']) < .26*abs(a['residual'])
    assert report['reports'][4]['max_residual'] > .72
    archive = json.loads(Path(station.__file__).with_name('results.json').read_text())
    np.testing.assert_allclose(archive['x_imag'],report['x_imag'],rtol=0,atol=2e-11)


def test_interval_jet_independent_of_range_fit():
    sample = np.linspace(2.9,3,12)
    a,v = enclosure.geometry_jet([enclosure.iv.mpf(float(x)) for x in sample])
    np.testing.assert_allclose([float(a.v.mid),float(v.v.mid)],station.geometry(sample),atol=1e-12,rtol=0)
    g = np.array([[float(x.mid) for x in jet.g] for jet in (a,v)])
    np.testing.assert_allclose(g,station.range_gradient(sample),atol=1e-12,rtol=0)


def test_native_exact_existence_certificate():
    data = certificate.verify()
    assert data['geometric_curvature_bound'] < 1
    assert data['necessary_log_error_floor_lower'] > .00254
    assert data['state_HS_error_upper'] == 2e-7
    assert data['global_survey']['boxes'] == 64
    assert data['global_survey']['curvature_upper_millionths'] < 720000
    assert float(data['signed_transfers']['right_adjuster_volume'].strip('[]').split(',')[1]) < -2e-5


def test_global_polynomials_match_independent_geometric_jet():
    from . import global_survey as global_proof

    ranges = [enclosure.iv.mpf(str(x)) for x in np.linspace(2.64, 2.98, 12)]
    log_ranges = [enclosure.iv.ln(x) for x in ranges]
    direction = np.array([enclosure.iv.mpf(i-5)/37 for i in range(12)], dtype=object)
    divided = direction / np.array(ranges, dtype=object)
    polynomials = global_proof.squared_measures()
    jets = []
    for measure in polynomials:
        value, gradient, second = global_proof.polynomial_jet(measure, log_ranges, np.zeros(12, dtype=object))
        root = enclosure.iv.sqrt(value)
        first = gradient @ divided / (2*root)
        curvature = (divided @ second @ divided) / (2*root) - first**2 / root
        jets.append((root, first, curvature))
    area, volume = enclosure.geometry_jet(ranges)
    for expected, group in ((area, jets[:4]), (volume, jets[4:])):
        for actual, oracle in zip(map(sum, zip(*group)), (expected.v, expected.g @ direction,
                                                          direction @ expected.h @ direction)):
            assert actual.a <= oracle.b and oracle.a <= actual.b


def test_correlated_boards_share_total_volume():
    from . import replicated

    data = replicated.measure()
    archived = json.loads(Path(replicated.__file__).with_suffix('.json').read_text())
    np.testing.assert_allclose(data['individual_volume_offsets'], archived['individual_volume_offsets'],
                               rtol=0, atol=1e-12)
    assert data['objective_gap_nats'] > data['lower_bound_nats'] > 0
    assert data['total_correlation_nats'] > 0
    pooled = data['pooled_volume_check']
    assert abs(pooled['volume_offset']) < 2e-14
    assert abs(pooled['separate_total_volume_offset']) > 1e-7
    assert pooled['objective_gap_nats'] > pooled['lower_bound_nats'] > 0
