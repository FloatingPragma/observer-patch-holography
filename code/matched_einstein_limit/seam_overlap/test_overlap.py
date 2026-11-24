import json
from fractions import Fraction

import numpy as np
from scipy.linalg import expm, logm

from . import measure, protected


def test_blank_write_and_sharp_coin_read_have_different_protected_support():
    seam = protected.seam
    law = seam.law
    writes = law._write_pattern()  # integer columns for sqrt(2) W
    balanced = np.eye(law.N_SRC, dtype=int)[law.SRC.index((1, 1))]
    zero_even = np.kron(seam.coin_basis(0, 1) @ [1, 1], balanced)
    np.testing.assert_array_equal(writes[:, law.SRC.index((1, 1))], zero_even)
    written = writes[:, law.SRC.index((2, 0))]
    coin = np.repeat([c for c, _ in law.LABELS], law.N_SRC)
    for parity in (1, -1):
        mode = np.kron(seam.coin_basis(2, parity) @ [1, -1], balanced)
        assert int(mode @ written) == 0
        norm = int(mode @ mode) * int(written @ written)
        weight = sum(Fraction(int(mode @ (written * (coin == c)))**2, norm)
                     for c in (0, 1))
        assert weight == Fraction(1, 4)
    # For odd total load the two written coin alternatives carry distinct endpoint loads.
    odd = writes[:, law.SRC.index((1, 0))].reshape(law.N_LAB, law.N_SRC)
    difference = odd[law.LABELS.index((0, 1))] - odd[law.LABELS.index((1, 1))]
    assert Fraction(int(difference @ difference), 2 * int((odd * odd).sum())) == Fraction(1, 2)


def test_operator_identity_for_coherent_states():
    # X has zero partial traces, and does not commute with h.
    sx = np.array([[0., 1.], [1., 0.]])
    sz = np.diag([1., -1.])
    h = np.kron(sz, np.eye(2)) + 2*np.kron(sx, sz) + np.kron(sz, sx)
    x = (np.kron(sx, sz) + np.kron(sx, sx))/16
    assert np.linalg.norm(x @ h - h @ x) > 0
    for state in (np.eye(4)/4 + x/2, np.eye(4)/4 - x/2):
        assert min(np.linalg.eigvalsh(state)) > 0
    sigma = expm(-measure.BETA*h)
    sigma /= sigma.trace()
    view = sigma.reshape(2, 2, 2, 2)
    marginals = np.kron(view.trace(axis1=1, axis2=3), view.trace(axis1=0, axis2=2))
    overlap = -logm(sigma) + logm(marginals)
    np.testing.assert_allclose(np.trace(x @ overlap)/measure.BETA, np.trace(x @ h), atol=1e-13)


def test_protected_columns_and_reference_selection():
    report = protected.audit()
    blind = report['protected_logical_registers']
    assert Fraction(blind['energy_contrast']) == Fraction(1, 4)
    assert blind['protected_input_columns'] == 36
    assert blind['max_integer_scalar_action_residual'] == 0
    region = report['regional']
    errors = [abs(row['rescaled_stress_slope']-report['rescaled_limit']) for row in region]
    assert all(new < old/12 for old, new in zip(errors, errors[1:]))
    p = protected.distribution(.5)
    product = np.outer(p.sum(axis=1), p.sum(axis=0))
    delta = protected.distribution(.4)-protected.distribution(-.4)
    np.testing.assert_allclose((delta*np.log(product)).sum(), 0, atol=1e-15)


def test_analytic_expanding_solution():
    for energy in (0, 4, 8, 12, 16):
        k, vacuum, v0 = measure.radar.COUPLING, measure.radar.VACUUM, measure.V0
        source = measure.canonical_source(energy)
        for t in (0, .35, .6, 1.1):
            z = np.sqrt(3*vacuum)*t/2
            volume = (np.sqrt(v0)*np.cosh(z) + np.sqrt(v0+k*energy/vacuum)*np.sinh(z))**2
            np.testing.assert_allclose(v0*np.prod(source(t)[:3]), volume, rtol=3e-12)


def test_independent_meter_and_changed_gravity():
    saved = json.loads((measure.HERE/'results.json').read_text())
    a, b = [saved['arms'][key] for key in ('interacting', 'protected')]
    for arm, exact in [(a, 8), (b, 0)]:
        np.testing.assert_allclose(arm['exact_overlap_energy'], exact, atol=1e-13)
        assert abs(arm['overlap_meter']-exact) < 3*arm['overlap_counting_se']
        table = np.array(arm['spectrum'])
        for mean, s in zip(arm['source_energies'], measure.SIGNS, strict=True):
            np.testing.assert_allclose(((1+s*measure.CYCLE)*table).sum()/4, mean)
    held = a['comparisons'][-1]
    k = held['ordinary_coupling']
    assert abs(k-1.6) < .000002
    for row in held['held']:
        assert abs(row['integrated_temporal']-8*k) < .00005
        np.testing.assert_allclose(row['integrated_axial_null'], 8*k, atol=.0005, rtol=0)
    changed = saved['selective_gravity_control']['comparisons'][-1]
    assert changed['ordinary_coupling'] == k
    for universal, selective in zip(held['held'], changed['held'], strict=True):
        assert universal['forecast'] == selective['forecast']
        assert selective['integrated_temporal'] == 0
        assert selective['integrated_axial_null'] == [0., 0., 0.]


def test_smaller_rims_at_equal_volume_and_converging_coupling():
    from . import rim

    saved = json.loads((measure.HERE/'rim.json').read_text())
    errors = []
    for row in saved['rows']:
        target = 4*np.pi*row['radius']**3/3
        a, b = row['high'], row['low']
        np.testing.assert_allclose([a['volume'], b['volume']], target, atol=1e-13, rtol=3e-13)
        assert a['area'] < b['area']
        residual = abs(row['inferred_coupling']-measure.radar.COUPLING)
        errors.append(residual)
    assert all(new < .26*old for old, new in zip(errors, errors[1:]))
    assert errors[-1] < .0015
    for e in (4, 8, 12):
        metric = rim.metric_at_preparation(e)
        observed = measure.canonical_source(e)
        for t in (0, .06, .21):
            scale = np.exp(metric(rim.balls.APPOINTMENT+t)[:3])
            np.testing.assert_allclose(scale, observed(t)[:3], rtol=2e-12)
