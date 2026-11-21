import json

import numpy as np
import pytest

from model.benchmarks.eight_slot_torus import law
from model.benchmarks.eight_slot_torus.electron_wall import walk as w
from model.benchmarks.eight_slot_torus.electron_wall.__main__ import HERE, digest, replay
from model.benchmarks.eight_slot_torus.electron_wall.gauge import flux_register
from model.benchmarks.eight_slot_torus.electron_wall.native import native_rows


def test_published_counts_all_settings():
    receipt = json.loads((HERE / "counts.json").read_text())
    assert receipt["sources"] == digest()
    _, result = replay(receipt)
    assert result["settings"] == 192


def test_reader_rejects_tampering():
    receipt = json.loads((HERE / "counts.json").read_text())
    receipt["counts"][next(iter(receipt["counts"]))][0] += 1
    with pytest.raises(ValueError, match="replay differs"):
        replay(receipt)


def noise(shape):
    rng = np.random.default_rng(682)
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def test_route_coin_and_endpoint_reference():
    shape = (4, 3, 2)
    u = np.exp(1j * noise((3, *shape)).real)
    state = noise((15, *shape, 2, 3, 2))
    chi = np.exp(1j * noise(shape).real)
    for q in (-1, 1):
        once = w.route(state, u, q)
        np.testing.assert_allclose(w.route(once, u, q), state, atol=1e-13)
        evolved = w.step(state, u, q)
        assert np.linalg.norm(evolved) == pytest.approx(np.linalg.norm(state))
        gauge = chi[None, ..., None, None, None] ** q
        transformed = np.array([np.roll(chi, -1, axis=i) * u[i] * chi.conj() for i in range(3)])
        np.testing.assert_allclose(w.step(state * gauge, transformed, q), evolved * gauge, atol=1e-13)
    with pytest.raises(ValueError):
        w.step(state, 2 * u, 1)


def test_all_carrier_and_dark_modes():
    length = 2
    dim = 15 * 4 * length
    # A twisted one-site torus realizes arbitrary Bloch momenta, including mixed signs.
    k = np.array([.37, -.81, .22])
    conn = np.exp(-1j * k)[:, None, None, None]
    shape = (15, 1, 1, 1, 2, length, 2)
    eye = np.eye(dim, dtype=complex)
    full = np.column_stack([w.step(e.reshape(shape), conn, 1).ravel() for e in eye])
    embed = np.kron(w.ROOT[:, None], np.eye(4 * length))
    d = w.symbol(k, length)
    np.testing.assert_allclose(embed.T @ full @ embed, d, atol=1e-13)
    np.testing.assert_allclose(full.conj().T @ full, eye, atol=1e-13)
    roots = np.linalg.eigvalsh(d)
    phases = np.arccos(roots)
    actual = np.linalg.eigvals(full)
    for sign in (1, -1):
        assert np.count_nonzero(abs(actual - sign) < 1e-10) == 26 * length
    mobile = actual[(abs(actual - 1) > 1e-10) & (abs(actual + 1) > 1e-10)]
    np.testing.assert_allclose(np.sort(np.angle(mobile)), np.sort(np.r_[-phases, phases]), atol=1e-12)


def test_wall_end_readout_in_time_and_all_private_modes():
    length = 4
    dim = 4 * length
    conn = w.links((1, 1, 1))
    site = np.eye(dim, dtype=complex)[0]
    signal = w.prepare(site.reshape(1, 1, 1, 2, length, 2))
    eigenvalues, vectors = np.linalg.eigh(w.symbol(np.zeros(3), length))
    weights = abs(vectors.conj().T @ site) ** 2
    unique = np.unique(eigenvalues.round(12))
    assert len(unique) == 2 * length
    for root in unique:
        assert weights[np.isclose(eigenvalues, root)].sum() > 1e-4
    rates = np.arccos(eigenvalues)
    for service in range(41):
        overlap = np.vdot(site, w.project(signal).ravel())
        assert overlap == pytest.approx(weights @ np.cos(service * rates), abs=1e-12)
        probs = w.pointer(overlap)
        assert probs.sum() == pytest.approx(1.)
        assert 2 * ((probs[0] - probs[2]) + 1j * (probs[1] - probs[3])) == pytest.approx(overlap)
        signal = w.step(signal, conn, 1)


def test_wall_controls_and_dispersion():
    for length in (2, 3, 4, 6, 8):
        props = w.rest_properties(length)
        mass = props["singular_mass"]
        assert mass < 1
        assert .75 < props["wall_endpoint_weight"] < 1
        assert np.linalg.svd(w.dimer(length, 1, True), compute_uv=False)[-1] == pytest.approx(1.)
        for corner in np.ndindex(2, 2, 2):
            k = np.pi * np.array(corner)
            plain = abs(np.linalg.eigvalsh(w.symbol(k, length, wilson=False))).min()
            assert 12 * plain == pytest.approx(mass)
            if any(corner):
                assert 12 * abs(np.linalg.eigvalsh(w.symbol(k, length))).min() >= 1 - 1e-12
        p = mass / 1000
        e = np.arcsin(abs(np.linalg.eigvalsh(w.symbol([p, 0., 0.], length))).min())
        curvature = (e - props["rest_rate"]) / p ** 2
        assert curvature == pytest.approx(props["curvature"], rel=5e-5)


def turn(state, axis):
    """Half-turn of the whole charged walk, including port labels and the spin frame."""
    slots = np.arange(15)
    for other in range(3):
        if other != axis:
            state = np.take(state, (-np.arange(state.shape[1 + other])) % state.shape[1 + other],
                            axis=1 + other)
            for start in (2 * other, 6 + 2 * other):
                slots[start:start + 2] = slots[start:start + 2][::-1]
    return state[slots] @ (-1j * w.PAULI[axis]).T


def test_realspace_rotations_and_spinor_loop():
    shape = (3, 3, 4)
    state = noise((15, *shape, 2, 3, 2))
    u = w.links(shape)
    for axis in (0, 1, 2):
        np.testing.assert_allclose(w.step(turn(state, axis), u, 1), turn(w.step(state, u, 1), axis),
                                   atol=1e-13)
    # Each inverse is minus its corresponding half turn.
    loop = turn(turn(turn(turn(state, 0), 1), 0), 1)
    np.testing.assert_allclose(loop, -state)
    overlap = np.vdot(state, loop) / np.vdot(state, state)
    np.testing.assert_allclose(w.pointer(overlap), [0, .25, .5, .25, 0], atol=1e-12)


def test_spin_charge_reversals_and_measured_moment():
    length = 3
    props = w.rest_properties(length)
    errors = []
    for period in (63, 251, 1001):
        shape = (period, 1, 1)
        p = 2 * np.pi / period
        a = .125
        def odd(charge, spin, spectator=False):
            args = dict(spectator=spectator)
            pos = w.read_grating(shape, length, charge, a, spin, **args)
            neg = w.read_grating(shape, length, charge, -a, spin, **args)
            return (pos - neg) / 2
        signal = odd(1, 0)
        np.testing.assert_allclose(signal, -odd(1, 1), atol=1e-13)
        np.testing.assert_allclose(signal, -odd(-1, 0), atol=1e-13)
        np.testing.assert_allclose(odd(1, 0, True), 0., atol=1e-13)
        _, r = w.wall_mode(length, p, 0)
        assert np.sin(p) > 0 and signal[0].imag < 0
        b = np.mean(np.sin(a * np.cos(2 * p * np.arange(period))) *
                    np.cos(2 * p * np.arange(period)))
        coefficient = abs(signal[0]) / (2 * p * b * np.sqrt(1 - (r / 12) ** 2))
        curvature = (np.arcsin(r / 12) - props["rest_rate"]) / p ** 2
        estimated = 2 * coefficient / curvature
        errors.append(abs(estimated - props["g"]))
    assert errors[2] < errors[1] < errors[0]
    assert errors[2] < .001


def test_private_blind_loophole():
    # Reflection exchanges two private harmonics and their private directional slots.
    period = 5
    def reflection(state):
        slots = [0, 1, 2, 3, 4, 5, 7, 6]
        return np.take(state, (-np.arange(period)) % period, axis=3)[..., slots]
    seed = noise((3, 3, 3, 1, 8))
    positive = seed * np.exp(2j * np.pi * np.arange(period) / period)[None, None, None, :, None]
    negative = reflection(positive)
    assert abs(np.vdot(positive, negative)) < 1e-11
    u = lambda state: law.step(state, law.LAWS["torus"])
    np.testing.assert_allclose(u(negative), reflection(u(positive)), atol=1e-13)
    np.testing.assert_allclose((abs(positive) ** 2).sum(axis=(3, 4)),
                               (abs(negative) ** 2).sum(axis=(3, 4)))
    # Distinguishing the harmonics by a private translation exposes opposite phases.
    for sign, state in ((1, positive), (-1, negative)):
        np.testing.assert_allclose(np.roll(state, 1, axis=3), np.exp(-sign * 2j * np.pi / period) * state)
    with pytest.raises(ValueError):
        native_rows(2)


def test_three_site_gauss_records():
    plus, leakage = flux_register(1, services=1)
    assert leakage < 1e-24
    # The plus and minus x sockets each have both kinetic and Wilson contributions.
    assert plus[1, 0, 0, 1] == pytest.approx(1 / 12)
    assert plus[0, 0, -1, 2] == pytest.approx(1 / 12)


def test_period_one_private_difference_is_an_invariant_local_line():
    seed = np.zeros((3, 3, 3, 1, 8), complex)
    seed[1, 1, 1, 0, 6:] = [1, -1]
    np.testing.assert_array_equal(law.coin(seed), -seed)
    np.testing.assert_array_equal(law.route(seed), -seed)
    np.testing.assert_array_equal(law.step(seed, law.LAWS["torus"]), seed)


def test_spinful_eight_port_control_has_scalar_offset_and_eight_cones():
    from model.benchmarks.eight_slot_torus.electron_wall.port_control import step

    shape = (4, 4, 4, 4)
    coords = np.indices(shape)
    conn = w.links(shape[:3])
    for k in [np.multiply(np.pi, key + (0,)) for key in np.ndindex(2, 2, 2)] + [
            np.array([1, -1, 0, 1]) * np.pi / 2]:
        wave = np.exp(1j * np.einsum("i,i...->...", k, coords)) / np.sqrt(np.prod(shape))
        basis = wave[..., None] * law.EQUAL
        reduced = np.empty((2, 2), complex)
        for spin in (0, 1):
            out = step(basis[..., None] * np.eye(2)[spin], conn, 1)
            reduced[:, spin] = np.einsum("...i,...ij->j", basis.conj(), out, optimize=True)
        expected = (np.cos(k[3]) * np.eye(2) + np.einsum("i,ijk->jk", np.sin(k[:3]), w.PAULI)) / 4
        np.testing.assert_allclose(reduced, expected, atol=1e-12)


def test_exchange_completions_against_full_two_particle_states():
    from model.benchmarks.eight_slot_torus.electron_wall.exchange import histories, site_pairs

    # Construct the two amplitudes explicitly and trace the detector sockets.
    shape = (5, 1, 1)
    states = []
    for sign in (-1, 1):
        p = sign * 2 * np.pi / shape[0]
        seed, _ = w.wall_mode(2, p, 0)
        state = w.step(w.prepare(w.wave(shape, p, seed)), w.links(shape), -1)
        states.append(state)
    records = site_pairs(*states)
    a, b = (s.swapaxes(0, 1).reshape(shape[0], -1) for s in states)
    direct = np.einsum("xa,yb->xayb", a, b)
    swapped = direct.transpose(2, 3, 0, 1)
    for name, sign in (("symmetric", 1), ("antisymmetric", -1)):
        prob = (abs(direct + sign * swapped) ** 2).sum(axis=(1, 3)) / 2
        unordered = np.diag(np.diag(prob)) + np.triu(prob + prob.T, 1)
        np.testing.assert_allclose(records[name], unordered, atol=1e-13)
    for reads in histories().values():
        assert np.trace(reads["symmetric"]) > np.trace(reads["tagged"]) > np.trace(reads["antisymmetric"])
        np.testing.assert_allclose((reads["symmetric"] + reads["antisymmetric"]) / 2, reads["tagged"])


def test_exchange_detector_resolves_transverse_sites():
    from model.benchmarks.eight_slot_torus.electron_wall.exchange import site_pairs

    a = np.zeros((15, 1, 2, 2, 2, 3, 2), complex)
    b = np.zeros_like(a)
    a[0, 0, 0, 1, 0, 0, 0] = 1
    b[0, 0, 1, 0, 0, 0, 0] = 1
    expected = np.zeros((4, 4))
    expected[1, 2] = 1
    for detector in site_pairs(a, b).values():
        np.testing.assert_array_equal(detector, expected)
    a.flat[0] = np.nan
    with pytest.raises(ValueError, match="finite amplitudes"):
        site_pairs(a, b)


def test_exact_root_and_anomaly_enclosures():
    from decimal import Decimal, localcontext
    from model.benchmarks.eight_slot_torus.electron_wall.wall_pole import continuants, enclose

    assert continuants(2) == ([1, -6, 1], [4, -4, 0])
    with localcontext() as ctx:
        ctx.prec = 80
        # At length two s(a)=sqrt(1+a^2)-1.
        mass = Decimal(2).sqrt() - 1
        t = mass / Decimal(2).sqrt()
        item = enclose(2)
        unit = 2 ** item["denominator_power"]
        for key, value in (("mass_squared", mass**2), ("magnetic_deficit", t/(1+t))):
            lo, hi = item[key]
            assert lo < unit * value < hi
    for length in (3, 4, 6, 8, 32, 64):
        row = enclose(length)
        more = enclose(length, 56)
        for key in ("mass_squared", "magnetic_deficit"):
            lo, hi = row[key]
            assert 0 < lo < hi < 2 ** row["denominator_power"]
            assert lo * 256 <= more[key][0] < more[key][1] <= hi * 256
