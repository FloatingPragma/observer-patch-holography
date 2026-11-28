import json

import numpy as np
from scipy import sparse
from scipy.linalg import eigh, expm

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    area_survey,
    contact_timing,
    full_repair,
    repair_exposure,
    weak_timing,
    exchange_geometry,
    poisson_precision,
    poisson_response,
    weak_source,
)


def test_uniform_poisson_derivative_and_permanent_mode_screening():
    # Six-cycle, delete one edge, write its two endpoints.
    n = 6
    adj = np.roll(np.eye(n), 1, axis=1) + np.roll(np.eye(n), -1, axis=1)
    lap = np.diag(adj.sum(axis=1)) - adj
    m = np.array([1, 1, 0, 0, 0, 0])
    inc = np.array([1, -1, 0, 0, 0, 0])
    h0 = lap / 2
    v = -np.outer(inc, inc) / 2 - poisson_response.K * np.diag(m)
    delta = poisson_response.weak_density(sparse.csr_matrix(lap), m)

    def occupation(x):
        return eigh(h0 + x * v)[1][:, 0] ** 2

    eps = 1e-6
    np.testing.assert_allclose(
        n * (occupation(eps) - occupation(-eps)) / (2 * eps), delta, atol=3e-8, rtol=0
    )
    live = lap - np.outer(inc, inc)
    e, modes = eigh(h0 + v)
    assert sum(e < 0) == 2
    for mode in modes[:, :2].T:
        p = mode**2
        assert max((live @ p)[2:]) < 0
        f = p + (np.diag(live.diagonal()) - live) @ p
        assert max((live @ f)[3:5]) < 0


def test_clock_average_agrees_with_literal_cycle_evolution():
    h = np.array([[0.8, -0.5], [-0.5, 0.55]])
    u = expm(-0.034j * h)
    v = np.array([0.6, 0.8])
    x = v.copy()
    rho = np.zeros((2, 2), complex)
    for _ in range(27):
        rho += np.outer(x, x.conj()) / 27
        x = u @ x
    e, q = eigh(h)
    coeff = q.T @ v
    expected = (
        q
        @ (
            np.outer(coeff, coeff)
            * weak_source.clock_average(-0.034 * (e[:, None] - e[None, :]), 27)
        )
        @ q.T
    )
    np.testing.assert_allclose(rho, expected, atol=5e-14)


def test_held_reading_follows_from_source_with_independent_detectors():
    root = poisson_response.HERE
    tape = json.loads((root / "mediated_release.json").read_bytes())
    standard = json.loads((root / "weak_source.json").read_bytes())
    raw = json.loads((root / "poisson_precision.json").read_bytes())
    b, keep, marks, y, h0, _, _, _, _ = weak_source.fixture(tape)
    h, i, j = poisson_precision.flight.detector(b, keep)
    meters = poisson_precision.meters(b, keep, marks, y, h0, h, i, j)
    original = poisson_precision.decode(raw, standard, tape, meters)
    # Corrupting every terminal record by directing it to the last shell
    # leaves the independently inferred source acceleration identical.
    for arm in raw["position"].values():
        for source in arm.values():
            source["later"] = [0] * (len(source["later"]) - 1) + [
                poisson_precision.TRIALS
            ]
    changed = poisson_precision.decode(raw, standard, tape, meters)
    for name in original:
        assert [r["from_source"] for r in changed[name]] == [
            r["from_source"] for r in original[name]
        ]
        assert all(r["from_source"][0] > 0 for r in original[name])
    assert all(r["terminal_interval"][0] > 0 for r in original["written"])


def test_native_swap_and_visible_ordinary_records():
    swap = exchange_geometry.contact_swap()
    np.testing.assert_array_equal(swap @ swap, np.eye(9))
    record = json.loads(
        (exchange_geometry.HERE / "exchange_geometry.json").read_bytes()
    )
    assert all(
        row["conditional_swap_visibility"][1] < 1
        for bank in record["decoded"]["arms"].values()
        for row in bank
    )
    assert record["decoded"]["contrast"][0]["residual"][0] > 0


def test_degree_class_area_weights_cannot_supply_square_law():
    from fractions import Fraction as Q

    data = area_survey.run()
    upper = data["charts"][-1]["comparisons"][-1]
    assert upper["excluded_for_both_types"]
    lo, hi = map(Q, upper["arbitrary_contact_areas"])
    assert (
        lo
        <= Q(upper["inner_to_outer_mean_area_ratio"])
        <= hi
        < Q(upper["square_law_lower"])
    )


def test_time_order_and_source_width_can_change_the_exterior_sign():
    from fractions import Fraction

    root = contact_timing.HERE
    tape = json.loads((root / "mediated_release.json").read_bytes())
    short = json.loads((root / "contact_timing.json").read_bytes())
    weak = json.loads((root / "weak_timing.json").read_bytes())
    thin = json.loads((root / "repair_exposure.json").read_bytes())
    full = json.loads((root / "full_repair.json").read_bytes())
    q = Fraction(weak["conditional_active_fraction"])
    assert q * short["active_ticks"] == Fraction(short["ticks"], 10000)
    assert weak_timing.run(short, tape, weak) == weak["decoded"]
    assert (
        weak["decoded"]["rare_recorded"]["2"][1]
        < 0
        < weak["decoded"]["recorded_short"]["2"][0]
    )
    assert repair_exposure.run(tape, short, thin) == thin["decoded"]
    assert full_repair.run(tape, thin, full) == full["decoded"]
    for n, arms in full["decoded"].items():
        for tag, bank in arms.items():
            for r in ("4", "6"):
                assert bank[r][1] < 0 < thin["decoded"][n][tag][r][0]


def test_passive_pair_source_and_three_load_release():
    from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
        pair_contact,
        pair_meter,
        pair_motion,
        pair_probe,
    )

    saved = json.loads((pair_contact.HERE / "pair_contact.json").read_bytes())
    measured = pair_contact.decode(saved)
    assert measured["full_support_nats"][0] > 0
    assert max(
        measured["from_source_fringe_shift"][0], measured["held_fringe_shift"][0]
    ) < min(measured["from_source_fringe_shift"][1], measured["held_fringe_shift"][1])
    frequencies = measured["ordinary_energy_bands"]
    assert frequencies["1,1"][0] <= -32 <= frequencies["1,1"][1]
    assert frequencies["-1,1"][0] <= 0 <= frequencies["-1,1"][1]
    assert len(pair_probe.check_generator()) == 6
    local = json.loads((pair_contact.HERE / "pair_meter.json").read_bytes())
    for row in local["readings"]:
        assert (
            pair_meter.estimate(local["reference"], local["clocks"], *row["bright"])
            == row["from_ordinary_fringe"]
        )
    # Off-diagonal phases survive propagation; compare all entries up to its global phase.
    h = sparse.csr_matrix([[1.5, -0.5, 0], [-0.5, 1, -0.5], [0, -0.5, 2.0]])
    psi = np.array([1, 1j, 0]) / np.sqrt(2)
    cases = [(0,), (-2,)]
    waves, cert = pair_motion.chebyshev(h, [1], cases, psi, [1, 4])
    center = sum(cert["spectral_bounds"]) / 2
    for t, bank in zip((1, 4), waves, strict=True):
        for j, (shift,) in enumerate(cases):
            exact = h.toarray()
            exact[1, 1] += shift
            np.testing.assert_allclose(
                bank[:, j] * np.exp(-1j * center * t),
                expm(-1j * t * exact) @ psi,
                rtol=0,
                atol=5e-13,
            )
