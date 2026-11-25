import json
from itertools import product

import numpy as np

from .. import measure as baths
from . import apparatus as app
from . import measure as m


def test_detector_only_forecast():
    receipt = json.loads((m.HERE / "results.json").read_text())
    with np.load(m.HERE / "counts.npz") as saved:
        off, slope = m.decode_calibration(saved["standard_return_ticks"])
        assert saved["standard_return_ticks"].shape == (5, 2, 31, 12, 12)
        assert np.all(saved["standard_return_ticks"].sum(axis=-1) == m.TRIALS)
        for name, row in receipt["held"].items():
            head, target = (saved[name + "_" + k] for k in ("parked", "release"))
            assert np.all(head.sum(axis=(1, 2)) == m.TRIALS)
            assert np.all(target.sum(axis=-1) == m.TRIALS)
            est, err = m.predict(head, off, slope)
            truth = target @ app.SCORES / m.TRIALS
            assert np.all(abs(est - truth) <= err + m.RANGE_ERROR)
            for n, pos in enumerate(est):
                np.testing.assert_allclose(
                    baths.station.geometry(pos)[1],
                    row["forecast_volume"][n],
                    atol=3e-14,
                )
                assert max(
                    row["forecast_volume_band"][n][0],
                    row["independent_volume_band"][n][0],
                ) <= min(
                    row["forecast_volume_band"][n][1],
                    row["independent_volume_band"][n][1],
                )
            for key in ("contrast_four", "reset"):
                assert np.all(saved[name + "_" + key].sum(axis=-1) == m.TRIALS)
        assert set(receipt["held"]) == {
            "uniform",
            "balanced_plus",
            "balanced_minus",
            "mean_matched_plus",
            "mean_matched_minus",
        } | {f"route_{i:02d}" for i in range(12)}
        assert (
            len(
                {
                    tuple(receipt["held"][f"route_{i:02d}"]["route_weights"])
                    for i in range(12)
                }
            )
            == 12
        )


def test_positive_preparations_and_blind_remote_control():
    with np.load(m.HERE / "counts.npz") as saved:
        off, slopes = m.decode_calibration(saved["standard_return_ticks"])
        variant, _ = m.control_sources(slopes)
        matched = m.mean_matched(off, slopes)
        for pair in (variant, matched):
            for x, _, _, _ in pair.values():
                assert np.linalg.eigvalsh(x).min() >= -1e-15
                np.testing.assert_allclose(
                    x.trace(axis1=1, axis2=2), m.ODDS, atol=1e-15
                )
        plus, minus = [v[0] for v in variant.values()]
        np.testing.assert_allclose(plus.sum(axis=0), minus.sum(axis=0), atol=1e-15)
        x = np.einsum(
            "aij,bji->ba",
            m.PAULI,
            matched["mean_matched_plus"][0] - matched["mean_matched_minus"][0],
        ).real
        np.testing.assert_allclose(np.einsum("bar,ba->r", slopes[0], x), 0, atol=4e-16)
        means = [saved[k + "_release"] @ app.SCORES / m.TRIALS for k in matched]
        assert np.max(abs(means[0][0] - means[1][0])) < 2 * m.RANGE_ERROR
        assert np.max(abs(means[0][3] - means[1][3])) > 10 * m.RANGE_ERROR


def test_interval_volume():
    rng = np.random.default_rng(123)
    center = np.array(
        [2.65, 2.89, 2.97, 2.71, 2.91, 2.98, 2.997, 2.925, 2.77, 2.99, 2.95, 2.76]
    )
    radius = np.linspace(1e-6, 4e-6, 12)
    low, hi = m.volume_band(center, radius)
    for sign in product((-1.0, 1.0), repeat=12):
        assert low < baths.station.geometry(center + radius * sign)[1] < hi
    for _ in range(50):
        assert (
            low
            < baths.station.geometry(center + rng.uniform(-1, 1, 12) * radius)[1]
            < hi
        )


def test_independent_integer_radar():
    # The independent integer program performs monitored launches, resets
    # their arrivals, and accumulates survivor partial traces at every tick.
    channels, effects = baths.load_native()
    archive = json.loads((m.HERE / "results.json").read_text())
    with np.load(m.HERE / "counts.npz") as saved:
        off, slopes = m.decode_calibration(saved["standard_return_ticks"])
        uniform = channels.mean(axis=1)
        stands = m.standards()
        raw = saved["standard_return_ticks"] @ app.SCORES / m.TRIALS
        for a, n in enumerate(app.DELAYS):
            for b in range(2):
                for k, state in enumerate(stands[b]):
                    independent = baths.station.pair(
                        effects, app.evolve(uniform, state, n)
                    )
                    assert np.max(abs(raw[a, b, k] - independent)) < m.RANGE_ERROR
        states = m.schedules(channels)
        states.update(m.control_sources(slopes)[0])
        states.update(m.mean_matched(off, slopes))
        for label, (state, _, _, _) in states.items():
            counts = saved[label + "_release"] @ app.SCORES / m.TRIALS
            for a, n in enumerate(app.DELAYS):
                ranges = baths.station.pair(effects, app.evolve(uniform, state, n))
                assert np.max(abs(counts[a] - ranges)) < m.RANGE_ERROR
                vol = baths.station.geometry(ranges)[1]
                low, hi = archive["held"][label]["forecast_volume_band"][a]
                assert low <= vol <= hi


def test_exact_fiber_on_independent_sparse_flight():
    from . import exact_memory

    archive = json.loads((m.HERE / "exact_memory.json").read_text())
    channels, clocks = app.instruments(8)
    effects = np.einsum("brtij,t->rbij", clocks, app.SCORES)
    replay = exact_memory.study(channels, effects)
    for k in ("right_inverse_defect_upper", "future_range_operator_error_upper"):
        assert replay[k] < 1e-5
    proof = archive["directions_on_exact_initial_range_fiber"]["volume"]
    for doc in (proof, replay["directions_on_exact_initial_range_fiber"]["volume"]):
        low, high = doc["finite_volume_difference_interval"]
        assert 2.2e-6 < low < high < 6.4e-6
        assert doc["minimum_source_eigenvalue_lower"] > 0.01794
    B = np.einsum("rbij,aji->rba", effects, m.PAULI).real.reshape(12, 30) / 4
    v = np.array(proof["direction_numerators"]) / archive["scale"]
    # Here the static effects and release are independently computed by sparse
    # exponentiation and terminal Kraus summation. Compare the actual endpoints.
    q = np.linalg.qr(B.T)[0]
    corrected = v - q @ (q.T @ v)
    perturb = np.einsum("ba,aij->bij", corrected.reshape(2, 15), m.PAULI) / 4
    ranges = []
    for sign in (-1, 1):
        incoming = m.ODDS[:, None, None] * np.eye(4) / 4 + sign * perturb / 32
        ranges.append(
            (
                app.census(incoming, clocks) @ app.SCORES,
                app.census(app.evolve(channels.mean(axis=1), incoming, 4), clocks)
                @ app.SCORES,
            )
        )
    np.testing.assert_allclose(ranges[1][0], ranges[0][0], rtol=0, atol=2e-14)
    diff = (
        m.baths.station.geometry(ranges[1][1])[1]
        - m.baths.station.geometry(ranges[0][1])[1]
    )
    low, high = proof["finite_volume_difference_interval"]
    assert low < diff < high
