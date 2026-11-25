import json
from itertools import combinations

import numpy as np
import pytest

from . import certificate, measure


@pytest.fixture(scope="module")
def native():
    return measure.load_native()


def stored():
    return json.loads((measure.HERE / "results.json").read_text())


def test_stationary_sources_and_held_local_records(native):
    channels, _ = native
    # Fresh sparse native radar kernels supply all 12 directed ranges independently
    # of the integer grid on which the matching route distributions were selected.
    effects = measure.station.radar_menus()[6]
    archive = stored()
    for mode, comparison in archive["comparisons"].items():
        low, high = comparison["endpoints"][-2:]
        for item in comparison["endpoints"]:
            weights = np.array(item["routes"])
            assert abs(weights.sum() - 1) < 1e-14 and weights.min() > 0.06
            states = []
            for ch, p in zip(channels, measure.ODDS, strict=True):
                val, vec = np.linalg.eig(np.einsum("d,dij->ij", weights, ch))
                k = np.argmin(abs(val - 1))
                assert abs(val[k] - 1) < 4e-11 and np.max(abs(np.delete(val, k))) < 0.92
                x = vec[:, k].reshape(4, 4)
                states.append(p * x / np.trace(x))
            states = np.array(states)
            ranges = measure.station.pair(effects, states)
            np.testing.assert_allclose(ranges, item["ranges"], rtol=0, atol=5e-11)
            labels = np.array(item["parity_blind_paulis"])
            np.testing.assert_allclose(
                2 * np.einsum("kij,bji->k", measure.BASIS, states).real,
                labels,
                atol=2e-10,
                rtol=0,
            )
            s = sum(
                -ev * np.log(ev) for r in states for ev in np.linalg.eigvalsh(r)
            ) + np.log(20)
            assert abs(s - item["closed_flight_entropy"]) < 2e-10
            assert item["pure_start_320_cycle_residual"] < 2e-10
            if mode == "six_edges":
                np.testing.assert_allclose(
                    item["edges"], archive["reference"]["edges"], atol=2e-12, rtol=0
                )
            for key in ("area", "volume"):
                assert abs(item[key] - archive["reference"][key]) < 2e-12
        assert (
            high["central_bridge_energy_current"] < low["central_bridge_energy_current"]
        )
        assert high["four_station_memory_entropy"] > low["four_station_memory_entropy"]
        clicks = comparison["independent_pauli_read"]
        rates = np.array(clicks["plus_counts"]) / clicks["trials_per_source"]
        width = clicks["simultaneous_99pct_hoeffding_halfwidth"]
        assert np.max(abs(rates - clicks["plus_probabilities"])) < width
        assert rates[0] - rates[1] > 2 * width


def test_reversed_raw_clocks_and_trivial_source_control(native):
    channels, effects = native
    state, _ = measure.prepare(channels, measure.UNIFORM)
    ranges = measure.station.pair(effects, state)
    reversed_ranges = measure.station.square(ranges).T[measure.station.ORDER]
    assert max(abs(ranges - reversed_ranges)) > 0.001
    # Reciprocity is not presumed. Reversing every route flips the fitted clock
    # gains; the resulting lengths are allowed to change. An exactly reciprocal
    # artificial input, on the other hand, must retain its six distinct edges.
    symmetric = np.array(
        [
            [1.0, 1.2, 1.3, 1.4],
            [1.2, 1.0, 1.5, 1.6],
            [1.3, 1.5, 1.0, 1.7],
            [1.4, 1.6, 1.7, 1.0],
        ]
    )
    np.testing.assert_allclose(
        measure.edges(symmetric[measure.station.ORDER]),
        [symmetric[a, b] for a, b in combinations(range(4), 2)],
        rtol=0,
        atol=2e-15,
    )
    dummy = np.repeat(channels.mean(axis=1, keepdims=True), 12, axis=1)
    arbitrary = np.arange(1, 13, dtype=float) / 78
    reused, _ = measure.prepare(dummy, arbitrary)
    np.testing.assert_allclose(reused, state, atol=2e-14, rtol=0)


def test_certified_entropy_and_native_energy_current(native):
    saved = json.loads((measure.HERE / "certificate.json").read_text())
    replay = certificate.verify(*native, stored())
    for mode in ("area_volume", "six_edges"):
        for key in (
            "entropy_increase_nats",
            "central_energy_current_change",
            "maximum_absolute_edge_difference",
            "absolute_area_difference_upper",
            "absolute_volume_difference_upper",
        ):
            low, high = (float(x) for x in replay[mode][key].strip("[]").split(","))
            archived = [float(x) for x in saved[mode][key].strip("[]").split(",")]
            np.testing.assert_allclose([low, high], archived, rtol=1e-8, atol=1e-14)
            if key == "entropy_increase_nats":
                assert low > 0
            if key == "central_energy_current_change":
                assert high < -0.0005


def independent_channels():
    """Propagate all four starts, subtract earlier hits, and sum terminal Kraus records."""
    from scipy.sparse.linalg import expm_multiply

    board = measure.station.iso.pairing.formation.federation(1)
    ports = measure.station.iso.pairing.formation.region(board)
    idx = measure.bridge.beam.indices(ports)
    start = np.zeros((4 * board.n, 16))
    start[idx] = np.eye(16)
    out = np.zeros((2, 12, 16, 16), complex)
    for b in (0, 1):
        for charge in (-1, 1):
            h = measure.native_precision.generator(board, ports, charge, bool(b))
            u = expm_multiply(
                -0.5j * h, start, start=0.0, stop=6.0, num=13, traceA=0.0
            ).reshape(13, -1, 4, 4, 4)
            local = u[:, ports].transpose(0, 1, 3, 2, 4)
            starts = {}
            survived = {}
            for i, j in measure.reuse.ROUTES:
                hits = measure.station.iso.ranging.first_hits(local, i, j)
                starts[i, j] = [np.kron(k, k.conj()) for k in hits]
                survivors = []
                for t in range(13):
                    v = u[t, :, :, i].copy()
                    for step in range(1, t + 1):
                        v -= u[t - step, :, :, j] @ hits[step]
                    particles = v.reshape(-1, 16)
                    gram = particles.T @ particles.conj()
                    survivors.append(
                        gram.reshape(4, 4, 4, 4).transpose(0, 2, 1, 3).reshape(16, 16)
                    )
                survived[i, j] = survivors
            for k, (i, j) in enumerate(measure.reuse.ROUTES):
                total = survived[i, j][12].copy()
                for step in range(1, 13):
                    arrival = starts[i, j][step]
                    total += survived[j, i][12 - step] @ arrival
                    for back in range(1, 13 - step):
                        total += starts[j, i][back] @ arrival
                out[b, k] += total / 2
    return out


def test_held_release_and_transfer_on_independent_propagator():
    from . import release

    receipt = stored()
    bank = independent_channels()
    effects = measure.station.radar_menus()[6]
    replay = release.study(bank, effects, receipt)
    archived = json.loads((measure.HERE / "release.json").read_text())
    for tested, saved in zip(replay["rows"], archived["rows"], strict=True):
        np.testing.assert_allclose(
            tested["two_endpoint_geometry"],
            saved["two_endpoint_geometry"],
            atol=2e-11,
            rtol=0,
        )
        assert tested["maximum_stationary_sham_one_cycle_drift"] < 1e-10
        assert (
            abs(
                tested["volume_second_difference_contrast"]
                - saved["volume_second_difference_contrast"]
            )
            < 4e-13
        )
    assert abs(replay["rows"][1]["observed_minus_forecast"]) < 1e-9
    for held in replay["rows"][2:]:
        ratio = (
            held["forecast_volume_second_difference_contrast"]
            / held["volume_second_difference_contrast"]
        )
        assert 23.8 < ratio < 24.2
