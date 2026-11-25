"""Check monitored trajectories, free source factorization and the archived protocols."""

from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.experiments.interval_echo_einstein_station.collar_pairing.ranging import (
    measure as legacy,
)
from . import live_radar as radar


def monitored_record(board, charge, sender, relay, state):
    """Evolve each surviving physical branch and project its terminal after each tick."""
    route = radar.source.formation.region(board)
    h = radar.source.system(board, route, charge)[0]
    scale = -1j * radar.STEP / 40
    u = scale * h
    trace = scale * h.diagonal().sum()
    source_rows, relay_rows = radar.station_rows(board, (sender, relay))
    live = np.zeros((h.shape[0], radar.TICKS + 1), complex)
    live[source_rows, 0] = state
    hits, returns, misses = [0.0], [0.0], [1.0]
    for t in range(1, radar.TICKS + 1):
        result = expm_multiply(u, live[:, :t], traceA=trace)
        at_relay = result[relay_rows, 0].copy()
        returned = result[source_rows, 1:].copy()
        hits.append(float(np.vdot(at_relay, at_relay).real))
        returns.append(float(np.vdot(returned, returned).real))
        result[relay_rows, 0] = 0.0
        result[source_rows, 1:] = 0.0
        live[:, :t] = result
        # Transfer only the branch just registered at the relay, preserving
        # the complete coherent eight-column strain and head state there.
        live[relay_rows, t] = at_relay
        misses.append(float(np.vdot(live[:, : t + 1], live[:, : t + 1]).real))
    return np.array(hits), np.array(returns), np.array(misses)


def test_retained_five_station_first_hit_histories_and_unsuccessful_trials():
    board = radar.source.formation.federation(0)
    stations = (*radar.source.formation.region(board), 5)
    state = np.array([1, 2j, -1.25, 1j, 0.4j, 0.5, -0.7j, 0.6])
    state /= np.linalg.norm(state)
    outward, returned = radar.directed_events(0, 1, stations=stations)
    # Mixed-sign phases probe incoming strain/head coherences and the endpoint
    # with index four exercises the independently enrolled extra station.
    for p, q in ((4, 1), (1, 4)):
        hits, back, miss = monitored_record(board, 1, stations[p], stations[q], state)
        np.testing.assert_allclose(
            np.einsum("i,tij,j->t", state.conj(), outward[p, q], state),
            hits,
            atol=3e-12,
            rtol=0,
        )
        np.testing.assert_allclose(
            np.einsum("i,tij,j->t", state.conj(), returned[p, q], state),
            back,
            atol=3e-12,
            rtol=0,
        )
        np.testing.assert_allclose(miss + back.cumsum(), 1, atol=8e-12, rtol=0)
        assert min(1 - hits.sum(), hits.sum() - back.sum()) >= -1e-12
    for cap in (8, 10, 12):
        ranges = radar.range_effects(0, 1, stations=stations, deadline=cap)
        assert ranges.shape == (10, 8, 8)
        np.testing.assert_allclose(
            ranges, ranges.conj().swapaxes(-1, -2), atol=1e-14, rtol=0
        )
        for row, (p, q) in zip(ranges, combinations(range(5), 2), strict=True):
            lens = []
            for i, j in ((p, q), (q, p)):
                bins = returned[i, j, : cap + 1]
                hits = outward[i, j, : cap + 1].sum(axis=0)
                missed = np.eye(8) - hits
                waiting = hits - bins.sum(axis=0)
                assert np.linalg.eigvalsh(missed).min() >= -1e-12
                assert np.linalg.eigvalsh(waiting).min() >= -1e-12
                lens.append(
                    (
                        np.einsum("t,tij->ij", np.arange(cap + 1), bins)
                        + cap * (missed + waiting)
                    )
                    / 4
                )
            np.testing.assert_allclose(row, (lens[0] + lens[1]) / 2, atol=2e-14, rtol=0)
        assert (
            np.linalg.eigvalsh(ranges).min() > 0
            and np.linalg.eigvalsh(ranges).max() <= cap / 4 + 2e-13
        )
    # At least one effect distinguishes joint states with both marginals fixed.
    a = ranges[0].reshape(2, 4, 2, 4)
    interaction = (
        ranges[0]
        - np.kron(np.eye(2) / 2, a.trace(axis1=0, axis2=2))
        - np.kron(a.trace(axis1=1, axis2=3), np.eye(4) / 4)
        + np.trace(ranges[0]) * np.eye(8) / 8
    )
    assert np.linalg.norm(interaction) > 1e-7


def test_zero_strain_coupling_reproduces_archived_two_head_radar():
    record = json.loads(Path(legacy.__file__).with_name("results.json").read_text())
    board = radar.source.formation.federation(1)
    route = radar.source.formation.region(board)
    states = dict(
        zip(
            ("formed", "independent"),
            radar.source.form(board, route, record["forming_load"]),
            strict=True,
        )
    )
    for charge in (-1, 1):
        outward, returned = radar.directed_events(1, charge, False)
        saved = next(
            row
            for row in record["runs"]
            if row["probe"] == charge and row["contrast"] == 8 and not row["left_odd"]
        )
        half = next(p for p in saved["panels"] if p["step"] == radar.STEP)
        for row in half["pairs"]:
            p, q = row["p"], row["q"]
            for tag, effect in (
                ("first_relay", outward[p, q]),
                ("first_return", returned[p, q]),
            ):
                # Free strain evolves under -X/2 throughout. Its final content
                # drops out of every event, including each signed head coherence.
                blocks = effect.reshape(-1, 2, 4, 2, 4)
                np.testing.assert_allclose(
                    blocks[:, 0, :, 0, :], blocks[:, 1, :, 1, :], atol=2e-13, rtol=0
                )
                np.testing.assert_allclose(blocks[:, 0, :, 1, :], 0, atol=2e-13, rtol=0)
                for name, rho in states.items():
                    np.testing.assert_allclose(
                        np.einsum("tij,ji->t", blocks[:, 0, :, 0, :], rho).real,
                        row["arms"][name][tag],
                        atol=8e-12,
                        rtol=0,
                    )


def test_integer_centers_bound_the_live_radar_and_decoupled_strain():
    from fractions import Fraction

    stations = (0, 5, 1)
    for coupled in (False, True):
        for k in (8, 10, 12):
            centers, denom, radius = radar.certified_ranges(0, -1, coupled, stations, k)
            assert centers.shape == (3, 2, 8, 8)
            assert 0 < radius < Fraction(2, 10**16)
            np.testing.assert_array_equal(centers[:, 0], centers[:, 0].swapaxes(1, 2))
            np.testing.assert_array_equal(centers[:, 1], -centers[:, 1].swapaxes(1, 2))
            nominal = (
                centers[:, 0].astype(float) + 1j * centers[:, 1].astype(float)
            ) / denom
            np.testing.assert_allclose(
                nominal,
                radar.range_effects(0, -1, coupled, stations, k),
                rtol=0,
                atol=4e-13,
            )
            if not coupled:
                for re, im in centers:
                    # Every exact off-diagonal strain block vanishes, and the two
                    # strain-diagonal blocks equal the same physical head effect.
                    # These integer Frobenius comparisons give operator upper tests.
                    off = sum(
                        x * x for part in (re[:4, 4:], im[:4, 4:]) for x in part.flat
                    )
                    difference = sum(
                        x * x
                        for part in (re[:4, :4] - re[4:, 4:], im[:4, :4] - im[4:, 4:])
                        for x in part.flat
                    )
                    assert Fraction(off, denom**2) < radius**2
                    assert Fraction(difference, denom**2) < (2 * radius) ** 2
