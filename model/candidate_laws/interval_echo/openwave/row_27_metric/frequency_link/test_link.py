"""Independent gate reduction, finite-amplitude records and integer-bin checks."""

import gzip
import json
import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.material_inertia.dark_collar import label
from . import measure, native, readout


def archive():
    return json.loads(gzip.decompress(measure.ARCHIVE.read_bytes()))


def raw(k, q, arm, indices):
    """Act actual marked and blank contacts on every column, before centering."""
    law = candidate.native_law
    h0, h1 = law.local_term_parts()
    local = h0 + k * h1
    cols = np.eye(99)[:, [law.SRC.index(w) for w in ((q, 0), (0, q))]]
    n = len(indices)
    h = np.zeros((n, n))
    for i in range(n if n == 10 else 1):
        rows = [i, (i + 1) % n]
        h[np.ix_(rows, rows)] += cols.T @ local @ cols
    mag, coin = {"sham": (0, 1), "coin_even": (2, 1), "coin_odd": (2, -1)}[arm]
    for j, position in enumerate(indices):
        for target in (2, 8):
            w = (1, q if position == target else 0)
            source = np.kron(label(mag, 1, coin), np.eye(9)[law.SRC.index(w)])
            out = local @ source
            e = source @ out
            np.testing.assert_allclose(out, e * source, atol=1e-12)
            h[j, j] += e
    return h


def test_all_native_columns_against_direct_two_stage_propagation():
    blocks, inserts, g = native.contact_certificate()
    for k, arm in measure.CONFIGS:
        for q in (1, 2):
            h, hands = native.histories(k, q, arm, blocks, inserts)
            hraw = raw(k, q, arm, list(range(10)))
            for i, b in enumerate((2, 4)):
                direct = raw(k, q, arm, [b, b + 1])
                np.testing.assert_allclose(
                    direct - direct.trace() * np.eye(2) / 2, hands[i], atol=1e-12
                )
                np.testing.assert_allclose(
                    hraw - hraw.trace() * np.eye(10) / 10, h, atol=1e-12
                )
                for parked in (False, True):
                    t = 4.27
                    expected = np.zeros(11)
                    for delay in readout.D:
                        clock = expm(-1j * (t - delay) * direct)[:, 0]
                        expected[0] += abs(clock[0]) ** 2 / 2
                        hflight = np.diag(np.diag(hraw)) if parked else hraw
                        orbital = expm(-1j * delay * hflight)[:, b + 1]
                        expected[1:] += abs(clock[1] * orbital) ** 2 / 2
                    actual = native.output(
                        hands[i], h, b + 1, np.array([t]), readout.D, parked
                    )
                    np.testing.assert_allclose(actual[0], expected, atol=2e-12)
    fg = federation(0)
    raw_routes = native.preparation_routes(g)
    assert len({v[-1] for v in raw_routes}) == 4
    memories = [{site, g["parked_unit"]} for site in g["source_contact_sites"]]
    for route in raw_routes:
        assert route[-1] not in g["perimeter"] and route[-1] != g["parked_unit"]
        assert all(
            b in fg.neighbours[a] and {a, b} not in memories
            for a, b in zip(route, route[1:])
        )
    assert native.formed_source()["even_echo_material_word"] == [1, 1]


def test_complete_public_ports_and_certified_frequencies():
    r = archive()
    assert (
        r["appointment_step"]["interval"][0]
        < 1 / 16
        < r["appointment_step"]["interval"][1]
    )
    for case in r["banks"]:
        primary = (
            case["tag"] == 1 and case["contrast"] == 8 and case["source"] == "coin_even"
        )
        for link in case["links"]:
            h = np.array(link["observation"]["counts"])
            pmf = np.array(link["observation"]["native_probability"])
            assert h.shape == (1024, 11)
            assert np.all(h.sum(axis=1) == readout.N)
            np.testing.assert_allclose(pmf.sum(axis=1), 1, atol=1e-13)
            idx = 1 + measure.OUT[link["receiver"]]
            y = h[:, idx] / readout.N
            np.testing.assert_array_less(abs(y - pmf[:, idx]), readout.errors(y))
            assert (
                link["independent_forecast_held"][
                    "independent_delayed_clock_max_error_fraction"
                ]
                < 0.4
            )
            if link["parked"]:
                assert not h[:, idx].any() and link["ratio"] is None
            else:
                expected = np.sqrt(1025) if link["sender"] == 0 else 1 / np.sqrt(1025)
                nominal = expected if primary else 1
                lo, hi = link["ratio"]["allowed"]
                assert lo < nominal < hi
                assert hi / lo < 1.03
        for c in case["local"]:
            assert c["fit"]["constant_excluded"]
    with pytest.raises(ValueError):
        measure.compare({"counts": [[2, 1]]}, {"counts": [[1, 2]]})


def test_frequency_conversion_only_reads_histograms():
    for f in (0.85, 22.371, 44.273):
        p = 0.002 * (1 + 0.6 * np.cos(f * readout.T + 0.62))
        counts = np.rint(readout.N * p).astype(int)
        fitted = readout.read(counts)
        assert abs(fitted["frequency"] - f) < 1e-4
        assert any(a < f < b for a, b in fitted["frequency_intervals"])
    stopped = readout.read(np.zeros(1024, dtype=int))
    assert stopped["frequency"] is None
    assert not stopped["constant_excluded"] and stopped["frequency_intervals"] == [
        [0.4, 45]
    ]
