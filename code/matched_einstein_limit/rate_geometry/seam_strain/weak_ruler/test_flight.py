"""Live eight-station first- and second-click regression and independent checks."""

from fractions import Fraction
from hashlib import sha256
from math import factorial
import json

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

from . import flight

ARCHIVE = flight.HERE / "flight.json"
CONTACT_HASH = "20e8058b9e1ebd7ab7f5a3ee65c1cc16adc25558f47e017862390ca465db52fd"


def archived():
    data = json.loads(ARCHIVE.read_text())
    # Pin every exact contact, leaving platform-specific propagator and
    # quadrature decimals to the independent numerical enclosures below.
    contacts = [
        data["stations"],
        [
            (
                data["onsets"][q]["order"],
                data["onsets"][q]["integer_block_frobenius_squared"],
            )
            for q in ("-1", "1")
        ],
    ]
    assert (
        sha256(json.dumps(contacts, separators=(",", ":")).encode()).hexdigest()
        == CONTACT_HASH
    )
    return data


def test_directed_first_and_second_click_keep_the_full_live_state():
    charge, divisor = -1, 16
    h, rows = flight.hardware(charge)
    gen = (-1j / (40 * divisor)) * h
    p = -np.expm1(-1 / divisor)

    def direct(a, b):
        wave = np.zeros((h.shape[0], 8), complex)
        wave[rows[a]] = np.eye(8)
        hits = [np.zeros((8, 8), complex)]
        for _ in range(8):
            wave = expm_multiply(gen, wave, traceA=gen.diagonal().sum())
            hits.append(np.sqrt(p) * wave[rows[b]])
            wave[rows[b]] *= np.sqrt(1 - p)
        return np.stack(hits)

    a, b = 1, 2
    f, g = direct(a, b), direct(b, a)
    regen = flight.discrete_hits(charge, divisor)
    np.testing.assert_allclose(regen[:9, a, b], f, rtol=0, atol=4e-14)
    np.testing.assert_allclose(regen[:9, b, a], g, rtol=0, atol=4e-14)
    _, return_prob = flight.densities(f, g, flight.states(), divisor, False)
    for t in range(2, 9):
        operators = [g[t - j] @ f[j] for j in range(1, t)]
        E = sum(k.conj().T @ k for k in operators)
        expected = np.trace(E @ flight.states(), axis1=1, axis2=2).real
        np.testing.assert_allclose(return_prob[:, t], expected, rtol=0, atol=8e-17)
    # Independent orientation check for the non-Hermitian continuous generator.
    projector = np.zeros(h.shape[0])
    projector[rows[b]] = -0.5
    generator = (-1j / 40) * h + diags(projector)
    initial = np.zeros((h.shape[0], 8))
    initial[rows[a]] = np.eye(8)
    dt = 7 / 32
    state = expm_multiply(
        dt * generator, initial, traceA=dt * generator.diagonal().sum()
    )
    np.testing.assert_allclose(
        flight.receiver_continuum(charge, b)[7, a], state[rows[b]], atol=4e-14, rtol=0
    )

    # Continuous relay probability at the final deadline can be checked without
    # time quadrature by tracing the complete unclicked state over the board.
    last = flight.DEADLINE * generator
    survivor = expm_multiply(last, initial, traceA=last.diagonal().sum())
    gain = np.eye(8) - survivor.conj().T @ survivor
    exact = np.trace(gain @ flight.states(), axis1=1, axis2=2).real
    errors = {}
    for panel in archived()["continuous_panels"]:
        if panel["charge"] == charge:
            p = np.array([panel["states"][k]["relay"][a][b] for k in flight.KEYS])
            errors[panel["divisor"]] = np.max(abs(exact - p))
    assert errors[32] < 1.2e-6 and 3.9 * errors[32] < errors[16] < 4.1 * errors[32]


def test_registered_probabilities_conserve_every_timeout_and_match_previous_traces():
    data = archived()
    assert data["stations"] == list(flight.STATIONS)
    off = ~np.eye(8, dtype=bool)
    old = json.loads((flight.HERE.parent / "live_cadence.json").read_text())["panels"]
    for panel in data["panels"] + data["continuous_panels"]:
        for name in flight.KEYS:
            d = {k: np.array(v) for k, v in panel["states"][name].items()}
            np.testing.assert_allclose(
                d["first_hist"].sum(axis=-1), d["relay"], atol=3e-15, rtol=0
            )
            np.testing.assert_allclose(
                d["return_hist"].sum(axis=-1), d["complete"], atol=3e-15, rtol=0
            )
            expected = off.astype(int) if name == "reference" else 0
            np.testing.assert_allclose(
                d["complete"] + d["relay_only"] + d["no_relay"],
                expected,
                rtol=0,
                atol=3e-15,
            )
            if name == "reference":
                for k in ("first_hist", "return_hist", "relay_only", "no_relay"):
                    assert float(d[k].min()) >= -2e-15
                assert np.all((d["complete"][off] > 0) & (d["complete"][off] < 1))
                assert np.all(
                    (d["capped_range"][off] > 0) & (d["capped_range"][off] < 3)
                )
    for panel in data["panels"]:
        d, q = panel["divisor"], panel["charge"]
        for old_panel in old:
            if old_panel["detector"] != "rate_one" or old_panel["step"] != 1 / d:
                continue
            a, b = map(flight.STATIONS.index, old_panel["ports"])
            prior = old_panel["by_charge"][str(q)]
            for name in ("reference", "ordinary", "decorrelated", "hidden"):
                for key in ("complete", "capped_range"):
                    m = np.array(panel["states"][name][key])
                    expected = prior[name][key]
                    if name in ("ordinary", "decorrelated"):
                        expected -= prior["reference"][key]
                    assert abs((m[a, b] + m[b, a]) / 2 - expected) < 3e-13


def test_positive_continuum_onsets_and_independent_semigroup_replay():
    data = archived()
    for q in (-1, 1):
        result = flight.onset(q)
        saved = data["onsets"][str(q)]
        order = np.array(result["order"])
        assert np.array_equal(order, saved["order"]) and np.array_equal(order, order.T)
        assert np.bincount(order.ravel()).tolist() == [8, 22, 16, 18]
        assert (
            result["integer_block_frobenius_squared"]
            == saved["integer_block_frobenius_squared"]
        )
        h, _ = flight.hardware(q)
        assert Fraction(int(abs(h).sum(axis=1).max()), 40) + Fraction(1, 2) <= 24
        stored_floor = json.loads(
            (flight.HERE.parent / "live_modular.json").read_text()
        )["reference_minimum_eigenvalue_lower"]
        assert Fraction(str(stored_floor)) >= Fraction(14, 10000)
        uniform = min(
            flight.early_return_floor(
                int(order[i, j]), saved["integer_block_frobenius_squared"][i][j]
            )
            for i, j in zip(*np.nonzero(order), strict=True)
        )
        assert Fraction(28, 10**103) < uniform < Fraction(29, 10**103)
        # Signed source densities have nominal round-off in their zero trace.
        # Compare each expectation at a tolerance per unit operator norm, even
        # when the observable is scalar and its exact derivative is zero.
        first_norm = np.zeros((8, 8))
        return_norm = np.zeros((8, 8))
        for i, j in zip(*np.nonzero(order), strict=True):
            m = int(order[i, j])
            intensity = (
                saved["integer_block_frobenius_squared"][i][j]
                / (40**m * factorial(m)) ** 2
            )
            first_norm[i, j] = intensity / (2 * m + 1)
            return_norm[i, j] = (
                intensity**2 * factorial(2 * m) ** 2 / factorial(4 * m + 2)
            )
        for name in flight.KEYS:
            for key, bound in (
                ("first_cumulative", first_norm),
                ("return_cumulative", return_norm),
            ):
                delta = np.array(result["states"][name][key]) - np.array(
                    saved["states"][name][key]
                )
                assert np.all(abs(delta) <= 2e-13 * bound)
        for key in ("first_cumulative", "return_cumulative"):
            assert np.min(np.array(saved["states"]["reference"][key])[order != 0]) > 0
    for continuous, q in ((False, -1), (True, 1)):
        panel = flight.measure(16, q, continuous)
        (saved,) = [
            p
            for p in data["continuous_panels" if continuous else "panels"]
            if p["divisor"] == 16 and p["charge"] == q
        ]
        for name in flight.KEYS:
            for key in (*flight.BINS, *flight.SCALARS):
                np.testing.assert_allclose(
                    panel["states"][name][key],
                    saved["states"][name][key],
                    rtol=1e-10,
                    atol=5e-14,
                )
    for key in ("capped_range", "relay", "return_hist"):
        result = flight.pooled(data, 32, key, continuous=True)
        for matrix in result.values():
            np.testing.assert_equal(matrix, matrix.swapaxes(0, 1))
