"""Check the local short-time sign, full native flight, raw reader and census denominators."""

from fractions import Fraction
from hashlib import sha256
import json

import numpy as np

from . import extended
from . import measure as m


def test_third_native_power_reaches_both_coins_before_any_remote_only_path():
    for level in range(3):
        b = m.source.federation(level)
        routes = list(m.routes(b, True))
        for r in (routes[0], routes[len(routes) // 2], routes[-1]):
            for q in (-1, 1):
                h = m.scaled_h(b, r, q) * 32
                i = m.initial_index(r, q)
                x = np.zeros(h.shape[0])
                x[4 * r[0] + i] = 1
                for _ in range(3):
                    assert not np.any(x.reshape(-1, 4)[:, [i ^ 1, i ^ 3]])
                    x = h @ x
                both = x.reshape(-1, 4)[:, i ^ 3]
                assert np.count_nonzero(both) == 1 and both[r[-1]] == -8
                assert not np.any(x.reshape(-1, 4)[:, i ^ 1])


def test_independent_native_exponential_and_reader_on_opposite_orientations():
    certificate_path = (
        m.HERE.parents[2]
        / "interval_echo_shell_flux/product_address"
        / "green_response/material_recoil/overlap_source/formation_certificate.json"
    )
    cert = {
        row["source_sign"]: row["covariance_billionths"]
        for row in json.loads(certificate_path.read_text())
    }
    for level in (0, 1):
        b = m.source.federation(level)
        r = m.source.region(b)
        for q in (-1, 1):
            i = m.initial_index(r, q)
            h = m.source.hamiltonian(b, r, q, 8)
            x = np.zeros(h.shape[0], complex)
            x[4 * r[0] + i] = 1
            flows = m.source.flow(h, x, m.TIMES)
            cheb = m.propagate(b, r, q)
            for t, state in enumerate(flows):
                rho = m.source.reduce(state)
                p = np.diag(rho).real
                np.testing.assert_allclose(p, cheb[t], rtol=0, atol=4e-13)
                perm = p[np.arange(4) ^ i]
                raw = m.covariance(p)
                np.testing.assert_allclose(
                    raw,
                    (1 if m.sign_class(i) == "alike" else -1) * m.covariance(perm),
                    atol=1e-13,
                )
                independent = np.kron(*m.source.marginals(rho))
                zz = np.trace(
                    np.kron(m.source.PAULI[2], m.source.PAULI[2]) @ (rho - independent)
                ).real
                np.testing.assert_allclose(4 * raw, zz, atol=1e-13)
            if level == 1:
                lo, hi = cert[q]
                assert lo < 1e9 * m.covariance(cheb[-1]) < hi


def test_placement_weights_census_controls_and_uniform_short_time_enclosures():
    p = m.HERE / "results.json"
    initial = json.loads(p.read_text())
    record = json.loads((m.HERE / "extended.json").read_text())
    assert record["first_archive_sha256"] == sha256(p.read_bytes()).hexdigest()
    assert (
        record["followup_sha256"]
        == sha256((m.HERE / "FOLLOWUP.md").read_bytes()).hexdigest()
    )
    assert (
        initial["protocol_sha256"]
        == sha256((m.HERE / "PROTOCOL.md").read_bytes()).hexdigest()
    )
    assert initial["error"]["covariance"] < 6e-8
    assert initial["controls"] == m.controls()
    for level in range(3):
        b = m.source.federation(level)
        sample = initial["boards"][str(b.cells)]
        rec = record["boards"][str(b.cells)]
        assert len(set(map(tuple, rec["minimum_routes"]))) == 3 * b.cells
        assert len(set(rec["standard_indices"])) == (3 * b.cells if level < 2 else 96)
        assert len(set(rec["survey_indices"])) == (1500 if level == 0 else 240)
        assert sum(sample["all_arm_placements"].values()) == 75 * b.cells
        bounds = rec["short_time"]
        assert bounds == extended.short_time_counts(b)
        assert Fraction(bounds["leading_flip_coefficient"]) == Fraction(16, 9)
        assert Fraction(bounds["all_arms_mean_over_t6"][1]) < 0
        for q in (-1, 1):
            a = rec["charges"][str(q)]
            assert len(a["census"]) == 3 * b.cells
            np.testing.assert_allclose(np.sum(a["standards"], axis=-1), 1, atol=3e-12)
            summary = a["summary"]
            for row in summary["original"]:
                assert (
                    row["positive"] + row["negative"] + row["unresolved"] == 3 * b.cells
                )
            assert (
                summary["original"][-1]["mean_radius_99"]
                == initial["error"]["covariance"]
            )
            for k in range(4):
                for t, row in enumerate(summary["standards"][str(k)]):
                    raw = np.asarray(a["standards"])[:, k, t, :]
                    np.testing.assert_allclose(
                        np.mean(m.covariance(raw)), row["mean"], atol=2e-14
                    )
            np.testing.assert_allclose(
                np.mean(m.covariance(np.asarray(a["standards"]).mean(axis=1)), axis=0),
                [row["mean"] for row in summary["independent_fair"]],
                atol=3e-14,
            )


def test_common_head_censuses_preserve_the_original_opportunities_and_sharp_inputs():
    source_path = m.HERE / "extended.json"
    held = json.loads(source_path.read_text())["boards"]
    final = json.loads((m.HERE / "common_heads.json").read_text())
    assert final["extended_sha256"] == sha256(source_path.read_bytes()).hexdigest()
    assert (
        final["protocol_sha256"]
        == sha256((m.HERE / "COMMON_HEADS.md").read_bytes()).hexdigest()
    )
    for n, width in (("20", 1500), ("320", 960)):
        original = held[n]
        roster = original["survey_routes"] if n == "20" else original["minimum_routes"]
        for sign in (-1, 1):
            old = original["charges"][str(sign)]
            src = old["survey_p"] if n == "20" else old["census"]
            known = {
                tuple(original["minimum_routes"][i]): row
                for i, row in zip(
                    original["standard_indices"], old["standards"], strict=True
                )
            }
            data = final["boards"][n]["charges"][str(sign)]
            extra = data["new_nontransmitting_inputs"]
            assert width == len(known) + len(extra) == len(roster)
            all_states = []
            for i, route in enumerate(roster):
                initial = m.initial_index(route, sign)
                if tuple(route) in known:
                    rows = known[tuple(route)]
                    assert str(i) not in extra
                else:
                    outside = iter(extra[str(i)])
                    rows = [src[i] if j == initial else next(outside) for j in range(4)]
                    assert next(outside, None) is None
                np.testing.assert_allclose(rows[initial], src[i], atol=2e-12)
                all_states.append(rows)
            all_states = np.asarray(all_states)
            for q in ("0", "1", "2", "3", "independent_fair"):
                density = (
                    all_states.mean(axis=1)
                    if q == "independent_fair"
                    else all_states[:, int(q)]
                )
                cov = m.covariance(density)
                calc = data["summary"][q]
                np.testing.assert_allclose(
                    cov.mean(axis=0), [row["mean"] for row in calc], atol=1e-14
                )
                assert all(row["n"] == width for row in calc)
