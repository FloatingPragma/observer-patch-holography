"""Predict the relative acceleration of two freely released ordinary packets from a local joint reading."""

import argparse
import json
import math
from fractions import Fraction as Q
from hashlib import sha256

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import laplacian, shortest_path
from scipy.sparse.linalg import expm_multiply

from . import sector_passage_chart as chart
from . import sector_passage_inversion as inversion
from . import sector_passage_release as release

HERE = chart.source.balls.HERE


def project(L, field):
    # Each sharp blank-to-material contact has off-diagonal amplitude -1/2
    # for either load. Each of two independent packets contributes one half
    # of the neighbor-minus-launch sum to the initial distance acceleration.
    v = L @ field
    return -(v + v.T) / 2


def exterior_resistance_identity(adj, L, changes, i, j, err):
    # Differences of the two network Green functions are harmonic in either
    # coordinate at clean launches. Thus the acceleration change is half the
    # sum of changed incident edge resistances at each launch.
    maximum = 0.0
    for field in changes:
        incident = np.sum(adj * field, axis=1) / 2
        residual = project(L, field)[i, j] - incident[i] - incident[j]
        maximum = max(maximum, float(max(abs(residual))))
    assert maximum <= 24 * err + 1e-11
    return maximum


def response(R, dR, root):
    if not root:
        return dR
    den = np.sqrt(np.maximum(0, R + dR)) + np.sqrt(np.maximum(0, R))
    np.fill_diagonal(den, 1)
    return dR / den


def finite_flight(board, physical, ends, field, theta, raw):
    # The two releases are independent on repeatable source preparations.
    # The meter compares their terminal separation on the two observed charts.
    local = np.searchsorted(physical, ends)
    before = field[tuple(local)]
    result = {}
    n = 2**28
    span = np.ptp(field)
    for load in (1, 2):
        J = release.doubled_generator(board, physical, raw, load)
        rng = np.random.default_rng(248406 + load)
        eye = np.eye(len(J))[:, local]
        cols = expm_multiply(-0.5j * theta * csr_matrix(J), eye)
        probabilities = abs(cols) ** 2
        joint = np.outer(*probabilities.T)
        np.testing.assert_allclose(joint.sum(), 1, atol=5e-13)
        counts = rng.multinomial(n, joint.ravel() / joint.sum()).reshape(joint.shape)
        observed = float(np.sum(counts * (field - before)) / n)
        width = span * math.sqrt(math.log(4 / 0.005) / (2 * n)) + 2e-10
        result[str(load)] = {
            "trials": n,
            "mean_change": observed,
            "interval": [observed - width, observed + width],
            "forecast": float(np.sum(joint * (field - before))),
            "marginal_counts": [counts.sum(axis=k).tolist() for k in (0, 1)],
        }
    return result


def run():
    receipt = json.loads((HERE / "sector_passage.json").read_text())
    w, board, ports, clear, A, uv, loads, periods = chart.apparatus()
    configs = list(periods)
    rs = []
    errs = []
    for times in periods.values():
        _, _, r, e, _ = chart.law(A, uv, times)
        rs.append(r)
        errs.append(e)
    rs = np.array(rs)
    _F, C, ce, _ = inversion.basis(ports, A, configs, rs, errs)
    group = np.array([configs.index(tuple(z > 0)) for z in loads])
    weights = np.bincount(group) / len(w)
    R = np.einsum("i,ijk->jk", weights, rs)
    arms = np.array(receipt["spatial_modular_source"]["reference"]["control_word_ids"])
    low = rs[group[arms[1]]].mean(0)
    # The matched preparation has identical single-home histograms. Compute
    # its entire distance change from the measured coincidence difference.
    dR = C[-1] / 2
    np.testing.assert_allclose(rs[group[arms[0]]].mean(0) - low, dR, atol=2e-12)
    L = csr_matrix(laplacian(A))
    dist = shortest_path(board.adjacency(), unweighted=True)
    i, j = np.triu_indices(len(ports), 1)
    keep = (clear[i] >= 2) & (clear[j] >= 2) & (dist[ports[i], ports[j]] >= 3)
    i, j = i[keep], j[keep]
    identity_residual = exterior_resistance_identity(A, L, rs - R, i, j, 2 * max(errs))
    assert np.all(R[~np.eye(len(R), dtype=bool)] >= 1 / 7)
    active = np.flatnonzero(clear > 0)
    features = receipt["spatial_modular_source"]["rates"]["8"][
        "local_feature_contrasts"
    ]
    f = np.asarray(features)
    x, y = (loads[:, [0, 3]] > 0).astype(int).T
    first = Q(int(y[x == 0].sum()), int((x == 0).sum()))
    last = Q(int(x[y == 0].sum()), int((y == 0).sum()))
    residual = f[-1] - float(first) * f[0] - float(last) * f[3]
    peer = receipt["spatial_modular_source"]
    covariance_band = peer["held_source"]["covariance_interval"]
    min_covariance_rate = peer["all_rates"]["minimum_covariance_gain_per_rate"]
    peer_dR = np.einsum("k,kij->ij", f, C[1:])
    peer_field_error = (
        ce * sum(abs(f))
        + 5 * abs(C[1:]).max() * receipt["rates"]["8"]["source_error_l1"]
        + 1e-12
    )
    raw = chart.source.thermal.pair.stock.kernels()[2]
    result = {}
    for name, root in (("commute", False), ("square_root", True)):
        field = response(low, dR, root)
        acc = project(L, field)[i, j]
        e = 12 * (ce + 8 * abs(dR).max() * max(errs)) + 1e-10
        assert np.max(abs(acc)) > e
        # At a common reference the variation of sqrt(R) multiplies every
        # coefficient by the same positive local derivative on each chord.
        denom = 2 * np.sqrt(np.maximum(R, 0)) if root else np.ones_like(R)
        np.fill_diagonal(denom, 1)
        tangent = np.column_stack([project(L, v / denom)[i, j] for v in C[1:]])
        picks = inversion.baselines(tangent)
        M = tangent[picks]
        V = np.linalg.inv(M)
        entry_error = 12 * (2 * ce + 10 * abs(C[1:]).max() * max(errs)) + 1e-10
        slack = np.max(abs(np.eye(5) - V @ M).sum(1) + abs(V).sum(1) * 5 * entry_error)
        assert slack < 1
        # Compare the exact finite controlled excursion with an actual fill
        # of the ordinary empty home at 20. This home lies in a separate block.
        fillbase = rs[group[loads[:, 1] == 0]].mean(0)
        ordacc = project(L, response(fillbase, C[2], root))[i, j]
        sel = int(np.argmin(acc))
        ends = ports[[i[sel], j[sel]]]
        uncert = (2 * ce + 8 * abs(dR).max() * max(errs)) + 1e-10
        assert uncert < 1e-7
        flights = finite_flight(
            board, ports[active], ends, field[np.ix_(active, active)], 0.25, raw
        )
        for r in flights.values():
            r["interval"][0] -= 2 * uncert
            r["interval"][1] += 2 * uncert
        result[name] = {
            "decreasing": int(sum(acc < -e)),
            "increasing": int(sum(acc > e)),
            "unresolved": int(sum(abs(acc) <= e)),
            "lower": float(acc.min()),
            "upper": float(acc.max()),
            "acceleration_error": e,
            "witness_ports": ends.tolist(),
            "witness_ordinary_20": float(ordacc[sel]),
            "ordinary_and_source_both_inward": int(sum((acc < -e) & (ordacc < -2 * e))),
            "tangent_baselines": ports[np.column_stack((i[picks], j[picks]))].tolist(),
            "tangent_matrix": M.tolist(),
            "invertibility_residual": slack,
            "material_at_quarter_tick": flights,
        }
        peer_acc = project(L, response(R, peer_dR, root))[i, j]
        # Off-diagonal R is at least 1/7 in every arm. The two square roots
        # each have Lipschitz constant below 1.33 on this domain.
        peer_error = 12 * (2 * peer_field_error + 4 * max(errs)) + 1e-10
        result[name]["whole_peer_at_one_eighth"] = {
            "witness": float(peer_acc[sel]),
            "error": peer_error,
            "more_inward": int(sum(peer_acc < -peer_error)),
            "more_outward": int(sum(peer_acc > peer_error)),
        }
        if not root:
            coefficient = 2 * acc
            # Linearity of this ruler makes the isolated covariance term
            # additive even when the ordinary marginal occupancies change.
            coef_range = coefficient[sel] + np.array([-2 * e, 2 * e])
            assert coef_range[1] < 0
            result[name]["connected_witness"] = {
                "unit_covariance_coefficient": coef_range.tolist(),
                "upper_acceleration_per_positive_rate": float(
                    coef_range[1] * min_covariance_rate
                ),
                "local_measurement_acceleration_interval": [
                    float(coef_range[0] * covariance_band[1]),
                    float(coef_range[1] * covariance_band[0]),
                ],
                "forecast_at_one_eighth": float(
                    coefficient[sel] * peer["rates"]["8"]["covariance_gain"]
                ),
            }
            donor_gains = (float(first), 0, 0, float(last))
            comparisons = {}
            for k, home in enumerate((8, 20, 32, 44)):
                donor = tangent[:, k] + donor_gains[k] * tangent[:, -1]
                donor_error = (1 + donor_gains[k]) * entry_error
                # A negative contrast means more inward acceleration than on
                # the corresponding preparation before the intervention.
                safety = np.minimum.reduce(
                    (-acc - e, -peer_acc - peer_error, -donor - donor_error)
                )
                best = int(np.argmax(safety))
                comparisons[str(home)] = {
                    "simultaneously_more_inward": int(sum(safety > 0)),
                    "minimum_donor_acceleration": float(donor.min()),
                    "resolved_inward_donor_pairs": int(sum(donor < -donor_error)),
                    "donor_error": donor_error,
                    "witness": None
                    if safety[best] <= 0
                    else {
                        "ports": ports[[i[best], j[best]]].tolist(),
                        "controlled_change": float(acc[best]),
                        "whole_peer_change": float(peer_acc[best]),
                        "donor_change": float(donor[best]),
                    },
                }
            result[name]["ordinary_and_whole_peer"] = comparisons
            assert all(
                v["resolved_inward_donor_pairs"] == 0 for v in comparisons.values()
            )
    return {
        "source_sha256": sha256(
            (HERE / "sector_passage.json").read_bytes()
        ).hexdigest(),
        "chords": len(i),
        "maximum_incident_resistance_identity_residual": identity_residual,
        "ordinary_joint_gains": [str(first), str(last)],
        "native_tangent_unexplained_joint": residual,
        "native_tangent_error": 4 * receipt["rates"]["8"]["source_error_l1"] + 1e-11,
        "charts": result,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    out = run()
    dst = HERE / "sector_passage_tide.json"
    if opt.write:
        dst.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        chart.source.balls.radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(dst.read_text()), out
        )
    print(
        out["chords"],
        [
            (k, z["decreasing"], z["increasing"], z["witness_ports"], z["lower"])
            for k, z in out["charts"].items()
        ],
    )
