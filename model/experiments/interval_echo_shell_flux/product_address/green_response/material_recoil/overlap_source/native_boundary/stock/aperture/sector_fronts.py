"""Observer first-arrival times, counts and receiving ports on the occupied board."""

import argparse
import json
from fractions import Fraction
from hashlib import sha256
from itertools import product

import numpy as np
from scipy.linalg import qr
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from sympy import Matrix

from . import sector_fronts_lab as lab
from . import sector_passage_chart as chart
from .sector_passage_tide import project

HERE = chart.source.balls.HERE


def prompt(adj, distances, occupied):
    """Count empty-chart geodesics still on time at each station and incoming port."""
    allowed = adj * (~occupied[:, None] & ~occupied)
    adjacency = csr_matrix(allowed)
    paths = np.diag((~occupied).astype(np.int64))
    for distance in range(1, int(distances.max()) + 1):
        layer = paths * (distances == distance - 1)
        assert layer.max() * adj.sum(1).max() < 2**62
        following = layer @ adjacency
        mask = distances == distance
        paths[mask] = following[mask]
    u, v = np.nonzero(adj)
    by_port = paths[:, u] * allowed[u, v] * (distances[:, u] + 1 == distances[:, v])
    return paths, by_port


def calibrate(adj, homes, exterior, periods):
    bare = shortest_path(csr_matrix(adj), unweighted=True).astype(int)
    u, v = np.nonzero(adj)
    extpair = np.logical_and.outer(exterior, exterior)
    extport = extpair[:, v] & (np.arange(len(adj))[:, None] != v)
    all_paths = []
    all_ports = []
    same_time = same_support = 0
    L = csr_matrix(np.diag(adj.sum(1)) - adj)
    neighborhood = (adj @ exterior.astype(int)) > 0
    padded_control = np.where(adj, 2, 0)
    np.testing.assert_array_equal(shortest_path(csr_matrix(padded_control)), 2 * bare)
    for config in product((0, 1), repeat=4):
        occ = np.zeros(len(adj), dtype=bool)
        occ[homes] = config
        weight = adj * (1 + occ[:, None].astype(int) + occ)
        padding = padded_control - weight
        assert np.all(padding >= 0)
        np.testing.assert_array_equal(weight + padding, padded_control)
        if config in periods:
            np.testing.assert_array_equal(weight[np.triu(adj, 1) > 0], periods[config])
        T = shortest_path(csr_matrix(weight)).astype(int)
        assert np.all(T[:, v] <= T[:, u] + weight[u, v])
        support = T[:, u] + weight[u, v] == T[:, v]
        # Every nontrivial distance has a predecessor one positive step earlier.
        has_predecessor = np.logical_or.reduceat(
            support[:, np.argsort(v)], np.cumsum(np.r_[0, np.bincount(v)[:-1]]), axis=1
        )
        assert np.array_equal(has_predecessor, ~np.eye(len(adj), dtype=bool))
        n, incoming = prompt(adj, bare, occ)
        np.testing.assert_array_equal(
            T[np.ix_(neighborhood, exterior)], bare[np.ix_(neighborhood, exterior)]
        )
        for ruler in (lambda t: t, np.sqrt):
            np.testing.assert_array_equal(
                project(L, ruler(T) - ruler(bare))[extpair], 0
            )
        np.testing.assert_array_equal(T[extpair], bare[extpair])
        np.testing.assert_array_equal(support[extport], incoming[extport] > 0)
        if all_ports:
            np.testing.assert_array_equal(
                incoming[extport] > 0, all_ports[0][extport] > 0
            )
        assert np.all(n[extpair] > 0)
        same_time += 1
        same_support += 1
        all_paths.append(n)
        all_ports.append(incoming)
    # This graph has ordinary geodesics avoiding all four homes for every
    # exterior pair, and for every port tied at the front on the bare chart.
    # Depositing at home 20 or 32 never removes an exterior bare geodesic.
    for k in range(16):
        np.testing.assert_array_equal(all_paths[k][extpair], all_paths[k & 9][extpair])
        np.testing.assert_array_equal(all_ports[k][extport], all_ports[k & 9][extport])
    return (
        bare,
        np.array(all_paths)[[0, 8, 1, 9]],
        np.array(all_ports)[[0, 8, 1, 9]],
        {
            "configurations": same_time,
            "identical_front_port_sets": same_support,
            "native_configurations": len(periods),
            "exterior_stations": int(exterior.sum()),
            "front_acceleration_change_both_rulers": 0,
        },
    )


def coefficients(read):
    base, x, y, both = read
    return np.array([x - base, y - base, both - x - y + base])


def cert_image(cols, numerators, denominators, coords):
    """Ordinary calibration selects rows. Rational arithmetic certifies their rank."""
    active = np.linalg.norm(cols, axis=1)
    if np.any(active == 0):
        return {"rank_upper": int(sum(active != 0))}
    pivots = qr(cols / active[:, None], pivoting=True, mode="economic")[2][:3]
    rows = []
    intercepts = []
    for ix in pivots:
        b, x, y, both = [
            Fraction(int(n[ix]), int(d[ix]))
            for n, d in zip(numerators, denominators, strict=True)
        ]
        rows.append([x - b, y - b, both - x - y + b])
        intercepts.append(str(b))
    mat = Matrix(rows)
    return {
        "rows": coords[pivots].tolist(),
        "ordinary_mixed_matrix": [[str(v) for v in row] for row in rows],
        "ordinary_intercepts": intercepts,
        "determinant": str(mat.det()),
        "certified_rank": mat.rank(),
    }


def run():
    sourcefile = HERE / "sector_passage.json"
    tidefile = HERE / "sector_passage_tide.json"
    tide = json.loads(tidefile.read_bytes())["charts"]["commute"]["connected_witness"][
        "unit_covariance_coefficient"
    ]
    data = json.loads(sourcefile.read_bytes())
    src = data["spatial_modular_source"]
    _, _, ports, clearance, adj, _, loads, periods = chart.apparatus()
    homes = np.searchsorted(ports, [8, 20, 32, 44])
    ext = clearance >= 2
    d, N, P, cert = calibrate(adj, homes, ext, periods)
    u, v = np.nonzero(adj)
    i, j = np.nonzero(np.triu(np.logical_and.outer(ext, ext), 1))
    # With positive costs, only a route that visits both homes can have a
    # connected response. Its mass cannot exceed its share in the launch beam.
    mixed = N[3] - N[1] - N[2] + N[0]
    overlap = mixed[i, j] / N[0, i, j]
    k = int(np.argmax(overlap))
    s, t = i[k], j[k]
    share = Fraction(int(mixed[s, t]), int(N[0, s, t]))
    covariance = np.array(src["held_source"]["covariance_interval"])
    rate = src["rates"]["8"]
    native = np.array(rate["local_feature_contrasts"])[[0, 3, 4]]
    moment = float(rate["covariance_gain"])
    images = {}
    launch, e = np.nonzero(
        ext[:, None] & ext[v] & (np.arange(len(adj))[:, None] != v) & (P[0] > 0)
    )
    camera = ports[np.column_stack((launch, u[e], v[e]))]
    num = P[:, launch, e]
    fixed = np.broadcast_to(N[0, launch, v[e]], num.shape)
    fastest = N[:, launch, v[e]]
    camera_witnesses = []
    for key, denom in (
        ("all_arrivals", fixed),
        ("prompt_per_launch", fixed),
        ("fastest_each_preparation", fastest),
    ):
        observed = num if key != "all_arrivals" else np.broadcast_to(num[0], num.shape)
        C = (
            coefficients(observed / denom)
            if key == "fastest_each_preparation"
            else coefficients(observed) / denom[0]
        )
        at = int(np.argmax(abs(C[-1])))
        cert_rank = cert_image(C, observed, denom, camera)
        dl = C[-1] * moment
        full = native @ C
        images[key] = {
            "channels": len(dl),
            "active_connected": int(sum(abs(C[-1]) > 1e-12)),
            "minimum_connected_at_eighth": float(dl.min()),
            "maximum_connected_at_eighth": float(dl.max()),
            "full_native_range_at_eighth": [float(full.min()), float(full.max())],
            "calibration": cert_rank,
            "witness": camera[at].tolist(),
            "witness_mixed": float(C[-1, at]),
            "witness_matched": float(C[-1, at] / 2),
        }
        if key == "fastest_each_preparation":
            camera_witnesses.append((int(launch[at]), int(v[e[at]]), int(u[e[at]])))
            for start, port, end in cert_rank["rows"]:
                camera_witnesses.append(
                    tuple(np.searchsorted(ports, [start, end, port]))
                )
    x = float(share)
    F = coefficients(N[:, i, j] / N[0, i, j])
    peer = native @ F
    ordinary = -F[:2]
    # The 20/32 ordinary responses, and both geometric distance changes,
    # vanish identically across exterior pairs for the front-time readout.
    image_acq = lab.measure(
        adj,
        ports,
        homes,
        loads,
        src["reference"]["control_word_ids"],
        d,
        [(s, t, None), *camera_witnesses],
        chart.boundary.exchange_geometry.contact_swap(),
    )
    inversion = lab.invert(
        images["fastest_each_preparation"]["calibration"], image_acq[2:], tide
    )
    # Remove the exchange part of the 99-state generator and keep its
    # native agreement part, alphabet and port reversal symmetry. A blank
    # link cannot deliver either probe under this law at the same clock.
    emission = chart.boundary.exchange_geometry.lab.lab.cold.emission
    agreement = emission.four_h(1) - emission.four_h(0)
    assert np.count_nonzero(agreement[9:, :9]) == 0
    np.testing.assert_array_equal(
        agreement[:9, :9], np.diag(np.diag(agreement[:9, :9]))
    )
    np.testing.assert_allclose(
        image_acq[0]["records"]["1/fixed"]["forecast"][:3], x * np.array([1, -2, 1]) / 2
    )
    np.testing.assert_allclose(
        image_acq[1]["records"]["1/fastest"]["forecast"][3],
        images["fastest_each_preparation"]["witness_matched"],
    )
    return {
        "inputs_sha256": {
            f.name: sha256(f.read_bytes()).hexdigest()
            for f in (sourcefile, tidefile, HERE / "SECTOR_FRONTS.md")
        },
        "calibration": cert,
        "exterior_pairs": len(i),
        "mixed_sensitive_pairs": int(sum(mixed[i, j] > 0)),
        "marginal_sensitive_pairs_8_44": (ordinary > 0).sum(1).tolist(),
        "pair": ports[[s, t]].tolist(),
        "bare_time": int(d[s, t]),
        "both_path_fraction": str(share),
        "ordinary_route_counts_00_10_01_11": N[:, s, t].tolist(),
        "mixed_prompt_change_matched": x / 2,
        "mixed_prompt_at_eighth": x * moment,
        "mixed_prompt_from_source_counts": (x * covariance).tolist(),
        "minimum_prompt_change_per_positive_rate": x
        * src["all_rates"]["minimum_covariance_gain_per_rate"],
        "extra_time_variance_from_source_counts": (8 * x * covariance).tolist(),
        "whole_native_prompt_change_at_pair": float(peer[k]),
        "whole_native_prompt_range": [float(peer.min()), float(peer.max())],
        "delay_histogram_offsets": [0, 2, 4],
        "matched_change_in_offset_probabilities": (
            x * np.array([1, -2, 1]) / 2
        ).tolist(),
        "changed_clock_maximum_effect": 0,
        "agreement_only_changed_law_blank_link_delivery_probability": [0, 0],
        "images": images,
        "receiver_trials": image_acq,
        "exterior_image_inversion": inversion,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    out = run()
    path = HERE / "sector_fronts.json"
    if parser.parse_args().write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        chart.source.balls.radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_text()), out
        )
    print(
        out["pair"],
        out["both_path_fraction"],
        out["images"]["fastest_each_preparation"]["calibration"],
    )
