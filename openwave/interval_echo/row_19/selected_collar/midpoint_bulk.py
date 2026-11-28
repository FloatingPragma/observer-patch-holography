"""Read the two surface words and their inferred midpoint on a calibrated angular screen."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json

import networkx as nx
import numpy as np
from scipy.spatial import cKDTree

from model.executable_model import cell, federation, carrier
from model.experiments.mobile_memory.calibrated import geometry
from . import resolved_shells as shells

HERE = shells.HERE
SEED = 28122811
RADII = np.array([0.5, 0.75, 15 / 16, 1.25, 1.5])
TOL = 2e-12


def base_chart():
    vertices, edges = carrier.icosahedron_from_coordinates()
    # Native blank exchange reads all one-hop edges. The first nonconstant
    # eigenspace determines the twelve reference dot products up to rotation.
    adj = cell().adjacency().toarray()
    vals, vec = np.linalg.eigh(5 * np.eye(12) - adj)
    triplet = vec[:, np.isclose(vals, vals[1], atol=1e-10)]
    assert triplet.shape == (12, 3)
    relabel = carrier.isomorphism(cell().seams, edges, 12)
    vertices = np.array([geometry._unit(x) for x in vertices])
    refs = vertices[[relabel[i] for i in range(12)]]
    np.testing.assert_allclose(
        4 * triplet @ triplet.T, refs @ refs.T, atol=1e-12, rtol=0
    )
    faces = [
        (a, b, c)
        for a, b, c in combinations(range(12), 3)
        if {(a, b), (a, c), (b, c)} <= set(edges)
    ]
    return vertices, faces


def sites(vertices, faces):
    corners = np.asarray(vertices)[faces]
    directions = corners.sum(axis=1)
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    a, b, c = corners.transpose(1, 0, 2)
    omega = 2 * np.arctan2(
        abs(np.einsum("ij,ij->i", a, np.cross(b, c))),
        1 + np.sum(a * b + b * c + c * a, axis=1),
    )
    assert abs(omega.sum() - 4 * np.pi) < 1e-10 and omega.min() > 0
    # A positive normalized combination of the vertices stays in the same
    # cap, so this also bounds every point of a triangle, not just its corners.
    mesh = np.linalg.norm(corners - directions[:, None], axis=2).max() + TOL
    return directions, omega / (4 * np.pi), float(mesh)


def ruler(graph, points):
    graph_dist = np.array(
        nx.floyd_warshall_numpy(graph, nodelist=list(range(len(points))))
    )
    chord = np.linalg.norm(points[:, None] - points, axis=2)
    pairs = np.triu_indices(len(points), 1)
    ratios = 7 * graph_dist[pairs] / chord[pairs]
    a, b = float(min(ratios)), float(max(ratios))
    return {
        "swap_appointments_per_chord_range": [a, b],
        "minimum_relative_error_single_speed": (b - a) / (a + b),
    }


def moments(x, center, weights, radii):
    # Independently enumerate pairs on the finite screen. Remove the diagonal
    # before dividing by the observed distinct-address probability.
    z = np.sum((x - center) ** 2, axis=1)
    keep = z[:, None] <= radii**2 + TOL
    end_weight = weights @ keep
    endpoints = end_weight**2 - weights**2 @ keep
    mid = np.zeros_like(endpoints)
    raw = np.zeros_like(mid)
    # Stable set of three-dimensional dot products in blocks.
    for start in range(0, len(x), 64):
        stop = min(len(x), start + 64)
        dist2 = np.maximum(
            (z[start:stop, None] + z[None, :]) / 2 - (1 - x[start:stop] @ x.T) / 2, 0
        )
        for j, r in enumerate(radii):
            prob = np.maximum(1 - dist2 / (r * r), 0)
            inside = dist2 <= r * r + TOL
            for k, w in enumerate(weights):
                mid[k, j] += w[start:stop] @ (prob @ w)
                raw[k, j] += w[start:stop] @ (inside @ w)
    diagonal = weights**2
    norm = 1 - diagonal.sum(axis=1)[:, None]
    mid -= diagonal @ np.maximum(1 - z[:, None] / radii**2, 0)
    raw -= diagonal @ keep
    return np.stack([endpoints / norm, raw / norm, mid / norm], axis=-1).clip(0, 1)


def native_read(p, kappa, word, rng):
    values, counts = [], []
    for chance in word:
        tape, value = shells.event_counts(chance, p, rng)
        counts.append(tape)
        values.append(value)
    return {
        "counts": counts,
        "charge": shells.bounds(
            kappa * (values[1] - values[0]) / shells.interval(word[1] - word[0])
        ),
        "finite_admission": float(p),
    }


def direct_addresses(x, w, center, radii, rng):
    n = 2**17
    ids = rng.choice(len(x), (2, n), p=w)
    ids = ids[:, ids[0] != ids[1]]
    a, b = x[ids]
    m = (a + b) / 2
    far = np.maximum(
        np.linalg.norm(a - center, axis=1), np.linalg.norm(b - center, axis=1)
    )
    m2 = np.sum((m - center) ** 2, axis=1)
    coin = rng.random(len(m2))
    return len(m2), [
        [
            int(np.sum(far**2 <= r * r + TOL)),
            int(np.sum(m2 <= r * r + TOL)),
            int(np.sum(coin < np.maximum(1 - m2 / (r * r), 0))),
        ]
        for r in radii
    ]


def separation_counts(x, w, radii, rng):
    ids = rng.choice(len(x), (2, 2**18), p=w)
    a, b = x[ids[:, ids[0] != ids[1]]]
    separation = np.sum((a - b) ** 2, axis=1)
    return len(a), [int(np.sum(separation >= 4 * (1 - r * r - TOL))) for r in radii]


def continuum(r, central):
    if central:
        return 0.0, r * r, r * r / 2

    # The angular integral for a shell of occupied midpoint radius t is
    # proportional to the difference of the two squared positive parts.
    def primitive(u):
        return r**4 * u - 2 * r * r * u**3 / 3 + u**5 / 5

    integrated = 2 * primitive(min(1, r)) - primitive(r)

    # The unweighted midpoint probability integrates the first power.
    def first(u):
        return r * r * u - u**3 / 3

    raw = (2 * first(min(1, r)) - first(r)) / 2
    return r**4 / 16, raw, integrated / (4 * r * r)


def run():
    shells.iv.dps = 35
    reference = shells.dist.read("resolved_shells")
    kappa = shells.iv.mpf(reference["conditional_entropy_slope"])
    word = [
        F(reference["word_event_probabilities"][arm]) for arm in ("fresh", "shared")
    ]
    rng = np.random.default_rng(SEED)
    separation_rng = np.random.default_rng(SEED + 1)
    vertices, faces = base_chart()
    graphs, cal = geometry._atlas([federation(k) for k in range(3)])
    control = ruler(graphs[0], cal[0])
    levels = []
    for level in range(1, 5):
        vertices, faces = geometry._subdivide(vertices, faces)
        x, angular, mesh = sites(vertices, faces)
        if level <= 2:
            distances, ids = cKDTree(x).query(cal[level])
            assert distances.max() < 1e-12 and len(set(ids.tolist())) == len(x)
        weights = np.array([np.full(len(x), 1 / len(x)), angular])
        rows = []
        for site in [0, 4**level // 2, 4**level - 1, None]:
            center = np.zeros(3) if site is None else x[site]
            radii = np.array([0.25, 0.5, 0.75]) if site is None else RADII
            pops = moments(x, center, weights, radii)
            for k, mode in enumerate(("uniform_station", "solid_angle")):
                enrolled, tallies = direct_addresses(x, weights[k], center, radii, rng)
                if site is None:
                    attempts, hits = separation_counts(
                        x, weights[k], radii, separation_rng
                    )
                for j, r in enumerate(radii):
                    expected = continuum(r, site is None)
                    reads = [native_read(p, kappa, word, rng) for p in pops[k, j]]
                    end, _, mid = [shells.iv.mpf(v["charge"]) for v in reads]
                    extra = {}
                    if site is None:
                        f = hits[j] / attempts
                        err = shells.reg.src.probability_error(np.asarray(f), attempts)
                        extra["separation_read"] = {
                            "starts": attempts,
                            "hits": hits[j],
                            "enclosed_fraction": [f - err, f + err],
                        }
                    elif r > 1:
                        extra["additional_central_mass_fraction"] = shells.bounds(
                            (end - 15 * r * mid / 32)
                            / (kappa * (15 * r / 32) * (1 - 1 / r**2))
                        )
                    rows.append(
                        {
                            **extra,
                            "center": site,
                            "radius": float(r),
                            "address_law": mode,
                            "distinct_starts": enrolled,
                            "independent_address_hits": tallies[j],
                            "endpoint": reads[0],
                            "midpoint": reads[1],
                            "weighted_midpoint": reads[2],
                            "round_limit_admissions": list(expected),
                            "endpoint_minus_converted_midpoint": shells.bounds(
                                end - 15 * r * mid / 32
                            ),
                        }
                    )
        levels.append(
            {
                "level": level,
                "cells": len(x),
                "native_incidence_verified": level <= 2,
                "angular_mesh_upper": mesh,
                "angular_self_pair": float(angular @ angular),
                "uniform_to_area_total_variation": float(
                    abs(weights[0] - angular).sum() / 2
                ),
                "observations": rows,
            }
        )
    return {
        "seed": SEED,
        "starts_per_word_arm": shells.STARTS,
        "conditional_charge": shells.bounds(kappa),
        "native_range_control": control,
        "levels": levels,
        "inputs": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in ("MIDPOINT_BULK.md", "resolved_shells.json")
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    z = run()
    path = HERE / "midpoint_bulk.json"
    if args.write:
        path.write_text(json.dumps(z, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), z)
    print(z["native_range_control"])
    for v in z["levels"]:
        print(
            v["level"],
            v["angular_mesh_upper"],
            v["observations"][-1]["weighted_midpoint"]["charge"],
        )
