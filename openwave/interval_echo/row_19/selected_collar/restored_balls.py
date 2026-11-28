"""Measure fixed-volume area and regional position entropy on the ordinary-calibrated source."""

from fractions import Fraction as F
from math import log, sqrt

from mpmath import iv
import networkx as nx
import numpy as np

from . import stationary_polarization as steady

exp = steady.exp
local = exp.local
N = 2**58
SEED = 28122915


def bounds(z):
    return [
        float(np.nextafter(float(z.a), -np.inf)),
        float(np.nextafter(float(z.b), np.inf)),
    ]


def proper_clock(rng):
    # Five active contacts and one sharp baseline mark at the gateway. Read
    # native exchange parity and the load contrast energy on sharp stations.
    law = local.survey.routes.volume.candidate.native_law
    values, vecs = np.linalg.eigh(law.seam_term(8))
    rows = []
    for neighbor in (0, 1):
        energy = []
        counts = []
        for q in (1, 2):
            reading = iv.mpf(0)
            sample = []
            for other in [neighbor, 0, 0, 0, 0]:
                v = law.SRC.index((q, other))
                prob = vecs[v] ** 2
                spec = np.rint(values).astype(int)
                dist = {
                    int(e): float(prob[spec == e].sum())
                    for e in set(spec)
                    if prob[spec == e].sum() > 1e-13
                }
                keys = sorted(dist)
                read = local.acquisition([dist[k] for k in keys], keys, rng)
                reading += iv.mpf(read["mean"])
                sample.append({"energy": keys, "counts": read["counts"]})
            energy.append(reading)
            counts.append(sample)
        rows.append({"counts": counts, "frequency": bounds(energy[1] - energy[0])})
    ratio = iv.mpf(rows[1]["frequency"]) / iv.mpf(rows[0]["frequency"])
    return rows, ratio


def analyze(p0, p1, ball1, ball2, neighbor, V, A, rate):
    diff = [y - x for x, y in zip(p0, p1, strict=True)]
    pk = sum(p0[i] for i in neighbor)
    dp = sum(diff[i] for i in neighbor)
    c2 = 1 + (rate**2 - 1) * pk
    c3 = 1 + (rate**3 - 1) * pk
    dc2 = (rate**2 - 1) * dp
    dc3 = (rate**3 - 1) * dp
    coin = -(V[0] + V[1]) * dc3 / (2 * (V[1] - V[0]) * c3)
    da = (A[0] + A[1]) * dc2 / 2 + coin * (A[1] - A[0]) * c2
    ent = []
    ds = []
    for ball in (ball1, ball2):
        other = 1 - sum(p0[i] for i in ball)
        assert other.a > 0
        ent.append(-other * iv.log(other) - sum(p0[i] * iv.log(p0[i]) for i in ball))
        ds.append(sum(diff[i] * iv.log(other / p0[i]) for i in ball))
    change = (ds[0] + ds[1]) / 2 + coin * (ent[1] - ent[0])
    return da, change, coin


def labeled(da, ds, coin):
    gain = -ds / da if 0 not in da else None
    return {
        "area_tangent": bounds(da),
        "entropy_tangent": bounds(ds),
        "public_probability_tangent": bounds(coin),
        "stationary_gain": bounds(gain) if gain is not None else None,
    }


def run():
    iv.dps = 45
    rng = np.random.default_rng(SEED)
    setup, _ = steady.symmetric_instrument(local.survey.bulk.federation(0))
    _, cells, ports, ab, sites, degree, edges, _ = setup
    graph = nx.Graph(edges)
    assert all(degree[g] == 6 and graph.degree(g) == 5 for g in ab)
    clocks, rate = proper_clock(rng)
    row = local.inputs("stationary_polarization")["mirror_matched"]
    err = float(F(row["combined_trace_error_upper"]))
    source = row["outcomes"]["unread"]["arms"]
    surveys = []
    acquisitions = {}
    for direction in ("forward", "reversed"):
        obs = {}
        for name in ("fresh", "shared"):
            den = 2**100
            p = (
                np.array(
                    source[direction]["reads"][name][
                        "stationary_position_floor_2pow100"
                    ],
                    float,
                )
                / den
            )
            # Stable binning at these counts, with the certified finite wait
            # correction carried through every single probability interval.
            p = np.array(p, float)
            count = rng.multinomial(2**50, p / p.sum(), size=N // 2**50).sum(axis=0)
            half = sqrt(log(2 * 1000 / 0.01) / (2 * N)) + err + 1e-14
            band = [iv.mpf([float(z / N - half), float(z / N + half)]) for z in count]
            obs[name] = band
            acquisitions[direction + "_" + name] = count.tolist()
        ordinary = [iv.mpf(int(j == sites[0])) for j in range(len(ports))]
        for c in ab:
            distances = nx.single_source_shortest_path_length(graph, c)
            balls = [
                {i for i, d in distances.items() if d <= r}
                for r in range(max(distances.values()))
            ]
            for r in range(len(balls) - 1):
                if r + 1 >= float(rate.a) * (r + 2):
                    continue
                inner, outer = balls[r : r + 2]
                V = [len(b) for b in (inner, outer)]
                A = [nx.cut_size(graph, b) for b in (inner, outer)]
                params = (inner, outer, list(graph[c]), V, A, rate)
                a, s, b = analyze(obs["fresh"], obs["shared"], *params)
                oa, os, ob = analyze(obs["fresh"], ordinary, *params)
                # The determinant vanishes for two interventions that obey a
                # common coefficient, even if that coefficient is unrestricted.
                determinant = a * os - s * oa
                surveys.append(
                    {
                        "direction": direction,
                        "center": ports[c],
                        "hop_radii": [r, r + 1],
                        "ball_sites": [sorted(inner), sorted(outer)],
                        "volumes": V,
                        "cuts": A,
                        "record_source": labeled(a, s, b),
                        "ordinary_standard": labeled(oa, os, ob),
                        "common_gain_determinant": bounds(determinant),
                    }
                )
    return {
        "seed": SEED,
        "starts": N,
        "ordinary_center_partner": ports[sites[0]],
        "ports": ports,
        "cells": list(cells),
        "clocks": clocks,
        "proper_clock_ratio": bounds(rate),
        "position_counts": acquisitions,
        "surveys": surveys,
        "inputs": exp.save_inputs(
            "RESTORED_BALLS.md", "stationary_polarization.json", "contact_census.json"
        ),
    }


if __name__ == "__main__":
    result = run()
    exp.emit("restored_balls", result)
    for a in result["surveys"]:
        print(
            a["direction"],
            a["center"],
            a["hop_radii"],
            a["record_source"]["stationary_gain"],
            a["common_gain_determinant"],
        )
