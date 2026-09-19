"""Measured order intervals and the paired record-metric family.

Production analysis only. verify.py independently authenticates and replays
the tapes and recomputes the counts with a different reachability algorithm.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import importlib.util
import json
from math import exp, lgamma, log
from pathlib import Path

import numpy as np

if __package__:
    from .experiment import ARCHIVE, HERE, ROOT, canonical, save_json, sha
else:
    from experiment import ARCHIVE, HERE, ROOT, canonical, save_json, sha


def dimension(f):
    if f == 1:
        return 1.0
    if not 0 < f < 1:
        return None
    def theory(d):
        return exp(lgamma(d + 1) + lgamma(d / 2) - log(2) - lgamma(1.5 * d))
    lo, hi = 1., 2.
    while theory(hi) > f:
        hi *= 2
    for _ in range(70):
        mid = (lo + hi) / 2
        if theory(mid) > f:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def interval(parents, bottom, top):
    """Ancestor intersection followed by exact topological bitsets.

    Bitsets are retired after their last child, bounding live storage. No
    Euclidean or layer-distance predicate enters the order calculation.
    """
    pending, past = [int(top)], set()
    while pending:
        v = pending.pop()
        if v < bottom or v in past:
            continue
        past.add(v)
        pending.extend(int(p) for p in parents[v] if p >= bottom)
    if bottom not in past:
        raise ValueError("incomparable interval tips")
    members = {int(bottom)}
    for v in sorted(past):
        if any(int(p) in members for p in parents[v]):
            members.add(v)
    ids = sorted(members)
    remaining = Counter(int(p) for v in ids for p in parents[v] if int(p) in members)
    live, heights, pairs, max_height = {}, {}, 0, 0
    for i, v in enumerate(ids):
        bits, height = 0, 1
        for p in set(map(int, parents[v])) & members:
            bits |= live[p]
            height = max(height, heights[p] + 1)
        pairs += bits.bit_count()
        live[v] = bits | (1 << i)
        heights[v] = height
        max_height = max(max_height, height)
        for p in parents[v]:
            p = int(p)
            if p in members:
                remaining[p] -= 1
                if remaining[p] == 0:
                    live.pop(p)
                    heights.pop(p)
    n = len(ids)
    fraction = Fraction(2 * pairs, n * (n - 1)) if n > 1 else None
    return np.asarray(ids, dtype=np.int32), {
        "bottom": int(bottom), "top": int(top), "events": n, "strict_pairs": pairs,
        "ordering_fraction": str(fraction) if fraction is not None else None,
        "ordering_fraction_float": float(fraction) if fraction is not None else None,
        "mm_dimension": dimension(float(fraction)) if fraction is not None else None,
        "height_in_events": max_height,
        "members_sha256": hashlib.sha256(np.asarray(ids, dtype="<i4").tobytes()).hexdigest()}


def spatial_summary(points):
    mean = points.mean(axis=0)
    delta = points - mean
    return {"centroid": mean.tolist(), "rms_radius": float(np.sqrt(np.mean(np.sum(delta**2, axis=1)))),
            "max_radius_from_centroid": float(np.linalg.norm(delta, axis=1).max()),
            "bounding_box_min": points.min(axis=0).tolist(), "bounding_box_max": points.max(axis=0).tolist()}


def trace_graph(directory):
    manifest = json.loads((directory / "trace.json").read_text())
    count = manifest["events"]
    parents = np.full((count, 2), -1, dtype=np.int32)
    owner = np.empty(count, dtype=np.int32)
    phase = np.empty(count, dtype=np.int16)
    positions, checkpoints = {}, {}
    for item in manifest["chunks"]:
        z = np.load(directory / item["file"])
        lo, end = item["first_event"], item["first_event"] + item["event_count"]
        if item["kind"] in ("intra", "glued"):
            parents[lo:end] = np.repeat(z["parents"], 2, axis=0)
            owner[lo:end] = z["endpoints"].reshape(-1) // 12
        else:
            parents[lo:end, 0] = z["parents"]
            owner[lo:end] = np.arange(item["event_count"]) // 12
        phase[lo:end] = item["phase"]
        positions[item["phase"]] = z["positions"]
        if item["kind"] == "glued":
            checkpoints[item["level"], item["sweep"]] = (item["phase"], z["port0_writers"])
    return manifest, parents, owner, phase, positions, checkpoints


def provenance(directory, output):
    manifest, parents, owner, phase, positions, checkpoints = trace_graph(directory)
    geometry = {level: dict(np.load(HERE / f"geometry/geometry_l{level}.npz")) for level in (3, 4, 5)}
    at, writers = checkpoints[3, 1]
    p = positions[at]
    anchors = {"writer_readback": int(np.argmin(np.sum((p-p.mean(axis=0))**2, axis=1))),
               "cell_centre": int(np.argmax(geometry[3]["centres"][:, 2]))}
    rows, member_arrays = [], {}
    for selection, anchor in anchors.items():
        child4 = int(geometry[4]["children"][anchor, 0])
        child5 = int(geometry[5]["children"][child4, 0])
        bottom = int(writers[anchor])
        previous = None
        for level, sweep, cell in ((3, 2, anchor), (3, 3, anchor), (4, 2, child4), (5, 2, child5)):
            end_phase, end_writers = checkpoints[level, sweep]
            top = int(end_writers[cell])
            ids, row = interval(parents, bottom, top)
            key = f"{selection}_l{level}_s{sweep}"
            member_arrays[key] = ids
            rb, cc, cells_by_level = [], [], {}
            for ph in np.unique(phase[ids]):
                selected = ids[phase[ids] == ph]
                item = manifest["chunks"][ph]
                lev = item["level"]
                rb.append(positions[ph][owner[selected]])
                cc.append(geometry[lev]["centres"][owner[selected]])
                cells_by_level.setdefault(lev, set()).update(map(int, owner[selected]))
            row.update({"selection": selection, "anchor_l3": anchor, "top_level": level, "top_sweep": sweep,
                        "member_array": key, "phase_gap": end_phase-int(phase[bottom]),
                        "interior": {"spatial": True, "reason": "closed S2 support; no spatial excision",
                                     "past_and_future_executed": int(phase[bottom]) > 0 and end_phase < len(manifest["chunks"])-1},
                        "placements": {"writer_readback": spatial_summary(np.concatenate(rb)),
                                       "cell_centre": spatial_summary(np.concatenate(cc))},
                        "cell_coverage": {str(l): len(s) for l, s in sorted(cells_by_level.items())},
                        "event_kinds": dict(Counter(manifest["chunks"][int(ph)]["kind"] for ph in phase[ids])),
                        "count_ratio_to_previous": None if previous is None else len(ids)/previous,
                        "placement_does_not_change_fixed_interval_order": True})
            previous = len(ids)
            rows.append(row)
            print(f"{key}: N={row['events']}, C={row['strict_pairs']}, d={row['mm_dimension']}", flush=True)
    np.savez_compressed(output / "interval_members.npz", **member_arrays)
    save_json(output / "provenance.json", rows)
    return rows


def record_metric_family(output, dim):
    source = ROOT / "evidence/source_net_causal_poset/build_causal_poset.py"
    spec = importlib.util.spec_from_file_location("paired_source_net", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    q, rounds = 13, 4
    values = module.orbit(q)
    a, b = module.axis_tables(values)
    sites = module.site_coordinates(q, dim)
    indptr, indices = module.site_graph(q, dim, values, a, b, sites)
    n = len(sites)
    payload = np.zeros((5, n), dtype=np.int64)
    payload[0] = np.arange(1, n+1)
    parents = [np.empty(0, dtype=np.int32) for _ in range(n)]
    for t in range(1, rounds+1):
        for site in range(n):
            reads = indices[indptr[site]:indptr[site+1]]
            payload[t, site] = 1 + sum(map(int, payload[t-1, reads]))
            parents.append((t-1)*n+reads)
    offsets = np.r_[0, np.cumsum([len(p) for p in parents])]
    flat = np.concatenate(parents)
    filename = "q13_reads.npz" if dim == 3 else f"q13_control_d{dim}_reads.npz"
    np.savez_compressed(output / filename, indptr=indptr, indices=indices,
                        parent_offsets=offsets, read_writers=flat, values=payload)
    centre = module.centre_axis_label(values)*sum(q**i for i in range(dim))
    historical = module.family_row(json.loads((source.parent / "source_net_causal_limit_receipt.json").read_text()), 13, dim)
    coordinates = np.asarray([float(x[0]+x[1]*(1+np.sqrt(5))/2) for x in values])[sites]
    rows = []
    previous = None
    for lag in range(1, 5):
        ids, row = interval(parents, centre, lag*n+centre)
        old = historical["vertical_intervals"][lag-1]
        assert row["events"] == old["inclusive_event_count"]
        assert row["strict_pairs"] == old["strict_pair_count"]
        # Exact containment: k^2/(4q) <= min(xi,1-xi)^2 in Q(phi).
        xi = values[module.centre_axis_label(values)]
        squares = [module.phi_square(xi), module.phi_square((1-xi[0], -xi[1]))]
        interior = all(module.phi_sign((4*q*x[0]-lag*lag, 4*q*x[1])) >= 0 for x in squares)
        assert interior == old["continuum_diamond_inside_cube"]
        row.update({"lag": lag, "interior": interior, "historical_counts_equal": True,
                    "counts_by_layer": np.bincount(ids//n, minlength=lag+1).tolist(),
                    "placement_in_units_of_L": spatial_summary(coordinates[ids % n]),
                    "count_ratio_to_previous": None if previous is None else len(ids)/previous})
        previous = len(ids)
        rows.append(row)
    result = {"q": 13, "dimension": dim, "sites": n, "rounds": 4,
              "reads": int(len(flat)), "events": len(parents), "centre": centre, "intervals": rows,
              "source_producer_sha256": sha(source), "trace_sha256": sha(output / filename),
              "law": "v(t+1,i)=1+sum(v(t,j) for metric neighbours j, including i)",
              "granularity": "one write per carrier per layer; distinct from two port writes per seam mean",
              "closed_interval_endpoint_convention": "included; r=1 is reported as d=1, unlike historical null"}
    return result


def q13(output):
    result = record_metric_family(output, 3)
    save_json(output / "q13.json", result)
    controls = {"schema": "oph.support_wiring.q13_controls.v1",
                "addendum_sha256": sha(HERE / "CONTROL_ADDENDUM.md"),
                "historical_receipt_sha256": sha(ROOT / "evidence/source_net_causal_poset/source_net_causal_limit_receipt.json"),
                "families": [record_metric_family(output, dim) for dim in (1, 2)],
                "flat_reference_ordering_fractions": {"1": "1/2", "2": "8/35"},
                "flat_reference_mm_dimensions": {"1": 2, "2": 3},
                "reference_is_acceptance_target": False}
    save_json(output / "q13_controls.json", controls)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=["provenance", "q13"])
    p.add_argument("--output", type=Path, default=ARCHIVE)
    a = p.parse_args()
    if a.mode == "provenance":
        provenance(a.output / "trace", a.output)
    else:
        q13(a.output)
