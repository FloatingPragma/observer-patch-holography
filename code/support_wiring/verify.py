"""Independent verification: imports neither experiment.py nor readouts.py.

Checks raw geometry, versioned execution and exact rational arithmetic before
using any parent edge. Interval pairs are counted by reverse descendant
bitsets, independently of the producer's forward ancestor bitsets. Kernels
are checked by n propagations and the full Gram Y^T Y, rather than 2n
propagations and a local return block.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.optimize import linear_sum_assignment

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT = ROOT / "evidence/support_wiring_776"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def encode(x):
    return (json.dumps(x, sort_keys=True, separators=(",", ":"), allow_nan=False)+"\n").encode()


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, "duplicate JSON key")
            result[key] = value
        return result
    def reject_constant(value):
        raise ValueError("nonfinite JSON constant: "+value)
    return json.loads(Path(path).read_text(), object_pairs_hook=unique, parse_constant=reject_constant)


def close(actual, expected, label, atol=2e-11):
    require(np.allclose(actual, expected, atol=atol, rtol=2e-11, equal_nan=False), label)


def integer_columns(a):
    require(a.ndim == 2 and a.shape[1] == 3 and a.dtype == np.uint64, "dyadic limb format")
    return np.asarray([int(low) | int(mid) << 64 | int(high) << 128 for low, mid, high in a], dtype=object)


def initial_values(n):
    return np.asarray([hashlib.sha256(b"oph776:20260919:" + str(i).encode()).digest()[0] & 7
                       for i in range(n)], dtype=object)


def schedule(g, wiring):
    # Explicit list of colours rather than the producer's endpoint-used sets.
    colors = []
    for edge in sorted(g["seams"].tolist()):
        for bucket in colors:
            if all(not (set(edge) & set(other)) for other in bucket):
                bucket.append(edge)
                break
        else:
            colors.append([edge])
    c = len(g["faces"])
    result = []
    for bucket in colors:
        e = np.asarray(bucket)
        edges = (12*np.arange(c)[:, None, None]+e).reshape(-1, 2)
        result.append(("intra", edges))
    if wiring != "isolated":
        x = g[wiring]
        result.append(("glued", np.asarray([12*x[:, 0]+x[:, 1], 12*x[:, 2]+x[:, 3]]).T))
    return result


def frame_from_graph(g):
    adjacency = np.zeros((12, 12))
    for a, b in g["seams"]:
        adjacency[a, b] = adjacency[b, a] = 1
    close(adjacency.sum(axis=1), 5, "carrier degree")
    lap = 5*np.eye(12)-adjacency
    # Spectral projector by its polynomial, not by the producer's eigensolver.
    lam = 5-math.sqrt(5)
    p = lap @ (lap-6*np.eye(12)) @ (lap-(5+math.sqrt(5))*np.eye(12))
    p /= lam*(lam-6)*(lam-5-math.sqrt(5))
    close(p@p, p, "rank-three projector")
    close(np.trace(p), 3, "projector rank")
    gram = 4*p
    import itertools
    for ix in itertools.combinations(range(12), 3):
        block = gram[np.ix_(ix, ix)]
        if np.linalg.det(block) > 1e-6:
            f = np.linalg.solve(np.linalg.cholesky(block), gram[list(ix)]).T
            break
    close(f@f.T, gram, "rank-three frame")
    return p, f


def verify_geometry():
    manifest = read_json(HERE / "geometry/geometry.json")
    require(manifest["revision"] == "7faa47b5cf00b42f6bf6b3e95ed4f7eb64f4239f", "source revision")
    old_path = ROOT / "code/source_routing/support_w12_l3.json"
    require(digest(old_path) == manifest["l3_original_sha256"], "original L3 pin")
    old = read_json(old_path)
    geometries, stats = {}, []
    for level in (3, 4, 5):
        path = HERE / f"geometry/geometry_l{level}.npz"
        require(digest(path) == manifest["levels"][str(level)]["sha256"], "geometry hash")
        g = dict(np.load(path))
        geometries[level] = g
        verts, faces = g["vertices"], g["faces"]
        c = 20*4**level
        require(faces.shape == (c, 3), "cell count")
        close(np.linalg.norm(verts, axis=1), 1, "unit vertices")
        centres = verts[faces].sum(axis=1)
        centres /= np.linalg.norm(centres, axis=1)[:, None]
        close(centres, g["centres"], "cell centres")
        close(g["areas"].sum(), 4*math.pi, "sphere area")
        va, vb, vc = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
        area = 2*np.arctan2(abs(np.sum(va*np.cross(vb, vc), axis=1)),
                           1+np.sum(va*vb+vb*vc+vc*va, axis=1))
        close(g["areas"], area, "individual spherical triangle areas")
        v_incidence, e_incidence = defaultdict(list), defaultdict(list)
        for cell, face in enumerate(faces):
            for j in range(3):
                v_incidence[int(face[j])].append(cell)
                e_incidence[tuple(sorted((int(face[j]), int(face[(j+1)%3]))))].append(cell)
        require(all(len(x) == 2 for x in e_incidence.values()), "closed sphere edges")
        require(Counter(map(len, v_incidence.values())) == {5: 12, 6: len(verts)-12}, "pentagonal vertices")
        if level == 3:
            require(np.array_equal(faces, old["faces"]) and np.array_equal(g["w12"], old["glued_pairs"]), "L3 reproduction")
        else:
            coarse = geometries[level-1]
            parents, children = g["parent"], g["children"]
            require(np.array_equal(np.sort(children.ravel()), np.arange(c)), "child partition")
            require(np.array_equal(parents[children], np.arange(c//4)[:, None].repeat(4, axis=1)), "parent transport")
            # Check the four ordered child triangles geometrically, not only labels.
            points = coarse["vertices"][coarse["faces"]]
            x, y, z = points[:, 0], points[:, 1], points[:, 2]
            xy, yz, zx = x+y, y+z, z+x
            for mid in (xy, yz, zx):
                mid /= np.linalg.norm(mid, axis=1)[:, None]
            expected = np.stack((np.stack((x, xy, zx), 1), np.stack((y, yz, xy), 1),
                                 np.stack((z, zx, yz), 1), np.stack((xy, yz, zx), 1)), 1)
            close(verts[faces[children]], expected, "geodesic child triangles")
            close(g["areas"][children].sum(axis=1), coarse["areas"], "area refinement")
            close(g["expectation_weights"], g["areas"]/coarse["areas"][parents], "conditional expectation")
            close(g["expectation_weights"][children].sum(axis=1), 1, "expectation left inverse")
        p, f = frame_from_graph(g)
        base = verts[:12]
        close(base@base.T, 4*p, "port labels agree with geometry template")
        antipode = np.argmin(base@base.T, axis=1)
        normals = centres
        reference = np.zeros_like(normals)
        reference[:, 2] = 1
        reference[abs(normals[:, 2]) > .9] = [1, 0, 0]
        tx = np.cross(reference, normals)
        tx /= np.linalg.norm(tx, axis=1)[:, None]
        ty = np.cross(normals, tx)
        frames = np.stack([tx, ty, normals], axis=2)
        for wiring, inc in (("w12", v_incidence), ("w3", e_incidence)):
            pairs = sorted({(a, b) for cells in inc.values() for a in cells for b in cells if a < b})
            wire = g[wiring]
            require(np.array_equal(wire[:, [0, 2]], pairs), "complete geometric neighbourhood")
            support = sparse.coo_matrix((np.ones(len(wire)), (wire[:, 0], wire[:, 2])), shape=(c, c)).tocsr()
            require(sparse.csgraph.connected_components(support, directed=False, return_labels=False) == 1, "connected glued support")
            slots = np.r_[wire[:, 0]*12+wire[:, 1], wire[:, 2]*12+wire[:, 3]]
            require(np.all((wire[:, [1, 3]] >= 0) & (wire[:, [1, 3]] < 12)), "port range")
            require(len(np.unique(slots)) == len(slots), "unique glued port slots")
            by_cell = [[] for _ in range(c)]
            for a, pa, b, pb in wire:
                by_cell[a].append((int(b), int(pa)))
                by_cell[b].append((int(a), int(pb)))
            # Independently check optimal objective at every cell; ties may be
            # resolved by the captured source's stable endpoint ordering.
            for cell, endpoints in enumerate(by_cell):
                other, port = np.asarray(endpoints).T
                direction = centres[other]-centres[cell]
                direction /= np.linalg.norm(direction, axis=1)[:, None]
                scores = ((direction@frames[cell])@base.T).astype(np.float32)
                r, col = linear_sum_assignment(-scores)
                observed = scores[np.arange(len(port)), port].sum(dtype=np.float64)
                optimal = scores[r, col].sum(dtype=np.float64)
                require(abs(observed-optimal) < 2e-6, "geometric assignment is not optimal")
            degree = np.bincount(wire[:, [0, 2]].ravel(), minlength=c)
            defects = np.asarray([any(len(v_incidence[int(v)]) == 5 for v in face) for face in faces])
            require(np.count_nonzero(defects) == 60, "sixty pentagonal cells")
            if wiring == "w12":
                require(np.all(degree == np.where(defects, 11, 12)), "pentagonal unused slots")
            else:
                require(np.all(degree == 3), "three-port degrees")
            va = np.einsum("nij,nj->ni", frames[wire[:, 0]], base[wire[:, 1]])
            vb = np.einsum("nij,nj->ni", frames[wire[:, 2]], base[wire[:, 3]])
            cosine = np.sum(va*vb, axis=1)
            stats.append({"level": level, "wiring": wiring, "carriers": c, "glued_seams": len(wire),
                          "degree_histogram": {str(k): v for k, v in sorted(Counter(map(int, degree)).items())},
                          "unused_ports": 12*c-2*len(wire), "defect_cells": int(defects.sum()),
                          "unused_on_defect_cells": int((12-degree[defects]).sum()),
                          "antipodal_label_pairs": int(np.count_nonzero(antipode[wire[:, 1]] == wire[:, 3])),
                          "global_direction_dot": {"min": float(cosine.min()), "mean": float(cosine.mean()), "max": float(cosine.max())}})
    return geometries, stats


def verify_trace(directory, geometries):
    manifest = read_json(directory / "trace.json")
    binding = {"schema": "oph.support_wiring.trace.v1", "specification_sha256": digest(HERE / "SPECIFICATION.md"),
               "geometry_sha256": digest(HERE / "geometry/geometry.json"), "sweeps_per_level": 4,
               "levels": [3, 4, 5], "numerator_bits": 192, "law": "exact_pair_mean",
               "events": "two simultaneous destination writes per mean action"}
    require(manifest["binding"] == binding, "trace contract binding")
    expected = []
    for level in (3, 4, 5):
        expected.append((level, 0, "prepare" if level == 3 else "copy", None))
        for sweep in range(1, 5):
            expected.extend((level, sweep, kind, edge) for kind, edge in schedule(geometries[level], "w12"))
    require(len(manifest["chunks"]) == len(expected), "complete scheduled phase count")
    total_events = sum(12*len(geometries[l]["faces"]) if edge is None else 2*len(edge)
                       for l, _, _, edge in expected)
    require(manifest["events"] == total_events, "complete write count")
    parents = np.full((total_events, 2), -1, dtype=np.int32)
    owner = np.empty(total_events, dtype=np.int32)
    phase = np.empty(total_events, dtype=np.int16)
    snapshots, checkpoints = {}, {}
    chain = hashlib.sha256(encode(binding)).hexdigest()
    current_id, exponent = 0, 0
    state = writers = None
    mean_actions, identity_actions = 0, 0
    for index, (item, (level, sweep, kind, edge)) in enumerate(zip(manifest["chunks"], expected)):
        g = geometries[level]
        path = directory / item["file"]
        require(item["file"] == f"phase_{index:03d}.npz", "chunk filename")
        require(item["sha256"] == digest(path), "trace chunk hash")
        require(item["previous"] == chain, "trace previous hash")
        rest = {k: v for k, v in item.items() if k != "chain"}
        chain = hashlib.sha256(encode(rest)).hexdigest()
        require(chain == item["chain"], "trace chain")
        require((item["level"], item["sweep"], item["kind"], item["phase"], item["first_event"])
                == (level, sweep, kind, index, current_id), "scheduled phase metadata")
        z = dict(np.load(path))
        values = integer_columns(z["values"])
        if kind == "prepare":
            state = initial_values(12*len(g["faces"]))
            require(np.array_equal(values, state), "preparation values")
            require(np.array_equal(z["parents"], np.full(len(state), -1)), "preparation has no reads")
            writers = np.arange(current_id, current_id+len(state), dtype=np.int32)
            owned = np.arange(len(state))//12
            count = len(state)
        elif kind == "copy":
            source = (g["parent"][:, None]*12+np.arange(12)).ravel()
            require(np.array_equal(z["parents"], writers[source]), "refinement current parent writers")
            require(np.array_equal(values, state[source]), "refinement componentwise copy")
            state = values
            count = len(state)
            parents[current_id:current_id+count, 0] = z["parents"]
            writers = np.arange(current_id, current_id+count, dtype=np.int32)
            owned = np.arange(count)//12
        else:
            require(np.array_equal(z["endpoints"], edge), "missing or changed seam attempts")
            a, b = edge.T
            require(np.array_equal(z["parents"], np.column_stack([writers[a], writers[b]])), "stale or unauthenticated read")
            # Equality in integers at the new common denominator, all attempts.
            require(np.array_equal(values, state[a]+state[b]), "noncanonical seam mean")
            identity_actions += int(np.count_nonzero(state[a] == state[b]))
            old_total, old_square = sum(state), sum(state*state)
            loss = 2*sum((state[a]-state[b])**2)
            state = state*2
            state[a] = state[b] = values
            exponent += 1
            require(sum(state) == 2*old_total, "exact total conservation")
            require(4*old_square-sum(state*state) == loss, "exact quadratic descent identity")
            count = 2*len(edge)
            new_ids = np.arange(current_id, current_id+count, dtype=np.int32).reshape(-1, 2)
            writers[a], writers[b] = new_ids[:, 0], new_ids[:, 1]
            parents[current_id:current_id+count] = np.repeat(z["parents"], 2, axis=0)
            owned = edge.ravel()//12
            mean_actions += len(edge)
        require(item["exponent"] == exponent and item["event_count"] == count, "dyadic exponent or event count")
        require(np.array_equal(z["port0_writers"], writers[::12]), "checkpoint writers")
        _, frame = frame_from_graph(g)
        positions = (np.asarray(state, dtype=float)/2.**exponent).reshape(-1, 12)@frame
        close(z["positions"], positions, "writer rank-three readback")
        require(np.all(parents[current_id:current_id+count] < current_id), "causal phase has no internal ordering")
        owner[current_id:current_id+count] = owned
        phase[current_id:current_id+count] = index
        snapshots[index] = positions
        if kind == "glued":
            checkpoints[level, sweep] = (index, writers[::12].copy())
        current_id += count
    require(chain == manifest["final_chain"] and current_id == total_events, "terminal custody")
    summary = {"events": total_events, "mean_actions": mean_actions, "identity_mean_actions": identity_actions,
               "final_exponent": exponent, "exact_total_and_descent_checks": True, "final_chain": chain}
    print(f"Authenticated {total_events:,} events; exact means, joins and readbacks pass", flush=True)
    return manifest, parents, owner, phase, snapshots, checkpoints, summary


def count_interval(parents, bottom, top):
    """Independent reverse sweep; count strict descendants, not ancestors."""
    visited, todo = set(), [int(top)]
    while todo:
        v = todo.pop()
        if v in visited or v < bottom:
            continue
        visited.add(v)
        todo.extend(int(p) for p in parents[v] if p >= bottom)
    children = defaultdict(list)
    for v in visited:
        for p in set(map(int, parents[v])):
            if p in visited:
                children[p].append(v)
    members, todo = set(), [int(bottom)]
    while todo:
        v = todo.pop()
        if v in members:
            continue
        members.add(v)
        todo.extend(children[v])
    require(top in members and bottom in visited, "interval tips must be comparable")
    ids = sorted(members)
    positions = {v: k for k, v in enumerate(ids)}
    descendants, lengths = {}, {}
    pairs = 0
    for v in reversed(ids):
        bits, height = 0, 1
        for child in children[v]:
            bits |= descendants[child] | (1 << positions[child])
            height = max(height, lengths[child]+1)
        descendants[v] = bits
        lengths[v] = height
        pairs += bits.bit_count()
    return np.asarray(ids, dtype=np.int32), pairs, lengths[bottom]


def check_order_row(row, ids, pairs, height):
    n = len(ids)
    require(row["events"] == n and row["strict_pairs"] == pairs, "interval counts")
    fraction = Fraction(2*pairs, n*(n-1)) if n > 1 else None
    require(row["ordering_fraction"] == (str(fraction) if fraction is not None else None), "ordering fraction")
    require(row["height_in_events"] == height, "interval height")
    require(row["members_sha256"] == hashlib.sha256(ids.astype("<i4").tobytes()).hexdigest(), "interval membership")
    if fraction is not None:
        close(row["ordering_fraction_float"], float(fraction), "fraction float")
    d = row["mm_dimension"]
    if fraction and 0 < fraction <= 1:
        require(d is not None and d >= 1, "MM inversion domain")
        theory = math.exp(math.lgamma(d+1)+math.lgamma(d/2)-math.log(2)-math.lgamma(1.5*d))
        close(theory, float(fraction), "MM inverse substitution")
    else:
        require(d is None, "undefined dimension")


def check_spatial(row, points):
    mean = np.sum(points, axis=0)/len(points)
    radii2 = np.sum((points-mean)**2, axis=1)
    close(row["centroid"], mean, "position centroid")
    close(row["rms_radius"], math.sqrt(float(radii2.sum()/len(points))), "position RMS")
    close(row["max_radius_from_centroid"], math.sqrt(float(max(radii2))), "position radius")
    close(row["bounding_box_min"], np.min(points, axis=0), "position min")
    close(row["bounding_box_max"], np.max(points, axis=0), "position max")


def verify_intervals(output, geometries, replay):
    manifest, parents, owner, phase, positions, checkpoints, _ = replay
    rows = read_json(output / "provenance.json")
    stored_ids = np.load(output / "interval_members.npz")
    ph, writer = checkpoints[3, 1]
    rb = positions[ph]
    anchors = {"writer_readback": int(np.argmin(np.linalg.norm(rb-rb.mean(axis=0), axis=1))),
               "cell_centre": int(np.argmin(np.linalg.norm(geometries[3]["centres"]-[0, 0, 1], axis=1)))}
    require(len(rows) == 8, "complete interval census")
    for j, selection in enumerate(("writer_readback", "cell_centre")):
        cell = anchors[selection]
        child = int(geometries[4]["children"][cell, 0])
        grandchild = int(geometries[5]["children"][child, 0])
        previous = None
        for k, (lev, sweep, dest) in enumerate(((3, 2, cell), (3, 3, cell), (4, 2, child), (5, 2, grandchild))):
            row = rows[4*j+k]
            end_phase, end_writer = checkpoints[lev, sweep]
            bottom, top = int(writer[cell]), int(end_writer[dest])
            require((row["selection"], row["anchor_l3"], row["top_level"], row["top_sweep"], row["bottom"], row["top"])
                    == (selection, cell, lev, sweep, bottom, top), "frozen tip selection")
            ids, pairs, height = count_interval(parents, bottom, top)
            check_order_row(row, ids, pairs, height)
            require(np.array_equal(ids, stored_ids[row["member_array"]]), "stored interval members")
            require(row["phase_gap"] == end_phase-int(phase[bottom]), "count growth phase clock")
            if previous is None:
                require(row["count_ratio_to_previous"] is None, "first count ratio")
            else:
                close(row["count_ratio_to_previous"], len(ids)/previous, "count growth ratio")
            previous = len(ids)
            rb_points, centre_points, cells_by_level = [], [], defaultdict(set)
            for event in ids:
                event_phase, cell_id = int(phase[event]), int(owner[event])
                level = manifest["chunks"][event_phase]["level"]
                cells_by_level[level].add(cell_id)
                rb_points.append(positions[event_phase][cell_id])
                centre_points.append(geometries[level]["centres"][cell_id])
            check_spatial(row["placements"]["writer_readback"], np.asarray(rb_points))
            check_spatial(row["placements"]["cell_centre"], np.asarray(centre_points))
            require(row["cell_coverage"] == {str(l): len(s) for l, s in cells_by_level.items()}, "cell coverage")
            require(row["event_kinds"] == dict(Counter(manifest["chunks"][int(ph)]["kind"] for ph in phase[ids])), "event kind counts")
            require(row["interior"]["spatial"] and row["interior"]["past_and_future_executed"]
                    and int(phase[bottom]) > 0 and end_phase < len(manifest["chunks"])-1, "interior tips")
            require(row["placement_does_not_change_fixed_interval_order"] is True, "placement/order separation")
            print(f"Verified interval {selection} L{lev} sweep {sweep}: {len(ids)} events", flush=True)
    return rows


def qphi_sign(a, b):
    """Exact sign of a+b*phi by comparing integers u^2 and 5v^2."""
    u, v = 2*np.asarray(a)+np.asarray(b), np.asarray(b)
    return np.where((u >= 0) & (v >= 0), np.sign(u+v),
                    np.where((u <= 0) & (v <= 0), np.sign(u+v), np.sign(u)*np.sign(u*u-5*v*v)))


def verify_metric_family(output, geometry, row, dim):
    path = output / ("q13_reads.npz" if dim == 3 else f"q13_control_d{dim}_reads.npz")
    require(row["trace_sha256"] == digest(path), "q13 trace hash")
    require(row["source_producer_sha256"] == digest(ROOT / "evidence/source_net_causal_poset/build_causal_poset.py"), "q13 definition pin")
    q, n = 13, 13**dim
    labels = np.indices((q,)*dim).reshape(dim, -1).T
    floors = np.asarray([(i+math.isqrt(5*i*i))//2 for i in range(q)])
    integer = -floors[labels]
    a = integer[:, None, :]-integer[None, :, :]
    b = labels[:, None, :]-labels[None, :, :]
    aa = np.sum(a*a+b*b, axis=2)
    bb = np.sum(2*a*b+b*b, axis=2)
    relation = qphi_sign(q*aa-1, q*bb) <= 0
    # Readback metric identity for the same source port records; full pair census.
    ma, mb, mc = np.pad(integer, ((0, 0), (0, 3-dim))).T
    ba, b_b, bc = np.pad(labels, ((0, 0), (0, 3-dim))).T
    record = np.column_stack([b_b-ma, b_b+ma, bc-mb, bc+mb, ba-mc, ba+mc])
    positive, negative = (0, 1, 4, 5, 8, 9), (3, 2, 7, 6, 11, 10)
    loads = np.zeros((n, 12))
    for col, (p, m) in enumerate(zip(positive, negative)):
        loads[:, p] = np.maximum(record[:, col], 0)
        loads[:, m] = np.maximum(-record[:, col], 0)
    _, frame = frame_from_graph(geometry)
    rb = loads@frame
    phi = (1+math.sqrt(5))/2
    scale2 = 4/(phi+2)
    for start in range(0, n, 64):
        distances = np.sum((rb[start:start+64, None, :]-rb[None, :, :])**2, axis=2)
        close(distances, scale2*(aa[start:start+64]+bb[start:start+64]*phi), "q13 readback metric identity", atol=2e-10)
    del a, b
    z = dict(np.load(path))
    ptr = np.r_[0, np.cumsum(relation.sum(axis=1))]
    ids = np.nonzero(relation)[1]
    # Authenticate the complete typed tape before indexing or doing arithmetic.
    # int(float_payload) would silently turn an invalid half-integer into a
    # valid value; unchecked trailing offsets could hide unexecuted reads.
    shapes = {"indptr": ((n+1,), np.int64), "indices": ((len(ids),), np.int32),
              "parent_offsets": ((5*n+1,), np.int64),
              "read_writers": ((4*len(ids),), np.int32), "values": ((5, n), np.int64)}
    require(set(z) == set(shapes), "q13 complete array census")
    for name, (shape, dtype) in shapes.items():
        require(z[name].shape == shape and z[name].dtype == dtype, "q13 exact array format: "+name)
    expected_offsets = np.r_[np.zeros(n+1, dtype=np.int64),
                             np.cumsum(np.tile(relation.sum(axis=1), 4))]
    require(np.array_equal(z["parent_offsets"], expected_offsets), "q13 complete read offsets")
    require(np.array_equal(z["indptr"], ptr) and np.array_equal(z["indices"], ids), "q13 exact metric reads")
    payload = z["values"]
    require(payload.shape == (5, n) and np.array_equal(payload[0], np.arange(1, n+1)), "q13 preparation")
    parents = [np.asarray([], dtype=np.int32) for _ in range(n)]
    for t in range(1, 5):
        for site in range(n):
            event = t*n+site
            expected = (t-1)*n+np.flatnonzero(relation[site])
            actual = z["read_writers"][z["parent_offsets"][event]:z["parent_offsets"][event+1]]
            require(np.array_equal(actual, expected), "q13 authenticated previous versions")
            total = 1+sum(int(payload.ravel()[x]) for x in actual)
            require(payload[t, site] == total, "q13 read law")
            parents.append(actual)
    require(np.array_equal(z["parent_offsets"][:n+1], np.zeros(n+1)), "q13 initial records do not read")
    require(z["parent_offsets"][-1] == len(z["read_writers"]) == row["reads"], "q13 read count")
    require((row["q"], row["dimension"], row["sites"], row["rounds"], row["events"]) == (13, dim, n, 4, 5*n), "q13 census")
    coords = integer+phi*labels
    centre_axis = min(range(q), key=lambda i: (abs(-floors[i]+phi*i-.5), i))
    centre = centre_axis*sum(q**i for i in range(dim))
    require(row["centre"] == centre and len(row["intervals"]) == 4, "q13 central tips")
    historical = read_json(ROOT / "evidence/source_net_causal_poset/source_net_causal_limit_receipt.json")
    historical_level = next(level for level in historical["levels"] if level["q"] == q)
    historical_family = next(f for f in historical_level["families"] if f["dimension"] == dim)
    previous = None
    for lag, item in enumerate(row["intervals"], 1):
        require((item["lag"], item["bottom"], item["top"]) == (lag, centre, lag*n+centre), "q13 interval tips")
        members, pairs, height = count_interval(parents, centre, lag*n+centre)
        check_order_row(item, members, pairs, height)
        old = historical_family["vertical_intervals"][lag-1]
        require(item["historical_counts_equal"] is True and len(members) == old["inclusive_event_count"]
                and pairs == old["strict_pair_count"], "q13 historical interval counts")
        require(item["counts_by_layer"] == np.bincount(members//n, minlength=lag+1).tolist(), "q13 layers")
        x, y = -int(floors[centre_axis]), centre_axis
        interior = all(qphi_sign(4*q*(u*u+v*v)-lag*lag, 4*q*(2*u*v+v*v)) >= 0
                       for u, v in ((x, y), (1-x, -y)))
        require(item["interior"] == bool(interior), "q13 interior/clipping flag")
        check_spatial(item["placement_in_units_of_L"], coords[members % n])
        if previous is not None:
            close(item["count_ratio_to_previous"], len(members)/previous, "q13 growth")
        else:
            require(item["count_ratio_to_previous"] is None, "q13 initial growth")
        previous = len(members)
    print(f"Verified q13 d={dim} metric identity, {len(z['read_writers']):,} versioned reads and four interval counts", flush=True)
    return row


def verify_q13(output, geometry):
    return verify_metric_family(output, geometry, read_json(output / "q13.json"), 3)


def verify_q13_controls(output, geometry):
    controls = read_json(output / "q13_controls.json")
    require(controls["schema"] == "oph.support_wiring.q13_controls.v1", "q13 controls schema")
    require(controls["addendum_sha256"] == digest(HERE / "CONTROL_ADDENDUM.md"), "q13 control declaration")
    require(controls["historical_receipt_sha256"] == digest(ROOT / "evidence/source_net_causal_poset/source_net_causal_limit_receipt.json"),
            "q13 control historical pin")
    require([row["dimension"] for row in controls["families"]] == [1, 2], "complete q13 control families")
    require(controls["flat_reference_ordering_fractions"] == {"1": "1/2", "2": "8/35"}
            and controls["flat_reference_mm_dimensions"] == {"1": 2, "2": 3}
            and controls["reference_is_acceptance_target"] is False, "q13 control comparison boundary")
    for dim, row in enumerate(controls["families"], 1):
        verify_metric_family(output, geometry, row, dim)
    return controls


def verify_controls(output, geometries, kernels=True):
    rows = read_json(output / "controls.json")
    matrices = np.load(output / "kernels.npz")
    require(len(rows) == 9, "all levels and wiring controls")
    for index, row in enumerate(rows):
        level, wiring = 3+index//3, ("isolated", "w3", "w12")[index%3]
        require((row["level"], row["wiring"]) == (level, wiring), "control ordering")
        g = geometries[level]
        c = len(g["faces"])
        phases = schedule(g, wiring)
        initial = np.asarray(initial_values(12*c), dtype=float)
        target = (np.repeat(initial.reshape(-1, 12).sum(axis=1)/12, 12)
                  if wiring == "isolated" else np.full(12*c, initial.sum()/(12*c)))
        endpoints = []
        require(len(row["confluence"]) == 2, "both finite schedules")
        for reverse, summary in enumerate(row["confluence"]):
            x = initial.copy()
            residuals = [float(np.dot(x-target, x-target))]
            for _ in range(4):
                for _, edges in (phases[::-1] if reverse else phases):
                    left, right = edges.T
                    # Independent delta-current implementation of the mean.
                    delta = (x[right]-x[left])/2
                    x[left] += delta
                    x[right] -= delta
                residuals.append(float(np.dot(x-target, x-target)))
            require(summary["reverse"] == bool(reverse), "control schedule direction")
            close(summary["squared_distance_by_sweep"], residuals, "finite confluence residuals", atol=2e-8)
            require(all(a >= b-1e-9 for a, b in zip(residuals, residuals[1:])), "quadratic descent")
            close(summary["max_distance_to_component_mean"], np.max(abs(x-target)), "finite control max residual")
            close(summary["sum_drift"], x.sum()-initial.sum(), "finite control total", atol=2e-8)
            drift = np.max(abs((x-initial).reshape(-1, 12).sum(axis=1))) if wiring == "isolated" else abs(x.sum()-initial.sum())
            close(summary["max_component_sum_drift"], drift, "component conservation", atol=2e-8)
            endpoints.append(x)
        close(row["finite_schedule_max_difference"], np.max(abs(endpoints[0]-endpoints[1])), "finite schedule dependence")
        degree = np.bincount(g["w12"][:, [0, 2]].ravel(), minlength=c)
        defect, regular = np.flatnonzero(degree == 11), np.flatnonzero(degree == 12)
        samples = [0] if wiring == "isolated" else sorted({int(defect[0]), int(defect[-1]), int(regular[0]), int(regular[-1])})
        require([(r["cell"], r["n"]) for r in row["kernels"]] == [(s, n) for s in samples for n in (1, 5, 30, 100, 300)], "frozen kernel sample census")
        projector, _ = frame_from_graph(g)
        q = np.eye(12)-np.ones((12, 12))/12
        if not kernels:
            continue
        ep = np.concatenate([x for _, x in phases])
        D = 2*len(ep)/c
        close(row["denominator"], D, "expectation normalization")
        left, right = ep.T
        adjacency = sparse.coo_matrix((np.ones(2*len(ep)), (np.r_[left, right], np.r_[right, left])), shape=(12*c, 12*c)).tocsr()
        operator = adjacency/D+sparse.diags(1-np.asarray(adjacency.sum(axis=1)).ravel()/D)
        if wiring == "isolated":
            operator = operator[:12, :12]
        for cell in samples:
            y = np.zeros((operator.shape[0], 12))
            y[12*cell:12*cell+12] = q
            for n in range(1, 301):
                y = operator@y
                norm = np.linalg.norm(y)
                if norm:
                    y /= norm
                if n not in (1, 5, 30, 100, 300):
                    continue
                # Full-field Gram equals the producer's local T^(2n) return.
                gram = y.T@y
                kernel = 12*gram/np.trace(gram)
                key = f"l{level}_{wiring}_c{cell}_n{n}"
                saved = next(r for r in row["kernels"] if r["cell"] == cell and r["n"] == n)
                require(saved["matrix"] == key, "kernel matrix name")
                close(matrices[key], kernel, "independent response Gram", atol=1e-8)
                close(saved["slow_share"], np.trace(projector@kernel)/12, "slow-band share", atol=2e-9)
                close(saved["eigenvalues"], np.linalg.eigvalsh(kernel)[::-1], "kernel eigenvalues", atol=1e-8)
            print(f"Verified kernels L{level} {wiring} cell {cell}", flush=True)
    return rows


def verify_inventory(output):
    manifest = read_json(output/"archive_manifest.json")
    require(manifest["schema"] == "oph.support_wiring.archive.v1", "archive schema")
    actual = {p.relative_to(output).as_posix() for p in output.rglob("*") if p.is_file() and p.name != "archive_manifest.json"}
    require(set(manifest["artifacts"]) == actual, "complete evidence inventory")
    for name, pin in manifest["artifacts"].items():
        require(".." not in Path(name).parts and not Path(name).is_absolute(), "safe archive path")
        path = output/name
        require(path.stat().st_size == pin["bytes"] and digest(path) == pin["sha256"], "archive artifact digest: "+name)
    expected_sources = {p.relative_to(ROOT).as_posix() for p in HERE.rglob("*")
                        if p.is_file() and "__pycache__" not in p.parts}
    expected_sources.add(".github/workflows/support-wiring.yml")
    require(set(manifest["sources"]) == expected_sources, "complete source inventory")
    for name, pin in manifest["sources"].items():
        require(digest(ROOT/name) == pin, "source digest: "+name)
    require(manifest["mirror"] == "evidence/source_net_causal_poset/support_wiring_receipt.json", "mirror location")
    if output.resolve() == DEFAULT.resolve():
        require((ROOT/manifest["mirror"]).read_bytes() == (output/"support_wiring_receipt.json").read_bytes(), "byte-exact causal-poset mirror")


def verify_receipt(output, stats, replay, intervals, q13, controls, q13_controls):
    r = read_json(output/"support_wiring_receipt.json")
    require(r["schema"] == "oph.support_wiring.paired_receipt.v1" and r["issue"] == 776, "paired receipt schema")
    require(r["specification_sha256"] == digest(HERE/"SPECIFICATION.md"), "receipt specification")
    require(r["geometry"] == read_json(HERE/"geometry/geometry.json"), "receipt geometry")
    require(len(r["wiring"]) == len(stats), "wiring census size")
    for reported, expected in zip(r["wiring"], stats):
        for key in expected:
            if key == "global_direction_dot":
                for moment in ("min", "mean", "max"):
                    close(reported[key][moment], expected[key][moment], "global antipodal statistic")
            else:
                require(reported[key] == expected[key], "wiring statistic: "+key)
    require(r["provenance_intervals"] == intervals, "receipt provenance readouts")
    require(r["record_metric_q13"] == q13, "receipt q13 readouts")
    require(r["record_metric_q13_controls"] == q13_controls, "receipt q13 control readouts")
    require(r["canonical_controls"] == controls, "receipt control readouts")
    manifest = replay[0]
    execution = {"law": "exact_pair_mean", "events": manifest["events"], "phases": len(manifest["chunks"]),
                 "mean_actions": replay[-1]["mean_actions"],
                 "glued_actions": sum(c["event_count"]//2 for c in manifest["chunks"] if c["kind"] == "glued"),
                 "preparation_writes": 12*1280, "refinement_writes": 12*(5120+20480),
                 "final_denominator_exponent": replay[-1]["final_exponent"], "final_chain": replay[-1]["final_chain"]}
    require(r["execution"] == execution, "receipt executed law and census")
    expected_boundary = {"declared_wiring_port_assignment_schedule_loads_and_joins": True,
                         "canonical_law_for_every_provenance_seam": True,
                         "physical_clock": False, "continuum_limit": False, "source_selection": False,
                         "M1_derived": False, "universal_finite_schedule_confluence": False,
                         "dimension_is_acceptance_target": False, "port_label_antipodes_imposed_on_gluing": False,
                         "join_is_commutative_cell_pullback_in_common_port_labels": True}
    require(r["claim_boundary"] == expected_boundary, "receipt nonclaims")
    h = r["historical_L6"]
    baseline = ROOT/"evidence/exact_federation_L6_canonical_20260909"
    require(h["archive"] == baseline.relative_to(ROOT).as_posix() and h["level"] == 6 and h["carriers"] == 81920, "historical regulator labels")
    require(set(h["sha256"]) == {p.name for p in baseline.iterdir() if p.is_file()}, "historical file census")
    for name, pin in h["sha256"].items():
        require(pin == digest(baseline/name), "historical input pin")
    kernel = read_json(baseline/"kernel_readout.json")
    require(h["isolated_kernel_readout"] == kernel["per_cell"][0]["isolated"], "historical isolated kernel readout")
    require(h["w3_kernel_summary"] == kernel["summary_glued"] and h["kernel_steps"] == kernel["steps"], "historical slow-band comparison")
    require(h["kernel_sample_count"] == len(kernel["cells"]), "historical sample count")
    require(h["isolated_n300_max_deviation_from_4P"] == kernel["isolated_max_abs_deviation_from_4_P_slow_at_n_300"], "historical isolated comparison")
    required_fields = {
        "canonical_w3": ("law", "sweeps", "attempts", "terminated", "budget", "max_abs_deviation_from_component_mean"),
        "historical_integer_control": ("law", "sweeps", "unique_quotient_hash_count", "quotient_hash_equals_expected_all")}
    for field, name in (("canonical_w3", "float_law_port_pair.json"), ("historical_integer_control", "integer_law_port_pair.json")):
        source = read_json(baseline/name)
        require(h[field] == {k: source[k] for k in required_fields[field]}, "complete historical law-labelled comparison")
    source = read_json(baseline/"float_law_isolated.json")
    require(h["canonical_isolated"] == {"law": source["law"], "schedules": len(source["schedules"]),
                                        "sweeps": source["sweeps"], "all_terminated": source["all_terminated"],
                                        "component_mean_within_1e-9_all": source["component_mean_within_1e-9_all"]}, "historical isolated confluence")
    print("Verified paired receipt, historical comparisons and causal-poset mirror", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=DEFAULT)
    parser.add_argument("--part", choices=["all", "trace", "kernels"], default="all")
    args = parser.parse_args()
    if args.part == "all":
        verify_inventory(args.archive)
    geometries, stats = verify_geometry()
    if args.part != "kernels":
        replay = verify_trace(args.archive / "trace", geometries)
        intervals = verify_intervals(args.archive, geometries, replay)
        q13 = verify_q13(args.archive, geometries[3])
        q13_controls = verify_q13_controls(args.archive, geometries[3])
    if args.part != "trace":
        controls = verify_controls(args.archive, geometries)
    if args.part == "all":
        verify_receipt(args.archive, stats, replay, intervals, q13, controls, q13_controls)
    print("SUPPORT_WIRING_VERIFIED " + args.part, flush=True)
