"""Equal-inventory free two-vacancy interferometry on disjoint winding guides."""

import argparse
import io
import json
from hashlib import sha256
from itertools import pairwise

import networkx as nx
import numpy as np
from scipy.linalg import eigh, expm

from . import geometry, harmonics, measure, preparation

HERE = measure.HERE
LOCK = "c3e68ac395ae8a7b397a73496617cfb5bac4c94d2c73f463d6effdaf1d5c720a"
T = 2.0 ** np.arange(-5, 14)
SHOTS = 2**22
EPS = np.sqrt(np.log(2 * 100000 / 0.01) / (2 * SHOTS))
SIGNS = np.array([1, -1])
FOUR = np.outer(SIGNS, SIGNS).ravel()


def second_cycle(graph, root, n):
    paths = nx.single_source_shortest_path(
        nx.bfs_tree(graph, root, sort_neighbors=sorted), root
    )
    records = []
    for a, b in graph.edges:
        if a not in paths or b not in paths:
            continue
        u, v = paths[a], paths[b]
        common = next(
            (i for i, p in enumerate(zip(u, v)) if p[0] != p[1]), min(len(u), len(v))
        )
        if len(u) + len(v) - 2 * common + 1 != n:
            continue
        cycle = u[: common - 1 : -1] + v[common - 1 :]
        i = cycle.index(min(cycle))
        ordered = cycle[i:] + cycle[:i]
        records.append(min(ordered, [ordered[0], *ordered[:0:-1]]))
    return min(records, default=None)


def stock():
    board = geometry.federation(2)
    graph = nx.Graph(board.seams)
    flights = json.loads((HERE / "geometry.json").read_text())["flights"]
    for A in sorted(
        (r["route"] for r in flights if r["level"] == 2), key=len, reverse=True
    ):
        rem = graph.copy()
        rem.remove_nodes_from(A)
        dist = nx.multi_source_dijkstra_path_length(graph, A)
        for root in sorted(rem, key=lambda x: (-dist[x], x)):
            B = second_cycle(rem, root, len(A))
            if B:
                break
        if B:
            break
    else:
        raise ValueError("No disjoint equal-length native BFS stock.")
    used = set(A) | set(B)
    assert len(used) == 2 * len(A)
    active = {tuple(sorted(edge)) for r in (A, B) for edge in pairwise(r + r[:1])}
    marks = {
        tuple(sorted((r[i], r[(i + 1) % len(r)])))
        for r in (A, B)
        for i in range(1, len(r), 2)
    }
    cuts = [e for e in board.seams if used.intersection(e) and e not in active]
    closed = graph.copy()
    closed.remove_edges_from(cuts)
    assert all(set(nx.node_connected_component(closed, r[0])) == set(r) for r in (A, B))
    distance = nx.multi_source_dijkstra_path_length(graph, used)
    d0 = max(graph, key=lambda x: (distance[x], -x))
    from_donor = nx.single_source_shortest_path_length(graph, d0)
    d1 = max(graph, key=lambda x: (min(distance[x], from_donor[x]), -x))
    assert min(distance[d0], distance[d1], nx.shortest_path_length(graph, d0, d1)) > 2
    access = graph.copy()
    access.remove_nodes_from([d0, d1])
    access.remove_edges_from(marks)
    parked = [min(graph.neighbors(d)) for d in (d0, d1)]
    assert all(distance[p] > 1 for p in parked)
    preparing = access.copy()
    preparing.remove_node(parked[1])
    assert nx.is_connected(preparing)
    ref = min(
        p
        for p in preparing
        if p not in used and p not in parked and board.degree[p] == 5
    )
    bypasses = []
    for r in (A, B):
        for i in range(1, len(r), 2):
            a, b = r[i], r[(i + 1) % len(r)]
            allowed = nx.restricted_view(preparing, [a], [])
            bypasses.append(nx.shortest_path_length(allowed, b, ref))
    # Both vacancies park beyond every shutter during equal-load writes. Open one
    # exterior entrance per guide, move its own vacancy in, and park it on the
    # first blank while writing that entrance closed. The other vacancy stays put.
    ingresses = []
    lengths = []
    shifts = []
    source_routes = []
    for seat, r in enumerate((A, B)):
        outside = nx.restricted_view(graph, used | {d0, d1, parked[1 - seat]}, [])
        distances = nx.single_source_shortest_path_length(outside, parked[seat])
        i, v = min(
            (i, int(v))
            for i in range(0, len(r), 2)
            for v in board.neighbours[r[i]]
            if v in distances and all(v != pair[0] for pair in ingresses)
        )
        shifts.append(i)
        source_routes.append(r[i:] + r[:i])
        ingresses.append([v, r[i]])
        lengths.append(distances[v] + 1)
    A, B = source_routes
    assert len(set(ingresses[0]) & set(ingresses[1])) == 0
    hA, ports = geometry.microscopic(A)
    hB, portsB = geometry.microscopic(B)
    np.testing.assert_array_equal(hA, hB)
    np.testing.assert_array_equal(ports, portsB)
    # On every exterior shutter a zero-load record removes the entire column on
    # all nine ordinary endpoint loads, including simultaneous nearby vacancies.
    table = measure.source.measure.blocks()
    dark = table[(0, 1)]
    assert np.max(abs(dark)) == 0
    return (
        dict(
            routes=[A, B],
            second_root=root,
            donors=[d0, d1],
            donor_clearance=[distance[d0], distance[d1]],
            off_route_phase=ref,
            exterior_cuts=len(cuts),
            prep_access_biconnected=nx.is_biconnected(access),
            held_endpoint_bypasses_checked=len(bypasses),
            longest_borrowed_phase_path=max(bypasses),
            sharp_stock_parking=parked,
            ingress_shutters=ingresses,
            even_rotations_to_reachable_ingress=shifts,
            vacancy_loading_path_lengths=lengths,
            two_vacancy_dimension=len(hA) ** 2,
        ),
        hA / 8,
        ports,
    )


def frame(i, L):
    v = np.zeros(2 * L)
    v[1:-1] = np.repeat(np.sin(np.pi * (L - i) * np.arange(1, L) / L), 2)
    mirror = v.reshape(-1, 2)[:, ::-1].ravel()
    cols = np.column_stack((v + mirror, v - mirror))
    return cols / np.linalg.norm(cols, axis=0)


def planes(h):
    L = len(h) // 2
    modes, rows = {}, []
    for i in range(1, L):
        F = frame(i, L)
        vals, rotation = eigh(F.T @ h @ F)
        waves = F @ rotation
        modes[i] = (waves.T, vals)
        s, a = F.T
        wedge = (np.outer(s, a) - np.outer(a, s)) / np.sqrt(2)
        per = np.trace(F.T @ h @ F)
        err = np.linalg.norm(h @ wedge + wedge @ h - per * wedge)
        product = np.outer(waves[:, 0], waves[:, 1])
        pres = np.linalg.norm(h @ product + product @ h - per * product)
        rows.append([float(per), float(err), float(pres)])
    return modes, np.array(rows)


def read_axes(native, ref, signal):
    terms = measure.source.contacts()
    assert abs(ref @ signal) < 1e-12
    initial = (ref + signal) / np.sqrt(2)
    projectors = np.array([initial, (ref - 1j * signal) / np.sqrt(2)])
    for state in [initial, *projectors]:
        np.testing.assert_allclose(
            preparation.compile_state(state, terms), state, atol=1e-12
        )
    e, v = native
    psi = v @ (np.exp(-1j * np.outer(e, T)) * (v.T @ initial)[:, None])
    return np.clip(abs(projectors.conj() @ psi) ** 2, 0, 1)


def joint(p, q, rng):
    counts = np.empty((2, 2, len(T), 4), dtype=int)
    for i in range(2):
        for j in range(2):
            for t, (a, b) in enumerate(zip(p[i], q[j], strict=True)):
                counts[i, j, t] = rng.multinomial(
                    SHOTS, [a * b, a * (1 - b), (1 - a) * b, (1 - a) * (1 - b)]
                )
    return counts


def unwrap(z, err):
    angles = np.angle(z)
    phase = angles[0]
    margin = np.pi - abs(phase) - err
    for angle in angles[1:]:
        val = angle + 2 * np.pi * round((2 * phase - angle) / (2 * np.pi))
        margin = min(margin, np.pi - abs(val - 2 * phase) - 3 * err)
        phase = val
    assert margin > 0
    return float(phase / T[-1]), float(margin)


def sum_read(counts):
    c = counts @ FOUR / SHOTS
    z = c[0, 0] - c[1, 1] + 1j * (c[0, 1] + c[1, 0])
    radius = 4 * np.sqrt(2) * EPS
    return *unwrap(z, np.arcsin(radius)), float(min(abs(z)))


def passage(raw):
    rng = np.random.default_rng(314144)
    result = {}
    N = harmonics.N_METER
    eps = np.sqrt(np.log(2 * 100000 / 0.01) / (2 * N))
    for name, (d, v) in measure.source.contacts().items():
        if name == "odd":
            continue
        local = np.array([[d, v], [v, d]]) / 8
        times = np.array([1 / 8, 3 / 16])
        p = np.array([abs(expm(-1j * t * local)[1, 0]) ** 2 for t in times])
        hits = rng.binomial(N, p)
        raw["meter_" + name] = np.column_stack((hits, N - hits))
        lo, hi = (
            np.arcsin(np.sqrt(np.clip(hits[:, None] / N + [-eps, eps], 0, 1)))
            / times[:, None]
        ).T
        assert max(lo) < min(hi)
        result[name] = [float(max(lo)), float(min(hi))]
    return result


def packet_reads(native, ports, rng, raw):
    n = max(ports) + 1
    cells = ports // 2
    mid = n // 4
    out = []
    t = np.array([0, 2, 4, 6])
    for q in [-1.6, -1.3, 1.3, 1.6]:
        k = np.pi - q
        z, g = measure.symbol(k)
        u = np.exp(-(((cells - mid) / 4.0) ** 2) + 1j * k * cells)
        u *= np.where(ports % 2 == 0, 1, z.conjugate() / g)
        u[(ports < 4) | (ports >= n - 4)] = 0
        u /= np.linalg.norm(u)
        ev, vec = native
        psi = (np.exp(-1j * np.outer(t, ev)) * (vec.T @ u)) @ vec.T
        p = geometry.physical(psi, ports, n)
        counts = np.array([rng.multinomial(SHOTS, row / row.sum()) for row in p])
        raw[f"packet_{q}"] = counts
        x = np.arange(n) // 2
        out.append(
            dict(
                q=q,
                sampled_speed=np.polyfit(t, (counts / SHOTS) @ x, 1)[0],
                exact_speed=np.polyfit(t, p @ x, 1)[0],
                bulk_point_speed=-2 * np.sin(q) / g,
                max_boundary_weight=max(p[:, :2].sum(1) + p[:, -2:].sum(1)),
            )
        )
    return out


def scalar(q):
    return 0.61 * (1 - np.cos(q)) - 0.09 * (1 - np.cos(2 * q))


def run():
    assert sha256((HERE / "COPIES.md").read_bytes()).hexdigest() == LOCK
    info, h, ports = stock()
    pos, checks = planes(h)
    _, c1 = planes(geometry.microscopic(info["routes"][0], 1)[0] / 8)
    np.testing.assert_allclose(checks[:, 0], 9, atol=1e-12)
    np.testing.assert_allclose(c1[:, 0], 2, atol=1e-12)
    assert max(checks[:, 1:].max(), c1[:, 1:].max()) < 1e-12
    native = eigh(h)
    L = len(h) // 2
    rest = (-1.0) ** (np.arange(len(h)) // 2) / np.sqrt(len(h))
    np.testing.assert_allclose(h @ rest, 8 * rest, atol=1e-13)
    rng = np.random.default_rng(313144)
    raw, rows = {}, []
    calibrated = passage(raw)
    for i in range(L):
        j = L // 2 if i != L // 2 else L // 3
        r_lower, r_upper = pos[j][0]
        wave = pos[i][0][1] if i else rest
        a, b = read_axes(native, r_upper, wave), read_axes(native, r_lower, wave)
        pair = joint(a, b, rng)
        raw[f"pair_{i}"] = pair
        full, margin, vis = sum_read(pair)
        if i:
            shape, raw[f"space_{i}"] = harmonics.spatial(rest, wave, rng, L)
            assert shape["accepted_indices"] == [i]
            k = np.pi * shape["best_index"] / L
            p = read_axes(native, rest, wave)
            y = rng.binomial(SHOTS, p)
            raw[f"recoil_{i}"] = np.stack([y, SHOTS - y], axis=-1)
            rec, _ = unwrap(
                2 * y[0] / SHOTS - 1 + 1j * (2 * y[1] / SHOTS - 1),
                np.arcsin(2 * np.sqrt(2) * EPS),
            )
        else:
            k, rec, shape = 0.0, 0.0, {"known_stationary_carrier": True}
        # A scalar upper band with matching lower partner and the same reference
        # rule is intentionally given access to the same two-copy clock.
        qs = np.pi * j / L
        energies = [
            3.5 + scalar(k) - (3.5 + scalar(qs)),
            3.5 + scalar(k) + (3.5 + scalar(qs)),
        ]
        af, bf = [(1 + np.array([np.cos(e * T), np.sin(e * T)])) / 2 for e in energies]
        foil = joint(af, bf, rng)
        raw[f"foil_{i}"] = foil
        rows.append(
            dict(
                index=i,
                reference_index=j,
                spatial=shape,
                q=k,
                recoil=rec,
                energy=full / 2,
                min_joint_visibility=vis,
                unwrap_margin=margin,
                foil_energy=sum_read(foil)[0] / 2,
            )
        )
    central = L // 2
    # A physically occupied upper/lower pair at distinct wavelengths is also
    # stationary but has a different energy. Record its read against two B's.
    up1 = pos[1][0][1]
    p = read_axes(native, up1, rest)
    controls = {}
    for name, v in [
        ("equal_side", pos[central][0][1]),
        ("mismatched_momenta", pos[central][0][0]),
    ]:
        y = joint(p, read_axes(native, v, rest), rng)
        raw[name] = y
        controls[name] = sum_read(y)[0] / 2
    flat = joint(p, np.full_like(p, 0.5), rng)
    raw["dephased"] = flat
    c = flat @ FOUR / SHOTS
    controls["dephased_max_contrast"] = max(
        abs(c[0, 0] - c[1, 1] + 1j * (c[0, 1] + c[1, 0]))
    )
    shifted = (native[0] + 1.75, native[1])
    shift_err = float(np.max(abs(read_axes(shifted, up1, rest) - p)))
    assert shift_err < 1e-10
    controls["uniform_shift_local_read_max_error"] = shift_err
    q, d, D, foil = np.array(
        [[r[f] for r in rows] for f in ("q", "recoil", "energy", "foil_energy")]
    )
    x = q[1:3] ** 2
    curvature = np.polyfit(x, 2 * d[1:3] / x, 1)[1]
    a, b = (np.mean(calibrated[k]) for k in ("blank", "even"))
    C2 = a * b
    inv = d * (2 * D[0] + d)
    held = np.arange(3, L)
    residual = inv - C2 * q * q
    design = np.column_stack([np.ones(L), 1 - np.cos(q), 1 - np.cos(2 * q)])

    def fit_error(spectrum):
        train, test = design[::2], design[1::2]
        return np.max(
            abs(
                test @ np.linalg.lstsq(train, spectrum[::2], rcond=None)[0]
                - spectrum[1::2]
            )
        )

    energy_rad = np.arcsin(4 * np.sqrt(2) * EPS) / (2 * T[-1])
    recoil_rad = np.arcsin(2 * np.sqrt(2) * EPS) / T[-1]
    x1, x2 = x
    curvature_radius = 2 * recoil_rad / (x2 - x1) * (x2 / x1 + x1 / x2)
    c2lo, c2hi = np.prod([calibrated[k] for k in ("blank", "even")], axis=0)
    linear = abs((inv - C2 * (2 * np.sin(q / 2)) ** 2)[held])
    counts_err = (
        2 * (abs(D[0]) + max(d) + energy_rad + recoil_rad) * recoil_rad
        + 2 * max(d) * energy_rad
    )
    worst_lattice = max(
        linear + counts_err + max(C2 - c2lo, c2hi - C2) * (2 * np.sin(q[held] / 2)) ** 2
    )
    pair_err = max(abs(D - D[0] - d))
    assert pair_err + 2 * energy_rad + recoil_rad < 0.00002 and worst_lattice < 0.003
    details = dict(
        protocol_sha256=LOCK,
        stock=info,
        shots_per_joint_bank=SHOTS,
        dyadic_times=T.tolist(),
        reference_traces=[min(checks[:, 0]), max(checks[:, 0])],
        max_product_residual=max(checks[:, 2]),
        max_wedge_residual=max(checks[:, 1]),
        contrast_one_trace_range=[min(c1[:, 0]), max(c1[:, 0])],
        contrast_one_wedge_residual=max(c1[:, 1]),
        paired_frequency_radius=energy_rad,
        recoil_frequency_radius=recoil_rad,
        meter_intervals=calibrated,
        rest=D[0],
        inverse_kinetic_mass=curvature,
        inverse_mass_count_radius=curvature_radius,
        kinetic_mass=1 / curvature,
        contact_speed=np.sqrt(C2),
        rest_and_curvature_speed=np.sqrt(D[0] * curvature),
        finite_grid_max_upper_secant=max(np.diff(d) / np.diff(q)),
        packet_centroids=packet_reads(native, ports, rng, raw),
        recoil_consistency_max_error=pair_err,
        held_lattice_invariant_max_error=max(linear),
        held_lattice_max_including_intervals=worst_lattice,
        held_continuum_max_error=max(abs(residual[held])),
        held_continuum_pass_indices=held[abs(residual[held]) < 0.003].tolist(),
        curvature_matched_one_hop_max_recoil_error=max(
            abs(d[held] - curvature * (1 - np.cos(q[held])))
        ),
        host_rest_lattice_invariant_max_error=max(
            abs((d * (16 + d) - C2 * 4 * np.sin(q / 2) ** 2)[held])
        ),
        native_scalar_two_range_odd_max_error=fit_error(D),
        planted_two_range_odd_max_error=fit_error(foil),
        controls=controls,
        rows=rows,
    )
    archive = io.BytesIO()
    np.savez_compressed(archive, **raw)
    blob = archive.getvalue()
    details["counts_sha256"] = sha256(blob).hexdigest()
    return json.dumps(details, sort_keys=True, indent=2) + "\n", blob


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    write = parser.parse_args().write
    out, blob = run()
    for name, data in (("copies.json", out.encode()), ("copies_counts.npz", blob)):
        path = HERE / name
        if write:
            path.write_bytes(data)
        elif path.read_bytes() != data:
            raise AssertionError(f"Copy replay mismatch: {name}")
    print("Verified passive two-copy vacancy interferometer.")
