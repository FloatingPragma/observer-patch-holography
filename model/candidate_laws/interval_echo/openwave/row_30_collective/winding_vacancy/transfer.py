"""Direct vacancy recoil and exchange against the independent winding-gateway receiver."""

import argparse
import hashlib
import io
import json
from itertools import combinations, pairwise, product

import networkx as nx
import numpy as np
from scipy.linalg import eigh, expm
from scipy.optimize import least_squares

from model.experiments.interval_echo_angular.charged_orientation.two_winding_loops import (
    measure as receiver,
)
from . import density, measure, phase, preparation

HERE = measure.HERE
PROTOCOL = "9ca0f17960b41dce60e3ec9fd0164be1c235ff76dda07a6140cc23b3e744bd81"
SCAN = np.arange(128) / 4
TIMES = np.r_[1 / 64, np.arange(128) / 8]


def embedded():
    route, _ = measure.source.material.geometry()
    board = measure.source.material.federation(2)
    graph = nx.Graph(board.seams)
    used = {v // 12 for v in route}

    def triangle(v):
        for a, b in combinations(board.neighbours[v], 2):
            if board.degree[a] == board.degree[b] == 5 and graph.has_edge(a, b):
                return [v, a, b]
        return []

    for u, v in board.external:
        if u // 12 in used or v // 12 in used:
            continue
        ref = triangle(u) + triangle(v)
        if len(ref) != 6:
            continue
        used_ports = set(route) | set(ref)
        access = graph.copy()
        access.remove_nodes_from(used_ports - {route[0], ref[3]})
        if not nx.has_path(access, route[0], ref[3]):
            continue
        path = nx.shortest_path(access, route[0], ref[3])
        break
    assert len(set(route + ref)) == len(route) + 6 and len(path) > 2
    edges = set(tuple(sorted((ref[a], ref[b]))) for a, b in receiver.EDGES)
    ring_edges = set(tuple(sorted(e)) for e in pairwise(route + route[:1]))
    ring_marks = set(
        tuple(sorted((route[i], route[(i + 1) % len(route)])))
        for i in range(1, len(route), 2)
    )
    ref_marks = set(tuple(sorted((ref[a], ref[b]))) for a, b in receiver.EDGES[:6])
    cuts = {
        e for e in board.seams if set(e) & used_ports and e not in edges | ring_edges
    }
    left, right = tuple(sorted(path[:2])), tuple(sorted(path[-2:]))
    assert left in cuts and right in cuts
    stops = nx.multi_source_dijkstra_path_length(graph, used_ports | set(path))
    donor = max(graph, key=lambda p: (stops[p], -p))
    assert stops[donor] > 2
    prepgraph = graph.copy()
    prepgraph.remove_node(donor)
    prepgraph.remove_edges_from(ring_marks | ref_marks)
    assert nx.is_biconnected(prepgraph)
    other = graph.copy()
    other.remove_edges_from(cuts)
    assert set(nx.node_connected_component(other, route[0])) == set(route)
    assert set(nx.node_connected_component(other, ref[0])) == set(ref)
    # During closure the vacancy parks on the ring's first blank or the receiving
    # gateway's blank bridge. Both outside entrances then have ordinary loads on both arms.
    parked = {route[1], ref[0]}
    assert not parked.intersection((*left, *right))
    assert all(edge in graph.edges for edge in edges | ring_edges)
    return dict(
        signal_ports=route,
        reference_ports=ref,
        interferometer_route=path,
        compensator=donor,
        prep_bypasses_biconnected=True,
        common_filled_ingress_shutters=[left, right],
        parking_ports=sorted(parked),
        exterior_shutters=len(cuts),
    )


def space_line(vals, cells, held):
    def fit(k):
        m = phase.spatial_basis(cells, k)
        c = np.linalg.lstsq(m, vals, rcond=None)[0]
        return np.linalg.norm(vals - m @ c)

    grid = np.linspace(0.00001, np.pi - 0.00001, 145)
    q = density.bounded(grid, fit)
    c = np.linalg.lstsq(phase.spatial_basis(cells, q), vals, rcond=None)[0]
    return q, phase.spatial_basis(held, q) @ c


def read_recoil(bank):
    cells = np.arange(3, 35)
    held = np.arange(35, 70)
    z = phase.quadratures(bank)
    y = z[:, 2 * cells]
    # Estimate two signed eigenphases on the entire principal interval from charge counts.
    u = np.linalg.svd(y, full_matrices=False)[0][:, :2]
    f = np.sort(
        np.angle(np.linalg.eigvals(np.linalg.lstsq(u[:-1], u[1:], rcond=None)[0]))
        / 0.25
    )

    def components(freq):
        d = np.exp(1j * np.outer(SCAN, freq))
        c = np.linalg.lstsq(d, y, rcond=None)[0]
        return c, d @ c - y

    def residual(freq):
        _, r = components(freq)
        return np.r_[r.real.ravel(), r.imag.ravel()]

    f = least_squares(residual, f, ftol=1e-11, gtol=1e-12).x
    c, _ = components(f)
    fitted = [space_line(amp, cells, held) for amp in c]
    forecast = np.exp(1j * np.outer(SCAN[1::2], f)) @ np.array([v[1] for v in fitted])
    return dict(
        upper_recoil=-f[0],
        lower_relative_rest=-f[1],
        upper_q=fitted[0][0],
        lower_q=fitted[1][0],
        held_profile_rms=float(
            np.sqrt(np.mean(abs(forecast - z[1::2, 2 * held]) ** 2))
        ),
    )


def acquire_recoil(native, raw):
    rng = np.random.default_rng(307144)
    _, b = phase.clocks(*native[1:])
    banks, rows = [], []
    for j in range(1, 73):
        frame = density.states(j)
        probe = frame @ np.array([1, 1j]) / np.sqrt(2)
        # Occupying both real eigenmodes with equal weight does not require their frequencies.
        bank = phase.read(native, b, probe, SCAN, 2**19, rng)
        banks.append(bank.astype(np.int32))
        rows.append(dict(j=j, **read_recoil(bank)))
    raw["recoil_counts"] = np.array(banks)
    q = np.array([row["upper_q"] for row in rows])
    energy = np.array([row["upper_recoil"] for row in rows])
    low = np.array([row["lower_relative_rest"] for row in rows])
    ql = np.array([row["lower_q"] for row in rows])
    x = 1 - np.cos(q)
    train = np.arange(1, 73) % 2 == 0
    forecasts = dict(
        native=np.sqrt(12.25 + 4 * x) - 3.5,
        continuum=np.sqrt(12.25 + 2 * q * q) - 3.5,
        host=np.sqrt(64 + (32 / 7) * q * q) - 8,
        curvature_matched_hopping=(4 / 7) * x,
    )
    fits = {}
    for name, m in [
        ("nearest", np.column_stack((np.ones(72), x))),
        ("next_nearest", np.column_stack((np.ones(72), x, 1 - np.cos(2 * q)))),
    ]:
        pars = np.linalg.lstsq(m[train], energy[train], rcond=None)[0]
        forecasts[name] = m @ pars
        fits[name] = pars.tolist()
    errors = {}
    for key, c in forecasts.items():
        err = c - energy
        errors[key] = dict(
            max_all=max(abs(err)),
            max_odd=max(abs(err[~train])),
            max_long=max(abs(err[q < 0.45])),
            n_held_above_gate=int(sum(abs(err[~train]) > 0.003)),
        )
    return dict(
        rows=rows,
        fits=fits,
        errors=errors,
        all_held_spatial_max=max(x["held_profile_rms"] for x in rows),
        max_lower_error=max(
            abs(low - (-3.5 - np.sqrt(12.25 + 8 * np.sin(ql / 2) ** 2)))
        ),
        q_range=[min(q), max(q)],
        long_count=int(sum(q < 0.45)),
    )


def probes(native):
    _, e, v = native
    h = (v * e) @ v.T
    r, b = phase.clocks(e, v)
    out = [("upper_rest", 8.0, b), ("upper_top", 9.0, r)]
    for j in (69, 64, 57, 49, 36, 21):
        frame = density.states(j)
        val, vec = eigh(frame.T @ h @ frame)
        for tag, i in [("lower", 0), ("upper", 1)]:
            wave = frame @ vec[:, i]
            np.testing.assert_allclose(h @ wave, val[i] * wave, atol=3e-13)
            out.append((f"{tag}_{j}", val[i], wave))
    return out


def ref_evolution(h, cols, step=None):
    if step is None:
        e, v = eigh(h)
        amp = np.einsum(
            "ai,ti,bi->tab",
            v[cols],
            np.exp(-1j * np.outer(TIMES, e)),
            v[cols],
            optimize=True,
        )
        return np.diagonal(amp, axis1=1, axis2=2), np.ones((len(TIMES), len(cols)))
    r = expm(-1j * step * h)[np.ix_(cols, cols)]
    powers = []
    for t in TIMES:
        n = round(t / step)
        assert abs(n * step - t) < 1e-13
        powers.append(np.linalg.matrix_power(r, n))
    p = np.array(powers)
    return np.diagonal(p, axis1=1, axis2=2), np.sum(abs(p) ** 2, axis=1)


def gateway(native, raw):
    words, ix, terms = receiver.build(True)
    sub = [ix[3, a + 8 * b] for a, b in product((0, 7), repeat=2)]
    h8 = sum(terms)
    outside = np.array([i for i in range(len(words)) if i not in sub])
    np.testing.assert_array_equal(h8[np.ix_(sub, sub)], 36 * np.eye(4))
    escape = h8[np.ix_(outside, sub)]
    np.testing.assert_array_equal(escape.T @ escape, 1040 * np.eye(4))
    h = h8 / 8
    poles = eigh(h, eigvals_only=True)
    rng = np.random.default_rng(308144)
    reference = dict(
        complete_levels=poles.tolist(),
        nearest_passive_pole_error=min(abs(poles - 4.5)),
        exact_gateway_centroid=4.5,
        exact_gateway_variance=16.25,
        reference_words=[words[i] for i in sub],
    )
    energies = probes(native)
    data = []
    for mode, step in [
        ("free", None),
        ("check_64", 1 / 64),
        ("check_256", 1 / 256),
        ("check_1024", 1 / 1024),
        ("dephased", 1 / 1024),
    ]:
        amp, success = ref_evolution(h, sub, step)
        draws, stats = [], []
        for label, E, _ in energies:
            cross = np.exp(-1j * TIMES[:, None] * E) * amp.conj()
            coherent = cross if mode != "dephased" else np.zeros_like(cross)
            outcomes = []
            for visibility in (coherent.real, coherent.imag):
                bright = (1 + abs(amp) ** 2 + 2 * visibility) / 4
                dark = (1 + abs(amp) ** 2 - 2 * visibility) / 4
                other = (success - abs(amp) ** 2) / 2
                fail = (1 - success) / 2
                p = np.stack([bright, dark, other, fail], axis=-1)
                assert p.min() > -2e-12
                p = np.maximum(p, 0)
                outcomes.append(p / p.sum(axis=-1, keepdims=True))
            p = np.array(outcomes)
            shots = density.samples(p, 2**20, rng)
            read = (
                shots[0, :, :, 0]
                - shots[0, :, :, 1]
                + 1j * (shots[1, :, :, 0] - shots[1, :, :, 1])
            ).mean(axis=1) / 2**20
            true = coherent.mean(axis=1)
            eligible = abs(true) > 0.1
            center = np.exp(1j * TIMES * (E - 4.5))
            stats.append(
                dict(
                    label=label,
                    signal_formation=E,
                    exposures_above_visibility=int(sum(eligible)),
                    largest_phase_error=float(
                        max(abs(np.angle((read * center)[eligible])))
                    )
                    if any(eligible)
                    else None,
                    exact_phase_error=float(
                        max(abs(np.angle((true * center)[eligible])))
                    )
                    if any(eligible)
                    else None,
                    short_read_slope=(
                        -np.angle(read[0]) * 64 if mode != "dephased" else None
                    ),
                    exact_short_slope=(
                        -np.angle(true[0]) * 64 if mode != "dephased" else None
                    ),
                    single_branch_expected=E - 4.5,
                )
            )
            draws.append(shots)
        raw["gateway_" + mode] = np.array(draws)
        data.append(
            dict(
                mode=mode,
                minimum_conditional_reference_survival=float(success.min()),
                minimum_reference_return_modulus=float(min(abs(amp.mean(axis=1)))),
                minimum_interferometric_visibility=float(min(abs(amp.mean(axis=1))))
                if mode != "dephased"
                else 0.0,
                phase_rows=stats,
            )
        )
    return dict(reference=reference, modes=data)


def run():
    assert hashlib.sha256((HERE / "TRANSFER.md").read_bytes()).hexdigest() == PROTOCOL
    chart = embedded()
    native = density.native(chart["signal_ports"], 8)
    # Verify that the signal eigenmodes used by the final interferometer are native-address preparable.
    terms = measure.source.contacts()
    for _, _, target in probes(native):
        np.testing.assert_allclose(
            preparation.compile_state(target, terms), target, atol=1e-12
        )
    raw = {}
    res = dict(
        protocol_sha256=PROTOCOL,
        physical_chart=chart,
        source_commit="fc4d85caa25158b01b140c91954443e820c3ad1f",
        recoil=acquire_recoil(native, raw),
        gateway=gateway(native, raw),
    )
    stream = io.BytesIO()
    np.savez_compressed(stream, **raw)
    data = stream.getvalue()
    res["counts_sha256"] = hashlib.sha256(data).hexdigest()
    return json.dumps(res, indent=2, sort_keys=True, default=int) + "\n", data


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    write = p.parse_args().write
    output, counts = run()
    for filename, blob in [
        ("transfer.json", output.encode()),
        ("transfer_counts.npz", counts),
    ]:
        path = HERE / filename
        if write:
            path.write_bytes(blob)
        else:
            assert path.read_bytes() == blob, f"replay mismatch: {filename}"
    print(
        "Single-charge transfers written."
        if write
        else "Single-charge transfers verified."
    )
