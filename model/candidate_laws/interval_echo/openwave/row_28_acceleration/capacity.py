"""Measure sharp-record cut capacity and independently audit two timed signal routes."""

import argparse
import heapq
import json

import numpy as np
from scipy.sparse.csgraph import connected_components
from model.candidate_laws.interval_echo import candidate

from .experiment import HERE, layout, geometry
from .signals import adjacency


def path_records(level, commit, start, emit, targets, passing_writes=False):
    n, edges, *_ = layout(level)
    incident = [[] for _ in range(n)]
    for (a, b), deadline in zip(edges, commit, strict=True):
        incident[a].append((int(b), deadline))
        incident[b].append((int(a), deadline))
    queue = [(emit, start)]
    parent = {start: None}
    stamps = {start: emit}
    while queue:
        t, p = heapq.heappop(queue)
        for v, deadline in incident[p]:
            if v not in stamps and (passing_writes or not deadline or t < deadline):
                stamps[v] = t + 1
                parent[v] = p
                heapq.heappush(queue, (t + 1, v))
    paths = []
    for v in targets:
        if v not in stamps:
            paths.append(None)
            continue
        way = [v]
        while way[-1] != start:
            way.append(parent[way[-1]])
        paths.append(dict(arrival=stamps[v], ports=way[::-1]))
    return paths


def cut_panel(level, seed, dilute, commit, rng):
    n, edges, _, stations, _ = layout(level)
    live = commit == 0
    _, labels = connected_components(adjacency(n, edges[live]), directed=False)
    ds = geometry(n, edges, stations, live)
    cuts = labels[edges[:, 0]] != labels[edges[:, 1]]
    assert not np.any(cuts & live)
    rows = []
    components = {}
    for i, p in enumerate(stations):
        group = int(labels[p])
        if group not in components:
            boundary = np.flatnonzero(cuts & np.any(labels[edges] == group, axis=1))
            components[group] = dict(
                seams=boundary.tolist(),
                coin_ones=rng.binomial(2**18, 0.5, len(boundary)).tolist(),
            )
        b = len(components[group]["seams"])
        rad = int(ds[i].max())
        rows.append(
            dict(
                origin=p,
                component=group,
                accessible=int(np.sum(ds[i] >= 0)),
                farthest=rad,
                cut_bits=b,
                radius_squared_over_nats=rad**2 / (b * np.log(2)) if b else None,
            )
        )
    return dict(
        level=level, seed=seed, dilute=dilute, stations=rows, components=components
    )


def run():
    cases = json.loads((HERE / "signals.json").read_text())["cases"]
    clocks = json.loads((HERE / "homogeneity.json").read_text())
    law = candidate.native_law
    for w in ((0, 2), (2, 0)):
        output = law.echo()[:, law.SRC.index(w)].reshape(-1, law.N_SRC)
        for c in (0, 1):
            amp = output[1 + law.LABELS.index((c, w[0] - w[1])), law.SRC.index((1, 1))]
            np.testing.assert_allclose(amp, 2**-0.5, atol=1e-14)
    panels = []
    traces = []
    rng = np.random.default_rng(280901)
    with np.load(HERE / "signals.npz") as sig, np.load(HERE / "records.npz") as survey:
        for p in cases:
            level, seed, dilute = p["level"], p["seed"], p["dilute"]
            commit = sig[f"{seed}_1_commits"] if dilute else survey[f"epoch_{seed}"]
            panels.append(cut_panel(level, seed, dilute, commit, rng))
            reports = [
                v
                for v in clocks
                if v["seed"] == seed
                and v["dilute"] == dilute
                and v["indices"] == [0, 1]
            ]
            witness = reports[0]["witness"]
            result = []
            for e, times in ((0, [witness["first"]] * 2), (4, witness["second"])):
                a = path_records(
                    level, commit, witness["origin"], e, witness["destinations"]
                )
                assert [w["arrival"] for w in a] == times
                b = path_records(
                    level, commit, witness["origin"], e, witness["destinations"], True
                )
                assert all(w["arrival"] == e + witness["initial_range"] for w in b)
                result.append(dict(emit=e, native=a, passing_write=b))
            traces.append(dict(level=level, seed=seed, dilute=dilute, traces=result))
    return dict(coin_shots=2**18, panels=panels, traces=traces)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    result = json.loads(json.dumps(run()))
    path = HERE / "capacity.json"
    if cli.parse_args().write:
        path.write_text(json.dumps(result, indent=1) + "\n")
    else:
        assert result == json.loads(path.read_text())
        print(
            "Verified reachable cut-record capacities and native witness trajectories."
        )


if __name__ == "__main__":
    main()
