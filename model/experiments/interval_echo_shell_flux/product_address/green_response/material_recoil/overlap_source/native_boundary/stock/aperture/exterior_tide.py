"""Initial material tides on equal-radius chords outside every written contact."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
import json

import numpy as np
from scipy.sparse import diags

from . import coherent_tide as coherent
from . import exterior_source

tidal = coherent.tidal
HERE = tidal.HERE
RADII = (2, 4, 6, 7, 8, 10, 12)


def fields(h, ends, cost, adj, pairs):
    i, j = pairs.T
    base, used = tidal.chart(h, ends, cost)
    lap = adj - diags(np.asarray(adj.sum(1)).ravel())

    def action(s):
        return np.asarray(lap @ s.T).T.astype(np.int64)

    def chord(w):
        return w[i, j] + w[j, i] - w[i, i] - w[j, j]

    d = base[i, j] + cost * used[i, j]
    s0, s1, s2 = [chord(action(w)) for w in (base**2, 2 * base * used, used)]
    curvature = (s0 + cost * s1 + cost**2 * s2) / d**2
    l0, l1 = action(base), action(used)
    pair = 2 * (l0[i, j] + l0[j, i] + cost * (l1[i, j] + l1[j, i])) / d
    # Both simultaneous probes already occupy this bond; its swap keeps the
    # unordered pair fixed. Remove the two nonexistent transitions to zero range.
    pair += 4 * np.asarray(adj[i, j]).ravel().astype(int)
    return d, curvature, pair


def summarize(values):
    n = len(values)
    if not n:
        return dict(chords=0)
    return dict(
        chords=n,
        inward=sum(v < 0 for v in values),
        outward=sum(v > 0 for v in values),
        silent=sum(v == 0 for v in values),
        mean=str(sum(values) / n),
        minimum=str(min(values)),
        maximum=str(max(values)),
    )


def run():
    _, w, *_ = tidal.lab.energy.pair.setup()
    board, _, marks, *_ = tidal.lab.energy.returned(w)
    adj, h, _ = tidal.flight.records.thermal_radar.chart(board, marks)
    apparatus = json.loads((HERE / "tidal_release.json").read_bytes())["schedules"][
        "forward"
    ]["all"]["1"]
    costs = {k: Q(v["chart_cost"]["fastest_relay"]) for k, v in apparatus.items()}
    distance = np.min(h[np.ravel(marks)], axis=0)
    shell = {}
    for r in range(1, int(h[219].max()) + 1):
        stations = np.flatnonzero((h[219] == r) & (distance > 0))
        pairs = np.array(list(combinations(stations, 2)), int)
        if len(pairs) == 0:
            shell[str(r)] = {"stations": len(stations), "chords": 0}
            continue
        calc = {
            name: fields(h, marks[-1], cost, adj, pairs) for name, cost in costs.items()
        }
        panel = {}
        for src, ref, ordinary in [
            ("source", "original", ("single", "ordinary")),
            (
                "positive",
                "reference",
                ("reference_plus_single", "reference_plus_double"),
            ),
        ]:
            measured = {}
            for bank, mask in [
                ("all", np.ones(len(pairs), bool)),
                ("equal_initial_range", calc[src][0] == calc[ref][0]),
                (
                    "distant_equal_range",
                    (calc[src][0] == calc[ref][0]) & (distance[pairs].min(1) >= 2),
                ),
            ]:
                contrast = {}
                for name in (src, *ordinary):
                    d, c, p = calc[name]
                    d0, c0, p0 = calc[ref]
                    contrast[name] = dict(
                        projected=summarize((c - c0)[mask]),
                        pair=summarize((p - p0)[mask]),
                        absolute_projected=summarize((d * c - d0 * c0)[mask]),
                    )
                response = contrast[src]["projected"]
                prices = {}
                if response["chords"]:
                    for donor in ordinary:
                        base = Q(contrast[donor]["projected"]["mean"])
                        prices[donor] = (
                            str(Q(response["mean"]) / base) if base < 0 else None
                        )
                measured[bank] = dict(contrasts=contrast, ordinary_mass_ratios=prices)
                if mask.any() and bank == "distant_equal_range":
                    c = calc[src][1] - calc[ref][1]
                    for tag, fun in [
                        ("largest_inward", np.argmin),
                        ("largest_outward", np.argmax),
                    ]:
                        idx = np.flatnonzero(mask)[fun(c[mask])]
                        measured[bank][tag] = dict(
                            pair=pairs[idx].tolist(),
                            change=str(c[idx]),
                            common_length=str(calc[src][0][idx]),
                        )
            panel[src] = measured
        shell[str(r)] = dict(
            stations=len(stations), chords=len(pairs), comparisons=panel
        )
    local = exterior_source.source_bank()
    example = shell["4"]["comparisons"]["source"]["distant_equal_range"][
        "largest_inward"
    ]["pair"]
    held = exterior_source.held_witness(
        h, marks[-1], adj, example, local["preparations"], costs
    )
    support = {}
    for name, reading in local["preparations"].items():
        cost = 1 / Q(reading["success"])
        d, proj, _ = fields(h, marks[-1], cost, adj, np.array([example]))
        support[name] = dict(
            cost=str(cost), initial_range=str(d[0]), projected_acceleration=str(proj[0])
        )
    assert (
        Q(support["source"]["projected_acceleration"])
        < Q(support["original"]["projected_acceleration"])
        < Q(support["low"]["projected_acceleration"])
    )
    assert (
        support["source"]["initial_range"]
        == support["original"]["initial_range"]
        == support["low"]["initial_range"]
    )
    return dict(
        same_charge_witness=support,
        local_source=local,
        held_four_hop=held,
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "EXTERIOR_TIDE.md",
                HERE / "coherent_tide.json",
                HERE / "tidal_release.json",
            )
        },
        registered_radii=list(RADII),
        chart_costs={k: str(p) for k, p in costs.items()},
        shells=shell,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    output = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    path = HERE / "exterior_tide.json"
    if args.write:
        path.write_text(output)
    else:
        assert path.read_text() == output
    print("Checked all registered exterior chords and both ordinary calibrations.")
