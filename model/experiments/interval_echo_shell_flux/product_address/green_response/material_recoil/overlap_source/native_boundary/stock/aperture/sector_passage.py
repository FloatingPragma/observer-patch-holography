"""Exterior mean passage, beacon choice and coherent material records on the same source."""

import argparse
import json
from hashlib import sha256

import numpy as np

from . import sector_passage_chart as chart
from . import sector_passage_inversion as inversion
from . import sector_passage_lab as lab
from . import sector_passage_release as release

HERE = chart.source.balls.HERE


def range_statistics(rr, eps, exterior, ports):
    u, v = np.triu_indices(len(exterior), 1)
    i, j = exterior[u], exterior[v]
    vals = rr[i, j]

    def witness(k):
        return {"ports": [int(ports[i[k]]), int(ports[j[k]])], "value": float(vals[k])}

    return {
        "pairs": len(vals),
        "positive": int(np.sum(vals > eps)),
        "negative": int(np.sum(vals < -eps)),
        "unresolved": int(np.sum(abs(vals) <= eps)),
        "error": eps,
        "maximum": witness(np.argmax(vals)),
        "minimum": witness(np.argmin(vals)),
    }


def run():
    words, board, ports, clearance, adj, uv, loads, periods = chart.apparatus()
    configs = list(periods)
    group = np.array([configs.index(tuple(z > 0)) for z in loads])
    weights = np.bincount(group, minlength=len(configs)) / len(words)
    prepared, initial, subset, full, mats = chart.source.sources()
    ab = np.searchsorted(ports, [25, 36])
    steps, rs, hs, costs, errs, herrs = [], [], [], [], [], []
    for ticks in periods.values():
        step, _, R, err, cost = chart.law(adj, uv, ticks)
        prob, he = chart.bearing(R, ab, err)
        assert (
            np.max(
                abs((step @ prob / 12 - prob)[np.setdiff1d(np.arange(len(ports)), ab)])
            )
            < 1e-12
        )
        steps.append(step)
        rs.append(R)
        hs.append(prob)
        costs.append(cost)
        errs.append(err)
        herrs.append(he)
    rs, hs, costs = map(np.array, (rs, hs, costs))
    # A nonnegative addition to a fixed background can only increase each
    # reciprocal resistance: it reduces some link conductances. A real
    # ordinary delivery from the excluded storage rail supplies the standard.
    ordinary = np.zeros(len(configs))
    empty = loads[loads[:, 1] == 0]
    for row in empty:
        delivered = row.copy()
        delivered[1] = 1
        ordinary[configs.index(tuple(row > 0))] -= 1 / len(empty)
        ordinary[configs.index(tuple(delivered > 0))] += 1 / len(empty)
    ordR = np.einsum("i,ijk->jk", ordinary, rs)
    ordH = ordinary @ hs
    assert ordR.min() >= -2 * max(errs)
    orig_cost = weights @ costs
    exterior = np.flatnonzero(clearance >= 2)
    # Every signed initial distribution has exact sum zero. Center the
    # numerical transfer tables before bounding its certified L1 residual.
    radius = rs - rs[0]
    bearings = hs - hs[0]
    links = costs - costs[0]
    archive = json.loads((HERE / "sector_balls.json").read_text())
    held = archive["held_native"]["8"]
    observed = np.bincount(
        group, weights=held["contrast_by_terminal_word"], minlength=len(configs)
    )
    rows = {}
    response = {}
    for divisor, (d, eps) in prepared.items():
        gain = np.bincount(group, weights=d, minlength=len(configs))
        rr = np.einsum("i,ijk->jk", gain, radius)
        h = gain @ bearings
        er = float(eps * abs(radius).max() + 2 * abs(gain).sum() * max(errs) + 1e-12)
        eh = float(eps * abs(bearings).max() + 2 * abs(gain).sum() * max(herrs) + 1e-13)
        serial = chart.serial(adj, orig_cost, gain @ links)
        spans = rs[:, ab[0], ab[1]]
        pooled = ((weights + gain) * spans) @ hs / ((weights + gain) @ spans)
        baseline = (weights * spans) @ hs / (weights @ spans)
        rows[str(divisor)] = {
            "serial": serial,
            "service_phase_ordinary_resistance": range_statistics(
                rr, er, exterior, ports
            ),
            "beacon_probability_change": h.tolist(),
            "beacon_error": eh,
            "pooled_range_bearing_change": (pooled - baseline).tolist(),
            "joint_occupied_gains": gain.tolist(),
            "source_error_l1": eps,
        }
        response[divisor] = (rr, h, gain, eh)
    rr, h, gain, eh = response[8]
    launch = int(exterior[np.argmax(abs(h[exterior]))])
    predicted = range_statistics(
        rr, rows["8"]["service_phase_ordinary_resistance"]["error"], exterior, ports
    )
    pairs = np.array(
        [
            np.searchsorted(ports, predicted[key]["ports"])
            for key in ("minimum", "maximum")
        ]
    )
    other_rates = {
        str(k): [float(v[0][i, j]) for i, j in pairs] for k, v in response.items()
    }
    times = release.run(
        board,
        ports,
        clearance,
        ab,
        h,
        ordH,
        eh,
        2 * max(herrs),
        chart.source.thermal.pair.stock.kernels()[2],
    )
    acquired = lab.source_and_receiver(
        steps,
        hs[:, launch],
        gain,
        observed,
        np.array(initial, float),
        held["copies_per_start"],
        launch,
        ab,
    )
    tensor = inversion.experiment(
        ports,
        exterior,
        adj,
        configs,
        rs,
        errs,
        steps,
        group,
        weights,
        prepared,
        (words, initial, subset, full, mats),
        ordinary,
    )
    # Average the physically observed squared ranges before the public observer
    # constructs a beacon coordinate on the unsorted data.
    before = np.bincount(group[loads[:, 1] == 0], minlength=len(configs)) / len(empty)

    def public_coordinate(f):
        return chart.bearing(np.einsum("i,ijk->jk", f, rs), ab, max(errs))[0]

    public = release.run(
        board,
        ports,
        clearance,
        ab,
        public_coordinate(weights + gain) - public_coordinate(weights),
        public_coordinate(before + ordinary) - public_coordinate(before),
        4 * max(herrs) + 2 * prepared[8][1] * max(spans) / min(spans),
        4 * max(herrs),
        chart.source.thermal.pair.stock.kernels()[2],
    )
    proof = {}
    heat = chart.source.thermal.apparatus()[-1][0].diagonal(axis1=1, axis2=2).T
    # Match both local charge readings, the full contact-energy triple and
    # the final written head. Test whether other home reports are needed.
    _, tags = np.unique(
        np.column_stack((words[:, 1], heat, words[:, -1])), axis=0, return_inverse=True
    )
    v = hs[group, launch]
    candidates = []
    for tag in np.unique(tags):
        ids = np.flatnonzero(tags == tag)
        candidates.append(
            (np.ptp(v[ids]), int(ids[np.argmin(v[ids])]), int(ids[np.argmax(v[ids])]))
        )
    gap, i, j = max(candidates)
    proof.update(
        identical_home_charge_and_thermometer={
            "words": [words[i].tolist(), words[j].tolist()],
            "contact_means_times_four": heat[i].tolist(),
            "difference_in_beacon_probability": float(gap),
            "error": 2 * max(herrs),
            "ordinary_word_indices": [i, j],
        }
    )
    return {
        "sha256": {
            p: sha256((HERE / p).read_bytes()).hexdigest()
            for p in (
                "SECTOR_PASSAGE.md",
                "sector_home.json",
                "sector_balls.json",
                "sector_calorimetry.json",
            )
        },
        "ports": ports.tolist(),
        "physical_clearance": clearance.tolist(),
        "written_seams_bypassed": 3,
        "service_patterns": [list(map(int, pattern)) for pattern in configs],
        "rates": rows,
        "pinned_pair_contrasts": other_rates,
        "ordinary_fill_resistance": range_statistics(
            ordR, 2 * max(errs), exterior, ports
        ),
        "beacons": ports[ab].tolist(),
        "selected_launch": int(ports[launch]),
        "two_clock_control": proof,
        "independent_destination": acquired,
        "material_release": times,
        "public_range_release": public,
        "spatial_modular_source": tensor,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    result = run()
    path = HERE / "sector_passage.json"
    if opt.write:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        chart.source.balls.radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_text()), result
        )
    print(result["rates"]["8"]["service_phase_ordinary_resistance"])
    print(result["independent_destination"]["simultaneous_interval"])
