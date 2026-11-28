"""Ordinary mass deposited inside a surveyed volume and the same two source readings."""

import argparse
import itertools
import json
from fractions import Fraction
from hashlib import sha256

import numpy as np
from scipy.sparse.csgraph import shortest_path

from . import sector_ball_source as source
from . import sector_boundary_clock as boundary
from . import sector_cadence, sector_home_charge
from . import sector_home_lab as lab
from . import sector_loaded_ball as ranging

HERE = source.balls.HERE


def deposit(graph, hops, words, homes):
    center, target = 13, 20
    donor = int(np.argmax(hops[center]))
    adjacency = graph.toarray().astype(bool)
    usable = adjacency.copy()
    other = [i for i in homes if i != target]
    usable[other] = False
    usable[:, other] = False
    distance, parent = shortest_path(
        usable, indices=donor, unweighted=True, return_predecessors=True
    )
    assert np.isfinite(distance[target])
    path = [target]
    while path[-1] != donor:
        path.append(int(parent[path[-1]]))
    path.reverse()
    swap = boundary.exchange_geometry.contact_swap().argmax(0)
    footprint = homes + [donor]
    admitted = words[words[:, homes.index(target)] == 0]
    inventories = np.unique(admitted[:, :4], axis=0)
    records = []
    for inventory in inventories:
        for load in (1, 2):
            ports = np.zeros(len(adjacency), dtype=int)
            ports[homes] = inventory
            ports[donor] = load
            initial = ports.copy()
            for u, v in itertools.pairwise(path):
                assert adjacency[u, v] and ports[u] == load and ports[v] == 0
                ports[u], ports[v] = divmod(int(swap[3 * ports[u] + ports[v]]), 3)
            assert ports[target] == load and np.array_equal(
                ports[other], initial[other]
            )
            assert np.array_equal(np.sort(ports), np.sort(initial))
            records.extend([initial[footprint], ports[footprint]])
    records = np.array(records)
    dist, idx, trips = ranging.routes(graph, records, footprint, center, 11)
    # Public marker memberships and occupied masses are independently recovered
    # on both sides of the delivery. Check both conventions for tied ranges.
    witnesses = {}
    for name, order in (
        ("ascending", np.arange(len(adjacency))),
        ("descending", np.arange(len(adjacency))[::-1]),
    ):
        rank = order[np.argsort(dist[:, order], axis=1, kind="stable")][idx]
        inside = np.argsort(rank, axis=1)[:, footprint] < 11
        assert inside[:, homes.index(target)].all() and not inside[:, -1].any()
        rest = (48 * records**2 * inside).sum(1)
        delta_rest = rest[1::2] - rest[::2]
        assert np.array_equal(delta_rest, 48 * records[1::2, homes.index(target)] ** 2)
        cuts = []
        u, v = np.nonzero(np.triu(adjacency, 1))
        for region in rank[:, :11]:
            mark = np.isin(np.arange(len(adjacency)), region)
            cuts.append(int(np.count_nonzero(mark[u] != mark[v])))
        witnesses[name] = {
            "enclosed_rest_changes": sorted(set(map(int, delta_rest))),
            "cut_changes": sorted(
                set(map(int, np.diff(np.array(cuts).reshape(-1, 2), axis=1).ravel()))
            ),
            "identical_memberships": bool(
                (np.sort(rank[::2, :11]) == np.sort(rank[1::2, :11])).all()
            ),
        }
    durations = boundary.rates(graph, records, footprint, [center], [target])[:, 0] * 2
    ticks = np.array(
        [
            np.unique(durations[records[:, homes.index(target)] == k]).item()
            for k in range(3)
        ]
    )
    assert np.array_equal(ticks, [2, 4, 4])
    return ticks, {
        "donor": donor,
        "path": path,
        "empty_input_words": len(admitted),
        "tested_inventories": records[::2].tolist(),
        "ordinary_round_trips": trips,
        "roundtrip_ticks": ticks.tolist(),
        "surveyed": witnesses,
    }


def rate_bounds(
    values,
    delta,
    subset,
    mats,
    names=("home_single_units", "four_home_single_units", "twice_bond_energy"),
):
    history = json.loads((HERE / "sector_retention.json").read_text())
    clock = source.experiment.ReturnClock(
        (mats[2] @ mats[1] @ mats[0]).T.tocsr(), subset
    )
    power, eta = sector_cadence.mixing(clock, len(delta))
    ticks = clock.grid * values
    contrast = np.asarray(delta, float)
    norm = np.abs(contrast).sum()
    err = np.zeros(values.shape[1])
    intervals = []
    for _ in range(history["last_age"]):
        ticks, step = clock.advance(ticks)
        err += step
        val = ticks / clock.grid
        mean = contrast @ val
        width = norm * err + 1e-10
        intervals.append(np.stack((mean - width, mean + width), axis=1))
    data = np.array(intervals)
    tail = ((np.ptp(val, axis=0) + 2 * err) * norm / 2 + 1e-10) * power / (
        1 - eta
    ) + 1e-9
    result = {}
    assert len(names) == values.shape[1]
    for i, name in enumerate(names):
        lo, hi = data[:, i].T
        result[name] = []
        for row in history["continuous_cadence"]["rate_bins"]:
            a, b = row["rate"]
            bound = sector_cadence.evaluate(1 - b, 1 - a, lo, hi, tail[i])
            result[name].append({"rate": [a, b], "per_rate": list(bound)})
    return result


def run():
    words = source.thermal.apparatus()[0]
    board, homes, seams, terms, _, total = source.balls.energy.returned(words)
    graph, hops, _ = source.balls.radar.chart(board, seams)
    ticks, calibration = deposit(graph, hops, words, homes)
    four_h = source.thermal.pair.stock.kernels()[2]
    lines = lab.spectral_lines(four_h)
    home = words[:, homes.index(20)]
    bond = terms[board.seams.index((13, 20))] / 4
    assert np.array_equal(lines.mean(1)[home], bond)
    single_units = words[:, :4] ** 2
    inventory = single_units.sum(1)
    counts = np.bincount(inventory)
    assert np.all(words[:, :4].sum(1) == 4)
    reference = Fraction(int(inventory.sum()), len(inventory))
    # The mean physical energy of every supported preparation is the
    # expectation of this compression, including coherences. The matrix
    # need not keep the support invariant under subsequent evolution.
    h4 = source.balls.energy.compressed_hamiltonian(words, board, homes, seams, total)
    upper = int((h4.diagonal() + np.abs(h4).sum(1) - np.abs(h4.diagonal())).max())
    assert upper == int(total.max())
    ordinary_energy = Fraction(int(total.sum()), 4 * len(words))
    energy_gain = Fraction(upper, 4) - ordinary_energy
    assert np.array_equal(6 * bond - 1.5 * (ticks[home] - 2), 48 * home**2)
    std = json.loads((HERE / "thermal_standard.json").read_text())
    assert [std[k]["uniform_energy"] for k in ("single", "double")] == [48, 192]
    prepared, delta, subset, full, mats = source.sources()
    rates = rate_bounds(
        np.column_stack((home**2, single_units.sum(1), (2 * bond).astype(np.int64))),
        delta,
        subset,
        mats,
    )
    lower = min(r["per_rate"][0] for r in rates["home_single_units"])
    assert lower > 0 and min(r["per_rate"][0] for r in rates["twice_bond_energy"]) > 0
    observables = np.column_stack((bond, ticks[home], 48 * single_units, total / 4))
    names = ["bond", "ticks", *[f"rest_{j}" for j in homes], "whole_board_energy"]
    forecast = {}
    for divisor, (d, eps) in prepared.items():
        centered = observables - (observables.min(0) + observables.max(0)) / 2
        means = d @ centered
        width = eps * abs(centered).max(0) + 1e-10
        forecast[str(divisor)] = dict(
            zip(
                names,
                np.column_stack((means - width, means + width)).tolist(),
                strict=True,
            )
        )
    measurements = sector_home_charge.reconstruct(
        json.loads((HERE / "sector_balls.json").read_text())
    )
    native = np.full(len(full), -1)
    native[subset] = np.arange(len(words))
    tables = source.experiment.acquisition.returns.ticket_tables(mats)
    held = lab.acquire(tables, subset, native, home, lines, ticks)
    return {
        "sha256": {
            p: sha256((HERE / p).read_bytes()).hexdigest()
            for p in (
                "SECTOR_HOME.md",
                "sector_retention.json",
                "sector_balls.json",
                "thermal_standard.json",
            )
        },
        "calibration": calibration,
        "modular_rest_inference": measurements,
        "clean_bond_spectra": lines.tolist(),
        "ordinary_rest_budget": {
            "reference_counts_in_single_units": counts.tolist(),
            "reference_mean_single_units": str(reference),
            "maximum_single_units_any_preparation": int(inventory.max()),
            "gain_per_reference_upper": str(int(inventory.max()) / reference - 1),
            "whole_board_reference_energy": str(ordinary_energy),
            "maximum_mean_board_energy": str(Fraction(upper, 4)),
            "sharp_word_attaining_mean_energy_bound": words[
                int(total.argmax())
            ].tolist(),
            "board_energy_gain_per_ordinary_rest_upper": str(
                energy_gain / (48 * reference)
            ),
        },
        "all_rates": rates,
        "minimum_home_rest_energy_per_reporting_probability": 48 * lower,
        "native_contrasts": forecast,
        "held_at_one_eighth": held,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    dst = HERE / "sector_home.json"
    if args.write:
        dst.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        source.balls.radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(dst.read_text()), out
        )
    print(
        out["calibration"],
        out["minimum_home_rest_energy_per_reporting_probability"],
        out["held_at_one_eighth"]["simultaneous_99_percent"],
    )
