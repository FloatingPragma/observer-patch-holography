"""Exchange two calibrated occupied residents across a remeasured material boundary."""

import argparse
import json
from hashlib import sha256

import numpy as np
from scipy.stats import beta

from . import ordinary_mass_crossing as ordinary

home = ordinary.home
HERE = ordinary.HERE


def regional_energy(stock, board, seams, word, residents, inside):
    """Eighths of an energy unit, including the half-cut share on each side."""
    diag = stock.diagonal().reshape(11, 9)
    energy = 64 * residents**2 * (6 - board.degree)
    for a, b in board.seams:
        mark = word[4 + seams.index((a, b))] if (a, b) in seams else 0
        price = diag[mark, 3 * residents[a] + residents[b]]
        energy[a] += price
        energy[b] += price
    return int(energy[inside].sum()), int(energy[~inside].sum())


def compare_all_links(graph, states):
    # Test feasibility before calling the reversible pulse sequence. A missing
    # parking neighbor gives the same unavailable appointment in both arms.
    u, v = np.nonzero(np.triu(graph.toarray(), 1))
    active = np.flatnonzero(states[0])
    occ = set(active)
    assert set(np.flatnonzero(states[1])) == occ
    ok = []
    for x, y in zip(u, v, strict=True):
        ends = sorted({x, y} & occ)
        free = [
            set(graph.indices[graph.indptr[k] : graph.indptr[k + 1]]) - occ - {x, y}
            for k in ends
        ]
        ok.append(all(free) and (len(free) < 2 or len(free[0] | free[1]) >= 2))
    ok = np.array(ok)
    measured = home.boundary.rates(graph, states[:, active], active, u[ok], v[ok])
    np.testing.assert_array_equal(measured[0], measured[1])
    return {
        "tested_links": int(sum(ok)),
        "unavailable": np.column_stack((u[~ok], v[~ok])).tolist(),
        "measured_period_counts": np.bincount(measured[0], minlength=4).tolist(),
    }


def bond_lines(raw, left, right):
    """Native spectral lines and their equal-parity weights on a sharp blank bond."""
    i, j = 3 * left + right, 3 * right + left
    assert set(np.flatnonzero(raw[:, i])) <= {i, j}
    assert set(np.flatnonzero(raw[:, j])) <= {i, j}
    assert raw[i, i] == raw[j, j] and raw[i, j] == raw[j, i]
    diagonal = raw[i, i] / 4
    splitting = abs(raw[i, j]) / 4 if i != j else 0
    return diagonal - splitting, diagonal + splitting


def read_changed_bonds(raw, board, pair, mask):
    """Estimate the in-region energy contrast from independent bond spectral records."""
    before, after = pair
    changed = before != after
    private = 8 * (6 - board.degree) * (after**2 - before**2)
    estimate = float(private[mask].sum())
    span = 0.0
    records = []
    n = 2**18
    rng = np.random.default_rng(28112080)
    for x, y in board.seams:
        weight = (int(mask[x]) + int(mask[y])) / 2
        if weight == 0 or not (changed[x] or changed[y]):
            continue
        read = []
        for sign, state in zip((-1, 1), pair, strict=True):
            lo, hi = bond_lines(raw, state[x], state[y])
            upper = int(rng.binomial(n, 0.5)) if lo != hi else 0
            estimate += sign * weight * (lo + (hi - lo) * upper / n)
            span += weight * (hi - lo)
            read.append({"lines": [float(lo), float(hi)], "counts": [n - upper, upper]})
        records.append({"edge": [int(x), int(y)], "weight": weight, "arms": read})
    # Two sharp arms per bond. All sampled spectral fractions share 99% coverage.
    width = span * np.sqrt(np.log(4 * len(records) / 0.01) / (2 * n))
    return {
        "trials_per_bond_arm": n,
        "private_energy_change": int(private[mask].sum()),
        "readings": records,
        "inside_change": estimate,
        "inside_change_interval": [estimate - width, estimate + width],
    }


def run():
    receipt = json.loads((HERE / "ordinary_mass_crossing.json").read_bytes())
    words = home.source.thermal.apparatus()[0]
    board, homes, seams, *_ = home.source.balls.energy.returned(words)
    graph, _, _ = home.source.balls.radar.chart(board, seams)
    graph = graph.tocsr()
    raw = home.source.thermal.pair.stock.kernels()[2]
    swap = home.boundary.exchange_geometry.contact_swap().argmax(0)
    initial = receipt["inventories"][0]["orders"]["ascending"]
    baseline = set(initial["markers"])
    rows = []
    calorimeter = None
    for p in initial["cut_edges"]:
        i, j = p if p[0] in baseline else p[::-1]
        if {i, j} & set(homes) or {i, j} & set(np.ravel(seams)):
            rows.append(
                {"inside_outside": [i, j], "skip": "resident or stored seam endpoint"}
            )
            continue
        allowed = words[words[:, homes.index(ordinary.TARGET)] == 0]
        row = {"inside_outside": [i, j], "backgrounds": []}
        for inv in np.unique(allowed[:, :4], axis=0):
            word = allowed[(allowed[:, :4] == inv).all(1)][0]
            before = np.zeros(board.n, dtype=int)
            before[homes] = inv
            before[[i, j]] = 1, 2
            after = before.copy()
            after[i], after[j] = divmod(int(swap[3 * before[i] + before[j]]), 3)
            assert (after[i], after[j]) == (2, 1)
            pair = np.array([before, after])
            sites = np.arange(board.n)
            dist, group, trips = home.ranging.routes(
                graph, pair, sites, ordinary.CENTER, ordinary.VOLUME
            )
            np.testing.assert_array_equal(dist[group[0]], dist[group[1]])
            result = {"load_inventory": inv.tolist(), "paths_tested": trips}
            admitted = True
            for k, order in (("ascending", sites), ("descending", sites[::-1])):
                ranked = order[np.argsort(dist[group[0], order], kind="stable")][
                    : ordinary.VOLUME
                ]
                if i not in ranked or j in ranked:
                    admitted = False
                    continue
                mask = np.isin(sites, ranked)
                energies = np.array(
                    [regional_energy(raw, board, seams, word, s, mask) for s in pair]
                )
                local, donor = np.diff(energies, axis=0)[0] / 8
                mass = 48 * (after**2 - before**2)
                assert mass[mask].sum() == 144 and mass.sum() == 0
                result[k] = {
                    "markers": np.sort(ranked).tolist(),
                    "inside_rest_gain": 144,
                    "outside_rest_gain": -144,
                    "inside_native_energy_gain": float(local),
                    "outside_native_energy_gain": float(donor),
                    "whole_board_native_energy_gain": float(local + donor),
                }
            result["both_orders"] = admitted
            receiver, donor = map(int, dist[group[0], [i, j]])
            horizon = int(np.count_nonzero(dist[group[0]] <= max(receiver, donor)))
            shell_dist, arm, trips = home.ranging.routes(
                graph, pair, sites, ordinary.CENTER, horizon
            )
            shells = shell_dist[arm] <= receiver
            np.testing.assert_array_equal(shells[0], shells[1])
            shell = {
                "receiver_range": receiver,
                "donor_range": donor,
                "within_radius": np.flatnonzero(shells[0]).tolist(),
                "paths_tested": trips,
            }
            if receiver < donor:
                mask = shells[0]
                u, v = np.nonzero(np.triu(graph.toarray(), 1))
                cut = mask[u] != mask[v]
                gains = (
                    np.diff(
                        [
                            regional_energy(raw, board, seams, word, state, mask)
                            for state in pair
                        ],
                        axis=0,
                    )[0]
                    / 8
                )
                shell.update(
                    cut_edges=np.column_stack((u[cut], v[cut])).tolist(),
                    inside_native_energy_gain=float(gains[0]),
                    outside_native_energy_gain=float(gains[1]),
                    whole_board_native_energy_gain=float(gains.sum()),
                    inside_rest_gain=int(
                        48 * (after[mask] ** 2 - before[mask] ** 2).sum()
                    ),
                    protected_timing=compare_all_links(graph, pair),
                )
                if calorimeter is None:
                    calorimeter = read_changed_bonds(raw, board, pair, mask)
                    lo, hi = calorimeter["inside_change_interval"]
                    assert lo <= gains[0] <= hi
            result["complete_shell"] = shell
            row["backgrounds"].append(result)
        rows.append(row)
    enrolled = [
        (x, z)
        for x in rows
        for z in x.get("backgrounds", [])
        if "protected_timing" in z["complete_shell"]
    ]
    if not enrolled:
        raise ValueError("No complete radar shell enrolls an original boundary edge")
    # Hold an admitted chart fixed for the two blinded arms. Read disjoint
    # standing-energy lines at both ends on independent copies of each ticket.
    rng = np.random.default_rng(28112079)
    N = 2**18
    lines = home.lab.spectral_lines(raw)
    tapes = []
    for p in (0.25, 0.5):
        flags = rng.random(N) < p
        q = np.stack((1 + flags, 2 - flags), axis=1)
        energies = lines[q, rng.integers(2, size=q.shape)]
        decode = 1 + (energies >= 32)
        assert np.array_equal(q, decode) and np.all(decode.sum(1) == 3)
        k = int(sum(flags))
        lim = [
            float(beta.ppf(0.0025, k, N - k + 1)),
            float(beta.ppf(0.9975, k + 1, N - k)),
        ]
        tapes.append(
            {
                "probability": p,
                "source_fraction_interval": lim,
                "joint_spectral_counts": np.bincount(
                    energies[:, 0] * 34 + energies[:, 1], minlength=34**2
                )
                .reshape(34, 34)[[8, 9, 32, 33]][:, [8, 9, 32, 33]]
                .tolist(),
                "agreed_swaps": k,
                "regional_entropy_interval": ordinary.entropy_band(*lim),
            }
        )

    def delta(key):
        return (np.array(tapes[1][key]) - np.array(tapes[0][key])[::-1]).tolist()

    return {
        "input_sha256": {
            p: sha256((HERE / p).read_bytes()).hexdigest()
            for p in (
                "OCCUPIED_MASS_EXCHANGE.md",
                "OCCUPIED_RADAR_SHELL.md",
                "ordinary_mass_crossing.json",
            )
        },
        "boundary_edges": rows,
        "enrolled_eleven_marker_backgrounds": sum(
            z["both_orders"] for x in rows for z in x.get("backgrounds", [])
        ),
        "enrolled_edge_backgrounds": len(enrolled),
        "regional_calorimeter": calorimeter,
        "shots_per_arm": N,
        "spectral_lines": [8, 9, 32, 33],
        "sampled_edge": enrolled[0][0]["inside_outside"],
        "sampled_background": enrolled[0][1]["load_inventory"],
        "source_arms": tapes,
        "enclosed_rest_increase_interval": (
            144 * np.array(delta("source_fraction_interval"))
        ).tolist(),
        "enclosed_record_entropy_increase_interval": delta("regional_entropy_interval"),
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    data = run()
    dest = HERE / "occupied_mass_exchange.json"
    if cli.parse_args().write:
        dest.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    else:
        ordinary.standard.exp.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(dest.read_bytes()), data
        )
    print(
        data["enrolled_edge_backgrounds"],
        data["sampled_edge"],
        data["enclosed_rest_increase_interval"],
        data["enclosed_record_entropy_increase_interval"],
    )
