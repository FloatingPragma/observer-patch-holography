"""Move calibrated ordinary loads into the same measured volume and compare the records."""

import argparse
import json
from hashlib import sha256
from itertools import pairwise
from math import log

import numpy as np
from scipy.special import entr
from scipy.stats import beta

from . import sector_home as home
from . import thermal_standard as standard

HERE = home.HERE
CENTER, TARGET, VOLUME = 13, 20, 11


def entropy_band(lo, hi):
    def h(p):
        return float(entr(p) + entr(1 - p))

    return [min(h(lo), h(hi)), h(np.clip(0.5, lo, hi))]


def energy_calibration(board, raw):
    adjacency = board.adjacency().toarray()
    L = np.diag(adjacency.sum(1)) - adjacency
    lam = np.linalg.eigvalsh(L)[1]
    dt = 1 / lam
    out = {}
    for q in (1, 2):
        g = standard.standing_generator(board, raw, q)
        # Extract the intercept directly from the uniform-mode action.
        gap = int(g[0].sum()) // 4
        np.testing.assert_array_equal(g - 4 * gap * np.eye(len(g), dtype=int), 2 * L)
        omega = lam / 2
        t = 2**24
        counts = np.random.default_rng(987013 + q).binomial(
            t, (1 + np.array([np.cos(omega * dt), np.sin(omega * dt)])) / 2
        )
        xy = 2 * counts / t - 1
        angle = np.arctan2(xy[1], xy[0])
        width = np.arcsin(2 * np.sqrt(log(800) / (t)) / np.linalg.norm(xy))
        slope = (angle + np.array([-width, width])) / (dt * lam)
        out[str(q)] = {
            "standing_gap_native": gap,
            "first_laplacian_eigenvalue": float(lam),
            "time": dt,
            "mode_phase_counts": counts.tolist(),
            "trials_per_quadrature": t,
            "kinetic_mass_at_p_squared_equal_lambda": (1 / (2 * slope[::-1])).tolist(),
        }
    return out


def run():
    words = home.source.thermal.apparatus()[0]
    board, homes, seams, *_ = home.source.balls.energy.returned(words)
    graph, hops, _ = home.source.balls.radar.chart(board, seams)
    ticks, check = home.deposit(graph, hops, words, homes)
    raw = home.source.thermal.pair.stock.kernels()[2]
    spectral = home.lab.spectral_lines(raw)
    scale = energy_calibration(board, raw)
    path = check["path"]
    donor = path[0]
    support = homes + [donor]
    h = homes.index(TARGET)
    residents = np.unique(words[words[:, h] == 0, :4], axis=0)
    adj = graph.toarray().astype(bool)
    u, v = np.nonzero(np.triu(adj, 1))
    lookup = home.boundary.exchange_geometry.contact_swap().argmax(0)
    inventory_rows = []
    for inv in residents:
        end = []
        during = []
        for q in (0, 1, 2):
            same = np.zeros(board.n, dtype=int)
            same[homes] = inv
            same[donor] = q
            history = [same.copy()]
            for i, j in pairwise(path):
                assert adj[i, j] and same[i] == q and same[j] == 0
                same[i], same[j] = divmod(int(lookup[3 * same[i] + same[j]]), 3)
                history.append(same.copy())
            hist = np.array(history)
            assert np.all((hist**2).sum(1) == (hist[0] ** 2).sum())
            end.extend([hist[0, support], hist[-1, support]])
            during.append(hist)
        end = np.array(end)
        dist, ids, trips = home.ranging.routes(graph, end, support, CENTER, VOLUME)
        sites = np.arange(board.n)
        readings = [
            home.ranging.routes(graph, hist, sites, CENTER, VOLUME) for hist in during
        ]
        sided = {}
        for name, order in (
            ("ascending", np.arange(board.n)),
            ("descending", np.arange(board.n)[::-1]),
        ):
            rank = order[np.argsort(dist[:, order], axis=1, kind="stable")][ids]
            region = np.sort(rank[:, :VOLUME], axis=1)
            np.testing.assert_array_equal(region, np.tile(region[0], (len(end), 1)))
            member = region[0]
            mask = np.isin(np.arange(board.n), member)
            assert mask[TARGET] and not mask[donor]
            boundary = mask[u] != mask[v]
            periods = home.boundary.rates(graph, end, support, u[boundary], v[boundary])
            np.testing.assert_array_equal(periods, np.tile(periods[0], (len(end), 1)))
            enrolled = []
            passage = []
            for q, hist in enumerate(during):
                ds, group, _ = readings[q]
                balls = order[np.argsort(ds[:, order], axis=1, kind="stable")][
                    group, :VOLUME
                ]
                steps = []
                for place, sharp, B in zip(path, hist, balls, strict=True):
                    active = np.isin(sites, B)
                    cut = active[u] != active[v]
                    duration = home.boundary.rates(
                        graph, [sharp], sites, u[cut], v[cut]
                    )[0]
                    steps.append(
                        {
                            "port": place,
                            "cut_count": int(sum(cut)),
                            "cut_service_per_slot": float(np.sum(1 / duration)),
                            "keeps_initial_markers": bool(
                                np.array_equal(np.sort(B), member)
                            ),
                            "rest_inside": int(48 * (sharp[active] ** 2).sum()),
                        }
                    )
                passage.append(steps)
                energy = 48 * (hist[:, mask] ** 2).sum(1)
                assert energy[-1] - energy[0] == 48 * q**2
                events = np.flatnonzero(np.diff(energy))
                enrolled.append(
                    [
                        [int(path[k]), int(path[k + 1]), int(energy[k + 1] - energy[k])]
                        for k in events
                    ]
                )
            sided[name] = {
                "markers": member.tolist(),
                "cut_edges": np.column_stack((u[boundary], v[boundary])).tolist(),
                "protected_cut_periods": periods[0].tolist(),
                "crossings_by_load": enrolled,
                "target_out_and_back_by_arm": (2 * dist[ids, TARGET]).tolist(),
                "moving_surveys": passage,
            }
        inventory_rows.append(
            {
                "fixed_homes": inv.tolist(),
                "returned_ranging_probes": trips,
                "orders": sided,
            }
        )
    # Two independent spectral readings on replicas sharing the actual
    # donor ticket. Each reads the orthogonal load class with disjoint lines.
    # The ticket is withheld from the two local observers.
    N = 2**18
    rng = np.random.default_rng(886154)
    records = {}
    for label, p, sharp in (
        ("blank", 0, 0),
        ("known_single", 0, 1),
        ("quarter_double", 0.25, None),
        ("half_double", 0.5, None),
        ("known_double", 0, 2),
    ):
        q = (
            np.full(N, sharp, dtype=int)
            if sharp is not None
            else 2 * rng.binomial(1, p, N)
        )
        scores = spectral[q[:, None], rng.integers(2, size=(N, 2))]
        reads = (scores > 0).astype(int) + (scores >= 32).astype(int)
        np.testing.assert_array_equal(reads, np.column_stack((q, q)))
        bank = np.bincount((reads[:, 0] * 3 + reads[:, 1]), minlength=9).reshape(3, 3)
        row = {
            "agreed_load_counts": bank.tolist(),
            "spectral_counts": [
                np.bincount(col, minlength=34).tolist() for col in scores.T
            ],
            "echo_counts": [int(np.count_nonzero(ticks[q] == k)) for k in (2, 4)],
            "reconstructed_mean_rest_energy": float(np.mean(48 * reads[:, 0] ** 2)),
        }
        if sharp is None:
            k = int(bank[2, 2])
            # Two two-sided exact binomial intervals and a union bound.
            lo = float(beta.ppf(0.0025, k, N - k + 1))
            hi = float(beta.ppf(0.9975, k + 1, N - k))
            row["fraction_interval"] = [lo, hi]
            row["entropy_interval"] = entropy_band(lo, hi)
        records[label] = row
    lo, hi = [records[name] for name in ("quarter_double", "half_double")]
    return {
        "inputs": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in (
                "ORDINARY_MASS_CROSSING.md",
                "sector_home.json",
                "thermal_standard.json",
            )
        },
        "rest_and_motion": scale,
        "route": path,
        "inventories": inventory_rows,
        "shots_per_arm": N,
        "readings": records,
        "half_minus_quarter_entropy": [
            hi["entropy_interval"][0] - lo["entropy_interval"][1],
            hi["entropy_interval"][1] - lo["entropy_interval"][0],
        ],
        "half_minus_quarter_energy": [
            192 * (hi["fraction_interval"][0] - lo["fraction_interval"][1]),
            192 * (hi["fraction_interval"][1] - lo["fraction_interval"][0]),
        ],
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    dst = HERE / "ordinary_mass_crossing.json"
    result = run()
    if cli.parse_args().write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        standard.exp.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(dst.read_bytes()), result
        )
    print(result["half_minus_quarter_energy"], result["half_minus_quarter_entropy"])
