"""Compare renewal time and full-generator release on every clean ordinary edge."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

import numpy as np

from . import clock_release as release

lab = release.lab
HERE = release.HERE


def run():
    _, w, cells, *_ = lab.energy.pair.setup()
    h = lab.energy.pair.stock.kernels()[2]
    board, home, marks, *_ = lab.energy.returned(w)
    _, dist, _ = release.thermal_radar.chart(board, marks)
    clean = set(range(board.n)) - set(home) - set(np.ravel(marks))
    edges = [e for e in board.seams if set(e) <= clean]
    shift = (lab.energy.least_change(cells) - 1) / len(w)
    assert np.all(shift @ w[:, :4] == 0)
    donors = {
        x: [j for j, v in enumerate(home) if v in board.neighbours[x]] for x in clean
    }
    assert max(map(len, donors.values())) == 1
    original = np.bincount(cells)
    assert sum(original) == len(w)
    means = [shift @ ((w[:, j] - 1) ** 2) for j in range(4)]
    ratios = {}
    checked = 0
    for t in (1, 2):
        delays = {}
        prices = {}
        for x in clean:
            nbr = board.neighbours[x]
            if not donors[x]:
                delays[x] = 0.0
                prices[x] = F()
                continue
            j = donors[x][0]
            # The addressed rail can be stored beyond this star before the
            # exposure. All probe contacts and both storage moves are blank.
            pairs = [
                (u, v)
                for u in nbr
                if u in clean
                for v in board.neighbours[u]
                if v in clean and v not in nbr and v != x
            ]
            assert pairs
            src = home[j]
            _, bank = release.clock.read_clock(h, nbr, t, source=src)
            periods = 1 / bank[w[:, j], 0]
            dt = float(np.asarray(shift, float) @ periods)
            delays[x] = float(np.log1p(dt / periods.mean()))
            prices[x] = F(3, 4) * means[j] if t == 1 else F()
        counter = dict(both_silent=0, clock_only=0, release_only=0, opposite=0, same=0)
        data = []
        for x, y in edges:
            # Enumerate all distinct source restrictions in the double
            # commutator's support. Outside generators commute at this order.
            inner = {x, y, *board.neighbours[x], *board.neighbours[y]}
            outer = inner.union(*(board.neighbours[v] for v in inner))
            keep = [v in outer for v in home] + [
                bool(inner.intersection(e)) for e in marks
            ]
            words, indices = np.unique(
                w * np.asarray(keep), axis=0, return_inverse=True
            )
            _, _, col = release.instrument(board, home, marks, h, (x, y), dist[home[1]])
            delta = F()
            for i, word in enumerate(words):
                p = sum(shift[indices == i])
                if not p:
                    continue
                actual = [col(word, x, t), col(word, y, t)]
                background = release.acceleration([col(word, None, t)], [1])[0]
                value = (release.acceleration(actual, [1, 1])[0] - background) / t
                delta += p * value
                checked += 1
            assert delta == prices[y] - prices[x]
            clock = delays[x] - delays[y]
            category = (
                ("opposite" if float(delta) * clock < 0 else "same")
                if delta and abs(clock) > 1e-12
                else (
                    "release_only"
                    if delta
                    else "clock_only"
                    if abs(clock) > 1e-12
                    else "both_silent"
                )
            )
            counter[category] += 1
            if category != "both_silent":
                data.append(
                    dict(
                        ports=[int(x), int(y)],
                        extra_log_period_contrast=clock,
                        extra_inward_release_acceleration=str(delta),
                    )
                )
        ratios[str(t)] = dict(classification=counter, nonzero_edges=data)
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "clock_release.json", HERE / "renewal_clock.json")
        },
        clean_ports=len(clean),
        edges=len(edges),
        nonzero_weight_local_restrictions_checked=checked,
        readers=ratios,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "clock_chart.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        lab.check(json.loads(dst.read_bytes()), result)
    print(result["edges"], [r["classification"] for r in result["readers"].values()])
