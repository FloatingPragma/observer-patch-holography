"""Complete the directed-bridge census and compare independently commissioned source heads."""

import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json

import numpy as np

from . import measure as m


def short_time_counts(board):
    t = F(1, 2**20)
    a = F(4, 3)
    remainder = 32**4 * t / 24
    hi = (a + remainder) ** 2
    lo = (a - remainder) ** 2 - (32 * t) ** 2 * hi

    def bounds(c):
        same, opp = c["alike"], c["unlike"]
        return [
            str((same * lo - opp * hi) / (same + opp)),
            str((same * hi - opp * lo) / (same + opp)),
        ]

    counts = m.geometry_counts(board)
    inherited = Counter(
        m.sign_class(m.initial_index(r, -1)) for r in m.routes(board, False)
    )
    return dict(
        time_max=str(t),
        leading_flip_coefficient=str(a * a),
        covariance_over_t6=[str(lo), str(hi)],
        all_arms_mean_over_t6=bounds(counts),
        minimum_arm_mean_over_t6=bounds(inherited),
    )


def summarize(rows, population, upper):
    tables = np.asarray(rows)
    n = len(rows)
    radius = (
        0
        if n == population
        else float(
            0.5 * np.sqrt((1 - (n - 1) / population) * np.log(2 * 120 / 0.01) / (2 * n))
        )
    )
    result = []
    for i in range(len(m.TIMES)):
        c = m.covariance(tables[:, i])
        result.append(
            dict(
                n=n,
                positive=int((c > upper).sum()),
                negative=int((c < -upper).sum()),
                unresolved=int((abs(c) <= upper).sum()),
                minimum=float(c.min()),
                maximum=float(c.max()),
                mean=float(c.mean()),
                pooled=float(m.covariance(tables[:, i].mean(axis=0))),
                mean_radius_99=upper + radius,
            )
        )
    return result


def study():
    first_path = m.HERE / "results.json"
    pilot = json.loads(first_path.read_text())
    upper = pilot["error"]["covariance"]
    boards = {}
    for level in range(3):
        board = m.source.federation(level)
        minimum = list(m.routes(board, False))
        arms = list(m.routes(board, True))
        old = pilot["boards"][str(board.cells)]
        if level < 2:
            standards = range(len(minimum))
        else:
            standards = old["sample"]
        output = {}
        for q in (-1, 1):
            previous = {row["index"]: row["p"] for row in old["records"][str(q)]}
            census = []
            common = []
            for i, route in enumerate(minimum):
                start = m.initial_index(route, q)
                h = None
                if i in previous:
                    unchanged = previous[i]
                else:
                    h = m.scaled_h(board, route, q)
                    unchanged = m.propagate_scaled(h, route, start).tolist()
                census.append(unchanged)
                if i in standards:
                    if h is None:
                        h = m.scaled_h(board, route, q)
                    common.append(
                        [
                            unchanged
                            if j == start
                            else m.propagate_scaled(h, route, j).tolist()
                            for j in range(4)
                        ]
                    )
            survey = (
                list(range(len(arms)))
                if level == 0
                else sorted(
                    map(
                        int,
                        np.random.Generator(np.random.PCG64(1828096 + level)).choice(
                            len(arms), 240, replace=False
                        ),
                    )
                )
            )
            reuse = {r: p for r, p in zip(minimum, census, strict=True)}
            arm_p = [
                reuse[arms[i]]
                if arms[i] in reuse
                else m.propagate(board, arms[i], q).tolist()
                for i in survey
            ]
            arms_class = [m.sign_class(m.initial_index(arms[i], q)) for i in survey]
            cl = [m.sign_class(m.initial_index(r, q)) for r in minimum]
            common_arr = np.array(common)
            summary = dict(
                original=summarize(census, len(minimum), upper),
                all_arms=summarize(arm_p, len(arms), upper),
                standards={
                    str(j): summarize(common_arr[:, j], len(minimum), upper)
                    for j in range(4)
                },
                independent_fair=summarize(
                    common_arr.mean(axis=1), len(minimum), upper
                ),
                sharp_classes={
                    k: summarize(
                        [x for x, y in zip(census, cl, strict=True) if y == k],
                        cl.count(k),
                        upper,
                    )
                    for k in ("alike", "unlike")
                },
                arm_sign_failures=[
                    sum(
                        m.covariance(np.asarray(x)[t]) * (1 if k == "alike" else -1)
                        <= upper
                        for x, k in zip(arm_p, arms_class, strict=True)
                    ).item()
                    for t in range(len(m.TIMES))
                ],
            )
            output[str(q)] = dict(
                census=census, standards=common, survey_p=arm_p, summary=summary
            )
        boards[str(board.cells)] = dict(
            minimum_routes=minimum,
            standard_indices=list(standards),
            survey_indices=survey,
            survey_routes=[arms[i] for i in survey],
            short_time=short_time_counts(board),
            charges=output,
        )
    return dict(
        first_archive_sha256=sha256(first_path.read_bytes()).hexdigest(),
        followup_sha256=sha256((m.HERE / "FOLLOWUP.md").read_bytes()).hexdigest(),
        boards=boards,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = study()
    path = m.HERE / "extended.json"
    if args.write:
        path.write_text(
            json.dumps(result, separators=(",", ":"), allow_nan=False) + "\n"
        )
    else:
        m.compare(json.loads(path.read_text()), result)
    for size, b in result["boards"].items():
        print(
            size,
            {
                q: {
                    k: s["summary"][k][-1]["mean"]
                    for k in ("original", "all_arms", "independent_fair")
                }
                for q, s in b["charges"].items()
            },
            flush=True,
        )
