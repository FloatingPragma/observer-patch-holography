"""Score the complete records without the flight model or a fitted receiver."""

from itertools import product
import json
from math import log, sqrt
from pathlib import Path

from scipy.stats import beta

HERE = Path(__file__).parent
BUDGET = 2**30
# Three services (signal, echo, no encoding), both laws, both input labels,
# three distances and four deadlines. Clopper-Pearson with two-sided union bound.
TAIL = 0.01 / (2 * 3 * 2 * 2 * 3 * 4)


def interval(x):
    lo = beta.ppf(TAIL, x, BUDGET - x + 1) if x else 0
    hi = beta.isf(TAIL, x + 1, BUDGET - x) if x < BUDGET else 1
    return float(lo), float(hi)


def score_bounds(zeros, ones):
    return ((1 + zeros[0] - ones[1]) / 2, (1 + zeros[1] - ones[0]) / 2)


def keys(n, echo):
    addresses = {f"local/{v}" for v in range(12)} | {"quiet"}
    expires = {f"timeout/{v}" for v in range(1, n + 1)}
    if not echo:
        return addresses | expires
    return {f"outward/{k}" for k in expires} | {
        f"{tag}/{k}" for tag in ("marked", "unmarked") for k in addresses | expires
    }


def evaluate(data):
    assert data["shots"] == BUDGET and data["port_pairs"] == [
        [82, 74],
        [72, 77],
        [182, 184],
    ]
    assert set(data["banks"]) == {"native", "half_exchange"}
    rows = []
    for law, near, n in product(
        ("native", "half_exchange"), (82, 72, 182), range(1, 5)
    ):
        bank = data["banks"][law][str(near)][str(n)]
        assert set(bank) == {"one_way", "sham", "echo"}
        for name, pair in bank.items():
            assert len(pair) == 2
            for row in pair:
                assert (
                    set(row) == keys(n, name == "echo") and sum(row.values()) == BUDGET
                )
                assert all(type(v) is int and v >= 0 for v in row.values())
        zero, one = [x["local/1"] for x in bank["one_way"]]
        e0, e1 = [x["marked/local/1"] for x in bank["echo"]]
        sham0, sham1 = [x["local/1"] for x in bank["sham"]]
        actual = score_bounds(interval(zero), interval(one))
        forecast = score_bounds(
            [sqrt(v) for v in interval(e0)], [sqrt(v) for v in interval(e1)]
        )
        sham = score_bounds(interval(sham0), interval(sham1))
        data0, data1 = bank["one_way"]
        empirical_optimum = sum(max(v, data1[k]) for k, v in data0.items()) / (
            2 * BUDGET
        )
        # A two-sided family for all 24 settings and at most 2**17 possible
        # complete-record decision tables, including the first failed tick.
        uniform = sqrt(log(2 * 24 * 2**17 / 0.01) / (4 * BUDGET))
        rows.append(
            dict(
                law=law,
                port=near,
                deadline=n,
                hits=(zero, one),
                echoes=(e0, e1),
                correct=BUDGET + zero - one,
                optimum_interval=(
                    empirical_optimum - uniform,
                    empirical_optimum + uniform,
                ),
                estimate=(1 + (zero - one) / BUDGET) / 2,
                actual=actual,
                forecast=forecast,
                sham=sham,
                far_marked=[
                    sum(v for k, v in r.items() if k.startswith("marked/"))
                    for r in bank["echo"]
                ],
                unmarked_return=[r["unmarked/local/1"] for r in bank["echo"]],
            )
        )
    return rows


if __name__ == "__main__":
    rows = evaluate(json.loads((HERE / "records.json").read_text()))
    for row in rows:
        a, e = row["actual"], row["forecast"]
        verdict = (
            "excluded"
            if max(a[0], e[0]) > min(a[1], e[1])
            else "inconclusive"
            if max(a[1] - a[0], e[1] - e[0]) > 0.01
            else "compatible"
        )
        print(
            row["law"],
            row["port"],
            row["deadline"],
            row["echoes"],
            f"direct {row['estimate']:.7f} [{a[0]:.7f},{a[1]:.7f}]",
            f"echo [{e[0]:.7f},{e[1]:.7f}]",
            verdict,
        )
    print(
        "All no-encoding intervals cover chance:",
        all(r["sham"][0] <= 0.5 <= r["sham"][1] for r in rows),
    )

    for law in ("native", "half_exchange"):
        subset = [r for r in rows if r["law"] == law]
        best = max(subset, key=lambda r: r["optimum_interval"][1])
        print(
            law,
            "highest complete-record optimum, joint 99% interval",
            best["port"],
            best["deadline"],
            best["optimum_interval"],
        )
