"""Rebuild the local decisions and range calibration from deposited counts alone."""

from collections import defaultdict
import json
from math import log, sqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
SITES = ("82", "72", "182")
ARMS = ("passive", "relay", "sham", "half_exchange")


def decisions(row, n, deadline):
    keys = row["keys"]
    assert len(keys) == len(set(keys))
    for key in keys:
        tag, word, address = key.split(":")
        assert tag == "unmarked" and 1 <= len(word) <= 4 and set(word) <= {"0", "1"}
        if address == "timeout":
            assert len(word) == 4 and word.count("1") < deadline
        else:
            assert address in {str(i) for i in range(13)}
            assert word.endswith("1") and word.count("1") == deadline
    for split in ("train", "test"):
        assert len(row[split]) == 2
        for counts in row[split]:
            assert len(counts) == len(keys) and sum(counts) == n and min(counts) >= 0
            assert all(type(c) is int for c in counts)
    for j in range(2):
        observed = row["test"][j]
        assert row["expired"][j] == sum(
            c for k, c in zip(keys, observed) if k.endswith(":timeout")
        )
        assert row["detected"][j] == sum(
            c
            for k, c in zip(keys, observed)
            if k.rsplit(":", 1)[1] not in ("12", "timeout")
        )
    scores = {}
    for instrument in ("address", "clock", "joint"):
        groups = defaultdict(list)
        for i, key in enumerate(keys):
            _, word, address = key.split(":")
            reading = {"address": address, "clock": word, "joint": (word, address)}
            groups[(address == "timeout", reading[instrument])].append(i)
        hits = 0
        for (expired, _), indices in groups.items():
            votes = [sum(counts[i] for i in indices) for counts in row["train"]]
            guess = int(not expired and votes[1] > votes[0])
            hits += sum(row["test"][guess][i] for i in indices)
        scores[instrument] = hits
    assert scores == row["correct"]
    return scores


def best_address_count(row):
    # There are 13 nonexpired address symbols including quiet. A timeout
    # must keep the default decision. For each other symbol the challenge
    # counts choose the winning bit; a union bound below covers this choice.
    by_address = defaultdict(lambda: [0, 0])
    for j, key in enumerate(row["keys"]):
        a = key.rsplit(":", 1)[1]
        for bit in (0, 1):
            by_address[a][bit] += row["test"][bit][j]
    return sum(
        counts[0] if a == "timeout" else max(counts) for a, counts in by_address.items()
    )


def run():
    data = json.loads((HERE / "records.json").read_text())
    n = data["shots_per_bit"]
    eps = sqrt(log(2 * 72 / 0.01) / (4 * n))
    assert abs(eps - data["simultaneous_half_width"]) < 1e-14
    assert data["threshold"] == 0.52
    passive = {}
    for k, site in enumerate(SITES, 1):
        sonar = data["addressed_ranging"][site]
        good, cut = sonar["service"], sonar["outward_cut"]
        assert good["remote_mark"] and good["return_tick"] == 2 * k
        assert good["tick_word"] == [1] * (2 * k)
        assert good["local_gateway"] == [0] * (2 * k - 1) + [1]
        assert not cut["remote_mark"] and cut["return_tick"] is None
        src = data["stations"][site]
        rates = {}
        for arm in ARMS:
            bits = src[arm]["reset"]["bit_addresses"]
            assert bits == ([1, 12] if k == 1 and arm == "relay" else [12, 12])
            rates[arm] = []
            for d in range(1, 4):
                row = src[arm][str(d)]
                scores = decisions(row, n, d)
                rates[arm].append(scores["address"] / (2 * n))
        passive[site] = rates["passive"]
        passive_range = []
        for d in range(1, 4):
            r = src["passive_ranging"][str(d)]
            assert sum(r["counts"]) == n
            count = sum(
                c
                for key, c in zip(r["keys"], r["counts"])
                if key.startswith("visited:")
                and key.rsplit(":", 1)[1] not in ("12", "timeout")
            )
            assert count == r["marked_returns"]
            passive_range.append(count)
        print(
            "Station",
            site,
            "round trip",
            good["return_tick"],
            "passive returns",
            passive_range,
        )
        for arm in ARMS:
            print(" ", arm, *[f"{100 * v:.5f}%" for v in rates[arm]])
    # Additional diagnostic after protocol definition: every such table,
    # every station, flight arm and deadline in a single one-sided family.
    uniform = sqrt(log(36 * 2**13 / 0.01) / (4 * n))
    for site in SITES:
        rows = data["stations"][site]["passive"]
        ceiling = [
            min(1, best_address_count(rows[str(d)]) / (2 * n) + uniform)
            for d in (1, 2, 3)
        ]
        print(
            site,
            "simultaneous 99% ceilings for every address decoder:",
            [round(float(x), 7) for x in ceiling],
        )
    calibration = passive[SITES[0]][0] - eps > data["threshold"]
    print("Near first tick clears threshold:", calibration)
    for site, by in zip(SITES[1:], (2, 3)):
        rates = passive[site][:by]
        if calibration:
            excluded = max(rates) + eps < data["threshold"]
            passed = max(rates) - eps > data["threshold"]
            print(
                site,
                "proportional deadline forecast:",
                "excluded" if excluded else "met" if passed else "unresolved",
            )
    print(
        "Simultaneous 99% sampling half-width:", f"{100 * eps:.5f}", "percentage points"
    )


if __name__ == "__main__":
    run()
