"""Balanced source mixtures at identical symmetric contact probabilities."""

import argparse
from itertools import combinations
import json


from . import bridge_radar as bridge, paired_relay as pair, paired_survey as radar

PROTOCOL = bridge.HERE / "bridge_range_protocol.md"
RECEIPT = bridge.HERE / "bridge_range.json"


def differences(values):
    return [values[i + 1] - values[i] for i in (0, 2, 4, 6)]


def balance(s, a, t):
    c = [
        s[(i + 1) % 3] * a[(i + 2) % 3] - s[(i + 2) % 3] * a[(i + 1) % 3]
        for i in range(3)
    ]
    det = sum(x * y for x, y in zip(c, t, strict=True))
    norm = sum(abs(x) for x in c)
    if 0 in det or 0 in norm:
        return None
    return dict(
        signed_weights=[radar.limits(x / norm) for x in c],
        delta_range=radar.limits(det / norm),
    )


def compute(data):
    # All quantiles used below were counted in the primary simultaneous event.
    if data["results"] != bridge.decode(data):
        raise ValueError("Invalid primary measurements")
    result = {}
    for key, x in data["radar"].items():
        _, courier, _ = key.split(",")
        coin = pair.checked(data["contact"][courier], (10, 3), bridge.SHOTS)
        s, a = [differences([bridge.interval(v) for v in coin[:, i]]) for i in (0, 1)]
        x = pair.checked(x, (10, 6, 11), bridge.SHOTS)
        lengths = [bridge.ranges(row) for row in x[:8]]
        names = {}
        for i, j in radar.ROUTES:
            names[f"{i}>{j}"] = [v[i, j] for v in lengths]
            if i < j:
                names[f"{i}-{j}"] = [(v[i, j] + v[j, i]) / 2 for v in lengths]
        take = {}
        for route, values in names.items():
            t = differences(values)
            for inds in combinations(range(4), 3):
                selected = [[row[i] for i in inds] for row in (s, a, t)]
                name = route + " " + ",".join(bridge.LABELS[i] for i in inds)
                take[name] = balance(*selected)
        result[key] = take
    return result


def run():
    parent = json.loads(bridge.RECEIPT.read_text())
    return dict(
        protocol=bridge.digest(PROTOCOL),
        parent=bridge.digest(bridge.RECEIPT),
        mixtures=compute(parent),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = run()
    if args.write:
        with RECEIPT.open("x") as f:
            json.dump(data, f, separators=(",", ":"), allow_nan=False)
            f.write("\n")
    elif data != json.loads(RECEIPT.read_text()):
        raise ValueError("Mixture certificate differs")
    for arm, rows in data["mixtures"].items():
        valid = {k: v for k, v in rows.items() if v is not None}
        print(arm, len(valid), sum("ordinary" not in k for k in valid))
        if valid:
            key = max(valid, key=lambda k: min(abs(a) for a in valid[k]["delta_range"]))
            print(key, valid[key])
