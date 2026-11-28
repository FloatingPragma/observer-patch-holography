"""Certify the sign for continuously variable two-head entropy priority."""

import argparse
from fractions import Fraction as Q
import json
import math

import numpy as np

from . import landmark_entropy as entropy
from . import landmark_preparation as prep
from . import collective_source as source

HERE = source.tidal.HERE


def run():
    records, _ = prep.prepare()
    words, probs, _ = source.records()
    orig = probs["original"]
    _, _, cells, _, _, _, ids, _ = source.tidal.lab.energy.pair.setup()
    quadratic = np.array(
        list(map(Q, records["preparations"]["minimum_change"]["group_weights"])), object
    )[ids]
    # The least-squares first-order condition makes h a linear combination
    # of the feature rows and a constant on the free components.
    h = (quadratic - orig) / orig
    assert max(abs(h)) <= 1 and orig @ h == 0
    partial = np.isin(words[:, -1], source.PARTIAL)
    assert np.all(h[partial] == 0)
    free_mass = sum(orig[~partial])
    blank = words[:, -1] == 0
    b0 = orig @ partial
    n = (words[:, :4] == 1).sum(1)
    # Four disjoint 12-port owners per chain bound (m-1)/V by 1/48.
    # The native pooled response is affine between these two score changes.
    scores = [blank.astype(object), blank + b0 * n / (4 * 48)]
    span = [Q(1), 1 + b0 / 48]
    slope = (quadratic - orig) @ blank
    features = []
    for x, width in zip(scores, span, strict=True):
        features.append(
            (
                orig[~partial] @ x[~partial],
                {
                    int(i): sum(x[ids == i]) / Q(int(sum(ids == i)))
                    for i in np.unique(ids[~partial])
                },
                width,
            )
        )
    edge = Q(1, 1024)
    floor = Q(1, 512)
    # For small s, exp(s*h)/Z differs from 1+s*h by at most 4*s²
    # relative to conditional uniform weight and is at least 1/8.
    # Therefore conditional chi²<=128*s⁴ and TV(p*,linear)<=16F*s².
    assert all(
        (quadratic - orig) @ x - 16 * free_mass * edge * width > floor
        for x, width in zip(scores, span, strict=True)
    )
    _, left_w, _ = entropy.solve(True, edge)
    intervals = []
    left = edge
    while left < 1:
        right = min(Q(1), left * Q(17, 16))
        mid = (left + right) / 2
        _, _, detail = entropy.solve(True, mid)
        _, right_w, _ = entropy.solve(True, right)
        q, sizes, free, gids, total = detail
        fids = np.flatnonzero(free)
        margins = []
        for baseline, average, width in features:
            lower = total * sum(
                lo * average[i] for i, (lo, _) in zip(fids, q, strict=True)
            )
            margins.append((lower - baseline - floor * right) / width)
        margin = min(margins)
        worst = Q(0)
        for p in (left_w, right_w):
            group = np.bincount(gids)  # group cardinalities checked below
            chi = Q(0)
            for i, size, (lo, hi) in zip(fids, sizes, q, strict=True):
                assert group[i] == size and lo > 0
                mass = sum(p[gids == i]) / total
                chi += max(abs(mass - lo), abs(mass - hi)) ** 2 / lo
            worst = max(worst, chi)
        # Chi² is convex along the affine feasible witness joining the
        # endpoints. One rational comparison covers every intermediate s.
        assert margin > 0 and margin**2 > total**2 * worst / 2, (
            float(left),
            float(right),
        )
        intervals.append(
            dict(
                left=str(left),
                right=str(right),
                squared_margin_over_needed=float(2 * margin**2 / (total**2 * worst)),
            )
        )
        left, left_w = right, right_w
    joint = np.array([orig @ (cells == j) for j in range(4)], object).reshape(2, 2)
    prod = np.outer(joint.sum(0), joint.sum(1)).T
    target = np.array(
        [probs["source"] @ (cells == j) for j in range(4)], object
    ).reshape(2, 2)
    assert np.array_equal(prod, target)

    def shannon(p):
        return -sum(float(v) * math.log(float(v)) for v in p.ravel() if v)

    return dict(
        initial_joint=list(map(str, joint.ravel())),
        independent_joint=list(map(str, target.ravel())),
        two_head_entropy_gain_nats=shannon(prod) - shannon(joint),
        initial_word_entropy_nats=math.log(len(words)),
        repaired_word_entropy_nats=shannon(left_w),
        small_interval_upper=str(edge),
        linear_blank_slope=str(slope),
        source_gain_lower_per_fraction=str(floor),
        interval_count=len(intervals),
        intervals=intervals,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    raw = json.dumps(out, indent=2, allow_nan=False) + "\n"
    path = HERE / "landmark_strength.json"
    if args.write:
        path.write_text(raw)
    else:
        assert path.read_text() == raw
    print(
        out["interval_count"],
        min(row["squared_margin_over_needed"] for row in out["intervals"]),
    )
