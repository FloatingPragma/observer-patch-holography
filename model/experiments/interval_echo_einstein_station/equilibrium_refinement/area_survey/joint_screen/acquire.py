"""Fresh same-trajectory joint hits, independent shams, source counts, and rim clocks."""

import json
from math import log

import numpy as np
from mpmath import iv

from . import measure

SHOTS = 2**60
BATCH = 2**48
L = log(2 * 4000 / 0.01)
TERMINAL_BIAS = 1e-10
CLOCK_BIAS = 1e-10


def draw(rng, p):
    p = np.asarray(p, dtype=float).copy()
    p[-1] = 1 - p[:-1].sum()
    assert p.min() >= 0
    return rng.multinomial(BATCH, p, size=SHOTS // BATCH).sum(0)


def band(x, bias=0.0):
    y = float(x)
    d = np.sqrt(2 * y * L) + 2 * L
    return iv.mpf(
        [max(0.0, (y - d) / SHOTS - bias - 2e-15), (y + d) / SHOTS + bias + 2e-15]
    )


def ends(v):
    return [
        float(np.nextafter(float(v.a), -np.inf)),
        float(np.nextafter(float(v.b), np.inf)),
    ]


def conditional(p):
    # Near an isotropic image this identity avoids cancellation in slope bounds.
    t = sum(p)
    diff = sum((p[i] - p[j]) ** 2 for i, j in ((0, 1), (0, 2), (1, 2)))
    return iv.sqrt(1 - diff / (2 * t**2))


def rims(s, p, rng):
    H = p @ s.starts
    R = H - p @ s.scores / 1024
    d = []
    counts = []
    for h, r in zip(H, R, strict=True):
        row = draw(rng, [r, h - r, 1 - h])
        counts.append(row.tolist())
        hband = band(int(row[0] + row[1]))
        rband = band(int(row[0]), h * CLOCK_BIAS)
        d.append(1024 * (1 - rband / hband))
    edges = [
        (d[measure.clock.RAYS.index((a, b))] + d[measure.clock.RAYS.index((b, a))]) / 2
        for a, b in measure.EDGES
    ]
    faces = []
    for pairs in measure.PAIRS:
        a, b, c = [edges[j] ** 2 for j in pairs]
        faces.append(iv.sqrt(4 * a * b - (a + b - c) ** 2) / 4)
    return counts, faces


def screens(s, im, p, rng):
    flux = (p[:, None] * s.occupied).T @ im.terminal
    h, joint, twin = im.face_probs(p)
    report = []
    coeff = []
    for i, face in enumerate(measure.FACES):
        two = joint[i]
        single = flux[i, 6:][list(face)] - np.array(
            [two[0] + two[1], two[0] + two[2], two[1] + two[2]]
        )
        none = h[i] - sum(two) - sum(single)
        both = h[i] ** 2
        sham = both * twin[i]
        real = draw(rng, [*two, *single, none, 1 - h[i]])
        control = draw(rng, [*sham, both - sum(sham), 1 - both])
        report.append(
            dict(two_resident=real.tolist(), two_independent=control.tolist())
        )
        c = [band(int(v), h[i] * TERMINAL_BIAS) for v in real[:3]]
        # Each conditional marginal for one randomly drawn resident has
        # sensitivity <= bias/2. Bound both marginals' product per eligible pair.
        y = [band(int(v), both * 2 * TERMINAL_BIAS) for v in control[:3]]
        sc, sy = conditional(c), conditional(y)
        coeff.append([sc, sc * sum(c), sy, sy * sum(y)])
    return report, np.array(coeff, dtype=object)


def run():
    iv.dps = 40
    s = measure.clock.survey(1)
    im = measure.instrument(s)
    u = s.uniform
    grad, _ = measure.clock.differentials(u, s.starts, s.scores)
    moves, ctrl, _ = measure.clock.survey_directions(s, grad)
    older = json.loads((measure.clock.HERE / "counts.json").read_text())
    step = older["source_probability_step"]
    old_tickets = older["held_lottery"]["source_word_tickets"]
    new = measure.dyadic_stock(s, im)
    stocks = {
        "baseline": u,
        "old_clock_lottery": np.array(old_tickets) / 2**50,
        "augmented_lottery": new / 2**50,
        "interior_reallocation": measure.interior_reallocation(s)[1]
        / (2 * len(s.words)),
    }
    for name, arms in older["counts"].items():
        for sign, info in arms.items():
            stocks[f"{name}/{sign}"] = (
                u + int(sign) * step * moves[name] + info["compensation"] * ctrl
            )
    labels = [tuple(x for x in w if x not in s.rim) for w in s.words]
    ids = {w: i for i, w in enumerate(sorted(set(labels)))}
    labels = [ids[w] for w in labels]
    rng = np.random.default_rng(1247087)
    records = {}
    obs = {}
    reference = None
    for name, p in stocks.items():
        counts, areas = rims(s, p, rng)
        joint, coeff = screens(s, im, p, rng)
        if reference is None:
            reference = areas
        images = np.r_[
            np.asarray(reference, dtype=object) @ coeff,
            np.asarray(areas, dtype=object) @ coeff,
        ]
        bulk = draw(rng, np.bincount(labels, weights=p))
        ent = -sum(band(n) * iv.log(band(n)) for n in bulk)
        records[name] = dict(
            complete_rim_clocks=counts,
            joint_faces=joint,
            bulk_word_counts=bulk.tolist(),
        )
        mu = sum(band(SHOTS - face["two_resident"][-1]) for face in joint)
        heralds = sum(band(row[0] + row[1]) for row in counts)
        cap = (heralds - 2 * mu) * iv.log(4) + (3 * mu - heralds) * iv.log(6) / 2
        obs[name] = (ent, images, cap)
    cal = {}
    for name in ("zero_to_one", "zero_to_two", "matched_energy_and_cut_load"):
        a, plus, cap_a = obs[f"{name}/1"]
        b, minus, cap_b = obs[f"{name}/-1"]
        de = (a - b) / (2 * step)
        di = (plus - minus) / (2 * step)
        if name == "zero_to_one":
            gains = [None if (v.a <= 0 <= v.b) else -de / v for v in di]
        cal[name] = dict(
            entropy=ends(de),
            independent_capacity=ends((cap_a - cap_b) / (2 * step)),
            slopes=[ends(x) for x in di],
            residuals=[
                None if g is None else ends(de + g * v)
                for g, v in zip(gains, di, strict=True)
            ],
        )
    differences = {}
    for name in ("old_clock_lottery", "augmented_lottery", "interior_reallocation"):
        a, x, cap_a = obs[name]
        b, y, cap_b = obs["baseline"]
        differences[name] = dict(
            entropy=ends(a - b),
            independent_capacity=ends(cap_a - cap_b),
            combined_finite_change=ends((a - b) + (cap_a - cap_b)),
            image=[ends(q) for q in x - y],
        )
    return dict(
        shots_per_setting=SHOTS,
        batch=BATCH,
        simultaneous_coverage=0.99,
        registered_mean_wait=measure.clock.MEAN_TIMEOUT,
        area_gains=[None if v is None else ends(v) for v in gains],
        held=cal,
        lotteries=differences,
        augmented_source_tickets=new.tolist(),
        total_tickets=2**50,
        counts=records,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
