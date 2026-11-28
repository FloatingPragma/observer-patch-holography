"""Test every radial embedding of the measured compact-mode profile on the outer SPARC rings."""

from fractions import Fraction as F
from functools import cmp_to_key
from hashlib import sha256
import json
from os.path import relpath

from .polarization_sky import HERE, experiment, sky

# Squared chords of the five noncentral shells on the unit dodecahedron.
CHORD2 = [
    (F(2), F(-2, 3)),
    (F(4, 3), F(0)),
    (F(8, 3), F(0)),
    (F(2), F(2, 3)),
    (F(4), F(0)),
]


def compare(u, v):
    a, b = u[0] - v[0], u[1] - v[1]

    def sign(value):
        return int(value > 0) - int(value < 0)

    if not a or not b or sign(a) == sign(b):
        return sign(a) or sign(b)
    return sign(a) * sign(a * a - 5 * b * b)


def embeddings(radii):
    scales = {}
    for r in radii:
        for j, (a, b) in enumerate(CHORD2, 1):
            factor = r * r / (a * a - 5 * b * b)
            scales.setdefault((factor * a, -factor * b), []).append((r, j))
    index = {r: 5 for r in radii}
    lower = None
    for limit in sorted(scales, key=cmp_to_key(compare)):
        yield lower, limit, index.copy()
        for r, j in scales[limit]:
            assert index[r] == j
            index[r] -= 1
        lower = limit
    yield lower, None, index


def program(rings, which, profile, zlo, zhi):
    gainmax = min(fast * r / profile[which[r]][0] for r, _, fast, _, _ in rings)
    caps = [zhi, zhi, F(7, 5) * zhi, gainmax]
    rows = []

    def add(c, b):
        a = [v * w for v, w in zip(c, caps, strict=True)]
        norm = max(F(1), abs(b), *map(abs, a))
        rows.append(([z / norm for z in a], b / norm))

    add([-1, 0, 0, 0], -zlo)
    for j, high in ((1, F(1)), (2, F(7, 5))):
        a = [F(1, 10), 0, 0, 0]
        a[j] = -1
        add(a, 0)
        a = [-high, 0, 0, 0]
        a[j] = 1
        add(a, 0)
    for r, slow, fast, bar, _ in rings:
        low, high = profile[which[r]]
        add([*bar, low / r], fast)
        add([*(-x for x in bar), -high / r], -slow)
    return rows, caps


def trial(rows):
    lp = sky.solve(rows)
    assert lp.success
    if lp.x[-1] < -1e-10:
        y = [F(float(-v)).limit_denominator(10**11) for v in lp.ineqlin.marginals]
        assert min(y) >= 0
        left = [
            sum(w * a[j] for w, (a, b) in zip(y, rows, strict=True)) for j in range(4)
        ]
        gap = sum(min(v, 0) for v in left) - sum(
            w * b for w, (a, b) in zip(y, rows, strict=True)
        )
        assert gap > 0
        return None, {
            "dual": [[i, str(w)] for i, w in enumerate(y) if w],
            "margin": str(gap),
        }
    x = [F(float(v)).limit_denominator(10**11) for v in lp.x[:4]]
    assert all(0 <= v <= 1 for v in x) and all(
        sum(v * w for v, w in zip(a, x, strict=True)) <= b for a, b in rows
    )
    return x, None


def compare_galaxy(rings, bounds, leading, finite):
    lo, hi = bounds
    outer = (max(F(0), F(sky.sample.lower(lo))), F(sky.sample.upper(hi)))
    inner = (max(F(0), F(sky.sample.upper(lo))), F(sky.sample.lower(hi)))
    rejected = []
    finite_rejected = []
    finite_possible = False
    radii = sorted({r[0] for r in rings})
    for k, (low, high, index) in enumerate(embeddings(radii)):
        rows, _ = program(rings, index, leading, *outer)
        x, certificate = trial(rows)
        if x is not None:
            rows, caps = program(rings, index, leading, *inner)
            x, _ = trial(rows)
            assert x is not None
            return {
                "leading_compatible": True,
                "scale_squared_interval": [
                    list(map(str, u)) if u is not None else None for u in (low, high)
                ],
                "shell_indices": [index[r] for r in radii],
                "parameters": [str(v * w) for v, w in zip(x, caps, strict=True)],
                "finite_envelope_excludes": False,
            }
        rejected.append(certificate)
        if not finite_possible:
            rows, _ = program(rings, index, finite, *outer)
            y, c = trial(rows)
            if y is None:
                finite_rejected.append(c)
            else:
                finite_possible = True
    return {
        "leading_compatible": False,
        "leading_rejections": rejected,
        "finite_envelope_excludes": not finite_possible,
        "finite_rejections": finite_rejected if not finite_possible else [],
    }


def run():
    import networkx as nx
    from mpmath import iv
    from . import local_trapping as trap

    sky.sample.iv.dps = 40
    tr = trap.local.inputs("local_trapping")
    source = tr["surveys"][0]
    c = list(map(F, source["exact_shell_t6"]))
    leading = [(sum(c[: i + 1]) / sum(c),) * 2 for i in range(len(c))]
    board = trap.local.survey.bulk.federation(0)
    graph = nx.Graph((a // 12, b // 12) for a, b in board.external)
    hops = nx.single_source_shortest_path_length(graph, 0)
    _, coords = trap.local.survey.bulk.geometry._atlas([board])
    for v, k in hops.items():
        if k:
            a, b = CHORD2[k - 1]
            difference = coords[0][v] - coords[0][0]
            assert abs(float(difference @ difference) - float(a + b * 5**0.5)) < 1e-12
    intervals = [[] for _ in c]
    remainder = F(8) ** 4 * F(tr["duration"]) / 24
    for v, co in enumerate(map(F, source["exact_cell_t6"])):
        root = iv.sqrt(sky.sample.box(144 * co))
        lower = root - sky.sample.box(remainder)
        assert lower.a > 0
        upper = root + sky.sample.box(remainder)
        intervals[hops[v]].append(
            (F(sky.sample.lower(lower**2 / 144)), F(sky.sample.upper(upper**2 / 144)))
        )
    finite = [
        tuple(
            sum(z[i] for group in intervals[: k + 1] for z in group) / sum(c)
            for i in (0, 1)
        )
        for k in range(len(c))
    ]
    frozen = sky.HERE / "outer_requirement.json"
    table, kin = [sky.sample.measure.DATA / f for f in ("table1.dat", "table2.dat")]
    meta = {line[:11].strip(): line for line in table.read_text().splitlines()}
    raw = kin.read_text().splitlines()
    result = {}
    for name, v in json.loads(frozen.read_bytes())["galaxies"].items():
        rings = sky.rings_in(raw, v["rings"])
        row = compare_galaxy(rings, sky.sample.slope_box(meta[name]), leading, finite)
        row["lines"] = [r[-1] for r in rings]
        result[name] = row
    files = [
        HERE / "LOCAL_TRAPPING.md",
        HERE / "local_trapping.json",
        frozen,
        table,
        kin,
    ]
    return {
        "cumulative_leading": [str(x[0]) for x in leading],
        "cumulative_finite_intervals": [list(map(str, z)) for z in finite],
        "summary": {
            "leading_compatible": sum(r["leading_compatible"] for r in result.values()),
            "finite_envelope_excludes": sum(
                r["finite_envelope_excludes"] for r in result.values()
            ),
            "galaxies": len(result),
        },
        "results": result,
        "inputs": {relpath(f, HERE): sha256(f.read_bytes()).hexdigest() for f in files},
    }


if __name__ == "__main__":
    z = run()
    experiment.emit("trapped_orbits", z)
    print(z["summary"])
