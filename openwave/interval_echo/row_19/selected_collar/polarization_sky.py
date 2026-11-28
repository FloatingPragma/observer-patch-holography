"""Test a decreasing excess enclosed mass and the composition dependence of the passive source."""

from fractions import Fraction as F
from hashlib import sha256
import json
from os.path import relpath

from . import contact_census as experiment
from .contact_lapse_sky import sky

HERE = experiment.HERE


def constraints(rings, zlo, zhi):
    caps = [zhi, zhi, F(7, 5) * zhi]
    rows = []

    def add(vector, b):
        v = [x * c for x, c in zip(vector, caps, strict=True)]
        size = max(F(1), abs(b), *map(abs, v))
        rows.append(([x / size for x in v], b / size))

    add([-1, 0, 0], -zlo)
    for i, high in ((1, F(1)), (2, F(7, 5))):
        low = [F(1, 10), 0, 0]
        low[i] = -1
        upper = [-high, 0, 0]
        upper[i] = 1
        add(low, 0)
        add(upper, 0)
    for r, _, upper, bar, _ in rings:
        add(bar, upper)
        for s, lower, _, other, _ in rings:
            if s >= r:
                add(
                    [r * x - s * y for x, y in zip(bar, other, strict=True)],
                    r * upper - s * lower,
                )
    return rows, caps


def certify(rings, box):
    lo, hi = box
    outer, _ = constraints(
        rings, max(F(0), F(sky.sample.lower(lo))), F(sky.sample.upper(hi))
    )
    res = sky.solve(outer)
    assert res.success
    if res.x[-1] < -1e-10:
        dual = [F(float(-u)).limit_denominator(10**11) for u in res.ineqlin.marginals]
        weights = [
            sum(y * a[j] for y, (a, b) in zip(dual, outer, strict=True))
            for j in range(3)
        ]
        gain = sum(min(x, 0) for x in weights) - sum(
            y * b for y, (a, b) in zip(dual, outer, strict=True)
        )
        assert min(dual) >= 0 and gain > 0
        return {
            "possible": False,
            "dual": [[k, str(v)] for k, v in enumerate(dual) if v],
            "rejection_margin": str(gain),
        }
    inner, caps = constraints(
        rings, max(F(0), F(sky.sample.upper(lo))), F(sky.sample.lower(hi))
    )
    trial = sky.solve(inner)
    assert trial.success
    x = [F(float(z)).limit_denominator(10**11) for z in trial.x[:3]]
    assert min(x) >= 0 and max(x) <= 1
    assert all(sum(v * w for v, w in zip(row, x, strict=True)) <= b for row, b in inner)
    parameters = [v * c for v, c in zip(x, caps, strict=True)]
    lows = {}
    for r, slow, _, bar, _ in rings:
        val = r * (slow - sum(a * b for a, b in zip(bar, parameters, strict=True)))
        lows[r] = max(F(0), val, lows.get(r, 0))
    enclosed = {r: max(z for s, z in lows.items() if s >= r) for r in lows}
    for r, slow, fast, bar, _ in rings:
        model = enclosed[r] / r + sum(
            a * b for a, b in zip(bar, parameters, strict=True)
        )
        assert slow <= model <= fast
    return {
        "possible": True,
        "parameters": list(map(str, parameters)),
        "mass_witness": [[str(r), str(enclosed[r])] for r in sorted(enclosed)],
    }


def run():
    sky.sample.iv.dps = 40
    source = sky.HERE / "outer_requirement.json"
    table, kin = [sky.sample.measure.DATA / f for f in ("table1.dat", "table2.dat")]
    meta = {line[:11].strip(): line for line in table.read_text().splitlines()}
    raw = kin.read_text().splitlines()
    rows = {}
    for name, entry in json.loads(source.read_bytes())["galaxies"].items():
        rings = sky.rings_in(raw, entry["rings"])
        row = certify(rings, sky.sample.slope_box(meta[name]))
        row["catalogue_lines"] = [r[-1] for r in rings]
        row["distinct_radii"] = len({r[0] for r in rings})
        rows[name] = row
    # Eta is normalized by the inferred Galactic dark acceleration in both
    # publications. The native composition of that reference is unknown.
    external = {
        "Be_Ti": {
            "source": "https://doi.org/10.1103/PhysRevLett.100.041101",
            "observable": "differential acceleration over Galactic dark acceleration",
            "eta": "-0.00004",
            "one_sigma": "0.00007",
            "absolute_envelope": "0.00018",
            "allowance": "absolute central value plus twice the quoted uncertainty",
        },
        "PSR_J1713+0747": {
            "source": "https://doi.org/10.1103/PhysRevLett.120.241104",
            "observable": "neutron-star and white-dwarf differential fall toward Galactic dark acceleration",
            "absolute_envelope": "0.004",
            "allowance": "rounded published pulsar limit",
        },
    }
    for r in external.values():
        cap = F(r["absolute_envelope"])
        r["max_source_fraction_times_composition_difference"] = str(cap * 65 / 48)
    files = [
        source,
        table,
        kin,
        HERE / "CONTACT_CENSUS.md",
        HERE / "MATERIAL_POLARIZATION.md",
        HERE / "material_polarization.json",
    ]
    return {
        "galaxies": rows,
        "summary": {
            "ring_count": sum(len(r["catalogue_lines"]) for r in rows.values()),
            "single_radius": sum(r["distinct_radii"] == 1 for r in rows.values()),
            "compatible": sum(r["possible"] for r in rows.values()),
            "excluded": [name for name, row in rows.items() if not row["possible"]],
        },
        "galactic_free_fall": external,
        "pure_response_min_max": ["48/65", "48/17"],
        "inputs": {relpath(f, HERE): sha256(f.read_bytes()).hexdigest() for f in files},
    }


if __name__ == "__main__":
    result = run()
    experiment.emit("polarization_sky", result)
    print(result["summary"])
