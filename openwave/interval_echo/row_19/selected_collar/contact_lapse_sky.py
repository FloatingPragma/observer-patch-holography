"""Infer allowed differences of a serial-record clock rate from outer SPARC orbits."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from os.path import relpath

from cosmology.rar_deep_regime.brightness_pairs import midpoint_halo as sky
from . import shared_contact as local

HERE = local.HERE
C = F("299792.458")


def rotation_bounds(line):
    s = sky.sample
    inc, err = F(line[30:34]), F(line[35:39])

    def sin2(v):
        return s.iv.sin(s.box(v) * s.iv.pi / 180) ** 2

    low = F(s.lower(sin2(inc) / sin2(min(F(90), inc + 2 * err))))
    high = F(s.upper(sin2(inc) / sin2(inc - 2 * err))) if inc > 2 * err else None
    d, e = F(line[15:21]), F(line[22:27])
    return low, high, max(0, (d - 2 * e) / d), (d + 2 * e) / d


def span(rings, bounds):
    vmin, vmax, dmin, dmax = bounds
    spans = {}
    for r, slow, fast, (g, d, b), _ in rings:
        ordinary = [g + d * x + b * y for x, y in sky.sample.RECT]
        low = r * max(F(0), vmin * slow - dmax * max(ordinary))
        high = r * (vmax * fast - dmin * min(ordinary)) if vmax is not None else None
        if r in spans:
            low = max(low, spans[r][0])
            if high is not None:
                high = min(high, spans[r][1])
        spans[r] = (low, high)
    radii = sorted(spans)
    if len(radii) < 2:
        return {"one_radius": True}
    lo = [max(spans[x][0] for x in radii[: i + 1]) for i in range(len(radii))]
    hi = [
        min(spans[x][1] for x in radii[i:]) if vmax is not None else None
        for i in range(len(radii))
    ]
    if vmax is not None and any(x > y for x, y in zip(lo, hi, strict=True)):
        return {"monotone_mass_bounds_conflict": True}
    weights = [
        (1 / a - 1 / b) / C**2 for a, b in zip(radii[:-1], radii[1:], strict=True)
    ]
    lower = sum(x * y for x, y in zip(weights, lo[:-1], strict=True))
    upper = (
        sum(x * y for x, y in zip(weights, hi[1:], strict=True))
        if vmax is not None
        else None
    )
    return {
        "radii_kpc": list(map(str, radii)),
        "cumulative_lower": list(map(str, lo)),
        "suffix_upper": [str(z) if z is not None else None for z in hi],
        "weak_log_lapse_span": [str(lower), str(upper) if upper is not None else None],
    }


def run():
    sky.sample.iv.dps = 40
    source = sky.HERE / "outer_requirement.json"
    names = json.loads(source.read_bytes())["galaxies"]
    files = [sky.sample.measure.DATA / f for f in ("table1.dat", "table2.dat")]
    meta = {line[:11].strip(): line for line in files[0].read_text().splitlines()}
    raw = files[1].read_text().splitlines()
    relay = local.inputs("shared_contact")["relay"]["serial"]
    fresh, shared = [F(relay[k]["probability"]) for k in ("fresh", "shared")]
    row = next(
        v
        for v in local.inputs("routed_contacts")["profiles"]
        if v["stations"] == 20 and v["stop_bits"] == 5
    )
    a, b = [
        (shared - fresh) * u / (5 - fresh * u)
        for u in map(F, row["exact_odd_on_held_edge"])
    ]
    B = sky.sample.box
    native = sky.sample.iv.log(B(1 - b) / B(1 - a))
    galaxies = {}
    for name, entry in names.items():
        d = span(sky.rings_in(raw, entry["rings"]), rotation_bounds(meta[name]))
        d["lines"] = [v["line"] for v in entry["rings"]]
        if "weak_log_lapse_span" in d:
            _, hi = d["weak_log_lapse_span"]
            if hi is not None:
                redshift = sky.sample.iv.expm1(B(F(hi)))
                cap = redshift / (B(a - b) + redshift * B(a))
                d["maximum_embedded_dose_upper"] = sky.sample.upper(cap)
        galaxies[name] = d
    used = [(name, z) for name, z in galaxies.items() if "weak_log_lapse_span" in z]
    caps = [
        (name, z["maximum_embedded_dose_upper"])
        for name, z in used
        if "maximum_embedded_dose_upper" in z
    ]
    return {
        "native_log_lapse_full_dose": [
            sky.sample.lower(native),
            sky.sample.upper(native),
        ],
        "galaxies": galaxies,
        "summary": {
            "enrolled": len(names),
            "multi_ring_spans": len(used),
            "positive_lower": sum(F(z["weak_log_lapse_span"][0]) > 0 for _, z in used),
            "finite_upper": len(caps),
            "full_dose_excluded": sum(v < 1 for _, v in caps),
            "most_restrictive_finite": list(min(caps, key=lambda v: v[1])),
            "weakest_finite": list(max(caps, key=lambda v: v[1])),
        },
        "inputs": {
            relpath(f, HERE): sha256(f.read_bytes()).hexdigest()
            for f in [
                *files,
                source,
                HERE / "CONTACT_LAPSE_SKY.md",
                HERE / "shared_contact.json",
                HERE / "routed_contacts.json",
            ]
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    z = run()
    path = HERE / "contact_lapse_sky.json"
    if parser.parse_args().write:
        path.write_text(json.dumps(z, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), z)
    print(z["native_log_lapse_full_dose"], z["summary"])
