"""Reduce the published outer H I orbit, component profiles and gas masses."""

import argparse
import hashlib
import io
import json
import math
import re
import subprocess
from pathlib import Path
from zipfile import ZipFile

import numpy as np

from ..udg_transfer.measure import GM_SUN, KPC

SOURCES = {
    "catalog": (
        "https://astroweb.case.edu/SPARC/ETG_Lelli2017.mrt",
        "bf48cf65ffd0b2e9cde960d915798da3cb795527cd10f4e0de71b4ef152c4376",
    ),
    "profiles": (
        "https://astroweb.case.edu/SPARC/Rotmod_ETG.zip",
        "e640cd0f387588b83d8755321a9f385bb674369b6c42f05e89cfafdccad7af34",
    ),
    "hi": (
        "https://arxiv.org/pdf/1509.05236",
        "9b0d3f541a41d4daa3e7ca624749fb5c69cd0409880f749e300a44f3cf85d620",
    ),
}


def header_value(raw, label):
    return float(re.search(r"^!" + label + r"\s*:\s*([\d.]+)", raw, re.M)[1])


def enclosed_disk(raw, radius):
    """Integrate tabulated piecewise-linear surface density over an inscribed cylinder."""
    if "!Z-Density law   : EXPONENTIAL" not in raw:
        raise ValueError("Unexpected vertical distribution")
    model_r, density, _ = np.loadtxt(io.StringIO(raw), comments="!", unpack=True)
    if model_r[0] != 0 or np.any(np.diff(model_r) <= 0) or np.any(density < 0):
        raise ValueError("Invalid disk table")
    limit = min(0.9 * radius, model_r[-1])
    r = np.r_[model_r[model_r < limit], limit]
    a, b = r[:-1], r[1:]
    middle = (a + b) / 2

    def integrand(x):
        return x * np.interp(x, model_r, density)

    integral = (
        (b - a) * (integrand(a) + 4 * integrand(middle) + integrand(b))
    ).sum() / 6
    height = header_value(raw, "Z Scaleheight")
    return (
        2
        * math.pi
        * 1e6
        * integral
        * (-math.expm1(-math.sqrt(1 - 0.9**2) * radius / height))
    )


def reduce(catalog, profiles, hi):
    for key, path in zip(SOURCES, (catalog, profiles, hi)):
        if hashlib.sha256(path.read_bytes()).hexdigest() != SOURCES[key][1]:
            raise ValueError(f"Unexpected primary source digest for {key}")
    text = subprocess.check_output(["pdftotext", "-layout", str(hi), "-"], text=True)
    hi_rows = {}
    for line in text.splitlines():
        t = line.split()
        if (
            len(t) == 12
            and t[0] in {"NGC", "UGC"}
            and t[1].isdigit()
            and t[3] in {"c", "v"}
        ):
            hi_rows[t[0] + t[1]] = (float(t[2]), float(t[6]), float(t[8]))
    rows = []
    with ZipFile(profiles) as archive:
        for line in catalog.read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            t = line.split()
            name = t[0]
            d, ed, inc, ei, lum, reff = (float(t[i]) for i in (1, 2, 4, 5, 6, 8))
            hi_dist, hi_speed, loghi = hi_rows[name]

            def read(part):
                return archive.read(f"{name}_{part}.dat").decode()

            rad, v, ev, vg, vd, vb = np.loadtxt(io.StringIO(read("rotmod")))[-1]
            disk = read("disk")
            bulge = read("bulge")
            md = 1e9 * header_value(disk, "Mass of disk")
            mb = 1e9 * header_value(bulge, "Mass of bulge")
            enclosed = min(md, enclosed_disk(disk, rad)) + min(
                mb, rad * KPC * vb**2 * 1e6 / GM_SUN
            )
            gas = 1.4 * 10**loghi * (d / hi_dist) ** 2
            logobs = math.log10(v**2 * 1e6 / (rad * KPC))
            logstar = math.log10(0.8 * (vd**2 + vb**2) * 1e6 / (rad * KPC))
            if (
                abs(v - hi_speed) > ev
                or vg != 0
                or abs((md + mb) / (lum * 1e9) - 1) > 0.02
            ):
                raise ValueError(
                    f"Component census or H I identification failed for {name}"
                )
            rows.append(
                dict(
                    name=name,
                    distance_mpc=d,
                    distance_error_mpc=ed,
                    inc_deg=inc,
                    inc_error_deg=ei,
                    radius_kpc=rad,
                    re_kpc=reff,
                    speed_kms=v,
                    speed_error_kms=ev,
                    vdisk_unit_kms=vd,
                    vbul_unit_kms=vb,
                    stellar_unit_msun=md + mb,
                    stellar_inside_unit_msun=enclosed,
                    helium_hi_msun=gas,
                    source_speed_hi_kms=hi_speed,
                    source_loghi=loghi,
                    source_distance_hi_mpc=hi_dist,
                    loggobs_check_dex=logobs - float(t[14]),
                    loggstar_check_dex=logstar - float(t[18]),
                )
            )
    if len(rows) != 16 or set(hi_rows) != {r["name"] for r in rows}:
        raise ValueError("Sample mismatch")
    return dict(
        sources={k: dict(url=v[0], sha256=v[1]) for k, v in SOURCES.items()},
        galaxies=rows,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in SOURCES:
        p.add_argument("--" + key, type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    result = reduce(args.catalog, args.profiles, args.hi)
    with args.output.open("x") as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write("\n")
