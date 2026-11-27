"""Compare orbital response in bright and faint annuli of the same galaxy."""

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import sqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
KPC_M = Fraction("3.0856775814913673e19")
CUTS = {"deep": Fraction("1.2e-11"), "wider": Fraction("1.2e-10")}


@dataclass(frozen=True)
class Ring:
    line: int
    radius: Fraction
    brightness: Fraction
    gbar: Fraction
    gas_fraction: Fraction


def photometric_sources(table1, table2):
    names = {
        r[:11].strip()
        for r in table1
        if 1 <= int(r[112:115]) <= 2 and Fraction(r[30:34]) >= 30
    }
    rings = defaultdict(list)
    for line, row in enumerate(table2, 1):
        name = row[:11].strip()
        if name not in names:
            continue
        radius, gas, disk, bulge, sbd, sbb = [
            Fraction(row[i:j])
            for i, j in [(19, 25), (39, 45), (46, 52), (53, 59), (60, 67), (68, 76)]
        ]
        bright = sbd + sbb
        vsq = gas * abs(gas) + disk * abs(disk) / 2 + 7 * bulge * abs(bulge) / 10
        if radius > 0 and vsq > 0 and bright > 0:
            rings[name].append(
                Ring(
                    line,
                    radius,
                    bright,
                    vsq * 10**6 / (radius * KPC_M),
                    gas * abs(gas) / vsq,
                )
            )
    return names, rings


def select(rings, cut):
    winner = None
    best = None
    for x, y in combinations((r for r in rings if r.gbar < cut), 2):
        faint, bright = sorted((x, y), key=lambda r: r.brightness)
        gap = bright.brightness / faint.brightness
        mismatch = max(x.gbar, y.gbar) / min(x.gbar, y.gbar)
        if gap < 2 or mismatch > Fraction(11, 10):
            continue
        key = (gap, -mismatch, -min(x.line, y.line), -max(x.line, y.line))
        if best is None or key > best:
            best = key
            winner = (faint, bright)
    return winner


def speed_interval(table, row):
    return Fraction(table[row.line - 1][26:32]), Fraction(table[row.line - 1][33:38])


def ordering(lo, hi, t):
    # Both x and sqrt(t) are nonnegative. Compare their exact squares.
    if lo > t and lo * lo > t:
        return "larger_faint"
    if hi is not None and hi < t and hi * hi < t:
        return "smaller_faint"
    return "overlap"


def decode(name, pair, table):
    faint, bright = pair
    speeds = [speed_interval(table, r) for r in pair]
    (vf, ef), (vb, eb) = speeds
    assert min(vf, vb) > 0 and min(ef, eb) >= 0
    ratio = faint.gbar / bright.gbar
    radial = bright.radius / faint.radius
    obs = radial * (vf / vb) ** 2
    domains = {}
    for n in (0, 1, 2):
        lo = radial * (max(0, vf - n * ef) / (vb + n * eb)) ** 2
        hi = radial * ((vf + n * ef) / (vb - n * eb)) ** 2 if vb > n * eb else None
        domains[str(n)] = dict(
            observed_range=[
                round(float(lo), 6),
                round(float(hi), 6) if hi is not None else None,
            ],
            sign=ordering(lo, hi, ratio),
        )

    def annotation(r, s):
        return dict(
            table2_line=r.line,
            radius_kpc=float(r.radius),
            light_Lsun_pc2=float(r.brightness),
            baryon_acceleration_1e10=round(float(r.gbar * 10**10), 7),
            signed_gas_fraction=round(float(r.gas_fraction), 5),
            speed_km_s=float(s[0]),
            quoted_error_km_s=float(s[1]),
        )

    return dict(
        galaxy=name,
        distance_Mpc=float(table[faint.line - 1][12:18]),
        faint=annotation(faint, speeds[0]),
        bright=annotation(bright, speeds[1]),
        photometric_contrast=round(float(bright.brightness / faint.brightness), 4),
        baryon_acceleration_ratio=round(float(ratio), 7),
        common_amplitude_shape_range=sorted(
            [round(float(ratio), 7), round(sqrt(float(ratio)), 7)]
        ),
        observed_acceleration_ratio=round(float(obs), 6),
        speed_error_boxes=domains,
    )


def run():
    paths = [DATA / "table1.dat", DATA / "table2.dat"]
    tables = [p.read_text().splitlines() for p in paths]
    names, rings = photometric_sources(*tables)
    # Complete both enrollments before extracting the first observed speed.
    selected = {
        cut: {
            n: pair
            for n in sorted(names)
            if (pair := select(rings[n], threshold)) is not None
        }
        for cut, threshold in CUTS.items()
    }
    output = {}
    for cut, selection in selected.items():
        rows = [decode(n, p, tables[1]) for n, p in selection.items()]
        totals = {}
        for k in ("0", "1", "2"):
            count = Counter(row["speed_error_boxes"][k]["sign"] for row in rows)
            totals[k] = {
                x: count[x] for x in ("larger_faint", "overlap", "smaller_faint")
            }
        output[cut] = dict(
            selected=len(rows),
            eligible_quality_inclination_without_pair=len(names) - len(rows),
            faint_ring_outside=sum(
                r["faint"]["radius_kpc"] > r["bright"]["radius_kpc"] for r in rows
            ),
            count_by_number_of_speed_errors=totals,
            pairs=rows,
        )
    return dict(
        data_sha256={p.name: sha256(p.read_bytes()).hexdigest() for p in paths},
        catalogue_galaxies=len(tables[0]),
        quality_inclination_galaxies=len(names),
        cuts=output,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    path = HERE / "results.json"
    d = run()
    if parser.parse_args().write:
        path.write_text(json.dumps(d, indent=2) + "\n")
    else:
        assert d == json.loads(path.read_text())
        print("Replayed every matched brightness pair.")
