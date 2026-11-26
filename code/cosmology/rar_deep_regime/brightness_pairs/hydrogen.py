"""Transfer independent total H I masses into the two-annulus orbital test."""

import argparse
from collections import Counter
import csv
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json

if __package__:
    from . import sensitivity as s
else:
    import sensitivity as s


@dataclass(frozen=True)
class Interval:
    lo: F = F(0)
    hi: F | None = None
    open_lo: bool = False
    open_hi: bool = False

    def cut(self, a, c, strict=False):
        """Intersect with a*x <= c (or < c)."""
        if not a:
            return self if c > 0 or (c == 0 and not strict) else None
        lo, hi, ol, oh = self.lo, self.hi, self.open_lo, self.open_hi
        b = c / a
        if a < 0 and b >= lo:
            ol = strict or (b == lo and ol)
            lo = b
        if a > 0 and (hi is None or b <= hi):
            oh = strict or (b == hi and oh)
            hi = b
        if hi is not None and (lo > hi or (lo == hi and (ol or oh))):
            return None
        return Interval(lo, hi, ol, oh)

    def report(self):
        return {
            "lower": float(self.lo),
            "upper": float(self.hi) if self.hi is not None else None,
            "endpoints_excluded": [self.open_lo, self.open_hi],
        }


def project(rows, disk, gas):
    """Eliminate disk M/L exactly from linear inequalities a*g + b*d <= c."""
    lower = [(0, disk[0], False)]
    upper = [(0, disk[1], False)]
    limits = []
    for a, b, c, strict in rows:
        if b:
            (lower if b < 0 else upper).append((-a / b, c / b, strict))
        else:
            limits.append((a, c, strict))
    for m, h, op in lower:
        for n, j, oq in upper:
            limits.append((m - n, j - h, op or oq))
    for a, c, strict in limits:
        gas = gas.cut(a, c, strict)
        if gas is None:
            return None
    return gas


def records(pair, table):
    f, h = [s.measured(table[pair[k]["table2_line"] - 1]) for k in ("faint", "bright")]
    assert not f[3][2] and not h[3][2]  # Both measured bulge velocities are zero.
    lo = h[0] / f[0] * (max(0, f[1] - 2 * f[2]) / (h[1] + 2 * h[2])) ** 2
    hi = (
        h[0] / f[0] * ((f[1] + 2 * f[2]) / (h[1] - 2 * h[2])) ** 2
        if h[1] > 2 * h[2]
        else None
    )
    return (
        f[3][:2],
        h[3][:2],
        min(lo, lo**2),
        max(hi, hi**2) if hi is not None else None,
    )


def constraints(data, matched=False, lower_shape=False, upper_shape=False):
    f, h, left, right = data

    def compare(x, t, y):
        return (x[0] - t * y[0], x[1] - t * y[1], 0, False)

    rows = [(-g, -d, 0, True) for g, d in (f, h)]
    if matched:
        rows.extend((g, d, s.DEEP, True) for g, d in (f, h))
        rows += [compare(f, s.MATCH, h), compare(h, s.MATCH, f)]
    if lower_shape:
        rows.append((left * h[0] - f[0], left * h[1] - f[1], 0, False))
    if upper_shape and right is not None:
        rows.append(compare(f, right, h))
    return rows


def shape(data, disk, gas, matched=False):
    def possible(low, high):
        return project(constraints(data, matched, low, high), disk, gas) is not None

    if not possible(False, False):
        return "no_eligible_force"
    if not possible(True, False):
        return "larger_faint"
    if not possible(False, True):
        return "smaller_faint"
    assert possible(True, True)
    return "overlap"


def power10_bounds(exponent):
    """Enclose 10^exponent on a 1e-10 rational grid; certify by integer powers."""
    q = F(exponent)
    with localcontext() as ctx:
        ctx.prec = 32
        scaled = 10**10 * (
            Decimal(10) ** (Decimal(q.numerator) / Decimal(q.denominator))
        )
        guess = int(scaled)
    lo, hi = F(guess - 1, 10**10), F(guess + 2, 10**10)
    target = F(10) ** q.numerator
    assert lo > 0 and lo**q.denominator <= target <= hi**q.denominator
    return lo, hi


def calibrations(row):
    x, y = (Decimal(row[k]) for k in ("log_alfalfa", "log_sparc"))
    width = (
        Decimal(1).scaleb(x.as_tuple().exponent)
        + Decimal(1).scaleb(y.as_tuple().exponent)
    ) / 2
    central = Interval(*power10_bounds(x - y))
    rec = Interval(power10_bounds(x - y - width)[0], power10_bounds(x - y + width)[1])
    return {
        "printed_centres": central,
        "print_rounding": rec,
        "rounding_and_mass_errors": Interval(rec.lo * F(9, 11), rec.hi * F(11, 9)),
    }


def single_dish(target):
    with (s.HERE / "ugc06667_flux.tsv").open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(rows) == 1 and rows[0]["galaxy"] == "U6667"
    r = rows[0]
    assert not r["flux_flag"] and r["reference_code"] == "1"
    wsrt, ew, dish, ed = [
        F(r[k]) for k in ("wsrt", "wsrt_error", "single_dish", "single_dish_error")
    ]
    t1 = (s.HERE.parent / "data/table1.dat").read_text().splitlines()
    row = next(row for row in t1 if row[:11].strip() == "UGC06667")
    mass, distance = F(row[86:93]) * 10**9, F(row[15:21])
    catalogue_flux = mass / (235600 * distance**2)
    lower_flux = dish - 2 * ed
    # Read the exact rational requirements directly from the same two annuli.
    domain = project(
        constraints(target, matched=True, upper_shape=True), s.DISK, Interval()
    )
    return {
        "galaxy": "UGC06667",
        "inclination_deg": float(F(row[30:34])),
        "dish_to_wsrt_central": float(dish / wsrt),
        "dish_to_wsrt_two_error_box": [
            float(lower_flux / (wsrt + 2 * ew)),
            float((dish + 2 * ed) / (wsrt - 2 * ew)),
        ],
        "catalogue_distance_Mpc": float(distance),
        "catalogue_HI_mass_Msun": float(mass),
        "catalogue_equivalent_flux_Jy_km_s": float(catalogue_flux),
        "smallest_dish_to_catalogue_two_error_ratio": float(
            lower_flux / catalogue_flux
        ),
        "largest_shape_allowed_flux_Jy_km_s": float(catalogue_flux * domain.hi),
        "shared_gas_mass_suffices_at_two_errors": lower_flux
        <= catalogue_flux * domain.hi,
    }


def run():
    pairs = json.loads((s.HERE / "results.json").read_text())["cuts"]["deep"]["pairs"]
    table = (s.HERE.parent / "data/table2.dat").read_text().splitlines()
    with (s.HERE / "hydrogen_masses.tsv").open() as f:
        masses = {r["galaxy"]: r for r in csv.DictReader(f, delimiter="\t")}
    matches, missing, thresholds = {}, [], {}
    for pair in pairs:
        name = pair["galaxy"]
        if name in masses:
            data = records(pair, table)
            matches[name] = {}
            for label, gas in calibrations(masses[name]).items():
                matches[name][label] = {
                    "gas": gas.report(),
                    "nominal": shape(data, (F(".5"), F(".5")), gas),
                    "varied": shape(data, s.DISK, gas),
                    "varied_requiring_original_cuts": shape(data, s.DISK, gas, True),
                }
        else:
            missing.append(name)
        if pair["speed_error_boxes"]["2"]["sign"] == "smaller_faint":
            data = records(pair, table)
            thresholds[name] = {}
            for label, matched in (("fixed_annuli", False), ("retaining_cuts", True)):
                limit = project(
                    constraints(data, matched, upper_shape=True), s.DISK, Interval()
                )
                thresholds[name][label] = limit.report() if limit else None
    return {
        "single_dish_target": single_dish(
            records(next(p for p in pairs if p["galaxy"] == "UGC06667"), table)
        ),
        "total_pairs": len(pairs),
        "published_matches": matches,
        "pairs_without_published_match": missing,
        "lower_shape_gas_domains": thresholds,
        "counts": {
            phase: {
                test: dict(Counter(v[phase][test] for v in matches.values()))
                for test in ("nominal", "varied", "varied_requiring_original_cuts")
            }
            for phase in (
                "printed_centres",
                "print_rounding",
                "rounding_and_mass_errors",
            )
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    receipt = s.HERE / "hydrogen.json"
    result = run()
    if parser.parse_args().write:
        receipt.write_text(json.dumps(result, indent=2) + "\n")
    else:
        assert json.loads(receipt.read_text()) == result
        print(
            "Replayed all published hydrogen matches and shared-normalization bounds."
        )
