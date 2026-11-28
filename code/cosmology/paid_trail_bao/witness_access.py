"""Invert orbital bands for constant and spatially variable access to an extra record."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json

from . import record_profiles as rp

HERE = rp.HERE


def rational(x):
    return Q(str(x))


class Bounds:
    """Exact nonnegative gain interval from inequalities coefficient * gain <= rhs."""

    def __init__(self):
        self.lower = Q(0)
        self.upper = None
        self.low_pair = "nonnegative gain"
        self.high_pair = None
        self.impossible = None

    def add(self, coefficient, rhs, label):
        if coefficient == 0:
            if rhs < 0:
                self.impossible = label
        elif coefficient < 0:
            value = rhs / coefficient
            if value > self.lower:
                self.lower, self.low_pair = value, label
        else:
            value = rhs / coefficient
            if self.upper is None or value < self.upper:
                self.upper, self.high_pair = value, label

    def archive(self):
        return dict(
            possible=self.impossible is None
            and (self.upper is None or self.lower <= self.upper),
            lower=str(self.lower),
            upper=None if self.upper is None else str(self.upper),
            lower_constraint=self.low_pair,
            upper_constraint=self.high_pair,
            contradiction=self.impossible,
        )


def constant_access(radius, low, high, x, y):
    # B + X*x_i + Y*y_i must fall in every band with B, X and Y nonnegative.
    # Eliminate B using all ordered lower/upper pairs (and 0 <= B).
    y_lower = [(Q(0), Q(0), "nonnegative third-head gain")]
    y_upper = []
    gain = Bounds()

    def halfspace(a, b, c, label):
        if b == 0:
            gain.add(a, c, label)
        else:
            (y_upper if b > 0 else y_lower).append((-a / b, c / b, label))

    for j, outer in enumerate(radius):
        halfspace(x[j], y[j], high[j], [0, outer])
        for i, inner in enumerate(radius):
            halfspace(x[j] - x[i], y[j] - y[i], high[j] - low[i], [inner, outer])
    for (a, b, i), (c, d, j) in product(y_lower, y_upper):
        gain.add(a - c, d - b, [i, j])
    out = gain.archive()
    if out["possible"]:
        g0 = gain.lower
        g1 = max(a * g0 + b for a, b, _ in y_lower)
        B = max(
            [
                Q(0),
                *(v - g0 * xx - g1 * yy for v, xx, yy in zip(low, x, y, strict=True)),
            ]
        )
        mass = [B + g0 * xx + g1 * yy for xx, yy in zip(x, y, strict=True)]
        assert all(
            lower <= v <= upper for lower, v, upper in zip(low, mass, high, strict=True)
        )
        out.update(
            native_gain=str(g0),
            third_head_gain=str(g1),
            first_mass_msun=str(B),
            third_head_access=None if g0 + g1 == 0 else float(g1 / (g0 + g1)),
        )
    return out


def variable_access(radius, low, high, shell_low, shell_high):
    gain = Bounds()
    for i, r in enumerate(radius):
        lo, hi = Q(0), Q(0)
        for j in range(i + 1, len(radius)):
            lo += shell_low[j]
            hi += shell_high[j]
            gain.add(lo, high[j] - low[i], [r, radius[j], "minimum shell"])
            gain.add(-hi, high[i] - low[j], [r, radius[j], "maximum shell"])
    out = gain.archive()
    if out["possible"]:
        g = gain.lower
        feasible = [(low[0], high[0])]
        for lower, upper, a, b in zip(
            low[1:], high[1:], shell_low[1:], shell_high[1:], strict=True
        ):
            old_l, old_u = feasible[-1]
            feasible.append((max(lower, old_l + g * a), min(upper, old_u + g * b)))
            assert feasible[-1][0] <= feasible[-1][1]
        mass = [feasible[-1][0]]
        for i in reversed(range(1, len(radius))):
            lower, upper = feasible[i - 1]
            v = max(lower, mass[-1] - g * shell_high[i])
            assert v <= min(upper, mass[-1] - g * shell_low[i])
            mass.append(v)
        out["feasible_msun"] = list(map(str, reversed(mass)))
    return out


def run():
    inputs, *_ = rp.galaxy.transfer.budget.load_inputs()
    name = "NGC2841"
    row = inputs[name]
    indices = [i for i, r in enumerate(row["radius_kpc"]) if 20 <= r <= 50]
    assert indices == list(range(indices[0], indices[-1] + 1))
    radius = [row["radius_kpc"][i] for i in indices]
    G = rational(rp.galaxy.transfer.budget.local.profiles.galaxies.G)
    low, high = [], []
    for i in indices:
        r, v, e, b = (
            rational(row[k][i])
            for k in ("radius_kpc", "observed_kms", "error_kms", "baryons_kms")
        )
        low.append(r * (max(0, v - 2 * e) ** 2 - b * b) / G)
        high.append(r * ((v + 2 * e) ** 2 - b * b) / G)
    assert min(low) > 0
    src = json.loads((HERE / "record_profile_maps.json").read_bytes())
    tests = {}
    for tag, scales in src["ngc2841_shell_access"].items():
        tests[tag] = {}
        for cap, shell in scales.items():
            profiles = src["enclosed"][tag][cap]

            def select(seq):
                return [rational(seq[i]) for i in indices]

            def relative(seq):
                v = select(seq)
                return [a - v[0] for a in v]

            a = constant_access(
                radius,
                low,
                high,
                relative(profiles["none"][name]),
                relative(profiles["third_head"][name]),
            )
            b = variable_access(
                radius, low, high, select(shell["lower"]), select(shell["upper"])
            )
            assert not a["possible"] or b["possible"]
            tests[tag][cap] = dict(constant=a, cellwise=b)
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "WITNESS_ACCESS.md",
                HERE / "record_profile_maps.json",
                HERE / "source_profile.json",
            )
        },
        radius_kpc=radius,
        mass_lower_msun=list(map(str, low)),
        mass_upper_msun=list(map(str, high)),
        tests=tests,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    path = HERE / "witness_access.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        rp.galaxy.check(json.loads(path.read_bytes()), out)
    print(
        {
            t: {
                c: (v["constant"]["possible"], v["cellwise"]["possible"])
                for c, v in caps.items()
            }
            for t, caps in out["tests"].items()
        }
    )
