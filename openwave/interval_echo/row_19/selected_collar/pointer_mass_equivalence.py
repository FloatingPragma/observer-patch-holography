"""Test all pointer phases and timing mixtures for positive ordinary response."""

import argparse
from hashlib import sha256
from itertools import product
import json

from . import ordinary_service_distribution as distributions

HERE, mat = distributions.HERE, distributions.mat


def negate(pair):
    return [-pair[1], -pair[0]]


def with_even(readings):
    result = dict(readings)
    for meter in mat.METERS:
        a, b = (mat.iv.mpf(readings[f"{meter}/{s}"]) for s in (-1, 1))
        result["even/" + meter] = mat.lab.six.bounds((a + b) / 2)
    return result


def comparison(q, ordinary, times, groups):
    # For some base with p_base > 0, divide by k*p_base. The RHS
    # becomes O_base and the other ordinary weights p_t/p_base are
    # nonnegative. The unknown normalized quadratures may have either sign.
    result = {name: {} for name in groups}
    basis = [(t, j) for t in times for j in range(2)]
    for signs in product((-1, 1), repeat=len(basis)):
        pattern = "".join("+" if s > 0 else "-" for s in signs)
        axes = {
            f"{t}/{'RI'[j]}": (
                q[t][j] if sign > 0 else {k: negate(v) for k, v in q[t][j].items()}
            )
            for (t, j), sign in zip(basis, signs, strict=True)
        }
        for base in times:
            cols = {
                **axes,
                **{
                    t + "/ordinary": {k: negate(v) for k, v in ordinary[t].items()}
                    for t in times
                    if t != base
                },
            }
            for group, labels in groups.items():
                result[group][pattern + "/" + base] = distributions.separate(
                    cols, ordinary[base], labels
                )
    return result


def run():
    mat.iv.dps = 35
    p = HERE / "charge_selected_flight.json"
    o = HERE / "ordinary_service_distribution.json"
    frames = json.loads(p.read_bytes())["windows"]
    controls = json.loads(o.read_bytes())["windows"]
    q = {
        t: [
            with_even({k: v[j] for k, v in row["meter_quadratures"].items()})
            for j in range(2)
        ]
        for t, row in frames.items()
    }
    ordinary = {t: with_even(r["held_ordinary"]) for t, r in controls.items()}
    meters = [f"{name}/{sign}" for name in mat.METERS for sign in (-1, 1)]
    clocks = [k for k in ordinary["1"] if k.startswith("clock/")]
    names = {
        "impulses": ["shift/-1", "shift/1"],
        "material": meters,
        "charge_even": ["even/" + k for k in mat.METERS],
        "clocks": clocks,
        "all": meters + clocks,
    }
    cases = {t: comparison(q, ordinary, [t], names) for t in frames}
    cases["mixture"] = comparison(q, ordinary, list(frames), names)
    totals = {
        t: {
            group: {
                state: sum(r["status"] == state for r in tests.values())
                for state in ("excluded", "compatible_relaxation", "unresolved")
            }
            for group, tests in case.items()
        }
        for t, case in cases.items()
    }
    return {
        "quadratures": q,
        "ordinary": ordinary,
        "groups": names,
        "counts": totals,
        "certificates": cases,
        "inputs": {
            f.name: sha256(f.read_bytes()).hexdigest()
            for f in (p, o, HERE / "CHARGE_SELECTED_FLIGHT.md")
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    r = run()
    path = HERE / "pointer_mass_equivalence.json"
    if args.write:
        path.write_text(json.dumps(r, indent=1, allow_nan=False) + "\n")
    else:
        assert json.loads(path.read_bytes()) == r
    for t, z in r["counts"].items():
        print(t, z)
