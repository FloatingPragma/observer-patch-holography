"""Clock and circumference records for a multi-radius spherical release."""

import argparse
from fractions import Fraction as F
import json
from math import hypot, sqrt
from pathlib import Path

from model.candidate_laws.interval_echo_published_records.spherical_dust.body import gravitational_mass
from model.candidate_laws.interval_echo_published_records.spherical_dust.exterior_reader import (
    release_enclosure, rk4, source_counts,
)


RADII = (2, 4, 6, 8, 16)
DURATION = F(1, 8)
TICK = F(1, 2**44)
ARCHIVE = Path(__file__).with_name("radial_release.json")


def carrier(mu, radius, changed):
    """Squared drift and its radial derivative for the Hamiltonian generator."""
    correction = slope = 0
    if changed and radius > 4:
        q = 4 / radius
        correction = q * (1 - q)**3
        slope = q * (1 - q)**2 * (2 - 5*q)
    return 2*mu*(1 + correction)/radius, -2*mu*(1 + slope)/radius**2


def release(mu, start, mass, changed, panels):
    if mu == 0:
        return start
    w, _ = carrier(mu, start, changed)
    p0 = mass * sqrt(w/(1-w))

    def field(_, state):
        r, p = state
        w, dw = carrier(mu, r, changed)
        v = -sqrt(w)
        gamma = hypot(mass, p)/mass
        return v*gamma + p/mass, -gamma*p*dw/(2*v)

    return rk4(field, 0, float(DURATION), panels, (start, p0))[0]


def mass_interval(start, end, duration):
    """Necessary GM interval for inverse-square attraction from rounded radii."""
    r0, r1 = start*TICK, end*TICK
    loss = r0-r1
    lower = max(0, (loss-TICK)*2/duration**2) * (r1-TICK/2)**2
    upper = max(0, (loss+TICK)*2/duration**2) * (r0+TICK/2)**2
    return lower, upper


def endpoint_bounds(mu, r, changed):
    """Rational enclosure for the exact massive worldline, including the control."""
    r = F(r)
    floor = r - F(1, 300)
    # |1 + q(1-q)^2(2-5q)| <= 4 on 0 <= q <= 1.
    loss = 4*mu*DURATION**2/(2*floor**2)
    if loss >= r - floor:
        raise ValueError("trajectory needs a larger enclosing annulus")
    _, gradient = carrier(mu, r, changed)
    # |d[(1+2q-9q^2+12q^3-5q^4)/r^2]/dr| <= 84/r^3.
    variation = (84 if changed and r > 4 else 2) * mu * loss/floor**3
    return tuple(r + (gradient/2 + sign*variation)*DURATION**2/2 for sign in (-1, 1))


def as_numbers(interval):
    return list(map(float, interval))


def survey(mu, changed=False):
    records = []
    intervals = []
    exact_intervals = []
    for r in RADII:
        start = int(r/TICK)
        bounds = endpoint_bounds(mu, r, changed)
        squared_drift, _ = carrier(mu, F(r), changed)
        # Allow any rounded synthetic ticket compatible with the rigorous enclosure.
        small_end, large_end = (round(bound/TICK) for bound in bounds)
        exact_lower = mass_interval(start, large_end, DURATION)[0]
        exact_upper = mass_interval(start, small_end, DURATION)[1]
        exact_intervals.append((exact_lower, exact_upper))
        tickets = []
        for m in (1, 3):
            end = round(release(float(mu), r, m, changed, 16)/float(TICK))
            fine = round(release(float(mu), r, m, changed, 32)/float(TICK))
            if abs(end-fine) > 2 or not small_end <= fine <= large_end:
                raise ArithmeticError("Hamiltonian integrations disagree with the trajectory enclosure")
            if not changed or r <= 4:
                # Rescale the analytically enclosed unit-radius ratio to this release.
                scale = F(r, 2)
                low, high = release_enclosure(mu/scale, DURATION/scale)
                if not low*scale-TICK <= fine*TICK <= high*scale+TICK:
                    raise ArithmeticError("ticket disagrees with the enclosed cycloid")
            tickets.append(fine)
            intervals.append(mass_interval(start, fine, DURATION))
        records.append({"start_radius_ticks": start, "end_radius_ticks_masses_1_3": tickets,
                        "enclosed_GM_from_spherical_constraint": str(r*squared_drift/2),
                        "GM_each_from_tickets": [as_numbers(i) for i in intervals[-2:]],
                        "GM_enclosure_over_all_admissible_tickets": as_numbers((exact_lower, exact_upper))})
    measured = max(lo for lo, _ in intervals), min(hi for _, hi in intervals)
    exact = max(lo for lo, _ in exact_intervals), min(hi for _, hi in exact_intervals)
    return {"radii": records, "largest_lower_smallest_upper": as_numbers(measured),
            "largest_certified_lower_smallest_certified_upper": as_numbers(exact),
            "one_GM_admitted_by_tickets": measured[0] <= measured[1]}


def run():
    counts = source_counts()
    ordinary, collar = counts["held"]
    if counts["ordinary"] != (ordinary, 0) or counts["blank"] != (0, 0):
        raise ArithmeticError("source preparations do not match the reader")
    sources = {"blank": 0, "ordinary": ordinary, "held": ordinary+collar}
    results = {name: survey(gravitational_mass(mass)) for name, mass in sources.items()}
    mu = gravitational_mass(sources["held"])
    results["changed_outer_profile"] = survey(mu, changed=True)
    results["doubled_coupling"] = survey(2*mu)
    for name, mass in {**sources, "doubled_coupling": 2*sources["held"]}.items():
        low, high = results[name]["largest_lower_smallest_upper"]
        if not low <= float(gravitational_mass(mass)) <= high:
            raise ArithmeticError("ordinary or held strength is outside the release records")
    control = results["changed_outer_profile"]
    lo, hi = control["largest_certified_lower_smallest_certified_upper"]
    if control["one_GM_admitted_by_tickets"] or lo <= hi:
        raise ArithmeticError("the changed profile admits one inverse-square strength")
    if control["radii"][:2] != results["held"]["radii"][:2]:
        raise ArithmeticError("the changed profile altered the inner experiment")
    return {"schema": "oph.row26.spherical_radius_release.v1", "duration": str(DURATION),
            "circumference_over_2pi_tick": str(TICK), "radii": RADII,
            "source_stock": {name: str(count) for name, count in sources.items()}, "experiments": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        ARCHIVE.write_text(result)
    elif json.loads(ARCHIVE.read_text()) != json.loads(result):
        raise ArithmeticError("saved release tickets differ from the replay")
    print("Verified five radii, two probe masses, two source stocks and both changed carriers.")


if __name__ == "__main__":
    main()
