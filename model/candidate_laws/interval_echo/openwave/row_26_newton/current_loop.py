"""Counter-propagating light on a stationary loop around a weak spherical shell."""

import argparse
from fractions import Fraction as F
import json
from math import pi, sqrt
from pathlib import Path

import numpy as np
from mpmath.ctx_iv import MPIntervalContext
from numpy.polynomial.legendre import leggauss


MASS, SHELL_RADIUS, G, SPEED = F(1, 8), F(1), F(1, 4), F(1, 32)
TICK = F(1, 2**40)
FILE = Path(__file__).with_suffix(".json")


def solve(radius, omega, current_couples):
    # Green integral over the measured shell. u is the cosine to the detector
    # radius. Azimuthal integration leaves j along the detector's tangent.
    a, r = float(SHELL_RADIUS), float(radius)
    u, weights = leggauss(32)
    distance = np.sqrt(r*r+a*a-2*r*a*u)
    scalar = float(G*MASS)*sum(weights/distance)/2
    tangent = float(G*MASS*omega)*a*sum(weights*u/distance)/2
    gtt = -1+2*scalar
    cross = -4*tangent*r if current_couples else 0.0
    angular = (1+2*scalar)*r*r
    if gtt >= 0:
        raise ArithmeticError("loop observer cannot be held stationary")
    # Angular null roots for dt/dphi in each direction, read on the same clock.
    discr = sqrt(cross*cross-gtt*angular)
    times = [2*pi*sqrt(-gtt)*(discr+hand*cross)/(-gtt) for hand in (1, -1)]
    return [round(t/float(TICK)) for t in times]


def certify_tickets(radius, omega, coupled, tickets):
    """Use the analytic monopole/dipole integrals to check each clock rounding."""
    iv = MPIntervalContext()
    iv.dps = 48
    def val(q):
        return iv.mpf(q.numerator)/q.denominator
    u = G*MASS/radius
    j = F(2, 3)*MASS*SHELL_RADIUS**2*omega
    lapse_squared = val(1-2*u)
    cross = val(-2*G*j/radius if coupled else F(0))
    angular = val((1+2*u)*radius**2)
    root = iv.sqrt(cross**2+lapse_squared*angular)
    for hand, ticket in zip((1, -1), tickets, strict=True):
        ticks = 2*iv.pi*(root+hand*cross)/(iv.sqrt(lapse_squared)*val(TICK))
        if not (ticks.a > ticket-0.5 and ticks.b < ticket+0.5):
            raise ArithmeticError("quadrature ticket lies outside its rounding cell")


def run():
    angular_momentum = F(2, 3)*MASS*SHELL_RADIUS**2*SPEED
    rings = []
    for r in (4, 8, 16):
        clocks = {}
        for arm, coupling in (("gem", True), ("scalar_only", False)):
            clocks[arm] = {}
            for s in (-1, 0, 1):
                tickets = solve(r, s*SPEED, coupling)
                certify_tickets(r, s*SPEED, coupling, tickets)
                clocks[arm][str(s)] = tickets
        plus, reverse = clocks["gem"]["1"], clocks["gem"]["-1"]
        odd = F((plus[0]-plus[1])-(reverse[0]-reverse[1]), 2)*TICK
        expected = -8*pi*float(G*angular_momentum)/(r*sqrt(1-2*float(G*MASS)/r))
        if abs(float(odd)-expected) > float(TICK):
            raise ArithmeticError("independent signed return prediction outside its rounding allowance")
        if not (clocks["gem"]["0"][0] == clocks["gem"]["0"][1]
                and all(a == b for a, b in clocks["scalar_only"].values())
                and plus == reverse[::-1]):
            raise ArithmeticError("source reversal or scalar control failed")
        # The controlled law keeps both static clock and gravitoelectric signals fixed.
        rings.append({"reference_radius": r, "returns_in_receiver_ticks": clocks,
                      "current_odd_proper_time_difference": float(odd),
                      "exact_stationary_metric_prediction": expected,
                      "linear_GEM_prediction": -8*pi*float(G*angular_momentum)/r})
    return {"schema": "oph.row26.gem_current_loop.v1", "tick": str(TICK),
            "shell_mass": str(MASS), "shell_radius": str(SHELL_RADIUS), "G": str(G),
            "rotation_rates": [str(-SPEED), "0", str(SPEED)], "positive_J": str(angular_momentum),
            "rings": rings}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        FILE.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    elif result != json.loads(FILE.read_text()):
        raise ArithmeticError("current-loop return records differ from the archive")
    print("Verified both directions, three source currents, two couplings and three loop radii.")
