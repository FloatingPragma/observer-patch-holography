"""Prepare and read two distinct loads on a connected Interval Echo federation."""
import argparse
from fractions import Fraction
import json
from math import pi

import numpy as np
from scipy.linalg import expm

from . import open_pair
from .corridor_release import monitored_source


ROUNDS = 3
ATTEMPTS = 8
EXPOSURES = (1, 7, 37)


def source_bound(rounds=ROUNDS, attempts=ATTEMPTS):
    """Unconditional dark weight for a finite number of pulses and sharp retries."""
    if any(type(n) is not int or n < 1 for n in (rounds, attempts)):
        raise ValueError("Expected positive integer preparation budgets")
    bright = Fraction(1, 9)**rounds
    failed = (1-bright)/2
    return (1-failed**attempts)/(1+bright), failed**attempts


def source():
    installation = open_pair.seal_from_sharp()
    exchange, contrast = open_pair.source_parts()
    dark = np.array([1,0,-1])
    if (not np.array_equal(2*exchange @ dark, dark)
            or not np.array_equal(2*exchange, [[1,-1,0],[-1,2,-1],[0,-1,1]])
            or not np.array_equal(contrast, 2*np.eye(3))):
        raise ArithmeticError("The native source does not have the certified bright and dark blocks")
    # Every addressed pulse on this three-state source has the same scalar
    # contrast 2*kappa. Dropping it keeps every branch probability and phase difference.
    pulse = expm((-2j*pi/3)*exchange)
    rho, exhausted, reach, first = monitored_source(pulse, ROUNDS, ATTEMPTS)
    floor, predicted_exhaustion = source_bound()
    if (abs(dark @ rho @ dark/2 - float(floor)) > 5e-14
            or abs(exhausted-float(predicted_exhaustion)) > 5e-14):
        raise ArithmeticError("Acquired source differs from the exact finite preparation")
    open_pair.dark_certificate()
    return rho, {"setup": installation, "first_vacancy_in_one_attempt": first,
                 "attempt_reach": reach, "exhaustion": exhausted, "guaranteed_joint_home": str(floor)}


def terminal(source_density, exhaustion, kappa, t):
    if not 0 <= exhaustion <= 1 or abs(source_density.trace().real+exhaustion-1) > 5e-13:
        raise ValueError("Inconsistent source acceptance and exhaustion")
    arms = {}
    for arm in ("coherent", "dephased", "stochastic_swap", "unfiltered"):
        raw = arm == "unfiltered"
        rho = np.diag([1.,0.,0.]) if raw else source_density
        bins = open_pair.read_release(kappa,t,rho,flight="coherent" if raw else arm)
        bins[-1] += 0 if raw else exhaustion
        if not np.isfinite(bins).all() or min(bins) < 0 or abs(bins.sum()-1) > 2e-10:
            raise ArithmeticError("Incomplete joint pair receiver")
        arms[arm] = bins.tolist()
    dark = np.array([1, 0, -1])
    protected = (dark @ source_density @ dark).real/2
    if sum(arms["coherent"][:-1])+2e-10 < protected:
        raise ArithmeticError("Native release lost protected weight")
    return arms


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--records", type=float, nargs=2, metavar=("KAPPA","EXPOSURE"))
    group.add_argument("--full", action="store_true", help="emit all enrolled appointments and arms")
    args = parser.parse_args()
    rho, report = source()
    if args.records or args.full:
        runs = [args.records] if args.records else [(k,t) for k in (0,.5,8) for t in EXPOSURES]
        report["home_seams_in_order"] = open_pair.geometry()[1]
        report["appointments"] = [dict(kappa=k, exposure=t, bins=terminal(rho,report["exhaustion"],k,t))
                                  for k,t in runs]
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
