"""Acquire the M10.3 seam records after bounded neutral-collar preparation."""
import argparse
from fractions import Fraction
import json
from math import exp, hypot, isfinite, pi, sin, sqrt

import numpy as np
from scipy.linalg import expm

from .corridor_pair import BASIS, LAW, LIVE, corridor_generator, exact_flight, receiver_bins


ROUNDS = 128
ATTEMPTS = 8
PULSE = Fraction(1, 2)


def uniform_floor(pulse=PULSE, rounds=ROUNDS, attempts=ATTEMPTS):
    """Joint lower bound at every real coupling and every terminal exposure."""
    if (not isinstance(pulse, Fraction) or not 0 < pulse <= 1 or type(rounds) is not int
            or rounds < 1 or type(attempts) is not int or attempts < 1):
        raise ValueError("Expected a positive rational pulse <= 1 and positive finite source budgets")
    # In the first half lobe sin(x)/x >= 2/pi and pi² < 10.
    # Beyond it the bright-to-quiet amplitude squared is <= 2*pulse²/pi².
    # The bright fraction of the accepted release is always <= 1/2.
    bright = (1-pulse**2/5)**rounds
    loss = max(bright/(1+bright), pulse**2/9)
    return (1-Fraction(1, 2)**attempts)*(1-loss)


def monitored_source(step, rounds, attempts):
    """First-vacancy distribution on one sharp attempt and the accepted ensemble."""
    state = np.array([1, 0, 0], dtype=complex)
    first = []
    for _ in range(rounds):
        state = step @ state
        first.append(float(abs(state[1])**2))
        state[1] = 0
    loss = sum(first)
    if abs(loss + np.vdot(state,state).real - 1) > 5e-12:
        raise ArithmeticError("Incomplete first-vacancy instrument")
    reaches = [loss**i for i in range(attempts)]
    rho = sum(reaches)*np.outer(state, state.conj())
    exhausted = loss**attempts
    if abs(rho.trace().real+exhausted-1) > 5e-12:
        raise ArithmeticError("Lost a source branch")
    return rho, exhausted, reaches, first


def verify_retry(kappa):
    """A first vacancy puts the same two loads at ac; exchange on bc restores ab."""
    ids = [LAW.SRC.index(word) for word in ((0,1),(1,0))]
    swap = expm(-1j*pi*LAW.seam_term(kappa)[np.ix_(ids, ids)])
    if not np.allclose(abs(swap[:,0])**2, [0,1], rtol=0, atol=1e-13):
        raise ArithmeticError("The retry does not restore the original sharp pair")


def flight_unitary(kappa, t):
    """Exactly unitary flight from the spectrum of the real symmetric three-word generator."""
    h = corridor_generator(kappa)
    if not np.array_equal(h, h.T):
        raise ArithmeticError("Expected a real symmetric corridor generator")
    energies, modes = np.linalg.eigh(h)
    return (modes*np.exp(-1j*t*energies)) @ modes.T


def acquire(kappa):
    verify_retry(kappa)
    return monitored_source(flight_unitary(kappa, float(PULSE)), ROUNDS, ATTEMPTS)


def analytic_home(kappa, t):
    """Joint home probability from the bright and dark solutions, with all retries."""
    omega = hypot(float(kappa)+.25, sqrt(.5))
    q = .5/omega**2
    r = (1-q*sin(float(PULSE)*omega)**2)**ROUNDS
    exhaustion = ((1-r)/2)**ATTEMPTS
    return (1-exhaustion)*(1-q*r/(1+r)*sin(t*omega)**2)


def terminal(kappa, t, source):
    if not isfinite(t) or t < 0:
        raise ValueError("Expected a finite nonnegative exposure")
    rho, exhaustion, _, _ = source
    flight = flight_unitary(kappa, t)
    ports = {}
    for arm, density, failure in (("coherent", rho, exhaustion),
                                  ("dephased", np.diag(np.diag(rho)), exhaustion),
                                  ("unfiltered", np.diag([1., 0., 0.]), 0.)):
        population = np.diag(flight @ density @ flight.conj().T).real
        if population.min() < -2e-13:
            raise ArithmeticError("Negative terminal population")
        bins = receiver_bins(np.maximum(population,0))
        bins[-1] += failure
        ports[arm] = bins
    # A rate-1/2 full swap on either live seam moves the vacancy on a three-site path.
    # Evolve its entire spatial distribution including unequal endpoint populations.
    x, _, z = np.diag(rho).real
    mass = x+z
    mid = mass*(-np.expm1(-1.5*t))/3
    difference = (x-z)*exp(-.5*t)
    swap = receiver_bins([(mass-mid+difference)/2, mid, (mass-mid-difference)/2])
    swap[-1] += exhaustion
    ports["stochastic_swap"] = swap
    if (any(min(bins) < 0 or abs(sum(bins)-1) > 5e-12 for bins in ports.values())
            or abs(sum(ports["coherent"][:-1])-analytic_home(kappa,t)) > 5e-12):
        raise ArithmeticError("Incomplete seam law or disagreement with the solvable flight")
    return ports


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--records", type=float, nargs=2, metavar=("KAPPA","EXPOSURE"),
                     help="emit source branches, exhaustion and all four complete M10.3 detector laws")
    args = cli.parse_args()
    installation = exact_flight()
    floor = uniform_floor()
    if floor <= Fraction(94,100):
        raise ArithmeticError("Finite source misses the continuous home gate")
    if args.records is not None:
        kappa, t = args.records
        source = acquire(kappa)
        bins = terminal(kappa,t,source)
        print(json.dumps(dict(seals=installation["neutral_edges"], pair_basis=BASIS,
                              reset_contact=LIVE[1], first_vacancy=source[3], attempt_reach=source[2],
                              source_exhaustion=source[1], bins={k:v.tolist() for k,v in bins.items()}),
                         allow_nan=False))
        return
    print(f"All real strengths, all free exposures: joint home probability >= {floor}")
    for kappa in (-8, -.25, 0, .5, 8, 128, sqrt((2*pi)**2-.5)-.25):
        source = acquire(kappa)
        bins = terminal(kappa,37.,source)
        print(f"kappa={kappa:.8g} at exposure 37: coherent home {sum(bins['coherent'][:-1]):.9f}, "
              f"swap home {sum(bins['stochastic_swap'][:-1]):.9f}, exhaustion {source[1]:.9g}")


if __name__ == "__main__":
    main()
