"""Certified local phase service on the vacuum and one-load sectors of three ports."""

from functools import lru_cache

from mpmath import iv

from . import circuit

# Each decimal is the specified exposure, in native exchange times. Obtained by
# fitting only a one-load input/output gate, without network or Bell scores.
TIMES = {
    -1: (
        "3.481580353843096",
        "0.08470382309068888",
        "2.527286074151944",
        "0.9180658384326942",
        "0.36171016832955605",
        "2.4420733790608753",
        "0.34317132800597755",
        "2.085258834868358",
        "5.52703910148635",
        "5.762677725283091",
    ),
    1: (
        "1.5627757964054876",
        "4.693716067189554",
        "4.853763433988885",
        "4.378840962048423",
        "0.7500571121174973",
        "1.798701647738304",
        "1.3437943721298375",
        "4.092046047648764",
        "1.7951511905472524",
        "2.159746203677181",
    ),
    2: (
        "0.5506019433871812",
        "0.4886399259440983",
        "5.90670950050828",
        "0.8417048133820738",
        "0.05215114937267621",
        "0.5097750139626672",
        "1.0877318036091592",
        "1.7567549770102098",
        "0.3780644892132115",
        "0.9821779803470796",
    ),
}


@lru_cache(maxsize=None)
def guards(port):
    board = circuit.carrier.federation(0)
    available = [
        p
        for p in board.neighbours[port]
        if port // 12 == p // 12 and p % 12 not in (0, 1, 2, 3, 5, 7)
    ]
    if len(available) < 2:
        raise ValueError("Each requested port needs two idle internal neighbours")
    return tuple(available[:2])


def native_phase(state, p, angle):
    """Replay in the sparse occupation circuit. Unrounded leakage is bounded below."""
    import numpy as np

    turns = round(angle / (np.pi / 2))
    g = guards(p)
    for j, decimal in enumerate(TIMES[turns]):
        state = circuit.contact(state, (p, g[j % 2]), circuit.pulse(float(decimal)))
    return state


def phase_bound():
    """Outward interval arithmetic, using exact decimals and even/odd seam eigenvalues."""
    prev = iv.dps
    try:
        iv.dps = 45
        upper = {}
        for turns, times in TIMES.items():
            v = [iv.mpc(1), iv.mpc(0), iv.mpc(0)]
            for n, t in enumerate(times):
                j = 1 + n % 2
                tau = iv.mpf(t)
                even = (v[0] + v[j]) * iv.exp(-8j * tau) / 2
                odd = (v[0] - v[j]) * iv.exp(-9j * tau) / 2
                v[0], v[j] = even + odd, even - odd
            target = {-1: -1j, 1: 1j, 2: -1}[turns]
            v[0] -= target
            norm = iv.sqrt(sum(abs(a) ** 2 for a in v))
            if not norm.b < iv.mpf("1e-13"):
                raise ArithmeticError("Native pulse-word gate exceeds its error budget")
            upper[turns] = float(norm.b)
        return upper
    finally:
        iv.dps = prev
