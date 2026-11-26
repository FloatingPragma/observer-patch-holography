"""Screen a common regional response using simultaneous interval data.

Run this package module with PYTHONPATH=.:code. The two unknown
constants are an overall positive slope and a free common offset. Bounds
are exact for the supplied rectangular intervals; correlated inputs can
only narrow that set. A compatible slope is a necessary finite control.
"""
from fractions import Fraction as F
from itertools import product
import json
import math
from pathlib import Path
import mpmath as mp

from .. import (
    experiment as survey,
    finite_rods,
)
from . import (
    experiment as counted,
)
from ..archive import source_densities


def slope_window(rows):
    """For each reported x/y rectangle, ask if one line y=kx+b hits every rectangle."""
    lower, upper = F(0), None
    lower_pair, upper_pair = None, None
    for (i, (x_min, x_max, y_min, y_max)), (j, other) in product(rows.items(), repeat=2):
        assert x_min <= x_max and y_min <= y_max
        delta_x, delta_y = x_max - other[0], y_min - other[3]
        if delta_x == 0:
            if delta_y > 0:
                return None, (i, j, float(delta_y))
            continue
        cutoff = delta_y / delta_x
        if delta_x > 0 and cutoff > lower:
            lower, lower_pair = cutoff, (i, j)
        elif delta_x < 0 and (upper is None or cutoff < upper):
            upper, upper_pair = cutoff, (i, j)
    if upper is not None and (upper < lower or upper <= 0):
        return None, (lower_pair, upper_pair)
    return (float(lower), float(upper) if upper is not None else None), None


def phase_interval(bank, blank, cfg):
    """Linear phase response inferred directly from amplified counts, including blank uncertainty."""
    with mp.workdps(70):
        def phase_limits(record):
            freq = mp.mpf(record["ones"]) / (record["ones"] + record["zeros"])
            radius = cfg["probability_radius"]
            width = mp.mpf(radius.numerator) / radius.denominator
            theta = mp.pi / 6
            lo, hi = [mp.mpf(1)/4 + mp.sin(theta + t)**2 / 2 for t in (-.125, .125)]
            if freq + width < lo or freq - width > hi:
                raise ValueError("count bank is incompatible with the enrolled phase range")
            def phase(x):
                x = max(lo, min(hi, x))
                return mp.asin(mp.sqrt(2*x - mp.mpf(1)/2)) - theta
            return phase(freq - width), phase(freq + width)
        low, high = phase_limits(bank)
        base_low, base_high = phase_limits(blank)
        gain = (mp.log(mp.mpf(5)/3) * mp.sqrt(3)/4
                * cfg["inverse_radius"]**4 / cfg["coherent_repeats"])
        return (F(math.nextafter(float((low-base_high)*gain), -math.inf)),
                F(math.nextafter(float((high-base_low)*gain), math.inf)))


def source_controls():
    receipt = json.loads((Path(counted.__file__).parent / "results.json").read_text(),
                         parse_float=F)
    ordinary, collar = (F(v).limit_denominator(33) for v in source_densities()["held"])
    stock = {"ordinary": (ordinary, ordinary), "held": (ordinary + collar, ordinary + collar),
             "half_defect": (ordinary + collar / 2, ordinary + collar),
             "blind_pulse": (ordinary + collar, ordinary + collar)}
    # Opposite axial directions share physical shapes on this static source.
    acquisition = list(dict.fromkeys(counted.FRAMES.values()))
    labels = list(counted.FRAMES.values())
    entropy_windows, phase_windows = [], []
    for level in receipt["levels"]:
        d = level["design"]["inverse_radius"]
        rows, phase_rows = {}, {}
        for name, (rod_density, matter_density) in stock.items():
            intervals = level["comparisons"][name]["ten_entropy_intervals_over_ell4"]
            for shape in acquisition:
                u = [F(v).limit_denominator(3) for v in survey._frames()[shape]]
                ell = F(1, d)
                rod = finite_rods.read_ball(finite_rods.angle_gain(ell) * rod_density, ell, u)
                # The 4*pi prefactor common to every regional work integral is absorbed in k.
                x = (rod.weighted_over_4pi * d**4 * matter_density
                     * (1 + u[1]**2 + u[2]**2))
                y_min, y_max = intervals[labels.index(shape)]
                rows[name, shape] = (x, x, y_min, y_max)
                banks = level["banks"][shape]
                phase_rows[name, shape] = (x, x, *phase_interval(banks[name], banks["blank"], level["design"]))
        ordinary_and_held = {key: v for key, v in rows.items() if key[0] in ("ordinary", "held")}
        allowed, _ = slope_window(ordinary_and_held)
        phase_allowed, _ = slope_window({k: phase_rows[k] for k in ordinary_and_held})
        assert allowed is not None and phase_allowed is not None
        entropy_windows.append(allowed)
        phase_windows.append(phase_allowed)
        # Isolate two acquisitions with the same input and different recorded outputs.
        shape = acquisition[0]
        blind = {key: rows[key] for key in (("held", shape), ("blind_pulse", shape))}
        rejected, witness = slope_window(blind)
        assert rejected is None and witness is not None and len(witness) == 3
        print(f"inverse radius {d}: entropy {allowed}; phase {phase_allowed}; "
              f"matched-input {shape} gap to blind control {witness[2]:.7f}")
    def intersection(windows):
        return max(lo for lo, _ in windows), min(hi for _, hi in windows)
    lo, hi = intersection(entropy_windows)
    assert lo > hi
    lo, hi = intersection(phase_windows)
    assert lo < 8 * math.pi**2 < hi
    print(f"three-level phase slope window [{lo:.8f}, {hi:.8f}]; 8*pi²={8*math.pi**2:.8f}")


if __name__ == "__main__":
    # A positive slope of two exists with unknown intercept; four corners
    # allow an exact calibration without a least-squares or tensor fit.
    assert slope_window({"a": (F(0), F(0), F(1), F(1)),
                         "b": (F(2), F(2), F(5), F(5))})[0] == (2., 2.)
    assert slope_window({"a": (F(0), F(0), F(0), F(0)),
                         "b": (F(1), F(1), F(-1), F(-1))})[0] is None
    source_controls()
