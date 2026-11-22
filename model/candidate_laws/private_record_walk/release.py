"""Coherent sharp-sign and odd-sign releases on the unweighted private graph."""
import argparse
import json
from math import comb, isqrt
from pathlib import Path

import numpy as np
from scipy import sparse

from model.executable_model import federation
from model.experiments.interval_echo_consensus.bounded_clock_flight import DetectorEvolution
from model.experiments.interval_echo_encounter.check_transport import checked_steps
from .probe import pair_source

HERE = Path(__file__).parent


def word_transform(width):
    """Amplitude for each individual sharp output word, grouped by its Hamming weight."""
    multiplicity = np.array([comb(width, m) for m in range(width+1)])
    transform = np.zeros((width+1, width+1), dtype=np.int64)
    for r in range(width+1):
        for m in range(width+1):
            transform[r, m] = sum((-1)**j * comb(r, j) * comb(width-r, m-j)
                                  for j in range(max(0, m+r-width), min(m, r)+1))
    if not np.array_equal(transform.T @ (multiplicity[:, None]*transform), 2**width*np.diag(multiplicity)):
        raise ArithmeticError("Incomplete sign-character transform")
    return transform / 2**width, multiplicity


def joint_read(sectors, home):
    """All home seams and exterior quiet, resolved by the final tagged word."""
    transform, counts = word_transform(len(sectors)-1)
    amplitudes = transform @ sectors
    home_mask = np.zeros(sectors.shape[1], dtype=bool)
    home_mask[home] = True
    ports = abs(amplitudes[:, home])**2
    quiet = (abs(amplitudes[:, ~home_mask])**2).sum(axis=1)
    return counts[:, None]*np.column_stack((ports, quiet))


def fourth_derivative_bound(h, launch):
    """Bound the fourth derivative of any projection probability with integer launch moments."""
    if (h.shape[0] != h.shape[1] or h.dtype.kind not in "fiu" or not np.isfinite(h.data).all()
            or np.any(h.data < -128) or np.any(h.data > 128)):
        raise ValueError("Expected bounded real material entries")
    twice = (2*h.astype(float)).tocsr()
    if (not np.array_equal(twice.data, np.rint(twice.data)) or (twice != twice.T).nnz):
        raise ValueError("Expected the bounded real symmetric half-integer generator")
    a = twice.astype(np.int64) - int(twice[launch, launch])*sparse.eye(h.shape[0], dtype=np.int64, format="csr")
    radius = int(abs(a).sum(axis=1).max())
    if radius**4 >= 2**63:
        raise ArithmeticError("Insufficient integer range for the fourth moment")
    v = np.zeros(h.shape[0], dtype=np.int64)
    v[launch] = 1
    bounds = [1]
    for k in range(1, 5):
        v = a @ v
        square = sum(int(n)**2 for n in v[v != 0])
        bounds.append((isqrt(square)+1)/2**k)
    return sum(comb(4, k)*bounds[k]*bounds[4-k] for k in range(5))


def certified_window(curve, slopes, fourth, grid, allowance):
    """Range of the cubic Hermite interpolant by subdivision of its Bernstein ordinates."""
    if (curve.ndim != 1 or curve.size < 2 or slopes.shape != curve.shape
            or not np.isfinite(curve).all() or not np.isfinite(slopes).all()
            or not np.isfinite(fourth) or fourth < 0 or type(grid) is not int or grid <= 0
            or not np.isfinite(allowance) or not 0 <= allowance < 1
            or curve.min() < -allowance or curve.max() > 1+allowance):
        raise ValueError("Expected finite detector samples, slopes and error bounds")
    h = 1/grid
    # Each interval's cubic has endpoint values and derivatives from the detector.
    coefficients = np.array([curve[:-1], curve[:-1]+h*slopes[:-1]/3,
                             curve[1:]-h*slopes[1:]/3, curve[1:]]).T
    for _ in range(3):
        a, b, c, d = coefficients.T
        ab, bc, cd = (a+b)/2, (b+c)/2, (c+d)/2
        left, right = (ab+bc)/2, (bc+cd)/2
        mid = (left+right)/2
        coefficients = np.vstack((np.column_stack((a, ab, left, mid)), np.column_stack((mid, right, cd, d))))
    remainder = fourth*h**4/384
    lower = max(0., float(coefficients.min())-remainder-allowance)
    upper = min(1., float(curve.min())+allowance)
    verdict = "certified" if lower > .94 else "refuted" if upper < .94 else "undecided"
    return {"grid_min": float(curve.min()), "min_grid_exposure": int(np.argmin(curve))/grid,
            "hermite_remainder": remainder, "continuous_lower": lower, "minimum_upper": upper, "gate": verdict}


def experiment():
    spec = json.loads((HERE/"protocol.json").read_text(encoding="utf-8"))
    board = federation(spec["board_level"])
    edge = tuple(spec["launch"])
    times = spec["exposures"]
    pair, odd, common, active, width, rest, bits = pair_source(board, edge)
    launch = pair.index(*edge)
    home_seams = [(a, b) for a, b in board.seams if a//12 == b//12 == edge[0]//12]
    home = [pair.index(*s) for s in home_seams]
    initial = np.zeros(pair.P, complex)
    initial[launch] = 1
    grid = 16
    ticks = range(1, grid*max(times)+1)
    appointments = {t: [] for t in times}
    curves = []
    slopes = []
    fourths = []
    max_rate = 0.
    budget = 0.
    for m in range(width+1):
        # m tagged signs are odd. Every other private sign is odd.
        h = odd - sparse.diags((width-m)*active, dtype=float, format="csr")
        home_rows = h[home, :]
        max_rate = max(max_rate, float(abs(h).sum(axis=1).max()))
        fourths.append(fourth_derivative_bound(h, launch))
        current = np.empty(len(ticks)+1)
        slope = np.empty_like(current)
        current[0], slope[0] = 1., 0.
        steps = checked_steps(h/grid, initial, ticks, degree=24, substeps=1, center=0)
        for expected_tick, (tick, state, error) in zip(ticks, steps, strict=True):
            if tick != expected_tick or not np.isfinite(error) or error < 0:
                raise ArithmeticError("Invalid propagation step or error enclosure")
            budget = max(budget, error)
            local = state[home]
            current[tick] = np.vdot(local, local).real
            slope[tick] = 2*np.vdot(local, home_rows @ state).imag
            if tick % grid == 0 and tick//grid in appointments:
                appointments[tick//grid].append(state.copy())
        curves.append(current)
        slopes.append(slope)
    if not np.isfinite(budget) or budget >= 1e-7:
        raise ArithmeticError("Propagation failed its probability error budget")
    # Nine-term Walsh sums, multiplicities and the final probability sums.
    # The derivation of this additional allowance is in the candidate specification.
    read_error = budget + 2e-10
    curves = np.array(curves)
    slopes = np.array(slopes)
    odd_home = curves[-1]
    weights = np.array([comb(width, m)/2**width for m in range(width+1)])
    sharp_home = weights @ curves
    allowance = read_error*(1+2*max_rate/(3*grid))
    windows = [certified_window(c, s, f, grid, allowance)
               for c, s, f in zip(curves, slopes, fourths, strict=True)]
    reads = {}
    for t, states in appointments.items():
        joint = joint_read(np.array(states), home)
        # The binomial mixture is the unread-sign marginal, not the joint instrument.
        np.testing.assert_allclose([joint.sum(), joint[:, :-1].sum()],
                                   [1, sharp_home[grid*t]], atol=1e-8, rtol=0)
        reads[t] = {"weight_and_home_or_quiet": joint.tolist(), "changed_word": float(joint[1:].sum()),
                    "sharp_home": float(joint[:, :-1].sum()), "odd_home": float(odd_home[grid*t])}
    # OR agrees on every single. On a pair it counts a shared record once.
    rivals = {}
    for name, h in (("either_occupied", odd+sparse.diags(common, dtype=float)),
                    ("all_even", pair.exchange())):
        engine = DetectorEvolution(h.tocsr(), times)
        rivals[name] = {"home_and_quiet": engine.read(launch, home).tolist(), "error": engine.probability_error}
    moving = pair.exchange().getcol(launch).tocoo()
    separating = sum(1 for row in moving.row if active[row])
    if separating != sum(len(board.neighbours[v])-1 for v in edge):
        raise ArithmeticError("The first changed-word paths disagree with the launch geometry")
    return {"launch": edge, "seams": home_seams, "tagged_signs": width, "private_signs": bits,
            "single_load_rest": rest, "material_states": pair.P, "sector_probability_error": budget,
            "joint_probability_error": read_error,
            "window": {"odd": windows[-1], "sharp": certified_window(sharp_home, weights @ slopes,
                                                                      weights @ fourths, grid, allowance),
                       "by_odd_count": windows,
                       "all_tagged_preparations_lower": min(w["continuous_lower"] for w in windows)},
            "first_changed_word_quartic": separating*width/64, "appointments": reads, "controls": rivals}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="compare all regenerated outcomes to results.json")
    args = parser.parse_args(argv)
    actual = json.loads(json.dumps(experiment(), allow_nan=False))
    if args.check:
        expected = json.loads((HERE/"results.json").read_text(encoding="utf-8"))
        check_receipt(actual, expected)
        print("All sharp-word, seam and quiet probabilities, controls and window bounds agree.")
    else:
        print(json.dumps(actual, indent=2, allow_nan=False))


def check_receipt(actual, stored):
    if type(actual) is not type(stored):
        raise ValueError("Receipt value has the wrong type")
    if isinstance(actual, dict):
        if actual.keys() != stored.keys():
            raise ValueError("Missing or extra receipt fields")
        for key in actual:
            check_receipt(actual[key], stored[key])
    elif isinstance(actual, list):
        for a, b in zip(actual, stored, strict=True):
            check_receipt(a, b)
    elif isinstance(actual, float):
        if not np.isfinite(stored) or abs(actual-stored) > 1e-11:
            raise ValueError(f"Expected {actual}, read {stored}")
    elif actual != stored:
        raise ValueError(f"Expected {actual}, read {stored}")
