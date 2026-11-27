"""Ten-frame interpolation and a finite comoving lattice of native endpoint memories."""

import argparse
from decimal import localcontext
from fractions import Fraction as F
from itertools import permutations, product
import json
import math
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo_partition_repair.published_strings.memories import Pair
from ..body import archive_bodies
from ..survey import locate, mass_tally, read_ball


HERE = Path(__file__).resolve().parent
REST = (F(1), F(0), F(0), F(0))
AXIS = (F(5, 3), F(4, 3), F(0), F(0))
FACE = (F(3), F(2), F(2), F(0))
UNENROLLED = (F(2), F(1), F(1), F(1))
DIRECTIONS = tuple(tuple(b for b in product((-1, 0, 1), repeat=3)
                         if sum(t*t for t in b) == k and next(t for t in b if t) > 0)
                   for k in (1, 2, 3))
FRAMES = (REST, *((AXIS[0], *(sign * AXIS[1] if j == i else F(0) for j in range(3)))
                 for sign in (1, -1) for i in range(3)),
          *((FACE[0], *(F(0) if j == i else F(2) for j in range(3))) for i in (2, 1, 0)))
CENTERS = ((F(0), F(0), F(0)), (F(7, 1021), F(-5, 1021), F(3, 1021)),
           (F(81, 1021), F(-43, 1021), F(31, 1021)))
LENGTH = F(1, 8)


def kernel_interval(frame, b):
    """A rational enclosure of gamma/sqrt(1-(v.n)^2) on an unoriented lattice line."""
    gamma, *p = frame
    squared = 1 - sum(x*y for x, y in zip(p, b))**2 / (gamma**2 * sum(y*y for y in b))
    scale = 10**24
    lo = math.isqrt(scale**2 * squared.denominator // squared.numerator)
    return gamma * lo / scale, gamma * (lo + 1) / scale


def transfer_row(frame):
    return [tuple(sum(kernel_interval(frame, b)[k] for b in group) / len(group) for k in (0, 1))
            for group in DIRECTIONS]


def multiply(x, y):
    values = [a*b for a in x for b in y]
    return min(values), max(values)


def determinant(matrix):
    low = high = 0
    for perm in permutations(range(3)):
        sign = sum(perm[i] > perm[j] for i in range(3) for j in range(i+1, 3)) % 2
        term = multiply(multiply(matrix[0][perm[0]], matrix[1][perm[1]]), matrix[2][perm[2]])
        if sign:
            low, high = low - term[1], high - term[0]
        else:
            low, high = low + term[0], high + term[1]
    return low, high


def calibration():
    frames = (REST, AXIS, FACE)
    matrix = [transfer_row(frame) for frame in frames]
    det = determinant(matrix)
    if det[0] <= 0:
        raise ValueError("cubic angular responses do not separate the three calibration classes")
    weights = []
    for k in range(3):
        replaced = [[(f[0]**2, f[0]**2) if j == k else entry for j, entry in enumerate(row)]
                    for f, row in zip(frames, matrix)]
        weights.append(multiply(determinant(replaced), (1/det[1], 1/det[0])))
    return det, weights


def relative_factor(frame, weights):
    products = [multiply(w, q) for w, q in zip(weights, transfer_row(frame))]
    return sum(p[0] for p in products), sum(p[1] for p in products)


def gate(cut_entropy):
    if not cut_entropy:
        return Pair(0)
    low, high = -0.02, 0.02
    if not Pair(low).cut_entropy_change <= cut_entropy <= Pair(high).cut_entropy_change:
        raise ValueError("material exposure outside the calibrated monotone endpoint branch")
    for _ in range(100):
        middle = (low + high) / 2
        if Pair(middle).cut_entropy_change < cut_entropy:
            low = middle
        else:
            high = middle
    return Pair((low + high) / 2)


def lattice_points(den, center, radius, frame):
    """Count the projected tangent ellipsoid using an exact integer membership test."""
    gamma, *p = frame
    unit = math.lcm(den, radius.denominator, *(x.denominator for x in center))
    scale = math.lcm(*(x.denominator for x in frame))
    g = int(scale * gamma)
    momenta = [int(scale * v) for v in p]
    axes = []
    for c, v in zip(center, p):
        extent = float(radius) * math.sqrt(1 + float(v)**2)
        axes.append(np.arange(math.floor(den*(float(c) - extent))-1,
                              math.ceil(den*(float(c) + extent))+2, dtype=np.int64))
    maximum = max(max(abs(int(a[0])*(unit//den) - int(unit*c)),
                      abs(int(a[-1])*(unit//den) - int(unit*c))) for a, c in zip(axes, center))
    bound = max(3*g*g, sum(abs(v) for v in momenta)**2) * maximum**2
    if max(bound, (int(radius*unit)*g)**2) >= 2**63:
        raise ValueError("finite replay exceeds its exact 64-bit membership range")
    points = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    offset = points*(unit//den) - [int(unit*c) for c in center]
    limit = int(radius*unit) * g
    dot = sum(momenta[j]*offset[:, j] for j in range(3))
    inside = g*g*(offset*offset).sum(axis=1) - dot*dot <= limit*limit
    return points[inside]


def complete_pairs(points, b, steps):
    """Sum joint inclusions over all short separations on one primitive lattice direction."""
    pivot = next(i for i, v in enumerate(b) if v)
    if np.ptp(points[:, pivot]) > steps:
        raise ValueError("the ball must fit inside the separation window")
    transverse = [points[:, j] - b[j]*points[:, pivot] for j in range(3) if j != pivot]
    x, y = (v-v.min() for v in transverse)
    keys = x*(int(y.max())+1) + y
    counts = np.bincount(keys)
    order = np.lexsort((points[:, pivot], keys))
    same_line = np.diff(keys[order]) == 0
    adjacent = int(np.count_nonzero(same_line & (np.diff(points[order, pivot]) == 1)))
    return int((counts*(counts-1)).sum()//2), adjacent


def read_bank(den, center, radius, frame, density, weights):
    """Aggregate an explicitly finite common array over the complete and cut native pairs."""
    points = lattice_points(den, center, radius, frame)
    bulk, divergence = [], []
    full_by_orbit = []
    for index, (dirs, weight) in enumerate(zip(DIRECTIONS, weights), 1):
        step = math.sqrt(index)/den
        short_steps = math.floor(float(LENGTH)/step)
        signed = -(4*math.pi/15)*density*weight*step/(den**3*len(dirs))
        # The first positive separation represents its Voronoi interval [0, 3 step/2].
        short, first, long = gate(signed), gate(1.5*signed), gate(-(short_steps+0.5)*signed)
        counts = [complete_pairs(points, b, short_steps) for b in dirs]
        full, adjacent = (sum(values) for values in zip(*counts))
        one = 2*len(points)*len(dirs)
        later_full, first_cut = full-adjacent, one-2*adjacent
        later_cut = one*(short_steps-1) - 2*later_full
        bulk.extend((short.cut_entropy_change*later_cut, first.cut_entropy_change*first_cut,
                     long.cut_entropy_change*one))
        divergence.extend((short.cut_relative_entropy*later_cut, short.whole_modular_change*later_full,
                           first.cut_relative_entropy*first_cut, first.whole_modular_change*adjacent,
                           long.cut_relative_entropy*one))
        full_by_orbit.append([full, adjacent])
    return {"address_count": len(points), "complete_short_pairs": full_by_orbit,
            "entropy_over_radius4": math.fsum(bulk)/float(radius)**4,
            "relative_entropy_over_radius4": math.fsum(divergence)/float(radius)**4}


def separate_blank_control(mass):
    """A single rest-frame lattice address enrolled by only one equal-volume survey."""
    projected = []
    with localcontext() as ctx:
        ctx.prec = 65
        for stock in (F(0), mass):
            ball = read_ball(stock, F(1), 64)
            _, x = locate(ball.sample_timestamps, ball.measured_h, ball.sample_tick)
            projected.append(float(x[2]))
    distance = F(156252, 10**7)
    if not projected[0] < distance < projected[1]:
        raise ValueError("the witnessed comoving address must cross the measured rim")
    registers = 2*sum(len(dirs)*(math.floor(float(LENGTH)*256/math.sqrt(i))+1)
                      for i, dirs in enumerate(DIRECTIONS, 1))
    cut_entropy = -(3*math.log(3/8)+5*math.log(5/8))/8
    whole_entropy = -(5*math.log(5/8)+3*math.log(1/8))/8
    return {"radius": "1/64", "mesh_denominator": 256, "witness_distance": str(distance),
            "flat_projected_radius": projected[0], "held_projected_radius": projected[1],
            "blank_endpoints_at_witness": registers,
            "blank_entropy_increase_lower": registers*min(cut_entropy, whole_entropy-cut_entropy)}


def measure():
    det, enclosures = calibration()
    weights = [float((a+b)/2) for a, b in enclosures]
    hold = relative_factor(UNENROLLED, enclosures)
    held = archive_bodies()["held"]
    tickets = [mass_tally(F(m).limit_denominator(33), 128) for m in (held.ordinary, held.collar)]
    shell_mass = sum(F(outer-inner, unit) for unit, inner, outer in tickets)
    density = float(3*shell_mass/(4*(F(1, 2)**3 - F(1, 4)**3))) / math.pi
    surveys = []
    for den in (256, 512, 1024):
        for center in CENTERS:
            for radius in (F(1, 64), F(3, 256)):
                rows = []
                for u in (*FRAMES, UNENROLLED):
                    reading = read_bank(den, center, radius, u, density, weights)
                    coefficient = float((sum(relative_factor(u, enclosures)))/2)
                    rows.append({"frame": list(map(str, u)), **reading,
                                 "continuum_factor": coefficient, "dust_factor": float(u[0]**2)})
                surveys.append({"mesh_denominator": den, "radius": str(radius),
                                "center": list(map(str, center)), "frames": rows})
    return {"schema": "oph.dust.homogeneous_pairs.v1", "determinant": float(sum(det)/2),
            "orbit_weights": weights, "unenrolled_factor": float(sum(hold)/2),
            "weighed_shell_tickets": [list(map(int, t)) for t in tickets],
            "held_rest_density": density, "surveys": surveys,
            "separately_enrolled_blank_control": separate_blank_control(F(64, 7)*shell_mass)}


def verify():
    exact = measure()
    saved = json.loads((HERE/"results.json").read_text())
    def check(value, expected):
        if isinstance(value, float):
            return math.isclose(value, expected, abs_tol=2e-12, rel_tol=2e-12)
        if isinstance(value, list):
            return len(value) == len(expected) and all(check(x,y) for x,y in zip(value,expected))
        if isinstance(value, dict):
            return value.keys() == expected.keys() and all(check(v,expected[k]) for k,v in value.items())
        return value == expected
    if not check(exact, saved):
        raise ValueError("homogeneous pair records changed")
    return exact


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    if parser.parse_args().write:
        (HERE/"results.json").write_text(json.dumps(measure(), indent=2) + "\n")
    else:
        verified = verify()
        print("Verified", len(verified["surveys"]), "overlapping surveys and the unenrolled boost.")
