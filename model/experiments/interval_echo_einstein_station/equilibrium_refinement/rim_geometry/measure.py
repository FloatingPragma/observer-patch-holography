"""Score equal-width and independently clocked apertures on the five measured rim stations."""
from fractions import Fraction
from itertools import combinations
import json
import math
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo.clock import Registers
from model.executable_model import federation
from ..clock_ball import CENTER, RADIUS, read_word_menu, shortest_chains
from ..measure import apparatus
from ...proper_ball_ruler.measure import pulse_calibration


SHOTS = 2**18
SEED = 710162


def arc_classes(count):
    """Pair incidence averaged over rotations of a consecutive run."""
    rows = []
    for k in range(count + 1):
        pairs = list(combinations(range(k), 2))
        near = sum(min(b - a, count - b + a) == 1 for a, b in pairs)
        rows.append(dict(visible=k, near=near, far=len(pairs) - near))
    return rows


def sharp_arc_bound():
    classes = arc_classes(5)
    target = math.log(5 / 2) / math.log(5)
    # All possible angular quadratures lie on the segments between consecutive rows.
    vertices = [(row["near"] / row["visible"], row["far"] / row["visible"])
                for row in classes[1:]]
    candidates = []
    for left, right in zip(vertices, vertices[1:]):
        slope = tuple(b - a for a, b in zip(left, right, strict=True))
        mix = [0., 1.]
        # Changes of the active absolute deviation and both sign crossings.
        for intercept, derivative in [(left[0] - target, slope[0]),
                                      (left[1] - target, slope[1]),
                                      (sum(left) - 2 * target, sum(slope)),
                                      (left[0] - left[1], slope[0] - slope[1])]:
            if derivative and 0 < -intercept / derivative < 1:
                mix.append(-intercept / derivative)
        candidates.extend(max(abs(x + t * dx - target)
                              for x, dx in zip(left, slope, strict=True)) for t in mix)
    bound = math.log(5) * min(candidates)
    analytic_bound = math.log(5) / 2 - 2 * math.log(5 / 2) / 3
    if not math.isclose(bound, analytic_bound, abs_tol=1e-14):
        raise ValueError("Angular minimax disagrees with the explicit bound")
    return dict(consecutive_classes=classes, required_overlap_per_single=target,
                sharp_minimum_worst_pair_error_nats=bound,
                attaining_uniform_arc_length=3 / (1.5 - target),
                fitted_full_rim_lottery_probability=target / (5 - 4 * target))


def receiver_probability(ranges, occupied):
    """Retention on each independent routed report uses a fresh half-exchange coin per leg."""
    miss = Fraction(1)
    for p in occupied:
        miss *= 1 - Fraction(1, 2**ranges[p])
    return 1 - miss


def sampled_panel_total(rows, receivers):
    trials = len(rows) * receivers * SHOTS
    point = sum(sum(row["hits"]) for row in rows) * receivers / trials
    # Two-sided Hoeffding for each of three independent averages, union bound 0.01.
    radius = receivers * math.sqrt(math.log(600) / (2 * trials))
    return point, (point - radius, point + radius)


def measure():
    car = federation(0)
    ports, contacts, border, words = apparatus(car, [0])
    rays = shortest_chains(ports, contacts, CENTER)
    interior = [p for p in ports if len(rays[p]) < RADIUS]
    rim = [p for p in ports if len(rays[p]) == RADIUS]
    surveys = {p: shortest_chains(ports, contacts, p) for p in [CENTER, *rim]}
    distances = {p: {q: len(path) for q, path in routes.items()}
                 for p, routes in surveys.items()}
    reg = Registers(car, seed=31709)
    flags = {w: read_word_menu(reg, border, ports, surveys, w) for w in words}
    if ports != list(range(len(ports))) or len(rim) != 5:
        raise ValueError("Measurement requires the original five-rim address table")
    forks = []
    for i in interior:
        if distances[CENTER][i] != 1:
            continue
        ends = [q for q in rim if distances[q][i] == 1]
        for q, r in combinations(ends, 2):
            forks.append(dict(inner=i, rim=[q, r], rim_separation=distances[q][r]))
    halfway = pulse_calibration()["quadrature_probabilities"][0][2]
    if not math.isclose(halfway, 0.5, rel_tol=0, abs_tol=2e-14):
        raise ValueError("Material clock half-exchange coin misses its calibration")
    rng = np.random.default_rng(SEED)
    rows = []
    admitted = set(interior + rim)
    for word in words:
        if not set(word) <= admitted:
            continue
        occupied = [q for q in rim if flags[word][q]]
        probabilities = [receiver_probability(distances[p], occupied) for p in rim]
        counts = [int(rng.binomial(SHOTS, float(x))) for x in probabilities]
        rows.append(dict(word=list(word), rim_occupancy=len(occupied),
                         occupied_rim_separation=distances[occupied[0]][occupied[1]]
                         if len(occupied) == 2 else None,
                         hit_probabilities=[str(x) for x in probabilities], hits=counts,
                         mean_lit_panels=str(sum(probabilities))))
    k_one = [row for row in rows if row["rim_occupancy"] == 1]
    single = {Fraction(row["mean_lit_panels"]) for row in k_one}
    if len(single) != 1:
        raise ValueError("Ordinary clocked aperture response depends on placement")
    gain = math.log(5) / float(single.pop())
    measured_single, single_interval = sampled_panel_total(k_one, len(rim))
    sampled_gain = math.log(5) / measured_single
    positive_control = {}
    for distance in [1, 2]:
        pairs = [row for row in rows if row["occupied_rim_separation"] == distance]
        measured_pair, pair_interval = sampled_panel_total(pairs, len(rim))
        positive_control[str(distance)] = dict(
            placements=len(pairs),
            predicted_pair_capacity_nats=gain * float(Fraction(pairs[0]["mean_lit_panels"])),
            sampled_pair_capacity_nats=sampled_gain * measured_pair,
            simultaneous_99_interval_nats=[math.log(5) * pair_interval[0] / single_interval[1],
                                            math.log(5) * pair_interval[1] / single_interval[0]])
    return dict(center=CENTER, rim=rim, available_geodesic_rows=[
        dict(origin=p, distances=[distances[p][q] for q in ports]) for p in surveys],
        geodesic_forks=forks, arc_width_test=sharp_arc_bound(),
        clocked_aperture=dict(single_coin_probability="1/2", simulation_seed=SEED,
            independent_starts_per_word_and_receiver=SHOTS, analytic_ordinary_gain=gain,
            sampled_ordinary_gain=sampled_gain,
            measured_pair_classes=positive_control, records=rows))


if __name__ == "__main__":
    Path(__file__).with_name("results.json").write_text(
        json.dumps(measure(), indent=2, allow_nan=False) + "\n")
