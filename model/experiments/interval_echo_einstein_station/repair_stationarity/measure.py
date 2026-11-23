"""Regional entropy and seam readings along a feasible overlap-repair variation."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from math import log
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import candidate, projection
from model.candidate_laws.interval_echo.clock import Registers
from model.executable_model import federation

HERE = Path(__file__).parent
OUTCOMES = list(product((0, 1), repeat=3))
SHOTS = 2**22
STEP = Fraction(1, 64)


def digit(owner, shape):
    return tuple(12 * owner + value for value in shape)


def branch(car, owners, ports, bits, order, mode, home=0):
    """All controls use hardware fair contacts, signed pointer reads and owned routing."""
    regs = Registers(car, seed=0)
    if max(abs(rate - .5) for rate in regs.rates.values()) > 1e-12:
        raise ValueError("The enrolled native contact is not fair")
    region = {home, *owners}
    regs.rng = projection.Branch([*bits[:3], *(b for i in order for b in bits[3+2*i:5+2*i])])
    shared = digit(home, (3, 2, 4))
    extra = [digit(home, d) for d in ((5, 11, 4), (6, 7, 10))]
    for owner in region:
        regs.loads[12 * owner + 2] = 2
    regs.initialize_digits(extra)
    for address in extra:
        regs.assign(address, regs.fair(home, {home}), {home})
    regs.fair(home, {home})
    retained = [digit(owner, (6, 7, 10)) for owner in owners]
    for i in order:
        owner = owners[i]
        pair = {home, owner}
        local = {owner}
        scratch = digit(owner, (0, 1, 10))
        regs.initialize_digits([scratch, retained[i]])
        common = regs.read_tick(shared, pair, at_hand=owner)
        source_a, source_b = [regs.read_tick(address, pair, at_hand=owner) for address in extra]
        regs.assign(scratch, regs.fair(owner, local), local)
        second = regs.fair(owner, local)
        independent_error = regs.read(scratch, local) & second
        if mode == "independent":
            mismatch = independent_error
        elif mode == "common" or i == 0:
            mismatch = source_a & source_b
        elif mode == "disjoint":
            mismatch = source_a & (1 - source_b)
        else:
            raise ValueError("Unknown correlation preparation")
        regs.assign(retained[i], common ^ mismatch, local)
        regs.assign(scratch, 0, local)
        regs.assign(digit(owner, (3, 2, 4)), 0, local)
    for address in extra:
        regs.assign(address, 0, {home})
    x = regs.read(shared, {home})
    y, z = [regs.read(address, {owner}) for owner, address in zip(owners, retained)]
    terminal = regs.loads.copy()
    # The two full routes exchange their endpoints, restoring the intermediate rails.
    routes = [regs.route(address[0], target, owner=region)
              for address, target in zip(retained, (12*home, 12*home+1))]
    for route in routes:
        for a, b in route:
            regs.swap(a, b, region)
    mismatch = int(regs.flag(12*home, 12*home+1, region))
    for route in reversed(routes):
        for a, b in reversed(route):
            regs.swap(a, b, region)
    if not np.array_equal(regs.loads, terminal) or mismatch != (y ^ z):
        raise ValueError("The native parity reader failed to restore the three-cell state")
    load = tuple(int(terminal[12*c:12*c+12].sum()) for c in sorted(region))
    local_word = tuple(int(p) for p in np.flatnonzero(terminal))
    cut = tuple(int(terminal[p]) for p in ports)
    return (x, y, z), mismatch, local_word, cut, load


def apparatus(level):
    car = federation(level)
    owners = sorted(b//12 for a, b in car.external if a//12 == 0)[:2]
    region = {0, *owners}
    cut = sorted(p for a, b in car.external if (a//12 in region) != (b//12 in region)
                 for p in (a, b))
    tables, detectors, cuts, loads, unique_words = {}, {}, {}, set(), {}
    for mode in ("independent", "common", "disjoint"):
        tables[mode], detectors[mode], cuts[mode] = Counter(), Counter(), Counter()
        for bits in product((0, 1), repeat=7):
            result = branch(car, owners, cut, bits, (0, 1), mode)
            if result != branch(car, owners, cut, bits, (1, 0), mode):
                raise ValueError("Observer orders give different material readouts")
            row, disagreement, word, boundary, material = result
            tables[mode][row] += 1
            detectors[mode][disagreement] += 1
            cuts[mode][boundary] += 1
            loads.add(material)
            if unique_words.setdefault(row, word) != word:
                raise ValueError("Terminal sharp words retain more than the three digits")
    if len(set(unique_words.values())) != len(OUTCOMES) or len(loads) != 1:
        raise ValueError("The record spectrum or volume proxy changed")
    if cuts["independent"] != cuts["common"] or cuts["common"] != cuts["disjoint"]:
        raise ValueError("The complete measured cut is sensitive to the selected mixture")
    independent = Counter(projection.sample(car, owners, bits, (0, 1))
                          for bits in product((0, 1), repeat=5))
    if any(4*independent[row] != tables["independent"][row] for row in OUTCOMES):
        raise ValueError("The independent control differs from the registered projection")
    optimum = [["3/8", "1/8", "1/8", "3/8"]] * 2
    if any(projection.accessible_tables(table) != optimum for table in tables.values()):
        raise ValueError("The A3-accessible tables differ")
    receipt = dict(cells=car.cells, regional_cells=sorted(region), cut_ports=cut,
                   occupied_region_loads=list(next(iter(loads))), complete_material_word_count=len(unique_words),
                   branches_per_mode_per_order=128, owner_orders=2, accessible_tables=optimum,
                   common_cut_distribution={"".join(map(str, k)): str(Fraction(n, 128))
                                            for k, n in cuts["independent"].items()},
                   joint_word_counts={mode: [table[row] for row in OUTCOMES] for mode, table in tables.items()},
                   routed_disagreement_counts={mode: [table[0], table[1]] for mode, table in detectors.items()})
    return tables, detectors, receipt


def mixed(tables, theta, outcomes):
    # P(both errors) is 1/4 in the common source and zero in the disjoint source.
    weight = 4 * theta
    return tuple((weight*tables["common"][row] + (1-weight)*tables["disjoint"][row])/128
                 for row in outcomes)


def entropy(p):
    return -sum(x*log(x) for x in p if x)


def tangent_report(tables, detectors, theta, seed):
    p = np.array(mixed(tables, theta, OUTCOMES), dtype=float)
    minus, plus = [np.array(mixed(tables, theta + sign*STEP, OUTCOMES), dtype=float)
                   for sign in (-1, 1)]
    tangent = (plus-minus)/(2*float(STEP))
    modular = -np.log(p)
    analytic = log(float((Fraction(1,4)-theta)**2/((Fraction(1,2)+theta)*theta)))
    first_law = float(modular @ tangent)
    np.testing.assert_allclose(first_law, analytic, atol=1e-12)
    # The second central moment is the delta-method variance for the plug-in entropy.
    rng = np.random.default_rng(seed)
    counts = [rng.multinomial(SHOTS, probabilities) for probabilities in (minus, plus)]
    entropy_slope = (entropy(counts[1]/SHOTS) - entropy(counts[0]/SHOTS))/(2*float(STEP))
    errors = [probs @ np.log(probs)**2 - (probs @ np.log(probs))**2 for probs in (minus, plus)]
    exact_step = (entropy(plus)-entropy(minus))/(2*float(STEP))
    qminus, qplus = [np.array(mixed(detectors, theta+sign*STEP, (0, 1)), dtype=float)
                      for sign in (-1, 1)]
    probe_counts = [rng.multinomial(SHOTS, probabilities) for probabilities in (qminus, qplus)]
    probe_slope = (probe_counts[1][1]-probe_counts[0][1])/(SHOTS*2*float(STEP))
    return dict(theta=str(theta), exact_entropy_slope=first_law, finite_entropy_slope=exact_step,
                measured_entropy_slope=float(entropy_slope),
                entropy_slope_standard_error=float(np.sqrt(sum(errors)/SHOTS)/(2*float(STEP))),
                joint_minus_counts=counts[0].tolist(), joint_plus_counts=counts[1].tolist(),
                exact_routed_mismatch_slope=float((qplus[1]-qminus[1])/(2*float(STEP))),
                measured_routed_mismatch_slope=float(probe_slope),
                routed_slope_standard_error=float(np.sqrt((qminus.prod()+qplus.prod())/SHOTS)/(2*float(STEP))),
                routed_minus_counts=probe_counts[0].tolist(), routed_plus_counts=probe_counts[1].tolist())


def run():
    campaigns, calibrations = [], []
    constraints = np.array([[1, x, x ^ y, x ^ z] for x, y, z in OUTCOMES], dtype=float)
    for level in range(3):
        tables, detectors, board = apparatus(level)
        for theta in (Fraction(1, 16), Fraction(1, 8)):
            result = tangent_report(tables, detectors, theta, [632573, level, theta.denominator])
            campaigns.append(dict(cells=board["cells"], **result))
        calibrations.append(board)
    independent = np.array([tables["independent"][row]/128 for row in OUTCOMES])
    log_meter = constraints @ np.array([log(32/9), 0., log(3), log(3)])
    np.testing.assert_allclose(-np.log(independent), log_meter, atol=1e-14)
    # Native stored energy on the temporary comparator seam for its sharp outcomes.
    term = candidate.native_law.seam_term(candidate.KAPPA)
    indices = [candidate.native_law.SRC.index(word) for word in ((0,0), (2,2), (0,2), (2,0))]
    energies = term.diagonal()[indices]
    np.testing.assert_allclose(energies, [0,0,32.5,32.5], atol=1e-13)
    return dict(outcome_order=["".join(map(str, row)) for row in OUTCOMES], preparations_per_setting=SHOTS,
                central_difference_step=str(STEP), feasible_joint_tangent_dimension=8-int(np.linalg.matrix_rank(constraints)),
                equal_rail_and_unequal_rail_seam_energy=[0., float(energies[-1])],
                exact_comparator_energy_slope_per_theta=-2*float(energies[-1]),
                apparatus=calibrations, acquisitions=campaigns)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(), indent=2)+"\n")
