"""Project ancestral-agreement measurements of two routed physical marks."""
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from math import log
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo.clock import Registers
from model.candidate_laws.interval_echo.projection import Branch
from model.executable_model import federation
from model.experiments.mobile_memory.refinement import _cell_seams, certify_step

MENUS = ("common_parent", "always_share", "independent")
HALF = Fraction(1, 2)
CHILDREN = range(4)


def source_step(r, carry, menu):
    """Take the same five physical readings on every branch and publish the next carry."""
    owner = {0}
    share, a0, a1, b0, b1 = (r.fair(0, owner) for _ in range(5))
    a, b = 2*a0 + a1, 2*b0 + b1
    common = r.read(carry, owner)
    if share and (menu == "always_share" or (menu == "common_parent" and common)):
        b = a
    r.assign(carry, int(common and a == b), owner)
    if r.read(carry, owner) != (common and a == b):
        raise ValueError("incorrect retained prefix-comparison bit")
    return a, b


CARRY = (5, 11, 4)
EMITTERS = (1, 7)


def prepared_registers(board):
    r = Registers(board, seed=0)
    if any(abs(p - float(HALF)) > 1e-12 for p in r.rates.values()):
        raise ValueError("source half-contact is not fair from both sharp inputs")
    r.loads[2] = 2
    r.initialize_digits([CARRY])
    r.loads[list(EMITTERS)] = 2
    return r


def native_kernel(board, equal_parent, menu):
    """Read all equiprobable branches of one decision from both comparison states."""
    if menu not in MENUS:
        raise ValueError("unrecognized source menu")
    output = Counter()
    r = prepared_registers(board)
    for bits in product(range(2), repeat=5):
        r.assign(CARRY, int(equal_parent), {0})
        r.rng = Branch(bits)
        output[source_step(r, CARRY, menu)] += 1
    return output


def mark_and_read(r, owners, left, right):
    """Emit two distinguishable marks, read every terminal port and undo the routing."""
    contacts = (r.route(EMITTERS[0], 12*left, owner=owners)
                + r.route(EMITTERS[1], 12*right+6, owner=owners))
    for a, b in contacts:
        r.swap(a, b, owners)
    readings = [[], []]
    for face in owners:
        vacancy = 12*face+10
        for marker, rail in enumerate((12*face, 12*face+6)):
            if r.flag(rail, vacancy, {face}):
                readings[marker].append(face)
    for a, b in reversed(contacts):
        r.swap(a, b, owners)
    if readings != [[left], [right]]:
        raise ValueError("local terminal records disagree with the emitted address pair")
    return readings[0][0], readings[1][0]


def complete_branch_records(board, depth, menu):
    """Carry the measured prefixes, flag and marks through one physical experiment."""
    r = prepared_registers(board)
    owners = range(4**depth)
    counts = Counter()
    for bits in product(range(2), repeat=5*depth):
        r.rng = Branch(bits)
        r.assign(CARRY, 1, {0})
        left = right = 0
        for index in range(1, depth+1):
            a, b = source_step(r, CARRY, menu)
            stride = 4**(depth-index)
            prepared = r.loads.copy()
            x, y = mark_and_read(r, owners, (4*left+a)*stride, (4*right+b)*stride)
            if not np.array_equal(r.loads, prepared):
                raise ValueError("mark readback failed to restore the other material registers")
            left, right = x//stride, y//stride
            if (left == right) != r.read(CARRY, {0}):
                raise ValueError("measured addresses disagree with the local carry")
        counts[left, right] += 1
    total = sum(counts.values())
    return Counter({k: Fraction(n, total) for k, n in counts.items()}), {
        "distinct_ordered_pair_readbacks": len(counts), "complete_histories": total,
        "half_contacts": r.gates["half_contacts"],
        "directed_intercell_swaps": r.gates["intercell_pi_contacts"]}


def extend(table, kernels):
    result = Counter()
    for (x, y), mass in table.items():
        for (a, b), count in kernels[x == y].items():
            result[4*x+a, 4*y+b] += mass * Fraction(count, 32)
    return result


def push(table):
    result = Counter()
    for (a, b), weight in table.items():
        result[a//4, b//4] += weight
    return result


def kl(table, reference):
    return sum(float(p)*log(float(p/reference[k])) for k, p in table.items() if p)


def stratum_optimum(depth):
    """Solve the projection directly from stratum masses and uniform counting."""
    size = 4**depth
    groups = {}
    for x, y in product(range(size), repeat=2):
        common = depth
        a, b = x, y
        while a != b:
            common -= 1
            a, b = a//4, b//4
        groups.setdefault(common, []).append((x, y))
    s = Fraction(5, 8)
    output = Counter()
    for common, pairs in groups.items():
        mass = s**common * (1 if common == depth else 1-s)
        for pair in pairs:
            output[pair] = mass/len(pairs)
    return output


def agreement(table):
    return sum(p for (x, y), p in table.items() if x == y)


def cover_weights(board, depth):
    """Stop at the first fair 'one', with the last resolution taking the residual weight."""
    r = prepared_registers(board)
    counts = Counter()
    for bits in product(range(2), repeat=depth-1):
        r.rng = Branch(bits)
        chosen = depth
        for index in range(1, depth):
            if r.fair(0, {0}):
                chosen = index
                break
        counts[chosen] += 1
    total = sum(counts.values())
    return [Fraction(counts[j], total) for j in range(1, depth+1)]


def run():
    boards = [federation(j) for j in range(3)]
    for coarse, fine in zip(boards, boards[1:]):
        certify_step(coarse.cells, _cell_seams(coarse), fine.cells, _cell_seams(fine))
    tables = {menu: Counter({(0, 0): Fraction(1)}) for menu in MENUS}
    receipt = []
    prior_weights = []
    s = Fraction(5, 8)
    one_step_cost = float(s)*log(2.5) + float(1-s)*log(.5)
    for depth, board in enumerate(boards[1:], 1):
        gate = {menu: {equal: native_kernel(board, equal, menu)
                       for equal in (True, False)} for menu in MENUS}
        fine = {menu: extend(tables[menu], gate[menu]) for menu in MENUS}
        for menu in MENUS:
            if push(fine[menu]) != tables[menu]:
                raise ValueError("source address tables violate exact pushforward")
        native, control, ref = (fine[key] for key in MENUS)
        if native != stratum_optimum(depth) or set(ref.values()) != {Fraction(1, 16**depth)}:
            raise ValueError("source or reference disagrees with the finite projection")
        if agreement(native) != s**depth or agreement(control) != s**depth:
            raise ValueError("source and changed wiring do not satisfy the same constraints")
        # The prior changed-wiring distribution is a concrete, non-optimal feasible table
        # from depth two onwards. Test the same lift on a visibly nonsymmetric competitor.
        coarse_competitor = Counter(tables["common_parent"])
        if depth > 1:
            delta = coarse_competitor[0, 1]/2
            coarse_competitor[0, 1] -= delta
            coarse_competitor[0, 2] += delta
        lifted = extend(coarse_competitor, gate["common_parent"])
        if push(lifted) != coarse_competitor or agreement(lifted) != s**depth:
            raise ValueError("canonical lift lost a coarse feasible table or constraint")
        offset = kl(lifted, ref) - kl(coarse_competitor, tables["independent"])
        gap = kl(control, ref) - kl(native, ref)
        predicted = one_step_cost * sum(1-float(s**k) for k in range(depth))
        if abs(offset-float(s**(depth-1))*one_step_cost) > 1e-12 or abs(gap-predicted) > 1e-12:
            raise ValueError("finite KL chain or changed-wiring gap failed")
        replay = {}
        for menu in MENUS:
            readbacks, replay[menu] = complete_branch_records(board, depth, menu)
            if readbacks != fine[menu]:
                raise ValueError("same-apparatus joint histories disagree with the projected tables")
        weights = cover_weights(board, depth)
        if prior_weights and weights[:-2] + [sum(weights[-2:])] != prior_weights:
            raise ValueError("observer-selection weights do not coarsen")
        prior_weights = weights
        receipt.append({
            "cover_weights_coarse_to_fine": [str(w) for w in weights],
            "depth": depth, "sphere_cells": board.cells,
            "marked_trials_by_menu": replay, "ancestor_equality": str(agreement(native)),
            "native_KL_nats": round(kl(native, ref), 12), "changed_wiring_KL_excess_nats": round(gap, 12),
            "feasible_lift_cost_nats": round(offset, 12),
            "half_contact_Born_branches_per_context": sum(gate["common_parent"][True].values()),
            "same_parent_diagonal_counts": [gate["common_parent"][True][a, a] for a in CHILDREN],
            "separated_parent_diagonal_counts": [gate["common_parent"][False][a, a] for a in CHILDREN],
            "uniform_reference_mass": str(next(iter(ref.values()))),
        })
        tables = fine
    return {"root_face": 0, "child_choices": 4, "source_share_or_fresh": "1/2",
            "conditional_child_agreement": str(s), "levels": receipt}


if __name__ == "__main__":
    result = run()
    if result != json.loads(Path(__file__).with_name("results.json").read_text()):
        raise SystemExit("the native readout differs from the registered receipt")
    print(json.dumps(result, indent=2))
