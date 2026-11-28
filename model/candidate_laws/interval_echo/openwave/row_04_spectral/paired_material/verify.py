"""Signed pointer, formation and free-flight comparisons on one ternary material."""

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import cache
from itertools import product
import json

import numpy as np

from .. import raw_verify as raw


LOCAL = tuple(product((0, 1), raw.LOADS))


@cache
def seam():
    echo, exchange, agreement = raw.local()[:3]
    full_words = tuple(product(raw.LABELS, raw.LOADS))
    # Count each species of sharp port load on BOTH ends of every native term.
    for value in (0, 2):
        count = np.array([pair.count(value) for _, pair in full_words])
        for part in (exchange, agreement):
            np.testing.assert_array_equal(count[:, None] * part, part * count)

    conjugate = [raw.ix(None if label is None else (1-label[0], -label[1]), (2-a, 2-b))
                 for label, (a, b) in full_words]
    for part in (echo, exchange, agreement):
        np.testing.assert_array_equal(part[np.ix_(conjugate, conjugate)], part)

    basis = np.eye(99, dtype=int)
    # An addressed echo can spend an installed neutral-load collar with |d|=2.
    for sharp in ((2, 0), (0, 2)):
        j = raw.ix(None, sharp)
        seed = sum(basis[:, raw.ix((c, sharp[0]-sharp[1]), (1, 1))] for c in (0, 1))
        np.testing.assert_array_equal(echo[:, j], 2 * seed)
        np.testing.assert_array_equal(echo @ seed, 2 * basis[:, j])

    # Each arrow is a physical even superposition of the two signed memories.
    arrows = np.stack([basis[:, raw.ix((c, -1), pair)] + basis[:, raw.ix((c, 1), pair)]
                       for c, pair in LOCAL], axis=1)
    np.testing.assert_array_equal(exchange @ arrows, 0 * arrows)
    twice_agreement = arrows.T @ agreement @ arrows // 4
    np.testing.assert_array_equal(agreement @ arrows, 2 * arrows @ twice_agreement)
    twice_flight = twice_agreement - 2 * np.eye(len(LOCAL), dtype=int)
    for coin in (0, 1):
        np.testing.assert_array_equal(twice_flight[:, LOCAL.index((coin, (1, 1)))], 0)

    flag = raw.projector4([label is None for label, _ in full_words])
    twice_hit = arrows.T @ flag @ arrows // 4
    np.testing.assert_array_equal(flag @ arrows, 2 * arrows @ twice_hit)
    for ref, n in product(range(3), repeat=2):
        j = raw.ix(None, (n, ref))
        np.testing.assert_array_equal(agreement[:, j], 4 * (n-ref)**2 * basis[:, j])
    for sign, n in product((-1, 0, 1), range(3)):
        j = raw.ix(None, (n, 1))
        pointer = raw.projector4([label is not None and label[1] == sign for label, _ in full_words])
        np.testing.assert_array_equal(pointer[:, j], 4 * int(sign == n-1) * basis[:, j])
    return twice_flight, twice_hit


@dataclass(frozen=True)
class Word:
    excess: int | None
    vacancy: int | None
    flips: frozenset

    def population(self, port):
        return 1 + int(port == self.excess) - int(port == self.vacancy)


def contact(word, edge, matrix):
    coin = int(edge not in word.flips)  # stock head at the earlier physical port
    pair = tuple(word.population(p) for p in edge)
    entries = matrix[:, LOCAL.index((coin, pair))]
    for index in np.flatnonzero(entries):
        end_coin, end_pair = LOCAL[index]
        if sorted(pair) != sorted(end_pair):
            raise ValueError("This contact left the indicated charge sector")
        pos = edge[end_pair.index(2)] if 2 in end_pair else word.excess
        neg = edge[end_pair.index(0)] if 0 in end_pair else word.vacancy
        flips = word.flips ^ {edge} if end_coin != coin else word.flips
        yield Word(pos, neg, flips), int(entries[index])


def free(ray, ref=1):
    edges, degree, near = raw.board(0)
    matrix = seam()[0]
    result = Counter()
    for word, amplitude in ray.items():
        ports = [p for p in (word.excess, word.vacancy) if p is not None]
        for idx in {e for p in ports for e in near[p].values()}:
            for out, value in contact(word, edges[idx], matrix):
                result[out] += amplitude * value
        for p in ports:
            cost = (word.population(p)-ref)**2 - (1-ref)**2
            result[word] += amplitude * 2 * (6-int(degree[p])) * cost
    return {word: val for word, val in result.items() if val}


def edges_in(path):
    return frozenset(tuple(sorted((a, b))) for a, b in zip(path, path[1:]))


def parent_and_tips(paths, sign):
    parent = paths[0][:-1]
    if paths[1][:-1] != parent:
        raise ValueError("Tips do not share a parent")
    return tuple(Word(path[-1] if sign == 1 else None, path[-1] if sign == -1 else None,
                      edges_in(path)) for path in (parent, *paths))


def prepare(words, sign):
    # The two edge flags and center |d|=1 flag are read after native echoes.
    ports = [w.excess if sign == 1 else w.vacancy for w in words]
    if raw.board(0)[1][ports[0]] != 5:
        raise ValueError("The five-query preparation needs a center degree reference")
    hits = []
    for tip in ports[1:]:
        edge = tuple(sorted((ports[0], tip)))
        hit = np.zeros((3, 3), dtype=int)
        for i, word in enumerate(words):
            for out, val in contact(word, edge, seam()[1]):
                hit[words.index(out), i] += val
        hits.append(hit)
    center = 2 * np.diag([1, 0, 0])  # signed pointer columns checked in seam()
    v = np.array([1, 0, 0])
    refusals = []
    for step, (effect, yes) in enumerate(zip([*hits, center, hits[0], center],
                                              [True, False, False, False, False], strict=True), 1):
        hit = effect @ v
        quiet = 2*v - hit
        v, rejected = (hit, quiet) if yes else (quiet, hit)
        refusals.append(Fraction(int(rejected @ rejected), 4**step))
    acceptance = Fraction(int(v @ v), 4**5)
    if not (acceptance == Fraction(1, 8) and sum(refusals) + acceptance == 1
            and v[0] == 0 and v[1] == -v[2]):
        raise ValueError("Signed preparation or refusal law differs")
    return acceptance, {words[1]: 1, words[2]: -1}


def join(a, b):
    if a.excess == b.vacancy or a.flips & b.flips:
        raise ValueError("The selected charge and vacancy routes overlap")
    return Word(a.excess, b.vacancy, a.flips | b.flips)


def eigen(ray, reference, value):
    expected = {word: value * amplitude for word, amplitude in ray.items()}
    if free(ray, ref=reference) != expected:
        raise ArithmeticError("Free flight mixes the prepared signed record")


def verify():
    _, degree, near = raw.board(0)
    # Same chart for both charges. Reverse motion starts at the largest address.
    root = len(degree)-1
    leaves = [p for p in sorted(near[root]) if degree[p] == 5 and all(q > p for q in near[p])]
    if len(leaves) < 2:
        raise ValueError("The vacancy source needs two terminal tips")
    neg_paths = tuple((root, tip) for tip in leaves[:2])
    negative = parent_and_tips(neg_paths, -1)
    p_neg, neg_ray = prepare(negative, -1)
    # At the original port 0 the two signs already give different transport.
    empty = frozenset()
    source_moves = [sum(w != source for w in free({source: 1}))
                    for source in (Word(0, None, empty), Word(None, 0, empty))]
    if source_moves != [5, 0]:
        raise ArithmeticError("Different fixed-chart source control")
    sources = []
    for column, paths in raw.spectroscopy():
        positive = parent_and_tips(paths, 1)
        if {p for path in paths for p in path} & {p for path in neg_paths for p in path}:
            raise ValueError("The signed preparation patches share an addressed port")
        p_pos, pos_ray = prepare(positive, 1)
        pair = {join(a, b): x*y for a, x in pos_ray.items() for b, y in neg_ray.items()}
        if {w.excess for w in pair} != {9, 11} or {w.vacancy for w in pair} != {228, 230}:
            raise ValueError("The terminal pointer addresses have changed")
        if tuple(sorted(x*x for x in pair.values())) != (1, 1, 1, 1):
            raise ValueError("The joint charge readings are not equiprobable")
        # Single-defect and actual joint generators, each relative to its own
        # all-one material background at the specified degree-reference load.
        for ref, pos, neg in ((1, 1, 1), (0, 5, -3), (2, -3, 5)):
            eigen(pos_ray, ref, pos)
            eigen(neg_ray, ref, neg)
            eigen(pair, ref, 2)
        changed = {join(a, b): y for a in pos_ray for b, y in neg_ray.items()}
        escaped = set(free(changed)) - set(changed)
        if not escaped or not any(w.excess == positive[0].excess for w in escaped):
            raise ArithmeticError("Changing the excess tip phase failed to open its return path")
        sources.append({"excess_donor_column": column, "joint_acceptance_given_sharp_parents": str(p_pos*p_neg)})
    return {"signed_material_ports": {"excess": [9, 11], "vacancy": [p[-1] for p in neg_paths]},
            "each_joint_tip_probability_given_acceptance": "1/4",
            "off_diagonal_moves_from_port_0_excess_vacancy": source_moves,
            "twice_vacuum_relative_agreement_eigenvalues": {
                "ordered_as": ["excess", "vacancy", "pair"],
                "reference_1": [1, 1, 2], "reference_0": [5, -3, 2], "reference_2": [-3, 5, 2]},
            "conserved_zero_and_two_counts_on_all_native_free_contact_columns": True, "sources": sources}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, allow_nan=False))
