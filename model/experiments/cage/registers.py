"""The two operators on a profile as words on the 432 registers of one cell.

The register form of the scalar channel (`..compiler`) writes "mode `i` is occupied" as a product `B_i` of `Z`
and "an object passes between `i` and `j`" as a word `A_ij`. For two modes that no link joins, `A_ij` is the
product of the link words along a path, with the phase that a closed loop of links would fix. All four seams of
the ring are joined through one tree: the two seams at the first port of the ring through that port's fuel, the
two seams at the opposite port through that port's fuel, and the two fuels through the diagonal seam between
the two ports. Every word uses paths of this one tree.

    Q_e = 4 - sum_i B_i + sum_{i<j} s_i s_j (i A_ij B_i - i A_ij B_j),        n(e) = Q_e / 8
    K_e = sum_i s_i (i A_if B_i - i A_if B_f),                                c(e)^+ f + f^+ c(e) = K_e / 4

All statements of this module are identities between polynomials of words with integer coefficients. No matrix
on the 432 registers is built.
"""
from __future__ import annotations

import itertools

from model.candidate_laws.scalar_channel.carrier import require
from model.candidate_laws.scalar_channel.compiler import anti, bits, build, local_cycles, mul, weight
from .construction import RING

UNITS = ((1, 0), (0, 1), (-1, 0), (0, -1))          # the powers of i as pairs (real, imaginary)


def add(poly, coefficient, word):
    """Add `coefficient` times a word to a polynomial. A coefficient is a pair of integers (real, imaginary)."""
    phase, x, z = word
    unit = UNITS[phase % 4]
    c = (coefficient[0] * unit[0] - coefficient[1] * unit[1], coefficient[0] * unit[1] + coefficient[1] * unit[0])
    old = poly.get((x, z), (0, 0))
    new = (old[0] + c[0], old[1] + c[1])
    if new == (0, 0):
        poly.pop((x, z), None)
    else:
        poly[(x, z)] = new


def product(left, right):
    out = {}
    for (xa, za), a in left.items():
        for (xb, zb), b in right.items():
            add(out, (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]), mul((0, xa, za), (0, xb, zb)))
    return out


def combination(*terms):
    """An integer combination of polynomials: terms are pairs (integer, polynomial)."""
    out = {}
    for factor, poly in terms:
        for (x, z), c in poly.items():
            add(out, (factor * c[0], factor * c[1]), (0, x, z))
    return out


def hermitian(poly):
    """A polynomial equals its adjoint: the adjoint of c X^x Z^z is conj(c) (-1)^(x.z) X^x Z^z."""
    return all(c == ((c[0], -c[1]) if (x & z).bit_count() % 2 == 0 else (-c[0], c[1])) for (x, z), c in poly.items())


def along(g, path):
    """The word of a passage along a path of links, with the phase that closed loops fix."""
    word = ((len(path) - 2) % 4, 0, 0)
    for a, b in zip(path, path[1:]):
        word = mul(word, g['oriented'](a, b))
    return word


def profile_number(g, seams, signs, fuels, diagonal):
    """Q_e = 8 n(e) for the profile with the given signs on four seam modes.

    `fuels[k]` is the fuel mode linked to seam `k`. Two different fuels are joined through the mode `diagonal`.
    """
    def route(a, b):
        if fuels[a] == fuels[b]:
            return [seams[a], fuels[a], seams[b]]
        return [seams[a], fuels[a], diagonal, fuels[b], seams[b]]
    Q = {}
    add(Q, (4, 0), (0, 0, 0))
    for s in seams:
        add(Q, (-1, 0), g['B'][s])
    for a, b in itertools.combinations(range(4), 2):
        A = along(g, route(a, b))
        sign = signs[a] * signs[b]
        add(Q, (0, sign), mul(A, g['B'][seams[a]]))
        add(Q, (0, -sign), mul(A, g['B'][seams[b]]))
    return Q


def converter(g, seams, signs, fuels, diagonal, fuel):
    """K_e = 4 (c(e)^+ f + f^+ c(e)) for the fuel mode `fuel`, on paths of the same tree."""
    K = {}
    for k, s in enumerate(seams):
        A = along(g, [s, fuel] if fuels[k] == fuel else [s, fuels[k], diagonal, fuel])
        add(K, (0, signs[k]), mul(A, g['B'][s]))
        add(K, (0, -signs[k]), mul(A, g['B'][fuel]))
    return K


def tree(g, ring=RING):
    """The four seam modes of a ring, their signs, their fuels and the diagonal seam, on the registers of one cell."""
    index = g['data']['index']
    first, opposite = ring[0], ring[2]
    require(tuple(sorted((first, opposite))) in index, 'the first and the third port of the ring are joined by a seam')
    seams, fuels = [], []
    for k in range(4):
        pair = tuple(sorted((ring[k], ring[(k + 1) % 4])))
        seams.append(index[pair])
        fuels.append(g['F'][first] if first in pair else g['F'][opposite])
    diagonal = index[tuple(sorted((first, opposite)))]
    links = [(s, f) for s, f in zip(seams, fuels)] + [(diagonal, g['F'][first]), (diagonal, g['F'][opposite])]
    require(all(tuple(sorted(link)) in g['ix'] for link in links), 'a link of the tree is no register')
    return seams, (1, -1, 1, -1), fuels, diagonal


def run():
    """The identities of the two operators on the 432 registers of one cell, as exact statements."""
    g = build()
    loops = [row['pauli'] for row in local_cycles(g)]
    seams, signs, fuels, diagonal = tree(g)
    fuel = g['F'][RING[0]]
    Q = profile_number(g, seams, signs, fuels, diagonal)
    K = converter(g, seams, signs, fuels, diagonal, fuel)
    identity = {(0, 0): (1, 0)}
    R = combination((1, identity), (-1, {g['B'][fuel][1:]: (1, 0)}))                 # R = 2 n_f
    numbers = combination(*[(1, {g['B'][i][1:]: (1, 0)}) for i in range(g['nsource'])])
    words_q = [(0, x, z) for x, z in Q]
    words_k = [(0, x, z) for x, z in K]
    stars = lambda words: sorted({g['owner'][k] for w in words for k in bits(w[1] | w[2])})
    most = lambda words: max(len({g['owner'][k] for k in bits(w[1] | w[2])}) for w in words)
    others = [k for k in g['F'] + g['R']]
    return dict(
        registers=g['M'], modes=g['n'], loops=len(loops),
        profile_number=dict(
            words=len(Q), largest_weight=max(weight(w) for w in words_q), hermitian=hermitian(Q),
            terms_left_in_Q_squared_minus_8_Q=len(combination((1, product(Q, Q)), (-8, Q))),
            pairs_with_the_loops=len(words_q) * len(loops),
            pairs_with_the_loops_that_do_not_commute=sum(anti(w, loop) for w in words_q for loop in loops),
            words_that_do_not_commute_with_a_fuel_or_reference_number=sum(anti(w, g['B'][k]) for w in words_q for k in others),
            terms_left_in_the_commutator_with_the_number_of_objects=len(combination((1, product(Q, numbers)), (-1, product(numbers, Q)))),
            stars_of_all_words=stars(words_q), most_stars_in_one_word=most(words_q)),
        converter=dict(
            words=len(K), largest_weight=max(weight(w) for w in words_k), hermitian=hermitian(K),
            # G_e^2 = n(e) + n_f - 2 n(e) n_f, written for K = 4 G_e, Q = 8 n(e), R = 2 n_f
            terms_left_in_K_squared_minus_2Q_minus_8R_plus_2QR=len(combination((1, product(K, K)), (-2, Q), (-8, R), (2, product(Q, R)))),
            # [G_e, n(e) + n_f] = 0
            terms_left_in_the_commutator_with_Q_plus_4R=len(combination((1, product(K, combination((1, Q), (4, R)))),
                                                                       (-1, product(combination((1, Q), (4, R)), K)))),
            pairs_with_the_loops=len(words_k) * len(loops),
            pairs_with_the_loops_that_do_not_commute=sum(anti(w, loop) for w in words_k for loop in loops),
            words_that_do_not_commute_with_another_fuel_or_a_reference_number=sum(anti(w, g['B'][k]) for w in words_k for k in others if k != fuel),
            stars_of_all_words=stars(words_k), most_stars_in_one_word=most(words_k)),
        guard='one product of the code checks of the ports %s stands in front of each whole operator' % (list(RING),))
