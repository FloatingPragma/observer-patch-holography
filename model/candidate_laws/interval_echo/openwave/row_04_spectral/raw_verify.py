"""Independent calculations for the spectral release. Inputs: carrier NPZs and donor JSON."""

from collections import Counter
from dataclasses import dataclass
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from functools import cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm


HERE = Path(__file__).resolve().parent
if not __debug__:
    raise RuntimeError("Run the independent oracle without -O so its assertion gates are active")
DONOR = HERE.parent / "row_30_selective"
NPZ = HERE.parents[3] / "executable_model" / "carriers"
BOARD_HASHES = (
    "99870398e3d516f899c600e9bf25ba53c7f513142ba8facbe1cedd94b23f915b",
    "c7d4ad87b09d9158a3d9be60c537d1aaecc69e951e6a9979e505e0162e134d4b",
    "e7273bbd0e814abf49f3a333ecfe58a542262aee7d2399db677d15c19fcac673",
)
CAPS = (1, 8, 16, 32, 96)
STRENGTHS = (0, .5, -.5, 8, -8)
TIMES = (0, 1, 7, 37, 150)
ARMS = ("minus", "bare_fuel_minus", "plus", "dephase_at_release", "random_swaps")
LOADS = tuple(product(range(3), repeat=2))
LABELS = (None, *product(range(2), range(-2, 3)))


def read_json(path):
    def unique(pairs):
        assert len(dict(pairs)) == len(pairs), "Repeated JSON key"
        return dict(pairs)

    return json.loads(path.read_bytes(), object_pairs_hook=unique, parse_constant=invalid)


def invalid(value):
    raise ValueError(f"Invalid numeric literal: {value}")


def ix(label, loads):
    return 9 * LABELS.index(label) + LOADS.index(loads)


@cache
def local():
    """Integer 2U in coordinates with physical ket = S x, S_written = 1/sqrt(2)."""
    unit = np.eye(99, dtype=int)
    writes = np.zeros((99, 9), dtype=int)
    reversal = unit.copy()
    reflect = []
    d2 = np.zeros(99, dtype=int)
    for j, (a, b) in enumerate(LOADS):
        lo, odd = divmod(a + b, 2)
        for coin in (0, 1):
            pair = (lo + odd * coin, lo + odd * (1 - coin))
            writes[ix((coin, a-b), pair), j] = 1
    for label, pair in product(LABELS, LOADS):
        if label is None:
            mirror = None
        else:
            c, d = label
            mirror = (1-c, -d)
            j = ix(label, pair)
            reversal[:, j] = unit[:, ix((c, -d), pair)]
            d2[j] = d*d
        reflect.append(ix(mirror, pair[::-1]))
    echo2 = 2*unit - writes @ writes.T
    echo2[:9, :9] = 0
    echo2[:, :9] = 2*writes
    echo2[:9, :] = writes.T
    np.testing.assert_array_equal(echo2 @ echo2, 4*unit)
    np.testing.assert_array_equal(echo2[np.ix_(reflect, reflect)], echo2)
    # 8 H_exchange and 4 H_agreement in the same rational coordinates.
    ex8 = echo2 @ (unit-reversal) @ echo2
    ag4 = echo2 @ (d2[:, None]*echo2)
    for h in (ex8, ag4):
        np.testing.assert_array_equal(h[np.ix_(reflect, reflect)], h)
    arrows = np.stack([unit[:, ix((c, -1), w)] + unit[:, ix((c, 1), w)]
                       for c, w in product((0, 1), ((0, 0), (1, 0), (0, 1)))], axis=1)
    # Written sign-even columns have unit norm in this metric.
    np.testing.assert_array_equal(ex8 @ arrows, 0*arrows)
    compressed2 = (arrows.T @ ag4 @ arrows) // 4
    np.testing.assert_array_equal(ag4 @ arrows, 2*arrows @ compressed2)
    return echo2, ex8, ag4, arrows, compressed2, reflect


def projector4(mask):
    echo2 = local()[0]
    return echo2 @ (np.asarray(mask)[:, None]*echo2)


@cache
def queries():
    """Contract each physical hit/quiet with its three distinct joint input histories."""
    arrows = local()[3]
    hit4 = projector4(np.arange(99) < 9)
    quiet4 = 4*np.eye(99, dtype=int) - hit4
    ancilla4 = projector4(np.repeat([t is not None and abs(t[1]) == 1 for t in LABELS], 9))
    output = []
    for tip in (1, 2):
        projected = np.zeros((3, 3), dtype=int)
        projected[0, 0] = projected[tip, tip] = projected[0, tip] = projected[tip, 0] = 2
        for c in (0, 1):
            # Head of coin 1 is the low endpoint, head of coin 0 the high one.
            cols = [3*c, 3*c, 3*c]
            cols[0], cols[tip] = 3*c+2-c, 3*(1-c)+1+c
            embedding = arrows[:, cols]
            for native, small in ((hit4, projected), (quiet4, 4*np.eye(3, dtype=int)-projected)):
                np.testing.assert_array_equal(native @ embedding, embedding @ small)
                # Coin resolved after the flag: four flag/coin branches on a sharp head.
                for old_coin in (0, 1):
                    mask = np.repeat([l is not None and l[0] == old_coin for l in LABELS], 9)
                    post = mask*(native @ embedding[:, 0])
                    index = 0 if c == old_coin else tip
                    assert np.array_equal(post, 2*embedding[:, index]) or np.array_equal(post, -2*embedding[:, index])
        output.append(projected // 2)
    # The center uses a blank degree reference. The two absent inputs
    # have distinct spectator records and the same vacant local column.
    E = np.eye(99, dtype=int)
    for pair in ((1, 0), (0, 1)):
        embed = E[:, [ix(None, pair), 0, 0]]
        P = np.diag([1, 0, 0])
        np.testing.assert_array_equal(ancilla4 @ embed, 4*embed @ P)
        np.testing.assert_array_equal((4*E-ancilla4) @ embed, 4*embed @ (np.eye(3, dtype=int)-P))
        np.testing.assert_array_equal(local()[2] @ embed, 4*embed @ P)
    return (*output, 2*P)  # Twice the three physical hit maps (A, B, C).


@cache
def attempt(sign, omit_last_center=False):
    q = queries()
    v = np.array([F(1), F(0), F(0)], dtype=object)
    rows = []
    for label, yes in list(zip((0, 1, 2, 0, 2), (sign == -1, False, False, False, False)))[:4 if omit_last_center else 5]:
        hit = (q[label] @ v) / 2
        kept, refused = (hit, v-hit) if yes else (v-hit, hit)
        # Reading both child coins resolves (unchanged, A flipped, B flipped).
        rows.append(tuple(refused**2))
        assert sum(kept**2) + sum(rows[-1]) == sum(v**2)
        v = kept
    return tuple(v), tuple(rows)


def joint_attempt(sign, mounts):
    """Apply the native selective maps to both full collar memories and the occupied port."""
    blank_hit4 = projector4(np.arange(99) < 9)
    query4 = []
    basis = np.zeros((3, 11, 11, 3), dtype=int)
    for word in range(3):
        for s, t in product((-1, 1), repeat=2):
            a = LABELS.index((mounts[0] ^ (word == 1), s))
            b = LABELS.index((mounts[1] ^ (word == 2), t))
            basis[word, a, b, word] = 1
    for tip, c in enumerate(mounts, start=1):
        occupied = (c, 1-c)
        pairs = [(0, 0)]*3
        pairs[0], pairs[tip] = occupied, occupied[::-1]
        ids = [ix(label, pair) for label, pair in product(LABELS, pairs)]
        np.testing.assert_array_equal(np.delete(blank_hit4[:, ids], ids, axis=0), 0)
        query4.append(blank_hit4[np.ix_(ids, ids)])
    ref4 = projector4(np.repeat([lab is not None and abs(lab[1]) == 1 for lab in LABELS], 9))
    center4 = np.array([ref4[ix(None, w), ix(None, w)] for w in ((1, 0), (0, 0), (0, 0))])
    v = basis[0]
    refusals = []
    for step, (q, yes) in enumerate(zip((0, 1, 2, 0, 2), (sign == -1, False, False, False, False)), 1):
        if q == 2:
            hit = v*center4
        else:
            perm = (q, 2, 1-q)
            block = v.transpose(perm).reshape(33, 11)
            hit = (query4[q] @ block).reshape(11, 3, 11).transpose(np.argsort(perm))
        quiet = 4*v-hit
        v, rejected = (hit, quiet) if yes else (quiet, hit)
        # After a refusal read both coins in place. Each ordered pair of coin
        # outcomes must leave exactly the full corresponding sharp word.
        row = [F(0)]*3
        for flipped in product((0, 1), repeat=2):
            keep = [np.array([lab is not None and lab[0] == c ^ f for lab in LABELS])
                    for c, f in zip(mounts, flipped)]
            post = rejected*keep[0][:, None, None]*keep[1][None, :, None]
            if flipped == (1, 1):
                assert not post.any()
                continue
            word = flipped[0]+2*flipped[1]
            amplitude = post[basis[word] != 0][0]
            np.testing.assert_array_equal(post, amplitude*basis[word])
            row[word] = F(int(amplitude), 4**step)**2
        refusals.append(tuple(row))
    coeff = [int(v[b != 0][0]) for b in basis]
    np.testing.assert_array_equal(v, np.einsum("i,ijkl->jkl", coeff, basis))
    return tuple(F(x, 4**5) for x in coeff), tuple(refusals)


@cache
def fuel(k, matched, coin):
    """Both measured 99-component post-pulse kets, from the full two-rail exponential."""
    _, ex8, ag4, arrows, _, mirror = local()
    scale = np.r_[np.ones(9), np.full(90, 1/np.sqrt(2))]
    h = (ex8/8 + k*ag4/4)*scale[:, None]/scale[None, :]
    T = np.zeros((99, 99))
    for pair in LOADS:
        for d in (-1, 1):
            T[ix((0, d), pair), ix(None, pair)] = 1/np.sqrt(2)
    if coin:
        T = T[np.ix_(mirror, mirror)]
    h2 = np.block([[h, 64*T], [64*T.T, h + (k if matched else 0)*np.eye(99)]])
    state = expm((-1j*np.pi/256)*h2)[:, 99].reshape(2, 99)
    target = scale * arrows[:, 3*coin]
    vac = np.eye(99)[0]
    for ket, ray in zip(state, (target, vac)):
        np.testing.assert_allclose(ket, (ray @ ket)*ray, atol=2e-14, rtol=0)
    mass = np.sum(abs(state)**2, axis=1)
    assert abs(sum(mass)-1) < 2e-14
    return float(mass[0])


@cache
def ledger(edges, depth, n, sign, p):
    """Forward absorbing instrument; Decimal exponents retain subnormal original-start yields."""
    with localcontext() as ctx:
        ctx.prec = 85
        v, rows = attempt(sign)
        as_d = lambda f: D(f.numerator)/D(f.denominator)
        wf, hf = (1-D(str(p)))**n, D(1)/D(2)**n
        # Each of the first E-1 empty writers has the same fresh fuel/read/reset.
        through = (1-wf)**(edges-1)
        f = {"writer_before_return": 1-through, "commissioning_hop": through*hf}
        through *= 1-hf
        f["last_write"] = through*wf
        through *= 1-wf
        f["route_to_parent"] = D(0)
        for _ in range(depth):
            f["route_to_parent"] += through*hf
            through *= 1-hf
        accept, f["failed_reset"], f["no_accepted_word"] = D(0), D(0), D(0)
        accepted = as_d(sum(x*x for x in v))
        at = [as_d(sum(t)) for t in zip(*rows)]
        for t in range(n):
            accept += through*accepted
            at_center, left, right = (through*x for x in at)
            if t+1 == n:
                f["no_accepted_word"] = at_center+left+right
            else:
                f["failed_reset"] += (left+right)*hf
                through = at_center + (left+right)*(1-hf)
        assert abs(accept + sum(f.values()) - 1) < D("1e-79")
        return accept, f


@cache
def board(level):
    path = NPZ / f"FEDERATION_LEVEL{level}.npz"
    assert sha256(path.read_bytes()).hexdigest() == BOARD_HASHES[level]
    with np.load(path, allow_pickle=False) as z:
        edges = tuple(tuple(sorted(map(int, row))) for row in z["edges"])
        degree = z["degree"]
    assert len(set(edges)) == len(edges) and all(0 <= a < b < len(degree) for a, b in edges)
    incidence = [dict() for _ in degree]
    for i, (a, b) in enumerate(edges):
        incidence[a][b] = incidence[b][a] = i
    np.testing.assert_array_equal([len(row) for row in incidence], degree)
    return edges, degree, incidence


def parking(level, m, f, wrong_mount=False):
    """Check the specified installation by writing and flipping physical coin bits."""
    edges, _, inc = board(level)
    assert m in inc[0] and f in inc[0] and f in inc[m] and max(m, f) < 12
    coin = np.full(len(edges), -1)
    here = hops = 0
    mf, rm, rf = inc[m][f], inc[0][m], inc[0][f]

    def write(i, c=1):
        assert here not in edges[i] and coin[i] == -1
        coin[i] = c

    def transfer(to, blank):
        nonlocal here, hops
        i = inc[here][to]
        if blank:
            assert coin[i] == -1
        else:
            assert coin[i] == int(here == edges[i][0]), "Return arrow points away from the load"
            coin[i] ^= 1
        here, hops = to, hops+1

    for i, e in enumerate(edges):
        if 0 not in e and i != mf:
            write(i)
    transfer(m, True)
    for i in inc[0].values():
        if i != rm:
            write(i, int(i != rf or wrong_mount))
    transfer(f, True)
    write(rm)
    transfer(0, False)
    write(mf)
    assert here == 0 and np.all(coin == 1)
    return len(coin)+hops


@dataclass
class Component:
    states: list
    words: list
    sites: np.ndarray
    h2: np.ndarray | None
    lap: np.ndarray | None


@cache
def component(level, order="ascending"):
    """Flood from physical outgoing columns on (occupied port, flipped physical seam set)."""
    edges, degree, inc = board(level)
    key = {"ascending": lambda p: p, "descending": lambda p: -p,
           "reverse_local": lambda p: (p//12, -p % 12)}[order]
    init = [int(key(a) < key(b)) for a, b in edges]
    reduced = local()[4] - 2*np.eye(6, dtype=int)  # twice the vacuum-subtracted material term
    assert not reduced[:, [0, 3]].any()
    states = [(0, frozenset())]
    index = {states[0]: 0}
    sites, words, columns = [], [], []
    for port, flips in states:
        occupancy = Counter({0: 1})
        vertices = {0}
        for e in flips:
            a, b = edges[e]
            tail, head = (a, b) if init[e] else (b, a)
            occupancy[tail] -= 1
            occupancy[head] += 1
            vertices.update((tail, head))
        assert {p: w for p, w in occupancy.items() if w} == {port: 1}
        sites.append(next(p for p in occupancy if occupancy[p] == 1))
        word = tuple(sorted(vertices, key=key))
        assert word[-1] == port and len(word) == len(flips)+1
        assert frozenset(inc[a][b] for a, b in zip(word, word[1:])) == flips
        words.append(word)
        col = Counter({index[port, flips]: 2*(6-int(degree[port]))})
        for other in sorted(inc[port]):
            e = inc[port][other]
            c = init[e] ^ (e in flips)
            slot = 1 if port < other else 2
            column = reduced[:, 3*c+slot]
            for j in np.flatnonzero(column):
                new_c, new_slot = divmod(int(j), 3)
                assert new_slot != 0
                dest = edges[e][new_slot-1]
                state = dest, (flips ^ {e} if c != new_c else flips)
                if state not in index:
                    index[state] = len(states)
                    states.append(state)
                col[index[state]] += int(column[j])
        columns.append(col)
    assert len(set(words)) == len(states) == len(set(s[1] for s in states))
    if order == "ascending":
        h2 = np.zeros((len(states), len(states)), dtype=int)
        for i, col in enumerate(columns):
            for j, val in col.items():
                h2[j, i] = val
        np.testing.assert_array_equal(h2.T, h2)
        # Geometric edges come from physically distinct unit transitions.
        adj = (h2 == -1).astype(int)
        np.fill_diagonal(adj, 0)
        assert adj.sum() == 2*(len(states)-1)
        lap = np.diag(adj.sum(axis=0))-adj
    else:
        h2 = lap = None
    return Component(states, words, np.array(sites), h2, lap)


def bands(eigenvalues):
    return [np.flatnonzero(abs(eigenvalues-x) < 1e-8) for x in
            eigenvalues[np.r_[True, np.diff(eigenvalues) > 1e-8]]]


@cache
def spectroscopy():
    """Return every enrolled plain spatial ray; check both charts against the donor."""
    tree = component(0)
    energy, v = np.linalg.eigh(4*tree.h2)
    k2, w = np.linalg.eigh(tree.lap)
    primary = read_json(DONOR/"results.json")["spectral_audit"]["path"]
    alternate = read_json(DONOR/"staggered.json")["eigenspaces"]
    parity = np.array([(-1)**len(flips) for _, flips in tree.states])
    groups = bands(k2)
    chosen = []
    for g, record, alt in zip(bands(energy), primary, alternate, strict=True):
        wave = v[:, g]
        support = np.linalg.norm(wave, axis=1) > 1e-7
        assert record["possible_ports"] == sorted(set(tree.sites[support]))
        assert record["multiplicity"] == len(g) and alt["indices"] == g.tolist()
        for signs, archive, score in ((1, record, "maximum_weight"),
                                      (parity[:, None], alt, "best_weight")):
            weights = [np.linalg.norm((w[:, b]*signs).T @ wave, ord=2)**2 for b in groups]
            best = groups[np.argmax(weights)][0]
            np.testing.assert_allclose([archive["energy"], archive["best_k2"], archive[score]],
                                       [energy[g[0]], k2[best], max(weights)], rtol=0, atol=5.1e-10)
            if score == "best_weight":
                np.testing.assert_allclose(weights[0], alt["rest_weight"], rtol=0, atol=5.1e-10)
                continue
            hits = [b for b, x in zip(groups, weights) if x >= .9]
            np.testing.assert_allclose([k2[b[0]] for b in hits], record["compatible_k2"], atol=5.1e-10, rtol=0)
            for b in hits:
                # The donor's chart fixes the basis by projecting successive coordinates.
                residual = w[:, b] @ w[:, b].T
                rays = []
                for i in range(len(k2)):
                    col = residual[:, i].copy()
                    if np.linalg.norm(col) < 1e-7:
                        continue
                    ray = col/np.linalg.norm(col)
                    residual -= np.outer(ray, ray)
                    rays.append(ray)
                for j, ray in zip(b, rays, strict=True):
                    tip = np.flatnonzero(abs(ray) > 1e-8)
                    np.testing.assert_allclose(ray[tip], [1/np.sqrt(2), -1/np.sqrt(2)], rtol=0, atol=1e-12)
                    chosen.append((int(j), tuple(tree.words[t] for t in tip)))
    return chosen


@cache
def historical_pulse(column, k):
    tree = component(0)
    program = read_json(DONOR/"analyzers.json")
    pulse = program["pulses"][column]
    ports = program["ports"]
    assert ports == sorted(p for p in set(tree.sites) if board(0)[1][p] == 5)
    angles = np.array(pulse["angles"])
    assert angles.shape == (program["layers"], len(ports))
    v = np.eye(len(tree.states), dtype=complex)[0]
    single = expm((-1j*k*program["step"]/2)*tree.h2)
    for layer in angles:
        v = single @ v
        for port, phase in zip(ports, layer):
            v[tree.sites == port] *= np.exp(-1j*k*phase/8)
    paths = dict(spectroscopy())[column]
    a, b = map(tree.words.index, paths)
    overlap = abs(v[a]-v[b])**2/2
    if k == 8:
        assert abs(1-overlap-pulse["infidelity"]) < 2e-10
    out = expm((-1j*k*37/2)*tree.h2) @ v
    return overlap, np.linalg.norm(out[np.isin(tree.sites, tree.sites[[a, b]])])**2


@cache
def propagation(level, k, t):
    tree = component(level)
    return expm((-1j*k*t/2)*tree.h2), expm((-abs(k)*t/2)*tree.lap)


def flight(level, paths, k, t):
    tree = component(level)
    a, b = map(tree.words.index, paths)
    two = np.eye(len(tree.words), dtype=int)[:, [a, b]]
    minus, plus = two @ np.array([1, -1]), two.sum(axis=1)
    np.testing.assert_array_equal(tree.h2 @ minus, minus)
    np.testing.assert_array_equal(tree.lap @ minus, minus)
    assert (tree.h2 @ plus)[tree.words.index(paths[0][:-1])] == -2
    U, M = propagation(level, k, t)
    u, v = (U @ two).T
    np.testing.assert_allclose(u-v, np.exp(-1j*k*t/2)*minus, rtol=0, atol=8e-12)
    out = (abs(u-v)**2/2, abs(u-v)**2/2, abs(u+v)**2/2,
           (abs(u)**2+abs(v)**2)/2, M @ (plus/2))
    for p in out:
        assert min(p) >= -2e-14 and abs(p.sum()-1) < 6e-12
    return out


if __name__ == "__main__":
    # JSON lines preserve tiny original-start probabilities as decimal strings.
    # Both signed words have the same resolved refusal law.
    assert attempt(-1)[1] == attempt(1)[1]
    for level, (column, paths), n in product(range(3), spectroscopy(), CAPS):
        edges = len(board(level)[0])
        for matched, strengths in ((True, (None,)), (False, (*STRENGTHS, -128, 128))):
            for k in strengths:
                p = .5 if matched else fuel(k, False, 0)
                acc, reject = ledger(edges, len(paths[0])-2, n, -1, p)
                print(json.dumps({"cells": 20*4**level, "column": column, "cap": n,
                                  "matched": matched, "strength": k, "accept": str(acc),
                                  "refused": {name: str(q) for name, q in reject.items()}}))
