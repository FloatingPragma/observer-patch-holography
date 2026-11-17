"""Exact CAR and link-polynomial utilities for the CCF-2 construction.

CAR words use (mode, creation). Link exponents describe full rotors. Finite
hard-wall operators must be multiplied as matrices, never by cancelling
opposite formal exponents. The exact CAR normal-ordering implementation is
part of this package and checked against complete finite matrices.
"""
from __future__ import annotations
from collections import defaultdict
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import sympy as sy

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
from . import gauging_compiler as compiler
normal = compiler.normal


def clean(poly):
    return {w: value for w, raw in poly.items() if (value := sy.expand(raw)) != 0}


def add(*polys):
    out = defaultdict(lambda: sy.Integer(0))
    for p in polys:
        for w, c in p.items():
            out[w] += c
    return clean(out)


def scale(poly, coefficient):
    coefficient = sy.sympify(coefficient)
    return clean({w: coefficient*c for w, c in poly.items()})


def mono(word=(), coefficient=1):
    return {w: sy.sympify(coefficient)*s for w, s in normal(tuple(word))}


def letter(i, creation=False):
    return mono(((int(i), bool(creation)),))


def multiply(left, right):
    return compiler.multiply(left, right)


mul = multiply
adjoint = compiler.dagger
dagger = adjoint


def number(i):
    return mono(((int(i), True), (int(i), False)))


def index(v, s, c, f):
    return 96*v + 48*s + 3*c + f


def substitute(poly, mapping):
    expanded_mapping = {}
    for key, value in mapping.items():
        if isinstance(key, int):
            expanded_mapping[(key, False)] = value
            expanded_mapping[(key, True)] = adjoint(value)
        else:
            expanded_mapping[key] = value
    out = {}
    for word, coefficient in poly.items():
        image = mono(coefficient=coefficient)
        for item in word:
            image = multiply(image, expanded_mapping.get(item, letter(*item)))
        out = add(out, image)
    return out


def reindex(poly, mapping):
    out = defaultdict(lambda: sy.Integer(0))
    for word, coefficient in poly.items():
        for w, sign in normal(tuple((mapping.get(i, i), c) for i, c in word)):
            out[w] += sign*coefficient
    return clean(out)


def vacuum_restrict(poly, allowed_modes):
    """Vacuum matrix element on other modes, NOT an invariant-sector claim."""
    allowed = set(allowed_modes)
    return {w: c for w, c in poly.items() if all(i in allowed for i, _ in w)}


def lift(poly, nedges):
    return dress(poly, (0,)*nedges)


def dress(poly, shift):
    return {(w, tuple(shift)): c for w, c in poly.items()}


gadd = add


def gadjoint(poly):
    out = defaultdict(lambda: sy.Integer(0))
    for (word, shift), coefficient in poly.items():
        for w, sign in normal(tuple((i, not c) for i, c in reversed(word))):
            out[(w, tuple(-x for x in shift))] += sign*sy.conjugate(coefficient)
    return clean(out)


def forget_links(poly):
    out = defaultdict(lambda: sy.Integer(0))
    for (word, _), coefficient in poly.items():
        out[word] += coefficient
    return clean(out)


@lru_cache(None)
def native_charge_coordinates():
    """The literal native-coordinate charge ordering of the retained source."""
    from .snapshot import snapshot
    source = snapshot()
    charges = tuple(source['mode_charge_native_thirds'][3*c] for c in range(16))
    assert len(charges) == 16
    return charges


def mode_metadata(nvertices):
    q = native_charge_coordinates()
    return ([q[c] for v in range(nvertices) for s in range(2)
             for c in range(16) for f in range(3)],
            [v for v in range(nvertices) for _ in range(96)])


def gauge(poly, nvertices, edges):
    charges, vertices = mode_metadata(nvertices)
    rows = compiler.compile_poly(poly, charges, vertices, nvertices, edges)
    return {(tuple((int(i), bool(c)) for i, c in row['word']),
             tuple(row['link_shift'])): sy.sympify(row['coefficient_exact_as_supplied'])
            for row in rows}


def check_gauged(poly, nvertices, edges, charges=None):
    if charges is None:
        charges, _ = mode_metadata(nvertices)
    assert clean(poly) == gadjoint(poly), 'non-Hermitian gauged polynomial'
    for (word, shift), coefficient in poly.items():
        assert coefficient != 0 and len(shift) == len(edges)
        s = [0]*nvertices
        for i, creation in word:
            assert 0 <= i < 96*nvertices
            s[i//96] += charges[i]*(1 if creation else -1)
        divergence = [0]*nvertices
        for value, (u, v) in zip(shift, edges):
            assert isinstance(value, (int, sy.Integer))
            divergence[u] += value
            divergence[v] -= value
        assert s == divergence, (word, s, divergence)
    return {'terms': len(poly), 'hermitian': True, 'all_local_gauss_identities': True}


def polynomial_digest(poly):
    payload = [(repr(k), str(v)) for k, v in sorted(poly.items())]
    return hashlib.sha256(json.dumps(payload, separators=(',', ':')).encode()).hexdigest()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def apply_word(word, mask):
    """Literal complete Fock action; return (new_mask, sign) or None."""
    sign = 1
    for i, creation in reversed(word):
        occupied = bool(mask & (1 << i))
        if occupied == creation:
            return None
        sign *= (-1)**((mask & ((1 << i)-1)).bit_count())
        mask ^= 1 << i
    return mask, sign
