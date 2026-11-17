"""Loads on the carrier with one address register, moved by the public local completion of the seam event.

A state is a dictionary {(configuration, record, address index): amplitude}. The configuration is the tuple of
occupations of the twelve ports, the record takes the four values of the completion and the address index runs over
a declared list of ordered seams. One tick is the repair after the address reflection, as in `Observer.step`; the
inverse tick is the reflection after the repair, as in `Observer.inverse_step`. Amplitudes are exact fractions, or
any number type with + and *.
"""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from model.coherent_seam_repair.common import local_repair_unitary

SOURCES = ('00', '10', '01', '11', '20', '02')   # occupation order of the public completion


@lru_cache(maxsize=1)
def completion():
    """The public 24 x 24 completion (source times four record values) as exact fractions."""
    matrix = local_repair_unitary()
    out = []
    for i in range(24):
        row = []
        for j in range(24):
            value = s.Rational(matrix[i, j])
            row.append(F(int(value.p), int(value.q)))
        out.append(tuple(row))
    return tuple(out)


@lru_cache(maxsize=None)
def images(cfg, rec, a, b):
    """Images of (configuration, record) under the completion on the ordered seam (a, b)."""
    na, nb = cfg[a], cfg[b]
    if na + nb > 2:
        raise ValueError('a seam with more than two loads lies outside the local completion')
    col = 4 * SOURCES.index(f'{na}{nb}') + rec
    out = []
    for row, line in enumerate(completion()):
        u = line[col]
        if u:
            source, r2 = divmod(row, 4)
            c2 = list(cfg)
            c2[a], c2[b] = int(SOURCES[source][0]), int(SOURCES[source][1])
            out.append((tuple(c2), r2, u))
    return tuple(out)


def add(d, key, value):
    v = d.get(key, 0) + value
    if v:
        d[key] = v
    else:
        d.pop(key, None)


def config(*ports):
    c = [0] * 12
    for p in ports:
        c[p] += 1
    return tuple(c)


def ports_of(cfg):
    return tuple(p for p in range(12) for _ in range(cfg[p]))


class Walk:
    """One address register over `addresses`; at address k the completion acts on the ordered seam `pairs[k]`.

    With pairs equal to the addresses this is the plain law of `Observer.repair`. A frame that uses one fixed
    ordered seam at both addresses of a seam is the oriented law.
    """

    def __init__(self, addresses, pairs=None):
        self.addresses = tuple(tuple(a) for a in addresses)
        self.pairs = tuple(tuple(p) for p in pairs) if pairs is not None else self.addresses
        self.size = len(self.addresses)

    def repair(self, state):
        out = {}
        for (cfg, rec, k), amp in state.items():
            a, b = self.pairs[k]
            for c2, r2, u in images(cfg, rec, a, b):
                add(out, (c2, r2, k), amp * u)
        return out

    def reflect(self, state):
        """2 A A^dag - 1 about the uniform address state, the form of `Observer.reflect`."""
        blocks = {}
        for (cfg, rec, k), amp in state.items():
            blocks.setdefault((cfg, rec), {})[k] = amp
        out = {}
        for (cfg, rec), block in blocks.items():
            total = 0
            for v in block.values():
                total = total + v
            mean = total * F(1, self.size)
            for k in range(self.size):
                v = 2 * mean - block.get(k, 0)
                if v:
                    out[(cfg, rec, k)] = v
        return out

    def step(self, state):
        return self.repair(self.reflect(state))

    def inverse_step(self, state):
        return self.reflect(self.repair(state))

    def uniform(self, amplitudes, rec=0):
        """A|v>: {configuration: amplitude} with a record value and amplitude one on every address (unnormalized)."""
        return {(cfg, rec, k): a for cfg, a in amplitudes.items() for k in range(self.size)}


def norm2(state):
    total = 0
    for v in state.values():
        total = total + (abs(v) ** 2 if isinstance(v, complex) else v * v)
    return total


def dot(u, v):
    if len(u) > len(v):
        u, v = v, u
    total = 0
    for key, c in u.items():
        if key in v:
            total = total + (c.conjugate() if isinstance(c, complex) else c) * v[key]
    return total


def restrict(state, predicate):
    return {key: v for key, v in state.items() if predicate(key)}


def combine(u, v, b=1):
    out = dict(u)
    for key, c in v.items():
        add(out, key, b * c)
    return out


def weights(state, n2, classify):
    out = {}
    for key, amp in state.items():
        label = classify(key)
        out[label] = out.get(label, 0) + amp * amp / n2
    return out


# ----------------------------------------------------------------------------------------- independent pairs of objects
def tensor(u, v):
    return {(ka, kb): ca * cb for ka, ca in u.items() for kb, cb in v.items()}


def on_first(pair, op):
    groups = {}
    for (ka, kb), c in pair.items():
        groups.setdefault(kb, {})[ka] = c
    out = {}
    for kb, state in groups.items():
        for ka, c in op(state).items():
            add(out, (ka, kb), c)
    return out


def on_second(pair, op):
    groups = {}
    for (ka, kb), c in pair.items():
        groups.setdefault(ka, {})[kb] = c
    out = {}
    for ka, state in groups.items():
        for kb, c in op(state).items():
            add(out, (ka, kb), c)
    return out


def on_both(pair, op_first, op_second=None):
    return on_second(on_first(pair, op_first), op_second or op_first)
