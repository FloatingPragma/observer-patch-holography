"""Exact full-Gauss Hamiltonians and their complete flux records."""
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
from math import isqrt

from model.interacting_extension.higgs_join.interacting_carrier import (
    HERE as SOURCE, ROOT, NativeCarrier, reduce_exact, diagonal_operators)


def fraction(value):
    if type(value) not in (int, F):
        raise TypeError('exact integer or Fraction control required')
    return F(value)


@dataclass(frozen=True)
class Controls:
    electric: F = F(1)
    magnetic: F = F(1)
    hopping: F = F(1)
    higgs: F = F(1)
    scalar: F = F(1)

    def __post_init__(self):
        for name in self.__dataclass_fields__:
            object.__setattr__(self, name, fraction(getattr(self, name)))

    def packet(self):
        return {k: str(getattr(self, k)) for k in self.__dataclass_fields__}


def add(out, key, value):
    value = out.get(key, F(0)) + value
    if value:
        out[key] = value
    else:
        out.pop(key, None)


def sqrt_upper(value, bits=16):
    """Outward dyadic square root, using integer comparisons only."""
    value = F(value)
    if value < 0:
        raise ValueError('nonnegative square required')
    scaled = value * (1 << (2 * bits))
    n = isqrt(scaled.numerator // scaled.denominator)
    if n*n*scaled.denominator < scaled.numerator:
        n += 1
    return F(n, 1 << bits)


def inputs():
    receipt = json.loads((SOURCE/'joined_operator_receipt.json').read_text())
    for name, key in (
        ('joined_basis.json', 'basis_sha256'),
        ('joined_exact_operators.json', 'exact_operator_sha256'),
    ):
        if hashlib.sha256((SOURCE/name).read_bytes()).hexdigest() != receipt[key]:
            raise ValueError('source exact packet checksum: '+name)
    if receipt['families'] != 3 or receipt['complete_Gauss_dimensions'] != {'0':438, '1':1596}:
        raise ValueError('complete three-family source required')
    return (json.loads((SOURCE/'joined_basis.json').read_text()),
            json.loads((SOURCE/'joined_exact_operators.json').read_text())['operators'])


class Sector:
    """Coordinates use the source's orthogonal, unnormalized exact basis."""
    def __init__(self, q, controls=Controls(), data=None):
        if type(q) is not int or q not in (0, 1):
            raise ValueError('source Qstar sector 0 or 1 required')
        self.q, self.controls = q, controls
        basis, operators = inputs() if data is None else data
        self.basis = basis[str(q)]
        self.labels = self.basis['labels']
        self.norms = tuple(F(x) for x in self.basis['raw_norms'])
        self.dimension = self.basis['dimension']
        if len(self.norms) != self.dimension or any(x <= 0 for x in self.norms):
            raise ValueError('positive complete Gram diagonal required')
        prefix = f'q{q}_'
        diagonal = operators[prefix+'diagonal']
        numerators = {}
        for i, (e, nh) in enumerate(zip(diagonal['electric_native_trace'],
                                       diagonal['N_H'], strict=True)):
            add(numerators, (i, i), self.norms[i]*(controls.electric*F(e)+controls.scalar*F(nh)))
        def hermitian(name, coefficient):
            for i, j, value in operators[prefix+name]:
                add(numerators, (i, j), coefficient*F(value))
                add(numerators, (j, i), coefficient*F(value))
        hermitian('loop', -controls.magnetic)
        for f in range(3):
            hermitian(f'D0_{f}_{f}', -controls.hopping)
            hermitian(f'D1_{f}_{f}', -controls.hopping)
            hermitian(f'W_{f}_{f}', controls.higgs)
        if any(v != numerators.get((j, i), 0) for (i, j), v in numerators.items()):
            raise ArithmeticError('weighted Hermitian source matrix required')
        self.numerators = numerators
        self.columns = defaultdict(list)
        row_bounds = [F(0)]*self.dimension
        for (i, j), value in numerators.items():
            self.columns[j].append((i, value/self.norms[i]))
            row_bounds[i] += sqrt_upper(value*value/(self.norms[i]*self.norms[j]))
        # A Hermitian matrix's spectral norm is at most its largest absolute
        # row sum. Every radical in this row sum is rounded upward exactly.
        self.norm_bound = max(row_bounds)
        seeds = [i for i, x in enumerate(self.labels)
                 if x['N_H'] == 0 and x['link0_block'] == x['link1_block'] == -1
                 and x['e_at1'] == [] and x['e_at0'] == ([] if q == 0 else [0])]
        if len(seeds) != 1 or self.norms[seeds[0]] != 1:
            raise ArithmeticError('unique normalized reference-compensated seed required')
        self.seed = seeds[0]
        self.records = defaultdict(list)
        for i, x in enumerate(self.labels):
            self.records[x['N_H'], x['link0_block'], x['link1_block']].append(i)
        self.records = dict(sorted(self.records.items()))
        if sum(map(len, self.records.values())) != self.dimension:
            raise ArithmeticError('complete record partition required')

    def apply(self, vector):
        out = {}
        for j, x in vector.items():
            for i, value in self.columns[j]:
                add(out, i, value*x)
        return out

    def powers(self, order):
        out = [{self.seed:F(1)}]
        for _ in range(order):
            out.append(self.apply(out[-1]))
        return out

    def inner(self, left, right, indices=None):
        return sum((self.norms[i]*x*right.get(i, 0) for i, x in left.items()
                    if indices is None or i in indices), F(0))

    def block(self, r, nh=None):
        return {i for i, x in enumerate(self.labels)
                if x['link0_block'] == x['link1_block'] == r
                and (nh is None or x['N_H'] == nh)}

    def moments(self, r, nh=None):
        p = self.powers(5)
        keep = self.block(r, nh)
        dot = lambda i, j: self.inner(p[i], p[j], keep)
        return {2:dot(1, 1), 4:dot(2, 2)/4-dot(1, 3)/3,
                6:dot(1, 5)/60-dot(2, 4)/24+dot(3, 3)/36}


def source_replay():
    """Replay every basis vector and every term used in the full Hamiltonian."""
    basis, operators = inputs()
    carrier = NativeCarrier(ROOT)
    colored = []
    for r, b in enumerate(carrier.data['blocks']):
        if b['color_casimir']:
            values = {carrier.data['casimir'][i, i] for i in b['indices']}
            if len(values) != 1:
                raise ArithmeticError('uniform native Casimir required')
            colored.append({'block':r, 'charge6':int(b['q6Y']),
                            'dimension':b['dimension'], 'casimir':str(values.pop())})
    if colored != [dict(block=0, charge6=-4, dimension=3, casimir='4'),
                   dict(block=2, charge6=1, dimension=6, casimir='25/8'),
                   dict(block=3, charge6=2, dimension=3, casimir='3')]:
        raise ArithmeticError('native colored representation changed')
    sectors = {}
    checked = 0
    for q in (0, 1):
        sector = carrier.physical_basis(3, q)
        sectors[q] = sector
        if sector['labels'] != basis[str(q)]['labels'] or sector['dimension'] != basis[str(q)]['dimension']:
            raise ArithmeticError('source detector labels or dimension differ')
        for actual, saved in zip(sector['vectors'], basis[str(q)]['basis'], strict=True):
            vector = {(h, int(mask), l0, l1):F(c) for h, mask, l0, l1, c in saved['entries']}
            if actual != vector:
                raise ArithmeticError('source basis differs')
        if list(map(str, sector['norms'])) != basis[str(q)]['raw_norms']:
            raise ArithmeticError('source Gram norms differ')
        _, diagonal = diagonal_operators(carrier, sector)
        if {k:list(map(str, values)) for k, values in diagonal.items()} != operators[f'q{q}_diagonal']:
            raise ArithmeticError('source electric or occupation diagonal differs')
        terms = [('loop', carrier.loop)]
        for f in range(3):
            for edge in (0, 1):
                terms.append((f'D{edge}_{f}_{f}',
                              lambda v, e=edge, f=f: carrier.hopping(v, e, f, f, 3)))
            terms.append((f'W_{f}_{f}', lambda v, f=f: carrier.higgs(v, f, f, 3)))
        for name, apply in terms:
            exact, count = reduce_exact(sector, sector, apply, f'color readout q{q} {name}')
            expected = {(i, j):F(v) for i, j, v in operators[f'q{q}_{name}']}
            if exact != expected:
                raise ArithmeticError('literal source operator differs: '+name)
            checked += count
    exact, count = reduce_exact(sectors[0], sectors[1],
                               lambda v: carrier.addition(v, 0, 0, 3), 'seed addition')
    expected = {(i, j):F(v) for i, j, v in operators['addition_0_0']}
    if exact != expected:
        raise ArithmeticError('source reference-neutral addition differs')
    s0, s1 = Sector(0, data=(basis, operators)), Sector(1, data=(basis, operators))
    image = {i:v/s1.norms[i] for (i, j), v in exact.items() if j == s0.seed}
    if image != {s1.seed:F(1)}:
        raise ArithmeticError('the two seeds are not the same addition pair')
    return {'dimensions':{'0':438, '1':1596}, 'colored_blocks':colored,
            'basis_vectors_rebuilt':2034, 'operator_columns_replayed':checked+count,
            'seed_pair_is_actual_reference_neutral_addition':True,
            'raw_source_imports_and_bytes_verified':True}
