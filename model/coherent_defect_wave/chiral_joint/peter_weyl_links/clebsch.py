"""Complete native-current Clebsch intertwiners, derived by exact linear algebra."""
from collections import defaultdict
from functools import lru_cache

import sympy as s

from ...physical_attachment.chiral import representation
from ..periodic_family.model import native_blocks
from ..oriented_links.model import labels, label_index


@lru_cache(maxsize=1)
def irreps():
    currents = representation()[2]
    result = {0: {'dimension': 1, 'currents': tuple(s.ImmutableSparseMatrix([[0]]) for _ in currents),
                  'native_indices': (), 'orientation': 0}}
    for block in native_blocks():
        charge = int(6 * currents[-1][block[0], block[0]])
        restricted = tuple(s.ImmutableSparseMatrix(t.extract(block, block)) for t in currents)
        result[charge] = {'dimension': len(block), 'currents': restricted,
                          'native_indices': block, 'orientation': 1}
        result[-charge] = {'dimension': len(block),
            'currents': tuple(s.ImmutableSparseMatrix(-t.T) for t in restricted),
            'native_indices': block, 'orientation': -1}
    assert len(result) == 11 and sum(row['dimension'] ** 2 for row in result.values()) == 119
    return result


@lru_cache(maxsize=1)
def native_charges():
    return tuple(q for q, row in irreps().items() if row['orientation'] == 1)


@lru_cache(maxsize=512)
def carrier_index(charge, a, b):
    if charge == 0:
        if a != 0 or b != 0:
            raise ValueError('vacuum representation has one row')
        return 0
    rep = irreps()[charge]; block = rep['native_indices']
    return label_index()[rep['orientation'], block[a], block[b]]


@lru_cache(maxsize=128)
def intertwiner(r, source, target):
    """Every isometric embedding V_target -> V_r tensor V_source.

    Exact Cartan weights remove necessarily zero entries. All twelve current
    equations are then solved on every remaining coefficient.
    """
    if r not in irreps() or source not in irreps() or target not in irreps():
        raise ValueError('retained irreducible representation required')
    rr, rs, rt = (irreps()[q] for q in (r, source, target))
    dr, ds, dt = (rep['dimension'] for rep in (rr, rs, rt))
    if r + source != target:
        return (), {'cartan_variables': 0, 'constraint_rank': 0, 'multiplicity': 0,
                    'hypercharge_rejected': True}
    product_currents = tuple(s.ImmutableSparseMatrix(
        s.kronecker_product(a, s.eye(ds)) + s.kronecker_product(s.eye(dr), b))
        for a, b in zip(rr['currents'], rs['currents'], strict=True))
    variables = tuple((row, col) for row in range(dr * ds) for col in range(dt)
        if all(product_currents[j][row, row] == rt['currents'][j][col, col]
               for j in (6, 7, 10, 11)))
    if not variables:
        return (), {'cartan_variables': 0, 'constraint_rank': 0, 'multiplicity': 0,
                    'hypercharge_rejected': False}
    equations = defaultdict(lambda: defaultdict(lambda: s.Integer(0)))
    for current, (left, right) in enumerate(zip(product_currents, rt['currents'], strict=True)):
        left_columns = defaultdict(list); right_rows = defaultdict(list)
        for (row, col), value in left.todok().items():
            left_columns[col].append((row, value))
        for (row, col), value in right.todok().items():
            right_rows[row].append((col, value))
        for variable, (row, col) in enumerate(variables):
            for out_row, value in left_columns[row]:
                equations[current, out_row, col][variable] += value
            for out_col, value in right_rows[col]:
                equations[current, row, out_col][variable] -= value
    entries = {}; count = 0
    for equation in equations.values():
        terms = {j: s.expand(value) for j, value in equation.items() if s.expand(value) != 0}
        if terms:
            for j, value in terms.items():
                entries[count, j] = value
            count += 1
    constraints = s.ImmutableSparseMatrix(count, len(variables), entries)
    nullspace = constraints.nullspace()
    embeddings = []
    for vector in nullspace:
        raw = s.ImmutableSparseMatrix(dr * ds, dt,
            {entry: value for entry, value in zip(variables, vector, strict=True) if value})
        for previous in embeddings:
            raw -= previous * s.simplify(s.trace(previous.H * raw) / dt)
        gram = raw.H * raw; norm = s.simplify(s.trace(gram) / dt)
        assert gram == norm * s.eye(dt) and norm.is_positive is True
        c = s.ImmutableSparseMatrix(raw / s.sqrt(norm))
        assert c.H * c == s.eye(dt)
        for previous in embeddings:
            assert previous.H * c == s.zeros(dt)
        for left, right in zip(product_currents, rt['currents'], strict=True):
            assert left * c == c * right
        embeddings.append(c)
    return tuple(embeddings), {'cartan_variables': len(variables),
        'constraint_rank': len(variables) - len(nullspace),
        'multiplicity': len(embeddings), 'hypercharge_rejected': False}


@lru_cache(maxsize=1)
def all_channels():
    channels = {}; census = []
    for r in native_charges():
        for source in irreps():
            target = r + source
            if target not in irreps():
                census.append({'field_charge': r, 'source_charge': source,
                               'target_charge': target, 'retained_target': False,
                               'multiplicity': 0})
                continue
            embeddings, checks = intertwiner(r, source, target)
            census.append({'field_charge': r, 'source_charge': source,
                           'target_charge': target, 'retained_target': True, **checks})
            if embeddings:
                channels[r, source, target] = embeddings
    return channels, tuple(census)


@lru_cache(maxsize=128)
def coordinate(r, a, b):
    """P M_(D^r_ab) P, with Haar coefficient sqrt(dim(source)/dim(target))."""
    if r not in native_charges() or a not in range(irreps()[r]['dimension']) or b not in range(irreps()[r]['dimension']):
        raise ValueError('one native coordinate matrix entry required')
    out = defaultdict(lambda: s.Integer(0))
    for (field, source, target), embeddings in all_channels()[0].items():
        if field != r:
            continue
        ds, dt = irreps()[source]['dimension'], irreps()[target]['dimension']
        factor = s.sqrt(s.Rational(ds, dt))
        for c in embeddings:
            left = [(m, i, c[a * ds + m, i]) for m in range(ds) for i in range(dt)
                    if c[a * ds + m, i]]
            right = [(n, j, c[b * ds + n, j]) for n in range(ds) for j in range(dt)
                     if c[b * ds + n, j]]
            for m, i, value in left:
                for n, j, other in right:
                    out[carrier_index(target, i, j), carrier_index(source, m, n)] += factor * value * s.conjugate(other)
    return s.ImmutableSparseMatrix(119, 119, {key: s.simplify(value) for key, value in out.items() if value})


@lru_cache(maxsize=256)
def link(a, b):
    for charge in native_charges():
        block = irreps()[charge]['native_indices']
        if a in block and b in block:
            return coordinate(charge, block.index(a), block.index(b))
    raise ValueError('native indices must belong to the same irreducible block')


@lru_cache(maxsize=5)
def transporter(charge):
    d = irreps()[charge]['dimension']; entries = {}
    for a in range(d):
        for b in range(d):
            for (row, col), value in coordinate(charge, a, b).todok().items():
                entries[119 * a + row, 119 * b + col] = value
    return s.ImmutableSparseMatrix(119 * d, 119 * d, entries)
