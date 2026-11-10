"""Exact full Gauss kernels for every uniform two-cell lepton sector.

Keys are (mask, link, hleft, hright), with six CAR modes per cell.
The basis is rational, has an exact coordinate left inverse, and never uses
a numerical rank tolerance. Raising maps retain their adjacent-weight image.
"""
from collections import defaultdict
from functools import lru_cache

import sympy as s

from ..model import active_terms, sequence
from ...physical_attachment.chiral import representation
from ...chiral_higgs.model import higgs_memory
from ...chiral_source.model import link_generator, link_labels


def _clean(values):
    return {key: s.expand(value) for key, value in values.items() if s.expand(value) != 0}


def _integer(value):
    value = s.simplify(value)
    if value.is_Integer is not True:
        raise ValueError('exact integral weight required')
    return int(value)


def boundary_charges(boundary):
    """Integer q means (0,q); a pair specifies both static endpoint charges."""
    charges = (0, boundary) if type(boundary) is int else boundary
    if not isinstance(charges, tuple) or len(charges) != 2 or any(type(q) is not int for q in charges):
        raise ValueError('boundary is an exact integer q or tuple (q_left,q_right)')
    return charges


@lru_cache(maxsize=1)
def current_data():
    indices, _, _, links = active_terms()
    native = representation()[2]
    _, scalar, _, _ = higgs_memory()
    positions = [0] + [link_labels().index((indices[a], indices[b])) for a, b in links[1:]]
    matter = tuple(s.ImmutableSparseMatrix(t.extract(indices, indices)) for t in native)
    edge = tuple(tuple(s.ImmutableSparseMatrix(link_generator(t, site).extract(positions, positions))
                       for t in native) for site in (0, 1))
    scalar = tuple(s.ImmutableSparseMatrix(t) for t in scalar)
    # Verify the actual representation convention instead of guessing an
    # entrywise raising operator from remembered Standard Model weights.
    for triple in (matter[8:11], scalar[8:11], edge[0][8:11], edge[1][8:11]):
        x, y, z = triple
        assert x == x.H and y == y.H and z == z.H
        assert x * y - y * x == s.I * z
        plus = x + s.I * y
        assert z * plus - plus * z == plus
        assert plus * plus.H - plus.H * plus == 2 * z
    for matrices in (matter, scalar, edge[0], edge[1]):
        assert all(t == s.zeros(t.rows) for t in matrices[:8])
    return matter, scalar, edge


@lru_cache(maxsize=1)
def _local_weight_states():
    matter, scalar, _ = current_data()
    y = [_integer(6 * matter[11][j, j]) for j in range(3)] * 2
    z = [_integer(2 * matter[10][j, j]) for j in range(3)] * 2
    sy = [_integer(6 * scalar[11][j, j]) for j in range(5)]
    sz = [_integer(2 * scalar[10][j, j]) for j in range(5)]
    number = (-1, 1, 1) * 2
    spin = (1, 1, 1, -1, -1, -1)
    local = defaultdict(list)
    for mask in range(64):
        occupied = [j for j in range(6) if mask >> j & 1]
        yy = sum(y[j] for j in occupied)
        zz = sum(z[j] for j in occupied)
        ell = sum(number[j] for j in occupied)
        ss = sum(spin[j] for j in occupied)
        for h in range(5):
            local[yy + sy[h], zz + sz[h]].append((mask, h, ell, ss))
    return {weight: tuple(values) for weight, values in local.items()}


@lru_cache(maxsize=None)
def _weight_blocks(boundary, left_two_T3=0, right_two_T3=0):
    left_charge, right_charge = boundary_charges(boundary)
    local = _local_weight_states()
    _, _, edge = current_data()
    blocks = defaultdict(list)
    for link in range(6):
        ly, ry = [_integer(6 * edge[site][11][link, link]) for site in (0, 1)]
        lz, rz = [_integer(2 * edge[site][10][link, link]) for site in (0, 1)]
        for ml, hl, ll, sl in local.get((-6 * left_charge - ly, left_two_T3 - lz), ()):
            for mr, hr, lr, sr in local.get((-6 * right_charge - ry, right_two_T3 - rz), ()):
                blocks[ll + lr, sl + sr].append((ml | (mr << 6), link, hl, hr))
    return {label: tuple(sorted(keys)) for label, keys in sorted(blocks.items())}


def block_labels(boundary):
    return tuple(_weight_blocks(boundary))


def diagonal_basis(boundary, L, twoSz):
    return _weight_blocks(boundary).get((L, twoSz), ())


def character_multiplicities(boundary):
    """Independent SU(2)xSU(2) singlet count from four actual weight spaces.

For each endpoint, multiplicity(j=0)=dim(T3=0)-dim(T3=1).
The product difference retains hypercharge, lepton number and source spin.
    """
    weights = {(a, b): _weight_blocks(boundary, a, b) for a in (0, 2) for b in (0, 2)}
    labels = set().union(*(set(value) for value in weights.values()))
    result = {}
    for label in sorted(labels):
        number = (len(weights[0, 0].get(label, ())) - len(weights[2, 0].get(label, ()))
                  - len(weights[0, 2].get(label, ())) + len(weights[2, 2].get(label, ())))
        assert number >= 0
        if number:
            result[label] = number
    return result


@lru_cache(maxsize=None)
def _operator_entries(site, current):
    matter, scalar, edge = current_data()
    if current == 'raising':
        matrices = (matter[8] + s.I * matter[9], scalar[8] + s.I * scalar[9],
                    edge[site][8] + s.I * edge[site][9])
    else:
        matrices = (matter[current], scalar[current], edge[site][current])
    return tuple(tuple((i, j, s.expand(v)) for (i, j), v in matrix.todok().items()
                       if s.expand(v) != 0) for matrix in matrices)


def _current_column(key, site, current, boundary=0):
    if site not in (0, 1):
        raise ValueError('site is 0 or 1')
    mask, link, hl, hr = key
    matter, scalar, edge = _operator_entries(site, current)
    result = defaultdict(lambda: s.Integer(0))
    offset = 6 * site
    for i, j, coefficient in matter:
        for spin in range(2):
            target = sequence(mask, ((offset + 3 * spin + i, True),
                                     (offset + 3 * spin + j, False)))
            if target is not None:
                result[target[0], link, hl, hr] += coefficient * target[1]
    for i, j, coefficient in scalar:
        if j == (hl if site == 0 else hr):
            result[mask, link, i if site == 0 else hl, hr if site == 0 else i] += coefficient
    for i, j, coefficient in edge:
        if j == link:
            result[mask, i, hl, hr] += coefficient
    if current == 11:
        result[key] += boundary_charges(boundary)[site]
    return _clean(result)


def local_current_column(key, site, current_index, boundary=0):
    """All twelve actual local Gauss currents, including the fixed boundary."""
    if current_index not in range(12):
        raise ValueError('current_index must be in range(12)')
    return _current_column(key, site, current_index, boundary)


def local_raising_column(key, site):
    return _current_column(key, site, 'raising')


@lru_cache(maxsize=None)
def raising_matrix(boundary, L, twoSz):
    """Return M and rows (site,targetkey), retaining both adjacent weights."""
    keys = diagonal_basis(boundary, L, twoSz)
    columns = []
    rows = set()
    for key in keys:
        column = {(site, target): coefficient for site in (0, 1)
                  for target, coefficient in local_raising_column(key, site).items()}
        columns.append(column)
        rows.update(column)
    rows = tuple(sorted(rows))
    row_index = {key: i for i, key in enumerate(rows)}
    entries = {(row_index[key], j): value for j, column in enumerate(columns)
               for key, value in column.items()}
    return s.ImmutableSparseMatrix(len(rows), len(keys), entries), rows


@lru_cache(maxsize=None)
def _physical_block(boundary, L, twoSz):
    m, rows = raising_matrix(boundary, L, twoSz)
    n = m.cols
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    row_columns = defaultdict(list)
    for (row, col), value in m.todok().items():
        row_columns[row].append((col, value))
    for entries in row_columns.values():
        first = entries[0][0]
        for col, _ in entries[1:]:
            parent[find(col)] = find(first)
    components = defaultdict(list)
    for col in range(n):
        components[find(col)].append(col)
    rows_by_component = defaultdict(list)
    for row, entries in row_columns.items():
        rows_by_component[find(entries[0][0])].append(row)
    entries = {}
    coordinate_rows = []
    certificates = []
    rank = 0
    for columns in sorted(components.values(), key=lambda values: values[0]):
        component_rows = rows_by_component[find(columns[0])]
        small = s.Matrix(m.extract(component_rows, columns))
        rref, pivots = small.rref()
        rank += len(pivots)
        free = [j for j in range(len(columns)) if j not in pivots]
        for free_column in free:
            output_column = len(coordinate_rows)
            coordinate_rows.append(columns[free_column])
            entries[columns[free_column], output_column] = s.Integer(1)
            for row, pivot in enumerate(pivots):
                value = -rref[row, free_column]
                if value:
                    entries[columns[pivot], output_column] = value
        certificates.append({'domain_columns': len(columns), 'image_rows': len(component_rows),
                             'exact_rank': len(pivots), 'nullity': len(free)})
    b = s.ImmutableSparseMatrix(n, len(coordinate_rows), entries)
    assert m * b == s.zeros(m.rows, b.cols)
    assert b.extract(coordinate_rows, range(b.cols)) == s.eye(b.cols)
    assert b.cols + rank == n
    gram = s.ImmutableSparseMatrix(b.H * b)
    return b, gram, tuple(coordinate_rows), {
        'boundary': boundary, 'boundary_hypercharges': list(boundary_charges(boundary)), 'L': L, 'twoSz': twoSz,
        'diagonal_dimension': n, 'raising_image_rows': len(rows),
        'exact_constraint_rank': rank, 'physical_dimension': b.cols,
        'constraint_components': certificates,
        'largest_component_domain': max((c['domain_columns'] for c in certificates), default=0),
        'exact_null_residual': True, 'coordinate_left_inverse_checked': True,
    }


def physical_basis(boundary, L, twoSz):
    """Return exact rational B and Gram B.H*B in the diagonal_basis order."""
    result = _physical_block(boundary, L, twoSz)
    return result[0], result[1]


def coordinate_rows(boundary, L, twoSz):
    """Selecting these rows is an exact left inverse of B on its image."""
    return _physical_block(boundary, L, twoSz)[2]


def physical_metadata(boundary, L, twoSz):
    return _physical_block(boundary, L, twoSz)[3]


def physical_dimension(boundary, L, twoSz):
    return physical_metadata(boundary, L, twoSz)['physical_dimension']


def complete_metadata(boundary):
    blocks = [physical_metadata(boundary, ell, spin) for ell, spin in block_labels(boundary)]
    character_counts = character_multiplicities(boundary)
    assert {(b['L'], b['twoSz']): b['physical_dimension'] for b in blocks if b['physical_dimension']} == character_counts
    return {'boundary': boundary, 'boundary_hypercharges': list(boundary_charges(boundary)), 'active_dimension': 614400,
            'diagonal_dimension': sum(block['diagonal_dimension'] for block in blocks),
            'physical_dimension': sum(block['physical_dimension'] for block in blocks),
            'blocks': blocks, 'rank_method': 'exact rational RREF of disjoint constraint components',
            'numerical_rank_tolerance_used': False,
            'independent_SU2_character_multiplicities_match': True}
