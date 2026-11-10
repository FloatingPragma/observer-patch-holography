"""Exact graph, equitable intertwiner and physical odd cyclic operator."""
from collections import Counter, defaultdict, deque
from functools import lru_cache
from pathlib import Path
import hashlib
import json

import numpy as np
from scipy.sparse import coo_matrix
import sympy as s

from ..family_cell.model import column, exact_ground
from ..family_cell.quotient import add, mul, divide, evaluate, sqrt_interval, outward

PARTITION = Path(__file__).with_name('equitable_partition.npz')


def seed():
    """Unnormalized J_up^dagger Omega; squared norm two, both terms retained."""
    return {(1 << 11, 1): 1, (1 << 10, 2): -1}


@lru_cache(maxsize=1)
def graph():
    source = seed()
    seen = set(source); front = set(source); depth = 0
    while front:
        more = {target for key in front for target in column(key)} - seen
        if not more:
            break
        seen.update(more); front = more; depth += 1
    keys = tuple(sorted(seen)); index = {key: j for j, key in enumerate(keys)}
    rows = []; cols = []; values = []
    for j, key in enumerate(keys):
        for target, value in column(key).items():
            if target not in index:
                raise ValueError('full local Hamiltonian leaves specified graph')
            rows.append(index[target]); cols.append(j); values.append(value)
    h = coo_matrix((np.array(values, dtype=np.int64), (rows, cols)),
                   shape=(len(keys), len(keys))).tocsr()
    h.sort_indices()
    assert (h - h.T).nnz == 0
    origin = index[next(iter(source))]
    phase = {origin: 1}; todo = deque((origin,))
    while todo:
        j = todo.popleft()
        for offset in range(h.indptr[j], h.indptr[j + 1]):
            i, value = int(h.indices[offset]), int(h.data[offset])
            if i == j:
                continue
            desired = -phase[j] * (1 if value > 0 else -1)
            if i in phase:
                if phase[i] != desired:
                    raise ValueError('frustrated graph does not satisfy this Perron certificate')
            else:
                phase[i] = desired; todo.append(i)
    if len(phase) != len(keys):
        raise ValueError('physical seed component is disconnected')
    phases = np.array([phase[j] for j in range(len(keys))], dtype=np.int64)
    for key, value in source.items():
        assert phase[index[key]] * value == 1
    payload = json.dumps({'keys': keys, 'rows': h.indptr.tolist(),
        'cols': h.indices.tolist(), 'values': h.data.tolist(), 'seed': sorted(source.items())},
        separators=(',', ':')).encode()
    return keys, h, phases, depth, hashlib.sha256(payload).hexdigest()


def generate_partition():
    keys, h, phase, _, digest = graph()
    source = seed(); mapping = {}; colors = []
    for i, key in enumerate(keys):
        signature = int(key[1] != 0), int(phase[i]) * source.get(key, 0)
        if signature not in mapping:
            mapping[signature] = len(mapping)
        colors.append(mapping[signature])
    for _ in range(len(keys)):
        mapping = {}; next_colors = []
        for i in range(h.shape[0]):
            totals = defaultdict(int)
            for offset in range(h.indptr[i], h.indptr[i + 1]):
                j = int(h.indices[offset])
                totals[colors[j]] += int(h.data[offset]) * int(phase[i]) * int(phase[j])
            signature = colors[i], tuple(sorted((c, v) for c, v in totals.items() if v))
            if signature not in mapping:
                mapping[signature] = len(mapping)
            next_colors.append(mapping[signature])
        if len(set(colors)) == len(mapping):
            colors = next_colors
            break
        colors = next_colors
    else:
        raise ValueError('partition failed finite stabilization')
    np.savez_compressed(PARTITION, colors=np.array(colors, dtype=np.int32),
                        source_digest=np.array(digest))


@lru_cache(maxsize=1)
def quotient():
    keys, h, phase, depth, digest = graph()
    with np.load(PARTITION, allow_pickle=False) as archive:
        if str(archive['source_digest']) != digest:
            raise ValueError('equitable partition belongs to a different full source')
        colors = archive['colors'].copy()
    if colors.dtype != np.int32 or colors.shape != (len(keys),):
        raise ValueError('complete integer partition array required')
    classes = sorted(set(map(int, colors)))
    if classes != list(range(len(classes))):
        raise ValueError('partition classes must be contiguous and nonempty')
    sizes = Counter(map(int, colors)); rows = {}; source_values = {}; source = seed()
    for i, key in enumerate(keys):
        totals = defaultdict(int)
        for offset in range(h.indptr[i], h.indptr[i + 1]):
            j = int(h.indices[offset])
            totals[int(colors[j])] += int(h.data[offset]) * int(phase[i]) * int(phase[j])
        row = {c: v for c, v in totals.items() if v}; c = int(colors[i])
        if c in rows and rows[c] != row:
            raise ValueError('full-source equitable intertwiner fails')
        rows[c] = row
        value = int(phase[i]) * source.get(key, 0)
        if c in source_values and source_values[c] != value:
            raise ValueError('physical two-component seed is outside the quotient')
        source_values[c] = value
    q = s.ImmutableSparseMatrix(len(classes), len(classes),
        {(i, j): value for i, row in rows.items() for j, value in row.items()})
    g = s.ImmutableSparseMatrix(s.diag(*(sizes[i] for i in classes)))
    v = s.ImmutableSparseMatrix([source_values[i] for i in classes])
    assert g * q == q.T * g and (v.T * g * v)[0] == 2
    blank = tuple(i for i in classes if q[i, i] == 0)
    excited = tuple(i for i in classes if q[i, i] == 1)
    assert len(blank) + len(excited) == len(classes)
    assert q.extract(blank, blank) == s.zeros(len(blank))
    assert q.extract(excited, excited) == s.eye(len(excited))
    assert all(v[i] == 0 for i in blank)
    return q, g, v, blank, excited, depth


def modular_cyclic(matrix, vector, prime=1000003):
    """Propose pivot rows and verify column independence over a finite field."""
    a = np.array(matrix.tolist(), dtype=np.int64) % prime
    current = np.array(vector, dtype=np.int64).ravel() % prime
    pivots = []; basis = []
    if matrix.rows * (prime - 1) ** 2 >= 2 ** 63:
        raise ValueError('modular matrix product would exceed int64')
    for degree in range(matrix.rows + 1):
        residual = current.copy()
        for pivot, previous in zip(pivots, basis, strict=True):
            residual = (residual - int(residual[pivot]) * previous) % prime
        support = np.flatnonzero(residual)
        if not len(support):
            return tuple(pivots)
        pivot = int(support[0])
        residual = residual * pow(int(residual[pivot]), -1, prime) % prime
        pivots.append(pivot); basis.append(residual)
        current = a @ current % prime
    raise ValueError('finite cyclic sequence failed to close')


@lru_cache(maxsize=1)
def cyclic():
    q, g, source, blank, excited, _ = quotient()
    a = q.extract(excited, blank) * q.extract(blank, excited)
    metric = g.extract(excited, excited); v = source.extract(excited, (0,))
    assert metric * a == a.T * metric
    pivots = modular_cyclic(a, v); degree = len(pivots)
    columns = [v]
    for _ in range(degree):
        columns.append(a * columns[-1])
    c = s.Matrix.hstack(*columns[:-1])
    minor = c.extract(pivots, range(degree))
    coefficients = minor.inv(method='DM') * columns[-1].extract(pivots, (0,))
    assert c * coefficients == columns[-1]
    companion = s.zeros(degree)
    for j in range(degree - 1):
        companion[j + 1, j] = 1
    companion[:, degree - 1] = coefficients
    assert a * c == c * companion
    gram = c.T * metric * c / 2
    assert gram * companion == companion.T * gram
    x = s.Symbol('x')
    p = s.Poly(x ** degree - sum(coefficients[j] * x ** j for j in range(degree)), x)
    moments = [(v.T * metric * col)[0] / 2 for col in columns]
    numerator = s.Poly(sum(sum(p.nth(j) * moments[j - k - 1]
        for j in range(k + 1, degree + 1)) * x ** k for k in range(degree)), x)
    assert s.gcd(p, numerator).degree() == 0 and moments[0] == 1
    energy = s.Symbol('E')
    hp = s.Poly(s.cancel(p.as_expr().subs(x, energy * (energy - 1)) / energy), energy)
    residual = s.zeros(q.rows, 1)
    for coefficient in hp.all_coeffs():
        residual = q * residual + coefficient * source
    assert residual == s.zeros(q.rows, 1)
    assert len(modular_cyclic(q, source)) == hp.degree()
    return p, numerator, hp, tuple(moments), pivots, gram


@lru_cache(maxsize=1)
def spectrum():
    q, g, source, blank, excited, depth = quotient()
    p, numerator, hp, moments, pivots, gram = cyclic()
    roots = []
    for factor, multiplicity in s.factor_list(p)[1]:
        assert multiplicity == 1
        for interval, count in s.polys.polytools.intervals(factor, eps=s.Rational(1, 10 ** 90)):
            assert count == 1 and interval[0] >= 0
            roots.append(tuple(map(s.Rational, interval)))
    roots.sort(); records = []
    for root in roots:
        weight = divide(evaluate(numerator, root), evaluate(p.diff(), root))
        assert 0 < weight[0] <= weight[1] <= 1
        if root == (0, 0):
            records.append({'energy_interval': ['1', '1'],
                'probability_interval': outward(weight, 80), 'squared_root_interval': ['0', '0']})
            continue
        discriminant = sqrt_interval(add((1, 1), mul((4, 4), root)))
        inverse = divide((1, 1), discriminant)
        for sign in (-1, 1):
            energy = mul((s.Rational(1, 2),) * 2,
                add((1, 1), mul((sign, sign), discriminant)))
            split = mul((s.Rational(1, 2),) * 2,
                add((1, 1), mul((sign, sign), inverse)))
            probability = mul(weight, split)
            records.append({'energy_interval': outward(energy, 48),
                'probability_interval': outward(probability, 80),
                'squared_root_interval': outward(root, 48)})
    records.sort(key=lambda row: s.Rational(row['energy_interval'][0]))
    assert len(records) == hp.degree() == 49
    assert all(s.Rational(row['probability_interval'][0]) > 0 for row in records)
    assert sum(s.Rational(row['probability_interval'][0]) for row in records) <= 1
    assert sum(s.Rational(row['probability_interval'][1]) for row in records) >= 1
    vacuum = exact_ground()['ground_energy_interval']
    odd = records[0]['energy_interval']
    difference = [str(s.Rational(odd[0]) - s.Rational(vacuum[1])),
                  str(s.Rational(odd[1]) - s.Rational(vacuum[0]))]
    assert s.Rational(difference[1]) < 0
    keys, h, _, _, digest = graph()
    return {
        'full_local_dimension': 5 * 2 ** 30,
        'full_seed_component_dimension': len(keys),
        'full_component_H_nonzeros': int(h.nnz),
        'source_graph_sha256': digest,
        'graph_radius_from_two_seed_components': depth,
        'equitable_dimension': q.rows,
        'scalar_blank_classes': len(blank), 'scalar_excited_classes': len(excited),
        'squared_physical_cyclic_dimension': p.degree(),
        'full_physical_cyclic_dimension': hp.degree(),
        'squared_minimal_polynomial': str(p.as_expr()),
        'squared_polynomial_factorization': str(s.factor(p.as_expr())),
        'physical_resolvent_numerator': str(numerator.as_expr()),
        'full_H_minimal_polynomial': str(hp.as_expr()),
        'squared_moments': list(map(str, moments)),
        'independence_modulus': 1000003, 'independent_pivot_rows': list(pivots),
        'exact_full_graph_and_cyclic_intertwiners': True,
        'full_H_polynomial_residual': '0 on the physical two-component seed',
        'positive_Perron_overlap_of_physical_seed': True,
        'complete_energy_records': records,
        'outside_cyclic_sector_record': {'probability_interval': ['0', '0'], 'retained': True},
        'odd_component_ground_energy_interval': odd,
        'ground_filter_success_interval': records[0]['probability_interval'],
        'ground_filter_failure_interval': [str(1 - s.Rational(records[0]['probability_interval'][1])),
                                            str(1 - s.Rational(records[0]['probability_interval'][0]))],
        'vacuum_component_ground_comparison_interval': vacuum,
        'odd_minus_vacuum_component_ground_interval': difference,
        'global_ground_across_all_local_sectors_proved': False,
        'full_component_excitation_gap_inferred_from_visible_line_spacing': False,
        'seed_is_Ddagger_on_interacting_ground': False,
        'reference_is_product_ground_or_particle_mass': False,
        'numerical_eigensolver_in_certificate': False,
    }
