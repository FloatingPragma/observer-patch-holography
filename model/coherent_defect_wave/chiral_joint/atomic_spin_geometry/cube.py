"""Full eight-spin spectrum and an exact on-site rotation of the bond axes."""
from functools import lru_cache
import sympy as s


def cube_edges():
    return tuple((x, x ^ bit, axis) for x in range(8)
                 for axis, bit in enumerate((4, 2, 1)) if x < (x ^ bit))


def pauli(site, axis):
    rows = {}
    for column in range(256):
        bit = (column >> site) & 1
        if axis == 2:
            rows[column, column] = 1 - 2 * bit
        else:
            rows[column ^ (1 << site), column] = 1 if axis == 0 else s.I * (1 - 2 * bit)
    return s.ImmutableSparseMatrix(256, 256, rows)


@lru_cache(maxsize=1)
def full_operators():
    spins = [[pauli(site, axis) for axis in range(3)] for site in range(8)]
    plain = s.SparseMatrix(256, 256, {}); twisted = s.SparseMatrix(256, 256, {})
    for left, right, axis in cube_edges():
        for j in range(3):
            term = spins[left][j] * spins[right][j]
            plain += term
            twisted += (1 if j == axis else -1) * term
    unitary = s.SparseMatrix.eye(256)
    for site in range(8):
        for axis, bit in enumerate((4, 2, 1)):
            if site & bit:
                unitary = s.I * spins[site][axis] * unitary
    unitary = s.ImmutableSparseMatrix(unitary)
    assert unitary.H * unitary == s.SparseMatrix.eye(256)
    assert unitary * plain * unitary.H == twisted
    total = [sum((spins[site][j] for site in range(8)), s.SparseMatrix(256, 256, {}))/2 for j in range(3)]
    assert all(plain * spin == spin * plain for spin in total)
    return s.ImmutableSparseMatrix(plain), s.ImmutableSparseMatrix(twisted), unitary


def magnetization_block(down):
    if type(down) is not int or not 0 <= down <= 8:
        raise ValueError('a complete magnetization sector from zero through eight is required')
    basis = tuple(mask for mask in range(256) if mask.bit_count() == down)
    positions = {mask: i for i, mask in enumerate(basis)}
    matrix = s.zeros(len(basis))
    for column, mask in enumerate(basis):
        matrix[column, column] = sum((1-2*((mask >> left)&1))*(1-2*((mask >> right)&1))
                                    for left, right, _ in cube_edges())
        for left, right, _ in cube_edges():
            if ((mask >> left) ^ (mask >> right)) & 1:
                matrix[positions[mask ^ (1 << left) ^ (1 << right)], column] += 2
    return basis, s.ImmutableMatrix(matrix)


@lru_cache(maxsize=1)
def spectrum():
    x = s.Symbol('x'); factors = {}; blocks = []
    for down in range(9):
        basis, matrix = magnetization_block(down)
        entries = []
        for factor, multiplicity in s.factor_list(matrix.charpoly(x).as_expr())[1]:
            multiplicity = int(multiplicity)
            polynomial = s.Poly(factor, x).monic()
            key = str(polynomial.as_expr())
            if key not in factors:
                factors[key] = {'polynomial': polynomial, 'sectors': {}}
            factors[key]['sectors'][down] = multiplicity
            entries.append({'factor_of_4H': key, 'multiplicity': multiplicity})
        assert sum(s.degree(s.sympify(row['factor_of_4H']), x)*row['multiplicity']
                   for row in entries) == len(basis)
        blocks.append({'down_spins': down, 'twice_magnetization': 8-2*down,
                       'dimension': len(basis), 'complete_characteristic_factors': entries})
    rows = []
    for key, data in factors.items():
        for root_index, (interval, mult) in enumerate(data['polynomial'].intervals(eps=s.Rational(1, 10**24))):
            assert mult == 1
            low, high = (v/4 for v in interval)
            rows.append({'factor_of_4H': key, 'root_index': root_index,
                         'energy_interval': [str(low), str(high)],
                         'full_multiplicity': sum(data['sectors'].values()),
                         'magnetization_multiplicities': {str(k): v for k, v in data['sectors'].items()}})
    rows.sort(key=lambda row: s.Rational(row['energy_interval'][0]))
    assert sum(row['full_multiplicity'] for row in rows) == 256
    for left, right in zip(rows, rows[1:]):
        assert s.Rational(left['energy_interval'][1]) < s.Rational(right['energy_interval'][0])
    ground, first = rows[:2]
    assert ground['full_multiplicity'] == 1 and ground['magnetization_multiplicities'] == {'4': 1}
    assert first['energy_interval'] == ['-4', '-4']
    assert first['full_multiplicity'] == 3 and first['magnetization_multiplicities'] == {'3': 1, '4': 1, '5': 1}
    gap = [-4-s.Rational(ground['energy_interval'][1]), -4-s.Rational(ground['energy_interval'][0])]
    assert gap[0] > s.Rational(4, 5)
    return {'Hamiltonian': 'H_cube=sum_12_undirected_edges S_x dot S_y',
            'full_spin_dimension': 256, 'complete_magnetization_blocks': blocks,
            'complete_spectrum': rows, 'ground_rank': 1, 'ground_total_spin': '0',
            'ground_energy_interval': ground['energy_interval'],
            'first_excited_rank': 3, 'first_excited_total_spin': '1',
            'first_excited_energy': '-4', 'gap_interval': [str(v) for v in gap],
            'source_second_order_operator': '-24 A_edge I + 2 J_edge H_cube',
            'spin_labels_are_rotated_atomic_spin_labels': True,
            'physical_relativistic_particle_or_GeV_mass': False}
