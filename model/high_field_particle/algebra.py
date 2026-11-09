"""Full-source finite checks of the high-field one-flip expansion.

The coefficient calculation uses exact matrix entries and all intermediate
states. It supplies perturbative coefficients, not a finite-coupling remainder
or a physical particle assignment. The infinite-volume band theorem is an
explicit external input in the owning proof.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import product

import numpy as np

from model.native_ordered_ground import source


def _graph(sites, edges):
    source.parameters(sites, edges, F(1), F(1))
    if sites > 8:
        raise ValueError('full source checks are limited to eight sites')
    adjacency = [[0] * sites for _ in range(sites)]
    for x, y in edges:
        adjacency[x][y] += 1
        adjacency[y][x] += 1
    colors = {}
    for start in range(sites):
        if start in colors:
            continue
        colors[start] = 1
        queue = [start]
        for x in queue:
            for y, multiplicity in enumerate(adjacency[x]):
                if multiplicity and y not in colors:
                    colors[y] = -colors[x]
                    queue.append(y)
    return tuple(tuple(row) for row in adjacency), tuple(colors[x] for x in range(sites))


def _walsh(sites):
    dimension = 1 << sites
    return np.array([[1 - 2 * ((x & y).bit_count() % 2)
                      for y in range(dimension)] for x in range(dimension)],
                    dtype=np.int64)


def transformed_source(sites, edges, coupling=F(1), field=F(1)):
    """Return the exact integer global-Hadamard source and its denominator."""
    _graph(sites, edges)
    native, scale = source.matrix(sites, edges, coupling, field)
    dimension = 1 << sites
    walsh = _walsh(sites)
    transformed = walsh @ native @ walsh
    if np.any(transformed % dimension):
        raise ArithmeticError('the source frame has a nonintegral numerator')
    transformed //= dimension
    expected = np.zeros_like(native)
    bond, selection = coupling * scale / 2, field * scale / 2
    if bond.denominator != 1 or selection.denominator != 1:
        raise ArithmeticError('source denominator does not represent its coefficients')
    for mask in range(dimension):
        expected[mask, mask] -= int(selection) * (sites - 2 * mask.bit_count())
        for x, y in edges:
            if ((mask >> x) & 1) == ((mask >> y) & 1):
                expected[mask, mask] -= int(bond)
                expected[mask ^ (1 << x) ^ (1 << y), mask] += int(bond)
    if not np.array_equal(transformed, expected):
        raise ArithmeticError('the native source disagrees with XX-YY-ZZ-I')
    transformed.flags.writeable = False
    return transformed, scale


def frame_audit(sites, edges, coupling=F(2, 7), field=F(3, 5)):
    """Check source frame, full parity/charge and each local Z kick exactly."""
    adjacency, signs = _graph(sites, edges)
    h, scale = transformed_source(sites, edges, coupling, field)
    dimension = 1 << sites
    indices = np.arange(dimension)
    parity = np.array([1 - 2 * (mask.bit_count() % 2) for mask in indices])
    twice_charge = np.array([sum(signs[x] * (1 - 2 * ((int(mask) >> x) & 1))
                                 for x in range(sites)) for mask in indices])
    if np.any(h * (parity[:, None] - parity[None, :])):
        raise ArithmeticError('full source parity is not conserved')
    if np.any(h * (twice_charge[:, None] - twice_charge[None, :])):
        raise ArithmeticError('full staggered source charge is not conserved')
    walsh = _walsh(sites)
    for x in range(sites):
        z_diagonal = 1 - 2 * ((indices >> x) & 1)
        if not np.array_equal(walsh * z_diagonal[None, :],
                              walsh[indices ^ (1 << x), :]):
            raise ArithmeticError('local source Z does not become a one-site X kick')
    return {'sites': sites, 'all_source_rows': dimension, 'edge_count': len(edges),
            'distinct_edges': sum(bool(adjacency[x][y]) for x in range(sites)
                                  for y in range(x + 1, sites)),
            'coupling': str(coupling), 'selection_field': str(field),
            'source_denominator': scale, 'exact_global_Hadamard_identity': True,
            'full_source_parity_conserved': True,
            'full_source_staggered_charge_conserved': True,
            'local_source_Z_equals_transformed_X': True,
            'one_flip_charge_relative_to_product_vacuum': [-s for s in signs],
            'one_flip_parity_relative_to_product_vacuum': -1,
            'finite_coupling_local_kick_is_exact_band_projection': False}


def _zero_matrix(rows, columns):
    return [[F(0) for _ in range(columns)] for _ in range(rows)]


def expansion(sites, edges):
    """Exact coefficients in H/B=Nflip+epsilon V, epsilon=J/B.

The source scalar -BN/2 is removed. Both returned relative coefficients
subtract the corresponding perturbative ground energy. Repeated edge entries
are physical parallel bonds and their amplitudes add before squaring.
"""
    adjacency, signs = _graph(sites, edges)
    h, scale = transformed_source(sites, edges)
    if scale != 2:
        raise ArithmeticError('unit source has an unexpected denominator')
    dimension = 1 << sites
    twice_v = h.copy()
    for mask in range(dimension):
        twice_v[mask, mask] += sites - 2 * mask.bit_count()
    one = tuple(1 << x for x in range(sites))
    one_set = set(one)
    ground_first = F(int(twice_v[0, 0]), 2)
    ground_second = sum((F(int(twice_v[mask, 0]) ** 2, -4 * mask.bit_count())
                         for mask in range(1, dimension)), F(0))
    absolute_first = [[F(int(twice_v[x, y]), 2) for y in one] for x in one]
    absolute_second = _zero_matrix(sites, sites)
    first_wave = _zero_matrix(dimension, sites)
    virtual = set()
    for mask in range(dimension):
        if mask in one_set:
            continue
        denominator = 1 - mask.bit_count()
        for x, column in enumerate(one):
            entry = int(twice_v[mask, column])
            first_wave[mask][x] = F(entry, 2 * denominator)
            if entry:
                virtual.add(mask)
            for y, row in enumerate(one):
                absolute_second[y][x] += F(int(twice_v[row, mask]) * entry,
                                           4 * denominator)
    if any(mask.bit_count() != 3 for mask in virtual):
        raise ArithmeticError('the one-flip source has a non-three-flip virtual state')
    first = [[value - (ground_first if x == y else 0)
              for y, value in enumerate(row)] for x, row in enumerate(absolute_first)]
    second = [[value - (ground_second if x == y else 0)
               for y, value in enumerate(row)] for x, row in enumerate(absolute_second)]
    degree = tuple(sum(row) for row in adjacency)
    squared_degree = tuple(sum(w * w for w in row) for row in adjacency)
    expected_first = [[F(degree[x], 2) if x == y else F(0)
                       for y in range(sites)] for x in range(sites)]
    expected_second = [[F(2 * squared_degree[x] * (x == y)
                           - sum(adjacency[x][z] * adjacency[z][y] for z in range(sites)), 8)
                        for y in range(sites)] for x in range(sites)]
    if first != expected_first or second != expected_second:
        raise ArithmeticError('full intermediate-state sum disagrees with the graph formula')
    if ground_first != -F(len(edges), 2) or ground_second != -F(sum(squared_degree), 16):
        raise ArithmeticError('full source ground subtraction disagrees with the graph formula')
    if any(second[x][y] and signs[x] != signs[y]
           for x in range(sites) for y in range(sites)):
        raise ArithmeticError('the effective matrix mixes opposite one-flip charges')

    # Intermediate normalization P W_j=0. This checks the full order-two
    # intertwining identity, including the five-flip rows of W_2.
    neighbors = [tuple((row, F(int(twice_v[row, col]), 2))
                       for row in np.flatnonzero(twice_v[:, col]))
                 for col in range(dimension)]
    v_wave = _zero_matrix(dimension, sites)
    for col in range(dimension):
        for row, value in neighbors[col]:
            for x in range(sites):
                v_wave[row][x] += value * first_wave[col][x]
    second_wave = _zero_matrix(dimension, sites)
    for mask in range(dimension):
        for x, column in enumerate(one):
            order_one = (mask.bit_count() - 1) * first_wave[mask][x] + F(int(twice_v[mask, column]), 2)
            if mask in one_set:
                order_one -= absolute_first[one.index(mask)][x]
            if order_one:
                raise ArithmeticError('first-order full-space intertwiner failed')
            value = v_wave[mask][x] - sum(first_wave[mask][y] * absolute_first[y][x]
                                        for y in range(sites))
            if mask in one_set:
                if value != absolute_second[one.index(mask)][x]:
                    raise ArithmeticError('second-order projected intertwiner failed')
            else:
                second_wave[mask][x] = value / (1 - mask.bit_count())
                if (mask.bit_count() - 1) * second_wave[mask][x] + value:
                    raise ArithmeticError('second-order complementary intertwiner failed')
    return {'sites': sites, 'edges': edges, 'degree': degree,
            'squared_multiplicity_degree': squared_degree,
            'ground_first': ground_first, 'ground_second': ground_second,
            'absolute_first': tuple(map(tuple, absolute_first)),
            'absolute_second': tuple(map(tuple, absolute_second)),
            'relative_first': tuple(map(tuple, first)),
            'relative_second': tuple(map(tuple, second)),
            'three_flip_intermediate_masks': tuple(sorted(virtual)),
            'all_complementary_rows_checked': dimension - sites,
            'full_order_two_intertwiner_exact': True,
            'second_wave_maximum_flip_count': max((mask.bit_count() for mask in range(dimension)
                                                   if any(second_wave[mask])), default=0)}


def periodic_edges(shape):
    """Positive-axis periodic bonds; length-two axes have explicit duplicates."""
    if (type(shape) is not tuple or not shape or
            any(type(n) is not int or n < 2 or n % 2 for n in shape)):
        raise ValueError('positive even periodic side lengths at least two required')
    volume = 1
    for n in shape:
        volume *= n
    if volume > 8:
        raise ValueError('finite periodic source check is limited to eight sites')
    vertices = tuple(product(*(range(n) for n in shape)))
    index = {x: i for i, x in enumerate(vertices)}
    edges = []
    for x in vertices:
        for axis, length in enumerate(shape):
            y = list(x)
            y[axis] = (y[axis] + 1) % length
            edges.append(tuple(sorted((index[x], index[tuple(y)]))))
    return tuple(edges)


def cubic_paths():
    """Exact two-edge paths in Z^3; no finite-torus or matrix extrapolation."""
    directions = tuple(tuple(sign if j == axis else 0 for j in range(3))
                       for axis in range(3) for sign in (-1, 1))
    walks = Counter(tuple(x[j] + y[j] for j in range(3))
                    for x in directions for y in directions)
    coefficients = {x: F((12 if x == (0, 0, 0) else 0) - count, 8)
                    for x, count in walks.items()}
    hessian = tuple(tuple(-sum(value * x[i] * x[j] for x, value in coefficients.items())
                          for j in range(3)) for i in range(3))
    if sum(walks.values()) != 36 or len(walks) != 19:
        raise ArithmeticError('cubic two-edge support count failed')
    if hessian != tuple(tuple(F(3 if i == j else 0) for j in range(3)) for i in range(3)):
        raise ArithmeticError('cubic second-order Hessian failed')
    return {'walks': dict(sorted(walks.items())),
            'second_order_coefficients': dict(sorted(coefficients.items())),
            'at_zero_momentum': sum(coefficients.values(), F(0)),
            'hessian_at_zero_momentum': hessian,
            'first_order_constant': F(3)}


def _json_exact(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_exact(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_exact(item) for item in value]
    return value


def certificate():
    """JSON-safe exact controls; no spectral fit or measured input."""
    graphs = [('isolated_site', 1, ()), ('bond', 2, ((0, 1),)),
              ('path', 3, ((0, 1), (1, 2))),
              ('periodic_cycle_6', 6, periodic_edges((6,))),
              ('open_cube', 8, source.cube()[0]),
              ('periodic_cube_declared_parallel_bonds', 8, periodic_edges((2, 2, 2)))]
    controls = []
    for name, sites, edges in graphs:
        controls.append({'name': name, 'frame': frame_audit(sites, edges),
                         'expansion': expansion(sites, edges)})
    return _json_exact({
        'finite_source_controls': controls,
        'bulk_unit_edge_cubic_paths': cubic_paths(),
        'dimensionless_parameter': 'epsilon=J/B',
        'removed_source_scalar': '-B*sites/2',
        'relative_first_formula': 'diag(degree)/2',
        'relative_second_formula': '(2*diag(sum_edge_multiplicity_squared)-A^2)/8',
        'cubic_energy_through_second_order':
            'B+3J+(J^2/(2B))*(3-(cos(k1)+cos(k2)+cos(k3))^2)',
        'cubic_hessian_through_second_order': '3*J^2/B times I3',
        'first_order_hopping_zero': True,
        'infinite_lattice_band_uses_external_Yarotsky_smallness_theorem': True,
        'explicit_Yarotsky_threshold_supplied': False,
        'finite_coupling_remainder_computed_by_this_certificate': False,
        'specified_positive_selection_field_required': True,
        'joint_n_inverse_fifth_field_trajectory_used': False,
        'finite_periodic_double_edges_identified_with_bulk_unit_edges': False,
        'native_ground_preparation_executed_here': False,
        'measured_particle_inputs_used': False,
        'physical_particle_assignment': False,
        'physical_mass_prediction': False,
        'relativistic_dispersion_claimed': False})
