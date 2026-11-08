"""Exact local energy continuity for the selected finite spin action.

The currents refer to Hamiltonian evolution. An intervening recorded
instrument has its own energy change. No Lorentz stress tensor is inferred.
"""
from itertools import combinations

import sympy as sp


PAULI = (
    sp.Matrix([[0, 1], [1, 0]]),
    sp.Matrix([[0, -sp.I], [sp.I, 0]]),
    sp.diag(1, -1),
)


def commutator(a, b):
    return a * b - b * a


def spin(size, site, axis):
    """Spin one-half with tensor factors in ascending site order."""
    if (type(size) is not int or not 1 <= size <= 6
            or type(site) is not int or not 0 <= site < size
            or type(axis) is not int or not 0 <= axis < 3):
        raise ValueError('require one to six sites and valid site/axis')
    factors = [PAULI[axis] / 2 if x == site else sp.eye(2)
               for x in range(size)]
    return sp.kronecker_product(*factors)


def source(size, edges, *, coupling=sp.Rational(1), field=sp.Rational(0),
           staggered=None):
    """H_B with exact rational inputs; each unoriented bond occurs once.

    The staggered selection field is admitted only on a bipartite graph.
    The zero-field identity also applies to other finite graphs.
    """
    if type(size) is not int or not 1 <= size <= 6:
        raise ValueError('finite exact verifier requires one to six sites')
    if (isinstance(coupling, bool) or isinstance(field, bool)
            or not isinstance(coupling, (int, sp.Rational))
            or not isinstance(field, (int, sp.Rational))
            or coupling <= 0 or field < 0):
        raise ValueError('require exact rational J>0 and B>=0')
    edges = tuple(edges)
    if any(type(e) is not tuple or len(e) != 2
           or any(type(x) is not int for x in e)
           or not 0 <= e[0] < e[1] < size for e in edges):
        raise ValueError('bonds require distinct ordered valid endpoints')
    if len(set(edges)) != len(edges):
        raise ValueError('count each undirected bond once')
    if staggered is None:
        if field:
            raise ValueError('selection field requires bipartite signs')
        staggered = (1,) * size
    if (len(staggered) != size
            or any(type(x) is not int or x not in (-1, 1) for x in staggered)
            or (field and any(staggered[x] == staggered[y] for x, y in edges))):
        raise ValueError('invalid selection-field bipartite signs')
    eye = sp.eye(2 ** size)
    spins = [[spin(size, x, a) for a in range(3)] for x in range(size)]
    bond = {(x, y): coupling * (sum((spins[x][a] * spins[y][a]
                                      for a in range(3)), sp.zeros(2 ** size))
                                  - eye / 4)
            for x, y in edges}
    density = [-field * staggered[x] * spins[x][0] for x in range(size)]
    for (x, y), h in bond.items():
        density[x] += h / 2
        density[y] += h / 2
    hamiltonian = sum(density, sp.zeros(2 ** size))
    currents = {(x, y): sp.I * commutator(density[x], density[y])
                for x, y in combinations(range(size), 2)}
    return {'hamiltonian': hamiltonian, 'density': density,
            'currents': currents, 'bonds': bond, 'spin': spins}


def verify_continuity(size, edges, **kwargs):
    data = source(size, edges, **kwargs)
    h, e, currents = (data[k] for k in ('hamiltonian', 'density', 'currents'))
    zero = sp.zeros(2 ** size)
    for x in range(size):
        outgoing = sum((currents[x, y] if x < y else -currents[y, x]
                        for y in range(size) if y != x), zero)
        if sp.I * commutator(h, e[x]) + outgoing != zero:
            raise ArithmeticError('local energy continuity failed')
    if any(j.H != j for j in currents.values()):
        raise ArithmeticError('energy current is not Hermitian')
    tests = tuple(sp.Rational((x + 1) ** 2, x + 2) for x in range(size))
    weighted = sum((f * density for f, density in zip(tests, e, strict=True)), zero)
    flux = sum(((tests[y] - tests[x]) * j for (x, y), j in currents.items()), zero)
    if sp.I * commutator(h, weighted) != flux:
        raise ArithmeticError('weighted energy identity failed')
    adjacency = {x: set() for x in range(size)}
    for x, y in edges:
        adjacency[x].add(y)
        adjacency[y].add(x)
    for (x, y), current in currents.items():
        if y not in adjacency[x] and not adjacency[x].intersection(adjacency[y]):
            if current != zero:
                raise ArithmeticError('energy current exceeds graph range two')
    return {'sites': size, 'bonds': len(edges),
            'nonzero_site_current_pairs': sum(j != zero for j in currents.values()),
            'continuity_exact': True, 'weighted_identity_exact': True,
            'graph_range_at_most_two': True}


def exact_checks():
    rows = [verify_continuity(3, ((0, 1), (1, 2))),
            verify_continuity(3, ((0, 1), (0, 2), (1, 2))),
            verify_continuity(4, ((0, 1), (0, 2), (1, 3), (2, 3))),
            verify_continuity(3, ((0, 1), (1, 2)),
                              field=sp.Rational(3, 8), staggered=(1, -1, 1))]
    data = source(3, ((0, 1), (1, 2)), coupling=sp.Rational(3, 2))
    spins, bonds = data['spin'], data['bonds']
    chirality = sum((sp.LeviCivita(a, b, c) * spins[0][a] * spins[1][b] * spins[2][c]
                    for a in range(3) for b in range(3) for c in range(3)), sp.zeros(8))
    if sp.I * commutator(bonds[0, 1], bonds[1, 2]) != sp.Rational(9, 4) * chirality:
        raise ArithmeticError('oriented bond chirality coefficient changed')
    return {'finite_graphs': rows, 'oriented_bond_chirality_exact': True,
            'continuum_stress_or_momentum_tensor_constructed': False,
            'apparatus_energy_or_Einstein_geometry_constructed': False}
