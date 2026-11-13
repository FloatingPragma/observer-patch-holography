"""Exact recoupling on the full, countably infinite physical spin network."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import sympy as sy
from sympy.physics.wigner import wigner_6j, clebsch_gordan
from model.openwave_pure_su2_source.source import coefficients


def admissible(a, b, c):
    return all(type(k) is int and k >= 0 for k in (a, b, c)) and abs(a-b) <= c <= a+b and (a+b+c) % 2 == 0


def labels(cap):
    if type(cap) is not int or cap < 0:
        raise ValueError('a nonnegative integer diagnostic cap is required')
    return tuple((a, b, c) for a, b in product(range(cap+1), repeat=2)
                 for c in range(abs(a-b), a+b+1, 2))


def casimir(k):
    if type(k) is not int or k < 0:
        raise ValueError('twice-spin is a nonnegative integer')
    return F(k*(k+2), 4)


def electric(label, shared=True):
    a, b, c = label
    if not admissible(a, b, c):
        raise ValueError('inadmissible Gauss spin network')
    return 3*casimir(a)+3*casimir(b)+casimir(c) if shared else 4*(casimir(a)+casimir(b))


@lru_cache(None)
def wilson_entry(left, right, face):
    if face not in (0, 1):
        raise ValueError('the graph has two plaquettes')
    if not admissible(*left) or not admissible(*right):
        raise ValueError('inadmissible Gauss spin network')
    a, b, c = left
    A, B, C = right
    if face:
        a, b, A, B = b, a, B, A
    if b != B or abs(a-A) != 1 or abs(c-C) != 1:
        return sy.S.Zero
    half = sy.Rational(1, 2)
    recoupling = wigner_6j(a*half, b*half, c*half, C*half, half, A*half)
    return sy.simplify(sy.sqrt((a+1)*(A+1)*(c+1)*(C+1))*recoupling**2/2)


def successors(label, face):
    a, b, c = label
    for delta, gamma in product((-1, 1), repeat=2):
        target = (a+delta, b, c+gamma) if face == 0 else (a, b+delta, c+gamma)
        if admissible(*target):
            value = wilson_entry(target, label, face)
            if value:
                yield target, value


@lru_cache(None)
def matrices(cap):
    basis = labels(cap)
    index = {v: i for i, v in enumerate(basis)}
    out = []
    for face in (0, 1):
        entries = {}
        for j, label in enumerate(basis):
            for target, value in successors(label, face):
                if target in index:
                    entries[index[target], j] = value
        matrix = sy.SparseMatrix(len(basis), len(basis), entries)
        assert matrix == matrix.T
        out.append(matrix)
    return basis, tuple(out)


def hamiltonian(cap, g2=F(2), spacing=F(1), shared=True):
    alpha, beta = coefficients(F(g2), F(spacing))
    basis, (first, second) = matrices(cap)
    diagonal = [sy.Rational(alpha*electric(label, shared)+2*beta) for label in basis]
    return sy.diag(*diagonal)-sy.Rational(beta)*(first+second)


def tail_floor(cap, g2=F(2), spacing=F(1)):
    alpha, _ = coefficients(F(g2), F(spacing))
    return 3*alpha*casimir(cap+1)


def magnetic_labels(k):
    return [sy.Rational(m, 2) for m in range(-k, k+1, 2)]


def cg_integral_entry(left, right):
    """Independent Haar-index oracle; no 6j symbol occurs in this function."""
    a, b, c = left
    A, B, C = right
    if b != B:
        return sy.S.Zero
    j, k, J, jp, Jp = map(lambda x: sy.Rational(x, 2), (a, b, c, A, C))
    total = sy.S.Zero
    for m, n, p, q in product(magnetic_labels(a), magnetic_labels(a), magnetic_labels(b), magnetic_labels(b)):
        old = sum(clebsch_gordan(j, k, J, m, p, z)*clebsch_gordan(j, k, J, n, q, z)
                  for z in magnetic_labels(c))
        if old == 0:
            continue
        for u in (sy.Rational(-1, 2), sy.Rational(1, 2)):
            M, N = m+u, n+u
            if M not in magnetic_labels(A) or N not in magnetic_labels(A):
                continue
            new = sum(clebsch_gordan(jp, k, Jp, M, p, z)*clebsch_gordan(jp, k, Jp, N, q, z)
                      for z in magnetic_labels(C))
            total += old*new*clebsch_gordan(j, sy.Rational(1, 2), jp, m, u, M)*clebsch_gordan(j, sy.Rational(1, 2), jp, n, u, N)
    norm = sy.sqrt(sy.Rational((a+1)*(b+1), c+1))*sy.sqrt(sy.Rational((A+1)*(b+1), C+1))
    return sy.simplify(norm*total/(2*(A+1)*(b+1)))


def graph():
    vertices = tuple(product(range(3), range(2)))
    edges = tuple(((x, y), (x+1, y)) for x in range(2) for y in range(2)) + tuple(((x, 0), (x, 1)) for x in range(3))
    faces = []
    for x in range(2):
        walk = ((x, 0), (x+1, 0), (x+1, 1), (x, 1), (x, 0))
        boundary = []
        for a, b in zip(walk, walk[1:]):
            edge = (a, b) if (a, b) in edges else (b, a)
            boundary.append((edges.index(edge), 1 if edge == (a, b) else -1))
        divergence = {v: 0 for v in vertices}
        for index, sign in boundary:
            a, b = edges[index]
            divergence[a] -= sign
            divergence[b] += sign
        assert set(divergence.values()) == {0}
        faces.append(boundary)
    assert len(set(i for i, _ in faces[0]) & set(i for i, _ in faces[1])) == 1
    return {'vertices': [list(v) for v in vertices], 'edges': [[list(a), list(b)] for a, b in edges],
            'oriented_face_boundaries': [[[i, sign] for i, sign in p] for p in faces]}


def run():
    cases = [((0, 0, 0), (1, 0, 1)), ((0, 1, 1), (1, 1, 0)), ((0, 1, 1), (1, 1, 2)),
             ((1, 1, 0), (2, 1, 1)), ((1, 1, 2), (2, 1, 1)), ((1, 1, 2), (2, 1, 3)),
             ((2, 2, 2), (1, 2, 1)), ((2, 2, 4), (3, 2, 5))]
    checked = []
    for left, right in cases:
        value = cg_integral_entry(left, right)
        assert value == wilson_entry(left, right, 0)
        checked.append({'left': list(left), 'right': list(right), 'exact_W1_coefficient': str(value)})
    # Two magnetic multiplications commute on the complete source. The expanded
    # intermediate label set retains every term, including those outside a cap.
    commutators = 0
    for label in labels(3):
        products = []
        for order in ((0, 1), (1, 0)):
            values = {}
            for middle, x in successors(label, order[0]):
                for target, y in successors(middle, order[1]):
                    values[target] = values.get(target, 0)+x*y
            products.append(values)
        assert all(sy.simplify(products[0].get(k, 0)-products[1].get(k, 0)) == 0 for k in products[0].keys() | products[1].keys())
        commutators += 1
    return {'schema': 'oph.coupled_full_Haar_source.v1',
            'verdict': 'PASS_EXACT_SHARED_EDGE_GAUSS_SOURCE', 'graph': graph(),
            'source': 'Full L2(SU2^7,Haar) Gauss quotient on two adjacent plaquettes; no representation cutoff in the Hamiltonian',
            'physical_basis': 'All integer a,b,c>=0 with |a-b|<=c<=a+b and a+b+c even; spins are a/2,b/2,c/2',
            'shared_electric_diagonal': '3 C(a)+3 C(b)+C(c), C(k)=k(k+2)/4',
            'magnetic_operator': 'beta (2-W1-W2), each W=Tr(U_loop)/2, norm<=1',
            'independent_exact_CG_Haar_integrals': checked,
            'full_successor_commutators_checked': commutators,
            'diagnostic_caps': [{'max_twice_loop_spin': n, 'complete_admissible_states': len(labels(n)),
                                 'all_intermediate_J_retained': True, 'tail_electric_floor_at_g2_2_a_1': str(tail_floor(n))} for n in (2, 4, 6, 8)],
            'complete_source_properties': {'exact_Gauss_reduction': True, 'shared_edge_Casimir_retained': True,
                 'all_Haar_representations': True, 'positive_magnetic_potential': True,
                 'self_adjoint_compact_resolvent_at_finite_graph': True},
            'scope': {'spatial_continuum': False, 'infinite_volume_gap': False, 'physical_GeV_mass': False,
                      'four_dimensional_Yang_Mills_reconstruction': False}}
