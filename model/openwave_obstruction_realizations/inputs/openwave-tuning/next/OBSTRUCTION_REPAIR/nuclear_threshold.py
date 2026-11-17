"""Independent highest-weight CAR reconstruction and two-baryon audit.

No producer matrix, singlet basis, or inertia function is imported. The
physical source is the frozen four-site signed-hopping colour Hamiltonian.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import lcm
from pathlib import Path
import hashlib
import json
import numpy as np
import sympy as sy

HERE = Path(__file__).resolve().parent


def hop(word, dest, source):
    if not (word >> source) & 1:
        return None, 0
    sign = (-1) ** ((word & ((1 << source) - 1)).bit_count())
    word ^= 1 << source
    if (word >> dest) & 1:
        return None, 0
    sign *= (-1) ** ((word & ((1 << dest) - 1)).bit_count())
    return word | (1 << dest), sign


def current(word, sites, a, b):
    result = {}
    for x in sites:
        z, c = hop(word, 3*x+a, 3*x+b)
        if z is not None:
            result[z] = result.get(z, 0) + c
    return {z: c for z, c in result.items() if c}


def build(k):
    choices = list(combinations(range(4), k))
    words = sorted(sum(1 << (3*x+c) for c, xs in enumerate(xs3) for x in xs)
                   for xs3 in product(choices, repeat=3))
    index = {w: i for i, w in enumerate(words)}
    # Independent route: zero Cartan weight and both simple raising currents.
    equations = {}
    for j, word in enumerate(words):
        for a, b in [(0, 1), (1, 2)]:
            for out, coefficient in current(word, range(4), a, b).items():
                row = equations.setdefault((a, b, out), [0]*len(words))
                row[j] += coefficient
    constraints = sy.Matrix(list(equations.values())) if equations else sy.zeros(0, len(words))
    vectors = constraints.nullspace()
    R = np.array([[int(z*lcm(*(int(x.q) for x in v))) for z in v]
                  for v in vectors], dtype=np.int64).T
    integer_constraints = np.array(list(equations.values()), dtype=np.int64).reshape(-1, len(words))
    assert not np.any(integer_constraints@R)
    G = R.T@R
    K = np.zeros((len(words), len(words)), dtype=np.int64)
    mass = np.diag([sum((1 if x % 2 == 0 else -1)*((w >> (3*x)) & 7).bit_count()
                           for x in range(4)) for w in words])
    electric6 = np.zeros_like(K)
    for j, word in enumerate(words):
        for x in range(3):
            for c in range(3):
                for left, right in [(x, x+1), (x+1, x)]:
                    out, coefficient = hop(word, 3*left+c, 3*right+c)
                    if out is not None:
                        K[index[out], j] -= coefficient
        for cut in [1, 2, 3]:
            n = (word & ((1 << (3*cut))-1)).bit_count()
            electric6[j, j] -= n*n
            for a, b in product(range(3), repeat=2):
                for mid, c1 in current(word, range(cut), b, a).items():
                    for out, c2 in current(mid, range(cut), a, b).items():
                        electric6[index[out], j] += 3*c1*c2
    components = []
    for X in [mass, K, electric6]:
        assert np.array_equal(X, X.T)
        assert not np.any(integer_constraints@(X@R))
        components.append(R.T@X@R)
    return G, components, len(words)


def inertia(matrix):
    # Fraction-free symmetric Bareiss elimination. A common positive scale
    # preserves inertia; successive leading determinants give LDL signs.
    denominator = lcm(*(x.denominator for row in matrix for x in row))
    a = [[int(x*denominator) for x in row] for row in matrix]
    negative = 0
    previous = 1
    for k in range(len(a)):
        pivot = a[k][k]
        assert pivot != 0, 'A zero pivot requires a different exact decomposition.'
        negative += pivot*previous < 0
        for i in range(k+1, len(a)):
            for j in range(i, len(a)):
                numerator = pivot*a[i][j]-a[i][k]*a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator//previous
                a[j][i] = a[i][j]
        previous = pivot
    return negative


def main():
    paths = [HERE/'hadrons'/name for name in
             ['STAGGERED_CHAIN4_RESULT.json', 'RAW_COEFFICIENT_INPUTS.json',
              'CHAIN4_MANIFEST.json', 'INDUCED_INTERACTION_RESULT.json']]
    pins = {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    receipt = json.loads(paths[0].read_text())
    inputs = json.loads(paths[1].read_text())
    mass = F(receipt['rational_computational_mass'])
    error = F(receipt['uniform_mass_operator_error'])
    unit = F(inputs['electron_anchor_GeV']) / F(receipt['actual_mass_exact'])
    assert error == 12*abs(mass-F(receipt['actual_mass_exact']))
    sectors = []
    for k in range(5):
        sectors.append(build(k))
        print('reconstructed', k, len(sectors[-1][0]), flush=True)
    assert [len(s[0]) for s in sectors] == [1, 20, 50, 20, 1]
    checks, thresholds = [], []
    for row in receipt['results']:
        kappa = F(row['kappa'])
        E = {}
        for record, (G, components, weight_dimension) in zip(row['complete_sector_records'], sectors):
            n = record['N']
            M, K, C = components
            lo, hi = map(F, record['ground_interval_exact'])
            E[n] = lo, hi
            if n in [0, 12]:
                assert all(not np.any(X) for X in components) and lo == hi == 0
                continue
            counts = []
            for bound in [lo+error, hi-error]:
                A = [[mass*int(M[i,j])+F(1,8)*int(K[i,j])+kappa*F(int(C[i,j]),6)
                      -bound*int(G[i,j]) for j in range(len(G))] for i in range(len(G))]
                counts.append(inertia(A))
            assert counts == [0, 1]
            checks.append(dict(kappa=str(kappa), N=n, independent_dimension=len(G),
                               Cartan_weight_dimension=weight_dimension, inertia=counts))
            print('verified', str(kappa), n, flush=True)
        # Relative to the SAME actual N6 vacuum: excess over twice the B1 gap.
        # E12+E6-2E9; no arbitrary energy offset can affect this combination.
        excess = E[12][0]+E[6][0]-2*E[9][1], E[12][1]+E[6][1]-2*E[9][0]
        assert excess[0] > 0
        anti_excess = E[0][0]+E[6][0]-2*E[3][1], E[0][1]+E[6][1]-2*E[3][0]
        assert anti_excess == excess
        thresholds.append(dict(kappa=str(kappa),
                               excess_over_twice_same_box_single_gap_exact=list(map(str, excess)),
                               conditional_GeV=list(map(float, [v*unit for v in excess])),
                               below_same_box_two_particle_reference=False,
                               infinite_volume_binding_test=False))
    # Domain-wall term is a repulsive density coupling plus a one-body term.
    identities = []
    for nx, ny in product([0, 1], repeat=2):
        assert int(nx != ny) == nx+ny-2*nx*ny
        identities.append(dict(nx=nx, ny=ny, domainwall=int(nx != ny)))
    for path in paths:
        assert hashlib.sha256(path.read_bytes()).hexdigest() == pins[str(path.relative_to(HERE))]
    out = dict(schema='oph.obstruction_repair.nuclear_threshold_independent.v1',
               independently_reconstructed_CAR_and_highest_weight_kernel=True,
               exact_ground_brackets_verified=len(checks), checks=checks, thresholds=thresholds,
               order2_density_coefficient='positive 9 t^2/(2 kappa)',
               order2_one_body_term='-9 t^2/(4 kappa) sum_edges(n_x+n_y)',
               order2_identity_checks=identities,
               no_claim_that_negative_virtual_shift_implies_binding=True,
               physical_spin_flavour_nucleus_identification=False,
               continuum_or_halo_lifetime_present=False,
               physical_nuclear_row_solved=False, pins=pins)
    (HERE/'NUCLEAR_THRESHOLD_REVIEW.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(thresholds, indent=2))


if __name__ == '__main__':
    main()
