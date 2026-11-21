"""Centre classes of the five-level code, Standard Model fields, fluxes, reactions and closed carriers.

Every quantity is an integer, a rational number or an element of Q(sqrt 5); no floating point enters.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import itertools

from .exact import icosahedron

# ------------------------------------------------------------------ 1. the 32 occupation patterns
LEVELS = ('v1', 'v2', 'v3', 's1', 's2')
CHARGES = (F(-1, 3),) * 3 + (F(1, 2),) * 2          # hypercharge of each level
COLOUR_OF_N3 = {0: '1', 1: '3', 2: '3b', 3: '1'}       # Lambda^n3 C^3
WEAK_OF_N2 = {0: '1', 1: '2', 2: '1'}                  # Lambda^n2 C^2


def _particle_name(n3, n2, s1):
    if n3 == 1:
        return {0: 'd_R', 2: 'u_R'}.get(n2) or ('u_L' if s1 else 'd_L')
    return {0: 'e_R', 2: 'nu_R'}.get(n2) or ('nu_L' if s1 else 'e_L')


MULTIPLET_OF = {'u_L': 'Q', 'd_L': 'Q', 'u_R': 'u_R', 'd_R': 'd_R',
                'nu_L': 'L', 'e_L': 'L', 'e_R': 'e_R', 'nu_R': 'nu_R'}


def patterns():
    rows = []
    for mask in range(32):
        bits = [mask >> i & 1 for i in range(5)]
        n3, n2 = sum(bits[:3]), sum(bits[3:])
        y = sum((c for c, b in zip(CHARGES, bits) if b), F(0))
        t3 = F(bits[3] - bits[4], 2)
        six_y = 6 * y
        if n3 % 2:
            name = _particle_name(n3, n2, bits[3])
            multiplet = MULTIPLET_OF[name]
        else:
            comp = [1 - b for b in bits]
            base = _particle_name(sum(comp[:3]), sum(comp[3:]), comp[3])
            name, multiplet = 'anti-' + base, 'anti-' + MULTIPLET_OF[base]
        rows.append(dict(mask=mask, occupied=[LEVELS[i] for i in range(5) if bits[i]],
                         n3=n3, n2=n2, Y=y, six_Y=six_y, T3=t3, Q=t3 + y,
                         k=int(six_y) % 6, colour=COLOUR_OF_N3[n3], weak=WEAK_OF_N2[n2],
                         name=name, multiplet=multiplet))
    return rows


# ------------------------------------------------------------------ Standard Model fields
# name, SU(3), SU(2), Y, Lorentz type (L, R Weyl; S scalar; V vector), left-handed conjugate name
SU3_TRIALITY = {'1': 0, '3': 1, '3b': 2, '8': 0}
SU2_TWICE_J = {'1': 0, '2': 1, '2b': 1, '3': 2}
CONJ3 = {'1': '1', '3': '3b', '3b': '3', '8': '8'}
CONJ2 = {'1': '1', '2': '2b', '2b': '2', '3': '3'}
CONJ_HAND = {'L': 'R', 'R': 'L', 'S': 'S', 'V': 'V'}
BASE_FIELDS = [
    ('Q', '3', '2', F(1, 6), 'L'), ('u_R', '3', '1', F(2, 3), 'R'), ('d_R', '3', '1', F(-1, 3), 'R'),
    ('L', '1', '2', F(-1, 2), 'L'), ('e_R', '1', '1', F(-1), 'R'), ('nu_R', '1', '1', F(0), 'R'),
    ('H', '1', '2', F(1, 2), 'S'),
    ('G', '8', '1', F(0), 'V'), ('W', '1', '3', F(0), 'V'), ('B', '1', '1', F(0), 'V')]
LEFT_CONJUGATE_NAME = {'anti-u_R': 'u^c', 'anti-d_R': 'd^c', 'anti-e_R': 'e^c', 'anti-nu_R': 'nu^c'}
SELF_CONJUGATE = {'G', 'W', 'B'}


def fields():
    out = []
    for name, c, w, y, h in BASE_FIELDS:
        out.append(dict(name=name, su3=c, su2=w, Y=y, lorentz=h))
        if name not in SELF_CONJUGATE:
            out.append(dict(name='anti-' + name, su3=CONJ3[c], su2=CONJ2[w], Y=-y, lorentz=CONJ_HAND[h]))
    for f in out:
        f['k'] = int(6 * f['Y']) % 6
        f['left_handed_name'] = LEFT_CONJUGATE_NAME.get(f['name'])
    return out


def conj(f):
    return dict(f, name=f['name'][5:] if f['name'].startswith('anti-') else 'anti-' + f['name'],
                su3=CONJ3[f['su3']], su2=CONJ2[f['su2']], Y=-f['Y'], lorentz=CONJ_HAND[f['lorentz']])


# ------------------------------------------------------------------ 2. the Z6 kernel
def kernel_phase_turns(triality, twice_j, y):
    """Phase (in turns) of (omega I3, -I2, exp(2 pi i Y)) on a component: triality/3 + 2j/2 + Y."""
    return F(triality, 3) + F(twice_j, 2) + y


# ------------------------------------------------------------------ exact invariant counting for SU(3) and SU(2)
def _unit(n, a, b):
    m = [[0] * n for _ in range(n)]
    m[a][b] = 1
    return m


def _neg_transpose(m):
    return [[-m[j][i] for j in range(len(m))] for i in range(len(m))]


def _commutator(a, b):
    n = len(a)
    ab = [[sum(a[i][k] * b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    ba = [[sum(b[i][k] * a[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    return [[ab[i][j] - ba[i][j] for j in range(n)] for i in range(n)]


def _adjoint(gens, basis, coords):
    """Matrices of ad(X) for each generator X in the given basis with integer coordinate map."""
    out = []
    for x in gens:
        cols = [coords(_commutator(x, b)) for b in basis]
        out.append([[cols[c][r] for c in range(len(basis))] for r in range(len(basis))])
    return out


def _sl3():
    gens = [_unit(3, 0, 1), _unit(3, 1, 2), _unit(3, 1, 0), _unit(3, 2, 1)]   # E1, E2, F1, F2
    off = [(i, j) for i in range(3) for j in range(3) if i != j]
    h1 = [[1, 0, 0], [0, -1, 0], [0, 0, 0]]
    h2 = [[0, 0, 0], [0, 1, 0], [0, 0, -1]]
    basis = [_unit(3, i, j) for i, j in off] + [h1, h2]

    def coords(m):
        assert m[0][0] + m[1][1] + m[2][2] == 0
        return [m[i][j] for i, j in off] + [m[0][0], -m[2][2]]
    return {'1': [[[0]]] * 4, '3': gens, '3b': [_neg_transpose(g) for g in gens],
            '8': _adjoint(gens, basis, coords)}


def _sl2():
    gens = [_unit(2, 0, 1), _unit(2, 1, 0)]                                   # E, F
    basis = [_unit(2, 0, 1), [[1, 0], [0, -1]], _unit(2, 1, 0)]

    def coords(m):
        return [m[0][1], m[0][0], m[1][0]]
    return {'1': [[[0]]] * 2, '2': gens, '2b': [_neg_transpose(g) for g in gens],
            '3': _adjoint(gens, basis, coords)}


REPS = {'su3': _sl3(), 'su2': _sl2()}


def _rank(rows, ncols):
    """Exact rank over Q of a sparse integer matrix given as dict rows."""
    rows = [dict(r) for r in rows if r]
    rank = 0
    for col in range(ncols):
        pivot = next((r for r in rows if r.get(col)), None)
        if pivot is None:
            continue
        rows.remove(pivot)
        pv = F(pivot[col])
        reduced = []
        for r in rows:
            if r.get(col):
                f = F(r[col]) / pv
                r = dict(r)
                for c, v in pivot.items():
                    r[c] = r.get(c, 0) - f * v
                    if r[c] == 0:
                        del r[c]
            if r:
                reduced.append(r)
        rows = reduced
        rank += 1
    return rank


@lru_cache(maxsize=None)
def invariants(group, reps):
    """Dimension of the invariant subspace of the tensor product of the named representations."""
    mats = [REPS[group][r] for r in reps]
    dims = [len(m[0]) for m in mats]
    total = 1
    for d in dims:
        total *= d
    index = list(itertools.product(*[range(d) for d in dims]))
    pos = {t: i for i, t in enumerate(index)}
    rows = []
    for g in range(len(mats[0])):
        # action of generator g on the tensor product: sum over factors
        cols = {}
        for t in index:
            col = pos[t]
            for f, m in enumerate(mats):
                gm = m[g]
                for r in range(dims[f]):
                    v = gm[r][t[f]]
                    if v:
                        target = pos[t[:f] + (r,) + t[f + 1:]]
                        cols.setdefault(target, {})
                        cols[target][col] = cols[target].get(col, 0) + v
        rows.extend(r for r in cols.values() if any(r.values()))
    return total - _rank(rows, total)


def gauge_singlets(a, b, c):
    """Number of SU(3) x SU(2) singlets in A (x) B (x) C*, and whether hypercharge balances."""
    n3 = invariants('su3', (a['su3'], b['su3'], CONJ3[c['su3']]))
    n2 = invariants('su2', (a['su2'], b['su2'], CONJ2[c['su2']]))
    return n3 * n2, a['Y'] + b['Y'] - c['Y'] == 0


def lorentz_kind(a, b, c):
    """Renormalizable vertex type of the all-incoming triple (A, B, conj C)."""
    kinds = Counter((a['lorentz'], b['lorentz'], CONJ_HAND[c['lorentz']]))
    if kinds['S'] == 1 and (kinds['L'] == 2 or kinds['R'] == 2):
        return 'Yukawa'
    if kinds['V'] == 1 and kinds['L'] == 1 and kinds['R'] == 1:
        return 'fermion gauge current'
    if kinds['V'] == 1 and kinds['S'] == 2:
        return 'scalar gauge current'
    if kinds['V'] == 3:
        return 'gauge self-coupling'
    if kinds['L'] + kinds['R'] == 3:
        return 'three-fermion gauge singlet'
    return 'no renormalizable Lorentz scalar'


# ------------------------------------------------------------------ 5. reactions
REACTION_TYPES = ('1+5->0', '1+1->2', '5+5->4', '1+4->5', '1+2->3', '4+5->3', '1+3->4', '3+5->2',
                  '2+2->4', '2+3->5', '2->1+1', '3->1+2', '1->2+5', '4->2+2', '4->1+3', '5->2+3')
PAIR_CHANNELS = ('1+5->0', '2+4->0', '3+3->0')
# vertex type claimed for each class-level reaction
CLAIMED_VERTEX = {
    '1+5->0': 'fermion gauge current', '2+4->0': 'fermion gauge current', '3+3->0': 'Yukawa',
    '1+2->3': 'Yukawa', '4+5->3': 'Yukawa', '1+3->4': 'Yukawa', '3+5->2': 'Yukawa', '2+3->5': 'Yukawa',
    '3->1+2': 'Yukawa', '4->1+3': 'Yukawa', '5->2+3': 'Yukawa',
    '1+1->2': 'three-fermion gauge singlet', '5+5->4': 'three-fermion gauge singlet',
    '1+4->5': 'three-fermion gauge singlet', '2+2->4': 'three-fermion gauge singlet',
    '2->1+1': 'three-fermion gauge singlet', '1->2+5': 'three-fermion gauge singlet',
    '4->2+2': 'three-fermion gauge singlet'}
READING = {
    '1+5->0': 'annihilation of a left-handed quark doublet with its antiparticle',
    '2+4->0': 'annihilation of a right-handed quark with its antiparticle',
    '3+3->0': 'lepton mass vertex (L with the Higgs into e_R or nu_R); also L Lbar and H Hbar annihilation',
    '1+2->3': 'quark Yukawa vertex: Q with d^c into the Higgs class',
    '4+5->3': 'conjugate of 1+2->3',
    '1+3->4': 'quark mass vertex: Q with the Higgs into u_R or d_R',
    '3+5->2': 'conjugate of 1+3->4',
    '2+3->5': 'conjugate quark mass vertex: u^c or d^c with the Higgs into anti-Q',
    '3->1+2': 'reverse of 1+2->3', '4->1+3': 'reverse of 1+3->4', '5->2+3': 'reverse of 2+3->5',
    '1+1->2': 'two quark doublets into the antitriplet singlet class (diquark)',
    '5+5->4': 'conjugate of 1+1->2', '1+4->5': 'quark doublet with a right-handed quark into the antitriplet doublet class',
    '2+2->4': 'two right-handed antiquarks into the triplet singlet class',
    '2->1+1': 'reverse of 1+1->2', '1->2+5': 'reverse of 5+2->1, conjugate of 1+4->5', '4->2+2': 'reverse of 2+2->4'}


def parse(reaction):
    lhs, rhs = reaction.split('->')
    a = [int(x) for x in lhs.split('+')]
    b = [int(x) for x in rhs.split('+')]
    if len(a) == 1:                    # split c -> a + b is the vertex a + b -> c read backwards
        a, b = b, a
    return a[0], a[1], b[0]


def conjugate_reaction(reaction):
    lhs, rhs = reaction.split('->')
    flip = lambda side: '+'.join(str((-int(x)) % 6) for x in side.split('+'))
    return flip(lhs) + '->' + flip(rhs)


def realizations(reaction, field_list):
    a, b, c = parse(reaction)
    by_k = {k: [f for f in field_list if f['k'] == k] for k in range(6)}
    out = []
    for fa, fb in itertools.product(by_k[a], by_k[b]):
        if a == b and fa['name'] > fb['name']:
            continue
        for fc in by_k[c]:
            n, y_ok = gauge_singlets(fa, fb, fc)
            if n and y_ok:
                colour = sorted((fa['su3'], fb['su3'], CONJ3[fc['su3']]))
                out.append(dict(fields=[fa['name'], fb['name'], fc['name']], singlets=n,
                                hypercharge_sum=str(fa['Y'] + fb['Y'] - fc['Y']),
                                colour_epsilon=colour in (['3', '3', '3'], ['3b', '3b', '3b']),
                                vertex=lorentz_kind(fa, fb, fc)))
    return out


# Standard Model Yukawa operators: (conjugated left field, scalar, right field)
YUKAWA_OPERATORS = {
    'Q^dag H d_R': ('anti-Q', 'H', 'd_R'),
    'Q^dag Htilde u_R': ('anti-Q', 'anti-H', 'u_R'),
    'L^dag H e_R': ('anti-L', 'H', 'e_R'),
    'L^dag Htilde nu_R': ('anti-L', 'anti-H', 'nu_R'),
}


def yukawa_checks(field_list):
    by_name = {f['name']: f for f in field_list}
    one = dict(name='1', su3='1', su2='1', Y=F(0), lorentz='S')
    rows = {}
    for op, names in YUKAWA_OPERATORS.items():
        a, h, b = (by_name[n] for n in names)
        # A (x) H (x) B contains a singlet: A + H -> conj B with C = conj(B)
        n, y_ok = gauge_singlets(a, h, conj(b))
        rows[op] = dict(hypercharge_sum=str(a['Y'] + h['Y'] + b['Y']), gauge_singlets=n,
                        classes=[a['k'], h['k'], b['k']], class_sum_mod_6=(a['k'] + h['k'] + b['k']) % 6,
                        fusion=f"{conj(a)['name']} + {conj(h)['name']} -> {b['name']}: "
                               f"{(-a['k']) % 6}+{(-h['k']) % 6}->{b['k']}")
        rows[op]['neutral'] = (y_ok and n > 0 and a['Y'] + h['Y'] + b['Y'] == 0 and rows[op]['class_sum_mod_6'] == 0)
    return rows


# ------------------------------------------------------------------ 4. Aharonov-Bohm visibility
def aharonov_bohm():
    phase = {f: {k: F(f * k, 6) % 1 for k in range(6)} for f in range(6)}
    seen_by = {f: [k for k in range(6) if phase[f][k] != 0] for f in range(6)}
    return phase, seen_by


# ------------------------------------------------------------------ 6. closed carriers
CARRIER_CONTENTS = ((2, 4), (1, 5), (3, 3), (2, 2, 4, 4), (4, 4, 4), (2, 2, 2), (2, 3, 3, 4), (1, 2, 3), (3, 4, 5))


def carrier_stokes():
    """Face fluxes of any link assignment on the closed icosahedral carrier sum to zero."""
    verts, edges, faces = icosahedron()
    incidence = {e: [] for e in edges}
    for fi, (a, b, c) in enumerate(faces):
        for u, v in ((a, b), (b, c), (c, a)):
            e, s = ((u, v), 1) if u < v else ((v, u), -1)
            incidence[e].append((fi, s))
    column_sums = {e: sum(s for _, s in inc) for e, inc in incidence.items()}
    return dict(ports=len(verts), seams=len(edges), faces=len(faces),
                faces_per_seam=sorted({len(v) for v in incidence.values()}),
                seam_column_sums=sorted(set(column_sums.values())))


def neutral_multisets(max_size):
    out = []
    for size in range(1, max_size + 1):
        for ms in itertools.combinations_with_replacement(range(1, 6), size):
            if sum(ms) % 6 == 0:
                out.append(ms)
    return out
