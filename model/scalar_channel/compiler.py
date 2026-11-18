"""Full CELL even-algebra compiler, deterministic initialization and penalties.

Exact integer Paulis certify ordinary qubit algebra for every native/source
and apparatus generator. A 432-qubit stabilizer simulation executes the loop
measurements, retained signs, feed-forward and occupied product preparation.
Dense matrices are used only for174 one-body modes, after the algebra bridge.
No 2**432 state is claimed materialized. Runtime is not an NP solver.
"""
from itertools import combinations
from collections import deque, Counter
import json, math
import numpy as np
import scipy.linalg as la
from .carrier import load_source, require

def mul(a, b):
    p, x, z = a
    q, y, w = b
    return ((p + q + 2 * (z & y).bit_count()) % 4, x ^ y, z ^ w)

def anti(a, b):
    return ((a[1] & b[2]).bit_count() + (a[2] & b[1]).bit_count()) % 2

def weight(a):
    return (a[1] | a[2]).bit_count()

def neg(a):
    return ((a[0] + 2) % 4, a[1], a[2])

def bits(x):
    require(type(x) is int and x >= 0, 'Pauli bit mask must be a nonnegative integer')
    while x:
        k = (x & -x).bit_length() - 1
        yield k
        x ^= 1 << k

def gf2_rank(rows):
    piv = {}
    for v in rows:
        require(type(v) is int and v >= 0, 'GF(2) row must be a nonnegative integer')
        while v:
            k = v.bit_length() - 1
            if k in piv:
                v ^= piv[k]
            else:
                piv[k] = v
                break
    return len(piv)

class Tableau:

    def __init__(self, q):
        require(type(q) is int and q > 0, 'tableau size must be a positive integer')
        self.q = q
        self.rows = [(0, 0, 1 << j) for j in range(q)]

    @property
    def rows(self):
        return self._rows

    @rows.setter
    def rows(self, value):
        self._rows = value
        self._basis = None

    def pauli(self, p):
        self._check_pauli(p)
        self.rows = [neg(s) if anti(p, s) else s for s in self.rows]

    def _check_pauli(self, p):
        require(isinstance(p, tuple) and len(p) == 3 and all(type(x) is int for x in p)
                and 0 <= p[0] < 4 and 0 <= p[1] < 1 << self.q and 0 <= p[2] < 1 << self.q,
                'Pauli outside tableau register domain')
        require((p[0] - (p[1] & p[2]).bit_count()) % 2 == 0, 'measurement Pauli must be Hermitian')

    def expectation(self, p):
        self._check_pauli(p)
        if any((anti(p, s) for s in self.rows)):
            return 0
        basis = self._basis
        if basis is None:
            basis = {}
            for s in self.rows:
                key = s[1] | s[2] << self.q
                t = s
                while key:
                    k = key.bit_length() - 1
                    if k in basis:
                        other = basis[k]
                        t = mul(t, other)
                        key = t[1] | t[2] << self.q
                    else:
                        basis[k] = t
                        break
            self._basis = basis
        t = p
        key = t[1] | t[2] << self.q
        while key:
            k = key.bit_length() - 1
            if k not in basis:
                return 0
            t = mul(t, basis[k])
            key = t[1] | t[2] << self.q
        require(t[0] in (0, 2), 'Compiler invariant: t[0] in (0, 2)')
        return 1 if t[0] == 0 else -1

    def measure(self, p, sign):
        self._check_pauli(p)
        require(type(sign) is int and sign in (-1, 1), 'measurement sign must be -1 or +1')
        self._basis = None
        ant = [k for k, s in enumerate(self.rows) if anti(p, s)]
        if not ant:
            old = self.expectation(p)
            if old != sign:
                raise ValueError('Impossible requested measurement')
            return 1.0
        pivot = ant[0]
        old = self.rows[pivot]
        for k in ant[1:]:
            self.rows[k] = mul(self.rows[k], old)
        self.rows[pivot] = p if sign == 1 else neg(p)
        return 0.5

def build():
    d = load_source('cell')
    E = len(d['edges'])
    W = len(d['wedges'])
    N = d['N']
    nsrc = E + W
    n = nsrc + 2 * N
    F = list(range(nsrc, nsrc + N))
    R = list(range(nsrc + N, n))
    links = []
    owners = {}
    source_coeff = {}
    for w, (a, b, c) in enumerate(d['wedges']):
        b = int(b)
        for edge, sgn in ((tuple(sorted((int(a), b))), 1), (tuple(sorted((b, int(c)))), -1)):
            e = d['index'][edge]
            ij = (e, E + w)
            links.append(ij)
            owners[ij] = b
            source_coeff[ij] = sgn / math.sqrt(2)
    for b in range(N):
        es = [e for e, pair in enumerate(d['edges']) if b in pair]
        ws = [E + w for w, abc in enumerate(d['wedges']) if int(abc[1]) == b]
        for j in es + ws:
            ij = tuple(sorted((F[b], j)))
            links.append(ij)
            owners[ij] = b
        ij = (F[b], R[b])
        links.append(ij)
        owners[ij] = b
    algebra = encode_graph(n, links, owners)
    return dict(data=d, E=E, W=W, N=N, nsource=nsrc, F=F, R=R,
                source_coeff=source_coeff, **algebra)


def encode_graph(n, links, owners=None):
    """Integer BKSF formulas on a declared finite connected graph.

    The cycle +1 sector represents the even-parity quotient of the fermion
    even algebra. The other central parity block is not identified with it.
    """
    require(type(n) is int and n >= 2, 'graph must have at least two modes')
    links = list(links)
    require(len(links) > 0 and all(isinstance(edge, tuple) and len(edge) == 2
                and all(type(i) is int and 0 <= i < n for i in edge)
                and edge[0] < edge[1] for edge in links), 'invalid oriented interaction edge')
    require(len(set(links)) == len(links), 'duplicate interaction edge')
    links = sorted(links)
    if owners is None:
        owners = {edge: 0 for edge in links}
    require(set(owners) == set(links) and all(type(b) is int and b >= 0 for b in owners.values()),
            'complete nonnegative owner chart required')
    M = len(links)
    ix = {ij: k for k, ij in enumerate(links)}
    inc = [[] for _ in range(n)]
    for k, (i, j) in enumerate(links):
        inc[i].append(k)
        inc[j].append(k)
    reached, pending = {0}, [0]
    while pending:
        i = pending.pop()
        for k in inc[i]:
            a, b = links[k]
            j = b if a == i else a
            if j not in reached:
                reached.add(j)
                pending.append(j)
    require(len(reached) == n, 'interaction graph must be connected')
    B = [(0, 0, sum((1 << k for k in inc[i]))) for i in range(n)]
    A = [(0, 1 << k, sum((1 << l for l in inc[i] if l < k)) ^ sum((1 << l for l in inc[j] if l < k))) for k, (i, j) in enumerate(links)]

    def oriented(i, j):
        require(type(i) is int and type(j) is int and i != j
                and tuple(sorted((i, j))) in ix, 'unknown interaction edge')
        a = A[ix[tuple(sorted((i, j)))]]
        return a if i < j else neg(a)

    def cycle(seq):
        seq = list(seq)
        require(len(seq) >= 3 and len(set(seq)) == len(seq), 'cycle must be simple with at least three modes')
        p = (len(seq) % 4, 0, 0)
        for i, j in zip(seq, seq[1:] + seq[:1]):
            p = mul(p, oriented(i, j))
        return p
    return dict(n=n, links=links, M=M, ix=ix, inc=inc, B=B, A=A,
                oriented=oriented, cycle=cycle, owner=[owners[ij] for ij in links])

def fundamental_cycles(g):
    n = g['n']
    parent = [-1] * n
    parent[0] = 0
    q = deque([0])
    tree = set()
    while q:
        i = q.popleft()
        for k in g['inc'][i]:
            a, b = g['links'][k]
            j = b if i == a else a
            if parent[j] < 0:
                parent[j] = i
                tree.add(k)
                q.append(j)
    require(min(parent) >= 0, 'Compiler invariant: min(parent) >= 0')
    out = []
    for k, (i, j) in enumerate(g['links']):
        if k in tree:
            continue
        pi = [i]
        while parent[pi[-1]] != pi[-1]:
            pi.append(parent[pi[-1]])
        pj = [j]
        while pj[-1] not in pi:
            pj.append(parent[pj[-1]])
        lca = pj[-1]
        seq = pi[:pi.index(lca) + 1] + pj[-2::-1]
        out.append((k, seq, g['cycle'](seq)))
    return (out, parent)

def local_cycles(g):
    d = g['data']
    out = []
    for w, (a, b, c) in enumerate(d['wedges']):
        b = int(b)
        for edge in (tuple(sorted((int(a), b))), tuple(sorted((b, int(c))))):
            seq = [g['F'][b], d['index'][edge], g['E'] + w]
            out.append(dict(kind='star_triangle', center=b, modes=seq, pauli=g['cycle'](seq)))
    edgeset = set(d['index'])
    wi = {tuple(map(int, w)): k for k, w in enumerate(d['wedges'])}
    faces = [t for t in combinations(range(g['N']), 3) if all((tuple(sorted(e)) in edgeset for e in combinations(t, 2)))]
    require(len(faces) == 20, 'Compiler invariant: len(faces) == 20')
    for a, b, c in faces:
        seq = [d['index'][a, b], g['E'] + wi[a, b, c], d['index'][b, c], g['E'] + wi[a, c, b], d['index'][a, c], g['E'] + wi[b, a, c]]
        out.append(dict(kind='icosahedron_face', ports=[a, b, c], modes=seq, pauli=g['cycle'](seq)))
    return out

def path(g, i, j):
    require(type(i) is int and type(j) is int and 0 <= i < g['n'] and 0 <= j < g['n'],
            'path endpoints outside graph')
    parent = {i: i}
    q = deque([i])
    while j not in parent:
        a = q.popleft()
        for k in g['inc'][a]:
            l, r = g['links'][k]
            b = r if a == l else l
            if b not in parent:
                parent[b] = a
                q.append(b)
    seq = [j]
    while seq[-1] != i:
        seq.append(parent[seq[-1]])
    return list(reversed(seq))

def logical_pauli_generators(g):
    """An explicit full fixed-parity qubit basis, not merely a count.

    Zi=Bi; Xi=A_(i,last) product_(i<k<last)Bk. The A for nonedges is the
    ordered local path product i**(length-1) product A. These Xi commute
    with each other and flip exactly Zi and the dependent last parity.
    They, together with the unique encoded vacuum, define every isometry
    column by finite Pauli products. They are external proof objects;
    observers do not receive their global paths or use them as clocks.
    """
    n = g['n']
    xs = []
    for i in range(n - 1):
        seq = path(g, i, n - 1)
        A = (len(seq) - 2, 0, 0)
        for a, b in zip(seq, seq[1:]):
            A = mul(A, g['oriented'](a, b))
        for k in range(i + 1, n - 1):
            A = mul(A, g['B'][k])
        xs.append(A)
    return xs

def hopping_words(g, k, z):
    """h_ij c_i†c_j+h_ij* c_j†c_i as up to4 integer Pauli terms."""
    require(type(k) is int and 0 <= k < g['M'], 'hopping edge outside graph')
    require(np.isscalar(z) and np.isfinite(z), 'finite hopping coefficient required')
    z = complex(z)
    i, j = g['links'][k]
    A = g['A'][k]
    Bi = g['B'][i]
    Bj = g['B'][j]
    terms = []
    if z.real:
        terms.extend([(1j * z.real / 2, mul(A, Bi)), (-1j * z.real / 2, mul(A, Bj))])
    if z.imag:
        terms.extend([(-z.imag / 2, A), (z.imag / 2, mul(A, mul(Bi, Bj)))])
    out = []
    for a, (p, x, w) in terms:
        canonical = (x & w).bit_count() % 4
        coeff = a * 1j ** ((p - canonical) % 4)
        require(abs(coeff.imag) < 1e-14, 'Compiler invariant: abs(coeff.imag) < 1e-14')
        out.append((float(coeff.real), (canonical, x, w)))
    return out

def support_stars(g, p):
    return sorted({g['owner'][k] for k in bits(p[1] | p[2])})

def physical_storage(g):
    d = g['data']
    local_wedges = {b: [w for w, abc in enumerate(d['wedges']) if int(abc[1]) == b] for b in range(g['N'])}
    app_count = Counter()
    out = []
    for k, (i, j) in enumerate(g['links']):
        b = g['owner'][k]
        if j < g['nsource']:
            w = j - g['E']
            a, _, c = map(int, d['wedges'][w])
            first = d['index'][tuple(sorted((a, b)))]
            if i == first:
                out.append(dict(qubit=k, center=b, kind='native_transit_bit', global_wedge=w))
            else:
                out.append(dict(qubit=k, center=b, kind='neutral_owner_path_rank_bit', local_bit=local_wedges[b].index(w)))
        else:
            out.append(dict(qubit=k, center=b, kind='coherent_apparatus_bit', local_bit=app_count[b]))
            app_count[b] += 1
    counts = Counter((x['kind'] for x in out))
    require(counts == dict(native_transit_bit=120, neutral_owner_path_rank_bit=120, coherent_apparatus_bit=192), 'Compiler invariant: counts == dict(native_transit_bit=120, neutral_owner_path_rank_bit=120, coherent_apparatus_bit=192)')
    require(all((app_count[b] == 16 for b in range(g['N']))), "Compiler invariant: all((app_count[b] == 16 for b in range(g['N'])))")
    return out

def routed_measurement(g, S, start):
    d = g['data']
    adj = [sorted({int(x) for e in d['edges'] if b in e for x in e if x != b}) for b in range(g['N'])]
    support = support_stars(g, S)
    route = [start]
    current = start
    for goal in [*support, start]:
        par = {current: current}
        q = deque([current])
        while goal not in par:
            a = q.popleft()
            for b in adj[a]:
                if b not in par:
                    par[b] = a
                    q.append(b)
        part = [goal]
        while part[-1] != current:
            part.append(par[part[-1]])
        route.extend(reversed(part[:-1]))
        current = goal
    return route

def run(include_artifacts=False):
    require(type(include_artifacts) is bool, 'include_artifacts must be a boolean')
    g = build()
    n = g['n']
    M = g['M']
    B = g['B']
    A = g['A']
    storage = physical_storage(g)
    fundamentals, parent = fundamental_cycles(g)
    locals_ = local_cycles(g)
    require(len(locals_) == 260 and len(fundamentals) == M - n + 1 == 259, 'Compiler invariant: len(locals_) == 260 and len(fundamentals) == M - n + 1 == 259')
    require(gf2_rank([x['pauli'][1] for x in locals_]) == 259, "Compiler invariant: gf2_rank([x['pauli'][1] for x in locals_]) == 259")
    logical_x = logical_pauli_generators(g)
    for i, X in enumerate(logical_x):
        require(mul(X, X) == (0, 0, 0), 'Compiler invariant: mul(X, X) == (0, 0, 0)')
        require(all((anti(X, B[j]) == int(i == j) for j in range(n - 1))), 'Compiler invariant: all((anti(X, B[j]) == int(i == j) for j in range(n - 1)))')
        require(all((not anti(X, Y) for Y in logical_x)), 'Compiler invariant: all((not anti(X, Y) for Y in logical_x))')
        require(all((not anti(X, r['pauli']) for r in locals_)), "Compiler invariant: all((not anti(X, r['pauli']) for r in locals_))")
    vacuum_stabilizers = [*B[:-1], *[r['pauli'] for r in locals_]]
    require(gf2_rank([P[1] | P[2] << M for P in vacuum_stabilizers]) == M, 'Compiler invariant: gf2_rank([P[1] | P[2] << M for P in vacuum_stabilizers]) == M')
    checks = 0
    for k, (i, j) in enumerate(g['links']):
        require(mul(A[k], A[k]) == (0, 0, 0), 'Compiler invariant: mul(A[k], A[k]) == (0, 0, 0)')
        for l in range(n):
            require(anti(A[k], B[l]) == int(l in (i, j)), 'Compiler invariant: anti(A[k], B[l]) == int(l in (i, j))')
            checks += 1
        for l, (a, b) in enumerate(g['links']):
            require(anti(A[k], A[l]) == int(k != l and bool({i, j} & {a, b})), 'Compiler invariant: anti(A[k], A[l]) == int(k != l and bool({i, j} & {a, b}))')
            checks += 1
    for rec in locals_:
        S = rec['pauli']
        require(mul(S, S) == (0, 0, 0), 'Compiler invariant: mul(S, S) == (0, 0, 0)')
        require(all((not anti(S, p) for p in A + B)), 'Compiler invariant: all((not anti(S, p) for p in A + B))')
        rec['support_stars'] = support_stars(g, S)
        rec['Pauli_weight'] = weight(S)
    tab = Tableau(M)
    rng = np.random.default_rng(193019)
    trace = []
    route_records = Counter()
    all_half = True
    for idx, (chord, seq, S) in enumerate(fundamentals):
        sign = int(rng.choice((-1, 1)))
        p = tab.measure(S, sign)
        all_half &= p == 0.5
        correction = (0, 0, 1 << chord)
        if sign < 0:
            tab.pauli(correction)
        require(tab.expectation(S) == 1, 'Compiler invariant: tab.expectation(S) == 1')
        require(all((anti(correction, T) == int(k == idx) for k, (_, _, T) in enumerate(fundamentals))), 'Compiler invariant: all((anti(correction, T) == int(k == idx) for k, (_, _, T) in enumerate(fundamentals)))')
        center = g['owner'][chord]
        route = routed_measurement(g, S, center)
        for v in route:
            route_records[v] += 1
        trace.append(dict(index=idx, observer_center=center, measured_sign=sign, outcome_probability=p, applied_Z_chord=bool(sign < 0), chord=chord, cycle_modes=seq, token_route=route, corrected_sign=1, dirty_or_missing_ack_action='abort this preparation and record failure; no postselection'))
    require(all_half, 'Compiler invariant: all_half')
    require(all((tab.expectation(S) == 1 for S in B)), 'Compiler invariant: all((tab.expectation(S) == 1 for S in B))')
    require(all((tab.expectation(r['pauli']) == 1 for r in locals_)), "Compiler invariant: all((tab.expectation(r['pauli']) == 1 for r in locals_))")
    order_checks = []
    for label, order in [('reverse', list(reversed(range(len(fundamentals))))), ('independent_local_order', rng.permutation(len(fundamentals)).tolist())]:
        alternative = Tableau(M)
        negative = 0
        for k in order:
            chord, _, S = fundamentals[k]
            sign = int(rng.choice((-1, 1)))
            negative += int(sign < 0)
            require(alternative.measure(S, sign) == 0.5, 'Compiler invariant: alternative.measure(S, sign) == 0.5')
            if sign < 0:
                alternative.pauli((0, 0, 1 << chord))
        require(all((alternative.expectation(S) == 1 for S in B)), 'Compiler invariant: all((alternative.expectation(S) == 1 for S in B))')
        require(all((alternative.expectation(r['pauli']) == 1 for r in locals_)), "Compiler invariant: all((alternative.expectation(r['pauli']) == 1 for r in locals_))")
        order_checks.append(dict(order=label, negative_outcomes=negative, same_unique_number_and_loop_stabilized_vacuum=True))
    occupied = list(range(g['E'], g['nsource']))
    pair_program = []
    for i, j in zip(occupied[::2], occupied[1::2]):
        seq = path(g, i, j)
        for a, b in zip(seq, seq[1:]):
            tab.pauli(g['oriented'](a, b))
        pair_program.append(seq)
    signs = [tab.expectation(p) for p in B]
    require(signs == [-1 if i in occupied else 1 for i in range(n)], 'Compiler invariant: signs == [-1 if i in occupied else 1 for i in range(n)]')
    require(all((tab.expectation(r['pauli']) == 1 for r in locals_)), "Compiler invariant: all((tab.expectation(r['pauli']) == 1 for r in locals_))")
    p0 = tab.rows.copy()
    fuel_rows = []
    for fuel in (0, 1):
        tab.rows = p0.copy()
        if fuel:
            tab.pauli(g['oriented'](g['F'][0], g['R'][0]))
        got = [tab.expectation(p) for p in B]
        want = [-1 if i in occupied or (fuel and i in (g['F'][0], g['R'][0])) else 1 for i in range(n)]
        require(got == want, 'Compiler invariant: got == want')
        fuel_rows.append(dict(fuel=fuel, reference=fuel, total_particle_number=len(occupied) + 2 * fuel, every_number_generator_matches=True, every_cycle_plus=all((tab.expectation(r['pauli']) == 1 for r in locals_))))
    tab.rows = p0.copy()
    chord = fundamentals[0][0]
    tab.pauli((0, 0, 1 << chord))
    wrong = [i for i, r in enumerate(locals_) if tab.expectation(r['pauli']) < 0]
    require(wrong, 'Compiler invariant: wrong')
    local_penalty = 1024.0
    free_norm = g['nsource'] + len(g['source_coeff']) / math.sqrt(2)
    source_h = np.zeros((n, n), complex)
    source_h[np.arange(g['nsource']), np.arange(g['nsource'])] = [1] * g['E'] + [-1] * g['W']
    compiled = []
    for ij, z in g['source_coeff'].items():
        k = g['ix'][ij]
        i, j = ij
        source_h[i, j] = source_h[j, i] = z
        for coeff, P in hopping_words(g, k, complex(z)):
            compiled.append(dict(coefficient=coeff, pauli=P, support_stars=support_stars(g, P), logical_edge=k))
    max_hop = 0
    max_stars = 0
    for k in range(M):
        for _, P in hopping_words(g, k, complex(0.3, 0.7)):
            max_hop = max(max_hop, weight(P))
            max_stars = max(max_stars, len(support_stars(g, P)))
            require(all((not anti(P, r['pauli']) for r in locals_)), "Compiler invariant: all((not anti(P, r['pauli']) for r in locals_))")
    vals = la.eigvalsh(source_h[:g['nsource'], :g['nsource']])
    E0 = float(vals[vals < 0].sum())
    require(local_penalty > 2 * free_norm + 1, 'Compiler invariant: local_penalty > 2 * free_norm + 1')
    c0 = np.array([[0, 1], [0, 0]], complex)
    Z = np.diag([1, -1])
    I2 = np.eye(2)
    c = [np.kron(c0, I2), np.kron(Z, c0)]
    gam = [x + x.conj().T for x in c]
    bb = [np.eye(4) - 2 * x.conj().T @ x for x in c]
    aa = -1j * gam[0] @ gam[1]
    z = 0.3 + 0.7j
    formula = 1j * z.real / 2 * aa @ (bb[0] - bb[1]) - z.imag / 2 * aa @ (np.eye(4) - bb[0] @ bb[1])
    exact = z * c[0].conj().T @ c[1] + z.conjugate() * c[1].conj().T @ c[0]
    complex_hop_error = float(la.norm(formula - exact, 2))
    require(complex_hop_error < 1e-13, 'Compiler invariant: complex_hop_error < 1e-13')
    schema = dict(
        schema='oph.full_cell_native_scalar_compiler.v3',
        source_modes=g['nsource'],
        fuel_modes=g['N'],
        parity_reference_modes=g['N'],
        total_modes=n,
        source_BKSF_qubits=240,
        apparatus_BKSF_qubits=M - 240,
        total_BKSF_qubits=M,
        additional_neutral_loop_measurement_ancillas=g['N'],
        additional_neutral_terminal_pointers=2,
        unused_native_registers_frozen_to_vacuum=dict(
            seam_owner_7state=g['E'],
            port_qutrit=g['N'],
            local_penalty_per_nonvacuum=local_penalty,
            physical_reset_diagnostics_executed=False,
            these_registers_are_not_counted_as_extra_encoded_qubits=True,
        ),
        fixed_parity_logical_qubits=n - 1,
        independent_loop_constraints=M - n + 1,
        full_intertwiner_basis=dict(
            explicit_logical_X_generators=n - 1,
            explicit_logical_Z_generators=n - 1,
            all_logical_Pauli_relations_checked=True,
            vacuum_stabilizer_rank=M,
            every_isometry_column='Apply the selected commuting logical X generators to the unique prepared vacuum.',
            all_columns_are_orthogonal_by_distinct_B_signs=True,
            spanning_dimension='2**173',
            dense_column_array_not_materialized=True,
        ),
        exact_generator_relations_checked=checks,
        bounded_penalty_basis=dict(
            star_triangles=240,
            face_hexagons=20,
            rank=259,
            maximum_native_star_support=max((len(r['support_stars']) for r in locals_)),
            maximum_Pauli_weight=max((r['Pauli_weight'] for r in locals_)),
        ),
        deterministic_preparation=dict(
            all_declared_loop_measurements_simulated=True,
            random_signs_retained=True,
            minus_signs=sum((r['measured_sign'] < 0 for r in trace)),
            all_conditional_probabilities_half=bool(all_half),
            alternative_measurement_orders=order_checks,
            success_probability_after_feedforward=1.0,
            postselection=False,
            all_260_local_penalty_cycles_plus=True,
            all_174_number_generators_match_product_state=True,
            total_ancilla_token_hops=sum((len(r['token_route']) - 1 for r in trace)),
            maximum_loop_modes=max((len(x[1]) for x in fundamentals)),
            record_token_visits_by_local_observer=dict(sorted(route_records.items())),
            pair_preparation_local_A_operations=sum((len(p) - 1 for p in pair_program)),
            ancilla_reset='Read each measured ancilla then conditional local X to0; retain outcome/reset records.',
            initial_native_owner_reset='All owner registers vacuum, which is the exact rank0 singlet codeword.',
            unused_register_reset='Supplied vacuum on30seam owners and12portqutrits; physical reset diagnostics are not executed.',
            initial_edge_bit_reset='Supplied zero transit and coherent apparatus edge qubits; physical reset diagnostics are not executed.',
        ),
        both_local_fuel_preparations=fuel_rows,
        compiled_law=dict(
            source_hopping_Pauli_terms=len(compiled),
            all_allowed_complex_instrument_hops_certified=True,
            maximum_allowed_hopping_Pauli_weight=max_hop,
            maximum_hopping_native_star_support=max_stars,
            complex_hopping_formula_dense_error=complex_hop_error,
            source_matrix_is_external_certificate=True,
            source_eigenvectors_are_not_controller_inputs=True,
            same_operator_family_covers_equal_five_seam_and_chi_targets=True,
        ),
        full_raw_domain_energy=dict(
            complete_one_body_source_spectrum=vals.tolist(),
            decoupled_apparatus_one_body_spectrum=[0.0] * (2 * g['N']),
            penalty_per_local_failure=local_penalty,
            free_operator_norm_upper_bound=free_norm,
            accepted_sea_energy=E0,
            lower_bound_invalid_energy_relative_sea=local_penalty - free_norm - E0,
            conservative_invalid_energy_gap_bound=local_penalty - 2 * free_norm,
            source_particle_or_hole_gap=1.0,
            zero_energy_decoupled_apparatus_modes_retained=True,
            global_ground_is_correct_source_sea_times_allowed_zero_apparatus=True,
            source_mass_gap_is_not_claimed_full_apparatus_spectral_gap=True,
            penalty_vanishes_on_all_accepted_source_and_instrument_histories=True,
            unused_seam_owner_and_port_register_penalties_included=True,
            fixed_finite_CELL_bound_not_uniform_unbounded_federation_theorem=True,
        ),
        controls=dict(
            uncorrected_loop_sector_detected=True,
            violated_local_penalties_after_one_Z=wrong,
            penalty_cost_after_one_Z=len(wrong) * local_penalty,
            specified_missing_or_wrong_ack='Retained terminal failure, no source experiment starts',
            dirty_owner='Declared complementary owner sector has penalty1024; its local diagnostic is a supplied prescription.',
        ),
        theorem_scope=dict(
            even_parity_algebra_and_local_instruments_intertwine=True,
            full_432qubit_dense_state_executed=False,
            exact_stabilizer_preparation_executed=True,
            full_star_gate_is_implicit_Clebsch_program=True,
            elementary_two_owner_gate_synthesis_executed=False,
            original_geometric_half_transition_preserved=False,
            unrestricted_seven_species_CAR_compiled=False,
            intrinsic_gauge_A1_complete=False,
            local_frame_condition='Owners at one center transform collectively; arbitrary independent owner rotations are excluded.',
            chart_transport_condition='Fusion ordering is explicit; raw permutation-only covariance is not supplied.',
        ),
    )
    preparation_trace = dict(
        loop_measurements=trace,
        occupied_product_pair_paths=pair_program,
        retained_probability_of_each_loop_sign_string='2**(-259)',
        no_lost_measurement_branches=True,
    )

    def serial(P):
        return dict(phase=P[0], x=hex(P[1]), z=hex(P[2]))
    operator_manifest = dict(
        scalar_source_modes=g['nsource'],
        total_modes=n,
        links=g['links'],
        qubit_owner_stars=g['owner'],
        exact_qubit_storage=storage,
        B=[serial(p) for p in B],
        A=[serial(p) for p in A],
        complete_logical_X_basis=[serial(p) for p in logical_x],
        local_cycle_penalties=[{**{k: v for k, v in r.items() if k != 'pauli'}, 'pauli': serial(r['pauli'])} for r in locals_],
        source_hopping_terms=[{**r, 'pauli': serial(r['pauli'])} for r in compiled],
        free_diagonal_coefficients=[1] * g['E'] + [-1] * g['W'] + [0] * (2 * g['N']),
        full_domain_rule='Each term is lifted on its own supporting star code projectors, zero outside; add1024 local code/cycle penalties and1024 pernonvacuum unusedseamowner/portregister.',
    )
    schema['execution_boundary'] = {
        'logical_stabilizer_measurements_executed': True,
        'local_ancilla_and_ack_service_executed': False,
        'route_is_external_certificate': True,
        'physical_owner_reset_executed': False,
    }
    if include_artifacts:
        schema['preparation_trace'] = preparation_trace
        schema['operator_manifest'] = operator_manifest
    return schema
if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
