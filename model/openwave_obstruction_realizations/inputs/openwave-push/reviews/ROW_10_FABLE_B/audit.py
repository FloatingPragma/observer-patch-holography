"""ROW_10_FABLE_B audit: the reviewer's own enumeration of the residual-symmetry mixing patterns of the carrier's
rotation group, its angle menus, a random-unitary null and an executed vacuum propagation.

No candidate code is imported.  Read as data: the carrier's twenty oriented faces (theory repository text, compared),
the candidate's patterns.json and angles.json (for the entry-by-entry comparison only).

Own construction: the 120 automorphisms of the seam graph by backtracking (the reviewer's code of the row-25
addendum), an own embedding of the ports on the vertices (0, +-1, +-phi) and cyclic, the 60 rotations of determinant
+1, the conjugate triplet by sqrt5 -> -sqrt5 on the identified entries, all subgroups by closure of pairs of
elements, lines by character projectors, squared moduli M = Tr(P Q), patterns up to row and column permutations.

    python audit.py      writes RESULT.json; stdout is RUN.log
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import hashlib
import json
import math
import re
import subprocess
import time
from itertools import permutations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
B = '/Users/muellerberndt/Projects/oph-meta/'
CAND = B + 'openwave-push/rowfull/ROW_10_symmetry/'
FACE_FILE = B + 'reverse-engineering-reality/evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py'
T0 = time.process_time()
S5 = math.sqrt(5.0)
PHI = (1 + S5) / 2
PERMS = list(permutations(range(3)))


def log(*a):
    print('[%6.1fs]' % (time.process_time() - T0), *a, flush=True)


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def date():
    return subprocess.run(['date'], capture_output=True, text=True).stdout.strip()


# ---- the carrier and its group ---------------------------------------------------------------------------------
FACES = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
         (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9), (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
_t = open(FACE_FILE).read()
_b = _t[_t.index('ORIENTED_BASE_FACES = ('):]
_b = _b[:_b.index('\n)\n')]
assert [tuple(int(x) for x in m) for m in re.findall(r'\((\d+), (\d+), (\d+)\)', _b)] == FACES
ADJ = {i: set() for i in range(12)}
for a, b, c in FACES:
    for x, y in ((a, b), (b, c), (c, a)):
        ADJ[x].add(y); ADJ[y].add(x)
assert all(len(ADJ[i]) == 5 for i in range(12))


def isomorphisms(adj_a, adj_b, first_only=False):
    order, seen = [0], {0}
    for v in order:
        for u in sorted(adj_a[v]):
            if u not in seen:
                seen.add(u); order.append(u)
    res = []

    def rec(k, img, used):
        if first_only and res:
            return
        if k == 12:
            res.append(tuple(img[i] for i in range(12))); return
        v = order[k]
        for w in range(12):
            if w not in used and all((img[u] in adj_b[w]) == (u in adj_a[v]) for u in order[:k]):
                img[v] = w; used.add(w); rec(k + 1, img, used); used.discard(w); del img[v]
    rec(0, {}, set())
    return res


AUT = isomorphisms(ADJ, ADJ)
VERT = []
for s1 in (1, -1):
    for s2 in (1, -1):
        VERT += [(0.0, s1, s2 * PHI), (s1, s2 * PHI, 0.0), (s2 * PHI, 0.0, s1)]
VERT = np.array(VERT, dtype=float)
ADJ_STD = {i: {j for j in range(12) if j != i and abs(((VERT[i] - VERT[j]) ** 2).sum() - 4) < 1e-9} for i in range(12)}
EMB = isomorphisms(ADJ, ADJ_STD, first_only=True)[0]
VEC = VERT[list(EMB)]                                      # port i sits at VEC[i]
NORM = 4 * (1 + PHI * PHI)


def rot(g):
    return sum(np.outer(VEC[g[i]], VEC[i]) for i in range(12)) / NORM


ALLM = [rot(g) for g in AUT]
PROPER = [i for i, R in enumerate(ALLM) if np.linalg.det(R) > 0]
G = [AUT[i] for i in PROPER]                               # 60 permutations
R3 = [ALLM[i] for i in PROPER]
IDX = {g: i for i, g in enumerate(G)}
MUL = [[IDX[tuple(g[h[k]] for k in range(12))] for h in G] for g in G]
E_ID = IDX[tuple(range(12))]


def ident(x):
    for q in range(-2, 3):
        p = round(4 * x - q * S5)
        if abs((p + q * S5) / 4 - x) < 1e-11 and abs(p) <= 4:
            return p, q
    raise ArithmeticError('entry outside (p + q sqrt5)/4: %r' % x)


R3P = []
for R in R3:
    Rp = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            p, q = ident(R[i, j])
            Rp[i, j] = (p - q * S5) / 4
    R3P.append(Rp)


def element_order(i):
    n, x = 1, i
    while x != E_ID:
        x = MUL[x][i]; n += 1
    return n


ORD = [element_order(i) for i in range(60)]


def closure(gens):
    S = {E_ID} | set(gens)
    frontier = list(S)
    while frontier:
        x = frontier.pop()
        for y in list(S):
            for z in (MUL[x][y], MUL[y][x]):
                if z not in S:
                    S.add(z); frontier.append(z)
    return frozenset(S)


SUBGROUPS = sorted({closure([i, j]) for i in range(60) for j in range(i, 60)}, key=lambda s: (len(s), sorted(s)))


def is_abelian(H):
    return all(MUL[a][b] == MUL[b][a] for a in H for b in H)


def projectors(H, rep):
    """Three rank-one projectors of a basis-fixing subgroup, or None.  Cyclic: by the character index k (k = 0 is
    the axis); Klein: the three half-turn axes."""
    comm = sum(abs(np.trace(rep[h])) ** 2 for h in H) / len(H)
    if not is_abelian(H) or abs(comm - 3) > 1e-9:
        return None
    gens = [h for h in H if ORD[h] == len(H)]
    if gens:
        a, n = min(gens), len(H)
        pw = [np.eye(3)]
        for _ in range(n - 1):
            pw.append(pw[-1] @ rep[a])
        out = []
        for k in range(n):
            P = sum(np.exp(-2j * np.pi * j * k / n) * pw[j] for j in range(n)) / n
            if abs(np.trace(P) - 1) < 1e-9:
                out.append(P)
        assert len(out) == 3
        return out
    out = [(np.eye(3) + rep[h]) / 2 + 0j for h in sorted(H) if h != E_ID]
    assert len(out) == 3 and all(abs(np.trace(P) - 1) < 1e-9 for P in out)
    return out


def unit(P):
    j = int(np.argmax(np.linalg.norm(P, axis=0)))
    return P[:, j] / np.linalg.norm(P[:, j])


def moduli(PA, PB):
    return np.array([[np.trace(p @ q).real for q in PB] for p in PA])


def canon(M):
    best = None
    for r in PERMS:
        for c in PERMS:
            t = tuple((np.round(M[np.ix_(r, c)], 9) + 0.0).ravel())
            if best is None or t < best:
                best = t
    return best


def find(classes, key, tol=1e-8):
    for k in classes:
        if max(abs(x - y) for x, y in zip(k, key)) < tol:
            return k
    return None


TYPE = {3: 'C3', 5: 'C5', 4: 'V4'}
REF = {'s12': (0.308, 0.012, 0.011), 's13': (0.02215, 0.00056, 0.00058), 's23': (0.470, 0.017, 0.013)}
RANGE3 = {'s12': (0.275, 0.345), 's13': (0.0203, 0.02388), 's23': (0.435, 0.585)}


def triple(M):
    s13 = M[0, 2]
    if s13 > 1 - 1e-12:
        return None
    return (M[0, 1] / (1 - s13), s13, M[1, 2] / (1 - s13))


def pulls(t):
    out = []
    for x, k in zip(t, ('s12', 's13', 's23')):
        c, up, dn = REF[k]
        out.append((x - c) / (up if x >= c else dn))
    return out


def chi2(t):
    return float(sum(p * p for p in pulls(t)))


def labelled(classes_reps):
    """All labelled moduli matrices (rows e, mu, tau; columns 1, 2, 3) of a set of patterns, and their triples."""
    mats, trip, undefined = [], [], 0
    for M in classes_reps:
        for r in PERMS:
            for c in PERMS:
                L = M[np.ix_(r, c)]
                key = tuple((np.round(L, 9) + 0.0).ravel())
                if find(mats, key) is None:
                    mats.append(key)
                    t = triple(L)
                    if t is None:
                        undefined += 1
                    else:
                        t = tuple(float(x) for x in t)
                        if find(trip, t) is None:
                            trip.append(t)
    return mats, trip, undefined


def main():
    pre = json.load(open(os.path.join(HERE, 'PREREG_AUDIT.json')))
    EX, TOL = pre['expected_if_the_candidate_is_exact'], pre['tolerances']
    out = {'schema': 'fable_b.row10.review.audit.v1', 'run_started': date(), 'audit_sha256': sha(os.path.join(HERE, 'audit.py')),
           'prereg_sha256': sha(os.path.join(HERE, 'PREREG_AUDIT.json'))}
    log('start', out['run_started'], out['audit_sha256'][:16], 'prereg', out['prereg_sha256'][:16])

    # ---- A. group --------------------------------------------------------------------------------------------
    hom = max(float(np.abs(R3[MUL[i][j]] - R3[i] @ R3[j]).max()) for i in range(60) for j in range(60))
    homp = max(float(np.abs(R3P[MUL[i][j]] - R3P[i] @ R3P[j]).max()) for i in range(60) for j in range(60))
    orth = max(float(np.abs(R @ R.T - np.eye(3)).max()) for R in R3 + R3P)
    by_order = {}
    for H in SUBGROUPS:
        by_order[len(H)] = by_order.get(len(H), 0) + 1
    el_orders = {}
    for o in ORD:
        el_orders[o] = el_orders.get(o, 0) + 1
    tr5 = sorted({round(float(np.trace(R3[i])), 9) for i in range(60) if ORD[i] == 5})
    tr5p = sorted({round(float(np.trace(R3P[i])), 9) for i in range(60) if ORD[i] == 5})
    same_class = all(abs(np.trace(R3[i]) + np.trace(R3P[i]) - 1) < 1e-9 for i in range(60) if ORD[i] == 5)
    comm_dims = {}
    for H in SUBGROUPS:
        comm_dims.setdefault(len(H), set()).add(round(sum(abs(np.trace(R3[h])) ** 2 for h in H) / len(H), 9))
    A = {'automorphisms_of_the_seam_graph': len(AUT), 'rotations_determinant_plus_1': len(G), 'homomorphism_defect_3': hom, 'homomorphism_defect_3_conjugate': homp,
         'orthogonality_defect': orth, 'element_orders': {str(k): v for k, v in sorted(el_orders.items())},
         'subgroups': len(SUBGROUPS), 'subgroups_by_order': {str(k): v for k, v in sorted(by_order.items())},
         'traces_of_fifth_turns_on_3': tr5, 'traces_of_fifth_turns_on_the_conjugate': tr5p, 'conjugate_exchanges_the_two_classes_of_fifth_turns': bool(same_class),
         'commutant_dimension_by_order_on_3': {str(k): sorted(v) for k, v in sorted(comm_dims.items())}}
    out['A_group'] = A
    log('A', json.dumps(A))

    # ---- B. patterns -----------------------------------------------------------------------------------------
    cand = json.load(open(CAND + 'patterns.json'))
    res_rep = {}
    for name, rep in (('3', R3), ("3'", R3P)):
        fix = [(H, TYPE[len(H)], projectors(H, rep)) for H in SUBGROUPS if projectors(H, rep) is not None]
        classes = {}
        for ia, (HA, tA, PA) in enumerate(fix):
            for ib, (HB, tB, PB) in enumerate(fix):
                M = moduli(PA, PB)
                key = canon(M)
                k = find(classes, key)
                if k is None:
                    k = key
                    classes[k] = {'pairs': {}, 'rep_pair': (ia, ib), 'gen_order': len(closure(list(HA) + list(HB)))}
                tp = tA + ',' + tB
                classes[k]['pairs'][tp] = classes[k]['pairs'].get(tp, 0) + 1
        # order-2 columns
        halves = [H for H in SUBGROUPS if len(H) == 2]
        cols = {}
        for HA, tA, PA in fix:
            for H2 in halves:
                h = [x for x in H2 if x != E_ID][0]
                axis = (np.eye(3) + rep[h]) / 2
                col = tuple(sorted(round(float(np.trace(p @ axis).real), 9) + 0.0 for p in PA))
                k = find(cols, col)
                if k is None:
                    k = col; cols[k] = {}
                cols[k][tA] = cols[k].get(tA, 0) + 1
        res_rep[name] = {'fix': fix, 'classes': classes, 'cols': cols}
    fix3, classes3 = res_rep['3']['fix'], res_rep['3']['classes']
    counts = {t: sum(1 for f in fix3 if f[1] == t) for t in ('C3', 'C5', 'V4')}
    # comparison with the candidate, entry by entry
    comp, used = [], set()
    for q in cand['patterns']:
        key = canon(np.array(q['float'], dtype=float))
        row = {'id': q['id']}
        for name in ('3', "3'"):
            k = find(res_rep[name]['classes'], key)
            if k is None:
                row['found_on_' + name] = False
                continue
            used.add((name, k))
            mine = res_rep[name]['classes'][k]['pairs']
            theirs = {tp.split('|')[1]: n for tp, n in q['type_pairs'].items() if tp.split('|')[0] == name}
            row['found_on_' + name] = True
            row['largest_entry_difference_' + name] = float(max(abs(x - y) for x, y in zip(k, key)))
            row['pair_counts_equal_' + name] = mine == theirs
            row['pairs_' + name] = mine
        # exact strings evaluated
        s = math.sqrt(30 + 6 * S5) / 2
        ex = np.array([[eval(e.replace('sqrt5', 'S5'), {'__builtins__': {}}, {'S5': S5, 's': s}) for e in r] for r in q['exact']], dtype=float)
        row['exact_strings_against_own_entries'] = float(max(abs(x - y) for x, y in zip(canon(ex), find(classes3, key) or key)))
        # transpose and conjugate identifications
        kT = find(classes3, canon(np.array(q['float'], dtype=float).T))
        row['transpose_is_candidate'] = q['transpose_is']; row['conjugate_is_candidate'] = q['conjugate_sqrt5_is']
        row['_kT'] = kT; row['_k'] = find(classes3, key)
        comp.append(row)
    idmap = {r['_k']: r['id'] for r in comp if r['_k'] is not None}
    for r in comp:
        r['transpose_is_own'] = idmap.get(r.pop('_kT'))
        r.pop('_k')
        r['transpose_agrees'] = r['transpose_is_own'] == r['transpose_is_candidate']
    # conjugation: the class of the same pair of subgroups on the conjugate triplet
    conj_map = {}
    fixp = res_rep["3'"]['fix']
    for k, v in classes3.items():
        ia, ib = v['rep_pair']
        Mp = moduli(fixp[ia][2], fixp[ib][2])
        kk = find(classes3, canon(Mp))
        conj_map[idmap.get(k)] = idmap.get(kk)
    for r in comp:
        r['conjugate_is_own'] = conj_map.get(r['id'])
        r['conjugate_agrees'] = r['conjugate_is_own'] == r['conjugate_is_candidate']
    missing_in_candidate = [list(k) for name in ('3', "3'") for k in res_rep[name]['classes'] if (name, k) not in used]
    # J of each pattern
    J = {}
    for k, v in classes3.items():
        ia, ib = v['rep_pair']
        ea, eb = [unit(P) for P in fix3[ia][2]], [unit(P) for P in fix3[ib][2]]
        U = np.array([[np.vdot(x, y) for y in eb] for x in ea])
        J[idmap.get(k)] = float(abs((U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0])).imag))
        v['U'] = U
    # columns
    cand_cols = {tuple(sorted(c['float'])): c for c in cand['c2_columns']}
    col_cmp = []
    for k, v in res_rep['3']['cols'].items():
        kk = find(cand_cols, k)
        theirs = {tp.split('G_e ')[1]: n for tp, n in cand_cols[kk]['ordered_pairs'].items() if tp.startswith('3|')} if kk else None
        col_cmp.append({'column': list(k), 'pairs': v, 'candidate_id': cand_cols[kk]['id'] if kk else None, 'pair_counts_equal': (theirs == v) if kk else False,
                        'same_on_the_conjugate_triplet': find(res_rep["3'"]['cols'], k) is not None and res_rep["3'"]['cols'][find(res_rep["3'"]['cols'], k)] == v})
    menus = {}
    for mname, ge in (('H', ('C3',)), ('X', ('C3', 'C5')), ('A', ('C3', 'C5', 'V4'))):
        ks = [k for k, v in classes3.items() if any(tp.split(',')[0] in ge for tp in v['pairs'])]
        npairs = sum(n for k in ks for tp, n in classes3[k]['pairs'].items() if tp.split(',')[0] in ge)
        ksp = [k for k, v in res_rep["3'"]['classes'].items() if any(tp.split(',')[0] in ge for tp in v['pairs'])]
        menus[mname] = {'ordered_pairs': npairs, 'patterns_on_3': len(ks), 'patterns_on_the_conjugate': len(ksp),
                        'same_set_on_both': all(find(ks, k) is not None for k in ksp) and len(ks) == len(ksp), 'ids': sorted(idmap.get(k) or '?' for k in ks), '_keys': ks}
    Bres = {'basis_fixing_subgroups': counts, 'ordered_pairs_per_triplet': len(fix3) ** 2, 'patterns_on_3': len(classes3), 'patterns_on_the_conjugate': len(res_rep["3'"]['classes']),
            'comparison_with_the_candidate': comp, 'own_patterns_absent_from_the_candidate': missing_in_candidate,
            'J_abs_by_pattern': J, 'candidate_J_abs': {q['id']: q['J_abs'] for q in cand['patterns']},
            'order_of_the_group_generated_by_the_pair': {idmap.get(k): v['gen_order'] for k, v in classes3.items()},
            'fixed_columns': col_cmp, 'fixed_columns_with_G_e_of_order_3': sum(1 for c in col_cmp if 'C3' in c['pairs']),
            'menus': {m: {k: v for k, v in d.items() if k != '_keys'} for m, d in menus.items()}}
    out['B_patterns'] = Bres
    log('B', json.dumps({k: Bres[k] for k in ('basis_fixing_subgroups', 'ordered_pairs_per_triplet', 'patterns_on_3', 'patterns_on_the_conjugate', 'own_patterns_absent_from_the_candidate', 'J_abs_by_pattern', 'order_of_the_group_generated_by_the_pair', 'menus', 'fixed_columns_with_G_e_of_order_3')}))
    for r in comp:
        log('B compare', json.dumps(r))
    for c in col_cmp:
        log('B column', json.dumps(c))

    # ---- C. angles -------------------------------------------------------------------------------------------
    ang = {}
    canda = json.load(open(CAND + 'angles.json'))
    for mname in ('H', 'X', 'A'):
        reps = [np.array(k).reshape(3, 3) for k in menus[mname]['_keys']]
        mats, trip, undef = labelled(reps)
        c2 = sorted((chi2(t), t) for t in trip)
        within3 = [t for t in trip if all(abs(p) <= 3 for p in pulls(t))]
        inrange = [t for t in trip if all(RANGE3[k][0] <= x <= RANGE3[k][1] for x, k in zip(t, ('s12', 's13', 's23')))]
        in30 = [t for t in trip if all(abs(x / REF[k][0] - 1) <= 0.3 for x, k in zip(t, ('s12', 's13', 's23')))]
        ang[mname] = {'labelled_matrices': len(mats), 'with_undefined_angles': undef, 'angle_triples': len(trip), 'bits': math.log2(len(trip)),
                      'within_3_sigma_on_all_three': len(within3), 'inside_the_three_3_sigma_ranges': len(inrange), 'within_30_percent_on_each': len(in30),
                      'nearest_chi2': c2[0][0], 'nearest_sin2': list(c2[0][1]), 'nearest_pulls': pulls(c2[0][1]), '_trip': trip}
    tH = ang['H']['_trip']
    cand_tr = [(t['sin2']['s12'], t['sin2']['s13'], t['sin2']['s23']) for t in canda['triples_menu_H']]
    ang['H']['own_triples_absent_from_the_candidate'] = sum(1 for t in tH if find(cand_tr, t, 1e-7) is None)
    ang['H']['candidate_triples_absent_here'] = sum(1 for t in cand_tr if find(tH, t, 1e-7) is None)
    per_angle = {}
    for i, k in enumerate(('s12', 's13', 's23')):
        vals = []                                          # distinct values by clustering within 1e-7 (the first run
        for x in sorted(t[i] for t in tH):                 # compared at fixed 9-digit rounding and split equal values)
            if not vals or x - vals[-1] > 1e-7:
                vals.append(float(x))
        near = min(vals, key=lambda x: abs(x - REF[k][0]))
        per_angle[k] = {'values': len(vals), 'nearest': near, 'pull': pulls([near if j == i else REF[kk][0] for j, kk in enumerate(('s12', 's13', 's23'))])[i],
                        'within_3_sigma': sum(1 for x in vals if abs(pulls([x if j == i else REF[kk][0] for j, kk in enumerate(('s12', 's13', 's23'))])[i]) <= 3)}
    ang['per_angle_menu_H'] = per_angle
    all_entries = sorted({round(float(x), 9) for k in classes3 for x in k} | {round(float(x), 9) for c in col_cmp for x in c['column']})
    ang['distinct_entries_of_all_patterns_and_columns'] = len(all_entries)
    ang['entries_between_0.0203_and_0.0239'] = [x for x in all_entries if 0.0203 <= x <= 0.0239]
    # menu R: G_e of order 3 with tau on the axis (k = 0) and (e, mu) on the pair in both orders; columns free
    tripR = []
    for HA, tA, PA in fix3:
        if tA != 'C3':
            continue
        for HB, tB, PB in fix3:
            M = moduli(PA, PB)
            for rows in ((1, 2, 0), (2, 1, 0)):
                for c in PERMS:
                    t = triple(M[np.ix_(rows, c)])
                    if t is not None and find(tripR, tuple(float(x) for x in t)) is None:
                        tripR.append(tuple(float(x) for x in t))
    c2R = sorted((chi2(t), t) for t in tripR)
    ang['R'] = {'angle_triples': len(tripR), 'bits': math.log2(len(tripR)), 'nearest_chi2': c2R[0][0], 'nearest_sin2': list(c2R[0][1]),
                'within_3_sigma_on_all_three': sum(1 for t in tripR if all(abs(p) <= 3 for p in pulls(t)))}
    out['C_angles'] = {m: {k: v for k, v in d.items() if k != '_trip'} if isinstance(d, dict) else d for m, d in ang.items()}
    log('C', json.dumps(out['C_angles']))

    # ---- D. real symmetric invariant operators (the Klein statement) ----------------------------------------
    rng = np.random.default_rng(77)
    gaps = {'C3': [], 'C5': [], 'V4': []}
    for H, tH_, PA in fix3:
        Smat = rng.normal(size=(3, 3)); Smat = Smat + Smat.T
        Sbar = sum(R3[h] @ Smat @ R3[h].T for h in H) / len(H)
        ev = np.linalg.eigvalsh(Sbar)
        gaps[tH_].append(float(np.min(np.diff(ev))))
    Hm = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); Hm = Hm + Hm.conj().T
    c3 = [f for f in fix3 if f[1] == 'C3'][0]
    Hbar = sum(R3[h] @ Hm @ R3[h].T for h in c3[0]) / 3
    out['D_real_symmetric_invariant_operators'] = {'largest_smallest_gap_cyclic': max(gaps['C3'] + gaps['C5']), 'smallest_gap_klein': min(gaps['V4']),
                                                   'complex_hermitian_invariant_operator_of_a_third_turn_smallest_gap': float(np.min(np.diff(np.linalg.eigvalsh(Hbar))))}
    log('D', json.dumps(out['D_real_symmetric_invariant_operators']))

    # ---- E. random-unitary null -------------------------------------------------------------------------------
    hn = pre['haar_null']
    rng = np.random.default_rng(hn['seed'])
    thr = ang['H']['nearest_chi2']
    below, s13sum, n_qr = 0, 0.0, hn['draws_qr']
    for _ in range(4):
        Zm = rng.normal(size=(n_qr // 4, 3, 3)) + 1j * rng.normal(size=(n_qr // 4, 3, 3))
        Q, Rr = np.linalg.qr(Zm)
        d = np.diagonal(Rr, axis1=1, axis2=2)
        Q = Q * (d / np.abs(d))[:, None, :]
        M2 = np.abs(Q) ** 2
        s13 = M2[:, 0, 2]
        t12, t23 = M2[:, 0, 1] / (1 - s13), M2[:, 1, 2] / (1 - s13)
        c = np.zeros(len(s13))
        for x, k in ((t12, 's12'), (s13, 's13'), (t23, 's23')):
            cc, up, dn = REF[k]
            c += ((x - cc) / np.where(x >= cc, up, dn)) ** 2
        below += int((c <= thr).sum()); s13sum += float(s13.sum())
    p_qr = below / n_qr
    nf = hn['draws_factorised']
    hits = 0
    for _ in range(4):
        u = rng.random((3, nf // 4))
        t12, t23, s13 = u[0], u[1], 1 - np.sqrt(u[2])
        c = np.zeros(nf // 4)
        for x, k in ((t12, 's12'), (s13, 's13'), (t23, 's23')):
            cc, up, dn = REF[k]
            c += ((x - cc) / np.where(x >= cc, up, dn)) ** 2
        hits += int((c <= thr).sum())
    p_f = hits / nf
    N = ang['H']['angle_triples']
    out['E_random_unitary_null'] = {'threshold_chi2': thr, 'draws_qr': n_qr, 'share_qr': p_qr, 'mean_sin2_theta13_qr': s13sum / n_qr,
                                    'draws_factorised': nf, 'share_factorised': p_f, 'standard_error_factorised': math.sqrt(p_f * (1 - p_f) / nf),
                                    'N': N, 'N_times_p': N * p_f, 'probability_that_one_of_N_random_draws_is_at_least_as_near': 1 - (1 - p_f) ** N,
                                    'candidate_share': canda['null']['share_chi2_at_or_below_nearest'], 'candidate_N_times_p': canda['null']['N_times_p_chi2']}
    log('E', json.dumps(out['E_random_unitary_null']))

    # ---- F. executed vacuum propagation on the nearest pattern -----------------------------------------------
    kP04 = [k for k, v in idmap.items() if v == 'P04'][0]
    U0 = classes3[kP04]['U']
    best = None
    for r in PERMS:
        for cidx in PERMS:
            Ul = U0[np.ix_(r, cidx)]
            t = triple(np.abs(Ul) ** 2)
            if t is not None and (best is None or chi2(t) < best[0]):
                best = (chi2(t), Ul, t)
    U = best[1]
    m_rep = np.array([0.0020642, 0.0636598, 0.9342760])
    r_rep = float((m_rep[1] ** 2 - m_rep[0] ** 2) / (m_rep[2] ** 2 - m_rep[0] ** 2))

    def propagate(Um, levels, x_end, steps=4000):
        Hf = 2 * Um @ np.diag(levels).astype(complex) @ Um.conj().T          # i dpsi/dx = Hf psi, flavour rows
        psi = np.eye(3, dtype=complex)                                        # column a = initial flavour a
        h = x_end / steps
        f = lambda y: -1j * (Hf @ y)
        for _ in range(steps):
            k1 = f(psi); k2 = f(psi + h / 2 * k1); k3 = f(psi + h / 2 * k2); k4 = f(psi + h * k3)
            psi = psi + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        return np.abs(psi) ** 2                                               # [final, initial]

    def closed(Um, levels, x):
        amp = (Um * np.exp(-2j * np.array(levels) * x)[None, :]) @ Um.conj().T
        return np.abs(amp) ** 2

    osc = {'assignment_sin2': list(best[2]), 'r_from_the_quoted_levels': r_rep}
    for name, r in (('levels_of_the_repository', r_rep), ('reference_splittings_input', 0.029805)):
        lv = [0.0, r, 1.0]
        Pn = propagate(U, lv, math.pi / 2)
        Pc = closed(U, lv, math.pi / 2)
        xs = np.linspace(0, 700, 70001)
        app = max(float(np.max(np.abs(np.einsum('bk,ak,xk->xba', U, U.conj(), np.exp(-2j * np.outer(xs, lv)))) ** 2 * (1 - np.eye(3))[None])), 0.0)
        osc[name] = {'P_ee': float(Pn[0, 0]), 'P_mu_to_e': float(Pn[0, 1]), 'P_mumu': float(Pn[1, 1]), 'P_mu_to_tau': float(Pn[2, 1]),
                     'integration_against_closed_form': float(np.abs(Pn - Pc).max()), 'unitarity_defect': float(np.abs(Pn.sum(axis=0) - 1).max()),
                     'largest_appearance_on_the_grid_to_700': app, 'distinct_frequencies': len({round(abs(a - b), 12) for i, a in enumerate(lv) for b in lv[i + 1:]}),
                     'slow_period_in_fast_periods': 1 / r}
    Pe = propagate(U, [0.3, 0.3, 0.3], math.pi / 2)
    Pd = propagate(U, [0.0, 0.0, 1.0], math.pi / 2)
    Pid = propagate(np.eye(3, dtype=complex), [0.0, r_rep, 1.0], math.pi / 2)
    osc['controls'] = {'three_equal_levels_largest_appearance': float((Pe * (1 - np.eye(3))).max()),
                       'two_equal_levels_distinct_frequencies': 1, 'two_equal_levels_largest_appearance_at_pi_over_2': float((Pd * (1 - np.eye(3))).max()),
                       'neutral_subgroup_equal_to_the_charged_one_largest_appearance': float((Pid * (1 - np.eye(3))).max())}
    out['F_oscillation'] = osc
    log('F', json.dumps(osc))

    # ---- G. checks against the written expectations ----------------------------------------------------------
    chk = {'rotations': A['rotations_determinant_plus_1'] == EX['rotations'], 'subgroups': A['subgroups'] == EX['subgroups'] and A['subgroups_by_order'] == EX['subgroups_by_order'],
           'basis_fixing': sum(counts.values()) == EX['basis_fixing'], 'pairs': Bres['ordered_pairs_per_triplet'] == EX['ordered_pairs_per_triplet'],
           'patterns_H': menus['H']['patterns_on_3'] == EX['patterns_menu_H'], 'patterns_all': Bres['patterns_on_3'] == EX['patterns_all'] and Bres['patterns_on_the_conjugate'] == EX['patterns_all'],
           'every_candidate_pattern_found_with_equal_entries_and_counts': all(r.get('found_on_3') and r.get("found_on_3'") and r['pair_counts_equal_3'] and r["pair_counts_equal_3'"] and r['largest_entry_difference_3'] <= TOL['entry_equality'] and r['exact_strings_against_own_entries'] <= TOL['entry_equality'] for r in comp),
           'no_own_pattern_absent_from_the_candidate': len(missing_in_candidate) == 0,
           'transposes_and_conjugates_agree': all(r['transpose_agrees'] and r['conjugate_agrees'] for r in comp),
           'fixed_columns': Bres['fixed_columns_with_G_e_of_order_3'] == EX['fixed_columns_with_G_e_of_order_3'] and all(c['pair_counts_equal'] for c in col_cmp),
           'menu_H_sizes': ang['H']['labelled_matrices'] == EX['menu_H_labelled_matrices'] and ang['H']['angle_triples'] == EX['menu_H_triples'],
           'menu_H_triples_equal_as_sets': ang['H']['own_triples_absent_from_the_candidate'] == 0 and ang['H']['candidate_triples_absent_here'] == 0,
           'menu_X_A_R_sizes': ang['X']['angle_triples'] == EX['menu_X_triples'] and ang['A']['angle_triples'] == EX['menu_A_triples'] and ang['R']['angle_triples'] == EX['menu_R_triples'],
           'nearest': abs(ang['H']['nearest_chi2'] - EX['nearest_chi2_menu_H']) <= TOL['chi2_abs'] and max(abs(a - b) for a, b in zip(ang['H']['nearest_sin2'], EX['nearest_sin2'])) <= 1e-5,
           'no_triple_within_3_sigma': ang['H']['within_3_sigma_on_all_three'] == 0 and ang['A']['within_3_sigma_on_all_three'] == 0,
           'haar_share': abs(p_f - canda['null']['share_chi2_at_or_below_nearest']) <= TOL['haar_probability_abs'],
           'oscillation': abs(osc['levels_of_the_repository']['P_ee'] - 0.936) <= TOL['oscillation_abs'] and abs(osc['levels_of_the_repository']['P_mu_to_tau'] - 0.677) <= TOL['oscillation_abs']}
    out['G_checks'] = chk
    out['finished'] = date(); out['cpu_seconds'] = round(time.process_time() - T0, 1)
    log('G', json.dumps(chk))
    json.dump(out, open(os.path.join(HERE, 'RESULT.json'), 'w'), indent=1, default=str)
    log('done', out['finished'], 'cpu', out['cpu_seconds'])


if __name__ == '__main__':
    main()
