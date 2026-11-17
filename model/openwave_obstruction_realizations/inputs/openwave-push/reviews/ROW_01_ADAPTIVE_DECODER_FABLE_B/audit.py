#!/usr/bin/env python3
"""RA1 audit: own recomputation for part A (adaptive phase decoder of row 1) and part B (neutral-start
record-only pointer join). Standard library only. This file imports no candidate and no producer code:
packages are read as data (JSON) or as text.

    audit.py own     own recomputation, writes RESULT.json
    audit.py final   reads RESULT.json of the `own` run back and embeds replay/COMPARE.json, the output of the
                     differential harness replay/compare_candidate.py (no recomputation in this step)
"""
import ast
import hashlib
import json
import math
import random
import subprocess
import sys
from collections import Counter, defaultdict, deque
from copy import deepcopy
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CONS = HERE / 'consumed'
SNAP = CONS / 'candidate_snapshot2'
FIRST = CONS / 'candidate'
PI = math.pi
TOL = 1e-9
MANIFEST_SHA = '85f9eefcf4c2292c1a2f47f522149cf1e28852705376d3abc2a65024243c7f3a'
CHECKS = []


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check(name, ok, detail=''):
    CHECKS.append({'check': name, 'holds': bool(ok), 'detail': str(detail)})
    print(('PASS  ' if ok else 'FAIL  ') + name + (' | ' + str(detail) if detail != '' else ''), flush=True)
    return bool(ok)


# ----------------------------------------------------------------------------- carrier, read as data
GEO = json.loads((CONS / 'row_completion_ROW_01/GEOMETRY.json').read_text())
SEAMS = [tuple(e) for e in GEO['seams']]
DUAL = [tuple(e) for e in GEO['dual_edges']]
LOOPS = [[(e, s) for e, s in loop] for loop in GEO['dual_port_boundaries']]
CYCLE = (0, 1, 8, 6, 2, 10)
SWEEP = list(range(7, 30)) + list(range(7))


def geometry_checks():
    ok = len(SEAMS) == 30 and len(set(SEAMS)) == 30 and all(0 <= a < b < 12 for a, b in SEAMS)
    deg = Counter(v for e in SEAMS for v in e)
    ok = ok and all(deg[v] == 5 for v in range(12))
    check('0.carrier: 30 seams on 12 ports, five per port', ok)
    ok = len(LOOPS) == 12
    use = defaultdict(list)
    for v, loop in enumerate(LOOPS):
        ok = ok and len(loop) == 5 and len({e for e, _ in loop}) == 5
        for i, (e, s) in enumerate(loop):
            a, b = SEAMS[e]
            ok = ok and ((s == 1 and v == a) or (s == -1 and v == b))
            head = DUAL[e][1] if s == 1 else DUAL[e][0]
            e2, s2 = loop[(i + 1) % 5]
            tail2 = DUAL[e2][0] if s2 == 1 else DUAL[e2][1]
            ok = ok and head == tail2
            use[e].append(s)
    ok = ok and all(sorted(use[e]) == [-1, 1] for e in range(30))
    check('0.carrier: every dual cell is a closed oriented pentagon, its signs are the incidence of the port, '
          'every dual edge bounds two cells with opposite signs', ok)
    ok = all(tuple(sorted((CYCLE[i], CYCLE[(i + 1) % 6]))) in set(SEAMS) for i in range(6))
    check('0.carrier: the pointer cycle 0,1,8,6,2,10 is a closed path of six seams', ok)
    tree = ast.parse((CONS / 'theory_repo/openwave_row1_integrity/strict_tape_validator.py').read_text())
    lit = None
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, 'id', None) == 'SEAMS' for t in node.targets):
            lit = ast.literal_eval(node.value)
    check('0.carrier: the seam list in the public validator (parsed as text) equals GEOMETRY.json',
          lit is not None and [tuple(e) for e in lit] == SEAMS)


# ----------------------------------------------------------------------------- the endpoint rule, own forms
def law(na, nb, c):
    d = na - nb
    num = d - (2 * c - 1) * (abs(d) % 2)
    assert num % 2 == 0
    return num // 2


def law_balanced(na, nb, c):
    tot = na + nb
    return na - (tot // 2 + c * (tot % 2))


# ============================================================================= PART A
def gv(v):
    return 0.31 * math.sin(v + 1)


def amp(a, b):
    return 0.23 * math.cos(a + b + 1)


def consts(e, deformed):
    if not deformed:
        return 0.0, 0.0, 0.0
    a, b = DUAL[e]
    return gv(a), gv(b), amp(a, b)


def declared_field(F, deformed):
    """The declared encoding of the public package: z_e(t) = exp(i phi_e(t)),
    phi_e(t) = g_a + (g_b - g_a) t + 2 pi F_e t + A_e sin(2 pi t)."""
    def field(e, t):
        ga, gb, A = consts(e, deformed)
        ph = ga + (gb - ga) * t + 2 * PI * F[e] * t + A * math.sin(2 * PI * t)
        return complex(math.cos(ph), math.sin(ph))
    return field


def rule(F):
    return {e: 8 * abs(F[e]) + 8 for e in range(30)}


def exact_degree(F, region):
    return sum(s * F[e] for v in region for e, s in LOOPS[v])


def divergence(F):
    q = [0] * 12
    for (a, b), j in zip(SEAMS, F):
        q[a] += j
        q[b] -= j
    return q


def principal(z0, z1):
    return math.atan2(z1.imag * z0.real - z1.real * z0.imag, z1.real * z0.real + z1.imag * z0.imag)


def own_decode(field, res, region, cap=None, guard=PI / 2):
    """Own decoder. The field is an argument; `res` is the number of subdivisions per edge."""
    edges = sorted({e for v in region for e, _ in LOOPS[v]})
    required = sum(res[e] + 1 for e in edges)
    if cap is not None and required > cap:
        return {'status': 'resource_stop', 'degree': None, 'required': required}
    total = 0.0
    max_inc = 0.0
    defect = 0.0
    for v in region:
        contour = []
        for e, s in LOOPS[v]:
            n = res[e]
            vals = [field(e, k / n) for k in range(n + 1)]
            if s == -1:
                vals.reverse()
            if contour:
                defect = max(defect, abs(contour[-1] - vals[0]))
                contour.extend(vals[1:])
            else:
                contour.extend(vals)
        defect = max(defect, abs(contour[0] - contour[-1]))
        if defect >= TOL:
            return {'status': 'open_curve', 'degree': None, 'required': required, 'defect': defect}
        incs = [principal(a, b) for a, b in zip(contour, contour[1:])]
        max_inc = max(max_inc, max(map(abs, incs)))
        if max_inc >= guard:
            return {'status': 'increment_guard', 'degree': None, 'required': required, 'max_increment': max_inc}
        total += math.fsum(incs) / (2 * PI)
    deg = round(total)
    if abs(total - deg) >= TOL:
        return {'status': 'non_integer', 'degree': None, 'required': required, 'float': total}
    return {'status': 'read', 'degree': deg, 'float': total, 'required': required,
            'max_increment': max_inc, 'defect': defect}


def part_a_bound():
    lo, hi = Fraction(333, 106), Fraction(355, 113)
    check('A1.rational bracket 333/106 < pi < 355/113', float(lo) < PI < float(hi))
    check('A1.analytic bound: (2 pi m + 0.62 + 0.46 pi)/(8m+8) < pi/4 for every m >= 0, because it is '
          '0.62 + 0.46 pi < 2 pi, exact in rationals',
          Fraction(62, 100) + Fraction(46, 100) * hi < 2 * lo)
    dgs = [abs(gv(b) - gv(a)) for a, b in DUAL]
    amps = [abs(amp(a, b)) for a, b in DUAL]
    check('A1.envelope: |g_b - g_a| <= 0.62 and |A_e| <= 0.23 on the 30 dual edges',
          max(dgs) <= 0.62 and max(amps) <= 0.23, f'max |dg| = {max(dgs):.6f}, max |A| = {max(amps):.6f}')

    def phase(e, F, t):
        ga, gb, A = consts(e, True)
        return ga + (gb - ga) * t + 2 * PI * F * t + A * math.sin(2 * PI * t)

    def closed(e, F):
        ga, gb, A = consts(e, True)
        n = 8 * abs(F) + 8
        return abs(gb - ga + 2 * PI * F) / n + abs(A) * math.sin(2 * PI / n)

    worst = 0.0
    for e in range(30):
        for F in (-128, -31, -4, -3, -1, 0, 1, 2, 3, 4, 7, 31, 128):
            n = 8 * abs(F) + 8
            direct = max(abs(phase(e, F, (k + 1) / n) - phase(e, F, k / n)) for k in range(n))
            worst = max(worst, abs(direct - closed(e, F)))
    check('A1.closed form of the largest increment, |g_b - g_a + 2 pi F|/n + |A| sin(2 pi/n), equals the '
          'direct evaluation (30 edges, 13 currents)', worst < 1e-12, f'max difference {worst:.3e}')
    table = []
    for m in (0, 1, 2, 3, 4, 8, 16, 128, 12494):
        n = 8 * m + 8
        best = max((closed(e, s * m), e, s) for e in range(30) for s in (1, -1))
        env = (0.62 + 2 * PI * m) / n + 0.23 * math.sin(2 * PI / n)
        cand = (2 * PI * m + 0.62 + 2 * PI * 0.23) / n
        table.append({'abs_F': m, 'n': n, 'worst_on_declared_geometry': best[0], 'edge': best[1], 'sign': best[2],
                      'worst_on_envelope': env, 'candidate_bound': cand, 'pi_over_4_minus_bound': PI / 4 - cand,
                      'closed_margin': (1.54 * PI - 0.62) / (8 * (m + 1))})
    check('A1.ordering at every tabulated |F|: declared geometry <= envelope <= candidate bound < pi/4',
          all(r['worst_on_declared_geometry'] <= r['worst_on_envelope'] + 1e-15
              and r['worst_on_envelope'] <= r['candidate_bound'] + 1e-15 and r['candidate_bound'] < PI / 4
              for r in table))
    check('A1.margin identity: pi/4 - bound = (1.54 pi - 0.62)/(8(|F|+1))',
          all(abs(r['pi_over_4_minus_bound'] - r['closed_margin']) < 1e-12 for r in table))
    envs = [(0.62 + 2 * PI * m) / (8 * m + 8) + 0.23 * math.sin(2 * PI / (8 * m + 8)) for m in range(20001)]
    check('A1.envelope worst case increases with |F| and stays below pi/4 (|F| = 0 .. 20000): supremum pi/4, '
          'not attained', all(x < y for x, y in zip(envs, envs[1:])) and envs[-1] < PI / 4,
          f'E(0) = {envs[0]:.6f}, E(4) = {envs[4]:.6f}, E(128) = {envs[128]:.6f}, E(20000) = {envs[-1]:.9f}, '
          f'pi/4 = {PI / 4:.9f}')
    cap_half = max(m for m in range(64) if (0.62 + 2 * PI * m + 0.46 * PI) / 64 < PI / 2)
    cap_pi = max(m for m in range(64) if (0.62 + 2 * PI * m + 0.46 * PI) / 64 < PI)
    return {'table': table, 'closed_form_max_difference': worst, 'max_abs_dg': max(dgs), 'max_abs_amp': max(amps),
            'fixed_64_subdivisions_certified_abs_F_with_guard_pi_over_2': cap_half,
            'fixed_64_subdivisions_certified_abs_F_for_exactness': cap_pi}


def part_a_witness():
    w = json.loads((CONS / 'background/G1_row1_strengthen/proposed_public/witness_four_turns.json').read_text())
    n = list(w['G'])
    F = [0] * 30
    for e, c in [tuple(x) for x in w['events']] + list(zip(SWEEP, map(int, w['sweep_coins']))):
        a, b = SEAMS[e]
        j = law(n[a], n[b], c)
        n[a] -= j
        n[b] += j
        F[e] += j
    region = w['region']
    exact = exact_degree(F, region)
    check('A2.witness: the tape ends with -4 turns on seam 1 and max |F_e| = 4',
          F[1] == -4 and max(map(abs, F)) == 4)
    check('A2.witness: B F = G - n, and the regional degree is sum_R (G - n)',
          divergence(F) == [g - x for g, x in zip(w['G'], n)] and exact == sum(w['G'][v] - n[v] for v in region),
          f'degree {exact}, final loads {n}')
    rec = json.loads((SNAP / 'ADAPTIVE_RESULT.json').read_text())['native_witness']
    check('A2.witness: own current, loads and degree equal the amendment\'s record',
          F == rec['current'] and n == rec['final_loads'] and exact == rec['exact_degree'])
    res = rule(F)
    r = own_decode(declared_field(F, True), res, region, cap=100000)
    edges = sorted({e for v in region for e, _ in LOOPS[v]})
    own_res = {str(e): res[e] for e in edges}
    check('A2.witness: own adaptive read is exact', r['status'] == 'read' and r['degree'] == exact,
          f"float {r.get('float')!r}, max increment {r.get('max_increment')!r}, samples {r['required']}")
    a = rec['amended_readout']
    check('A2.witness: resolutions, sample count, largest increment and degree equal the amendment\'s output',
          own_res == a['resolutions'] and r['required'] == a['required_samples']
          and abs(r['max_increment'] - a['max_principal_increment']) < 1e-12
          and abs(r['float'] - a['phase_degree']) < 1e-12,
          f"own {r['max_increment']!r} against {a['max_principal_increment']!r}")
    r64 = own_decode(declared_field(F, True), {e: 64 for e in range(30)}, region)
    check('A2.witness: 64 subdivisions per edge, a density fixed without the tape, read the same degree once '
          'the cap |F_e| <= 3 is lifted', r64['status'] == 'read' and r64['degree'] == exact,
          f"max increment {r64.get('max_increment')!r}")
    return {'current': F, 'final_loads': n, 'exact_degree': exact, 'own_read': r, 'own_resolutions': own_res,
            'fixed64_read': r64}


def a_cases():
    rng = random.Random(2026111401)
    cases = []
    for ampl in (0, 1, 2, 3, 4, 7, 31, 128):
        for _ in range(12):
            F = [rng.randint(-ampl, ampl) for _ in range(30)]
            if ampl:
                F[rng.randrange(30)] = ampl * rng.choice((-1, 1))
            region = [v for v in range(12) if rng.getrandbits(1)]
            cases.append((F, region))
    return cases


def part_a_sample(cases):
    st = {'reads': 0, 'exact': 0, 'max_error': 0.0, 'max_increment': 0.0, 'below_pi_over_4': 0,
          'nonzero_degrees': 0, 'degree_values': set()}
    fixed = Counter()
    for F, region in cases:
        q = exact_degree(F, region)
        assert q == sum(divergence(F)[v] for v in region)
        for deformed in (False, True):
            r = own_decode(declared_field(F, deformed), rule(F), region)
            st['reads'] += 1
            st['exact'] += r['status'] == 'read' and r['degree'] == q
            st['nonzero_degrees'] += q != 0
            st['degree_values'].add(q)
            if r['status'] == 'read':
                st['max_error'] = max(st['max_error'], abs(r['float'] - q))
                st['max_increment'] = max(st['max_increment'], r['max_increment'])
                st['below_pi_over_4'] += r['max_increment'] < PI / 4
        if max(map(abs, F)) <= 15:
            r = own_decode(declared_field(F, True), {e: 64 for e in range(30)}, region)
            fixed['cases'] += 1
            fixed['exact'] += r['status'] == 'read' and r['degree'] == q
    st['degree_values'] = [min(st['degree_values']), max(st['degree_values'])]
    check('A2.sample: every own adaptive read equals sum_R (B F)_v (signed currents up to |F_e| = 128, both '
          'deformation flags)', st['exact'] == st['reads'],
          f"{st['exact']} of {st['reads']}, max error {st['max_error']:.3e}, {st['nonzero_degrees']} reads with "
          f"nonzero degree, degrees from {st['degree_values'][0]} to {st['degree_values'][1]}")
    check('A1.sample: every sampled increment is below pi/4', st['below_pi_over_4'] == st['reads'],
          f"largest {st['max_increment']:.6f} against pi/4 = {PI / 4:.6f}")
    check('A3.fixed density: 64 subdivisions per edge read every sample tape with |F_e| <= 15 exactly',
          fixed['exact'] == fixed['cases'], f"{fixed['exact']} of {fixed['cases']}")
    st['fixed64'] = dict(fixed)
    return st


def part_a_mutated(cases):
    """The record F sets the resolution; the field is built from another current F'."""
    st = Counter()
    for F, region in cases:
        bnd = [e for e, (a, b) in enumerate(SEAMS) if (a in region) != (b in region)]
        if not bnd:
            st['no_boundary_edge'] += 1
            continue
        res = rule(F)
        qt = exact_degree(F, region)
        Fm = [-j for j in F]
        qf = exact_degree(Fm, region)
        if qf != qt:
            r = own_decode(declared_field(Fm, True), res, region)
            st['M1.cases'] += 1
            st['M1.reads_the_field'] += r.get('degree') == qf
            st['M1.reads_the_tape'] += r.get('degree') == qt
        e = bnd[0]
        sg = 1 if F[e] >= 0 else -1
        Fm = list(F)
        Fm[e] += sg
        qf = exact_degree(Fm, region)
        r = own_decode(declared_field(Fm, True), res, region)
        st['M2.cases'] += 1
        st['M2.reads_the_field'] += r.get('degree') == qf
        st['M2.reads_the_tape'] += r.get('degree') == qt
        Fm = list(F)
        Fm[e] = sg * (3 * abs(F[e]) + 3)
        r = own_decode(declared_field(Fm, True), res, region)
        st['M3.cases'] += 1
        st['M3.rejected_without_degree'] += r['status'] == 'increment_guard' and r['degree'] is None
        Fm = list(F)
        Fm[e] += res[e]
        qf = exact_degree(Fm, region)
        r = own_decode(declared_field(Fm, True), res, region)
        st['M4.cases'] += 1
        st['M4.reads_the_field'] += r.get('degree') == qf
        st['M4.reads_the_tape'] += r.get('degree') == qt
    check('A3.M1 field of the reversed tape, resolution from the tape: the read is the field\'s degree',
          st['M1.reads_the_field'] == st['M1.cases'] and st['M1.reads_the_tape'] == 0 and st['M1.cases'] > 0,
          f"field {st['M1.reads_the_field']} of {st['M1.cases']}, tape {st['M1.reads_the_tape']}")
    check('A3.M2 one more turn on a boundary edge: the read is the field\'s degree',
          st['M2.reads_the_field'] == st['M2.cases'] and st['M2.reads_the_tape'] == 0,
          f"field {st['M2.reads_the_field']} of {st['M2.cases']}, tape {st['M2.reads_the_tape']}")
    check('A4.M3 3|F_e|+3 turns on a boundary edge (increments between pi/2 and pi): rejected, no degree',
          st['M3.rejected_without_degree'] == st['M3.cases'], f"{st['M3.rejected_without_degree']} of {st['M3.cases']}")
    check('A3.M4 n_e more turns on a boundary edge: the samples are those of the tape\'s field and the read is '
          'the tape\'s degree (limit of every sampled reading)',
          st['M4.reads_the_tape'] == st['M4.cases'] and st['M4.reads_the_field'] == 0,
          f"tape {st['M4.reads_the_tape']} of {st['M4.cases']}, field {st['M4.reads_the_field']}")
    return dict(st)


def part_a_rejections():
    out = {}
    F = [0] * 30
    F[1] = 0.5
    r = own_decode(declared_field(F, True), {e: 8 for e in range(30)}, [0])
    out['half_turn_field'] = r
    check('A4.own decoder: a half turn on a boundary edge leaves the curve open (defect 2) and gives no degree',
          r['status'] == 'open_curve' and r['degree'] is None and abs(r['defect'] - 2.0) < 1e-9, f"defect {r.get('defect')}")
    F = [0] * 30
    F[1] = 12494
    r = own_decode(declared_field(F, True), rule(F), [0], cap=100000)
    out['cap_inside'] = {k: r.get(k) for k in ('status', 'degree', 'float', 'required', 'max_increment')}
    F2 = list(F)
    F2[1] = 12495
    r2 = own_decode(declared_field(F2, True), rule(F2), [0], cap=100000)
    out['cap_outside'] = r2
    check('A4.resource cap 100000 on one dual cell: |F_e| = 12494 needs 99997 samples and reads exactly, '
          '12495 needs 100005 and stops', r['status'] == 'read' and r['degree'] == exact_degree(F, [0])
          and r['required'] == 99997 and r2['status'] == 'resource_stop' and r2['required'] == 100005,
          f"error at 12494: {abs(r['float'] - r['degree']):.3e}, largest increment {r['max_increment']:.6f}")
    return out


# ============================================================================= PART B
def configs(q):
    out = [tuple([q] * 12)]
    for p in range(12):
        for m in range(12):
            if p != m:
                n = [q] * 12
                n[p] += 1
                n[m] -= 1
                out.append(tuple(n))
    return out


def in_sector(n, q):
    dev = [x - q for x in n]
    return sum(dev) == 0 and all(d in (-1, 0, 1) for d in dev) and dev.count(1) <= 1 and dev.count(-1) <= 1


def part_b_sector():
    out = {}
    check('B2.rule: the two own forms of the endpoint rule agree (450 branches)',
          all(law(a, b, c) == law_balanced(a, b, c) for a in range(-7, 8) for b in range(-7, 8) for c in (0, 1)))
    rng = random.Random(7)
    ok = True
    for _ in range(3000):
        na, nb, c, q = rng.randint(-50, 50), rng.randint(-50, 50), rng.getrandbits(1), rng.randint(-10**6, 10**6)
        ok = ok and law(na + q, nb + q, c) == law(na, nb, c)
    check('B2.translation covariance: the current depends on the difference and the coin only (3000 random '
          'pairs and shifts)', ok)
    per = {}
    types = Counter()
    for q in (1, 7, 100000, 2, 3, 12):
        cf = configs(q)
        assert len(cf) == 133 and len(set(cf)) == 133
        nbr = closed = unit = cov = 0
        for n in cf:
            for a, b in SEAMS:
                for c in (0, 1):
                    j = law(n[a], n[b], c)
                    o = list(n)
                    o[a] -= j
                    o[b] += j
                    nbr += 1
                    closed += in_sector(o, q)
                    unit += abs(j) <= 1
                    cov += law(n[a] - q, n[b] - q, c) == j
                    if q == 1:
                        d = n[a] - n[b]
                        types['no difference, no current' if d == 0 else 'adjacent pair annihilates' if abs(d) == 2
                              else 'a unit moves' if j else 'a unit holds'] += 1
        per[q] = {'branches': nbr, 'stay_in_sector': closed, 'current_at_most_one': unit, 'covariant': cov}
    theirs = sum(per[q]['branches'] for q in (1, 7, 100000))
    check('B3.count: 3 fillings x 133 configurations x 30 seams x 2 coins = 23,940 branches', theirs == 23940)
    check('B2.single-pair sector is closed under every event and every current is -1, 0 or +1 (six fillings, '
          '47,880 branches)', all(v['branches'] == v['stay_in_sector'] == v['current_at_most_one'] == v['covariant'] == 7980
                                  for v in per.values()), json.dumps(dict(types)))
    start = set()
    for a, b in SEAMS:
        for s in (1, -1):
            n = [1] * 12
            n[a] -= s
            n[b] += s
            start.add(tuple(n))
    seen = set(start)
    todo = list(start)
    while todo:
        n = todo.pop()
        for a, b in SEAMS:
            for c in (0, 1):
                j = law(n[a], n[b], c)
                o = list(n)
                o[a] -= j
                o[b] += j
                o = tuple(o)
                if o not in seen:
                    seen.add(o)
                    todo.append(o)
    check('B2.reachability: the 60 prepared pair states reach exactly the 133 configurations',
          len(start) == 60 and seen == set(configs(1)))
    uni = tuple([1] * 12)
    check('B2.the uniform state is absorbing: no event creates a pair',
          all(law(uni[a], uni[b], c) == 0 for a, b in SEAMS for c in (0, 1)))
    lift = {0: 0, 1: 1, 5: -1}
    cuts = 0
    ok = True
    for n in configs(1):
        dev = [x - 1 for x in n]
        sums = [0] * 4096
        for mask in range(1, 4096):
            low = mask & -mask
            sums[mask] = sums[mask ^ low] + dev[low.bit_length() - 1]
        for Q in sums:
            ok = ok and Q in (-1, 0, 1) and lift.get((-Q) % 6) == -Q
            cuts += 1
    check('B2.bounded regional charge: every region of every configuration has charge -1, 0 or +1 and the '
          'counter modulo six decodes it', ok, f'{cuts} cuts')
    # controls: laws under which the gates fail
    half = swap = 0
    for n in configs(1):
        for a, b in SEAMS:
            d = n[a] - n[b]
            half += Fraction(d, 2).denominator != 1
            swap += abs(d) > 1
    check('B2.control: under a half-current law 2 x {} of 7,980 branches leave the integer sector, and under '
          'a swap law (j = d) {} seam visits carry a current of two: both gates can fail'.format(half, swap),
          half > 0 and swap > 0)
    two = [2, 2, 0, 0] + [1] * 8
    check('B2.domain edge: two prepared pairs give a regional charge of two, outside the decoder',
          not in_sector(two, 1) and sum(x - 1 for x in two[:2]) == 2)
    out.update(per_filling={str(k): v for k, v in per.items()}, branch_types_at_filling_one=dict(types),
               regional_cuts=cuts, half_law_fractional_branches=2 * half, swap_law_double_currents=swap)
    return out


class Reject(Exception):
    pass


def need(ok, msg):
    if not ok:
        raise Reject(msg)


def is_int(x):
    return type(x) is int


def own_source(q, pair_edge, sign, prefix, sweep_coins, order=SWEEP, hidden=None):
    n = [q] * 12
    a, b = SEAMS[pair_edge]
    n[a] -= sign
    n[b] += sign
    recs = [{'kind': 'prepare_pair', 'seam': pair_edge, 'current': sign}]

    def ev(kind, e, c):
        a, b = SEAMS[e]
        d = n[a] - n[b]
        j = law(n[a], n[b], c)
        n[a] -= j
        n[b] += j
        r = {'kind': kind, 'seam': e, 'coin': c, 'current': j}
        if kind == 'E':
            r['contrast'] = d
        recs.append(r)

    def disturb(a, b):
        if n[a] < 1:
            a, b = b, a
        n[a] -= 1
        n[b] += 1

    for e, c in prefix:
        ev('C', e, c)
    if hidden and hidden[0] == 'pre':
        disturb(hidden[1], hidden[2])
    for e, c in zip(order, sweep_coins):
        ev('E', e, c)
    if hidden and hidden[0] == 'post':
        disturb(hidden[1], hidden[2])
    return {'records': recs, 'sweep_order': list(order)}, list(n)


def own_observe(q, tape, region):
    need(is_int(q) and q >= 1, 'filling')
    need(isinstance(region, (list, tuple)) and all(is_int(v) and 0 <= v < 12 for v in region)
         and len(set(region)) == len(region), 'region')
    need(isinstance(tape, dict) and isinstance(tape.get('records'), list), 'records')
    recs = tape['records']
    order = tape.get('sweep_order')
    need(len(recs) >= 31, 'length')
    need(isinstance(order, list) and all(is_int(e) for e in order) and sorted(order) == list(range(30)), 'order')
    n = [q] * 12
    cnt = [0] * 12
    first_e = len(recs) - 30
    for k, r in enumerate(recs):
        need(isinstance(r, dict), 'record')
        e, j = r.get('seam'), r.get('current')
        need(is_int(e) and 0 <= e < 30 and is_int(j), 'seam or current')
        a, b = SEAMS[e]
        if k == 0:
            need(r.get('kind') == 'prepare_pair' and j in (-1, 1), 'preparation')
        else:
            c = r.get('coin')
            need(is_int(c) and c in (0, 1), 'coin')
            d = n[a] - n[b]
            need(j == law(n[a], n[b], c), 'current differs from the rule')
            if k >= first_e:
                need(r.get('kind') == 'E' and e == order[k - first_e], 'closing order')
                need(is_int(r.get('contrast')) and r['contrast'] == d, 'closing contrast')
            else:
                need(r.get('kind') == 'C', 'kind')
        n[a] -= j
        n[b] += j
        cnt[a] = (cnt[a] + j) % 6
        cnt[b] = (cnt[b] - j) % 6
        need(in_sector(n, q), 'sector')
    res = sum(cnt[v] for v in region) % 6
    need(res in (0, 1, 5), 'counter')
    deg = {0: 0, 1: 1, 5: -1}[res]
    need(deg == sum(q - n[v] for v in region), 'counter against reconstruction')
    return {'degree': deg, 'charge': -deg, 'attempts': len(recs) - 1, 'loads': n}


def closing_only_charge(tape, region):
    """Regional charge from the thirty closing records alone: no filling, no preparation record, no prefix."""
    delta = [0] * 12
    eqs = []
    for r in tape['records'][-30:]:
        a, b = SEAMS[r['seam']]
        eqs.append((a, b, r['contrast'] - delta[a] + delta[b]))
        delta[a] -= r['current']
        delta[b] += r['current']
    m = {0: 0}
    todo = [0]
    while todo:
        v = todo.pop()
        for a, b, w in eqs:
            if a == v and b not in m:
                m[b] = m[a] - w
                todo.append(b)
            elif b == v and a not in m:
                m[a] = m[b] + w
                todo.append(a)
    if len(m) != 12 or any(m[a] - m[b] != w for a, b, w in eqs):
        return None
    n = [m[v] + delta[v] for v in range(12)]
    mean = Fraction(sum(n), 12)
    return sum((n[v] - mean for v in region), Fraction(0))


# ----------------------------------------------------------------------------- pointer, own state machine
def route(start, finish=0):
    prev = {start: None}
    dq = deque([start])
    while dq:
        v = dq.popleft()
        if v == finish:
            break
        for a, b in SEAMS:
            if v == a:
                w = b
            elif v == b:
                w = a
            else:
                continue
            if w not in prev:
                prev[w] = v
                dq.append(w)
    path = []
    cur = finish
    while cur != start:
        path.append(cur)
        cur = prev[cur]
    return tuple(reversed(path))


ROUTE = {v: route(v) for v in range(12)}


def p_request(st, w):
    if st[0] == 'D':
        return SEAMS[st[1]], 'discover', None, st[3]
    if st[0] == 'R':
        known, dest = st[3], st[1][0]
        return (min(known, dest), max(known, dest)), 'reset_route', None, known
    cursor, known = st[1], st[4]
    direction = 1 if w > 0 else -1
    dest = CYCLE[(cursor + direction) % 6]
    return (min(known, dest), max(known, dest)), 'measure', direction, known


def after_route(w):
    return ('END', 'success', 0) if w == 0 else ('M', 0, 0, 0, 0)


def p_step(st, coin, w, R):
    (a, b), stage, direction, pos = p_request(st, w)
    j = law(1 if pos == a else 0, 1 if pos == b else 0, coin)
    rec = {'stage': stage, 'seam': [a, b], 'coin': coin, 'current': j, 'requested_cycle_direction': direction}
    if st[0] == 'D':
        ei, t = st[1], st[2]
        if j:
            known = b if j > 0 else a
            nxt = ('R', ROUTE[known], 0, known) if ROUTE[known] else after_route(w)
        else:
            t += 1
            if t == R:
                ei, t = ei + 1, 0
            nxt = ('END', 'discovery_timeout', 0) if ei == 30 else ('D', ei, t, pos)
    elif st[0] == 'R':
        path, t, known = st[1], st[2], st[3]
        if j:
            known, path = path[0], path[1:]
            nxt = ('R', path, 0, known) if path else after_route(w)
        else:
            t += 1
            nxt = ('END', 'routing_timeout', 0) if t == R else ('R', path, t, known)
    else:
        cursor, steps, t, known = st[1], st[2], st[3], st[4]
        if j:
            cursor = (cursor + direction) % 6
            steps += 1
            known = CYCLE[cursor]
            nxt = ('END', 'success', direction * steps) if steps == 6 * abs(w) else ('M', cursor, steps, 0, known)
        else:
            t += 1
            nxt = ('END', 'measurement_timeout', direction * steps) if t == R else ('M', cursor, steps, t, known)
    return nxt, rec


def own_pointer_tape(x, w, R, coins):
    st = ('D', 0, 0, x)
    recs = []
    pos = x
    it = iter(coins)
    while st[0] != 'END':
        st, rec = p_step(st, next(it), w, R)
        recs.append(rec)
        if rec['current'] == 1:
            pos = rec['seam'][1]
        elif rec['current'] == -1:
            pos = rec['seam'][0]
    ok = st[1] == 'success'
    return {'success': ok, 'reason': 'complete_closed_pointer_word' if ok else st[1], 'records': recs,
            'final_loads': [int(v == pos) for v in range(12)], 'retained_signed_cycle_steps': st[2],
            'winding': st[2] // 6 if ok else None, 'fine_history_probability': str(Fraction(1, 2 ** len(recs)))}


def pointer_compatible(tape, x, w, R):
    recs = tape['records']
    st = ('D', 0, 0, x)
    for r in recs:
        if st[0] == 'END' or not isinstance(r, dict):
            return False
        (a, b), stage, direction, _ = p_request(st, w)
        c = r.get('coin')
        seam = r.get('seam')
        if not (is_int(c) and c in (0, 1)) or not isinstance(seam, (list, tuple)) or tuple(seam) != (a, b):
            return False
        if r.get('stage') != stage or r.get('requested_cycle_direction') != direction:
            return False
        st, rec = p_step(st, c, w, R)
        if not is_int(r.get('current')) or r['current'] != rec['current']:
            return False
    if st[0] != 'END':
        return False
    outcome, signed = st[1], st[2]
    if tape.get('success') is True:
        if outcome != 'success' or tape.get('reason') != 'complete_closed_pointer_word':
            return False
        if not is_int(tape.get('winding')) or tape['winding'] != w or signed != 6 * w:
            return False
    elif tape.get('success') is False:
        if outcome == 'success' or tape.get('reason') != outcome or tape.get('winding') is not None:
            return False
    else:
        return False
    if tape.get('retained_signed_cycle_steps') != signed:
        return False
    return tape.get('fine_history_probability') == str(Fraction(1, 2 ** len(recs)))


def own_admit(q, tape, region, pointer, R=2):
    try:
        src = own_observe(q, tape, region)
        need(isinstance(pointer, dict) and isinstance(pointer.get('records'), list), 'pointer tape')
        cols = [x for x in range(12) if pointer_compatible(pointer, x, src['degree'], R)]
        need(cols, 'no compatible pointer column')
        return {'admitted': True, 'readout': src['charge'] if pointer['success'] else None,
                'degree': src['degree'], 'columns': cols, 'N': src['attempts'] + len(pointer['records'])}
    except (Reject, TypeError, KeyError, IndexError, ValueError, AssertionError) as err:
        return {'admitted': False, 'reason': str(err)}


def vkey(v):
    return (True, v['readout'], v['degree'], v['N'], tuple(v['columns'])) if v['admitted'] else (False,)


def pointer_law(x, w, R):
    dist = {('D', 0, 0, x): Fraction(1)}
    out = defaultdict(Fraction)
    while dist:
        nxt = defaultdict(Fraction)
        for st, p in dist.items():
            for c in (0, 1):
                s2, _ = p_step(st, c, w, R)
                if s2[0] == 'END':
                    out[s2[1]] += p / 2
                else:
                    nxt[s2] += p / 2
        dist = nxt
    return out


def pointer_formula(x, w, R):
    h = Fraction(1, 2 ** R)
    s = Fraction(0)
    k = 0
    for a, b in SEAMS:
        if x in (a, b):
            y = b if x == a else a
            s += h ** k * (1 - h) ** (1 + len(ROUTE[y]) + 6 * abs(w))
            k += 1
    return s


def sector_dp(n0, seams):
    dist = {tuple(n0): Fraction(1)}
    for e in seams:
        a, b = SEAMS[e]
        nxt = defaultdict(Fraction)
        for n, p in dist.items():
            for c in (0, 1):
                j = law(n[a], n[b], c)
                o = list(n)
                o[a] -= j
                o[b] += j
                nxt[tuple(o)] += p / 2
        dist = nxt
    return dist


def part_b_laws():
    out = {}
    n_laws = 0
    ok = True
    laws = {}
    for R in (1, 2, 3):
        for x in range(12):
            for w in (-2, -1, 0, 1, 2):
                L = pointer_law(x, w, R)
                n_laws += 1
                ok = ok and sum(L.values()) == 1 and L['success'] == pointer_formula(x, w, R) \
                    and L['discovery_timeout'] == Fraction(1, 2 ** (5 * R))
                laws[(R, x, w)] = L
    check('B3.complete pointer law by own exact state recursion: success plus the three timeouts sum to one '
          'and success equals the closed formula (retry bounds 1, 2, 3; 12 columns; windings -2 .. 2)', ok,
          f'{n_laws} of {n_laws}')
    r2 = {w: [float(laws[(2, x, w)]['success']) for x in range(12)] for w in (0, 1)}
    out['success_probability_retry_bound_2'] = {'winding_0': [min(r2[0]), max(r2[0])],
                                                'winding_plus_or_minus_1': [min(r2[1]), max(r2[1])]}
    surv = []
    sums_ok = True
    for e, (a, b) in enumerate(SEAMS):
        for s in (1, -1):
            n = [1] * 12
            n[a] -= s
            n[b] += s
            d = sector_dp(n, SWEEP)
            sums_ok = sums_ok and sum(d.values()) == 1 and all(in_sector(k, 1) for k in d)
            surv.append((1 - d.get(tuple([1] * 12), Fraction(0)), e, s))
    lo, hi = min(surv), max(surv)
    mean = sum(p for p, _, _ in surv) / 60
    out['pair_survival_through_the_closing_sweep_empty_prefix'] = {
        'min': str(lo[0]), 'min_float': float(lo[0]), 'min_at': [lo[1], lo[2]], 'max': str(hi[0]),
        'max_float': float(hi[0]), 'max_at': [hi[1], hi[2]], 'mean_float': float(mean),
        'preparations_annihilated_with_certainty': sum(p == 0 for p, _, _ in surv)}
    check('B4.exact law of the sector through the closing sweep, 60 preparations with no prefix event: every '
          'law sums to one inside the sector', sums_ok,
          'the pair survives the sweep with probability between {:.4f} and {:.4f}, mean {:.4f}; {} preparations '
          'annihilate with certainty'.format(float(lo[0]), float(hi[0]), float(mean), sum(p == 0 for p, _, _ in surv)))
    return out, laws


def source_mutations(tape):
    out = []
    recs = tape['records']

    def mut(label, fn):
        t = deepcopy(tape)
        fn(t)
        out.append((label, t))

    def setf(k, key, f):
        return lambda t: t['records'][k].__setitem__(key, f(t['records'][k][key]))

    for s in range(30):
        if s != recs[0]['seam']:
            mut('prep.seam', setf(0, 'seam', lambda _, s=s: s))
    mut('prep.sign', setf(0, 'current', lambda j: -j))
    mut('prep.double', setf(0, 'current', lambda j: 2 * j))
    mut('prep.kind', setf(0, 'kind', lambda _: 'C'))
    for k in range(1, len(recs)):
        kind = recs[k]['kind']
        mut(kind + '.coin', setf(k, 'coin', lambda c: 1 - c))
        mut(kind + '.current+1', setf(k, 'current', lambda j: j + 1))
        mut(kind + '.current-1', setf(k, 'current', lambda j: j - 1))
        mut(kind + '.seam', setf(k, 'seam', lambda e: (e + 1) % 30))
        mut(kind + '.kind', setf(k, 'kind', lambda x: 'E' if x == 'C' else 'C'))
        if kind == 'E':
            mut('E.contrast+1', setf(k, 'contrast', lambda d: d + 1))
            mut('E.contrast-1', setf(k, 'contrast', lambda d: d - 1))
        mut(kind + '.delete', lambda t, k=k: t['records'].pop(k))
        mut(kind + '.duplicate', lambda t, k=k: t['records'].insert(k, deepcopy(t['records'][k])))
    mut('order.swap', lambda t: t['sweep_order'].__setitem__(slice(0, 2), t['sweep_order'][1::-1]))
    return out


def pointer_mutations(ptr):
    out = []

    def mut(label, fn):
        t = deepcopy(ptr)
        fn(t)
        out.append((label, t))

    def setf(k, key, f):
        return lambda t: t['records'][k].__setitem__(key, f(t['records'][k][key]))

    for k in range(len(ptr['records'])):
        mut('pointer.coin', setf(k, 'coin', lambda c: 1 - c))
        mut('pointer.current', setf(k, 'current', lambda j: j + 1))
        mut('pointer.stage', setf(k, 'stage', lambda s: 'measure' if s != 'measure' else 'discover'))
        mut('pointer.seam', setf(k, 'seam', lambda s: list(SEAMS[(SEAMS.index(tuple(s)) + 1) % 30])))
    mut('pointer.delete_last', lambda t: t['records'].pop())
    mut('pointer.duplicate_last', lambda t: t['records'].append(deepcopy(t['records'][-1])))
    mut('pointer.success_flag', lambda t: t.__setitem__('success', not t['success']))
    mut('pointer.reason', lambda t: t.__setitem__('reason', 'routing_timeout' if t['reason'] != 'routing_timeout'
                                                  else 'measurement_timeout'))
    mut('pointer.steps', lambda t: t.__setitem__('retained_signed_cycle_steps', t['retained_signed_cycle_steps'] + 6))
    return out


def valid_source_history(q, t):
    try:
        recs = t['records']
        prefix = [(r['seam'], r['coin']) for r in recs[1:-30]]
        coins = [r['coin'] for r in recs[-30:]]
        mine, _ = own_source(q, recs[0]['seam'], recs[0]['current'], prefix, coins, order=t['sweep_order'])
        return mine == t
    except Exception:
        return False


def valid_pointer_history(t, x, w, R):
    try:
        mine = own_pointer_tape(x, w, R, [r['coin'] for r in t['records']])
        return mine == t
    except Exception:
        return False


def load_examples():
    return json.loads((SNAP / 'NEUTRAL_RECORD_EXAMPLES.json').read_text())


def part_b_tapes(laws):
    ex = load_examples()
    out = {'tapes': len(ex)}
    st = Counter()
    outcomes = Counter()
    exp_s = Fraction(0)
    var = Fraction(0)
    exp_reason = defaultdict(Fraction)
    for x in ex:
        q, tape, region, ptr, v = x['filling'], x['source_records'], x['region'], x['pointer'], x['verdict']
        recs = tape['records']
        prefix = [(r['seam'], r['coin']) for r in recs[1:-30]]
        coins = [r['coin'] for r in recs[-30:]]
        mine, loads = own_source(q, recs[0]['seam'], recs[0]['current'], prefix, coins, order=tape['sweep_order'])
        st['source_tape_reproduced'] += mine == tape
        st['terminal_loads_reproduced'] += loads == x['audit_only']['final_loads']
        obs = own_observe(q, tape, region)
        st['charge_from_records_equals_charge_of_terminal_loads'] += obs['charge'] == sum(loads[u] - q for u in region)
        st['closing_records_alone_give_the_charge'] += closing_only_charge(tape, region) == obs['charge']
        own = own_admit(q, tape, region, ptr, 2)
        st['verdict_equal'] += (own['admitted'] is True and v['admitted'] is True and own['readout'] == v['readout']
                                and own['columns'] == v['compatible_pointer_initials']
                                and own['degree'] == v['source_degree'] and own['N'] == v['native_attempts']
                                and ptr['success'] == v['success'])
        wt = str(Fraction(1, 2 ** own['N']))
        st['weights_equal'] += v['fine_record_effect'] == [wt if u in own['columns'] else '0' for u in range(12)]
        st['audit_column_is_the_only_compatible_one'] += own['columns'] == [x['initial_pointer_for_audit_only']]
        st['pointer_tape_reproduced'] += own_pointer_tape(x['initial_pointer_for_audit_only'], own['degree'], 2,
                                                          [r['coin'] for r in ptr['records']]) == ptr
        outcomes[(ptr['reason'], v['readout'])] += 1
        L = laws[(2, x['initial_pointer_for_audit_only'], own['degree'])]
        exp_s += L['success']
        var += L['success'] * (1 - L['success'])
        for k, p in L.items():
            exp_reason[k] += p
        st['failures_admitted_without_readout'] += (not ptr['success']) and v['admitted'] and v['readout'] is None
        st['failures'] += not ptr['success']
    N = len(ex)
    check('B3.count: 288 recorded tapes', N == 288)
    for key in ('source_tape_reproduced', 'terminal_loads_reproduced',
                'charge_from_records_equals_charge_of_terminal_loads', 'closing_records_alone_give_the_charge',
                'verdict_equal', 'weights_equal', 'audit_column_is_the_only_compatible_one', 'pointer_tape_reproduced'):
        check('B3.288 tapes, own replay: ' + key.replace('_', ' '), st[key] == N, f'{st[key]} of {N}')
    check('B3.failures are kept: every timeout tape is admitted with no readout and the weight 2^-N of its '
          'full record', st['failures_admitted_without_readout'] == st['failures'] and st['failures'] > 0,
          f"{st['failures']} failures of {N}")
    obs_s = N - st['failures']
    z = (obs_s - float(exp_s)) / math.sqrt(float(var))
    check('B3.no filtering: successes among the 288 tapes against the exact law at retry bound 2',
          abs(z) < 3, f'observed {obs_s}, expected {float(exp_s):.2f}, standard deviation {math.sqrt(float(var)):.2f}, z = {z:+.2f}')
    out['outcomes'] = {f'{k[0]}|readout={k[1]}': n for k, n in sorted(outcomes.items(), key=str)}
    out['expected_outcomes_exact_law'] = {k: float(p) for k, p in exp_reason.items()}
    out['observed_successes'] = obs_s
    out['expected_successes'] = float(exp_s)
    out['z'] = z
    out['counts'] = dict(st)
    out['distinct_source_histories'] = len({json.dumps(x['source_records'], sort_keys=True) for x in ex})
    out['terminal_state_of_the_distinct_histories'] = sorted(
        {(x['case'], 'pair annihilated' if len(set(x['audit_only']['final_loads'])) == 1 else 'pair present')
         for x in ex})
    out['nonzero_readouts'] = sum(n for k, n in outcomes.items() if k[1] not in (0, None))
    return out, ex


def part_b_mutations(ex):
    out = {}
    cat = Counter()
    by_label = defaultdict(Counter)
    bad = 0
    firsts = {}
    for x in ex:
        firsts.setdefault((x['case'], tuple(x['region'])), x)
    for (case, region), x in sorted(firsts.items()):
        q, tape, ptr = x['filling'], x['source_records'], x['pointer']
        base = vkey(own_admit(q, tape, list(region), ptr, 2))
        assert base[0]
        for label, t in source_mutations(tape):
            v = vkey(own_admit(q, t, list(region), ptr, 2))
            if not v[0]:
                c = 'rejected'
            elif v == base:
                c = 'identical verdict'
                if not valid_source_history(q, t):
                    bad += 1
            elif v[1:3] != base[1:3]:
                c = 'admitted, other charge'
            else:
                c = 'admitted, same charge, other weight'
            cat[c] += 1
            by_label[label][c] += 1
    out['source_record_mutations'] = dict(cat)
    out['source_record_mutations_by_field'] = {k: dict(v) for k, v in sorted(by_label.items())}
    check('B1.source records: every single-field mutation is rejected, or changes the charge or the weight, '
          'or is another history the rule itself can produce from the same preparation (same charge, same '
          'length)', bad == 0 and cat['rejected'] > 0, json.dumps(dict(cat)))
    out['identical_verdicts_that_are_not_histories_of_the_rule'] = bad
    catp = Counter()
    by_label_p = defaultdict(Counter)
    badp = 0
    for x in ex[::3]:
        q, tape, region, ptr = x['filling'], x['source_records'], x['region'], x['pointer']
        base = own_admit(q, tape, region, ptr, 2)
        bk = vkey(base)
        src = own_observe(q, tape, region)
        for label, t in pointer_mutations(ptr):
            cols = [u for u in range(12) if pointer_compatible(t, u, src['degree'], 2)]
            v = (True, src['charge'] if t['success'] else None, src['degree'],
                 src['attempts'] + len(t['records']), tuple(cols)) if cols else (False,)
            if not v[0]:
                c = 'rejected'
            elif v == bk:
                c = 'identical verdict'
                if not valid_pointer_history(t, base['columns'][0], base['degree'], 2):
                    badp += 1
            else:
                c = 'admitted, other verdict'
            catp[c] += 1
            by_label_p[label][c] += 1
    out['pointer_record_mutations'] = dict(catp)
    out['pointer_record_mutations_by_field'] = {k: dict(v) for k, v in sorted(by_label_p.items())}
    check('B1.pointer records (96 of the 288 tapes): every single-field mutation is rejected or is another tape '
          'the pointer itself can produce (a coin on an attempt with no load difference)', badp == 0 and catp['rejected'] > 0
          and catp['admitted, other verdict'] == 0, json.dumps(dict(catp)))
    # terminal snapshots and the filling
    same = same_q = 0
    for x in ex:
        q, tape, region, ptr = x['filling'], x['source_records'], x['region'], x['pointer']
        base = vkey(own_admit(q, tape, region, ptr, 2))
        t2 = deepcopy(tape)
        t2['final_loads'] = [5] * 12
        p2 = deepcopy(ptr)
        p2['final_loads'] = [0] * 11 + [1] if ptr['final_loads'][11] == 0 else [1] + [0] * 11
        same += vkey(own_admit(q, t2, region, p2, 2)) == base
        same_q += all(vkey(own_admit(q2, tape, region, ptr, 2)) == base for q2 in (1, 2, 9, 10 ** 6) if q2 != q)
    check('B1.terminal snapshots: altered source and pointer terminal loads leave the verdict unchanged',
          same == len(ex), f'{same} of {len(ex)}')
    check('B1.the value of the filling does not enter the verdict: every other positive integer gives the same '
          'admission, charge, columns and weight', same_q == len(ex), f'{same_q} of {len(ex)}')
    out['verdict_unchanged_under_terminal_snapshot_mutation'] = same
    out['verdict_unchanged_under_other_filling'] = same_q
    return out


def own_histories():
    rng = random.Random(2026111402)
    hist = []
    for _ in range(400):
        q = rng.randint(1, 9)
        pe = rng.randrange(30)
        sign = rng.choice((-1, 1))
        prefix = [(rng.randrange(30), rng.getrandbits(1)) for _ in range(rng.randint(0, 12))]
        coins = [rng.getrandbits(1) for _ in range(30)]
        region = [v for v in range(12) if rng.getrandbits(1)]
        x = rng.randrange(12)
        pcoins = [rng.getrandbits(1) for _ in range(200)]
        ab = rng.sample(range(12), 2)
        hist.append({'q': q, 'pair_edge': pe, 'sign': sign, 'prefix': prefix, 'coins': coins, 'region': region,
                     'x': x, 'pcoins': pcoins, 'ab': ab})
    return hist


def part_b_own_sample():
    st = Counter()
    values = Counter()
    for h in own_histories():
        q, region = h['q'], h['region']
        tape, loads = own_source(q, h['pair_edge'], h['sign'], h['prefix'], h['coins'])
        obs = own_observe(q, tape, region)
        actual = sum(loads[v] - q for v in region)
        st['histories'] += 1
        st['charge_equals_terminal_charge'] += obs['charge'] == actual
        st['closing_records_alone_give_the_charge'] += closing_only_charge(tape, region) == actual
        ptr = own_pointer_tape(h['x'], obs['degree'], 2, h['pcoins'])
        v = own_admit(q, tape, region, ptr, 2)
        st['admitted'] += v['admitted'] and v['readout'] == (actual if ptr['success'] else None)
        values[(ptr['reason'], v.get('readout'))] += 1
        a, b = h['ab']
        t_pre, l_pre = own_source(q, h['pair_edge'], h['sign'], h['prefix'], h['coins'], hidden=('pre', a, b))
        st['unrecorded_transfer_before_the_sweep_rejected'] += not own_admit(q, t_pre, region, ptr, 2)['admitted']
        t_post, l_post = own_source(q, h['pair_edge'], h['sign'], h['prefix'], h['coins'], hidden=('post', a, b))
        v_post = own_admit(q, t_post, region, ptr, 2)
        st['unrecorded_transfer_after_the_sweep_same_records'] += t_post == tape
        st['unrecorded_transfer_after_the_sweep_same_verdict'] += vkey(v_post) == vkey(v)
        st['unrecorded_transfer_after_the_sweep_charge_differs_from_terminal'] += \
            obs['charge'] != sum(l_post[u] - q for u in region)
        wrong = 1 if obs['degree'] != 1 else -1
        pw = own_pointer_tape(h['x'], wrong, 2, h['pcoins'])
        if pw['success']:
            st['wrong_pointer_command_completed'] += 1
            st['wrong_pointer_command_rejected'] += not own_admit(q, tape, region, pw, 2)['admitted']
    N = st['histories']
    check('B1.own 400 histories: charge from the records equals the charge of the terminal loads',
          st['charge_equals_terminal_charge'] == N and st['closing_records_alone_give_the_charge'] == N
          and st['admitted'] == N, json.dumps({f'{k[0]}|{k[1]}': n for k, n in sorted(values.items(), key=str)}))
    check('B1.a transfer with no record before the sweep changes the closing records and is rejected',
          st['unrecorded_transfer_before_the_sweep_rejected'] == N, f"{st['unrecorded_transfer_before_the_sweep_rejected']} of {N}")
    check('B1.a transfer with no record after the sweep changes terminal loads and no record: same verdict',
          st['unrecorded_transfer_after_the_sweep_same_records'] == N
          and st['unrecorded_transfer_after_the_sweep_same_verdict'] == N,
          f"the admitted charge differs from the terminal charge in "
          f"{st['unrecorded_transfer_after_the_sweep_charge_differs_from_terminal']} of {N}")
    check('B1.a completed pointer tape of another command is rejected',
          st['wrong_pointer_command_rejected'] == st['wrong_pointer_command_completed'] > 0,
          f"{st['wrong_pointer_command_rejected']} of {st['wrong_pointer_command_completed']}")
    return {'counts': dict(st), 'outcomes': {f'{k[0]}|readout={k[1]}': n for k, n in sorted(values.items(), key=str)}}


def manifest_checks():
    out = {}
    man_path = SNAP / 'MANIFEST.json'
    out['manifest_sha256'] = sha(man_path)
    check('M.manifest SHA-256 equals the announced freeze', out['manifest_sha256'] == MANIFEST_SHA, out['manifest_sha256'])
    man = json.loads(man_path.read_text())
    pins = man['files_sha256']
    okp = {name: sha(SNAP / name) == h for name, h in pins.items()}
    check('M.the ten pinned files of the manifest are the bytes reviewed (copy of 16:40:50)', all(okp.values()),
          f'{sum(okp.values())} of {len(okp)}')
    changed = {}
    for name, h in pins.items():
        p = FIRST / name
        if p.exists():
            changed[name] = 'unchanged' if sha(p) == h else 'changed before the freeze'
        else:
            changed[name] = 'not in the copy of 16:33:37'
    out['first_reading_against_freeze'] = changed
    live = {}
    for rel, h in man['source_and_review_sha256'].items():
        if rel.endswith('hadrons/row1_review/RESULT.json'):
            live[rel] = 'opposing evidence: hashed after the verdict, see MANIFEST.json of this review'
            continue
        p = ROOT / rel
        live[rel] = 'equal' if p.exists() and sha(p) == h else 'DIFFERENT'
    out['source_pins_against_live_files'] = live
    check('M.source pins of the manifest equal the live files (14 of 15; the opposing evidence is read last)',
          all(v != 'DIFFERENT' for v in live.values()), f"{sum(v == 'equal' for v in live.values())} equal")
    pub = ROOT / 'reverse-engineering-reality/code/openwave_charge_quantization'
    snap_m = json.loads((pub / 'SNAPSHOT_MANIFEST.json').read_text())['files_sha256']
    bad = [rel for rel, h in snap_m.items() if sha(pub / 'workspace' / rel) != h]
    out['public_package'] = {'PUBLIC_MANIFEST.json': sha(pub / 'PUBLIC_MANIFEST.json'),
                             'SNAPSHOT_MANIFEST.json': sha(pub / 'SNAPSHOT_MANIFEST.json'),
                             'snapshot_files': len(snap_m), 'snapshot_files_different': len(bad)}
    check('M.public package unchanged: its two manifests have the hashes of the filing text and all {} '
          'snapshot files equal their pins'.format(len(snap_m)),
          not bad and out['public_package']['PUBLIC_MANIFEST.json'].startswith('d293ba59')
          and out['public_package']['SNAPSHOT_MANIFEST.json'].startswith('86c990ef'))
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'own'
    started = subprocess.run(['date', '+%Y-%m-%dT%H:%M:%S%z'], capture_output=True, text=True).stdout.strip()
    print('RA1 audit, mode', mode, 'started', started, 'python', sys.version.split()[0], flush=True)
    if mode == 'final':
        # no recomputation: the result of the `own` run is read back and the differential harness is embedded
        result = json.loads((HERE / 'RESULT.json').read_text())
        assert result['mode'] == 'own' and result['all_checks_hold'] in (True, False)
        cmp_path = HERE / 'replay' / 'COMPARE.json'
        cmp = json.loads(cmp_path.read_text())
        CHECKS.extend(result['checks'])
        for c in cmp['checks']:
            check('C.' + c['check'], c['holds'], c['detail'])
        result.update(mode='own, then final (comparison with the candidate copy embedded)', final_started=started,
                      comparison_with_candidate_copy=cmp, comparison_file_sha256=sha(cmp_path), checks=CHECKS,
                      checks_total=len(CHECKS), checks_holding=sum(c['holds'] for c in CHECKS),
                      all_checks_hold=all(c['holds'] for c in CHECKS))
        (HERE / 'RESULT.json').write_text(json.dumps(result, indent=1, default=str) + '\n')
        print(f"checks holding: {result['checks_holding']} of {result['checks_total']}", flush=True)
        return 0 if result['all_checks_hold'] else 1
    result = {'reviewer': 'Fable-B/ROW_01_ADAPTIVE_DECODER_FABLE_B', 'mode': mode, 'started': started,
              'python': sys.version.split()[0]}
    result['manifest'] = manifest_checks()
    geometry_checks()
    print('--- part A', flush=True)
    A = {'bound': part_a_bound(), 'witness': part_a_witness()}
    cases = a_cases()
    A['sample'] = part_a_sample(cases)
    A['mutated_field'] = part_a_mutated(cases)
    A['rejections'] = part_a_rejections()
    result['part_a'] = A
    print('--- part B', flush=True)
    B = {'sector': part_b_sector()}
    B['laws'], laws = part_b_laws()
    B['recorded_tapes'], ex = part_b_tapes(laws)
    B['mutations'] = part_b_mutations(ex)
    B['own_sample'] = part_b_own_sample()
    result['part_b'] = B
    result['checks'] = CHECKS
    result['checks_total'] = len(CHECKS)
    result['checks_holding'] = sum(c['holds'] for c in CHECKS)
    result['all_checks_hold'] = all(c['holds'] for c in CHECKS)
    result['finished'] = subprocess.run(['date', '+%Y-%m-%dT%H:%M:%S%z'], capture_output=True, text=True).stdout.strip()
    (HERE / 'RESULT.json').write_text(json.dumps(result, indent=1, default=str) + '\n')
    print(f"checks holding: {result['checks_holding']} of {result['checks_total']}; finished {result['finished']}", flush=True)
    return 0 if result['all_checks_hold'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
