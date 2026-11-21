"""Spin and statistics: the central involution of the binary icosahedral cover and the gradings of the source.

The cover acts on the five-level code by rho(q) = R(q) (+) S(sigma q): the rotation of the carrier on the three
vector levels and the SU(2) matrix of the Galois-conjugate quaternion on the two spinor levels. Both factors
are computed in exact Q(sqrt 5) arithmetic. The gradings are solved over GF(2).
"""
import itertools

from .exact import (QONE, QMINUS, ZERO, ONE, GENERATORS, binary_icosahedral, rotation, qmul, qgalois, qnorm,
                    det3, apply3, is_identity3, icosahedron)
from .model import patterns


# ------------------------------------------------------------------ the cover and its centre
def cover():
    group = binary_icosahedral()
    gset = set(group)
    rots = {q: rotation(q) for q in group}
    centre = [q for q in group if all(qmul(q, g) == qmul(g, q) for g in GENERATORS)]
    trivial_rotation = [q for q in group if is_identity3(rots[q])]
    galois_image = {qgalois(q) for q in group}
    galois_hom = all(qgalois(qmul(g, h)) == qmul(qgalois(g), qgalois(h)) for g in GENERATORS for h in group)
    galois_closed = all(qmul(a, b) in galois_image for a in (qgalois(GENERATORS[0]), qgalois(GENERATORS[1]))
                        for b in galois_image)
    # rho(q) = R(q) (+) S(sigma q); S(p) is left multiplication by p on H = C^2, so S(-1) = -I2
    z = QMINUS
    vector_part = [1, 1, 1] if is_identity3(rots[z]) else None
    spinor_part = [-1, -1] if qgalois(z) == QMINUS else None
    rho_z_diagonal = vector_part + spinor_part if vector_part and spinor_part else None
    # faithfulness: rho(q) = 1 needs R(q) = 1 (q = +-1) and sigma(q) = 1 (q = 1)
    kernel = [q for q in trivial_rotation if qgalois(q) == QONE]
    # the face placements of defects: rotations permute the faces of the carrier
    verts, _, faces = icosahedron()
    vindex = {v: i for i, v in enumerate(verts)}
    face_index = {frozenset(f): i for i, f in enumerate(faces)}
    face_kernel, face_closed = [], True
    for q in group:
        img = [vindex.get(apply3(rots[q], v)) for v in verts]
        if None in img:
            face_closed = False
            continue
        perm = [face_index.get(frozenset(img[x] for x in f)) for f in faces]
        if None in perm:
            face_closed = False
            continue
        if perm == list(range(len(faces))):
            face_kernel.append(q)
    return dict(
        order=len(group), unit_quaternions=all(qnorm(q) == ONE for q in group),
        closed=all(qmul(a, b) in gset for a in GENERATORS for b in group),
        rotations_special_orthogonal=all(det3(m) == ONE and all(
            sum((m[k][i] * m[k][j] for k in range(3)), ZERO) == (ONE if i == j else ZERO)
            for i in range(3) for j in range(3)) for m in rots.values()),
        centre_size=len(centre), centre_is_plus_minus_one=set(centre) == {QONE, QMINUS},
        elements_with_trivial_rotation=len(trivial_rotation),
        galois_conjugation_is_homomorphism=galois_hom, galois_image_order=len(galois_image),
        galois_image_closed=galois_closed,
        representation_kernel_size=len(kernel),
        rho_z_diagonal=rho_z_diagonal,
        rotations_permute_faces=face_closed, elements_fixing_every_face=len(face_kernel),
        z_fixes_every_face=QMINUS in face_kernel)


def central_involutions_of_G():
    """Centre of S(U(3) x U(2)) is {(a I3, b I2): a^3 b^2 = 1}; an involution has a, b in {+1, -1}."""
    return sorted([a] * 3 + [b] * 2 for a in (1, -1) for b in (1, -1) if a ** 3 * b ** 2 == 1)


def z_on_patterns(z_diag):
    rows = []
    for p in patterns():
        val = 1
        for i in range(5):
            if p['mask'] >> i & 1:
                val *= z_diag[i]
        six_y = int(p['six_Y'])
        rows.append(dict(mask=p['mask'], z=val, minus_one_to_n2=(-1) ** p['n2'], minus_one_to_6Y=(-1) ** (six_y % 2),
                         class3_holonomy=-1 if (3 * p['k']) % 6 == 3 else 1))
    return rows


def exp_i_pi_6Y_on_code():
    """exp(i pi 6Y) on the code, 6Y = diag(-2, -2, -2, 3, 3): the sign (-1)^(6Y) of each level."""
    return [(-1) ** (six_y % 2) for six_y in (-2, -2, -2, 3, 3)]


# ------------------------------------------------------------------ gradings over GF(2)
S1, S2 = 3, 4
PATS = [dict(mask=m, n3=bin(m & 0b00111).count('1'), n2=bin(m & 0b11000).count('1')) for m in range(32)]
for _p in PATS:
    _p['deg'] = _p['n3'] + _p['n2']
    _p['k'] = (3 * _p['n2'] - 2 * _p['n3']) % 6
BLOCKS = sorted({(p['n3'], p['n2']) for p in PATS})
VARS = [('q', m) for m in range(32)] + [('h', 1), ('h', 2)] + [('d', k) for k in range(1, 6)] + [('U', b) for b in BLOCKS]
IX = {v: i for i, v in enumerate(VARS)}
NV = len(VARS)
# one family of left-handed patterns: Q, u^c, d^c, L, e^c, nu^c (masks over v1 v2 v3 s1 s2)
FAMILY = {'Q': 0b01001, 'u^c': 0b00110, 'd^c': 0b11110, 'L': 0b01111, 'e^c': 0b11000, 'nu^c': 0b00000}
YUKAWA = (('Q', 'u^c', 1), ('Q', 'd^c', 2), ('L', 'e^c', 2), ('L', 'nu^c', 1))
FUSIONS = ((1, 1, 2), (5, 5, 4), (1, 4, 5), (2, 5, 1), (1, 2, 3), (4, 5, 3), (3, 5, 2), (1, 3, 4), (2, 2, 4),
           (4, 4, 2), (2, 3, 5), (3, 4, 1))
SPLITS = ((2, 1, 1), (4, 5, 5), (3, 1, 2), (3, 4, 5), (1, 2, 5), (5, 1, 4), (4, 2, 2), (2, 4, 4), (4, 1, 3),
          (2, 3, 5), (5, 2, 3), (1, 3, 4))
PAIRS = ((1, 5), (2, 4), (3, 3))


def _row(*vs, rhs=0):
    r = [0] * (NV + 1)
    for v in vs:
        if v is not None:
            r[IX[v]] ^= 1
    r[NV] = rhs
    return r


def _d(k):
    return None if k % 6 == 0 else ('d', k % 6)


def vertex_rows():
    rows, log = [], {}
    n = 0
    for p in PATS:            # currents: move one level inside the vector block or inside the spinor block
        a = p['mask']
        for i, j in list(itertools.permutations(range(3), 2)) + [(S1, S2), (S2, S1)]:
            if a >> i & 1 and not a >> j & 1:
                rows.append(_row(('q', a), ('q', a ^ (1 << i) ^ (1 << j))))
                n += 1
    log['currents'] = n
    n = 0
    for p in PATS:            # toggle of a spinor level with the register quantum of that level
        a = p['mask']
        for k, lev in ((1, S1), (2, S2)):
            if not a >> lev & 1:
                rows.append(_row(('q', a | (1 << lev)), ('q', a), ('h', k)))
                n += 1
    log['toggle'] = n
    n = 0
    for p in PATS:            # gauge-covariant hopping c*_x U c_y inside one block
        for p2 in PATS:
            if (p['n3'], p['n2']) == (p2['n3'], p2['n2']):
                rows.append(_row(('q', p['mask']), ('q', p2['mask']), ('U', (p['n3'], p['n2']))))
                n += 1
    log['hopping'] = n
    for f1, f2, hk in YUKAWA:
        rows.append(_row(('q', FAMILY[f1]), ('q', FAMILY[f2]), ('h', hk)))
    log['Yukawa'] = len(YUKAWA)
    for p in PATS:            # CPT: complement pattern, conjugate class
        rows.append(_row(('q', p['mask']), ('q', 31 ^ p['mask'])))
    for k in range(1, 6):
        rows.append(_row(_d(k), _d(6 - k)))
    log['CPT'] = 37
    ch = []
    for a, b, c in FUSIONS:
        ch.append(_row(_d(a), _d(b), _d(c)))
    for a, b, c in SPLITS:
        ch.append(_row(_d(a), _d(b), _d(c)))
    for a, b in PAIRS:
        ch.append(_row(_d(a), _d(b)))
    log['catalogue_channels'] = len(ch)
    return rows, ch, log


def gf2_solve(rows):
    a = [list(r) for r in rows]
    piv, r = [], 0
    for c in range(NV):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        for i in range(len(a)):
            if i != r and a[i][c]:
                a[i] = [x ^ y for x, y in zip(a[i], a[r])]
        piv.append(c)
        r += 1
    consistent = not any(row[NV] and not any(row[:NV]) for row in a)
    x = [0] * NV
    for i, c in enumerate(piv):
        x[c] = a[i][NV]
    basis = []
    for f in (c for c in range(NV) if c not in piv):
        v = [0] * NV
        v[f] = 1
        for i, c in enumerate(piv):
            if a[i][f]:
                v[c] = 1
        basis.append(v)
    return consistent, x, basis


def satisfies(rows, v):
    return all(sum(r[i] & v[i] for i in range(NV)) % 2 == r[NV] for r in rows)


def describe(v):
    q = {(p['n3'], p['n2']): v[IX[('q', p['mask'])]] for p in PATS}
    return dict(quanta_by_block={f'n3={a},n2={b}': q[(a, b)] for a, b in BLOCKS},
                quanta_constant_on_blocks=all(v[IX[('q', p['mask'])]] == q[(p['n3'], p['n2'])] for p in PATS),
                register=[v[IX[('h', 1)]], v[IX[('h', 2)]]],
                defects={str(k): v[IX[('d', k)]] for k in range(1, 6)},
                links=sorted({v[IX[('U', b)]] for b in BLOCKS}))


def assignment(rule):
    """A Z2 value on every object for a named rule."""
    v = [0] * NV
    for p in PATS:
        i = IX[('q', p['mask'])]
        v[i] = {'rotation': 1, 'six_Y': p['n2'] % 2, 'degree': p['deg'] % 2, 'class_k': 0}[rule]
    for k in (1, 2):
        v[IX[('h', k)]] = 1 if rule == 'six_Y' else 0
    for k in range(1, 6):
        v[IX[('d', k)]] = (k % 2) if rule in ('six_Y', 'class_k') else 0
    return v


def gradings():
    rows, ch, log = vertex_rows()
    # class-only gradings of the defects kept by every channel
    didx = [IX[('d', k)] for k in range(1, 6)]
    defect_solutions = []
    for bits in itertools.product((0, 1), repeat=5):
        v = [0] * NV
        for i, b in zip(didx, bits):
            v[i] = b
        if satisfies(ch, v):
            defect_solutions.append(list(bits))
    ok, _, basis = gf2_solve(rows + ch)
    decomposition = []
    for b in basis:
        d = describe(b)
        qb = d['quanta_by_block']
        alpha = qb['n3=0,n2=0']
        beta = (qb['n3=0,n2=1'] + alpha) % 2
        delta = (qb['n3=1,n2=0'] + alpha) % 2
        linear = all(qb[f'n3={a},n2={c}'] == (alpha + beta * c + delta * (a in (1, 2))) % 2 for a, c in BLOCKS)
        decomposition.append(dict(quanta_number=alpha, weak_centre_with_register=beta, coloured_quanta=delta,
                                  linear_in_these=linear, register=d['register'], defects=d['defects'],
                                  links=d['links'],
                                  trivial_on_physical_states=(alpha == 0 and delta == 0)))
    spin = assignment('rotation')
    named = {r: satisfies(rows + ch, assignment(r)) for r in ('rotation', 'six_Y', 'degree', 'class_k')}
    equal_to_spin = {r: assignment(r) == spin for r in ('six_Y', 'degree', 'class_k')}
    return dict(vertex_rows=log, consistent=ok, dimension=len(basis), number_of_gradings=2 ** len(basis),
                basis=decomposition,
                physical_generators=sorted({(r['quanta_number'], r['coloured_quanta'])
                                            for r in decomposition if not r['trivial_on_physical_states']}),
                defect_class_gradings=defect_solutions,
                rotation_grading=describe(spin), rotation_grading_is_a_grading=named['rotation'],
                rules_that_are_gradings=named, rules_equal_to_rotation=equal_to_spin)


def composites():
    """Rotation grading on composites: sum of the grading bits of the constituents."""
    spin = assignment('rotation')
    q = lambda m: spin[IX[('q', m)]]
    return dict(
        meson_quark_antiquark=sorted({(q(a) + q(31 ^ a)) % 2 for a in range(32)}),
        baryon_three_quanta=sorted({(q(a) + q(b) + q(c)) % 2 for a, b, c in itertools.product(range(32), repeat=3)}),
        pure_flux_contents=sorted({sum(spin[IX[('d', k)]] for k in ms) % 2 for ms in
                                   ((2, 4), (1, 5), (3, 3), (2, 2, 4, 4), (4, 4, 4), (2, 2, 2), (2, 3, 3, 4),
                                    (1, 2, 3), (3, 4, 5))}))


def links_two_sided():
    """z is constant on each block (n3, n2), so its two-sided action on End(block) is +1."""
    zd = [1, 1, 1, -1, -1]
    out = {}
    for b in BLOCKS:
        vals = set()
        for p in PATS:
            if (p['n3'], p['n2']) == b:
                v = 1
                for i in range(5):
                    if p['mask'] >> i & 1:
                        v *= zd[i]
                vals.add(v)
        out[f'n3={b[0]},n2={b[1]}'] = dict(z_values=sorted(vals), two_sided=sorted({x * y for x in vals for y in vals}))
    return out


def port_flux_parity():
    """Premise: vector levels carry load 2, spinor levels load -3, seen at flux 11 * load mod 20.

    Checked: the central lift S = -I2 multiplies every port coherent state by -1, so it acts by (-1)^k at flux
    k; the sectors of the levels then give the signs of rho(z)."""
    sectors = {'vector': (11 * 2) % 20, 'spinor': (11 * -3) % 20}
    signs = [(-1) ** sectors['vector']] * 3 + [(-1) ** sectors['spinor']] * 2
    pattern_rule = all(((11 * (2 * p['n3'] - 3 * p['n2'])) % 20) % 2 == p['n2'] % 2 for p in PATS)
    return dict(sectors=sectors, signs=signs, pattern_sector_parity_equals_n2=pattern_rule)
