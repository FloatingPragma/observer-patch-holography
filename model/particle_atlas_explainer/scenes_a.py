"""Opening and the four steps inside the machine."""
import math

from .canvas import (AMBER, BG, BLUE, CARD, CARD2, EDGE, FAINT, GOLD, GREEN, MONO, MUTED, ORANGE,
                     RED, SANS, SERIF, SOFT, TEAL, TEXT, VIOLET, n, text_width)
from .parts import (CW, ROOT3, SX, W, X0, X1, camera, card, carrier, header, hexkeep, lattice,
                    loop_triangle, mesh, mix, mobius, status, tri_fill, vortex)

STEPS = (('1', 'Machine', TEAL), ('2', 'Twist', TEAL), ('3', 'Protection', TEAL),
         ('4', 'Anatomy', TEAL), ('5', 'Motion', AMBER), ('6', 'Position', AMBER),
         ('7', 'Velocity', AMBER), ('8', 'Mass', AMBER), ('9', 'Fields', VIOLET))


def _centroid(pts, tri):
    return (sum(pts[k][0] for k in tri) / 3, sum(pts[k][1] for k in tri) / 3)


def _nearest(pts, tris, x, y):
    return min(tris, key=lambda t: (_centroid(pts, t)[0] - x) ** 2 + (_centroid(pts, t)[1] - y) ** 2)


def record_card(c, x, y, w, tick, det, op=1.0, h=46, checks=True):
    c.open(opacity=n(op))
    c.rect(x, y, w, h, '#0b1322', EDGE, 1.2, 11)
    c.circle(x + 22, y + h / 2, 5, AMBER)
    label = f'tick {tick}  detector {det}  count 1'.replace(' ', '\u00a0')
    c.text(x + 40, y + h / 2 + 5, label, 14, SOFT, family=MONO)
    if checks:
        for k, name in enumerate('AB'):
            bx = x + w - 64 + k * 30
            c.text(bx, y + h / 2 + 5, name, 13, GREEN, 700)
            c.tick(bx + 17, y + h / 2, GREEN, 6, 1.8)
    c.close()


def hero(c, data):
    c.caps(X0, 92, 'OPH particle atlas  ·  V3', TEAL, 15, ls=3.2)
    c.add(f'<text x="{X0}" y="178" font-family="{SANS}" font-size="80" font-weight="800" '
          f'letter-spacing="-2" fill="{TEXT}">From a <tspan fill="{AMBER}">twist</tspan> to a '
          f'<tspan fill="{VIOLET}">particle</tspan></text>')
    c.text(X0, 228, 'The proposed path from a protected source defect to an observer’s particle record',
           25, SOFT, maxw=CW)
    y = c.para(X0, 284, 1250,
               'The OPH simulation runs small observers that compare records and repair mismatches. '
               'Its hypothesis connects matter to protected {a|defects}. This illustrated path joins '
               'finite examples with the proposed physical interpretation. Each caption separates '
               'computed records from the proof needed to identify {v|position}, '
               '{v|velocity} and {v|mass} on one source.', 21, MUTED, 1.5)

    # --- triptych: program, records, description
    py, ph = 410, 292
    c.define('heroPanel', '<linearGradient id="heroPanel" x1="0" x2="1" y1="0" y2="0">'
                          '<stop offset="0" stop-color="#0c1d27"/><stop offset=".5" stop-color="#0d1424"/>'
                          '<stop offset="1" stop-color="#191533"/></linearGradient>')
    c.rect(X0, py, CW, ph, 'url(#heroPanel)', EDGE, 1.4, 26)

    # zone A: the mesh with one twist that has moved
    ax, ay = 345, py + ph / 2
    c.define('mAg', '<radialGradient id="mAg" cx="50%" cy="50%" r="62%">'
                    '<stop offset=".5" stop-color="#fff"/><stop offset="1" stop-color="#000"/></radialGradient>')
    c.define('mA', f'<mask id="mA" maskUnits="userSpaceOnUse" x="{X0}" y="{py}" width="450" height="{ph}">'
                   f'<rect x="{X0}" y="{py}" width="450" height="{ph}" fill="url(#mAg)"/></mask>')
    pts, edges, tris = lattice(ax, ay, 46, lambda x, y: abs(x - ax) < 232 and abs(y - ay) < 150)
    c.open(mask='url(#mA)')
    mesh(c, pts, edges, TEAL, 1.3, 0.42, 2.6, 0.85)
    c.close()
    trail = []
    for k in range(4):
        t = _nearest(pts, tris, ax + 62 - k * 23, ay - 6)
        if t not in trail:
            trail.append(t)
    for k, t in reversed(list(enumerate(trail))):
        if k == 0:
            tri_fill(c, pts, t, AMBER, 0.85, glow=58, cls='pulse')
        else:
            tri_fill(c, pts, t, AMBER, 0.34 - 0.09 * k, stroke=False)
    hx, hy = _centroid(pts, trail[0])
    c.arrow(hx + 34, hy, hx + 78, hy, AMBER, 2.2, 9, 0.9)

    # zone B: records
    bx = 614
    for k, (tick, det, op) in enumerate((('040', '06', 0.45), ('041', '07', 0.7), ('042', '08', 1.0))):
        record_card(c, bx, py + 66 + k * 58, 372, tick, det, op)

    # zone C: the smooth description
    zx, zy = 1040, py
    c.define('mCg', '<radialGradient id="mCg" cx="50%" cy="50%" r="62%">'
                    '<stop offset=".45" stop-color="#fff"/><stop offset="1" stop-color="#000"/></radialGradient>')
    c.define('mC', f'<mask id="mC" maskUnits="userSpaceOnUse" x="{zx}" y="{zy}" width="480" height="{ph}">'
                   f'<rect x="{zx}" y="{zy}" width="480" height="{ph}" fill="url(#mCg)"/></mask>')
    c.open(mask='url(#mC)')
    for i in range(8):
        yy = zy + 20 + i * 38
        c.curve([(zx + 10 + 460 * t / 30, yy + 11 * math.sin(t / 30 * math.pi + i * 0.45)) for t in range(31)],
                VIOLET, 1.1, 0.3)
    for j in range(10):
        xx = zx + 18 + j * 50
        c.curve([(xx + 13 * math.sin(t / 20 * math.pi + j * 0.35), zy + 8 + (ph - 16) * t / 20)
                 for t in range(21)], VIOLET, 1.1, 0.3)
    c.close()
    p0, p1, p2 = (zx + 70, zy + 232), (zx + 250, zy + 236), (zx + 356, zy + 104)
    prev = None
    for k in range(33):
        t = k / 32
        q = ((1 - t) ** 2 * p0[0] + 2 * t * (1 - t) * p1[0] + t * t * p2[0],
             (1 - t) ** 2 * p0[1] + 2 * t * (1 - t) * p1[1] + t * t * p2[1])
        if prev:
            c.line(prev[0], prev[1], q[0], q[1], mix(VIOLET, '#ffffff', t * 0.5), 1.2 + 2.2 * t,
                   round(0.12 + 0.8 * t, 2))
        prev = q
    c.glow(p2[0], p2[1], 56, VIOLET, 0.75, cls='pulse')
    c.glow(p2[0], p2[1], 22, AMBER, 0.95)
    c.circle(p2[0], p2[1], 8.5, AMBER, '#fff', 2)
    ang = math.atan2(p2[1] - p1[1], p2[0] - p1[0])
    c.arrow(p2[0] + 16 * math.cos(ang), p2[1] + 16 * math.sin(ang),
            p2[0] + 84 * math.cos(ang), p2[1] + 84 * math.sin(ang), TEXT, 2.4, 11)
    c.text(p2[0] + 96 * math.cos(ang) + 6, p2[1] + 96 * math.sin(ang) + 4, 'v', 26, TEXT, 400,
           family=SERIF, italic=True)
    c.line(p2[0], p2[1] + 16, p2[0], zy + ph - 34, VIOLET, 1.4, 0.8, '4 5')
    c.text(p2[0], zy + ph - 12, 'x', 26, TEXT, 400, 'middle', family=SERIF, italic=True)
    c.text(p2[0] - 30, p2[1] + 6, 'm', 26, TEXT, 400, 'end', family=SERIF, italic=True)

    for xx in (588, 1012):
        c.circle(xx, py + ph / 2, 17, CARD, EDGE, 1.3)
        c.arrow(xx - 7, py + ph / 2, xx + 8, py + ph / 2, SOFT, 2.2, 8)

    cy = py + ph + 38
    for x, col, lab, body in (
            (X0 + 8, TEAL, 'What the program runs', 'Observers, seams and repair records. Protection depends on the allowed moves.'),
            (bx + 4, AMBER, 'What gets recorded', 'Clicks that observers compare and agree on.'),
            (zx + 8, VIOLET, 'The physical target', 'A controlled particle description with position, velocity and mass.')):
        c.caps(x, cy, lab, col, 14)
        c.para(x, cy + 28, 372, body, 17, SOFT, 1.4)

    # --- chain strip
    sy = cy + 106
    gap, pw, phh = 10, (CW - 8 * 10) / 9, 48
    for a, b, lab, col in ((0, 3, 'Inside the machine', TEAL), (4, 7, 'The hand-over', AMBER),
                           (8, 8, 'The description', VIOLET)):
        xa, xb = X0 + a * (pw + gap), X0 + b * (pw + gap) + pw
        c.line(xa + 4, sy, xb - 4, sy, col, 2, 0.75)
        c.caps(xa + 4, sy - 12, lab, col, 12.5, ls=2)
    for i, (num, lab, col) in enumerate(STEPS):
        x = X0 + i * (pw + gap)
        c.rect(x, sy + 14, pw, phh, CARD, EDGE, 1.3, phh / 2)
        c.circle(x + 26, sy + 14 + phh / 2, 14, col)
        c.text(x + 26, sy + 14 + phh / 2 + 5.5, num, 15.5, BG, 800, 'middle')
        c.text(x + 50, sy + 14 + phh / 2 + 6, lab, 17, TEXT, 600, maxw=pw - 58)
    ly = sy + 104
    c.circle(X0 + 12, ly - 5, 6.5, TEAL)
    c.text(X0 + 30, ly, 'Computed in the model: a result in the stated construction', 16, MUTED)
    c.circle(X0 + 520, ly - 5, 6.5, 'none', AMBER, 1.8, dash='3.2 2.6')
    c.text(X0 + 538, ly, 'To derive: a link the program has to prove', 16, MUTED)
    c.glow(X0 + 922, ly - 5, 16, AMBER, 0.9)
    c.circle(X0 + 922, ly - 5, 5.5, AMBER)
    c.text(X0 + 942, ly, 'Amber: the illustrative defect followed through the chain', 16, MUTED)
    return ly + 30


# ---------------------------------------------------------------- step 1
def step1(c, data, y):
    yb = header(c, y, 1, 'Inside the machine', 'The machine compares, repairs and records',
                'The illustrated starting point is a loop of three moves.')
    h = 404
    card(c, X0, yb, CW, h)
    tx, tw = X0 + 36, 430
    yy = c.para(tx, yb + 50, tw,
                'An {t|observer} here is the smallest kind: a bookkeeper that holds a few records. '
                'Two neighbors overlap on a {t|seam}, where each can check the other.', 19.5, SOFT)
    yy = c.para(tx, yy + 12, tw,
                'At every tick, neighbors {b|compare}. A mismatch triggers a local {b|repair}. '
                'A reading that survives every comparison is kept as a {b|record}.', 19.5, SOFT)
    c.para(tx, yy + 12, tw,
           'The basic patch wires 12 observers with 30 seams, the corners and edges of an '
           'icosahedron. The physical target is to recover space, time and matter from '
           'the records made inside such a program.', 19.5, SOFT)

    # settle strip
    labels = ('a mismatch appears', 'repairs spread', 'settled')
    for k in range(3):
        cx, cy = 690 + k * 150, yb + 132
        pts, edges, tris = lattice(cx, cy, 27, hexkeep(cx, cy, 54.5))
        col = GREEN if k == 2 else TEAL
        mesh(c, pts, edges, col, 1.4, 0.55 if k < 2 else 0.8, 2.6)
        if k == 0:
            a, b = pts[(0, 0)], pts[(1, 0)]
            c.glow((a[0] + b[0]) / 2, (a[1] + b[1]) / 2, 30, RED, 0.6, cls='pulse')
            c.line(a[0], a[1], b[0], b[1], RED, 4)
        elif k == 1:
            for e in (((0, 0), (1, 0)),):
                c.line(pts[e[0]][0], pts[e[0]][1], pts[e[1]][0], pts[e[1]][1], GREEN, 3)
            for e in (((1, 0), (1, 1)), ((0, 0), (0, -1)), ((1, 0), (2, -1)), ((0, 0), (-1, 1))):
                if e[0] in pts and e[1] in pts:
                    c.line(pts[e[0]][0], pts[e[0]][1], pts[e[1]][0], pts[e[1]][1], RED, 2.6, 0.75, '5 5')
        else:
            c.circle(cx, cy, 15, BG, GREEN, 1.8)
            c.tick(cx, cy, GREEN, 8, 2.4)
        c.text(cx, yb + 226, labels[k], 15.5, MUTED, 400, 'middle')
        if k < 2:
            c.arrow(cx + 62, cy, cx + 88, cy, FAINT, 2, 8)
    # the loop of three moves
    ly = yb + 290
    c.rect(620, ly - 30, 440, 86, CARD2, EDGE, 1.2, 16)
    c.caps(840, ly - 6, 'one tick', FAINT, 12, 'middle')
    xs = (642, 786, 930)
    for i, (lab, col) in enumerate((('compare', TEAL), ('repair', RED), ('record', GREEN))):
        w = c.chip(xs[i], ly + 10, lab, col, 13, 13, 30)
        if i < 2:
            c.arrow(xs[i] + w + 8, ly + 25, xs[i + 1] - 8, ly + 25, FAINT, 1.8, 7)

    # the carrier with callouts
    kx, ky = 1296, yb + 184
    proj, depth = carrier(c, data, kx, ky, 80)
    faces = [f for i, f in enumerate(data['faces']) if depth[i] > 0]
    verts = sorted({v for f in faces for v in f})
    edges = sorted({tuple(sorted(e)) for f in faces for e in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0]))})
    v = max(verts, key=lambda i: proj[i][0] - proj[i][1])
    c.line(proj[v][0] + 6, proj[v][1] - 7, kx + 138, yb + 50, FAINT, 1.2)
    c.text(kx + 142, yb + 46, 'observer', 15.5, TEXT, 600, maxw=70)
    mid = lambda e: ((proj[e[0]][0] + proj[e[1]][0]) / 2, (proj[e[0]][1] + proj[e[1]][1]) / 2)
    e = min(edges, key=lambda e: mid(e)[0] + mid(e)[1])
    c.line(proj[e[0]][0], proj[e[0]][1], proj[e[1]][0], proj[e[1]][1], TEXT, 2.6)
    c.line(mid(e)[0] - 5, mid(e)[1] - 5, kx - 168, yb + 52, FAINT, 1.2)
    c.text(kx - 172, yb + 50, 'seam', 15.5, TEXT, 600, 'end')
    cen = lambda f: (sum(proj[k][0] for k in f) / 3, sum(proj[k][1] for k in f) / 3)
    f = max(faces, key=lambda f: cen(f)[0] * 0.6 + cen(f)[1])
    c.poly([proj[k][:2] for k in f], TEAL, TEXT, 2.2, 0.32)
    c.line(cen(f)[0] + 8, cen(f)[1] + 8, kx + 132, yb + 318, FAINT, 1.2)
    c.text(kx + 136, yb + 324, 'smallest', 15.5, TEXT, 600)
    c.text(kx + 136, yb + 343, 'loop', 15.5, TEXT, 600)
    c.text(kx, yb + h - 20, '12 observers  ·  30 seams  ·  20 triangles', 15.5, MUTED, 400, 'middle')
    return status(c, X0, yb + h + 16, CW,
                  'The repair law runs as code. Every outcome is kept as a record, holds and '
                  'failures included. One program, state, reference and clock are fixed before '
                  'anything is measured.',
                  'The carrier, the repair law, the scheduler and the detectors are design '
                  'choices of the construction. Running them needs no proof that they are the '
                  'only possible ones.', 'Declared')


# ---------------------------------------------------------------- step 2
def analogy(c, x, y, w, h, title, body):
    card(c, x, y, w, h)
    c.rect(x + 14, y + 14, w - 28, 132, CARD2, None, 0, 14)
    c.text(x + 28, y + 182, title, 20, TEXT, 700)
    c.para(x + 28, y + 212, w - 56, body, 17, MUTED, 1.45)


def step2(c, data, y):
    yb = header(c, y, 2, 'Inside the machine', 'A twist that no single seam holds',
                'A topological defect is a mismatch that belongs to a whole loop.')
    h = 476
    card(c, X0, yb, CW, h)
    tx = X0 + 36
    yy = c.para(tx, yb + 54, 318,
                'Alice, Bob and Charlie each hold a coin. Every seam carries one rule for the two '
                'coins it joins: {g|SAME} or {a|DIFFERENT}.', 19.5, SOFT)
    yy = c.para(tx, yy + 12, 318,
                'Multiply the rules around the loop, counting SAME as +1 and DIFFERENT as −1.',
                19.5, SOFT)
    ky = yb + h - 110
    c.circle(tx + 12, ky, 11, GOLD)
    c.text(tx + 12, ky + 4.5, 'H', 12, BG, 700, 'middle')
    c.circle(tx + 46, ky, 11, CARD, GOLD, 2)
    c.text(tx + 46, ky + 4.5, 'T', 12, GOLD, 700, 'middle')
    c.text(tx + 68, ky + 5, 'heads or tails', 15.5, MUTED)
    c.line(tx, ky + 32, tx + 56, ky + 32, GREEN, 3.2)
    c.text(tx + 68, ky + 37, 'rule met', 15.5, MUTED)
    c.line(tx, ky + 62, tx + 56, ky + 62, AMBER, 3.2, dash='7 6')
    c.text(tx + 68, ky + 67, 'rule broken', 15.5, MUTED)

    for cx, rules, title, col, prod, formula, cap in (
            (690, (1, 1, 1), 'Flat loop', GREEN, '+1', '+1 · +1 · +1 = +1',
             'All three rules can be met at once.'),
            (1010, (1, 1, -1), 'Twisted loop', AMBER, '−1', '+1 · +1 · −1 = −1',
             'Any two rules can be met. All three cannot.')):
        c.caps(cx, yb + 50, title, col, 14, 'middle')
        loop_triangle(c, cx, yb + 220, 196, rules, 'HHH', product=prod)
        c.text(cx, yb + 384, formula.replace(' ', ' '), 18, col, 700, 'middle', family=MONO)
        c.para(cx, yb + 416, 290, cap, 16, MUTED, 1.4, anchor='middle')

    # the mismatch moves and never leaves
    mx = 1352
    c.line(1196, yb + 40, 1196, yb + h - 40, EDGE, 1.2)
    c.caps(mx, yb + 50, 'repair moves it', AMBER, 13, 'middle')
    for k, coins in enumerate(('HHH', 'THH', 'TTH')):
        cy = yb + 112 + k * 110
        loop_triangle(c, mx - 26, cy, 78, (1, 1, -1), coins, small=True, show_rules=False)
        if k < 2:
            c.arrow(mx + 52, cy + 26, mx + 52, cy + 78, FAINT, 1.8, 7)
            c.text(mx + 64, cy + 57, 'repair', 13.5, FAINT)
    c.para(mx, yb + h - 54, 290, 'One seam is broken in every frame. The count stays at one.',
           16, MUTED, 1.4, anchor='middle')

    # definition bar
    dy = yb + h + 16
    c.rect(X0, dy, CW, 66, '#1c1708', AMBER, 1.3, 16)
    c.rect(X0, dy + 10, 6, 46, AMBER, None, 0, 3)
    c.para(X0 + 34, dy + 41, CW - 68,
           '{a|Topological defect:} a mismatch owned by the loop as a whole. It can sit on any '
           'seam of the loop. It cannot sit on none.', 20.5, TEXT)

    # everyday twists
    ay, aw, ah, gap = dy + 86, (CW - 48) / 3, 290, 24
    x1, x2, x3 = X0, X0 + aw + gap, X0 + 2 * (aw + gap)
    analogy(c, x1, ay, aw, ah, 'A wrinkle in a fitted carpet',
            'Press it flat here and it rises next door. The extra length belongs to the whole carpet.')
    fy = ay + 112

    def bump(cx0, wd, ht):
        return ''.join(f' L{n(cx0 - wd + 2 * wd * t / 28)},{n(fy - ht * math.exp(-((t / 28 - .5) / .2) ** 2))}'
                       for t in range(29))
    c.rect(x1 + 40, fy + 6, aw - 80, 5, '#1a2540', None, 0, 2)
    c.path(f'M{n(x1 + 40)},{fy}' + bump(x1 + 150, 64, 50) + f' L{n(x1 + aw - 40)},{fy}', TEAL, 3.6)
    c.path(f'M{n(x1 + 226)},{fy}' + bump(x1 + 300, 64, 50), AMBER, 2.4, op=0.9, dash='6 6')
    c.arrow(x1 + 150, ay + 30, x1 + 150, fy - 58, TEXT, 2.2, 9)
    c.text(x1 + 162, ay + 46, 'press', 14.5, MUTED)
    c.qarrow(x1 + 192, fy - 40, x1 + 226, fy - 78, x1 + 262, fy - 44, AMBER, 1.8, 8)

    analogy(c, x2, ay, aw, ah, 'A Möbius band',
            'Every short piece is an ordinary strip. The half turn belongs to the whole band, '
            'and only a cut removes it.')
    mobius(c, x2 + aw / 2, ay + 66, 90, 25)

    analogy(c, x3, ay, aw, ah, 'A whorl',
            'Comb the arrows smooth wherever you like. The point they circle stays: the full '
            'turn is counted on a ring around it.')
    c.add(f'<clipPath id="clipW"><rect x="{n(x3 + 14)}" y="{n(ay + 14)}" width="{n(aw - 28)}" '
          f'height="132" rx="14"/></clipPath>')
    c.open(**{'clip-path': 'url(#clipW)'})
    vortex(c, x3 + aw / 2, ay + 80, 230, 31, length=15, hot=60, fade=False)
    c.close()
    c.circle(x3 + aw / 2, ay + 80, 47, 'none', TEXT, 1.6, 0.9, '5 5')

    # the real carrier
    ry, rh = ay + ah + 24, 268
    card(c, X0, ry, CW, rh)
    carrier(c, data, X0 + 176, ry + rh / 2, 60, twisted=(0, 1), string=((0, 5),))
    c.text(X0 + 356, ry + 58, 'The same thing on the real carrier', 24, TEXT, 700)
    c.para(X0 + 356, ry + 94, 540,
           'Each of the 20 triangles is a smallest loop. Multiply the seam signs around a '
           'triangle: +1 is flat, −1 carries a twist. One flipped seam record ({r|red}) twists '
           'the two triangles that share it ({a|amber}). Further repair moves the broken seam '
           'around and leaves the twists in place.', 18, SOFT, 1.5)
    sx, sw = X0 + 946, 132
    for k, (big, small, col) in enumerate((('1 of 30', 'seam records flipped', RED),
                                           ('2,048', 'starting states, each run in 16 repair orders', TEAL),
                                           ('1', 'unit of disagreement left, every time', AMBER))):
        bx = sx + k * (sw + 14)
        c.rect(bx, ry + 34, sw, rh - 68, CARD2, EDGE, 1.2, 16)
        c.text(bx + sw / 2, ry + 98, big, 34, col, 800, 'middle', maxw=sw - 10)
        c.para(bx + sw / 2, ry + 132, sw - 22, small, 15, MUTED, 1.38, anchor='middle')
    return status(c, X0, ry + rh + 16, CW,
                  'A loop whose signs multiply to −1 has no fully settled state. That is a '
                  'finite theorem. The twelve-observer patch has 19 independent loops for such '
                  'signs to sit in.',
                  'The plain load records of the carrier have zero loop sums. The twists live in '
                  'added transport data: an ordering memory or signs on the seams.', 'Declared')


# ---------------------------------------------------------------- step 3
def sign_triangle(c, cx, cy, size, signs, title, hot=None):
    hh = size * ROOT3 / 2
    P = [(cx, cy - hh * 0.62), (cx + size / 2, cy + hh * 0.38), (cx - size / 2, cy + hh * 0.38)]
    for i in range(3):
        a, b = P[i], P[(i + 1) % 3]
        col = TEAL if signs[i] == 1 else AMBER
        c.line(a[0], a[1], b[0], b[1], col, 2.8)
        mx, my = (a[0] + b[0]) / 2, (a[1] + b[1]) / 2
        c.rect(mx - 19, my - 12, 38, 24, BG, col, 1.4, 12)
        c.text(mx, my + 5.5, '+1' if signs[i] == 1 else '−1', 14.5, col, 700, 'middle')
    for i, (x, y) in enumerate(P):
        if hot == i:
            c.circle(x, y, 20, 'none', VIOLET, 2.2, dash='4 4')
        c.circle(x, y, 12, CARD, GOLD if hot == i else MUTED, 2)
        c.text(x, y + 5, 'ABC'[i], 13, GOLD if hot == i else SOFT, 700, 'middle')
    c.text(cx, cy + hh * 0.38 + 46, title, 15, MUTED, 400, 'middle')
    prod = signs[0] * signs[1] * signs[2]
    c.text(cx, cy + 12, '−1' if prod < 0 else '+1', 22, AMBER, 800, 'middle')


def step3(c, data, y):
    yb = header(c, y, 3, 'Inside the machine', 'Why nothing local can erase it',
                'The twist is read on a boundary. Whatever happens inside the boundary leaves '
                'that reading alone.')
    gap = 24
    w = (CW - 2 * gap) / 3
    h = 462
    ill = 252
    titles = ('Relabeling cannot reach it', 'Repair moves it and never removes it',
              'Only a partner cancels it')
    bodies = (
        'An observer may swap its own names for heads and tails. Both of its seams in the loop '
        'flip together, so the product around the loop stays −1.',
        'Multiply the signs along the edge of a region and you learn whether it holds a twist. '
        'A repair inside the region never touches the edge. The twist wanders, and the edge '
        'reads −1 as before.',
        'Twists are created in pairs and paid for with work. When the two meet again they cancel '
        'and hand the work back. Alone, a twist has nothing to cancel against.')
    for k in range(3):
        x = X0 + k * (w + gap)
        card(c, x, yb, w, h)
        c.rect(x + 14, yb + 14, w - 28, ill, CARD2, None, 0, 14)
        c.text(x + 28, yb + ill + 56, titles[k], 20, TEXT, 700, maxw=w - 50)
        c.para(x + 28, yb + ill + 88, w - 56, bodies[k], 17, MUTED, 1.45)

    # (a) relabeling
    x = X0
    sign_triangle(c, x + 124, yb + 136, 134, (1, 1, -1), 'before')
    sign_triangle(c, x + w - 124, yb + 136, 134, (-1, -1, -1), 'after B relabels', hot=1)
    c.arrow(x + w / 2 - 20, yb + 136, x + w / 2 + 20, yb + 136, FAINT, 2, 8)

    # (b) the boundary reading
    x = X0 + w + gap
    cx, cy = x + w / 2 - 44, yb + 14 + ill / 2
    a = 58
    pts, edges, tris = lattice(cx, cy, a, hexkeep(cx, cy, 2 * a + 1))
    mesh(c, pts, edges, TEAL, 1.2, 0.3, 2.4, 0.6)
    inner = [t for t in tris if all(math.hypot(pts[k][0] - cx, pts[k][1] - cy) < a + 1 for k in t)]
    inner.sort(key=lambda t: math.atan2(_centroid(pts, t)[1] - cy, _centroid(pts, t)[0] - cx))
    ghost = _centroid(pts, inner[3])
    tri_fill(c, pts, inner[3], AMBER, 0.22, stroke=False)
    gx, gy = tri_fill(c, pts, inner[0], AMBER, 0.85, glow=50, cls='pulse')
    c.qarrow(ghost[0], ghost[1] + 4, cx, cy + 30, gx + 2, gy + 6, TEXT, 1.6, 7, 0.9)
    ring = [(cx + a * math.cos(math.radians(60 * i)), cy + a * math.sin(math.radians(60 * i))) for i in range(6)]
    c.poly(ring, 'none', TEXT, 2.8)
    for i in range(6):
        p, q = ring[i], ring[(i + 1) % 6]
        mx, my = (p[0] + q[0]) / 2, (p[1] + q[1]) / 2
        d = math.hypot(mx - cx, my - cy)
        ox, oy = (mx - cx) / d, (my - cy) / d
        neg = i == 4
        c.circle(mx + ox * 17, my + oy * 17, 11, BG, AMBER if neg else TEAL, 1.5)
        c.text(mx + ox * 17, my + oy * 17 + 5.2, '−' if neg else '+', 15.5, AMBER if neg else TEAL, 700, 'middle')
    lx = x + w - 34
    c.text(lx, yb + 100, 'the edge of', 14.5, TEXT, 600, 'end')
    c.text(lx, yb + 119, 'the region', 14.5, TEXT, 600, 'end')
    c.text(lx, yb + 144, 'reads −1', 15.5, AMBER, 700, 'end')
    c.text(lx, yb + 176, 'before and', 14, MUTED, 400, 'end')
    c.text(lx, yb + 194, 'after the move', 14, MUTED, 400, 'end')
    c.line(lx - 86, yb + 110, ring[0][0] + 4, ring[0][1] - 22, FAINT, 1.2)

    # (c) pairs and the energy ledger
    x = X0 + 2 * (w + gap)
    stages = data['path']
    picks = ((stages[0], 'flat'), (stages[1], 'pair created'), (stages[3], 'moved apart'))
    for k, (st, lab) in enumerate(picks):
        cx, cy = x + 84 + k * 142, yb + 86
        carrier(c, data, cx, cy, 30, twisted=tuple(st['defect_faces']), node_r=2.6, show_back=False)
        c.text(cx, yb + 164, lab, 14, SOFT, 600, 'middle')
        bw = 100
        work, cost = st['work_energy_in_J'], st['defect_energy_in_J']
        c.rect(cx - bw / 2, yb + 178, bw, 12, '#1a2540', None, 0, 6)
        if work:
            c.rect(cx - bw / 2, yb + 178, bw * work / 2, 12, BLUE, None, 0, 6)
        if cost:
            c.rect(cx + bw / 2 - bw * cost / 2, yb + 178, bw * cost / 2, 12, AMBER, None, 0, 6)
        c.text(cx, yb + 210, f'work {work}J', 13, BLUE, 600, 'middle')
        c.text(cx, yb + 228, f'twists {cost}J', 13, AMBER, 600, 'middle')
        if k < 2:
            c.arrow(cx + 60, cy, cx + 82, cy, FAINT, 1.8, 7)
    c.caps(x + w / 2, yb + 253, 'energy ledger: 2J at every stage', MUTED, 11.5, 'middle', ls=1.6)
    return status(c, X0, yb + h + 16, CW,
                  'Every outcome K of the repair law commutes with the charge Q: [K, Q] = 0. '
                  'No order of repairs switches the sector. Pair creation, transport and '
                  'annihilation keep the energy ledger at 2J. Here J is a model energy unit.',
                  'Protection covers the stated operations. Resetting the memory that carries '
                  'the charge, or reading it in the wrong basis, removes the charge.',
                  'Scope')


# ---------------------------------------------------------------- step 4
def step4(c, data, y):
    yb = header(c, y, 4, 'Inside the machine', 'Anatomy of a defect',
                'A small core, calm surroundings, and one number that every loop around the core reads.')
    h = 468
    card(c, X0, yb, CW, h)
    cx, cy, R = X0 + 248, yb + h / 2, 200
    c.add(f'<clipPath id="clipV"><circle cx="{n(cx)}" cy="{n(cy)}" r="{R}"/></clipPath>')
    c.circle(cx, cy, R, CARD2)
    c.open(**{'clip-path': 'url(#clipV)'})
    vortex(c, cx, cy, R + 20, 34, length=17, hot=96, cls='pulse')
    c.close()
    c.circle(cx, cy, 112, 'none', TEXT, 2, 0.95, '7 6')
    for k in range(8):
        ang = math.radians(45 * k)
        c.circle(cx + 112 * math.cos(ang), cy + 112 * math.sin(ang), 3.4, TEXT)
    lx = X0 + 492
    for ty, (px, py2), lab, body, col in (
            (yb + 72, (cx + 6, cy - 6), 'Core', 'The mismatch sits here. Repair meets it at every tick.', AMBER),
            (yb + 200, (cx + 112 * math.cos(-0.5), cy + 112 * math.sin(-0.5)), 'Boundary loop',
             'Shrink it or stretch it. While it encloses the core it reads the same full turn.', TEXT),
            (yb + 346, (cx + 172, cy + 84), 'Settled surroundings',
             'Away from the core, neighbors agree almost perfectly.', TEAL)):
        c.line(px, py2, lx - 12, ty - 5, FAINT, 1.2)
        c.circle(px, py2, 3.2, col)
        c.caps(lx, ty, lab, col, 13.5)
        c.para(lx, ty + 26, 236, body, 16.5, MUTED, 1.42)

    # what the number means
    dx = X0 + 792
    c.line(dx - 28, yb + 36, dx - 28, yb + h - 36, EDGE, 1.2)
    rows = (('Charge', 'The number the boundary loop reads. It changes only when something crosses '
                       'the loop, so it is conserved, and the charges of two defects add.'),
            ('Spin', 'Some source channels acquire a minus sign after a full turn. This is a '
                     'spinorial symmetry label. Identifying relativistic spin also requires '
                     'the field transformation law.'),
            ('Antiparticle', 'Charge conjugation pairs opposite source charges with equal energies. '
                             'Their interpretation as particles and antiparticles requires the '
                             'same physical field and interaction records.'))
    ys = (yb + 54, yb + 188, yb + 346)
    for (lab, body), ty in zip(rows, ys):
        c.text(dx + 96, ty + 6, lab, 21, TEXT, 700)
        c.para(dx + 96, ty + 36, 480, body, 17, MUTED, 1.45)
    ix = dx + 38
    c.circle(ix, ys[0] + 34, 27, 'none', TEXT, 1.8, dash='5 5')
    c.glow(ix, ys[0] + 34, 20, AMBER, 0.8)
    c.circle(ix, ys[0] + 34, 5, AMBER)
    c.text(ix, ys[0] + 84, '+1', 15, AMBER, 700, 'middle', family=MONO)
    by = ys[1] + 36
    for s in (-1, 1):
        c.path(f'M{ix - 30},{by} C{ix - 12},{by - 30 * s} {ix + 12},{by + 30 * s} {ix + 30},{by}',
               VIOLET if s > 0 else TEAL, 3)
    c.text(ix, by + 48, '2 turns', 14, VIOLET, 700, 'middle')
    ay = ys[2] + 30
    for s, col, sign in ((-1, AMBER, '+'), (1, BLUE, '−')):
        c.circle(ix + s * 15, ay, 16, mix(CARD, col, 0.18), col, 1.8)
        c.text(ix + s * 15, ay + 5.5, sign, 16, col, 700, 'middle')
    c.text(ix, ay + 44, 'conjugate', 14, MUTED, 600, 'middle')
    return status(c, X0, yb + h + 16, CW,
                  'A declared texture carries cores of winding +1 and −1. Every allowed repair '
                  'outcome preserves both, and the pair supports four localized field modes.',
                  'That the bare repair law forms such a core unaided. The map from the protected '
                  'memory charge of step 3 to that core, and from its exchange sign to the spin '
                  'and electric charge an observer measures.')
