"""Reusable illustrated parts: carrier, lattice patches, loops and fields."""
import math

from .canvas import (AMBER, BG, BLUE, CARD, CARD2, EDGE, FAINT, GOLD, GREEN, MONO, MUTED, ORANGE,
                     RED, SERIF, SOFT, TEAL, TEXT, VIOLET, n, text_width)

W = 1600
X0, X1 = 120, 1520
CW = X1 - X0
SX = 60
ROOT3 = math.sqrt(3)


def mix(a, b, t):
    """Blend two #rrggbb colors."""
    t = max(0.0, min(1.0, t))
    pa = [int(a[i:i + 2], 16) for i in (1, 3, 5)]
    pb = [int(b[i:i + 2], 16) for i in (1, 3, 5)]
    return '#' + ''.join(f'{round(x + (y - x) * t):02x}' for x, y in zip(pa, pb))


def card(c, x, y, w, h, fill=CARD, stroke=EDGE, r=22, op=None):
    c.rect(x, y, w, h, fill, stroke, 1.4, r, op)


def header(c, y, num, act, title, lede=None, color=TEAL, width=1240):
    """Numbered section heading. Registers a node on the thread."""
    c.nodes.append((y + 36, str(num), color))
    c.caps(X0, y + 10, f'{act}  ·  step {num}', color, 14)
    c.text(X0, y + 60, title, 40, TEXT, 700, ls=-0.4, maxw=CW)
    if lede:
        return c.para(X0, y + 100, width, lede, 22, MUTED, 1.4) + 6
    return y + 96


def status(c, x, y, w, computed, derive, label='To derive'):
    """Two-part scope line under a step: what is computed, what is owed."""
    half = (w - 20) / 2
    size, lh = 15.5, 1.42
    tw = half - 62
    # measure line counts with a dry canvas pass
    lines = max(_count(c, computed, tw, size), _count(c, derive, tw, size))
    h = 46 + lines * size * lh
    c.rect(x, y, w, h, CARD2, EDGE, 1.2, 16)
    c.line(x + w / 2, y + 14, x + w / 2, y + h - 14, EDGE, 1.2)
    for i, (lab, body, col) in enumerate((('Computed in the model', computed, TEAL), (label, derive, AMBER))):
        bx = x + 22 + i * (half + 20)
        if i == 0:
            c.circle(bx + 9, y + 24, 6.5, col)
        else:
            c.circle(bx + 9, y + 24, 6.5, 'none', col, 1.8, dash='3.2 2.6')
        c.caps(bx + 26, y + 29, lab, col, 12.5, ls=1.9)
        c.para(bx + 26, y + 29 + size * lh, tw, body, size, MUTED, lh)
    return y + h


def _count(c, markup, w, size):
    from .canvas import Canvas
    probe = Canvas()
    y = probe.para(0, 0, w, markup, size, lh=1.0)
    return round(y / size)


# ---------------------------------------------------------------- carrier
def camera(forward=(0.2, 1.0, 0.3), up=(0.0, 0.0, 1.0)):
    f = _unit(forward)
    r = _unit(_cross(up, f))
    u = _cross(f, r)
    return r, u, f


def _unit(v):
    m = math.sqrt(sum(a * a for a in v))
    return tuple(a / m for a in v)


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def carrier(c, data, cx, cy, scale, twisted=(), cam=None, labels=False, node_r=5.2,
            string=(), dim=1.0, show_back=True):
    """The certified twelve-port carrier, drawn from its exact coordinates.

    Returns projected vertex positions so callers can attach callouts.
    """
    r, u, f = cam or camera()
    verts = [tuple(v['coordinates_display']) for v in data['vertices']]
    proj = [(cx + scale * _dot(p, r), cy - scale * _dot(p, u), _dot(p, f)) for p in verts]
    faces, edges = data['faces'], data['edges']
    depth = [sum(proj[k][2] for k in face) / 3 for face in faces]
    front_edge, front_vertex = set(), set()
    for i, face in enumerate(faces):
        if depth[i] > 0:
            front_vertex.update(face)
            for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
                front_edge.add((min(a, b), max(a, b)))
    if show_back:
        for a, b in edges:
            if (a, b) not in front_edge:
                c.line(proj[a][0], proj[a][1], proj[b][0], proj[b][1], TEAL, 1.1, 0.16 * dim)
        for i, face in enumerate(faces):
            if depth[i] <= 0 and i in twisted:
                c.poly([proj[k][:2] for k in face], AMBER, op=0.16)
    for i in sorted(range(len(faces)), key=lambda j: depth[j]):
        if depth[i] <= 0:
            continue
        pts = [proj[k][:2] for k in faces[i]]
        if i in twisted:
            gx, gy = sum(p[0] for p in pts) / 3, sum(p[1] for p in pts) / 3
            c.glow(gx, gy, scale * 0.95, AMBER, 0.42)
            c.poly(pts, AMBER, AMBER, 1.6, 0.5)
        else:
            shade = 0.05 + 0.07 * depth[i] / 1.5
            c.poly(pts, TEAL, op=round(shade * dim, 3))
    for a, b in edges:
        if (a, b) in front_edge:
            hot = (a, b) in string or (b, a) in string
            c.line(proj[a][0], proj[a][1], proj[b][0], proj[b][1], RED if hot else TEAL,
                   3 if hot else 1.7, (0.95 if hot else 0.7) * dim)
    for i, (x, y, z) in enumerate(proj):
        if i in front_vertex:
            c.circle(x, y, node_r, BG, TEAL, 1.8, dim)
            c.circle(x, y, node_r * 0.38, TEAL, op=dim)
        elif show_back:
            c.circle(x, y, node_r * 0.55, TEAL, op=0.22 * dim)
    return proj, depth


# ---------------------------------------------------------------- lattice
def lattice(cx, cy, a, keep):
    """Triangular lattice clipped by keep(x, y). Returns points, edges, triangles."""
    pts = {}
    span = 14
    for q in range(-span, span + 1):
        for r in range(-span, span + 1):
            x, y = cx + a * (q + r / 2), cy + a * ROOT3 / 2 * r
            if keep(x, y):
                pts[(q, r)] = (x, y)
    edges, tris = [], []
    for (q, r) in pts:
        for dq, dr in ((1, 0), (0, 1), (1, -1)):
            if (q + dq, r + dr) in pts:
                edges.append(((q, r), (q + dq, r + dr)))
        if (q + 1, r) in pts and (q, r + 1) in pts:
            tris.append(((q, r), (q + 1, r), (q, r + 1)))
        if (q + 1, r) in pts and (q + 1, r + 1) in pts and (q, r + 1) in pts:
            tris.append(((q + 1, r), (q + 1, r + 1), (q, r + 1)))
    return pts, edges, tris


def hexkeep(cx, cy, radius):
    def keep(x, y):
        dx, dy = abs(x - cx), abs(y - cy)
        return dy <= radius * ROOT3 / 2 + 0.01 and ROOT3 * dx + dy <= radius * ROOT3 + 0.01
    return keep


def mesh(c, pts, edges, color=TEAL, sw=1.4, op=0.45, node=3.0, node_op=0.9):
    for a, b in edges:
        c.line(pts[a][0], pts[a][1], pts[b][0], pts[b][1], color, sw, op)
    if node:
        for x, y in pts.values():
            c.circle(x, y, node, BG, color, 1.3, node_op)


def tri_fill(c, pts, tri, color=AMBER, op=0.55, glow=0.0, stroke=True, cls=None):
    p = [pts[k] for k in tri]
    gx, gy = sum(a for a, _ in p) / 3, sum(b for _, b in p) / 3
    if glow:
        c.glow(gx, gy, glow, color, 0.5, cls=cls)
    c.poly(p, color, color if stroke else None, 1.6, op, cls=cls)
    return gx, gy


# ---------------------------------------------------------------- loops
def coin(c, x, y, face, name, r=21, color=TEXT):
    """An observer holding a coin: heads is filled, tails is hollow."""
    c.circle(x, y, r + 5, BG, EDGE, 1.2)
    if face == 'H':
        c.circle(x, y, r, GOLD)
        c.text(x, y + 6.5, 'H', 18, BG, 700, 'middle')
    else:
        c.circle(x, y, r, CARD, GOLD, 2.2)
        c.text(x, y + 6.5, 'T', 18, GOLD, 700, 'middle')
    return name


def loop_triangle(c, cx, cy, size, rules, coins, names=('Alice', 'Bob', 'Charlie'), small=False,
                  show_rules=True, product=None):
    """Three observers on a loop. rules[i] joins corner i and i+1: +1 same, -1 different."""
    h = size * ROOT3 / 2
    P = [(cx, cy - h * 0.62), (cx + size / 2, cy + h * 0.38), (cx - size / 2, cy + h * 0.38)]
    bad = []
    for i in range(3):
        a, b = P[i], P[(i + 1) % 3]
        same = coins[i] == coins[(i + 1) % 3]
        ok = same if rules[i] == 1 else not same
        col = GREEN if ok else AMBER
        if not ok:
            bad.append(i)
            c.line(a[0], a[1], b[0], b[1], AMBER, 12 if not small else 8, 0.2)
        c.line(a[0], a[1], b[0], b[1], col, 3.2 if not small else 2.4, None,
               None if ok else '7 6')
        mx, my = (a[0] + b[0]) / 2, (a[1] + b[1]) / 2
        ox, oy = mx - cx, my - (cy - h * 0.62 + h * 2 / 3)
        m = math.hypot(ox, oy) or 1
        ox, oy = ox / m, oy / m
        if show_rules:
            lab = 'SAME' if rules[i] == 1 else 'DIFFERENT'
            tx, ty = mx + ox * (34 if i != 1 else 20), my + oy * 26 + 5
            anchor = 'start' if ox > 0.3 else ('end' if ox < -0.3 else 'middle')
            c.text(tx, ty + (14 if i == 1 else 0), lab, 13, col, 700, anchor, ls=1.4)
        if not small:
            c.circle(mx, my, 12, BG, col, 1.6)
            (c.tick if ok else c.cross)(mx, my, col, 7 if ok else 6, 2)
    for i, (x, y) in enumerate(P):
        coin(c, x, y, coins[i], names[i], 21 if not small else 12)
        if not small:
            dy = -36 if i == 0 else 48
            c.text(x, y + dy, names[i], 15, MUTED, 400, 'middle')
    if product is not None:
        c.text(cx, cy + h * 0.02 + 8, product, 24, GREEN if product == '+1' else AMBER, 700, 'middle')
    return P, bad


# ---------------------------------------------------------------- fields
def vortex(c, cx, cy, radius, step=34, winding=1, phase=math.pi / 2, core=True, length=17,
           hot=70.0, fade=True, cls=None):
    """Arrow field circling a core. Arrows near the core disagree most."""
    k = int(radius // step) + 1
    for i in range(-k, k + 1):
        for j in range(-k, k + 1):
            x, y = cx + i * step, cy + j * step
            d = math.hypot(x - cx, y - cy)
            if d > radius or d < step * 0.55:
                continue
            ang = winding * math.atan2(y - cy, x - cx) + phase
            col = mix(AMBER, TEAL, (d - step * 0.6) / hot)
            op = 1.0 if not fade else max(0.25, 1 - (d / radius) ** 2.4 * 0.8)
            c.dart(x, y, ang, length, col, 1.7, round(op, 2))
    if core:
        c.glow(cx, cy, step * 1.25, AMBER, 0.75, cls=cls)
        c.circle(cx, cy, 6.5, AMBER)


def wave(c, x, y, w, amp, cycles, color, sw=2.4, envelope=None, op=None, phase=0.0, steps=140):
    pts = []
    for i in range(steps + 1):
        t = i / steps
        e = 1.0 if envelope is None else envelope(t)
        pts.append((x + w * t, y - amp * e * math.sin(2 * math.pi * cycles * t + phase)))
    c.curve(pts, color, sw, op)


def bell(t, mu=0.5, sigma=0.16):
    return math.exp(-((t - mu) / sigma) ** 2 / 2)


def lcg(seed):
    """Tiny deterministic generator for dot scatters."""
    state = seed & 0x7fffffff
    while True:
        state = (1103515245 * state + 12345) & 0x7fffffff
        yield state / 0x7fffffff


def mobius(c, cx, cy, R, half, segs=120, elev=math.radians(34), spin=math.pi / 2):
    """A Mobius band from its parametrization: shaded quads in depth order."""
    se, ce = math.sin(elev), math.cos(elev)

    def space(u, v):
        a = u + spin
        r = R + v * math.cos(u / 2)
        return (r * math.cos(a), r * math.sin(a), v * math.sin(u / 2))

    def screen(p):
        return (cx + p[0], cy - (p[1] * se + p[2] * ce), -p[1] * ce + p[2] * se)

    light = _unit((-0.35, -0.55, 0.76))
    quads = []
    strips = 3
    for i in range(segs):
        u0, u1 = 2 * math.pi * i / segs, 2 * math.pi * (i + 1) / segs
        um = (u0 + u1) / 2
        for j in range(strips):
            v0, v1 = -half + 2 * half * j / strips, -half + 2 * half * (j + 1) / strips
            corners = [space(u0, v0), space(u0, v1), space(u1, v1), space(u1, v0)]
            du = tuple(b - a for a, b in zip(space(um - 0.01, (v0 + v1) / 2), space(um + 0.01, (v0 + v1) / 2)))
            dv = tuple(b - a for a, b in zip(space(um, v0), space(um, v1)))
            normal = _unit(_cross(du, dv))
            shade = abs(_dot(normal, light))
            pts = [screen(p) for p in corners]
            z = sum(p[2] for p in pts) / 4
            rim = []
            if j == 0:
                rim.append((pts[0], pts[3]))
            if j == strips - 1:
                rim.append((pts[1], pts[2]))
            quads.append((z, [p[:2] for p in pts], shade, i / segs, rim))
    for z, pts, shade, t, rim in sorted(quads, key=lambda q: q[0]):
        base = mix(TEAL, VIOLET, 0.5 - 0.5 * math.cos(2 * math.pi * t))
        col = mix('#0b1524', base, 0.22 + 0.74 * shade)
        c.poly(pts, col, col, 0.9)
        for a, b in rim:
            c.line(a[0], a[1], b[0], b[1], mix(base, '#ffffff', 0.45), 1.5, 0.95)
