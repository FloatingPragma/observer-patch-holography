"""Draws executable_spec.svg: the hardware and the software of A-contrast-memory-spherical-v1 on one sheet.

Everything geometric is read from carrier files: the wiring of one cell from inputs/CELL.npz, the federations
from ../contrast_memory_spherical_v1/carriers/FEDERATION_LEVEL{0,1,2}.npz. The formulas are those of CONSTRUCTION_MANUAL.md.

    python make_diagram.py        # writes executable_spec.svg; needs numpy and networkx
"""
import collections
import itertools
import json
import math
from pathlib import Path

import networkx as nx
import numpy as np

HERE = Path(__file__).resolve().parent
INK, MUTED, FAINT, PANEL = '#1f2328', '#59636e', '#d0d7de', '#f6f8fa'
BLUE, ORANGE, PURPLE, GREEN, AMBER, RED = '#0969da', '#bc4c00', '#8250df', '#1a7f37', '#9a6700', '#cf222e'
SANS = 'Helvetica, Arial, sans-serif'
MONO = 'Menlo, Consolas, DejaVu Sans Mono, monospace'
out = []


def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def text(x, y, s, size=17, fill=INK, anchor='start', weight='normal', family=SANS, style=''):
    out.append('<text xml:space="preserve" x="%.1f" y="%.1f" font-family="%s" font-size="%d" fill="%s" text-anchor="%s" font-weight="%s"%s>%s</text>' % (
        x, y, family, size, fill, anchor, weight, (' font-style="%s"' % style) if style else '', esc(s)))


def mono(x, y, s, size=16, fill=INK, weight='normal', anchor='start'):
    text(x, y, s, size, fill, anchor, weight, MONO)


def lines(x, y, rows, size=17, gap=None, fill=INK, family=SANS, weight='normal'):
    gap = gap or size * 1.45
    for k, r in enumerate(rows):
        text(x, y + k * gap, r, size, fill, 'start', weight, family)
    return y + len(rows) * gap


def line(x1, y1, x2, y2, stroke=MUTED, w=1.5, dash=None, opacity=1.0, cap='round'):
    out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="%.2f" stroke-linecap="%s"%s%s/>' % (
        x1, y1, x2, y2, stroke, w, cap, (' stroke-dasharray="%s"' % dash) if dash else '', (' opacity="%.2f"' % opacity) if opacity < 1 else ''))


def circle(x, y, r, fill='#ffffff', stroke=INK, w=1.5, opacity=1.0):
    out.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" stroke="%s" stroke-width="%.2f"%s/>' % (
        x, y, r, fill, stroke, w, (' opacity="%.2f"' % opacity) if opacity < 1 else ''))


def rect(x, y, w, h, fill='none', stroke=FAINT, sw=1.5, rx=10, opacity=1.0):
    out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="%.1f" fill="%s" stroke="%s" stroke-width="%.2f"%s/>' % (
        x, y, w, h, rx, fill, stroke, sw, (' opacity="%.2f"' % opacity) if opacity < 1 else ''))


def poly(points, fill, stroke='none', w=1.0, opacity=1.0):
    out.append('<polygon points="%s" fill="%s" stroke="%s" stroke-width="%.2f" opacity="%.2f"/>' % (
        ' '.join('%.1f,%.1f' % p for p in points), fill, stroke, w, opacity))


def path(d, stroke=INK, w=1.5, fill='none', dash=None, marker=None):
    out.append('<path d="%s" fill="%s" stroke="%s" stroke-width="%.2f" stroke-linecap="round"%s%s/>' % (
        d, fill, stroke, w, (' stroke-dasharray="%s"' % dash) if dash else '', (' marker-end="url(#%s)"' % marker) if marker else ''))


def arrow(x1, y1, x2, y2, stroke=INK, w=1.8, both=False):
    out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="%.2f" marker-end="url(#ah-%s)"%s/>' % (
        x1, y1, x2, y2, stroke, w, stroke[1:], (' marker-start="url(#as-%s)"' % stroke[1:]) if both else ''))


PANELS = []


def panel(x, y, w, h, tag, title, name=None):
    PANELS.append([name or tag.lower().replace(' ', '_'), x, y, w, h, len(out)])
    rect(x, y, w, h, '#ffffff', FAINT, 1.5, 14)
    out.append('<path d="M%.1f,%.1f h%.1f a14,14 0 0 1 14,14 v34 h-%.1f v-34 a14,14 0 0 1 14,-14 z" fill="%s" stroke="%s" stroke-width="1.5"/>' % (
        x + 14, y, w - 28, w, PANEL, FAINT))
    text(x + 22, y + 32, tag, 15, MUTED, weight='bold')
    text(x + 22 + 9.6 * len(tag) + 14, y + 32, title, 21, INK, weight='bold')


# ---------------------------------------------------------------- geometry

def canonical_icosahedron():
    p = (1 + 5 ** 0.5) / 2
    V = []
    for a, b in itertools.product((-1, 1), (-p, p)):
        V += [(0, a, b), (a, b, 0), (b, 0, a)]
    V = np.array(V, float)
    V /= np.linalg.norm(V, axis=1)[:, None]
    d = np.linalg.norm(V[:, None] - V[None], axis=2)
    adj = np.abs(d - d[d > 1e-9].min()) < 1e-6
    return V, [(i, j) for i in range(12) for j in range(i + 1, 12) if adj[i, j]]


def cell():
    z = np.load(HERE / 'inputs' / 'CELL.npz', allow_pickle=False)
    seams = [tuple(map(int, e)) for e in z['edges']]
    V, E = canonical_icosahedron()
    iso = nx.vf2pp_isomorphism(nx.Graph(seams), nx.Graph(E))          # port label -> canonical vertex
    X = np.array([V[iso[p]] for p in range(12)])
    return seams, X


def rot(ax, deg):
    a = math.radians(deg)
    c, s = math.cos(a), math.sin(a)
    return {'x': np.array([[1, 0, 0], [0, c, -s], [0, s, c]]), 'y': np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]]),
            'z': np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])}[ax]


def align(v, target):
    """a rotation that takes the unit vector v to the unit vector target"""
    v, t = v / np.linalg.norm(v), np.asarray(target, float) / np.linalg.norm(target)
    c = float(v @ t)
    if c > 1 - 1e-12:
        return np.eye(3)
    if c < -1 + 1e-12:
        a = np.cross(v, [1.0, 0.0, 0.0] if abs(v[0]) < 0.9 else [0.0, 1.0, 0.0])      # a half turn about an axis across v
        a /= np.linalg.norm(a)
        return 2 * np.outer(a, a) - np.eye(3)
    k = np.cross(v, t)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + K + K @ K * (1 / (1 + c))


def draw_cell(X, seams, cx, cy, R, M, labels=True, star=None, face=None, wedge=None, ext=(), pr=17, fs=15, tag=None, signed=None):
    P = X @ M.T
    xy = lambda p: (cx + R * P[p, 0], cy - R * P[p, 1])
    order = sorted(range(len(seams)), key=lambda k: P[seams[k][0], 2] + P[seams[k][1], 2])
    signed = signed or {}
    if face:
        poly([xy(p) for p in face], '#ddf4ff', BLUE, 1.0, 0.85)
    for k in order:
        a, b = seams[k]
        back = (P[a, 2] + P[b, 2]) / 2 < -0.05
        hot = star is not None and star in (a, b)
        if (a, b) in signed:
            line(*xy(a), *xy(b), stroke=BLUE if signed[(a, b)] > 0 else RED, w=6.5)
            continue
        line(*xy(a), *xy(b), stroke=BLUE if hot else '#6e7781', w=3.6 if hot else (1.3 if back else 2.2),
             dash='5 6' if back and not hot else None, opacity=0.55 if back and not hot else 1.0)
    if wedge:
        a, b, c = wedge
        pa, pb, pc = np.array(xy(a)), np.array(xy(b)), np.array(xy(c))
        u, v = (pa - pb) / np.linalg.norm(pa - pb), (pc - pb) / np.linalg.norm(pc - pb)
        r = 0.34 * R
        s, e = pb + r * u, pb + r * v
        path('M%.1f,%.1f Q%.1f,%.1f %.1f,%.1f' % (s[0], s[1], *(pb + 1.25 * r * (u + v) / max(np.linalg.norm(u + v), 0.3)), e[0], e[1]), PURPLE, 3.2)
        circle(*(pb + 0.92 * r * (u + v) / max(np.linalg.norm(u + v), 0.3)), 6, PURPLE, PURPLE)
    for p in sorted(range(12), key=lambda q: P[q, 2]):
        back = P[p, 2] < -0.05
        is_ext = p in ext
        circle(*xy(p), pr, '#fff8c5' if is_ext else ('#f6f8fa' if back else '#ffffff'), ORANGE if is_ext else (MUTED if back else BLUE), 2.2 if not back else 1.4)
        if labels:
            text(xy(p)[0], xy(p)[1] + fs * 0.36, str(p), fs, MUTED if back and not is_ext else INK, 'middle', 'bold')
    if tag:
        text(cx, cy + R + 44, tag, 16, MUTED, 'middle')
    return xy


def federation(level):
    z = np.load(HERE.parent / 'contrast_memory_spherical_v1' / 'carriers' / f'FEDERATION_LEVEL{level}.npz', allow_pickle=False)
    edges = sorted({tuple(sorted(map(int, e))) for e in z['edges']})
    cells = len(z['degree']) // 12
    ext = [(a, b) for a, b in edges if a // 12 != b // 12]
    seams0 = [(a, b) for a, b in edges if a < 12 and b < 12]
    D = dict(nx.all_pairs_shortest_path_length(nx.Graph(seams0)))
    att = collections.defaultdict(list)
    for a, b in ext:
        att[a // 12].append(a % 12)
        att[b // 12].append(b % 12)
    sym = {c: sorted(D[x][y] for x, y in itertools.combinations(att[c], 2)) == [2, 2, 2] for c in range(cells)}
    # positions: the triangles of the geodesic tower, matched to the carrier's graph of cells
    V, E = canonical_icosahedron()
    F = [f for f in itertools.combinations(range(12), 3) if all(tuple(sorted(p)) in set(E) for p in itertools.combinations(f, 2))]
    V = [tuple(v) for v in V]
    for _ in range(level):
        mid, F2 = {}, []

        def m(a, b):
            k = (min(a, b), max(a, b))
            if k not in mid:
                x = (np.array(V[a]) + np.array(V[b])) / 2
                mid[k] = len(V)
                V.append(tuple(x / np.linalg.norm(x)))
            return mid[k]
        for a, b, c in F:
            ab, bc, ca = m(a, b), m(b, c), m(c, a)
            F2 += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        F = F2
    V = np.array(V)
    owner, dual = {}, []
    for f, tri in enumerate(F):
        for x, y in itertools.combinations(tri, 2):
            k = (min(x, y), max(x, y))
            if k in owner:
                dual.append((owner[k], f))
            else:
                owner[k] = f
    iso = nx.vf2pp_isomorphism(nx.Graph((a // 12, b // 12) for a, b in ext), nx.Graph(dual))
    cent = np.array([V[list(F[iso[c]])].sum(axis=0) for c in range(cells)])
    cent /= np.linalg.norm(cent, axis=1)[:, None]
    return cells, len(edges), len(ext), [(a // 12, b // 12) for a, b in ext], cent, sym, att, ext


def draw_sphere(level, cx, cy, R, M):
    cells, n_seams, n_ext, pairs, cent, sym, _, _ = federation(level)
    P = cent @ M.T
    circle(cx, cy, R * 1.035, '#fbfcfd', FAINT, 1.2)
    for a, b in sorted(pairs, key=lambda e: P[e[0], 2] + P[e[1], 2]):
        back = (P[a, 2] + P[b, 2]) / 2 < 0
        line(cx + R * P[a, 0], cy - R * P[a, 1], cx + R * P[b, 0], cy - R * P[b, 1], ORANGE, 0.8 if back else (2.0 if level < 2 else 1.3),
             opacity=0.16 if back else 0.95)
    r = {0: 11, 1: 6.5, 2: 3.6}[level]
    for c in sorted(range(cells), key=lambda q: P[q, 2]):
        if P[c, 2] < 0:
            continue
        circle(cx + R * P[c, 0], cy - R * P[c, 1], r, BLUE if sym[c] else '#ffffff', BLUE, 1.6 if level < 2 else 1.1)
    n_sym = sum(sym.values())
    text(cx, cy + R + 40, '%d cells' % cells, 20, INK, 'middle', 'bold')
    text(cx, cy + R + 64, '%d ports, %d seams' % (12 * cells, n_seams), 16, MUTED, 'middle')
    text(cx, cy + R + 86, '%d of them join two cells' % n_ext, 16, MUTED, 'middle')
    text(cx, cy + R + 108, 'junctions: %d symmetric, %d other' % (n_sym, cells - n_sym), 16, MUTED, 'middle')


# ---------------------------------------------------------------- the sheet

def idea_strip(x0, y0, w, h):
    panel(x0, y0, w, h, 'THE IDEA', 'A Universe that computes its own history', 'idea')
    cx, cy, r = x0 + 232, y0 + 200, 78
    out.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="none" stroke="%s" stroke-width="2.4" stroke-dasharray="2 9"/>' % (cx, cy, r, MUTED))
    nodes = [(-90, ['observers hold records'], BLUE), (0, ['they recover', 'the rules'], GREEN), (90, ['they build the machine'], AMBER), (180, ['its run is', 'this history'], PURPLE)]
    for ang, lab, col in nodes:
        a = math.radians(ang)
        px, py = cx + r * math.cos(a), cy + r * math.sin(a)
        circle(px, py, 11, col, col)
        a2 = math.radians(ang + 38)
        a3 = math.radians(ang + 52)
        arrow(cx + r * math.cos(a2), cy + r * math.sin(a2), cx + r * math.cos(a3), cy + r * math.sin(a3), MUTED, 2.2)
        if ang == -90:
            text(px, py - 24, lab[0], 15, col, 'middle', 'bold')
        elif ang == 90:
            text(px, py + 36, lab[0], 15, col, 'middle', 'bold')
        else:
            for k, part in enumerate(lab):
                text(px + (20 if ang == 0 else -20), py - 4 + 19 * k, part, 15, col, 'start' if ang == 0 else 'end', 'bold')
    text(cx, cy + 5, 'no outside', 15, MUTED, 'middle')
    c1, c2, c3 = x0 + 500, x0 + 1120, x0 + 1800
    text(c1, y0 + 92, 'A brain that thinks itself into existence', 19, INK, weight='bold')
    lines(c1, y0 + 122, ['The Universe is a network of very simple observers. Each holds',
                         'a few records and compares them with its neighbours. What',
                         'survives every comparison is a public fact. The observers inside',
                         'work out the rules and build the machine whose run is this same',
                         'history, its builders included. Nothing stands outside the loop.'], 16, 24, INK)
    text(c1, y0 + 262, 'closure:  R(D(U)) ≅ U', 16, MUTED, family=MONO)
    text(c2, y0 + 92, 'Why this problem can be solved', 19, INK, weight='bold')
    lines(c2, y0 + 122, ['A consistent history is a fixed point: every record agrees with every',
                         'record it overlaps. A search through all histories grows exponentially.',
                         'The observers do not search. Each repairs what it shares with one',
                         'neighbour, and every accepted repair lowers the total disagreement by',
                         'a fixed step. The work is bounded by the disagreement at the start:',
                         'polynomial in the number of ports.'], 16, 24, INK)
    lines(c2, y0 + 276, ['Proved for the native repair of the theory, with loads bounded by the alphabet.',
                         'For the two candidate laws tested here it is the target, and work in progress.'], 14, 19, MUTED)
    text(c3, y0 + 92, 'How large it is', 19, INK, weight='bold')
    lines(c3, y0 + 122, ['The carriers of the tests: 20, 80 and 320 cells.',
                         'The screen of the Universe we inhabit: about',
                         '8 × 10¹²² cells in the count of the theory.',
                         'One kind of cell, one wiring, one law on every seam,',
                         'repeated that many times.'], 16, 24, INK)


def main():
    Wd = 2400
    out.append('<<HEADER>>')
    text(40, 56, 'The executable construction of A-contrast-memory-spherical-v1: hardware and software', 34, INK, weight='bold')
    text(40, 88, 'Hardware: cells of twelve ports, their seams, the sphere of cells and the registers.  Software: two laws under test, a cage of the second, the instruments of an observer, its clock and its records.', 18, MUTED)
    idea_strip(40, 112, 2320, 330)

    seams, X = cell()
    nb = {p: sorted([b if a == p else a for a, b in seams if p in (a, b)]) for p in range(12)}

    # ---------------- panel 1: the cell
    x0, y0, w, h = 40, 472, 1140, 900
    panel(x0, y0, w, h, 'HARDWARE 1', 'The cell: an echosahedron', 'cell')
    M = rot('x', -62) @ rot('z', 18) @ align(X[0], (0, 0, 1))
    xy = draw_cell(X, seams, x0 + 330, y0 + 480, 270, M, star=0, face=(2, 4, 8), wedge=(1, 0, 2))
    px, py = xy(0)
    text(px + 26, py - 20, 'the star of port 0: its five seams', 16, BLUE, weight='bold')
    p0, p3, p4 = np.array(xy(0)), np.array(xy(1)), np.array(xy(2))
    u, v = (p3 - p0) / np.linalg.norm(p3 - p0), (p4 - p0) / np.linalg.norm(p4 - p0)
    dot = p0 + 0.92 * 0.34 * 270 * (u + v) / max(np.linalg.norm(u + v), 0.3)
    text(x0 + 36, y0 + 100, 'a wedge: two seams at one port', 16, PURPLE, weight='bold')
    text(x0 + 36, y0 + 122, 'here (1, 0, 2); ten at every port', 16, PURPLE)
    path('M%.1f,%.1f L%.1f,%.1f' % (x0 + 150, y0 + 134, dot[0] - 8, dot[1] - 8), PURPLE, 1.2, dash='3 4')
    fx = np.mean([xy(p)[0] for p in (2, 4, 8)])
    fy = np.mean([xy(p)[1] for p in (2, 4, 8)])
    text(fx, fy + 6, 'a face', 16, BLUE, 'middle', 'bold')
    tx = x0 + 680
    rows = [('12 ports', 'the corners; each is a small observer'), ('30 seams', 'the edges; two ports compare here'),
            ('20 faces', 'triangles of three seams'), ('120 wedges', 'pairs of seams at a port, ten each'),
            ('antipode', 'port p faces port 11 − p'), ('symmetry', '60 rotations; 120 with the inversion')]
    for k, (a, b) in enumerate(rows):
        mono(tx, y0 + 100 + 27 * k, a, 16, INK, 'bold')
        text(tx + 122, y0 + 100 + 27 * k, b, 16, INK)
    ty = y0 + 100 + 27 * len(rows) + 22
    text(tx, ty, 'Wiring: each port and its five neighbours', 17, INK, weight='bold')
    ty = lines(tx, ty + 30, ['%2d :  %s' % (p, '  '.join('%2d' % q for q in nb[p])) for p in range(12)], 16, 24, INK, MONO)
    text(tx, ty + 18, 'Why this shape', 17, INK, weight='bold')
    lines(tx, ty + 44, ['Twelve ports on an icosahedron: one kind of port, one',
                        'kind of seam, the largest symmetry twelve points allow.',
                        'Its six axes span three directions, and the slowest',
                        'response of the cell has rank three.'], 16, 23, MUTED)
    lines(x0 + 36, y0 + h - 62, ['Blue lines: the five seams of port 0. Dashed lines lie on the far side.',
                                 'The labels 0 to 11 are those of the carrier files; every table of the manual uses them.'], 16, 23, MUTED)

    # ---------------- panel 2: registers
    x0, y0, w, h = 1220, 472, 1140, 900
    panel(x0, y0, w, h, 'HARDWARE 2', 'What sits on a port, a seam and a wedge', 'registers')
    text(x0 + 30, y0 + 92, 'Law 1: loads with a seam memory', 19, INK, weight='bold')
    ax, bx, cy1 = x0 + 150, x0 + 590, y0 + 262
    line(ax, cy1, bx, cy1, '#6e7781', 5)
    for px_, lab in ((ax, 'a'), (bx, 'b')):
        circle(px_, cy1, 44, '#ffffff', BLUE, 3)
        text(px_, cy1 + 8, lab, 24, INK, 'middle', 'bold')
        rect(px_ - 92, cy1 - 126, 184, 54, '#ddf4ff', BLUE, 1.3, 8)
        text(px_, cy1 - 103, 'port register', 15, INK, 'middle', 'bold')
        text(px_, cy1 - 83, 'load n = 0, 1 or 2', 15, INK, 'middle')
        line(px_, cy1 - 72, px_, cy1 - 44, BLUE, 1.3, dash='3 3')
    mx = (ax + bx) / 2
    text(mx, cy1 - 14, 'seam e = {a, b}', 16, MUTED, 'middle')
    gx, gy, cw, ch = mx - 100, cy1 + 44, 50, 34
    rect(gx - 82, gy + 17, 64, ch - 2, '#ffffff', GREEN, 1.6, 6)
    text(gx - 50, gy + 39, 'blank', 14, GREEN, 'middle', 'bold')
    for r_ in range(2):
        for k, xval in enumerate(range(-2, 3)):
            rect(gx + k * cw, gy + r_ * ch, cw - 4, ch - 4, '#dafbe1', GREEN, 1.2, 5)
            text(gx + k * cw + (cw - 4) / 2, gy + r_ * ch + 21, ('%+d' % xval) if xval else '0', 14, INK, 'middle')
        text(gx + 5 * cw + 8, gy + r_ * ch + 21, 'c = %d' % r_, 14, MUTED)
    line(mx, cy1 + 4, mx, gy - 2, GREEN, 1.3, dash='3 3')
    text(x0 + 30, gy + 2 * ch + 28, 'private memory of the seam: blank, or one of ten records (c, x)', 16, GREEN, weight='bold')
    text(x0 + 30, gy + 2 * ch + 51, 'x: the difference of the two loads that was written.  c: the coin of an odd total.', 15, MUTED)
    text(x0 + 800, cy1 - 110, 'also on the hardware', 16, INK, weight='bold')
    lines(x0 + 800, cy1 - 84, ['a reference register at zero on', 'every port with five seams, so that', 'every port counts six;', 'an emitter with its fuel;', 'a reader of eleven bits per seam'], 15, 22, MUTED)

    text(x0 + 30, y0 + 500, 'Law 2: one mode on every seam and on every wedge', 19, INK, weight='bold')
    pa, pb, pc = (x0 + 120, y0 + 800), (x0 + 400, y0 + 600), (x0 + 680, y0 + 800)
    line(*pa, *pb, '#6e7781', 5)
    line(*pb, *pc, '#6e7781', 5)
    for (px_, py_), lab in ((pa, 'a'), (pb, 'b'), (pc, 'c')):
        circle(px_, py_, 40, '#ffffff', BLUE, 3)
        text(px_, py_ + 8, lab, 22, INK, 'middle', 'bold')
    s1 = ((pa[0] + pb[0]) / 2, (pa[1] + pb[1]) / 2)
    s2 = ((pb[0] + pc[0]) / 2, (pb[1] + pc[1]) / 2)
    for sx, sy in (s1, s2):
        rect(sx - 19, sy - 19, 38, 38, '#ffffff', INK, 2.4, 5)
        text(sx, sy + 6, '+m', 15, INK, 'middle', 'bold')
    text(s1[0] - 34, s1[1] - 30, 'seam mode (a, b)', 15, INK, 'end', 'bold')
    text(s2[0] + 34, s2[1] - 30, 'seam mode (b, c)', 15, INK, 'start', 'bold')
    wx_, wy_ = pb[0], pb[1] + 150
    poly([(wx_, wy_ - 28), (wx_ + 28, wy_), (wx_, wy_ + 28), (wx_ - 28, wy_)], '#fbefff', PURPLE, 2.4, 1.0)
    text(wx_, wy_ + 6, '−m', 15, PURPLE, 'middle', 'bold')
    arrow(s1[0] + 24, s1[1] + 8, wx_ - 34, wy_ - 6, PURPLE, 2, True)
    arrow(s2[0] - 24, s2[1] + 8, wx_ + 34, wy_ - 6, PURPLE, 2, True)
    mono((s1[0] + wx_) / 2 - 6, (s1[1] + wy_) / 2 - 16, '+1/√2', 15, PURPLE, 'bold', 'middle')
    mono((s2[0] + wx_) / 2 + 6, (s2[1] + wy_) / 2 - 16, '−1/√2', 15, PURPLE, 'bold', 'middle')
    text(wx_, wy_ + 58, 'wedge mode (a, b, c):', 15, PURPLE, 'middle', 'bold')
    text(wx_, wy_ + 78, 'the object in transit at port b', 15, PURPLE, 'middle')
    fx_, fy_ = x0 + 600, y0 + 556
    rect(fx_, fy_, 250, 100, '#fff8c5', AMBER, 1.4, 8)
    text(fx_ + 14, fy_ + 24, 'apparatus at every port', 15, INK, weight='bold')
    lines(fx_ + 14, fy_ + 46, ['a fuel mode at level 0', 'a parity reference', 'two pointers, a record memory'], 15, 20, INK)
    line(pb[0] + 40, pb[1] - 6, fx_, fy_ + 40, AMBER, 1.6, dash='4 4')
    text(x0 + 890, y0 + 580, 'one cell', 16, INK, weight='bold')
    lines(x0 + 890, y0 + 606, [' 30 seam modes', '120 wedge modes', ' 12 fuels', ' 12 references', '174 modes, written', 'on 432 two-state', 'registers'], 16, 23, INK, MONO)
    text(x0 + 30, y0 + h - 24, 'Each mode holds one object at most. An object on a seam moves through a wedge to the next seam of the same port.', 15, MUTED)

    # ---------------- panel 3: federation
    x0, y0, w, h = 40, 1412, 2320, 900
    panel(x0, y0, w, h, 'HARDWARE 3', 'The federation: cells on a sphere, one seam between neighbouring cells', 'federation')
    Ms = rot('x', -68) @ rot('z', 12)
    for level, cx in ((0, x0 + 250), (1, x0 + 740), (2, x0 + 1230)):
        draw_sphere(level, cx, y0 + 330, 205, Ms)
    text(x0 + 40, y0 + 690, 'Every dot is one whole cell of twelve ports. Every orange line is one external seam: it joins one port of a cell to one port of a neighbour.', 17, INK)
    text(x0 + 40, y0 + 716, 'A cell sits on a triangle of the sphere and has three neighbours. Each refinement cuts every triangle into four: 20, 80, 320 cells.', 17, INK)
    text(x0 + 40, y0 + 742, 'Filled dots: the three external ports of the cell lie symmetrically (pairwise two seams apart). Open dots: the gluing rule chose an unsymmetric triple.', 17, MUTED)
    lines(x0 + 40, y0 + 790, ['Gluing rule of the builder (a declared convention): the external seam towards a neighbour leaves through the port whose direction, in the',
                              'local frame of the cell, points most nearly at that neighbour; two seams never share a port. A port therefore has five or six seams.'], 17, 26, INK)
    lines(x0 + 40, y0 + 856, ['Why a sphere: the theory places observers on a closed screen. Why refinement: a result has to hold as the cells get finer.'], 17, 26, MUTED)
    cells, _, _, pairs, cent, sym, att, ext = federation(0)
    a_port, b_port = next((a, b) for a, b in ext if {a // 12, b // 12} == {0, 4})
    ca, cb = (a_port // 12, b_port // 12)
    pa_, pb_ = a_port % 12, b_port % 12
    zx, zy, zr = x0 + 1640, y0 + 310, 150
    Ma = rot('x', -22) @ rot('y', 10) @ align(X[pa_], (1, 0, 0))
    Mb = rot('x', -22) @ rot('y', -10) @ align(X[pb_], (-1, 0, 0))
    xa = draw_cell(X, seams, zx, zy, zr, Ma, ext=att[ca], pr=13, fs=12)
    xb = draw_cell(X, seams, zx + 430, zy, zr, Mb, ext=att[cb], pr=13, fs=12)
    line(*xa(pa_), *xb(pb_), ORANGE, 5)
    for xf, c_, own, cx_ in ((xa, ca, pa_, zx), (xb, cb, pb_, zx + 430)):
        for p in att[c_]:
            if p != own:
                px_, py_ = xf(p)
                dx, dy = px_ - cx_, py_ - zy
                n = math.hypot(dx, dy) or 1
                line(px_ + 14 * dx / n, py_ + 14 * dy / n, px_ + 66 * dx / n, py_ + 66 * dy / n, ORANGE, 4, dash='2 7')
    for xf, p in ((xa, pa_), (xb, pb_)):
        circle(*xf(p), 13, '#fff8c5', ORANGE, 2.4)
        text(xf(p)[0], xf(p)[1] + 4.5, str(p), 12, INK, 'middle', 'bold')
    text(zx, zy + zr + 46, 'cell %d' % ca, 18, INK, 'middle', 'bold')
    text(zx + 430, zy + zr + 46, 'cell %d' % cb, 18, INK, 'middle', 'bold')
    text(zx + 215, zy - zr - 34, 'The wiring between two cells of the 20-cell carrier', 18, INK, 'middle', 'bold')
    text(zx + 215, zy + zr + 80, 'one external seam: port %d of cell %d to port %d of cell %d' % (pa_, ca, pb_, cb), 16, ORANGE, 'middle', 'bold')
    text(zx + 215, zy + zr + 104, 'yellow ports carry the external seams of each cell; dotted stubs lead to the other neighbours', 15, MUTED, 'middle')
    mono(x0 + 1470, y0 + 620, 'cells   ports    seams = inside + between    wedges', 16, INK, 'bold')
    for k, level in enumerate((0, 1, 2)):
        z = np.load(HERE.parent / 'contrast_memory_spherical_v1' / 'carriers' / f'FEDERATION_LEVEL{level}.npz', allow_pickle=False)
        deg = z['degree']
        n = len(deg) // 12
        ns = len({tuple(sorted(map(int, e))) for e in z['edges']})
        mono(x0 + 1470, y0 + 648 + 26 * k, '%5d  %6d   %6d = %6d + %7d   %7d' % (n, 12 * n, ns, 30 * n, ns - 30 * n, int(sum(d * (d - 1) // 2 for d in deg))), 16)

    # ---------------- panel 4: law 1
    x0, y0, w, h = 40, 2352, 1140, 1080
    panel(x0, y0, w, h, 'SOFTWARE 1', 'Law 1: loads move, neighbours pay for disagreement', 'law1')
    y = y0 + 92
    mono(x0 + 30, y, 'H = ½ Σ_seams (1 − SWAP_e)  +  κ C', 20, INK, 'bold')
    mono(x0 + 30, y + 32, 'C = Σ_seams (n_a − n_b)² + Σ_ports (6 − d_v) n_v²  =  6 Σ_ports n_v² − 2 Σ_seams n_a n_b', 16)
    text(x0 + 30, y + 66, 'exchange: a load crosses its seam', 16, BLUE, weight='bold')
    text(x0 + 330, y + 66, 'agreement energy: a seam costs κ times the squared difference of its two loads', 16, RED, weight='bold')
    ex, ey = x0 + 70, y + 124
    for k, filled in enumerate((True, False)):
        circle(ex + 110 * k, ey, 24, '#ffffff', BLUE, 2.4)
        if filled:
            circle(ex + 110 * k, ey, 10, INK, INK)
    line(ex + 24, ey, ex + 86, ey, '#6e7781', 3)
    arrow(ex + 34, ey - 34, ex + 78, ey - 34, BLUE, 2, True)
    text(ex + 55, ey + 52, 'exchange', 15, MUTED, 'middle')
    cx0 = x0 + 350
    for k, (lab, cost, conf) in enumerate((('two loads on one seam', '10 κ', (1, 1, 0, 0)), ('two loads apart', '12 κ', (1, 0, 0, 1)), ('two loads on one port', '24 κ', (2, 0, 0, 0)))):
        bx_ = cx0 + 262 * k
        for j in range(4):
            if j < 3:
                line(bx_ + 46 * j + 16, ey, bx_ + 46 * (j + 1) - 16, ey, '#6e7781', 2.4)
            circle(bx_ + 46 * j, ey, 16, '#ffffff', BLUE, 2)
            if conf[j] == 1:
                circle(bx_ + 46 * j, ey, 7, INK, INK)
            if conf[j] == 2:
                circle(bx_ + 46 * j - 6, ey, 5.5, INK, INK)
                circle(bx_ + 46 * j + 6, ey, 5.5, INK, INK)
        text(bx_ + 69, ey + 46, lab, 15, MUTED, 'middle')
        mono(bx_ + 69, ey + 70, cost, 17, RED, 'bold', 'middle')
    text(x0 + 30, ey + 110, 'Two loads side by side pay less than two loads apart: they stay together and move as one object: a particle.', 16, INK)
    y = ey + 154
    text(x0 + 30, y, 'The same law as a program on one seam: the echo', 18, INK, weight='bold')
    bw, bh, by_ = 322, 132, y + 18
    steps = [('1  write', ['the seam takes the two loads to', 'their mean and notes the', 'difference x and the coin c', 'in its private memory'], GREEN),
             ('2  phase', ['the memory turns by', '(1 − R)/2 + κ x²,', 'where R flips the sign of x'], RED),
             ('3  unwrite', ['the same step run backwards:', 'the loads come back, exchanged', 'in part; the memory is blank'], GREEN)]
    for k, (ttl, body, col) in enumerate(steps):
        bx_ = x0 + 30 + (bw + 42) * k
        rect(bx_, by_, bw, bh, '#ffffff', col, 1.8, 10)
        text(bx_ + 14, by_ + 28, ttl, 17, col, weight='bold')
        lines(bx_ + 14, by_ + 54, body, 15, 21, INK)
        if k < 2:
            arrow(bx_ + bw + 6, by_ + bh / 2, bx_ + bw + 36, by_ + bh / 2, INK, 2)
    y = by_ + bh + 40
    mono(x0 + 30, y, 'h_e = U (G₀ + κ G₁) U       U = [[0, Wᵀ], [W, 1 − W Wᵀ]]       G₀ = (1 − R)/2,   G₁ = x²', 16)
    mono(x0 + 30, y + 28, 'W |a, b⟩ = ( |c=0, x=a−b⟩ |y₀⟩ + |c=1, x=a−b⟩ |y₁⟩ ) / √2        y_c: the two loads at their mean', 16)
    mono(x0 + 30, y + 56, 'with the memory blank:   h_e = (1 − SWAP_e)/2 + κ (n_a − n_b)²', 16, INK, 'bold')
    y += 100
    text(x0 + 30, y, 'What the law keeps, exactly', 18, INK, weight='bold')
    y = lines(x0 + 30, y + 28, ['the number of loads, and how many ports hold one load and how many hold two',
                                'for every seam: whether its memory is blank, and the size |x| of a record it holds'], 16, 24, INK)
    text(x0 + 30, y + 14, 'What it gives', 18, INK, weight='bold')
    gy2 = y + 62
    # the pair
    for j in range(4):
        if j < 3:
            line(x0 + 60 + 50 * j + 17, gy2, x0 + 60 + 50 * (j + 1) - 17, gy2, BLUE if j == 0 else '#6e7781', 5 if j == 0 else 2.4)
        circle(x0 + 60 + 50 * j, gy2, 17, '#ffffff', BLUE, 2.2)
        if j < 2:
            circle(x0 + 60 + 50 * j, gy2, 7.5, INK, INK)
    lines(x0 + 250, gy2 - 8, ['Particle stability: a pair written on a seam is read as a pair with probability',
                              'at least 1 − 40/(2κ − 6)² for κ > 3, at every size of the carrier and at every time.'], 16, 23, INK)
    gy3 = gy2 + 78
    line(x0 + 77, gy3, x0 + 143, gy3, '#6e7781', 2.4)
    for j in range(2):
        circle(x0 + 60 + 100 * j, gy3, 17, '#ffffff', BLUE, 2.2)
    rect(x0 + 91, gy3 - 13, 38, 26, '#dafbe1', GREEN, 2, 5)
    text(x0 + 110, gy3 + 5, '±1', 14, INK, 'middle', 'bold')
    lines(x0 + 250, gy3 - 8, ['An unrecovered seam: a memory that holds a record while its ports are empty',
                              'stays that way. It holds no load and lies κ above a blank seam.'], 16, 23, INK)
    text(x0 + 30, gy3 + 62, 'As events: every seam acts once per round by exp(−i τ h_e); each port observer counts its own rounds.', 16, INK)
    text(x0 + 30, y0 + h - 24, 'Supplied: the coupling κ (8 and 128 in the runs), the emitter that writes a pair, the reader. Candidate: A-contrast-memory-spherical-v1.', 15, MUTED)

    # ---------------- panel 5: law 2
    x0, y0, w, h = 1220, 2352, 1140, 1080
    panel(x0, y0, w, h, 'SOFTWARE 2', 'Law 2: an object hops between a seam and a wedge', 'law2')
    y = y0 + 92
    mono(x0 + 30, y, 'D = [[ m, B ], [ Bᵀ, −m ]]       D² = m² + K       K = B Bᵀ', 20, INK, 'bold')
    mono(x0 + 30, y + 32, 'B: seam (a,b) to wedge (a,b,c) with +1/√2, seam (b,c) to the same wedge with −1/√2', 16)
    mono(x0 + 30, y + 58, 'K = half the Laplacian of the graph of seams.   Frequencies: ω² = m² + λ, λ a level of K', 16)
    lx, ly = x0 + 50, y + 150
    line(lx, ly, lx + 230, ly, INK, 3)
    text(lx + 242, ly + 6, '+m   the object on a seam', 16, INK)
    line(lx, ly + 150, lx + 230, ly + 150, PURPLE, 3)
    text(lx + 242, ly + 156, '−m   the object in transit', 16, PURPLE)
    arrow(lx + 115, ly + 12, lx + 115, ly + 138, PURPLE, 2.2, True)
    mono(lx + 128, ly + 80, 'hop ±1/√2', 15, PURPLE)
    lines(lx, ly + 200, ['Two shelves and a hop between them:', 'a first-order law whose square is', 'the wave relation ω² = m² + λ.'], 15, 21, MUTED)
    sx, sy0 = x0 + 560, y + 100
    text(sx, sy0, 'the levels of one cell (m = 1)', 16, INK, weight='bold')
    scale = 56
    mid = sy0 + 190
    up = [(1.0, 1), (1.543362, 3), (2.0, 5), (2.148961, 3), (2.449490, 18)]
    down = [(1.0, 91)] + up[1:]
    rect(sx - 10, mid + 0.75 * scale, 276, 2.2 * scale, '#fbefff', 'none', 0, 6, 0.9)
    line(sx - 10, mid, sx + 520, mid, FAINT, 1.2, dash='4 5')
    text(sx + 526, mid + 5, '0', 14, MUTED)
    for sign, levels, col in ((-1, up, INK), (+1, down, PURPLE)):
        slot = None
        for val, mult in levels:
            ly_ = mid + sign * val * scale
            slot = ly_ if slot is None else (min(ly_, slot - 19) if sign < 0 else max(ly_, slot + 19))
            line(sx, ly_, sx + 100, ly_, col, 2.6)
            line(sx + 100, ly_, sx + 124, slot, col, 1.0)
            mono(sx + 130, slot + 5, '%s%.4f  ×%d' % ('+' if sign < 0 else '−', val, mult), 14, col)
    text(sx + 290, mid + 1.45 * scale, 'filled: the reference', 15, PURPLE, weight='bold')
    text(sx + 290, mid + 1.45 * scale + 21, '120 objects, the calm state', 15, PURPLE)
    text(sx + 290, mid - 1.75 * scale, 'one object more:', 15, INK, weight='bold')
    text(sx + 290, mid - 1.75 * scale + 21, 'a particle, at least m above', 15, INK)
    y = mid + 3.0 * scale + 30
    text(x0 + 30, y, 'The 150 modes of one cell', 18, INK, weight='bold')
    bx0, bwid, by2 = x0 + 30, 1080, y + 16
    parts = [(280, '23  waves over the twelve ports', '#ddf4ff', BLUE), (270, '36  flat branch: the cages', '#f6f8fa', MUTED), (530, '91  circulations at −m: the local clocks', '#fbefff', PURPLE)]
    xx = bx0
    for wseg, lab, fill, col in parts:
        rect(xx, by2, wseg - 4, 40, fill, col, 1.6, 6)
        text(xx + 12, by2 + 26, lab, 15, INK, weight='bold')
        xx += wseg
    y = by2 + 76
    text(x0 + 30, y, 'Many objects with exclusion', 18, INK, weight='bold')
    mono(x0 + 30, y + 30, 'H₀ = Σ_seams n_e − Σ_wedges n_w + (1/√2) Σ_wedges [ (c†_ab − c†_bc) c_w + h.c. ]', 16)
    y += 70
    text(x0 + 30, y, 'On ordinary registers (one cell)', 18, INK, weight='bold')
    mono(x0 + 30, y + 30, 'H = H₀ + 1024 Σ_ports (1 − P_b) + 512 Σ_loops (P_S − C_l) + 1024 Σ_unused (1 − |0⟩⟨0|)', 16)
    y = lines(x0 + 30, y + 58, ['432 two-state registers: 120 transit bits, 120 bits of the port codes, 192 apparatus bits;',
                                '259 loop constraints keep the registers a faithful picture of objects that exclude each other;',
                                'a broken code or loop costs 1024, far above every level of the law.'], 16, 24, INK)
    text(x0 + 30, y + 14, 'What it gives', 18, INK, weight='bold')
    lines(x0 + 30, y + 42, ['The lines of an added object are read at 1.0000, 1.5435, 2.0000, 2.1490 of the observer’s own',
                            'reference, predicted from spatial comparisons before the records exist (one cell).',
                            'A created object costs at least m, holds no charge of the twelve currents, and its number is kept.',
                            'On the sphere: round seam terms D = [[m, G B T], [T Bᵀ G, −m]] give the levels of a round sphere.',
                            'A cage in the flat branch swings after a pulse and is kept; see the next panel.'], 16, 24, INK)
    text(x0 + 30, y0 + h - 24, 'Supplied: m = 1, pulses, ramps, apparatus bits. Candidates: A-scalar-CAR-pair-transit-v2, A-native-scalar-overlap-v3, the round sibling.', 15, MUTED)

    # ---------------- panel 6: a cage in the flat branch of law 2
    x0, y0, w, h = 40, 3472, 2320, 620
    panel(x0, y0, w, h, 'SOFTWARE 3', 'A cage in the flat branch of law 2: a ring of four seams, its swing and its persistence', 'cage')
    cage = json.loads((HERE.parent / 'scalar_channel' / 'cage' / 'results.json').read_text(encoding='utf-8'))
    ring = (0, 2, 4, 6)
    sign = {}
    for k in range(4):
        a, b = ring[k], ring[(k + 1) % 4]
        sign[(min(a, b), max(a, b))] = (-1) ** k
    Mc = align(X[0] + X[4], (0, 0, 1))
    d = (X @ Mc.T)[0] - (X @ Mc.T)[4]
    Mc = rot('z', math.degrees(math.atan2(d[0], d[1]))) @ Mc
    ccx, ccy, cR = x0 + 300, y0 + 330, 200
    xyc = draw_cell(X, seams, ccx, ccy, cR, Mc, signed=sign, pr=15, fs=14)
    for (a, b), sg in sign.items():
        pa, pb = np.array(xyc(a)), np.array(xyc(b))
        mid = (pa + pb) / 2
        out_dir = mid - np.array([ccx, ccy])
        out_dir = out_dir / (np.linalg.norm(out_dir) or 1)
        circle(mid[0] + 24 * out_dir[0], mid[1] + 24 * out_dir[1], 13, '#ffffff', BLUE if sg > 0 else RED, 1.6)
        text(mid[0] + 24 * out_dir[0], mid[1] + 24 * out_dir[1] + 7, '+' if sg > 0 else '−', 20, BLUE if sg > 0 else RED, 'middle', 'bold')
    text(ccx, y0 + 86, 'the ring 0, 2, 4, 6 with the signs + − + −', 16, INK, 'middle', 'bold')
    text(ccx, y0 + h - 44, 'At each port of the ring a plus meets a minus:', 15, MUTED, 'middle')
    text(ccx, y0 + h - 24, 'the port passes nothing on to its other seams.', 15, MUTED, 'middle')

    tx = x0 + 600
    text(tx, y0 + 92, 'Why the law keeps it', 18, INK, weight='bold')
    mono(tx, y0 + 124, 'e  = ½ ( |0-2⟩ − |2-4⟩ + |4-6⟩ − |6-0⟩ )', 16)
    mono(tx, y0 + 150, 'e′ = Bᵀ e / √5     on the 28 wedges at the four ports', 16)
    mono(tx, y0 + 184, '2K = diag(d_a + d_b) − QᵀQ,    Q e = 0   ⇒   K e = 5 e', 16)
    mono(tx, y0 + 210, 'D (e, e′) = (e, e′) [[ m, √5 ], [ √5, −m ]]', 16, INK, 'bold')
    lines(tx, y0 + 244, ['Q adds the seams at a port without signs. The law moves e to e′',
                         'and back and to no other place: 4 seams and 28 wedges.',
                         '18 such profiles in a cell; its 30 rings span them.',
                         'On a federation: rings on ports with five seams.'], 15, 21, MUTED)
    lx, ly = tx + 20, y0 + 372
    line(lx, ly, lx + 190, ly, INK, 3)
    text(lx + 202, ly + 5, '+E = +2.4495   upper mode, empty in the reference', 15, INK)
    line(lx, ly + 150, lx + 190, ly + 150, PURPLE, 3)
    circle(lx + 95, ly + 150, 8, PURPLE, PURPLE)
    text(lx + 202, ly + 155, '−E = −2.4495   lower mode, filled in the reference', 15, PURPLE)
    arrow(lx + 95, ly + 136, lx + 95, ly + 12, AMBER, 2.4)
    lines(lx + 112, ly + 66, ['lifting the object: a neutral excitation,', '2E = 4.8990'], 15, 21, AMBER, weight='bold')
    text(tx, y0 + h - 44, 'Kept exactly: [H, P_X] = 0 for the excitation, at every reading.', 15, INK)
    text(tx, y0 + h - 24, 'One added object in the upper mode costs E; the swing is at 2E.', 15, MUTED)

    # the swing P_e(t) = 1 - (5/6) sin^2(sqrt6 t)
    gx, gy, gw, gh = x0 + 1240, y0 + 130, 500, 250
    text(gx - 40, y0 + 92, 'The swing after a finite pulse', 18, INK, weight='bold')
    rect(gx, gy, gw, gh, '#fbfcfd', FAINT, 1.2, 4)
    Eg = math.sqrt(6)
    tmax = 4.0
    pts = []
    for k in range(401):
        t = tmax * k / 400
        pv = 1 - (5 / 6) * math.sin(Eg * t) ** 2
        pts.append('%.1f,%.1f' % (gx + gw * t / tmax, gy + gh * (1 - pv)))
    out.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="2.6"/>' % (' '.join(pts), BLUE))
    for val, lab in ((1.0, '1'), (1 / 6, '1/6')):
        yy = gy + gh * (1 - val)
        line(gx, yy, gx + gw, yy, FAINT, 1.0, dash='3 5')
        text(gx - 8, yy + 5, lab, 14, MUTED, 'end')
    per = math.pi / Eg
    for k in range(1, 4):
        xx = gx + gw * k * per / tmax
        line(xx, gy + gh, xx, gy + gh + 7, MUTED, 1.4)
    text(gx + gw * per / tmax, gy + gh + 24, 'π/E', 14, MUTED, 'middle')
    text(gx + gw * 2 * per / tmax, gy + gh + 24, '2π/E', 14, MUTED, 'middle')
    text(gx + gw * 3 * per / tmax, gy + gh + 24, '3π/E', 14, MUTED, 'middle')
    text(gx, gy + gh + 24, '0', 14, MUTED, 'middle')
    text(gx + gw / 2, gy + gh + 46, 'the observer’s own phase', 14, MUTED, 'middle')
    mono(gx - 40, gy + gh + 80, 'P_e = 1 − (σ²/E²) sin²(E t):  a constant and one line at 2E', 15)
    lines(gx - 40, gy + gh + 108, ['The pulse H − 2E n(e), for π / (2 √(2E(E − m))), puts the object of the',
                                  'lower mode on e and does the work 3.456. Then no drive: the record',
                                  'starts from rest and swings. The reference alone does not swing.'], 15, 21, MUTED)

    rx = x0 + 1810
    text(rx, y0 + 92, 'What observers read', 18, INK, weight='bold')
    text(rx, y0 + 124, 'The clock, from counts of one cell', 16, BLUE, weight='bold')
    lines(rx, y0 + 148, ['the line of the seam record, against the', 'observer’s own reference:'], 15, 20, INK)
    mono(rx + 16, y0 + 196, '%.6f ± %.6f' % (cage['clock']['steps_of_equal_length']['line']['frequency'], cage['clock']['steps_of_equal_length']['line']['standard_error']), 16, INK, 'bold')
    lines(rx, y0 + 222, ['the level difference of the pair, read', 'separately with a fuel dial:'], 15, 20, INK)
    mono(rx + 16, y0 + 270, '%.6f ± %.6f' % (cage['dial']['level_difference'], cage['dial']['level_difference_standard_error']), 16, INK, 'bold')
    mono(rx + 16, y0 + 294, 'the law: 2√6 = 4.898979', 15, MUTED)
    text(rx, y0 + 332, 'Persistence: the excitation at five waits', 16, GREEN, weight='bold')
    mono(rx + 16, y0 + 358, 'sharp pulse      %.5f' % cage['object']['excitation_in_the_state']['sharp'], 15)
    mono(rx + 16, y0 + 380, 'coherent pulse   %.5f' % cage['object']['excitation_in_the_state']['coherent'], 15)
    mono(rx + 16, y0 + 402, 'no pulse         %.6f' % cage['object']['excitation_in_the_state']['none'], 15)
    text(rx, y0 + 440, 'Under a changed law', 16, RED, weight='bold')
    lines(rx, y0 + 464, ['(seam hopping round every port, same', 'alphabet, routing and symmetry): no single',
                         'line, and %.2f of the cage left at 64.' % cage['object']['changed_law']['weight_of_the_lower_mode_left_in_the_pair_at_1_4_16_64'][3]], 15, 20, INK)
    text(rx, y0 + h - 44, 'Supplied: the potential and the converter on e,', 15, MUTED)
    text(rx, y0 + h - 24, 'one guard over four port codes, the fuel dial.', 15, MUTED)

    # ---------------- panel 7: observers and a run
    x0, y0, w, h = 40, 3472 + 620 + 40, 2320, 500
    panel(x0, y0, w, h, 'SOFTWARE 4', 'Observers, instruments, clocks: one run from preparation to record', 'run')
    steps = [('1  prepare', ['measure every loop, keep the', 'sign, correct it; fill the', 'transits; ramp each star to', 'the calm reference'], MUTED),
             ('2  write', ['a fuel at level 0 is coupled', 'to a local target for a short', 'pulse while the law keeps', 'running'], AMBER),
             ('3  let it run', ['the law acts alone; the', 'observer counts cycles of its', 'own circulation, three', 'transits round its port'], PURPLE),
             ('4  read', ['a second pulse with one of', 'four phase settings; the fuel', 'is copied into two pointers', 'after the coupling is off'], AMBER),
             ('5  agree', ['two observers on an overlap', 'compare complete reports;', 'identical reports are admitted,', 'every other outcome is kept'], GREEN),
             ('6  record', ['three copies of the admitted', 'record; the order and count', 'of its own records are the', 'observer’s clock'], BLUE)]
    bw, bh, by_ = 340, 170, y0 + 76
    for k, (ttl, body, col) in enumerate(steps):
        bx_ = x0 + 28 + (bw + 44.8) * k
        rect(bx_, by_, bw, bh, '#ffffff', col, 2, 12)
        text(bx_ + 16, by_ + 32, ttl, 19, col, weight='bold')
        lines(bx_ + 16, by_ + 62, body, 16, 24, INK)
        if k < 5:
            arrow(bx_ + bw + 6, by_ + bh / 2, bx_ + bw + 39, by_ + bh / 2, INK, 2.2)
    y = by_ + bh + 44
    mono(x0 + 30, y, 'clock:   χ_b = ( |a,b,c⟩ + |c,b,f⟩ − |a,b,f⟩ ) / √3,   D χ_b = −m χ_b          reading:   p(fuel full) − p(fuel empty) = |amplitude of one added object|²', 16)
    mono(x0 + 30, y + 28, 'four settings θ = 0, π/2, π, 3π/2 turn the differences into one complex number per count q; its spectrum over q shows the frequencies in units of the observer’s reference', 15, MUTED)
    y += 76
    text(x0 + 30, y, 'Observers first.', 18, INK, weight='bold')
    lines(x0 + 180, y, ['Ports, seams and cells are the hardware of a computation. They are no space. What is compared with nature is what the observers record and agree on.',
                        'No observer reads the state of the program, its evolution parameter or a common clock; two observers compare clocks through the records on their overlap.'], 17, 26, INK)
    text(x0 + 30, y0 + h - 22, 'Both laws are candidates under test. Law 1 is specified in model/contrast_memory_spherical_v1/SPEC.md; the package for law 2 is model/scalar_channel, with the cage in model/scalar_channel/cage.', 15, MUTED)

    Ht = 3472 + 620 + 40 + 500 + 40
    defs = '<defs>' + ''.join(
        '<marker id="ah-%s" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="%s"/></marker>'
        '<marker id="as-%s" viewBox="0 0 10 10" refX="1" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M10,0 L0,5 L10,10 z" fill="%s"/></marker>' % (c[1:], c, c[1:], c)
        for c in (INK, BLUE, ORANGE, PURPLE, GREEN, AMBER, MUTED)) + '</defs>'
    body = out[1:]
    head_ = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" role="img" aria-label="%s">'
    sheet = [head_ % (Wd, Ht, Wd, Ht, 'The executable construction of A-contrast-memory-spherical-v1: hardware and software'), defs,
             '<rect x="0" y="0" width="%d" height="%d" fill="#ffffff"/>' % (Wd, Ht)] + body + ['</svg>']
    (HERE / 'executable_spec.svg').write_text('\n'.join(sheet) + '\n', encoding='utf-8')
    print('wrote executable_spec.svg', Wd, Ht)
    # every panel as a figure of its own, for the pages of the manual
    bounds = PANELS + [[None, 0, 0, 0, 0, len(out)]]
    for k, (name, x, y, w, h, start) in enumerate(PANELS):
        end = bounds[k + 1][5]
        m = 8
        fig = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="%d %d %d %d" width="%d" height="%d" role="img" aria-label="%s">' % (
            x - m, y - m, w + 2 * m, h + 2 * m, w + 2 * m, h + 2 * m, name), defs,
            '<rect x="%d" y="%d" width="%d" height="%d" fill="#ffffff"/>' % (x - m, y - m, w + 2 * m, h + 2 * m)] + out[start:end] + ['</svg>']
        (HERE / ('fig_%s.svg' % name)).write_text('\n'.join(fig) + '\n', encoding='utf-8')
        print('wrote fig_%s.svg' % name)


if __name__ == '__main__':
    main()
