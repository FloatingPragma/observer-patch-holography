"""Drawing primitives and carrier geometry for the figures of the executable model in the papers.

The wiring of one cell is read from model/executable_model/carriers/CELL.npz and the federations from
FEDERATION_LEVEL{0,1,2}.npz in the same folder. make_paper_figures.py draws the figures.
"""
import collections
import itertools
import json
import math
from pathlib import Path

import networkx as nx
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CARRIERS = ROOT / 'model' / 'executable_model' / 'carriers'
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
    z = np.load(CARRIERS / 'CELL.npz', allow_pickle=False)
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
    z = np.load(CARRIERS / f'FEDERATION_LEVEL{level}.npz', allow_pickle=False)
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
