"""Compact 16:9 social card of the explainer: six steps from a twist to a particle."""
import argparse
import math
from pathlib import Path

from . import build
from .canvas import (AMBER, BG, CARD, CARD2, EDGE, FAINT, GOLD, GREEN, MONO, MUTED, SANS, SERIF,
                     SOFT, TEAL, TEXT, VIOLET, Canvas, n)
from .parts import bell, hexkeep, lattice, lcg, loop_triangle, mesh, tri_fill, vortex
from .scenes_a import _centroid

W, H = 1600, 900
M = 56


def panel_twist(c, x, y, w, h):
    loop_triangle(c, x + w / 2, y + h / 2 + 12, 132, (1, 1, -1), 'HHH', small=True, show_rules=False,
                  product='−1')


def panel_protect(c, x, y, w, h):
    cx, cy, a = x + w / 2, y + h / 2 - 10, 42
    pts, edges, tris = lattice(cx, cy, a, hexkeep(cx, cy, 2 * a + 1))
    mesh(c, pts, edges, TEAL, 1.1, 0.3, 2, 0.6)
    inner = [t for t in tris if all(math.hypot(pts[k][0] - cx, pts[k][1] - cy) < a + 1 for k in t)]
    inner.sort(key=lambda t: math.atan2(_centroid(pts, t)[1] - cy, _centroid(pts, t)[0] - cx))
    g0 = _centroid(pts, inner[3])
    tri_fill(c, pts, inner[3], AMBER, 0.22, stroke=False)
    gx, gy = tri_fill(c, pts, inner[0], AMBER, 0.85, glow=40)
    c.qarrow(g0[0], g0[1] + 3, cx, cy + 24, gx + 2, gy + 5, TEXT, 1.5, 6, 0.9)
    ring = [(cx + a * math.cos(math.radians(60 * i)), cy + a * math.sin(math.radians(60 * i))) for i in range(6)]
    c.poly(ring, 'none', TEXT, 2.4)
    c.text(cx, y + h - 6, 'edge reads −1', 17, AMBER, 700, 'middle')


def panel_move(c, x, y, w, h):
    cell, rows, cols = 22, 6, 8
    gx, gy = x + (w - cols * cell) / 2 + 12, y + 30
    for r in range(rows + 1):
        c.line(gx, gy + r * cell, gx + cols * cell, gy + r * cell, EDGE, 1)
    for q in range(cols + 1):
        c.line(gx + q * cell, gy, gx + q * cell, gy + rows * cell, EDGE, 1)
    apex = (gx + 1.5 * cell, gy + 0.5 * cell)
    c.line(apex[0], apex[1], apex[0] + (rows - 0.5) * cell, gy + rows * cell, TEAL, 1.4, 0.8, '5 5')
    path = [1 + t // 2 for t in range(rows)]
    for t, q in enumerate(path):
        xx, yy = gx + q * cell, gy + t * cell
        if t == rows - 1:
            c.glow(xx + cell / 2, yy + cell / 2, 24, AMBER, 0.8)
        c.rect(xx + 2.5, yy + 2.5, cell - 5, cell - 5, AMBER, None, 0, 4, round(0.3 + 0.7 * t / (rows - 1), 2))
    c.text(gx - 8, gy + 12, 'tick', 14, MUTED, anchor='end')
    c.text(gx + cols * cell, gy - 6, 'place →', 14, MUTED, anchor='end')


def panel_position(c, x, y, w, h):
    nd, dw = 9, 22
    dx0, base = x + (w - nd * dw) / 2, y + h - 66
    hump = [(dx0 + nd * dw * t / 60, base - 10 - 62 * bell(t / 60, 0.5, 0.17)) for t in range(61)]
    c.path('M' + ' L'.join(f'{n(a)},{n(b)}' for a, b in hump) + f' L{n(hump[-1][0])},{n(base - 10)} '
           f'L{n(hump[0][0])},{n(base - 10)} Z', fill='url(#hump)')
    c.curve(hump, AMBER, 2)
    for i in range(nd):
        lit = i == 5
        if lit:
            c.glow(dx0 + i * dw + dw / 2, base + 8, 26, AMBER, 0.85)
        c.rect(dx0 + i * dw + 2, base - 2, dw - 4, 20, AMBER if lit else '#17233b', '#fff' if lit else EDGE,
               1.4 if lit else 1, 4)
    c.text(x + w / 2, base + 50, 'one click: here, now', 16, SOFT, 600, 'middle')


def panel_velocity(c, x, y, w, h):
    ox, oy, aw, ah = x + 40, y + 14, w - 72, h - 52
    c.line(ox, oy, ox, oy + ah, FAINT, 1.2)
    c.line(ox, oy + ah, ox + aw, oy + ah, FAINT, 1.2)
    recs = [(ox + 22 + k * (aw - 40) / 3, oy + ah - 18 - k * (ah - 36) / 3) for k in range(4)]
    c.line(recs[0][0] - 14, recs[0][1] + 13, recs[-1][0] + 14, recs[-1][1] - 13, VIOLET, 1.8, 0.9, '5 5')
    for a, b in recs:
        c.glow(a, b, 14, AMBER, 0.8)
        c.circle(a, b, 5, AMBER, '#fff', 1.4)
    c.text(ox + aw, oy + ah + 22, 'place', 14, MUTED, anchor='end')
    c.text(ox - 6, oy + 10, 'tick', 14, MUTED, anchor='end')


def panel_mass(c, x, y, w, h):
    cx, cy = x + w / 2, y + h / 2 - 2
    for r, op in ((70, 0.16), (50, 0.3)):
        c.circle(cx, cy, r, 'none', AMBER, 1.4, op, '4 5')
    c.glow(cx, cy, 44, AMBER, 0.9)
    c.circle(cx, cy, 9, AMBER, '#fff', 1.8)
    for k in range(8):
        ang = math.radians(45 * k + 22)
        c.arrow(cx + 66 * math.cos(ang), cy + 66 * math.sin(ang), cx + 24 * math.cos(ang),
                cy + 24 * math.sin(ang), TEAL, 1.4, 6, 0.8)
    c.text(cx, y + h - 6, 'E = mc² at rest', 17, TEXT, 700, 'middle', family=MONO)


STEPS = (
    ('A twist', 'A mismatch owned by a whole loop. Any seam can hold it.', TEAL, panel_twist),
    ('Protected', 'Repairs move it around. None can remove it.', TEAL, panel_protect),
    ('It moves', 'Re-established one place over each tick. Nothing travels.', AMBER, panel_move),
    ('Position', 'A detector writes one record of where.', AMBER, panel_position),
    ('Velocity', 'Records line up. Places per tick is the speed.', AMBER, panel_velocity),
    ('Mass', 'Holding the core costs repair at every tick. That cost is the mass.', AMBER, panel_mass),
)


def render():
    c = Canvas()
    c.define('hump', f'<linearGradient id="hump" x1="0" x2="0" y1="0" y2="1"><stop offset="0" '
                     f'stop-color="{AMBER}" stop-opacity=".5"/><stop offset="1" stop-color="{AMBER}" '
                     f'stop-opacity="0"/></linearGradient>')
    c.rect(0, 0, W, H, BG)
    c.define('dots', '<pattern id="dots" width="32" height="32" patternUnits="userSpaceOnUse">'
                     '<circle cx="2" cy="2" r="1.1" fill="#19233a"/></pattern>')
    c.rect(0, 0, W, H, 'url(#dots)')
    c.glow(200, 120, 520, TEAL, 0.16)
    c.glow(1420, 160, 560, VIOLET, 0.2)

    c.caps(M, 62, 'OPH particle atlas', TEAL, 16, ls=3.2)
    c.text(W - M, 62, 'github.com/FloatingPragma/observer-patch-holography', 16, MUTED, anchor='end')
    c.add(f'<text x="{M}" y="136" font-family="{SANS}" font-size="66" font-weight="800" '
          f'letter-spacing="-1.6" fill="{TEXT}">From a <tspan fill="{AMBER}">twist</tspan> to a '
          f'<tspan fill="{VIOLET}">particle</tspan></text>')
    c.text(M, 184, 'In the OPH simulation, matter is a topological defect: a twist that no repair can remove.',
           26, SOFT, maxw=W - 2 * M)
    c.text(M, 220, 'Here is how observers come to see it as a particle with a position, a velocity and a mass.',
           26, SOFT, maxw=W - 2 * M)

    gap, top, ch = 16, 292, 424
    cw = (W - 2 * M - 5 * gap) / 6
    zones = ((0, 1, 'In the simulation', TEAL), (2, 5, 'What observers measure', AMBER))
    for a, b, lab, col in zones:
        xa, xb = M + a * (cw + gap), M + b * (cw + gap) + cw
        c.line(xa + 2, top - 14, xb - 2, top - 14, col, 2.4, 0.8)
        c.caps(xa + 2, top - 26, lab, col, 14)
    for i, (title, body, col, draw) in enumerate(STEPS):
        x = M + i * (cw + gap)
        c.rect(x, top, cw, ch, CARD, EDGE, 1.4, 20)
        c.circle(x + 32, top + 36, 17, col)
        c.text(x + 32, top + 43, str(i + 1), 19, BG, 800, 'middle')
        c.text(x + 60, top + 45, title, 26, TEXT, 700, maxw=cw - 68)
        c.rect(x + 12, top + 70, cw - 24, 200, CARD2, None, 0, 14)
        draw(c, x + 12, top + 70, cw - 24, 200)
        c.para(x + 20, top + 306, cw - 36, body, 20.5, SOFT, 1.4)
        if i < 5:
            mx = x + cw + gap / 2
            c.circle(mx, top + 170, 13, BG, EDGE, 1.2)
            c.arrow(mx - 5, top + 170, mx + 6, top + 170, SOFT, 1.8, 6)

    by = top + ch + 20
    c.rect(M, by, W - 2 * M, 100, '#151129', VIOLET, 1.3, 18)
    c.circle(M + 40, by + 50, 20, VIOLET)
    c.text(M + 40, by + 57, '7', 20, BG, 800, 'middle')
    c.text(M + 76, by + 44, 'Fields', 26, TEXT, 700)
    c.text(M + 76, by + 76, 'Observers compress the agreed records into smooth fields. Their smallest '
                            'excitations are what physics calls particles.', 20, SOFT, maxw=W - 2 * M - 110)

    fy = H - 26
    c.circle(M + 8, fy - 6, 6.5, TEAL)
    c.text(M + 24, fy, 'Computed in the finite model: the twist, its protection, position and momentum records',
           17, MUTED)
    c.circle(M + 780, fy - 6, 6.5, 'none', AMBER, 1.8, dash='3.2 2.6')
    c.text(M + 796, fy, 'Work in progress: free motion of the defect and its rest mass', 17, MUTED)
    if c.overflow:
        raise ValueError('text overflow: ' + '; '.join(f'{d} {s[:60]!r}' for d, s in c.overflow))
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
            f'role="img" aria-label="From a twist to a particle: OPH particle atlas card">'
            f'<defs>{"".join(c.defs)}</defs>' + ''.join(c.items) + '</svg>\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('out', type=Path, help='directory for the SVG card')
    a = p.parse_args()
    build.load()
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / 'oph_particle_atlas_v3_card.svg').write_text(render())
    print(a.out / 'oph_particle_atlas_v3_card.svg')


if __name__ == '__main__':
    main()
