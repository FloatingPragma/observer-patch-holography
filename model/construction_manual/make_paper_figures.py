"""Draws the figures of the executable model for the papers: hardware, particle stability, the cage, one run.

The geometry and the drawing primitives are those of make_diagram.py, so the figures show the same cells and
federations as the sheet of the construction manual. The probabilities of the stability figure are computed from
the law of ../contrast_memory_spherical_v1 on the federation of 20 cells; the numbers of the cage figure are read from
../scalar_channel/cage/results.json.

    python model/construction_manual/make_paper_figures.py      # writes paper/figures/model_*.svg and .pdf

The PDFs are rendered with rsvg-convert and normalized with Ghostscript, as the book's diagrams are; both tools
must be on the PATH. The papers include the committed PDFs.
"""
from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT))

import make_diagram as md  # noqa: E402
from make_diagram import (BLUE, FAINT, GREEN, INK, MUTED, ORANGE, PURPLE, RED, AMBER, align, arrow,  # noqa: E402
                          circle, line, path, poly, rect, rot, text)

OUT = ROOT / 'paper' / 'figures'
LIGHT_BLUE, LIGHT_PURPLE, LIGHT_GREEN, LIGHT_RED, LIGHT_AMBER = '#ddf4ff', '#fbefff', '#dafbe1', '#ffebe9', '#fff8c5'


def label(x, y, s, size=24, fill=INK, anchor='start', weight='normal'):
    text(x, y, s, size, fill, anchor, weight)


def tag(x, y, s):
    text(x, y, s, 26, INK, 'start', 'bold')


def write(name, width, height):
    defs = '<defs>' + ''.join(
        '<marker id="ah-%s" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto">'
        '<path d="M0,0 L10,5 L0,10 z" fill="%s"/></marker>'
        '<marker id="as-%s" viewBox="0 0 10 10" refX="1" refY="5" markerWidth="7" markerHeight="7" orient="auto">'
        '<path d="M10,0 L0,5 L10,10 z" fill="%s"/></marker>' % (c[1:], c, c[1:], c)
        for c in (INK, BLUE, ORANGE, PURPLE, GREEN, AMBER, MUTED, RED)) + '</defs>'
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d">' % (width, height, width, height),
           defs, '<rect x="0" y="0" width="%d" height="%d" fill="#ffffff"/>' % (width, height)] + md.out + ['</svg>']
    OUT.mkdir(parents=True, exist_ok=True)
    source = OUT / f'{name}.svg'
    source.write_text('\n'.join(svg) + '\n', encoding='utf-8')
    with tempfile.TemporaryDirectory() as tmp:
        raw = Path(tmp) / 'raw.pdf'
        subprocess.run(['rsvg-convert', '-f', 'pdf', '-o', str(raw), str(source)], check=True)
        subprocess.run(['gs', '-q', '-dSAFER', '-dBATCH', '-dNOPAUSE', '-sDEVICE=pdfwrite', '-dCompatibilityLevel=1.5',
                        f'-sOutputFile={OUT / (name + ".pdf")}', str(raw)], check=True)
    md.out.clear()
    print('wrote', source.relative_to(ROOT), 'and its PDF')


# ---------------------------------------------------------------- hardware

def sphere(level, cx, cy, R, M):
    cells, n_seams, n_ext, pairs, cent, sym, _, _ = md.federation(level)
    P = cent @ M.T
    circle(cx, cy, R * 1.035, '#fbfcfd', FAINT, 1.4)
    for a, b in sorted(pairs, key=lambda e: P[e[0], 2] + P[e[1], 2]):
        back = (P[a, 2] + P[b, 2]) / 2 < 0
        line(cx + R * P[a, 0], cy - R * P[a, 1], cx + R * P[b, 0], cy - R * P[b, 1], ORANGE,
             0.9 if back else (2.4 if level < 2 else 1.4), opacity=0.16 if back else 0.95)
    r = {0: 11, 1: 6.5, 2: 3.4}[level]
    for c in sorted(range(cells), key=lambda q: P[q, 2]):
        if P[c, 2] >= 0:
            circle(cx + R * P[c, 0], cy - R * P[c, 1], r, BLUE, BLUE, 1.2)
    return cells, n_seams, n_ext


def hardware():
    seams, X = md.cell()
    W, H = 1400, 610
    tag(20, 40, '(a) one cell')
    M = rot('x', -62) @ rot('z', 18) @ align(X[0], (0, 0, 1))
    md.draw_cell(X, seams, 250, 265, 190, M, face=(2, 4, 8), wedge=(1, 0, 2), pr=19, fs=18)
    lines = ['12 ports, 30 seams, 20 faces', 'a port is an observer; a seam joins two', 'ports and keeps a private memory;',
             'a wedge (purple) is two seams at a port']
    for k, s in enumerate(lines):
        label(250, 500 + 26 * k, s, 21, MUTED if k else INK, 'middle', 'bold' if k == 0 else 'normal')
    tag(540, 40, '(b) federations of 20, 80 and 320 cells')
    Ms = rot('x', -68) @ rot('z', 12)
    for k, level in enumerate((0, 1, 2)):
        cx = 655 + 285 * k
        cells, n_seams, n_ext = sphere(level, cx, 290, 125, Ms)
        label(cx, 465, '%d cells' % cells, 24, INK, 'middle', 'bold')
        label(cx, 494, '{:,} ports'.format(12 * cells), 21, MUTED, 'middle')
        label(cx, 520, '{:,} seams'.format(n_seams), 21, MUTED, 'middle')
    label(940, 562, 'a dot is one cell, an orange line one seam between two cells;', 20, MUTED, 'middle')
    label(940, 588, 'each cell meets its three neighbours through one seam each', 20, MUTED, 'middle')
    write('model_hardware', W, H)


# ---------------------------------------------------------------- particle stability

def stability_series():
    from model.contrast_memory_spherical_v1 import carrier as cr
    from model.contrast_memory_spherical_v1 import stability as st
    times = np.linspace(0.0, 32.0, 161)
    ps = st.sector(0)
    series = {}
    for kind in ('law', 'exchange alone', 'agreement energy reversed'):
        lines, *_ = st.emitted_run(ps, st.law(ps, 8.0, kind), cr.EMITTER_SEAM[0], times=times)
        series[kind] = ([r['no_pair'] for r in lines], [r['pair_on_another_seam'] for r in lines])
    return times, series


def plot_frame(x0, y0, w, h, ymax, yticks, xlabel, ylabel, log=False):
    rect(x0, y0, w, h, 'none', FAINT, 1.4, 0)
    for v in yticks:
        y = log_y(y0, h, v) if log else y0 + h - h * v / ymax
        line(x0, y, x0 + w, y, FAINT, 1.0, dash='4 6')
        label(x0 - 12, y + 7, ('%g' % v), 20, MUTED, 'end')
    for t in (0, 8, 16, 24, 32):
        x = x0 + w * t / 32
        line(x, y0 + h, x, y0 + h + 8, MUTED, 1.4)
        label(x, y0 + h + 32, str(t), 20, MUTED, 'middle')
    label(x0 + w / 2, y0 + h + 62, xlabel, 21, MUTED, 'middle')
    if ylabel:
        label(x0, y0 - 16, ylabel, 22, INK, 'start', 'bold')


def log_y(y0, h, v, lo=-4.0, hi=0.0):
    return y0 + h - h * (math.log10(max(v, 10 ** lo)) - lo) / (hi - lo)


def curve(x0, y0, w, h, ymax, times, values, colour, width=3.2, dash=None, log=False):
    ys = [log_y(y0, h, v) if log else y0 + h - h * min(v, ymax) / ymax for v in values]
    pts = ' '.join('%.1f,%.1f' % (x0 + w * t / 32, y) for t, y in zip(times, ys))
    md.out.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="%.1f" stroke-linejoin="round"%s/>' % (
        pts, colour, width, (' stroke-dasharray="%s"' % dash) if dash else ''))


def stability():
    times, series = stability_series()
    results = json.loads((ROOT / 'model' / 'contrast_memory_spherical_v1' / 'RESULTS.json').read_text(encoding='utf-8'))
    bound = results['particle_stability']['law']['20 cells']['by_kappa']['8.0']['bound_with_carrier_constants']
    W, H = 1400, 560
    x0, y0, w, h = 110, 70, 520, 380
    tag(20, 36, '(a) probability of the outcome "no pair"')
    plot_frame(x0, y0, w, h, 1.0, (0.0001, 0.001, 0.01, 0.1, 1), 'evolution parameter of the calculation', None, log=True)
    yb = log_y(y0, h, bound)
    line(x0, yb, x0 + w, yb, RED, 2.2, dash='10 7')
    label(x0 + w - 6, yb - 10, 'bound of the theorem, constants of the carrier', 19, RED, 'end')
    curve(x0, y0, w, h, 1.0, times, series['exchange alone'][0], ORANGE, log=True)
    curve(x0, y0, w, h, 1.0, times, series['agreement energy reversed'][0], PURPLE, 2.6, '3 5', log=True)
    curve(x0, y0, w, h, 1.0, times, series['law'][0], BLUE, 3.6, log=True)
    x1 = 840
    tag(x1 - 110, 36, '(b) probability of the pair on another seam')
    plot_frame(x1, y0, w, h, 1.0, (0.2, 0.4, 0.6, 0.8, 1.0), 'evolution parameter of the calculation', None)
    curve(x1, y0, w, h, 1.0, times, series['exchange alone'][1], ORANGE)
    curve(x1, y0, w, h, 1.0, times, series['agreement energy reversed'][1], PURPLE, 2.6, '3 5')
    curve(x1, y0, w, h, 1.0, times, series['law'][1], BLUE, 3.6)
    ly = H - 24
    for k, (colour, name, dash) in enumerate(((BLUE, 'the contrast-memory law, κ = 8', None),
                                              (ORANGE, 'the exchange alone, κ = 0', None),
                                              (PURPLE, 'the agreement energy reversed, κ = −8', '3 5'))):
        lx = 70 + 440 * k
        line(lx, ly - 7, lx + 46, ly - 7, colour, 3.4, dash=dash)
        label(lx + 58, ly, name, 20, INK)
    write('model_stability', W, H)
    return series


# ---------------------------------------------------------------- the cage

def cage():
    seams, X = md.cell()
    results = json.loads((ROOT / 'model' / 'scalar_channel' / 'cage' / 'results.json').read_text(encoding='utf-8'))
    clock, dial, persistence = results['clock'], results['dial'], results['persistence']
    W, H = 1400, 575
    tag(20, 36, '(a) a ring of four seams')
    ring = (0, 2, 4, 6)
    sign = {}
    for k in range(4):
        a, b = ring[k], ring[(k + 1) % 4]
        sign[(min(a, b), max(a, b))] = (-1) ** k
    Mc = align(X[0] + X[4], (0, 0, 1))
    d = (X @ Mc.T)[0] - (X @ Mc.T)[4]
    Mc = rot('z', math.degrees(math.atan2(d[0], d[1]))) @ Mc
    ccx, ccy, cR = 210, 285, 165
    xyc = md.draw_cell(X, seams, ccx, ccy, cR, Mc, signed=sign, pr=17, fs=16)
    for (a, b), sg in sign.items():
        pa, pb = np.array(xyc(a)), np.array(xyc(b))
        mid = (pa + pb) / 2
        u = mid - np.array([ccx, ccy])
        u = u / (np.linalg.norm(u) or 1)
        circle(mid[0] + 26 * u[0], mid[1] + 26 * u[1], 15, '#ffffff', BLUE if sg > 0 else RED, 1.8)
        label(mid[0] + 26 * u[0], mid[1] + 26 * u[1] + 8, '+' if sg > 0 else '−', 24, BLUE if sg > 0 else RED, 'middle', 'bold')
    label(ccx, 500, 'at every port of the ring a plus meets a minus,', 20, MUTED, 'middle')
    label(ccx, 525, 'so the port passes nothing to its other seams', 20, MUTED, 'middle')

    # (b) levels
    bx = 450
    tag(bx, 36, '(b) the levels of the pair')
    E = math.sqrt(6)
    top, bottom = 150, 400
    line(bx + 20, top, bx + 230, top, INK, 3.2)
    label(bx + 240, top + 8, '+E = +√6', 22, INK)
    label(bx + 240, top + 34, 'upper mode, empty', 19, MUTED)
    line(bx + 20, bottom, bx + 230, bottom, PURPLE, 3.2)
    circle(bx + 125, bottom, 10, PURPLE, PURPLE)
    label(bx + 240, bottom + 8, '−E = −√6', 22, PURPLE)
    label(bx + 240, bottom + 34, 'lower mode, filled', 19, MUTED)
    arrow(bx + 125, bottom - 16, bx + 125, top + 10, AMBER, 3.0)
    label(bx + 140, (top + bottom) / 2 - 6, 'excitation', 21, AMBER, 'start', 'bold')
    label(bx + 140, (top + bottom) / 2 + 22, '2E = 4.898979', 21, AMBER)
    label(bx + 20, 480, 'the law moves the profile of the ring to its', 19, MUTED)
    label(bx + 20, 504, 'wedges and back, and to no other place;', 19, MUTED)
    label(bx + 20, 528, 'the excitation is kept at every reading', 19, MUTED)

    # (c) the swing
    sx, sy, sw, sh = 905, 95, 450, 300
    tag(sx - 40, 36, '(c) the seam record after a pulse')
    rect(sx, sy, sw, sh, 'none', FAINT, 1.4, 0)
    t = np.linspace(0.0, 3 * math.pi / E, 400)
    occ = 1 - (5 / 6) * np.sin(E * t) ** 2
    pts = ' '.join('%.1f,%.1f' % (sx + sw * x / t[-1], sy + sh - sh * y) for x, y in zip(t, occ))
    md.out.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="3.4"/>' % (pts, BLUE))
    for v, s in ((1.0, '1'), (1 / 6, '1/6')):
        y = sy + sh - sh * v
        line(sx, y, sx + sw, y, FAINT, 1.0, dash='4 6')
        label(sx - 10, y + 7, s, 20, MUTED, 'end')
    for k, s in enumerate(('0', 'π/E', '2π/E', '3π/E')):
        label(sx + sw * k / 3, sy + sh + 30, s, 20, MUTED, 'middle')
    label(sx + sw / 2, sy + sh + 60, "the observer's own phase", 21, MUTED, 'middle')
    f1 = clock['steps_of_equal_length']['line']
    label(sx - 40, 500, 'line read from counts:', 20, MUTED)
    label(sx + sw, 500, '%.6f ± %.6f' % (f1['frequency'], f1['standard_error']), 20, INK, 'end', 'bold')
    label(sx - 40, 527, 'level difference, fuel dial:', 20, MUTED)
    label(sx + sw, 527, '%.6f ± %.6f' % (dial['level_difference'], dial['level_difference_standard_error']), 20, INK, 'end', 'bold')
    label(sx - 40, 554, 'the law:', 20, MUTED)
    label(sx + sw, 554, '2√6 = %.6f' % (2 * E), 20, INK, 'end', 'bold')
    write('model_cage', W, H)
    return persistence


# ---------------------------------------------------------------- one run and two readings

def box(x, y, w, h, colour, fill, title, rows):
    rect(x, y, w, h, fill, colour, 2.2, 10)
    label(x + 14, y + 32, title, 23, colour, 'start', 'bold')
    for k, s in enumerate(rows):
        label(x + 14, y + 62 + 25 * k, s, 19, INK)


def run_figure():
    W, H = 1400, 735
    tag(20, 36, '(a) one run, from preparation to record')
    steps = (('prepare', INK, '#f6f8fa', ('blank memories, or', 'the calm reference', 'of the second law')),
             ('write', AMBER, LIGHT_AMBER, ('a local emitter', 'turns its fuel', 'into an object')),
             ('let it run', PURPLE, LIGHT_PURPLE, ('the law acts alone;', 'each observer', 'counts its events')),
             ('read', ORANGE, '#fff1e5', ('the two observers', 'of a seam copy one', 'answer to pointers')),
             ('agree', GREEN, LIGHT_GREEN, ('a record needs two', 'agreeing copies;', 'every outcome kept')),
             ('record', BLUE, LIGHT_BLUE, ('the agreed record;', 'its count is the', "observer's clock")))
    bw, gap = 196, 34
    for k, (title, colour, fill, rows) in enumerate(steps):
        x = 20 + k * (bw + gap)
        box(x, 60, bw, 150, colour, fill, title, rows)
        if k < len(steps) - 1:
            arrow(x + bw + 4, 135, x + bw + gap - 4, 135, INK, 2.4)

    tag(20, 272, '(b) two readings of the same pair on one cell')
    step = 18.0

    def leaves(x0, y, fill, stroke):
        for k in range(30):
            circle(x0 + step * k, y, 6, fill, stroke, 1.4)

    # the presence reading: thirty seam words, an OR tree of depth five, one public answer
    px, py = 70, 300
    label(px - 20, py + 22, 'the presence reading', 23, BLUE, 'start', 'bold')
    xs = [px + step * k for k in range(30)]
    y = py + 62
    leaves(px, y, LIGHT_BLUE, BLUE)
    label(px + step * 29 + 22, y + 7, 'thirty seam words', 19, MUTED)
    for level in range(5):
        parents = [(xs[2 * k] + xs[min(2 * k + 1, len(xs) - 1)]) / 2 for k in range((len(xs) + 1) // 2)]
        y2 = y + 42
        for k, xp in enumerate(parents):
            for c in (2 * k, 2 * k + 1):
                if c < len(xs):
                    line(xs[c], y + 6, xp, y2 - 7, BLUE, 1.3, opacity=0.7)
        for xp in parents:
            circle(xp, y2, 7 if len(parents) > 1 else 12, '#ffffff' if len(parents) > 1 else BLUE, BLUE, 1.6)
        xs, y = parents, y2
    label(xs[0] + 24, y + 7, 'OR of all thirty: one public answer', 19, INK, 'start', 'bold')
    label(px - 20, y + 60, 'the tree and the readers then run backwards: the seam words are unwritten', 19, MUTED)
    label(px - 20, y + 90, 'a positive answer acts as  X ↦ P X P', 21, BLUE, 'start', 'bold')

    # the reader with seam records
    rx = 790
    label(rx - 20, py + 22, 'the reader with seam records', 23, ORANGE, 'start', 'bold')
    leaves(rx, py + 62, '#fff1e5', ORANGE)
    for k in range(30):
        x = rx + step * k
        line(x, py + 70, x, py + 190, ORANGE, 1.2, opacity=0.6)
        rect(x - 6, py + 192, 12, 24, '#fff1e5', ORANGE, 1.2, 2)
    label(rx - 20, py + 250, 'one public record for every seam:', 19, INK, 'start', 'bold')
    label(rx - 20, py + 276, 'the pair is placed on one seam', 19, MUTED)
    label(rx - 20, py + 302 + 6, 'a positive answer acts as  X ↦ Σₑ Pₑ X Pₑ', 21, ORANGE, 'start', 'bold')
    label(W / 2, H - 16, 'Both find a pair with the same probability; the records taken afterwards follow different laws.', 21, INK, 'middle')
    write('model_run', W, H)


def main():
    if not (shutil.which('rsvg-convert') and shutil.which('gs')):
        raise SystemExit('make_paper_figures: rsvg-convert and gs must be on the PATH')
    hardware()
    stability()
    cage()
    run_figure()


if __name__ == '__main__':
    main()
