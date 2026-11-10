"""The hand-over to the observer, the field description and the roster."""
from fractions import Fraction
import math

from .canvas import (AMBER, BG, BLUE, CARD, CARD2, EDGE, FAINT, GOLD, GREEN, MONO, MUTED, ORANGE,
                     RED, SANS, SERIF, SOFT, TEAL, TEXT, VIOLET, n, text_width)
from .parts import (CW, ROOT3, SX, W, X0, X1, bell, card, header, lcg, mix, status, wave)
from .scenes_a import record_card


# ---------------------------------------------------------------- step 5
def step5(c, data, y):
    yb = header(c, y, 5, 'The hand-over', 'It moves, and nothing travels',
                'The pattern is re-established one place over. Different records carry it at every step.',
                AMBER)
    h = 484
    card(c, X0, yb, CW, h)
    half = 640
    lx, rx = X0 + 30, X1 - 30 - half
    for x, lab, col in ((lx, 'In the machine', TEAL), (rx, 'For an observer', VIOLET)):
        c.rect(x, yb + 26, half, 322, CARD2, None, 0, 16)
        c.caps(x + 22, yb + 56, lab, col, 13.5)
    mid = (lx + half + rx) / 2
    c.circle(mid, yb + 187, 19, CARD, EDGE, 1.3)
    c.arrow(mid - 8, yb + 187, mid + 9, yb + 187, SOFT, 2.2, 8)

    # machine: a staircase of hand-overs inside the cone
    cell, rows, cols = 33, 7, 15
    gx, gy = lx + 100, yb + 96
    for r in range(rows + 1):
        c.line(gx, gy + r * cell, gx + cols * cell, gy + r * cell, EDGE, 1, 0.9)
    for q in range(cols + 1):
        c.line(gx + q * cell, gy, gx + q * cell, gy + rows * cell, EDGE, 1, 0.9)
    start = 3
    apex = (gx + (start + 0.5) * cell, gy + 0.5 * cell)
    cone = [apex, (apex[0] + (rows - 0.5) * cell, gy + rows * cell), (apex[0] - (start + 0.5) * cell, gy + rows * cell),
            (gx, apex[1] + (start + 0.5) * cell)]
    c.poly([apex, (apex[0] + (rows - 0.5) * cell, gy + rows * cell),
            (gx, gy + rows * cell), (gx, apex[1] + (start + 0.5) * cell)], TEAL, op=0.07)
    c.line(apex[0], apex[1], apex[0] + (rows - 0.5) * cell, gy + rows * cell, TEAL, 1.6, 0.85, '6 6')
    c.line(apex[0], apex[1], gx, apex[1] + (start + 0.5) * cell, TEAL, 1.6, 0.85, '6 6')
    path = [start + t // 2 for t in range(rows)]
    for t, q in enumerate(path):
        x, yy = gx + q * cell, gy + t * cell
        last = t == rows - 1
        if last:
            c.glow(x + cell / 2, yy + cell / 2, 34, AMBER, 0.7)
        c.rect(x + 3, yy + 3, cell - 6, cell - 6, AMBER, None, 0, 6, round(0.3 + 0.7 * t / (rows - 1), 2))
        c.add(f'<rect x="{n(x + 1.5)}" y="{n(yy + 1.5)}" width="{cell - 3}" height="{cell - 3}" rx="7" '
              f'fill="none" stroke="#fff" stroke-width="2.2" opacity="0" class="hop" '
              f'style="animation-delay:{n(t * 0.8)}s"/>')
        if t and path[t - 1] != q:
            c.qarrow(x - cell / 2 + 4, yy - cell / 2 + 2, x + 2, yy - cell / 2 - 8, x + cell / 2 - 5,
                     yy + 5, TEXT, 1.5, 6, 0.9)
    c.arrow(gx, gy - 14, gx + 104, gy - 14, FAINT, 1.4, 7)
    c.text(gx + 112, gy - 9, 'places', 14, MUTED)
    c.arrow(gx - 16, gy, gx - 16, gy + 96, FAINT, 1.4, 7)
    c.text(gx - 22, gy + 116, 'ticks', 14, MUTED, anchor='end')
    c.text(gx + cols * cell - 6, gy + rows * cell + 22, 'fastest hand-over: one place per tick', 14,
           TEAL, 600, 'end')

    # observer: one smooth worldline inside the same cone
    ox, oy = rx + 100, yb + 96
    wd, ht = cols * cell, rows * cell
    c.line(ox, oy, ox + wd, oy, FAINT, 1.2, 0.7)
    c.line(ox, oy, ox, oy + ht, FAINT, 1.2, 0.7)
    a2 = (ox + (start + 0.5) * cell, oy + 0.5 * cell)
    c.poly([a2, (a2[0] + (rows - 0.5) * cell, oy + ht), (ox, oy + ht), (ox, a2[1] + (start + 0.5) * cell)],
           VIOLET, op=0.09)
    c.line(a2[0], a2[1], a2[0] + (rows - 0.5) * cell, oy + ht, VIOLET, 1.6, 0.85, '6 6')
    c.line(a2[0], a2[1], ox, a2[1] + (start + 0.5) * cell, VIOLET, 1.6, 0.85, '6 6')
    end = (a2[0] + (rows - 1) * cell / 2, a2[1] + (rows - 1) * cell)
    c.define('track', f'<linearGradient id="track" gradientUnits="userSpaceOnUse" x1="{n(a2[0])}" '
                      f'y1="{n(a2[1])}" x2="{n(end[0])}" y2="{n(end[1])}"><stop offset="0" stop-color="{AMBER}" '
                      f'stop-opacity=".25"/><stop offset="1" stop-color="{AMBER}"/></linearGradient>')
    c.line(a2[0], a2[1], end[0], end[1], 'url(#track)', 4)
    c.add(f'<g class="glide" style="--dx:{n(a2[0] - end[0])}px;--dy:{n(a2[1] - end[1])}px">')
    c.glow(end[0], end[1], 36, VIOLET, 0.8)
    c.glow(end[0], end[1], 17, AMBER, 0.95)
    c.circle(end[0], end[1], 7.5, AMBER, '#fff', 2)
    c.close()
    c.arrow(ox, oy - 14, ox + 104, oy - 14, FAINT, 1.4, 7)
    c.text(ox + 112, oy - 9, 'space', 14, MUTED)
    c.arrow(ox - 16, oy, ox - 16, oy + 96, FAINT, 1.4, 7)
    c.text(ox - 22, oy + 116, 'time', 14, MUTED, anchor='end')
    c.text(a2[0] + (rows - 0.5) * cell + 2, oy + ht + 22, 'candidate light cone', 14, VIOLET, 600, 'end')
    c.text(end[0] - 26, end[1] - 4, 'one object,', 14.5, TEXT, 600, 'end')
    c.text(end[0] - 26, end[1] + 15, 'one smooth path', 14.5, TEXT, 600, 'end')
    c.text(ox + wd, oy + 22, 'no signal reaches', 14, MUTED, 400, 'end')
    c.text(ox + wd, oy + 40, 'outside the cone', 14, MUTED, 400, 'end')

    c.para(lx + 4, yb + 384, half - 24,
           'Each tick, repair resolves the twist where it sits and finds it again next door. '
           'A wave crosses a lake the same way: the water goes up and down, and the shape advances.',
           17.5, SOFT, 1.46)
    c.para(rx + 4, yb + 384, half - 24,
           'A relativistic interpretation needs a common limiting speed for all observers. '
           'A finite causal bound alone does not establish this identification with light.',
           17.5, SOFT, 1.46)
    return status(c, X0, yb + h + 16, CW,
                  'Record-controlled moves carry a twist across the carrier, hop by hop. Influence '
                  'spreads inside a cone, and the bound holds uniformly as the lattice is refined.',
                  'Free motion of a defect under the repair law alone, with its own momentum and '
                  'energy law, and Lorentz symmetry of the limit.')


# ---------------------------------------------------------------- steps 6 and 7
def step67(c, data, y):
    gap = 24
    w = (CW - gap) / 2
    h = 668
    c.nodes.append((y + 36, '6', AMBER))
    c.nodes.append((y + 96, '7', AMBER))
    xa, xb = X0, X0 + w + gap
    for x, num, title in ((xa, 6, 'Position is where a record is written'),
                          (xb, 7, 'Velocity is the slope of the records')):
        card(c, x, y, w, h)
        c.caps(x + 32, y + 44, f'The hand-over  ·  step {num}', AMBER, 13.5)
        c.text(x + 32, y + 86, title, 29, TEXT, 700, ls=-0.3, maxw=w - 60)

    # --- 6: amplitude spread, detectors, one click, many runs
    x = xa
    px, pw, py = x + 28, w - 56, y + 112
    c.rect(px, py, pw, 258, CARD2, None, 0, 16)
    nd, dw = 13, 30
    dx0 = px + 96
    base = py + 110
    hump = [(dx0 + nd * dw * t / 80, base - 12 - 60 * bell(t / 80, 0.5, 0.17)) for t in range(81)]
    c.define('hump', f'<linearGradient id="hump" x1="0" x2="0" y1="0" y2="1"><stop offset="0" '
                     f'stop-color="{AMBER}" stop-opacity=".5"/><stop offset="1" stop-color="{AMBER}" '
                     f'stop-opacity="0"/></linearGradient>')
    c.path('M' + ' L'.join(f'{n(a)},{n(b)}' for a, b in hump) + f' L{n(hump[-1][0])},{n(base - 12)} '
           f'L{n(hump[0][0])},{n(base - 12)} Z', fill='url(#hump)')
    c.curve(hump, AMBER, 2.2)
    c.text(px + 20, py + 30, 'before the record: spread over places', 14, AMBER, 600)
    hit = 7
    rng = lcg(41)
    counts = [round(8.4 * bell((i + 0.5) / nd, 0.5, 0.17) + (next(rng) - 0.5) * 1.2) for i in range(nd)]
    for i in range(nd):
        xx = dx0 + i * dw
        lit = i == hit
        if lit:
            c.glow(xx + dw / 2, base + 9, 30, AMBER, 0.8, cls='pulse')
        c.rect(xx + 2.5, base - 4, dw - 5, 26, AMBER if lit else '#17233b', '#fff' if lit else EDGE,
               1.6 if lit else 1.1, 6)
        for k in range(max(0, counts[i])):
            c.circle(xx + dw / 2, base + 38 + k * 10.5, 3.5, AMBER, op=0.6)
    c.text(dx0 - 12, base + 14, 'detectors', 13.5, MUTED, anchor='end')
    c.text(dx0 - 12, base + 46, 'many runs', 13.5, MUTED, anchor='end')
    rx0, rw0 = px + pw - 132, 116
    c.rect(rx0, py + 74, rw0, 104, '#0b1322', AMBER, 1.3, 12)
    c.caps(rx0 + 14, py + 98, 'one record', AMBER, 11.5, ls=1.6)
    for k, line in enumerate(('tick 041', 'detector 07', 'count 1')):
        c.text(rx0 + 14, py + 122 + k * 19, line.replace(' ', ' '), 13.5, SOFT, family=MONO)
    c.line(dx0 + hit * dw + dw - 2, base + 9, rx0 - 2, base + 9, AMBER, 1.3, 0.9, '4 4')
    c.text(rx0 + rw0 / 2, py + 204, 'one run,', 13.5, MUTED, 400, 'middle')
    c.text(rx0 + rw0 / 2, py + 222, 'one click', 13.5, MUTED, 400, 'middle')
    yy = c.para(x + 32, y + 408, w - 64,
                'In this quantum preparation the amplitude is spread over places. A position '
                'detector writes a definite record: this detector, this tick, a count of one. '
                'The drawing illustrates the readout, rather than a hidden classical path.',
                17, SOFT, 1.46)
    yy = c.para(x + 32, yy + 8, w - 64,
                'The detector’s law fixes the odds of each click: the squared amplitude. Over many '
                'runs the clicks pile up under the hump. Observers who compare notes find the same click.',
                17, SOFT, 1.46)
    # how a detector senses the twist (certified two-sector record law)
    law = data['attachment']['same_time_scalar_records']
    iy = y + h - 84
    c.rect(x + 28, iy, w - 56, 62, CARD2, EDGE, 1.1, 14)
    c.caps(x + 46, iy + 26, 'A detector coupled', TEAL, 11.5, ls=1.6)
    c.caps(x + 46, iy + 45, 'to the twist fires in', TEAL, 11.5, ls=1.6)
    for k, (key, lab) in enumerate((('1', 'twist absent'), ('-1', 'twist present'))):
        frac = Fraction(law[key]['occupied'])
        bx = x + 282 + k * 196
        for j in range(frac.denominator):
            on = j < frac.numerator
            c.circle(bx + j * 19, iy + 22, 6.5, AMBER if on else 'none', None if on else FAINT, 1.4)
        c.text(bx - 7, iy + 48, f'{frac.numerator} of {frac.denominator} runs, {lab}', 13.5, SOFT)

    # --- 7: the slope of the records and the wavelength
    x = xb
    px, py = x + 28, y + 112
    c.rect(px, py, pw, 258, CARD2, None, 0, 16)
    ox, oy, aw, ah = px + 64, py + 28, 250, 196
    c.line(ox, oy, ox, oy + ah, FAINT, 1.3)
    c.line(ox, oy + ah, ox + aw, oy + ah, FAINT, 1.3)
    c.text(ox + aw, oy + ah + 22, 'place', 13.5, MUTED, anchor='end')
    c.text(ox - 10, oy + 10, 'tick', 13.5, MUTED, anchor='end')
    recs = [(ox + 34 + k * 52, oy + ah - 30 - k * 40) for k in range(4)]
    c.line(recs[0][0] - 26, recs[0][1] + 20, recs[-1][0] + 28, recs[-1][1] - 21.5, VIOLET, 2, 0.9, '6 5')
    for k, (a, b) in enumerate(recs):
        c.glow(a, b, 16, AMBER, 0.8)
        c.circle(a, b, 6, AMBER, '#fff', 1.6)
    a, b = recs[1], recs[2]
    c.line(a[0], a[1], b[0], a[1], TEXT, 1.4, 0.9)
    c.line(b[0], a[1], b[0], b[1], TEXT, 1.4, 0.9)
    c.text((a[0] + b[0]) / 2, a[1] + 18, 'places', 12.5, SOFT, anchor='middle')
    c.text(b[0] + 8, (a[1] + b[1]) / 2 + 4, 'ticks', 12.5, SOFT)
    c.text(ox + 16, oy + 12, 'velocity = places per tick', 14, VIOLET, 600)
    # momentum as wavelength
    wx = px + 372
    c.line(wx - 24, py + 24, wx - 24, py + 234, EDGE, 1.2)
    ww = pw - 372 - 22
    c.text(wx, py + 44, 'slow phase variation', 14, TEAL, 600)
    wave(c, wx, py + 84, ww, 18, 2.5, TEAL, 2.4)
    c.text(wx, py + 144, 'fast phase variation', 14, AMBER, 600)
    wave(c, wx, py + 184, ww, 18, 7, AMBER, 2.4)
    c.text(wx, py + 234, 'phase turns per place', 13.5, MUTED)
    yy = c.para(x + 32, y + 408, w - 64,
                'Successive records line up. Places per tick along that line is the velocity.',
                17, SOFT, 1.46)
    yy = c.para(x + 32, yy + 8, w - 64,
                'A Fourier detector reads how phase turns from place to place. Calling that label '
                'conserved momentum requires a translation symmetry of the same field dynamics.',
                17, SOFT, 1.46)
    c.para(x + 32, yy + 8, w - 64,
           'Place and momentum are two detector settings. A sharp answer to one leaves a broad '
           'spread of answers to the other.', 17, SOFT, 1.46)
    # sharp place, broad momentum
    iy = y + h - 84
    c.rect(x + 28, iy, w - 56, 62, CARD2, EDGE, 1.1, 14)
    c.caps(x + 46, iy + 25, 'Sharp place', TEAL, 11.5, ls=1.6)
    c.caps(x + 46, iy + 46, 'broad momentum', AMBER, 11.5, ls=1.6)
    bx = x + 226
    for i in range(9):
        hh = 30 * bell(i / 8, 0.5, 0.09)
        c.rect(bx + i * 14, iy + 46 - hh, 10, max(hh, 2), TEAL, None, 0, 2)
    c.arrow(bx + 142, iy + 31, bx + 176, iy + 31, FAINT, 1.6, 7)
    for i in range(9):
        hh = 30 * bell(i / 8, 0.5, 0.42)
        c.rect(bx + 196 + i * 14, iy + 46 - hh, 10, max(hh, 2), AMBER, None, 0, 2)
    c.text(x + w - 46, iy + 37, 'same state', 13.5, MUTED, anchor='end')
    return status(c, X0, y + h + 16, CW,
                  'Complete position records and conserved crystal-momentum records of one finite '
                  'field, with exact Born weights and a finite uncertainty relation. The selected '
                  'source has a complete spatial Fourier readout. A detector coupled to the '
                  'protected charge has a different click law in each sector.',
                  'The position detector of a moving defect in the selected program. Translation '
                  'symmetry of the same Hamiltonian and an energy-dispersion law that attach its '
                  'Fourier labels to physical momentum.')


# ---------------------------------------------------------------- step 8
def step8(c, data, y):
    yb = header(c, y, 8, 'The hand-over', 'Mass is the price of standing still',
                'The proposed link: turn a standing repair cost into a measured rest energy.', AMBER)
    h = 452
    card(c, X0, yb, CW, h)
    # left: the core that keeps costing
    cx, cy = X0 + 150, yb + 152
    for r, op in ((104, 0.14), (78, 0.24), (52, 0.4)):
        c.circle(cx, cy, r, 'none', AMBER, 1.6, op, '5 6')
    c.glow(cx, cy, 62, AMBER, 0.85, cls='pulse')
    c.circle(cx, cy, 11, AMBER, '#fff', 2)
    for k in range(8):
        ang = math.radians(45 * k + 22)
        c.arrow(cx + 96 * math.cos(ang), cy + 96 * math.sin(ang), cx + 34 * math.cos(ang),
                cy + 34 * math.sin(ang), TEAL, 1.5, 7, 0.75)
    c.text(cx, yb + 286, 'repair arrives at every tick', 14.5, TEAL, 600, 'middle')
    c.para(X0 + 300, yb + 62, 250,
           'Repair meets the mismatch at the core, resolves it, and finds it one seam over.',
           17, SOFT, 1.46)
    c.para(X0 + 300, yb + 178, 250,
           'To identify {a|rest energy}, derive an energy gap at zero physical momentum on '
           'the observer’s clock.',
           17, SOFT, 1.46)
    c.rect(X0 + 32, yb + 318, 520, 98, CARD2, EDGE, 1.1, 16)
    c.text(X0 + 292, yb + 360, 'mass = rest energy / c²'.replace(' ', '\u00a0'), 23, TEXT, 700, 'middle',
           family=MONO)
    c.text(X0 + 292, yb + 392, 'the physical target, once rest energy is identified', 15, MUTED, 400, 'middle')

    # middle: energy and momentum
    gx, gy, gw, gh = X0 + 652, yb + 48, 340, 262
    c.line(gx, gy, gx, gy + gh, FAINT, 1.4)
    c.line(gx, gy + gh, gx + gw, gy + gh, FAINT, 1.4)
    c.text(gx - 12, gy + 12, 'energy', 14, MUTED, anchor='end')
    c.text(gx + gw, gy + gh + 24, 'momentum', 14, MUTED, anchor='end')
    m, pmax = 0.62, 1.2
    sx, sy = gw / pmax, gh / 1.42
    hyp = [(gx + sx * p, gy + gh - sy * math.sqrt(p * p + m * m)) for p in (pmax * k / 60 for k in range(61))]
    c.path('M' + ' L'.join(f'{n(a)},{n(b)}' for a, b in hyp) + f' L{n(gx + gw)},{n(gy + gh - sy * pmax)} '
           f'L{n(gx)},{n(gy + gh)} Z', fill=AMBER, op=0.08)
    c.curve(hyp, AMBER, 3.2)
    c.line(gx, gy + gh, gx + sx * pmax, gy + gh - sy * pmax, TEAL, 2.4, None, '7 6')
    iy = gy + gh - sy * m
    c.glow(gx, iy, 22, AMBER, 0.9)
    c.circle(gx, iy, 6, AMBER, '#fff', 1.8)
    c.text(gx + 16, iy - 58, 'rest energy mc²', 15, AMBER, 700)
    c.text(gx + 16, iy - 39, 'the cost at zero', 13.5, MUTED)
    c.text(gx + 16, iy - 22, 'momentum', 13.5, MUTED)
    c.text(gx + 126, gy + 16, 'relativistic target law', 14.5, AMBER, 700)
    c.text(gx + gw - 6, gy + gh - 54, 'massless target', 14.5, TEAL, 700, 'end')
    c.text(gx + gw - 6, gy + gh - 35, 'no cost at rest', 13.5, MUTED, 400, 'end')
    c.rect(X0 + 620, yb + 346, 404, 70, CARD2, EDGE, 1.1, 16)
    c.text(X0 + 822, yb + 390, 'E² = (pc)² + (mc²)²'.replace(' ', ' '), 23, TEXT, 700, 'middle',
           family=MONO)

    # right: with and without a core
    rx, rw = X0 + 1090, CW - 1090 - 32
    c.line(rx - 26, yb + 36, rx - 26, yb + h - 36, EDGE, 1.2)
    c.caps(rx, yb + 58, 'The photon target', TEAL, 13)
    wave(c, rx, yb + 92, rw, 13, 5, TEAL, 2.2, envelope=lambda t: bell(t, 0.5, 0.26))
    c.para(rx, yb + 136, rw,
           'For a massless relativistic mode, E = c|p| and the wave speed is c. A source ripple '
           'must satisfy this law and the photon’s other properties before it is called light.',
           16.5, MUTED, 1.45)
    c.caps(rx, yb + 290, 'Measured comparison', AMBER, 13)
    masses = data['masses']
    ratio = float(masses['top']['value']) / float(masses['electron']['value'])
    c.para(rx, yb + 318, rw,
           f'Measured masses span a wide range: a top quark has about '
           f'{round(ratio, -3):,.0f} electron masses. These comparison values do not set '
           f'the source energies.', 16.5, MUTED, 1.45)
    return status(c, X0, yb + h + 16, CW,
                  'The law E² = (pc)² + (mc²)² holds for the declared free field action in its '
                  'stated limit. The interacting primary source has a rigorously isolated '
                  'gauge transition near 63 source energy units in a rare empty-matter sector.',
                  'An excitation that survives with a fixed energy band, its momentum shell, '
                  'and the rest mass of a defect computed from its repair cost. No mass in GeV '
                  'is predicted here.')


# ---------------------------------------------------------------- step 9
def step9(c, data, y):
    yb = header(c, y, 9, 'The description', 'Observers compress the records into fields',
                'Observers inside hold records. A field is how they summarize what they agree on.',
                VIOLET)
    h = 472
    card(c, X0, yb, CW, h)
    pw, ph, gap = 296, 214, 62
    x0 = X0 + 28
    labels = (('Records', 'finite clicks, each one committed', AMBER),
              ('Agreement', 'observers compare their copies', GREEN),
              ('Field', 'one smooth summary, with error bounds', VIOLET),
              ('Equations', 'predict the next records', TEXT))
    xs = [x0 + k * (pw + gap) for k in range(4)]
    for k, (lab, sub, col) in enumerate(labels):
        c.rect(xs[k], yb + 28, pw, ph, CARD2, None, 0, 16)
        c.circle(xs[k] + 22, yb + ph + 58, 12, 'none', col, 1.6)
        c.text(xs[k] + 22, yb + ph + 63, str(k + 1), 13.5, col, 700, 'middle')
        c.caps(xs[k] + 44, yb + ph + 63, lab, col, 13.5)
        c.text(xs[k] + 4, yb + ph + 90, sub, 15, MUTED, maxw=pw + 40)
        if k < 3:
            c.arrow(xs[k] + pw + 14, yb + 28 + ph / 2, xs[k] + pw + gap - 14, yb + 28 + ph / 2, SOFT, 2.2, 9)

    def profile(t):
        return 0.5 + 0.34 * math.sin(2 * math.pi * (t * 1.5 + 0.08)) * (0.55 + 0.45 * math.sin(math.pi * t))

    # 1 records
    rng = lcg(7)
    dots = []
    while len(dots) < 150:
        t, v = next(rng), next(rng)
        if v < profile(t) ** 2 * 1.5:
            dots.append((t, next(rng)))
    for t, v in dots:
        c.circle(xs[0] + 22 + t * (pw - 44), yb + 50 + v * (ph - 44), 2.7, AMBER, op=0.85)
    # 2 agreement
    ax, ay = xs[1] + pw / 2, yb + 28 + ph / 2
    for s, col, name in ((-1, TEAL, 'A'), (1, BLUE, 'B')):
        c.circle(ax + s * 52, ay, 84, col, col, 1.6, 0.1)
        c.circle(ax + s * 52, ay, 84, 'none', col, 1.6, 0.9)
        c.text(ax + s * 112, ay - 62, name, 17, col, 700, 'middle')
    for k, (dx, dy) in enumerate(((0, -40), (-10, -8), (12, 18), (-2, 46), (10, -22))):
        c.circle(ax + dx, ay + dy, 4.2, GREEN)
    for s in (-1, 1):
        for dx, dy in ((70, -26), (88, 22), (58, 52), (100, -8)):
            c.circle(ax + s * dx, ay + dy, 3, AMBER, op=0.55)
    c.circle(ax, ay + 78, 12, BG, GREEN, 1.6)
    c.text(ax, ay + 83.5, '✓', 15, GREEN, 700, 'middle')
    # 3 field
    fx, fy, fw, fh = xs[2] + 20, yb + 46, pw - 40, ph - 40
    c.define('fieldfill', f'<linearGradient id="fieldfill" x1="0" x2="0" y1="0" y2="1"><stop offset="0" '
                          f'stop-color="{VIOLET}" stop-opacity=".45"/><stop offset="1" stop-color="{VIOLET}" '
                          f'stop-opacity="0"/></linearGradient>')
    curve = [(fx + fw * k / 90, fy + fh - 14 - (fh - 40) * profile(k / 90) ** 2 * 1.25) for k in range(91)]
    c.path('M' + ' L'.join(f'{n(a)},{n(b)}' for a, b in curve) + f' L{n(fx + fw)},{n(fy + fh - 8)} '
           f'L{n(fx)},{n(fy + fh - 8)} Z', fill='url(#fieldfill)')
    c.curve(curve, VIOLET, 3)
    for t, v in dots[:46]:
        yy = fy + fh - 14 - (fh - 40) * profile(t) ** 2 * 1.25
        c.circle(fx + t * fw, yy + (v - 0.5) * 20, 2.2, AMBER, op=0.55)
    # 4 equations
    for k, (eq, lab) in enumerate((('(iγ·∂ − m) ψ = 0', 'matter'), ('∂·F = J', 'light'),
                                   ('E² = (pc)² + (mc²)²', 'energy and mass'))):
        ey = yb + 74 + k * 62
        c.text(xs[3] + 24, ey, eq, 19.5, TEXT, 700, family=MONO, maxw=pw - 40)
        c.text(xs[3] + 24, ey + 21, lab, 13.5, MUTED)
    # the loop back
    by = yb + ph + 112
    c.path(f'M{n(xs[3] + pw / 2)},{by} v 16 q 0,12 -12,12 H{n(xs[0] + pw / 2 + 12)} q -12,0 -12,-12 v -6',
           VIOLET, 1.8, op=0.9, dash='6 6')
    c.head(xs[0] + pw / 2, by + 2, -math.pi / 2, VIOLET, 10)
    c.text((xs[0] + xs[3] + pw) / 2, by + 50, 'predict, measure again, compare', 14.5, VIOLET, 600, 'middle')
    ty = by + 86
    c.para(X0 + 32, ty, 650,
           'When a controlled field description predicts the agreed records, its excitations '
           'can be studied as particle candidates. The program must establish that description '
           'and its errors on the same source and state.', 17, SOFT, 1.46)
    c.para(X0 + 730, ty, 640,
           'The program outputs finite events and records. The continuum belongs to the '
           'description, with a stated domain and error bounds. No extra continuum universe is '
           'produced underneath.', 17, SOFT, 1.46)
    return status(c, X0, yb + h + 16, CW,
                  'The selected program has one observable algebra, a stationary limiting state, '
                  'positive energy and convergent finite detector histories. Record copies and '
                  'observer dictionaries glue into one shared event description.',
                  'A nontrivial interacting relativistic field in that limit. Agreement between '
                  'observers fixes neither the Born weights nor a relativistic field limit.')


# ---------------------------------------------------------------- dictionary
ROWS = (('a twist that every loop around it reads', 'a particle'),
        ('the number on the boundary loop', 'charge'),
        ('a frame twist that needs two full turns', 'spin ½'),
        ('the standing repair cost of the core', 'rest energy, mass'),
        ('a pattern re-established one place over', 'motion, momentum'),
        ('the fastest hand-over, one place per tick', 'the speed of light'),
        ('a charge-conjugate source channel', 'antiparticle candidate'),
        ('a ripple with no core', 'photon, gluon: massless'),
        ('a committed record', 'a measurement outcome'),
        ('the shortest summary of agreed records', 'a quantum field'))


def dictionary(c, data, y):
    c.caps(X0, y + 10, 'The dictionary', VIOLET, 14)
    c.text(X0, y + 60, 'One object, two vocabularies', 40, TEXT, 700, ls=-0.4)
    yb = c.para(X0, y + 100, 1300, 'The left column is what the program holds. The right column '
                                   'is its proposed physical reading. These are interpretation '
                                   'targets, not ten established source-to-particle identifications.',
                                   22, MUTED, 1.4) + 6
    rh, gap = 52, 24
    w = (CW - gap) / 2
    h = 5 * rh + 74
    split = 412
    for col in range(2):
        x = X0 + col * (w + gap)
        card(c, x, yb, w, h)
        c.caps(x + 52, yb + 40, 'In the machine', TEAL, 12.5)
        c.caps(x + split + 44, yb + 40, 'Proposed reading', VIOLET, 12.5)
        for k in range(5):
            left, right = ROWS[col * 5 + k]
            ry = yb + 58 + k * rh
            c.line(x + 24, ry, x + w - 24, ry, EDGE, 1.1)
            c.circle(x + 36, ry + rh / 2, 4.5, TEAL)
            c.text(x + 52, ry + rh / 2 + 6, left, 16.5, SOFT, maxw=split - 62)
            c.arrow(x + split, ry + rh / 2, x + split + 30, ry + rh / 2, FAINT, 1.5, 7)
            c.text(x + split + 44, ry + rh / 2 + 6, right, 16.5, VIOLET, 700, maxw=w - split - 60)
    return yb + h


# ---------------------------------------------------------------- roster
def _mass(row):
    if row is None or row.get('value') is None:
        return None
    v = float(row['value'])
    if v == 0:
        return None
    unit = 'GeV'
    if v < 1:
        v, unit = v * 1000, 'MeV'
    sig = 4 if (unit == 'GeV' or v >= 100) else 3
    dec = max(0, sig - 1 - math.floor(math.log10(v)))
    return f'{v:.{dec}f} {unit}'


def tile(c, x, y, sym, name, mass, charge, col, w=106, h=92):
    c.rect(x, y, w, h, CARD2, col, 1.2, 14)
    c.text(x + w - 9, y + 18, charge, 11.5, MUTED, 400, 'end', family=MONO)
    c.text(x + w / 2, y + 41, sym, 26, col, 700, 'middle')
    c.text(x + w / 2, y + 60, name, 12.5, SOFT, 400, 'middle', maxw=w - 6)
    c.text(x + w / 2, y + 79, mass.replace(' ', ' '), 12, MUTED, 400, 'middle', family=MONO,
           maxw=w - 4)


def roster(c, data, y):
    c.caps(X0, y + 10, 'The roster', VIOLET, 14)
    c.text(X0, y + 60, 'The named particles to identify', 40, TEXT, 700, ls=-0.4)
    yb = c.para(X0, y + 100, 1300,
                'The Standard Model is the comparison target. Masses below are measured inputs, '
                'not source predictions. The small number is the electric charge.',
                22, MUTED, 1.4) + 6
    h = 392
    card(c, X0, yb, CW, h)
    m = data['masses']
    tw, th, g = 106, 92, 10
    groups = (
        ('Quarks', GOLD, 'Physical targets: spin ½, color charge and fractional electric charge. '
                         'Source weights suggest candidates; flavor is unassigned.',
         (('u', 'up', 'up', '+2/3'), ('c', 'charm', 'charm', '+2/3'), ('t', 'top', 'top', '+2/3'),
          ('d', 'down', 'down', '−1/3'), ('s', 'strange', 'strange', '−1/3'),
          ('b', 'bottom', 'bottom', '−1/3')), 3),
        ('Leptons', AMBER, 'Spin ½, no color. Flavor is unassigned in the source; '
                          'neutrino flavor labels are not mass eigenstates.',
         (('e', 'electron', 'electron', '−1'), ('μ', 'muon', 'muon', '−1'), ('τ', 'tau', 'tau', '−1'),
          ('νe', 'e neutrino', 'electron_neutrino', '0'), ('νμ', 'μ neutrino', 'muon_neutrino', '0'),
          ('ντ', 'τ neutrino', 'tau_neutrino', '0')), 3),
        ('Force carriers', TEAL, 'Physical targets: spin 1. The Higgs mechanism gives W and Z '
                                 'mass while preserving the massless photon.',
         (('γ', 'photon', 'photon', '0'), ('g', 'gluon', 'gluon', '0'), ('W', 'W boson', 'W_plus', '±1'),
          ('Z', 'Z boson', 'Z', '0')), 2),
        ('Background', VIOLET, 'The Higgs target has spin 0. Its couplings and vacuum response '
                              'must follow from the same source.',
         (('H', 'Higgs', 'Higgs', '0'),), 1))
    widths = [max(per * tw + (per - 1) * g, 232) for *_, per in groups]
    space = (CW - 72 - sum(widths)) / (len(groups) - 1)
    x = X0 + 36
    for (title, col, note, items, per), gw in zip(groups, widths):
        c.caps(x, yb + 46, title, col, 13.5)
        for k, (sym, name, key, charge) in enumerate(items):
            mass = _mass(m.get(key))
            if mass is None:
                mass = 'not fixed' if 'neutrino' in key else 'massless'
            tile(c, x + (k % per) * (tw + g), yb + 62 + (k // per) * (th + g), sym, name, mass, charge, col)
        c.para(x, yb + 62 + 2 * (th + g) + 26, gw, note, 15, MUTED, 1.42)
        x += gw + space
    return status(c, X0, yb + h + 16, CW,
                  'Under stated premises the twelve-port response carries the algebra '
                  'u(1) ⊕ su(2) ⊕ su(3), and its fifteen matter states carry the charges of one '
                  'Standard Model family.',
                  'Which defect is which particle. Spin, statistics, three families and the '
                  'measured masses on one construction. No row above is an identified OPH defect.')


# ---------------------------------------------------------------- evidence
def evidence(c, data, y):
    candidates, gauge = data['candidates'], data['gauge']
    collision = candidates['collision_status']
    contact = collision['separate_contact_extension']
    spectrum = gauge['interacting_source_spectrum']
    ratio = candidates['finite_sector_ratio_control']['result']
    dressed = data['identity_tests']['calculated_source_results']['charged']['original_compensated_dressed_energy_window']
    dressed_probability = dressed['original_dressed_recorded_simple_strict_lower']
    dressed_interval = '[' + ', '.join(dressed['reported_transfer_interval']) + ']'
    c.caps(X0, y + 10, 'Evidence behind the picture', TEAL, 14)
    c.text(X0, y + 60, 'Which links are computed?', 40, TEXT, 700, ls=-0.4)
    yb = c.para(X0, y + 100, 1320,
                'The picture is a proposed physical interpretation. These linked receipts tell '
                'us what the specified finite programs actually establish, including their '
                'source domains and detector failures.', 22, MUTED, 1.4) + 6
    gap, h = 24, 344
    w = (CW - gap) / 2
    low, high = spectrum['gap_interval_decimal_display_only']
    a, b = collision['adjacent_fused_energy_interval']
    adjacent = f'[{a}, {b}]'
    simple_ratio = float(Fraction(ratio['all_three_family_simple_ratio_upper']))
    pulse_count = contact['native_pulses']
    signal = contact['executed_interacting_echo']['pointer_probability_interval']
    signal_text = f'{float(Fraction(signal[0])):.4f} to {float(Fraction(signal[1])):.4f}'
    lower = contact['exact_ideal_word_signal_certificate']['ideal_word_interaction_signal_lower']
    cards = (
        ('Candidate map', AMBER,
         f'{len(candidates["target_rows"])} named targets; no identified species',
         'The fifteen internal gauge states match one fermion family’s weights. The positive '
         'charge-six component has charged-lepton-conjugate weights. Three lepton flavors are '
         'not supplied by copying that row. Neutral three-body source singlets are composite '
         'hypotheses, not electron or neutrino identifications.',
         'Protected labels, source excitations and named targets are separate tables.',
         'Inspect the candidate correspondence',
         '../../model/standard_model_source_map/candidates/CANDIDATES.md'),
        ('An isolated source transition', TEAL,
         f'{float(low):.6f} < energy gap < {float(high):.6f}',
         'This exact interacting gauge line uses source energy units. Its certified record '
         'comes from a rare empty-matter sector. A complete loop menu couples to six scalar, '
         'vector and diagonal-tensor cube multiplets. Even scalar and vector probes cannot '
         'create the isolated odd-scalar line from an even-scalar input.',
         'Full-complement bounds control these intervals. Static fusion is not a moving collision.',
         'Read the six-multiplet detector proof',
         '../../model/particle_candidate_bosons/PROOF.md'),
        ('An interaction experiment was executed', TEAL,
         f'{pulse_count:,} native pulses; all outcomes retained',
         f'A four-mode contact echo gives pointer probability {signal_text}; the free control '
         'gives zero within the stated numerical tolerance. Its exact ideal-word signal is '
         f'greater than {float(Fraction(lower)):.7f}. This is a separate supplied contact '
         'action and preparation. A 32-row moving-packet control is also executed.',
         'Primary-source moving collisions and the separate fine-grid packet run are unexecuted.',
         'Inspect the complete collision controls',
         '../../model/standard_model_source_map/candidates/CANDIDATES.md#finite-energies-ratios-and-collisions'),
        ('The mass comparison is a test', RED,
         'Source-derived GeV mass table: empty',
         f'The three computed neutral threshold families have any common-scale ratio below '
         f'{simple_ratio:.2f} at the smallest source member. They cannot be the electron, muon '
         'and tau hierarchy, whose muon/electron ratio is about 206.77. This excludes that '
         'assignment, not other source sectors or capacities.',
         'Required: physical-vacuum channels, a momentum/energy law, controlled capacity errors '
         'and one declared energy unit. Unit conversion alone cannot change mass ratios.',
         'Read the finite ratio exclusion',
         '../../model/standard_model_source_map/candidates/FINITE_SECTOR_RATIOS.md'),
        ('Testing the charged family', AMBER,
         'A dressed charged energy channel',
         'Pairing source components gives a scalar with Higgs-conjugate weights. Its '
         'dressed charged channel overlaps the bare-hole channel. A complete detector '
         f'records pair-energy transfers in {dressed_interval} source units with original '
         f'probability above {dressed_probability}.',
         'The compensating endpoint is retained. A Higgs phase, physical flavor and '
         'single-particle rest mass require their own identification.',
         'Inspect the charged-family calculation',
         '../../model/particle_candidate_charged/PROOF.md'),
        ('Testing neutral composites', TEAL,
         'Six source identities; seven records',
         'A conserved occupation counter records the change caused by an intervention. '
         'Its values distinguish the three composite families and their conjugates, plus '
         'failure. It cancels the unknown initial value and preserves the original energy '
         'and reference histories.',
         'These are source identities. Neutrino or neutron names need physical currents '
         'and propagation. Multiple excitations can alias outside this menu.',
         'Inspect all candidate identity tests',
         '../../model/particle_candidate_confirmation/RESULTS.md'))
    for i, (title, col, headline, body, scope, label, href) in enumerate(cards):
        x, yy = X0 + (i % 2) * (w + gap), yb + (i // 2) * (h + gap)
        card(c, x, yy, w, h)
        c.circle(x + 30, yy + 31, 6, col)
        c.caps(x + 46, yy + 36, title, col, 12.5)
        c.text(x + 28, yy + 77, headline, 22, TEXT, 700, maxw=w - 56)
        end = c.para(x + 28, yy + 110, w - 56, body, 17, SOFT, 1.42)
        c.para(x + 28, end + 10, w - 56, scope, 15, MUTED, 1.42)
        c.link(x + 28, yy + h - 25, label, href, 15)
    return yb + 3 * h + 2 * gap


# ---------------------------------------------------------------- standing
NODE_LABELS = (('Program, reference', 'and clock'), ('State and', 'boundary charges'),
               ('Local records', 'and causality'), ('One state,', 'one limit'),
               ('Relativistic', 'interacting field'), ('Matter, species', 'and defects'),
               ('Spectrum,', 'decays'), ('Independent', 'comparison'),
               ('Yang–Mills', 'gap'), ('Geometry', 'and gravity'))


def standing(c, data, y):
    c.caps(X0, y + 10, 'Where the program stands', TEAL, 14)
    c.text(X0, y + 60, 'Ten requirements, one construction', 40, TEXT, 700, ls=-0.4)
    yb = c.para(X0, y + 100, 1260,
                'Every claim about fields and particles runs through one dependency map on one '
                'program, state, reference and clock.', 22, MUTED, 1.4) + 6
    h = 372
    card(c, X0, yb, CW, h)
    nodes = data['nodes']
    step = 162
    x0, ny = X0 + (CW - 7 * step) / 2, yb + 92
    pos = [(x0 + k * step, ny) for k in range(8)] + [(x0 + 4 * step - 84, ny + 172), (x0 + 4 * step + 84, ny + 172)]
    done = [node['status'] == 'finite_constructed' for node in nodes]
    for k in range(7):
        solid = done[k] and done[k + 1]
        c.arrow(pos[k][0] + 30, ny, pos[k + 1][0] - 30, ny, TEAL if solid else AMBER, 2, 8,
                None if solid else 0.85, None if solid else '5 5')
    for k in (8, 9):
        c.arrow(pos[4][0] + (-16 if k == 8 else 16), ny + 80, pos[k][0] + (12 if k == 8 else -12),
                pos[k][1] - 30, FAINT, 1.6, 7, None, '4 5')
    for k, (px, py) in enumerate(pos):
        col = TEAL if done[k] else (AMBER if k < 8 else MUTED)
        if done[k]:
            c.glow(px, py, 40, TEAL, 0.5)
            c.circle(px, py, 23, TEAL)
            c.text(px, py + 7, str(k + 1), 20, BG, 800, 'middle')
        else:
            c.circle(px, py, 23, CARD2, col, 2, dash='4.5 3.5')
            c.text(px, py + 7, str(k + 1), 20, col, 800, 'middle')
        for j, line in enumerate(NODE_LABELS[k]):
            c.text(px, py + 46 + j * 18, line, 14, SOFT if k < 8 else MUTED, 600 if j == 0 else 400,
                   'middle', maxw=step - 4)
    # where this atlas attaches
    ax = pos[5][0]
    c.text(ax + 84, ny + 130, 'the defect chain of this atlas enters here:', 14, AMBER, 600)
    c.text(ax + 84, ny + 149, 'creation, protection and field attachment', 14, AMBER, 600)
    c.qarrow(ax + 76, ny + 132, ax + 24, ny + 130, ax + 10, ny + 92, AMBER, 1.5, 7)
    ly = yb + h - 92
    c.circle(X0 + 46, ly - 5, 8, TEAL)
    c.text(X0 + 66, ly, 'finite construction, computed', 15, SOFT)
    c.circle(X0 + 46, ly + 23, 8, 'none', AMBER, 1.8, dash='3.5 3')
    c.text(X0 + 66, ly + 28, 'physical proof required', 15, SOFT)
    c.circle(X0 + 46, ly + 51, 8, 'none', MUTED, 1.8, dash='3.5 3')
    c.text(X0 + 66, ly + 56, 'separate branch', 15, SOFT)
    c.link(X1 - 36, yb + h - 28, 'The dependency map with its acceptance conditions',
           '../../model/particle_joint_reconstruction/DEPENDENCIES.md', 15, anchor='end')
    return yb + h


# ---------------------------------------------------------------- footer
LINKS = (('Book, chapter 6: the leftover', '../../book/07-chapter-06.md'),
         ('Book, chapter 28: what a particle is', '../../book/29-chapter-28.md'),
         ('Book, chapter 35: what anything is made of', '../../book/36-chapter-35.md'),
         ('Particle paper', '../../paper/deriving_the_particle_zoo_from_observer_consistency.pdf'),
         ('Protected charge: every repair outcome', '../../model/defect_protection/README.md'),
         ('Creating and transporting a twist', '../../model/particle_observer_outputs/HOLONOMY_TRANSPORT.md'),
         ('A protected label read by a detector', '../../model/particle_defect_attachment/PROOF.md'),
         ('Position and crystal-momentum records', '../../model/particle_joint_translation_readout/README.md'),
         ('One source, state, reference and clock', '../../model/particle_common_contract/README.md'),
         ('Candidate map and collision evidence', '../../model/standard_model_source_map/candidates/CANDIDATES.md'),
         ('Candidate-by-candidate identity tests', '../../model/particle_candidate_confirmation/RESULTS.md'),
         ('Isolated gauge transition', '../../model/particle_primary_reduced_sector/PROOF.md'),
         ('Source charge and spin records', '../../model/particle_primary_regional_fermion/PROOF.md'),
         ('Equal work need not preserve finite energy', '../../model/particle_primary_equal_work_control/PROOF.md'),
         ('The executable-model hypothesis', '../OPH_EXECUTABLE_HYPOTHESIS.md'),
         ('Labelled source families and chronology', '../../model/particle_candidate_confirmation/CATALOGUE.md'))


def footer(c, data, y):
    c.line(X0, y, X1, y, EDGE, 1.2)
    c.caps(X0, y + 40, 'Read further', MUTED, 13)
    colw = CW / 3
    for k, (lab, href) in enumerate(LINKS):
        c.link(X0 + (k % 3) * colw, y + 72 + (k // 3) * 28, lab, href, 15.5)
    yy = y + 72 + math.ceil(len(LINKS) / 3) * 28 + 16
    c.para(X0, yy, CW,
           'Carrier drawings use the certified icosahedral geometry and the recorded creation and '
           'transport run. Reference masses: CODATA 2022 and PDG 2026. They are comparison data '
           'and no output of the program.', 14.5, FAINT, 1.45)
    return yy + 30
