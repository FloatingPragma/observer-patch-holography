"""Deterministic SVG primitives with metric text wrapping for the explainer."""
from html import escape
import math
import re

SANS = "'Helvetica Neue',Helvetica,Arial,'Liberation Sans',sans-serif"
MONO = "ui-monospace,'SF Mono',Menlo,Consolas,'Liberation Mono',monospace"
SERIF = "Georgia,'Times New Roman',serif"

# Palette. Teal marks the program, amber the defect, violet the observers'
# description. Red is an ordinary mismatch that repair removes.
BG = '#090d16'
CARD = '#111a2b'
CARD2 = '#0d1422'
EDGE = '#22304d'
TEXT = '#e9eef8'
SOFT = '#c2cde2'
MUTED = '#93a3c0'
FAINT = '#62718f'
TEAL = '#3fe0d0'
AMBER = '#ffb44a'
ORANGE = '#ff8a3c'
VIOLET = '#a995ff'
BLUE = '#6fb1ff'
GREEN = '#5fe3a1'
RED = '#ff6480'
GOLD = '#ffd77a'

STYLE = {'a': AMBER, 't': TEAL, 'v': VIOLET, 'b': TEXT, 'g': GREEN,
         'r': RED, 'o': GOLD, 'm': MUTED, 's': SOFT}

# Helvetica advance widths per 1000 em, ASCII 32..126.
_REG = [278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278, 278,
        556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 278, 278, 584, 584, 584, 556, 1015,
        667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833, 722, 778, 667, 778,
        722, 667, 611, 722, 667, 944, 667, 667, 611, 278, 278, 278, 469, 556, 333,
        556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833, 556, 556, 556, 556,
        333, 500, 278, 556, 500, 722, 500, 500, 500, 334, 260, 334, 584]
_BOLD = [278, 333, 474, 556, 556, 889, 722, 238, 333, 333, 389, 584, 278, 333, 278, 278,
         556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 333, 333, 584, 584, 584, 611, 975,
         722, 722, 722, 722, 667, 611, 778, 722, 278, 556, 722, 611, 833, 722, 778, 667, 778,
         722, 667, 611, 722, 667, 944, 667, 667, 611, 333, 278, 333, 584, 556, 333,
         556, 611, 556, 611, 556, 333, 611, 611, 278, 278, 556, 278, 889, 611, 611, 611, 611,
         389, 556, 333, 611, 556, 778, 556, 556, 500, 389, 280, 389, 584]
_EXTRA = {'−': 584, '–': 556, '×': 584, '·': 278, '→': 1000, '←': 1000, '↔': 1000,
          '⇄': 1000, '²': 333, '½': 834, '±': 584, '≈': 584, '≠': 584, '’': 222,
          '‘': 222, '“': 333, '”': 333, '…': 1000, 'γ': 520, 'μ': 576, 'τ': 440,
          'ν': 500, 'ψ': 720, '∂': 494, 'π': 560, '°': 400, '✓': 780, '✗': 780,
          '↑': 500, '↓': 500, 'ö': 556, '⊕': 760, '⊗': 760, '≥': 584, '≤': 584}
SAFETY = 1.035
assert len(_REG) == 95 and len(_BOLD) == 95


def text_width(s, size, bold=False, ls=0.0, mono=False):
    """Conservative rendered width of one line in user units."""
    if mono:
        return len(s) * (size * 0.602 + ls)
    table = _BOLD if bold else _REG
    total = 0
    for ch in s:
        o = ord(ch)
        total += table[o - 32] if 32 <= o <= 126 else _EXTRA.get(ch, 620)
    return total * size / 1000 * SAFETY + ls * len(s)


def n(v):
    """Short stable number formatting."""
    s = f'{v:.2f}'.rstrip('0').rstrip('.')
    return '0' if s == '-0' else s


def _runs(markup):
    """Split '{a|styled} plain' markup into (style, text) runs."""
    out = []
    for part in re.split(r'(\{[a-z]\|[^}]*\})', markup):
        if not part:
            continue
        if part.startswith('{') and part[2:3] == '|':
            out.append((part[1], part[3:-1]))
        else:
            out.append(('', part))
    return out


class Canvas:
    def __init__(self):
        self.items = []
        self.defs = []
        self._ids = set()
        self.overflow = []
        self.text_count = 0

    # ---- low level -------------------------------------------------
    def add(self, s):
        self.items.append(s)

    def define(self, key, body):
        if key not in self._ids:
            self._ids.add(key)
            self.defs.append(body)
        return key

    @staticmethod
    def _attrs(**kw):
        out = []
        for k, v in kw.items():
            if v is None or v is False:
                continue
            out.append(f'{k.replace("_", "-")}="{v}"')
        return (' ' + ' '.join(out)) if out else ''

    def open(self, **kw):
        self.add(f'<g{self._attrs(**kw)}>')

    def close(self):
        self.add('</g>')

    def rect(self, x, y, w, h, fill='none', stroke=None, sw=1.5, r=0, op=None, dash=None, cls=None):
        self.add(f'<rect x="{n(x)}" y="{n(y)}" width="{n(w)}" height="{n(h)}"'
                 + (f' rx="{n(r)}"' if r else '')
                 + self._attrs(fill=fill, stroke=stroke, stroke_width=n(sw) if stroke else None,
                               opacity=op, stroke_dasharray=dash, **{'class': cls}) + '/>')

    def circle(self, x, y, r, fill='none', stroke=None, sw=1.5, op=None, dash=None, cls=None):
        self.add(f'<circle cx="{n(x)}" cy="{n(y)}" r="{n(r)}"'
                 + self._attrs(fill=fill, stroke=stroke, stroke_width=n(sw) if stroke else None,
                               opacity=op, stroke_dasharray=dash, **{'class': cls}) + '/>')

    def ellipse(self, x, y, rx, ry, fill='none', stroke=None, sw=1.5, op=None, dash=None):
        self.add(f'<ellipse cx="{n(x)}" cy="{n(y)}" rx="{n(rx)}" ry="{n(ry)}"'
                 + self._attrs(fill=fill, stroke=stroke, stroke_width=n(sw) if stroke else None,
                               opacity=op, stroke_dasharray=dash) + '/>')

    def path(self, d, stroke=None, sw=2, fill='none', op=None, dash=None, cap='round', cls=None, extra=''):
        self.add(f'<path d="{d}"'
                 + self._attrs(fill=fill, stroke=stroke, stroke_width=n(sw) if stroke else None,
                               opacity=op, stroke_dasharray=dash,
                               stroke_linecap=cap if stroke else None,
                               stroke_linejoin='round' if stroke else None, **{'class': cls})
                 + (' ' + extra if extra else '') + '/>')

    def line(self, x1, y1, x2, y2, stroke, sw=2, op=None, dash=None, cap='round', cls=None):
        self.add(f'<line x1="{n(x1)}" y1="{n(y1)}" x2="{n(x2)}" y2="{n(y2)}"'
                 + self._attrs(stroke=stroke, stroke_width=n(sw), opacity=op,
                               stroke_dasharray=dash, stroke_linecap=cap, **{'class': cls}) + '/>')

    def poly(self, pts, fill='none', stroke=None, sw=1.5, op=None, cls=None, join='round'):
        d = ' '.join(f'{n(a)},{n(b)}' for a, b in pts)
        self.add(f'<polygon points="{d}"'
                 + self._attrs(fill=fill, stroke=stroke, stroke_width=n(sw) if stroke else None,
                               opacity=op, stroke_linejoin=join if stroke else None,
                               **{'class': cls}) + '/>')

    def curve(self, pts, stroke, sw=2, op=None, dash=None, fill='none', cls=None):
        d = 'M' + ' L'.join(f'{n(a)},{n(b)}' for a, b in pts)
        self.path(d, stroke, sw, fill, op, dash, cls=cls)

    # ---- composite marks -------------------------------------------
    def glow(self, x, y, r, color, op=0.55, cls=None):
        key = 'rg' + color.lstrip('#')
        self.define(key, f'<radialGradient id="{key}"><stop offset="0" stop-color="{color}" stop-opacity="1"/>'
                         f'<stop offset=".45" stop-color="{color}" stop-opacity=".35"/>'
                         f'<stop offset="1" stop-color="{color}" stop-opacity="0"/></radialGradient>')
        self.circle(x, y, r, f'url(#{key})', op=op, cls=cls)

    def head(self, x, y, ang, color, size=9, op=None):
        """Arrowhead with its tip at (x, y), pointing along ang (radians)."""
        a, b = ang + math.radians(152), ang - math.radians(152)
        self.poly([(x, y), (x + size * math.cos(a), y + size * math.sin(a)),
                   (x + size * math.cos(b), y + size * math.sin(b))], color, op=op)

    def arrow(self, x1, y1, x2, y2, color, sw=2, size=9, op=None, dash=None):
        ang = math.atan2(y2 - y1, x2 - x1)
        back = size * 0.75
        self.line(x1, y1, x2 - back * math.cos(ang), y2 - back * math.sin(ang), color, sw, op, dash)
        self.head(x2, y2, ang, color, size, op)

    def qarrow(self, x1, y1, cx, cy, x2, y2, color, sw=2, size=9, op=None, dash=None):
        """Quadratic curve with an arrowhead at the end."""
        ang = math.atan2(y2 - cy, x2 - cx)
        ex, ey = x2 - size * 0.7 * math.cos(ang), y2 - size * 0.7 * math.sin(ang)
        self.path(f'M{n(x1)},{n(y1)} Q{n(cx)},{n(cy)} {n(ex)},{n(ey)}', color, sw, op=op, dash=dash)
        self.head(x2, y2, ang, color, size, op)

    def tick(self, x, y, color, size=7, sw=2.2):
        """A check mark drawn as a path, centred on (x, y)."""
        self.path(f'M{n(x - size * 0.62)},{n(y + size * 0.05)} L{n(x - size * 0.15)},{n(y + size * 0.52)} '
                  f'L{n(x + size * 0.68)},{n(y - size * 0.5)}', color, sw)

    def cross(self, x, y, color, size=6, sw=2.2):
        self.line(x - size * 0.6, y - size * 0.6, x + size * 0.6, y + size * 0.6, color, sw)
        self.line(x - size * 0.6, y + size * 0.6, x + size * 0.6, y - size * 0.6, color, sw)

    def dart(self, x, y, ang, length, color, sw=1.6, op=None):
        """Small direction mark used in arrow fields."""
        dx, dy = math.cos(ang) * length / 2, math.sin(ang) * length / 2
        self.line(x - dx, y - dy, x + dx * 0.5, y + dy * 0.5, color, sw, op)
        self.head(x + dx, y + dy, ang, color, max(4.2, length * 0.36), op)

    # ---- text ------------------------------------------------------
    def text(self, x, y, s, size=20, fill=TEXT, weight=400, anchor='start', family=SANS,
             ls=0.0, italic=False, op=None, maxw=None, cls=None):
        s = str(s)
        self.text_count += 1
        if maxw is not None:
            w = text_width(s, size, weight >= 600, ls, family is MONO)
            if w > maxw + 0.5:
                self.overflow.append((round(w - maxw, 1), s))
        self.add(f'<text x="{n(x)}" y="{n(y)}" font-family="{family}" font-size="{n(size)}"'
                 + self._attrs(fill=fill, font_weight=weight if weight != 400 else None,
                               text_anchor=anchor if anchor != 'start' else None,
                               letter_spacing=n(ls) if ls else None,
                               font_style='italic' if italic else None, opacity=op,
                               **{'class': cls})
                 + f'>{escape(s)}</text>')

    def caps(self, x, y, s, fill=TEAL, size=14, anchor='start', ls=2.2, weight=700, maxw=None, op=None):
        self.text(x, y, s.upper(), size, fill, weight, anchor, ls=ls, maxw=maxw, op=op)

    def para(self, x, y, w, markup, size=20, fill=SOFT, lh=1.5, weight=400, anchor='start',
             bold_weight=700, max_lines=None):
        """Wrapped rich text. Returns the baseline y after the last line."""
        words = []
        for style, chunk in _runs(markup):
            for i, piece in enumerate(re.split(r'(\s+)', chunk)):
                if not piece:
                    continue
                if piece.isspace():
                    words.append((None, ' '))
                else:
                    words.append((style, piece))
        lines, cur, cur_w = [], [], 0.0
        space_w = text_width(' ', size, weight >= 600)
        pending_space = False
        for style, piece in words:
            if style is None:
                pending_space = bool(cur)
                continue
            bold = bool(style) or weight >= 600
            pw = text_width(piece, size, bold)
            add = pw + (space_w if pending_space else 0)
            if cur and pending_space and cur_w + add > w:
                lines.append(cur)
                cur, cur_w = [(style, piece)], pw
            else:
                if pending_space:
                    cur.append((cur[-1][0] if cur[-1][0] == style else '', ' '))
                cur.append((style, piece))
                cur_w += add
            pending_space = False
        if cur:
            lines.append(cur)
        if max_lines is not None and len(lines) > max_lines:
            self.overflow.append((len(lines) - max_lines, 'lines: ' + markup[:60]))
        step = size * lh
        for line in lines:
            merged = []
            for style, piece in line:
                if merged and merged[-1][0] == style:
                    merged[-1][1] += piece
                else:
                    merged.append([style, piece])
            body = ''
            for style, piece in merged:
                if style:
                    body += (f'<tspan fill="{STYLE[style]}" font-weight="{bold_weight}">'
                             f'{escape(piece)}</tspan>')
                else:
                    body += escape(piece)
            self.text_count += 1
            self.add(f'<text x="{n(x)}" y="{n(y)}" font-family="{SANS}" font-size="{n(size)}"'
                     + self._attrs(fill=fill, font_weight=weight if weight != 400 else None,
                                   text_anchor=anchor if anchor != 'start' else None)
                     + f' xml:space="preserve">{body}</text>')
            y += step
        return y

    def link(self, x, y, label, href, size=16, fill=BLUE, anchor='start'):
        self.add(f'<a href="{escape(href, quote=True)}" target="_blank">')
        self.text(x, y, label + ' ↗', size, fill, anchor=anchor)
        self.add('</a>')

    def chip(self, x, y, label, color, size=13, pad=11, h=26, filled=False, dash=None, anchor='start'):
        """Small rounded label. Returns its width."""
        w = text_width(label.upper(), size, True, 1.6) + 2 * pad
        if anchor == 'middle':
            x -= w / 2
        elif anchor == 'end':
            x -= w
        self.rect(x, y, w, h, color if filled else 'none', None if filled else color, 1.3, h / 2,
                  dash=dash, op=None)
        self.text(x + w / 2 + 0.8, y + h / 2 + size * 0.36, label.upper(), size,
                  BG if filled else color, 700, 'middle', ls=1.6)
        return w
