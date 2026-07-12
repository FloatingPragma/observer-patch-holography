#!/usr/bin/env python3
"""Generate a social-card SVG for the quark target-ancestry audit."""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, getcontext
from pathlib import Path
from xml.sax.saxutils import escape


getcontext().prec = 80

ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = ROOT.parent
PUBLIC_THEOREM_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_public_exact_yukawa_end_to_end_theorem.json"
REFERENCE_JSON = ROOT / "particles" / "data" / "particle_reference_values.json"
DEFAULT_OUTPUT = REPO_ROOT / "assets" / "quark_pdg_2025_exact_social.svg"

WIDTH = 1600
HEIGHT = 800

FONT_SANS = "'IBM Plex Sans','Avenir Next','Segoe UI',sans-serif"
FONT_MONO = "'IBM Plex Mono','SFMono-Regular','Consolas',monospace"

QUARK_ORDER = ["u", "d", "s", "c", "b", "t"]
REFERENCE_KEYS = {
    "u": "up_quark",
    "d": "down_quark",
    "s": "strange_quark",
    "c": "charm_quark",
    "b": "bottom_quark",
    "t": "top_quark",
}
ROW_COLORS = {
    "u": "#5eead4",
    "d": "#60a5fa",
    "s": "#c084fc",
    "c": "#f59e0b",
    "b": "#f97316",
    "t": "#fb7185",
}
COMPARISON_CHARTS = {
    "u": "MSbar at 2 GeV",
    "d": "MSbar at 2 GeV",
    "s": "MSbar at 2 GeV",
    "c": "MSbar at self-scale",
    "b": "MSbar at self-scale",
    "t": "separate pole extraction",
}


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"), parse_float=Decimal)


def _fmt_decimal(value: Decimal) -> str:
    return format(value, "f")


def _fmt_int(value: int | float) -> str:
    return str(value)


def _text(
    x: float,
    y: float,
    content: str,
    *,
    size: int,
    fill: str,
    weight: int = 400,
    anchor: str = "start",
    family: str = FONT_SANS,
) -> str:
    return (
        f'<text x="{x}" y="{y}" fill="{fill}" font-size="{size}" text-anchor="{anchor}" '
        f'font-weight="{weight}" font-family="{family}">{escape(content)}</text>'
    )


def _rect(x: float, y: float, w: float, h: float, *, fill: str, rx: float = 0, stroke: str | None = None, stroke_width: float = 0, opacity: float | None = None) -> str:
    parts = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}"']
    if rx:
        parts.append(f' rx="{rx}"')
    if stroke:
        parts.append(f' stroke="{stroke}" stroke-width="{stroke_width}"')
    if opacity is not None:
        parts.append(f' opacity="{opacity}"')
    parts.append("/>")
    return "".join(parts)


def _line(x1: float, y1: float, x2: float, y2: float, *, stroke: str, width: float, opacity: float = 1.0) -> str:
    return (
        f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke}" '
        f'stroke-width="{width}" opacity="{opacity}"/>'
    )


def build_svg() -> str:
    public_theorem = _load_json(PUBLIC_THEOREM_JSON)
    reference_bundle = _load_json(REFERENCE_JSON)
    exact_masses = public_theorem["public_exact_outputs"]["exact_running_values_gev"]
    references = reference_bundle["entries"]

    rows: list[dict[str, str]] = []
    max_abs_diff = Decimal("0")
    for symbol in QUARK_ORDER:
        oph_mass = Decimal(exact_masses[symbol])
        pdg_mass = Decimal(references[REFERENCE_KEYS[symbol]]["value_gev"])
        delta = oph_mass - pdg_mass
        if abs(delta) > max_abs_diff:
            max_abs_diff = abs(delta)
        rows.append(
            {
                "symbol": symbol,
                "oph_mass": _fmt_decimal(oph_mass),
                "pdg_mass": _fmt_decimal(pdg_mass),
                "delta_mass": _fmt_decimal(delta),
                "comparison_chart": COMPARISON_CHARTS[symbol],
            }
        )

    title = "Quark Target-Ancestry Audit: Not a Prediction"
    subtitle = "source-only spread fiber = (R>0)^2 | two moduli remain free | numeric public rows withheld"
    badge = "THEOREM-GRADE OBSTRUCTION"
    footer = "Stored rows copy a mixed-convention target packet. GeV matrices are mass textures, not physical Yukawas."

    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}" role="img" aria-labelledby="title desc">',
        f"<title id=\"title\">{escape(title)}</title>",
        (
            "<desc id=\"desc\">A target-ancestry audit showing that stored quark coordinates reproduce their "
            "input target packet. The source equations leave two independent spread moduli, the packet mixes "
            "comparison charts, and no physical dimensionless Yukawa matrix is certified.</desc>"
        ),
        "<defs>",
        '  <linearGradient id="bg" x1="0" y1="0" x2="1" y2="1"><stop offset="0%" stop-color="#07111c"/><stop offset="100%" stop-color="#0d1526"/></linearGradient>',
        '  <linearGradient id="panel" x1="0" y1="0" x2="1" y2="1"><stop offset="0%" stop-color="#0d1627"/><stop offset="100%" stop-color="#0b1220"/></linearGradient>',
        '  <linearGradient id="badge" x1="0" y1="0" x2="1" y2="0"><stop offset="0%" stop-color="#14b8a6"/><stop offset="100%" stop-color="#38bdf8"/></linearGradient>',
        "</defs>",
        _rect(0, 0, WIDTH, HEIGHT, fill="url(#bg)"),
        _text(48, 82, title, size=36, fill="#f8fbff", weight=800),
        _text(48, 118, subtitle, size=20, fill="#b9cae3", family=FONT_MONO),
        _rect(1170, 40, 374, 48, fill="url(#badge)", rx=24, opacity=0.98),
        _text(1357, 71, badge, size=18, fill="#06111d", weight=800, anchor="middle", family=FONT_MONO),
        _rect(40, 160, 930, 560, fill="url(#panel)", rx=28, stroke="#24364c", stroke_width=1.2),
        _rect(1000, 160, 560, 560, fill="url(#panel)", rx=28, stroke="#24364c", stroke_width=1.2),
        _rect(56, 182, 898, 52, fill="#0f1b2f", rx=18, stroke="#274565", stroke_width=1.0),
        _rect(1016, 182, 528, 52, fill="#0f1b2f", rx=18, stroke="#274565", stroke_width=1.0),
        _text(84, 216, "q", size=18, fill="#f8fbff", weight=700, family=FONT_MONO),
        _text(150, 216, "stored target audit [GeV]", size=18, fill="#f8fbff", weight=700, family=FONT_MONO),
        _text(420, 216, "input target [GeV]", size=18, fill="#f8fbff", weight=700, family=FONT_MONO),
        _text(690, 216, "copy residual [GeV]", size=18, fill="#f8fbff", weight=700, family=FONT_MONO),
        _text(1048, 216, "comparison chart", size=18, fill="#f8fbff", weight=700, family=FONT_MONO),
    ]

    start_y = 252
    row_h = 74
    for i, row in enumerate(rows):
        y = start_y + i * row_h
        fill = "#0d1728" if i % 2 == 0 else "#0a1322"
        svg.append(_rect(56, y, 898, 56, fill=fill, rx=16, stroke="#19283a", stroke_width=0.8))
        svg.append(_rect(1016, y, 528, 56, fill=fill, rx=16, stroke="#19283a", stroke_width=0.8))
        svg.append(_rect(72, y + 10, 46, 36, fill=ROW_COLORS[row["symbol"]], rx=18, opacity=0.96))
        svg.append(_text(95, y + 35, row["symbol"], size=24, fill="#08111b", weight=800, anchor="middle", family=FONT_MONO))
        svg.append(_text(150, y + 37, row["oph_mass"], size=18, fill="#e9f2ff", family=FONT_MONO))
        svg.append(_text(420, y + 37, row["pdg_mass"], size=18, fill="#dbeafe", family=FONT_MONO))
        svg.append(_text(690, y + 37, row["delta_mass"], size=18, fill="#d1fae5", family=FONT_MONO))
        svg.append(_rect(1032, y + 10, 46, 36, fill=ROW_COLORS[row["symbol"]], rx=18, opacity=0.96))
        svg.append(_text(1055, y + 35, row["symbol"], size=24, fill="#08111b", weight=800, anchor="middle", family=FONT_MONO))
        svg.append(_text(1106, y + 37, row["comparison_chart"], size=17, fill="#fde68a", family=FONT_MONO))

    table_bottom = start_y + len(rows) * row_h
    svg.extend(
        [
            _line(48, 744, 1552, 744, stroke="#23364a", width=1.0, opacity=0.9),
            _text(48, 776, footer, size=16, fill="#b8c8df", family=FONT_MONO),
        ]
    )

    svg.append("</svg>")
    return "\n".join(svg)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate the quark target-ancestry audit social-card SVG.")
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    args = parser.parse_args()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(build_svg() + "\n", encoding="utf-8")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
