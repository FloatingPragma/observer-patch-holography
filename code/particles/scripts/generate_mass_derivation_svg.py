#!/usr/bin/env python3
"""Generate a detailed dark-theme SVG for the implemented `/particles` pipeline."""

from __future__ import annotations

import argparse
import json
import pathlib
import textwrap
from collections import Counter
from datetime import datetime, timezone
from typing import Any, Dict, List, Sequence, Tuple
from xml.sax.saxutils import escape


ROOT = pathlib.Path(__file__).resolve().parents[2]
RESULTS_JSON = ROOT / "particles" / "results_status.json"
EXACT_NONHADRON_JSON = ROOT / "particles" / "exact_nonhadron_masses.json"
DEFAULT_OUTPUT = ROOT / "particles" / "particle_mass_derivation_graph.svg"
PUBLIC_PIXEL_DISPLAY = "1.630968"
PUBLIC_CAPACITY_DISPLAY = "3.5323546e122"

WIDTH = 2560
MARGIN_X = 72
PANEL_GAP_X = 36
PANEL_GAP_Y = 48
PANEL_W = (WIDTH - 2 * MARGIN_X - 2 * PANEL_GAP_X) / 3.0
HADRON_W = PANEL_W * 2.0 + PANEL_GAP_X

FONT_FAMILY = "'IBM Plex Sans','Avenir Next','Segoe UI',sans-serif"
MONO_FAMILY = "'IBM Plex Mono','SFMono-Regular','Consolas',monospace"

COLORS = {
    "bg": "#050b15",
    "bg_alt": "#0a1424",
    "panel": "#0c1729",
    "panel_alt": "#101e35",
    "panel_shell": "#0a1322",
    "panel_border": "#24415f",
    "panel_border_soft": "#20344d",
    "ink": "#eef4ff",
    "muted": "#a5b6d3",
    "subtle": "#7890b3",
    "line": "#7aa8ff",
    "line_soft": "#496b96",
    "axiom_fill": "#11203a",
    "axiom_stroke": "#60a5fa",
    "input_fill": "#162540",
    "input_stroke": "#7dd3fc",
    "logic_fill": "#10293a",
    "logic_stroke": "#2dd4bf",
    "status_fill": "#1a2340",
    "status_stroke": "#94a3b8",
    "prediction_fill": "#13251e",
    "prediction_stroke": "#4ade80",
    "task_fill": "#2a1816",
    "task_stroke": "#fb923c",
    "task_pending": "#f97316",
    "task_in_progress": "#fb923c",
    "task_complete": "#22c55e",
    "task_deferred": "#64748b",
    "task_out_of_scope": "#ef4444",
    "task_blocked": "#ef4444",
    "footer_fill": "#0f1d33",
    "footer_stroke": "#2e4c75",
    "green_note_fill": "#11281d",
    "green_note_stroke": "#22c55e",
    "green_note_text": "#d6ffe5",
    "green_note_body": "#b6e6c6",
}

STATUS_BAR = {
    "structural": "#38bdf8",
    "conditional_carrier_mode": "#f59e0b",
    "calibration": "#06b6d4",
    "candidate_trunk_compare_only": "#06b6d4",
    "secondary_quantitative": "#8b5cf6",
    "declared_surface_theorem": "#8b5cf6",
    "target_anchored_witness_no_go_boundary": "#f59e0b",
    "selected_class_theorem": "#14b8a6",
    "source_spread_nonidentifiability_obstruction": "#f59e0b",
    "continuation": "#f59e0b",
    "simulation_dependent": "#ef4444",
}

STATUS_TEXT = {
    "structural": "structural",
    "conditional_carrier_mode": "classical mode / quantum gate open",
    "calibration": "electroweak closure",
    "candidate_trunk_compare_only": "EW frontier",
    "secondary_quantitative": "secondary",
    "declared_surface_theorem": "declared-surface theorem",
    "target_anchored_witness_no_go_boundary": "witness boundary",
    "selected_class_theorem": "legacy selected-class audit",
    "source_spread_nonidentifiability_obstruction": "reciprocal-ray rejected",
    "continuation": "continuation",
    "simulation_dependent": "out of scope",
}

PARTICLE_INFO: Dict[str, Dict[str, str]] = {
    "photon": {"symbol": "gamma mode", "plain": "Two transverse classical Maxwell modes on the declared unbroken/deconfined action branch; quantum-particle gate open."},
    "gluon": {"symbol": "g mode", "plain": "Perturbative transverse Yang-Mills modes before confinement; no asymptotic colored particle is claimed."},
    "graviton": {"symbol": "TT mode", "plain": "Two classical TT modes of pure Einstein linearization; no graviton Hilbert space or particle pole is constructed."},
    "w_boson": {"symbol": "W", "plain": "Charged weak-force boson used in beta-decay-type processes."},
    "z_boson": {"symbol": "Z", "plain": "Neutral weak-force boson from the same electroweak sector as the W."},
    "higgs": {"symbol": "H", "plain": "Higgs boson tied to the Standard Model mass-giving field."},
    "top_quark": {"symbol": "t", "plain": "Heaviest known quark. Its conditional Higgs/top companion coordinate is separate from the rejected reciprocal-ray quark spectrum; no nonzero mass is printed here."},
    "electron": {"symbol": "e", "plain": "Light charged matter particle found in atoms."},
    "muon": {"symbol": "mu", "plain": "Heavier unstable cousin of the electron."},
    "tau": {"symbol": "tau", "plain": "Heaviest charged lepton; a short-lived electron cousin."},
    "up_quark": {"symbol": "u", "plain": "Elementary quark that helps build protons and neutrons."},
    "down_quark": {"symbol": "d", "plain": "Elementary quark that helps build protons and neutrons."},
    "strange_quark": {"symbol": "s", "plain": "Heavier quark seen inside unstable hadrons."},
    "charm_quark": {"symbol": "c", "plain": "Heavier quark seen inside unstable hadrons."},
    "bottom_quark": {"symbol": "b", "plain": "Very heavy quark used in flavor-physics tests."},
    "electron_neutrino": {"symbol": "nu_e", "plain": "Extremely light neutral lepton of the electron family."},
    "muon_neutrino": {"symbol": "nu_mu", "plain": "Extremely light neutral lepton of the muon family."},
    "tau_neutrino": {"symbol": "nu_tau", "plain": "Extremely light neutral lepton of the tau family."},
    "proton": {"symbol": "p", "plain": "Stable positively charged hadron found in atomic nuclei."},
    "neutron": {"symbol": "n", "plain": "Neutral hadron found in atomic nuclei."},
    "neutral_pion": {"symbol": "pi0 proxy", "plain": "Light meson placeholder row withheld from the source-only table because the hadron backend is absent."},
    "rho_770_0": {"symbol": "rho(770)0 proxy", "plain": "Vector-meson placeholder row withheld from the source-only table because the hadron backend is absent."},
}

GROUP_ROW_TEXT = {
    "Bosons": "This row sits on the force-carrier / Higgs side of the spectrum.",
    "Leptons": "This row sits on the lepton and neutrino side of the spectrum.",
    "Quarks": "This row sits on the quark flavor side of the spectrum.",
    "Hadrons": "This row sits on the bound-state QCD side of the spectrum.",
}

STATUS_EXPLAINER = {
    "structural": "exact structural theorem surface",
    "conditional_carrier_mode": "declared classical-action mode with the quantum-particle gate open",
    "calibration": "implemented P-driven electroweak quantitative-closure surface",
    "candidate_trunk_compare_only": "electroweak massive-boson frontier without public W/Z rows",
    "secondary_quantitative": "quantitative secondary branch with a separate proof package",
    "declared_surface_theorem": "theorem surface on declared running and matching conventions",
    "target_anchored_witness_no_go_boundary": "target-anchored witness with a corpus-limited no-go boundary",
    "selected_class_theorem": "legacy selected-class target audit with numeric witness withheld",
    "source_spread_nonidentifiability_obstruction": "common-scale six-Yukawa audit rejects the reciprocal-ray subfamily; numeric rows withheld",
    "continuation": "declared continuation or witness surface outside theorem-grade public output",
    "simulation_dependent": "source-backend-absent surface with empirical closure policy",
}

STATUS_NEXT_STEP = {
    "structural": "This row belongs to a structural theorem surface; particle interpretation requires the explicit carrier gates.",
    "conditional_carrier_mode": "This row records a zero hard quadratic parameter only; quantization, phase, spectrum, and positive-residue pole receipts remain required for particle promotion.",
    "calibration": "This row belongs to the implemented P-driven electroweak closure surface.",
    "candidate_trunk_compare_only": "This row belongs to the electroweak massive-boson frontier.",
    "secondary_quantitative": "This row belongs to a quantitative secondary branch with its own proof surface.",
    "declared_surface_theorem": "This row belongs to a theorem surface on declared running and matching conventions.",
    "target_anchored_witness_no_go_boundary": "This row belongs to a target-anchored witness with a corpus-limited no-go boundary.",
    "selected_class_theorem": "This legacy status denotes a selected-class target audit; numeric witness values are withheld from public prediction outputs.",
    "source_spread_nonidentifiability_obstruction": "The generic physical interface has six dimensionless Yukawa coordinates at one scale and scheme. The reciprocal-ray restriction fails the common-scale data, and a source-derived flavor-orbit selector is absent.",
    "continuation": "This row belongs to a declared continuation or witness surface outside theorem-grade public output.",
    "simulation_dependent": "This source-only row requires a backend bundle and publication-grade systematics. Empirical closure values use a separate e+e- payload class.",
}

PARTICLE_TITLE = {
    "photon": "Photon",
    "gluon": "Gluons",
    "graviton": "Graviton",
    "w_boson": "W Boson",
    "z_boson": "Z Boson",
    "higgs": "Higgs Boson",
    "top_quark": "Top Quark",
    "electron": "Electron",
    "muon": "Muon",
    "tau": "Tau Lepton",
    "up_quark": "Up Quark",
    "down_quark": "Down Quark",
    "strange_quark": "Strange Quark",
    "charm_quark": "Charm Quark",
    "bottom_quark": "Bottom Quark",
    "electron_neutrino": "Electron Neutrino",
    "muon_neutrino": "Muon Neutrino",
    "tau_neutrino": "Tau Neutrino",
    "proton": "Proton",
    "neutron": "Neutron",
    "neutral_pion": "Neutral Pion Proxy",
    "rho_770_0": "Rho(770)0 Proxy",
}

LANES: List[Dict[str, Any]] = [
    {
        "key": "structural",
        "title": "Conditional Classical Carrier Modes",
        "summary": "Maxwell, perturbative Yang-Mills, and pure-Einstein quadratic branches have massless classical modes under their stated action and phase assumptions.",
        "takeaway": "The displayed zero is a hard quadratic action parameter, not a 0 GeV quantum-particle prediction.",
        "logic": (
            "After separately selecting a positive Maxwell/Yang-Mills kinetic action or the pure Einstein-Hilbert "
            "kinetic branch, gauge fixing exposes two transverse modes per gauge generator or two TT tensor modes. "
            "The abstract compact group alone supplies neither the connection nor its kinetic term. Higgs/Stueckelberg "
            "phases, confinement, media, and extended gravitational field content remain outside the selected branch."
        ),
        "frontier_text": "Quantum-particle promotion requires a constructed quantization, positive physical Hilbert/Hamiltonian spectrum, positive-residue massless pole, and a stable deconfined/asymptotic state. The confining QCD phase does not pass the gluon-particle gate.",
        "prediction_surface": "No public particle mass is emitted. The separate receipt records only conditional classical/perturbative mode gates and their zero hard quadratic parameters.",
        "particles": ["photon", "gluon", "graviton"],
    },
    {
        "key": "d10",
        "status": "candidate_trunk_compare_only",
        "title": "Electroweak Massive-Boson Audit",
        "summary": "The W/Z lane separates a selected-carrier chart, a conditional five-premise value law, and a measured-reference inverse adapter.",
        "takeaway": "The older rounded boundary aliases belong to the conditional value-law candidate, not the inverse adapter. No nonzero source-only physical mass or complex-pole row is emitted.",
        "logic": (
            "From P the calculation builds the electroweak running family and the source basis that feeds the mass chart. "
            "The selected-carrier chart and the conditional five-premise value law are distinct objects. The finite-carrier "
            "certificate behind the latter is absent. The reference-fitted inverse comparison adapter consumes the reference "
            "boundary. The electromagnetic side is a separate endpoint-closure problem."
        ),
        "frontier_text": "A source-only W/Z result requires the finite quotient-transport certificate, one independently fixed physical scale and scheme, a single source branch, prospective provenance, and complex-pole and uncertainty receipts.",
        "prediction_surface": "No nonzero source-only physical W/Z mass is emitted. The chart, conditional value law, and inverse adapter remain on their separate audit surfaces.",
        "particles": [],
    },
    {
        "key": "d11",
        "title": "Higgs/Top Conditional Surface",
        "summary": "The Higgs/top lane emits a conditional exact candidate on the declared running, matching, and threshold surface. The one-scalar seed is a lower-rank fixed-ray branch.",
        "takeaway": "The Higgs coordinate and companion top coordinate are retained as conditional audit objects. This diagram prints no nonzero particle mass.",
        "logic": (
            "Take the electroweak substrate, impose the critical-surface condition, then use the synchronized "
            "transport core on the declared running, matching, and threshold surface. The fixed-ray certificate "
            "is a lower-rank one-scalar branch. The split theorem emits a shared Higgs/top scalar together "
            "with source-only residual selectors for the top and Higgs channels. The pair remains conditional on the "
            "declared electroweak running and matching surface and is not a source-only complex-pole result. The quark "
            "spectrum is assessed separately on a common-scale dimensionless Yukawa interface."
        ),
        "frontier_text": "Claim boundary: the Higgs/top pair is a conditional declared-surface candidate. The one-scalar seed is a lower-rank fixed-ray branch, while the inverse slice and direct-top conversion are comparison objects.",
        "prediction_surface": "No nonzero source-only physical Higgs or top mass is emitted. The conditional split remains on its declared audit surface.",
        "particles": ["higgs"],
    },
    {
        "key": "leptons",
        "title": "Charged Leptons",
        "summary": "The charged-lepton lane has exact icosahedral face geometry, a conditional affine response law, and an engineered finite digital model that proves the stipulated schema is nonempty.",
        "takeaway": "The digital witness closes fixed-cutoff satisfiability and the event-to-central-record bridge. It does not select the physical charged response or emit a nonzero source-only mass.",
        "logic": (
            "The finite witness instantiates a bounded observer-like self-reading patch: local register state, declared path "
            "interfaces, accepted/rejected central readback, recovery maps, and a public evidence archive. Eight register "
            "graphs account for 6467 matrix units, eight declared paths, 60 A5 rotation charts, and the central record dilation. "
            "The construction hard-codes its dimensions, automaton, signs, clock, and response; its dependency record omits "
            "target-informed ancestry. Per-event response maps do not form a checked global mutually exclusive channel, and "
            "the exported refinement is inert ancillary stabilization rather than physical cofinal refinement. A conditional "
            "nature/pole transport lemma assumes the physical Yukawa response and singularity readout it transports; its "
            "zero-self-energy kernel is a free witness rather than the interacting charged kernel."
        ),
        "frontier_text": "Physical promotion requires source-derived dynamical selection, documented target-free ancestry, identification with the charged-generation space, a global response channel, an interacting infrared-complete kernel, determinant and pole attachment, and cofinal physical refinement.",
        "prediction_surface": "No nonzero source-only charged-lepton mass is emitted. The finite digital schema witness and central-record construction remain separate from target-informed mass diagnostics.",
        "particles": ["electron", "muon", "tau"],
    },
    {
        "key": "quarks",
        "title": "Quarks",
        "summary": "The generic quark interface is six dimensionless Yukawa coordinates at one declared scale and scheme. The reciprocal-ray subfamily fails the common-scale data.",
        "takeaway": "The reciprocal product differs from one at every tested scale. Even after four endpoints are granted, the held-out charm/strange discrepancy remains about one fifth.",
        "logic": (
            "The six-scalar interface is (mu_u, sigma_u, rho_u, mu_d, sigma_d, rho_d), with an exact inverse to each "
            "ordered positive Yukawa pair. The reciprocal condition rho_u rho_d = 1 fails from the electroweak scale "
            "through the highest tested scale. Sub-percent residuals fit a mixed-convention target-informed mass "
            "chart, and a lower-order ablation fits that same chart better. Flavor-singlet source scalars cannot emit a "
            "nonzero fixed bifundamental Yukawa representative. The completed 64k simulator runs measure no quark/Higgs "
            "operator and are calibration-null. Independent and archived direct-carrier assays remain rank one in the Haar "
            "basis and miss both required transport shapes; complex orientation and refinement exports are absent."
        ),
        "frontier_text": "A physical result requires a source-derived flavor-orbit selector, a gauge-invariant chiral Higgs attachment, one common-scale six-Yukawa output, full coupled running and matching, and blind covariant validation. All eight prediction-contract gates fail.",
        "prediction_surface": "No nonzero source-only physical quark mass is emitted. The reciprocal-ray candidate is rejected; formula fits and simulator outputs remain comparison diagnostics.",
        "particles": ["up_quark", "down_quark", "strange_quark", "charm_quark", "bottom_quark", "top_quark"],
    },
    {
        "key": "neutrinos",
        "title": "Neutrinos",
        "summary": "The neutrino lane proves an exact no-go for the isotropic ansatz and retains the weighted-cycle construction only as a rejected, target-informed comparison candidate.",
        "takeaway": "No physical PMNS matrix, ordering, Majorana phases, or nonzero absolute neutrino masses are derived from source-closed OPH inputs.",
        "logic": (
            "The lane derives the light-neutrino scale from the electroweak core, builds the family-response tensor, the "
            "Majorana holonomy lift, the pullback metric, the forward Majorana matrix, and the ascending-splitting "
            "bundle. The old selector-law branch is S3-isotropic and is ruled out for the physical "
            "atmospheric scale by the exact general cap `max |Delta m^2| <= 8 a rho + 4 rho^2`. The repaired branch uses the "
            "same-label scalar certificate together with the flavor cocycle invariants `gamma`, `eps`, and "
            "`gamma_half` to define the positive load segment between `chi = 1 + eps` and `1 + gamma_half`. On that "
            "one-dimensional affine segment, the balanced selector and the least-distortion selector coincide exactly at "
            "the midpoint, so `D_nu = (chi + 1 + gamma_half) / 2` and the repaired edge law becomes "
            "`w_e = q_e^(1 + gamma + eps / D_nu)`. That construction is not source-closed: its family-transport kernel "
            "and line lift remain candidate inputs, while the stored charged-lepton left basis is open and nearly "
            "degenerate. The historical shared-basis cancellation therefore does not derive a physical PMNS matrix. "
            "After correcting the Majorana Takagi convention, the conditional intrinsic point also lies far outside the "
            "oscillation data; the earlier weighted point fails the official NuFIT 6.1 correlated "
            "`(sin^2 theta_23, delta_CP)` profile beyond the two-parameter 3-sigma threshold. The two-parameter exact adapter and `(B_nu, C_nu)` coordinates are "
            "target-fit diagnostics on this rejected candidate and cannot feed back into theorem state. A future physical "
            "result requires a source-closed neutrino operator, a nondegenerate charged-lepton basis, and an independently "
            "derived mass-eigenstate labeling rule fixed before comparison with data."
        ),
        "frontier_text": "Physical gates: derive a source-closed neutrino operator, a stable charged-lepton left basis, and a neutrino mass-eigenstate label/order rule without oscillation-target feedback. The weighted-cycle and bridge objects remain rejected diagnostics.",
        "prediction_surface": "No nonzero source-only neutrino mass or physical mixing row is emitted. The exact isotropic no-go is retained; weighted-cycle mixing, bridge coordinates, and absolute masses are rejected or comparison diagnostics.",
        "particles": ["electron_neutrino", "muon_neutrino", "tau_neutrino"],
    },
    {
        "key": "hadrons",
        "title": "Hadrons",
        "summary": "The source-only hadron lane lacks a production OPH backend. A separate measured e+e- hadron payload supplies an empirical fine-structure endpoint.",
        "takeaway": "Empirical closure gives alpha inverse 136.2636589431 on [136.1583832834, 136.3690975240], with a same-scheme anchor gap [0.6485541111, 0.8547920666].",
        "logic": (
            "The hadron path steps down from the electroweak and local quark masses into the strong scale, seeds the "
            "unquenched ensemble family, realizes deterministic cfg/source payload identifiers, attaches a fixed "
            "RHMC/HMC schedule shell and conditional execution receipt, then builds the stable-channel "
            "sequence population/evaluation shells and aggregates them into the ground-state readout surface. Numerical "
            "hadron masses require one production backend export bundle with publication-complete manifest provenance, real `pi_iso`, `N_iso_direct`, and `N_iso_exchange` arrays, the executed runtime receipt `(N_therm, N_sep)`, and then declared production continuum/volume/chiral/statistical systematics. "
            "The surrogate execution bridge is only a diagnostic proof that the schema closes; it is not a promotable hadron output surface."
        ),
        "frontier_text": "Source-only scope requires real correlator arrays and declared statistical, continuum, volume, and chiral systematics. The empirical endpoint is data-driven and does not close the source-only hadron backend.",
        "prediction_surface": "Four fine-structure coordinates stay separate: certified source root 136.9948351774, source-root-plus-public-width diagnostic 137.0359595136, empirical endpoint 136.2636589431, and measured endpoint 137.035999177. No nonzero source-only hadron mass is emitted.",
        "particles": ["proton", "neutron", "neutral_pion", "rho_770_0"],
    },
]


def load_results(path: pathlib.Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_exact_nonhadron_bundle(path: pathlib.Path) -> Dict[str, Dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return {entry["particle_id"]: entry for entry in payload.get("entries", [])}


def wrap_text(text: str, width: int) -> List[str]:
    normalized = " ".join(text.split())
    return textwrap.wrap(normalized, width=width, break_long_words=False, break_on_hyphens=False)


def wrap_identifier(text: str, width: int) -> List[str]:
    chunks: List[str] = []
    token = ""
    for char in text:
        token += char
        if char in ".-_/":
            chunks.append(token)
            token = ""
    if token:
        chunks.append(token)

    lines: List[str] = []
    line = ""
    for chunk in chunks:
        if len(line) + len(chunk) <= width:
            line += chunk
            continue
        if line:
            lines.append(line)
            line = ""
        while len(chunk) > width:
            lines.append(chunk[:width])
            chunk = chunk[width:]
        line = chunk
    if line:
        lines.append(line)
    return lines


def char_capacity(width_px: float, font_size: int, *, mono: bool = False) -> int:
    ratio = 0.62 if mono else 0.56
    return max(12, int(width_px / (font_size * ratio)))


def wrap_paragraphs(
    paragraphs: Sequence[str],
    width_px: float,
    font_size: int,
    *,
    family: str = FONT_FAMILY,
) -> List[str]:
    mono = family == MONO_FAMILY
    width = char_capacity(width_px, font_size, mono=mono)
    lines: List[str] = []
    for index, paragraph in enumerate(paragraphs):
        if not paragraph:
            continue
        wrapped = wrap_text(paragraph, width)
        lines.extend(wrapped or [""])
    return lines


def render_wrapped_text(
    x: float,
    y: float,
    lines: Sequence[str],
    *,
    font_size: int,
    fill: str,
    family: str = FONT_FAMILY,
    weight: int | str = 400,
    line_height: int = 18,
) -> str:
    tspans = []
    for index, line in enumerate(lines):
        dy = 0 if index == 0 else line_height
        tspans.append(f'<tspan x="{x:.1f}" dy="{dy}">{escape(line or " ")}</tspan>')
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="{family}" font-size="{font_size}" '
        f'font-weight="{weight}" fill="{fill}">' + "".join(tspans) + "</text>"
    )


def box_layout(
    *,
    title: str,
    body: Sequence[str],
    w: float,
    title_size: int,
    body_size: int,
    body_family: str = FONT_FAMILY,
    prewrapped_body_lines: Sequence[str] | None = None,
    accent: bool = False,
) -> Dict[str, Any]:
    title_lines = wrap_text(title, char_capacity(w - 36, title_size))
    if prewrapped_body_lines is not None:
        body_lines = list(prewrapped_body_lines)
    else:
        body_lines = wrap_paragraphs(body, w - 36, body_size, family=body_family)

    title_lh = max(22, int(round(title_size * 1.2)))
    body_lh = max(18, int(round(body_size * 1.35)))
    title_y_offset = 42 if accent else 28
    body_gap = 14 if body_lines else 0
    body_y_offset = title_y_offset + len(title_lines) * title_lh + body_gap
    total_h = body_y_offset + len(body_lines) * body_lh + 24

    return {
        "title_lines": title_lines,
        "body_lines": body_lines,
        "title_lh": title_lh,
        "body_lh": body_lh,
        "title_y_offset": title_y_offset,
        "body_y_offset": body_y_offset,
        "height": total_h,
    }


def estimate_box_height(
    *,
    title: str,
    body: Sequence[str],
    w: float,
    title_size: int = 18,
    body_size: int = 14,
    body_family: str = FONT_FAMILY,
    prewrapped_body_lines: Sequence[str] | None = None,
    accent: bool = False,
) -> float:
    return box_layout(
        title=title,
        body=body,
        w=w,
        title_size=title_size,
        body_size=body_size,
        body_family=body_family,
        prewrapped_body_lines=prewrapped_body_lines,
        accent=accent,
    )["height"]


def draw_chip(
    x: float,
    y: float,
    text: str,
    *,
    fill: str,
    text_fill: str,
    stroke: str | None = None,
    font_size: int = 12,
) -> str:
    width = max(74.0, len(text) * (font_size * 0.63) + 22.0)
    stroke_attr = f' stroke="{stroke}" stroke-width="1.4"' if stroke else ""
    return (
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{width:.1f}" height="28" rx="14" fill="{fill}"{stroke_attr}/>'
        + render_wrapped_text(
            x + 12,
            y + 19,
            [text],
            font_size=font_size,
            fill=text_fill,
            weight=700,
            line_height=14,
        )
    )


def draw_section_label(x: float, y: float, text: str, *, fill: str, stroke: str, text_fill: str) -> str:
    width = max(170.0, len(text) * 7.2 + 26.0)
    return (
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{width:.1f}" height="28" rx="14" fill="{fill}" stroke="{stroke}" stroke-width="1.6"/>'
        + render_wrapped_text(
            x + 14,
            y + 19,
            [text],
            font_size=13,
            fill=text_fill,
            weight=700,
            line_height=15,
        )
    )


def draw_box(
    *,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    body: Sequence[str],
    fill: str,
    stroke: str,
    title_size: int = 18,
    body_size: int = 14,
    title_fill: str = COLORS["ink"],
    body_fill: str = COLORS["muted"],
    accent: str | None = None,
    badge: str | None = None,
    badge_fill: str | None = None,
    badge_text_fill: str = "#07101d",
    body_family: str = FONT_FAMILY,
    prewrapped_body_lines: Sequence[str] | None = None,
    corner: int = 18,
) -> str:
    layout = box_layout(
        title=title,
        body=body,
        w=w,
        title_size=title_size,
        body_size=body_size,
        body_family=body_family,
        prewrapped_body_lines=prewrapped_body_lines,
        accent=accent is not None,
    )
    parts = [
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{corner}" fill="{fill}" stroke="{stroke}" stroke-width="2"/>'
    ]
    if accent:
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="18" rx="{corner}" fill="{accent}"/>')
        if badge:
            parts.append(
                render_wrapped_text(
                    x + 12,
                    y + 13,
                    [badge],
                    font_size=11,
                    fill="#ffffff",
                    weight=700,
                    line_height=13,
                )
            )
    if badge and badge_fill:
        badge_width = max(96.0, len(badge) * 7.0 + 24.0)
        parts.append(
            f'<rect x="{x+w-badge_width-14:.1f}" y="{y+14:.1f}" width="{badge_width:.1f}" height="28" rx="14" '
            f'fill="{badge_fill}" stroke="none"/>'
        )
        parts.append(
            render_wrapped_text(
                x + w - badge_width - 2,
                y + 33,
                [badge],
                font_size=12,
                fill=badge_text_fill,
                weight=700,
                line_height=14,
            )
        )
    parts.append(
        render_wrapped_text(
            x + 18,
            y + layout["title_y_offset"],
            layout["title_lines"],
            font_size=title_size,
            fill=title_fill,
            weight=700,
            line_height=layout["title_lh"],
        )
    )
    if layout["body_lines"]:
        parts.append(
            render_wrapped_text(
                x + 18,
                y + layout["body_y_offset"],
                layout["body_lines"],
                font_size=body_size,
                fill=body_fill,
                family=body_family,
                line_height=layout["body_lh"],
            )
        )
    return "".join(parts)


def draw_polyline(
    points: Sequence[Tuple[float, float]],
    *,
    color: str,
    width: float = 2.2,
    dashed: bool = False,
    arrow: bool = False,
) -> str:
    dash = ' stroke-dasharray="9 7"' if dashed else ""
    marker = ' marker-end="url(#arrow)"' if arrow else ""
    path = " L ".join(f"{x:.1f} {y:.1f}" for x, y in points)
    return f'<path d="M {path}" fill="none" stroke="{color}" stroke-width="{width}"{dash}{marker}/>'


def draw_vertical_arrow(
    x: float,
    y1: float,
    y2: float,
    *,
    color: str = COLORS["line"],
    width: float = 2.2,
    dashed: bool = False,
) -> str:
    if y2 <= y1:
        return ""
    return draw_polyline([(x, y1), (x, y2)], color=color, width=width, dashed=dashed, arrow=True)


def _format_exact_output(entry: Dict[str, Any]) -> str:
    del entry
    return "withheld; no nonzero source-only physical mass is emitted"


def public_exact_surface(row: Dict[str, Any], exact_entry: Dict[str, Any]) -> str:
    particle_id = row["particle_id"]
    if particle_id in {"w_boson", "z_boson"}:
        return "compare-only exact sidecar"
    if particle_id == "higgs":
        return "conditional declared-surface Higgs/top candidate"
    if particle_id in {"electron", "muon", "tau"}:
        return "target-anchored same-family witness withheld from public prediction output"
    if particle_id.endswith("_quark"):
        return "mixed-convention target-audit coordinate withheld from public prediction output"
    if particle_id.endswith("_neutrino"):
        return "scale-free weighted-cycle branch with compare-only absolute attachment"
    return STATUS_EXPLAINER.get(row["status"], "tracked particle row")


def public_exact_caveat(row: Dict[str, Any], exact_entry: Dict[str, Any]) -> str:
    particle_id = row["particle_id"]
    if particle_id in {"w_boson", "z_boson"}:
        return "Target-containing inverse-adapter row. It is distinct from the selected carrier chart and conditional value law; physical promotion requires the source, scale, scheme, and complex-pole receipts."
    if particle_id == "higgs":
        return "Conditional Higgs candidate on the declared electroweak running and matching surface. It is not a source-only complex-pole prediction."
    if particle_id in {"electron", "muon", "tau"}:
        return "Target-anchored charged witness. The finite digital model proves schema satisfiability and a central-record bridge, while physical selection, determinant attachment, cofinal refinement, and the pole scheme remain open."
    if particle_id.endswith("_quark"):
        return "Target-audit coordinate only. The reciprocal-ray subfamily fails the generic common-scale six-Yukawa interface, and the historic mixed-convention chart is not a physical postdiction."
    if particle_id.endswith("_neutrino"):
        return "Scale-free weighted-cycle branch with a compare-only absolute attachment candidate; the absolute mass normalization is not promoted."
    return STATUS_NEXT_STEP.get(row["status"], "")


def particle_body(row: Dict[str, Any], exact_entry: Dict[str, Any] | None = None) -> List[str]:
    info = PARTICLE_INFO[row["particle_id"]]
    lines = [
        f"Symbol: {info['symbol']}. What it is: {info['plain']}",
        f"Why this row matters: {GROUP_ROW_TEXT[row['group']]}",
    ]
    if exact_entry is not None:
        lines.extend(
            [
                f"Surface class: {public_exact_surface(row, exact_entry)}.",
                "OPH output: no nonzero source-only physical mass emitted.",
                "Reference mass: withheld from this diagram ledger.",
                f"Caveat: {public_exact_caveat(row, exact_entry)}",
            ]
        )
    else:
        status = row["status"]
        lines.extend(
            [
                f"Surface class: {STATUS_TEXT[status]} meaning {STATUS_EXPLAINER[status]}.",
                "Reference mass: withheld from this diagram ledger.",
            ]
        )
        if status in {"structural", "conditional_carrier_mode"}:
            lines.insert(-1, "OPH output: zero hard quadratic parameter on the declared classical-action branch; no quantum-particle pole is claimed.")
        else:
            lines.insert(-1, "OPH output: no nonzero source-only physical mass emitted.")
        lines.append(f"Caveat: {STATUS_NEXT_STEP[status]}")
    return lines


def particle_card_height(row: Dict[str, Any], w: float, exact_entry: Dict[str, Any] | None = None) -> float:
    return estimate_box_height(
        title=PARTICLE_TITLE[row["particle_id"]],
        body=particle_body(row, exact_entry),
        w=w,
        title_size=18,
        body_size=15,
        accent=True,
    )


def draw_particle_card(
    row: Dict[str, Any],
    x: float,
    y: float,
    w: float,
    exact_entry: Dict[str, Any] | None = None,
) -> Tuple[str, float]:
    badge_text = "not promoted" if exact_entry is not None and not exact_entry.get("promotable", False) else STATUS_TEXT[row["status"]]
    badge_color = COLORS["task_pending"] if exact_entry is not None and not exact_entry.get("promotable", False) else (COLORS["task_complete"] if exact_entry is not None else STATUS_BAR[row["status"]])
    height = particle_card_height(row, w, exact_entry)
    return (
        draw_box(
            x=x,
            y=y,
            w=w,
            h=height,
            title=PARTICLE_TITLE[row["particle_id"]],
            body=particle_body(row, exact_entry),
            fill=COLORS["prediction_fill"],
            stroke=COLORS["prediction_stroke"],
            title_size=18,
            body_size=15,
            accent=badge_color,
            badge=badge_text,
            body_fill="#d4ffe5",
        ),
        height,
    )


def particle_grid(count: int) -> int:
    return 2 if count >= 4 else 1


def lane_particle_ids(lane: Dict[str, Any], rows_by_id: Dict[str, Dict[str, Any]]) -> List[str]:
    return [particle_id for particle_id in lane["particles"] if particle_id in rows_by_id]


def measure_particle_section(lane: Dict[str, Any], rows_by_id: Dict[str, Dict[str, Any]], w: float) -> float:
    particle_ids = lane_particle_ids(lane, rows_by_id)
    if not particle_ids:
        return estimate_box_height(
            title="Public particle rows hidden",
            body=["This lane is off the public particle table because no closure-grade predictor is declared on its public surface."],
            w=w,
            title_size=18,
            body_size=15,
        )
    cols = particle_grid(len(particle_ids))
    gap = 14.0
    card_w = w if cols == 1 else (w - gap) / 2.0
    heights = [particle_card_height(rows_by_id[particle_id], card_w) for particle_id in particle_ids]
    if cols == 1:
        return sum(heights) + gap * max(0, len(heights) - 1)
    total = 0.0
    for index in range(0, len(heights), 2):
        pair = heights[index : index + 2]
        total += max(pair)
        if index + 2 < len(heights):
            total += gap
    return total


def draw_particle_section(
    lane: Dict[str, Any],
    rows_by_id: Dict[str, Dict[str, Any]],
    exact_by_id: Dict[str, Dict[str, Any]],
    x: float,
    y: float,
    w: float,
) -> Tuple[str, float]:
    particle_ids = lane_particle_ids(lane, rows_by_id)
    if not particle_ids:
        height = estimate_box_height(
            title="Public particle rows hidden",
            body=["This lane is off the public particle table because no closure-grade predictor is declared on its public surface."],
            w=w,
            title_size=18,
            body_size=15,
        )
        return (
            draw_box(
                x=x,
                y=y,
                w=w,
                h=height,
                title="Public particle rows hidden",
                body=["This lane is off the public particle table because no closure-grade predictor is declared on its public surface."],
                fill=COLORS["green_note_fill"],
                stroke=COLORS["green_note_stroke"],
                title_size=18,
                body_size=15,
                title_fill=COLORS["green_note_text"],
                body_fill=COLORS["green_note_body"],
            ),
            height,
        )
    cols = particle_grid(len(particle_ids))
    gap = 14.0
    card_w = w if cols == 1 else (w - gap) / 2.0
    parts: List[str] = []

    if cols == 1:
        cursor = y
        for particle_id in particle_ids:
            exact_entry = exact_by_id.get(particle_id)
            card, height = draw_particle_card(rows_by_id[particle_id], x, cursor, card_w, exact_entry)
            parts.append(card)
            cursor += height + gap
        return "".join(parts), cursor - y - gap

    cursor = y
    particles = [rows_by_id[particle_id] for particle_id in particle_ids]
    for index in range(0, len(particles), 2):
        left_row = particles[index]
        left_exact = exact_by_id.get(left_row["particle_id"])
        left_h = particle_card_height(left_row, card_w, left_exact)
        left_card, _ = draw_particle_card(left_row, x, cursor, card_w, left_exact)
        parts.append(left_card)

        row_height = left_h
        if index + 1 < len(particles):
            right_row = particles[index + 1]
            right_exact = exact_by_id.get(right_row["particle_id"])
            right_h = particle_card_height(right_row, card_w, right_exact)
            right_card, _ = draw_particle_card(right_row, x + card_w + gap, cursor, card_w, right_exact)
            parts.append(right_card)
            row_height = max(left_h, right_h)

        cursor += row_height + gap
    return "".join(parts), cursor - y - gap


def lane_status(lane: Dict[str, Any], rows_by_id: Dict[str, Dict[str, Any]]) -> str:
    particle_ids = lane_particle_ids(lane, rows_by_id)
    if particle_ids:
        return rows_by_id[particle_ids[0]]["status"]
    if "status" in lane:
        return lane["status"]
    return "simulation_dependent"


def lane_panel_height(
    lane: Dict[str, Any],
    rows_by_id: Dict[str, Dict[str, Any]],
    exact_by_id: Dict[str, Dict[str, Any]],
    w: float,
) -> float:
    inner_w = w - 36
    summary_lines = wrap_text(lane["summary"], char_capacity(inner_w - 18, 16))
    takeaway_lines = wrap_text(lane["takeaway"], char_capacity(inner_w - 18, 15))
    header_h = 118 + len(summary_lines) * 20 + len(takeaway_lines) * 18
    label_h = 28
    label_gap = 10
    section_gap = 20

    theorem_h = estimate_box_height(
        title=lane["summary"],
        body=[lane["logic"]],
        w=inner_w,
        title_size=20,
        body_size=15,
    )
    task_total_h = estimate_box_height(
        title="Claim frontier",
        body=[lane["frontier_text"]],
        w=inner_w,
        title_size=18,
        body_size=15,
    )
    prediction_h = estimate_box_height(
        title="Prediction surface",
        body=[lane["prediction_surface"]],
        w=inner_w,
        title_size=18,
        body_size=15,
    )
    particle_h = measure_particle_section(lane, rows_by_id, inner_w)

    return (
        header_h
        + label_h
        + label_gap
        + theorem_h
        + section_gap
        + label_h
        + label_gap
        + task_total_h
        + section_gap
        + label_h
        + label_gap
        + prediction_h
        + section_gap
        + label_h
        + label_gap
        + particle_h
        + 22
    )


def draw_lane_panel(
    lane: Dict[str, Any],
    rows_by_id: Dict[str, Dict[str, Any]],
    exact_by_id: Dict[str, Dict[str, Any]],
    x: float,
    y: float,
    w: float,
) -> Tuple[str, float]:
    inner_x = x + 18
    inner_w = w - 36
    status = lane_status(lane, rows_by_id)
    status_color = STATUS_BAR[status]
    summary_lines = wrap_text(lane["summary"], char_capacity(inner_w - 18, 16))
    takeaway_lines = wrap_text(lane["takeaway"], char_capacity(inner_w - 18, 15))
    header_h = 118 + len(summary_lines) * 20 + len(takeaway_lines) * 18
    label_h = 28
    label_gap = 10
    section_gap = 20

    theorem_h = estimate_box_height(
        title=lane["summary"],
        body=[lane["logic"]],
        w=inner_w,
        title_size=20,
        body_size=15,
    )
    task_total_h = estimate_box_height(
        title="Claim frontier",
        body=[lane["frontier_text"]],
        w=inner_w,
        title_size=18,
        body_size=15,
    )
    prediction_h = estimate_box_height(
        title="Prediction surface",
        body=[lane["prediction_surface"]],
        w=inner_w,
        title_size=18,
        body_size=15,
    )
    particle_h = measure_particle_section(lane, rows_by_id, inner_w)

    panel_h = (
        header_h
        + label_h
        + label_gap
        + theorem_h
        + section_gap
        + label_h
        + label_gap
        + task_total_h
        + section_gap
        + label_h
        + label_gap
        + prediction_h
        + section_gap
        + label_h
        + label_gap
        + particle_h
        + 22
    )

    parts: List[str] = [
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{panel_h:.1f}" rx="24" fill="{COLORS["panel_shell"]}" stroke="{COLORS["panel_border"]}" stroke-width="2"/>',
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="10" rx="24" fill="{status_color}"/>',
        f'<rect x="{x:.1f}" y="{y+10:.1f}" width="{w:.1f}" height="10" fill="{status_color}" opacity="0.22"/>',
    ]

    parts.append(
        render_wrapped_text(
            inner_x,
            y + 42,
            [lane["title"]],
            font_size=28,
            fill=COLORS["ink"],
            weight=700,
            line_height=30,
        )
    )
    parts.append(
        render_wrapped_text(
            inner_x,
            y + 74,
            summary_lines,
            font_size=16,
            fill=COLORS["muted"],
            line_height=20,
        )
    )
    summary_bottom_y = y + 74 + len(summary_lines) * 20
    parts.append(
        render_wrapped_text(
            inner_x,
            summary_bottom_y + 14,
            takeaway_lines,
            font_size=15,
            fill="#c7f0ff",
            line_height=18,
        )
    )
    visible_particles = lane_particle_ids(lane, rows_by_id)
    parts.append(
        render_wrapped_text(
            inner_x,
            summary_bottom_y + 18 + len(takeaway_lines) * 18 + 12,
            [
                (
                    f"{len(visible_particles)} tracked carrier-mode gate rows shown in this lane"
                    if lane.get("key") == "structural"
                    else f"{len(visible_particles)} tracked particle rows shown in this lane"
                )
            ],
            font_size=14,
            fill=COLORS["subtle"],
            line_height=16,
        )
    )
    parts.append(
        draw_chip(
            x + w - 136,
            y + 18,
            STATUS_TEXT[status],
            fill=status_color,
            text_fill="#08111b",
        )
    )

    cursor = y + header_h

    parts.append(
        draw_section_label(
            inner_x,
            cursor,
            "Implemented derivation / audit",
            fill="#0d2232",
            stroke=COLORS["logic_stroke"],
            text_fill=COLORS["ink"],
        )
    )
    cursor += label_h + label_gap
    theorem_y = cursor
    parts.append(
        draw_box(
            x=inner_x,
            y=theorem_y,
            w=inner_w,
            h=theorem_h,
            title=lane["summary"],
            body=[lane["logic"]],
            fill=COLORS["logic_fill"],
            stroke=COLORS["logic_stroke"],
            title_size=20,
            body_size=15,
        )
    )
    cursor += theorem_h + section_gap

    parts.append(
        draw_section_label(
            inner_x,
            cursor,
            "Named frontier objects",
            fill="#241611",
            stroke=COLORS["task_stroke"],
            text_fill=COLORS["ink"],
        )
    )
    cursor += label_h + label_gap

    no_task_h = estimate_box_height(
        title="Claim frontier",
        body=[lane["frontier_text"]],
        w=inner_w,
        title_size=18,
        body_size=15,
    )
    parts.append(draw_vertical_arrow(x + w / 2.0, theorem_y + theorem_h + 4, cursor - 6, color=COLORS["green_note_stroke"], dashed=True))
    parts.append(
        draw_box(
            x=inner_x,
            y=cursor,
            w=inner_w,
            h=no_task_h,
            title="Claim frontier",
            body=[lane["frontier_text"]],
            fill=COLORS["green_note_fill"],
            stroke=COLORS["green_note_stroke"],
            title_size=18,
            body_size=15,
            title_fill=COLORS["green_note_text"],
            body_fill=COLORS["green_note_body"],
        )
    )
    task_bottom = cursor + no_task_h
    cursor += no_task_h

    cursor += section_gap
    parts.append(
        draw_section_label(
            inner_x,
            cursor,
            "Claim / prediction surface",
            fill="#161f34",
            stroke=COLORS["status_stroke"],
            text_fill=COLORS["ink"],
        )
    )
    cursor += label_h + label_gap
    prediction_y = cursor
    parts.append(
        draw_vertical_arrow(x + w / 2.0, task_bottom + 4, prediction_y - 6, color=COLORS["line_soft"]))
    parts.append(
        draw_box(
            x=inner_x,
            y=prediction_y,
            w=inner_w,
            h=prediction_h,
            title="Claim / prediction surface",
            body=[lane["prediction_surface"]],
            fill=COLORS["status_fill"],
            stroke=COLORS["status_stroke"],
            title_size=18,
            body_size=15,
        )
    )
    cursor += prediction_h + section_gap

    parts.append(
        draw_section_label(
            inner_x,
            cursor,
            (
                "Tracked carrier-mode gates"
                if lane.get("key") == "structural"
                else "Tracked particle outputs"
            ),
            fill="#122219",
            stroke=COLORS["prediction_stroke"],
            text_fill=COLORS["ink"],
        )
    )
    cursor += label_h + label_gap
    particle_y = cursor
    parts.append(
        draw_vertical_arrow(x + w / 2.0, prediction_y + prediction_h + 4, particle_y - 6, color=COLORS["prediction_stroke"]))
    particle_markup, particle_total_h = draw_particle_section(lane, rows_by_id, exact_by_id, inner_x, particle_y, inner_w)
    parts.append(particle_markup)

    return "".join(parts), panel_h


def build_svg(results: Dict[str, Any], exact_by_id: Dict[str, Dict[str, Any]]) -> str:
    rows_by_id = {row["particle_id"]: row for row in results["rows"]}
    generated_utc = results["generated_utc"]
    built_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    counts = Counter(row["status"] for row in results["rows"])
    closedish = (
        counts["structural"]
        + counts["calibration"]
        + counts["secondary_quantitative"]
        + counts["selected_class_theorem"]
        + counts["source_spread_nonidentifiability_obstruction"]
    )
    total_rows = len(results["rows"])

    parts: List[str] = []
    current_y = 48.0

    # Header sizes
    intro_lines = [
        "A top-down map of the OPH particle program: start at the fixed-point closures, move through the sector lanes, and finish at the particle rows printed by the declared code surfaces.",
        "Each carrier must instantiate a bounded observer-like self-reading patch with local state, boundaries or ports, readback, records, repair moves, and a public evidence bundle.",
        "Blue and green regions show implemented derivation surfaces. Orange cards mark named theorem or computation burdens on the public lanes.",
    ]
    intro_height = 3 * 24

    legend_items = [
        ("Implemented theorem / technique", COLORS["logic_stroke"]),
        ("Named frontier object", COLORS["task_stroke"]),
        ("Particle end node", COLORS["prediction_stroke"]),
        ("Connector / dependency flow", COLORS["line"]),
    ]
    legend_w = 420.0
    legend_h = 152.0
    legend_x = WIDTH - MARGIN_X - legend_w
    legend_y = current_y

    header_bottom = max(current_y + 112 + intro_height, legend_y + legend_h)
    current_y = header_bottom + 30

    # Fixed-point closure section
    input_label_y = current_y
    current_y += 28

    input_gap = 24.0
    axiom_w = 620.0
    p_w = 360.0
    other_w = WIDTH - 2 * MARGIN_X - axiom_w - p_w - 2 * input_gap

    inputs = results["inputs"]
    input_specs = [
        {
            "x": MARGIN_X,
            "w": axiom_w,
            "title": "Five OPH Axioms",
            "body": ["This chart treats the OPH axioms as the common starting point. They provide the conceptual constraints upstream of every lane below."],
            "fill": COLORS["axiom_fill"],
            "stroke": COLORS["axiom_stroke"],
        },
        {
            "x": MARGIN_X + axiom_w + input_gap,
            "w": p_w,
            "title": "Local Fixed Point: P",
            "body": [f"P = {PUBLIC_PIXEL_DISPLAY}. This scalar is the local pixel ratio selected on the outer/inner closure surface and shared by the electroweak, flavor, and hadron lanes."],
            "fill": COLORS["input_fill"],
            "stroke": COLORS["input_stroke"],
        },
        {
            "x": MARGIN_X + axiom_w + input_gap + p_w + input_gap,
            "w": other_w,
            "title": "Global Fixed Point: N_CRC",
            "body": [f"N_CRC = {PUBLIC_CAPACITY_DISPLAY} is the conditional EW bridge coordinate. Physical capacity and the Lambda comparison remain open until F and CP-1 to CP-3 close and the joint posterior is propagated. Legacy neutrino side settings remain separate from the rejected weighted-cycle comparison lane. loops = {inputs['loops']}; hadron_profile = {inputs['hadron_profile']}."],
            "fill": COLORS["input_fill"],
            "stroke": COLORS["input_stroke"],
        },
    ]
    input_h = max(
        estimate_box_height(title=spec["title"], body=spec["body"], w=spec["w"], title_size=22, body_size=16)
        for spec in input_specs
    )
    input_y = current_y
    for spec in input_specs:
        parts.append(
            draw_box(
                x=spec["x"],
                y=input_y,
                w=spec["w"],
                h=input_h,
                title=spec["title"],
                body=spec["body"],
                fill=spec["fill"],
                stroke=spec["stroke"],
                title_size=22,
                body_size=16,
            )
        )
    current_y += input_h + 34

    # Scaffold section
    scaffold_label_y = current_y
    current_y += 28
    scaffold_y = current_y
    scaffold_w = WIDTH - 2 * MARGIN_X
    scaffold_body = [
        "Start with the declared P maps and the conditional N_CRC bridge. The physical Thomson map, F, and CP-1 to CP-3 remain open; no unconditional capacity closure is claimed.",
        "Then read each lane from top to bottom: implemented theorem content, named frontier objects, prediction surface, and the particle rows shown on the public table.",
        f"The badge reports {closedish} of {total_rows} tracked rows above continuation / simulation status. The quark count preserves a legacy machine status; the reader-facing lane records the common-scale reciprocal-ray rejection.",
        "The broader geometric-premise boundary sits above the particle lanes. Three cap-pair extraction witnesses are explicit. The open geometric clause concerns a common floor for the finitely many modular-transport marginals, followed by ordered null cut-pair rigidity.",
    ]
    scaffold_h = estimate_box_height(
        title="How to read the mass derivation chart",
        body=scaffold_body,
        w=scaffold_w,
        title_size=24,
        body_size=16,
    )
    parts.append(
        draw_box(
            x=MARGIN_X,
            y=scaffold_y,
            w=scaffold_w,
            h=scaffold_h,
            title="How to read the mass derivation chart",
            body=scaffold_body,
            fill=COLORS["panel_alt"],
            stroke="#31517c",
            title_size=24,
            body_size=16,
            badge=f"{closedish} / {total_rows} higher-closure rows",
            badge_fill="#7dd3fc",
            badge_text_fill="#07101d",
        )
    )
    current_y += scaffold_h + 38

    # Panels section header
    lanes_label_y = current_y
    current_y += 92

    row1 = LANES[:3]
    row2 = LANES[3:6]
    row3 = [LANES[6]]

    row1_xs = [MARGIN_X + index * (PANEL_W + PANEL_GAP_X) for index in range(3)]
    row2_xs = row1_xs
    hadron_x = (WIDTH - HADRON_W) / 2.0

    row1_heights = [lane_panel_height(lane, rows_by_id, exact_by_id, PANEL_W) for lane in row1]
    row1_y = current_y
    row1_h = max(row1_heights)

    row2_heights = [lane_panel_height(lane, rows_by_id, exact_by_id, PANEL_W) for lane in row2]
    row2_y = row1_y + row1_h + PANEL_GAP_Y
    row2_h = max(row2_heights)

    row3_h = lane_panel_height(row3[0], rows_by_id, exact_by_id, HADRON_W)
    row3_y = row2_y + row2_h + PANEL_GAP_Y

    total_panels_bottom = row3_y + row3_h

    footer_label_y = total_panels_bottom + 40
    footer_y = footer_label_y + 28
    footer_gap = 20.0
    footer_w = (WIDTH - 2 * MARGIN_X - 2 * footer_gap) / 3.0
    footer_items = [
        (
            "Status colors",
            [
                "classical carrier = zero hard quadratic parameter on a declared branch; quantum-particle gate open",
                "electroweak frontier = selected chart, conditional value law, and measured-reference inverse adapter remain separate",
                "secondary = quantitative branch built on a declared electroweak layer",
                "quark boundary = common-scale six-Yukawa data reject the reciprocal ray; numeric rows withheld",
                "continuation = declared continuation or witness surface",
                "simulation = execution-bound lane with backend and systematics prerequisites",
            ],
        ),
        (
            "Plain-English terms",
            [
                "Gauge boson = force carrier.",
                "Lepton = electron-like matter particle; neutrinos are neutral leptons.",
                "Quark = an elementary constituent used to build hadrons such as protons and neutrons.",
                "RG running / matching = translating couplings and masses between energy scales and schemes.",
                "Quenched = a simplified QCD simulation that omits some vacuum-quark effects.",
            ],
        ),
        (
            "Bottom line",
            [
                f"Status table generated: {generated_utc}. SVG generated: {built_utc}.",
                "This chart maps the implemented derivation pipeline and its claim boundaries. It emits no nonzero source-only physical particle mass. The full particle zoo does not sit on one uniform theorem chain from P plus axioms.",
            ],
        ),
    ]
    footer_h = max(
        estimate_box_height(title=title, body=body, w=footer_w, title_size=19, body_size=15) for title, body in footer_items
    )
    total_height = footer_y + footer_h + 56

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{total_height:.0f}" viewBox="0 0 {WIDTH} {total_height:.0f}">',
        "<defs>",
        '<linearGradient id="bg-grad" x1="0" y1="0" x2="0" y2="1">',
        f'<stop offset="0%" stop-color="{COLORS["bg_alt"]}"/>',
        f'<stop offset="100%" stop-color="{COLORS["bg"]}"/>',
        "</linearGradient>",
        '<marker id="arrow" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto">',
        f'<path d="M 0 0 L 12 6 L 0 12 z" fill="{COLORS["line"]}"/>',
        "</marker>",
        "</defs>",
        f'<rect width="{WIDTH}" height="{total_height:.0f}" fill="url(#bg-grad)"/>',
        f'<circle cx="{WIDTH-120}" cy="160" r="360" fill="#0ea5e9" opacity="0.08"/>',
        f'<circle cx="220" cy="{int(total_height-420)}" r="340" fill="#8b5cf6" opacity="0.06"/>',
        f'<circle cx="{WIDTH-220}" cy="{int(total_height*0.62)}" r="280" fill="#22c55e" opacity="0.05"/>',
    ]

    parts.append(
        render_wrapped_text(
            MARGIN_X,
            52,
            ["OPH Particle-Mass Pipeline"],
            font_size=50,
            fill=COLORS["ink"],
            weight=700,
            line_height=52,
        )
    )
    parts.append(
        render_wrapped_text(
            MARGIN_X,
            104,
            intro_lines,
            font_size=18,
            fill=COLORS["muted"],
            line_height=24,
        )
    )

    parts.append(
        f'<rect x="{legend_x:.1f}" y="{legend_y:.1f}" width="{legend_w:.1f}" height="{legend_h:.1f}" rx="20" fill="{COLORS["panel"]}" stroke="{COLORS["panel_border"]}" stroke-width="2"/>'
    )
    parts.append(
        render_wrapped_text(
            legend_x + 20,
            legend_y + 28,
            ["Legend"],
            font_size=18,
            fill=COLORS["ink"],
            weight=700,
            line_height=20,
        )
    )
    for index, (label, color) in enumerate(legend_items):
        yy = legend_y + 58 + index * 22
        parts.append(f'<rect x="{legend_x+20:.1f}" y="{yy-11:.1f}" width="14" height="14" rx="4" fill="{color}"/>')
        parts.append(
            render_wrapped_text(
                legend_x + 42,
                yy,
                [label],
                font_size=14,
                fill=COLORS["muted"],
                line_height=16,
            )
        )

    parts.append(render_wrapped_text(MARGIN_X, input_label_y + 20, ["1. Fixed-point closures"], font_size=24, fill=COLORS["ink"], weight=700, line_height=26))
    for spec in input_specs:
        parts.append(
            draw_box(
                x=spec["x"],
                y=input_y,
                w=spec["w"],
                h=input_h,
                title=spec["title"],
                body=spec["body"],
                fill=spec["fill"],
                stroke=spec["stroke"],
                title_size=22,
                body_size=16,
            )
        )

    parts.append(render_wrapped_text(MARGIN_X, scaffold_label_y + 20, ["2. Branching logic"], font_size=24, fill=COLORS["ink"], weight=700, line_height=26))
    parts.append(
        draw_box(
            x=MARGIN_X,
            y=scaffold_y,
            w=scaffold_w,
            h=scaffold_h,
            title="How to read the mass derivation chart",
            body=scaffold_body,
            fill=COLORS["panel_alt"],
            stroke="#31517c",
            title_size=24,
            body_size=16,
            badge=f"{closedish} / {total_rows} higher-closure rows",
            badge_fill="#7dd3fc",
            badge_text_fill="#07101d",
        )
    )

    parts.append(render_wrapped_text(MARGIN_X, lanes_label_y + 20, ["3. Sector-specific derivation lanes"], font_size=24, fill=COLORS["ink"], weight=700, line_height=26))

    # Shared branch connectors
    trunk_x = WIDTH / 2.0
    row1_centers = [x + PANEL_W / 2.0 for x in row1_xs]
    row2_centers = [x + PANEL_W / 2.0 for x in row2_xs]
    hadron_center = hadron_x + HADRON_W / 2.0
    bus1_y = row1_y - 22
    bus2_y = row2_y - 22
    bus3_y = row3_y - 22

    scaffold_center_y = scaffold_y + scaffold_h
    input_centers = [spec["x"] + spec["w"] / 2.0 for spec in input_specs]
    input_bus_y = input_y + input_h + 18
    parts.append(draw_polyline([(input_centers[0], input_bus_y), (input_centers[-1], input_bus_y)], color=COLORS["line"], width=2.2))
    for center in input_centers:
        parts.append(draw_polyline([(center, input_y + input_h), (center, input_bus_y)], color=COLORS["line"], width=2.2))
    parts.append(draw_polyline([(trunk_x, input_bus_y), (trunk_x, scaffold_y - 8)], color=COLORS["line"], width=2.4, arrow=True))
    parts.append(draw_polyline([(trunk_x, scaffold_center_y + 10), (trunk_x, bus1_y)], color=COLORS["line"], width=2.4))
    parts.append(draw_polyline([(row1_centers[0], bus1_y), (row1_centers[-1], bus1_y)], color=COLORS["line"], width=2.4))
    for center in row1_centers:
        parts.append(draw_vertical_arrow(center, bus1_y, row1_y - 6, color=COLORS["line"]))

    parts.append(draw_polyline([(trunk_x, bus1_y), (trunk_x, bus2_y)], color=COLORS["line"], width=2.2))
    parts.append(draw_polyline([(row2_centers[0], bus2_y), (row2_centers[-1], bus2_y)], color=COLORS["line"], width=2.2))
    for center in row2_centers:
        parts.append(draw_vertical_arrow(center, bus2_y, row2_y - 6, color=COLORS["line"]))

    parts.append(draw_polyline([(trunk_x, bus2_y), (trunk_x, bus3_y)], color=COLORS["line"], width=2.2))
    parts.append(draw_vertical_arrow(hadron_center, bus3_y, row3_y - 6, color=COLORS["line"]))

    # Panels
    for lane, x in zip(row1, row1_xs):
        panel_markup, _ = draw_lane_panel(lane, rows_by_id, exact_by_id, x, row1_y, PANEL_W)
        parts.append(panel_markup)
    for lane, x in zip(row2, row2_xs):
        panel_markup, _ = draw_lane_panel(lane, rows_by_id, exact_by_id, x, row2_y, PANEL_W)
        parts.append(panel_markup)
    hadron_markup, _ = draw_lane_panel(row3[0], rows_by_id, exact_by_id, hadron_x, row3_y, HADRON_W)
    parts.append(hadron_markup)

    # Footer
    parts.append(render_wrapped_text(MARGIN_X, footer_label_y + 20, ["4. Glossary / status key"], font_size=24, fill=COLORS["ink"], weight=700, line_height=26))
    for index, (title, body) in enumerate(footer_items):
        parts.append(
            draw_box(
                x=MARGIN_X + index * (footer_w + footer_gap),
                y=footer_y,
                w=footer_w,
                h=footer_h,
                title=title,
                body=body,
                fill=COLORS["footer_fill"],
                stroke=COLORS["footer_stroke"],
                title_size=19,
                body_size=15,
            )
        )

    parts.append("</svg>")
    return "\n".join(parts) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate the particle mass derivation SVG.")
    parser.add_argument("--results-json", default=str(RESULTS_JSON), help="Path to the results JSON.")
    parser.add_argument("--exact-nonhadron-json", default=str(EXACT_NONHADRON_JSON), help="Path to the exact non-hadron bundle JSON.")
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT), help="Output SVG path.")
    args = parser.parse_args()

    results = load_results(pathlib.Path(args.results_json))
    exact_by_id = load_exact_nonhadron_bundle(pathlib.Path(args.exact_nonhadron_json))
    svg = build_svg(results=results, exact_by_id=exact_by_id)

    output = pathlib.Path(args.output)
    output.write_text(svg, encoding="utf-8")
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
