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

import yaml


ROOT = pathlib.Path(__file__).resolve().parents[2]
RESULTS_JSON = ROOT / "particles" / "results_status.json"
TASK_TRACKER_YAML = ROOT / "particles" / "task_tracker.yaml"
EXACT_NONHADRON_JSON = ROOT / "particles" / "exact_nonhadron_masses.json"
DEFAULT_OUTPUT = ROOT / "particles" / "particle_mass_derivation_graph.svg"

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
    "calibration": "#06b6d4",
    "secondary_quantitative": "#8b5cf6",
    "continuation": "#f59e0b",
    "simulation_dependent": "#ef4444",
}

STATUS_TEXT = {
    "structural": "structural",
    "calibration": "calibration",
    "secondary_quantitative": "secondary",
    "continuation": "continuation",
    "simulation_dependent": "simulation",
}

PARTICLE_INFO: Dict[str, Dict[str, str]] = {
    "photon": {"symbol": "gamma", "plain": "Particle of light; the electromagnetic force carrier."},
    "gluon": {"symbol": "g", "plain": "Strong-force carrier. Free gluons are not seen because they are confined."},
    "graviton": {"symbol": "graviton", "plain": "Massless spin-2 quantum of the OPH dynamical-metric branch."},
    "w_boson": {"symbol": "W", "plain": "Charged weak-force boson used in beta-decay-type processes."},
    "z_boson": {"symbol": "Z", "plain": "Neutral weak-force boson from the same electroweak sector as the W."},
    "higgs": {"symbol": "H", "plain": "Higgs boson tied to the Standard Model mass-giving field."},
    "top_quark": {"symbol": "t", "plain": "Heaviest known quark. In the live table the main row comes from D11."},
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
    "neutral_pion": {"symbol": "pi0 proxy", "plain": "Light meson placeholder row kept suppressed until the active hadron pipeline closes."},
    "rho_770_0": {"symbol": "rho(770)0 proxy", "plain": "Vector-meson placeholder row kept suppressed until the active hadron pipeline closes."},
}

GROUP_ROW_TEXT = {
    "Bosons": "This row sits on the force-carrier / Higgs side of the spectrum.",
    "Leptons": "This row sits on the lepton and neutrino side of the spectrum.",
    "Quarks": "This row sits on the quark flavor side of the spectrum.",
    "Hadrons": "This row sits on the bound-state QCD side of the spectrum.",
}

STATUS_EXPLAINER = {
    "structural": "exact on the structural theorem surface",
    "calibration": "closed after fixing the shared calibration scale P",
    "secondary_quantitative": "quantitative and usable, but built on a later secondary branch",
    "continuation": "a live continuation output, not the final theorem-grade endpoint",
    "simulation_dependent": "execution-bound and not a paper-only derived output",
}

STATUS_NEXT_STEP = {
    "structural": "This row is effectively finished; the remaining work is presentation and audit.",
    "calibration": "This row is closed on the active theorem surface; the remaining work is mostly audit and synchronization.",
    "secondary_quantitative": "This row is quantitatively live, but the surrounding proof package still needs hardening.",
    "continuation": "Treat this as a real waypoint rather than a final answer: the theorem chain still has one or more open objects.",
    "simulation_dependent": "The bottleneck here is real computation and systematics rather than one more symbolic derivation.",
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
        "title": "Structural Massless Carriers",
        "summary": "Massless photon, gluon, and graviton sectors from the structural gauge and gravity branches.",
        "takeaway": "These are the cleanest rows in the whole poster: OPH says these carrier sectors stay massless, so the outputs are exact zeros rather than fitted near-misses.",
        "logic": (
            "Use the realized gauge/content branch for photon and gluons, and the dynamical-metric "
            "Einstein branch for the graviton. In plain English: OPH says the electromagnetic, color, "
            "and spin-2 carrier sectors remain massless at the structural level because their mass terms "
            "would break required gauge or diffeomorphism redundancies."
        ),
        "tasks_text": "No open blocker is attached to these structural massless rows in the ledger.",
        "prediction_surface": "Structural massless carrier output surface.",
        "particles": ["photon", "gluon", "graviton"],
        "tasks": [],
    },
    {
        "key": "d10",
        "title": "D10 Electroweak Calibration",
        "summary": "The D10 lane implements the single-P running family, the reduced two-scalar carrier, the exact current-carrier chart, the historical freeze-once validation surface, and the promoted target-free source-only electroweak repair theorem.",
        "takeaway": "This is the most finished numerical lane: once the shared scale P is fixed on the declared surface, the active theorem emits W and Z forward without reusing the target masses.",
        "logic": (
            "From P the code builds M_U, solves alpha_U from the pixel constraint, gets the electroweak "
            "scale v, runs couplings to mZ_run, emits the running-family observables, reduces them to the "
            "two-scalar `(sigma_EW, eta_EW)` carrier, and applies the carrier selector J_pop_EW. The emitted "
            "mass pair of the selected carrier remains explicit on disk, and the neutral running pair closes through "
            "the derived source-normalized hypercharge readout `chi_Y_EW`. The emitted fiberwise tree law isolates "
            "the charged-leg coordinate `tau2_tree_exact`, and the exact current-carrier mass chart closes on "
            "`(tau2_tree_exact, delta_n_tree_exact)`. On top of that chart, the promoted theorem "
            "`EWTargetFreeRepairValueLaw_D10` emits the same repaired chart directly from the D10 source basis "
            "through `lambda_EW = eta_source^2 / (4 * beta_EW)`, then emits one repaired coupling pair and one "
            "coherent electroweak quintet. The older freeze-once coherent repair law is retained only as compare-only "
            "validation and agrees with the target-free theorem to machine scale."
        ),
        "tasks_text": "No open electroweak mass-side blocker remains on the active D10 calibration surface. The remaining work is audit, presentation, and synchronization of the promoted target-free theorem across the OPH papers.",
        "prediction_surface": "Target-free source-only D10 repaired electroweak quintet, with the exact current-carrier chart and the freeze-once repair law retained as separate historical validation objects.",
        "particles": ["w_boson", "z_boson"],
        "tasks": [
            "particles.calibration.02-separate-p-resolution-from-d10-transport-mismatch",
            "particles.calibration.04-push-single-p-exact-ew-closure-branch",
        ],
    },
    {
        "key": "d11",
        "title": "D11 Higgs / Top Branch",
        "summary": "The D11 lane maps the D10 gauge core into a compact shared seed and then reads out Higgs and top masses through the D11 Jacobian surface.",
        "takeaway": "This is a compact heavy-particle lane: one small seed controls the Higgs and top outputs rather than a large new family-specific construction.",
        "logic": (
            "Take the D10 substrate, impose the critical-surface condition, then use the synchronized "
            "transport core to emit Higgs/top outputs. The reduced readback closes on the one-scalar "
            "seed `sigma_D11_HT`, and the forward certificate records `pi_y = pi_lambda`, `eta_HT = 0`, "
            "and `w_HT = 0` on the live branch. The public `m_H` and `m_t` rows use that promoted seed."
        ),
        "tasks_text": "No open D11 mass-side blocker remains on the live forward seed path.",
        "prediction_surface": "D11 exact promoted forward-seed plus Jacobian readout surface.",
        "particles": ["higgs", "top_quark"],
        "tasks": ["papers.compact.e.28-make-the-higgs-top-critical-surface-fully-rigorous"],
    },
    {
        "key": "leptons",
        "title": "Charged Leptons",
        "summary": "The charged-lepton lane contains a closed exact same-family readout on the fixed ordered carrier, while the theorem-grade absolute lane narrows to one branch-generator clause and one post-promotion affine scalar.",
        "takeaway": "One exact charged chain is already complete on its declared scope: shared excitation dictionary -> ordered charged carrier -> exact centered readback -> quadratic readout theorem -> exact same-family triple. The remaining theorem burden is the affine descent through `C_hat_e^{cand}` and then `mu_phys(Y_e)`.",
        "logic": (
            "The lane starts from the ordered charged package, proves that the realized support is a one-dimensional "
            "linear subray, exposes the canonical quadratic support-extension direction, maps that into the charged "
            "excitation gaps, closes the two-scalar support-extension law shell, isolates the eta source-readback "
            "primitive as a weighted midpoint-defect invariant, and then uses the endpoint-ratio breaker for sigma. "
            "On the fixed ordered three-point family, the exact same-family quadratic theorem already closes the readout to one exact charged triple on that scope. "
            "The live scalar order beneath the public charged theorem rows is eta first and sigma second on the same carrier. "
            "At theorem level, though, the deeper exact waiting set starts with the latent candidate "
            "`C_hat_e^{cand}`. Promoting it is blocked by the branch-generator splitting theorem and its "
            "commutator clause, and the local corpus proves neither exact vanishing nor uniform quadratic smallness of that descended commutator yet. The strongest new route is only an extension candidate: a centered Schur-type "
            "`P->Q->P` transfer theorem would promote the proxy bridge if a refinement-uniform middle-factor bound can be certified, with the current local proxy gap staying safe for about `M < 119.56`. On the absolute side the equalizer route is no-go under common-shift "
            "symmetry, so the future slot is one affine-covariant absolute anchor `A_ch`. Conditional on theorem-grade centered promotion, the exact post-promotion burden descends further to the physical affine scalar `mu_phys(Y_e)`, from which the uncentered lift, determinant-line section, and affine anchor follow canonically."
        ),
        "tasks_text": "Open task: search for the full charged closure chain first by proving `oph_generation_bundle_branch_generator_splitting`, reduced to `compression_descendant_commutator_vanishes_or_is_uniformly_quadratic_small_after_central_split`, and then deriving the descended physical affine scalar `mu_phys(Y_e)`. If that still fails, keep the layered frontier explicit rather than falling back to older fit language.",
        "prediction_surface": "Charged theorem surface with exact centered readback, an exact same-family sidecar triple on the fixed ordered carrier, a common-shift no-go for the absolute scale, and a layered post-promotion route through `mu_phys(Y_e)` to the uncentered determinant-line lift and affine anchor.",
        "particles": ["electron", "muon", "tau"],
        "tasks": [
            "papers.compact.e.29-derive-the-yukawa-excitation-dictionary",
            "papers.compact.e.30-replace-koide-assisted-lepton-fitting-with-a-theorem",
        ],
    },
    {
        "key": "quarks",
        "title": "Quarks",
        "summary": "The quark lane contains a closed exact same-family readout on the selected D12 sheet, while the theorem-grade physical branch still waits on one intrinsic ray law or a genuine branch-changing theorem.",
        "takeaway": "One exact quark chain is already complete on its declared scope: shared excitation dictionary -> D12 mass-side continuation -> selected `sigma_ref` sheet -> quadratic readout theorem -> exact same-family sextet. The remaining theorem gap is the physical branch, not the exact same-family readout.",
        "logic": (
            "The local quark path takes the shared flavor data, emits the quark sector mean split, assembles the "
            "quark descent, builds the forward Yukawa matrices, and fixes the pure-B source-readback shell "
            "`[-beta, 0, +beta]` on `B_ord = [-1, 0, 1]`. The even ordered-family surface is fixed by the mean split "
            "and diagonal gap machinery. The active builder still waits on the emitted pure-B payload pair "
            "`source_readback_u_log_per_side_and_source_readback_d_log_per_side`, but the March 28, 2026 "
            "consolidation against the OPH tier ledger shows that a nonzero light-quark pure-B selector is not "
            "available at recovered-core tier. The broader repair frontier is therefore a D12 light-quark "
            "isospin-breaking selector / overlap-defect scalar. On the current same-family continuation branch the "
            "mass side sharpens further to the emitted same-family ray `D12_ud_mass_ray`, with `Delta_ud_overlap = ray_modulus / 5` and `eta_Q_centered = -((1-x2^2)/27) * ray_modulus`. "
            "On the fixed ordered three-point family, the exact same-family quadratic theorem already closes the readout to one exact quark sextet on that scope. "
            "The CKM/CP side is no longer the open burden: the forward Yukawa step already emits the same-label "
            "transport unitary `V_CKM^fwd = U_u^dagger U_d`, and its principal logarithm gives the gauge-fixed "
            "generator. But the current D12 sheet is now a strict no-go for the physical CKM shell: same-sheet "
            "rephasing cannot change the CKM invariants, and the emitted angles are too small. The selector-side frontier is no longer open on the local solver surface: the same-label left-handed orbit already closes to the singleton `sigma_ref`, so the selected branch is fixed but remains the current wrong-branch D12 sheet. The only finite local scan on disk is a same-sheet 4001-point `Delta_ud_overlap` scan against `reference_targets`, and that route is comparison-only rather than a valid physical repair. A smaller finite local basis orbit is already extractable from the current forward Yukawa surface, but its nontrivial elements use right-singular-basis substitutions and are therefore diagnostic-only rather than admissible physical selectors. The exact next quark object is therefore the one-scalar law `quark_d12_t1_value_law` on the already-emitted mass ray, with `intrinsic_scale_law_D12` retained as the derived wrapper, unless a genuinely new branch-changing theorem is found outside the current selected sheet."
        ),
        "tasks_text": "Open task: search for a real branch-changing theorem that escapes the selected wrong-branch `sigma_ref` sheet and reaches the physical CKM shell. If no such route exists on the current corpus, target only `quark_d12_t1_value_law` on the emitted `D12_ud_mass_ray`, with `intrinsic_scale_law_D12` retained as the derived wrapper, and keep diagnostic orbit improvements out of theorem state.",
        "prediction_surface": "Local forward quark Yukawa surface on the public table, plus an exact same-family sextet on the selected D12 sheet, with CKM/CP transport closed on the selected wrong-branch branch, the light-quark mass side reduced to the emitted one-parameter ray `D12_ud_mass_ray`, and the next exact theorem-side burden reduced to `quark_d12_t1_value_law`, with `intrinsic_scale_law_D12` retained as the derived wrapper.",
        "particles": ["up_quark", "down_quark", "strange_quark", "charm_quark", "bottom_quark"],
        "tasks": [
            "papers.compact.e.29-derive-the-yukawa-excitation-dictionary",
            "particles.microphysics.01-bridge-edge-statistics-to-yukawa-dictionary",
            "papers.compact.e.31-derive-quark-textures-ckm-mixing-and-cp-from-overlap-defects",
        ],
    },
    {
        "key": "neutrinos",
        "title": "Neutrinos",
        "summary": "The neutrino lane closes the weighted-cycle PMNS/hierarchy branch and carries an exact two-parameter adapter on the explicit positive selector segment, while the theorem-grade absolute spectrum still waits on `C_nu`.",
        "takeaway": "One exact neutrino sidecar is already complete on the positive selector segment: the emitted scale-free branch plus one selector coordinate and one rescaling hit the representative central splittings exactly. The remaining theorem object is the reduced correction invariant `C_nu`.",
        "logic": (
            "The lane derives m_star = v^2 / mu_u from the D10 core, builds the family-response tensor, the "
            "Majorana holonomy lift, the pullback metric, the forward Majorana matrix, and the splitting/ordering "
            "bundle. The old selector-law branch is S3-isotropic and is now explicitly ruled out for the physical "
            "atmospheric scale by the exact cap `max |Delta m^2| <= 8 a rho`. The repaired live branch instead uses the "
            "same-label scalar certificate together with the flavor cocycle invariants `gamma`, `eps`, and "
            "`gamma_half` to define the positive load segment between `chi = 1 + eps` and `1 + gamma_half`. On that "
            "one-dimensional affine segment, the balanced selector and the least-distortion selector coincide exactly at "
            "the midpoint, so `D_nu = (chi + 1 + gamma_half) / 2` and the repaired edge law becomes "
            "`w_e = q_e^(1 + gamma + eps / D_nu)`. That repaired branch "
            "lands in the physical PMNS window and the correct splitting hierarchy. The normalized same-label overlap-defect weight section `qbar_e` is already closed below the absolute attachment problem. "
            "No hidden discrete branch remains on that repaired lane; once the normalizer is fixed, the residual quotient is exactly the positive rescaling orbit. "
            "On that already-emitted positive selector segment, the two-parameter exact adapter moves one selector coordinate and one positive rescaling to hit the representative central splittings exactly. "
            "After factoring out the already-emitted proxy `P_nu = I_nu^0.5 * ratio_hat^0.5 * sum_defect^-1`, the exact remaining theorem object is the reduced correction invariant `C_nu`, with the paper-facing amplitude `B_nu` and then the final positive normalization scalar induced only after that correction law is fixed."
        ),
        "tasks_text": "Open task: search for the full closure chain by deriving the reduced bridge-correction invariant `C_nu` above the emitted proxy and then collapsing the remaining positive normalization orbit. If closure still fails, keep `C_nu` as the canonical smallest exact missing object and keep corridor fits compare-only.",
        "prediction_surface": "Weighted-cycle neutrino branch with PMNS/hierarchy closure, a closed normalized overlap-defect weight section, an exact two-parameter sidecar adapter on the positive selector segment, and one reduced correction invariant `C_nu` still separating the emitted proxy from a theorem-grade absolute spectrum.",
        "particles": ["electron_neutrino", "muon_neutrino", "tau_neutrino"],
        "tasks": [
            "papers.compact.e.32-derive-neutrino-masses-from-screen-capacity-as-a-theorem",
            "particles.neutrino.11-derive-majorana-phase-pullback-metric",
            "particles.e.32b-close-neutrino-splittings-and-mixing",
        ],
    },
    {
        "key": "hadrons",
        "title": "Hadrons",
        "summary": "The hadron lane is deferred and execution-contract-frozen: the theorem-side schema, receipt, and local writeback path are explicit, but promotable masses still require one production backend export bundle and then real production systematics.",
        "takeaway": "This lane is not waiting on one clever symbolic trick. It is waiting on actual nonperturbative computation, which is why it is deferred in the active closure program.",
        "logic": (
            "The hadron path steps down from D10 and the local quark masses into Lambda_MSbar^(3), seeds the "
            "unquenched ensemble family, realizes deterministic cfg/source payload identifiers, attaches a fixed "
            "RHMC/HMC schedule shell and conditional execution receipt, then builds the stable-channel "
            "sequence population/evaluation shells and aggregates them into the ground-state readout surface. Numerical "
            "hadron masses require one production backend export bundle with publication-complete manifest provenance, real `pi_iso`, `N_iso_direct`, and `N_iso_exchange` arrays, the executed runtime receipt `(N_therm, N_sep)`, and then declared production continuum/volume/chiral/statistical systematics. "
            "The surrogate execution bridge is only a diagnostic proof that the schema closes; it is not a promotable hadron output surface."
        ),
        "tasks_text": "Deferred task: if this lane is resumed later, land the production backend export bundle and then publish the production systematics for the seeded 2+1 QED-off stable channels. It is not part of the active exact-spectrum closure program.",
        "prediction_surface": "Execution-contract-frozen stable-channel hadron shell; public hadron rows remain hidden because this lane still requires a real backend bundle and real production computation.",
        "particles": ["proton", "neutron", "neutral_pion", "rho_770_0"],
        "tasks": [
            "papers.compact.e.33-close-the-nonperturbative-qcd-hadron-branch",
            "particles.e.33a-unquench-hadron-branch-and-publish-systematics",
        ],
    },
]


def load_results(path: pathlib.Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_task_tracker(path: pathlib.Path) -> Dict[str, Dict[str, Any]]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    tasks: Dict[str, Dict[str, Any]] = {}
    for stage in payload.get("stages", []):
        for task in stage.get("tasks", []):
            tasks[task["id"]] = task
    return tasks


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


def task_badge_color(status: str) -> str:
    if status == "completed":
        return COLORS["task_complete"]
    if status == "in_progress":
        return COLORS["task_in_progress"]
    if status == "deferred":
        return COLORS["task_deferred"]
    if status == "blocked":
        return COLORS["task_blocked"]
    return COLORS["task_pending"]


def task_card_height(task: Dict[str, Any], w: float) -> float:
    task_id_lines = wrap_identifier(task["id"], char_capacity(w - 36, 14, mono=True))
    return estimate_box_height(
        title=task["title"],
        body=[],
        w=w,
        title_size=18,
        body_size=14,
        body_family=MONO_FAMILY,
        prewrapped_body_lines=["Task ID:"] + task_id_lines,
        accent=True,
    )


def draw_task_card(task: Dict[str, Any], x: float, y: float, w: float) -> Tuple[str, float]:
    task_id_lines = wrap_identifier(task["id"], char_capacity(w - 36, 14, mono=True))
    height = task_card_height(task, w)
    return (
        draw_box(
            x=x,
            y=y,
            w=w,
            h=height,
            title=task["title"],
            body=[],
            fill=COLORS["task_fill"],
            stroke=COLORS["task_stroke"],
            title_size=18,
            body_size=14,
            accent=task_badge_color(task["status"]),
            badge=task["status"].replace("_", " "),
            body_fill="#ffd2bf",
            body_family=MONO_FAMILY,
            prewrapped_body_lines=["Task ID:"] + task_id_lines,
        ),
        height,
    )


def _format_exact_output(entry: Dict[str, Any]) -> str:
    if "mass_gev" in entry:
        return f"{entry['mass_gev']} GeV"
    if "mass_eV" in entry:
        return f"{entry['mass_eV']} eV"
    return "n/a"


def particle_body(row: Dict[str, Any], exact_entry: Dict[str, Any] | None = None) -> List[str]:
    info = PARTICLE_INFO[row["particle_id"]]
    lines = [
        f"Symbol: {info['symbol']}. What it is: {info['plain']}",
        f"Why this row matters: {GROUP_ROW_TEXT[row['group']]}",
    ]
    if exact_entry is not None:
        lines.extend(
            [
                f"Surface class: exact {exact_entry['exact_kind']} on {exact_entry['scope']}.",
                f"OPH exact output: {_format_exact_output(exact_entry)}.",
                f"Reference: {row['reference_display']}.",
                f"Caveat: {exact_entry['note']}",
            ]
        )
    else:
        status = row["status"]
        lines.extend(
            [
                f"Surface class: {STATUS_TEXT[status]} meaning {STATUS_EXPLAINER[status]}.",
                f"OPH output: {row['prediction_display_gev']} GeV.",
                f"Reference: {row['reference_display']}.",
            ]
        )
        if row["delta_display"] != "n/a":
            lines.append(f"Gap vs reference: {row['delta_display']}.")
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
    badge_text = "exact sidecar" if exact_entry is not None and not exact_entry.get("promotable", False) else STATUS_TEXT[row["status"]]
    badge_color = COLORS["task_complete"] if exact_entry is not None else STATUS_BAR[row["status"]]
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
            body=["This lane stays off the public particle table until a closure-grade predictor exists."],
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
            body=["This lane stays off the public particle table until a closure-grade predictor exists."],
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
                body=["This lane stays off the public particle table until a closure-grade predictor exists."],
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
    return "simulation_dependent"


def lane_panel_height(
    lane: Dict[str, Any],
    rows_by_id: Dict[str, Dict[str, Any]],
    exact_by_id: Dict[str, Dict[str, Any]],
    tasks: Dict[str, Dict[str, Any]],
    w: float,
) -> float:
    inner_w = w - 36
    summary_lines = wrap_text(lane["summary"], char_capacity(inner_w - 18, 16))
    takeaway_lines = wrap_text(f"In one sentence: {lane['takeaway']}", char_capacity(inner_w - 18, 15))
    header_h = 118 + len(summary_lines) * 20 + len(takeaway_lines) * 18
    label_h = 28
    label_gap = 10
    section_gap = 20
    task_gap = 28

    theorem_h = estimate_box_height(
        title=lane["summary"],
        body=[lane["logic"]],
        w=inner_w,
        title_size=20,
        body_size=15,
    )
    if lane["tasks"]:
        task_total_h = sum(task_card_height(tasks[task_id], inner_w) for task_id in lane["tasks"]) + task_gap * (len(lane["tasks"]) - 1)
    else:
        task_total_h = estimate_box_height(
            title="No open blocker in ledger",
            body=[lane["tasks_text"]],
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
    tasks: Dict[str, Dict[str, Any]],
    x: float,
    y: float,
    w: float,
) -> Tuple[str, float]:
    inner_x = x + 18
    inner_w = w - 36
    status = lane_status(lane, rows_by_id)
    status_color = STATUS_BAR[status]
    summary_lines = wrap_text(lane["summary"], char_capacity(inner_w - 18, 16))
    takeaway_lines = wrap_text(f"In one sentence: {lane['takeaway']}", char_capacity(inner_w - 18, 15))
    header_h = 118 + len(summary_lines) * 20 + len(takeaway_lines) * 18
    label_h = 28
    label_gap = 10
    section_gap = 20
    task_gap = 28

    theorem_h = estimate_box_height(
        title=lane["summary"],
        body=[lane["logic"]],
        w=inner_w,
        title_size=20,
        body_size=15,
    )
    if lane["tasks"]:
        task_heights = [task_card_height(tasks[task_id], inner_w) for task_id in lane["tasks"]]
        task_total_h = sum(task_heights) + task_gap * (len(task_heights) - 1)
    else:
        task_heights = []
        task_total_h = estimate_box_height(
            title="No open blocker in ledger",
            body=[lane["tasks_text"]],
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
            [f"{len(visible_particles)} public outputs shown in this lane"],
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
            "Implemented theorem / technique",
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
            "Missing proof / closure tasks",
            fill="#241611",
            stroke=COLORS["task_stroke"],
            text_fill=COLORS["ink"],
        )
    )
    cursor += label_h + label_gap

    if lane["tasks"]:
        first_task_y = cursor
        parts.append(draw_vertical_arrow(x + w / 2.0, theorem_y + theorem_h + 4, first_task_y - 6, color=COLORS["task_stroke"], dashed=True))
        previous_bottom = None
        for index, task_id in enumerate(lane["tasks"]):
            task_y = cursor
            task_markup, task_h = draw_task_card(tasks[task_id], inner_x, task_y, inner_w)
            parts.append(task_markup)
            if previous_bottom is not None:
                parts.append(draw_vertical_arrow(x + w / 2.0, previous_bottom + 4, task_y - 6, color=COLORS["task_stroke"], dashed=True))
            previous_bottom = task_y + task_h
            cursor += task_h + task_gap
        cursor -= task_gap
        task_bottom = previous_bottom or first_task_y
    else:
        no_task_h = estimate_box_height(
            title="No open blocker in ledger",
            body=[lane["tasks_text"]],
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
                title="No open blocker in ledger",
                body=[lane["tasks_text"]],
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
            "Prediction surface",
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
            title="Prediction surface",
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
            "Tracked particle outputs",
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


def build_svg(results: Dict[str, Any], tasks: Dict[str, Dict[str, Any]], exact_by_id: Dict[str, Dict[str, Any]]) -> str:
    rows_by_id = {row["particle_id"]: row for row in results["rows"]}
    generated_utc = results["generated_utc"]
    built_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    counts = Counter(row["status"] for row in results["rows"])
    closedish = counts["structural"] + counts["calibration"] + counts["secondary_quantitative"]
    total_rows = len(results["rows"])

    parts: List[str] = []
    current_y = 48.0

    # Header sizes
    intro_lines = [
        "A top-down map of the current OPH particle program: start at the common inputs, move through the sector lanes, and finish at the particle rows the code can currently print.",
        "Blue and green regions show implemented derivation surfaces. Orange cards mark the exact theorem or computation burden that still blocks closure.",
        "This poster is designed as both a research snapshot and a reader-facing explainer, so each lane includes plain-language takeaways as well as the sharper technical status.",
    ]
    intro_height = 3 * 24

    legend_items = [
        ("Implemented theorem / technique", COLORS["logic_stroke"]),
        ("Open proof / closure task", COLORS["task_stroke"]),
        ("Particle end node", COLORS["prediction_stroke"]),
        ("Connector / dependency flow", COLORS["line"]),
    ]
    legend_w = 420.0
    legend_h = 152.0
    legend_x = WIDTH - MARGIN_X - legend_w
    legend_y = current_y

    header_bottom = max(current_y + 112 + intro_height, legend_y + legend_h)
    current_y = header_bottom + 30

    # Inputs section
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
            "title": "Declared External Input: P",
            "body": [f"P = {inputs['P']}. This scalar is the main numerical upstream input for the electroweak, flavor, and hadron branches."],
            "fill": COLORS["input_fill"],
            "stroke": COLORS["input_stroke"],
        },
        {
            "x": MARGIN_X + axiom_w + input_gap + p_w + input_gap,
            "w": other_w,
            "title": "Other Declared Input Surface",
            "body": [f"log_dim_H = {inputs['log_dim_H']} feeds the neutrino estimate lane. loops = {inputs['loops']} and hadron_profile = {inputs['hadron_profile']} are the settings used by this report."],
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
        "Start with the OPH axioms plus the declared live scalar P and the extra input surface used by the neutrino estimate lane.",
        "Then read each lane from top to bottom: what the lane already does, what is still missing, what prediction surface it produces, and which particle rows are currently publishable.",
        f"The badge reports {closedish} of {total_rows} tracked rows above continuation / simulation status. In plain terms: those are the rows that are already beyond the merely exploratory stage.",
        "The broader UV/BW premise boundary still sits above the particle lanes. Three cap-pair extraction witnesses are already explicit on disk; the remaining local gap is now sharper: on the local-Gibbs pullback branch the Markov/recovery side is already latent, so the live missing clause inside the fixed-local-collar Markov/faithfulness datum is the eventual lower spectral bound feeding the faithful modular-defect term, after which ordered null cut-pair rigidity remains.",
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

    row1_heights = [lane_panel_height(lane, rows_by_id, exact_by_id, tasks, PANEL_W) for lane in row1]
    row1_y = current_y
    row1_h = max(row1_heights)

    row2_heights = [lane_panel_height(lane, rows_by_id, exact_by_id, tasks, PANEL_W) for lane in row2]
    row2_y = row1_y + row1_h + PANEL_GAP_Y
    row2_h = max(row2_heights)

    row3_h = lane_panel_height(row3[0], rows_by_id, exact_by_id, tasks, HADRON_W)
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
                "structural = massless or exact structural rows",
                "calibration = produced by the implemented P-driven electroweak chain",
                "secondary = quantitative branch built on an earlier calibrated layer",
                "continuation = useful live number, but the theorem chain is still open",
                "simulation = execution-bound lane; real computation still dominates",
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
                "This chart maps the implemented derivation pipeline, its closure tasks, and the emitted public rows. It is a status poster, not a claim that the full particle zoo already sits on one finished theorem chain from P plus axioms.",
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

    parts.append(render_wrapped_text(MARGIN_X, input_label_y + 20, ["1. Inputs"], font_size=24, fill=COLORS["ink"], weight=700, line_height=26))
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
        panel_markup, _ = draw_lane_panel(lane, rows_by_id, exact_by_id, tasks, x, row1_y, PANEL_W)
        parts.append(panel_markup)
    for lane, x in zip(row2, row2_xs):
        panel_markup, _ = draw_lane_panel(lane, rows_by_id, exact_by_id, tasks, x, row2_y, PANEL_W)
        parts.append(panel_markup)
    hadron_markup, _ = draw_lane_panel(row3[0], rows_by_id, exact_by_id, tasks, hadron_x, row3_y, HADRON_W)
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
    parser.add_argument("--task-tracker-yaml", default=str(TASK_TRACKER_YAML), help="Path to the task tracker YAML.")
    parser.add_argument("--exact-nonhadron-json", default=str(EXACT_NONHADRON_JSON), help="Path to the exact non-hadron bundle JSON.")
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT), help="Output SVG path.")
    args = parser.parse_args()

    results = load_results(pathlib.Path(args.results_json))
    tasks = load_task_tracker(pathlib.Path(args.task_tracker_yaml))
    exact_by_id = load_exact_nonhadron_bundle(pathlib.Path(args.exact_nonhadron_json))
    svg = build_svg(results=results, tasks=tasks, exact_by_id=exact_by_id)

    output = pathlib.Path(args.output)
    output.write_text(svg, encoding="utf-8")
    print(f"saved: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
