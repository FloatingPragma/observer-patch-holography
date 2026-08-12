#!/usr/bin/env python3
"""Build the public OPH open-problem ledger from live GitHub issues."""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPO = "FloatingPragma/observer-patch-holography"
DEFAULT_JSON_OUT = ROOT / "tracking" / "open_issues" / "open_problem_ledger.json"
DEFAULT_MD_OUT = ROOT / "OPEN_PROBLEMS.md"


V2_TITLE_RE = re.compile(r"^\[([A-H]\d+)\]\s+")
V2_FIELD_NAMES = (
    "Deliverables",
    "Boundary",
    "Depends on",
    "Wave",
    "Suggested dependency",
    "Absorbs",
    "Plan",
)
V2_FIELD_MARKER = "|".join(re.escape(name) for name in V2_FIELD_NAMES)
# Contract fields may begin a paragraph or follow a sentence terminator when
# several compact fields share one paragraph.  Requiring that boundary keeps
# prose headings such as ``Audit boundary:`` and ``Remaining deliverables:``
# from being parsed as the exact ``Boundary:`` or ``Deliverables:`` fields.
V2_FIELD_PREFIX = r"(?:\A|(?<=[.;!?])[ \t]+|(?<=[\r\n])[ \t]*)"
V2_TRACKS: dict[str, dict[str, Any]] = {
    "track:foundations": {"code": "A", "slug": "foundations", "wave": 1},
    "track:observer-laws": {"code": "B", "slug": "observer-laws", "wave": 1},
    "track:geometry": {"code": "C", "slug": "geometry", "wave": 2},
    "track:time": {"code": "D", "slug": "time", "wave": 2},
    "track:covariant-net": {"code": "E", "slug": "covariant-net", "wave": 3},
    "track:gravity": {"code": "F", "slug": "gravity", "wave": 3},
    "track:custody": {"code": "G", "slug": "custody", "wave": "standing"},
    "track:quantitative": {"code": "H", "slug": "quantitative", "wave": 4},
}
V2_SIZE_LABELS = {"size:S", "size:M", "size:L"}
V2_STATE_LABELS = {"optional", "parked", "standing", "blocked"}
V2_WAVE_TITLES = {
    1: "Observer laws",
    2: "Geometry and time",
    3: "Covariant net and gravity",
    4: "Quantitative outputs",
}

V3_TITLE_RE = re.compile(r"^\[V3-([A-Z]{1,2}[0-9]+)\]\s+")
POLICY_FIELDS = {
    "phase",
    "claim_level",
    "blocker",
    "closure",
    "falsification",
    "chrome_policy",
}
PLACEHOLDER_POLICY = {
    "phase": "unclassified",
    "claim_level": "open",
    "blocker": "Classify blocker from the live issue body.",
    "closure": "Add exact closure criterion to this ledger.",
    "falsification": "Add exact falsification criterion to this ledger.",
    "chrome_policy": (
        "Do not launch workers until the issue has a concrete local packet."
    ),
}


ISSUE_POLICY: dict[int, dict[str, str]] = {
    28: {
        "phase": "compact-structure",
        "claim_level": "theorem gap",
        "blocker": "Hypercharge lattice and exact Z6 quotient proof packet.",
        "closure": "Standard Model gauge paper and DAG surfaces cite the same theorem-grade lattice/quotient derivation.",
        "falsification": "A required SM charge assignment cannot be represented on the proposed OPH lattice.",
        "chrome_policy": "Use only for proof audit after a local lattice packet exists.",
    },
    32: {
        "phase": "particle-rg",
        "claim_level": "exact finite representation-index frontier; physical RG packet open",
        "blocker": (
            "The finite matter packet fixes per-copy representation indices and "
            "the imported one-loop QFT functional gives a parametric beta law. "
            "Physical family and scalar multiplicities, the complete field census, "
            "ordered EFT intervals, thresholds, decoupling, scheme maps, Jacobians, "
            "term masks, deterministic remainders, and interval effects remain open."
        ),
        "closure": (
            "One target-clean source packet emits the complete same-branch RG, "
            "threshold, finite matching, scheme, Jacobian, mask, remainder, and "
            "interval data with independent replay and explicit imported-QFT provenance."
        ),
        "falsification": (
            "A hidden threshold, scalar or family choice, scheme freedom, target-fitted "
            "remainder, or unbounded truncation is needed to recover a quantitative row."
        ),
        "chrome_policy": "Use only to audit a concrete same-branch RG and matching packet.",
    },
    49: {
        "phase": "observer-monograph",
        "claim_level": "claim-hygiene gap",
        "blocker": "Proton stability and proton-spin claims need theorem-grade status or downgrade wording.",
        "closure": "Observer paper and summaries agree on exact theorem, conditional result, or open-continuation status.",
        "falsification": "A retained proton claim depends on unsupported QCD or baryon-structure assumptions.",
        "chrome_policy": "Use for independent claim audit if local wording remains ambiguous.",
    },
    55: {
        "phase": "observer-monograph",
        "claim_level": "decision gate",
        "blocker": "Critical-superstring lift needs either a theorem packet or an explicit non-claim decision.",
        "closure": "Paper and book consistently mark the lift as proved, conditional, or not part of the OPH core.",
        "falsification": "The lift is required for a core OPH claim but cannot be derived or cleanly removed.",
        "chrome_policy": "Use for proof/claim audit after local classification.",
    },
    58: {
        "phase": "observer-monograph",
        "claim_level": "theorem package gap",
        "blocker": "Observer continuation/backup theorem package.",
        "closure": "Continuation, recovery, and backup statements are packaged with hypotheses and proof boundary.",
        "falsification": "Observer continuation requires nonlocal state access outside OPH overlap rules.",
        "chrome_policy": "Use for proof audit after local theorem packet.",
    },
    59: {
        "phase": "observer-monograph",
        "claim_level": "audit gap",
        "blocker": "Line-by-line audit of the Additional Problem Closures table.",
        "closure": "Every table row links to a theorem, conditional claim, or open issue.",
        "falsification": "A table row asserts closure without a traceable theorem or artifact.",
        "chrome_policy": "Use for independent table audit after local row mapping.",
    },
    60: {
        "phase": "observer-monograph",
        "claim_level": "final audit gate",
        "blocker": "Final proof, citation, and reproducibility audit.",
        "closure": "The monograph passes local build, citation, theorem-status, and reproducibility checks.",
        "falsification": "A core claim cannot be linked to a proof or reproducible artifact.",
        "chrome_policy": "Use for final independent audit after local checks are green.",
    },
    66: {
        "phase": "reality-paper",
        "claim_level": "theorem gap",
        "blocker": "Extension from abelian cycle holonomy to the full OPH defect hierarchy.",
        "closure": "Reality paper states the full hierarchy construction with hypotheses and proof boundary.",
        "falsification": "Nonabelian or higher-defect holonomy is incompatible with the proposed reconciliation law.",
        "chrome_policy": "Use for proof search only after local hierarchy packet.",
    },
    70: {
        "phase": "reality-paper",
        "claim_level": "theorem gap",
        "blocker": "Proof that coarse-graining commutes with reconciliation.",
        "closure": "Coarse-graining/reconciliation commutation theorem is stated and synchronized across reality surfaces.",
        "falsification": "A counterexample shows reconciliation depends on refinement order.",
        "chrome_policy": "Use for proof audit after local counterexample search.",
    },
    112: {
        "phase": "archive-speculative",
        "claim_level": "personal/open hypothesis",
        "blocker": "No public theorem package or falsifiable closure artifact is defined.",
        "closure": "Either define a falsifiable OPH theorem target or archive as non-core speculation.",
        "falsification": "The hypothesis cannot be connected to OPH overlap consistency without extra metaphysics.",
        "chrome_policy": "Do not spend Chrome workers until a precise theorem target exists.",
    },
    113: {
        "phase": "observer-monograph",
        "claim_level": "theorem construction gap",
        "blocker": "Closure map and invariant sector construction.",
        "closure": "Observer paper gives a reproducible closure-map/invariant-sector theorem package.",
        "falsification": "The closure map is not well-defined or is not invariant under allowed patch refinements.",
        "chrome_policy": "Use for proof audit after local construction.",
    },
    153: {
        "phase": "hardware-gated-hadrons",
        "claim_level": "out of local scope",
        "blocker": "Working OPH hadron backend on suitable hardware such as GLORB/Echosahedron.",
        "closure": "Production backend output and continuum/volume/chiral/statistical systematics are published.",
        "falsification": "Surrogate local hadron artifacts are required as if they were production QCD outputs.",
        "chrome_policy": "Do not use Chrome workers for backend execution.",
    },
    155: {
        "phase": "continuation",
        "claim_level": "open branch",
        "blocker": "Strong-CP branch theorem or explicit continuation boundary.",
        "closure": "Strong-CP mechanism is derived, falsified, or downgraded consistently across public surfaces.",
        "falsification": "The proposed OPH branch leaves theta_QCD unconstrained while claiming closure.",
        "chrome_policy": "Use only after a concrete branch packet exists.",
    },
    157: {
        "phase": "hardware-gated-hadrons",
        "claim_level": "out of local scope",
        "blocker": "Nonperturbative QCD/hadron backend and systematics.",
        "closure": "The particle and gauge papers state the backend output and budgets or mark hadrons out-of-scope.",
        "falsification": "A specialist paper promotes hadron masses without production backend evidence.",
        "chrome_policy": "Do not use Chrome workers for backend execution.",
    },
    199: {
        "phase": "quark-global-classification",
        "claim_level": "selected-fiber descent closed; source spread non-identifiable; global classification separate",
        "blocker": "Class-uniform public quark-frame descent and sigma classification.",
        "closure": "Global quark frame classes are classified or every claim remains explicitly selected-fiber only.",
        "falsification": "Another admissible public quark class breaks the selected-fiber descent statement.",
        "chrome_policy": "Use for proof audit after local classification packet.",
    },
    201: {
        "phase": "charged-lepton-source",
        "claim_level": "source theorem gap",
        "blocker": "Sector-isolated charged determinant trace-lift attachment and normalization.",
        "closure": "The determinant character lands on the physical charged determinant line without target readback.",
        "falsification": "The additive determinant normalization remains underdetermined on the OPH source data.",
        "chrome_policy": "Use for proof search after local trace-lift packet.",
    },
    207: {
        "phase": "top-codomain-bridge",
        "claim_level": "constructive conversion contract",
        "blocker": "Source-side extraction response kernel from Q007TP4 coordinate to Q007TP codomain.",
        "closure": "Converted top coordinate and uncertainty budget close without using Q007TP as calibration input.",
        "falsification": "The direct-top codomain cannot be mapped from the current theorem coordinate without a free shift.",
        "chrome_policy": "Use only to audit a proposed response kernel.",
    },
    212: {
        "phase": "quark-off-canonical",
        "claim_level": "candidate-only",
        "blocker": "Target-free off-canonical flavor transport law and P-to-sigma evaluator.",
        "closure": "Off-canonical quark masses move by theorem-grade P-driven transport, not default-universe anchoring.",
        "falsification": "Off-canonical motion requires anchoring to the canonical PDG target surface.",
        "chrome_policy": "Use for proof audit after local off-canonical transport packet.",
    },
    223: {
        "phase": "thomson-endpoint",
        "claim_level": "constructive contract",
        "blocker": "Ward-projected source endpoint including hadronic spectral measure and interval certificate.",
        "closure": "alpha_Th(P) is emitted from source objects with certified transport/error bounds.",
        "falsification": "Measured alpha(0) or a free screened ansatz is required to close the endpoint.",
        "chrome_policy": "Use only after the source endpoint packet exists.",
    },
    224: {
        "phase": "p-root-adoption",
        "claim_level": "blocked on certified root",
        "blocker": "Certified P root after #223 and #32 close.",
        "closure": "Live particle consumers read the certified trunk root and legacy P paths are non-default.",
        "falsification": "Derived P cannot replace legacy/candidate paths without weakening prediction rows.",
        "chrome_policy": "Do not use until endpoint and RG interval gates close.",
    },
    225: {
        "phase": "publication-sync",
        "claim_level": "blocked on #223/#224",
        "blocker": "Certified derived P closure values and live consumer adoption.",
        "closure": "Observers and particle papers quote the same certified P, alpha, and dependent quantities.",
        "falsification": "Paper values must be published before the code root is certified.",
        "chrome_policy": "Use for claim-hygiene audit after code root certification.",
    },
    231: {
        "phase": "axiom-status",
        "claim_level": "independence decision gate",
        "blocker": "Interface independence proof or downgrade against the three-axiom basis.",
        "closure": "Basis language matches the proved independence/dependence status of each classified interface.",
        "falsification": "Interface independence is asserted while a dependency proof or countermodel is missing.",
        "chrome_policy": "Use for proof/countermodel audit after local packet.",
    },
    232: {
        "phase": "bw-cap-pair",
        "claim_level": "theorem gap",
        "blocker": "Transported BW/geometric cap-pair extraction and ordered cut-pair rigidity.",
        "closure": "BW cap-pair extraction and ordered cut-pair rigidity are emitted as theorem-grade artifacts.",
        "falsification": "The extracted cap/cut structure depends on arbitrary collar choices.",
        "chrome_policy": "Use for proof audit after local BW packet.",
    },
    233: {
        "phase": "retired-claim-audit",
        "claim_level": "withdrawn premise",
        "blocker": "No theorem gap remains: economy selection is not part of the three-axiom basis.",
        "closure": "All active surfaces classify finite economy scores as diagnostics that select no physics.",
        "falsification": "An active theorem continues to infer physical uniqueness from an economy score.",
        "chrome_policy": "Do not launch a worker; enforce the source guard.",
    },
    234: {
        "phase": "particle-provenance",
        "claim_level": "audit gap",
        "blocker": "Blind-prediction provenance and convention-sensitivity ledger.",
        "closure": "Every quantitative row records input use, blind/compare status, provenance, and convention sensitivity.",
        "falsification": "A public quantitative row hides target leakage or unquantified convention dependence.",
        "chrome_policy": "Use for audit only after local provenance artifact exists.",
    },
    235: {
        "phase": "p-closure-theorem",
        "claim_level": "residual theorem gap",
        "blocker": "Root monotonicity/uniqueness and exact endpoint boundary.",
        "closure": "P-closure has a scalar root theorem and endpoint status boundary synchronized across surfaces.",
        "falsification": "The claimed root is nonunique or endpoint-dependent in a way the papers do not state.",
        "chrome_policy": "Use for proof audit after local monotonicity packet.",
    },
    236: {
        "phase": "publication-ledger",
        "claim_level": "ledger artifact",
        "blocker": "Keep this ledger synchronized across README, papers, book, and public summaries.",
        "closure": "Public OPH surfaces point to the same open-problem ledger and do not imply hidden closure.",
        "falsification": "A public surface claims a branch is closed while this ledger or its issue says open.",
        "chrome_policy": "Not needed unless a downstream claim-hygiene audit is ambiguous.",
    },
    237: {
        "phase": "screen-microphysics",
        "claim_level": "benchmark suite gap",
        "blocker": "Explicit existence-program benchmark suite for screen microphysics reference architecture.",
        "closure": "Benchmarks specify model class, expected artifacts, pass/fail criteria, and reproducibility path.",
        "falsification": "The reference architecture cannot be exercised by any finite benchmark suite.",
        "chrome_policy": "Use for architecture audit after local benchmark suite skeleton.",
    },
    311: {
        "phase": "particle-ontology",
        "claim_level": (
            "bounded classical spectral match; physical particle criterion open"
        ),
        "blocker": (
            "On the exact twelve-vertex support, K_k = 5 I - A_k gives an "
            "explicit classical vector-spring realization in the positive "
            "hex-lattice metric, with an exact edge-Hessian proof, from which "
            "the twisted adjacency spectral family is recoverable. A separate "
            "issue-634 local-domain receipt matches its own six-sector and "
            "scalar spectra. No identity bridge joins those finite domains. "
            "Neither result covers a complete interface or extended source "
            "domain, and no measurement semantics, asymptotic dynamics, pole, "
            "or cofinal refinement family is emitted."
        ),
        "closure": (
            "An extended source-defined quantum measurement, state-update, dynamics, "
            "or equivalent physical spectral criterion excludes the relevant "
            "same-domain vector-spring completion and controls composition, "
            "asymptotic states, and refinement."
        ),
        "falsification": (
            "The proposed particle discriminator remains a function only of "
            "the exact-support spectral family or of the separate matched "
            "local-domain sector/scalar projection, without an extended "
            "measurement or state-update interface."
        ),
        "chrome_policy": "Use only to audit a concrete dynamics, pole, or quantization packet.",
    },
    334: {
        "phase": "gravity-scale-composition",
        "claim_level": "source clock plus gravity readout composition gap",
        "blocker": (
            "Issue #633 proves "
            "PHYSICAL_UNITS_NOT_EVALUABLE_ON_DECLARED_SERIALIZED_INTERFACE. This "
            "bounded verdict is not a complete-domain theorem. The extended matter "
            "and scalar/Yukawa domains in #569 and #630 must select a physical "
            "transition. A source-derived gravity length or radius, Einstein "
            "attachment, common provenance graph, and outward interval composition "
            "into G_SI remain open."
        ),
        "closure": (
            "A target-clean physical transition from the extended source domain and "
            "a source gravity readout compose on one branch into an independently "
            "replayed G_SI interval."
        ),
        "falsification": (
            "The composition imports measured G, Planck units, a cosmological target, "
            "or incompatible clock and gravity provenance."
        ),
        "chrome_policy": "Use only after #569 and #630 supply a physical transition and a concrete gravity-scale packet exists.",
    },
    547: {
        "phase": "n-dependent-electroweak-balance",
        "claim_level": "physical N-balance law open",
        "blocker": (
            "The local carrier belongs to #631 and direct record capacity belongs "
            "to #505. Their source-branch composition does not yet derive the "
            "N-dependent electroweak balance law."
        ),
        "closure": (
            "Completed #505 and #631 packets are bound on one target-clean branch "
            "and derive the balance equation with refinement, ancestry, stability, "
            "and deletion controls."
        ),
        "falsification": (
            "A nonzero balance residual is absorbed by redefining the carrier or "
            "capacity, or the equation is imposed rather than derived."
        ),
        "chrome_policy": "Use only after #505 and #631 emit their positive packets.",
    },
    505: {
        "phase": "direct-n-capacity-closure",
        "claim_level": "direct source-only record-capacity integration target",
        "blocker": (
            "The fixed-D=24 source checkpoint packet is complete, but #551 has "
            "not emitted one capacity-indexed source family with a regulator-stable "
            "unique slack zero or a complete stated-domain non-identifiability theorem."
        ),
        "closure": (
            "Consume the #551 verdict and certify robust direct closure "
            "F_set,r,0(D_star)={D_star}, with capacity extension, refinement, "
            "whole-fiber, ambiguity, and target-exclusion controls kept distinct."
        ),
        "falsification": (
            "The constructor reads an expected answer, a second admissible closure "
            "survives, capacity extension creates confusability, or the claimed "
            "cosmic value is imported from cosmology or the electroweak bridge."
        ),
        "chrome_policy": (
            "Use only to audit the composed #551 packet and fixed-D primitives; "
            "do not substitute target-guided equation search."
        ),
    },
    551: {
        "phase": "direct-n-capacity-closure",
        "claim_level": "capacity-indexed selector theorem-or-no-go",
        "blocker": (
            "The source corpus supplies a fixed-D checkpoint construction, not one "
            "target-clean rule D -> packet_r(D) -> M_0,r(D) with the exact "
            "subleading slack and regulator control."
        ),
        "closure": (
            "Prove existence, nontriviality, uniqueness, and cutoff independence "
            "of one slack zero, or give a complete same-antecedent counterfamily "
            "proving non-identifiability over the declared source class."
        ),
        "falsification": (
            "A second zero, regulator-dependent zero, separately authored per-D "
            "fixture, target input, or admissible same-source counterfamily survives "
            "a claimed positive selector theorem."
        ),
        "chrome_policy": (
            "Use for independent theorem or countermodel audit after a local "
            "capacity-indexed packet exists; local exact construction comes first."
        ),
    },
    522: {
        "phase": "cosmology-source-clock",
        "claim_level": "source generator and operational clock gap",
        "blocker": (
            "Exact collar counts do not select reserve weights, P/24, a weighted half, "
            "a generator coefficient, or an event-to-log-thickness clock scale."
        ),
        "closure": (
            "One target-free source DAG emits the reserve weights, full-collar derivative, "
            "weighted orientation identity, and operational clock map consumed by the simulator."
        ),
        "falsification": (
            "A second source-compatible weighting, generator, or clock scale survives every "
            "declared finite record and acceptance check."
        ),
        "chrome_policy": "Use only after a source generator and clock packet exists locally.",
    },
    623: {
        "phase": "physical-scalar-attachment",
        "claim_level": "physical scalar pole, vacuum, and multiplicity-discrimination gap",
        "blocker": (
            "The selected charged-response and pole-residue simulator artifacts "
            "expose no explicit scalar input, but that two-artifact audit is not "
            "an exhaustive producer inventory or a theorem on physical scalar "
            "completions. The #616 empty, duplicate, and inert rows are "
            "countermodels only to the enumerated grammar-visible checks."
        ),
        "closure": (
            "Land a source-derived scalar carrier with a pole or equivalent "
            "spectral receipt, potential/vacuum attachment, and discrimination "
            "against the retained grammar-scope countermodels; alternatively "
            "prove an exhaustive completion-indexed constant-observation theorem "
            "over a promotion-aware producer inventory."
        ),
        "falsification": (
            "The proposed negative theorem omits a scalar-sensitive pole, gauge "
            "response, potential/vacuum, multiplicity-rank, or Yukawa producer, "
            "or a proposed physical attachment imports a Higgs target."
        ),
        "chrome_policy": "Use only after a scalar-sensitive source packet or exhaustive producer inventory exists.",
    },
    624: {
        "phase": "noncentral-seam-classification",
        "claim_level": "full-schema classification gap with exact finite subresults",
        "blocker": (
            "Z7 has no faithful homomorphism into the order-two register lane or "
            "the existing order-six seam-class lane, but these two lanes are not "
            "proved exhaustive over all complete A1-A3 coefficient constructions. "
            "The noncontractible Z2->1 witness supplies a sector only through a "
            "separately supplied matter character."
        ),
        "closure": (
            "Prove an exhaustive classification theorem reducing every complete "
            "A1-A3 coefficient construction to the pinned lanes, or exhibit and "
            "classify every additional lane with its operational, refinement, "
            "meaning, optimizer, sector, and matter-character data."
        ),
        "falsification": (
            "The proposed product register splits a primitive central port, an A2 "
            "or refinement diagram fails, the A3 product optimizer is not proved, "
            "or the higher witness has trivial kernel and cokernel."
        ),
        "chrome_policy": "Use only to audit a complete full-schema lift and nontrivial 2-type.",
    },
    625: {
        "phase": "load-readback-classification",
        "claim_level": "full-schema independence or forcing gap",
        "blocker": (
            "The covariant half-count family and alternative quadratic form pass "
            "reduced carrier checks, but no complete A1 operational/refinement, "
            "A2 naturality, and A3 optimizer lift is proved."
        ),
        "closure": (
            "Either derive the integer counting normalization and discrete cost, "
            "or construct a complete A1-A3 alternative-readback model preserving "
            "every operational, refinement, meaning, and optimizer clause."
        ),
        "falsification": (
            "A purported alternative changes the feasible family or optimizer, "
            "breaks accepted-data naturality, or fails a complete refinement map."
        ),
        "chrome_policy": "Use only to audit a complete full-schema load/readback model.",
    },
    626: {
        "phase": "boolean-scatter-settling",
        "claim_level": "revised constructive theorem after exact negative baseline",
        "blocker": (
            "The canonical synchronous scatter-project rule has a satisfiable depth-two "
            "period-two orbit because downstream projections can rewrite parent wires."
        ),
        "closure": (
            "A source-pinned directional rule or isolation gadget passes the #328 regression, "
            "a universal compiler intertwiner, and explicit resource and settling bounds."
        ),
        "falsification": (
            "A satisfiable compiled circuit retains a nonzero cycle or the proposed "
            "intertwiner fails on a reachable production state."
        ),
        "chrome_policy": "Use for proof audit only after the revised production implementation exists.",
    },
    627: {
        "phase": "physical-routed-seam-selection",
        "claim_level": "physical grammar and matter-action selection gap",
        "blocker": (
            "No source packet selects the physical seam mechanism, character, or "
            "2-representation. The pure-hypercharge character menus have fixed "
            "dimensions 15/7, 15/3, and 15/7/3/1, while the canonical diagonal "
            "global-form kernel fixes all fifteen realized local matter states. "
            "These actions must not be conflated."
        ),
        "closure": (
            "A target-free source packet selects or classifies the seam grammar, refinement "
            "transport, and representation or 2-representation on realized matter."
        ),
        "falsification": (
            "A different complete #624 coefficient lane, character, flux action, "
            "or 2-representation survives the source selector, or the selected "
            "action contradicts diagonal-kernel single-valuedness or refinement."
        ),
        "chrome_policy": "Use only after a finite source selector and matter-action packet exists.",
    },
    628: {
        "phase": "operational-load-readback",
        "claim_level": "named-realization source mechanism gap",
        "blocker": (
            "The #625 complete-schema classification forces the counting readback up "
            "to one unit scale and excludes the adjacency-weighted cost from the "
            "operational cone. The physical mechanism that realizes the unit and the "
            "discrete cost from a source-defined process is the gap."
        ),
        "closure": (
            "A finite self-reading patch mechanism emits integer record differences and "
            "the discrete mismatch cost from one source-defined repair and readback process."
        ),
        "falsification": (
            "A rational normalization or inequivalent discrete cost preserves every "
            "operational clause and visible naturality check."
        ),
        "chrome_policy": "Use only to audit a concrete operational counting and cost mechanism.",
    },
    377: {
        "phase": "quark-source-spread",
        "claim_level": "closed locally as theorem-grade two-modulus non-identifiability obstruction",
        "blocker": "Current source identities leave an exact free (R_{>0})^2 action on the up/down spread pair.",
        "closure": "Accepted sharper-obstruction mode is met by a target-free dependency DAG, exact fiber classification, and countermodels with different mass readouts.",
        "falsification": "A theorem-grade source observable built from the stated corpus selects both positive spreads while preserving the declared identities.",
        "chrome_policy": "No workers needed unless a new action-breaking source observable is proposed.",
    },
    378: {
        "phase": "quark-global-classification",
        "claim_level": "separate global classification no-go boundary",
        "blocker": "No class-uniform selector or classification of every admissible public quark frame is emitted.",
        "closure": "Retain the existing corpus-limited global no-go or supply a theorem covering every admissible frame class.",
        "falsification": "An unclassified admissible frame invalidates any global uniqueness wording.",
        "chrome_policy": "No workers needed until a new global classifier packet exists.",
    },
    379: {
        "phase": "quark-source-spread",
        "claim_level": "closed locally as dependent non-identifiability obstruction",
        "blocker": "The generation-bundle path begins from a hand-written template and does not fix the two endpoint spans.",
        "closure": "Accepted sharper-obstruction mode is met by proving that generator splitting fixes profile rays only and leaves two positive moduli.",
        "falsification": "A source-derived generator with certified sector attachment emits both absolute spans without target ancestry.",
        "chrome_policy": "No workers needed until a source-derived generator replaces the template.",
    },
    380: {
        "phase": "quark-source-spread",
        "claim_level": "closed locally as dependent readback obstruction",
        "blocker": "Edge statistics leave two free coefficients in sigma_u=S13+c_u delta and sigma_d=S23+c_d delta.",
        "closure": "Accepted sharper-obstruction mode is met by the target-free two-parameter countermodel family.",
        "falsification": "A source identity fixes both coefficients independently without fitting quark rows.",
        "chrome_policy": "No workers needed until a new coefficient-fixing identity is proposed.",
    },
    381: {
        "phase": "quark-rg-scheme",
        "claim_level": "closed locally as finite-renormalization non-identifiability obstruction",
        "blocker": "Source physics does not choose an external running-mass chart; the light target rows use MSbar at 2 GeV.",
        "closure": "Accepted sharper-obstruction mode is met by the scheme-orbit theorem and explicit comparison-chart declaration.",
        "falsification": "A source-emitted operational observable selects a unique light-quark scheme/scale representative.",
        "chrome_policy": "No workers needed until an operational scheme selector is proposed.",
    },
    382: {
        "phase": "quark-rg-scheme",
        "claim_level": "closed locally as threshold/scheme non-identifiability obstruction",
        "blocker": "Charm and bottom use self-scale MSbar coordinates while top uses a separate pole extraction; no common-scale threshold map is emitted.",
        "closure": "Accepted sharper-obstruction mode is met by the heavy-coordinate scheme theorem and the dimensionful mass-texture audit.",
        "falsification": "A source RG trajectory, threshold map, top conversion, and same-scheme running Higgs normalization emit dimensionless physical Yukawas.",
        "chrome_policy": "No workers needed until a concrete RG/threshold packet exists.",
    },
    608: {
        "phase": "de-sitter-physical-attachment",
        "claim_level": "physical identification gap",
        "blocker": "The finite capacity, shock operator, and exact graph spectrum lack source-derived gauge, kinetic, response-coefficient, and physical-scale attachments.",
        "closure": "One source-derived packet fixes the de Sitter gauge mode, Laplacian normalization, ledger-to-source dictionary, response coefficient, and physical length scale.",
        "falsification": "The required source-derived attachments give the wrong shock sign or fail to reproduce the physical de Sitter mode spectrum.",
        "chrome_policy": "Use for independent audit only after a local attachment packet exists.",
    },
    546: {
        "phase": "charged-lepton-physical-landing",
        "claim_level": "target-anchored empirical closure with source landing open",
        "blocker": (
            "The charged ratios and higher-order remainder in the interval lane "
            "have measured-target ancestry. The 100-decimal-digit certificate "
            "proves outward arithmetic enclosure conditional on those inputs; it "
            "does not emit a source-only Yukawa operator, determinant character, "
            "physical scale, or scheme."
        ),
        "closure": (
            "A source-only physical charged Yukawa landing emits the centered "
            "operator, determinant line, labels, scale, and scheme with a "
            "no-target-leak dependency receipt."
        ),
        "falsification": (
            "A claimed source-only landing depends on measured charged masses, "
            "current-family checksums, a target-calibrated remainder, or an "
            "unfixed determinant normalization."
        ),
        "chrome_policy": (
            "Use the interval artifacts only as empirical-closure diagnostics; "
            "do not present them as prospective charged-lepton predictions."
        ),
    },
    569: {
        "phase": "physical-family-attachment",
        "claim_level": (
            "exact finite rank-45 tensor candidate with separate spin and "
            "local-operator packets; physical attachment open"
        ),
        "blocker": (
            "The family_band_attachment certificate proves the exact selection: among "
            "single complete faithful in-window multiplicity objects inside the screen "
            "coefficient space, the #625 operational cost order has the rank-three band "
            "as unique strict minimizer, selecting a finite rank-three candidate "
            "(Lean Screen/A5FamilyBand.lean, kernel-decided). The #599 simulator "
            "response artifact realizes the comparison clause, and the simulator "
            "pole-residue artifact realizes the realization clause for the response "
            "resolvent: the propagated dynamics has four pole clusters at the band "
            "costs, the rank-three frame residue sits at the lowest positive "
            "generator frequency "
            "(faithful, equivariant, Galois partner at the maximal pole), and the "
            "unitary channel conserves mode norms. The finite screen assembly has "
            "complex rank forty-five with the exact rank-15 generation factor "
            "imported. A conditional finite tensor construction recomputes the "
            "fifteen states and anomaly cancellation, fixes the diagonal Z6 action, "
            "and gives a nondegenerate chirality grading on the separate issue-314 "
            "twelve-port packet. The issue-634 local signed operator tensors with "
            "the rank-45 identity by a declared algebraic extension, which "
            "conditionally copies its exact positive dimensionless gap. This "
            "tensor extension is not source-selected. No source, domain, or "
            "transport bridge attaches the issue-314 "
            "spin packet to that local operator domain. "
            "This is not a physical matter-spectrum measurement. The "
            "#617 copy-count invisibility for external completions is unchanged. "
            "The bounded #627 classification also leaves the physical seam "
            "mechanism, character or 2-representation, and its line/flux action "
            "unselected."
        ),
        "closure": (
            "The matter-pole identification, physical Spin/locality bridge, third "
            "physical persistence leg, laboratory current attachment, exclusion of "
            "extra light sectors, and source-derived seam mechanism with its "
            "character or 2-representation identify the finite rank-45 "
            "response-resolvent realization with the physical matter-pole residue "
            "object and synchronize its line/flux action, or a countermodel separates "
            "them."
        ),
        "falsification": (
            "A physical pole-residue object realized in the screen with a faithful "
            "in-window multiplicity object that is not the rank-three band, or a "
            "physical comparison order that violates the #625 operational cone, "
            "or a selected seam action conflicts with the canonical diagonal "
            "kernel or measured line/flux transport."
        ),
        "chrome_policy": "Use the certificate for the exact selection only; never cite it as an unconditional three-family derivation.",
    },
    671: {
        "phase": "thermodynamic-completion",
        "claim_level": (
            "finite four-law package attained: state-side Axiom-3 projection is "
            "the Gibbs exponential family by the information-projection "
            "Pythagorean identity, transition-side projection onto the "
            "repaired-visible fibre is weighted conditional resampling, the "
            "kernel is stochastic, idempotent, reversible, stationary, and "
            "fixes fibre-measurable charges, relative entropy to the common "
            "reference contracts under repair, the exact first-law split "
            "carries its bilinear cross term, and the fixed-regulator third "
            "law gives the excited-mass bound with the ground-degeneracy "
            "entropy limit; Lean modules and the exact-rational certificate "
            "carry every statement"
        ),
        "blocker": (
            "Five typed receipts stay open. The weighted local objective needs "
            "its declared global representation for both optimizer objects. "
            "The state and transition optimizers need one common "
            "source-derived quotient weight, and the transition matrix may "
            "not be manufactured from a desired equilibrium output. The "
            "source-derived collar transition matrix must pass the "
            "equal-fibre-row matrix receipt or equal the conditional "
            "resampling kernel. One modular charge must be identified with "
            "physical energy and a calibrated clock beyond the attained "
            "finite central-interface split. A continuum third law needs "
            "refinement-uniform gap and multiplicity control. The "
            "strict-descent settlement normalizer carries no entropy "
            "inequality and may not be substituted for the repair kernel."
        ),
        "closure": (
            "Certify each of the five receipts, reject one with a bounded "
            "countermodel, or type it not evaluable at its first missing "
            "producer. A positive THERMO-REALIZATION verdict consumes the "
            "equal-fibre-row matrix receipt of the consensus construction "
            "and may not construct the transition matrix from the target "
            "formula."
        ),
        "falsification": (
            "A certified source collar matrix that fails the "
            "equal-fibre-row receipt while the physical reading is "
            "claimed, or distinct source references for the state and "
            "transition optimizers, voids the physical reading of the "
            "package; the finite theorems fall only with a failed exact "
            "statement on a finite instance."
        ),
        "chrome_policy": (
            "No comparison data enters this issue. Workers handle Lean "
            "extension, certificate replay, and source-matrix receipt "
            "review only; no laboratory temperature or calorimetric "
            "value is consumed."
        ),
    },
    655: {
        "phase": "physical-branch-bridge",
        "claim_level": (
            "frozen prospective branch prediction FZ-11 and exact conditional "
            "continuous-R3 scalar translation adapter; faithful proper-carrier "
            "isometric action and finite recharting cocycle attained; source "
            "selection of the vertex operator, physical attachment, time evolution, "
            "sector identification, and exclusivity open"
        ),
        "blocker": (
            "The response-selected record completion carries a faithful isometric "
            "action of all sixty proper carrier maps, and the declared finite "
            "recharting maps obey their exact cocycle. These maps act on one local "
            "completion and do not prove cofinal scale gluing. The finite source "
            "packet has an exact relabelling into the primitive twelve-port frame, "
            "and the declared paired-difference operator has the frozen cosine symbol "
            "on an auxiliary continuous R3 field. A1 through A3 do not source-select "
            "that vertex operator, derive a physical propagation "
            "sector or time equation, attach the auxiliary field to laboratory "
            "readout, or prove exclusivity. A3 fixes the normalized port-state weights, "
            "and complete directed support on an admitted antipodal frame derives "
            "reciprocity, but no theorem identifies state weights with kinetic hop "
            "rates or excludes generated multi-range terms. The exact normalized "
            "one-range/two-range family changes C4 while preserving the other declared "
            "kinematic symmetries. The physical carrier scale is the separate #664 "
            "gate. FZ-11 failure therefore scores the named branch and not the framework."
        ),
        "closure": (
            "Prove that the axioms force the real, reciprocal, "
            "direct-first-hop primitive-port cosine branch uniquely, including the "
            "state-to-hop action map, for one physical propagation sector, or record a bounded "
            "non-forcing/no-go verdict; a positive proof escalates the "
            "FZ-11 failure scope to the framework."
        ),
        "falsification": (
            "Under the frozen FZ-11 rule at five standard deviations: an "
            "isolated intrinsic positive C4; an isolated intrinsic "
            "anisotropic coefficient at ranks one through five; linked "
            "B0, B6, or the rigid rotated I6 vector excluded after "
            "resolved negative C4; or exclusion of the complete branch "
            "manifold with B0/C4^2 = 10/21 and B6/C4^2 = 32/315."
        ),
        "chrome_policy": (
            "Use only for the bridge/exclusivity theorem or its bounded "
            "no-go; the frozen prediction content may not change."
        ),
    },
    664: {
        "phase": "branch-specific-carrier-scale",
        "claim_level": (
            "exact algebraic carrier-rescaling theorem, conditional port-shell metric "
            "candidate, normalized support-to-response-frame isometry, three-orbit "
            "ambiguity theorem, exact dimensionless seam-current norm, and an exact "
            "joint physical-metric attachment theorem under three named premises; "
            "source selection, dimensionful scale, and a same-action positive lower "
            "bound remain open for both frozen branches"
        ),
        "blocker": (
            "Each frozen coefficient ray leaves its branch action length free. The "
            "exact rescaling counterfamily preserves the dimensionless ray while "
            "varying that length, and neither A1 through A3, the capacity root P, nor "
            "the N closure selects the same-action geometric factor in "
            "a_b^2 = kappa_b P ell_*^2. The D6 seam translation has exact squared "
            "response-pullback norm 2-2/sqrt(5) after the pinned labeled source "
            "identification. Its conditional physical formula additionally assumes "
            "12 a_cell = 4 pi s^2, a_vertex^2 = s^2, and "
            "a_seam(e)^2 = s^2 ||u_e||_G^2. No source theorem inhabits or uniquely "
            "selects those three premises. Port, face, and edge equal-area "
            "conventions give distinct coefficients. A scale or lower bound for one "
            "action does not transfer to the other. Without a source-certified "
            "same-action lower bound, a null result cannot exclude that branch "
            "because its signal can shrink with its action length."
        ),
        "closure": (
            "Derive and replay a quotient-visible carrier-to-cell metric map that "
            "selects or rejects the named port-dual attachment against the face and "
            "edge alternatives; bind the surviving map separately to the exact action "
            "admitted by #655 or #666 and emit a physical interval or lower bound for "
            "that same action by a source-native theorem or independently calibrated "
            "attachment. Otherwise prove bounded nonidentifiability for each branch. "
            "Close only after both frozen branches receive a typed outcome."
        ),
        "falsification": (
            "Two source-admissible physical metric attachments with different "
            "kappa_b, or an exact source rescaling that changes one action length while "
            "preserving every declared observable premise, refute unique scale "
            "selection on that class. After a same-action certified lower bound, a "
            "preregistered sufficiently powered null can reject only the attached "
            "branch."
        ),
        "chrome_policy": (
            "No target comparison data. Use only for the source metric, dimensional "
            "attachment, lower-bound theorem, and independent calibration audit."
        ),
    },
    666: {
        "phase": "source-seam-edge-discriminator",
        "claim_level": (
            "exact source-current D6 image, canonical internal Markov/Dirichlet "
            "action, exact q <= 1 symbol remainder, conditional rank-two transverse "
            "oscillator with an exact zero mode, edge-orbit coefficient ray, data-free "
            "frequency/velocity map, exhaustive synthetic finite-sample stress test, "
            "and frozen FZ-12 conditional prediction; a physical Maxwell sector, "
            "clock, frame, scale, readout, and comparison eligibility remain open"
        ),
        "blocker": (
            "The complete thirty-seam source current has exact D6 image. Under the "
            "pullback of the response-selected Gram metric, rather than the usual "
            "six-dimensional lattice metric, it has dense rank-three completion. "
            "The proper-carrier action is faithful. Record translations are exact "
            "isometries on the embedded D6 lattice and extend uniquely to its metric "
            "completion. A2-natural feasible/objective data plus an A3 unique "
            "minimizer force the source-counting convolution only under the separately "
            "declared complete-edge counting premise, which supplies weight 1/60. Its "
            "local generator satisfies the positive maximum principle, the exact "
            "carre-du-champ identity, and plane-wave diagonalization. A conditional "
            "rank-two oscillator gives omega_aux^2 = lambda_L and an exact zero mode "
            "without constructing Maxwell gauge structure, physical time, or a photon. "
            "The exact internal edge character is "
            "Lambda_int(k) = (6/a_edge^2) lambda_L(k). A separately assumed "
            "physical dilation gives Lambda_a(k) = (6/a^2) "
            "lambda_L((a/a_edge)k), and identifying Lambda_a with physical "
            "omega_phys^2 is another premise. The physical metric candidate in #664 "
            "also remains premise-typed. Exact data-free series arithmetic and a "
            "certified q <= 1 remainder map the frozen ray into formal radial and "
            "transverse velocity coefficients. The frozen edge relations are "
            "B0/C4^2 = 10/21 and B6/C4^2 = -2/63. The finite nonlinear "
            "repair kernel is state dependent and does not descend through the "
            "signed quotient. The internal convolution is not identified with a "
            "homogeneous physical position or field action, same-operator sector, "
            "clock, cofinal refinement, dimensionful scale, frame, analytic "
            "wave packet, or readout. A source-selected electroweak neutral channel "
            "would additionally require a certified (1,2)_(1/2) scalar and a nonzero "
            "neutral lower-component vacuum."
        ),
        "closure": (
            "Identify the exact internal homogeneous action with physical position "
            "and attach it to one source-selected physical sector with a physical "
            "clock, coherent frame and boost law, cofinal refinement, wave-packet "
            "readout, and nuisance contract, or record the precise bounded "
            "nonattachment. Direct time-of-flight uses do not require an electron or "
            "positron sector. Pair-production or decay-threshold uses require "
            "independent electron and positron dispersion and interaction kinematics. "
            "A nonzero-signal comparison may fit the amplitude after the applicable "
            "bridge closes. A branch-falsifying null also requires the same-action "
            "lower bound owned by #664. FZ-12 remains immutable in either case."
        ),
        "falsification": (
            "After the physical observable and covariance gates are frozen, a "
            "nonzero signal can test the scale-free linked FZ-12 manifold while "
            "fitting its amplitude. The leading quadratic term alone is generic and "
            "does not distinguish FZ-12; the full ray requires linked higher-order or "
            "angular sensitivity. Exclusion by a null rejects the attached branch "
            "only if #664 supplies a same-action physical positive lower bound and "
            "the registered test has power over the complete admitted manifold. A "
            "branch failure is OPH-wide only after forcedness and exclusivity are "
            "proved."
        ),
        "chrome_policy": (
            "Do not access comparison data. Use workers only for the source action, "
            "physical attachment, scale, and independent proof audit."
        ),
    },
    667: {
        "phase": "fz12-physical-comparison",
        "claim_level": (
            "frozen source-native edge ray with an exact data-free frequency and "
            "velocity map, certified q <= 1 remainder, and exhaustive synthetic "
            "coverage/power packet; the Fermi direct-timing contract is unarmed, the "
            "exposed Auger threshold projection is diagnostic only, and no physical "
            "FZ-12 comparison is eligible"
        ),
        "blocker": (
            "Issue #666 has not identified the internal edge-current action with "
            "physical position, a photon or other field sector, physical clock, "
            "coherent frame and boost transport, cofinal gluing, an analytic "
            "validity domain, wave-packet readout, or a nuisance model. Published "
            "Fermi-LAT quadratic timing reach corresponds to roughly 4.2e8 Planck "
            "lengths and cannot test a Planck-scale FZ-12 amplitude. Its leading "
            "quadratic statistic is also generic. A photon threshold use needs "
            "electron and positron dispersion plus interaction and composition "
            "kinematics; #670 owns that near-Planck preflight. The exposed Pierre "
            "Auger limit is conditional on its proton-rich source scenario and "
            "supplies no OPH verdict."
        ),
        "closure": (
            "Run the cheapest target-free power envelope before activating any data "
            "comparison. Continue only on a channel whose published or prospective "
            "exposure can distinguish a linked higher-order or angular FZ-12 feature "
            "after the applicable physical bridge closes. Otherwise close this issue "
            "as not activated and retain its exact source packet. Signal mode may fit "
            "the amplitude after the bridge closes. Powered-null mode additionally "
            "requires the same-action lower bound owned by #664."
        ),
        "falsification": (
            "A resolved leading coefficient with the wrong sign rejects the attached "
            "negative-C4 branch, while agreement in that coefficient alone is generic. "
            "A full FZ-12 signal requires at least one linked higher-order coefficient "
            "or angular feature while fitting its amplitude. A null rejects the branch "
            "only with a same-action positive lower bound and registered power over "
            "the complete exact-symbol manifold. An outcome becomes OPH-wide only "
            "after forcedness and exclusivity are proved."
        ),
        "chrome_policy": (
            "The declared Auger projection may be audited as exposed retrospective "
            "data with no score. Do not inspect a prospective comparison payload "
            "until the physical bridge, freeze, nuisance, and power contracts pass."
        ),
    },
    659: {
        "phase": "conditional-angular-discriminator",
        "claim_level": (
            "exact finite selection rule, registered levels-zero-through-five "
            "source declaration, and conditional nonzero continuum stiffness; "
            "complete source selection and physical transfer open"
        ),
        "blocker": (
            "The finite registered ladder emits unit-counting events only through "
            "level five and does not derive their identity from canonical A1-A3 "
            "or provide a complete-tower repair semigroup. The continuum packet "
            "also lacks a physical equivariant/isometric multiplicity-one "
            "transfer, radial-copy exclusion, covariance or operational-response "
            "identification, and screen-to-observable readout."
        ),
        "closure": (
            "Extend the registered ladder to a complete source law or freeze the "
            "equal-seam law as an explicit auxiliary premise; construct an "
            "operational stiffness response or a uniformly coercive inverse; then "
            "prove one physical angular transfer or record a bounded typed exit. "
            "A positive chain may populate only the existing frozen "
            "a5_angular_rules slot in #647."
        ),
        "falsification": (
            "A powered physical comparison is ineligible until the source and "
            "transfer gates close. After a frozen physical map, exclusion of the "
            "declared nonzero normalized rank-six interval rejects that complete "
            "mapped branch without licensing a replacement statistic."
        ),
        "chrome_policy": (
            "Do not access comparison data. Use workers only for source replay, "
            "continuum certification, and physical-transfer review until #647 "
            "declares the existing slot eligible."
        ),
    },
    658: {
        "phase": "native-volume-readout-bridge",
        "claim_level": (
            "exact conditional determinant and positive-volume-ratio kernel; "
            "native physical curvature attachment open"
        ),
        "blocker": (
            "The source does not emit a quotient-visible uniform-density cut or "
            "identify its positive collar readout with the physical spatial-volume "
            "ratio. Physical rechart covariance, conserved stress, adiabaticity, "
            "shear and gradient control, freeze-out transfer, isocurvature, and "
            "phase coherence remain open."
        ),
        "closure": (
            "Construct and replay the native positive-volume readout and the "
            "uniform-density, stress, adiabatic, and transfer premises needed to "
            "invoke the Lean q=zeta implication, or emit a bounded counterexample "
            "or exact not-evaluable boundary."
        ),
        "falsification": (
            "A source-attached positive readout whose volume ratio disagrees with "
            "the curvature scalar on the declared uniform-density class rejects "
            "the bridge; determinant algebra alone is not a physical test."
        ),
        "chrome_policy": (
            "Use only for source and transfer audits. Do not access cosmological "
            "comparison data before the complete bridge enters the frozen registry."
        ),
    },
    663: {
        "phase": "repair-law-adoption",
        "claim_level": (
            "conditional finite-word KL schedule theorem, exact signed source-load "
            "quotient, D6 seam image, and response-metric completion attained; "
            "A1-R/A2-R source and adoption obligations open"
        ),
        "blocker": (
            "The source loads map exactly onto the signed Z6 module, and conservative "
            "seam boundaries have the even-sum D6 image. Equipped with the pullback "
            "response-Gram metric, these modules share the rank-three completion. "
            "The conditional mean descends, while the nonlinear pathwise repair kernel "
            "does not. The basis does not derive a complete quotient-deduplicated "
            "primitive grammar, source counting as the A3 reference, full temporal "
            "PMF-simplex feasibility, adoption of the completion as physical position, "
            "a common refinement semigroup, or the coupled state-generator fixed point."
        ),
        "closure": (
            "Prove every A1-R/A2-R source, temporal, transition, refinement, and "
            "countermodel obligation before adopting the strengthening, or retain "
            "the surviving choice as an explicit branch premise or candidate extra "
            "axiom."
        ),
        "falsification": (
            "Two inequivalent source-admissible grammars, schedules, transition "
            "laws, or refinement completions satisfying the same proposed clauses "
            "refute uniqueness on that declared class."
        ),
        "chrome_policy": (
            "No physical target data. Use only for theorem, certificate, source, "
            "and basis-wide countermodel audits."
        ),
    },
    662: {
        "phase": "fz11-sealed-comparison",
        "claim_level": (
            "dormant one-shot physical comparison; signal mode requires the #655 "
            "physical bridge, while powered-null mode also requires #664"
        ),
        "blocker": (
            "Issue #655 has not proved that the frozen primitive twelve-port "
            "propagation branch is physically attached. A nonzero signal could fit "
            "the amplitude after that bridge closes. Issue #664 has not supplied "
            "the same-action lower bound needed to make a null decisive. No target "
            "payload may be selected, opened, or scored before the applicable "
            "mode-specific dependencies and this issue's freeze and power gates pass."
        ),
        "closure": (
            "After a positive #655 bridge, freeze one release, likelihood and "
            "covariance, nuisance treatment, exposure class, thresholds, and power; "
            "publish the commitment and score the unchanged FZ-11 branch exactly "
            "once. Require #664 only for powered-null activation. If the physical "
            "bridge has no positive exit, close this issue as not activated."
        ),
        "falsification": (
            "At the frozen power threshold, an isolated intrinsic positive C4, "
            "an intrinsic anisotropic coefficient at ranks one through five, or "
            "exclusion of the complete rigid FZ-11 branch manifold rejects the "
            "physically attached branch."
        ),
        "chrome_policy": (
            "Do not launch a comparison worker or inspect the target payload "
            "before #655 closes positively, the mode-specific #664 gate is resolved, "
            "and the sealed protocol hash is public."
        ),
    },
    594: {
        "phase": "physical-wz-source-to-pole",
        "claim_level": "OPH-native source-to-pole packet open",
        "blocker": (
            "Issue #646 has not emitted the minimal scale-free source contract, "
            "and the finite parents do not emit one target-clean action, complete "
            "W/Z-coupled census and Yukawa packet, FJ map, full RG/matching law, "
            "Lorentzian/Spin quantum-EFT transfer, or hermetic common-digest "
            "production replay. The physical-unit row is not evaluable on the "
            "declared clock interface."
        ),
        "closure": (
            "Issue #646 first emits the minimal scale-free source-input or "
            "non-identifiability verdict. Issues #569, #630, #631, #632, #32, "
            "and #635 then emit final positive or negative source verdicts. Any "
            "positive native row is substituted into the validated #593 consumer "
            "without algorithm changes and replays one recomputed subject digest; "
            "both output rows must end in final positive, rigorous negative, or "
            "NOT_EVALUABLE verdicts."
        ),
        "falsification": (
            "Any source choice, branch, tolerance, error budget, or producer depends "
            "on measured W, Z, Higgs, top, G_F, weak mixing, or a calibrated proxy."
        ),
        "chrome_policy": "Use only to audit a complete local source packet and replay bundle.",
    },
    630: {
        "phase": "electroweak-source-action",
        "claim_level": "bounded non-promoting frontier; positive scalar/Yukawa action open",
        "blocker": (
            "Finite current, matter, global-form, family, and negative-control "
            "receipts classify conditional operator spaces but select no physical "
            "scalar, normalized potential, vacuum, complete Yukawa matrices, or FJ map."
        ),
        "closure": (
            "One target-clean local source action emits the scalar carrier, kinetic "
            "normalization, stable potential and vacuum, complete Yu/Yd/Ye, and a "
            "symbolic v_chart-to-v_F map with scheme and uncertainty."
        ),
        "falsification": (
            "The positive action cannot distinguish the registered scalar and "
            "Higgs/top countermodels or imports a measured particle target."
        ),
        "chrome_policy": "Use only to audit a positive source-action packet; retain all negative controls.",
    },
    631: {
        "phase": "local-electroweak-carrier",
        "claim_level": "finite line theorem closed; physical common-load attachment open",
        "blocker": (
            "The finite screen and weak order-unit lines admit a unique positive "
            "unital map, but no source receipt identifies them as one physical load "
            "or derives its normalization."
        ),
        "closure": (
            "The positive #630 action supplies the weak/scalar line and one "
            "refinement-natural target-clean receipt identifies and normalizes the "
            "physical local screen/electroweak carrier."
        ),
        "falsification": (
            "A second admissible physical carrier, broken refinement map, wrong "
            "normalization, or target/cosmological input survives the controls."
        ),
        "chrome_policy": "Use only after the #630 positive action and local carrier receipt exist.",
    },
    503: {
        "phase": "direct-n-physical-carrier",
        "claim_level": "conditional Einstein theorem complete; inhabited source carrier open",
        "blocker": (
            "No single target-clean source-derived tower instantiates every "
            "geometry, modular, event, stress, entropy, vacuum, scale, and "
            "refinement premise. A finite direct N cannot supply that tower."
        ),
        "closure": (
            "Emit one inhabited common-domain Einstein/de Sitter tower with "
            "independent source provenance and certified tails, or the final "
            "bounded negative or not-evaluable verdict."
        ),
        "falsification": (
            "The tower reads N, Lambda, a cosmological target, or a conclusion "
            "coordinate while constructing any source premise."
        ),
        "chrome_policy": (
            "Use only to audit a complete local source tower; #589, not #503, "
            "owns any downstream cosmological comparison."
        ),
    },
    589: {
        "phase": "direct-n-horizon-attachment",
        "claim_level": "finite N and horizon-area identity open",
        "blocker": (
            "Issue #505 has not emitted a positive finite N together with an "
            "inhabited compatible #503 Einstein/de Sitter carrier, and no "
            "refinement-natural horizon-record order-unit identity is proved."
        ),
        "closure": (
            "Prove the target-clean horizon-record identity and emit the physical "
            "Lambda-times-Planck-area relation, or close bounded no-go or "
            "NOT_EVALUABLE_NO_HORIZON_RECORD_ATTACHMENT."
        ),
        "falsification": (
            "A finite N is compared with Lambda before the #503 carrier and "
            "#589 identity exist, or the public cosmological value is described "
            "as an unexposed prospective target."
        ),
        "chrome_policy": (
            "Use for attachment audit only after positive #505 and #503 packets; "
            "any present-day Lambda comparison is retrospective."
        ),
    },
    639: {
        "phase": "exposure-typed-dimensionless-forecast",
        "claim_level": (
            "C1 static draft controls complete; C3 invariant registry, "
            "exposure typing, custody, and scoring open"
        ),
        "blocker": (
            "The first-ranked direct N candidate is ineligible: #551 and #505 "
            "remain open at the universal source-membership and direct-capacity "
            "antecedent, while #589 closed not evaluable. The pointer advances "
            "on the complete frozen #647 "
            "registry. No physically typed registry row has passed exposure "
            "classification, durable custody, minimum-power review, and "
            "single-use unsealing, and the scoring surface awaits its "
            "separately reviewed unsealing change with custody and exposure "
            "typing."
        ),
        "closure": (
            "Score at most one eligible attached dimensionless candidate under "
            "its prospective, blind-postdiction, or exposed-retrospective class, "
            "or close with the complete bounded no-go or not-evaluable ledger. "
            "Exploratory rows are catalog-only."
        ),
        "falsification": (
            "An exposed target is called blind, producer strata become sequential "
            "second chances, a finite N bypasses #503/#589, minimum power or "
            "multiplicity is omitted, a post-unseal defect is repaired in place, "
            "or a composite failure is attributed to one auxiliary branch without "
            "an independent premise-isolation test."
        ),
        "chrome_policy": (
            "Use workers for adversarial ancestry and checker audits on the "
            "#647 inventory. Do not launch comparison evaluation until one "
            "attached packet passes the exposure, ranking, power, and freeze gates."
        ),
    },
    641: {
        "phase": "baryon-violation-boundary",
        "claim_level": "dimension-six operator census open; minimal X/Y channel absence closed",
        "blocker": (
            "The selected product adjoint excludes the ordinary minimal simple-GUT "
            "X/Y gauge-exchange generator, but no complete gauge-invariant baryon-"
            "violating operator basis through dimension six is classified on the "
            "selected matter and scalar boundary."
        ),
        "closure": (
            "Enumerate every gauge-invariant baryon-violating operator through "
            "dimension six, classify its OPH source status and mediators, and "
            "separate exact absences from open coefficients and ultraviolet channels."
        ),
        "falsification": (
            "A gauge-invariant operator or admitted mediator is omitted, an open "
            "coefficient is treated as zero, or the minimal X/Y channel result is "
            "promoted to general proton stability or a lifetime bound."
        ),
        "chrome_policy": (
            "Use only to audit a local symbolic operator basis and its completeness "
            "proof; no hadron or proton-lifetime computation is required."
        ),
    },
    642: {
        "phase": "global-form-discriminator",
        "claim_level": "finite Z6 packet closed; laboratory discriminator open",
        "blocker": (
            "The exact finite character, cocharacter, and line-class arithmetic "
            "has no completed continuum or laboratory attachment, particle mass, "
            "production rate, or discovery-level comparison rule."
        ),
        "closure": (
            "Emit the exact charge and line-operator incompatibility packet, "
            "physicalization verdict, frozen decision rule, independent replay, "
            "and explicit nondetection boundary."
        ),
        "falsification": (
            "A representation or line class is omitted, conventions are mixed, "
            "a source rate is invented, or experimental nondetection is counted "
            "as positive evidence without a predicted tested region."
        ),
        "chrome_policy": (
            "Use only to audit the exact finite discriminator and a read-only "
            "experimental crosswalk."
        ),
    },
    643: {
        "phase": "a5-angular-discriminator",
        "claim_level": (
            "finite branching and first invariant exact; frame lock and "
            "screen-to-sky observable open"
        ),
        "blocker": (
            "The pinned FZ-02 theorem packet does not separately prove the "
            "registered level-three/level-six frame lock. The A5 action may be "
            "internal, and no source-forced nonzero sky statistic or physical "
            "screen-to-sky map exists."
        ),
        "closure": (
            "Derive or remove the frame-lock clause, decide internal versus "
            "spatial action, and emit one nonzero rotation-controlled statistic "
            "with frozen foreground, covariance, trials, data manifest, and "
            "comparison contract for #639, or close with the exact no-go."
        ),
        "falsification": (
            "Known low-multipole anomalies select the statistic or frame, a free "
            "amplitude is promoted, or A5 invariance shared by isotropy is treated "
            "as a differentiating prediction."
        ),
        "chrome_policy": (
            "Use only to audit the frozen statistic and comparison contract. "
            "Issue #639 alone may access and score the comparison payload."
        ),
    },
    644: {
        "phase": "structural-specificity-audit",
        "claim_level": "bounded alternative-carrier specificity score open",
        "blocker": (
            "No carrier and response-law null ensemble is frozen, so the rarity "
            "of the finite gauge, quotient, matter, and family hit tuple is not "
            "computable."
        ),
        "closure": (
            "Freeze an exhaustive bounded grammar and equal algorithm budget, "
            "run every model, independently replay the full hit table, and emit "
            "the finite specificity or exact non-computability verdict."
        ),
        "falsification": (
            "Models are added or removed after scoring, isomorphic duplicates "
            "remain, target-aware pruning occurs, or the OPH carrier receives a "
            "larger selector or repair budget."
        ),
        "chrome_policy": (
            "Use only to audit a frozen local ensemble and independently "
            "recomputed full result table."
        ),
    },
    645: {
        "phase": "observer-overlap-interferometer",
        "claim_level": "source-to-readout and normalized cross-spectrum open",
        "blocker": (
            "No source map identifies a finite overlap observable with optical "
            "path or phase. Auto-spectrum detector noise also prevents the "
            "normalized measured cross-spectrum from being amplitude-free."
        ),
        "closure": (
            "Derive one nonconstant geometry-dependent source cross-spectrum, "
            "the source-to-readout coupling, signal/background decomposition, "
            "calibrated likelihood, and frozen comparison contract for #639, or "
            "close immediately with the declared no-map verdict."
        ),
        "falsification": (
            "An excluded holographic-noise template is relabelled as OPH, "
            "instrument noise is omitted, or a free amplitude or shape absorbs "
            "the comparison."
        ),
        "chrome_policy": (
            "Use only to audit a local source-to-readout packet and frozen "
            "comparison contract. Issue #639 alone may score public Holometer data."
        ),
    },
    646: {
        "phase": "electroweak-invariant-search",
        "claim_level": "scale-free W/Z invariant or no-go open",
        "blocker": (
            "The dependence of pole, width, residue, and asymmetry combinations "
            "on clock, vacuum scale, normalization, scheme, threshold, scalar, "
            "Yukawa, family, and continuum directions is not classified."
        ),
        "closure": (
            "Emit an exact or interval-certified invariant multi-output vector "
            "with its minimal source-input contract, or prove every frozen "
            "candidate remains non-identifiable."
        ),
        "falsification": (
            "A cancellation holds only at one point, measured W/Z values select "
            "the combination, weak-angle schemes are mixed, or a surviving open "
            "direction is omitted."
        ),
        "chrome_policy": (
            "Use only to audit the local dependency/Jacobian packet and "
            "target-free invariant proof."
        ),
    },
    647: {
        "phase": "oph-invariant-mining",
        "claim_level": "systematic OPH-only observable mine open",
        "blocker": (
            "The source-feature, typed nuisance, source-admissible completion, "
            "candidate-specific baseline, end-to-end physicalization, and exposed-"
            "data registries are not complete, and no candidate has entered #639."
        ),
        "closure": (
            "Freeze deterministic ranking and the bounded grammar before "
            "candidate generation, certify each OPH model-image inclusion "
            "globally and each baseline non-inclusion by counterexample, repeat "
            "the nuisance audit after physicalization, then emit exposure-typed "
            "candidates or the complete negative verdict."
        ),
        "falsification": (
            "A discrete branch is removed by a Jacobian, an exposed value "
            "influences mining, ranking remains discretionary, historical "
            "multiplicity or minimum power is omitted, direct N re-enters the "
            "fallback registry, or retrospective agreement is promoted to "
            "prospective survival. A composite failure cannot identify one "
            "auxiliary branch without an independent premise-isolation test."
        ),
        "chrome_policy": (
            "Use only for proof and ancestry audits against the frozen local "
            "grammar; #639 owns all comparison access and one-shot scoring."
        ),
    },
    632: {
        "phase": "electroweak-field-census",
        "claim_level": "complete W/Z-coupled census gap",
        "blocker": (
            "The finite representation inventory does not prove that every omitted "
            "field or operator has an identically zero W/Z vertex across the selected "
            "source action and EFT intervals."
        ),
        "closure": (
            "A source-complete field and operator census is emitted with exact "
            "zero-vertex decoupling certificates for every excluded object and "
            "mutation tests for sterile or source-invisible direct sums."
        ),
        "falsification": (
            "An omitted sector contributes to a W/Z self-energy, counterterm, "
            "threshold, mixing, or Yukawa vertex on an admissible completion."
        ),
        "chrome_policy": "Use only to audit a complete census and zero-vertex packet.",
    },
    633: {
        "phase": "operational-source-clock",
        "claim_level": (
            "positive dimensionless gap retained; physical units not evaluable "
            "on the declared serialized interface"
        ),
        "blocker": (
            "The bounded serialized-interface audit finds no emitted physical-unit "
            "field, and two attained producer runs ignore the named SI attachment "
            "channel. This is not a transitive source-closure or complete-domain "
            "non-identifiability theorem."
        ),
        "closure": (
            "The receipt certifies "
            "PHYSICAL_UNITS_NOT_EVALUABLE_ON_DECLARED_SERIALIZED_INTERFACE "
            "while retaining the positive dimensionless gap. Any physical "
            "transition on an extended matter domain is owned separately by "
            "#569 and #630."
        ),
        "falsification": (
            "A field in the declared serialized interface carries a physical-unit "
            "attachment, or the attained producer demonstrably consumes the named "
            "SI channel."
        ),
        "chrome_policy": (
            "Treat the closure as a bounded interface verdict; never cite it as "
            "a complete-domain clock no-go."
        ),
    },
    314: {
        "phase": "standard-model-source-realization",
        "claim_level": "conditional finite matter fixture exact; physical source selection open",
        "blocker": (
            "The exterior, anomaly, and finite Spin calculations consume a "
            "declared matrix-current fixture. No same-source producer selects "
            "the physical current action or matter representation."
        ),
        "closure": (
            "A source-passed current from #566 acts on a source-selected matter "
            "carrier and the same finite Spin implementers, with refinement "
            "commuting squares and adversarial controls."
        ),
        "falsification": (
            "The declared exterior/Spin theorem fails internally, or every "
            "target-free source matter producer selects an incompatible action."
        ),
        "chrome_policy": "Use only to audit a concrete same-source current/matter packet.",
    },
    566: {
        "phase": "standard-model-source-realization",
        "claim_level": "abstract A1/A2 Lie type forced; executable source current open",
        "blocker": (
            "The source artifact derives R=-J and response signs but contains no "
            "ordered current tomography, exact bracket reconstruction, or "
            "same-current closed overlap holonomy."
        ),
        "closure": (
            "One target-free source packet reconstructs twelve generators, their "
            "closed bracket, all proper-carrier projective implementers, "
            "response-component membership, and refinement intertwiners."
        ),
        "falsification": (
            "A model satisfying the explicit A1/A2 response premises has another "
            "compact Lie type, or the source producer cannot realize the forced "
            "type without importing a named current model."
        ),
        "chrome_policy": "Use only to audit raw ordered histories and an executable verifier.",
    },
    567: {
        "phase": "standard-model-source-realization",
        "claim_level": "declared-table Z6 kernel exact; physical global form open",
        "blocker": (
            "The order-six axis class uses declared diagonal and zero-sum "
            "relations. Character completeness and a same-source "
            "loop-to-kernel identity are absent."
        ),
        "closure": (
            "A source-derived complete character/relation lattice and genuine "
            "line category identify the carrier loop with the common physical "
            "current and matter kernel, with alternative global forms rejected."
        ),
        "falsification": (
            "Another source-admissible relation or transparent central action "
            "survives, or the loop map fails on an admitted character."
        ),
        "chrome_policy": "Use only to audit a concrete source character and loop packet.",
    },
    634: {
        "phase": "lorentzian-spin-local-source",
        "claim_level": (
            "finite causal and local-operator domain attained; continuum "
            "Lorentzian/Spin quantum-EFT promotion open"
        ),
        "blocker": (
            "The finite event complex, typed sections, seam topology, and local "
            "operators are inhabited. Negative cone margins, one Euclidean "
            "neighborhood fit, and the missing cofinal refinement limit block "
            "continuum promotion."
        ),
        "closure": (
            "The bounded finite-domain receipt remains closed. A separate "
            "regulator-controlled transfer must establish the continuum "
            "Lorentzian/Spin quantum-EFT domain."
        ),
        "falsification": (
            "The finite receipt fails replay or silently promotes its fitted "
            "finite coordinates to a continuum manifold."
        ),
        "chrome_policy": (
            "Use only to audit the bounded local-domain construction and exact "
            "bundle replay; #635 owns continuum transfer."
        ),
    },
    670: {
        "phase": "fz12-uhe-threshold-preflight",
        "claim_level": (
            "target-free near-Planck sensitivity work order; the exact photon-seam "
            "ray is frozen, while the electron/positron sector, interaction law, "
            "source and shower nuisances, and prospective power envelope are open"
        ),
        "blocker": (
            "The exposed Pierre Auger analysis constrains one alternative proton-rich "
            "Lorentz-violation scenario near the Planck scale, but its reference "
            "source scenarios do not supply the same bound and its shower response "
            "was not recomputed for the OPH ray. It is an exposed method fixture, not "
            "an OPH score. A valid threshold map requires independently sourced "
            "electron and positron dispersion, the conservation and composition law, "
            "pair-production and photon-decay interactions, source composition and "
            "backgrounds, atmospheric-shower propagation, detector classification, "
            "and a controlled exact-symbol remainder."
        ),
        "closure": (
            "Before opening any future data rows, derive and independently replay the "
            "joint photon-lepton threshold equations, record every source, interaction, "
            "composition, shower, detector, and nuisance premise, and compute a "
            "target-free AugerPrime power envelope over the complete attached branch. "
            "Arm a separately frozen prospective comparison only if that envelope has "
            "registered near-Planck sensitivity to a linked FZ-12 feature; otherwise "
            "close as not activated with the exact source packet retained."
        ),
        "falsification": (
            "Only a prospectively frozen, sufficiently powered outcome may reject the "
            "physically attached threshold branch. The result becomes OPH-wide only "
            "if source selection and exclusivity are proved. The existing exposed "
            "Auger rows cannot retrospectively falsify or confirm OPH."
        ),
        "chrome_policy": (
            "Do not inspect future event rows. The published Auger methods and exposed "
            "aggregate limits may be used only to construct and audit the target-free "
            "power and nuisance contract."
        ),
    },
}


# V3 is adequacy-first: these rows classify the issue contracts without
# inferring scientific closure from GitHub state.  Add each new ``[V3-*]``
# lane here before regenerating the ledger.  The builder and offline validator
# fail closed when a V3-labelled issue has no entry, which makes later lanes
# (for example #740+) a one-row declarative addition rather than a parser
# heuristic.
V3_ISSUE_POLICY: dict[int, dict[str, str]] = {
    740: {
        "phase": "v3-common-world-integration",
        "claim_level": "program-level same-architecture compatibility gate",
        "blocker": (
            "The lane results are carried by separate conditional witnesses; no "
            "artifact yet proves that one AV-n jointly inhabits their carriers, "
            "records, clocks, fields, matter, actions, repair semantics, scales, "
            "and physical readouts."
        ),
        "closure": (
            "One machine-checked common-world manifest joins every row included in "
            "a program-level adequacy claim on one inhabited AV-n, with typed maps, "
            "exact premise and empirical-input ancestry, receipts, and audit pointers."
        ),
        "falsification": (
            "Incompatible witness types, clocks, calibrations, architecture versions, "
            "or an owed required observation row block the common-world claim."
        ),
        "chrome_policy": (
            "Use an independent worker only to audit a sealed compatibility packet; "
            "dependency closure is not positive evidence and cannot fill a missing map."
        ),
    },
    741: {
        "phase": "v3-infrastructure-architecture-versioning",
        "claim_level": "immutable architecture, promotion, and audit index",
        "blocker": (
            "Promotions and audit assertions are not uniformly keyed to an immutable "
            "axiom/protocol/premise snapshot, reviewed commit, reviewer identity, "
            "repair history, and artifact hashes."
        ),
        "closure": (
            "The initial AV-n manifest, append-only decision history, semantic diff, "
            "promotion invalidation/replay rules, prediction-lineage guard, and "
            "machine-readable audit index land with mutation gates and independent audit."
        ),
        "falsification": (
            "A semantic architecture change silently preserves a promotion or frozen "
            "prediction, or an audit/promotion cannot be traced to its exact version "
            "and artifacts, invalidates the custody infrastructure."
        ),
        "chrome_policy": (
            "Workers may audit a committed version/diff/index packet; they cannot "
            "supply scientific evidence or decide an architecture change implicitly."
        ),
    },
    742: {
        "phase": "v3-l11-cosmology-astrophysics",
        "claim_level": "conditional cosmology and astrophysics adequacy lane",
        "blocker": (
            "No premise-complete common-world composition yet joins expansion, "
            "perturbations, structure, compact-object/gravitational propagation, "
            "dark-sector, ultra-high-energy, and birefringence targets to explicit "
            "clocks, matter, calibration, and empirical inputs."
        ),
        "closure": (
            "Every scoped target has an evidence-backed conservative status and "
            "atomic premise ancestry; qualifying compositions share one AV-n and "
            "eligible instruments and comparisons are preregistered under #737/#738."
        ),
        "falsification": (
            "A fitted background is promoted as structure formation, missing physics "
            "is relabeled as derived dark matter, or inspected data are called a "
            "prediction without a fresh freeze."
        ),
        "chrome_policy": (
            "Audit sealed, target-clean finite packets only; do not inspect comparison "
            "data to choose a cosmological mechanism or nuisance model."
        ),
    },
    743: {
        "phase": "v3-l12-interacting-qft-rg-scattering",
        "claim_level": "interacting-QFT, RG, and scattering composition gap",
        "blocker": (
            "The finite/free surfaces do not select one interaction, regulator, "
            "renormalization prescription, scale map, interacting state, asymptotic "
            "or finite operational scattering interface, and controlled error budget."
        ),
        "closure": (
            "At least one non-free same-AV-n benchmark composes end to end with "
            "declared regulator/scheme dependence, RG flow, unitarity/locality and "
            "scattering claims, error control, receipts, mutation gates, and audit."
        ),
        "falsification": (
            "An interaction symbol on a free theory, unrelated theorem re-exports, "
            "or an uncontrolled finite-volume receipt is used as continuum QFT, RG, "
            "or an S-matrix."
        ),
        "chrome_policy": (
            "Use workers for independent audit of a concrete regulator-explicit "
            "packet only; no target coupling or subtraction-point fitting."
        ),
    },
    744: {
        "phase": "v3-l13-qcd-hadron-nuclear-atomic",
        "claim_level": "layered strong-interaction-to-atomic adequacy lane",
        "blocker": (
            "Perturbative QCD, confinement, hadron, nuclear, and atomic targets lack "
            "a layered contract that separates source calculations from empirical "
            "transport and states the nonperturbative and laptop-scale boundaries."
        ),
        "closure": (
            "The perturbative and atomic benchmarks close within declared errors, "
            "each nonperturbative layer has an exact evidence-backed status/blocker, "
            "and no spectral or transport payload is relabeled as generated theory."
        ),
        "falsification": (
            "Empirical Ward projection is promoted to a hadron spectrum, one layer's "
            "success promotes another, or a fit/measurement is called a source-only "
            "prediction."
        ),
        "chrome_policy": (
            "Workers may audit finite perturbative/atomic or classification packets; "
            "do not use them as a substitute for nonperturbative compute."
        ),
    },
    745: {
        "phase": "v3-l14-electroweak-weak-flavor-neutrino",
        "claim_level": "conditional electroweak, weak, flavor, and neutrino composition",
        "blocker": (
            "No one action/version jointly supplies the scalar potential and broken "
            "vacuum, chiral weak currents, physical W/Z/H sheet, Yukawa matrices, "
            "mixing and CP structure, neutrino mechanism, calibration, and readout."
        ),
        "closure": (
            "One AV-n action inhabits the scalar, breaking, weak-current, flavor, "
            "and neutrino objects; all numbers and premises are classified, loop/RG "
            "imports are scoped, and ledger surfaces, receipts, gates, and audit agree."
        ),
        "falsification": (
            "A term table is called an assembled action, chart-local poles are called "
            "particles, fitted Yukawa/mixing inputs are called derived, or tree-level "
            "receipts are promoted to loop phenomenology."
        ),
        "chrome_policy": (
            "Audit sealed action/correspondence packets only; no target mass, angle, "
            "phase, or neutrino hierarchy may select the construction."
        ),
    },
    726: {
        "phase": "v3-infrastructure-observation-ledger",
        "claim_level": "standing adequacy-accounting and promotion gate",
        "blocker": (
            "The observation ledger must remain exhaustive and keep each physics "
            "target bound to one adequacy rung, owning lane, complete premise list, "
            "and evidence bundle; the ledger itself is not evidence."
        ),
        "closure": (
            "Every declared physics target has one machine-readable row and one "
            "generated Markdown row, deterministic parity checks pass, structural "
            "and emergent promotions cite audited lane receipts, and predictive "
            "rows point to frozen custody contracts."
        ),
        "falsification": (
            "An omitted target, unnamed consumed premise, unsupported rung "
            "promotion, moving numerical target, or source/generated drift invalidates "
            "the ledger state."
        ),
        "chrome_policy": (
            "Use an independent worker only to audit a committed ledger and its "
            "receipts; never let a worker promote a row from prose alone."
        ),
    },
    727: {
        "phase": "v3-infrastructure-premise-register",
        "claim_level": "standing conditionality and disposition gate",
        "blocker": (
            "Every premise consumed by a composition must be registered with an "
            "exact statement, type, consuming lanes, disposition, and the relevant "
            "source-selection no-go or provenance evidence."
        ),
        "closure": (
            "The machine register and generated surface are exhaustive and in "
            "parity, every consumer cites a registered row, numerical ancestry "
            "distinguishes P and N from derived values and flagged imports, and "
            "each disposition change is an audited recorded event."
        ),
        "falsification": (
            "A hidden premise, silent disposition change, unflagged numerical or "
            "empirical input, missing consumer edge, or claim that registration "
            "makes a premise true invalidates the register contract."
        ),
        "chrome_policy": (
            "Use workers only for independent premise-ancestry or no-go audits of "
            "a committed packet; architecture decisions remain explicit local edits."
        ),
    },
    728: {
        "phase": "v3-l1-spacetime-adequacy",
        "claim_level": (
            "conditional spacetime composition with a preregistered emergent-signature gate"
        ),
        "blocker": (
            "Dimensionality, Lorentz kinematics, causal order, event-manifold data, "
            "and the measured (1,3) signature must be composed with all stable-"
            "causality, open-image, affine, source, and finite-versus-continuum "
            "premises explicit; the emergent signature requires its independent "
            "preregistered replication and controls."
        ),
        "closure": (
            "One audited composed statement exists per spacetime target with the "
            "registered premise list and exact rung, the signature instrument has a "
            "frozen decision rule and reported controls, and any failed event or "
            "source premise is recorded rather than hidden."
        ),
        "falsification": (
            "A signature verdict read from the update rule, omitted replication "
            "control, unregistered premise, or continuum claim supported only by a "
            "finite receipt rejects the claimed promotion."
        ),
        "chrome_policy": (
            "Local composition and instrument receipts first; workers may audit a "
            "frozen replication packet without changing its analysis or controls."
        ),
    },
    729: {
        "phase": "v3-l2-gravitation-adequacy",
        "claim_level": (
            "conditional effective-gravity composition with calibrated and predictive boundaries"
        ),
        "blocker": (
            "The Einstein relation, finite Newtonian regime, horizon thermality, "
            "gravitational-wave stance, and N-linked cosmological row require one "
            "common typed tower, registered stress/entropy/scale premises, the "
            "declared repair law, and frozen comparison custody."
        ),
        "closure": (
            "The registered Einstein composition, finite inverse-square limit, and "
            "horizon temperature/entropy package are composed without an independent "
            "temperature input; gravitational-wave and cosmological rows retain their "
            "frozen or constants-lane contracts."
        ),
        "falsification": (
            "Target-theory curvature or stress imported as source data, calibration "
            "bypassed, a continuum claim inferred from finite receipts, or Lambda "
            "tuned outside the N ancestry invalidates the gravity composition."
        ),
        "chrome_policy": (
            "Use workers only to audit a concrete same-tower composition, finite "
            "Newtonian receipt, or horizon packet; no post-hoc gravitational fit."
        ),
    },
    730: {
        "phase": "v3-l3-quantum-adequacy",
        "claim_level": "conditional finite quantum and structural-field-theory composition",
        "blocker": (
            "Born and Lueders rules, dynamics, Tsirelson, no-signalling, phase-"
            "complete tomography, and regional-net structure must be composed from "
            "the registered algebra-state, effect-additivity, phase-instrument, "
            "carrier, and correlation premises with the finite boundary explicit."
        ),
        "closure": (
            "Each quantum target has one audited composed statement and premise "
            "list; the phase instrument supplies real outcomes rather than fabricated "
            "data; the locally covariant net/covariance/time-slice result is stated at "
            "its proved finite level; affinity discharge remains on #739."
        ),
        "falsification": (
            "An unregistered quantum representation, fabricated phase outcome, "
            "post-hoc statistic, or continuum/vacuum/local-field language beyond the "
            "finite receipts invalidates the promotion."
        ),
        "chrome_policy": (
            "Workers may audit a committed algebra, instrument, or finite-net packet; "
            "they may not supply missing quantum outcomes or promote continuum QFT."
        ),
    },
    731: {
        "phase": "v3-l4-mechanics-adequacy",
        "claim_level": "exact finite mechanics composition under registered premises",
        "blocker": (
            "Stationary action, the Euler-Lagrange/Hamilton bridge, and Noether "
            "conservation require the registered path reference and real enrichment, "
            "with source selection, saddle scope, units, clocks, and amplitudes "
            "kept outside the finite theorem."
        ),
        "closure": (
            "One mechanics surface composes the transition-derived action, committed "
            "mode/minimizer boundary, Legendre bridge, and Noether witnesses over the "
            "named register rows; reference/enrichment/saddle discharge remains #739."
        ),
        "falsification": (
            "A target-chosen action, consequence-selected enrichment, tuned reference "
            "or multiplier, or stationary-phase/continuum claim beyond the committed "
            "finite scope rejects the composition."
        ),
        "chrome_policy": (
            "Use workers only for proof or countermodel audit of a committed finite "
            "mechanics surface; do not search target dynamics to choose an enrichment."
        ),
    },
    732: {
        "phase": "v3-l5-thermodynamics-adequacy",
        "claim_level": "conditional finite four-law package under the declared repair law",
        "blocker": (
            "The four laws, fluctuation identities, Landauer bound, and repair arrow "
            "consume one common faithful reference/repaired-fibre law plus the named "
            "receipts; laboratory meaning additionally consumes explicit energy and "
            "clock calibration."
        ),
        "closure": (
            "One audited surface composes the finite four-law and fluctuation package "
            "under the registered repair law, keeps calibration explicit, and opens "
            "an emergent run only after an architecture export is recorded and the "
            "instrument is preregistered."
        ),
        "falsification": (
            "A stationary law manufactured from the desired equilibrium, distinct "
            "state/transition references, fitted coupling or calibration, or post-hoc "
            "run labeled validation invalidates the physical reading."
        ),
        "chrome_policy": (
            "Use workers for exact finite theorem, source-reference, or preregistered "
            "instrument audits only; no laboratory target fitting."
        ),
    },
    733: {
        "phase": "v3-l6-electromagnetism-adequacy",
        "claim_level": (
            "conditional finite Maxwell-shaped composition with frozen dispersion branches"
        ),
        "blocker": (
            "Maxwell-shaped equations, transverse massless propagation, Gauss/Coulomb "
            "structure, and light-signal emergence require the registered counting, "
            "oscillator, physical-frequency, sector, frame, exclusivity, and "
            "calibration bridges; dispersion verdicts remain frozen."
        ),
        "closure": (
            "One finite U(1) composition supplies the equations, rank-two transverse "
            "propagation, zero mode, and Gauss/Coulomb receipts under named premises; "
            "frozen dispersion rows keep custody, and the emergent light-signal run "
            "uses a preregistered instrument."
        ),
        "falsification": (
            "A forbidden hand-inserted kinetic term, frequency tuned to c, comparison "
            "outside the frozen rule, or continuum/laboratory Maxwell claim without "
            "the registered bridge invalidates the promotion."
        ),
        "chrome_policy": (
            "Audit committed finite-field, propagation, or frozen-branch packets only; "
            "never inspect comparison data to choose the field bridge."
        ),
    },
    734: {
        "phase": "v3-l7-standard-model-structure",
        "claim_level": "conditional Standard Model gauge, matter, family, and global-form composition",
        "blocker": (
            "The gauge algebra, one-generation matter content, chirality, family band, "
            "and global form must be composed over explicit selection, grammar, "
            "loop-to-kernel, and physical-chirality premises; source dynamics and "
            "character identities remain discharge work."
        ),
        "closure": (
            "One audited structure surface states exactly which gauge, matter, "
            "chirality, family, and global-form rows are attained or conditional, "
            "names every premise, and leaves coupling values and action assembly to "
            "their owning lanes."
        ),
        "falsification": (
            "Target-named generators, anomaly-informed table extension, a family "
            "count read from observation, or an unregistered global-form identity "
            "invalidates the structure claim."
        ),
        "chrome_policy": (
            "Use workers only for independent replay of a committed classification, "
            "matter, family, or kernel packet; no target-guided grammar changes."
        ),
    },
    735: {
        "phase": "v3-l8-standard-model-lagrangian",
        "claim_level": "term-by-term conditional effective-action correspondence",
        "blocker": (
            "Gauge and matter kinetic terms, scalar/Yukawa and electroweak-breaking "
            "stances, proton-stability boundary, and the assembled action must be "
            "classified term by term as derived, registered premise, partial, or "
            "absent; coupling values belong to #736."
        ),
        "closure": (
            "A generated correspondence table and composed action cover every "
            "textbook Standard Model Lagrangian term, cite the exact receipts and "
            "complete premise register, preserve absent/partial rows, and make no "
            "renormalized or laboratory claim."
        ),
        "falsification": (
            "A term transcribed from the target action, scalar structure tuned for "
            "symmetry breaking, an unsupported row upgraded, or a coupling number "
            "introduced outside #736 invalidates the assembly."
        ),
        "chrome_policy": (
            "Use workers only for a term-by-term audit of the committed correspondence "
            "and action; do not fill absent terms from the textbook target."
        ),
    },
    736: {
        "phase": "v3-l9-masses-and-constants",
        "claim_level": "conditional quantitative determinations with frozen predictive subsets",
        "blocker": (
            "Fine structure, N/Lambda, lepton, quark/mixing, electroweak, hadron, "
            "neutrino, and calibration rows require complete numerical ancestry: P "
            "and N are the only proposed fundamental free numerical parameters, while "
            "empirical numerical payloads and every "
            "selection or physical bridge remain explicitly registered."
        ),
        "closure": (
            "Each quantitative row has a certified enclosure or exact status, full "
            "input ancestry and falsifier; P and N routes state every selection and "
            "same-quantity bridge; target-anchored or predictive subsets retain their "
            "compare-only or frozen-custody classification."
        ),
        "falsification": (
            "An unflagged numerical input, target-anchored quantity called a "
            "prediction, post-comparison interval widening, correction selected by "
            "proximity, or failed registered kill band invalidates the row."
        ),
        "chrome_policy": (
            "Use workers for ancestry, interval, or frozen-rule audit of committed "
            "packets only; no target access during source construction."
        ),
    },
    737: {
        "phase": "v3-l10-emergent-instruments",
        "claim_level": "preregistered simulator-instrument and emergent-rung gate",
        "blocker": (
            "Each emergent target requires an instrument with its observable, "
            "analysis, controls, decision rule, architecture digest, and owning "
            "ledger row frozen before deterministic runs; architecture changes must "
            "first be recorded in the premise register."
        ),
        "closure": (
            "Every instrument binds exactly one observation row and composition lane, "
            "replays from committed artifacts, reports positive, negative, and control "
            "outcomes equally, and promotes or blocks only the rung named by its "
            "preregistered rule."
        ),
        "falsification": (
            "Post-hoc analysis, parameter revision after data, silent reruns, omitted "
            "controls, irreproducible configuration, or a laboratory claim from a "
            "laptop simulation invalidates the instrument verdict."
        ),
        "chrome_policy": (
            "Workers may independently replay or audit a frozen laptop-scale packet; "
            "they may not revise the instrument after seeing its run."
        ),
    },
    738: {
        "phase": "v3-standing-custody",
        "claim_level": "standing custody, comparison, and falsification enforcement",
        "blocker": (
            "This issue intentionally never closes: it continuously guards frozen "
            "targets, exposure classes, compare-only ancestry, independent audits, "
            "the tracking cascade, and immutable registration artifacts."
        ),
        "closure": (
            "A clean standing state means every prediction has a pre-exposure target, "
            "rule, band, and precision floor; every promotion survives independent "
            "audit; negative findings remain visible; and issue, ledger, scoreboard, "
            "and commit ordering stays synchronized."
        ),
        "falsification": (
            "Retrospective work called prediction, rewritten frozen bytes, hidden "
            "negative outcomes, promotion without audit, or a ledger change without "
            "receipts is a custody failure."
        ),
        "chrome_policy": (
            "Use workers only for adversarial audit of sealed artifacts; never use a "
            "worker to manufacture evidence or inspect embargoed comparison data."
        ),
    },
    739: {
        "phase": "v3-deferred-premise-discharge",
        "claim_level": "deferred necessity and axiomatization queue; adequacy remains conditional",
        "blocker": (
            "Every premise-register row marked remove or axiomatize retains its "
            "source-selection, representation, affinity, clock, bracket, matter, "
            "character, repair-export, or spectral-information obligation until an "
            "audited derivation or recorded architecture decision retires it."
        ),
        "closure": (
            "For each queued row, either derive it from the current basis with a "
            "producer/verifier packet and bounded countermodels, record and propagate "
            "an explicit axiomatization decision, or retain an honest scoped no-go; "
            "the queue closes only when no remove/axiomatize obligation remains."
        ),
        "falsification": (
            "Silent promotion, deleting a bounded negative, treating deferral as "
            "derivation, or claiming necessity while any consumed premise remains "
            "undischarged invalidates the necessity claim."
        ),
        "chrome_policy": (
            "Local source theorem or countermodel packet first; workers may audit one "
            "named discharge attempt but cannot decide axiomatization implicitly."
        ),
    },
}


# Explicitly classify the one parked non-V3 continuation proposal that remains
# open beside the V3 program.  It is deliberately non-promoting.
ISSUE_POLICY[715] = {
    "phase": "parked-external-extension-review",
    "claim_level": "parked external algebraic proposal; no physical promotion",
    "blocker": (
        "The external inductive-limit and admissible-instrument definitions, finite "
        "Python witnesses, and Lean-style sketch have not been pinned, independently "
        "replayed, or integrated with the canonical OPH algebra and repair interfaces."
    ),
    "closure": (
        "Import a hash-pinned review packet, replay the finite witnesses, formalize "
        "or reject the restriction-lift and W-invariance lemmas, and either absorb "
        "the definitions with explicit non-physical boundaries or archive the proposal."
    ),
    "falsification": (
        "A finite compatible-family counterexample to restriction lift, an admitted "
        "repair leaving W, or dependence on an extra axiom rejects the proposed layer; "
        "none of its finite witnesses closes a physical bridge."
    ),
    "chrome_policy": (
        "Do not launch a worker until the external revision is hash-pinned into a "
        "local review packet; then use one only for independent theorem/replay audit."
    ),
}


CLOSED_OUT_OF_SCOPE_ISSUES: dict[int, dict[str, str]] = {
    153: {
        "title": "hadron backend/systematics",
        "url": "https://github.com/FloatingPragma/observer-patch-holography/issues/153",
        "closed_status": "closed_not_planned_out_of_scope_computationally_blocked",
        "closure_note": (
            "Closed as out-of-scope, not solved. The current pipeline emits no hadron predictions "
            "until a working OPH hadron backend on suitable hardware publishes production output and systematics."
        ),
    },
    157: {
        "title": "nonperturbative QCD hadron branch",
        "url": "https://github.com/FloatingPragma/observer-patch-holography/issues/157",
        "closed_status": "closed_not_planned_out_of_scope_computationally_blocked",
        "closure_note": (
            "Closed as out-of-scope, not solved. Specialist-paper hadron claims remain suppressed until "
            "a working OPH backend emits the required nonperturbative QCD/hadron data."
        ),
    },
}


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _ascii(text: str) -> str:
    return text.encode("ascii", "ignore").decode("ascii")


def _run_gh() -> list[dict[str, Any]]:
    completed = subprocess.run(
        [
            "gh",
            "issue",
            "list",
            "-R",
            REPO,
            "--state",
            "open",
            "--limit",
            "500",
            "--json",
            "number,title,labels,updatedAt,url,body,milestone",
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    return json.loads(completed.stdout)


def _normalise_contract_text(text: str) -> str:
    return " ".join(text.split()).strip()


def _contract_field(body: str, name: str) -> str | None:
    """Read one compact V2 field even when several fields share a paragraph."""

    pattern = re.compile(
        rf"(?is){V2_FIELD_PREFIX}{re.escape(name)}:\s*(.*?)"
        rf"(?=(?:{V2_FIELD_PREFIX}(?:{V2_FIELD_MARKER}):)|\Z)"
    )
    matches = [_normalise_contract_text(match) for match in pattern.findall(body)]
    matches = [match for match in matches if match]
    if not matches:
        return None
    if len(matches) != 1:
        raise ValueError(f"V2 contract has {len(matches)} {name!r} fields")
    return matches[0]


def _v2_scope(body: str) -> str:
    """Return the first issue-scope paragraph, excluding administrative fields."""

    for paragraph in re.split(r"\n\s*\n", body):
        text = _normalise_contract_text(paragraph)
        if not text:
            continue
        marker_match = re.search(
            rf"(?is){V2_FIELD_PREFIX}(?:{V2_FIELD_MARKER}):",
            text,
        )
        if marker_match is not None:
            text = text[: marker_match.start()].strip()
        if text:
            return text
    raise ValueError("V2 contract has no issue-scope paragraph")


def _v2_policy(issue: dict[str, Any]) -> dict[str, str] | None:
    """Derive administrative ledger metadata from one live V2 issue contract.

    The derivation deliberately does not synthesize a scientific falsifier or
    infer that a dependency is discharged. It exposes the issue's own scope,
    deliverables, boundary, and dependency text and supplies only status-hygiene
    rules around them.
    """

    title = str(issue.get("title") or "")
    title_match = V2_TITLE_RE.match(title)
    labels = sorted(
        {
            str(label.get("name") or "").strip()
            for label in issue.get("labels") or []
            if str(label.get("name") or "").strip()
        }
    )
    track_labels = [label for label in labels if label.startswith("track:")]
    is_v2 = title_match is not None or bool(track_labels)
    if not is_v2:
        return None

    number = int(issue.get("number") or 0)
    where = f"V2 issue #{number}"
    if title_match is None:
        raise ValueError(f"{where} title lacks a task code such as [A1]")
    if len(track_labels) != 1:
        raise ValueError(
            f"{where} must carry exactly one track:* label; got {track_labels}"
        )
    track_label = track_labels[0]
    track = V2_TRACKS.get(track_label)
    if track is None:
        raise ValueError(f"{where} carries unknown track label {track_label!r}")
    task_code = title_match.group(1)
    if task_code[0] != track["code"]:
        raise ValueError(
            f"{where} task code {task_code} conflicts with {track_label}"
        )

    size_labels = [label for label in labels if label.startswith("size:")]
    if any(label not in V2_SIZE_LABELS for label in size_labels):
        raise ValueError(f"{where} carries an unknown size label: {size_labels}")
    state_labels = [label for label in labels if label in V2_STATE_LABELS]
    if len(state_labels) > 1:
        raise ValueError(f"{where} carries conflicting state labels: {state_labels}")
    standing = state_labels == ["standing"]
    if standing:
        if size_labels:
            raise ValueError(f"{where} standing lane must not carry a size label")
    elif len(size_labels) != 1:
        raise ValueError(f"{where} must carry exactly one V2 size label")

    body = str(issue.get("body") or "").strip()
    if not body:
        raise ValueError(f"{where} has no live issue-body contract")
    depends = _contract_field(body, "Depends on")
    wave = _contract_field(body, "Wave")
    if depends is None:
        raise ValueError(f"{where} has no Depends on: contract")
    if wave is None:
        raise ValueError(f"{where} has no Wave: contract")

    normalized_wave = wave.removesuffix(".").strip()
    if normalized_wave.casefold() == "standing":
        declared_wave: int | str = "standing"
    else:
        wave_match = re.match(r"^V2-W([1-9]\d*)\b", normalized_wave)
        if wave_match is None:
            raise ValueError(f"{where} has malformed Wave: contract {wave!r}")
        declared_wave = int(wave_match.group(1))
        suffix = normalized_wave[wave_match.end() :].strip()
        if suffix not in {"", V2_WAVE_TITLES.get(declared_wave)}:
            raise ValueError(f"{where} has malformed Wave: contract {wave!r}")
    if declared_wave != track["wave"]:
        raise ValueError(
            f"{where} declares wave {declared_wave!r}, expected {track['wave']!r}"
        )
    if standing != (declared_wave == "standing"):
        raise ValueError(f"{where} standing label and Wave: contract disagree")
    milestone = issue.get("milestone") or {}
    milestone_title = str(milestone.get("title") or "")
    if isinstance(declared_wave, int):
        milestone_match = re.match(r"^V2-W([1-9]\d*)\b", milestone_title)
        if milestone_match is None:
            raise ValueError(f"{where} has no V2 wave milestone")
        if int(milestone_match.group(1)) != declared_wave:
            raise ValueError(
                f"{where} Wave: contract conflicts with milestone {milestone_title!r}"
            )

    scope = _v2_scope(body)
    deliverables = _contract_field(body, "Deliverables")
    boundary = _contract_field(body, "Boundary")
    depends_none = bool(re.match(r"(?i)^none\b", depends))
    if standing:
        if task_code == "G2":
            phase = "v2-standing-discriminator-production"
            claim_level = "standing target-clean discriminator production"
            blocker = f"Standing producer scope with no issue prerequisite: {scope}"
            closure = (
                "Complete a declared target-clean bridge, obstruction, power, or "
                "exposure contract without rewriting frozen payloads or unsealing "
                "a comparison; issue state alone is not a scientific verdict."
            )
            chrome_policy = (
                "Use only to audit a concrete target-clean bridge or power packet; "
                "never use a worker to manufacture comparison evidence."
            )
        else:
            phase = f"v2-standing-{track['slug']}"
            claim_level = "standing custody/comparison control"
            blocker = f"Standing scope with no issue prerequisite: {scope}"
            closure = (
                "Keep the live custody, ledger, and validator scope current under "
                "its frozen eligibility rules; issue state alone is not a scientific "
                "verdict."
            )
            chrome_policy = (
                "Use only to audit a concrete sealed custody or validator packet; "
                "never use a worker to manufacture comparison evidence."
            )
    else:
        phase = f"v2-w{declared_wave}-{track['slug']}"
        if state_labels == ["optional"]:
            claim_level = "optional non-blocking investigation"
        elif state_labels == ["parked"]:
            claim_level = "parked stage-gated construction"
        elif state_labels == ["blocked"]:
            claim_level = "blocked construction"
        else:
            claim_level = "active issue-scoped construction"

        if state_labels == ["parked"] and depends_none:
            dependency_text = "Parked by V2 staging; no hard issue prerequisite."
        elif state_labels == ["parked"]:
            dependency_text = f"Parked; live prerequisite contract: {depends}"
        elif state_labels == ["blocked"] and depends_none:
            dependency_text = "Blocked by its live status; no issue prerequisite is named."
        elif state_labels == ["blocked"]:
            dependency_text = f"Blocked; live prerequisite contract: {depends}"
        elif depends_none:
            dependency_text = "No hard issue prerequisite."
        else:
            dependency_text = f"Live prerequisite contract: {depends}"
        blocker = f"{dependency_text} Open scope: {scope}"

        if deliverables is not None:
            closure = (
                "Verify the live deliverables at their declared status without "
                f"scientific promotion from issue state: {deliverables}"
            )
        else:
            closure = (
                "Complete the declared scope with theorem/certificate evidence "
                "or record an explicit scoped countermodel, no-go, or typed "
                f"not-evaluable exit. Declared scope: {scope}"
            )

        if state_labels == ["parked"]:
            chrome_policy = (
                "Do not launch workers while parked or while a listed prerequisite "
                "is open; begin from a concrete local issue packet."
            )
        elif state_labels == ["optional"]:
            chrome_policy = (
                "Use only to audit a concrete optional criteria/theorem packet; "
                "this lane blocks no other issue."
            )
        elif state_labels == ["blocked"]:
            chrome_policy = (
                "Do not launch workers until the recorded blocker changes and a "
                "concrete local issue packet exists."
            )
        else:
            chrome_policy = (
                "Local theorem or certificate work first; use an independent "
                "worker only to audit a concrete issue-scoped packet."
            )

    if boundary is not None:
        falsification = (
            "Scope-control failure if a result is promoted past this live boundary: "
            f"{boundary}"
        )
    else:
        falsification = (
            "Closure is invalid if the delivered artifact does not cover the live "
            "scope, a hard prerequisite remains unresolved, or theorem status "
            "exceeds its evidence. Scientific falsifiers remain in the canonical "
            "claim-local ledgers."
        )

    return {
        "phase": phase,
        "claim_level": claim_level,
        "blocker": blocker,
        "closure": closure,
        "falsification": falsification,
        "chrome_policy": chrome_policy,
    }


def _v3_policy(issue: dict[str, Any]) -> dict[str, str] | None:
    """Return the explicit V3 contract, failing closed on unregistered lanes."""

    title = str(issue.get("title") or "")
    labels = {
        str(label.get("name") or "").strip().casefold()
        for label in issue.get("labels") or []
        if str(label.get("name") or "").strip()
    }
    title_match = V3_TITLE_RE.match(title)
    if title_match is None and "v3" not in labels:
        return None

    number = int(issue.get("number") or 0)
    where = f"V3 issue #{number}"
    if title_match is None:
        raise ValueError(f"{where} title lacks a task code such as [V3-L1]")
    policy = V3_ISSUE_POLICY.get(number)
    if policy is None:
        raise ValueError(
            f"{where} has no explicit ledger policy; add it to V3_ISSUE_POLICY"
        )
    if set(policy) != POLICY_FIELDS:
        raise ValueError(
            f"{where} policy fields differ: "
            f"missing={sorted(POLICY_FIELDS - set(policy))}, "
            f"extra={sorted(set(policy) - POLICY_FIELDS)}"
        )
    for field, value in policy.items():
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{where} policy has empty {field}")
        if value == PLACEHOLDER_POLICY.get(field):
            raise ValueError(f"{where} policy retains placeholder {field}")
    if not policy["phase"].startswith("v3-"):
        raise ValueError(f"{where} policy phase must start with 'v3-'")
    return dict(policy)


def _fallback_policy(issue: dict[str, Any]) -> dict[str, str]:
    return dict(PLACEHOLDER_POLICY)


def build_ledger(issues: list[dict[str, Any]]) -> dict[str, Any]:
    rows = []
    open_numbers = {int(issue["number"]) for issue in issues}
    for issue in sorted(issues, key=lambda item: item["number"]):
        number = int(issue["number"])
        policy = _v3_policy(issue)
        if policy is None:
            policy = _v2_policy(issue)
        if policy is None:
            policy = ISSUE_POLICY.get(number, _fallback_policy(issue))
        labels = [label["name"] for label in issue.get("labels", [])]
        rows.append(
            {
                "number": number,
                "title": _ascii(issue["title"]),
                "url": issue["url"],
                "labels": labels,
                "updated_at": issue["updatedAt"],
                **policy,
            }
        )
    closed_out_of_scope_records = []
    for number, record in sorted(CLOSED_OUT_OF_SCOPE_ISSUES.items()):
        if number in open_numbers:
            continue
        policy = ISSUE_POLICY[number]
        closed_out_of_scope_records.append(
            {
                "number": number,
                "title": record["title"],
                "url": record["url"],
                "closed_status": record["closed_status"],
                "closure_note": record["closure_note"],
                **policy,
            }
        )
    return {
        "artifact": "oph_open_problem_ledger",
        "generated_utc": _now_utc(),
        "repo": REPO,
        "open_issue_count": len(rows),
        "closed_out_of_scope_count": len(closed_out_of_scope_records),
        "worker_policy": {
            "chrome_pro_workers_default": "local_first",
            "max_parallel_workers": 6,
            "launch_condition": "only after a concrete theorem, audit, or implementation packet exists",
            "obstruction_only_result_allowed": False,
        },
        "closed_out_of_scope_records": closed_out_of_scope_records,
        "rows": rows,
    }


def render_markdown(ledger: dict[str, Any]) -> str:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in ledger["rows"]:
        grouped[row["phase"]].append(row)

    lines = [
        "# OPH Open Problem Ledger",
        "",
        f"Generated: `{ledger['generated_utc']}` from live GitHub issues in `{ledger['repo']}`.",
        "",
        "This is the public boundary between closed OPH claims, conditional claims, and open work. "
        "Dedicated GitHub issues remain canonical for task state; this ledger records the current "
        "claim level, missing artifact, closure criterion, falsification route, and Chrome Pro worker policy.",
        "",
        "Worker policy: local artifacts first; up to six Chrome Pro workers may be used only after a "
        "concrete theorem, audit, or implementation packet exists. Obstruction-only worker output is not accepted.",
        "",
        f"Open issue count: `{ledger['open_issue_count']}`",
        "",
    ]
    closed_records = ledger.get("closed_out_of_scope_records", [])
    if closed_records:
        lines.extend(
            [
                "## Closed Out-Of-Scope Records",
                "",
                "These are not solved derivations. They are closed as non-current-scope tasks because "
                "the required computation depends on a working OPH hadron backend and suitable hardware.",
                "",
                "| Issue | Status | Claim level | Closure note | Reopen criterion | Chrome policy |",
                "| --- | --- | --- | --- | --- | --- |",
            ]
        )
        for row in closed_records:
            lines.append(
                f"| [#{row['number']}]({row['url']}) {row['title']} | `{row['closed_status']}` | "
                f"`{row['claim_level']}` | {row['closure_note']} | {row['closure']} | {row['chrome_policy']} |"
            )
        lines.append("")
    for phase in sorted(grouped):
        lines.extend([f"## {phase}", ""])
        lines.extend(
            [
                "| Issue | Claim level | Missing artifact / blocker | Closure criterion | Falsification route | Chrome policy |",
                "| --- | --- | --- | --- | --- | --- |",
            ]
        )
        for row in grouped[phase]:
            lines.append(
                f"| [#{row['number']}]({row['url']}) {row['title']} | `{row['claim_level']}` | "
                f"{row['blocker']} | {row['closure']} | {row['falsification']} | {row['chrome_policy']} |"
            )
        lines.append("")
    return "\n".join(lines).rstrip()


def validate_committed_ledger(json_path: Path, markdown_path: Path) -> list[str]:
    """Validate the tracked ledger offline, without GitHub credentials."""
    problems: list[str] = []
    try:
        ledger = json.loads(json_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return [f"cannot read ledger JSON: {exc}"]
    if ledger.get("artifact") != "oph_open_problem_ledger":
        problems.append("unexpected ledger artifact identifier")
    if ledger.get("repo") != REPO:
        problems.append("unexpected ledger repository")

    rows = ledger.get("rows")
    if not isinstance(rows, list):
        return problems + ["ledger rows must be a list"]
    required = {
        "number",
        "title",
        "url",
        "labels",
        "updated_at",
        "phase",
        "claim_level",
        "blocker",
        "closure",
        "falsification",
        "chrome_policy",
    }
    numbers: list[int] = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            problems.append(f"row {index} is not an object")
            continue
        missing = sorted(required - set(row))
        if missing:
            problems.append(f"row {index} is missing fields: {', '.join(missing)}")
        number = row.get("number")
        if not isinstance(number, int):
            problems.append(f"row {index} has a non-integer issue number")
        else:
            numbers.append(number)
        labels = row.get("labels")
        if isinstance(labels, list):
            track_labels = [
                label
                for label in labels
                if isinstance(label, str) and label.startswith("track:")
            ]
            normalized_labels = {
                label.casefold() for label in labels if isinstance(label, str)
            }
        else:
            track_labels = []
            normalized_labels = set()
        is_v3 = (
            V3_TITLE_RE.match(str(row.get("title") or "")) is not None
            or "v3" in normalized_labels
        )
        if is_v3:
            policy = V3_ISSUE_POLICY.get(number) if isinstance(number, int) else None
            if policy is None:
                problems.append(
                    f"row {index} has no explicit V3 policy in V3_ISSUE_POLICY"
                )
            phase = str(row.get("phase") or "")
            if not phase.startswith("v3-"):
                problems.append(f"row {index} has an unclassified V3 phase")
            for field in POLICY_FIELDS:
                value = row.get(field)
                if not isinstance(value, str) or not value.strip():
                    problems.append(f"row {index} has an empty V3 {field}")
                elif value == PLACEHOLDER_POLICY.get(field):
                    problems.append(
                        f"row {index} retains placeholder V3 {field}"
                    )
                elif policy is not None and value != policy[field]:
                    problems.append(
                        f"row {index} V3 {field} differs from explicit policy"
                    )
        if track_labels:
            phase = str(row.get("phase") or "")
            if not phase.startswith("v2-"):
                problems.append(f"row {index} has an unclassified V2 phase")
            if len(track_labels) == 1 and track_labels[0] in V2_TRACKS:
                track = V2_TRACKS[track_labels[0]]
                if track["wave"] == "standing":
                    if int(row.get("number") or 0) == 704:
                        expected_phase = "v2-standing-discriminator-production"
                    else:
                        expected_phase = f"v2-standing-{track['slug']}"
                else:
                    expected_phase = f"v2-w{track['wave']}-{track['slug']}"
                if phase != expected_phase:
                    problems.append(
                        f"row {index} V2 phase {phase!r} does not match "
                        f"{track_labels[0]}"
                    )
            placeholder_values = {
                field: PLACEHOLDER_POLICY[field]
                for field in ("blocker", "closure", "falsification")
            }
            for field, placeholder in placeholder_values.items():
                value = row.get(field)
                if not isinstance(value, str) or not value.strip():
                    problems.append(f"row {index} has an empty V2 {field}")
                elif value == placeholder:
                    problems.append(f"row {index} retains placeholder V2 {field}")
    if len(numbers) != len(set(numbers)):
        problems.append("ledger has duplicate issue numbers")
    if numbers != sorted(numbers):
        problems.append("ledger rows are not sorted by issue number")
    if ledger.get("open_issue_count") != len(rows):
        problems.append("open_issue_count does not equal the number of rows")

    closed = ledger.get("closed_out_of_scope_records")
    if not isinstance(closed, list):
        problems.append("closed_out_of_scope_records must be a list")
    elif ledger.get("closed_out_of_scope_count") != len(closed):
        problems.append(
            "closed_out_of_scope_count does not equal the number of records"
        )

    # Some release branches publish only the canonical JSON snapshot. If the
    # optional human-readable mirror exists, require exact synchronization;
    # its absence is not a hidden-input requirement.
    if markdown_path.exists():
        try:
            expected_markdown = render_markdown(ledger) + "\n"
            actual_markdown = markdown_path.read_text(encoding="utf-8")
        except (OSError, KeyError, TypeError) as exc:
            problems.append(f"cannot render/read ledger Markdown: {exc}")
        else:
            if actual_markdown != expected_markdown:
                problems.append(
                    "OPEN_PROBLEMS.md is out of sync with the ledger JSON"
                )
    return problems


def compare_committed_to_live(
    json_path: Path,
    live_issues: list[dict[str, Any]],
) -> list[str]:
    """Fail when the committed gate snapshot differs from live GitHub state.

    ``generated_utc`` is build metadata and is deliberately ignored. Every
    scientific field -- including open issue membership, title, labels,
    ``updated_at``, closure policy, and closed out-of-scope records -- must
    match the ledger rebuilt from the live issue list.
    """
    try:
        committed = json.loads(json_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return [f"cannot read committed ledger JSON: {exc}"]
    live = build_ledger(live_issues)
    committed.pop("generated_utc", None)
    live.pop("generated_utc", None)
    if committed == live:
        return []

    committed_open = {
        row.get("number")
        for row in committed.get("rows", [])
        if isinstance(row, dict)
    }
    live_open = {
        row.get("number")
        for row in live.get("rows", [])
        if isinstance(row, dict)
    }
    problems = []
    if committed_open != live_open:
        problems.append(
            "open issue membership differs from live GitHub: "
            f"missing={sorted(live_open - committed_open)}, "
            f"stale={sorted(committed_open - live_open)}"
        )
    else:
        problems.append(
            "committed ledger rows or policy fields differ from live GitHub"
        )
    return problems


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the OPH open-problem ledger from GitHub issues.")
    parser.add_argument("--json-out", default=str(DEFAULT_JSON_OUT))
    parser.add_argument(
        "--markdown-out",
        default=None,
        help=(
            "optional path for a rendered Markdown mirror; the committed "
            "surface is the JSON snapshot, and no mirror is written unless "
            "this flag is given"
        ),
    )
    parser.add_argument("--print-json", action="store_true")
    parser.add_argument(
        "--check",
        action="store_true",
        help="validate the committed JSON and Markdown offline; do not call GitHub",
    )
    parser.add_argument(
        "--check-live",
        action="store_true",
        help=(
            "compare the committed JSON with a fresh public GitHub issue query; "
            "requires gh authentication and ignores generated_utc only"
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    json_out = Path(args.json_out)
    markdown_out = Path(args.markdown_out) if args.markdown_out else DEFAULT_MD_OUT
    write_markdown = args.markdown_out is not None
    if args.check and args.check_live:
        raise SystemExit("choose at most one of --check and --check-live")
    if args.check:
        problems = validate_committed_ledger(json_out, markdown_out)
        if problems:
            print("open-problem ledger FAILED:")
            for problem in problems:
                print(f"  - {problem}")
            return 1
        print(
            "open-problem ledger OK: committed JSON is valid "
            "(Markdown mirror checked when present)"
        )
        return 0
    if args.check_live:
        problems = compare_committed_to_live(json_out, _run_gh())
        if problems:
            print("open-problem ledger LIVE CHECK FAILED:")
            for problem in problems:
                print(f"  - {problem}")
            return 1
        print("open-problem ledger OK: committed gate snapshot matches live GitHub")
        return 0

    ledger = build_ledger(_run_gh())
    json_text = json.dumps(ledger, indent=2, sort_keys=True, ensure_ascii=True) + "\n"

    json_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.write_text(json_text, encoding="utf-8")

    if write_markdown:
        markdown_out.write_text(render_markdown(ledger) + "\n", encoding="utf-8")

    if args.print_json:
        print(json_text, end="")
    else:
        print(f"saved: {json_out}")
        if write_markdown:
            print(f"saved: {markdown_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
