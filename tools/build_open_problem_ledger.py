#!/usr/bin/env python3
"""Build the public OPH open-problem ledger from live GitHub issues."""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPO = "FloatingPragma/observer-patch-holography"
DEFAULT_JSON_OUT = ROOT / "tracking" / "open_issues" / "open_problem_ledger.json"
DEFAULT_MD_OUT = ROOT / "OPEN_PROBLEMS.md"


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
        "claim_level": "physical particle criterion open after finite-interface no-go",
        "blocker": (
            "The exact finite flux chains, twisted spectra, fusion rows, and one-step "
            "refinement data have hash-identical classical-lattice and quantum-Hilbert "
            "interpretations. No asymptotic dynamics, pole, or all-depth intertwiner is emitted."
        ),
        "closure": (
            "A source-defined quantum dynamics or equivalent physical spectral criterion "
            "excludes the same-flux classical model and controls composition, asymptotic "
            "states, and refinement."
        ),
        "falsification": (
            "The proposed particle discriminator remains a function only of the finite "
            "interface shared by the classical and quantum countermodels."
        ),
        "chrome_policy": "Use only to audit a concrete dynamics, pole, or quantization packet.",
    },
    334: {
        "phase": "gravity-scale-composition",
        "claim_level": "source clock plus gravity readout composition gap",
        "blocker": (
            "Issue #633 proves PHYSICAL_UNITS_NOT_EVALUABLE on the complete declared "
            "finite domain. The extended matter and scalar/Yukawa domains in #569 "
            "and #630 must select a physical transition. A source-derived gravity "
            "length or radius, Einstein attachment, common provenance graph, and "
            "outward interval composition into G_SI remain open."
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
            "exact finite rank-45 chirality, Spin-data, locality, and gap "
            "attachment; physical and continuum attachment open"
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
            "imported. The finite issue-634 matter attachment recomputes the fifteen "
            "states and anomaly cancellation, fixes the diagonal Z6 action, gives a "
            "nondegenerate chirality grading, acts on the full rank-45 fiber with "
            "support-radius-one signed operators, inherits the exact positive "
            "dimensionless gap, and records the finite Spin data and lift ambiguity. "
            "This is not a physical matter-spectrum measurement. The "
            "#617 copy-count invisibility for external completions is unchanged. "
            "The bounded #627 classification also leaves the physical seam "
            "mechanism, character or 2-representation, and its line/flux action "
            "unselected."
        ),
        "closure": (
            "The matter-pole identification, continuum Spin/locality limit, third "
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
    594: {
        "phase": "physical-wz-source-to-pole",
        "claim_level": "OPH-native source-to-pole packet open",
        "blocker": (
            "The finite parents and bounded frontiers do not emit one target-clean "
            "action, complete W/Z-coupled census and Yukawa packet, FJ map, full "
            "RG/matching law, source point or law, operational clock, and hermetic "
            "common-digest production replay."
        ),
        "closure": (
            "The positive packets from #634, #569, #630, #631, #632, and #32, "
            "including the physical transition selected on the extended matter "
            "domain, are substituted into the validated #593 consumer without "
            "algorithm changes and independently replayed under one recomputed "
            "subject digest."
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
        "claim_level": "positive dimensionless gap closed; physical SI clock non-identifiable on the declared domain",
        "blocker": (
            "The declared source domain emits an exact positive dimensionless gap "
            "and admits source-indistinguishable SI attachments. It therefore cannot "
            "select a physical reference transition or SI chart."
        ),
        "closure": (
            "The source classification and same-domain countermodels certify "
            "PHYSICAL_UNITS_NOT_EVALUABLE while retaining the positive dimensionless "
            "gap. Any matter-anchored extension is owned separately by #569 and #630."
        ),
        "falsification": (
            "The declared source domain uniquely selects an SI attachment, or the "
            "same-domain alternatives cease to be source-indistinguishable."
        ),
        "chrome_policy": "Closed on the complete negative exit; use successor issues for extended source domains.",
    },
    634: {
        "phase": "lorentzian-spin-local-source",
        "claim_level": "inhabited local action domain open",
        "blocker": (
            "Conditional event, frame, cone, and spin interfaces do not yet construct "
            "one inhabited refinement-compatible Lorentzian/Spin local action domain."
        ),
        "closure": (
            "One source-derived local domain binds events, Lorentzian frames, Spin "
            "lifts, locality, Hamiltonian typing, refinement, and public receipts "
            "with independent replay."
        ),
        "falsification": (
            "The proposed domain is uninhabited, loses locality or Spin/Lorentzian "
            "compatibility under refinement, or requires a physical target."
        ),
        "chrome_policy": "Use only to audit a concrete inhabited local-action construction.",
    },
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
            "number,title,labels,updatedAt,url",
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    return json.loads(completed.stdout)


def _fallback_policy(issue: dict[str, Any]) -> dict[str, str]:
    return {
        "phase": "unclassified",
        "claim_level": "open",
        "blocker": "Classify blocker from the live issue body.",
        "closure": "Add exact closure criterion to this ledger.",
        "falsification": "Add exact falsification criterion to this ledger.",
        "chrome_policy": "Do not launch workers until the issue has a concrete local packet.",
    }


def build_ledger(issues: list[dict[str, Any]]) -> dict[str, Any]:
    rows = []
    open_numbers = {int(issue["number"]) for issue in issues}
    for issue in sorted(issues, key=lambda item: item["number"]):
        policy = ISSUE_POLICY.get(int(issue["number"]), _fallback_policy(issue))
        labels = [label["name"] for label in issue.get("labels", [])]
        rows.append(
            {
                "number": int(issue["number"]),
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
