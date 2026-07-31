#!/usr/bin/env python3
"""Three-axiom drift guard (AR1).

Enforces the canonical basis in ``claims/axiom_registry.yaml`` across active
surfaces:

1. the registry declares exactly three core axioms and a valid status enum;
2. no active surface states or counts a five-axiom basis, uses ``OPH5``,
   presents Axiom 4 as a core recovery principle, or presents Axiom 5 /
   Minimal Admissible Realization as a core economy principle;
3. designated entry surfaces carry the three-axiom basis;
4. mathematical uses of the group :math:`A_5` pass untouched.

The allowlist distinguishes archives, dated audits, provenance records,
bibliography keys, Lean kernel-axiom audits, migration documentation, and
deliberate stale-text test fixtures from active reader-facing prose.
``--pdf`` additionally extracts text from designated release PDFs and scans
it. Exit is nonzero on any violation.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REGISTRY = ROOT / "claims" / "axiom_registry.yaml"
# Active surface roots scanned for stale-basis tokens.
ACTIVE_GLOBS = [
    "README.md",
    "README_FR.md",
    "paper/*.tex",
    "paper/tex_fragments/*.tex",
    "flagship/*.tex",
    "extra/compact_proof_of_oph.tex",
    "cosmology/*.tex",
    "book/*.md",
    "docs/*.md",
    "claims/*.yaml",
    "claims/*.md",
    "claims/*.json",
    "Lean/Screen/*.lean",
    "Lean/ObserverPatchHolography/**/*.lean",
    "assets/prediction-chain.svg",
    "assets/book_diagrams/*.svg",
]

# Path substrings excluded as archives, provenance, or non-active records.
ALLOWLIST_PATHS = [
    "docs/THEORY_SYNC_AUDIT_",
    "docs/AXIOM_REFERENCE.md",       # defines the retired principles
    "docs/STYLE_GUIDE.md",
    "claims/axiom_registry.yaml",    # names the retired principles
    "claims/frozen_prediction_register.json",  # frozen custody bytes
    "docs/FROZEN_PREDICTION_LADDER.md",        # rendered from frozen rows
    "tools/test_check_axiom_consistency",
]

# Regexes whose ACTIVE use is a violation. Case-insensitive where noted.
STALE_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("five-axiom count", re.compile(r"\bfive\s+axioms\b", re.IGNORECASE)),
    ("five-axiom count (fr)", re.compile(r"\bcinq\s+axiomes\b", re.IGNORECASE)),
    ("five-axiom basis", re.compile(r"\bfive-axiom\b", re.IGNORECASE)),
    ("A1--A5 range", re.compile(r"\bA1\s*(?:--|-|–|—)\s*A5\b")),
    ("OPH5 umbrella", re.compile(r"\bOPH5\b")),
    ("Axiom 4 as principle", re.compile(r"\bAxiom[~ ]*4\b")),
    ("Axiom 5 as principle", re.compile(r"\bAxiom[~ ]*5\b")),
    ("axiom five", re.compile(r"\baxiom five\b", re.IGNORECASE)),
    ("economy axiom", re.compile(r"\beconomy axiom\b", re.IGNORECASE)),
    (
        "retired recovery principle",
        re.compile(r"\bRecoverable Generalized Entropy axiom\b", re.IGNORECASE),
    ),
    (
        "retired A3 refinement clause",
        re.compile(r"\brefinement[- ]closure clause of (?:the third OPH axiom|Axiom 3)\b", re.IGNORECASE),
    ),
    (
        "retired combined A3 branch",
        re.compile(r"\bAxiom-?3 MaxEnt/refinement branch\b", re.IGNORECASE),
    ),
    (
        "retired economy selector",
        re.compile(r"\bMAR is used\b", re.IGNORECASE),
    ),
    (
        "retired selector name",
        re.compile(r"\bminimal\s+admissib(?:le|ility)\b", re.IGNORECASE),
    ),
    (
        "retired selector history",
        re.compile(
            r"\beconomy\s+(?:rule|selector|axiom)\b|"
            r"\bformer\s+least-value\s+selection\b|"
            r"\bwithdrawn\s+from\s+its\s+declared\s+branch\b",
            re.IGNORECASE,
        ),
    ),
    (
        "MAR as core principle",
        re.compile(r"Minimal Admissible Realization economy rule as an axiom"),
    ),
]

# Line-level allowances: a stale token on a line matching one of these is not
# a violation (bibliography keys, group names, kernel-axiom audits, quoted
# retired terminology inside migration/status prose).
LINE_ALLOW = [
    re.compile(r"\\cite\{oph5\}|\\bibitem\{oph5\}"),
    re.compile(r"alternatingGroup|A5OPH|a5_closure|A5_WZ_PHYSICAL_STATUS"),
    re.compile(r"#print axioms|sorryAx|propext|Quot\.sound|Classical\.choice"),
    re.compile(r"former_A4|former_A5|retired_principles|retired as|withdrawn"),
    re.compile(r"no longer|formerly|superseded"),  # explicit supersession notes
    re.compile(r'^\|\s*Five axioms \(also "OPH5"\)\s*\|'),
    re.compile(r"^\|\s*Recoverable generalized entropy axiom \(the recovery bundle\)\s*\|"),
]

ENTRY_SURFACES = {
    "README.md": "three core axioms",
    "README_FR.md": "trois axiomes",
    "docs/AXIOM_REFERENCE.md": "exactly three core axioms",
}

REQUIRED_ENUM = [
    "axiom_forced",
    "exact_named_realization",
    "discovery_only",
    "conditional_open_interface",
    "independence_limited",
    "physical_identification",
    "withdrawn",
]

PDF_SURFACES = [
    "flagship/from_observer_consensus_to_standard_physics.pdf",
    "extra/compact_proof_of_oph.pdf",
    "paper/observers_are_all_you_need.pdf",
    "paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf",
    "paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf",
    "paper/deriving_the_particle_zoo_from_observer_consistency.pdf",
    "paper/reality_as_consensus_protocol.pdf",
    "paper/screen_microphysics_and_observer_synchronization.pdf",
    "book/reverse-engineering-reality-book.pdf",
]


def _strip_comments(text: str) -> str:
    lines = [line for line in text.splitlines() if not line.lstrip().startswith("#")]
    return "\n".join(lines)


def registry_checks(errors: list[str]) -> dict:
    registry = json.loads(_strip_comments(REGISTRY.read_text(encoding="utf-8")))
    if registry.get("core_axiom_count") != 3:
        errors.append("registry: core_axiom_count must be exactly 3")
    axioms = registry.get("axioms", [])
    if [a.get("id") for a in axioms] != ["A1", "A2", "A3"]:
        errors.append("registry: axioms must be exactly A1, A2, A3 in order")
    for axiom in axioms:
        for field in ("key", "informal", "formal_concise", "reference_anchor",
                      "constrains", "does_not_imply"):
            if not axiom.get(field):
                errors.append(f"registry: axiom {axiom.get('id')} missing {field}")
    if len(axioms) == 3:
        a1 = axioms[0].get("formal_concise", "")
        for needle in (
            "finite-dimensional complex Hilbert response space",
            "faithful real-linear map",
            "commutator-closed",
            "complete for the declared quotient-visible infinitesimal port response",
            "positive definite",
        ):
            if needle not in a1:
                errors.append(f"registry: A1 formal_concise missing response contract: {needle}")
        a2 = axioms[1].get("formal_concise", "")
        for needle in (
            "proper carrier automorphism",
            "one projective implementer",
            "PU(H)",
            "endogenous",
            "unitary centralizer",
        ):
            if needle not in a2:
                errors.append(f"registry: A2 formal_concise missing transport contract: {needle}")
        for forbidden in ("2I", "A5", "SU(2)", "SU(3)", "su(2)", "su(3)"):
            if forbidden in a1 or forbidden in a2:
                errors.append(
                    "registry: A1/A2 formal clauses contain a forbidden "
                    f"target group name: {forbidden}"
                )
        a3 = axioms[2].get("formal_concise", "")
        for needle in (
            "restriction map is injective",
            "strictly positive exact weights",
            "Umegaki relative entropies",
            "identity-proportional",
        ):
            if needle not in a3:
                errors.append(f"registry: A3 formal_concise missing exact contract: {needle}")
    if registry.get("status_enum") != REQUIRED_ENUM:
        errors.append("registry: status_enum does not match the canonical seven values")
    for row in registry.get("premise_interfaces", []):
        if row.get("class") not in REQUIRED_ENUM:
            errors.append(f"registry: interface {row.get('id')} has invalid class")
        if row.get("class") == "physical_identification" and row.get(
            "attachment_state"
        ) not in ("pending", "supported", "rejected"):
            errors.append(
                f"registry: physical identification {row.get('id')} needs attachment_state"
            )
    return registry


def path_allowed(path: Path) -> bool:
    rel = path.relative_to(ROOT).as_posix()
    return any(token in rel for token in ALLOWLIST_PATHS)


def scan_text(rel: str, text: str, errors: list[str]) -> None:
    for number, line in enumerate(text.splitlines(), start=1):
        if any(pattern.search(line) for pattern in LINE_ALLOW):
            continue
        for label, pattern in STALE_PATTERNS:
            if pattern.search(line):
                errors.append(f"{rel}:{number}: stale token ({label}): {line.strip()[:120]}")
    for label, pattern in STALE_PATTERNS:
        for match in pattern.finditer(text):
            if "\n" in match.group(0):
                errors.append(
                    f"{rel}: cross-line stale token ({label}): "
                    f"{' '.join(match.group(0).split())[:120]}"
                )


def scan_surfaces(errors: list[str]) -> None:
    seen: set[Path] = set()
    for glob in ACTIVE_GLOBS:
        for path in sorted(ROOT.glob(glob)):
            if not path.is_file() or path in seen or path_allowed(path):
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8", errors="replace")
            scan_text(path.relative_to(ROOT).as_posix(), text, errors)


def entry_surface_checks(errors: list[str]) -> None:
    for rel, needle in ENTRY_SURFACES.items():
        path = ROOT / rel
        if not path.is_file():
            errors.append(f"{rel}: designated entry surface missing")
            continue
        text = path.read_text(encoding="utf-8", errors="replace").lower()
        if needle.lower() not in text:
            errors.append(f"{rel}: entry surface does not state the three-axiom basis ('{needle}')")


def pdf_checks(errors: list[str]) -> None:
    try:
        from pypdf import PdfReader
    except ImportError:
        errors.append("--pdf requested but pypdf is unavailable")
        return
    for rel in PDF_SURFACES:
        path = ROOT / rel
        if not path.is_file():
            errors.append(f"{rel}: designated PDF missing")
            continue
        try:
            text = "\n".join((page.extract_text() or "") for page in PdfReader(str(path)).pages)
            text = text.translate(str.maketrans({"ﬁ": "fi", "ﬂ": "fl", "ﬀ": "ff"}))
        except Exception as exc:  # noqa: BLE001
            errors.append(f"{rel}: PDF extraction failed: {exc}")
            continue
        for label, pattern in STALE_PATTERNS:
            if pattern.search(text):
                errors.append(f"{rel}: extracted text carries stale token ({label})")


def write_inventory() -> None:
    rows = []
    seen: set[Path] = set()
    for glob in ACTIVE_GLOBS:
        for path in sorted(ROOT.glob(glob)):
            if not path.is_file() or path in seen:
                continue
            seen.add(path)
            rows.append(
                {
                    "path": path.relative_to(ROOT).as_posix(),
                    "glob": glob,
                    "allowlisted": path_allowed(path),
                }
            )
    payload = {
        "schema": "oph.active_surface_inventory.v1",
        "generator": "tools/check_axiom_consistency.py --inventory",
        "surface_count": len(rows),
        "entry_surfaces": sorted(ENTRY_SURFACES),
        "pdf_surfaces": PDF_SURFACES,
        "surfaces": rows,
    }
    out = ROOT / "claims" / "active_surface_inventory.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out.relative_to(ROOT)} ({len(rows)} surfaces)")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--pdf", action="store_true", help="also scan designated release PDFs")
    parser.add_argument("--list-only", action="store_true", help="print violations without failing")
    parser.add_argument("--inventory", action="store_true", help="write claims/active_surface_inventory.json")
    args = parser.parse_args(argv)

    if args.inventory:
        write_inventory()
        return 0

    errors: list[str] = []
    registry_checks(errors)
    scan_surfaces(errors)
    entry_surface_checks(errors)
    if args.pdf:
        pdf_checks(errors)

    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        print(f"axiom consistency: {len(errors)} violation(s)", file=sys.stderr)
        return 0 if args.list_only else 1
    print("axiom consistency OK: three-axiom basis is coherent across active surfaces")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
