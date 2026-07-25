#!/usr/bin/env python3
"""Shared validator and renderer for public quantitative claim blocks.

The policy is intentionally narrower than a raw numeral scan. Dates, issue
numbers, theorem counts, paths, representation dimensions, and equations are
not automatically physical claims. Public tables that compare an OPH value
with an external value are governed: they must be generated from the canonical
manifest, whose rows resolve to a registry claim, registry class, emitting
script, artifact pointer, and any external comparison fixture.
"""

from __future__ import annotations

import json
import re
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_RELATIVE = Path("claims/public_surface_quantitative_claims.json")
REGISTRY_RELATIVE = Path("claims/claim_registry.yaml")

BLOCK_START = "<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->"
BLOCK_END = "<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->"

CLASS_VOCABULARY = {
    "physical_establishment",
    "empirical_implementation",
    "emitted_artifact",
    "branch_entry",
    "conditional_implication",
    "declared_structure",
}
OPH_VALUE_CLASSES = {"physical_establishment", "empirical_implementation"}
ROLES = {
    "diagnostic",
    "conditional_branch",
    "rejected_candidate",
    "target_anchored_backsolve",
    "structural_non_discriminating",
    "oph_result",
}

PUBLIC_SURFACE_GLOBS = ("README.md", "README_FR.md", "docs/**/*.md", "extra/*.md")
EXTERNAL_HEADER_TERMS = re.compile(
    r"\b(?:PDG|NIST|CODATA|external|measurement|measured|experimental|reference"
    r"|mesure|mesurée|expérimental|référence|comparaison)\b",
    re.IGNORECASE,
)
OPH_HEADER_TERM = re.compile(r"(?:^|[\s|])OPH(?:[\s|]|$)", re.IGNORECASE)
NUMERIC_TOKEN = re.compile(
    r"(?<![A-Za-z0-9_])[-+]?(?:\d+(?:[.,]\d*)?|[.,]\d+)(?:[eE][-+]?\d+)?"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def pointer_get(document: Any, pointer: str) -> Any:
    """Resolve one RFC-6901 JSON pointer."""
    if pointer == "":
        return document
    if not isinstance(pointer, str) or not pointer.startswith("/"):
        raise ValueError(f"invalid JSON pointer {pointer!r}")
    value = document
    for raw_token in pointer[1:].split("/"):
        token = raw_token.replace("~1", "/").replace("~0", "~")
        if isinstance(value, list):
            try:
                value = value[int(token)]
            except (ValueError, IndexError) as exc:
                raise KeyError(pointer) from exc
        elif isinstance(value, dict) and token in value:
            value = value[token]
        else:
            raise KeyError(pointer)
    return value


def as_decimal(value: Any, label: str) -> Decimal:
    if isinstance(value, bool) or value is None:
        raise ValueError(f"{label} is not numeric: {value!r}")
    try:
        number = Decimal(str(value))
    except (InvalidOperation, ValueError) as exc:
        raise ValueError(f"{label} is not numeric: {value!r}") from exc
    if not number.is_finite():
        raise ValueError(f"{label} is not finite: {value!r}")
    return number


def _load_artifact(
    root: Path,
    relative: str,
    cache: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    if relative not in cache:
        cache[relative] = load_json(root / relative)
    return cache[relative]


def _registry_by_id(root: Path) -> tuple[dict[str, dict[str, Any]], int]:
    registry = load_json(root / REGISTRY_RELATIVE)
    claims = registry.get("claims", [])
    by_id = {
        claim["claim_id"]: claim
        for claim in claims
        if isinstance(claim, dict) and isinstance(claim.get("claim_id"), str)
    }
    physical_count = sum(
        claim.get("claim_class") == "physical_establishment" for claim in claims
    )
    return by_id, physical_count


def _validate_localized_text(
    value: Any,
    label: str,
    issues: list[str],
) -> None:
    if not isinstance(value, dict) or set(value) != {"en", "fr"}:
        issues.append(f"{label}: expected exact en/fr translations")
        return
    if not all(isinstance(value[locale], str) and value[locale] for locale in ("en", "fr")):
        issues.append(f"{label}: translations must be nonempty strings")


def validate_manifest(
    root: Path = ROOT,
    manifest: dict[str, Any] | None = None,
) -> tuple[list[str], dict[str, Any], dict[str, dict[str, Any]], int]:
    """Validate manifest-to-registry and manifest-to-producer resolution."""
    issues: list[str] = []
    manifest = manifest or load_json(root / MANIFEST_RELATIVE)
    registry, physical_count = _registry_by_id(root)
    cache: dict[str, dict[str, Any]] = {}

    if manifest.get("schema") != "oph.public_surface_quantitative_claims.v1":
        issues.append("manifest: unsupported or missing schema")

    surfaces = manifest.get("surfaces")
    expected_surfaces = {
        ("README.md", "en"),
        ("README_FR.md", "fr"),
    }
    actual_surfaces: set[tuple[str, str]] = set()
    if not isinstance(surfaces, list):
        issues.append("manifest: surfaces must be a list")
        surfaces = []
    for surface in surfaces:
        if not isinstance(surface, dict):
            issues.append("manifest: surface rows must be objects")
            continue
        path = surface.get("path")
        locale = surface.get("locale")
        if isinstance(path, str) and isinstance(locale, str):
            actual_surfaces.add((path, locale))
        if not isinstance(path, str) or not (root / path).is_file():
            issues.append(f"manifest: public surface does not exist: {path!r}")
    if actual_surfaces != expected_surfaces:
        issues.append(
            "manifest: surfaces must be exactly README.md/en and README_FR.md/fr"
        )

    rows = manifest.get("rows")
    if not isinstance(rows, list) or not rows:
        issues.append("manifest: rows must be a nonempty list")
        rows = []
    row_ids = [row.get("row_id") for row in rows if isinstance(row, dict)]
    duplicates = sorted(
        row_id
        for row_id in set(row_ids)
        if isinstance(row_id, str) and row_ids.count(row_id) > 1
    )
    if duplicates:
        issues.append(f"manifest: duplicate row IDs: {duplicates}")

    required_roles = {
        "diagnostic",
        "rejected_candidate",
        "target_anchored_backsolve",
        "structural_non_discriminating",
    }
    seen_roles: set[str] = set()

    for index, row in enumerate(rows):
        label = f"manifest row {index}"
        if not isinstance(row, dict):
            issues.append(f"{label}: row must be an object")
            continue
        row_id = row.get("row_id")
        if not isinstance(row_id, str) or not row_id:
            issues.append(f"{label}: row_id must be a nonempty string")
            row_id = f"index-{index}"
        label = f"row {row_id}"

        claim_id = row.get("claim_id")
        claim_class = row.get("claim_class")
        claim = registry.get(claim_id)
        if claim is None:
            issues.append(f"{label}: unknown registry claim ID {claim_id!r}")
        elif claim_class != claim.get("claim_class"):
            issues.append(
                f"{label}: declared class {claim_class!r} does not match registry "
                f"class {claim.get('claim_class')!r}"
            )
        if claim_class not in CLASS_VOCABULARY:
            issues.append(f"{label}: class {claim_class!r} is outside scoreboard vocabulary")

        role = row.get("role")
        if role not in ROLES:
            issues.append(f"{label}: unsupported role {role!r}")
        else:
            seen_roles.add(role)
        if role == "oph_result" and claim_class not in OPH_VALUE_CLASSES:
            issues.append(
                f"{label}: an OPH result value requires physical_establishment or "
                "empirical_implementation class"
            )

        _validate_localized_text(row.get("label"), f"{label} label", issues)
        _validate_localized_text(row.get("status"), f"{label} status", issues)
        if role == "target_anchored_backsolve":
            status = row.get("status", {})
            if (
                not isinstance(status, dict)
                or "never a prediction" not in status.get("en", "").lower()
                or "jamais une prédiction" not in status.get("fr", "").lower()
            ):
                issues.append(
                    f"{label}: target-anchored back-solves must say they are never "
                    "a prediction in both locales"
                )
        if role == "rejected_candidate":
            status = row.get("status", {})
            if (
                not isinstance(status, dict)
                or "rejected" not in status.get("en", "").lower()
                or "rejet" not in status.get("fr", "").lower()
            ):
                issues.append(
                    f"{label}: rejected candidates must state the rejection "
                    "disposition in both locales"
                )

        producer = row.get("producer")
        if not isinstance(producer, dict):
            issues.append(f"{label}: producer must be an object")
            continue
        script_rel = producer.get("script")
        artifact_rel = producer.get("artifact")
        if not isinstance(script_rel, str) or not script_rel.endswith(".py"):
            issues.append(f"{label}: producer script must be a .py path")
        elif not (root / script_rel).is_file():
            issues.append(f"{label}: emitting script does not exist: {script_rel}")
        if not isinstance(artifact_rel, str) or not (root / artifact_rel).is_file():
            issues.append(f"{label}: producer artifact does not exist: {artifact_rel!r}")
            continue
        try:
            artifact = _load_artifact(root, artifact_rel, cache)
        except (OSError, json.JSONDecodeError) as exc:
            issues.append(f"{label}: cannot load producer artifact: {exc}")
            continue
        artifact_id = producer.get("artifact_id")
        if artifact.get("artifact") != artifact_id:
            issues.append(
                f"{label}: artifact ID {artifact.get('artifact')!r} does not match "
                f"declared {artifact_id!r}"
            )
        if isinstance(script_rel, str) and (root / script_rel).is_file():
            script_text = (root / script_rel).read_text(encoding="utf-8", errors="ignore")
            if isinstance(artifact_id, str) and artifact_id not in script_text:
                issues.append(
                    f"{label}: emitting script does not name artifact ID {artifact_id!r}"
                )
        try:
            producer_value = as_decimal(
                pointer_get(artifact, producer.get("value_pointer")),
                f"{label} producer value",
            )
        except (KeyError, TypeError, ValueError) as exc:
            issues.append(f"{label}: unresolved producer value: {exc}")
            producer_value = None

        guards = producer.get("guards", {})
        if not isinstance(guards, dict) or not guards:
            issues.append(f"{label}: producer guards must be a nonempty object")
        else:
            for pointer, expected in guards.items():
                try:
                    actual = pointer_get(artifact, pointer)
                except (KeyError, TypeError, ValueError):
                    issues.append(f"{label}: unresolved producer guard {pointer!r}")
                    continue
                if actual != expected:
                    issues.append(
                        f"{label}: producer guard {pointer} is {actual!r}, "
                        f"expected {expected!r}"
                    )

        reference = row.get("reference")
        reference_value: Decimal | None = None
        uncertainty: Decimal | None = None
        if reference is not None:
            if not isinstance(reference, dict):
                issues.append(f"{label}: reference must be an object")
            else:
                reference_rel = reference.get("artifact")
                if not isinstance(reference_rel, str) or not (root / reference_rel).is_file():
                    issues.append(
                        f"{label}: reference artifact does not exist: {reference_rel!r}"
                    )
                else:
                    try:
                        reference_artifact = _load_artifact(root, reference_rel, cache)
                        reference_value = as_decimal(
                            pointer_get(reference_artifact, reference.get("value_pointer")),
                            f"{label} reference value",
                        )
                        if "uncertainty_pointer" in reference:
                            uncertainty = as_decimal(
                                pointer_get(
                                    reference_artifact,
                                    reference["uncertainty_pointer"],
                                ),
                                f"{label} reference uncertainty",
                            )
                            if uncertainty <= 0:
                                issues.append(
                                    f"{label}: reference uncertainty must be positive"
                                )
                    except (
                        OSError,
                        json.JSONDecodeError,
                        KeyError,
                        TypeError,
                        ValueError,
                    ) as exc:
                        issues.append(f"{label}: unresolved reference: {exc}")
                _validate_localized_text(reference.get("label"), f"{label} reference label", issues)

        supporting = row.get("supporting_artifacts", [])
        if not isinstance(supporting, list):
            issues.append(f"{label}: supporting_artifacts must be a list")
        else:
            for support_index, support in enumerate(supporting):
                support_label = f"{label} supporting artifact {support_index}"
                if not isinstance(support, dict):
                    issues.append(f"{support_label}: entry must be an object")
                    continue
                support_rel = support.get("artifact")
                if not isinstance(support_rel, str) or not (root / support_rel).is_file():
                    issues.append(
                        f"{support_label}: artifact does not exist: {support_rel!r}"
                    )
                    continue
                try:
                    support_artifact = _load_artifact(root, support_rel, cache)
                except (OSError, json.JSONDecodeError) as exc:
                    issues.append(f"{support_label}: cannot load artifact: {exc}")
                    continue
                support_id = support.get("artifact_id")
                if support_artifact.get("artifact") != support_id:
                    issues.append(
                        f"{support_label}: artifact ID "
                        f"{support_artifact.get('artifact')!r} does not match "
                        f"declared {support_id!r}"
                    )
                for pointer, expected in support.get("guards", {}).items():
                    try:
                        actual = pointer_get(support_artifact, pointer)
                    except (KeyError, TypeError, ValueError):
                        issues.append(
                            f"{support_label}: unresolved guard {pointer!r}"
                        )
                        continue
                    if actual != expected:
                        issues.append(
                            f"{support_label}: guard {pointer} is {actual!r}, "
                            f"expected {expected!r}"
                        )
                for pointer, fragment in support.get("contains", {}).items():
                    try:
                        actual = pointer_get(support_artifact, pointer)
                    except (KeyError, TypeError, ValueError):
                        issues.append(
                            f"{support_label}: unresolved contains guard {pointer!r}"
                        )
                        continue
                    if not isinstance(actual, str) or not isinstance(fragment, str):
                        issues.append(
                            f"{support_label}: contains guards require strings"
                        )
                    elif fragment not in actual:
                        issues.append(
                            f"{support_label}: {pointer} does not contain "
                            f"{fragment!r}"
                        )

        if producer_value is not None and reference_value is not None:
            display = row.get("display", {})
            scale = as_decimal(display.get("scale", "1"), f"{label} display scale")
            value_display = producer_value * scale
            reference_display = reference_value * scale
            uncertainty_display = uncertainty * scale if uncertainty is not None else None
            if uncertainty_display is not None:
                indistinguishable = abs(value_display - reference_display) <= uncertainty_display
            else:
                decimals = int(display.get("reference_decimals", 0))
                rounding_cell = Decimal(1).scaleb(-decimals)
                indistinguishable = (
                    abs(value_display - reference_display) < rounding_cell / Decimal(2)
                )
            if indistinguishable and role != "target_anchored_backsolve":
                issues.append(
                    f"{label}: self-comparison agrees within the quoted or displayed "
                    "external precision without a target-anchored back-solve label"
                )

    missing_roles = sorted(required_roles - seen_roles)
    if missing_roles:
        issues.append(f"manifest: required quantitative roles are absent: {missing_roles}")

    return issues, manifest, registry, physical_count


def _fixed(number: Decimal, decimals: int, locale: str) -> str:
    quantum = Decimal(1).scaleb(-decimals)
    rendered = format(number.quantize(quantum, rounding=ROUND_HALF_UP), f".{decimals}f")
    if locale == "fr":
        rendered = rendered.replace(".", ",")
    return rendered


def _scientific(number: Decimal, significant: int, locale: str) -> str:
    if number == 0:
        mantissa, exponent = Decimal(0), 0
    else:
        exponent = number.copy_abs().adjusted()
        mantissa = number / (Decimal(10) ** exponent)
    decimals = max(0, significant - 1)
    mantissa_text = _fixed(mantissa, decimals, locale)
    return f"${mantissa_text}\\times10^{{{exponent}}}$"


def _parenthetical_uncertainty(
    value: Decimal,
    uncertainty: Decimal,
    decimals: int,
    locale: str,
) -> str:
    value_text = _fixed(value, decimals, locale)
    digits = int((uncertainty * (Decimal(10) ** decimals)).to_integral_value())
    return f"{value_text}({digits})"


def _row_values(
    root: Path,
    row: dict[str, Any],
    cache: dict[str, dict[str, Any]],
) -> tuple[Decimal, Decimal | None, Decimal | None]:
    producer = row["producer"]
    artifact = _load_artifact(root, producer["artifact"], cache)
    value = as_decimal(
        pointer_get(artifact, producer["value_pointer"]),
        f"{row['row_id']} producer value",
    )
    reference = row.get("reference")
    reference_value = None
    uncertainty = None
    if reference:
        ref_artifact = _load_artifact(root, reference["artifact"], cache)
        reference_value = as_decimal(
            pointer_get(ref_artifact, reference["value_pointer"]),
            f"{row['row_id']} reference value",
        )
        if "uncertainty_pointer" in reference:
            uncertainty = as_decimal(
                pointer_get(ref_artifact, reference["uncertainty_pointer"]),
                f"{row['row_id']} reference uncertainty",
            )
    return value, reference_value, uncertainty


def render_section(
    root: Path,
    manifest: dict[str, Any],
    registry: dict[str, dict[str, Any]],
    physical_count: int,
    locale: str,
) -> str:
    cache: dict[str, dict[str, Any]] = {}
    if locale == "en":
        heading = "## Quantitative Claim Status"
        intro = (
            f"The canonical claim registry contains **{physical_count} claims of class "
            "`physical_establishment`**. The rows below are generated from "
            "[machine-readable claim annotations]"
            "(claims/public_surface_quantitative_claims.json). Rejected candidates, "
            "diagnostics, target-anchored back-solves, and structural zeros retain "
            "their separate roles."
        )
        headers = (
            "Quantity",
            "Registered branch value",
            "External comparison",
            "Status",
        )
        no_comparison = "none"
        class_word = "class"
        footer = (
            "`m_u` and `m_c` are absent because the Clebsch lane does not emit "
            "them. `c` is absent because its SI value is definitional. The "
            "geometric $G$ identity is shown, but an SI value is omitted because "
            "the operational clock packet is not source-closed."
        )
    else:
        heading = "## Statut des énoncés quantitatifs"
        intro = (
            f"Le registre canonique contient **{physical_count} énoncé de classe "
            "`physical_establishment`**. Les lignes ci-dessous sont produites à "
            "partir d’[annotations lisibles par machine]"
            "(claims/public_surface_quantitative_claims.json). Les candidats "
            "rejetés, les diagnostics, les résolutions à rebours ancrées sur une "
            "cible et les zéros structurels conservent des rôles distincts."
        )
        headers = (
            "Quantité",
            "Valeur de la branche enregistrée",
            "Comparaison externe",
            "Statut",
        )
        no_comparison = "aucune"
        class_word = "classe"
        footer = (
            "`m_u` et `m_c` sont absents, car la branche de Clebsch ne les émet "
            "pas. `c` est absent, car sa valeur dans le Système international est "
            "une définition. L’identité géométrique pour $G$ est affichée, mais "
            "aucune valeur dans le SI ne l’est, car le paquet d’horloge "
            "opérationnel n’est pas clos à partir de la source."
        )

    lines = [
        BLOCK_START,
        "<!-- Generated by tools/build_public_quantitative_section.py; do not edit. -->",
        heading,
        "",
        intro,
        "",
        f"| {headers[0]} | {headers[1]} | {headers[2]} | {headers[3]} |",
        "| --- | ---: | ---: | --- |",
    ]

    for row in manifest["rows"]:
        value, reference, uncertainty = _row_values(root, row, cache)
        display = row["display"]
        scale = as_decimal(display.get("scale", "1"), f"{row['row_id']} scale")
        value_scaled = value * scale
        value_text = _fixed(value_scaled, int(display["value_decimals"]), locale)
        unit = display.get("unit", "")
        value_cell = f"`{value_text}{f' {unit}' if unit else ''}`"

        if reference is None:
            reference_cell = no_comparison
        else:
            reference_scaled = reference * scale
            style = display.get("reference_style", "plain")
            if style == "parenthetical_uncertainty":
                if uncertainty is None:
                    raise ValueError(
                        f"{row['row_id']}: parenthetical uncertainty is missing"
                    )
                reference_text = _parenthetical_uncertainty(
                    reference_scaled,
                    uncertainty * scale,
                    int(display["reference_decimals"]),
                    locale,
                )
            elif style == "plus_minus":
                if uncertainty is None:
                    raise ValueError(f"{row['row_id']}: plus-minus uncertainty is missing")
                reference_text = (
                    _fixed(reference_scaled, int(display["reference_decimals"]), locale)
                    + " ± "
                    + _fixed(
                        uncertainty * scale,
                        int(display["uncertainty_decimals"]),
                        locale,
                    )
                )
            elif style == "plain":
                reference_text = _fixed(
                    reference_scaled,
                    int(display["reference_decimals"]),
                    locale,
                )
            else:
                raise ValueError(f"{row['row_id']}: unknown reference style {style!r}")
            reference_label = row["reference"]["label"][locale]
            reference_cell = (
                f"`{reference_text}{f' {unit}' if unit else ''}` "
                f"({reference_label})"
            )

        substitutions: dict[str, str] = {}
        if reference is not None:
            relative_gap = abs(value - reference) / abs(reference)
            substitutions["percent_gap"] = (
                _fixed(relative_gap * Decimal(100), 1, locale) + "%"
            )
            substitutions["relative_gap_scientific"] = _scientific(
                relative_gap, 2, locale
            )
            if uncertainty is not None:
                substitutions["sigma_scientific"] = _scientific(
                    abs(value - reference) / uncertainty,
                    2,
                    locale,
                )
        status = row["status"][locale].format(**substitutions)
        claim_class = registry[row["claim_id"]]["claim_class"]
        status = f"{status} ({class_word}: `{claim_class}`)"
        lines.append(
            f"| {row['label'][locale]} | {value_cell} | {reference_cell} | {status} |"
        )

    lines.extend(["", footer, BLOCK_END])
    return "\n".join(lines)


def replace_generated_block(text: str, block: str, path: str) -> str:
    if text.count(BLOCK_START) != 1 or text.count(BLOCK_END) != 1:
        raise ValueError(
            f"{path}: expected exactly one {BLOCK_START!r}/{BLOCK_END!r} marker pair"
        )
    before, remainder = text.split(BLOCK_START, 1)
    _, after = remainder.split(BLOCK_END, 1)
    return before + block + after


def expected_surface_texts(
    root: Path = ROOT,
    manifest: dict[str, Any] | None = None,
) -> tuple[list[str], dict[str, str]]:
    issues, manifest, registry, physical_count = validate_manifest(root, manifest)
    if issues:
        return issues, {}
    outputs: dict[str, str] = {}
    for surface in manifest["surfaces"]:
        path = surface["path"]
        text = (root / path).read_text(encoding="utf-8")
        block = render_section(
            root,
            manifest,
            registry,
            physical_count,
            surface["locale"],
        )
        try:
            outputs[path] = replace_generated_block(text, block, path)
        except ValueError as exc:
            issues.append(str(exc))
    return issues, outputs


def iter_public_surfaces(root: Path) -> list[Path]:
    paths: set[Path] = set()
    for pattern in PUBLIC_SURFACE_GLOBS:
        paths.update(path for path in root.glob(pattern) if path.is_file())
    return sorted(paths)


def _table_bounds(lines: list[str], header_index: int) -> tuple[int, int]:
    start = header_index
    while start > 0 and lines[start - 1].lstrip().startswith("|"):
        start -= 1
    end = header_index
    while end + 1 < len(lines) and lines[end + 1].lstrip().startswith("|"):
        end += 1
    return start, end


def unmanaged_oph_comparison_tables(root: Path = ROOT) -> list[str]:
    """Reject hand-written OPH-versus-external numeric tables.

    The generated block is the annotation mechanism. Numerals elsewhere are
    outside this comparison-table gate, so dates, issue IDs, theorem counts,
    equations, and paths cannot become false positives.
    """
    issues: list[str] = []
    for path in iter_public_surfaces(root):
        text = path.read_text(encoding="utf-8", errors="ignore")
        block_start = text.find(BLOCK_START)
        block_end = text.find(BLOCK_END)
        lines = text.splitlines()
        offsets: list[int] = []
        offset = 0
        for line in lines:
            offsets.append(offset)
            offset += len(line) + 1
        seen_bounds: set[tuple[int, int]] = set()
        for index, line in enumerate(lines):
            if not line.lstrip().startswith("|"):
                continue
            if not OPH_HEADER_TERM.search(line) or not EXTERNAL_HEADER_TERMS.search(line):
                continue
            bounds = _table_bounds(lines, index)
            if bounds in seen_bounds:
                continue
            seen_bounds.add(bounds)
            table = "\n".join(lines[bounds[0] : bounds[1] + 1])
            if not NUMERIC_TOKEN.search(table):
                continue
            table_offset = offsets[bounds[0]]
            inside_generated = (
                block_start >= 0
                and block_end >= 0
                and block_start <= table_offset <= block_end
            )
            if not inside_generated:
                relative = path.relative_to(root)
                issues.append(
                    f"{relative}:{index + 1}: unmanaged numeric OPH-versus-external "
                    "comparison table; add a registry-resolved row to "
                    f"{MANIFEST_RELATIVE}"
                )
    return issues


def check_repository(root: Path = ROOT) -> list[str]:
    issues, outputs = expected_surface_texts(root)
    if not issues:
        for relative, expected in outputs.items():
            actual = (root / relative).read_text(encoding="utf-8")
            if actual != expected:
                issues.append(
                    f"{relative}: generated quantitative claim block is stale; run "
                    "python tools/build_public_quantitative_section.py"
                )
    issues.extend(unmanaged_oph_comparison_tables(root))
    return issues
