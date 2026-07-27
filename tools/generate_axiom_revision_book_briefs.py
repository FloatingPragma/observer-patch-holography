#!/usr/bin/env python3
"""Generate the internal book-preservation ledger and per-source revision briefs."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BOOK = ROOT / "book"
OUT = BOOK / "axiom_revision_briefs"
DEDICATION = "*For my wife Noon, and for Douglas Adams.*"

SOURCES = [
    BOOK / "prologue.md",
    *(BOOK / f"chapter-{index:02d}-{slug}.md" for index, slug in [
        (1, "consistency"),
        (2, "lineage"),
        (3, "screen"),
        (4, "entropy"),
        (5, "algebra"),
        (6, "overlap"),
        (7, "recovery"),
        (8, "holography"),
        (9, "entanglement"),
        (10, "error-correction"),
        (11, "maxent"),
        (12, "symmetry"),
        (13, "desitter"),
        (14, "standard-model"),
        (15, "relativity"),
        (16, "matter"),
        (17, "darwin"),
        (18, "strangeloop"),
        (19, "synthesis"),
        (20, "metaphysics"),
    ]),
    BOOK / "epilogue.md",
    BOOK / "appendix-concept-glossary.md",
    BOOK / "appendix-extended-interludes.md",
]


def sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def paragraphs(text: str) -> list[tuple[int, int, str]]:
    lines = text.splitlines()
    blocks: list[tuple[int, int, str]] = []
    start: int | None = None
    buffer: list[str] = []
    for line_no, line in enumerate(lines, start=1):
        if line.strip():
            if start is None:
                start = line_no
            buffer.append(line)
            continue
        if buffer and start is not None:
            blocks.append((start, line_no - 1, "\n".join(buffer)))
        start = None
        buffer = []
    if buffer and start is not None:
        blocks.append((start, len(lines), "\n".join(buffer)))
    return blocks


def narrative_entries(path: Path, text: str) -> list[dict[str, object]]:
    blocks = paragraphs(text)
    prose = [
        block
        for block in blocks
        if not block[2].lstrip().startswith(("#", "```", "$$", "|"))
    ]
    selected: list[tuple[str, tuple[int, int, str]]] = []
    if prose:
        selected.append(("opening", prose[0]))
    if len(prose) > 1:
        selected.append(("closing", prose[-1]))
    for block in prose:
        lower = block[2].lower()
        if DEDICATION.lower() in lower:
            selected.append(("dedication", block))
        elif any(
            marker in lower
            for marker in (
                "douglas adams",
                "hitchhiker",
                "debugg",
                "hacker",
                "joke",
                "imagine",
                "story",
                "coffee",
            )
        ):
            selected.append(("voice", block))

    unique: dict[tuple[int, int], dict[str, object]] = {}
    for kind, (start, end, body) in selected:
        key = (start, end)
        entry = unique.setdefault(
            key,
            {
                "kind": kind,
                "line_start": start,
                "line_end": end,
                "sha256": sha256(body),
                "narrative_purpose": (
                    "Protected reader-facing voice, scene, analogy, opening, or closing."
                ),
                "owner_only_deletion": kind == "dedication",
            },
        )
        if kind == "dedication":
            entry["kind"] = "dedication"
            entry["owner_only_deletion"] = True
    return list(unique.values())


def brief_text(entry: dict[str, object]) -> str:
    path = str(entry["source"])
    title = str(entry["title"])
    protected = entry["protected_passages"]
    protected_lines = "\n".join(
        f"- `{item['kind']}` lines {item['line_start']}-{item['line_end']}, "
        f"SHA-256 `{item['sha256']}`, owner-only deletion "
        f"`{str(item['owner_only_deletion']).lower()}`"
        for item in protected
    )
    if not protected_lines:
        protected_lines = "- No dedicated passage marker; preserve the full source voice and rhythm."

    return f"""# Axiom revision brief: {title}

Source: `{path}`

## Reader purpose

Preserve this source's place in the book-wide narrative while presenting the
three-axiom basis and the finite-to-physical claim boundary without internal
tracking language.

## Protected narrative material

{protected_lines}

The source-file SHA-256 at this review point is
`{entry['source_sha256']}`. Technical repairs may change the full-file hash;
the passage hashes above are the preservation checks.

## Canonical claim contract

- Core premises are exactly A1 observer-screen architecture, A2 agreement on
  accepted operational meaning, and A3 conditional maximum randomness on a
  fixed feasible family.
- Cross-scale optimizer compatibility, recovery, generalized-entropy
  stationarity, focusing, physical scale, vacuum, and economy-like choices
  are theorem targets, named interfaces, diagnostics, or withdrawn claims.
- A finite exact result keeps the scope supplied by its explicit premises.
- A physical interpretation requires its named source attachment and
  identification map.

## Retained and removed inferences

Retain exact finite mathematics, working scenes, analogies, jokes, chapter
rhythm, and externally supported physics. Remove any inference that treats
recovery as a fourth axiom, economy selection as a fifth axiom, or A3 as a
cross-scale optimizer rule. Countermodels and negatively closed hypotheses
belong in technical ledgers unless they are needed to protect the reader from
a false inference.

## Axiom use and open interfaces

Use the informal form before formal notation. Spell out the carrier
architecture when A1 is introduced. State the accepted-data diagrams for A2.
State the fixed feasible family, exact reference, state-determining cover,
positive exact weights, and relative-entropy objective when A3 is formalized.
Name every additional branch premise locally.

## Diagrams, examples, and cross-references

Check every included diagram against the three-axiom registry. Keep internal
claim identifiers, receipt names, schema fields, repository paths, and release
history out of reader-facing prose.

## Adversarial reader questions

1. Which conclusion follows from A1-A3 alone?
2. Which conclusion consumes a named branch premise or physical attachment?
3. Does any selector enter without an independently defined source rule and a
   same-architecture negative control?
4. Would the paragraph remain true for every state allowed by the exact A3
   feasible family?

## Style and signoff

- [x] Claim-registry dependency checked by the migration audit.
- [x] Preservation entries generated from the reviewed source.
- [x] Machine style guide applicable to the source.
- [ ] Owner voice and narrative signoff.
"""


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    entries: list[dict[str, object]] = []
    for source in SOURCES:
        text = source.read_text(encoding="utf-8")
        title_match = re.search(r"^#\s+(.+)$", text, re.MULTILINE)
        title = title_match.group(1) if title_match else source.stem
        entry: dict[str, object] = {
            "source": source.relative_to(ROOT).as_posix(),
            "title": title,
            "source_sha256": sha256(text),
            "protected_passages": narrative_entries(source, text),
            "claim_registry_signoff": True,
            "machine_style_signoff": True,
            "owner_voice_signoff": False,
        }
        entries.append(entry)
        (OUT / f"{source.stem}.md").write_text(brief_text(entry), encoding="utf-8")

    dedication_count = sum(
        item["kind"] == "dedication"
        for entry in entries
        for item in entry["protected_passages"]
    )
    if dedication_count != 1 or DEDICATION not in (BOOK / "prologue.md").read_text(
        encoding="utf-8"
    ):
        raise SystemExit("protected dedication is missing or duplicated")

    manifest = {
        "schema": "oph.book.axiom_revision_preservation.v1",
        "canonical_basis": ["A1", "A2", "A3"],
        "protected_dedication": DEDICATION,
        "source_count": len(entries),
        "owner_voice_signoff": False,
        "sources": entries,
    }
    (OUT / "preservation_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    ledger_lines = [
        "# Book axiom-revision preservation ledger",
        "",
        f"Protected dedication: {DEDICATION}",
        "",
        "| Source | Protected passages | Claim audit | Machine style | Owner voice |",
        "|---|---:|---|---|---|",
    ]
    for entry in entries:
        ledger_lines.append(
            f"| `{entry['source']}` | {len(entry['protected_passages'])} | pass | pass | pending |"
        )
    ledger_lines.extend(
        [
            "",
            "Owner voice approval is intentionally not inferred from machine checks. "
            "The technical migration may merge without publishing; owner approval is "
            "required before the publication gate.",
            "",
        ]
    )
    (OUT / "PRESERVATION_LEDGER.md").write_text(
        "\n".join(ledger_lines), encoding="utf-8"
    )
    print(f"generated {len(entries)} briefs and preservation ledger")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
