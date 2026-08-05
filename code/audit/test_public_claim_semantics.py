"""Regression gates for high-risk public epistemic-status wording."""

from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def _highlight_blocks(text: str, heading: str, terminator: str) -> list[str]:
    start = text.index(heading)
    body = text[start : text.index(terminator, start)]
    blocks: list[str] = []
    for number in range(1, 8):
        block_start = body.index(f"{number}. **")
        block_end = (
            body.index(f"{number + 1}. **", block_start)
            if number < 7
            else len(body)
        )
        blocks.append(body[block_start:block_end])
    return blocks


def test_fz10_public_wording_matches_the_frozen_decision_rule() -> None:
    register = json.loads(
        (ROOT / "claims/frozen_prediction_register.json").read_text(encoding="utf-8")
    )
    row = next(item for item in register["rows"] if item["id"] == "FZ-10")
    rule = row["kill_band"].lower()
    assert "more than three standard uncertainties" in rule
    assert "sigma <= 0.045 mev" in rule
    assert "within two standard uncertainties" in rule

    surfaces = {
        "compact": ROOT / "extra/compact_proof_of_oph.tex",
        "flagship": ROOT / "flagship/from_observer_consensus_to_standard_physics.tex",
        "particle paper": ROOT / "paper/deriving_the_particle_zoo_from_observer_consistency.tex",
        "Koide paper": ROOT / "extra/koide_identity_from_positive_c3_face_circulants.tex",
        "standard-model chapter": ROOT / "book/chapter-14-standard-model.md",
        "matter chapter": ROOT / "book/chapter-16-matter.md",
        "book synthesis": ROOT / "book/chapter-19-synthesis.md",
        "English README": ROOT / "README.md",
        "French README": ROOT / "README_FR.md",
        "English glossary": ROOT / "book/appendix-concept-glossary.md",
    }
    forbidden = re.compile(
        r"(?:outside|lands? outside).{0,40}(?:72.{0,5}ev|window|enclosure)"
        r".{0,40}(?:refutes?|kills?)",
        flags=re.IGNORECASE | re.DOTALL,
    )
    for name, path in surfaces.items():
        text = path.read_text(encoding="utf-8")
        assert not forbidden.search(text), f"{name} restores the obsolete window kill rule"

    compact = surfaces["compact"].read_text(encoding="utf-8").lower()
    readme = surfaces["English README"].read_text(encoding="utf-8").lower()
    assert "more than three" in compact and "standard uncertainties" in compact
    assert re.search(r"target-informed\s+conditional\s+postdiction", readme)


def test_readme_highlights_remain_short_and_reader_facing() -> None:
    surfaces = [
        (
            ROOT / "README.md",
            "## Seven Reproducible Physics Receipts",
            "A separate signed-graph theorem",
        ),
        (
            ROOT / "README_FR.md",
            "## Sept reçus de physique reproductibles",
            "Un théorème distinct de graphe signé",
        ),
    ]
    for path, heading, terminator in surfaces:
        text = path.read_text(encoding="utf-8")
        for number, block in enumerate(
            _highlight_blocks(text, heading, terminator), start=1
        ):
            prose = re.sub(r"\[([^]]+)\]\([^)]+\)", r"\1", block)
            words = re.findall(r"\b[\wÀ-ÿ]+(?:[-’'][\wÀ-ÿ]+)*\b", prose)
            assert len(words) <= 100, (
                f"{path.name} receipt {number} has {len(words)} words; "
                "proof inventories belong in the linked technical surfaces"
            )
            assert len(re.findall(r"\]\(", block)) <= 3, (
                f"{path.name} receipt {number} links more than three destinations"
            )

    english = (ROOT / "README.md").read_text(encoding="utf-8")
    french = (ROOT / "README_FR.md").read_text(encoding="utf-8")
    assert "**The four laws of thermodynamics" in english
    assert "conditional finite theorem package" not in english.lower()
    assert "72-eV-wide interval" in english
    assert "lands within 72 eV" not in english
    assert "intervalle large de 72 eV" in french
    assert "tombe à 72 eV" not in french
    assert re.search(r"separately\s+specified matter structure", english)
    assert "Two declared twelve-port wave" in english


def test_active_v2_owner_display_separates_historical_milestones() -> None:
    ladder = (ROOT / "docs/FROZEN_PREDICTION_LADDER.md").read_text(encoding="utf-8")
    assert "historical milestone C4" in ladder
    assert "historical milestone C3-V" in ladder
    assert not re.search(r"active V2 [^|]+ \(C(?:3|4)", ladder)


def test_screen_paper_separates_same_radius_and_two_radius_controls() -> None:
    paper = (
        ROOT / "paper/screen_microphysics_and_observer_synchronization.tex"
    ).read_text(encoding="utf-8")
    normalized = " ".join(paper.split())
    assert "same-radius positive-weight member" in normalized
    assert "isotropic through \\(k^8\\)" in normalized
    assert "face radius \\(2\\) and weight \\(27/1600\\)" in normalized
    assert "raw eighth-moment \\(\\mathcal I_6\\) multiple \\(-256/125\\)" in normalized
    stale = re.compile(
        r"same-radius.{0,500}leading anisotropy.{0,100}k\^8.{0,100}-256/125",
        flags=re.IGNORECASE,
    )
    assert not stale.search(normalized)
