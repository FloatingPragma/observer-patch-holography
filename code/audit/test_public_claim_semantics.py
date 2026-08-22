"""Regression gates for high-risk public epistemic-status wording."""

from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def _highlight_blocks(text: str, heading: str, terminator: str) -> list[str]:
    """Return the numbered receipt sections between the heading and terminator.

    Receipts are discovered by their numbered subheadings, so adding or
    removing one changes the returned list rather than breaking the parse.
    """
    start = text.index(heading)
    body = text[start : text.index(terminator, start)]
    starts = [match.start() for match in re.finditer(r"(?m)^### \d+\. ", body)]
    bounds = starts + [len(body)]
    return [body[bounds[i] : bounds[i + 1]] for i in range(len(starts))]


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

    readme = " ".join(surfaces["English README"].read_text(encoding="utf-8").lower().split())
    french = " ".join(surfaces["French README"].read_text(encoding="utf-8").lower().split())
    # The reader-facing case states the Koide premise and presents capacity
    # as a high-level scientific program. Detailed negative/comparison status
    # belongs in the technical surfaces, not the marketing README files.
    assert "stated balance premise" in readme
    assert "the capacity program asks whether" in readme
    assert "so neither is a prediction" not in readme
    assert "prémisse d’équilibre déclarée" in french
    assert "le programme de capacité demande si" in french
    assert "aucun des deux n’est une prédiction" not in french


def test_readme_case_section_remains_short_and_reader_facing() -> None:
    # The case section is a fixed set of bold-led bullets; each stays short,
    # reader-facing, and lightly linked, so the front door never regrows a
    # receipt inventory. The French renderings run longer for the same
    # content, so each surface carries its own word budget.
    surfaces = [
        (
            ROOT / "README.md",
            "## One Architecture, All Of Physics",
            120,
        ),
        (
            ROOT / "README_FR.md",
            "## Une seule architecture, toute la physique",
            140,
        ),
    ]
    for path, heading, word_budget in surfaces:
        text = path.read_text(encoding="utf-8")
        start = text.index(heading) + len(heading)
        body = text[start : text.index("\n## ", start)]
        bullets = re.split(r"(?m)^- ", body)[1:]
        assert 5 <= len(bullets) <= 8, f"{path.name} lists {len(bullets)} case bullets"
        for number, bullet in enumerate(bullets, start=1):
            assert bullet.startswith("**"), (
                f"{path.name} case bullet {number} drops its bold lead"
            )
            prose = re.sub(r"\[([^]]+)\]\([^)]+\)", r"\1", bullet)
            words = re.findall(r"\b[\wÀ-ÿ]+(?:[-’'][\wÀ-ÿ]+)*\b", prose)
            assert len(words) <= word_budget, (
                f"{path.name} case bullet {number} has {len(words)} words; "
                "detail belongs in the linked technical surfaces"
            )
            assert len(re.findall(r"\]\(", bullet)) <= 3, (
                f"{path.name} case bullet {number} links more than three destinations"
            )

    # Flattened so the checks survive the line wrapping of the README prose.
    english = " ".join((ROOT / "README.md").read_text(encoding="utf-8").split())
    french = " ".join((ROOT / "README_FR.md").read_text(encoding="utf-8").split())
    assert "lands within 72 eV" not in english
    assert "tombe à 72 eV" not in french
    assert "carries diagnostic status" in english
    assert "statut diagnostique" in french


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
