"""Tree-file-level stale-absence guards for the E4 status matrix.

Companion to ``Lean/QFT/InheritanceMatrixGuards.lean``.  The Lean module
scans the loaded environment (pinned Mathlib plus the import closure of
``QFT.InheritanceMatrix``) on every build; this test covers what that
scan structurally cannot see:

* tree files outside the import closure (a future stronger
  ``Lean/QFT/TimeSlice.lean`` appearing anywhere, a spin structure or
  Cauchy embedding landing in an unimported module);
* comment-level citations, above all the Tomita TODO in
  ``Mathlib/Analysis/InnerProductSpace/StandardSubspace.lean`` that row 1
  cites verbatim; the TODO being resolved is exactly the "blocker
  starts moving" event the guard exists to catch;
* the row-4 presence citations (finite response/Hodge precursors and the
  B15/B16 partial artifacts),
  which live outside Lean altogether.

Fail-closed contract: an unreadable or missing citation path is a
FAILURE, never a pass or a skip.  If this test cannot evaluate a cited
absence (because pinned Mathlib is not provisioned, a file moved, or
the guard tables drifted) that is a finding.

What this test does NOT cover: it checks that named probes are absent
from pinned files.  It does not prove any matrix row's claim that a
target is unstatable by any route; absence of a name is not absence of
a concept, and a structure landing under an unprobed name passes.
"""

from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
LEAN_DIR = REPO_ROOT / "Lean"
MATHLIB_SRC = LEAN_DIR / ".lake" / "packages" / "mathlib" / "Mathlib"
GUARD_LEAN = LEAN_DIR / "QFT" / "InheritanceMatrixGuards.lean"
MATRIX_LEAN = LEAN_DIR / "QFT" / "InheritanceMatrix.lean"

# Mirrors of the Lean guard tables.  test_probe_tables_in_sync fails if
# either side drifts, so neither table can be edited alone.
MATHLIB_ABSENCE_PROBES = {
    1: ["tomita", "modularconjugation", "modularoperator", "kms"],
    6: ["lorentz", "globallyhyperbolic", "globalhyperbolic",
        "cauchysurface", "cauchyembedding", "pseudoriemannian",
        "semiriemannian"],
    7: ["haag", "ruelle", "scattering", "waveoperator"],
}
TREE_ABSENCE_PROBES = {
    1: ["tomita"],
    2: ["spinstructure"],
    3: ["quasilocal", "vacuumrepresentation", "doplicherroberts"],
    6: ["cauchyembedding"],
    7: ["energymomentum", "spectrumcondition", "spectraladapter"],
}

# Row-1 citation: the TODO line this exit quotes.  Its disappearance
# means the blocker moved, which is a red for the row as written.
TOMITA_TODO_FILE = MATHLIB_SRC / "Analysis" / "InnerProductSpace" / "StandardSubspace.lean"
TOMITA_TODO_TEXT = "Define the Tomita conjugation"

# Row-2 citation: SpinGroup exists but is unconnected.
SPIN_GROUP_FILE = MATHLIB_SRC / "LinearAlgebra" / "CliffordAlgebra" / "SpinGroup.lean"
SPIN_GROUP_DECL = "def spinGroup"

# Row-6 re-review trigger: a dedicated stronger TimeSlice module landing
# may change the relative-Cauchy exit even though the current finite
# PathTimeSliceInterface is already cited as present.
TIME_SLICE_FILE = LEAN_DIR / "QFT" / "TimeSlice.lean"

# Row-6 presence citation (post V3.10 wording): the typed time-indexed
# net-evolution interface exists and the row says so; its disappearance
# makes the corrected citation stale.
ROW6_PRESENT_PATHS = [
    LEAN_DIR / "QFT" / "PathTimeSliceInterface.lean",
]

# Row-4 presence citations (post-4012dea5 wording): these exist, and the
# row says so; their disappearance makes the citation stale.
ROW4_PRESENT_PATHS = [
    LEAN_DIR / "Screen" / "A5ResponseWordAlgebra.lean",
    LEAN_DIR / "Screen" / "PositionSpaceMaxwellAction.lean",
    REPO_ROOT / "code" / "b15_matter_freeze" / "matter_class_freeze_v1.json",
    REPO_ROOT / "code" / "b16_lattices" / "lattices_v1.json",
]

# A Lean declaration line, permissive on modifiers/attributes; group 1
# is the declared name.  Matching names (not prose) keeps citation text
# and doc comments from tripping the probes.
DECL_RE = re.compile(
    r"^\s*(?:@\[[^\]]*\]\s*)?"
    r"(?:noncomputable\s+|private\s+|protected\s+|partial\s+|unsafe\s+|scoped\s+)*"
    r"(?:def|theorem|lemma|structure|class|abbrev|inductive|instance|opaque|axiom)\s+"
    r"([A-Za-z0-9_.'«»₀-ₜᵢ-ᵪ!?]+)"
)

# The tree scan matches declaration names, never prose or string
# literals, so the two files that talk ABOUT the probes are scanned
# like any other file; a probe-matching declaration hiding inside the
# guard layer itself is a red, not a blind spot.
TREE_SCAN_EXCLUDE: set = set()


def _require(path: Path, findings: list[str], what: str) -> bool:
    """Fail-closed existence gate: a missing citation path is a finding."""
    if not path.exists():
        findings.append(
            f"UNREADABLE CITATION: {what} at {path} does not exist. "
            "The guard cannot evaluate the cited claim; that is a finding, not a pass."
        )
        return False
    return True


def _decl_names(path: Path) -> list[str]:
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:  # unreadable file = finding, surfaced by caller
        raise AssertionError(f"UNREADABLE CITATION PATH {path}: {exc}") from exc
    return [m.group(1) for line in text.splitlines() if (m := DECL_RE.match(line))]


def _scan(files: list[Path], probes: dict[int, list[str]], origin: str) -> list[str]:
    findings: list[str] = []
    for path in files:
        low_name = path.name.lower()
        names = None
        for row, row_probes in probes.items():
            for probe in row_probes:
                if probe in low_name:
                    findings.append(
                        f"row {row}: {origin} FILE name {path} matches absence "
                        f"probe '{probe}': cited absence has ended"
                    )
                if names is None:
                    names = _decl_names(path)
                for name in names:
                    if probe in name.lower():
                        findings.append(
                            f"row {row}: {origin} declaration '{name}' in {path} "
                            f"matches absence probe '{probe}': cited absence has ended"
                        )
    return findings


def test_mathlib_provisioned_fail_closed():
    """No pinned Mathlib source, no evaluation: that is a failure."""
    assert MATHLIB_SRC.is_dir(), (
        f"Pinned Mathlib source not found at {MATHLIB_SRC}. The stale-absence "
        "guards cannot evaluate their cited claims without it; run the Lean "
        "provisioning first. An unevaluable guard is a finding, never a green."
    )


def test_scanned_mathlib_is_the_pinned_revision():
    """The scanned checkout must be the manifest-pinned Mathlib revision.

    The other tests trust whatever sits at Lean/.lake/packages/mathlib.
    In CI the build step reconciles that checkout to the manifest, but a
    local run against a stale checkout would evaluate the absence claims
    against the wrong revision and report them as if they held for the
    pin. A checkout that cannot be resolved is a finding, never a pass.
    """
    manifest = json.loads((LEAN_DIR / "lake-manifest.json").read_text(encoding="utf-8"))
    pinned = next(
        p["rev"] for p in manifest["packages"] if p["name"] == "mathlib"
    )
    result = subprocess.run(
        ["git", "-C", str(MATHLIB_SRC.parent), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"cannot resolve the Mathlib checkout revision at {MATHLIB_SRC.parent}: "
        f"{result.stderr.strip()}. An unverifiable checkout is a finding."
    )
    checked_out = result.stdout.strip()
    assert checked_out == pinned, (
        f"Mathlib checkout {checked_out} differs from the manifest pin {pinned}: "
        "the absence claims would be evaluated against the wrong revision"
    )


def test_mathlib_absence_probes():
    """Rows 1, 6, 7: the cited pinned-Mathlib absences hold."""
    findings: list[str] = []
    if not _require(MATHLIB_SRC, findings, "pinned Mathlib source tree"):
        raise AssertionError("\n".join(findings))
    files = sorted(MATHLIB_SRC.rglob("*.lean"))
    # Fail closed on scan vacuity: pinned Mathlib has 7871 source files.
    if len(files) < 7000:
        findings.append(
            f"scan integrity: only {len(files)} Mathlib source files seen "
            "(floor 7000): this is not the library the guard claims to scan"
        )
    findings += _scan(files, MATHLIB_ABSENCE_PROBES, "Mathlib")
    assert not findings, "\n".join(findings)


def test_row1_tomita_todo_open():
    """Row 1 cites the StandardSubspace TODO verbatim; the TODO being
    resolved or reworded is the blocker moving, which is a red."""
    findings: list[str] = []
    if _require(TOMITA_TODO_FILE, findings, "row-1 cited Tomita-TODO module"):
        if TOMITA_TODO_TEXT not in TOMITA_TODO_FILE.read_text(encoding="utf-8"):
            findings.append(
                f"row 1: the cited TODO line '{TOMITA_TODO_TEXT}' is gone from "
                f"{TOMITA_TODO_FILE}: the cited blocker has moved; re-review row 1"
            )
    assert not findings, "\n".join(findings)


def test_row2_spin_group_present():
    """Row 2 cites SpinGroup.lean as present-but-unconnected."""
    findings: list[str] = []
    if _require(SPIN_GROUP_FILE, findings, "row-2 cited SpinGroup module"):
        if SPIN_GROUP_DECL not in SPIN_GROUP_FILE.read_text(encoding="utf-8"):
            findings.append(
                f"row 2: '{SPIN_GROUP_DECL}' is absent from "
                f"{SPIN_GROUP_FILE}: the citation is stale"
            )
    assert not findings, "\n".join(findings)


def test_tree_absence_probes():
    """Rows 1, 2, 3, 6, 7: the cited tree absences hold, over the
    whole Lean tree, imported or not."""
    findings: list[str] = []
    if not _require(LEAN_DIR, findings, "project Lean tree"):
        raise AssertionError("\n".join(findings))
    files = [
        p for p in sorted(LEAN_DIR.rglob("*.lean"))
        if ".lake" not in p.parts and p.resolve() not in TREE_SCAN_EXCLUDE
    ]
    if len(files) < 100:
        findings.append(
            f"scan integrity: only {len(files)} project Lean files seen "
            "(floor 100): this is not the tree the guard claims to scan"
        )
    findings += _scan(files, TREE_ABSENCE_PROBES, "project")
    if TIME_SLICE_FILE.exists():
        findings.append(
            f"row 6: {TIME_SLICE_FILE} exists: a dedicated stronger "
            "time-slice module has landed; re-review the relative-Cauchy exit"
        )
    assert not findings, "\n".join(findings)


def test_row4_presence_citations():
    """Row 4: the named finite precursors and partial artifacts exist; their
    disappearance makes the corrected citation stale."""
    findings: list[str] = []
    for path in ROW4_PRESENT_PATHS:
        _require(path, findings, "row-4 cited artifact")
    assert not findings, "\n".join(findings)


def test_row6_presence_citation():
    """Row 6 (post-V3.10 wording): the typed time-indexed interface the
    corrected row cites as present exists; its disappearance makes the
    citation stale."""
    findings: list[str] = []
    for path in ROW6_PRESENT_PATHS:
        _require(path, findings, "row-6 cited artifact")
    if _require(MATRIX_LEAN, findings, "row-6 matrix citation"):
        matrix = MATRIX_LEAN.read_text(encoding="utf-8")
        for required in (
            "QFT/PathTimeSliceInterface.lean supplies a finite",
            "e4Row6_interface_anchor",
            "relative-Cauchy evolution",
            "stress-response",
        ):
            if required not in matrix:
                findings.append(
                    f"row 6: corrected matrix citation is missing {required!r}"
                )
        for stale in (
            "QFT/TimeSlice.lean is a declared unstarted E3 deliverable",
            "QFT/TimeSlice.lean \\\n+       is a declared unstarted E3 deliverable",
        ):
            if stale in matrix:
                findings.append(
                    f"row 6: stale pre-interface absence wording remains: {stale!r}"
                )
    assert not findings, "\n".join(findings)


def test_probe_tables_in_sync():
    """Single-source-of-truth check: every probe listed here appears as
    a string literal in the Lean guard, and vice versa, so the two
    layers cannot silently drift apart."""
    findings: list[str] = []
    if not _require(GUARD_LEAN, findings, "Lean guard module"):
        raise AssertionError("\n".join(findings))
    guard_text = GUARD_LEAN.read_text(encoding="utf-8")
    lean_probes = set(re.findall(r'"([a-z]+)"', guard_text))
    here = {p for probes in (*MATHLIB_ABSENCE_PROBES.values(),
                             *TREE_ABSENCE_PROBES.values()) for p in probes}
    for probe in sorted(here):
        if probe not in lean_probes:
            findings.append(
                f"probe '{probe}' is checked here but absent from "
                f"{GUARD_LEAN}: the two guard layers have drifted"
            )
    for probe in sorted(lean_probes - here):
        findings.append(
            f"lowercase string literal '{probe}' in {GUARD_LEAN} is not a "
            "probe known to this test: either a probe is present in Lean "
            "only (the layers have drifted) or a stray lowercase literal "
            "entered the guard; reconcile the tables"
        )
    assert not findings, "\n".join(findings)
