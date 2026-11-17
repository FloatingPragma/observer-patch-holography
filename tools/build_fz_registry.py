"""Build and validate the frozen-prediction ladder surface (issue #607).

The machine-readable register is ``claims/frozen_prediction_register.json``.
This tool validates it fail-closed and renders
``docs/FROZEN_PREDICTION_LADDER.md``; ``--check`` fails when the committed
page differs from the render, and the mandatory suite runs that check.

Fail-closed rules: frozen rows carry a parseable, non-future UTC freeze time,
custody, a typed attestation state, content hash, kill band, and comparison
protocol. Historical owner numbers remain unchanged as scientific provenance;
current execution ownership lives in GitHub and the workspace DAG. Pending rows
carry an open owning issue and a milestone; that pointer may
be updated before a freeze. Retrospective results occupy a separate collection;
their former reservations cannot also occur as ladder rows. The issue-506
record is checked against a fresh replay of its canonical producer as well as
its recomputed payload digest. Committed custody contracts bind the source-side
FZ-02 receipt and Lean module. The custody directories are under
``evidence/custody/`` with two explicitly pinned metadata-normalized derivatives;
the verifier compares those against their historical originals. Every clone verifies each
manifest artifact, detached OpenTimestamps digest, attestation class, and the
append-only FZ-02 custody and scientific errata. If that directory is removed,
the tool reports ``external_custody_not_present`` explicitly rather than
claiming that the artifact set was verified. FZ-11 additionally resolves its
original and repaired source commits, hashes the historical receipt and Lean
blobs, and requires a direct-parent repair. A custody root that is itself the original custody Git checkout additionally
gets strict history checks for its append-only proof and decision-rule repairs;
the vendored copy carries the bytes and timestamps, not that history. The explicit
``--verify-fz11-lean`` gate re-elaborates the repaired proof and accepts exactly
five standard-axiom reports with no ``sorryAx``; it never skips when requested.
FZ-12 is bound independently to its exact source commit, root custody commit,
canonical source receipt, declared physical-promotion gates, and frozen snapshot.
The October 4 FZ-13/FZ-15 additions retain that historical bootstrap and bind
their complete conditional targets and protocols to separate stamped snapshots;
they do not alter historical frozen rows or establish their physical premises.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import strict_json
except ModuleNotFoundError:  # Support importlib loading from repository tests.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import strict_json

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_bitcoin_upgrade_receipt import RECEIPT as BITCOIN_UPGRADE_RECEIPT
from verify_bitcoin_upgrade_receipt import verify as verify_bitcoin_upgrade_receipt
from closed_lanes import is_open, live_owner

ROOT = Path(__file__).resolve().parents[1]
REGISTER_PATH = ROOT / "claims" / "frozen_prediction_register.json"
PENDING_TEST_OWNER = live_owner(1025)
COMMON_MODEL_OWNER = live_owner(1026)
FROZEN_ROWS_BOOTSTRAP_REVISION = "68663e9e52a3931c322676a127dd0af144a01de3"
FROZEN_ROWS_BOOTSTRAP_GIT_BLOB_SHA1 = "731ff48c6005a0d2b3930c1d33a6f00c2dd3e345"
METADATA_NORMALIZATION_PATH = ROOT / "evidence/custody/metadata_normalization_2026-10-11.json"
METADATA_NORMALIZATION_SHA256 = "29056c76c475620748d6518bc68d96ae9d55685d2f97075dfe5329c7f4e7947d"


def custody_metadata_normalization() -> dict:
    """Load the fixed, metadata-only normalization, never a target exemption."""
    raw = METADATA_NORMALIZATION_PATH.read_bytes()
    if hashlib.sha256(raw).hexdigest() != METADATA_NORMALIZATION_SHA256:
        fail("custody metadata normalization bytes differ from fixed pin")
    return strict_json.loads(raw.decode("utf-8"))


def verify_custody_artifact_digest(path: Path, expected_hash: str) -> bool:
    """Bind a displayed copy to its timestamped original without restamping it.

    Only two exact normalized byte strings are admitted. Their original bytes
    must exist at the pinned public source revision and reproduce the original
    digest. All scientific JSON fields and all other erratum text are compared.
    """
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() == expected_hash:
        return True
    metadata = custody_metadata_normalization()
    # A caller may supply a copied custody root (including mutation tests).
    matches = [
        (relative, pin) for relative, pin in metadata["artifacts"].items()
        if path.name == Path(relative).name
    ]
    if len(matches) != 1:
        return False
    relative, pin = matches[0]
    if (pin["original_sha256"] != expected_hash or
            hashlib.sha256(raw).hexdigest() != pin["normalized_sha256"]):
        return False
    original = read_commit_blob(ROOT, metadata["source_revision"], relative,
                                "custody metadata normalization")
    if hashlib.sha256(original).hexdigest() != expected_hash:
        return False
    if path.suffix == ".json":
        source = strict_json.loads(original.decode("utf-8"))
        source["pinned_source"]["file"] = Path(source["pinned_source"]["file"]).name
        return source == strict_json.loads(raw.decode("utf-8"))
    normalized, count = re.subn(
        r"(?m)^   the [^\n]+ commit$", "   the custody commit",
        original.decode("utf-8"),
    )
    return count == 1 and normalized.encode("utf-8") == raw

# Append-only October 4 extension. The V3 row schema is unchanged; the original
# bootstrap still binds every historical frozen/void row. New conditional
# cosmology rows are bound separately to their externally stamped snapshots.
COSMOLOGY_FREEZE_SCHEMA = "oph.conditional_cosmology_freeze.v1"
COSMOLOGY_FREEZE_PINS = {
    "FZ-13": {
        "snapshot_sha256": "60ae0827d652a5434046b5ce8ceb536cf75bdd0ced6b19e8e1005c3ed3d5fb89",
        "manifest_sha256": "d95e68b703661bad49f99839244bece3c4d449447e20e72ed4113a047c0d2f5d",
    },
    "FZ-15": {
        "snapshot_sha256": "150f9f6f771bbde5ec0c3548756a83d8ab3499d9ff8c0e0364388883ccead1d8",
        "manifest_sha256": "f2ab98661a5b14ae781991a0c269dcd12f7376761ecd3cac381d5246ede9d149",
    },
}
SURFACE_PATH = ROOT / "docs" / "FROZEN_PREDICTION_LADDER.md"
FZ02_RECEIPT_PATH = (
    ROOT
    / "code"
    / "a5_closure"
    / "receipts"
    / "a5_angular_multiplet_reference.receipt.json"
)
FZ11_RECEIPT_REL = (
    "code/a5_fingerprint/runtime/"
    "spin_six_primitive_port_prediction_receipt.json"
)
FZ11_RECEIPT_PATH = ROOT / FZ11_RECEIPT_REL
FZ11_LEAN_REL = "Lean/Screen/A5PrimitivePortPrediction.lean"
FZ11_SOURCE_COMMIT = "66176656dc1143f9ec50ba1a6e409c403545857f"
FZ11_CUSTODY_COMMIT = "97202365784ad5bfc96c482b95b429c396afb5bf"
FZ10_SOURCE_COMMIT = "bbc2cebedf190af351fba1eb6c8887fb340eb6bb"
FZ10_CUSTODY_COMMIT = "df097d8fe7c38d008a1ba7827f7d573286ae2012"
FZ10_REGISTRATION = "FZ-10 Koide conditional tau window"
FZ10_CUSTODY_COMMIT_UTC = "2026-07-28T09:01:08Z"
FZ11_LEAN_REPAIR_SOURCE_COMMIT = "05771b773ffaef4be10ae67a72e51cc17e3a38fb"
FZ11_LEAN_REPAIR_CUSTODY_COMMIT = "c35679ab93e4c121fd2d0d51d6e253f9c5bc6257"
FZ11_DECISION_RULE_CUSTODY_COMMIT = "8cc5261653e37cbca0e6017fcc95a9fe7f649963"
FZ11_SCHEMA = "oph.spin_six_primitive_port_prediction.v1"
FZ11_STATUS = (
    "FROZEN_PROSPECTIVE_PRIMITIVE_TWELVE_PORT_BRANCH_PREDICTION__"
    "PHYSICAL_COMPARISON_UNARMED"
)
FZ11_AXIOM_DECLARATIONS = (
    "OPH.A5PrimitivePortPrediction.invariant_port_weights_are_equal",
    "OPH.A5PrimitivePortPrediction.b6_over_c4_squared",
    "OPH.A5PrimitivePortPrediction.b0_over_c4_squared",
    "OPH.A5PrimitivePortPrediction.b6_over_b0",
    "OPH.A5PrimitivePortPrediction.binary_refinement",
)
FZ11_ALLOWED_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}
FZ11_ROW_CONTENT = (
    "Primitive twelve-port propagation prediction: under the named real, "
    "reciprocal, finite-range cosine kinetic branch, the complete primitive "
    "port orbit is the sole hop support and no independent isotropic k^4 or "
    "k^6 term is present through the displayed order. Proper-carrier covariance "
    "forces equal weights and continuum normalization gives "
    "omega^2=(1/(2a^2)) sum_i[1-cos(a k u_i.n)]. Writing its expansion as "
    "k^2+C4 k^4+B0 k^6+B6 k^6 I6+O(k^8), the branch predicts C4=-a^2/20, "
    "B0=a^4/840, B6=2a^4/7875 and therefore B6/C4^2=32/315, "
    "B0/C4^2=10/21 and B6/B0=16/75. Intrinsic anisotropic ranks j=1 through "
    "5 vanish, while j=6 has the unique rotated I6 shape. After C4 fixes the "
    "scale, the only angular nuisance is one shared global orientation class "
    "valued in the three-dimensional quotient SO(3)/A5. It is not one scalar "
    "parameter and does not permit a deformation of I6. Spin six denotes "
    "angular rank, not particle spin. The scalar or polarization-independent "
    "sector bridge, coherent frame transport and exclusivity are declared "
    "physical premises; issue #655 owns their derivation or rejection."
)
FZ11_ROW_COMPARISON_PROTOCOL = (
    "The prospective primitive-port coefficient manifold and append-only "
    "corrected decision rule are frozen. Physical comparison is unarmed. An "
    "eligible comparison must test the full real, reciprocal, finite-range "
    "cosine branch with the complete primitive port orbit as sole hop support "
    "and no independent isotropic k^4 or k^6 term through the displayed order. "
    "Its dataset-specific contract must fix one post-freeze data release, the "
    "joint likelihood or full covariance, physical sector and carrier frame, "
    "one shared orientation class in the three-dimensional quotient SO(3)/A5, "
    "the boost law, source, medium, gravitational and instrumental nuisance "
    "models, trials accounting, sensitivity floor, calibrated joint coverage, "
    "and isolation of the carrier contribution before exposure. A photon test "
    "also requires equal action on both transverse polarizations. The 2026-07-17 "
    "WMAP/CMB template campaign and every data product it inspected are "
    "explicitly excluded."
)
FZ11_ROW_KILL_BAND = (
    "Scope premise: the real, reciprocal, finite-range cosine kinetic branch "
    "uses the complete primitive port orbit as its sole hop support and has no "
    "independent isotropic k^4 or k^6 term through the displayed order. "
    "FZ11-R01 (FAIL): an isolated intrinsic C4 is positive at five or more "
    "standard deviations. FZ11-R02 (FAIL): an isolated intrinsic anisotropic "
    "coefficient at angular rank one through five is nonzero at five or more "
    "standard deviations. FZ11-R03 (FAIL): for resolved negative C4 with "
    "adequate sixth-order sensitivity, the linked B0, B6, or rigid rotated I6 "
    "vector is excluded after fitting one shared global orientation class in "
    "the three-dimensional quotient SO(3)/A5 at five or more standard "
    "deviations. FZ11-R04 (FAIL): the calibrated "
    "joint likelihood excludes the complete C4<0, B0/C4^2=10/21, "
    "B6/C4^2=32/315, rotated-I6 branch manifold at five or more standard "
    "deviations. FZ11-R05 (SUPPORT): the zero-coefficient minimal locally "
    "Lorentz-invariant Standard Model plus General Relativity baseline is "
    "excluded at five or more standard deviations, the complete linked branch "
    "manifold agrees within two standard deviations, named systematic "
    "alternatives are rejected, "
    "and an independent eligible release replicates it. Every null, "
    "underpowered, incomplete-covariance, unresolved-frame, polarization-split, "
    "or non-isolated outcome is INCONCLUSIVE. FAIL rejects the primitive "
    "twelve-port physical propagation branch; it is OPH-wide only if issue #655 "
    "proves that branch forced and exclusive."
)
FZ12_RECEIPT_REL = (
    "code/a5_fingerprint/runtime/"
    "seam_current_edge_prediction_receipt.json"
)
FZ12_RECEIPT_PATH = ROOT / FZ12_RECEIPT_REL
FZ12_CARRIER_LEAN_REL = "Lean/Screen/SeamCurrentCarrierQuotient.lean"
FZ12_MOMENT_LEAN_REL = "Lean/Screen/SeamCurrentEdge30Moment.lean"
FZ12_RAY_LEAN_REL = "Lean/Screen/A5OrbitRaySeparation.lean"
FZ12_SOURCE_COMMIT = "bc5595f8dbb2d2886e2a64ddf447f69fbb00eb3f"
FZ12_CUSTODY_COMMIT = "54b450af0bb5bd0fee4842f5c5f654d08d6baa2d"
FZ12_FROZEN_UTC = "2026-08-02T11:52:27Z"
FZ12_DECISION_RULE_CUSTODY_COMMIT = "25da61a800226e0232336ccc86de8dec7d6b51c6"
FZ12_DECISION_RULE_UTC = "2026-08-02T12:37:49Z"
FZ12_DECISION_RULE_MANIFEST = "fz12_decision_rule_manifest_2026-08-02.json"
FZ12_DECISION_RULE_JSON = "fz12_decision_rule_v2_2026-08-02.json"
FZ12_DECISION_RULE_NOTE = "FZ12_DECISION_RULE_CLARIFICATION_2026-08-02.md"
FZ12_SCHEMA = "oph.seam_current_edge_prediction_candidate.v1"
FZ12_STATUS = (
    "EXACT_SOURCE_NATIVE_EDGE_RAY__PROSPECTIVE_PHYSICAL_BRANCH_UNARMED__"
    "PHYSICAL_PRODUCER_OPEN"
)
FZ12_FREEZE_SCHEMA = "oph.fz12.seam_current_edge_prediction.freeze.v1"
FZ12_FREEZE_STATUS = (
    "FROZEN_CONDITIONAL_SOURCE_SEAM_EDGE_BRANCH__PHYSICAL_PRODUCER_OPEN__"
    "COMPARISON_INELIGIBLE"
)
FZ12_CUSTODY_PATH = "falsification/frozen_targets/fz12_2026-08-02"
FZ12_TARGET_FILE = "frozen_target_seam_current_edge_prediction_2026-08-02.md"
FZ12_PREDICTION_FILE = "seam_current_edge_prediction_frozen_2026-08-02.json"
FZ12_MANIFEST_FILE = "registration_manifest_2026-08-02.json"
FZ12_IN_REPO_ARTIFACTS = {
    FZ12_RECEIPT_REL,
    FZ12_CARRIER_LEAN_REL,
    FZ12_MOMENT_LEAN_REL,
    FZ12_RAY_LEAN_REL,
}
FZ04_VERDICT_REL = (
    "code/particles/alpha_hvp_audit/outputs/alpha_hvp_class_verdict.json"
)
FZ04_VERDICT_PATH = ROOT / FZ04_VERDICT_REL
FZ04_BUILDER_PATH = (
    ROOT
    / "code"
    / "particles"
    / "alpha_hvp_audit"
    / "build_alpha_hvp_verdict.py"
)
# Byte-exact copies of the frozen-target custody directories, with their
# OpenTimestamps proofs, live in this repository so any clone can verify them.
DEFAULT_CUSTODY_ROOT = ROOT / "evidence" / "custody"
# The official OpenTimestamps client imports python-bitcoinlib, which loads
# OpenSSL through ctypes.util.find_library("ssl"). That lookup returns None on
# Windows, so the client cannot start there; POSIX platforms verify custody and
# Windows reports the explicit state below.
CUSTODY_PARSER_UNSUPPORTED = sys.platform == "win32"

SCHEMA = "oph.frozen_prediction_register.v3"
STATUSES = {
    "frozen_attested",
    "frozen_stamped_upgrade_pending",
    "standing_frozen",
    "registered_pending_freeze",
    "resource_deferred",
    "superseded_void",
}
FROZEN_STATUSES = {
    "frozen_attested",
    "frozen_stamped_upgrade_pending",
    "standing_frozen",
}

ROW_KEYS = {
    "id",
    "content",
    "status",
    "frozen_utc",
    "custody",
    "attestation",
    "content_sha256",
    "kill_band",
    "comparison_protocol",
    "owning_issue",
    "milestone",
}
REGISTER_KEYS = {
    "schema",
    "issue",
    "generated_surface",
    "policy",
    "external_custody_contracts",
    "retrospective_results",
    "rows",
}
RETROSPECTIVE_RESULT_KEYS = {
    "id",
    "former_ladder_reservation",
    "content",
    "status",
    "payload_path",
    "payload_sha256",
    "comparison_protocol",
    "evidential_boundary",
    "owning_issue",
    "milestone",
}
COMMON_CONTRACT_KEYS = {
    "rows",
    "custody_path",
    "registration_manifest",
    "registration_manifest_sha256",
    "attestation_state",
    "artifact_sha256",
    "in_repo_artifact_sha256",
}
FZ02_CONTRACT_EXTRA_KEYS = {
    "custody_commit",
    "custody_commit_utc",
    "source_commit",
    "custody_erratum",
    "custody_erratum_sha256",
    "scientific_erratum",
    "scientific_erratum_sha256",
    "target_file",
    "target_block_sha256",
    "target_payload_sha256",
}
FZ11_CONTRACT_EXTRA_KEYS = {
    "source_commit",
    "custody_commit",
    "frozen_utc",
    "target_file",
    "prediction_file",
    "original_in_repo_artifact_sha256",
    "lean_repair_source_commit",
    "lean_repair_custody_commit",
    "lean_repair_utc",
    "lean_repair_manifest",
    "lean_repair_manifest_sha256",
    "lean_repair_artifact_sha256",
    "lean_repair_attestation_state",
    "decision_rule_custody_commit",
    "decision_rule_utc",
    "decision_rule_manifest",
    "decision_rule_manifest_sha256",
    "decision_rule_artifact_sha256",
    "decision_rule_attestation_state",
}
FZ12_CONTRACT_EXTRA_KEYS = {
    "source_commit",
    "custody_commit",
    "frozen_utc",
    "target_file",
    "prediction_file",
    "decision_rule_custody_commit",
    "decision_rule_utc",
    "decision_rule_manifest",
    "decision_rule_manifest_sha256",
    "decision_rule_artifact_sha256",
    "decision_rule_attestation_state",
}
ATTESTATION_STATES = {"calendar_pending", "bitcoin_attested"}
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
GIT_COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")

FZ04_SCOPE = {
    "comparison_timing": "retrospective",
    "prospective_freeze": False,
    "independent_hvp_implementation_supplied": False,
    "empirical_input_promoted_to_source_output": False,
    "physical_alpha_prediction_emitted": False,
}
FZ04_CLAIM_BOUNDARY = (
    "The multi-class independent alpha/HVP test is not evaluable. One "
    "byte-pinned KNT19 accounting row is compatible under a secondary "
    "arithmetic replay. Raw-dispersive, independent-code, and lattice-HVP "
    "classes lack frozen repository ingests. The result is retrospective and "
    "supplies neither a prospective freeze nor a physical OPH alpha prediction."
)

# DetachedTimestampFile header followed by the SHA-256 operation tag and the
# 32-byte digest of the paired file. Reading this prefix does not require the
# optional opentimestamps Python package or network access.
OTS_DETACHED_HEADER = bytes.fromhex(
    "004f70656e54696d657374616d707300" "0050726f6f6600bf89e2e884e8929401"
)
OTS_SHA256_TAG = b"\x08"
OTS_PENDING_ATTESTATION_TAG = bytes.fromhex("83dfe30d2ef90c8e")
OTS_BITCOIN_ATTESTATION_TAG = bytes.fromhex("0588960d73d71901")
_OTS_CLI_USABLE: bool | None = None
_OTS_CLI_PATH: str | None = None
_OTS_OFFICIAL_PARSED_DIGESTS: set[str] = set()


def fail(message: str) -> None:
    raise SystemExit(f"frozen-prediction register: {message}")


def load_json(path: Path) -> dict:
    try:
        return strict_json.load(path)
    except FileNotFoundError:
        fail(f"missing input {path.relative_to(ROOT)}")
    except (json.JSONDecodeError, strict_json.DuplicateKeyError) as error:
        fail(f"invalid JSON in {path.relative_to(ROOT)}: {error}")
    raise AssertionError("unreachable")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode("utf-8")


def _run_git(
    repo_root: Path, arguments: list[str], where: str
) -> subprocess.CompletedProcess[bytes]:
    try:
        return subprocess.run(
            ["git", *arguments],
            cwd=repo_root,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except OSError as error:
        fail(f"{where}: cannot execute git: {error}")
    raise AssertionError("unreachable")


def git_checkout_root(path: Path) -> Path | None:
    """Return the exact Git checkout rooted at ``path``, when one exists.

    A parent repository does not make an arbitrary copied custody directory a
    source of historical evidence. Worktrees are accepted because Git, rather
    than a literal ``.git`` directory test, resolves their checkout root.
    """

    git_metadata_present = (path / ".git").exists()
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            cwd=path,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except OSError as error:
        if git_metadata_present:
            fail(f"Git metadata exists at {path}, but git cannot run: {error}")
        return None
    if result.returncode != 0:
        if git_metadata_present:
            fail(
                f"Git metadata exists at {path}, but the checkout cannot be "
                "resolved"
            )
        return None
    try:
        resolved = Path(result.stdout.decode("utf-8").strip()).resolve()
    except UnicodeDecodeError:
        return None
    requested = path.resolve()
    if resolved == requested:
        return resolved
    if git_metadata_present:
        fail(f"Git metadata at {path} resolves to unexpected root {resolved}")
    return None


def verify_commit_blobs(
    repo_root: Path,
    commit: str,
    expected_hashes: dict[str, str],
    where: str,
) -> dict[str, str]:
    """Resolve an ancestral commit and hash the requested blobs from Git."""

    if GIT_COMMIT_RE.fullmatch(commit) is None:
        fail(f"{where}: commit must be a full lowercase Git object id")
    resolved = _run_git(repo_root, ["cat-file", "-e", f"{commit}^{{commit}}"], where)
    if resolved.returncode != 0:
        fail(f"{where}: commit does not resolve: {commit}")
    ancestral = _run_git(repo_root, ["merge-base", "--is-ancestor", commit, "HEAD"], where)
    if ancestral.returncode != 0:
        fail(f"{where}: commit is not an ancestor of HEAD: {commit}")

    verified: dict[str, str] = {}
    for relative_path, expected_hash in expected_hashes.items():
        path = Path(relative_path)
        if path.is_absolute() or ".." in path.parts or not relative_path:
            fail(f"{where}: historical blob path must be repository-relative")
        require_sha256(expected_hash, f"{where}[{relative_path!r}]")
        blob = read_commit_blob(repo_root, commit, relative_path, where)
        actual_hash = sha256_bytes(blob)
        if actual_hash != expected_hash:
            fail(
                f"{where}: historical blob hash mismatch for {relative_path} "
                f"at {commit}: {actual_hash} != {expected_hash}"
            )
        verified[relative_path] = actual_hash
    return verified


def read_commit_blob(
    repo_root: Path, commit: str, relative_path: str, where: str
) -> bytes:
    """Read one exact historical blob after validating its relative path."""

    path = Path(relative_path)
    if path.is_absolute() or ".." in path.parts or not relative_path:
        fail(f"{where}: historical blob path must be repository-relative")
    blob = _run_git(repo_root, ["show", f"{commit}:{relative_path}"], where)
    if blob.returncode != 0:
        fail(f"{where}: commit lacks {relative_path}: {commit}")
    return blob.stdout


def verify_direct_parent(
    repo_root: Path, child: str, expected_parent: str, where: str
) -> None:
    """Require a single-parent correction commit directly atop its frozen parent."""

    result = _run_git(repo_root, ["show", "-s", "--format=%P", child], where)
    if result.returncode != 0:
        fail(f"{where}: cannot read correction-commit parents")
    try:
        parents = result.stdout.decode("ascii").strip().split()
    except UnicodeDecodeError:
        fail(f"{where}: correction-commit parent list is not ASCII")
    if parents != [expected_parent]:
        fail(
            f"{where}: {child} must be the direct single-parent child of "
            f"{expected_parent}, got {parents}"
        )


def verify_fz11_source_history(
    contract: dict[str, Any], repo_root: Path = ROOT
) -> dict[str, Any]:
    """Bind the FZ-11 original and corrected bytes to their source commits.

    Exported source trees without Git metadata retain byte-level validation but
    are explicitly classified as lacking history. Any actual checkout must
    resolve both declared commits, their blobs, and the direct-parent repair.
    """

    checkout = git_checkout_root(repo_root)
    if checkout is None:
        return {
            "state": "git_history_not_present",
            "repo_root": str(repo_root),
        }

    original_hashes = contract["original_in_repo_artifact_sha256"]
    repaired_hashes = contract["in_repo_artifact_sha256"]
    original = verify_commit_blobs(
        checkout,
        contract["source_commit"],
        original_hashes,
        "FZ-11 original source commit",
    )
    repaired = verify_commit_blobs(
        checkout,
        contract["lean_repair_source_commit"],
        repaired_hashes,
        "FZ-11 repaired source commit",
    )
    verify_direct_parent(
        checkout,
        contract["lean_repair_source_commit"],
        contract["source_commit"],
        "FZ-11 source repair ancestry",
    )
    if original[FZ11_RECEIPT_REL] != repaired[FZ11_RECEIPT_REL]:
        fail("FZ-11 source history shows prediction-receipt drift during Lean repair")
    if original[FZ11_LEAN_REL] == repaired[FZ11_LEAN_REL]:
        fail("FZ-11 source history does not contain distinct original and repaired Lean bytes")
    return {
        "state": "verified",
        "repo_root": str(checkout),
        "source_commit": contract["source_commit"],
        "repair_commit": contract["lean_repair_source_commit"],
        "original_blobs": original,
        "repaired_blobs": repaired,
        "repair_is_direct_child": True,
    }


def verify_fz12_source_history(
    contract: dict[str, Any], repo_root: Path = ROOT
) -> dict[str, Any]:
    """Bind the FZ-12 receipt and its complete input closure to one commit."""

    checkout = git_checkout_root(repo_root)
    if checkout is None:
        return {
            "state": "git_history_not_present",
            "repo_root": str(repo_root),
        }
    verified = verify_commit_blobs(
        checkout,
        contract["source_commit"],
        contract["in_repo_artifact_sha256"],
        "FZ-12 source commit",
    )
    receipt_raw = read_commit_blob(
        checkout,
        contract["source_commit"],
        FZ12_RECEIPT_REL,
        "FZ-12 historical receipt",
    )
    try:
        receipt = json.loads(receipt_raw)
    except json.JSONDecodeError as error:
        fail(f"FZ-12 historical receipt is invalid JSON: {error}")
    if receipt_raw != canonical_json_bytes(receipt):
        fail("FZ-12 historical receipt must use canonical JSON bytes")
    pins = receipt.get("parent_pins")
    if not isinstance(pins, list) or len(pins) != 4:
        fail("FZ-12 historical receipt must retain exactly four parent pins")
    historical_parent_blobs: dict[str, str] = {}
    for index, pin in enumerate(pins):
        if not isinstance(pin, dict) or set(pin) != {"bytes", "path", "role", "sha256"}:
            fail(f"FZ-12 historical parent pin {index} has invalid structure")
        relative_path = pin.get("path")
        digest = pin.get("sha256")
        byte_count = pin.get("bytes")
        if not isinstance(relative_path, str) or not isinstance(digest, str):
            fail(f"FZ-12 historical parent pin {index} lacks a typed path or hash")
        require_sha256(
            digest.removeprefix("sha256:"),
            f"FZ-12 historical parent pin {relative_path!r}",
        )
        blob = read_commit_blob(
            checkout,
            contract["source_commit"],
            relative_path,
            "FZ-12 historical parent closure",
        )
        actual_digest = "sha256:" + sha256_bytes(blob)
        if actual_digest != digest or len(blob) != byte_count:
            fail(
                "FZ-12 historical parent pin mismatch at source commit: "
                f"{relative_path}"
            )
        historical_parent_blobs[relative_path] = actual_digest.removeprefix("sha256:")
    return {
        "state": "verified",
        "repo_root": str(checkout),
        "source_commit": contract["source_commit"],
        "blobs": verified,
        "historical_parent_blobs": historical_parent_blobs,
    }


def verify_fz12_custody_history(
    contract: dict[str, Any], custody_root: Path
) -> dict[str, Any]:
    """Bind the original FZ-12 package and append-only rule clarification."""

    checkout = git_checkout_root(custody_root)
    if checkout is None:
        return {
            "state": "git_history_not_present",
            "repo_root": str(custody_root),
        }
    custody_rel = Path(contract["custody_path"])
    original_expected: dict[str, str] = {
        (custody_rel / contract["registration_manifest"]).as_posix(): contract[
            "registration_manifest_sha256"
        ]
    }
    for filename, digest in contract["artifact_sha256"].items():
        original_expected[(custody_rel / filename).as_posix()] = digest
    decision_expected = dict(original_expected)
    decision_expected[
        (custody_rel / contract["decision_rule_manifest"]).as_posix()
    ] = contract["decision_rule_manifest_sha256"]
    for filename, digest in contract["decision_rule_artifact_sha256"].items():
        decision_expected[(custody_rel / filename).as_posix()] = digest
    original = verify_commit_blobs(
        checkout,
        contract["custody_commit"],
        original_expected,
        "FZ-12 root custody commit",
    )
    clarified = verify_commit_blobs(
        checkout,
        contract["decision_rule_custody_commit"],
        decision_expected,
        "FZ-12 decision-rule root custody commit",
    )
    original_proof_count = verify_commit_ots_bindings(
        checkout,
        contract["custody_commit"],
        original_expected,
        contract["attestation_state"],
        "FZ-12 root custody proofs",
    )
    decision_rule_proof_count = verify_commit_ots_bindings(
        checkout,
        contract["decision_rule_custody_commit"],
        decision_expected,
        contract["decision_rule_attestation_state"],
        "FZ-12 decision-rule root custody proofs",
    )
    verify_direct_parent(
        checkout,
        contract["decision_rule_custody_commit"],
        contract["custody_commit"],
        "FZ-12 decision-rule custody ancestry",
    )
    for path, digest in original.items():
        if clarified.get(path) != digest:
            fail(
                "FZ-12 decision-rule clarification changed original custody blob "
                f"{path}"
            )
    return {
        "state": "verified",
        "repo_root": str(checkout),
        "custody_commit": contract["custody_commit"],
        "decision_rule_commit": contract["decision_rule_custody_commit"],
        "original_blob_count": len(original),
        "decision_rule_blob_count": len(clarified),
        "original_proof_count": original_proof_count,
        "decision_rule_proof_count": decision_rule_proof_count,
        "decision_rule_is_direct_child": True,
    }


def verify_fz11_custody_history(
    contract: dict[str, Any], custody_root: Path, directory: Path
) -> dict[str, Any]:
    """Verify the root custody commits when the coordinated checkout exists."""

    checkout = git_checkout_root(custody_root)
    if checkout is None:
        return {
            "state": "git_history_not_present",
            "repo_root": str(custody_root),
        }

    custody_rel = Path(contract["custody_path"])
    original_expected: dict[str, str] = {
        (custody_rel / contract["registration_manifest"]).as_posix(): contract[
            "registration_manifest_sha256"
        ]
    }
    for filename, digest in contract["artifact_sha256"].items():
        original_expected[(custody_rel / filename).as_posix()] = digest

    repair_expected = dict(original_expected)
    repair_expected[
        (custody_rel / contract["lean_repair_manifest"]).as_posix()
    ] = contract["lean_repair_manifest_sha256"]
    for filename, digest in contract["lean_repair_artifact_sha256"].items():
        repair_expected[(custody_rel / filename).as_posix()] = digest
    decision_expected = dict(repair_expected)
    decision_expected[
        (custody_rel / contract["decision_rule_manifest"]).as_posix()
    ] = contract["decision_rule_manifest_sha256"]
    for filename, digest in contract["decision_rule_artifact_sha256"].items():
        decision_expected[(custody_rel / filename).as_posix()] = digest
    original = verify_commit_blobs(
        checkout,
        contract["custody_commit"],
        original_expected,
        "FZ-11 original root custody commit",
    )
    repaired = verify_commit_blobs(
        checkout,
        contract["lean_repair_custody_commit"],
        repair_expected,
        "FZ-11 repaired root custody commit",
    )
    decision_corrected = verify_commit_blobs(
        checkout,
        contract["decision_rule_custody_commit"],
        decision_expected,
        "FZ-11 decision-rule root custody commit",
    )
    original_proof_count = verify_commit_ots_bindings(
        checkout,
        contract["custody_commit"],
        original_expected,
        "calendar_pending",
        "FZ-11 original root custody proofs",
    )
    repair_proof_count = verify_commit_ots_bindings(
        checkout,
        contract["lean_repair_custody_commit"],
        repair_expected,
        "calendar_pending",
        "FZ-11 repaired root custody proofs",
    )
    decision_rule_proof_count = verify_commit_ots_bindings(
        checkout,
        contract["decision_rule_custody_commit"],
        decision_expected,
        "calendar_pending",
        "FZ-11 decision-rule root custody proofs",
    )
    verify_direct_parent(
        checkout,
        contract["lean_repair_custody_commit"],
        contract["custody_commit"],
        "FZ-11 root custody repair ancestry",
    )
    verify_direct_parent(
        checkout,
        contract["decision_rule_custody_commit"],
        contract["lean_repair_custody_commit"],
        "FZ-11 decision-rule custody ancestry",
    )
    for path, digest in original.items():
        if repaired.get(path) != digest:
            fail(f"FZ-11 root custody repair changed frozen historical blob {path}")
    for path, digest in repaired.items():
        if decision_corrected.get(path) != digest:
            fail(
                "FZ-11 decision-rule correction changed prior historical blob "
                f"{path}"
            )
    return {
        "state": "verified",
        "repo_root": str(checkout),
        "custody_commit": contract["custody_commit"],
        "repair_commit": contract["lean_repair_custody_commit"],
        "decision_rule_commit": contract["decision_rule_custody_commit"],
        "original_blob_count": len(original),
        "repair_blob_count": len(repaired),
        "decision_rule_blob_count": len(decision_corrected),
        "original_proof_count": original_proof_count,
        "repair_proof_count": repair_proof_count,
        "decision_rule_proof_count": decision_rule_proof_count,
        "repair_is_direct_child": True,
        "decision_rule_is_direct_child": True,
    }


def parse_fz11_axiom_reports(output: str) -> dict[str, list[str]]:
    """Parse the exact five ``#print axioms`` reports emitted by FZ-11."""

    if "sorryAx" in output:
        fail("FZ-11 Lean replay reported sorryAx")
    report_start = re.compile(r"^'([^']+)' depends on axioms: ", re.MULTILINE)
    starts = list(report_start.finditer(output))
    if len(starts) != len(FZ11_AXIOM_DECLARATIONS):
        fail(
            "FZ-11 Lean replay must emit exactly five #print axioms reports; "
            f"got {len(starts)}"
        )

    reports: dict[str, list[str]] = {}
    for index, match in enumerate(starts):
        end = starts[index + 1].start() if index + 1 < len(starts) else len(output)
        body = output[match.end() : end].strip()
        if not body.startswith("[") or not body.endswith("]"):
            fail(f"FZ-11 Lean replay has a malformed axiom report for {match.group(1)}")
        name = match.group(1)
        if name in reports:
            fail(f"FZ-11 Lean replay repeated the axiom report for {name}")
        axioms = [part.strip() for part in body[1:-1].split(",") if part.strip()]
        reports[name] = axioms

    if tuple(reports) != FZ11_AXIOM_DECLARATIONS:
        fail(
            "FZ-11 Lean replay reported the wrong declarations or order: "
            f"{list(reports)}"
        )
    for name, axioms in reports.items():
        if set(axioms) != FZ11_ALLOWED_AXIOMS or len(axioms) != len(
            FZ11_ALLOWED_AXIOMS
        ):
            fail(f"FZ-11 Lean replay reported unexpected axioms for {name}: {axioms}")
    return reports


MATHLIB_OLEAN_CANDIDATES = (
    ROOT / "Lean" / ".lake" / "packages" / "mathlib" / ".lake" / "build" / "lib"
    / "lean" / "Mathlib.olean",
    ROOT / "Lean" / ".lake" / "packages" / "mathlib" / ".lake" / "build" / "lib"
    / "Mathlib.olean",
)


def mathlib_oleans_built() -> bool:
    """True when compiled Mathlib oleans exist under the Lean workspace."""

    if any(path.is_file() for path in MATHLIB_OLEAN_CANDIDATES):
        return True
    lib = ROOT / "Lean" / ".lake" / "packages" / "mathlib" / ".lake" / "build" / "lib"
    return lib.is_dir() and next(lib.glob("**/Mathlib.olean"), None) is not None


def fz11_lean_replay_available() -> bool:
    """True when `lake` is on PATH and Mathlib is built, so the replay can elaborate.

    The optional pytest replay uses this; the explicit `--verify-fz11-lean` gate
    never skips.
    """

    return shutil.which("lake") is not None and mathlib_oleans_built()


def verify_fz11_lean_replay() -> dict[str, Any]:
    """Re-elaborate FZ-11 and verify its five public axiom reports fail-closed."""

    lake = shutil.which("lake")
    if lake is None:
        fail(
            "FZ-11 Lean replay requested, but lake is unavailable; the explicit "
            "replay gate never skips"
        )
    try:
        replay = subprocess.run(
            [lake, "env", "lean", "Screen/A5PrimitivePortPrediction.lean"],
            cwd=ROOT / "Lean",
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    except OSError as error:
        fail(f"FZ-11 Lean replay could not execute lake: {error}")
    if replay.returncode != 0:
        tail = replay.stdout[-4000:]
        fail(f"FZ-11 Lean replay exited {replay.returncode}:\n{tail}")
    reports = parse_fz11_axiom_reports(replay.stdout)
    return {
        "state": "verified",
        "lean_exit_code": replay.returncode,
        "report_count": len(reports),
        "declarations": list(reports),
        "sorry_ax_present": False,
    }


def rebuild_issue506_verdict() -> dict[str, Any]:
    """Replay the canonical producer instead of trusting its stored digest."""

    spec = importlib.util.spec_from_file_location(
        "_oph_alpha_hvp_verdict_builder", FZ04_BUILDER_PATH
    )
    if spec is None or spec.loader is None:
        fail("cannot import the issue-506 canonical verdict producer")
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
        rebuilt = module.build_verdict()
    except Exception as error:
        fail(f"issue-506 canonical verdict replay failed: {error}")
    if not isinstance(rebuilt, dict):
        fail("issue-506 canonical verdict producer did not return an object")
    return rebuilt


def require_sha256(value: Any, where: str) -> str:
    if not isinstance(value, str) or not SHA256_RE.fullmatch(value):
        fail(f"{where} must be a lowercase SHA-256 hex digest")
    return value


def parse_utc(value: Any, where: str, *, reject_future: bool = True) -> datetime:
    if not isinstance(value, str) or not value.endswith("Z"):
        fail(f"{where} must be an ISO-8601 UTC timestamp ending in Z")
    try:
        parsed = datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError:
        fail(f"{where} is not a valid ISO-8601 UTC timestamp")
    if parsed.utcoffset() != timezone.utc.utcoffset(parsed):
        fail(f"{where} must be UTC")
    if reject_future and parsed > datetime.now(timezone.utc):
        fail(f"{where} cannot be in the future")
    return parsed


def validate_hash_mapping(mapping: Any, where: str) -> dict[str, str]:
    if not isinstance(mapping, dict):
        fail(f"{where} must be an object")
    for path, digest in mapping.items():
        if not isinstance(path, str) or not path or Path(path).is_absolute():
            fail(f"{where} keys must be nonempty relative paths")
        require_sha256(digest, f"{where}[{path!r}]")
    return mapping


def validate_custody_contracts(
    register: dict, rows_by_id: dict[str, dict]
) -> dict[str, dict]:
    contracts = register.get("external_custody_contracts")
    expected_contracts = {
        "FZ-01", "FZ-02", "FZ-10", "FZ-11", "FZ-12", "FZ-13", "FZ-15"
    }
    if not isinstance(contracts, dict) or set(contracts) != expected_contracts:
        fail(
            "external_custody_contracts must contain exactly FZ-01, FZ-02, "
            "FZ-10, FZ-11, FZ-12, FZ-13, and FZ-15"
        )

    claimed_rows: set[str] = set()
    for contract_id, contract in contracts.items():
        if contract_id == "FZ-02":
            expected_keys = COMMON_CONTRACT_KEYS | FZ02_CONTRACT_EXTRA_KEYS
        elif contract_id == "FZ-11":
            expected_keys = COMMON_CONTRACT_KEYS | FZ11_CONTRACT_EXTRA_KEYS
        elif contract_id == "FZ-12":
            expected_keys = COMMON_CONTRACT_KEYS | FZ12_CONTRACT_EXTRA_KEYS
        else:
            expected_keys = COMMON_CONTRACT_KEYS
        if not isinstance(contract, dict) or set(contract) != expected_keys:
            fail(f"external_custody_contracts[{contract_id}]: keys mismatch")
        row_ids = contract["rows"]
        if (
            not isinstance(row_ids, list)
            or not row_ids
            or any(not isinstance(row_id, str) for row_id in row_ids)
        ):
            fail(f"external_custody_contracts[{contract_id}].rows must be nonempty")
        for row_id in row_ids:
            if row_id not in rows_by_id:
                fail(
                    f"external_custody_contracts[{contract_id}] names unknown row {row_id}"
                )
            if row_id in claimed_rows:
                fail(f"custody row {row_id} is claimed by more than one contract")
            claimed_rows.add(row_id)

        custody_path = contract["custody_path"]
        manifest_name = contract["registration_manifest"]
        if (
            not isinstance(custody_path, str)
            or Path(custody_path).is_absolute()
            or not isinstance(manifest_name, str)
            or Path(manifest_name).name != manifest_name
        ):
            fail(f"external_custody_contracts[{contract_id}] paths must be relative")
        require_sha256(
            contract["registration_manifest_sha256"],
            f"external_custody_contracts[{contract_id}].registration_manifest_sha256",
        )
        artifacts = validate_hash_mapping(
            contract["artifact_sha256"],
            f"external_custody_contracts[{contract_id}].artifact_sha256",
        )
        in_repo = validate_hash_mapping(
            contract["in_repo_artifact_sha256"],
            f"external_custody_contracts[{contract_id}].in_repo_artifact_sha256",
        )
        state = contract["attestation_state"]
        if state not in ATTESTATION_STATES:
            fail(
                f"external_custody_contracts[{contract_id}] has unknown attestation state"
            )

        for relative_path, expected_hash in in_repo.items():
            source_path = ROOT / relative_path
            if not source_path.is_file():
                fail(f"{contract_id} missing in-repo custody artifact {relative_path}")
            if sha256_file(source_path) != expected_hash:
                fail(
                    f"{contract_id} in-repo custody artifact hash mismatch: {relative_path}"
                )

        for row_id in row_ids:
            row = rows_by_id[row_id]
            attestation = str(row["attestation"] or "").lower()
            if state == "calendar_pending":
                if row["status"] != "frozen_stamped_upgrade_pending":
                    fail(f"{row_id} calendar-pending custody requires pending status")
                if "pending" not in attestation or "bitcoin" not in attestation:
                    fail(f"{row_id} must state that its Bitcoin upgrade is pending")
                if (
                    "complete bitcoin" in attestation
                    or "bitcoin-attested" in attestation
                ):
                    fail(f"{row_id} must not claim a completed Bitcoin attestation")
            elif row["status"] == "superseded_void":
                if row["frozen_utc"] is not None:
                    fail(f"{row_id} superseded-void custody cannot retain a freeze time")
                if "bitcoin" not in attestation:
                    fail(f"{row_id} must name the attestation on its historical bytes")
            elif row["status"] not in {"frozen_attested", "standing_frozen"}:
                fail(f"{row_id} Bitcoin custody requires an attested frozen status")
            elif "bitcoin" not in attestation:
                fail(f"{row_id} must name its Bitcoin attestation")

        if contract_id == "FZ-02":
            for key in (
                "registration_manifest_sha256",
                "custody_erratum_sha256",
                "scientific_erratum_sha256",
                "target_block_sha256",
                "target_payload_sha256",
            ):
                require_sha256(
                    contract[key], f"external_custody_contracts[FZ-02].{key}"
                )
            for key in ("custody_commit", "source_commit"):
                if not isinstance(contract[key], str) or not GIT_COMMIT_RE.fullmatch(
                    contract[key]
                ):
                    fail(
                        f"external_custody_contracts[FZ-02].{key} must be a full commit"
                    )
            custody_time = parse_utc(
                contract["custody_commit_utc"],
                "external_custody_contracts[FZ-02].custody_commit_utc",
            )
            row_time = parse_utc(rows_by_id["FZ-02"]["frozen_utc"], "FZ-02.frozen_utc")
            if row_time != custody_time:
                fail("FZ-02 frozen_utc must equal the corrected custody commit time")
            if contract["custody_commit"] != "1e7d7c73dadeef9aa10ec60061a85cee8426c5b1":
                fail(
                    "FZ-02 custody commit must equal the append-only correction record"
                )
            if contract["source_commit"] != "091658ce585c107a260e7b980352be904d2419b2":
                fail(
                    "FZ-02 source commit must contain the frozen receipt and Lean module"
                )
            if contract["target_file"] not in artifacts:
                fail("FZ-02 target_file must be present in its artifact hash contract")
            for key in ("custody_erratum", "scientific_erratum", "target_file"):
                value = contract[key]
                if not isinstance(value, str) or Path(value).name != value:
                    fail(f"external_custody_contracts[FZ-02].{key} must be a file name")

        if contract_id == "FZ-11":
            for key in (
                "source_commit",
                "custody_commit",
                "lean_repair_source_commit",
                "lean_repair_custody_commit",
                "decision_rule_custody_commit",
            ):
                if not isinstance(contract[key], str) or not GIT_COMMIT_RE.fullmatch(
                    contract[key]
                ):
                    fail(
                        f"external_custody_contracts[FZ-11].{key} must be a full commit"
                    )
            if contract["source_commit"] != FZ11_SOURCE_COMMIT:
                fail("FZ-11 source commit must contain the frozen producer bytes")
            if contract["custody_commit"] != FZ11_CUSTODY_COMMIT:
                fail("FZ-11 custody commit must contain the stamped freeze package")
            if contract["lean_repair_source_commit"] != FZ11_LEAN_REPAIR_SOURCE_COMMIT:
                fail("FZ-11 Lean repair source commit must contain the corrected proof")
            if contract["lean_repair_custody_commit"] != FZ11_LEAN_REPAIR_CUSTODY_COMMIT:
                fail("FZ-11 Lean repair custody commit must contain the append-only erratum")
            if (
                contract["decision_rule_custody_commit"]
                != FZ11_DECISION_RULE_CUSTODY_COMMIT
            ):
                fail(
                    "FZ-11 decision-rule custody commit must contain the "
                    "append-only correction"
                )
            contract_time = parse_utc(
                contract["frozen_utc"],
                "external_custody_contracts[FZ-11].frozen_utc",
            )
            row_time = parse_utc(rows_by_id["FZ-11"]["frozen_utc"], "FZ-11.frozen_utc")
            if row_time != contract_time:
                fail("FZ-11 row and custody contract must carry one freeze time")
            repair_time = parse_utc(
                contract["lean_repair_utc"],
                "external_custody_contracts[FZ-11].lean_repair_utc",
            )
            if repair_time <= contract_time:
                fail("FZ-11 Lean repair time must follow the original prediction freeze")
            decision_rule_time = parse_utc(
                contract["decision_rule_utc"],
                "external_custody_contracts[FZ-11].decision_rule_utc",
            )
            if decision_rule_time <= repair_time:
                fail("FZ-11 decision-rule correction must follow the Lean repair")
            for key in ("target_file", "prediction_file"):
                value = contract[key]
                if not isinstance(value, str) or Path(value).name != value:
                    fail(f"external_custody_contracts[FZ-11].{key} must be a file name")
                if value not in artifacts:
                    fail(f"FZ-11 {key} must be present in its artifact hash contract")
            if set(in_repo) != {FZ11_RECEIPT_REL, FZ11_LEAN_REL}:
                fail("FZ-11 must pin exactly its prediction receipt and Lean proof")
            original_in_repo = validate_hash_mapping(
                contract["original_in_repo_artifact_sha256"],
                "external_custody_contracts[FZ-11].original_in_repo_artifact_sha256",
            )
            if set(original_in_repo) != {FZ11_RECEIPT_REL, FZ11_LEAN_REL}:
                fail("FZ-11 original manifest must retain its two source pins")
            if original_in_repo[FZ11_RECEIPT_REL] != in_repo[FZ11_RECEIPT_REL]:
                fail("FZ-11 Lean repair must not alter the frozen prediction receipt")
            if original_in_repo[FZ11_LEAN_REL] == in_repo[FZ11_LEAN_REL]:
                fail("FZ-11 Lean erratum must bind distinct original and repaired proof bytes")
            require_sha256(
                contract["lean_repair_manifest_sha256"],
                "external_custody_contracts[FZ-11].lean_repair_manifest_sha256",
            )
            repair_artifacts = validate_hash_mapping(
                contract["lean_repair_artifact_sha256"],
                "external_custody_contracts[FZ-11].lean_repair_artifact_sha256",
            )
            if set(repair_artifacts) != {
                "A5PrimitivePortPrediction_repaired.lean",
                "LEAN_PROOF_ERRATUM_2026-07-31.md",
            }:
                fail("FZ-11 Lean repair artifact contract is incomplete")
            if (
                not isinstance(contract["lean_repair_manifest"], str)
                or Path(contract["lean_repair_manifest"]).name
                != contract["lean_repair_manifest"]
            ):
                fail("FZ-11 lean_repair_manifest must be a file name")
            if contract["lean_repair_attestation_state"] not in ATTESTATION_STATES:
                fail("FZ-11 Lean repair has an unknown attestation state")
            require_sha256(
                contract["decision_rule_manifest_sha256"],
                "external_custody_contracts[FZ-11].decision_rule_manifest_sha256",
            )
            decision_artifacts = validate_hash_mapping(
                contract["decision_rule_artifact_sha256"],
                "external_custody_contracts[FZ-11].decision_rule_artifact_sha256",
            )
            if set(decision_artifacts) != {
                "FZ11_DECISION_RULE_ERRATUM_2026-07-31.md",
                "fz11_decision_rule_v2_2026-07-31.json",
            }:
                fail("FZ-11 decision-rule artifact contract is incomplete")
            if (
                not isinstance(contract["decision_rule_manifest"], str)
                or Path(contract["decision_rule_manifest"]).name
                != contract["decision_rule_manifest"]
            ):
                fail("FZ-11 decision_rule_manifest must be a file name")
            if contract["decision_rule_attestation_state"] not in ATTESTATION_STATES:
                fail("FZ-11 decision-rule correction has an unknown attestation state")

        if contract_id == "FZ-12":
            for key in (
                "source_commit",
                "custody_commit",
                "decision_rule_custody_commit",
            ):
                if not isinstance(contract[key], str) or not GIT_COMMIT_RE.fullmatch(
                    contract[key]
                ):
                    fail(
                        f"external_custody_contracts[FZ-12].{key} must be a full commit"
                    )
            if contract["source_commit"] != FZ12_SOURCE_COMMIT:
                fail("FZ-12 source commit must contain the frozen producer bytes")
            if contract["custody_commit"] != FZ12_CUSTODY_COMMIT:
                fail("FZ-12 custody commit must contain the stamped freeze package")
            if (
                contract["decision_rule_custody_commit"]
                != FZ12_DECISION_RULE_CUSTODY_COMMIT
            ):
                fail(
                    "FZ-12 decision-rule custody commit must contain the "
                    "append-only clarification"
                )
            contract_time = parse_utc(
                contract["frozen_utc"],
                "external_custody_contracts[FZ-12].frozen_utc",
            )
            row_time = parse_utc(rows_by_id["FZ-12"]["frozen_utc"], "FZ-12.frozen_utc")
            if (
                contract["frozen_utc"] != FZ12_FROZEN_UTC
                or rows_by_id["FZ-12"]["frozen_utc"] != FZ12_FROZEN_UTC
                or row_time != contract_time
            ):
                fail("FZ-12 row and custody contract must carry the fixed freeze time")
            decision_rule_time = parse_utc(
                contract["decision_rule_utc"],
                "external_custody_contracts[FZ-12].decision_rule_utc",
            )
            if (
                decision_rule_time <= contract_time
                or contract["decision_rule_utc"] != FZ12_DECISION_RULE_UTC
            ):
                fail("FZ-12 decision-rule clarification must follow the original freeze")
            if contract["rows"] != ["FZ-12"]:
                fail("FZ-12 custody contract must bind exactly row FZ-12")
            if contract["custody_path"] != FZ12_CUSTODY_PATH:
                fail("FZ-12 custody path drifted from its stamped package")
            if contract["registration_manifest"] != FZ12_MANIFEST_FILE:
                fail("FZ-12 registration manifest file name drifted")
            if contract["target_file"] != FZ12_TARGET_FILE:
                fail("FZ-12 target file name drifted")
            if contract["prediction_file"] != FZ12_PREDICTION_FILE:
                fail("FZ-12 prediction file name drifted")
            if set(artifacts) != {FZ12_TARGET_FILE, FZ12_PREDICTION_FILE}:
                fail("FZ-12 custody contract must pin exactly its target and snapshot")
            if set(in_repo) != FZ12_IN_REPO_ARTIFACTS:
                fail(
                    "FZ-12 must pin exactly its prediction receipt and three Lean proofs"
                )
            require_sha256(
                contract["decision_rule_manifest_sha256"],
                "external_custody_contracts[FZ-12].decision_rule_manifest_sha256",
            )
            decision_artifacts = validate_hash_mapping(
                contract["decision_rule_artifact_sha256"],
                "external_custody_contracts[FZ-12].decision_rule_artifact_sha256",
            )
            if set(decision_artifacts) != {
                FZ12_DECISION_RULE_NOTE,
                FZ12_DECISION_RULE_JSON,
            }:
                fail("FZ-12 decision-rule artifact contract is incomplete")
            if contract["decision_rule_manifest"] != FZ12_DECISION_RULE_MANIFEST:
                fail("FZ-12 decision_rule_manifest file name drifted")
            if contract["decision_rule_attestation_state"] not in ATTESTATION_STATES:
                fail("FZ-12 decision-rule clarification has an unknown attestation state")
            row = rows_by_id["FZ-12"]
            if (
                row.get("status") != "frozen_stamped_upgrade_pending"
                or row.get("owning_issue") != 666
                or row.get("milestone") != "C3-V"
                or FZ12_SOURCE_COMMIT not in str(row.get("custody", ""))
                or FZ12_CUSTODY_COMMIT not in str(row.get("custody", ""))
                or FZ12_DECISION_RULE_CUSTODY_COMMIT
                not in str(row.get("custody", ""))
                or FZ12_CUSTODY_PATH not in str(row.get("custody", ""))
            ):
                fail("FZ-12 row ownership, status, or custody binding drifted")

    return contracts


def validate_fz11_prediction(
    rows_by_id: dict[str, dict], contract: dict[str, Any]
) -> None:
    """Bind FZ-11 to one canonical, prospective branch prediction."""

    raw = FZ11_RECEIPT_PATH.read_bytes()
    try:
        receipt = json.loads(raw)
    except json.JSONDecodeError as error:
        fail(f"FZ-11 prediction receipt is invalid JSON: {error}")
    if raw != canonical_json_bytes(receipt):
        fail("FZ-11 prediction receipt must use canonical JSON bytes")
    if receipt.get("schema") != FZ11_SCHEMA or receipt.get("status") != FZ11_STATUS:
        fail("FZ-11 prediction receipt schema or status drifted")
    body = {key: value for key, value in receipt.items() if key != "receipt_sha256"}
    expected_self_digest = "sha256:" + sha256_bytes(canonical_json_bytes(body))
    if receipt.get("receipt_sha256") != expected_self_digest:
        fail("FZ-11 prediction receipt self-digest mismatch")

    row = rows_by_id["FZ-11"]
    raw_digest = sha256_bytes(raw)
    if row["content_sha256"] != raw_digest:
        fail("FZ-11 content hash must equal the canonical prediction receipt bytes")
    if contract["in_repo_artifact_sha256"].get(FZ11_RECEIPT_REL) != raw_digest:
        fail("FZ-11 custody contract does not bind the canonical prediction receipt")
    if receipt.get("frozen_utc") != row["frozen_utc"]:
        fail("FZ-11 receipt and register row must carry one freeze time")
    if receipt.get("issue") != 655 or row.get("owning_issue") != 655:
        fail("FZ-11 must retain historical physical-bridge issue #655")
    if FZ11_DECISION_RULE_CUSTODY_COMMIT not in str(row.get("custody", "")):
        fail("FZ-11 row does not name the append-only decision-rule custody")
    exact_row_fields = {
        "content": FZ11_ROW_CONTENT,
        "comparison_protocol": FZ11_ROW_COMPARISON_PROTOCOL,
        "kill_band": FZ11_ROW_KILL_BAND,
    }
    for field, expected in exact_row_fields.items():
        if row.get(field) != expected:
            fail(f"FZ-11 {field} must equal its exact bounded registry contract")

    scope = receipt.get("prediction_scope", {})
    if (
        scope.get("name") != "primitive twelve-port scalar propagation branch"
        or scope.get("type") != "prospective conditional physical-branch prediction"
        or "not particle spin" not in str(scope.get("spin_language", ""))
        or "#655" not in str(scope.get("stronger_derivation_open", ""))
    ):
        fail("FZ-11 prediction scope drifted")
    premises = receipt.get("branch_premises", {})
    if set(premises) != {
        "support",
        "equal_weights",
        "normalization",
        "finite_scale",
        "physical_sector",
        "frame",
        "exclusivity",
    }:
        fail("FZ-11 must retain every named physical branch premise")

    exact = receipt.get("exact_prediction", {})
    if exact.get("coefficients") != {
        "C4_over_a2": "-1/20",
        "B0_over_a4": "1/840",
        "B6_over_a4": "2/7875",
    }:
        fail("FZ-11 primitive-port coefficients drifted")
    if exact.get("scale_free_relations") != {
        "B6_over_C4_squared": "32/315",
        "B0_over_C4_squared": "10/21",
        "B6_over_B0": "16/75",
    }:
        fail("FZ-11 scale-free coefficient manifold drifted")
    if exact.get("signs") != {
        "C4": "negative",
        "B0": "positive",
        "B6": "positive",
    }:
        fail("FZ-11 coefficient signs drifted")
    if "j=1,2,3,4,5" not in str(exact.get("harmonic_nulls", "")):
        fail("FZ-11 must retain the rank-one-through-five angular nulls")
    if "SO(3)/A5" not in str(exact.get("fit_freedom_after_C4", "")):
        fail("FZ-11 must retain only the declared orientation nuisance")

    boundary = receipt.get("exposure_and_custody_boundary", {})
    if (
        boundary.get("comparison_permitted") is not False
        or boundary.get("new_comparison_data_read") is not False
        or boundary.get("comparison_state")
        != "UNARMED_PENDING_DATASET_SPECIFIC_PREREGISTRATION"
        or "WMAP" not in str(boundary.get("excluded_data_class", ""))
    ):
        fail("FZ-11 exposure boundary or unarmed comparison state drifted")
    decision = receipt.get("prospective_decision_rule", {})
    if (
        set(decision) != {
            "eligible_data",
            "trigger",
            "fail",
            "support",
            "inconclusive",
            "scope_of_failure",
        }
        or "five standard deviations" not in str(decision.get("fail", ""))
        or "primitive twelve-port physical propagation branch"
        not in str(decision.get("scope_of_failure", ""))
        or "#655" not in str(decision.get("scope_of_failure", ""))
    ):
        fail("FZ-11 prospective decision rule drifted")
    baseline = receipt.get("baseline_contrast", {})
    if (
        baseline.get("baseline_prediction") != "C4 = B0 = B6 = 0 for intrinsic vacuum propagation"
        or "minimal locally Lorentz-invariant Standard Model plus General Relativity"
        not in str(baseline.get("baseline", ""))
        or "imitate" not in str(baseline.get("nonuniqueness", ""))
    ):
        fail("FZ-11 baseline contrast drifted")


def validate_fz12_prediction(
    rows_by_id: dict[str, dict], contract: dict[str, Any]
) -> None:
    """Bind FZ-12 to its exact source ray and keep comparison fail-closed."""

    raw = FZ12_RECEIPT_PATH.read_bytes()
    try:
        receipt = json.loads(raw)
    except json.JSONDecodeError as error:
        fail(f"FZ-12 prediction receipt is invalid JSON: {error}")
    if raw != canonical_json_bytes(receipt):
        fail("FZ-12 prediction receipt must use canonical JSON bytes")
    expected_keys = {
        "baseline_contrast",
        "conditional_physical_candidate",
        "exact_source_result",
        "exposure_and_custody_boundary",
        "fz11_separation",
        "issue",
        "parent_pins",
        "physical_premises",
        "producer_scope",
        "promotion_gates",
        "prospective_decision_rule",
        "receipt_sha256",
        "schema",
        "status",
    }
    if set(receipt) != expected_keys:
        fail("FZ-12 prediction receipt keys drifted")
    if receipt.get("schema") != FZ12_SCHEMA or receipt.get("status") != FZ12_STATUS:
        fail("FZ-12 prediction receipt schema or status drifted")
    body = {key: value for key, value in receipt.items() if key != "receipt_sha256"}
    expected_self_digest = "sha256:" + sha256_bytes(canonical_json_bytes(body))
    if receipt.get("receipt_sha256") != expected_self_digest:
        fail("FZ-12 prediction receipt self-digest mismatch")

    row = rows_by_id["FZ-12"]
    raw_digest = sha256_bytes(raw)
    if row.get("content_sha256") != raw_digest:
        fail("FZ-12 content hash must equal the canonical prediction receipt bytes")
    if contract["in_repo_artifact_sha256"].get(FZ12_RECEIPT_REL) != raw_digest:
        fail("FZ-12 custody contract does not bind the canonical prediction receipt")
    if receipt.get("issue") != 666 or row.get("owning_issue") != 666:
        fail("FZ-12 must remain owned by physical-producer issue #666")
    if (
        "exact D6 image" not in row.get("content", "")
        or "B0/C4^2=10/21" not in row.get("content", "")
        or "B6/C4^2=-2/63" not in row.get("content", "")
        or "B6/B0=-1/15" not in row.get("content", "")
        or "without proving the physical field action" not in row.get("content", "")
        or "opposite rank-six sign" not in row.get("content", "")
    ):
        fail("FZ-12 registry content drifted from its bounded edge-ray claim")
    if (
        "Physical comparison is ineligible and unarmed"
        not in row.get("comparison_protocol", "")
        or "source-derived homogeneous position action"
        not in row.get("comparison_protocol", "")
        or "post-custody dataset-specific contract"
        not in row.get("comparison_protocol", "")
        or "#664" not in row.get("comparison_protocol", "")
        or "FZ-11" not in row.get("comparison_protocol", "")
    ):
        fail("FZ-12 registry comparison protocol is no longer fail-closed")
    if (
        any(f"FZ12-R0{index}" not in row.get("kill_band", "") for index in range(1, 6))
        or "append-only clarification" not in row.get("kill_band", "")
        or "including every a >= a_min" not in row.get("kill_band", "")
        or "Without that lower bound and power calculation" not in row.get(
            "kill_band", ""
        )
        or "OPH-wide only after a forced and exclusive physical bridge theorem"
        not in row.get("kill_band", "")
    ):
        fail("FZ-12 registry decision rule or scope of failure drifted")

    source = receipt.get("exact_source_result", {})
    geometry = source.get("finite_geometry_replay", {})
    if set(source) != {
        "finite_geometry_replay",
        "metric_completion",
        "scope_boundary",
        "seam_current_image",
    }:
        fail("FZ-12 exact source-result keys drifted")
    if (
        source.get("seam_current_image")
        != "D6 = {z in Z^6 : the coordinate sum is even}; the residual cokernel is one parity bit"
        or "same response-selected Euclidean three-carrier"
        not in str(source.get("metric_completion", ""))
        or "do not identify seam currents with physical motion"
        not in str(source.get("scope_boundary", ""))
        or geometry.get("source_ports") != 12
        or geometry.get("source_seams") != 30
        or geometry.get("signed_edge_directions") != 30
        or geometry.get("unoriented_axes") != 15
        or geometry.get("directed_seam_labels") != 60
        or geometry.get("port_degree") != 5
        or geometry.get("seam_difference_norm_squared") != "4+0*sqrt5"
        or geometry.get("edge_midpoint_norm_squared") != "6+2*sqrt5"
        or geometry.get("even_moments_on_unit_sphere")
        != {
            "sum_w_dot_n_squared": "10",
            "sum_w_dot_n_fourth": "6",
            "sum_w_dot_n_sixth": "30/7 - (2/7) I6(n)",
        }
    ):
        fail("FZ-12 exact finite source geometry or moments drifted")

    candidate = receipt.get("conditional_physical_candidate", {})
    if candidate.get("coefficients") != {
        "C4_over_a2": "-1/20",
        "B0_over_a4": "1/840",
        "B6_over_a4": "-1/12600",
    }:
        fail("FZ-12 edge coefficients drifted")
    if candidate.get("scale_free_relations") != {
        "B0_over_C4_squared": "10/21",
        "B6_over_C4_squared": "-2/63",
        "B6_over_B0": "-1/15",
    }:
        fail("FZ-12 scale-free coefficient manifold drifted")
    if candidate.get("signs") != {
        "B0": "positive",
        "B6": "negative",
        "C4": "negative",
    }:
        fail("FZ-12 coefficient signs drifted")
    if (
        candidate.get("operator")
        != "Λ_a(k,n) = (1/(5 a^2)) sum_{j=1}^{30} [1 - cos(a k w_j.n)]"
        or "spatial kinetic eigenvalue" not in str(candidate.get("symbol_name", ""))
        or "not automatically a frequency squared"
        not in str(candidate.get("symbol_name", ""))
        or "ranks one through five vanish"
        not in str(candidate.get("harmonic_nulls", ""))
        or "SO(3)/A5" not in str(candidate.get("fit_freedom_after_C4", ""))
    ):
        fail("FZ-12 operator, angular nulls, or orientation nuisance drifted")

    expected_premises = {
        "cofinal_gluing",
        "complete_support",
        "continuum_normalization",
        "equal_weights",
        "exclusivity",
        "finite_scale",
        "frame_and_boost",
        "homogeneous_translation_action",
        "physical_sector",
        "positive_scale_lower_bound",
        "readout_and_nuisance",
        "seam_as_displacement",
    }
    premises = receipt.get("physical_premises", {})
    if set(premises) != expected_premises:
        fail("FZ-12 must retain every named physical premise")
    if (
        "a >= a_min > 0" not in str(premises.get("positive_scale_lower_bound", ""))
        or "no experimental power against a null"
        not in str(premises.get("positive_scale_lower_bound", ""))
    ):
        fail("FZ-12 must retain its source-derived lower-bound gate")

    exposure = receipt.get("exposure_and_custody_boundary", {})
    if set(exposure) != {
        "candidate_registration",
        "comparison_data_read",
        "comparison_inputs",
        "comparison_permitted",
        "comparison_state",
        "excluded_data_class",
        "public_measurement_read",
        "source_inputs_only",
        "target_values_read",
    }:
        fail("FZ-12 exposure-boundary keys drifted")
    if (
        exposure.get("comparison_data_read") is not False
        or exposure.get("comparison_inputs") != []
        or exposure.get("comparison_permitted") is not False
        or exposure.get("public_measurement_read") is not False
        or exposure.get("target_values_read") is not False
        or exposure.get("comparison_state")
        != "INELIGIBLE_UNARMED_PHYSICAL_PRODUCER_OPEN"
        or "not the frozen prediction register"
        not in str(exposure.get("candidate_registration", ""))
        or "WMAP" not in str(exposure.get("excluded_data_class", ""))
        or "FZ-11" not in str(exposure.get("excluded_data_class", ""))
    ):
        fail("FZ-12 target exposure or comparison boundary drifted")

    producer = receipt.get("producer_scope", {})
    if (
        producer.get("frozen_prediction_registered") is not False
        or producer.get("physical_producer_closed") is not False
        or producer.get("type")
        != "target-clean prospective conditional physical-branch candidate"
        or "only if every physical promotion gate"
        not in str(producer.get("statement", ""))
    ):
        fail("FZ-12 source receipt must remain a pre-freeze physical candidate")

    promotion = receipt.get("promotion_gates", {})
    expected_gates = [
        ("seam-as-displacement identification", "#666"),
        ("homogeneous translation action", "#663/#666"),
        ("sole complete edge support and equal weights", "#655/#666"),
        ("physical scalar or polarization-independent sector", "#655/#666"),
        ("cofinal scale and observer gluing", "#663"),
        ("finite physical carrier scale", "#664"),
        ("source-derived positive physical scale lower bound", "#664"),
        ("frame, boost, readout, nuisance, and exclusivity contract", "#666"),
        ("dataset-specific post-custody preregistration", "#639"),
    ]
    actual_gates = promotion.get("gates")
    if (
        set(promotion) != {
            "all_discharged",
            "comparison_eligible",
            "fail_closed_rule",
            "gates",
            "physical_producer_closed",
        }
        or promotion.get("all_discharged") is not False
        or promotion.get("comparison_eligible") is not False
        or promotion.get("physical_producer_closed") is not False
        or not isinstance(actual_gates, list)
        or [
            (gate.get("gate"), gate.get("owner"), gate.get("status"))
            for gate in actual_gates
            if isinstance(gate, dict)
        ]
        != [(name, owner, "OPEN") for name, owner in expected_gates]
        or any(set(gate) != {"gate", "owner", "status"} for gate in actual_gates)
        or "no comparison may be armed"
        not in str(promotion.get("fail_closed_rule", ""))
    ):
        fail("FZ-12 physical promotion gates drifted or were promoted without proof")

    decision = receipt.get("prospective_decision_rule", {})
    if (
        set(decision) != {
            "eligibility",
            "fail",
            "inconclusive",
            "no_null_verdict",
            "scope_of_failure",
            "support",
            "trigger",
        }
        or "all physical promotion gates" not in str(decision.get("eligibility", ""))
        or "negative-B6 edge relation" not in str(decision.get("fail", ""))
        or "source-derived lower bound"
        not in str(decision.get("no_null_verdict", ""))
        or "forced and exclusive" not in str(decision.get("scope_of_failure", ""))
    ):
        fail("FZ-12 prospective decision rule drifted")

    separation = receipt.get("fz11_separation", {})
    if separation != {
        "edge_B6_over_C4_squared": "-2/63",
        "fz11_bytes_modified": False,
        "fz11_prediction_receipt_read": False,
        "fz11_register_id": "FZ-11",
        "opposite_rank_six_sign": True,
        "relationship": (
            "FZ-11 freezes the conditional primitive-vertex ray; this packet "
            "records the distinct source-native edge ray and cannot amend, "
            "reinterpret, or score FZ-11"
        ),
        "supersedes_fz11": False,
        "vertex_B6_over_C4_squared": "32/315",
    }:
        fail("FZ-12 must remain distinct from and must not supersede FZ-11")
    fz11 = rows_by_id.get("FZ-11", {})
    if (
        fz11.get("status") not in FROZEN_STATUSES
        or "B6/C4^2=32/315" not in str(fz11.get("content", ""))
        or fz11.get("content_sha256") == row.get("content_sha256")
    ):
        fail("FZ-11 and FZ-12 must remain separately frozen coefficient rays")

    pins = receipt.get("parent_pins")
    if not isinstance(pins, list) or len(pins) != 4:
        fail("FZ-12 must retain exactly four source parent pins")
    pins_by_path = {
        pin.get("path"): pin for pin in pins if isinstance(pin, dict)
    }
    expected_pin_paths = {
        "code/a5_fingerprint/runtime/a5_multipole_fixed_point_receipt.json",
        FZ12_CARRIER_LEAN_REL,
        FZ12_MOMENT_LEAN_REL,
        FZ12_RAY_LEAN_REL,
    }
    if set(pins_by_path) != expected_pin_paths:
        fail("FZ-12 source parent-pin set drifted")
    for relative_path, pin in pins_by_path.items():
        if set(pin) != {"bytes", "path", "role", "sha256"}:
            fail(f"FZ-12 parent pin {relative_path} has unexpected keys")
        source_path = ROOT / relative_path
        expected_hash = "sha256:" + sha256_file(source_path)
        if pin.get("sha256") != expected_hash or pin.get("bytes") != source_path.stat().st_size:
            fail(f"FZ-12 parent pin differs from live source bytes: {relative_path}")
        if relative_path in FZ12_IN_REPO_ARTIFACTS and (
            contract["in_repo_artifact_sha256"].get(relative_path)
            != expected_hash.removeprefix("sha256:")
        ):
            fail(f"FZ-12 custody contract does not bind parent pin {relative_path}")


def validate_retrospective_results(register: dict) -> set[str]:
    results = register.get("retrospective_results")
    if not isinstance(results, list) or not results:
        fail("retrospective_results must be a nonempty list")

    seen_ids: set[str] = set()
    former_reservations: set[str] = set()
    for index, result in enumerate(results):
        where = f"retrospective_results[{index}] ({result.get('id')})"
        if not isinstance(result, dict) or set(result) != RETROSPECTIVE_RESULT_KEYS:
            fail(f"{where}: keys mismatch")
        result_id = result["id"]
        if (
            not isinstance(result_id, str)
            or not result_id
            or result_id in seen_ids
        ):
            fail(f"{where}: id must be a unique nonempty string")
        seen_ids.add(result_id)
        former = result["former_ladder_reservation"]
        if (
            not isinstance(former, str)
            or re.fullmatch(r"FZ-\d{2}", former) is None
            or former in former_reservations
        ):
            fail(f"{where}: former_ladder_reservation must be a unique FZ id")
        former_reservations.add(former)
        if result["status"] != "retrospective_not_evaluable":
            fail(f"{where}: unsupported retrospective status")
        for key in (
            "content",
            "payload_path",
            "comparison_protocol",
            "evidential_boundary",
            "milestone",
        ):
            if not isinstance(result[key], str) or not result[key].strip():
                fail(f"{where}: {key} must be nonempty")
        payload_path = Path(result["payload_path"])
        if payload_path.is_absolute() or ".." in payload_path.parts:
            fail(f"{where}: payload_path must be a repository-relative path")
        require_sha256(result["payload_sha256"], f"{where}.payload_sha256")
        if not isinstance(result["owning_issue"], int):
            fail(f"{where}: owning_issue must identify the closed source issue")

    fz04 = results[0]
    if (
        len(results) != 1
        or fz04["id"] != "RR-506-ALPHA-HVP"
        or fz04["former_ladder_reservation"] != "FZ-04"
        or fz04["owning_issue"] != 506
        or fz04["payload_path"] != FZ04_VERDICT_REL
    ):
        fail("the issue-506 retrospective result binding is malformed")

    verdict = load_json(FZ04_VERDICT_PATH)
    if (
        verdict.get("schema") != "oph.alpha_hvp_class_verdict.v2"
        or verdict.get("issue") != 506
        or verdict.get("row_class")
        != "retrospective_empirical_same_scheme_accounting_audit"
        or verdict.get("verdict")
        != "MULTI_CLASS_NOT_EVALUABLE__ONE_RECORDED_ACCOUNTING_REPLAY_COMPATIBLE"
    ):
        fail("the issue-506 payload has the wrong retrospective verdict identity")
    if verdict.get("scope") != FZ04_SCOPE:
        fail("the issue-506 payload scope differs from the bounded retrospective scope")
    if verdict.get("claim_boundary") != FZ04_CLAIM_BOUNDARY:
        fail("the issue-506 payload claim boundary differs from the bounded statement")

    payload_without_digest = {
        key: value for key, value in verdict.items() if key != "verdict_sha256"
    }
    computed_hash = sha256_bytes(canonical_json_bytes(payload_without_digest))
    reported_hash = verdict.get("verdict_sha256")
    if reported_hash != f"sha256:{computed_hash}":
        fail(
            "the issue-506 payload self-digest does not equal its canonical "
            f"content hash {computed_hash}"
        )
    if fz04["payload_sha256"] != computed_hash:
        fail(
            "the issue-506 retrospective payload hash does not equal the "
            f"canonical content hash {computed_hash}"
        )
    rebuilt_verdict = rebuild_issue506_verdict()
    if verdict != rebuilt_verdict:
        fail(
            "the issue-506 stored payload does not equal the canonical producer replay"
        )
    return former_reservations


def validate_frozen_rows_bootstrap(rows: list[dict]) -> None:
    """Bind frozen and void rows to their immutable V3-bootstrap payload.

    Pending and resource-deferred rows remain editable until a freeze occurs.
    The October 4 append-only extension admits exactly FZ-13 and FZ-15 under
    separate immutable custody pins. It does not repin or rewrite the original
    bootstrap. A new freeze must obtain its own reviewed certificate; changing
    status or content cannot silently add it to the historical population.
    """

    where = "frozen-row V3 bootstrap"
    relative = REGISTER_PATH.relative_to(ROOT).as_posix()
    resolved = _run_git(
        ROOT,
        ["rev-parse", f"{FROZEN_ROWS_BOOTSTRAP_REVISION}:{relative}"],
        where,
    )
    if resolved.returncode != 0:
        fail(f"{where}: pinned register Git object is unavailable")
    actual_blob = resolved.stdout.decode("ascii", errors="replace").strip()
    if actual_blob != FROZEN_ROWS_BOOTSTRAP_GIT_BLOB_SHA1:
        fail(f"{where}: pinned register Git blob changed")
    raw = read_commit_blob(ROOT, FROZEN_ROWS_BOOTSTRAP_REVISION, relative, where)
    try:
        historical = strict_json.loads(raw.decode("utf-8"))
    except (
        UnicodeDecodeError,
        json.JSONDecodeError,
        strict_json.DuplicateKeyError,
    ) as error:
        fail(f"{where}: pinned register is invalid strict JSON: {error}")
    historical_rows = historical.get("rows") if isinstance(historical, dict) else None
    if not isinstance(historical_rows, list):
        fail(f"{where}: pinned register has no row collection")
    immutable_statuses = FROZEN_STATUSES | {"superseded_void"}
    historical_immutable = [
        row for row in historical_rows if row.get("status") in immutable_statuses
    ]
    # The fixed normalization changes display custody metadata only. Every
    # scientific field continues to compare against the original Git blob.
    metadata_rows = custody_metadata_normalization()["rows"]
    for row in historical_immutable:
        replacements = metadata_rows.get(row["id"], {})
        if not set(replacements) <= {"custody", "attestation"}:
            fail("custody normalization cannot change scientific row fields")
        row.update(replacements)
    historical_ids = {row["id"] for row in historical_immutable}
    current_immutable = [row for row in rows if row.get("id") in historical_ids]
    if historical_immutable != current_immutable:
        fail(
            f"{where}: frozen or void rows differ from the immutable bootstrap; "
            "the append-only extension cannot rewrite historical frozen rows"
        )
    added_immutable = {
        row["id"] for row in rows
        if row.get("status") in immutable_statuses and row["id"] not in historical_ids
    }
    if added_immutable != set(COSMOLOGY_FREEZE_PINS):
        fail("frozen-row extension must contain exactly the registered FZ-13 and FZ-15 freezes")


def validate_cosmology_freezes(register: dict, rows_by_id: dict[str, dict]) -> None:
    """Bind the new row and complete protocol to the accepted custody bytes.

    A content hash alone cannot protect a comparison rule. The full snapshot
    contains the complete register row and structured protocol, and its fixed
    digest prevents a changed rule from being concealed by rebuilding custody.
    The independently fixed manifest digest also protects the human target.
    """

    for row_id, pins in COSMOLOGY_FREEZE_PINS.items():
        where = f"{row_id} conditional cosmology freeze"
        contract = register["external_custody_contracts"][row_id]
        expected_path = f"falsification/frozen_targets/{row_id.lower().replace('-', '')}_2026-10-04"
        if (
            contract["rows"] != [row_id]
            or contract["custody_path"] != expected_path
            or contract["registration_manifest"] != "registration_manifest_2026-10-04.json"
        ):
            fail(f"{where}: custody identity drifted")
        if contract["registration_manifest_sha256"] != pins["manifest_sha256"]:
            fail(f"{where}: immutable manifest pin drifted")
        if contract["artifact_sha256"].get("frozen_prediction.json") != pins["snapshot_sha256"]:
            fail(f"{where}: immutable snapshot pin drifted")
        snapshot_path = DEFAULT_CUSTODY_ROOT / expected_path / "frozen_prediction.json"
        if not snapshot_path.is_file() or sha256_file(snapshot_path) != pins["snapshot_sha256"]:
            fail(f"{where}: immutable snapshot bytes drifted")
        snapshot = load_json(snapshot_path)
        if set(snapshot) != {"schema", "row", "protocol"} or snapshot["schema"] != COSMOLOGY_FREEZE_SCHEMA:
            fail(f"{where}: snapshot schema drifted")
        if snapshot["row"] != rows_by_id[row_id]:
            fail(f"{where}: register row differs from immutable snapshot")
        row = rows_by_id[row_id]
        if row["content_sha256"] != sha256_bytes(row["content"].encode("utf-8")):
            fail(f"{where}: content hash does not bind the exact content text")
        if not isinstance(snapshot["protocol"], dict) or not snapshot["protocol"]:
            fail(f"{where}: complete structured protocol is required")
        if "OWNER SLOT" in json.dumps(snapshot["protocol"]):
            fail(f"{where}: unresolved protocol owner slot")
        verify_cosmology_source_artifacts(snapshot["protocol"], where)


def verify_cosmology_source_artifacts(protocol: dict, where: str) -> None:
    """Verify the freeze's historical source bytes without pinning live files."""

    commit = protocol.get("source_commit")
    if not isinstance(commit, str) or re.fullmatch(r"[0-9a-f]{40}", commit) is None:
        fail(f"{where}: source_commit must be a full historical commit")
    pins = validate_hash_mapping(protocol.get("source_artifacts"), f"{where}.source_artifacts")
    if not pins:
        fail(f"{where}: source artifacts are required")
    for relative, digest in pins.items():
        if ".." in Path(relative).parts:
            fail(f"{where}: source artifact path must stay inside the repository")
        payload = read_commit_blob(ROOT, commit, relative, where)
        if sha256_bytes(payload) != digest:
            fail(f"{where}: historical source artifact hash mismatch: {relative}")


def validate(register: dict) -> list[dict]:
    if set(register) != REGISTER_KEYS:
        fail("top-level keys mismatch")
    if register.get("schema") != SCHEMA:
        fail(f"schema must equal {SCHEMA}")
    if register.get("issue") != 607:
        fail("issue must equal 607")
    rows = register.get("rows")
    if not isinstance(rows, list) or not rows:
        fail("rows must be a nonempty list")

    seen_ids = [row.get("id") for row in rows]
    if (
        any(not isinstance(row_id, str) for row_id in seen_ids)
        or len(set(seen_ids)) != len(seen_ids)
        or any(re.fullmatch(r"FZ-\d{2}", row_id) is None for row_id in seen_ids)
        or seen_ids != sorted(seen_ids)
    ):
        fail("rows must carry unique ascending FZ identifiers")

    former_reservations = validate_retrospective_results(register)
    overlap = set(seen_ids) & former_reservations
    if overlap:
        fail(
            "retrospective reservations must not appear as prospective ladder "
            f"rows: {sorted(overlap)}"
        )
    allocated = sorted(set(seen_ids) | former_reservations)
    expected_allocated = [f"FZ-{index:02d}" for index in range(1, 16)]
    if allocated != expected_allocated:
        fail(
            "ladder rows and explicitly retired reservations must account for "
            f"{expected_allocated}, got {allocated}"
        )

    rows_by_id: dict[str, dict] = {}
    for index, row in enumerate(rows):
        where = f"rows[{index}] ({row.get('id')})"
        if set(row) != ROW_KEYS:
            fail(f"{where}: keys mismatch")
        if row["status"] not in STATUSES:
            fail(f"{where}: unknown status {row['status']}")
        owning_issue = row["owning_issue"]
        if owning_issue is not None and (
            not isinstance(owning_issue, int)
            or isinstance(owning_issue, bool)
            or owning_issue <= 0
        ):
            fail(f"{where}: owning_issue must be null or a positive integer")
        for key in ("content", "kill_band", "comparison_protocol"):
            if not isinstance(row[key], str) or not row[key].strip():
                fail(f"{where}: {key} must be nonempty")
        if row["status"] in FROZEN_STATUSES:
            for key in ("frozen_utc", "custody", "attestation", "content_sha256"):
                if not isinstance(row[key], str) or not row[key].strip():
                    fail(f"{where}: a frozen row requires {key}")
            parse_utc(row["frozen_utc"], f"{where}.frozen_utc")
        elif row["status"] == "registered_pending_freeze":
            if owning_issue is None:
                fail(f"{where}: a pending row requires a current owning issue")
            if not is_open(owning_issue):
                fail(
                    f"{where}: a pending row's owning issue #{owning_issue} must be"
                    " one of the open issues in tools/closed_lanes.py"
                )
            if not isinstance(row["milestone"], str) or not row["milestone"].strip():
                fail(f"{where}: a pending row requires a milestone")
        elif row["status"] == "resource_deferred":
            if row["owning_issue"] is not None:
                fail(f"{where}: a resource-deferred row cannot retain an open owner")
            if not isinstance(row["milestone"], str) or not row["milestone"].strip():
                fail(f"{where}: a resource-deferred row requires a disposition")
        else:
            if row["status"] != "superseded_void":
                raise AssertionError("unreachable status branch")
            if row["frozen_utc"] is not None:
                fail(f"{where}: a superseded-void row cannot retain a freeze time")
            if row["owning_issue"] is not None:
                fail(f"{where}: a superseded-void row cannot retain an open owner")
            if row["milestone"] != "superseded":
                fail(f"{where}: a superseded-void row requires milestone superseded")
            if row["comparison_protocol"].lower().startswith("none;") is False:
                fail(f"{where}: a superseded-void row must refuse comparison")
            if row["kill_band"].lower().startswith("none;") is False:
                fail(f"{where}: a superseded-void row must carry no kill band")
        rows_by_id[row["id"]] = row

    fz01 = rows_by_id.get("FZ-01", {})
    if (
        "four retained frozen targets" not in str(fz01.get("content", ""))
        or "ringdown row is excluded from scoring"
        not in str(fz01.get("comparison_protocol", ""))
    ):
        fail("FZ-01 must exclude the superseded ringdown row from its retained targets")
    fz06 = rows_by_id.get("FZ-06", {})
    if (
        fz06.get("status") != "superseded_void"
        or fz06.get("frozen_utc") is not None
        or "VOID/MISATTRIBUTED" not in str(fz06.get("content", ""))
        or "no frozen numeric prediction" not in str(fz06.get("content", ""))
    ):
        fail("FZ-06 must retain the alpha=4 record only as superseded void history")

    contracts = validate_custody_contracts(register, rows_by_id)
    verify_fz11_source_history(contracts["FZ-11"])
    validate_fz11_prediction(rows_by_id, contracts["FZ-11"])
    verify_fz12_source_history(contracts["FZ-12"])
    validate_fz12_prediction(rows_by_id, contracts["FZ-12"])

    fz02 = rows_by_id.get("FZ-02", {})
    if (
        "frame-lock clause is not established" not in fz02.get("content", "")
        or "ineligible pending issue #643" not in fz02.get("content", "")
        or "no frame-lock verdict may be issued" not in fz02.get(
            "comparison_protocol", ""
        )
        or "FZ02-R03a" not in fz02.get("kill_band", "")
        or "FZ02-R03b (INELIGIBLE; #643)" not in fz02.get("kill_band", "")
    ):
        fail(
            "FZ-02 must keep the unsupported frame-lock clause ineligible "
            "pending issue #643"
        )

    fz03 = rows_by_id.get("FZ-03", {})
    if (
        fz03.get("status") != "registered_pending_freeze"
        or fz03.get("owning_issue") != PENDING_TEST_OWNER
    ):
        fail(f"FZ-03 must remain pending under the open test owner issue #{PENDING_TEST_OWNER}")

    fz05 = rows_by_id.get("FZ-05", {})
    if (
        "#736" not in fz05.get("content", "")
        or "#729" not in fz05.get("content", "")
        or "retrospective" not in fz05.get("comparison_protocol", "").lower()
        or not all(
            issue in fz05.get("kill_band", "")
            for issue in ("#729", "#736", f"#{COMMON_MODEL_OWNER}", f"#{PENDING_TEST_OWNER}")
        )
        or "historical lane" not in fz05.get("content", "")
        or fz05.get("owning_issue") != PENDING_TEST_OWNER
    ):
        fail(
            "FZ-05 must keep the capacity, common-tower, comparison, and "
            "retrospective exposure boundaries"
        )
    receipt = load_json(FZ02_RECEIPT_PATH)
    live_hash = receipt.get("receipt_sha256")
    if fz02["content_sha256"] != live_hash:
        fail(
            "the FZ-02 content hash does not equal the live angular-multiplet "
            f"receipt hash {live_hash}"
        )
    validate_frozen_rows_bootstrap(rows)
    validate_cosmology_freezes(register, rows_by_id)
    return rows


def ots_binding_bytes(data: bytes, where: str) -> dict[str, Any]:
    prefix = OTS_DETACHED_HEADER + OTS_SHA256_TAG
    if not data.startswith(prefix) or len(data) < len(prefix) + 32:
        fail(f"{where}: unsupported or malformed detached OpenTimestamps proof")
    digest_start = len(prefix)
    return {
        "file_sha256": data[digest_start : digest_start + 32].hex(),
        "has_pending_calendar": OTS_PENDING_ATTESTATION_TAG in data,
        "has_bitcoin_attestation": OTS_BITCOIN_ATTESTATION_TAG in data,
    }


def ots_binding(path: Path) -> dict[str, Any]:
    return ots_binding_bytes(path.read_bytes(), str(path))


def verify_ots_with_official_cli(path: Path) -> str:
    """Parse a detached stamp with ``ots info`` when the client is usable.

    ``ots info`` is entirely local and does not contact calendars or Bitcoin
    nodes. The manual digest and attestation checks remain mandatory because
    the command's human-readable output is not a stable machine interface.
    """

    global _OTS_CLI_PATH, _OTS_CLI_USABLE

    proof_digest = sha256_file(path)
    if proof_digest in _OTS_OFFICIAL_PARSED_DIGESTS:
        return "verified_by_ots_info"

    if _OTS_CLI_USABLE is False:
        return "official_cli_not_available"
    candidates: list[Path] = []
    if _OTS_CLI_PATH is not None:
        candidates.append(Path(_OTS_CLI_PATH))
    discovered = shutil.which("ots")
    if discovered is not None:
        candidates.append(Path(discovered))
    candidates.append(Path(sys.executable).with_name("ots"))
    pyenv_root = Path(os.environ.get("PYENV_ROOT", Path.home() / ".pyenv"))
    candidates.append(pyenv_root / "versions" / "sherlock2" / "bin" / "ots")
    seen: set[Path] = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        if not candidate.is_file() or not os.access(candidate, os.X_OK):
            continue
        executable = str(candidate)
        try:
            parsed = subprocess.run(
                [executable, "info", str(path)],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        except OSError:
            continue
        combined = f"{parsed.stdout}\n{parsed.stderr}".lower()
        if parsed.returncode == 0:
            _OTS_CLI_PATH = executable
            _OTS_CLI_USABLE = True
            _OTS_OFFICIAL_PARSED_DIGESTS.add(proof_digest)
            return "verified_by_ots_info"
        # A pyenv/asdf shim can resolve even when the selected environment does
        # not contain the client. Continue to a concrete environment candidate.
        if parsed.returncode in {126, 127} and (
            "command not found" in combined
            or "no such file" in combined
            or "not installed" in combined
        ):
            continue
        fail(
            f"official OpenTimestamps parser rejected {path}: "
            f"exit {parsed.returncode}: {combined[-1000:]}"
        )
    _OTS_CLI_USABLE = False
    return "official_cli_not_available"


def verify_ots_binding(
    path: Path, expected_hash: str, expected_state: str
) -> dict[str, Any]:
    if not path.is_file():
        fail(f"missing detached OpenTimestamps proof {path}")
    official_parser = verify_ots_with_official_cli(path)
    if official_parser != "verified_by_ots_info":
        fail(
            "official OpenTimestamps parser is required for external custody "
            f"verification: {path}"
        )
    proof = ots_binding(path)
    if proof["file_sha256"] != expected_hash:
        fail(f"OpenTimestamps digest mismatch for {path}")
    if expected_state == "calendar_pending":
        if not proof["has_pending_calendar"] or proof["has_bitcoin_attestation"]:
            fail(f"{path} must be calendar-pending and not yet Bitcoin-attested")
    elif expected_state == "bitcoin_attested":
        if not proof["has_bitcoin_attestation"]:
            fail(f"{path} does not contain a Bitcoin block attestation")
    else:
        fail(f"unknown expected OpenTimestamps state {expected_state}")
    return {**proof, "official_parser": official_parser}


def verify_historical_ots_binding(
    payload: bytes,
    expected_hash: str,
    expected_state: str,
    where: str,
) -> dict[str, Any]:
    """Parse one proof exactly as stored in a historical Git commit.

    Historical proof bytes are intentionally independent from the current
    working-copy proof, which may receive an in-place Bitcoin attestation
    upgrade. The artifact digest and attestation state remain checked at the
    historical commit.
    """

    with tempfile.TemporaryDirectory(prefix="oph-fz11-ots-") as temp_dir:
        proof_path = Path(temp_dir) / "historical.ots"
        proof_path.write_bytes(payload)
        official_parser = verify_ots_with_official_cli(proof_path)
        if official_parser != "verified_by_ots_info":
            fail(
                "official OpenTimestamps parser is required for historical "
                f"custody verification: {where}"
            )
    proof = ots_binding_bytes(payload, where)
    if proof["file_sha256"] != expected_hash:
        fail(f"OpenTimestamps digest mismatch for {where}")
    if expected_state == "calendar_pending":
        if not proof["has_pending_calendar"] or proof["has_bitcoin_attestation"]:
            fail(f"{where} must be calendar-pending and not Bitcoin-attested")
    elif expected_state == "bitcoin_attested":
        if not proof["has_bitcoin_attestation"]:
            fail(f"{where} does not contain a Bitcoin block attestation")
    else:
        fail(f"unknown expected OpenTimestamps state {expected_state}")
    return {**proof, "official_parser": official_parser}


def verify_commit_ots_bindings(
    repo_root: Path,
    commit: str,
    expected_artifacts: dict[str, str],
    expected_state: str,
    where: str,
) -> int:
    """Verify historical detached proofs without pinning them to live bytes."""

    for artifact_path, digest in expected_artifacts.items():
        proof_path = artifact_path + ".ots"
        payload = read_commit_blob(repo_root, commit, proof_path, where)
        verify_historical_ots_binding(
            payload,
            digest,
            expected_state,
            f"{where}:{proof_path}",
        )
    return len(expected_artifacts)


def fenced_target_hashes(path: Path) -> tuple[str, str]:
    payload = path.read_bytes()
    start_marker = b"<!-- FZ02-TARGET-BEGIN -->"
    end_marker = b"<!-- FZ02-TARGET-END -->"
    try:
        start = payload.index(start_marker)
        start_line_end = payload.index(b"\n", start) + 1
        end = payload.index(end_marker, start_line_end)
    except ValueError:
        fail(f"{path}: missing FZ-02 target fence")
    end_with_marker = end + len(end_marker)
    if payload[end_with_marker : end_with_marker + 2] == b"\r\n":
        end_with_marker += 2
    elif payload[end_with_marker : end_with_marker + 1] == b"\n":
        end_with_marker += 1
    return (
        sha256_bytes(payload[start:end_with_marker]),
        sha256_bytes(payload[start_line_end:end]),
    )


def verify_external_custody(
    register: dict,
    custody_root: Path = DEFAULT_CUSTODY_ROOT,
) -> dict[str, Any]:
    contracts = register["external_custody_contracts"]
    custody_base = custody_root / "falsification" / "frozen_targets"
    unverified_state = None
    if not custody_base.is_dir():
        unverified_state = "external_custody_not_present"
    elif CUSTODY_PARSER_UNSUPPORTED:
        unverified_state = "custody_parser_unsupported_on_platform"
    if unverified_state is not None:
        return {
            "state": unverified_state,
            "custody_root": str(custody_root),
            "fz11_custody_git_history": {"state": unverified_state},
            "fz12_custody_git_history": {"state": unverified_state},
            "contracts": {
                contract_id: {
                    "verification": unverified_state,
                    "attestation_state": contract["attestation_state"],
                }
                for contract_id, contract in contracts.items()
            },
        }

    results: dict[str, dict[str, str]] = {}
    fz11_custody_git_history: dict[str, Any] = {
        "state": "external_custody_not_present"
    }
    fz12_custody_git_history: dict[str, Any] = {
        "state": "external_custody_not_present"
    }
    for contract_id, contract in contracts.items():
        directory = custody_root / contract["custody_path"]
        if not directory.is_dir():
            fail(f"{contract_id} external custody directory is missing: {directory}")
        manifest_path = directory / contract["registration_manifest"]
        if not manifest_path.is_file():
            fail(f"{contract_id} registration manifest is missing")
        if sha256_file(manifest_path) != contract["registration_manifest_sha256"]:
            fail(f"{contract_id} registration manifest hash mismatch")
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as error:
            fail(f"{contract_id} registration manifest is invalid JSON: {error}")
        if manifest.get("artifacts") != contract["artifact_sha256"]:
            fail(
                f"{contract_id} registration manifest differs from the in-repo contract"
            )
        if contract_id == "FZ-10":
            if (
                manifest.get("registration") != FZ10_REGISTRATION
                or manifest.get("source_commit") != FZ10_SOURCE_COMMIT
            ):
                fail("FZ-10 registration manifest bindings drifted")
        if contract_id == "FZ-11":
            if set(manifest) != {
                "registration",
                "frozen_utc",
                "source_commit",
                "policy",
                "artifacts",
                "in_repo_artifacts",
            }:
                fail("FZ-11 registration manifest keys drifted")
            if (
                manifest.get("source_commit") != contract["source_commit"]
                or manifest.get("frozen_utc") != contract["frozen_utc"]
                or manifest.get("in_repo_artifacts")
                != contract["original_in_repo_artifact_sha256"]
                or "prediction frozen before" not in str(manifest.get("policy", ""))
            ):
                fail("FZ-11 registration manifest bindings drifted")
        if contract_id == "FZ-12":
            if set(manifest) != {
                "registration",
                "frozen_utc",
                "source_commit",
                "policy",
                "artifacts",
                "in_repo_artifacts",
            }:
                fail("FZ-12 registration manifest keys drifted")
            if (
                manifest.get("registration")
                != "FZ-12 source-seam edge propagation branch"
                or manifest.get("source_commit") != contract["source_commit"]
                or manifest.get("frozen_utc") != contract["frozen_utc"]
                or manifest.get("in_repo_artifacts")
                != contract["in_repo_artifact_sha256"]
                or "before any new target or comparison data are examined"
                not in str(manifest.get("policy", ""))
                or "FZ-11" not in str(manifest.get("policy", ""))
                or "physical producer" not in str(manifest.get("policy", ""))
                or "unarmed" not in str(manifest.get("policy", ""))
            ):
                fail("FZ-12 registration manifest bindings drifted")

        state = contract["attestation_state"]
        for filename, expected_hash in contract["artifact_sha256"].items():
            artifact_path = directory / filename
            if not artifact_path.is_file():
                fail(f"{contract_id} custody artifact is missing: {filename}")
            if not verify_custody_artifact_digest(artifact_path, expected_hash):
                fail(f"{contract_id} custody artifact hash mismatch: {filename}")
            verify_ots_binding(Path(str(artifact_path) + ".ots"), expected_hash, state)
        verify_ots_binding(
            Path(str(manifest_path) + ".ots"),
            contract["registration_manifest_sha256"],
            state,
        )

        if contract_id == "FZ-02":
            erratum_path = directory / contract["custody_erratum"]
            if not erratum_path.is_file():
                fail("FZ-02 append-only custody erratum is missing")
            if not verify_custody_artifact_digest(erratum_path, contract["custody_erratum_sha256"]):
                fail("FZ-02 custody erratum hash mismatch")
            scientific_erratum_path = directory / contract["scientific_erratum"]
            if not scientific_erratum_path.is_file():
                fail("FZ-02 append-only scientific erratum is missing")
            if (
                sha256_file(scientific_erratum_path)
                != contract["scientific_erratum_sha256"]
            ):
                fail("FZ-02 scientific erratum hash mismatch")
            target_path = directory / contract["target_file"]
            block_hash, payload_hash = fenced_target_hashes(target_path)
            if block_hash != contract["target_block_sha256"]:
                fail("FZ-02 fenced target block hash mismatch")
            if payload_hash != contract["target_payload_sha256"]:
                fail("FZ-02 fenced target payload hash mismatch")

        if contract_id == "FZ-11":
            prediction_path = directory / contract["prediction_file"]
            try:
                prediction = json.loads(prediction_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError as error:
                fail(f"FZ-11 frozen prediction snapshot is invalid JSON: {error}")
            if set(prediction) != {
                "schema",
                "frozen_utc",
                "source_commit",
                "source_receipt",
                "source_receipt_sha256",
                "status",
                "operator",
                "coefficients",
                "scale_free_relations",
                "harmonic_nulls",
                "shape",
                "comparison_state",
                "new_comparison_data_read",
                "excluded_prior_data",
                "scope_of_failure",
            }:
                fail("FZ-11 frozen prediction snapshot keys drifted")
            if (
                prediction.get("schema")
                != "oph.fz11.primitive_port_prediction.freeze.v1"
                or prediction.get("frozen_utc") != contract["frozen_utc"]
                or prediction.get("source_commit") != contract["source_commit"]
                or prediction.get("source_receipt") != FZ11_RECEIPT_REL
                or prediction.get("source_receipt_sha256")
                != contract["in_repo_artifact_sha256"][FZ11_RECEIPT_REL]
                or prediction.get("status") != FZ11_STATUS
                or prediction.get("coefficients")
                != {
                    "C4_over_a2": "-1/20",
                    "B0_over_a4": "1/840",
                    "B6_over_a4": "2/7875",
                }
                or prediction.get("scale_free_relations")
                != {
                    "B6_over_C4_squared": "32/315",
                    "B0_over_C4_squared": "10/21",
                    "B6_over_B0": "16/75",
                }
                or prediction.get("harmonic_nulls") != [1, 2, 3, 4, 5]
                or prediction.get("new_comparison_data_read") is not False
                or prediction.get("comparison_state")
                != "UNARMED_PENDING_DATASET_SPECIFIC_PREREGISTRATION"
            ):
                fail("FZ-11 frozen prediction snapshot content drifted")

            repair_manifest_path = directory / contract["lean_repair_manifest"]
            if not repair_manifest_path.is_file():
                fail("FZ-11 Lean repair manifest is missing")
            if (
                sha256_file(repair_manifest_path)
                != contract["lean_repair_manifest_sha256"]
            ):
                fail("FZ-11 Lean repair manifest hash mismatch")
            try:
                repair_manifest = json.loads(
                    repair_manifest_path.read_text(encoding="utf-8")
                )
            except json.JSONDecodeError as error:
                fail(f"FZ-11 Lean repair manifest is invalid JSON: {error}")
            if set(repair_manifest) != {
                "schema",
                "repair_utc",
                "original_prediction_frozen_utc",
                "original_source_commit",
                "repair_source_commit",
                "original_registration_manifest",
                "original_registration_manifest_sha256",
                "prediction_receipt_sha256",
                "prediction_bytes_changed",
                "new_comparison_data_read",
                "artifacts",
                "live_in_repo_artifact",
                "live_in_repo_artifact_sha256",
                "verification",
            }:
                fail("FZ-11 Lean repair manifest keys drifted")
            if (
                repair_manifest.get("schema")
                != "oph.fz11_lean_repair_manifest.v1"
                or repair_manifest.get("repair_utc")
                != contract["lean_repair_utc"]
                or repair_manifest.get("original_prediction_frozen_utc")
                != contract["frozen_utc"]
                or repair_manifest.get("original_source_commit")
                != contract["source_commit"]
                or repair_manifest.get("repair_source_commit")
                != contract["lean_repair_source_commit"]
                or repair_manifest.get("original_registration_manifest")
                != contract["registration_manifest"]
                or repair_manifest.get("original_registration_manifest_sha256")
                != contract["registration_manifest_sha256"]
                or repair_manifest.get("prediction_receipt_sha256")
                != contract["in_repo_artifact_sha256"][FZ11_RECEIPT_REL]
                or repair_manifest.get("prediction_bytes_changed") is not False
                or repair_manifest.get("new_comparison_data_read") is not False
                or repair_manifest.get("artifacts")
                != contract["lean_repair_artifact_sha256"]
                or repair_manifest.get("live_in_repo_artifact") != FZ11_LEAN_REL
                or repair_manifest.get("live_in_repo_artifact_sha256")
                != contract["in_repo_artifact_sha256"][FZ11_LEAN_REL]
                or repair_manifest.get("verification")
                != {
                    "lean_exit_code": 0,
                    "sorry_ax_present": False,
                    "public_theorem_count": 5,
                }
            ):
                fail("FZ-11 Lean repair manifest bindings drifted")
            repair_state = contract["lean_repair_attestation_state"]
            for filename, expected_hash in contract[
                "lean_repair_artifact_sha256"
            ].items():
                repair_artifact = directory / filename
                if not repair_artifact.is_file():
                    fail(f"FZ-11 Lean repair artifact is missing: {filename}")
                if sha256_file(repair_artifact) != expected_hash:
                    fail(f"FZ-11 Lean repair artifact hash mismatch: {filename}")
                verify_ots_binding(
                    Path(str(repair_artifact) + ".ots"),
                    expected_hash,
                    repair_state,
                )
            verify_ots_binding(
                Path(str(repair_manifest_path) + ".ots"),
                contract["lean_repair_manifest_sha256"],
                repair_state,
            )
            repaired_copy = directory / "A5PrimitivePortPrediction_repaired.lean"
            if sha256_file(repaired_copy) != contract["in_repo_artifact_sha256"][FZ11_LEAN_REL]:
                fail("FZ-11 repaired Lean custody copy differs from the live proof")

            decision_manifest_path = directory / contract["decision_rule_manifest"]
            if not decision_manifest_path.is_file():
                fail("FZ-11 decision-rule correction manifest is missing")
            if (
                sha256_file(decision_manifest_path)
                != contract["decision_rule_manifest_sha256"]
            ):
                fail("FZ-11 decision-rule correction manifest hash mismatch")
            try:
                decision_manifest = json.loads(
                    decision_manifest_path.read_text(encoding="utf-8")
                )
            except json.JSONDecodeError as error:
                fail(f"FZ-11 decision-rule manifest is invalid JSON: {error}")
            if set(decision_manifest) != {
                "schema",
                "created_utc",
                "original_prediction_frozen_utc",
                "original_source_commit",
                "original_custody_commit",
                "prediction_bytes_changed",
                "comparison_data_read",
                "original_bindings",
                "append_only_artifacts",
            }:
                fail("FZ-11 decision-rule correction manifest keys drifted")
            expected_original_bindings = {
                contract["target_file"]: contract["artifact_sha256"][
                    contract["target_file"]
                ],
                contract["prediction_file"]: contract["artifact_sha256"][
                    contract["prediction_file"]
                ],
                contract["registration_manifest"]: contract[
                    "registration_manifest_sha256"
                ],
                "source_prediction_receipt": contract["in_repo_artifact_sha256"][
                    FZ11_RECEIPT_REL
                ],
            }
            if (
                decision_manifest.get("schema")
                != "oph.fz11.decision_rule_erratum_manifest.v1"
                or decision_manifest.get("created_utc")
                != contract["decision_rule_utc"]
                or decision_manifest.get("original_prediction_frozen_utc")
                != contract["frozen_utc"]
                or decision_manifest.get("original_source_commit")
                != contract["source_commit"]
                or decision_manifest.get("original_custody_commit")
                != contract["custody_commit"]
                or decision_manifest.get("prediction_bytes_changed") is not False
                or decision_manifest.get("comparison_data_read") is not False
                or decision_manifest.get("original_bindings")
                != expected_original_bindings
                or decision_manifest.get("append_only_artifacts")
                != contract["decision_rule_artifact_sha256"]
            ):
                fail("FZ-11 decision-rule correction manifest bindings drifted")

            decision_state = contract["decision_rule_attestation_state"]
            for filename, expected_hash in contract[
                "decision_rule_artifact_sha256"
            ].items():
                decision_artifact = directory / filename
                if not decision_artifact.is_file():
                    fail(f"FZ-11 decision-rule artifact is missing: {filename}")
                if sha256_file(decision_artifact) != expected_hash:
                    fail(f"FZ-11 decision-rule artifact hash mismatch: {filename}")
                verify_ots_binding(
                    Path(str(decision_artifact) + ".ots"),
                    expected_hash,
                    decision_state,
                )
            verify_ots_binding(
                Path(str(decision_manifest_path) + ".ots"),
                contract["decision_rule_manifest_sha256"],
                decision_state,
            )

            decision_path = directory / "fz11_decision_rule_v2_2026-07-31.json"
            try:
                decision_rule = json.loads(decision_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError as error:
                fail(f"FZ-11 corrected decision rule is invalid JSON: {error}")
            expected_fail_rules = [
                "isolated intrinsic C4 is positive",
                "an isolated intrinsic anisotropic coefficient at j=1,2,3,4,5 is nonzero",
                (
                    "for resolved negative C4 with adequate sixth-order sensitivity, "
                    "linked B0, B6, or the rigid rotated I6 vector is excluded after "
                    "the fixed orientation profile"
                ),
                (
                    "the calibrated joint likelihood excludes the complete C4<0, "
                    "B0/C4^2=10/21, B6/C4^2=32/315, rotated-I6 branch manifold"
                ),
            ]
            if set(decision_rule) != {
                "schema",
                "created_utc",
                "original_target",
                "original_target_sha256",
                "prediction_receipt_sha256",
                "prediction_bytes_changed",
                "comparison_data_read",
                "operator_scope_clarification",
                "admissibility",
                "fail_at_five_sigma",
                "support",
                "inconclusive",
                "scope_of_failure",
                "excluded_exposure",
            }:
                fail("FZ-11 corrected decision-rule keys drifted")
            if (
                decision_rule.get("schema") != "oph.fz11.decision_rule_erratum.v1"
                or decision_rule.get("created_utc") != contract["decision_rule_utc"]
                or decision_rule.get("original_target") != contract["target_file"]
                or decision_rule.get("original_target_sha256")
                != contract["artifact_sha256"][contract["target_file"]]
                or decision_rule.get("prediction_receipt_sha256")
                != contract["in_repo_artifact_sha256"][FZ11_RECEIPT_REL]
                or decision_rule.get("prediction_bytes_changed") is not False
                or decision_rule.get("comparison_data_read") is not False
                or "sole hop support"
                not in str(decision_rule.get("operator_scope_clarification", ""))
                or "no independent isotropic fourth- or sixth-order term"
                not in str(decision_rule.get("operator_scope_clarification", ""))
                or "joint likelihood or full covariance"
                not in str(decision_rule.get("admissibility", ""))
                or "equal transverse-polarization action"
                not in str(decision_rule.get("admissibility", ""))
                or decision_rule.get("fail_at_five_sigma") != expected_fail_rules
                or "complete linked branch manifold within two sigma"
                not in str(decision_rule.get("support", ""))
                or "null or underpowered result"
                not in str(decision_rule.get("inconclusive", ""))
                or decision_rule.get("scope_of_failure")
                != (
                    "primitive twelve-port physical propagation branch; OPH-wide "
                    "only after issue #655 proves the branch forced and exclusive"
                )
                or "WMAP" not in str(decision_rule.get("excluded_exposure", ""))
            ):
                fail("FZ-11 corrected decision-rule bindings drifted")
            fz11_custody_git_history = verify_fz11_custody_history(
                contract, custody_root, directory
            )

        if contract_id == "FZ-12":
            prediction_path = directory / contract["prediction_file"]
            try:
                prediction = json.loads(prediction_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError as error:
                fail(f"FZ-12 frozen prediction snapshot is invalid JSON: {error}")
            if set(prediction) != {
                "schema",
                "frozen_utc",
                "source_commit",
                "source_receipt",
                "source_receipt_sha256",
                "status",
                "operator",
                "coefficients",
                "scale_free_relations",
                "harmonic_nulls",
                "shape",
                "comparison_state",
                "new_target_or_comparison_data_read",
                "positive_scale_lower_bound_required_for_null",
                "excluded_prior_data",
                "fz11_superseded",
                "scope_of_failure",
            }:
                fail("FZ-12 frozen prediction snapshot keys drifted")
            if (
                prediction.get("schema") != FZ12_FREEZE_SCHEMA
                or prediction.get("frozen_utc") != contract["frozen_utc"]
                or prediction.get("source_commit") != contract["source_commit"]
                or prediction.get("source_receipt") != FZ12_RECEIPT_REL
                or prediction.get("source_receipt_sha256")
                != contract["in_repo_artifact_sha256"][FZ12_RECEIPT_REL]
                or prediction.get("status") != FZ12_FREEZE_STATUS
                or prediction.get("operator")
                != "Lambda_a(k,n) = (1/(5 a^2)) sum_{j=1}^{30} [1 - cos(a k w_j.n)]"
                or prediction.get("coefficients")
                != {
                    "C4_over_a2": "-1/20",
                    "B0_over_a4": "1/840",
                    "B6_over_a4": "-1/12600",
                }
                or prediction.get("scale_free_relations")
                != {
                    "B0_over_C4_squared": "10/21",
                    "B6_over_C4_squared": "-2/63",
                    "B6_over_B0": "-1/15",
                }
                or prediction.get("harmonic_nulls") != [1, 2, 3, 4, 5]
                or "SO(3)/A5" not in str(prediction.get("shape", ""))
                or prediction.get("comparison_state")
                != "INELIGIBLE_UNARMED_PHYSICAL_PRODUCER_OPEN"
                or prediction.get("new_target_or_comparison_data_read") is not False
                or prediction.get("positive_scale_lower_bound_required_for_null")
                is not True
                or "WMAP" not in str(prediction.get("excluded_prior_data", ""))
                or "FZ-11" not in str(prediction.get("excluded_prior_data", ""))
                or prediction.get("fz11_superseded") is not False
                or prediction.get("scope_of_failure")
                != (
                    "seam-current edge physical propagation branch; OPH-wide "
                    "only after a forced and exclusive physical bridge theorem"
                )
            ):
                fail("FZ-12 frozen prediction snapshot content drifted")

            target_path = directory / contract["target_file"]
            try:
                target_text = target_path.read_text(encoding="utf-8")
            except UnicodeDecodeError as error:
                fail(f"FZ-12 frozen target is not UTF-8: {error}")
            required_target_fragments = (
                "# FZ-12 frozen target: source-seam edge propagation branch",
                f"Freeze time: {contract['frozen_utc']}",
                f"Source commit: `{contract['source_commit']}`",
                (
                    "Canonical receipt SHA-256: `"
                    f"{contract['in_repo_artifact_sha256'][FZ12_RECEIPT_REL]}`"
                ),
                "M6 = (30/7) r^6 - (2/7) I6",
                "Lambda_a(k,n) = (1/(5 a^2))",
                "B0/C4^2 = 10/21",
                "B6/C4^2 = -2/63",
                "B6/B0 = -1/15",
                "FZ12-R01 FAIL",
                "FZ12-R02 FAIL",
                "FZ12-R03 FAIL",
                "FZ12-R04 FAIL",
                "FZ12-R05 SUPPORT",
                "A null cannot",
                "source-derived lower bound on `a`",
                "FZ-11 remains the immutable primitive-vertex branch",
                "No new target or comparison data were read",
            )
            if any(
                fragment not in target_text for fragment in required_target_fragments
            ):
                fail("FZ-12 frozen target statement or decision rule drifted")

            decision_manifest_path = directory / contract["decision_rule_manifest"]
            if not decision_manifest_path.is_file():
                fail("FZ-12 decision-rule clarification manifest is missing")
            if (
                sha256_file(decision_manifest_path)
                != contract["decision_rule_manifest_sha256"]
            ):
                fail("FZ-12 decision-rule clarification manifest hash mismatch")
            try:
                decision_manifest = json.loads(
                    decision_manifest_path.read_text(encoding="utf-8")
                )
            except json.JSONDecodeError as error:
                fail(f"FZ-12 decision-rule clarification manifest is invalid: {error}")
            expected_original_bindings = {
                contract["target_file"]: contract["artifact_sha256"][
                    contract["target_file"]
                ],
                contract["prediction_file"]: contract["artifact_sha256"][
                    contract["prediction_file"]
                ],
                contract["registration_manifest"]: contract[
                    "registration_manifest_sha256"
                ],
                "source_prediction_receipt": contract["in_repo_artifact_sha256"][
                    FZ12_RECEIPT_REL
                ],
            }
            if (
                set(decision_manifest)
                != {
                    "schema",
                    "created_utc",
                    "original_prediction_frozen_utc",
                    "original_source_commit",
                    "original_custody_commit",
                    "prediction_bytes_changed",
                    "source_receipt_bytes_changed",
                    "comparison_data_read",
                    "reason",
                    "original_bindings",
                    "append_only_artifacts",
                }
                or decision_manifest.get("schema")
                != "oph.fz12.decision_rule_clarification_manifest.v1"
                or decision_manifest.get("created_utc")
                != contract["decision_rule_utc"]
                or decision_manifest.get("original_prediction_frozen_utc")
                != contract["frozen_utc"]
                or decision_manifest.get("original_source_commit")
                != contract["source_commit"]
                or decision_manifest.get("original_custody_commit")
                != contract["custody_commit"]
                or decision_manifest.get("prediction_bytes_changed") is not False
                or decision_manifest.get("source_receipt_bytes_changed") is not False
                or decision_manifest.get("comparison_data_read") is not False
                or "null is inconclusive" not in str(decision_manifest.get("reason", ""))
                or decision_manifest.get("original_bindings")
                != expected_original_bindings
                or decision_manifest.get("append_only_artifacts")
                != contract["decision_rule_artifact_sha256"]
            ):
                fail("FZ-12 decision-rule clarification manifest bindings drifted")

            decision_state = contract["decision_rule_attestation_state"]
            for filename, expected_hash in contract[
                "decision_rule_artifact_sha256"
            ].items():
                artifact_path = directory / filename
                if not artifact_path.is_file():
                    fail(f"FZ-12 decision-rule artifact is missing: {filename}")
                if sha256_file(artifact_path) != expected_hash:
                    fail(f"FZ-12 decision-rule artifact hash mismatch: {filename}")
                verify_ots_binding(
                    Path(str(artifact_path) + ".ots"), expected_hash, decision_state
                )
            verify_ots_binding(
                Path(str(decision_manifest_path) + ".ots"),
                contract["decision_rule_manifest_sha256"],
                decision_state,
            )

            decision_path = directory / FZ12_DECISION_RULE_JSON
            try:
                decision_rule = json.loads(decision_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError as error:
                fail(f"FZ-12 clarified decision rule is invalid JSON: {error}")
            expected_rule_keys = {
                "schema",
                "created_utc",
                "original_target",
                "original_target_sha256",
                "original_snapshot",
                "original_snapshot_sha256",
                "source_receipt",
                "source_receipt_sha256",
                "prediction_bytes_changed",
                "comparison_data_read",
                "prediction_manifold_changed",
                "admissibility",
                "fail_at_five_sigma",
                "support",
                "inconclusive",
                "null_rule",
                "scope_of_failure",
                "excluded_exposure",
            }
            failures = decision_rule.get("fail_at_five_sigma")
            if (
                set(decision_rule) != expected_rule_keys
                or decision_rule.get("schema")
                != "oph.fz12.decision_rule_clarification.v1"
                or decision_rule.get("created_utc") != contract["decision_rule_utc"]
                or decision_rule.get("original_target") != contract["target_file"]
                or decision_rule.get("original_target_sha256")
                != contract["artifact_sha256"][contract["target_file"]]
                or decision_rule.get("original_snapshot")
                != contract["prediction_file"]
                or decision_rule.get("original_snapshot_sha256")
                != contract["artifact_sha256"][contract["prediction_file"]]
                or decision_rule.get("source_receipt") != FZ12_RECEIPT_REL
                or decision_rule.get("source_receipt_sha256")
                != contract["in_repo_artifact_sha256"][FZ12_RECEIPT_REL]
                or decision_rule.get("prediction_bytes_changed") is not False
                or decision_rule.get("comparison_data_read") is not False
                or decision_rule.get("prediction_manifold_changed") is not False
                or not isinstance(failures, list)
                or len(failures) != 4
                or not str(failures[0]).startswith("FZ12-R01")
                or "no negative-C4 trigger" not in str(failures[0])
                or not str(failures[1]).startswith("FZ12-R02")
                or "no negative-C4 trigger" not in str(failures[1])
                or not str(failures[2]).startswith("FZ12-R03")
                or "resolved negative C4" not in str(failures[2])
                or not str(failures[3]).startswith("FZ12-R04")
                or "same-action, same-sector source theorem"
                not in str(decision_rule.get("null_rule", ""))
                or "including all a >= a_min"
                not in str(decision_rule.get("null_rule", ""))
                or "five standard deviations or more"
                not in str(decision_rule.get("null_rule", ""))
                or "Without that lower bound and power calculation"
                not in str(decision_rule.get("null_rule", ""))
                or "OPH-wide only after a separate theorem"
                not in str(decision_rule.get("scope_of_failure", ""))
                or "WMAP" not in str(decision_rule.get("excluded_exposure", ""))
            ):
                fail("FZ-12 clarified decision-rule bindings drifted")

            fz12_custody_git_history = verify_fz12_custody_history(
                contract, custody_root
            )

        results[contract_id] = {
            "verification": "verified",
            "attestation_state": state,
        }
    return {
        "state": "verified",
        "custody_root": str(custody_root),
        "fz11_custody_git_history": fz11_custody_git_history,
        "fz12_custody_git_history": fz12_custody_git_history,
        "contracts": results,
    }


def render(register: dict, rows: list[dict], bitcoin_upgrades: dict | None = None) -> str:
    if bitcoin_upgrades is None:
        bitcoin_upgrades = load_json(BITCOIN_UPGRADE_RECEIPT)
    upgraded_rows = set(bitcoin_upgrades["packages"])
    lines: list[str] = []
    lines.append("# The frozen-prediction ladder")
    lines.append("")
    lines.append(
        "Generated by `tools/build_fz_registry.py` from"
        " `claims/frozen_prediction_register.json`; edit the JSON, then"
        " regenerate. The standing register was established under issue #607."
    )
    lines.append("")
    lines.append(
        "The table below shows the registered targets and the current proof-file"
        " custody. Expand a record for its full content and decision rule."
    )
    lines.append("")
    lines.append("<details>")
    lines.append("<summary>Scope, exposure and historical classification</summary>")
    lines.append("")
    lines.append(register["policy"])
    lines.append("")
    lines.append(
        "For frozen rows, the Content, Kill band, and Comparison protocol text"
        " preserves the exact registration-time payload. Any present-tense issue"
        " reference inside that immutable text is historical provenance, not a"
        " live project-status assertion."
    )
    lines.append("")
    lines.append(
        "A freeze fixes a target against data released after its frozen time."
        " FZ-01's retained targets were frozen on 2026-07-17 at 07:18 UTC, so"
        " they are forward tests against later releases. Comparisons of the"
        " same formulas with data published before that time (Planck 2018,"
        " SPT-3G D1, ACT DR6 and DESI DR2 for the scalar tilt and the w-law)"
        " are postdictions and are recorded as such in the postdiction ledgers."
        " The WMAP ILC search dated 2026-07-17 tested an A5 angular template at"
        " l >= 6. FZ-11 and FZ-12 carry that template's coefficients, so their"
        " protocols exclude the search's data products; FZ-01's targets are not"
        " A5 angular templates and the search examined none of them."
    )
    lines.append("")
    lines.append(
        "Finite completion-lane theorem packages are not automatically"
        " prediction rungs. A1, A3, A4, B1, B2, B3, B4, B5, B6, B7, B8, B9, B11/B13, C1, C2, D1, and E1 emit no row here:"
        " each lacks"
        " a prospectively registered physical observable, attachment, and"
        " decision rule, so their exact finite results remain in the claim"
        " and postdiction ledgers only. A4 is a bounded conditional finite "
        "endpoint with source and cross-regulator realization open. B17's "
        "seven-clause operational observer receipt is fixed-regulator and has "
        "no source or laboratory attachment. B6 is an "
        "exact finite holonomy/character-phase packet with no physical "
        "connection or interference readout. C1 is an "
        "intrinsic Herm2 Lorentz module and Einstein-coordinate bridge with no "
        "physical soldering. C2 is a bounded algebraic event-frame soldering "
        "contract whose source atlas, physical cone, causal, and clock receipts "
        "remain open. D1 proves that record order supplies only order data, "
        "together with a non-affine regrading no-go and "
        "conditional affine clock/proper-time comparison only after supplied "
        "event and calibration data; it has no source physical clock, "
        "observable, or prospective decision rule. The observer-net work (historical lane #728, closed): "
        "its five-module continuation conditionally proves post-hoc source-operator "
        "generation, exact factors on a declared Cartesian carrier, a "
        "conditional-expectation diamond with coverage, and constant-tower "
        "transport. The source does not select that carrier, region map, or slot "
        "split, and the exact off-diagonal control shows the two slot expectations "
        "are not jointly injective; a source correlation/descent receipt and a "
        "nonconstant realization are not supplied. The mechanics work (historical lane #731, closed): supplied counting "
        "and trivial data conditionally realize only the uniform transition "
        "kernel, not those inputs, an initial law, or the complete reference. Its "
        "stationary-maximum control excludes universal Gibbs mode/minimizer "
        "recovery, not saddles or stationary phase, and its same-history theorem "
        "admits distinct regular real Legendre enrichments. None of these packets "
        "selects a physical continuation. B9's exact independent-sign pinching average "
        "and totalized-logarithm support countermodel fix algebraic interfaces, "
        "not a physical information observable. B11/B13's exact frame-rank, "
        "continuous-binary, finite-unsharp-web, real-source-tomography, and "
        "phase-free real-closure no-gos locate an algebraic commutator phase "
        "lift. The conjugation-orbit theorem and post-hoc repair-count diagnostic "
        "add no rung: the statistic and its arbitrary completion pairing were not "
        "preregistered and are ineligible as validation. The branch has no "
        "source-produced public quantum instrument or prospective observable; "
        "the complex effect, operational-additivity, and instrument continuation, "
        f"scoped under the historical lane #730, is owned by #{COMMON_MODEL_OWNER}. The B12 "
        "source-realization result likewise emits no row: in "
        "addition to its nonreversible H-theorem probe, an exact spectral and "
        "empirical-denominator obstruction rules out the hoped-for common "
        "reference through two exact direct mechanisms on the current artifact. "
        "The stochastic, nonlinear, dilated, and enriched replacement "
        "source/collar/tail routes, scoped under the historical premise program "
        f"#739, are owned by #{COMMON_MODEL_OWNER}; energy-clock calibration "
        "is empirical import PR-15. No observable "
        "or prospective decision rule is fixed. B14's "
        "oriented-face nearest-compact comparison is robust in the three added "
        "coordinate norms but remains metric- and repair-rule-conditional. E9's invariant-form "
        "dimension drop, constructed multifactor binding, and relative-coefficient "
        "identifiability are structural, source-unselected results and likewise "
        "do not create prediction rows. The fixed-cutoff unitary-power no-go "
        "rules out only direct convergence of a nonidentity exact update; its "
        "relative-evolution control leaves scattering open, and it supplies no "
        "prospective observable or prediction row. Each registered target remains "
        "bound to its own frozen content and decision rule."
    )
    lines.append("")
    lines.append(
        "Owner numbers and issue references inside immutable frozen content"
        " record the ownership in force when the row was registered and are"
        " explicitly historical; words such as `current` inside those immutable"
        " bytes also describe the registration-time state. A pending row's"
        " owning issue is open; for every historical lane the open issue that"
        " owns its remaining work is named in `tools/closed_lanes.py`, and"
        " current issue state lives in GitHub."
    )
    lines.append("")
    lines.append("</details>")
    lines.append("")
    lines.append("## Register at a glance")
    lines.append("")
    lines.append(
        "The short labels below are for navigation. Each record retains its"
        " complete registered content and decision rule below the table."
    )
    lines.append("")
    lines.append("| Freeze | Topic | Current proof custody | Frozen (UTC) |")
    lines.append("| --- | --- | --- | --- |")
    active_rows = [row for row in rows if row["status"] != "superseded_void"]
    superseded_rows = [row for row in rows if row["status"] == "superseded_void"]
    topic_labels = {
        "FZ-01": "Four retained cosmology targets",
        "FZ-02": "A5 angular invariants",
        "FZ-03": "Neutrino ordering and mixing",
        "FZ-05": "Direct record-count closure",
        "FZ-07": "Dark-sector RAR/BTFR scatter",
        "FZ-08": "DESI dark-energy evolution",
        "FZ-09": "Vacuum birefringence",
        "FZ-10": "Conditional tau-mass window",
        "FZ-11": "Primitive-port dispersion",
        "FZ-12": "Seam-edge dispersion",
        "FZ-13": "Fixed global capacity (conditional)",
        "FZ-14": "Integer-division Kerr comb",
        "FZ-15": "Edge-center scalar tilt (conditional)",
    }
    if {row["id"] for row in active_rows} != set(topic_labels):
        fail("ladder navigation labels must cover exactly the active rows")
    state_labels = {
        "frozen_attested": "Bitcoin proof stored",
        "frozen_stamped_upgrade_pending": "Calendar proof stored",
        "registered_pending_freeze": "Pending registration",
        "resource_deferred": "Deferred",
    }
    for row in active_rows:
        frozen = row["frozen_utc"] or "—"
        proof_state = (
            "Bitcoin proof added (explorer-assisted)"
            if row["id"] in upgraded_rows
            else state_labels[row["status"]]
        )
        lines.append(
            f"| [{row['id']}](#{row['id'].lower()}) | {topic_labels[row['id']]} |"
            f" {proof_state} | {frozen} |"
        )
    lines.append("")
    lines.append(
        "“Calendar proof stored” describes the committed `.ots` file without"
        " a Bitcoin block attestation. “Bitcoin proof added (explorer-assisted)”"
        " means an upgraded copy is committed with its original file binding"
        " and Bitcoin block header checked, and two public explorers reported"
        " the block on the best chain. Independent full-node verification has"
        " not been run for those upgrades."
    )
    lines.append("")
    lines.append("## Registered records")
    lines.append("")
    for row in active_rows:
        if row["owning_issue"] is None:
            owner = row["milestone"]
        elif row["status"] == "registered_pending_freeze":
            owner = (
                f"[#{row['owning_issue']}]"
                "(https://github.com/FloatingPragma/"
                "observer-patch-holography/issues/"
                f"{row['owning_issue']}); milestone {row['milestone']}"
            )
        else:
            owner = (
                f"historical [#{row['owning_issue']}]"
                "(https://github.com/FloatingPragma/"
                "observer-patch-holography/issues/"
                f"{row['owning_issue']}); historical milestone {row['milestone']}"
            )
        if row["frozen_utc"]:
            frozen = row["frozen_utc"]
        elif row["status"] == "resource_deferred":
            frozen = "not registered"
        elif row["status"] == "superseded_void":
            frozen = "not a valid freeze"
        else:
            frozen = "to freeze"
        lines.append(f"### {row['id']}")
        lines.append("")
        lines.append(f"- **Registered state:** `{row['status']}`")
        lines.append(f"- **Frozen (UTC):** {frozen}")
        lines.append(f"- **Registration-era owner:** {owner}")
        lines.append("")
        lines.append("<details>")
        lines.append("<summary>Full content and decision rule</summary>")
        lines.append("")
        lines.append("**Content**")
        lines.append("")
        lines.append(row["content"])
        lines.append("")
        lines.append("**Kill band**")
        lines.append("")
        lines.append(row["kill_band"])
        lines.append("")
        lines.append("</details>")
        lines.append("")
    lines.append("")
    lines.append("## Current structural landing and prediction boundary")
    lines.append("")
    lines.append(
        "The local face-curvature action and separate scalar Coulomb minimum"
        " strengthen the finite seam-carrier mathematics but do not arm or"
        " promote FZ-12: they supply no temporal propagation, physical"
        " position/scale/frame, source-to-laboratory field attachment,"
        " calibrated readout, or same-operator lower-bound comparison. The"
        " triple-observer quantum coupling and ambient-cofinal thermodynamic"
        " interface likewise emit no eligible empirical target. No new frozen"
        " row is minted by these structural theorems, and all immutable frozen"
        " bytes remain unchanged."
    )
    lines.append("")
    lines.append("## Executable model and exposed-data correspondence")
    lines.append("")
    lines.append(
        "The [native model audit](../model/PHYSICAL_AUDIT.md) records exact"
        " single-carrier signal/clock laws, repair-action decompositions and"
        " bounded calibration contracts. These are finite mathematical and"
        " software results under declared laws and preparations. They do not"
        " identify a physical primordial field, a physical Lorentzian action"
        " or an observed clock, and require no federation for their local scope."
        " Their validity is independent of whether an additional physical"
        " correspondence succeeds. The active modular-clock construction and"
        " its direct-gap comparison with published optical ratios are conditional"
        " model evidence and a retrospective attachment audit; neither creates"
        " a prospective freeze."
    )
    lines.append("")
    lines.append(
        "The [ACT DR6 ordered-chain reanalysis](../code/cosmology/tilt_posterior/README.md)"
        " evaluates the fixed edge-center tilt and preserves the historical"
        " clock alternative as a separate diagnostic. The P-ACT-LB data were"
        " released in 2025, before the FZ-01 and FZ-15 freezes. The result is"
        " retrospective tension under baseline Lambda-CDM; its posterior tail"
        " positions are neither p-values nor frozen verdicts. The ordinary"
        " convergence checkpoint and per-chain sensitivity do not certify"
        " the FZ-15 tail-ESS requirement. This reanalysis neither substitutes"
        " an eligible dataset nor changes any registered content, decision rule,"
        " custody record or status. The [public-data audit](research/PUBLIC_DATA_PREDICTION_AUDIT.md)"
        " states the additional physical hypotheses and comparison eligibility."
    )
    lines.append("")
    lines.append("## Superseded records outside the ladder")
    lines.append("")
    lines.append(
        "These identifiers preserve attested historical bytes and their current"
        " scientific disposition. They are not predictions, freezes, or scoring"
        " surfaces."
    )
    lines.append("")
    lines.append("| Record | Content | Status | Comparison authority |")
    lines.append("| --- | --- | --- | --- |")
    for row in superseded_rows:
        lines.append(
            f"| {row['id']} | {row['content']} | {row['status']} |"
            f" {row['comparison_protocol']} |"
        )
    lines.append("")
    lines.append("## Retrospective results outside the ladder")
    lines.append("")
    lines.append(
        "These records were evaluated after their comparison inputs were known."
        " They are not freezes, ladder rungs, predictions, or evidence from a"
        " prospective test. A former reservation remains visible only to make"
        " the bookkeeping transition traceable."
    )
    lines.append("")
    lines.append(
        "| Record | Former reservation | Result | Status | Source | Payload hash |"
    )
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for result in register["retrospective_results"]:
        owner = (
            f"[#{result['owning_issue']}](https://github.com/FloatingPragma/"
            "observer-patch-holography/issues/"
            f"{result['owning_issue']}) ({result['milestone']})"
        )
        lines.append(
            f"| {result['id']} | {result['former_ladder_reservation']} |"
            f" {result['content']} | {result['status']} | {owner} |"
            f" `{result['payload_sha256']}` |"
        )
    lines.append("")
    for result in register["retrospective_results"]:
        lines.append(
            f"- **{result['id']} protocol**: {result['comparison_protocol']}"
        )
        lines.append(
            f"  Evidential boundary: {result['evidential_boundary']}"
        )
        lines.append(f"  Payload: `{result['payload_path']}`.")
    lines.append("")
    lines.append("## Custody and attestation")
    lines.append("")
    lines.append(
        "The custody statements below preserve their scientific registration content;"
        " display metadata is normalized by the fixed October 11 receipt."
        " Later proof upgrades appear in the dated receipt below; they do not"
        " change a target, decision rule, freeze time, or scientific verdict."
    )
    lines.append("")
    for row in rows:
        if row["status"] in FROZEN_STATUSES:
            attestation = str(row["attestation"]).rstrip(".")
            lines.append(f"- **{row['id']}**: {row['custody']} {attestation}.")
            lines.append(
                f"  Content hash: `{row['content_sha256']}`."
                f" Comparison protocol: {row['comparison_protocol']}"
            )
    lines.append("")
    lines.append("## Later Bitcoin proof upgrades")
    lines.append("")
    lines.append(
        "The [2026-10-04 upgrade receipt]"
        "(../evidence/custody/bitcoin_upgrades_2026-10-04/upgrade_audit.json)"
        " records the original and upgraded proof hashes, file digests, earliest"
        " attesting blocks, block headers and explorer observations. Every"
        " archived `.ots` file remains unchanged. The official OpenTimestamps"
        " client parsed all upgraded paths; local checks recomputed their file"
        " bindings, block-header hashes, Merkle roots and proof of work."
        " Blockstream and mempool.space reported the same blocks on the best"
        " chain. This is explorer-assisted evidence, not independent Bitcoin"
        " full-node verification."
    )
    lines.append("")
    lines.append("| Freeze | Upgraded proofs | Manifest anchor (Bitcoin height) |")
    lines.append("| --- | ---: | ---: |")
    for row_id in sorted(upgraded_rows):
        summary = bitcoin_upgrades["packages"][row_id]
        height = summary["manifest_anchor_height"]
        block_hash = bitcoin_upgrades["blocks"][str(height)]["block_hash"]
        lines.append(
            f"| {row_id} | {summary['proof_count']} |"
            f" [{height}](https://blockstream.info/block/{block_hash}) |"
        )
    lines.append("")
    lines.append("## Custody verification contracts")
    lines.append("")
    lines.append("The states in this table are the original registration requirements.")
    lines.append("")
    lines.append("| Contract | Rows | Registration-time attestation | External custody |")
    lines.append("| --- | --- | --- | --- |")
    state_labels = {
        "bitcoin_attested": "Bitcoin block attestation present",
        "calendar_pending": "calendar commitments present; Bitcoin upgrade pending",
    }
    for contract_id, contract in register["external_custody_contracts"].items():
        lines.append(
            f"| {contract_id} | {', '.join(contract['rows'])} |"
            f" {state_labels[contract['attestation_state']]} |"
            f" `{contract['custody_path']}` |"
        )
    lines.append("")
    fz02_contract = register["external_custody_contracts"]["FZ-02"]
    lines.append(
        "FZ-02 is bound to original custody commit"
        f" `{fz02_contract['custody_commit']}` at"
        f" `{fz02_contract['custody_commit_utc']}` and source commit"
        f" `{fz02_contract['source_commit']}`. Its append-only custody erratum"
        " corrects"
        " the original timestamp, source-commit, and whole-file-versus-fenced-"
        "block hash metadata. Its append-only scientific erratum marks the"
        " unsupported level-six/level-three frame-lock clause ineligible"
        " pending issue #643. Neither erratum modifies any stamped artifact."
    )
    lines.append("")
    lines.append(
        "FZ-10 is bound to source commit"
        f" `{FZ10_SOURCE_COMMIT}`, the commit its vendored registration"
        " manifest names for the byte copy of the producing artifact, and to"
        f" original custody commit `{FZ10_CUSTODY_COMMIT}` at"
        f" `{FZ10_CUSTODY_COMMIT_UTC}`; the row's frozen time is the payload"
        " generation time recorded in the frozen target statement. The custody"
        " verifier checks the vendored manifest against that source commit."
    )
    lines.append("")
    lines.append(
        "FZ-01's kill band names per-target decision policies. The four policies"
        " are vendored with their registration manifest and OpenTimestamps proofs"
        " under `evidence/custody/falsification/frozen_targets/fz01_2026-07-17/`,"
        " and their numeric thresholds are these. Capacity band: N in"
        " [3.30836, 3.31137] x 10^122, equivalently Lambda l_P^2 in"
        " [2.84618, 2.84878] x 10^-122, scored in ln N against the band edges on"
        " the frozen Planck 2018 TT,TE,EE+lowE+lensing Omega_Lambda h^2 posterior;"
        " the band dies at z > 3 once the posterior width reaches"
        " sigma_lnN <= 0.01, and the uncorrected bridge value 3.532131543 x 10^122"
        " is carried as the null alternative. Scalar tilt: the two candidates"
        " n_s = 1 - P/48 = 0.9660214956 and n_s = 1 - e alpha sqrt(pi) ="
        " 0.96484114303 (menu size two) are scored per candidate with"
        " z = |n_s_post - n_s_candidate| / sigma_post, COMPATIBLE at z <= 2, FAIL"
        " at z > 3, INCONCLUSIVE between, in the frozen scorer order ACT DR6 with"
        " Planck, SPT-3G with Planck, CMB-S4; the lane fails when both candidates"
        " fail under one scorer. Acceleration scale: the presence reading"
        " 1.184737388 x 10^-10 m s^-2 (lambda_c = 1 - P/24) and the Poisson"
        " reading 1.179018696 x 10^-10 m s^-2 (lambda_c = exp(-P/24)) are scored"
        " per reading with the same z thresholds against a held-out sub-percent"
        " a_0 determination, SPARC and SPARC-derived recalibrations being"
        " ineligible; the branch fails when both readings fail under one scorer."
        " Dark-energy w-law: the point (w_0, w_a) = (-1, 0) on the frozen"
        " DESI DR3 + CMB + SNe posterior kills the fixed-N lane when it is"
        " excluded at more than 3 sigma with a thawing posterior mean, kills the"
        " fixed-capacity branch outright when evolving w is confirmed at more"
        " than 5 sigma, passes inside the 1 sigma region, and carries no verdict"
        " between 1 and 3 sigma. No FZ-01 verdict is registered at current"
        " precision."
    )
    lines.append("")
    lines.append(
        "Append-only note, 2026-10-19: DOE and NSF ended their support for the"
        " CMB-S4 project on 9 July 2025 ([CMB-S4](https://cmb-s4.org/)). FZ-01's"
        " frozen custody files, decision rule and scorer order are unchanged."
        " The terminal CMB-S4 pass has no decision data unless a CMB-S4 baseline"
        " product is released; substituting another dataset changes the"
        " combination, which the frozen policy allows only through a new"
        " registration before the data opening."
    )
    lines.append("")
    lines.append(
        "Append-only note, 2026-10-26: the five FZ-01 target files under"
        " `evidence/custody/falsification/frozen_targets/fz01_2026-07-17/` open"
        " with a drafting-stage header that calls them unregistered drafts."
        " Those headers are anchored bytes from the drafting stage; the"
        " registration manifest, this register and that directory's"
        " `PROVENANCE.md` supersede them. No frozen byte changes."
    )
    lines.append("")
    lines.append(
        "Append-only exposure note, 2026-10-26, FZ-11 and FZ-12: both rows"
        " share the quadratic coefficient delta_gamma2 = -a^2/20. With the"
        " natural scale a = sqrt(P) l_P it is -5.47 x 10^-58 eV^-2, which is"
        " 5.47 times the photon-only Pierre Auger 2022 bound"
        " delta_gamma2 > -1 x 10^-58 eV^-2; that bound requires a < 0.546 l_P."
        " The bound is scenario- and interaction-dependent and predates both"
        " freezes, so the exposure is retrospective and decides nothing under"
        " the frozen rules; the rows, decision rules and eligible data are"
        " unchanged. Sources:"
        " `code/a5_fingerprint/runtime/fz12_auger_threshold_diagnostic_receipt.json`,"
        " `docs/research/DISPERSION_FEASIBILITY.md` and"
        " `docs/research/PUBLIC_DATA_PREDICTION_AUDIT.md`."
    )
    lines.append("")
    lines.append(
        "Append-only disclosure note, 2026-10-26, FZ-13: the row was"
        " registered on 2026-10-04, after the DESI DR2 combinations with"
        " Pantheon+, Union3 and DESY5 supernovae were public; the cosmology"
        " ledger's seen-data diagnostics for (w_0, w_a) = (-1, 0) are 2.60,"
        " 3.38 and 3.93 standardized units for those three combinations. The"
        " registration names Pantheon+ and states no reason for that choice;"
        " its eligible decision data are the first DESI DR3 release named in"
        " its protocol."
    )
    lines.append("")
    lines.append(
        "Status tokens are registration-time states. `frozen_attested`: frozen"
        " with a Bitcoin attestation at registration. `frozen_stamped_upgrade_pending`:"
        " frozen with calendar commitments at registration; a later Bitcoin"
        " upgrade is recorded in the custody tree and shown in the current-custody"
        " column. `registered_pending_freeze`: registered, not yet frozen."
        " `resource_deferred`: deferred, with no open owner. `superseded_void`:"
        " retained only as void history."
    )
    lines.append("")
    lines.append(
        "FZ-10 retains its existing armed physical comparison. The October 4"
        " FZ-13 and FZ-15 registrations add specified conditional branch tests:"
        " their fixed-capacity and reserve-generator premises are hypotheses,"
        " not newly established OPH derivations. Eligibility is determined by"
        " each complete frozen protocol, and no comparison has been scored."
        " FZ-11 and FZ-12 carry exact five-sigma kill bands and an unarmed, ineligible"
        " physical comparison: they are frozen predictions whose arming requires"
        " the source-derived attachments their comparison protocols name."
    )
    lines.append("")
    lines.append(
        "The custody directories named in the frozen rows are under"
        " `evidence/custody/`, where the custody paths below resolve."
        " Two metadata-normalized copies are identified in"
        " `metadata_normalization_2026-10-11.json`; their scientific content is"
        " checked against the pinned historical public artifacts. Detached proofs"
        " for these copies bind the original bytes, not the normalized copies."
        " The validator checks the committed in-repo hash contracts, recomputes"
        " every manifest and artifact hash with that explicit distinction, checks each"
        " detached `.ots` digest, and distinguishes pending calendar commitments"
        " from Bitcoin block attestations. The original custody commits are"
        " recorded as provenance; the timestamps anchor the bytes independently"
        " of any Git history."
        " The local structural check does not contact a Bitcoin node; independent"
        " chain verification remains the job of `ots verify` after an upgrade."
    )
    lines.append("")
    lines.append(
        "Pending rows freeze at their milestones, before their comparison data"
    )
    lines.append(
        "is examined; validation requires each pending row to carry a milestone"
        " and an owning issue that is open; the open issues are listed in"
        " `tools/closed_lanes.py`. Pending ownership may be updated before a"
        " freeze."
    )
    lines.append(
        "Issue numbers in frozen rows are registration-time history; the open"
        " issue that owns each historical lane's remaining work is named in"
        " `tools/closed_lanes.py`."
    )
    lines.append(
        "Retrospective results are validated and rendered in their separate"
        " section. Their former reservations do not occur in the ladder table."
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--verify-fz11-lean",
        action="store_true",
        help=(
            "re-elaborate the FZ-11 Lean module and fail unless exactly five "
            "admission-free axiom reports are emitted"
        ),
    )
    args = parser.parse_args()
    register = load_json(REGISTER_PATH)
    rows = validate(register)
    source_history = verify_fz11_source_history(
        register["external_custody_contracts"]["FZ-11"]
    )
    fz12_source_history = verify_fz12_source_history(
        register["external_custody_contracts"]["FZ-12"]
    )
    lean_replay = verify_fz11_lean_replay() if args.verify_fz11_lean else None
    custody = verify_external_custody(register)
    bitcoin_upgrades = verify_bitcoin_upgrade_receipt()
    surface = render(register, rows, bitcoin_upgrades)
    if args.check:
        committed = (
            SURFACE_PATH.read_text(encoding="utf-8") if SURFACE_PATH.is_file() else ""
        )
        if committed != surface:
            print(
                "frozen-prediction register: docs/FROZEN_PREDICTION_LADDER.md is "
                "stale; run python tools/build_fz_registry.py",
                file=sys.stderr,
            )
            return 1
        print(
            "frozen-prediction register: external custody "
            f"{custody['state']} at {custody['custody_root']}"
        )
        print(
            "frozen-prediction register: FZ-11 source history "
            f"{source_history['state']}"
        )
        print(
            "frozen-prediction register: FZ-12 source history "
            f"{fz12_source_history['state']}"
        )
        if lean_replay is not None:
            print(
                "frozen-prediction register: FZ-11 Lean replay verified "
                f"({lean_replay['report_count']} axiom reports)"
            )
        print("frozen-prediction register: surface is current")
        return 0
    SURFACE_PATH.write_text(surface, encoding="utf-8", newline="\n")
    print(
        "frozen-prediction register: external custody "
        f"{custody['state']} at {custody['custody_root']}"
    )
    print(
        "frozen-prediction register: FZ-11 source history "
        f"{source_history['state']}"
    )
    print(
        "frozen-prediction register: FZ-12 source history "
        f"{fz12_source_history['state']}"
    )
    if lean_replay is not None:
        print(
            "frozen-prediction register: FZ-11 Lean replay verified "
            f"({lean_replay['report_count']} axiom reports)"
        )
    print(f"frozen-prediction register: wrote {SURFACE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
