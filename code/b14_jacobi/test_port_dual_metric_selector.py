"""Semantic mutation tests for the B14 port-dual metric-selector verifier.

The committed certificate must verify, and each semantic tamper class
must be rejected: a selection flip, a distance edit, a gap edit, a
sector-scale edit, a weight edit, a cone-flag flip, a controls edit, a
pinned-hash edit, a schema edit, a pinned-input-set shrinkage, a
smuggled status field, and a boundary or declared-rules weakening even
when the tamper recomputes the canonical self-hash.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CERT = HERE / "port_dual_metric_selector.certificate.json"
VERIFIER = HERE / "verify_port_dual_metric_selector.py"


def run_verifier(path: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(VERIFIER), str(path)],
        capture_output=True,
        text=True,
    )


def rehash(cert: dict) -> None:
    body = {k: v for k, v in cert.items() if k != "certificate_sha256"}
    payload = json.dumps(body, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True).encode("ascii")
    cert["certificate_sha256"] = "sha256:" + hashlib.sha256(payload).hexdigest()


def mutate(tmp_path: Path, transform, recompute_hash: bool = False) -> Path:
    cert = json.loads(CERT.read_text())
    transform(cert)
    if recompute_hash:
        rehash(cert)
    out = tmp_path / "mutated.json"
    out.write_text(json.dumps(cert, indent=1, sort_keys=True) + "\n")
    return out


def test_committed_certificate_verifies():
    result = run_verifier(CERT)
    assert result.returncode == 0, result.stderr


def test_selection_flip_rejected(tmp_path):
    def transform(cert):
        cert["selection"] = "F"

    assert run_verifier(mutate(tmp_path, transform)).returncode != 0


def test_distance_edit_rejected(tmp_path):
    def transform(cert):
        cert["squared_distances_at_measure_point"]["G"] = ["540", "0"]

    assert run_verifier(mutate(tmp_path, transform)).returncode != 0


def test_gap_edit_rejected(tmp_path):
    def transform(cert):
        cert["gaps"]["F_minus_G"] = ["0", "0"]

    assert run_verifier(mutate(tmp_path, transform)).returncode != 0


def test_sector_scale_edit_rejected(tmp_path):
    def transform(cert):
        cert["sector_scales"]["three_plus"] = ["1/6", "0"]

    assert run_verifier(mutate(tmp_path, transform)).returncode != 0


def test_weight_edit_rejected(tmp_path):
    def transform(cert):
        cert["port_dual_weight"] = "1/6"

    assert run_verifier(mutate(tmp_path, transform)).returncode != 0


def test_cone_flag_flip_rejected(tmp_path):
    def transform(cert):
        cert["sector_balanced"] = False

    assert run_verifier(mutate(tmp_path, transform)).returncode != 0


def test_controls_edit_rejected(tmp_path):
    def transform(cert):
        cert["controls"]["skewed_face_weighting_rejected"] = False

    assert run_verifier(mutate(tmp_path, transform)).returncode != 0


def test_pinned_hash_edit_rejected(tmp_path):
    def transform(cert):
        cert["pinned_inputs"]["invariant_metric_phase.certificate.json"] = "0" * 64

    assert run_verifier(mutate(tmp_path, transform)).returncode != 0


def test_schema_edit_rejected(tmp_path):
    def transform(cert):
        cert["schema"] = "oph.b14_port_dual_metric_selector.v2"

    assert run_verifier(mutate(tmp_path, transform, True)).returncode != 0


def test_pinned_input_shrinkage_rejected(tmp_path):
    def transform(cert):
        del cert["pinned_inputs"]["invariant_metric_phase.certificate.json"]

    assert run_verifier(mutate(tmp_path, transform, True)).returncode != 0


def test_boundary_weakening_rejected_despite_rehash(tmp_path):
    def transform(cert):
        cert["boundary"] = (
            "The source uniquely selects the compact family G; "
            "no conditionality applies."
        )

    assert run_verifier(mutate(tmp_path, transform, True)).returncode != 0


def test_declared_rules_weakening_rejected_despite_rehash(tmp_path):
    def transform(cert):
        cert["declared_rules"]["face_weighting"] = (
            "derived from the source axioms"
        )

    assert run_verifier(mutate(tmp_path, transform, True)).returncode != 0


def test_smuggled_field_rejected_despite_rehash(tmp_path):
    def transform(cert):
        cert["source_selection_status"] = "SOURCE-SELECTED UNCONDITIONALLY"

    assert run_verifier(mutate(tmp_path, transform, True)).returncode != 0


def test_stale_self_hash_rejected(tmp_path):
    def transform(cert):
        cert["scaling_receipt"] = "tampered"

    # No rehash: the canonical self-hash must catch the edit.
    assert run_verifier(mutate(tmp_path, transform)).returncode != 0
