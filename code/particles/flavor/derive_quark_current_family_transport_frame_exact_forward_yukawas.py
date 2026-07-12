#!/usr/bin/env python3
"""Emit target-anchored quark mass textures on the declared transport frame.

The legacy artifact name is retained for compatibility. Its entries have units
of GeV and are assembled from a mixed-convention target packet, so they are not
physical dimensionless Yukawa matrices.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
STRENGTHENED_THEOREM_JSON = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_strengthened_physical_sigma_lift_theorem.json"
)
COMPLETION_JSON = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_exact_pdg_completion.json"
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "flavor" / "quark_current_family_transport_frame_exact_forward_yukawas.json"
)


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _decode_complex_matrix(payload: Any) -> np.ndarray:
    return np.asarray(payload["real"], dtype=float) + 1j * np.asarray(payload["imag"], dtype=float)


def _encode_complex_matrix(matrix: np.ndarray) -> dict[str, Any]:
    return {"real": np.real(matrix).tolist(), "imag": np.imag(matrix).tolist()}


def _load_forward_builder():
    module_path = ROOT / "particles" / "flavor" / "build_forward_yukawas.py"
    spec = importlib.util.spec_from_file_location("build_forward_yukawas", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"unable to load {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.build_artifact


def build_artifact(strengthened_theorem: dict[str, Any], completion: dict[str, Any]) -> dict[str, Any]:
    u_left = _decode_complex_matrix(strengthened_theorem["restricted_sigma_ud_phys_element"]["U_u_left"])
    d_left = _decode_complex_matrix(strengthened_theorem["restricted_sigma_ud_phys_element"]["U_d_left"])
    masses = completion["exact_running_values_gev"]
    masses_u = np.diag([float(masses["u"]), float(masses["c"]), float(masses["t"])])
    masses_d = np.diag([float(masses["d"]), float(masses["s"]), float(masses["b"])])
    y_u = u_left @ masses_u
    y_d = d_left @ masses_d

    payload = {
        "Y_u": _encode_complex_matrix(y_u),
        "Y_d": _encode_complex_matrix(y_d),
    }
    build_forward = _load_forward_builder()
    forward = build_forward(payload)
    forward["artifact"] = "oph_quark_current_family_transport_frame_exact_forward_yukawas"
    forward["generated_utc"] = _timestamp()
    forward["scope"] = "current_family_common_refinement_transport_frame_only"
    forward["proof_status"] = "closed_target_anchored_mixed_convention_mass_texture_audit"
    forward["public_promotion_allowed"] = False
    forward["matrix_classification"] = "mixed_scheme_GeV_mass_texture_matrices"
    forward["stored_entry_dimension"] = "GeV"
    forward["physical_yukawa_certified"] = False
    forward["single_common_scale_running_sextet"] = False
    forward["promotion_blockers"] = ["QUARK_COMMON_SCALE_DIMENSIONLESS_YUKAWA_CERTIFICATE"]
    forward["supporting_artifacts"] = {
        "strengthened_physical_sigma_lift_theorem": strengthened_theorem["artifact"],
        "exact_pdg_completion": completion["artifact"],
    }
    forward["metadata"] = {
        "note": (
            "Target-anchored mass textures on the declared current-family/common-refinement transport-frame "
            "carrier. Their singular values reproduce the chosen mixed-convention mass coordinates. They are "
            "not one common-scale running sextet and not physical dimensionless Yukawa matrices."
        ),
        "right_basis_gauge": "identity",
    }
    return forward


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the target-audit quark mass textures on the declared transport-frame carrier.")
    parser.add_argument("--strengthened-theorem", default=str(STRENGTHENED_THEOREM_JSON))
    parser.add_argument("--completion", default=str(COMPLETION_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_artifact(
        _load_json(Path(args.strengthened_theorem)),
        _load_json(Path(args.completion)),
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
