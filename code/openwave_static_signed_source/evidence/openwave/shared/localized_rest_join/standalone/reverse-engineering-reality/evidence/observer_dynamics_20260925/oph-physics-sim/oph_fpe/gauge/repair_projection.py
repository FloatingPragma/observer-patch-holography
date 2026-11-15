from __future__ import annotations

from typing import Any

import numpy as np

from oph_fpe.claims import RECOVERED_CORE, with_claim_metadata


def exact_repair_projection_receipt(
    repair_operator: np.ndarray,
    projection_operator: np.ndarray | None = None,
    *,
    tolerance: float = 1e-10,
) -> dict[str, Any]:
    repair = np.asarray(repair_operator, dtype=np.complex128)
    projection = repair if projection_operator is None else np.asarray(projection_operator, dtype=np.complex128)
    if repair.shape != projection.shape or repair.ndim != 2 or repair.shape[0] != repair.shape[1]:
        raise ValueError("repair and projection operators must be square matrices with matching shape")
    tol = float(tolerance)
    repair_idempotence = _relative_norm(repair @ repair - repair, repair)
    projection_idempotence = _relative_norm(projection @ projection - projection, projection)
    projection_self_adjoint = _relative_norm(projection - projection.conj().T, projection)
    equality = _relative_norm(repair - projection, projection)
    receipt = (
        repair_idempotence <= tol
        and projection_idempotence <= tol
        and projection_self_adjoint <= tol
        and equality <= tol
    )
    report = {
        "mode": "finite_exact_repair_equals_projection",
        "EXACT_REPAIR_EQUALS_PROJECTION_RECEIPT": bool(receipt),
        "receipt": bool(receipt),
        "operator_dimension": int(repair.shape[0]),
        "tolerance": tol,
        "repair_idempotence_error": repair_idempotence,
        "projection_idempotence_error": projection_idempotence,
        "projection_self_adjoint_error": projection_self_adjoint,
        "repair_projection_relative_error": equality,
        "claim_boundary": (
            "finite support-visible collar check that the declared exact repair map equals the "
            "orthogonal/conditional-expectation projection; not a 4D Yang-Mills mass-gap proof"
        ),
    }
    return with_claim_metadata(
        report,
        claim_level=RECOVERED_CORE,
        receipt="EXACT_REPAIR_EQUALS_PROJECTION_RECEIPT",
    )


def _relative_norm(delta: np.ndarray, reference: np.ndarray) -> float:
    return float(np.linalg.norm(delta, ord="fro") / (np.linalg.norm(reference, ord="fro") + 1e-12))
