"""Independent matrix-exponential controls; never imports the producer.

Geometry is reconstructed with exact Q(sqrt(5)) arithmetic. Numeric evolution
is computed from the six-by-six generator and an augmented twelve-by-twelve
matrix exponential, rather than the producer's trigonometric formulas.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "runtime/seam_maxwell_continuum_receipt.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load_receipt(path: Path) -> dict:
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, f"duplicate JSON key: {key}")
            result[key] = value
        return result
    def reject_constant(value):
        raise ValueError(f"non-finite JSON value: {value}")
    return json.loads(path.read_text(), object_pairs_hook=unique,
                      parse_constant=reject_constant)


def exact_geometry() -> np.ndarray:
    phi = (1 + sp.sqrt(5)) / 2
    vertices = [sp.Matrix(v) for v in
                [(0, -1, -phi), (-1, -phi, 0), (-phi, 0, -1), (1, -phi, 0),
                 (0, 1, -phi), (-phi, 0, 1), (phi, 0, -1), (0, -1, phi),
                 (-1, phi, 0), (phi, 0, 1), (1, phi, 0), (0, 1, phi)]]
    directions = []
    for x, y in itertools.combinations(vertices, 2):
        d = (y - x) / 2
        if sp.simplify(d.dot(d) - 1) == 0:
            directions.append(d)
    require(len(directions) == 30, "exact seam census")
    v = sp.Matrix(sp.symbols("x y z"))
    r2 = v.dot(v)
    for order, target in [(2, 10 * r2), (4, 6 * r2 ** 2)]:
        residual = sp.Poly(sp.expand(sum(d.dot(v) ** order for d in directions) - target), *v)
        require(all(sp.simplify(c) == 0 for c in residual.coeffs()), "exact moment identity")
    return np.array([list(d) for d in directions], dtype=float)


def generator(k: np.ndarray, omega: float) -> tuple[np.ndarray, np.ndarray]:
    r = np.linalg.norm(k)
    if r == 0:
        return np.zeros((6, 6), complex), np.zeros((6, 6))
    # Construct cross-product columns via its action on coordinate vectors.
    curl = 1j * np.column_stack([np.cross(k, e) for e in np.eye(3)]) * omega / r
    g = np.block([[np.zeros((3, 3)), curl], [-curl, np.zeros((3, 3))]])
    p = np.kron(np.eye(2), np.eye(3) - np.outer(k, k) / (r * r))
    return g, p


def verify(receipt: dict, *, geometry: np.ndarray | None = None) -> int:
    require(isinstance(receipt, dict) and set(receipt) == {
        "schema", "status", "source", "implementation", "universal_bounds",
        "proof_scope", "assumptions", "nonclaims", "momenta", "scales", "times", "cases"}, "root schema")
    require(receipt.get("schema") == "oph.seam_maxwell_continuum.v1", "schema")
    require(receipt.get("status") ==
            "ANALYTIC_CONTINUUM_ON_SUPPLIED_DOMAIN__PHYSICAL_ATTACHMENT_OPEN", "status")
    source = receipt.get("source", {})
    require(source == {
        "path": "code/a5_fingerprint/runtime/seam_current_edge_prediction_receipt.json",
        "sha256": "0b8f0f7573f556ef0f47158fe07eca002c5b35790231d1ef2b75518057d12915"}, "source contract")
    require(hashlib.sha256((ROOT / source["path"]).read_bytes()).hexdigest() == source["sha256"],
            "immutable source")
    paths = {"code/electromagnetism/seam_maxwell_continuum.py",
             "code/electromagnetism/verify_seam_maxwell_continuum.py",
             "Lean/Screen/SeamMaxwellContinuum.lean"}
    require(set(receipt.get("implementation", {})) == paths, "implementation census")
    for path, expected in receipt["implementation"].items():
        require(hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == expected, "implementation digest")
    require(receipt.get("universal_bounds") == {
        "symbol": "0 <= |k|^2 - Lambda_a(k) <= a^2 |k|^4 / 20",
        "frequency": "0 <= |k| - sqrt(Lambda_a(k)) <= a^2 |k|^3 / 20",
        "propagator": "||U_a(t,k)-U_0(t,k)|| <= min(2, |t| a^2 |k|^3 / 20)",
        "constant_forcing": "||B_a(t,k)-B_0(t,k)|| <= t^2 a^2 |k|^3 / 40"}, "bound contract")
    require(receipt.get("proof_scope") == {
        "lean": "global symbol, frequency, cosine and sine bounds",
        "paper": "unitary real L2 assembly, strong limit, H3 and Duhamel bounds",
        "numeric": "finite off-axis, zero-mode, ultraviolet and signed-time controls only"}, "proof scope")
    require(receipt.get("assumptions") == [
        "supplied common Euclidean R3 domain and Lebesgue measure",
        "complete equal-weight source seam symbol at each positive scale",
        "declared reversible transverse curl-pair evolution and common time",
        "same initial field and same supplied forcing across scales"], "assumptions")
    require(receipt.get("nonclaims") == [
        "source-selected spacetime or operational clock",
        "finite-scale local curl or exact finite propagation cone",
        "source-produced charged current or Gauss-sector attachment",
        "laboratory identification, units, measurement or prediction verdict"], "nonclaims")
    ks = [(0, 0, 0), (1, 0, 0), (1, 2, -2), (2, -1, 3), (0.125, 0.25, 0.5), (12, -8, 4)]
    scales, times = [1.0, 0.5, 0.25, 0.125], [-2.0, 0.0, 1 / 3, 1.0, 3.0]
    require(receipt.get("momenta") == [list(k) for k in ks], "momentum census")
    require(receipt.get("scales") == scales and receipt.get("times") == times, "scale/time census")
    expected_cases = list(itertools.product(range(len(ks)), scales, times))
    rows = receipt.get("cases", [])
    require(len(rows) == len(expected_cases), "case census")
    directions = exact_geometry() if geometry is None else geometry
    for row, (ki, a, t) in zip(rows, expected_cases, strict=True):
        require(set(row) == {"momentum_index", "scale", "time", "symbol", "propagator_error", "forcing_error"}, "row schema")
        require(type(row["momentum_index"]) is int and
                (row["momentum_index"], row["scale"], row["time"]) == (ki, a, t), "case identity")
        k = np.array(ks[ki], float)
        r = float(np.linalg.norm(k))
        lam = float(np.sum(1 - np.cos(a * directions @ k)) / (5 * a * a))
        tol = 2e-9 * max(1.0, r * r)
        require(np.isfinite(lam) and lam >= -tol, "symbol positivity")
        require(-tol <= r * r - lam <= a * a * r ** 4 / 20 + tol, "symbol bound")
        omega = np.sqrt(max(0.0, lam))
        require(-tol <= r - omega <= a * a * r ** 3 / 20 + tol, "frequency bound")
        ga, p = generator(k, omega)
        g0, _ = generator(k, r)
        ua, u0 = expm(t * ga), expm(t * g0)
        block_a = np.block([[ga, np.eye(6)], [np.zeros((6, 12))]])
        block_0 = np.block([[g0, np.eye(6)], [np.zeros((6, 12))]])
        ba, b0 = expm(t * block_a)[:6, 6:], expm(t * block_0)[:6, 6:]
        err_u, err_b = float(np.linalg.norm(ua - u0, 2)), float(np.linalg.norm(ba - b0, 2))
        for name, value in [("symbol", lam), ("propagator_error", err_u), ("forcing_error", err_b)]:
            require(isinstance(row[name], (int, float)) and not isinstance(row[name], bool)
                    and np.isfinite(row[name]) and abs(row[name] - value) < tol, f"independent {name}")
        require(err_u <= min(2.0, abs(t) * a * a * r ** 3 / 20) + tol, "propagator bound")
        require(err_b <= t * t * a * a * r ** 3 / 40 + tol, "forcing bound")
        require(np.linalg.norm(ua.conj().T @ ua - np.eye(6)) < tol, "unitarity")
        require(np.linalg.norm(ua @ p - p @ ua) < tol, "constraint preservation")
        gn, _ = generator(-k, omega)
        require(np.linalg.norm(expm(t * gn) - ua.conj()) < tol, "opposite-momentum reality")
    return len(rows)


if __name__ == "__main__":
    print(f"independent seam-Maxwell continuum controls: {verify(load_receipt(OUTPUT))} cases")
