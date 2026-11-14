"""Reproduce bounded controls for the analytic seam-Maxwell continuum theorem.

This is a deterministic mathematical regression, not observed physical data
or a simulation of source-selected spacetime. The universal estimates are
proved in Lean and the L2 assembly is proved in the owning paper.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "runtime/seam_maxwell_continuum_receipt.json"
SOURCE = ROOT / "code/a5_fingerprint/runtime/seam_current_edge_prediction_receipt.json"
SOURCE_SHA = "0b8f0f7573f556ef0f47158fe07eca002c5b35790231d1ef2b75518057d12915"
LEAN = ROOT / "Lean/Screen/SeamMaxwellContinuum.lean"
SCHEMA = "oph.seam_maxwell_continuum.v1"
STATUS = "ANALYTIC_CONTINUUM_ON_SUPPLIED_DOMAIN__PHYSICAL_ATTACHMENT_OPEN"
MOMENTA = [(0, 0, 0), (1, 0, 0), (1, 2, -2), (2, -1, 3),
           (0.125, 0.25, 0.5), (12, -8, 4)]
SCALES = [1.0, 0.5, 0.25, 0.125]
TIMES = [-2.0, 0.0, 1 / 3, 1.0, 3.0]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value: dict) -> str:
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def seam_directions() -> np.ndarray:
    phi = (1 + np.sqrt(5.0)) / 2
    vertices = np.array([(0, s, t * phi) for s in (-1, 1) for t in (-1, 1)] +
                        [(t * phi, 0, s) for s in (-1, 1) for t in (-1, 1)] +
                        [(s, t * phi, 0) for s in (-1, 1) for t in (-1, 1)])
    directions = [(vertices[j] - vertices[i]) / 2
                  for i, j in itertools.combinations(range(12), 2)
                  if abs(np.linalg.norm(vertices[j] - vertices[i]) - 2) < 1e-12]
    assert len(directions) == 30
    return np.array(directions)


def symbol(a: float, k: np.ndarray) -> float:
    if not np.isfinite(a) or a <= 0 or not np.isfinite(k).all():
        raise ValueError("positive finite scale and finite momentum required")
    # Stable evaluation of the complete cosine symbol, including near k=0.
    return float(2 * np.sum(np.sin(a * (seam_directions() @ k) / 2) ** 2) / (5 * a * a))


def blocks(k: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    r = float(np.linalg.norm(k))
    if r == 0:
        return r, np.zeros((6, 6)), np.zeros((6, 6), dtype=complex)
    x, y, z = k / r
    cross = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    transverse = np.eye(3) - np.outer(k, k) / (r * r)
    p = np.kron(np.eye(2), transverse)
    j = np.block([[np.zeros((3, 3)), 1j * cross],
                  [-1j * cross, np.zeros((3, 3))]])
    return r, p, j


def propagator(omega: float, t: float, p: np.ndarray, j: np.ndarray) -> np.ndarray:
    return np.eye(6) - 2 * np.sin(omega * t / 2) ** 2 * p + np.sin(omega * t) * j


def constant_forcing(omega: float, t: float, p: np.ndarray, j: np.ndarray) -> np.ndarray:
    """Integral from 0 to t of the same propagator, including omega=0."""
    return (t * np.eye(6) + (t * np.sinc(omega * t / np.pi) - t) * p +
            omega * t * t / 2 * np.sinc(omega * t / (2 * np.pi)) ** 2 * j)


def build() -> dict:
    if digest(SOURCE) != SOURCE_SHA:
        raise ValueError("frozen source symbol drift")
    rows = []
    for ki, values in enumerate(MOMENTA):
        k = np.asarray(values, dtype=float)
        r, p, j = blocks(k)
        for a in SCALES:
            lam = symbol(a, k)
            omega = np.sqrt(lam)
            for t in TIMES:
                ua = propagator(omega, t, p, j)
                u0 = propagator(r, t, p, j)
                ba = constant_forcing(omega, t, p, j)
                b0 = constant_forcing(r, t, p, j)
                def rounded(x: float) -> float:
                    return float(f"{x:.10g}")
                rows.append({"momentum_index": ki, "scale": a, "time": t,
                             "symbol": rounded(lam),
                             "propagator_error": rounded(np.linalg.norm(ua - u0, 2)),
                             "forcing_error": rounded(np.linalg.norm(ba - b0, 2))})
    return {
        "schema": SCHEMA, "status": STATUS,
        "source": {"path": SOURCE.relative_to(ROOT).as_posix(), "sha256": SOURCE_SHA},
        "implementation": {p.relative_to(ROOT).as_posix(): digest(p) for p in
                           [Path(__file__).resolve(), HERE / "verify_seam_maxwell_continuum.py", LEAN]},
        "universal_bounds": {
            "symbol": "0 <= |k|^2 - Lambda_a(k) <= a^2 |k|^4 / 20",
            "frequency": "0 <= |k| - sqrt(Lambda_a(k)) <= a^2 |k|^3 / 20",
            "propagator": "||U_a(t,k)-U_0(t,k)|| <= min(2, |t| a^2 |k|^3 / 20)",
            "constant_forcing": "||B_a(t,k)-B_0(t,k)|| <= t^2 a^2 |k|^3 / 40",
        },
        "proof_scope": {
            "lean": "global symbol, frequency, cosine and sine bounds",
            "paper": "unitary real L2 assembly, strong limit, H3 and Duhamel bounds",
            "numeric": "finite off-axis, zero-mode, ultraviolet and signed-time controls only",
        },
        "assumptions": ["supplied common Euclidean R3 domain and Lebesgue measure",
                        "complete equal-weight source seam symbol at each positive scale",
                        "declared reversible transverse curl-pair evolution and common time",
                        "same initial field and same supplied forcing across scales"],
        "nonclaims": ["source-selected spacetime or operational clock",
                      "finite-scale local curl or exact finite propagation cone",
                      "source-produced charged current or Gauss-sector attachment",
                      "laboratory identification, units, measurement or prediction verdict"],
        "momenta": [list(k) for k in MOMENTA], "scales": SCALES, "times": TIMES,
        "cases": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = canonical(build())
    if args.check:
        if OUTPUT.read_text() != result:
            raise SystemExit("seam-Maxwell continuum receipt is stale")
        print("seam-Maxwell continuum receipt parity OK")
    else:
        OUTPUT.parent.mkdir(exist_ok=True)
        OUTPUT.write_text(result)
        print(f"wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
