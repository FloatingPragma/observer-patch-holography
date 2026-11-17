"""Independent paper-claim verification; not a candidate experiment or G6."""
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
IDS = ["OPH-NATIVE-" + suffix for suffix in (
    "DYNAMIC-UNIT-SELECTION", "LOCALIZED-LOWEST-MODE", "TRANSPORTED-CLOCK-LIMIT",
    "PLANAR-ROTATION-HOLONOMY", "SELECTED-STATE-JOIN")]
PATHS = ("code/native_selected_records/PROOF.md",
         "paper/tex_fragments/DYNAMIC_UNIT_SELECTION_AND_LOCAL_PHASES.tex",
         "claims/claim_registry.yaml", "claims/native_selected_records.json")
bindings = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in PATHS}
registry = json.loads((ROOT / "claims/claim_registry.yaml").read_text())
projection = {"claims": [c for c in registry["claims"] if c["claim_id"] in IDS],
              "claim_context": {k: registry["claim_context"][k] for k in IDS}}
assert len(projection["claims"]) == 5
mirror = json.loads((ROOT / "claims/native_selected_records.json").read_text())["claims"]
assert mirror == projection["claims"]
projection_sha = hashlib.sha256(json.dumps(
    projection, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()).hexdigest()


def repair(a, b, coin):
    low = (a + b) // 2
    high = a + b - low
    x, y = (low, high) if coin == 0 else (high, low)
    return x, y, a - x


probes = badprobe = 0
for a in range(-7, 8):
    for b in range(-7, 8):
        for coin in (0, 1):
            x, y, j = repair(a, b, coin)
            xx, yy, same = repair(x, y, coin)
            _, _, opposite = repair(xx, yy, 1 - coin)
            assert same == 0 and 2 * j + opposite == a - b
            badprobe += j + opposite != a - b
            probes += 1

variational = []
for n in (2, 3, 12, 24, 36):
    w, lam, epsilon = F(2 * n), F(2 * n + 2), F(1, 2 * (n + 2))
    threshold = lam + 2 * w * (2 * n + 3) * F(n - 1, n)
    mr = w * (1 - epsilon) * F(n + 2, n)
    zr = threshold + 1 - 2 * w * (1 / epsilon - 1) * F(n - 1, n)
    assert mr > lam and zr > lam
    variational.append({"N": n, "m_margin": str(mr - lam), "z_margin": str(zr - lam)})

constant = math.sqrt(3 / 343)
assert constant < .1 and 2 * math.asin(math.sqrt(.3)) < 1.2
assert math.sqrt((1 / 3) ** 2 * (.2 / .8) ** 3) <= 1 / 24 + 1e-15
edges = [(0, 1), (0, 1), (0, 2), (0, 2), (1, 2)]
R = np.eye(33)
for p in range(3):
    for e, (a, b) in enumerate(edges):
        for coin in (0, 1):
            q = [int(p == i) for i in range(3)]
            q[a], q[b], j = repair(q[a], q[b], coin)
            y, cc = q.index(1), coin ^ int(j != 0)
            i, o = 11 * p + 2 * e + coin, 11 * y + 2 * e + cc
            R[:, i] = 0
            R[o, i] = 1
assert np.array_equal(R @ R, np.eye(33))
omega = np.ones(3) / math.sqrt(3)
f = np.array([0., 1., -1.]) / math.sqrt(2)
g = np.array([-2., 1., 1.]) / math.sqrt(6)
clock, badclock = [], 0
for k in range(7):
    epsilon = 2. ** (-k)
    av = np.r_[np.full(10, epsilon / math.sqrt(10)), math.sqrt(1 - epsilon ** 2)]
    A = np.kron(np.eye(3), av[:, None])
    W = R @ (2 * A @ A.T - np.eye(33))
    columns, phases = [A @ omega], [0.]
    for vec, lam in ((f, 4), (g, 6)):
        p = 1 - epsilon ** 2 * lam / 10
        u = A @ vec
        w = (R @ u - p * u) / math.sqrt(1 - p * p)
        columns.append((u - 1j * w) / math.sqrt(2))
        phases.append(2 * math.asin(epsilon * math.sqrt(lam / 20)))
    J = np.column_stack(columns)
    residual = float(np.linalg.norm(W @ J - J @ np.diag(np.exp(1j * np.array(phases)))))
    assert residual < 1e-11 and np.linalg.norm(J.conj().T @ J - np.eye(3)) < 1e-11
    badclock += int(np.linalg.norm(W @ J - J @ np.diag(np.exp(2j * np.array(phases)))) > 1e-4)
    errors = []
    for t in (0., .13, 1., 3.7):
        actual = J.conj().T @ np.linalg.matrix_power(W, math.floor(t / epsilon)) @ J
        expected = np.diag(np.exp(1j * t * np.sqrt(np.array([0., 4., 6.]) / 5)))
        error = float(np.linalg.norm(actual - expected, 2))
        bound = t / 10 * epsilon ** 2 + 1.2 * epsilon
        assert error <= bound + 1e-10
        errors.append({"t": t, "error": error, "bound": bound})
    clock.append({"k": k, "chart_residual": residual, "floor_bounds": errors})

probability = F(1, 50) * F(3, 5) ** 16 + F(3, 200) * F(2, 5) ** 16
fidelity = 1 / (1 + F(3, 4) * F(2, 3) ** 16)
intervals = badinterval = 0
for omega in (.3, .89630114120727, 1.0988975913, 2.4):
    for eta in (.001, .01, .2):
        for q in (0., eta / 2, eta):
            for complement in (-1., 0., 1.):
                r = (1 - q) * math.cos(omega) + q * complement
                low = max(-1., (r - eta) / (1 - eta))
                high = min(1., (r + eta) / (1 - eta))
                assert low - 1e-14 <= math.cos(omega) <= high + 1e-14
                badinterval += abs(math.acos(max(-1, min(1, r))) - omega) > 1e-6
                intervals += 1

result = {
    "kind": "Independent mathematical audit, not a candidate experiment",
    "bindings": bindings, "canonical_projection_sha256": projection_sha,
    "canonical_projection": projection,
    "signed_probe_cases": probes, "wrong_probe_formula_rejected_cases": badprobe,
    "variational_coefficients": variational, "clock_derivative_constant": constant,
    "native_clock_checks": clock, "doubled_phase_mutant_rejected_refinements": badclock,
    "selection_probability_n8": float(probability), "selection_fidelity_n8": float(fidelity),
    "interval_cases": intervals, "uncertified_point_estimate_disagreements": badinterval,
    "scope": "Independent algebra/finite checks; no full instrument replay or PDF-layout audit.",
    "external_G6": False,
}
(Path(__file__).parent / "PAPER_REVIEW_CHECK_V2.json").write_text(
    json.dumps(result, indent=2, ensure_ascii=False) + "\n")
print(json.dumps({k: result[k] for k in ("signed_probe_cases", "wrong_probe_formula_rejected_cases",
    "doubled_phase_mutant_rejected_refinements", "interval_cases", "canonical_projection_sha256")}, indent=2))
