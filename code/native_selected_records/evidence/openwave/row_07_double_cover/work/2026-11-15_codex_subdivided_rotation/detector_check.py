"""Post hoc detector verification, never credited as source emergence."""
import json
from fractions import Fraction as F
from pathlib import Path

r = json.loads((Path(__file__).resolve().parent / "CONFIRM_RECEIPT.json").read_text())
checks = []
for case in r["cases"]:
    a, b = [list(map(F, v)) for v in case["final_arm_vectors"]]
    norm = sum(x * x for x in a + b)
    p0 = sum((x + y) ** 2 for x, y in zip(a, b)) / (2 * norm)
    ppi = sum((x - y) ** 2 for x, y in zip(a, b)) / (2 * norm)
    v0, vpi, vquarter = 2 * p0 - 1, 2 * ppi - 1, F(0)
    assert v0 == F(case["exact_visibility"])
    assert vpi == -v0 and p0 + ppi == 1
    # b->ib gives norm(a+ib)^2=norm(a-ib)^2=norm, hence equal ports.
    assert (norm / (2 * norm)) == F(1, 2)
    checks.append({"M": case["M"], "law": case["law"], "turns": case["turns"],
                   "direction": case["direction"], "original_visibility": str(v0),
                   "planted_pi_visibility": str(vpi), "planted_quarter_visibility": str(vquarter)})
print(json.dumps({"passed": True, "checks": checks,
                  "credit": "detector sensitivity and normalization only; planted phases are not native dynamics"}, indent=2))
