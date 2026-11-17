"""Full frozen matrix enumeration, not Monte Carlo histories."""
import json
from fractions import Fraction as F
from model import CALLS, density
from process import committed_visibility, execute, preparation

cases = [execute(h=4, turns=t, law=law, direction=d, start=s)
         for s in range(3) for d in (1, -1) for t in (1, 2)
         for law in ("native", "swap")]
cases += [execute(h=6, turns=t, direction=d) for d in (1, -1) for t in (1, 2)]
inputs = {"basis0": (1, 0, 0), "basis1": (0, 1, 0), "basis2": (0, 0, 1),
          "uniform": (1, 1, 1), "contrast": (1, -1, 0),
          "large_uniform": (10**18 + 1, 10**18 - 1, 10**18)}
preps = {law: {str(s): {name: preparation(law, density(v), s, True)[1]
                       for name, v in inputs.items()} for s in range(3)}
         for law in ("native", "swap")}
coherence = {str(t): committed_visibility(4, t) for t in (1, 2)}
failures = []
for c in cases:
    v, r, t = c["visibility"]["float"], c["ray_overlap"], c["turns"]
    good = ((v < -0.99 and r < -0.999) if t == 1 else (v > 0.97 and r > 0.999))
    if c["law"] == "swap":
        good = v > 0 and F(c["raw_visibility"]["exact"]) == 1
    if not good:
        failures.append({k: c[k] for k in ("h", "turns", "law", "direction", "start", "visibility")})
for law, starts in preps.items():
    for tests in starts.values():
        for result in tests.values():
            p = F(result["success_probability"])
            assert p >= F(1, 72) if law == "native" else p == F(1, 9)
print(json.dumps({"candidate": "A-native-incidence-rotation-v1", "cases": cases,
                  "preparation_challenges": preps, "committed_record_control": coherence,
                  "gate_failures": failures, "literal_endpoint_calls": CALLS,
                  "execution": "complete finite branches and exact rational products; no sampled histories"}, indent=2))
if failures:
    raise SystemExit(1)
