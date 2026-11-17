"""Canonical serialization, fixed experiment specification and source custody."""
import hashlib
import re
from pathlib import Path

from source_temporal_acceptance.codec import canonical, digest, equal, load
import pin_scope

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SUPPORT = "code/source_routing/support_w12_l3.json"


def cases():
    result = []
    for name, ports, horizons, demands, weights in (
        ("chain", [0, 1, 2, 3], range(4, 9), [[1, 0], [0, 1]], [1, 1, 1]),
        ("captured", [0, 1, 14, 23, 45], range(6, 9), [[1, 0], [0, 1]], [1]*4),
        ("aggregate", [0, 1, 14, 23, 45], (4, 7), [[1, 1]], [1]*4),
        ("tilted", [0, 1, 14, 23, 45], (7, 8), [[1, 0], [0, 1]], [1, 2, 3, 5]),
    ):
        for horizon in horizons:
            result.append({"name": name, "ports": ports, "horizon": horizon,
                           "targets": demands, "weights": weights,
                           "preparation_writes": 4 if name == "chain" else 15360})
    return result


def pins():
    names = {SUPPORT, "docs/AXIOM_REFERENCE.md", "claims/axiom_registry.yaml",
             "Lean/lean-toolchain", "Lean/lakefile.lean", "Lean/lake-manifest.json",
             "code/source_checkpoint_selection/CONTRACT.md"}
    for package in ("source_checkpoint_selection", "source_temporal_acceptance"):
        names.update(str(p.relative_to(ROOT)).replace("\\", "/")
                     for p in (ROOT/"code"/package).glob("*.py"))
    pending = ["Geometry.SourceCheckpointAxiomAudit"]
    while pending:
        relative = "Lean/"+pending.pop().replace(".", "/")+".lean"
        path = ROOT/relative
        if relative in names or not path.is_file():
            continue
        names.add(relative)
        for line in re.findall(r"^import (.+)$", path.read_text(encoding="utf-8"), re.M):
            pending.extend(line.split())
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(pin_scope.bound(names))}


SCOPE = {
    "optimizer": "transition distribution on the complete finite move-word simplex",
    "cover": "one full-history classical algebra with weight one",
    "reference": "product of the specified positive move weights, normalized on all words",
    "constraint": "all specified public records must be recoverable at the supplied deadline",
    "sampling": "initial receiver sample and one sample after every native mean",
    "policy": "exact continuation-mass transition kernel; no failed-run rejection sampling",
    "execution": "every selected word executes its full horizon, including after first publication",
    "controls": "two independently variable real preparation coordinates; common baseline three",
    "host_work": "planning, retained samples, exact decoding and random choice are supplied host interfaces",
    "canonical_axiom_selection": False,
    "metric_radius_selection": False,
    "physical_clock_or_precision": False,
}
