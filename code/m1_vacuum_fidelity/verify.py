"""Fail-closed receipt reader with independent full-catalog reconstruction."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "m1_vacuum_fidelity"

from .check import need, same, verify_evidence

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SOURCE_FILES = [f"code/m1_vacuum_fidelity/{p}" for p in
                ("__init__.py", "model.py", "check.py", "numerics.py", "verify.py", "build.py",
                 "test_fidelity.py", "test_scalar_observables.py",
                 "README.md", "CONTRACT.md", "DERIVATION.md")]
SOURCE_FILES += ["Lean/Geometry/M1VacuumFidelity.lean", "Lean/Geometry/M1VacuumFidelityAxiomAudit.lean",
                 "paper/tex_fragments/M1_VACUUM_FIDELITY.tex",
                 "paper/tex_fragments/SOURCE_SCALAR_QUANTUM.tex",
                 "paper/tex_fragments/RECORD_GLUING_PRINCIPLE.tex",
                 "Lean/Geometry/RecordGluingPrinciple.lean",
                 ".github/workflows/m1-vacuum-fidelity.yml", ".github/workflows/lean-ci.yml",
                 "requirements.txt", ".gitattributes"]
PARENTS = ("OPH-QFT-EFFECTIVE-QUANTUM-COMPARISON", "OPH-RG-INFORMATION-RESOURCE-BOUND",
           "OPH-RG-LOCAL-CLOCK-REDUCTION")


def pins(root=ROOT):
    return {p: hashlib.sha256((root/p).read_bytes()).hexdigest() for p in SOURCE_FILES}


def no_duplicates(pairs):
    result = {}
    for key, value in pairs:
        need(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_constant(value):
    raise ValueError(f"nonfinite JSON constant: {value}")


def strict_load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=no_duplicates,
                      parse_constant=reject_constant)


def parent_claims(path=ROOT/"claims/claim_registry.yaml"):
    registry = strict_load(path)
    result = {}
    for name in PARENTS:
        matches = [row for row in registry["claims"] if row["claim_id"] == name]
        need(len(matches) == 1, f"parent claim must occur exactly once: {name}")
        payload = json.dumps(matches[0], sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        result[name] = hashlib.sha256(payload.encode()).hexdigest()
    return result


def verify(path=HERE/"receipt.json"):
    packet = strict_load(path)
    need(type(packet) is dict and set(packet) == {"schema", "sources", "parent_claims", "evidence"},
         "receipt envelope is incomplete or has extra keys")
    same(packet["schema"], "oph-m1-vacuum-fidelity-v1", "unknown schema")
    same(packet["sources"], pins(), "source pins differ")
    same(packet["parent_claims"], parent_claims(), "parent claims differ")
    verify_evidence(packet["evidence"])
    return packet


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", nargs="?", type=Path, default=HERE/"receipt.json")
    args = parser.parse_args()
    verify(args.path)
    print("Verified 4 complete spectra, 512 observations, 24 full spatial covariance checks, "
          "702 forward and 702 inverse scalar writes, 6 complete signal rays, "
          "exact cubic algebra and resource powers")
