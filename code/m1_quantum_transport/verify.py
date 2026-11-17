"""Strict source custody and complete independently reconstructed quantum evidence."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "m1_quantum_transport"

from .check import need, same, verify_evidence
import pin_scope

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SOURCE_FILES = [f"code/m1_quantum_transport/{name}" for name in
                ("__init__.py", "model.py", "check.py", "tetra_model.py", "tetra_check.py",
                 "stencil_model.py", "stencil_check.py",
                 "verify.py", "build.py", "test_transport.py",
                 "README.md", "CONTRACT.md", "DERIVATION.md")]
SOURCE_FILES += ["Lean/Geometry/M1QuantumTransport.lean", "Lean/Geometry/M1QuantumTransportAxiomAudit.lean",
                 "paper/tex_fragments/M1_QUANTUM_TRANSPORT.tex",
                 "paper/tex_fragments/RECORD_GLUING_PRINCIPLE.tex",
                 "code/source_selection_model/RECORD_GLUING.md", "code/rg_principle/DERIVATION.md",
                 "Lean/Geometry/RecordGluingPrinciple.lean", ".github/workflows/lean-ci.yml",
                 ".github/workflows/m1-quantum-transport.yml", "requirements.txt", ".gitattributes"]
SOURCE_FILES = pin_scope.bound(SOURCE_FILES)
PARENTS = ("OPH-RG-LOCAL-CLOCK-REDUCTION", "OPH-RG-INFORMATION-RESOURCE-BOUND")


def pins(root=ROOT):
    return {name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in SOURCE_FILES}


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
        rows = [r for r in registry["claims"] if r["claim_id"] == name]
        need(len(rows) == 1, f"parent claim must occur exactly once: {name}")
        payload = json.dumps(rows[0], sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        result[name] = hashlib.sha256(payload.encode()).hexdigest()
    return result


def verify(path=HERE/"receipt.json"):
    packet = strict_load(path)
    need(type(packet) is dict and set(packet) == {"schema", "sources", "parent_claims", "evidence"},
         "receipt envelope is incomplete or has extra fields")
    same(packet["schema"], "oph-m1-quantum-transport-v1", "unknown receipt schema")
    same(pin_scope.bound_pins(packet["sources"]), pins(), "source pins differ")
    same(packet["parent_claims"], parent_claims(), "parent claims differ")
    verify_evidence(packet["evidence"])
    return packet


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", nargs="?", type=Path, default=HERE/"receipt.json")
    args = parser.parse_args()
    verify(args.path)
    print("Verified charged quantum flights, minimal tetrahedral transport, complete spectra, "
          "all Fock sectors, finite-stencil velocity bounds and exact-delay controls")
