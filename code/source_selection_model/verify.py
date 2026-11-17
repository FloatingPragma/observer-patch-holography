"""Recompute every finite certificate and bind the analytic/proof source bytes."""
from pathlib import Path
import argparse
import hashlib
import json
import re

from . import channels, finite_model, geometry, grammar, verify_response
import pin_scope

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def pins():
    names = {str(p.relative_to(ROOT)).replace("\\", "/") for p in HERE.glob("*.py")}
    names.update({"code/source_selection_model/DERIVATION.md",
                  "code/source_selection_model/RECORD_GLUING.md",
                  "paper/tex_fragments/SOURCE_RECORD_GLUING.tex",
                  "code/source_selection_model/CONTRACT.md",
                  "code/source_selection_model/README.md",
                  "code/source_selection_model/response.json",
                  "code/source_selection_model/records.json",
                  "claims/axiom_registry.yaml", "docs/AXIOM_REFERENCE.md",
                  "code/a5_closure/port_current_inner_certificate.py",
                  "code/a5_closure/echosahedral_selector_certificate.py",
                  "code/a5_closure/response_grammar_completeness_certificate.py",
                  ".github/workflows/source-selection-model.yml",
                  ".github/workflows/lean-ci.yml", "Lean/Geometry.lean",
                  "Lean/lean-toolchain", "Lean/lake-manifest.json"})
    pending = ["Geometry.SourceSelectionLocalityAxiomAudit"]
    while pending:
        relative = "Lean/" + pending.pop().replace(".", "/") + ".lean"
        if relative in names or not (ROOT/relative).is_file():
            continue
        names.add(relative)
        for line in re.findall(r"^import (.+)$", (ROOT/relative).read_text(encoding="utf-8"), re.M):
            pending.extend(line.split())
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(pin_scope.bound(names))}


def evidence(response_path=None, record_path=None):
    support = geometry.tower()
    return {"schema": "oph.source-selection-model-finite-evidence.v1",
            "verification_scope": "finite identities and executions; full model proof is analytic",
            "response": verify_response.verify(verify_response.strict_load(response_path or HERE/"response.json")),
            "algebra": finite_model.algebra_controls(),
            "constraint_grammar": grammar.controls(),
            "central_channels": channels.controls(),
            "support": geometry.check_tower(support),
            "records": finite_model.record_controls(verify_response.strict_load(record_path or HERE/"records.json")),
            "refinement": finite_model.refinement_controls(support), "source_pins": pins()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-receipt", action="store_true")
    parser.add_argument("--receipt", type=Path, default=HERE/"receipt.json")
    parser.add_argument("--response", type=Path)
    parser.add_argument("--records", type=Path)
    args = parser.parse_args()
    actual = evidence(args.response, args.records)
    if args.write_receipt:
        args.receipt.write_text(json.dumps(actual, indent=2)+"\n", encoding="utf-8", newline="\n")
    else:
        retained = verify_response.strict_load(args.receipt)
        retained["source_pins"] = pin_scope.bound_pins(retained["source_pins"])
        # JSON comparison retains type distinctions (True is not an integer).
        if json.dumps(retained, sort_keys=True) != json.dumps(actual, sort_keys=True):
            raise ValueError("retained receipt or source custody disagrees with independent replay")
    print("SOURCE_SELECTION_MODEL_FINITE_EVIDENCE_VERIFIED")
    print(json.dumps({key: value for key, value in actual.items()
                      if key not in ("schema", "source_pins")}, sort_keys=True))


if __name__ == "__main__":
    main()
