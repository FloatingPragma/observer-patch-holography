"""Live reference execution followed by independent semantic verification."""
from pathlib import Path
import argparse
import hashlib
import json
import re

from . import check_process, process
from source_selection_model.verify_response import strict_load
import pin_scope

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def pins():
    paths = {str(path.relative_to(ROOT)).replace("\\", "/") for path in HERE.glob("*.py")}
    paths.update({"code/source_record_gluing/capture.json",
                  "code/source_record_gluing/README.md", "Lean/Geometry.lean",
                  "paper/tex_fragments/SOURCE_RECORD_GLUING.tex",
                  "code/source_selection_model/RECORD_GLUING.md",
                  "code/source_selection_model/CONTRACT.md",
                  "code/source_selection_model/verify_response.py",
                  "code/a5_closure/port_current_inner_certificate.py",
                  "code/a5_closure/echosahedral_selector_certificate.py",
                  ".github/workflows/source-selection-model.yml", ".github/workflows/lean-ci.yml",
                  "Lean/lean-toolchain", "Lean/lake-manifest.json"})
    pending = ["Geometry.SourceRecordGluingAxiomAudit", "Geometry.GoldenSourceVolumeLimit",
               "Geometry.SourcePublicationAxiomAudit"]
    while pending:
        path = "Lean/" + pending.pop().replace(".", "/") + ".lean"
        if path in paths or not (ROOT/path).is_file():
            continue
        paths.add(path)
        for line in re.findall(r"^import (.+)$", (ROOT/path).read_text(encoding="utf-8"), re.M):
            pending.extend(line.split())
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(pin_scope.bound(paths))}


def verify_capture(packet):
    require = check_process.require
    require(type(packet) is dict and set(packet) == {"schema", "scope", "cases"}, "capture keys")
    require(packet["schema"] == "oph.record-gluing-reference.v2" and packet["scope"] ==
            "generated process theory of proposed RG; no native implementation asserted", "capture scope")
    require(type(packet["cases"]) is list, "case list")
    expected = {(q, -1, None) for q in (2, 3, 5)} | {
        (q, stage, i) for q in (2, 3, 5) for stage in (-1, 0) for i in range(q**3)}
    cases = {}
    for row in packet["cases"]:
        require(type(row) is dict and set(row) == {"q", "layers", "intervention", "intervention_stage", "stream_sha256",
                                                  "operation_counts", "checkpoint_sha256"}, "case keys")
        require(type(row["q"]) is int and type(row["layers"]) is int and row["layers"] == 2,
                "case parameters")
        require(row["intervention"] is None or type(row["intervention"]) is int, "case intervention")
        require(type(row["intervention_stage"]) is int, "case intervention stage")
        key = row["q"], row["intervention_stage"], row["intervention"]
        require(key in expected and key not in cases, "case coverage/uniqueness")
        cases[key] = row
    require(set(cases) == expected, "complete intervention denominator")
    summaries = []
    total_events = 0
    for q in (2, 3, 5):
        reference = None
        # Independent path multiplicities are built from the actually consumed
        # writers, not from the emitter's geometric menu or declared parents.
        coefficients = {}
        probes = [(-1, None)]+[(stage, i) for stage in (-1, 0) for i in range(q**3)]
        for stage, intervention in probes:
            result = check_process.check(q, 2, intervention,
                                         process.execute(q, 2, intervention, stage), stage)
            row = cases[q, stage, intervention]
            require(row["stream_sha256"] == result["sha256"], "retained event commitment")
            require(json.dumps(row["operation_counts"], sort_keys=True) ==
                    json.dumps(result["operation_counts"], sort_keys=True), "retained resource census")
            require(row["checkpoint_sha256"] == hashlib.sha256(json.dumps(
                result["layers"], separators=(",", ":")).encode()).hexdigest(), "retained checkpoints")
            total_events += result["events"]
            if intervention is None:
                reference = result
                for i in range(q**3):
                    coefficients[-1, i] = {(-1, i): 1}
                for (layer, target), parents in sorted(result["parents"].items()):
                    counts = {}
                    for parent in parents:
                        for source, value in coefficients[parent].items():
                            counts[source] = counts.get(source, 0)+value
                    # An independent perturbation at this commit is a new
                    # basis direction, even when its unperturbed value is
                    # functionally identical to another stored version.
                    counts[layer, target] = 1
                    coefficients[layer, target] = counts
            else:
                require(result["parents"] == reference["parents"], "intervention changed routing")
                for layer in range(2):
                    for target in range(q**3):
                        require(result["layers"][layer][target]-reference["layers"][layer][target] ==
                                coefficients[layer, target].get((stage, intervention), 0),
                                "complete consumed-writer intervention propagation")
        summaries.append({"q": q, "population": q**3, "attempts": 2*q**3+1,
                          "initial_interventions": q**3, "commit_interventions": q**3,
                          "events_per_attempt": reference["events"],
                          "reads_per_layer": reference["reads_per_layer"],
                          "minimum_squared_clock_margin": reference["minimum_squared_clock_margin"]})
    return {"schema": "oph.record-gluing-verification.v2",
            "verification_scope": "reference primitive execution and finite identities; RG native realization is proposed",
            "families": summaries, "total_reference_events": total_events, "source_pins": pins()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--capture", type=Path, default=HERE/"capture.json")
    parser.add_argument("--receipt", type=Path, default=HERE/"receipt.json")
    parser.add_argument("--write-receipt", action="store_true")
    args = parser.parse_args()
    result = verify_capture(strict_load(args.capture))
    if args.write_receipt:
        args.receipt.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8", newline="\n")
    else:
        retained = strict_load(args.receipt)
        retained["source_pins"] = pin_scope.bound_pins(retained["source_pins"])
        if json.dumps(retained, sort_keys=True) != json.dumps(result, sort_keys=True):
            raise ValueError("reference receipt or source custody disagrees with replay")
    print(json.dumps({k: v for k, v in result.items() if k != "source_pins"}, sort_keys=True))


if __name__ == "__main__":
    main()
