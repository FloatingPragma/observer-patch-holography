"""Reproduce the conditional Standard Model derivation and its spectrum."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from .contract import theorem_contract
from .defect_constraints import constraint_certificate
from .defects import SOURCE_PATHS as DEFECT_SOURCES, defect_catalogue
from .matching import example_certificate
from .spectrum import build_catalogue


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE / "results.json"
PARENTS = (
    "model/common_source/native_matter_attachment_results.json",
    "model/common_source/calibrated_matter_results.json",
    "model/defect_field_benchmark/results.json",
    "model/topological_defect_spectrum/results.json",
    "model/bounded_defect_index/results.json",
    "model/spin_texture_energy/results.json",
)
DIRECT_SOURCES = (
    "Lean/Screen/SMStructureComposition.lean",
    "Lean/Screen/SMStructureAdequacySurface.lean",
    "Lean/Screen/ElectroweakBreakingComposition.lean",
    "Lean/Screen/WeylYukawaConventions.lean",
    "Lean/Screen/LocalGaugeJetAction.lean",
    "Lean/Screen/ExteriorSelection.lean",
    "model/matched_einstein_limit/PROOF.md",
    "tracking/premise_register.json",
    "tools/strict_json.py",
)


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"


def source_bindings():
    """Keep inherited custody while computing this package's finite result."""
    inherited = {}
    for name in PARENTS:
        parent = strict_json.load(ROOT / name)
        bindings = parent.get("source_bindings", parent.get("source_sha256"))
        if type(bindings) is not dict or not bindings:
            raise ValueError("parent calculation lacks source bindings")
        for path, digest in bindings.items():
            target = (ROOT / path).resolve()
            if not target.is_relative_to(ROOT.resolve()):
                raise ValueError("source binding escapes the repository")
            if hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                raise ValueError("stale inherited source: " + path)
            if path in inherited and inherited[path] != digest:
                raise ValueError("inconsistent inherited source: " + path)
            inherited[path] = digest
    paths = set(inherited) | set(PARENTS) | set(DIRECT_SOURCES) | set(DEFECT_SOURCES)
    contract = theorem_contract()
    for section in ("theorem_chain", "conditional_interfaces"):
        for row in contract[section]:
            paths.update(ref["path"] for ref in row["references"])
    paths.update(str(path.relative_to(ROOT)) for path in PACKAGE.iterdir()
                 if path.is_file() and path.suffix in {".py", ".md"})
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in sorted(paths)}


def report():
    bindings = source_bindings()
    result = {
        "schema": "oph.conditional_standard_model.v1",
        "all_pass": True,
        "pass_scope": "finite algebra, conditional theorem references and full report reproducibility",
        "theorem_contract": theorem_contract(),
        "catalogue": build_catalogue(),
        "spectral_transfer_control": example_certificate(),
        "additive_defect_charge_constraints": constraint_certificate(),
        "defect_catalogue": defect_catalogue(),
        "canonical_physical_premises_discharged": [],
        "integrated_gate_accepted": False,
        "independently_predicted_measured_mass_count": 0,
        "source_selection_derived": False,
        "source_bindings": bindings,
    }
    if source_bindings() != bindings:
        raise ValueError("source files changed during calculation")
    return result


def validate(value):
    if canonical(value) != canonical(report()):
        raise ValueError("conditional Standard Model receipt differs from full fresh calculation")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    operation = parser.add_mutually_exclusive_group(required=True)
    operation.add_argument("--write", nargs="?", const=str(OUTPUT))
    operation.add_argument("--check", nargs="?", const=str(OUTPUT))
    args = parser.parse_args()
    if args.write:
        path = Path(args.write)
        path.write_text(canonical(report()), encoding="utf-8")
        print("wrote " + str(path))
    else:
        validate(strict_json.load(Path(args.check)))
        print("conditional Standard Model derivation and spectrum replay passed")


if __name__ == "__main__":
    main()
