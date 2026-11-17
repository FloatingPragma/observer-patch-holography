"""Fresh complete-source admission, public closure and frozen science comparison."""
import copy
import hashlib
import json
from pathlib import Path

from model.particle_common_contract.provenance import import_closure
from .model import calculate, source, PARENT, PARENT_SHA, FINITE_CLOCK_PROOF

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
REFERENCE = HERE / "research_reference.json"
REFERENCE_SHA = "1efaae14f6da07373d4e2cbeec8b1b4904cf745a5edb544515e4fc288237c972"
IGNORED_COMPARISON_FIELDS = {"schema", "input_sha256", "local_sha256"}


def validate_pins(pins):
    for name, expected in pins.items():
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("source custody requires repository-relative paths")
        path = ROOT / relative
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError("source custody mismatch: " + name)


def scientific_fields(value):
    return {key: item for key, item in value.items() if key not in IGNORED_COMPARISON_FIELDS}


def compare_research(value):
    if hashlib.sha256(REFERENCE.read_bytes()).hexdigest() != REFERENCE_SHA:
        raise ValueError("the immutable independent research receipt changed")
    if scientific_fields(value) != scientific_fields(json.loads(REFERENCE.read_text())):
        raise ValueError("public apparatus science differs from the frozen calculation")
    return {"reference_sha256": REFERENCE_SHA, "all_scientific_fields_equal": True,
            "ignored_only_provenance_fields": sorted(IGNORED_COMPARISON_FIELDS),
            "fresh_calculation_was_replayed": True}


def input_hashes():
    parent = source()
    inherited = dict(parent["input_sha256"])
    inherited[PARENT] = PARENT_SHA
    validate_pins(inherited)
    seeds = {str(path.relative_to(ROOT)) for path in HERE.iterdir() if path.suffix in (".py", ".md")}
    seeds.update(inherited)
    seeds.update({str(REFERENCE.relative_to(ROOT)), FINITE_CLOCK_PROOF, "model/requirements.txt"})
    names = import_closure(seeds)
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sorted(names)}


def certificate():
    pins = input_hashes()
    value = copy.deepcopy(calculate())
    comparison = compare_research(value)
    value["schema"] = "oph.affine_gauss_charged_clock.v1"
    value["input_sha256"] = pins
    value["independent_research_comparison"] = comparison
    return value


def write():
    value = certificate()
    (HERE / "results.json").write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    return value
