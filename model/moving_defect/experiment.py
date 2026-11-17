"""Reproduce the exact moving-code packet and freshly gate source execution."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from tools import strict_json
from . import transport
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "results.json"
SOURCES = (
    "model/moving_defect/__init__.py", "model/moving_defect/transport.py",
    "model/moving_defect/experiment.py", "model/moving_defect/test_transport.py",
    "model/moving_defect/README.md", "model/moving_defect/campaign.json",
    "model/particle_screening/screening.py", "model/particle_screening/README.md",
    "model/screened_recovery/protocol.py", "model/screened_recovery/README.md",
    "model/screened_transport/transport.py", "model/exchange_sector/sector.py",
    "model/common_source/record_code.py", "model/common_source/ENCODED_MEMORY.md",
    "model/common_source/native_routing.py", "model/common_source/slow_geometry.py",
    "model/observer/native.py", "model/observer/record_recovery.py",
    "model/observer/adaptive_record_recovery.py",
    "model/observer/source_discrimination/circuit.py",
    "model/observer/source_discrimination/experiment.py",
    "model/observer/source_scalar/compiler.py",
    "model/local_repair_family/local_bias.py", "model/local_repair_family/mixing.py",
    "model/local_repair_family/record_clock.py", "model/local_repair_family/record_clock_campaign.json",
    "model/local_repair_family/record_clock_results.json",
    "model/local_repair_family/inference.py",
)
CAP = 1e-8


def contract():
    value = strict_json.load(HERE / "campaign.json")
    if (value["baseline_head"] != "67aa7b6f4a904a1a90ea6d77bf2f66c741f2ba0c"
            or value["scope"] != "FINITE_ACTIVE_MOVING_MEMORY_AND_OBSTRUCTION"
            or value["fixed_decoration"] != list(transport.parameters().port_weights)
            or value["source_replay_fixed_error_cap"] != "1/100000000"
            or value["whole_program_frozen"] is not False
            or value["integrated_gate_accepted"] is not False
            or type(value["physical_comparison_count_added"]) is not int
            or value["physical_comparison_count_added"] != 0):
        raise ValueError("moving code campaign authority changed")
    return value


def admit(diagnostics):
    names = ("all_sixty_encoder_max_operator_error", "complete_path_max_branch_operator_error",
             "complete_path_max_total_mass_error")
    if set(diagnostics) != set(names)|{"changing_core_C_E_path"}:
        raise ArithmeticError("source diagnostics missing or extra")
    result = {}
    for name in names:
        value = diagnostics[name]
        if type(value) is not float or not math.isfinite(value) or not 0 <= value <= CAP:
            raise ArithmeticError("source execution failed the fixed numerical admission cap")
        result[name] = {"kind": "fresh_numerical_admission", "lower": "0", "upper": "1/100000000"}
    expected = []
    weight = transport.F(1)
    for core, law, seam, coin, raw in transport.witness_path():
        _, target, _ = transport.branch(core, law, seam, coin, raw)
        weight *= transport.branch_probability(core, law, seam, coin, raw)
        expected.append({"core_before": list(core), "law": law, "seam": list(seam),
                         "coin": coin, "raw": raw, "core_after": list(target),
                         "status": "continue" if target in transport.LIVE else "opposite_core_stop",
                         "retained_raw_slots": 2*(1 << transport.circuit.RECORD_BITS[law]),
                         "cumulative_flag_weight_given_selected_seams": str(weight)})
    compare(expected, diagnostics["changing_core_C_E_path"])
    result["changing_core_C_E_path"] = expected
    return result


def science():
    uniform = transport.stopping_certificate(transport.local.Parameters())
    if uniform["initial_expected_stopping_attempts"] != "150":
        raise ArithmeticError("independent fair-coin control changed")
    return {"service": transport.exact_service_certificate(),
            "stopping": transport.stopping_certificate(),
            "uniform_control": {"initial_expected_stopping_attempts": "150",
                                "edge_distance2_stop_lumped_rows": [["29/30", "1/30", "0"], ["1/30", "19/20", "1/60"]]},
            "exterior": transport.exterior_certificate(),
            "fixed_reference_modular_compatibility": transport.clock_compatibility(),
            "source_layout": transport.compiled_layout_certificate(),
            "resources": transport.resource_certificate()}


def run():
    value = contract()
    return {"schema": "oph.moving_defect.results.v1", "campaign": value,
            "science": science(), "source_replay": admit(transport.source_replay_diagnostics()),
            "source_sha256": bound_pins({name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES})}


def compare(expected, actual, path="root"):
    if type(expected) is not type(actual):
        raise ValueError(f"receipt type differs at {path}")
    if isinstance(expected, dict):
        if expected.keys() != actual.keys():
            raise ValueError(f"receipt keys differ at {path}")
        for key in expected:
            compare(expected[key], actual[key], f"{path}.{key}")
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise ValueError(f"receipt length differs at {path}")
        for index, (a,b) in enumerate(zip(expected, actual, strict=True)):
            compare(a, b, f"{path}[{index}]")
    elif expected != actual:
        raise ValueError(f"receipt value differs at {path}")


def verify(path=OUTPUT):
    compare(run(), strict_json.load(Path(path)))
    return {"verified": True, "scope": "finite_active_moving_memory", "owner_issue": 1139}


def main():
    parser = argparse.ArgumentParser()
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--write", action="store_true")
    choice.add_argument("--check", action="store_true")
    choice.add_argument("--diagnostics", action="store_true")
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(json.dumps(run(), sort_keys=True, indent=2)+"\n")
        print(json.dumps({"written": str(OUTPUT.relative_to(ROOT))}))
    elif args.check:
        print(json.dumps(verify()))
    else:
        print(json.dumps(transport.source_replay_diagnostics(), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
