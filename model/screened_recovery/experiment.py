"""Source-bound receipt for the screened q=1 active recovery service."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

from tools import strict_json
from ..particle_screening import screening
from . import protocol


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "screened_recovery_results.json"
PREDECESSOR_FILES = (
    "model/particle_screening/screening.py",
    "model/particle_screening/screening_results.json",
    "model/common_source/record_code.py",
    "model/common_source/native_routing.py",
    "model/observer/record_recovery.py",
    "model/observer/source_discrimination/circuit.py",
    "model/observer/source_discrimination/experiment.py",
    "model/observer/source_scalar/compiler.py",
    "model/exchange_sector/sector.py",
)
PACKET_FILES = (
    "model/screened_recovery/__init__.py",
    "model/screened_recovery/protocol.py",
    "model/screened_recovery/experiment.py",
    "model/screened_recovery/test_protocol.py",
    "model/screened_recovery/README.md",
    "model/screened_recovery/campaign.json",
)


def _sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def science_result():
    """Machine-stable exact outputs; numerical diagnostics live beside them."""
    return {
        "schema": "oph.screened_recovery.science.v1",
        "scope": "FINITE_COMPONENT_DIAGNOSTICS",
        "service": protocol.exact_service_certificate(),
        "compiled_layout": protocol.compiled_layout_certificate(),
        "resources": protocol.resource_report(),
        "cp03_accepted": False,
        "autonomous_particle_identified": False,
        "rest_energy_identified": False,
        "physical_lifetime_identified": False,
    }


def campaign_binding(science):
    campaign = strict_json.load(HERE/"campaign.json")
    expected = {
        "schema": "oph.screened_recovery.campaign.v1",
        "campaign_id": "screened-recovery-20261028-v1",
        "predecessor_git_commit": "7cef80f58f51f06e70ccda3d4373344ca3993dc6",
        "scope": "FINITE_COMPONENT_DIAGNOSTICS",
        "whole_program_frozen": False,
        "boundary_status": "FROZEN_DIAGNOSTIC_CONTRACT",
        "integrated_gate_accepted": False,
        "native_family_changed": False,
        "physical_comparison_count_added": 0,
    }
    for key, value in expected.items():
        if type(campaign.get(key)) is not type(value) or campaign[key] != value:
            raise ValueError("screened-recovery campaign scope changed")
    source = campaign.get("predecessor_source_sha256")
    if (type(source) is not dict or set(source) != set(PREDECESSOR_FILES)
            or any(type(digest) is not str or digest != _sha(path)
                   for path, digest in source.items())):
        raise ValueError("screened-recovery predecessor source changed")
    result_hash = hashlib.sha256(_canonical(science)).hexdigest()
    if campaign.get("science_result_sha256") != result_hash:
        raise ValueError("screened-recovery exact science result changed")
    return {"campaign_id": campaign["campaign_id"],
            "campaign_sha256": _sha("model/screened_recovery/campaign.json"),
            "science_result_sha256": result_hash}


def source_bindings():
    bindings = screening.source_bindings()
    for path in PREDECESSOR_FILES+PACKET_FILES:
        bindings[path] = _sha(path)
    bindings["model/particle_screening/screening_results.json"] = _sha(
        "model/particle_screening/screening_results.json")
    return dict(sorted(bindings.items()))


def execution_diagnostics():
    """Fresh backend observations; never historical claims in the receipt."""
    science = science_result()
    return {
        **science,
        "source_replay": protocol.source_replay_report(),
        "source_composition": protocol.source_composition_report(),
        **campaign_binding(science),
        "source_sha256": source_bindings(),
    }


def _admission_intervals():
    """Fixed numerical acceptance ceilings, not outward arithmetic proofs."""
    intervals = {}
    def error(path, upper='1/100000000'):
        intervals[path] = ('0', upper)
    for field in ('forward_full512_macro_operator_error_plus_polar_change',
                  'inverse_full512_macro_operator_error_plus_polar_change',
                  'source_X_conjugation_operator_error', 'fair_coin_preparation_operator_error',
                  'encoder_operator_error_bound', 'decoder_operator_error_bound'):
        error(('source_composition', field))
    for count in ('1', '2', '3', '4'):
        for field in ('requested_inverse_adjoint_error', 'requested_inverse_uniform_to_zero_error',
                      'requested_word_uniform_column_error', 'requested_word_unitarity_residual',
                      'uniform_target_unitary_completion_error_upper'):
            error(('source_composition', 'mixing_angle_synthesis_diagnostics', count, field))
    for law in ('C', 'E'):
        for field in ('maximum_one_corrected_attempt_log_error', 'maximum_one_corrected_attempt_operator_error',
                      'complete_local_native_dilation_error'):
            error(('source_composition', 'laws', law, field))
        for attempts in ('1', '2', '10', '100'):
            error(('source_composition', 'laws', law, 'encoded_t_attempt_operator_error_upper', attempts),
                  '1/1000000')
        error(('source_replay', 'all_positive_raw_branch_source_recoveries', law,
               'maximum_full512_source_recovery_error'))
    for field in ('encoder_operator_error', 'maximum_source_replay_error'):
        error(('source_replay', field))
    for step, probability in enumerate((F(1), F(1, 4))):
        prefix = ('source_replay', 'corrected_two_attempt_CE_continuation', step)
        for field in ('post_native_column_gram_error', 'corrected_against_exact_code_operator_error'):
            error((*prefix, field))
        intervals[(*prefix, 'raw_instrument_total_probability_for_two_input_columns')] = (
            str(probability-F(1, 10**8)), str(probability+F(1, 10**8)))
    return intervals


def _admit(diagnostics):
    intervals = _admission_intervals()
    threshold = ('source_replay', 'maximum_allowed_numerical_execution_error')
    observed = set()
    def visit(value, path=()):
        if path in intervals:
            lower, upper = intervals[path]
            if (type(value) is not float or not math.isfinite(value)
                    or not F(lower) <= F(value) <= F(upper)):
                raise ArithmeticError(f'fresh numerical diagnostic exceeds fixed admission at {path}')
            observed.add(path)
            return {'kind': 'fresh_numerical_admission', 'lower': lower, 'upper': upper}
        if path == threshold:
            if type(value) is not float or value != 1e-8:
                raise ArithmeticError('the existing source replay gate changed')
            observed.add(path)
            return '1/100000000'
        if type(value) is dict:
            return {key: visit(item, (*path, key)) for key, item in value.items()}
        if type(value) is list:
            return [visit(item, (*path, index)) for index, item in enumerate(value)]
        if type(value) is float:
            raise ArithmeticError(f'unregistered backend numerical diagnostic at {path}')
        return value
    result = visit(diagnostics)
    if observed != set(intervals) | {threshold}:
        raise ArithmeticError('a required fresh numerical admission check is missing')
    return result


def run():
    result = _admit(execution_diagnostics())
    result['receipt_schema'] = 'oph.screened_recovery.receipt.v2'
    result['numerical_admission_semantics'] = (
        'Fixed acceptance intervals checked against fresh backend diagnostics; '
        'not observed historical residuals, outward interval proofs, processor-failure '
        'density certificates, hardware error rates or classical-reset fault bounds. '
        'The composition observations concern successful processor blocks and all native raw controls.')
    return result


def _compare(expected, actual, path="root"):
    if type(expected) is not type(actual):
        raise ValueError(f"receipt type mismatch at {path}")
    if isinstance(expected, dict):
        if set(expected) != set(actual):
            raise ValueError(f"receipt key mismatch at {path}")
        for key in expected:
            _compare(expected[key], actual[key], f"{path}.{key}")
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise ValueError(f"receipt length mismatch at {path}")
        for i, (left, right) in enumerate(zip(expected, actual, strict=True)):
            _compare(left, right, f"{path}[{i}]")
    elif isinstance(expected, float):
        raise ValueError(f"raw backend floats do not belong in the canonical receipt at {path}")
    elif expected != actual:
        raise ValueError(f"receipt mismatch at {path}")


def verify(path=OUTPUT):
    _compare(strict_json.load(Path(path)), run())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    options = parser.add_mutually_exclusive_group()
    options.add_argument("--write", nargs="?", const=OUTPUT, type=Path)
    options.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    options.add_argument('--diagnostics', action='store_true')
    args = parser.parse_args()
    result = execution_diagnostics() if args.diagnostics else run()
    if args.write:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n",
                              encoding='utf-8', newline='\n')
    elif args.check:
        _compare(strict_json.load(args.check), result)
    else:
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
