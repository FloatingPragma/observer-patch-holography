"""Shared source-bound diagnostic contract; no physical gate promotion."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
REQUIRED_SOURCES = frozenset((
    "model/completion/campaign.json", "model/completion/observer.py",
    "model/completion/family.py",
    "model/native_modular_reference/history_transport.py",
    "model/native_modular_reference/history_clock.py",
    "model/native_modular_reference/local_preparation.py",
    "model/observer/adaptive_record_recovery.py", "model/observer/native.py",
    "model/observer/quantum_mixing.py", "model/observer/record_recovery.py",
    "model/observer_time_fields/spatial_dynamics.py",
    "model/relative_observer_geometry/intrinsic_geometry.py",
))

# These two canonical digests were checked against the byte-pinned inputs at
# 3d282b3d59763fa1a5171ad9a4fbc09fc0abf35a. They allow only renewal of the
# retained control's documentation custody, never a new scientific result.
_PREDECESSOR = "model/completion/campaign.json"
_PREDECESSOR_SHA256 = "6b12ab2491a801473e79917b0684fa6ead988056a1a8ee205e30c178271424ad"
_PREDECESSOR_CANONICAL = "866940997b82c46f0cea42de0bb653a0dcb4c07f66057868a5a50bfb0b4f8099"
_CONTROL = "model/common_source/native_defect_decision_results.json"
_CONTROL_SHA256 = "6f33c6225924af280c8f92da3df6e2f004637e0896c967ff124030fc3359c1f5"
_CONTROL_CANONICAL = "0dc04c091af1dd535a8b8f8fab5b8364ba944fc46ba7275a1b84bc267d7678d7"
_DOCUMENT = "docs/research/OPH_FIELD_RECONSTRUCTION.md"
_DOCUMENT_SHA256 = "c5514b54a45fd2411ac49c6bfc332e1148481d66d84d5195aa2fedbf0263ca5a"


def source_bindings(paths):
    result = {}
    for name in paths:
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            raise ValueError("existing repository source required")
        result[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return dict(sorted(result.items()))


def _canonical_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def _metadata_only_predecessor_compatible(expected, current):
    """Validate one documented custody renewal without rewriting frozen pins.

    The predecessor may change only its first retained-control digest. That
    control may change only its existing documentation digest. Both current
    references must resolve, the six executable dependencies stay byte-exact,
    and all scientific values, types, source paths and scope flags stay exact.
    Canonical comparison permits JSON whitespace/key-order differences only.
    """
    if (expected.keys() != current.keys()
            or expected.get(_PREDECESSOR) != _PREDECESSOR_SHA256
            or {p for p in expected if expected[p] != current[p]} != {_PREDECESSOR}):
        return False
    try:
        predecessor = strict_json.load(ROOT / _PREDECESSOR)
        reference = predecessor["retained_controls"][0]
        if (reference["path"] != _CONTROL
                or reference["sha256"] != source_bindings((_CONTROL,))[_CONTROL]):
            return False
        reference["sha256"] = _CONTROL_SHA256
        if _canonical_digest(predecessor) != _PREDECESSOR_CANONICAL:
            return False
        control = strict_json.load(ROOT / _CONTROL)
        sources = control["source_sha256"]
        if (type(sources) is not dict or _DOCUMENT not in sources
                or sources != source_bindings(sources)):
            return False
        sources[_DOCUMENT] = _DOCUMENT_SHA256
        return _canonical_digest(control) == _CONTROL_CANONICAL
    except (KeyError, IndexError, TypeError, ValueError, OSError):
        return False


def load_contract(path=None):
    value = strict_json.load(HERE / "campaign.json" if path is None else path)
    if value["schema"] != "oph.clock_defect_energy.diagnostic_contract.v1":
        raise ValueError("unknown diagnostic contract")
    for name in ("integrated_gate_accepted", "whole_program_frozen", "native_family_changed"):
        if value[name] is not False:
            raise ValueError("diagnostic contract cannot promote a physical or source claim")
    if (type(value["physical_comparison_count_added"]) is not int
            or value["physical_comparison_count_added"] != 0):
        raise ValueError("this diagnostic has no empirical comparison")
    if (value["scope"] != "FINITE_COMPONENT_DIAGNOSTICS"
            or value["boundary_status"] != "FROZEN_DIAGNOSTIC_CONTRACT"
            or value["predecessor"]["merged_into_main_claimed"] is not False):
        raise ValueError("diagnostic scope or pending dependency changed")
    if type(value["source_sha256"]) is not dict or set(value["source_sha256"]) != REQUIRED_SOURCES:
        raise ValueError("complete predecessor source bindings required")
    current = source_bindings(value["source_sha256"])
    if (current != value["source_sha256"]
            and not _metadata_only_predecessor_compatible(value["source_sha256"], current)):
        raise ValueError("predecessor source changed; review a new diagnostic contract")
    return value


def check_receipt(path, expected):
    """Reject duplicate keys and type-changing edits (including false for zero)."""
    actual = strict_json.load(path)
    canonical = lambda value: json.dumps(value, sort_keys=True, allow_nan=False)
    if canonical(actual) != canonical(expected):
        raise ValueError("diagnostic receipt differs from exact replay")


def binding():
    value = load_contract()
    return {"campaign_id": value["campaign_id"],
            "contract_sha256": hashlib.sha256((HERE / "campaign.json").read_bytes()).hexdigest(),
            "predecessor_git_commit": value["predecessor"]["git_commit"]}
