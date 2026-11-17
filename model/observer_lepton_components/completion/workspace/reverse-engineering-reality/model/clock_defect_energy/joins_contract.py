"""Successor contract for explicit finite theorem applications and services."""
from __future__ import annotations

import hashlib

from tools import strict_json
from .contract import HERE, source_bindings, check_receipt, load_contract as load_diagnostic

REQUIRED_SOURCES = frozenset((
    "model/clock_defect_energy/campaign.json",
    "model/native_modular_reference/history_transport.py",
    "model/native_modular_reference/history_clock.py",
    "model/native_modular_reference/common_caps.py",
    "model/native_modular_reference/full_flow.py",
    "model/observer/record_recovery.py", "model/observer/source_scalar/compiler.py",
    "model/observer_time_fields/spatial_dynamics.py",
    "model/observer_time_fields/continuum_dynamics.py",
    "Lean/Geometry/ProperTimeInternalAction.lean",
    "Lean/Geometry/InternalEnergyInertia.lean",
    "Lean/Geometry/MassShellKinematics.lean",
    "Lean/Geometry/SourceActionMeasure.lean",
    "Lean/Thermodynamics/GibbsReferenceEnergyIdentification.lean",
))


def load_contract():
    predecessor = load_diagnostic()
    value = strict_json.load(HERE/"joins_campaign.json")
    if (value["schema"] != "oph.clock_defect_energy.theorem_joins.v1"
            or value["scope"] != "FINITE_COMPONENT_DIAGNOSTICS"
            or value["predecessor_campaign"] != predecessor["campaign_id"]):
        raise ValueError("unknown finite theorem join contract")
    for name in ("integrated_gate_accepted", "whole_program_frozen", "native_family_changed"):
        if value[name] is not False:
            raise ValueError("theorem joins do not promote the native family or physical gates")
    if (type(value["physical_comparison_count_added"]) is not int
            or value["physical_comparison_count_added"] != 0
            or value["service_inventory_extended"] is not True):
        raise ValueError("wrong empirical or service scope")
    sources = value["source_sha256"]
    if (type(sources) is not dict or set(sources) != REQUIRED_SOURCES
            or source_bindings(sources) != sources):
        raise ValueError("finite theorem join source custody changed")
    return value


def binding():
    value = load_contract()
    return {"campaign_id": value["campaign_id"],
            "contract_sha256": hashlib.sha256((HERE/"joins_campaign.json").read_bytes()).hexdigest(),
            "predecessor_campaign": value["predecessor_campaign"],
            "pending_pr_head": value["pending_pr_head"]}
