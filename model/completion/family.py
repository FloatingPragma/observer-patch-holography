"""Exact parameter boundary for the first finite C/E completion campaign.

The tagged mixture is an instrument, not a coherent superposition of laws.
All probabilities are rational declarations. An event attempt is not duration.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path


CAMPAIGN_ID = "OPH-CE-FINITE-MEMORY-2026-10-28-v1"
ROOT = Path(__file__).resolve().parents[2]
CONTRACT = Path(__file__).with_name("campaign.json")
MAX_DENOMINATOR = 256
POSITIVE_FLOOR = Fraction(1, 16)


def _probability(value: Fraction, name: str, lower: Fraction = Fraction(0)) -> None:
    if type(value) is not Fraction:
        raise ValueError(f"{name} must be an exact Fraction, not a float or Boolean")
    if not lower <= value <= 1 or value.denominator > MAX_DENOMINATOR:
        raise ValueError(f"{name} outside the frozen rational campaign boundary")


@dataclass(frozen=True)
class FamilyParameters:
    c_weight: Fraction = Fraction(1, 2)
    repair_probability: Fraction = Fraction(1)
    exploration_probability: Fraction = Fraction(1, 2)
    memory_capacity: int = 1
    overflow: str = "export_oldest"

    def __post_init__(self) -> None:
        self.validate()

    def validate(self) -> None:
        _probability(self.c_weight, "c_weight")
        _probability(self.repair_probability, "repair_probability", POSITIVE_FLOOR)
        _probability(self.exploration_probability, "exploration_probability", POSITIVE_FLOOR)
        if type(self.memory_capacity) is not int or not 1 <= self.memory_capacity <= 8:
            raise ValueError("memory_capacity must be an integer in 1..8")
        if self.overflow not in ("export_oldest", "halt"):
            raise ValueError("overflow must export the old record or explicitly halt")

    def as_dict(self) -> dict:
        return {
            "c_weight": str(self.c_weight),
            "repair_probability": str(self.repair_probability),
            "exploration_probability": str(self.exploration_probability),
            "memory_capacity": self.memory_capacity,
            "overflow": self.overflow,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "FamilyParameters":
        if not isinstance(data, dict) or set(data) != {
            "c_weight", "repair_probability", "exploration_probability", "memory_capacity", "overflow"
        }:
            raise ValueError("complete exact parameter mapping required")
        probabilities = {}
        for key in ("c_weight", "repair_probability", "exploration_probability"):
            value = data[key]
            if type(value) is not str:
                raise ValueError("serialized probabilities must be rational strings")
            try:
                probabilities[key] = Fraction(value)
            except (ValueError, ZeroDivisionError) as error:
                raise ValueError("invalid serialized rational probability") from error
        return cls(**probabilities, memory_capacity=data["memory_capacity"], overflow=data["overflow"])


def load_contract(path: Path = CONTRACT, *, check_sources: bool = True) -> dict:
    contract = json.loads(path.read_text())
    if contract.get("schema") != "oph.finite_completion_campaign.v1" or contract.get("campaign_id") != CAMPAIGN_ID:
        raise ValueError("unknown finite campaign")
    if contract.get("scope") != "FINITE_COMPONENT_FAMILY" or contract.get("whole_program_frozen") is not False:
        raise ValueError("finite campaign must not promote the full program boundary")
    if contract.get("boundary_status") != "FROZEN_FINITE_CAMPAIGN":
        raise ValueError("the finite campaign boundary must be explicit")
    if contract.get("physical_comparison_count_added") != 0:
        raise ValueError("this campaign adds no physical comparison")
    expected = {
        "probability_max_denominator": MAX_DENOMINATOR,
        "repair_probability_min": str(POSITIVE_FLOOR),
        "exploration_probability_min": str(POSITIVE_FLOOR),
        "memory_capacity_range": [1, 8],
        "controller_state_limit": 256,
        "laws": ["C", "E"],
        "overflow_policies": ["export_oldest", "halt"],
    }
    if contract.get("parameter_boundary") != expected:
        raise ValueError("parameter boundary differs from executable family")
    members = contract.get("initial_members", [])
    initial = {
        "C": FamilyParameters(c_weight=Fraction(1)),
        "E": FamilyParameters(c_weight=Fraction(0)),
        "mixed_delayed": FamilyParameters(repair_probability=Fraction(3, 4)),
    }
    if len(members) != 3 or {row["name"] for row in members} != set(initial):
        raise ValueError("C and E endpoints plus mixed control are required")
    for member in members:
        if FamilyParameters.from_dict(member["parameters"]) != initial[member["name"]]:
            raise ValueError("a retained initial member differs from its frozen parameters")
    if check_sources:
        for group in ("source_bindings", "retained_controls"):
            rows = contract.get(group)
            if not rows:
                raise ValueError(f"missing {group}")
            for row in rows:
                source = (ROOT / row["path"]).resolve()
                if not source.is_relative_to(ROOT) or not source.is_file():
                    raise ValueError("missing or escaping source binding")
    return contract


def binding() -> dict:
    load_contract()
    return {"campaign_id": CAMPAIGN_ID,
            "path": "model/completion/campaign.json",
            "sha256": hashlib.sha256(CONTRACT.read_bytes()).hexdigest()}


if __name__ == "__main__":
    print(json.dumps(binding(), indent=2))
