"""Bind the conserved-sector action theorem to complete native C/E kernels.

The sparse kernel is assembled from every binary state, core, seam and coin.
The Python checks establish the finite model premises; Lean proves their
consequences for arbitrary positive sector factors. No physical energy is set.
"""
from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..clock_defect_energy import reference_balance
from ..exchange_sector import sector
from ..observer import record_recovery as records
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "sector_action_results.json"
BASE_SOURCES = (
    "model/observer/record_recovery.py",
    "model/observer/native.py",
    "model/observer/quantum_mixing.py",
    "model/common_source/slow_geometry.py",
    "model/common_source/native_join.py",
    "model/repair_scheduler/coins.py",
    "model/completion/observer.py",
    "model/completion/family.py",
    "model/exchange_sector/sector.py",
    "model/clock_defect_energy/continuation.py",
    "model/clock_defect_energy/defect_history.py",
    "model/clock_defect_energy/contract.py",
    "model/clock_defect_energy/joins_contract.py",
    "model/clock_defect_energy/reference_balance.py",
    "model/clock_defect_energy/reference_balance_results.json",
    "model/active_code_energy/energy_bridge_results.json",
    "Lean/Geometry/SourceActionMeasure.lean",
    "evidence/native_geometric_source_20260925/spec.json",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py",
)
PACKET_SOURCES = (
    "Lean/Geometry/ConservedSectorAction.lean",
    "model/conserved_sector_action/__init__.py",
    "model/conserved_sector_action/sector_action.py",
    "model/conserved_sector_action/test_sector_action.py",
    "model/conserved_sector_action/README.md",
    "model/conserved_sector_action/campaign.json",
)
# The campaign predates later documentation and receipt-source updates. Keep
# its exact predecessor receipt rather than rebinding that frozen input to
# a newer live receipt. All other campaign inputs remain live and exact.
PREDECESSOR_ENERGY_SOURCE = "model/active_code_energy/energy_bridge_results.json"
PREDECESSOR_ENERGY_CAPTURE = "model/conserved_sector_action/predecessor/energy_bridge_results.json"
PREDECESSOR_ENERGY_SHA256 = "d59512d3d37c2f6a5ef80e3b36aa20a90a3bb8d365049dc8f33e535d2320bbfd"


def digest(path: str) -> str:
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def predecessor_bindings() -> dict[str, str]:
    """Resolve only the campaign's named historical receipt to its capture."""
    if digest(PREDECESSOR_ENERGY_CAPTURE) != PREDECESSOR_ENERGY_SHA256:
        raise ValueError("sector action predecessor capture changed")
    return bound_pins({name: PREDECESSOR_ENERGY_SHA256 if name == PREDECESSOR_ENERGY_SOURCE
            else digest(name) for name in BASE_SOURCES})


def source_bindings() -> dict[str, str]:
    """Bind current code and actual input files, including the frozen capture."""
    current = bound_pins({name: digest(name) for name in (*BASE_SOURCES, *PACKET_SOURCES)
               if name != PREDECESSOR_ENERGY_SOURCE})
    current[PREDECESSOR_ENERGY_CAPTURE] = digest(PREDECESSOR_ENERGY_CAPTURE)
    return dict(sorted({**sector.source_bindings(), **current}.items()))


def native_kernel(law: str) -> tuple[list[dict[int, F]], tuple[int, ...]]:
    """Every raw basis outcome, including all zero-current holds."""
    if law not in ("C", "E"):
        raise ValueError("expected C or E")
    basis, cores = sector.BASIS, sector.CORES
    source_index = {x: i for i, x in enumerate(basis)}
    core_index = {c: i for i, c in enumerate(cores)}
    size = len(cores)
    rows, charges = [], []
    seams = records.seams()
    denominator = 2 * len(seams)
    for x in basis:
        outcomes = [(y, record) for seam in seams for coin in (0, 1)
                    for y, record in [records.recorded_step(law, x, seam, coin)]]
        for core in cores:
            counts = Counter(source_index[y] * size + core_index[
                sector.moved_core(core, law, record)] for y, record in outcomes)
            rows.append({j: F(n, denominator) for j, n in sorted(counts.items())})
            charges.append(sum(x[p] for p in core))
    return rows, tuple(charges)


def certify_kernel(rows: list[dict[int, F]], charges: tuple[int, ...],
                   reference: tuple[F, ...]) -> dict:
    """Check the exact support, normalization and sectorwise reciprocity premises."""
    size = len(charges)
    if not size or len(rows) != size or len(reference) != size:
        raise ValueError("state domain mismatch")
    if min(reference) <= 0 or sum(reference, F()) != 1:
        raise ValueError("reference is not faithful and normalized")
    incoming = [F() for _ in rows]
    movements = 0
    for i, row in enumerate(rows):
        if sum(row.values(), F()) != 1 or any(p <= 0 for p in row.values()):
            raise ValueError("kernel is not a positive sparse stochastic row")
        for j, p in row.items():
            if not 0 <= j < size:
                raise ValueError("transition escapes state domain")
            if charges[i] != charges[j]:
                raise ValueError("transition changes the declared sector")
            if reference[i] * p != reference[j] * rows[j].get(i, F()):
                raise ValueError("reference reciprocity fails")
            incoming[j] += reference[i] * p
            movements += i != j
    if tuple(incoming) != reference:
        raise ValueError("reference is not stationary")
    components = {}
    for q in sorted(set(charges)):
        domain = {i for i, a in enumerate(charges) if a == q}
        root = min(domain)
        seen, todo = {root}, deque([root])
        while todo:
            for j in rows[todo.popleft()]:
                if j not in seen:
                    seen.add(j)
                    todo.append(j)
        if seen != domain:
            raise ValueError("a sector is not connected")
        components[str(q)] = len(seen)
    return {
        "states": size,
        "directed_movement_entries": movements,
        "positive_entries_including_holds": sum(map(len, rows)),
        "minimum_hold_probability": str(min(row.get(i, F()) for i, row in enumerate(rows))),
        "connected_sector_sizes": components,
        "reference_sector_masses": {str(q): str(sum(
            (p for p, a in zip(reference, charges, strict=True) if a == q), F()))
            for q in sorted(set(charges))},
        "support_preserves_sector": True,
        "faithful_reference_is_normalized_stationary_and_reversible": True,
    }


def report() -> dict:
    campaign = strict_json.load(HERE / "campaign.json")
    required = {
        "schema": "oph.conserved_sector_action.campaign.v1",
        "campaign_id": "conserved-sector-action-20261028-v1",
        "predecessor_git_commit": "7cef80f58f51f06e70ccda3d4373344ca3993dc6",
        "scope": "FINITE_COMPONENT_DIAGNOSTICS",
        "boundary_status": "FROZEN_DIAGNOSTIC_CONTRACT",
        "whole_program_frozen": False,
        "native_family_changed": False,
        "physical_comparison_count_added": 0,
        "integrated_gate_accepted": False,
    }
    if any(type(campaign.get(k)) is not type(v) or campaign[k] != v
           for k, v in required.items()):
        raise ValueError("sector action campaign changed")
    if bound_pins(campaign.get("source_sha256")) != predecessor_bindings():
        raise ValueError("sector action predecessor binding changed")
    kernels = {law: native_kernel(law) for law in ("C", "E")}
    if kernels["C"] != kernels["E"]:
        raise ValueError("C/E joint population kernels differ")
    laws = {}
    for law, (rows, charges) in kernels.items():
        reference = tuple(sector.WEIGHTS[q] / len(sector.CORES) for q in charges)
        laws[law] = certify_kernel(rows, charges, reference)
    quantum = reference_balance.co_moving()
    if not quantum["charge_projectors_preserved_by_each_complete_record_block"]:
        raise ValueError("complete quantum record support changes q")
    return {
        "schema": "oph.conserved_sector_action.results.v1",
        "campaign_id": campaign["campaign_id"],
        "source_sha256": bound_pins(source_bindings()),
        "predecessor_evidence": {
            "frozen_source_path": PREDECESSOR_ENERGY_SOURCE,
            "captured_path": PREDECESSOR_ENERGY_CAPTURE,
            "sha256": PREDECESSOR_ENERGY_SHA256,
            "predecessor_git_commit": campaign["predecessor_git_commit"],
            "role": "historical campaign input; distinct from the current energy receipt",
        },
        "laws": laws,
        "complete_raw_quantum_source_core_pair_checks": quantum["raw_kraus_source_core_pair_checks"],
        "theorem_module": "OPH.ConservedSectorAction",
        "theorem_premises_checked_by": "exact Python rational and discrete enumeration; not a Lean encoding of the full native kernel",
        "theorem_consequences": {
            "arbitrary_positive_sector_reweighting": "pi_f(x,c)=f(q)*pi(x,c)/Z is faithful, normalized and stationary for every positive f",
            "conditional_reference": "conditioning pi_f on a fixed q cancels f and Z; the same conditional density has the same instrument records",
            "sector_potential_gradient": "every U(q) has zero native transition-gradient form",
            "reciprocity_uniqueness": "SourceActionMeasure.ratio_along_walk fixes the reference ratio inside each connected sector; global normalization leaves two independent sector probabilities",
        },
        "scope": "reduced diagonal source/core population law under independent uniform seams and fair coins; full raw record support preserves q",
        "full_archive_stationarity_claimed": False,
        "physical_energy_or_mass_selected": False,
        "physical_comparison_count_added": 0,
        "integrated_gate_accepted": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", nargs="?", const=OUTPUT, type=Path)
    mode.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    args = parser.parse_args()
    result = report()
    if args.write:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    elif strict_json.load(args.check) != result:
        raise ValueError("conserved-sector theorem receipt mismatch")
    print("conserved-sector action theorem premises verified for C/E")


if __name__ == "__main__":
    main()
