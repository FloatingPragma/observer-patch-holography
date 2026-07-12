#!/usr/bin/env python3
"""Realized cyclic cap-net repair run for the #503 nonemptiness gate.

This module executes a genuine finite transactional OPH repair run, real
seeded conflicts, transactional majority repair, runtime-verified termination,
confluence, and schedule independence (the node-D1 clauses checked at run
time, as in the #238 rooted-tree export), on a refinement tower of cyclic
cap nets, and exports the repaired record layer as a consensus artifact for
the realized-branch receipt evaluation.

Provenance is stated exactly: the overlap net of each stage is CHOSEN to be a
sphere triangulation (icosahedron plus edge subdivisions). Per the
dimension-selection boundary of the compact paper
(Corollary `cor:dimension-selection-status`), this is explicit branch
selection, allowed as a named branch input, never claimed as a consensus
output. What the run genuinely witnesses is the joint realizability of the
node-D1 repair clauses with the spherical-incidence, mesh, and
refinement-naturality receipts on one tower: the repair layer does not know
the topology, and the incidence complex is recomputed from the repaired
normal-form records alone.

NOT witnessed here: the modular cross-ratio receipt and the geometric
2pi-KMS receipt, which require modular-clock instrumentation of the MaxEnt
states on this tower and are pending; and the E/MI/scale/UC/VR families.

Run:
    python3 code/geometry/cyclic_cap_net_run.py
writes code/geometry/runs/cyclic_cap_net_run_domain.json.
"""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from quotient_cap_readout import (  # noqa: E402
    RepairSystem,
    classify_surface,
    euler_characteristic,
    icosahedron,
    incidence_complex,
    is_closed_surface,
    is_connected,
    normal_form_records,
    orient,
    refinement_is_simplicial,
    refinement_subdivide,
    spherical_incidence_receipt,
)

ARTIFACT_PATH = HERE / "runs" / "cyclic_cap_net_run_domain.json"
N_SCHEDULES = 24
N_CONFLICT_SEEDS = 5


def build_tower(stages: int = 3) -> list[list[tuple[int, ...]]]:
    """Record layers of the refinement tower: icosahedron plus subdivisions."""
    records, _ = icosahedron()
    tower = [list(records)]
    for _ in range(stages - 1):
        finer, _ = refinement_subdivide(tower[-1])
        tower.append(finer)
    return tower


def run_repair_stage(records: list[tuple[int, ...]], stage: int,
                     rng: random.Random) -> dict:
    """Execute the transactional repair with runtime D1 verification.

    Several independent conflict seeds are repaired under many random
    schedules; the stage passes only if every (seed, schedule) pair
    terminates in the same conflict-free quotient normal form for that seed.
    """
    verification = {
        "n_records": len(records),
        "n_conflict_seeds": N_CONFLICT_SEEDS,
        "n_schedules_per_seed": N_SCHEDULES,
        "terminating": True,
        "conflict_free_normal_form": True,
        "schedule_independent": True,
    }
    reference_nf = None
    for seed_record in range(min(N_CONFLICT_SEEDS, len(records))):
        system = RepairSystem(records, seed_record=seed_record)
        reference = system.repair()
        # conflict-free check
        probe = RepairSystem(records, seed_record=seed_record)
        probe.state = dict(reference)
        if probe.conflicted_records():
            verification["conflict_free_normal_form"] = False
        for _ in range(N_SCHEDULES):
            schedule = list(range(len(records)))
            rng.shuffle(schedule)
            if system.repair(schedule) != reference:
                verification["schedule_independent"] = False
        if reference_nf is None:
            reference_nf = reference
    verification["stage"] = stage
    return verification


def mesh_receipt(records: list[tuple[int, ...]]) -> dict:
    """Mesh nondegeneracy modulus: caps are vertex stars; theta_r is the
    largest star size as a fraction of the record layer, which must decrease
    along refinement (cofinal nondegenerate cap mesh)."""
    k = incidence_complex(records)
    star_sizes = []
    for v in k.vertices:
        star_sizes.append(sum(1 for t in k.triangles if v in t))
    return {
        "max_star_size": max(star_sizes),
        "n_triangles": len(k.triangles),
        "theta": max(star_sizes) / len(k.triangles),
    }


def run_tower(stages: int = 3, seed: int = 20260712) -> dict:
    rng = random.Random(seed)
    tower = build_tower(stages)
    stage_reports = []
    prev_records = None
    prev_complex = None
    for stage, records in enumerate(tower):
        d1 = run_repair_stage(records, stage, rng)
        # the incidence readout consumes ONLY the repaired normal-form records
        system = RepairSystem(records, seed_record=0)
        nf = system.repair()
        nf_records = normal_form_records(system, nf)
        k = incidence_complex(nf_records)
        report = {
            "d1_verification": d1,
            "incidence": {
                "n_vertices": len(k.vertices),
                "n_edges": len(k.edges),
                "n_triangles": len(k.triangles),
                "euler_characteristic": euler_characteristic(k),
                "connected": is_connected(k),
                "closed_surface": is_closed_surface(k),
                "orientable_coherent": orient(k) is not None,
                "spherical_incidence_receipt": spherical_incidence_receipt(k),
                "surface_classification": classify_surface(k),
            },
            "mesh_receipt": mesh_receipt(nf_records),
        }
        if prev_records is not None:
            _, projection = refinement_subdivide(prev_records)
            report["refinement_natural"] = refinement_is_simplicial(
                k, prev_complex, projection
            )
        stage_reports.append(report)
        prev_records = nf_records
        prev_complex = k
    thetas = [r["mesh_receipt"]["theta"] for r in stage_reports]
    return {
        "artifact": "oph_cyclic_cap_net_repair_run",
        "object_id": "CyclicCapNetRepairTower_Issue503",
        "issue": 503,
        "provenance": (
            "consensus_product_on_selected_branch: the repair, termination, "
            "confluence, schedule-independence, and incidence readout are "
            "genuine run outputs; the overlap net is CHOSEN spherical "
            "(explicit branch selection per cor:dimension-selection-status), "
            "never claimed as a consensus output"
        ),
        "stages": stage_reports,
        "receipts_witnessed": {
            "d1_repair_clauses": all(
                r["d1_verification"]["terminating"]
                and r["d1_verification"]["conflict_free_normal_form"]
                and r["d1_verification"]["schedule_independent"]
                for r in stage_reports
            ),
            "spherical_incidence_all_stages": all(
                r["incidence"]["spherical_incidence_receipt"]
                for r in stage_reports
            ),
            "mesh_modulus_decreasing": all(
                thetas[i] > thetas[i + 1] for i in range(len(thetas) - 1)
            ),
            "refinement_naturality": all(
                r.get("refinement_natural", True) for r in stage_reports
            ),
        },
        "receipts_pending": [
            "modular_cross_ratio_receipt (needs modular-clock instrumentation "
            "of the MaxEnt states on this tower)",
            "geometric_2pi_kms_receipt (same instrumentation)",
            "Cyc/NTI/weak-additivity/MI null-net receipts",
            "E1-E4 event receipts",
            "UC/VR/scale physical-identification receipts",
        ],
    }


def main() -> None:
    report = run_tower()
    ARTIFACT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(ARTIFACT_PATH, "w") as f:
        json.dump(report, f, indent=2)
        f.write("\n")
    w = report["receipts_witnessed"]
    print(json.dumps(w, indent=2))
    print("pending:", len(report["receipts_pending"]), "receipt families")


if __name__ == "__main__":
    main()
