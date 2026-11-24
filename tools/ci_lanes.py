#!/usr/bin/env python3
"""Lane table for the push CI (ci.yml) and the nightly replays (nightly.yml).

Each lane is one independent replay: shell commands run from the repository
root with the lane's environment. A lane belongs to one tier ("push" runs on
every push and pull request; "nightly" is reserved for the schedule and holds no lane) and one install
profile (the requirements file its job installs). The runner packs the lanes
of one tier and profile into a fixed number of buckets by their measured
minutes, so a workflow matrix job only names its bucket:

    python tools/ci_lanes.py run --tier push --profile root --bucket 0 --buckets 8
    python tools/ci_lanes.py run --lane m1_quantum_transport
    python tools/ci_lanes.py list --tier push --profile root --buckets 8

Every lane in a bucket runs even after an earlier one fails; the bucket exits
non-zero when any lane failed and prints a summary with durations.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

BLAS1 = {"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}
CODE = {"PYTHONPATH": "code"}


@dataclass(frozen=True)
class Lane:
    name: str
    tier: str
    profile: str
    minutes: float
    commands: tuple[str, ...]
    env: dict[str, str] = field(default_factory=dict)


def lane(name, tier, profile, minutes, *commands, env=None):
    return Lane(name, tier, profile, minutes, tuple(commands), dict(env or {}))


def pytest_lane(name, minutes, path, *, tier="push", profile="root", env=None, flags="-q"):
    return lane(name, tier, profile, minutes, f"python -m pytest {flags} {path}", env=env)


LANES: tuple[Lane, ...] = (
    # ------------------------------------------------------------------
    # Push tier, root requirements: code/ replays and their controls.
    pytest_lane("capacity_readback", 1, "code/capacity_readback -W error", env={"OPENBLAS_NUM_THREADS": "1"}),
    lane("causal_poset_archive", "push", "root", 2,
         "python evidence/source_net_causal_poset/verify_causal_poset_archive.py", env=BLAS1),
    pytest_lane("collar_gap", 1, "code/yang_mills/test_collar_gap_certificate.py "
                "code/yang_mills/tests/test_finite_collar_gap_certificate.py "
                "code/yang_mills/tests/test_collar_certificate_accounting.py tools/test_strict_json.py",
                flags="-q -W error"),
    lane("consensus", "push", "root", 1.5, "python -W error -m pytest -q code/consensus", env=BLAS1),
    lane("cosmology_accounting", "push", "root", 0.5,
         "python -m pytest -q code/cosmology/fixed_capacity_wlaw code/cosmology/postdiction_ledger -W error",
         "python code/cosmology/postdiction_ledger/verify_cosmology_postdiction_ledger_independent.py"),
    lane("finite_gauge_transfer", "push", "root", 2,
         "python -W error -m pytest -q code/yang_mills/tests/test_z2_finite_transfer_receipt.py "
         "code/yang_mills/tests/test_z2_finite_transfer_verifier.py code/yang_mills/tests/test_z2_transfer_spectrum.py "
         "code/yang_mills/tests/test_z2_kogut_precision.py code/yang_mills/tests/test_z2_full_space_controls.py "
         "code/yang_mills/tests/test_kogut_susskind_fiber_rate_instances.py code/yang_mills/tests/test_su2_cubic_gap.py "
         "code/yang_mills/tests/test_loop_vacuum_controls.py",
         "python -W error -m pytest -q code/edge_sectors tools/test_post_r2029_audit_surfaces.py", env=BLAS1),
    pytest_lane("geometry", 2, "code/geometry", env={"OPENBLAS_NUM_THREADS": "1"}),
    lane("k7_gauge_kinetic_cross", "push", "root", 1,
         "python code/angular_sprint/k7_gauge_kinetic_cross_control.py --check",
         "python code/angular_sprint/verify_k7_gauge_kinetic_cross_control_independent.py "
         "--receipt code/angular_sprint/runtime/k7_gauge_kinetic_cross_control.json",
         "python code/angular_sprint/kinetic_form_selection_certificate.py --verify",
         "python code/angular_sprint/verify_kinetic_form_selection_independent.py "
         "--receipt code/angular_sprint/runtime/kinetic_form_selection_receipt.json"),
    pytest_lane("m1_interfaces", 1, "code/m1_interfaces"),
    pytest_lane("m1_necessity", 1.5, "code/m1_necessity"),
    pytest_lane("m1_operational_clocks", 1, "code/m1_operational_clocks"),
    pytest_lane("m1_quantum_transport", 2.5, "code/m1_quantum_transport"),
    pytest_lane("m1_source_realization", 1.5, "code/m1_source_realization"),
    pytest_lane("m1_vacuum_fidelity", 1, "code/m1_vacuum_fidelity", flags="-q -W error"),
    lane("mckay_golden_field", "push", "root", 1.5,
         'python -B code/a5_closure/mckay_golden_field_certificate.py --output "$LANE_TMP/mckay-golden-raw.json"',
         'python -B code/a5_closure/verify_mckay_golden_field_independent.py --input "$LANE_TMP/mckay-golden-raw.json" '
         '--output "$LANE_TMP/mckay-golden-verified.json"',
         "python -B code/a5_closure/test_mckay_golden_field_certificate.py",
         "python -B code/a5_closure/test_mckay_golden_field_cross_control.py",
         "python -B code/a5_closure/verify_mckay_golden_field_cross_control.py"),
    lane("neutral_packet", "push", "root", 2,
         "python -m pytest -q -W error code/electromagnetism/test_neutral_packet_projection.py "
         "code/electromagnetism/test_neutral_packet_observables.py code/electromagnetism/test_whitney_quantum_packet.py",
         "python -m pytest -q -W error code/electromagnetism/test_whitney_interacting_quantum.py",
         env={"OPENBLAS_NUM_THREADS": "1"}),
    lane("private_links", "push", "root", 0.5, "python tools/check_private_links.py"),
    pytest_lane("quantum_information", 2.5,
                "code/quantum_information code/collar_alignment code/maxent",
                env={**CODE, "OPENBLAS_NUM_THREADS": "1"}),
    lane("radial_inverse", "push", "root", 1,
         "python -W error -m pytest -q code/cosmology/test_oph_radial_lift_330.py", env=BLAS1),
    pytest_lane("rg_principle", 1, "code/rg_principle"),
    lane("ringdown_arithmetic", "push", "root", 1, "python -W error -m pytest -q code/gravitation/ringdown_comb"),
    lane("source_selection_model", "push", "root", 6.5,
         "python -m source_selection_model.verify", "python -m source_record_gluing.verify",
         "python -m pytest -q code/source_selection_model code/source_record_gluing", env=CODE),
    lane("certificate_smoke", "push", "root", 3,
         "python tools/run_mandatory_suite.py --certificate-smoke-only",
         "python -m pytest -q code/a5_closure/tests/test_source_repair_generator_certificate.py"),

    # ------------------------------------------------------------------
    # Push tier, model requirements (model/requirements.txt): programs of the
    # main candidate law and the theory certificates that use the same pins.
    lane("support_refinement", "push", "model", 2,
         "python -m pytest -q -W error code/support_refinement",
         "python -m support_refinement --check code/support_refinement/results.json", env={**CODE, **BLAS1}),
    lane("matched_einstein_limit", "push", "model", 1,
         "python -m pytest -q -W error code/matched_einstein_limit", env={**CODE, **BLAS1}),
    lane("defect_particle_classes", "push", "model", 1,
         "python -m defect_particle_classes --check",
         "python -m pytest -q -W error code/defect_particle_classes", env={**CODE, **BLAS1}),
    lane("maxwell_certificates", "push", "model", 1,
         "python -m pytest -q -W error code/maxwell_local_selection code/maxwell_weak_field", env={**CODE, **BLAS1}),
    lane("selection_certificates", "push", "model", 1,
         "python -m pytest -q -W error code/carrier_selection_reduction code/finite_instrument_calibration",
         env={**CODE, **BLAS1}),
)


def select(tier: str | None, profile: str | None) -> list[Lane]:
    return [x for x in LANES if (tier is None or x.tier == tier) and (profile is None or x.profile == profile)]


def pack(lanes: list[Lane], buckets: int) -> list[list[Lane]]:
    """Longest-processing-time packing; deterministic for a fixed lane table."""
    if buckets < 1:
        raise SystemExit("--buckets must be positive")
    bins: list[list[Lane]] = [[] for _ in range(buckets)]
    loads = [0.0] * buckets
    for item in sorted(lanes, key=lambda x: (-x.minutes, x.name)):
        index = min(range(buckets), key=lambda i: (loads[i], i))
        bins[index].append(item)
        loads[index] += item.minutes
    return bins


def run_lane(item: Lane, tmp_root: Path) -> tuple[bool, float]:
    tmp = tmp_root / item.name
    tmp.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update(item.env)
    env["LANE_TMP"] = str(tmp)
    start = time.monotonic()
    ok = True
    print(f"::group::lane {item.name}", flush=True)
    for command in item.commands:
        print(f"$ {command}", flush=True)
        result = subprocess.run(["bash", "-c", "set -euo pipefail\n" + command], cwd=ROOT, env=env)
        if result.returncode != 0:
            print(f"::error::lane {item.name} failed (exit {result.returncode}): {command}", flush=True)
            ok = False
            break
    print("::endgroup::", flush=True)
    return ok, time.monotonic() - start


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", choices=("run", "list"))
    parser.add_argument("--tier", choices=("push", "nightly"))
    parser.add_argument("--profile", choices=("root", "model"))
    parser.add_argument("--bucket", type=int)
    parser.add_argument("--buckets", type=int, default=1)
    parser.add_argument("--lane", action="append", default=[])
    args = parser.parse_args()

    names = [x.name for x in LANES]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        raise SystemExit("duplicate lane names: " + ", ".join(duplicates))

    if args.lane:
        by_name = {x.name: x for x in LANES}
        unknown = [n for n in args.lane if n not in by_name]
        if unknown:
            raise SystemExit("unknown lanes: " + ", ".join(unknown))
        groups = [[by_name[n] for n in args.lane]]
    else:
        lanes = select(args.tier, args.profile)
        if not lanes:
            raise SystemExit("no lanes match the selection")
        groups = pack(lanes, args.buckets)
        if args.bucket is not None:
            if not 0 <= args.bucket < args.buckets:
                raise SystemExit("--bucket must satisfy 0 <= bucket < --buckets")
            groups = [groups[args.bucket]]

    if args.action == "list":
        for index, group in enumerate(groups):
            load = sum(x.minutes for x in group)
            print(f"bucket {index}: {load:.1f} min: " + ", ".join(x.name for x in group))
        return

    tmp_root = Path(os.environ.get("RUNNER_TEMP") or ROOT / "temp") / "ci-lanes"
    results = []
    for group in groups:
        for item in group:
            ok, seconds = run_lane(item, tmp_root)
            results.append((item.name, ok, seconds))
    print("\nlane summary:")
    for name, ok, seconds in results:
        print(f"  {'PASS' if ok else 'FAIL'} {seconds:7.1f}s {name}")
    failed = [name for name, ok, _ in results if not ok]
    if failed:
        raise SystemExit("failed lanes: " + ", ".join(failed))


if __name__ == "__main__":
    main()
