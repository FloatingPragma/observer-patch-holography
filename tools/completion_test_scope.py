"""Explicit current-law and independent-component test domains.

The obsolete CP adoption/issue-routing wrappers do not validate CCF-LT-2.
Their component mathematics remains in its owning suites. This selector never
changes a test, input, assertion or pytest collection hook. Every top-level
completion test file must have an explicit domain before this command runs.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

# These validate the current CCF-LT-2 programme, graph and premise domains.
CURRENT = (
    'test_assumption_closure.py',
    'test_theorem_map.py',
)

# Numerical, operator, observer, source-custody and scope controls retain their own laws.
COMPONENTS = (
    'test_acceptance_audit.py',
    'test_adaptive_reference.py',
    'test_defect_algebra.py',
    'test_family.py',
    'test_finite_record_dynamics.py',
    'test_four_lane_joins.py',
    'test_full_field_joins.py',
    'test_modular_access.py',
    'test_observer.py',
    'test_observer_reference.py',
    'test_particle_premise_joins.py',
    'test_source_flavor_comparison.py',
    'test_stationary_reference.py',
)

# Each wrapper projects component evidence into the replaced CP-00..CP-13
# registry or finite_attachment_map; mathematical tests run in owning packages.
COMPONENT_ROUTING = (
    'test_compact_source_campaign.py',
    'test_defect_geometry_campaign.py',
    'test_defect_selection_campaign.py',
    'test_dynamics_response_campaign.py',
    'test_finite_attachment_map.py',
    'test_finite_observable_campaign.py',
    'test_native_observation_campaign.py',
    'test_particle_consistency_campaign.py',
    'test_particle_derivation_campaign.py',
    'test_radial_source_campaign.py',
    'test_source_premise_campaign.py',
    'test_spectrum_einstein_campaign.py',
)

# Keep independent arithmetic/resource checks without re-adopting its signed class.
PARTIAL_COMPONENTS = {
    "test_program_class_v2.py": (
        "test_actual_continuing_observer_counts_fit_declared_finite_request_budget",
        "test_exact_comparison_rule_retains_inconclusive_and_calibrated_cases",
    ),
}


def paths(root=ROOT):
    root = Path(root)
    groups = (CURRENT, COMPONENTS, COMPONENT_ROUTING, tuple(PARTIAL_COMPONENTS))
    assigned = [name for group in groups for name in group]
    discovered = {p.name for p in (root / "model/completion").glob("test_*.py")}
    if len(assigned) != len(set(assigned)) or set(assigned) != discovered:
        raise ValueError("every completion test file needs one explicit source domain")
    selected = ["tools/test_completion_program.py"]
    selected += ["model/completion/" + name for name in (*CURRENT, *COMPONENTS)]
    selected += [f"model/completion/{name}::{case}"
                 for name, cases in PARTIAL_COMPONENTS.items() for case in cases]
    selected.append("model/clock_defect_energy")
    return selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        parser.error("shard index must be in 0..shard-count-1")
    selected = paths()
    if not args.run:
        print("\n".join(selected))
        return
    command = [sys.executable, "-m", "pytest", "-q", "--durations=10", "-W", "error",
               "-p", "model.pytest_shards",
               "--model-shard-index", str(args.shard_index),
               "--model-shard-count", str(args.shard_count), *selected]
    raise SystemExit(subprocess.call(command, cwd=ROOT))


if __name__ == "__main__":
    main()
