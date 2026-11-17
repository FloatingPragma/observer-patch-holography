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
    'test_field_current_measurement.py',
    'test_field_current_measurement_consumers.py',
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

# Existing authenticated custody jobs execute these exact signed review inputs.
RETAINED_CUSTODY = (
    'test_custody_renewal.py',
    'test_field_custody.py',
    'test_field_navigation.py',
    'test_source_records_custody.py',
    'test_source_records_navigation.py',
)

# The existing field-current collector executes every original adapter test.
ROUTED_CONSUMER = (
    'test_native_measurement_parent_adapter.py',
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
    groups = (CURRENT, COMPONENTS, COMPONENT_ROUTING, RETAINED_CUSTODY,
              ROUTED_CONSUMER, tuple(PARTIAL_COMPONENTS))
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


ISOLATED_CUSTODY_TEST = (
    "model/charged_lepton_response/test_packet.py::"
    "test_generated_results_and_runtime_custody"
)


def model_paths(root=ROOT):
    """Cover every model test outside completion without altering collection."""
    root = Path(root)
    model = root / "model"
    # Group by the same top-level packages traversed by `pytest model`.
    # Existing pytest ignore rules still govern descendants of each package.
    outside = {
        "model/" + path.name for path in model.iterdir()
        if path.name not in {"completion", "__pycache__"} and not path.name.startswith(".")
        and (path.is_dir() or (path.name.startswith("test_") and path.suffix == ".py"))
    }
    selected_completion = [name for name in paths(root)
                           if name.startswith("model/completion/")]
    return sorted(outside) + selected_completion


def command(*, root=ROOT, model_suite=False, shard_index=0, shard_count=1, deselect=()):
    if not 0 <= shard_index < shard_count:
        raise ValueError("shard index must be in 0..shard-count-1")
    if (len(deselect) != len(set(deselect))
            or any(name != ISOLATED_CUSTODY_TEST for name in deselect)
            or (deselect and not model_suite)):
        raise ValueError("only the separately executed charged-lepton custody test may be deselected")
    selected = model_paths(root) if model_suite else paths(root)
    argv = [sys.executable, "-m", "pytest", "-q", "--durations=30", "-W", "error",
            "-p", "tools.completion_contract_inputs", "-p", "model.pytest_shards",
            "--model-shard-index", str(shard_index), "--model-shard-count", str(shard_count)]
    for name in deselect:
        argv += ["--deselect", name]
    return argv + selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--model-suite", action="store_true")
    parser.add_argument("--deselect", action="append", default=[])
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    args = parser.parse_args()
    try:
        argv = command(model_suite=args.model_suite, shard_index=args.shard_index,
                       shard_count=args.shard_count, deselect=args.deselect)
    except ValueError as exc:
        parser.error(str(exc))
    if not args.run:
        selected = model_paths() if args.model_suite else paths()
        print("\n".join(selected))
        return
    raise SystemExit(subprocess.call(argv, cwd=ROOT))


if __name__ == "__main__":
    main()
