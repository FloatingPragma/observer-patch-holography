"""Semantic mutation tests for the B20 random-scan preflight validator.

The committed certificate must validate, and each semantic tamper class
must be rejected offline (no run directory needed): a subset promoted
to passing, a component count altered, a fixed-space dimension altered
(either scheduler), a non-idempotence witness altered, a visit-count
edit, a state-label edit, a verdict flip, a reference edit, and an
emptied extra-field inventory.  The connectivity and elimination
routines are additionally cross-checked against independent oracles on
synthetic instances, because the committed data exercises only the
one-component/dimension-one case.
"""

from __future__ import annotations

import importlib.util
import json
import os
import random
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
CERT = HERE / "runtime" / "b20_preflight_certificate.json"
VALIDATOR = HERE / "validate_random_scan.py"
PRODUCER = HERE / "preflight_random_scan.py"


def load(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_validator(path: Path, offline: bool = True) -> subprocess.CompletedProcess:
    # Mutation tests run with the run directory redirected to a
    # nonexistent path, proving each tamper class is rejected by the
    # offline algebra layer alone.
    env = dict(os.environ)
    if offline:
        env["B20_RUN_DIR"] = "/nonexistent-b20-run-dir"
    else:
        env.pop("B20_RUN_DIR", None)
    return subprocess.run(
        [sys.executable, str(VALIDATOR), str(path)],
        capture_output=True,
        text=True,
        env=env,
    )


def mutate(tmp_path: Path, transform) -> Path:
    cert = json.loads(CERT.read_text())
    transform(cert)
    out = tmp_path / "mutated.json"
    out.write_text(json.dumps(cert, indent=1, sort_keys=True) + "\n")
    return out


def test_committed_certificate_validates_offline():
    result = run_validator(CERT, offline=True)
    assert result.returncode == 0, result.stderr
    assert "custody layer skipped" in result.stdout


def test_committed_certificate_validates_with_custody():
    run_dir = Path(__file__).resolve().parents[3] / "oph-physics-sim" / \
        "runs" / "b12_prereg_16k_20260806"
    if not run_dir.is_dir():
        import pytest

        pytest.skip("pinned run directory not present")
    result = run_validator(CERT, offline=False)
    assert result.returncode == 0, result.stderr
    assert "custody layer verified" in result.stdout


def test_promoted_pass_flag_rejected(tmp_path):
    def transform(cert):
        cert["subset_results"][0]["passes_all_requirements"] = True
        cert["subset_results"][0]["nonconstant_protected_observable"] = True

    result = run_validator(mutate(tmp_path, transform))
    assert result.returncode != 0


def test_component_count_tamper_rejected(tmp_path):
    def transform(cert):
        cert["subset_results"][0]["join_component_count"] = 2

    result = run_validator(mutate(tmp_path, transform))
    assert result.returncode != 0


def test_dimension_tamper_rejected(tmp_path):
    def transform(cert):
        cert["subset_results"][0]["fixed_space_dimension"] = 2

    result = run_validator(mutate(tmp_path, transform))
    assert result.returncode != 0


def test_second_scheduler_tamper_rejected(tmp_path):
    def transform(cert):
        entry = cert["subset_results"][0]["second_scheduler"]
        entry["fixed_space_dimension"] = 2

    result = run_validator(mutate(tmp_path, transform))
    assert result.returncode != 0


def test_witness_tamper_rejected(tmp_path):
    def transform(cert):
        cert["subset_results"][0]["non_idempotence_witness"]["mix_entry"] = "0"

    result = run_validator(mutate(tmp_path, transform))
    assert result.returncode != 0


def test_verdict_flip_rejected(tmp_path):
    def transform(cert):
        cert["verdict"] = "positive"
        cert["designated_subset"] = ["checkpoint_class", "stable_flag"]

    result = run_validator(mutate(tmp_path, transform))
    assert result.returncode != 0


def test_reference_edit_rejected(tmp_path):
    def transform(cert):
        cert["reference"][0] = "1/2"

    result = run_validator(mutate(tmp_path, transform))
    assert result.returncode != 0


def test_extra_field_inventory_tamper_rejected(tmp_path):
    def transform(cert):
        cert["extra_step_fields_realized_values"] = {}

    result = run_validator(mutate(tmp_path, transform))
    assert result.returncode != 0


def test_visit_count_tamper_rejected(tmp_path):
    def transform(cert):
        cert["visit_counts"][0] += 32

    result = run_validator(mutate(tmp_path, transform))
    assert result.returncode != 0


def test_state_label_tamper_rejected(tmp_path):
    def transform(cert):
        label = json.loads(cert["state_labels"][0])
        label[-1][1] = 99  # move one state into a new repair-load block
        cert["state_labels"][0] = json.dumps(label)

    result = run_validator(mutate(tmp_path, transform))
    assert result.returncode != 0


def test_arena2_mass_tamper_rejected(tmp_path):
    def transform(cert):
        cert["arena2_states"][0][2] += 1

    result = run_validator(mutate(tmp_path, transform))
    assert result.returncode != 0


def test_grammar_prose_tamper_rejected(tmp_path):
    def transform(cert):
        cert["no_go_grammar"] = "all schedulers and references"

    result = run_validator(mutate(tmp_path, transform))
    assert result.returncode != 0


def test_next_route_foreclosure_tamper_rejected(tmp_path):
    def transform(cert):
        cert["recorded_next_route"] = "no viable route remains"

    result = run_validator(mutate(tmp_path, transform))
    assert result.returncode != 0


def test_offline_layer_runs_without_run_directory():
    # The offline algebra layer must fully verify the vendored objects
    # even when the pinned run directory is absent, and the
    # --require-custody flag must then fail closed.
    module = load(VALIDATOR)
    module.RUN_DIR = Path("/nonexistent-b20-run-dir")
    cert = json.loads(CERT.read_text())
    module.offline_layer(cert)  # must not raise

    argv = [str(VALIDATOR), str(CERT), "--require-custody"]
    old_argv = sys.argv
    sys.argv = argv
    try:
        try:
            module.main()
        except SystemExit as exc:
            assert exc.code == 1
        else:
            raise AssertionError("--require-custody must fail without run dir")
    finally:
        sys.argv = old_argv


def union_find_oracle(fibres, n):
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for f in fibres:
        first = {}
        for x in range(n):
            if f[x] in first:
                ra, rb = find(x), find(first[f[x]])
                if ra != rb:
                    parent[ra] = rb
            else:
                first[f[x]] = x
    return len({find(x) for x in range(n)})


def test_component_routine_matches_union_find_oracle():
    # The committed data yields one component everywhere, so this guards
    # the connectivity routine against a degenerate always-one failure.
    module = load(VALIDATOR)
    assert module.components_multi(
        [[0, 0, 1, 2, 2, 3], [0, 1, 1, 2, 3, 3]], 6) == 2
    assert module.components_multi(
        [list(range(6)), list(range(6))], 6) == 6
    rng = random.Random(7)
    for _ in range(200):
        n = rng.randrange(2, 40)
        k = rng.randrange(2, 5)
        fibres = [
            [rng.randrange(1, 8) for _ in range(n)] for _ in range(k)
        ]
        assert module.components_multi(fibres, n) == union_find_oracle(
            fibres, n
        )


def test_fixed_dimension_routines_on_synthetic_instances():
    # Both elimination routines (producer and validator use different
    # pivot orders) must find dimension >= 2 on a genuinely disconnected
    # instance and agree with the component count on random instances.
    validator = load(VALIDATOR)
    producer = load(PRODUCER)
    rng = random.Random(11)
    for _ in range(25):
        n = rng.randrange(2, 12)
        k = rng.randrange(2, 4)
        fibres = [
            [rng.randrange(0, 4) for _ in range(n)] for _ in range(k)
        ]
        pi = [Fraction(rng.randrange(1, 9)) for _ in range(n)]
        total = sum(pi, Fraction(0))
        pi = [p / total for p in pi]
        kernels = [producer.heat_bath(pi, f) for f in fibres]
        mix = [
            [
                sum((ker[x][y] for ker in kernels), Fraction(0)) / k
                for y in range(n)
            ]
            for x in range(n)
        ]
        comp = union_find_oracle(fibres, n)
        assert producer.fixed_space_dimension(mix, n) == comp
        assert validator.fixed_dimension(mix, n) == comp
    # A hand-built disconnected instance: two islands, dimension two.
    fibres = [[0, 0, 1, 1], [2, 2, 3, 3]]
    pi = [Fraction(1, 4)] * 4
    kernels = [producer.heat_bath(pi, f) for f in fibres]
    mix = [
        [(kernels[0][x][y] + kernels[1][x][y]) / 2 for y in range(4)]
        for x in range(4)
    ]
    assert producer.fixed_space_dimension(mix, 4) == 2
    assert validator.fixed_dimension(mix, 4) == 2
