from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest


PROGRAMS = Path(__file__).resolve().parents[1] / "programs"
MODULE_PATH = PROGRAMS / "record_gated_cayley_runtime.py"
SPEC = importlib.util.spec_from_file_location("record_gated_cayley_runtime", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
runtime = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = runtime
SPEC.loader.exec_module(runtime)


class FakeCircuit:
    def __init__(self, blob: bytes = b"logical") -> None:
        self.blob = blob


def fake_qpy_dump(circuits, handle) -> None:
    handle.write(b"|".join(circuit.blob for circuit in circuits))


def make_bundle_files(tmp_path: Path):
    opaque_id = "opaque-circuit-1"
    heldout_opaque_id = "opaque-circuit-2"
    public_role = "role-opaque-development"
    heldout_public_role = "role-opaque-heldout"
    public_backend = "backend-opaque-fez"
    heldout_public_backend = "backend-opaque-kingston"
    public_layout = "layout-opaque-fez"
    heldout_public_layout = "layout-opaque-kingston"
    qpy_sha = hashlib.sha256(b"logical").hexdigest()
    resources = {
        "account_quota_seconds": 600,
        "account_reserve_seconds": 180,
        "estimated_qpu_seconds_ceiling": 420,
        "max_execution_time_seconds": 420,
        "max_pubs_per_job": 300,
        "shots_per_circuit": 256,
        "backend_slot_count": 2,
        "logical_qubits_per_circuit": 4,
        "classical_bits_per_circuit": 7,
        "catalog_entries_per_backend_slot": 1,
        "total_circuit_instances": 2,
        "total_planned_shots": 512,
        "transpile_optimization_level": 1,
        "transpiler_seed": 509,
        "resilience_level": 0,
        "dynamical_decoupling": False,
        "gate_twirling": False,
        "measurement_twirling": False,
        "common_duration_padding": True,
    }
    manifest = {
        "schema_version": runtime.MANIFEST_SCHEMA,
        "protocol_id": "opaque-protocol",
        "backend_slots": [
            {"role": public_role, "backend": public_backend, "layout": public_layout},
            {
                "role": heldout_public_role,
                "backend": heldout_public_backend,
                "layout": heldout_public_layout,
            },
        ],
        "circuits": [
            {
                "opaque_id": opaque_id,
                "qpy_sha256": qpy_sha,
                "shots": 256,
                "backend_role": public_role,
            },
            {
                "opaque_id": heldout_opaque_id,
                "qpy_sha256": qpy_sha,
                "shots": 256,
                "backend_role": heldout_public_role,
            },
        ],
        "resources": resources,
        "thresholds": {},
        "nulls": [],
        "exclusions": [],
    }
    catalog_sha = runtime.sha256_json(manifest)
    manifest["catalog_precommitment_sha256"] = catalog_sha
    analysis = {
        "schema_version": "analysis.v1",
        "blind_manifest_commitment": catalog_sha,
        "analysis": "frozen",
    }
    analysis_sha = runtime.sha256_json(analysis)
    analysis["analysis_lock_sha256"] = analysis_sha
    manifest["analysis_lock_sha256"] = analysis_sha

    core_sha = runtime.sha256_json(manifest)
    descriptor = {
        "family": "cayley",
        "backend_role": "development",
        "backend_role_opaque_id": public_role,
        "shots": 256,
        "qpy_sha256": qpy_sha,
        "circuit_name": opaque_id,
        "circuit_metadata": {},
        "parameters": {},
    }
    heldout_descriptor = dict(descriptor)
    heldout_descriptor.update(
        {
            "backend_role": "heldout",
            "backend_role_opaque_id": heldout_public_role,
            "circuit_name": heldout_opaque_id,
        }
    )
    payload = {
        "public_manifest_core_sha256": core_sha,
        "protocol_id": "opaque-protocol",
        "backend_slots": [
            {
                "role": "development",
                "role_opaque_id": public_role,
                "backend": "ibm_fez",
                "backend_opaque_id": public_backend,
                "layout": [10, 18, 94, 124],
                "layout_opaque_id": public_layout,
                "properties_last_update": "2026-07-11T10:15:57+07:00",
            },
            {
                "role": "heldout",
                "role_opaque_id": heldout_public_role,
                "backend": "ibm_kingston",
                "backend_opaque_id": heldout_public_backend,
                "layout": [21, 47, 50, 107],
                "layout_opaque_id": heldout_public_layout,
                "properties_last_update": "2026-07-11T09:49:05+07:00",
            },
        ],
        "mappings": {},
        "circuits": {opaque_id: descriptor, heldout_opaque_id: heldout_descriptor},
    }
    secret_hex = "11" * 32
    commitment = runtime.sha256_bytes(
        runtime.COMMITMENT_DOMAIN
        + bytes.fromhex(secret_hex)
        + b"\0"
        + runtime.canonical_json_bytes(payload)
    )
    manifest["secret_commitment_sha256"] = commitment
    manifest["manifest_sha256"] = runtime.sha256_json(
        runtime._without(manifest, "manifest_sha256")
    )
    reveal = {
        "schema_version": runtime.REVEAL_SCHEMA,
        "secret_hex": secret_hex,
        "sealed_payload": payload,
    }
    manifest_path = tmp_path / "manifest.json"
    reveal_path = tmp_path / "reveal.json"
    analysis_path = tmp_path / "analysis.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    reveal_path.write_text(json.dumps(reveal), encoding="utf-8")
    analysis_path.write_text(json.dumps(analysis), encoding="utf-8")
    return manifest_path, reveal_path, analysis_path


def test_operator_bundle_verifies_full_digest_dag_without_exposing_secret(tmp_path: Path) -> None:
    manifest_path, reveal_path, analysis_path = make_bundle_files(tmp_path)
    bundle = runtime.verify_operator_bundle(
        manifest_path,
        reveal_path,
        analysis_path,
        circuit_rebuilder=lambda _opaque_id, _descriptor: FakeCircuit(),
        ideal_validator=lambda rows: {
            "passed": True,
            "circuits_checked": len(rows),
            "families": {"cayley": len(rows)},
        },
        qpy_dump=fake_qpy_dump,
    )
    assert bundle.backend_slots[0]["role"] == "development"
    assert bundle.backend_slots[0]["backend"] == "ibm_fez"
    assert bundle.circuit_families == {
        "opaque-circuit-1": "cayley",
        "opaque-circuit-2": "cayley",
    }
    assert "11" * 32 not in json.dumps(dataclasses_to_json(bundle))


def dataclasses_to_json(bundle) -> dict:
    return {
        "manifest_sha256": bundle.manifest_sha256,
        "analysis_lock_sha256": bundle.analysis_lock_sha256,
        "backend_slots": bundle.backend_slots,
        "ideal_validation": bundle.ideal_validation,
    }


def test_operator_bundle_rejects_tampered_analysis_lock(tmp_path: Path) -> None:
    manifest_path, reveal_path, analysis_path = make_bundle_files(tmp_path)
    document = json.loads(analysis_path.read_text(encoding="utf-8"))
    document["analysis"] = "changed"
    analysis_path.write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(runtime.RuntimeSafetyError, match="analysis-lock"):
        runtime.verify_operator_bundle(
            manifest_path,
            reveal_path,
            analysis_path,
            circuit_rebuilder=lambda _opaque_id, _descriptor: FakeCircuit(),
            ideal_validator=lambda rows: {"passed": True, "circuits_checked": len(rows)},
            qpy_dump=fake_qpy_dump,
        )


def test_submission_requires_flag_and_both_literal_hashes() -> None:
    bundle = SimpleNamespace(manifest_sha256="aa" * 32, analysis_lock_sha256="bb" * 32)
    with pytest.raises(runtime.RuntimeSafetyError, match="confirm-submit"):
        runtime.validate_submission_confirmation(
            confirm_submit=False,
            confirmed_manifest_sha256="aa" * 32,
            confirmed_analysis_lock_sha256="bb" * 32,
            bundle=bundle,
        )
    with pytest.raises(runtime.RuntimeSafetyError, match="manifest digest"):
        runtime.validate_submission_confirmation(
            confirm_submit=True,
            confirmed_manifest_sha256="cc" * 32,
            confirmed_analysis_lock_sha256="bb" * 32,
            bundle=bundle,
        )
    runtime.validate_submission_confirmation(
        confirm_submit=True,
        confirmed_manifest_sha256="aa" * 32,
        confirmed_analysis_lock_sha256="bb" * 32,
        bundle=bundle,
    )


def test_qpu_guard_preserves_preregistered_reserve() -> None:
    resources = {
        "account_quota_seconds": 600,
        "account_reserve_seconds": 180,
        "estimated_qpu_seconds_ceiling": 420,
    }
    usage = {
        "usage_limit_seconds": 600,
        "usage_remaining_seconds": 600,
        "usage_limit_reached": False,
    }
    assert runtime.guard_qpu_usage(usage, resources, 419)["headroom_after_estimate_seconds"] == 181
    with pytest.raises(runtime.RuntimeSafetyError, match="ceiling"):
        runtime.guard_qpu_usage(usage, resources, 421)
    usage["usage_remaining_seconds"] = 300
    with pytest.raises(runtime.RuntimeSafetyError, match="reserve"):
        runtime.guard_qpu_usage(usage, resources, 121)


class FakeTarget(dict):
    dt = 1e-9
    operation_names = ("x", "measure", "reset", "if_else")


class FakeQubit:
    def __init__(self, index: int) -> None:
        self.index = index


class FakeInstruction:
    def __init__(self, operation, qubits) -> None:
        self.operation = operation
        self.qubits = tuple(qubits)


class FakeDurationCircuit:
    def __init__(self, data, qubits) -> None:
        self.data = data
        self.qubits = qubits

    def find_bit(self, qubit):
        return SimpleNamespace(index=qubit.index)


class FakeBackend:
    name = "ibm_fez"
    num_qubits = 156
    supported_instructions = ("if_else", "measure", "reset")
    default_rep_delay = 0.00025

    def __init__(self) -> None:
        props = SimpleNamespace(duration=2e-6, error=0.0)
        self.target = FakeTarget(x={(0,): props})

    def status(self):
        return SimpleNamespace(operational=True)

    def configuration(self):
        return SimpleNamespace(
            simulator=False,
            n_qubits=156,
            supported_instructions=self.supported_instructions,
        )


def test_recursive_duration_estimator_bounds_durationless_if_else() -> None:
    q0 = FakeQubit(0)
    branch = FakeDurationCircuit(
        [FakeInstruction(SimpleNamespace(name="x", blocks=()), [q0])], [q0]
    )
    conditional = SimpleNamespace(name="if_else", blocks=(branch,))
    circuit = FakeDurationCircuit(
        [
            FakeInstruction(SimpleNamespace(name="x", blocks=()), [q0]),
            FakeInstruction(conditional, [q0]),
        ],
        [q0],
    )
    duration = runtime.estimate_compiled_duration(circuit, FakeBackend())
    assert duration == pytest.approx(104e-6)


def test_dynamic_backend_validation_rejects_missing_control_flow() -> None:
    backend = FakeBackend()
    runtime.validate_dynamic_backend(backend, "ibm_fez", [0, 1, 2, 3])
    backend.supported_instructions = ("measure", "reset")
    backend.target.operation_names = ("x", "measure", "reset")
    with pytest.raises(runtime.RuntimeSafetyError, match="dynamic-circuit"):
        runtime.validate_dynamic_backend(backend, "ibm_fez", [0, 1, 2, 3])


def test_hash_chained_event_journal_detects_mutation(tmp_path: Path) -> None:
    path = tmp_path / "events.ndjson"
    first = runtime.append_event(path, {"event_type": "one"})
    second = runtime.append_event(path, {"event_type": "two"})
    assert second["previous_event_sha256"] == first["event_sha256"]
    rows = path.read_text(encoding="utf-8").splitlines()
    tampered = json.loads(rows[-1])
    tampered["event_type"] = "changed"
    rows[-1] = json.dumps(tampered)
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    with pytest.raises(runtime.RuntimeSafetyError, match="hash chain"):
        runtime.append_event(path, {"event_type": "three"})


class FakeSamplerOptions:
    def __init__(self) -> None:
        self.dynamical_decoupling = SimpleNamespace(enable=None)
        self.twirling = SimpleNamespace(enable_gates=None, enable_measure=None)
        self.environment = SimpleNamespace(job_tags=None)
        self.max_execution_time = None


class FakeJob:
    def __init__(self, job_id: str, status: str = "DONE") -> None:
        self._job_id = job_id
        self._status = status

    def job_id(self):
        return self._job_id

    def status(self):
        return self._status

    @property
    def usage_estimation(self):
        return {"quantum_seconds": 2.0}

    def metrics(self):
        return {"usage": {"quantum_seconds": 2.0}}


class FakeService:
    def __init__(self) -> None:
        self.jobs_by_id = {}

    def usage(self):
        return {
            "usage_limit_seconds": 600,
            "usage_remaining_seconds": 600,
            "usage_consumed_seconds": 0,
            "usage_limit_reached": False,
        }

    def job(self, job_id):
        return self.jobs_by_id[job_id]


def test_submit_is_resumable_and_submits_only_one_completed_group_at_a_time(tmp_path: Path) -> None:
    service = FakeService()
    counter = {"value": 0}

    class FakeSampler:
        def __init__(self, mode, options) -> None:
            assert mode == "backend-object"
            assert options.dynamical_decoupling.enable is False
            assert options.twirling.enable_gates is False
            assert options.twirling.enable_measure is False

        def run(self, circuits, shots):
            counter["value"] += 1
            job = FakeJob(f"job-{counter['value']}")
            service.jobs_by_id[job.job_id()] = job
            return job

    bindings = SimpleNamespace(SamplerV2=FakeSampler, SamplerOptions=FakeSamplerOptions)
    circuit = runtime.PreparedCircuit(
        opaque_id="opaque",
        family="cayley",
        shots=256,
        backend_role="public-role",
        logical_qpy_sha256="aa" * 32,
        compiled_qpy_sha256="bb" * 32,
        compiled_duration_seconds=0.001,
        logical_circuit=object(),
        compiled_circuit=object(),
    )

    def group(number):
        return runtime.PreparedGroup(
            group_id=f"{number:064x}",
            backend_role="development",
            backend_role_opaque_id="public-role",
            family=f"family-{number}",
            backend_name="ibm_fez",
            physical_layout=(10, 18, 94, 124),
            shots=256,
            estimated_qpu_seconds=3.0,
            circuits=(circuit,),
            logical_qpy=b"logical",
            compiled_qpy=b"compiled",
        )

    prepared = runtime.PreparedRun(
        plan={
            "selected_backend_role": "development",
            "plan_sha256": "cc" * 32,
        },
        groups=(group(1), group(2)),
        backends={"public-role": "backend-object"},
    )
    bundle = SimpleNamespace(
        manifest={
            "resources": {
                "account_quota_seconds": 600,
                "account_reserve_seconds": 180,
                "estimated_qpu_seconds_ceiling": 420,
                "max_execution_time_seconds": 420,
            }
        },
        manifest_sha256="dd" * 32,
        analysis_lock_sha256="ee" * 32,
    )
    journal = tmp_path / "runtime_development_submission_events.ndjson"
    first = runtime.submit_prepared_run(
        prepared, bundle, service, bindings, tmp_path, submission_journal=journal
    )
    assert [event["job_id"] for event in first] == ["job-1"]
    second = runtime.submit_prepared_run(
        prepared, bundle, service, bindings, tmp_path, submission_journal=journal
    )
    assert [event["job_id"] for event in second] == ["job-2"]
    third = runtime.submit_prepared_run(
        prepared, bundle, service, bindings, tmp_path, submission_journal=journal
    )
    assert third == []
    assert counter["value"] == 2


def test_sampler_joint_counts_and_each_register_are_preserved() -> None:
    class Bits:
        def __init__(self, counts):
            self._counts = counts

        def get_counts(self):
            return self._counts

    data = SimpleNamespace(
        keys=lambda: ["heated", "decision", "final"],
        heated=Bits({"010": 3}),
        decision=Bits({"1": 3}),
        final=Bits({"001": 3}),
    )
    pub = SimpleNamespace(
        data=data,
        join_data=lambda: Bits({"0011010": 3}),
        metadata={"shots": 3},
    )
    result = runtime.extract_pub_result(pub, "opaque")
    assert result["joint_counts"] == {"0011010": 3}
    assert result["register_counts"]["heated"] == {"010": 3}
