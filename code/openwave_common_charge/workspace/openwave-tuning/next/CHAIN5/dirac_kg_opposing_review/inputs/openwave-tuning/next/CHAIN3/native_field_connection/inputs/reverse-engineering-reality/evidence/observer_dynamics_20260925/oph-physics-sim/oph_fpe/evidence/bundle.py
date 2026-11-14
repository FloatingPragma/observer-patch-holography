from __future__ import annotations

import csv
from importlib import metadata
import json
import math
import platform
from pathlib import Path
from typing import Any

import yaml

from oph_fpe.dynamics.repair import RepairEvent
from oph_fpe.evidence.cross_repo_artifacts import repository_provenance
from oph_fpe.evidence.hashes import CANONICAL_HASH_SCHEMA, stable_json_hash


def strict_jsonable(value: Any) -> Any:
    """Replace non-finite floats with None so emitted JSON stays strict-parser safe."""

    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {key: strict_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [strict_jsonable(item) for item in value]
    return value


class RunBundle:
    def __init__(self, root: Path, run_id: str):
        self.path = root / run_id
        self.path.mkdir(parents=True, exist_ok=False)
        (self.path / "controls").mkdir()
        (self.path / "plots").mkdir()

    def write_config(self, config: dict[str, Any]) -> None:
        (self.path / "config.yml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        (self.path / "seed_material.json").write_text(
            json.dumps(
                {
                    "config_hash": stable_json_hash(config),
                    "hash_schema": CANONICAL_HASH_SCHEMA,
                    "seed": config.get("seed"),
                },
                indent=2,
            ),
            encoding="utf-8",
        )

    def write_manifest(self, manifest: dict[str, Any]) -> None:
        manifest = dict(manifest)
        source_provenance = _source_provenance()
        manifest["git_commit"] = source_provenance["simulator"]["commit"]
        manifest["source_provenance"] = source_provenance
        manifest["dependency_provenance"] = _dependency_provenance()
        (self.path / "manifest.json").write_text(
            json.dumps(strict_jsonable(manifest), indent=2, default=str), encoding="utf-8"
        )
        (self.path / "git_commit.txt").write_text(str(manifest["git_commit"]) + "\n", encoding="utf-8")

    def write_json(
        self,
        name: str,
        data: Any,
        *,
        compact: bool = False,
        max_bytes: int | None = None,
    ) -> dict[str, Any]:
        dump_options = {"separators": (",", ":")} if compact else {"indent": 2}
        serialized = json.dumps(strict_jsonable(data), default=str, **dump_options)
        byte_count = len(serialized.encode("utf-8"))
        if max_bytes is not None and byte_count >= int(max_bytes):
            raise ValueError(
                f"{name} is {byte_count} bytes; hard gate requires < {int(max_bytes)} bytes"
            )
        (self.path / name).write_text(serialized, encoding="utf-8")
        return {
            "path": name,
            "byte_count": byte_count,
            "compact_json": bool(compact),
            "hard_maximum_bytes_exclusive": int(max_bytes) if max_bytes is not None else None,
            "under_hard_limit": bool(max_bytes is None or byte_count < int(max_bytes)),
        }

    def write_jsonl(self, name: str, rows: list[dict[str, Any]]) -> None:
        with (self.path / name).open("w", encoding="utf-8") as handle:
            for row in rows:
                handle.write(json.dumps(strict_jsonable(row), default=str) + "\n")

    def write_mismatch_trace(self, rows: list[dict[str, Any]]) -> None:
        _write_csv(self.path / "mismatch_trace.csv", rows)

    def write_repair_events(self, events: list[RepairEvent]) -> None:
        _write_csv(
            self.path / "repair_events.csv",
            [
                {
                    "cycle": event.cycle,
                    "node": event.node,
                    "beta": event.beta,
                    "phi_before": event.phi_before,
                    "phi_after": event.phi_after,
                    "delta_phi": event.delta_phi,
                    "accepted": event.accepted,
                    "reason": event.reason,
                }
                for event in events
            ],
        )


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _source_provenance() -> dict[str, Any]:
    simulator_root = Path(__file__).resolve().parents[2]
    research_root = simulator_root.parent / "reverse-engineering-reality"
    result: dict[str, Any] = {"simulator": repository_provenance(simulator_root)}
    if research_root.is_dir():
        research = repository_provenance(research_root)
        release_file = research_root / "paper" / "release_info.tex"
        research["paper_release_file_sha256"] = (
            stable_json_hash({"text": release_file.read_text(encoding="utf-8")})
            if release_file.is_file()
            else None
        )
        result["research"] = research
    return result


def _dependency_provenance() -> dict[str, Any]:
    packages = {}
    for name in ("oph-fpe", "numpy", "scipy", "networkx", "jsonschema", "pyyaml"):
        try:
            packages[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            packages[name] = None
    pyproject = Path(__file__).resolve().parents[2] / "pyproject.toml"
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": packages,
        "pyproject_sha256": (
            stable_json_hash({"text": pyproject.read_text(encoding="utf-8")}) if pyproject.is_file() else None
        ),
    }
