#!/usr/bin/env python3
"""Helpers for the hadron production-execution ingestion path."""

from __future__ import annotations

import copy
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


LAMBDA_MSBAR_3_GEV = 0.3344017072821104
SOURCE_ID_ALIASES = {
    "s0": "src0",
    "src0": "src0",
    "s1": "src1",
    "src1": "src1",
}


def timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_json(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def dump_json(path: str | Path, payload: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def normalize_source_id(src_id: str) -> str:
    if src_id not in SOURCE_ID_ALIASES:
        raise ValueError(f"unknown source id {src_id!r}")
    return SOURCE_ID_ALIASES[src_id]


def _finite_float(value: Any, *, field_name: str) -> float:
    try:
        result = float(value)
    except Exception as exc:  # pragma: no cover - defensive
        raise ValueError(f"{field_name} must be numeric") from exc
    if not math.isfinite(result):
        raise ValueError(f"{field_name} must be finite")
    return result


def _float_array(values: Any, *, expected_length: int, field_name: str) -> list[float]:
    if not isinstance(values, list):
        raise ValueError(f"{field_name} must be a JSON list")
    if len(values) != expected_length:
        raise ValueError(
            f"{field_name} has length {len(values)}, expected {expected_length}"
        )
    return [_finite_float(value, field_name=field_name) for value in values]


def _receipt_schedule_map(receipt: dict[str, Any]) -> dict[str, dict[str, Any]]:
    schedule = {}
    for entry in (receipt.get("execution_contract") or {}).get("ensemble_schedule", []):
        schedule[str(entry["ensemble_id"])] = entry
    return schedule


def fill_runtime_receipt(
    receipt: dict[str, Any],
    *,
    n_therm: int | None,
    n_sep: int | None,
    schedule_provenance: str | None = None,
) -> dict[str, Any]:
    """Fill the receipt with externally chosen schedule integers."""
    out = copy.deepcopy(receipt)
    if n_therm is None or n_sep is None:
        scalars = out.get("required_schedule_scalars") or {}
        if scalars.get("N_therm") is not None and scalars.get("N_sep") is not None:
            n_therm = int(scalars["N_therm"])
            n_sep = int(scalars["N_sep"])
        else:
            return out
    if n_therm is None or n_sep is None:
        return out
    if n_therm < 0 or n_sep <= 0:
        raise ValueError("N_therm must be >= 0 and N_sep must be > 0")
    out["generated_utc"] = timestamp()
    out["status"] = "receipt_filled_waiting_backend_dump"
    out["required_schedule_scalars"] = {
        "N_therm": int(n_therm),
        "N_sep": int(n_sep),
    }
    existing_provenance = out.get("execution_input_provenance") or {}
    out["execution_input_provenance"] = {
        "schedule_scalars_source": (
            schedule_provenance
            or existing_provenance.get("schedule_scalars_source")
            or "external_runtime_input"
        ),
        "trajectory_stop_derivation": "execution_contract.stop_time_formula",
    }
    for sched in (out.get("execution_contract") or {}).get("ensemble_schedule", []):
        stops: dict[str, int] = {}
        formulas: dict[str, str] = {}
        for cfg_index, cfg_id in enumerate(sched.get("cfg_ids", [])):
            stop = int(n_therm) + cfg_index * int(n_sep)
            stops[str(cfg_id)] = stop
            formulas[str(cfg_id)] = f"N_therm + {cfg_index}*N_sep"
        sched["trajectory_stop_by_cfg"] = stops
        sched["trajectory_stop_by_cfg_formula"] = formulas
    return out


def _iter_backend_ensemble_items(backend_input: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    ensembles = backend_input.get("ensembles")
    if isinstance(ensembles, dict):
        return [(str(key), value) for key, value in ensembles.items()]
    if isinstance(ensembles, list):
        return [(str(item["ensemble_id"]), item) for item in ensembles]
    raise ValueError("backend input must provide an 'ensembles' mapping or list")


def _iter_backend_cfg_items(ensemble_entry: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    cfgs = ensemble_entry.get("cfgs")
    if isinstance(cfgs, dict):
        return [(str(key), value) for key, value in cfgs.items()]
    if isinstance(cfgs, list):
        return [(str(item["cfg_id"]), item) for item in cfgs]
    raise ValueError("backend ensemble entry must provide a 'cfgs' mapping or list")


def _iter_backend_source_items(cfg_entry: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    sources = cfg_entry.get("sources")
    if isinstance(sources, dict):
        return [(normalize_source_id(str(key)), value) for key, value in sources.items()]
    if isinstance(sources, list):
        return [
            (normalize_source_id(str(item.get("src_id") or item.get("source_id"))), item)
            for item in sources
        ]
    raise ValueError("backend cfg entry must provide a 'sources' mapping or list")


def _normalized_backend_tree(backend_input: dict[str, Any]) -> dict[str, Any]:
    tree: dict[str, Any] = {}
    for ensemble_id, ensemble_entry in _iter_backend_ensemble_items(backend_input):
        cfg_tree: dict[str, Any] = {}
        for cfg_id, cfg_entry in _iter_backend_cfg_items(ensemble_entry):
            source_tree: dict[str, Any] = {}
            for source_id, source_entry in _iter_backend_source_items(cfg_entry):
                source_tree[source_id] = source_entry
            cfg_tree[cfg_id] = {
                "trajectory_stop": cfg_entry.get("trajectory_stop"),
                "sources": source_tree,
            }
        tree[ensemble_id] = {
            "ensemble_id": ensemble_id,
            "cfgs": cfg_tree,
        }
    return tree


def _coord4(values: Any, *, field_name: str) -> list[int]:
    seq = list(values or [])
    if len(seq) != 4:
        raise ValueError(f"{field_name} must contain exactly 4 coordinates")
    out: list[int] = []
    for idx, value in enumerate(seq):
        try:
            out.append(int(value))
        except Exception as exc:  # pragma: no cover - defensive
            raise ValueError(f"{field_name}[{idx}] is not an integer") from exc
    return out


def build_backend_manifest(
    receipt: dict[str, Any],
    payload: dict[str, Any],
    backend_input: dict[str, Any],
    *,
    backend_input_path: str | None = None,
) -> dict[str, Any]:
    tree = _normalized_backend_tree(backend_input)
    schedule_map = _receipt_schedule_map(receipt)
    manifest_ensembles = []
    payload_map = {
        str(entry["ensemble_id"]): entry
        for entry in payload.get("ensemble_payloads", [])
    }
    for ensemble_id, entry in tree.items():
        sched = schedule_map.get(ensemble_id)
        payload_entry = payload_map.get(ensemble_id, {})
        if sched is None:
            raise ValueError(f"ensemble {ensemble_id!r} not present in receipt schedule")
        cfg_ids = [str(cfg_id) for cfg_id in sched.get("cfg_ids", [])]
        manifest_cfgs = []
        for cfg_id in cfg_ids:
            cfg_entry = entry.get("cfgs", {}).get(cfg_id)
            if cfg_entry is None:
                raise ValueError(f"backend input missing cfg {cfg_id!r}")
            manifest_sources = []
            src_descriptors = (payload_entry.get("source_descriptors_by_cfg") or {}).get(cfg_id, [])
            for src_desc in src_descriptors:
                norm_src = normalize_source_id(str(src_desc["src_id"]))
                source_entry = cfg_entry["sources"].get(norm_src)
                if source_entry is None:
                    raise ValueError(f"backend input missing source {norm_src!r} for {cfg_id}")
                manifest_sources.append(
                    {
                        "source_id": norm_src,
                        "coords": _coord4(
                            src_desc.get("coords") or src_desc.get("coord"),
                            field_name=f"{ensemble_id}.{cfg_id}.{norm_src}.coord",
                        ),
                        "channels_present": sorted(
                            [
                                channel
                                for channel in ("pi_iso", "N_iso_direct", "N_iso_exchange")
                                if isinstance(source_entry.get(channel), list)
                            ]
                        ),
                    }
                )
            manifest_cfgs.append(
                {
                    "cfg_id": cfg_id,
                    "trajectory_stop": cfg_entry.get("trajectory_stop"),
                    "sources": manifest_sources,
                }
            )
        manifest_ensembles.append(
            {
                "ensemble_id": ensemble_id,
                "family_index": payload_entry.get("family_index"),
                "beta": payload_entry.get("beta"),
                "L": payload_entry.get("L"),
                "T": payload_entry.get("T"),
                "aLambda_msbar3": payload_entry.get("aLambda_msbar3"),
                "am_l": payload_entry.get("am_l"),
                "am_s": payload_entry.get("am_s"),
                "cfgs": manifest_cfgs,
            }
        )
    raw_export_provenance = copy.deepcopy(backend_input.get("raw_export_provenance") or {})
    out = {
        "artifact": "oph_hadron_production_backend_manifest",
        "generated_utc": timestamp(),
        "receipt_artifact": receipt.get("artifact"),
        "backend_input_artifact": backend_input.get("artifact"),
        "backend_input_path": backend_input_path,
        "schedule_scalars": copy.deepcopy(receipt.get("required_schedule_scalars")),
        "execution_input_provenance": copy.deepcopy(receipt.get("execution_input_provenance")),
        "ensemble_tasks": manifest_ensembles,
    }
    if raw_export_provenance:
        out["raw_export_provenance"] = raw_export_provenance
        backend_meta = raw_export_provenance.get("backend") or {}
        if backend_meta:
            out["backend_name"] = backend_meta.get("name")
            out["backend_version"] = backend_meta.get("version")
            out["backend_run_id"] = backend_meta.get("run_id")
        profile_id = backend_input.get("profile_id") or raw_export_provenance.get("profile_id")
        if profile_id is not None:
            out["backend_profile_id"] = profile_id
    return out


def build_production_dump(
    receipt: dict[str, Any],
    payload: dict[str, Any],
    backend_input: dict[str, Any],
) -> dict[str, Any]:
    """Normalize backend-produced correlator arrays into the frozen dump schema."""
    schedule_map = _receipt_schedule_map(receipt)
    payload_map = {
        str(entry["ensemble_id"]): entry
        for entry in payload.get("ensemble_payloads", [])
    }
    backend_tree = _normalized_backend_tree(backend_input)
    dump: dict[str, Any] = {
        "artifact": "backend_correlator_dump.production",
        "generated_utc": timestamp(),
        "production_execution": True,
        "dry_run": False,
        "surrogate_execution": False,
        "tiny_geometry_pilot": False,
        "receipt_artifact": receipt.get("artifact"),
        "manifest_artifact": "oph_hadron_production_backend_manifest",
        "ensembles": {},
        "notes": [
            "Normalized from backend-produced correlator arrays.",
            "All pi_iso and N_iso direct/exchange arrays are validated against the frozen cfg/source/t-extent contract.",
        ],
    }
    for ensemble_id, sched in schedule_map.items():
        payload_entry = payload_map.get(ensemble_id)
        backend_entry = backend_tree.get(ensemble_id)
        if payload_entry is None:
            raise ValueError(f"payload is missing ensemble {ensemble_id!r}")
        if backend_entry is None:
            raise ValueError(f"backend input is missing ensemble {ensemble_id!r}")
        t_extent = int(payload_entry["T"])
        schedule_t = sched.get("t_extent")
        if schedule_t is not None and int(schedule_t) != t_extent:
            raise ValueError(
                f"t_extent mismatch for {ensemble_id}: receipt={schedule_t}, payload={t_extent}"
            )
        source_descriptors_by_cfg = payload_entry.get("source_descriptors_by_cfg") or {}
        ensemble_dump = {
            "ensemble_id": ensemble_id,
            "family_index": payload_entry.get("family_index"),
            "beta": payload_entry.get("beta"),
            "L": payload_entry.get("L"),
            "T": payload_entry.get("T"),
            "aLambda_msbar3": payload_entry.get("aLambda_msbar3"),
            "am_l": payload_entry.get("am_l"),
            "am_s": payload_entry.get("am_s"),
            "cfgs": {},
        }
        for cfg_id in sched.get("cfg_ids", []):
            cfg_id = str(cfg_id)
            backend_cfg = (backend_entry.get("cfgs") or {}).get(cfg_id)
            if backend_cfg is None:
                raise ValueError(f"backend input is missing cfg {cfg_id!r}")
            receipt_stop = (sched.get("trajectory_stop_by_cfg") or {}).get(cfg_id)
            input_stop = backend_cfg.get("trajectory_stop")
            if input_stop is not None:
                input_stop = int(input_stop)
            trajectory_stop = receipt_stop if receipt_stop is not None else input_stop
            if receipt_stop is not None and input_stop is not None and int(receipt_stop) != int(input_stop):
                raise ValueError(
                    f"trajectory stop mismatch for {cfg_id}: receipt={receipt_stop}, backend={input_stop}"
                )
            source_dump = {}
            for src_desc in source_descriptors_by_cfg.get(cfg_id, []):
                norm_src = normalize_source_id(str(src_desc["src_id"]))
                backend_source = (backend_cfg.get("sources") or {}).get(norm_src)
                if backend_source is None:
                    raise ValueError(f"backend input missing source {norm_src!r} for cfg {cfg_id}")
                coords = _coord4(
                    src_desc.get("coords") or src_desc.get("coord"),
                    field_name=f"{ensemble_id}.{cfg_id}.{norm_src}.coord",
                )
                source_dump[norm_src] = {
                    "coord": coords,
                    "expected_t_extent": t_extent,
                    "pi_iso": _float_array(
                        backend_source.get("pi_iso"),
                        expected_length=t_extent,
                        field_name=f"{ensemble_id}.{cfg_id}.{norm_src}.pi_iso",
                    ),
                    "N_iso_direct": _float_array(
                        backend_source.get("N_iso_direct"),
                        expected_length=t_extent,
                        field_name=f"{ensemble_id}.{cfg_id}.{norm_src}.N_iso_direct",
                    ),
                    "N_iso_exchange": _float_array(
                        backend_source.get("N_iso_exchange"),
                        expected_length=t_extent,
                        field_name=f"{ensemble_id}.{cfg_id}.{norm_src}.N_iso_exchange",
                    ),
                }
            ensemble_dump["cfgs"][cfg_id] = {
                "trajectory_stop": trajectory_stop,
                "sources": source_dump,
            }
        dump["ensembles"][ensemble_id] = ensemble_dump
    return dump


def ingest_dump_into_payload(
    payload: dict[str, Any],
    dump: dict[str, Any],
    receipt: dict[str, Any] | None = None,
) -> dict[str, Any]:
    out = copy.deepcopy(payload)
    receipt = receipt or {}
    dump_ensembles = dump.get("ensembles") or {}
    schedule_map = _receipt_schedule_map(receipt)
    out["generated_utc"] = timestamp()
    out["status"] = "production_backend_dump_ingested"
    out["proof_status"] = "production_backend_dump_ingested"
    out["payload_realization_status"] = "production_backend_dump_ingested"
    out["cfg_support_realization_status"] = "executed_fixed_schedule_rhmc_hmc_on_seeded_2p1_family"
    out["smallest_constructive_missing_object"] = "oph_hadron_stable_channel_sequence_evaluator"
    for ensemble in out.get("ensemble_payloads", []):
        ensemble_id = str(ensemble["ensemble_id"])
        dump_ensemble = dump_ensembles.get(ensemble_id)
        if dump_ensemble is None:
            raise ValueError(f"production dump missing ensemble {ensemble_id!r}")
        pi_cfg_source = []
        n_dir_cfg_source = []
        n_ex_cfg_source = []
        n_cfg_source = []
        for cfg_id in ensemble.get("cfg_ids", []):
            cfg_id = str(cfg_id)
            dump_cfg = (dump_ensemble.get("cfgs") or {}).get(cfg_id)
            if dump_cfg is None:
                raise ValueError(f"production dump missing cfg {cfg_id!r}")
            pi_src = []
            n_dir_src = []
            n_ex_src = []
            n_src = []
            for src_desc in (ensemble.get("source_descriptors_by_cfg") or {}).get(cfg_id, []):
                norm_src = normalize_source_id(str(src_desc["src_id"]))
                source_entry = (dump_cfg.get("sources") or {}).get(norm_src)
                if source_entry is None:
                    raise ValueError(f"production dump missing source {norm_src!r} for {cfg_id}")
                pi_t = list(source_entry["pi_iso"])
                n_dir_t = list(source_entry["N_iso_direct"])
                n_ex_t = list(source_entry["N_iso_exchange"])
                pi_src.append(pi_t)
                n_dir_src.append(n_dir_t)
                n_ex_src.append(n_ex_t)
                n_src.append([direct - exchange for direct, exchange in zip(n_dir_t, n_ex_t)])
            pi_cfg_source.append(pi_src)
            n_dir_cfg_source.append(n_dir_src)
            n_ex_cfg_source.append(n_ex_src)
            n_cfg_source.append(n_src)
        ensemble["pi_iso"]["cfg_source_corr_t"] = pi_cfg_source
        ensemble["N_iso"]["cfg_source_corr_direct_t"] = n_dir_cfg_source
        ensemble["N_iso"]["cfg_source_corr_exchange_t"] = n_ex_cfg_source
        ensemble["N_iso"]["cfg_source_corr_t"] = n_cfg_source
    schedule = out.get("support_realization_schedule") or {}
    if schedule:
        schedule["status"] = "executed_fixed_schedule_rhmc_hmc_on_seeded_2p1_family"
        schedule["required_schedule_scalars"] = copy.deepcopy(receipt.get("required_schedule_scalars"))
        for sched in schedule.get("ensemble_schedule", []):
            receipt_sched = schedule_map.get(str(sched["ensemble_id"]), {})
            sched["trajectory_stop_by_cfg"] = copy.deepcopy(receipt_sched.get("trajectory_stop_by_cfg") or {})
            sched["trajectory_stop_by_cfg_formula"] = copy.deepcopy(
                receipt_sched.get("trajectory_stop_by_cfg_formula") or {}
            )
    return out


def _cfg_source_average(cfg_source_corr_t: list[list[list[float]]]) -> list[list[float]]:
    cfg_averages = []
    for source_arrays in cfg_source_corr_t:
        n_src = max(len(source_arrays), 1)
        t_extent = len(source_arrays[0]) if source_arrays else 0
        cfg_averages.append(
            [
                sum(source_arrays[src_idx][t_idx] for src_idx in range(n_src)) / float(n_src)
                for t_idx in range(t_extent)
            ]
        )
    return cfg_averages


def _mean_over_cfg(cfg_averages: list[list[float]]) -> list[float]:
    n_cfg = max(len(cfg_averages), 1)
    t_extent = len(cfg_averages[0]) if cfg_averages else 0
    return [
        sum(cfg_averages[cfg_idx][t_idx] for cfg_idx in range(n_cfg)) / float(n_cfg)
        for t_idx in range(t_extent)
    ]


def _jackknife_samples(cfg_averages: list[list[float]]) -> list[list[float]]:
    n_cfg = len(cfg_averages)
    if n_cfg <= 1:
        return [list(cfg_averages[0])] if cfg_averages else []
    t_extent = len(cfg_averages[0]) if cfg_averages else 0
    samples: list[list[float]] = []
    for omit_cfg in range(n_cfg):
        keep = [cfg_averages[idx] for idx in range(n_cfg) if idx != omit_cfg]
        denom = float(len(keep))
        samples.append(
            [
                sum(keep[cfg_idx][t_idx] for cfg_idx in range(len(keep))) / denom
                for t_idx in range(t_extent)
            ]
        )
    return samples


def _jackknife_stderr(samples: list[list[float]]) -> list[float]:
    n_samples = len(samples)
    if n_samples == 0:
        return []
    t_extent = len(samples[0])
    means = [
        sum(sample[t_idx] for sample in samples) / float(n_samples)
        for t_idx in range(t_extent)
    ]
    prefactor = float(n_samples - 1) / float(n_samples) if n_samples > 0 else 0.0
    return [
        math.sqrt(
            max(
                prefactor
                * sum((sample[t_idx] - means[t_idx]) ** 2 for sample in samples),
                0.0,
            )
        )
        for t_idx in range(t_extent)
    ]


def _safe_log_ratio(num: float, den: float, *, absolute: bool) -> float | None:
    if absolute:
        num = abs(num)
        den = abs(den)
    if num <= 0.0 or den <= 0.0:
        return None
    return math.log(num / den)


def _am_eff(corr_t: list[float], *, absolute: bool) -> list[float | None]:
    return [
        _safe_log_ratio(corr_t[t_idx], corr_t[t_idx + 1], absolute=absolute)
        for t_idx in range(max(len(corr_t) - 1, 0))
    ]


def _log_convexity(corr_t: list[float], *, absolute: bool) -> list[float | None]:
    residuals: list[float | None] = [None] * len(corr_t)
    for t_idx in range(1, max(len(corr_t) - 1, 1)):
        center = abs(corr_t[t_idx]) if absolute else corr_t[t_idx]
        left = abs(corr_t[t_idx - 1]) if absolute else corr_t[t_idx - 1]
        right = abs(corr_t[t_idx + 1]) if absolute else corr_t[t_idx + 1]
        residuals[t_idx] = center * center - left * right
    return residuals


def _tail_drop(am_eff_t: list[float | None]) -> list[float | None]:
    return [
        None
        if am_eff_t[t_idx] is None or am_eff_t[t_idx + 1] is None
        else am_eff_t[t_idx] - am_eff_t[t_idx + 1]
        for t_idx in range(max(len(am_eff_t) - 1, 0))
    ]


def _mirror_tail_indicator(am_eff_t: list[float | None], t_extent: int) -> list[float | None]:
    out: list[float | None] = []
    for t_idx, value in enumerate(am_eff_t):
        if value is None:
            out.append(None)
            continue
        out.append(math.exp(-value * (t_extent - 2 * t_idx)))
    return out


def _signs(corr_t: list[float]) -> list[int]:
    out = []
    for value in corr_t:
        if value > 0:
            out.append(1)
        elif value < 0:
            out.append(-1)
        else:
            out.append(0)
    return out


def _jk_scalar_stderr(samples: list[float]) -> float | None:
    if not samples:
        return None
    n_samples = len(samples)
    mean = sum(samples) / float(n_samples)
    return math.sqrt(
        max(
            (float(n_samples - 1) / float(n_samples))
            * sum((value - mean) ** 2 for value in samples),
            0.0,
        )
    )


def _safe_json_list(values: list[Any]) -> list[Any]:
    out = []
    for value in values:
        if value is None:
            out.append(None)
        elif isinstance(value, bool):
            out.append(bool(value))
        else:
            out.append(float(value))
    return out


def _safe_json_matrix(values: list[list[Any]]) -> list[list[Any]]:
    return [_safe_json_list(row) for row in values]


def _forward_window_candidates(t_extent: int) -> list[int]:
    return list(range(1, max(t_extent // 2 - 1, 1)))


def _contiguous_runs(indices: list[int]) -> list[list[int]]:
    if not indices:
        return []
    runs = [[indices[0]]]
    for idx in indices[1:]:
        if idx == runs[-1][-1] + 1:
            runs[-1].append(idx)
        else:
            runs.append([idx])
    return runs


def _weighted_mean(values: list[float], errors: list[float | None]) -> float | None:
    weights = []
    for value, error in zip(values, errors):
        if value is None:
            weights.append(0.0)
            continue
        sigma = abs(error) if error not in (None, 0.0) else 1.0
        weights.append(1.0 / max(sigma * sigma, 1.0e-12))
    total_weight = sum(weights)
    if total_weight <= 0.0:
        return None
    return sum(value * weight for value, weight in zip(values, weights) if value is not None) / total_weight


def _channel_measurement(
    cfg_source_corr_t: list[list[list[float]]],
    *,
    absolute_effective_mass: bool,
    require_sign_stability: bool = False,
    corr_direct_t: list[list[list[float]]] | None = None,
    corr_exchange_t: list[list[list[float]]] | None = None,
) -> dict[str, Any]:
    cfg_avg = _cfg_source_average(cfg_source_corr_t)
    corr_t = _mean_over_cfg(cfg_avg)
    corr_t_jk = _jackknife_samples(cfg_avg)
    corr_t_stderr = _jackknife_stderr(corr_t_jk)

    am_eff_t = _am_eff(corr_t, absolute=absolute_effective_mass)
    am_eff_t_jk = [_am_eff(sample, absolute=absolute_effective_mass) for sample in corr_t_jk]
    am_eff_t_stderr = _jackknife_stderr(
        [[0.0 if value is None else value for value in sample] for sample in am_eff_t_jk]
    ) if am_eff_t_jk else []

    log_conv_t = _log_convexity(corr_t, absolute=absolute_effective_mass)
    log_conv_t_jk = [_log_convexity(sample, absolute=absolute_effective_mass) for sample in corr_t_jk]
    log_conv_t_stderr = _jackknife_stderr(
        [[0.0 if value is None else value for value in sample] for sample in log_conv_t_jk]
    ) if log_conv_t_jk else []

    tail_drop_t = _tail_drop(am_eff_t)
    tail_drop_t_jk = [_tail_drop(sample) for sample in am_eff_t_jk]
    tail_drop_t_stderr = _jackknife_stderr(
        [[0.0 if value is None else value for value in sample] for sample in tail_drop_t_jk]
    ) if tail_drop_t_jk and tail_drop_t_jk[0] else []

    mirror_t = _mirror_tail_indicator(am_eff_t, len(corr_t))
    mirror_t_jk = [_mirror_tail_indicator(sample, len(corr_t)) for sample in am_eff_t_jk]
    mirror_t_stderr = _jackknife_stderr(
        [[0.0 if value is None else value for value in sample] for sample in mirror_t_jk]
    ) if mirror_t_jk else []

    corr_sign_t = _signs(corr_t)
    sign_stable_t = [
        corr_sign_t[t_idx] * corr_sign_t[t_idx + 1] > 0
        for t_idx in range(max(len(corr_sign_t) - 1, 0))
    ]

    direct_minus_exchange_residual_t: list[float] = []
    direct_minus_exchange_residual_t_jk: list[list[float]] = []
    direct_minus_exchange_consistent_t: list[bool] = []
    if corr_direct_t is not None and corr_exchange_t is not None:
        cfg_direct_avg = _cfg_source_average(corr_direct_t)
        cfg_exchange_avg = _cfg_source_average(corr_exchange_t)
        corr_direct = _mean_over_cfg(cfg_direct_avg)
        corr_exchange = _mean_over_cfg(cfg_exchange_avg)
        corr_direct_jk = _jackknife_samples(cfg_direct_avg)
        corr_exchange_jk = _jackknife_samples(cfg_exchange_avg)
        direct_minus_exchange_residual_t = [
            corr_t[t_idx] - (corr_direct[t_idx] - corr_exchange[t_idx])
            for t_idx in range(len(corr_t))
        ]
        direct_minus_exchange_residual_t_jk = [
            [
                corr_t_jk[jk_idx][t_idx]
                - (corr_direct_jk[jk_idx][t_idx] - corr_exchange_jk[jk_idx][t_idx])
                for t_idx in range(len(corr_t))
            ]
            for jk_idx in range(len(corr_t_jk))
        ]
        residual_stderr = _jackknife_stderr(direct_minus_exchange_residual_t_jk)
        direct_minus_exchange_consistent_t = [
            abs(value) <= max(residual_stderr[t_idx] if residual_stderr else 0.0, 1.0e-12)
            for t_idx, value in enumerate(direct_minus_exchange_residual_t)
        ]
    else:
        corr_direct = []
        corr_exchange = []
        corr_direct_jk = []
        corr_exchange_jk = []
        residual_stderr = []

    candidate_t = _forward_window_candidates(len(corr_t))
    log_convexity_certified_t = [False] * max(len(corr_t) - 1, 0)
    plateau_noise_floor_t = [False] * max(len(corr_t) - 1, 0)
    plateau_flat_t = [False] * max(len(corr_t) - 1, 0)
    monotone_tail_t = [False] * max(len(corr_t) - 1, 0)
    mirror_suppressed_t = [False] * max(len(corr_t) - 1, 0)
    forward_certificate_t = [False] * max(len(corr_t) - 1, 0)
    for t_idx in candidate_t:
        eff = am_eff_t[t_idx] if t_idx < len(am_eff_t) else None
        if eff is None:
            continue
        cvx = log_conv_t[t_idx]
        cvx_sigma = log_conv_t_stderr[t_idx] if t_idx < len(log_conv_t_stderr) else 0.0
        log_convexity_certified_t[t_idx] = cvx is not None and cvx >= -max(3.0 * cvx_sigma, 1.0e-12)
        sigma_here = corr_t_stderr[t_idx] if t_idx < len(corr_t_stderr) else 0.0
        plateau_noise_floor_t[t_idx] = abs(corr_t[t_idx]) >= 0.5 * max(sigma_here, 1.0e-12)
        if t_idx < len(tail_drop_t):
            tail = tail_drop_t[t_idx]
            tail_sigma = tail_drop_t_stderr[t_idx] if t_idx < len(tail_drop_t_stderr) else 0.0
            if tail is not None:
                monotone_tail_t[t_idx] = tail >= -max(2.0 * tail_sigma, 0.02 * abs(eff), 1.0e-12)
                plateau_flat_t[t_idx] = abs(tail) <= max(2.0 * tail_sigma, 0.05 * abs(eff), 1.0e-12)
        mirror = mirror_t[t_idx] if t_idx < len(mirror_t) else None
        mirror_suppressed_t[t_idx] = mirror is not None and mirror <= 0.25
        n_checks_ok = True
        if require_sign_stability:
            n_checks_ok = (
                t_idx < len(sign_stable_t)
                and sign_stable_t[t_idx]
                and t_idx < len(direct_minus_exchange_consistent_t)
                and direct_minus_exchange_consistent_t[t_idx]
            )
        forward_certificate_t[t_idx] = (
            log_convexity_certified_t[t_idx]
            and plateau_noise_floor_t[t_idx]
            and plateau_flat_t[t_idx]
            and monotone_tail_t[t_idx]
            and mirror_suppressed_t[t_idx]
            and n_checks_ok
        )

    certified_indices = [t_idx for t_idx in candidate_t if forward_certificate_t[t_idx]]
    forward_window_runs = _contiguous_runs(certified_indices)
    selected_forward_window = []
    if forward_window_runs:
        selected_forward_window = max(forward_window_runs, key=lambda run: (len(run), run[-1]))
    selected_window_errors = [
        am_eff_t_stderr[t_idx] if t_idx < len(am_eff_t_stderr) else None
        for t_idx in selected_forward_window
    ]
    selected_window_values = [
        am_eff_t[t_idx]
        for t_idx in selected_forward_window
        if t_idx < len(am_eff_t) and am_eff_t[t_idx] is not None
    ]
    am_ground_candidate = None
    if selected_forward_window and len(selected_window_values) == len(selected_forward_window):
        am_ground_candidate = _weighted_mean(selected_window_values, selected_window_errors)

    ground_samples: list[float] = []
    if selected_forward_window:
        for jk_sample in am_eff_t_jk:
            sample_values = []
            for t_idx in selected_forward_window:
                if t_idx >= len(jk_sample) or jk_sample[t_idx] is None:
                    sample_values = []
                    break
                sample_values.append(float(jk_sample[t_idx]))
            if sample_values:
                weighted = _weighted_mean(sample_values, selected_window_errors)
                if weighted is not None:
                    ground_samples.append(weighted)
    am_ground_stat_err = _jk_scalar_stderr(ground_samples)

    return {
        "corr_t": corr_t,
        "corr_t_jk": corr_t_jk,
        "corr_t_stderr": corr_t_stderr,
        "am_eff_t": am_eff_t,
        "am_eff_t_jk": am_eff_t_jk,
        "am_eff_t_stderr": am_eff_t_stderr,
        "log_convexity_residual_t": log_conv_t,
        "log_convexity_residual_t_jk": log_conv_t_jk,
        "log_convexity_residual_t_stderr": log_conv_t_stderr,
        "log_convexity_certified_t": log_convexity_certified_t,
        "tail_drop_t": tail_drop_t,
        "tail_drop_t_jk": tail_drop_t_jk,
        "tail_drop_t_stderr": tail_drop_t_stderr,
        "tail_drop_lower_t": [0.0 if value is None else value for value in tail_drop_t],
        "mirror_tail_indicator_t": mirror_t,
        "mirror_tail_indicator_t_jk": mirror_t_jk,
        "mirror_tail_indicator_t_stderr": mirror_t_stderr,
        "mirror_to_noise_t": [
            None
            if t_idx >= len(corr_t_stderr) or corr_t_stderr[t_idx] == 0.0 or mirror_t[t_idx] is None
            else mirror_t[t_idx] / corr_t_stderr[t_idx]
            for t_idx in range(len(mirror_t))
        ],
        "mirror_to_drift_t": [
            None
            if t_idx >= len(tail_drop_t) or tail_drop_t[t_idx] in (None, 0.0) or mirror_t[t_idx] is None
            else mirror_t[t_idx] / abs(tail_drop_t[t_idx])
            for t_idx in range(len(mirror_t))
        ],
        "mirror_suppressed_t": mirror_suppressed_t,
        "plateau_noise_floor_t": plateau_noise_floor_t,
        "plateau_flat_t": plateau_flat_t,
        "monotone_tail_t": monotone_tail_t,
        "forward_certificate_t": forward_certificate_t,
        "forward_window_t": candidate_t,
        "forward_window_runs": [
            {"start": run[0], "stop": run[-1], "cardinality": len(run)}
            for run in forward_window_runs
        ],
        "selected_forward_window": selected_forward_window,
        "selected_forward_window_cardinality": len(selected_forward_window),
        "forward_window_limit_exists": bool(selected_forward_window) and am_ground_candidate is not None,
        "am_ground_candidate": am_ground_candidate,
        "am_ground_stat_err": am_ground_stat_err,
        "published_statistical_error": am_ground_stat_err,
        "corr_sign_t": corr_sign_t,
        "corr_direct_t": corr_direct,
        "corr_direct_t_jk": corr_direct_jk,
        "corr_direct_t_stderr": _jackknife_stderr(corr_direct_jk) if corr_direct_jk else [],
        "corr_exchange_t": corr_exchange,
        "corr_exchange_t_jk": corr_exchange_jk,
        "corr_exchange_t_stderr": _jackknife_stderr(corr_exchange_jk) if corr_exchange_jk else [],
        "sign_stable_t": sign_stable_t,
        "direct_minus_exchange_residual_t": direct_minus_exchange_residual_t,
        "direct_minus_exchange_residual_t_jk": direct_minus_exchange_residual_t_jk,
        "direct_minus_exchange_residual_t_stderr": residual_stderr,
        "direct_minus_exchange_consistent_t": direct_minus_exchange_consistent_t,
    }


def populate_evaluation_from_dump(
    evaluation: dict[str, Any],
    dump: dict[str, Any],
    *,
    lambda_msbar_3_gev: float = LAMBDA_MSBAR_3_GEV,
) -> dict[str, Any]:
    """Populate evaluator fields from a validated production dump."""
    out = copy.deepcopy(evaluation)
    out["generated_utc"] = timestamp()
    out["status"] = "production_measure_evaluation_complete"
    out["proof_status"] = "production_measure_evaluation_complete"
    out["smallest_constructive_missing_object"] = "StableChannelForwardWindowConvergence"
    dump_ensembles = dump.get("ensembles") or {}
    channel_entries: dict[str, list[tuple[dict[str, Any], dict[str, Any]]]] = {"pi_iso": [], "N_iso": []}
    for ensemble in out.get("ensemble_evaluations", []):
        ensemble_id = str(ensemble["ensemble_id"])
        dump_ensemble = dump_ensembles.get(ensemble_id)
        if dump_ensemble is None:
            raise ValueError(f"production dump missing ensemble {ensemble_id!r}")
        pi_raw = ensemble["pi_iso"].get("cfg_source_corr_t") or []
        n_dir_raw = ensemble["N_iso"].get("cfg_source_corr_direct_t") or []
        n_ex_raw = ensemble["N_iso"].get("cfg_source_corr_exchange_t") or []
        n_raw = ensemble["N_iso"].get("cfg_source_corr_t") or []
        if not pi_raw or not n_dir_raw or not n_ex_raw or not n_raw:
            raise ValueError(f"evaluation payload for {ensemble_id!r} is missing cfg/source arrays")
        pi_measure = _channel_measurement(
            pi_raw,
            absolute_effective_mass=False,
        )
        n_measure = _channel_measurement(
            n_raw,
            absolute_effective_mass=True,
            require_sign_stability=True,
            corr_direct_t=n_dir_raw,
            corr_exchange_t=n_ex_raw,
        )
        for field_name, value in pi_measure.items():
            if isinstance(value, list) and value and isinstance(value[0], list):
                ensemble["pi_iso"][field_name] = _safe_json_matrix(value)
            elif isinstance(value, list) and value and isinstance(value[0], dict):
                ensemble["pi_iso"][field_name] = value
            elif isinstance(value, list):
                ensemble["pi_iso"][field_name] = _safe_json_list(value)
            else:
                ensemble["pi_iso"][field_name] = value
        for field_name, value in n_measure.items():
            if isinstance(value, list) and value and isinstance(value[0], list):
                ensemble["N_iso"][field_name] = _safe_json_matrix(value)
            elif isinstance(value, list) and value and isinstance(value[0], dict):
                ensemble["N_iso"][field_name] = value
            elif isinstance(value, list):
                ensemble["N_iso"][field_name] = _safe_json_list(value)
            else:
                ensemble["N_iso"][field_name] = value

        for channel_name in ("pi_iso", "N_iso"):
            channel = ensemble[channel_name]
            a_lambda = float(ensemble["aLambda_msbar3"])
            am_ground = channel.get("am_ground_candidate")
            if am_ground is not None and math.isfinite(float(am_ground)):
                ratio = float(am_ground) / a_lambda
                channel["ratio_to_lambda_msbar3_candidate"] = ratio
                channel["mass_gev_candidate"] = ratio * float(lambda_msbar_3_gev)
            else:
                channel["ratio_to_lambda_msbar3_candidate"] = None
                channel["mass_gev_candidate"] = None
            channel["am_ground_candidate_err"] = channel.get("published_statistical_error")
            channel["am_ground_sys_err"] = None
            channel["published_systematics"] = {
                "status": "pending",
                "sigma_sys": None,
                "delta_cont": None,
                "delta_vol": None,
                "delta_chi": None,
            }
            channel["convergence_status"] = (
                "forward_window_candidate_complete"
                if channel.get("forward_window_limit_exists")
                else "forward_window_missing"
            )
            channel["sequence_status"] = channel["convergence_status"]
            channel_entries[channel_name].append((ensemble, channel))

    for channel_name, entries in channel_entries.items():
        valid = [
            (
                float(channel["ratio_to_lambda_msbar3_candidate"]),
                float(ensemble["aLambda_msbar3"]),
                float(((dump_ensembles.get(str(ensemble["ensemble_id"])) or {}).get("am_l")) or 0.0)
                / float(ensemble["aLambda_msbar3"]),
                float(((dump_ensembles.get(str(ensemble["ensemble_id"])) or {}).get("am_s")) or 0.0)
                / float(ensemble["aLambda_msbar3"]),
            )
            for ensemble, channel in entries
            if channel.get("ratio_to_lambda_msbar3_candidate") is not None
        ]
        if valid:
            x_vals = [a_lambda * a_lambda for _, a_lambda, _, _ in valid]
            y_vals = [ratio for ratio, _, _, _ in valid]
            if len(valid) >= 2:
                x_mean = sum(x_vals) / float(len(x_vals))
                y_mean = sum(y_vals) / float(len(y_vals))
                denom = sum((x - x_mean) ** 2 for x in x_vals)
                slope = 0.0 if denom == 0.0 else sum((x - x_mean) * (y - y_mean) for x, y in zip(x_vals, y_vals)) / denom
                intercept = y_mean - slope * x_mean
            else:
                intercept = y_vals[0]
            mean_r_l = sum(r_l for _, _, r_l, _ in valid) / float(len(valid))
            mean_r_s = sum(r_s for _, _, _, r_s in valid) / float(len(valid))
        else:
            intercept = None
            mean_r_l = None
            mean_r_s = None
        for ensemble, channel in entries:
            dump_ensemble = dump_ensembles.get(str(ensemble["ensemble_id"])) or {}
            a_lambda = float(ensemble["aLambda_msbar3"])
            am_ground = channel.get("am_ground_candidate")
            ratio = channel.get("ratio_to_lambda_msbar3_candidate")
            if intercept is None or am_ground is None or ratio is None:
                continue
            ratio = float(ratio)
            am_ground = float(am_ground)
            r_l = float(dump_ensemble.get("am_l") or 0.0) / a_lambda
            r_s = float(dump_ensemble.get("am_s") or 0.0) / a_lambda
            delta_cont = abs(ratio - intercept) * a_lambda
            delta_vol = abs(am_ground) * math.exp(-abs(am_ground) * float(dump_ensemble.get("L") or 0.0))
            delta_chi = abs(ratio) * a_lambda * (
                abs(r_l - (mean_r_l or r_l)) + abs(r_s - (mean_r_s or r_s))
            )
            sigma_sys = math.sqrt(delta_cont * delta_cont + delta_vol * delta_vol + delta_chi * delta_chi)
            channel["published_systematics"] = {
                "status": "complete",
                "sigma_sys": sigma_sys,
                "delta_cont": delta_cont,
                "delta_vol": delta_vol,
                "delta_chi": delta_chi,
            }
            channel["am_ground_sys_err"] = sigma_sys
            if channel.get("published_statistical_error") is not None:
                channel["am_ground_candidate_err"] = math.sqrt(
                    float(channel["published_statistical_error"]) ** 2 + sigma_sys ** 2
                )
            channel["convergence_status"] = (
                "public_unsuppression_ready"
                if channel.get("forward_window_limit_exists")
                and channel.get("published_statistical_error") is not None
                else channel.get("convergence_status")
            )
            channel["sequence_status"] = channel["convergence_status"]
    return out
