"""Acquire and inspect the successor OW-30 detector's three enrolled boards.

The verified launcher authenticates every binding and supplies the selectors.
"""

from hashlib import sha256
import json
from pathlib import Path
import re
from tempfile import NamedTemporaryFile, TemporaryDirectory
from threading import RLock
from unittest.mock import patch

import numpy as np

from .. import coupling_meter, experiment, formation, holdout
from ..registration import verify as legacy
from model.executable_model import federation

EXIT_CODES = {"passed": 0, "rejected": 2, "incomplete": 3, "invalid": 4}
HERE = Path(__file__).resolve().parent
PILOT_SEEDS = (holdout.SEED,) * 3
CANDIDATE_FILE = HERE / "candidate.json"
REPEAT_SEEDS = (374621, 614357, 801173)
VALIDATION_SEEDS = (785029, 514823, 690271)
LOCAL_SELECTORS = {
    "pilot": PILOT_SEEDS,
    "repeat": REPEAT_SEEDS,
    "validation": VALIDATION_SEEDS,
}
BINDING_KEYS = {
    "mode",
    "anchor",
    "source_commit",
    "candidate_pin",
    "candidate_sha256",
    "protocol_sha256",
}
GATE_KEYS = {
    "changed_reader",
    "meter",
    "commissioned",
    "lattice",
    "closure",
    "live",
    "changed",
    "transport",
    "quiet",
    "shape",
    "continuum",
}
# The reused acquisition functions reference module-level names. Concurrent
# Python callers must not observe another run's temporary selector or admission.
_ACQUISITION = RLock()


class UnmeasurableStrand(Exception):
    """The enrolled geometry cannot support the predeclared measurement."""


def valid_inputs(binding, seeds):
    if not isinstance(binding, dict) or set(binding) != BINDING_KEYS:
        return False
    if (
        not isinstance(seeds, (tuple, list))
        or len(seeds) != 3
        or any(type(s) is not int for s in seeds)
    ):
        return False

    def hex_value(key, width):
        value = binding[key]
        return (
            isinstance(value, str)
            and re.fullmatch(f"[a-f0-9]{{{width}}}", value) is not None
        )

    if not (
        hex_value("source_commit", 40)
        and binding["candidate_pin"] == binding["source_commit"]
        and hex_value("protocol_sha256", 64)
        and hex_value("candidate_sha256", 64)
    ):
        return False
    if binding["mode"] in LOCAL_SELECTORS:
        return (
            binding["anchor"] is None
            and tuple(seeds) == LOCAL_SELECTORS[binding["mode"]]
        )
    return (
        binding["mode"] == "public"
        and hex_value("anchor", 40)
        and len(set(seeds)) == 3
        and all(2**63 <= s < 2**64 for s in seeds)
    )


def matching_bytes(binding):
    try:
        return binding["protocol_sha256"] == legacy.file_hash(
            HERE / "TASK.md"
        ) and binding["candidate_sha256"] == legacy.file_hash(CANDIDATE_FILE)
    except OSError:
        return False


def write_json(path, value):
    with NamedTemporaryFile("w", dir=path.parent, delete=False) as stream:
        try:
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write("\n")
        except BaseException:
            Path(stream.name).unlink(missing_ok=True)
            raise
    Path(stream.name).replace(path)


def read_json(path):
    return json.loads(path.read_text())


def same_readings(a, b):
    if type(a) is not type(b):
        return False
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(same_readings(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(same_readings(x, y) for x, y in zip(a, b))
    if isinstance(a, float):
        return bool(np.isfinite(a) and np.isfinite(b) and abs(a - b) <= 1e-8)
    return a == b


def enrollment(level, seed):
    routes = experiment.grow(federation(level), seed)
    route = max((p for closed, p in routes if not closed), key=len, default=[])
    cells = sorted({port // 12 for dimer in route for port in dimer})
    reason = None
    if not route:
        reason = "no open component"
    elif len(cells) < 2:
        reason = "enrolled component lies in a single cell"
    elif len(route) + 1 < 4:
        reason = "fewer than four distinct launch profiles exist"
    return {
        "level": level,
        "seed": seed,
        "route": [list(p) for p in route],
        "cells": cells,
        "ineligible": reason,
    }


def expected_counts(n):
    m = n + 1
    last = min(8, max(4, (m + 1) // 4))
    # (1,b), (2,3), and withheld (2,b) for b >= 4.
    pairs = [(a, b) for a in (1, 2) for b in range(a + 1, last + 1)]
    specs = {f"clock_{hand}": ((experiment.SAMPLES,), 2**20) for hand in ("cos", "sin")}
    for a, b in pairs:
        for tag in ("sum", "changed", "scan0", "scan1"):
            shape = (
                (experiment.SAMPLES, m)
                if tag.startswith("scan")
                else (experiment.SAMPLES,)
            )
            specs[f"p{a}_{b}_{tag}"] = (shape, experiment.SHOTS)
    for arm in ("native", "control"):
        for time in (0, 2, 4, 8):
            specs[f"{arm}_flight_{time}"] = ((2, m), 2**20)
    return specs


def audit_counts(folder, specs):
    found = {p.stem for p in folder.glob("*.npy")}
    missing = sorted(set(specs) - found)
    invalid = {k: "unexpected count array" for k in sorted(found - set(specs))}
    words = {}
    for name in sorted(found.intersection(specs)):
        shape, bound = specs[name]
        path = folder / (name + ".npy")
        try:
            if path.is_symlink():
                raise ValueError("symlink in release")
            array = np.load(path, allow_pickle=False, mmap_mode="r")
            if (
                array.dtype.str != "<u8"
                or array.shape != shape
                or not array.flags.c_contiguous
            ):
                raise ValueError("incorrect type, orientation or shape")
            if 0 in shape or int(array.max()) > bound:
                raise ValueError("count exceeds trial budget")
            # NPY payload is contiguous; reject unparsed trailing data as well.
            if array.offset + array.nbytes != path.stat().st_size:
                raise ValueError("unexpected trailing bytes")
            words[name] = dict(
                shape=list(shape),
                dtype="<u8",
                trials=bound,
                sha256=sha256(memoryview(array).cast("B")).hexdigest(),
            )
        except (ValueError, OSError, EOFError, TypeError) as err:
            invalid[name] = str(err)
    return dict(expected=len(specs), files=words, missing=missing, invalid=invalid)


def worst(statuses):
    return max(statuses, key=EXIT_CODES.__getitem__)


def measured_meter():
    path = experiment.HERE / "coupling_meter.json"
    bench = read_json(path)
    values, _ = coupling_meter.measure(
        np.asarray(bench["all_counts"]),
        np.asarray(bench["hand_calibration"]) / bench["calibration_trials"],
    )
    return dict(
        source_sha256=legacy.file_hash(path),
        h=float(values[0]),
        j=float(values[1]),
        m=float(values[2]),
        c_squared=float(values[3]),
    )


def even_hold_dispersion(physical, spatial, meter):
    """Use time-counter frequencies and even-scan k, with coefficients from the bench."""
    modes = physical["observer"]["modes"]
    omega = np.asarray(spatial["omega_from_pair_counts"])
    k = np.asarray(spatial["even_trained_k"])
    if not np.allclose(omega, [r["omega"] for r in modes], atol=1e-8, rtol=0):
        raise ValueError(
            "reported branch frequencies disagree with raw temporal counters"
        )
    m, c2 = meter["m"], meter["c_squared"]
    q2 = 2 - 2 * np.cos(k)
    lattice = np.sqrt(m * m + c2 * q2)
    ordinary = np.sqrt(m * m + c2 * k * k)
    tangent = m + c2 * q2 / (2 * m)
    low = (k > 0) & (k <= 0.50)
    low_count = int(np.count_nonzero(low))
    low_span = float(np.ptp(k[low])) if low_count else 0.0
    low_max = float(np.max(np.abs(ordinary[low] - omega[low]))) if low_count else None
    lattice_max = float(np.max(np.abs(lattice - omega)))

    def anchored_rms(forecast):
        residual = (forecast[1:] - forecast[0]) - (omega[1:] - omega[0])
        return float(np.sqrt(np.mean(residual**2)))

    lattice_rms, tangent_rms = anchored_rms(lattice), anchored_rms(tangent)
    rows = []
    for i, source in enumerate(modes):
        rows.append(
            dict(
                launch=source["spatial_launch"],
                k=float(k[i]),
                omega=float(omega[i]),
                full_recurrence_diagnostic_k=source["k"],
                low_k=bool(low[i]),
                lattice=float(lattice[i]),
                ordinary=float(ordinary[i]),
                matched_tangent=float(tangent[i]),
                lattice_residual=float(omega[i] - lattice[i]),
                ordinary_residual=float(omega[i] - ordinary[i]),
                tangent_residual=float(omega[i] - tangent[i]),
            )
        )
    return dict(
        spatial_ruler="even_hold_fit",
        meter=meter,
        modes=rows,
        low_k_lines=low_count,
        low_k_span=low_span,
        low_k_max=low_max,
        lattice_max=lattice_max,
        shape_rms_lattice=lattice_rms,
        shape_rms_tangent=tangent_rms,
        gates=dict(
            lattice=bool(np.all(omega > 0) and lattice_max < 0.0048),
            continuum=bool(low_count >= 2 and low_span >= 0.02 and low_max < 0.0048),
            shape=bool(lattice_rms < 0.5 * tangent_rms),
        ),
    )


def board_receipt(plan, physical, folder, meter):
    from .spatial import score as spatial_score

    evidence = {"enrollment": plan}
    if plan["ineligible"]:
        # The source guard throws this sentinel before spectroscopy, so a
        # hardware outage cannot masquerade as a measured enrollment failure.
        valid = (
            physical
            and physical.get("error") == "UnmeasurableStrand"
            and physical.get("detail") == plan["ineligible"]
            and (folder / "failure.json").is_file()
        )
        evidence["counts"] = audit_counts(folder, {})
        evidence["status"] = (
            "invalid"
            if evidence["counts"]["invalid"]
            else "rejected"
            if valid
            else "incomplete"
        )
        return evidence
    inventory = audit_counts(folder, expected_counts(len(plan["route"])))
    evidence["counts"] = inventory
    if inventory["invalid"]:
        evidence["status"] = "invalid"
    elif inventory["missing"] or not physical or "observer" not in physical:
        evidence["status"] = "incomplete"
    else:
        recorded = physical["observer"]["count_sha256"]
        actual = {
            k: v["sha256"] for k, v in inventory["files"].items() if "_flight_" not in k
        }
        if recorded != actual or physical["route"] != plan["route"]:
            evidence.update(
                status="invalid", reason="count digests or enrolled route disagree"
            )
            return evidence
        if set(physical["gates"]) != GATE_KEYS or any(
            type(x) is not bool for x in physical["gates"].values()
        ):
            return dict(
                evidence, status="invalid", reason="missing or invalid physical gate"
            )
        for label, key in (("native", "transport"), ("control", "changed_transport")):
            flight = physical[key]
            if [f["t"] for f in flight["flights"]] != [0, 2, 4, 8] or flight[
                "trials_per_arm"
            ] != 2**20:
                return dict(
                    evidence,
                    status="invalid",
                    reason="invalid transport schedule or trial count",
                )
            for reading in flight["flights"]:
                arms = np.load(
                    folder / f"{label}_flight_{reading['t']}.npy", allow_pickle=False
                )
                difference = arms[0].astype(np.int64) - arms[1].astype(np.int64)
                if not np.array_equal(reading["differences"], difference):
                    return dict(
                        evidence,
                        status="invalid",
                        reason="terminal differences disagree with raw counts",
                    )
        spatial = spatial_score(physical, folder)
        gates = dict(physical["gates"])
        dispersion = None
        if "even_trained_k" in spatial:
            dispersion = even_hold_dispersion(physical, spatial, meter)
            gates.update(dispersion["gates"])
        else:
            gates.update(lattice=False, continuum=False, shape=False)
        if spatial["status"] == "missing_records":
            outcome = "incomplete"
        elif spatial["status"] == "invalid_records":
            outcome = "invalid"
        else:
            accepted = (
                spatial["records_complete"]
                and spatial["profile_pass"]
                and spatial["launch_label_rejected"]
                and all(gates.values())
            )
            outcome = "passed" if accepted else "rejected"
        evidence.update(
            status=outcome,
            native_gates=gates,
            legacy_recurrence_gates=physical["gates"],
            legacy_recurrence_score=physical["score"],
            spatial=spatial,
            dispersion=dispersion,
        )
    if physical and "error" in physical:
        evidence["acquisition_error"] = {k: physical[k] for k in ("error", "detail")}
    return evidence


def measure_release(dest, draw):
    records = dest / "records"
    result_path = records / "result.json"
    raw = read_json(result_path) if result_path.exists() else {}
    boards = raw.get("boards", [])
    raw_draw = (
        read_json(records / "draw.json") if (records / "draw.json").exists() else None
    )
    anchor = draw["binding"]["anchor"] or "local-" + draw["binding"]["mode"]
    terms_hash = legacy.file_hash(HERE / "TASK.md")
    expected_raw = dict(
        registration=anchor,
        source_commit=draw["binding"]["source_commit"],
        adopted_pin=draw["binding"]["candidate_pin"],
        protocol_sha256=terms_hash,
        seeds=draw["seeds"],
    )
    metadata = sorted(records.rglob("*.json"))
    receipt = dict(
        draw=draw,
        boards=[],
        metadata_sha256={
            p.relative_to(records).as_posix(): legacy.file_hash(p) for p in metadata
        },
    )
    if not matching_bytes(draw["binding"]):
        return dict(
            receipt,
            status="invalid",
            reason="frozen terms or candidate declaration differ from the binding",
        )
    if raw_draw is not None and raw_draw != expected_raw:
        return dict(
            receipt, status="invalid", reason="raw draw differs from bound inputs"
        )
    if raw and (raw.get("registration") != anchor or len(boards) != 3):
        return dict(
            receipt, status="invalid", reason="incorrect board count or registration"
        )
    meter = measured_meter()
    if raw and not np.allclose(
        [raw["m"], raw["c_squared"]],
        [meter["m"], meter["c_squared"]],
        atol=1e-8,
        rtol=0,
    ):
        return dict(
            receipt,
            status="invalid",
            reason="meter forecasts differ from the disjoint bench counts",
        )
    for level, seed in enumerate(draw["seeds"]):
        plan = enrollment(level, seed)
        board = boards[level] if boards else None
        if board is not None and (board["level"] != level or board["seed"] != seed):
            return dict(receipt, status="invalid", reason="raw board identity differs")
        receipt["boards"].append(
            board_receipt(plan, board, records / str(level), meter)
        )
    allowed = {records / file for file in ("draw.json", "result.json")}
    expected_dirs = {records / str(i) for i in range(3)}
    for i, row in enumerate(receipt["boards"]):
        folder = records / str(i)
        if not boards:
            allowed.add(folder / "failure.json")
        elif "error" in boards[i]:
            failure = folder / "failure.json"
            if failure.exists() and read_json(failure) != boards[i]:
                return dict(
                    receipt, status="invalid", reason="board failure receipts disagree"
                )
            if not failure.exists():
                row["status"] = "incomplete"
            allowed.add(failure)
        allowed.update(folder / (k + ".npy") for k in row["counts"]["files"])
    if any(
        p.is_symlink()
        or (p.is_dir() and p not in expected_dirs)
        or (p.is_file() and p not in allowed)
        for p in records.rglob("*")
    ):
        return dict(
            receipt,
            status="invalid",
            reason="unexpected file or directory in raw release",
        )
    status = worst(b["status"] for b in receipt["boards"])
    if raw_draw is None or not result_path.exists():
        status = worst((status, "incomplete"))
    return dict(receipt, status=status)


def acquire(binding, seeds, dest):
    """Issue exactly one attempt; a retry of an existing release only inspects it."""
    if not valid_inputs(binding, seeds):
        return "invalid"
    if not matching_bytes(binding):
        return "invalid"
    dest = Path(dest)
    try:
        dest.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        return inspect_release(binding, seeds, dest)
    draw = dict(binding=dict(binding), seeds=list(seeds))
    write_json(dest / "draw.json", draw)
    write_json(
        dest / "manifest.json", dict(draw=draw, status="incomplete", stage="acquiring")
    )
    try:
        # Reuse the count-generating body and every original physical gate. The
        # guarded driver supplies the new draw and verified source binding.
        with _ACQUISITION:
            original_prepare = formation.prepare

            def prepare(level, seed):
                plan = enrollment(level, seed)
                if plan["ineligible"]:
                    raise UnmeasurableStrand(plan["ineligible"])
                return original_prepare(level, seed)

            with (
                patch.object(legacy, "HERE", HERE),
                patch.object(legacy, "settings", return_value=list(seeds)),
                patch.object(legacy, "PIN", binding["candidate_pin"]),
                patch.object(formation, "prepare", prepare),
            ):
                anchor = binding["anchor"] or "local-" + binding["mode"]
                legacy.run_bank(anchor, dest / "records", binding["source_commit"])
        report = measure_release(dest, draw)
    except Exception as err:
        status = (
            "invalid"
            if isinstance(err, (AssertionError, ValueError, KeyError))
            else "incomplete"
        )
        try:
            report = measure_release(dest, draw)
        except Exception:
            report = dict(draw=draw, status=status)
        report.update(
            status=worst((status, report["status"])),
            error=type(err).__name__,
            detail=str(err),
        )
    write_json(dest / "manifest.json", report)
    return report["status"]


def inspect_release(binding, seeds, release):
    """Verify the saved inventory and spatial decision without reacquiring."""
    if not valid_inputs(binding, seeds) or not matching_bytes(binding):
        return "invalid"
    release = Path(release)
    try:
        draw = read_json(release / "draw.json")
        if draw != dict(binding=binding, seeds=list(seeds)):
            return "invalid"
        stored = read_json(release / "manifest.json")
        if stored.get("draw") != draw:
            return "invalid"
        if stored.get("status") == "invalid":
            return "invalid"
        if "boards" not in stored:
            return (
                stored["status"]
                if stored["status"] in ("incomplete", "invalid")
                else "invalid"
            )
        current = measure_release(release, draw)
        if current["status"] in ("invalid", "incomplete"):
            return current["status"]
        return current["status"] if same_readings(current, stored) else "invalid"
    except FileNotFoundError:
        return "incomplete"
    except (ValueError, KeyError, TypeError, IndexError, OSError):
        return "invalid"


def replay(binding, seeds, release):
    """Retain the release and compare it with one fresh run on its original draw."""
    status = inspect_release(binding, seeds, release)
    if status in ("invalid", "incomplete"):
        return status
    release = Path(release)
    with TemporaryDirectory() as tmp:
        repeated = Path(tmp) / "repeat"
        observed = acquire(binding, seeds, repeated)
        if observed in ("invalid", "incomplete"):
            return observed
        if observed != status:
            return "invalid"
        try:
            a, b = (read_json(p / "manifest.json") for p in (release, repeated))
            # Raw count payload hashes, integer types and shapes are exact.
            # JSON metadata also carries calibrated floats, compared to 1e-8.
            if a["metadata_sha256"].keys() != b["metadata_sha256"].keys():
                return "invalid"
            for name in a["metadata_sha256"]:
                original, regenerated = (
                    read_json(p / "records" / name) for p in (release, repeated)
                )
                if not same_readings(original, regenerated):
                    return "invalid"
            if not same_readings(a["boards"], b["boards"]):
                return "invalid"
        except (ValueError, TypeError, AssertionError, OSError):
            return "invalid"
    return status


def pilot_summary(release):
    manifest = read_json(Path(release) / "manifest.json")
    if (
        manifest["draw"]["binding"]["mode"] != "pilot"
        or tuple(manifest["draw"]["seeds"]) != PILOT_SEEDS
    ):
        raise ValueError("The golden receipt requires the fixed pilot selectors")
    boards = []
    for board in manifest["boards"]:
        counts = board["counts"]["files"]
        inventory = json.dumps(counts, sort_keys=True, separators=(",", ":")).encode()
        boards.append(
            dict(
                level=board["enrollment"]["level"],
                dimers=len(board["enrollment"]["route"]),
                status=board["status"],
                raw_arrays=len(counts),
                inventory_sha256=sha256(inventory).hexdigest(),
                failed_native_gates=[
                    k for k, ok in board.get("native_gates", {}).items() if not ok
                ],
                spatial=board.get("spatial"),
                dispersion=board.get("dispersion"),
            )
        )
    return dict(
        scope="Exposed v3 internal pilot",
        seeds=manifest["draw"]["seeds"],
        status=manifest["status"],
        boards=boards,
    )
