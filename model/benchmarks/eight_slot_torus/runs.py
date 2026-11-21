"""Portable run files: one gzip-compressed JSON document per simulation.

A run can be computed on any machine and replayed in the web app. The
format is versioned; arrays are stored as base64 of zlib-compressed bytes.
"""
from __future__ import annotations

import base64
import gzip
import json
import platform
import time
import zlib
from pathlib import Path

import numpy as np

FORMAT = "eight-slot-torus-run/1"


def _pack(array: np.ndarray, dtype: str) -> dict:
    data = np.ascontiguousarray(array, dtype=dtype)
    return {"dtype": dtype, "shape": list(data.shape),
            "zlib_base64": base64.b64encode(zlib.compress(data.tobytes(), 9)).decode()}


def unpack(block: dict) -> np.ndarray:
    raw = zlib.decompress(base64.b64decode(block["zlib_base64"]))
    return np.frombuffer(raw, dtype=block["dtype"]).reshape(block["shape"])


def _law(law) -> dict:
    return {"key": law.key, "title": law.title,
            "equal_phase": [float(np.real(law.equal_phase)), float(np.imag(law.equal_phase))],
            "encounter_phase": float(law.encounter_phase), "access": law.access,
            "status": law.status, "note": law.note,
            "rail_contact_phases": [float(v) for v in law.rail_contact_phases],
            "rail_extra_phases": [float(v) for v in law.rail_extra_phases]}


def _provenance(seconds: float) -> dict:
    return {"created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "python": platform.python_version(), "numpy": np.__version__,
            "machine": platform.machine(), "seconds": round(seconds, 2)}


def continuum_document(result: dict, convergence: list[dict] | None = None) -> dict:
    """A field3d run: the law's scalar field at every update with the continuum it is compared with."""
    size = result["size"]
    scale = result["norm"]
    frames = np.stack(result["frames"]) / scale
    quantized = np.clip(np.round(frames * 127), -127, 127)
    return {
        "format": FORMAT,
        "kind": "field3d",
        "title": "Scalar field released at rest, law and continuum",
        "law": _law(result["law"]),
        "carrier": {"size": size, "chart_spacing": 1 / size, "private_order": result["private_order"],
                    "sector_angle": 2 * np.pi * result["private_order"] / size},
        "preparation": {"modes": [[a, list(k)] for a, k in result["modes"]], "released_at_rest": True},
        "continuum": {"equation": "phi_tt = (1/4) Laplacian phi - (Q^2/4) phi",
                      "c": 0.5, "Q": 2 * np.pi * result["private_order"],
                      "rest_angular_frequency": np.pi * result["private_order"]},
        "updates": result["updates"],
        "chart_time": result["chart_time"],
        "frames": {**_pack(quantized, "int8"), "scale": 1 / 127,
                   "note": "phi / max|phi(0)| quantized to 1/127; the error table uses full precision"},
        "stations": result["stations"].tolist(),
        "station_law": np.round(result["station_law"] / scale, 6).tolist(),
        "station_continuum": np.round(result["station_continuum"] / scale, 6).tolist(),
        "errors": [[round(t, 6), e, i] for t, e, i in result["errors"]],
        "max_error": result["max_error"],
        "convergence": convergence or [],
        "provenance": _provenance(result["seconds"]),
    }


def journeys_document(journeys: list[dict], seconds: float) -> dict:
    from .law import LAWS
    return {"format": FORMAT, "kind": "clock-journeys",
            "title": "Two clocks on the diagonal and axial routes, exact expectations",
            "law": _law(LAWS["torus"]), "journeys": journeys,
            "provenance": _provenance(seconds)}


def write(document: dict, path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = json.dumps(document, separators=(",", ":"), allow_nan=False).encode()
    path.write_bytes(gzip.compress(data, mtime=0))
    return path


def read(path: Path) -> dict:
    document = json.loads(gzip.decompress(Path(path).read_bytes()))
    if document.get("format") != FORMAT:
        raise ValueError(f"unknown run format {document.get('format')!r}")
    return document


def journey_document(result: dict, title: str, step: str, summary: str) -> dict:
    """An animated journey: traveller and stationary twin densities, clock advances, forecast."""
    density = result["density"]  # (F, 2, M, M): traveller, twin
    scale = density.max(axis=(0, 2, 3), keepdims=True)
    scaled = np.sqrt(np.clip(density / scale, 0, 1))
    return {
        "format": FORMAT, "kind": "journey", "title": title, "step": step, "summary": summary,
        "law": _law(result["law"]),
        "carrier": {"size": result["size"], "stride": result["stride"]},
        "route": result["route"], "ticks": result["ticks"], "centre": np.round(result["centre"], 2).tolist(),
        "speed": result["speed"], "reunion": result["reunion"], "forecast": result["forecast"],
        "reference_speed": 1 / 2.002604166666667,
        "frames": {**_pack(np.round(scaled * 255), "uint8"),
                   "note": "sqrt of probability over (x, y), sampled every stride addresses; "
                           "channel 0 traveller, channel 1 stationary twin"},
        "clocks": [{**c, "rest": np.round(c["rest"], 7).tolist(), "moving": np.round(c["moving"], 7).tolist()}
                   for c in result["clocks"]],
        "provenance": _provenance(result["seconds"]),
    }


def race_document(result: dict) -> dict:
    strip = result["strip"]
    logp = np.log10(np.maximum(strip, 1e-300))
    encoded = np.clip(np.round((logp + 60) / 60 * 255), 0, 255)
    return {
        "format": FORMAT, "kind": "race", "step": "Signals",
        "title": "A fresh message on the free wave",
        "summary": "Field-only access: a compact launch at A, read by a detector at B. The law's strict "
                   "front moves one address per update and carries only exponentially small precursors; "
                   "the message arrives with the smooth cone at c = 1/2.",
        "law": _law(result["law"]),
        "carrier": {"size": result["size"]}, "source": result["source"], "receiver": result["receiver"],
        "distance": result["distance"], "stride": result["stride"], "deadlines": result["deadlines"],
        "strip": {**_pack(encoded, "uint8"), "note": "log10 probability per address, mapped from [-60, 0] to [0, 255]"},
        "detector": [float(f"{v:.6e}") for v in result["detector"]],
        "first_nonzero": result["first_nonzero"], "largest_before_early": result["largest_before_early"],
        "largest_before_smooth": result["largest_before_smooth"], "received_by_late": result["received_by_late"],
        "provenance": _provenance(result["seconds"]),
    }


def scan_document(title: str, step: str, summary: str, x_label: str, y_label: str, series: list,
                  curves: list, law, seconds: float, table: list) -> dict:
    return {"format": FORMAT, "kind": "scan", "title": title, "step": step, "summary": summary,
            "law": _law(law), "x_label": x_label, "y_label": y_label, "series": series,
            "curves": curves, "table": table, "provenance": _provenance(seconds)}


def formation_document(result: dict, control: dict) -> dict:
    joint = result["joint"]
    scaled = np.sqrt(np.clip(joint / joint.max(), 0, 1))
    return {
        "format": FORMAT, "kind": "formation", "step": "Matter",
        "title": "A pair formed in an encounter",
        "summary": "Two constituents launched toward each other on a strip. Where they meet head-on, the "
                   "capture vertex turns the collision into a bound pair plus a marker that leaves on the "
                   "junction's own wire. The pair stays together under the contact law, and as it matures "
                   "it opens its hop and travels on.",
        "law": _law(result["law"]),
        "carrier": {"size": result["size"]}, "updates": result["updates"], "orders": result["orders"],
        "angle": result["angle"], "contact_phase": result["contact_phase"],
        "joint": {**_pack(np.round(scaled * 255), "uint8"),
                  "note": "sqrt of the probability of the unmarked pair over (z_A, z_B)"},
        "bound": np.round(result["bound"], 9).tolist(),
        "collected": np.round(result["collected"] / max(result["clicks"][-1], 1e-300), 9).tolist(),
        "captured": [round(v, 9) for v in result["captured"]],
        "clicks": [float(f"{v:.6e}") for v in result["clicks"]],
        "captured_without_vertex": control["captured"][-1],
        "markers": [[z, n, float(f"{w:.4e}")] for z, n, w in result["markers"]],
        "birth_centre": result["birth_centre"],
        "forward_fraction": result["forward_fraction"],
        "forward_fraction_collected": result["forward_fraction_collected"],
        "provenance": _provenance(result["seconds"]),
    }


def matter_wave_document(result: dict, frames: np.ndarray, stride: int) -> dict:
    """The formed pair's matter-wave relation, with a resting and a moving packet."""
    amplitude = np.abs(frames)
    scale = amplitude.max()
    encoded = np.stack([np.round(np.sqrt(amplitude / scale) * 255),
                        np.round((np.angle(frames) + np.pi) / (2 * np.pi) * 255) % 256], axis=2)
    return {
        "format": FORMAT, "kind": "dispersion", "step": "Matter waves",
        "title": "The matter-wave relation of a collision-formed pair",
        "summary": "The frequency relation of a formed pair follows from the mature step of the law. For each of the fourteen registered order-three directions at "
                   "chart periods 64, 128 and 256 the frequency comes from the native step, and the "
                   "observer's test asks whether the kinetic increment (omega^2 - m^2) / (c^2 |k|^2) "
                   "lies within 3.5% of one, with m the rest frequency and c^2 the scalar standard of "
                   "the unbound walk. The changed law moves half the spatial weight of the formed "
                   "vector into the private slots.",
        "law": _law(result["law"]),
        "phi": result["phi"], "rest": result["rest"], "coefficient": result["coefficient"],
        "c2": result["c2"], "tolerance": 0.035, "identity_error": result["identity_error"],
        "rows": [{**r, "omega": round(r["omega"], 9), "kinetic": round(r["kinetic"], 7)} for r in result["rows"]],
        "recorded": result["recorded"],
        "stride": stride, "updates": (frames.shape[0] - 1) * stride,
        "frames": {**_pack(encoded, "uint8"),
                   "note": "axis 1: resting, moving packet; axis 2: sqrt amplitude, phase of the S-reading"},
        "provenance": _provenance(result["seconds"]),
    }


def pair_persistence_document(result: dict) -> dict:
    """The protected circulation of a bound pair and the independent-transport control."""
    scale = max(result["pair_density"].max(), result["control_density"].max())
    encode = lambda a: np.round(np.sqrt(np.clip(a / scale, 0, 1)) * 255)
    return {
        "format": FORMAT, "kind": "persistence", "step": "Pair persistence",
        "title": "A bound pair that persists in free evolution",
        "summary": "A signed circulation of a bound pair around one plaquette, read by observers at every "
                   "update. The law keeps it whole, and a returning marker cannot release it. The control "
                   "gives the two constituents independent transport.",
        "law": _law(result["law"]),
        "phi": result["phi"], "size": result["size"], "updates": result["updates"], "ages": result["ages"],
        "eigenvalue": result["eigenvalue"], "eigen_error": result["eigen_error"],
        "profile_change": result["profile_change"], "capture_projection": result["capture_projection"],
        "contact_matched": [round(v, 12) for v in result["contact_matched"]],
        "control_matched": [round(v, 9) for v in result["control_matched"]],
        "frames": {**_pack(np.stack([encode(result["pair_density"]), encode(result["control_density"])], axis=1), "uint8"),
                   "note": "sqrt of the probability of constituent A over (x, y); channel 0 the law, channel 1 the control"},
        "recorded": result["recorded"],
        "provenance": _provenance(result["seconds"]),
    }
