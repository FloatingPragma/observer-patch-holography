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

FORMAT = "oph-simulator-run/1"


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
            "status": law.status, "note": law.note}


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
            "title": "Two clocks on the step-3 routes, exact expectations",
            "law": _law(LAWS["reflected-routing"]), "journeys": journeys,
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
