from __future__ import annotations

import dataclasses
import hashlib
import json
import math
import struct
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np


_TYPE_KEY = "__oph_hash_type__"
CANONICAL_HASH_SCHEMA = "oph_typed_canonical_json_v1"


def canonicalize_for_hash(data: Any) -> Any:
    """Return a typed, JSON-safe representation for evidence hashing.

    ``json.dumps(..., default=str)`` is attractive but unsafe for receipts: it
    makes otherwise unsupported objects hash according to an implementation's
    display representation, and ordinary JSON numbers lose distinctions such
    as integer versus binary64 state.  This encoder is deliberately explicit.
    Floating-point values are encoded by their normalized IEEE-754 binary64
    bytes, NumPy containers retain shape/dtype metadata, and unordered
    containers are sorted by their canonical byte representation.

    Negative zero is normalized to positive zero because the simulator treats
    them as the same numeric state.  All NaN payloads are likewise canonical;
    infinities retain their sign.
    """

    if data is None or isinstance(data, (bool, str)):
        return data
    if isinstance(data, (int, np.integer)) and not isinstance(data, (bool, np.bool_)):
        return int(data)
    if isinstance(data, (float, np.floating)):
        return {_TYPE_KEY: "binary64", "hex": _canonical_float64_hex(float(data))}
    if isinstance(data, (complex, np.complexfloating)):
        value = complex(data)
        return {
            _TYPE_KEY: "complex128",
            "realHex": _canonical_float64_hex(float(value.real)),
            "imagHex": _canonical_float64_hex(float(value.imag)),
        }
    if isinstance(data, Path):
        return {_TYPE_KEY: "path", "value": data.as_posix()}
    if isinstance(data, (bytes, bytearray, memoryview)):
        return {_TYPE_KEY: "bytes", "hex": bytes(data).hex()}
    if isinstance(data, np.ndarray):
        array = np.asarray(data)
        return {
            _TYPE_KEY: "ndarray",
            "dtype": _canonical_dtype_name(array.dtype),
            "shape": [int(value) for value in array.shape],
            "values": canonicalize_for_hash(array.tolist()),
        }
    if isinstance(data, np.generic):
        return canonicalize_for_hash(data.item())
    if dataclasses.is_dataclass(data) and not isinstance(data, type):
        return {
            _TYPE_KEY: "dataclass",
            "class": f"{type(data).__module__}.{type(data).__qualname__}",
            "fields": canonicalize_for_hash(dataclasses.asdict(data)),
        }
    if isinstance(data, Mapping):
        rows = [
            [canonicalize_for_hash(key), canonicalize_for_hash(value)]
            for key, value in data.items()
        ]
        rows.sort(key=lambda row: canonical_json_bytes(row[0]))
        return {_TYPE_KEY: "mapping", "items": rows}
    if isinstance(data, (set, frozenset)):
        values = [canonicalize_for_hash(value) for value in data]
        values.sort(key=canonical_json_bytes)
        return {_TYPE_KEY: "set", "items": values}
    if isinstance(data, tuple):
        return {_TYPE_KEY: "tuple", "items": [canonicalize_for_hash(value) for value in data]}
    if isinstance(data, Sequence):
        return [canonicalize_for_hash(value) for value in data]
    raise TypeError(
        "unsupported value in canonical evidence hash: "
        f"{type(data).__module__}.{type(data).__qualname__}"
    )


def canonical_json_bytes(data: Any, *, already_canonical: bool = False) -> bytes:
    """Serialize *data* to the canonical byte stream used by receipts."""

    payload = data if already_canonical else canonicalize_for_hash(data)
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def stable_json_hash(data: Any) -> str:
    return "sha256:" + hashlib.sha256(canonical_json_bytes(data)).hexdigest()


def _canonical_float64_hex(value: float) -> str:
    if math.isnan(value):
        return "7ff8000000000000"
    if value == 0.0:
        value = 0.0
    return struct.pack(">d", value).hex()


def _canonical_dtype_name(dtype: np.dtype[Any]) -> Any:
    normalized = np.dtype(dtype)
    if normalized.fields:
        return canonicalize_for_hash(normalized.descr)
    if normalized.kind in {"b", "i", "u", "f", "c", "S", "U", "V"}:
        return f"{normalized.kind}{normalized.itemsize}"
    return normalized.name
