#!/usr/bin/env python3
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from qiskit_ibm_runtime import QiskitRuntimeService


def load_api_key(credentials_file: Path) -> str:
    text = credentials_file.read_text().strip()
    match = re.fullmatch(r"IBM cloud API key:\s*(\S+)", text)
    api_key = match.group(1) if match else text
    if not api_key or any(character.isspace() for character in api_key):
        raise ValueError(
            f"IBM cloud API key not found in {credentials_file}; expected one raw token "
            "or 'IBM cloud API key: <token>'"
        )
    return api_key


def get_service(credentials_file: Path) -> QiskitRuntimeService:
    return QiskitRuntimeService(
        channel="ibm_quantum_platform",
        token=load_api_key(credentials_file),
    )


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
