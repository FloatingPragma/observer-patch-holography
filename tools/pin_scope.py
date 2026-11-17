"""tools/ and model/ entry point for code/pin_scope.py (one implementation)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "_oph_pin_scope", Path(__file__).resolve().parents[1] / "code" / "pin_scope.py"
)
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)

exempt = _MODULE.exempt
bound = _MODULE.bound
bound_pins = _MODULE.bound_pins
repo_sha256 = _MODULE.repo_sha256
scoped = _MODULE.scoped
