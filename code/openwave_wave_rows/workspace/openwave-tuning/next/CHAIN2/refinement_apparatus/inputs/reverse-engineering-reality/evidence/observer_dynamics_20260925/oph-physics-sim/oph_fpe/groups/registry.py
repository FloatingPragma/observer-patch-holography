from __future__ import annotations

from oph_fpe.groups.base import FiniteGroup
from oph_fpe.groups.clock import build_clock
from oph_fpe.groups.s3 import build_s3
from oph_fpe.groups.z2 import build_z2


def get_group(name: str) -> FiniteGroup:
    normalized = name.strip().upper()
    if normalized == "Z2":
        return build_z2()
    if normalized == "S3":
        return build_s3()
    if normalized.startswith("C") and normalized[1:].isdigit():
        return build_clock(int(normalized[1:]))
    raise ValueError(f"unknown finite group: {name}")
