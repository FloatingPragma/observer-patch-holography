"""Read the immutable native charge representation and gamma convention."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .._source import verify_sources

ROOT = Path(__file__).resolve().parents[3]
FILES = (
    "model/observer_particles/results/representation_identification_operators.json",
    "model/observer_particles/results/representation_identification.json",
    "model/common_source/overlap_regulator.py",
)


def load(repo=ROOT):
    """Verify consumed source bytes before reading the native integer charges."""
    repo = Path(repo).resolve()
    provenance = verify_sources(repo, FILES)
    raw = json.loads((repo / FILES[0]).read_text())
    owner = json.loads((repo / FILES[1]).read_text())
    if hashlib.sha256((repo / FILES[0]).read_bytes()).hexdigest() != owner["artifacts"][Path(FILES[0]).name]:
        raise ValueError("Native representation receipt does not bind the operator packet")
    q = raw["even_identification_operators"]["q"]
    if q["shape"] != [15, 15] or len(q["entries"]) != 15:
        raise ValueError("Expected the complete native M15 integer charge matrix")
    entries = {i: int(value) for i, j, value in q["entries"] if i == j and str(int(value)) == value}
    if set(entries) != set(range(15)):
        raise ValueError("Native charge is not the complete diagonal integer operator")
    charges = tuple(entries[i] for i in range(15))
    blocks = [{"q6Y": charge, "indices": tuple(i for i, q in enumerate(charges) if q == charge),
               "dimension": charges.count(charge)} for charge in sorted(set(charges))]
    if {b["q6Y"]: b["dimension"] for b in blocks} != {-4: 3, -3: 2, 1: 6, 2: 3, 6: 1}:
        raise ValueError("Native charge representation changed")
    # Module execution belongs to this repository. An alternate source path is
    # usable only when it agrees with these immutable consumed source bytes.
    from model.common_source.overlap_regulator import gamma_matrices
    if repo != ROOT:
        verify_sources(ROOT, FILES)
    return {"repo": repo, "charges": charges, "blocks": blocks,
            "gamma": gamma_matrices(), "provenance": provenance}
