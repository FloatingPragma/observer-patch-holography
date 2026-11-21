"""The registry of candidate laws: every package in model/candidate_laws that has a candidate.py."""
from __future__ import annotations

import importlib
from pathlib import Path

CANDIDATES = Path(__file__).resolve().parents[1] / 'candidate_laws'


def names() -> list[str]:
    return sorted(p.parent.name for p in CANDIDATES.glob('*/candidate.py'))


def load(name: str):
    """The candidate.py module of one candidate law."""
    if name not in names():
        raise SystemExit(f'unknown candidate law {name!r}; candidates: {", ".join(names())}')
    return importlib.import_module(f'model.candidate_laws.{name}.candidate')


def main_candidate() -> str:
    mains = [n for n in names() if load(n).ROLE == 'main']
    if len(mains) != 1:
        raise SystemExit(f'exactly one candidate must have ROLE "main", found {mains}')
    return mains[0]
