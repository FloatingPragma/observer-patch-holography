"""Replay with one math thread from the RER root as described in README.md."""

import json
from pathlib import Path

from .experiment import run
from .passive import run_passive


def serialize(x):
    if isinstance(x, dict):
        return {k: serialize(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [serialize(v) for v in x]
    if isinstance(x, float):
        return round(x, 11)
    return x


if __name__ == "__main__":
    out = Path(__file__).with_name("results.json")
    packet = run()
    packet["passive_boards"] = run_passive()
    out.write_text(json.dumps(serialize(packet), indent=1, allow_nan=False) + "\n")
