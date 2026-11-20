"""Replay the calibrated interacting source."""

import json

from .experiment import run


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
