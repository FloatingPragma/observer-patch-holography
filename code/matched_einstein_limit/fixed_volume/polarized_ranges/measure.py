"""Reproduce the independent curved and native reciprocal range tables."""
import json
from . import native, smooth


if __name__ == "__main__":
    print(json.dumps(dict(smooth=smooth.run(), native=native.run()), indent=2, allow_nan=False))
