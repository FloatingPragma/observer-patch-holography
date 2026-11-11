"""Replay the exact finite repair-to-quotient bridge."""
import json

from .bridge import report


if __name__ == "__main__":
    print(json.dumps(report(), indent=2, sort_keys=True))
