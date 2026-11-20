"""Check the candidate laws.

Runs the exactness and regression checks of the simulator that executes the three
candidate laws and their components (``python -m simulator check``).
"""
from __future__ import annotations

import sys


def main(argv=None):
    from simulator.__main__ import main as simulator
    return simulator(["check"] if argv is None else argv)


if __name__ == "__main__":
    sys.exit(main())
