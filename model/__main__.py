"""Check the OPH law.

Runs the exactness and regression checks of the simulator that executes the OPH law
and the controls of its experiments (``python -m simulator check``).
"""
from __future__ import annotations

import sys


def main(argv=None):
    from simulator.__main__ import main as simulator
    return simulator(["check"] if argv is None else argv)


if __name__ == "__main__":
    sys.exit(main())
