"""Run the polynomial-sv model and compare its results with the stored result file.

Arguments are passed to ``model.polynomial_sv.verify``: ``--full`` recomputes every
result, ``--output PATH`` writes a fresh result file, ``--verify PATH`` compares with
an explicit result file.
"""
from __future__ import annotations

import sys


def main(argv=None):
    from .polynomial_sv.verify import main as verify
    return verify(sys.argv[1:] if argv is None else argv)


if __name__ == '__main__':
    main()
