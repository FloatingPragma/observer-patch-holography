"""Run the A-contrast-memory-spherical-v1 model and compare its results with the stored result file.

Arguments are passed to ``model.contrast_memory_spherical_v1.verify``: ``--full`` recomputes every
result, ``--output PATH`` writes a fresh result file, ``--verify PATH`` compares with
an explicit result file.
"""
from __future__ import annotations

import sys


def main(argv=None):
    from .contrast_memory_spherical_v1.verify import main as verify
    return verify(sys.argv[1:] if argv is None else argv)


if __name__ == '__main__':
    main()
