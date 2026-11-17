"""Low-level finite-group tables with no simulator package dependencies.

This module deliberately lives above ``defects`` and ``bulk`` in the import
graph.  Gauge-covariant overlap is a core dynamics primitive; importing it must
not initialize the defect-analysis or bulk-extraction packages.
"""

from __future__ import annotations

from itertools import permutations

import numpy as np


S3_ELEMENTS: tuple[tuple[int, int, int], ...] = tuple(permutations((0, 1, 2)))
S3_INDEX = {element: index for index, element in enumerate(S3_ELEMENTS)}


def _compose(left: tuple[int, int, int], right: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(left[right[index]] for index in range(3))


S3_MUL = np.array(
    [[S3_INDEX[_compose(left, right)] for right in S3_ELEMENTS] for left in S3_ELEMENTS],
    dtype=np.int16,
)
S3_INV = np.array(
    [S3_INDEX[tuple(element.index(index) for index in range(3))] for element in S3_ELEMENTS],
    dtype=np.int16,
)
S3_CLASS = np.array(
    [
        0
        if element == (0, 1, 2)
        else 1
        if sum(1 for index, image in enumerate(element) if index != image) == 2
        else 2
        for element in S3_ELEMENTS
    ],
    dtype=np.int16,
)


__all__ = ["S3_CLASS", "S3_ELEMENTS", "S3_INDEX", "S3_INV", "S3_MUL"]
