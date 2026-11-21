"""The optional loader must reject an unpinned alternative without an external checkout."""

import numpy as np
import pytest

from .carrier import federation
from .refinement import load_refinement


def test_altered_board_rejected_even_with_same_cell_graph(tmp_path):
    seams = np.asarray(federation(0).external, np.int32)
    # Preserve every cell adjacency, but replace one physical gateway label.
    before = seams[0, 0]
    taken = set(seams.ravel())
    seams[0, 0] = next(
        p for p in range(before // 12 * 12, before // 12 * 12 + 12) if p not in taken
    )
    path = tmp_path / "FEDERATION_LEVEL0.npz"
    np.savez_compressed(path, level=0, cells=20, external=seams)
    with pytest.raises(ValueError, match="checksum"):
        load_refinement(path, level=0)
