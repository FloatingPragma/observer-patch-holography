"""Read an independent board-seam midpoint distance for every eligible pair."""

import json
from pathlib import Path

from scipy.sparse.csgraph import shortest_path

from . import material as mat


def audit():
    result = []
    for lev in range(3):
        graph = mat.federation(lev).adjacency()
        ports = mat.route(lev)
        distances = shortest_path(graph, directed=False, unweighted=True, indices=ports)
        lengths = {}
        for length in mat.LENGTHS:
            positions = []
            for a in range(1, 16 - length):
                left = 2 * a - 1
                right = 2 * (a + length) - 1
                distance = (
                    int(distances[left : left + 2, ports[right : right + 2]].min()) + 1
                )
                positions.append(
                    {
                        "first": a,
                        "left_seam": ports[left : left + 2],
                        "right_seam": ports[right : right + 2],
                        "ambient": distance,
                    }
                )
            lengths[length] = {
                "material": 2 * length,
                "positions": positions,
                "middle_ambient": positions[(len(positions) - 1) // 2]["ambient"],
            }
        result.append(lengths)
    return result


if __name__ == "__main__":
    Path(__file__).with_name("ruler.json").write_text(
        json.dumps(audit(), indent=1) + "\n"
    )
