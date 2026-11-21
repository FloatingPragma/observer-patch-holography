"""Read an explicitly supplied Openwave spherical refinement at its deposited byte identity."""

from dataclasses import dataclass
from hashlib import sha256
from io import BytesIO
from pathlib import Path
from typing import ClassVar

import numpy as np

from .carrier import CELL_SEAMS, federation

DEPOSIT_REVISION = "b189afa2420912c02e2ab5cee9808490ee6dd4d9"
# SHA-256 of simulation/carriers/FEDERATION_LEVEL{level}.npz at DEPOSIT_REVISION.
DEPOSIT_SHA256 = (
    "bb937e80df03327ec1b2001145f5218733865f370abd51b003fbb7c991a87659",
    "39bae3eacd040f3ed7f8f2d59974327a94b7cb2d31452f20ebc16f19962ea324",
    "6e701ed57c17b83db5217e8a881155a4361a36a11b15ef76449ea165c3ba0988",
    "8c25ed8277cf50f33ad0f53c5c3df7e13e93df23b555d28d3f010deadd695d39",
    "4b14f8f52abc0ce856340e58dd632ee69002d04b7f8c5325c095fc81dfa0f27c",
    "fe58adce37f3f12b59fc61cdd1b219c57501c645ba5df5f3c80bf63302767dae",
    "215cdfd9d692dc43cde020bd2cb50e6a734dddbe45f3f66ce33482d6ebec7c0e",
    "578bc497a10c989ce558768f89f95cac99071900b9838926e380d06a425c8700",
)


@dataclass(frozen=True, eq=False)
class Refinement:
    cell_seams: ClassVar = CELL_SEAMS
    level: int
    cells: int
    sha256: str
    external: np.ndarray


def load_refinement(path: Path, *, level: int) -> Refinement:
    """Keep native port labels and seam order. Every cell uses the shared CELL_SEAMS incidence."""
    if type(level) is not int or not 0 <= level < len(DEPOSIT_SHA256):
        raise ValueError(
            "Expected an integer deposited refinement level from 0 through 7"
        )
    content = Path(path).read_bytes()
    signature = sha256(content).hexdigest()
    if signature != DEPOSIT_SHA256[level]:
        raise ValueError(
            f"Level {level} refinement checksum differs from the deposited board"
        )
    # Load from the same immutable bytes that passed the identity check.
    with np.load(BytesIO(content), allow_pickle=False) as stored:
        cells = int(stored["cells"])
        seams = stored["external"]
        if int(stored["level"]) != level or cells != 20 * 4**level:
            raise ValueError("Inconsistent deposited subdivision level or cell count")
    if (
        seams.shape != (3 * cells // 2, 2)
        or not np.issubdtype(seams.dtype, np.integer)
        or np.any(seams < 0)
        or np.any(seams >= 12 * cells)
    ):
        raise ValueError("Invalid external port addresses")
    owners = seams // 12
    if (
        np.any(owners[:, 0] == owners[:, 1])
        or np.unique(seams).size != seams.size
        or not np.all(np.bincount(owners.ravel(), minlength=cells) == 3)
    ):
        raise ValueError("External seams must match three distinct ports of every cell")
    if level <= 2:
        small = federation(level)
        endpoints = set(map(tuple, np.sort(seams, axis=1)))
        if small.cells != cells or set(small.external) != endpoints:
            raise ValueError(
                "The supplied refinement and shared carrier have different labelled seams"
            )
    # An immutable bytes backing also prevents callers from re-enabling writes.
    external = np.frombuffer(seams.tobytes(), dtype=seams.dtype).reshape(seams.shape)
    return Refinement(level=level, cells=cells, sha256=signature, external=external)
