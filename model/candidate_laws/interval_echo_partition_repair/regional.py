"""Finite local interface repairs and the observer's regional word."""
from dataclasses import dataclass
from fractions import Fraction
from typing import ClassVar, NamedTuple


class RegionRead(NamedTuple):
    word: tuple[int, ...]
    open_boundary: int
    volume: int


@dataclass(frozen=True)
class PartitionRepair:
    """An interface owns one legal two-end word; commands are two fair bits."""

    bond_words: ClassVar[tuple[tuple[int, int], ...]] = ((0, 0), (1, 1), (1, 2), (2, 2), (2, 1))

    @classmethod
    def advance_bond(cls, index: int, command: int) -> int:
        if index not in range(len(cls.bond_words)) or command not in range(4):
            raise ValueError("invalid interface state or two-bit command")
        return (index + (-1, 0, 0, 1)[command]) % len(cls.bond_words)

    @staticmethod
    def advance_volume(value: int, command: int) -> int:
        if value not in (0, 1) or command not in range(4):
            raise ValueError("invalid volume bit or two-bit command")
        if value == 0:
            return int(command == 3)
        return int(command >= 2)

    @classmethod
    def bond_kernel(cls):
        return tuple(tuple(Fraction(sum(cls.advance_bond(i, c) == j for c in range(4)), 4)
                           for j in range(len(cls.bond_words))) for i in range(len(cls.bond_words)))

    @classmethod
    def read_region(cls, edge_words, volume_words, edges, region):
        """The local endpoint alphabet and each cell's own visible volume bit.

        edge_words contains legal word indices in the fixed order of edges.
        Edge orientation and cell numbers specify the register addresses.
        Empty and whole regions are supported; a cut contributes only its
        inside endpoint. The two counters are independent digital apertures.
        """
        region = frozenset(region)
        if len(edge_words) != len(edges) or not region <= set(range(len(volume_words))):
            raise ValueError("unregistered edges or cells")
        if any(v not in (0, 1) for v in volume_words):
            raise ValueError("invalid volume word")
        visible = []
        area = 0
        for index, (left, right) in zip(edge_words, edges, strict=True):
            if index not in range(len(cls.bond_words)) or left == right:
                raise ValueError("invalid interface")
            if not 0 <= left < len(volume_words) or not 0 <= right < len(volume_words):
                raise ValueError("interface leaves the cell register")
            a, b = cls.bond_words[index]
            if left in region:
                visible.append(a)
            if right in region:
                visible.append(b)
            if (left in region) != (right in region):
                area += (a if left in region else b) != 0
        local_volumes = tuple(volume_words[node] for node in sorted(region))
        return RegionRead(tuple(visible)+local_volumes, area, sum(local_volumes))
