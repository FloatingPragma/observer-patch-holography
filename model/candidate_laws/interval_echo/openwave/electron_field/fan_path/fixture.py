"""Extend the archived four-face charge preparation along actual federation contacts."""

from contextlib import contextmanager
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import networkx as nx
import numpy as np
from scipy import sparse

FIELD = Path(__file__).resolve().parents[1]
ROOT = FIELD.parents[4]


@contextmanager
def archived_fixture():
    """Use the checked local ring, two-cell generator and doublet decoder."""
    manifest = dict(
        line.split("  ", 1)[::-1]
        for line in (FIELD / "records.sha256").read_text().splitlines()
    )
    with TemporaryDirectory() as work, ZipFile(FIELD / "records.zip") as archive:
        assert set(archive.namelist()) == set(manifest)
        for name, digest in manifest.items():
            assert sha256(archive.read(name)).hexdigest() == digest
        directory = Path(work)
        archive.extractall(directory)
        package = directory / "experiments" / "interval_echo_collective_spin"
        saved_path = list(sys.path)
        # Each of these archived, top-level imports also installs the next
        # experiment's path. Reusing one from an expired extraction skips
        # that installation and can silently mix two copies of the fixture.
        imports = (
            "search", "winding", "laplacian", "threshold", "run",
            "geometry", "contact", "current_algebra", "operators",
        )
        saved_modules = {
            name: sys.modules.pop(name) for name in imports if name in sys.modules
        }
        sys.path[:0] = [str(package), str(ROOT)]
        try:
            import search

            witness = json.loads((package / "results.json").read_text())["headline"][
                "best_worst_qubit_success"
            ]
            yield search, witness
        finally:
            sys.path[:] = saved_path
            for name in imports:
                sys.modules.pop(name, None)
            sys.modules.update(saved_modules)


def terminal_fan(board, hub, degrees=None):
    rim = nx.Graph(sorted(set(face) - {hub}) for face in board.faces() if hub in face)
    assert len(rim) == len(rim.edges) == 5
    paths = []
    for a, b in list(rim.edges):
        rim.remove_edge(a, b)
        path = nx.shortest_path(rim, a, b)
        paths.extend((path, path[::-1]))
        rim.add_edge(a, b)
    if degrees is not None:
        paths = [p for p in paths if [int(board.degree[v]) for v in p] == degrees]
    return (hub, *min(paths))


@lru_cache(maxsize=1)
def asymmetric_paths(federation):
    """Geometry-only census, including boards with no qualifying routes."""
    census = []
    for level in range(3):
        board = federation(level)
        external = {v: u for a, b in board.external for u, v in ((a, b), (b, a))}
        rim = {
            v: [p for p in board.neighbours[v] if p in external and p != external[v]]
            for v in external
        }
        assert all(len(options) <= 1 for options in rim.values())
        paths = [(a, b) for a, b in external.items() if rim[a]]
        counts = []
        for _ in (2, 3):
            paths = [
                (*p, v, external[v])
                for p in paths
                for v in rim[p[-1]]
                if external[v] // 12 not in {q // 12 for q in p}
            ]
            eligible = [p for p in paths if rim[p[-1]]]
            counts.append(len(eligible))
        census.append((board.cells, *counts))
        if eligible:
            path = min(eligible)
            assert path == (98, 241, 240, 1011, 1010, 1045)
            links = tuple(zip(path[::2], path[1::2], strict=True))
            hubs = (path[0], *path[1::2])
            cells = tuple(terminal_fan(board, h, [5, 5, 5, 5, 6]) for h in hubs)
            assert all(
                a in cells[i] and b == cells[i + 1][0] for i, (a, b) in enumerate(links)
            )
            return board, cells, links, tuple(census)
    raise ValueError("No asymmetric three-bridge route on the declared boards")


def physical_path(board, witness, n):
    actual = witness["actual_ports"]
    left, right = tuple(actual["left"]), tuple(actual["right"])
    assert (actual["cells"], left[0], right[0]) == (20, 14, 24)
    neighbours = {v: u for a, b in board.external for u, v in ((a, b), (b, a))}
    assert neighbours[left[0]] == right[0]
    assert [v for v in left[1:] if v in neighbours] == [22]
    assert [v for v in right[1:] if v in neighbours] == [34]
    far_left = terminal_fan(board, neighbours[22])
    far_right = terminal_fan(board, neighbours[34])
    paths = {
        2: ((left, right), ((14, 24),)),
        3: ((left, right, far_right), ((14, 24), (34, 116))),
        4: ((far_left, left, right, far_right), ((68, 22), (14, 24), (34, 116))),
    }
    cells, bridges = paths[n]
    assert len({v for fan in cells for v in fan}) == 6 * n
    assert len({v[0] // 12 for v in cells}) == n
    for (before, after), (a, b) in zip(zip(cells, cells[1:]), bridges, strict=True):
        assert a in before and b in after and neighbours[a] == b
    actual_faces = {frozenset(face) for face in board.faces()}
    for cell in cells:
        assert all(
            frozenset((cell[0], u, v)) in actual_faces
            for u, v in zip(cell[1:], cell[2:])
        )
    return cells, bridges


@dataclass
class PortWords:
    cell: int
    vertex: int
    words: tuple
    rows: slice

    @property
    def shape(self):
        return tuple(len(words) for words in self.words)


class FanPath:
    def __init__(self, archived, size, sign, asymmetric=False):
        search, witness = archived
        assert sign in (-1, 1)
        self.sign = sign
        self.size = size
        if asymmetric:
            self.board, cells, bridges, self.census = asymmetric_paths(
                search.federation
            )
            self.cells, self.bridges = cells[:size], bridges[: size - 1]
        else:
            self.board = search.federation(0)
            self.cells, self.bridges = physical_path(self.board, witness, size)
            self.census = ()
        self.field = search.fan(4)
        self.heads = tuple(witness["release"]["head_sectors"][0])
        if sign == -1:
            self.heads = tuple(
                sum(i in e for e in self.field.edges) - h
                for i, h in enumerate(self.heads)
            )
        self.vacuum, self.code = search.vacuum_code(self.field, self.heads)
        self.entry = [
            0,
            *(self.cells[i + 1].index(b) for i, (_, b) in enumerate(self.bridges)),
        ]
        self.exit = [
            *(self.cells[i].index(a) for i, (a, _) in enumerate(self.bridges)),
            0,
        ]
        self.words_at = [
            tuple(
                self.field.words[
                    tuple(
                        h + sign * (int(v == x) - int(v == 0))
                        for v, h in enumerate(self.heads)
                    )
                ]
            )
            for x in range(6)
        ]
        self.basis = []
        self.blocks = {}
        owners = []
        for c, cell in enumerate(self.cells):
            for x, port in enumerate(cell):
                words = tuple(self.words_at[ref] for ref in self.references(c, x))
                start = len(self.basis)
                self.basis.extend((port, w) for w in product(*words))
                self.blocks[c, x] = PortWords(
                    c, x, words, slice(start, len(self.basis))
                )
                owners.extend([c] * (len(self.basis) - start))
        self.owner = np.array(owners)
        self.ids = dict(zip(self.basis, range(len(self.basis)), strict=True))
        assert len(self.ids) == len(self.basis)
        self.parked, self.bridge_term = self.generators()
        self.h = self.parked + self.bridge_term
        if size == 2:
            original = search.JoinedSpace(self.field, (self.heads, self.heads), sign)
            degrees = [int(self.board.degree[p]) for c in self.cells for p in c]
            np.testing.assert_array_equal(
                self.h.toarray(), original.hamiltonian(degrees)
            )
            for c in (0, 1):
                rows = self.blocks[c, 0].rows
                np.testing.assert_array_equal(
                    original.portal[c].ravel(), np.arange(rows.start, rows.stop)
                )

    def references(self, charge_cell, vertex):
        labels = (
            list(self.exit[:charge_cell]) + [vertex] + self.entry[charge_cell + 1 :]
        )
        assert len(labels) == self.size
        return labels

    def generators(self):
        size = len(self.basis)
        local = []
        crossing = []
        potential = np.zeros(size)
        for (c, x), block in self.blocks.items():
            port = self.cells[c][x]
            for i in range(block.rows.start, block.rows.stop):
                _, w = self.basis[i]
                potential[i] = 8 * (6 - int(self.board.degree[port]))
                for (a, b), mask in self.field.bits.items():
                    if x not in (a, b) or ((x == a) == bool(w[c] & mask)) != (
                        self.sign == 1
                    ):
                        continue
                    target = self.cells[c][b if x == a else a]
                    j = self.ids[target, w[:c] + (w[c] ^ mask,) + w[c + 1 :]]
                    if i < j:
                        local.append((i, j, 4))
                for f, word in enumerate(w):
                    for mask, cw in zip(
                        self.field.masks, self.field.clockwise, strict=True
                    ):
                        if word & mask in (cw, cw ^ mask):
                            j = self.ids[port, w[:f] + (word ^ mask,) + w[f + 1 :]]
                            if i < j:
                                local.append((i, j, 0.5))
                if c < self.size - 1 and x == self.exit[c]:
                    crossing.append((i, self.ids[self.bridges[c][1], w], 4))

        def laplacian(edges):
            i, j, weight = (np.asarray(v) for v in zip(*edges, strict=True))
            matrix = sparse.coo_matrix(
                (-weight, (i, j)), shape=(size, size), dtype=float
            ).tocsr()
            matrix = matrix + matrix.T
            return matrix - sparse.diags(np.asarray(matrix.sum(axis=0)).ravel())

        return laplacian(local) + sparse.diags(potential), laplacian(crossing)

    def source(self, c=0):
        block = self.blocks[c, 0]
        initial = np.zeros((len(self.basis), 2))
        for phase in (0, 1):
            row = np.ones(1)
            for i, words in enumerate(block.words):
                term = (
                    self.code[:, phase]
                    if i == c
                    else np.ones(len(words)) / np.sqrt(len(words))
                )
                row = np.kron(row, term)
            initial[block.rows, phase] = row
        np.testing.assert_allclose(initial.T @ initial, np.eye(2), atol=2e-15)
        return initial

    def decode(self, columns, c=None):
        """Return vectors indexed by unobserved head records, output code, input."""
        if c is None:
            c = self.size - 1
        block = self.blocks[c, 0]
        amplitudes = columns[block.rows].reshape(*block.shape, columns.shape[1])
        output = np.tensordot(self.code.T, amplitudes, axes=(1, c))
        return np.moveaxis(output, 0, -2).reshape(-1, 2, columns.shape[1])
