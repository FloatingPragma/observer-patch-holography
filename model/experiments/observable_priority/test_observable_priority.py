"""Regression checks for rank-blind recovery, physical separation and certified replay."""

from ast import literal_eval
from pathlib import Path
import re
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from .__main__ import compare_candidates
from .cells import cell_inputs
from .closure import certificate, validate
from .equipment import install
from .observed import Screen, reconstruct
from .priority import authenticate, observe, recover_document
from .source import ROOT, board_hash, canonical, read


class PriorityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.screen, cls.native = install(0, "native")
        cls.document = observe(0, cls.screen, cls.native)
        cls.relation = authenticate(0, cls.screen, cls.native, cls.document)

    def compare(self, left, right, **kwargs):
        options = dict(level=0, cell=0, left=left, right=right, left_hits=None, right_hits=None)
        return compare_candidates(SimpleNamespace(**(options | kwargs)))

    def test_recovery_accepts_any_row_order(self):
        shuffled = dict(self.document, rows=self.document["rows"][::-1])
        recovered = recover_document(0, self.screen, shuffled)
        self.assertEqual(recovered, self.relation)

    def test_installed_lean_charts_match_executable_federations(self):
        lean = ROOT / "Lean/ObserverPatchHolography/ObservablePriority"

        def literal(path, name):
            # Read the whole constant so an appended or computed entry cannot
            # silently disappear from the comparison. Rows contain nested tuples.
            match = re.search(rf"^def {name}\b[^\n]*?:=\s*(.*?)(?=^(?:def|theorem|end)\b|\Z)",
                              path.read_text(encoding="utf-8"), re.MULTILINE | re.DOTALL)
            self.assertIsNotNone(match, f"{path.name} is missing literal {name}")
            return literal_eval(match[1].strip().translate(str.maketrans("⟨⟩", "[]")))

        chart = literal(lean / "PriorityCellData.lean", "originalChart")
        for level in range(3):
            screen, service = install(level, "native")
            cells = cell_inputs(screen, range(len(screen.edges)))
            blocks = sorted(lean.glob(f"Screen{len(cells)}Block*.lean"),
                            key=lambda path: int(path.stem.split("Block")[1]))
            rows = [row for block in blocks for row in literal(block, "rows")]
            self.assertEqual(len(rows), len(cells))
            for cell, (inner, boundary) in enumerate(cells):
                with self.subTest(level=level, cell=cell):
                    inner = sorted(inner, key=service.rank.__getitem__)
                    self.assertEqual(chart, [(screen.edges[e][0] % 12, screen.edges[e][1] % 12)
                                             for e in inner])
                    self.assertEqual([service.rank[e] for e in inner], list(range(1, len(chart) + 1)))
                    self.assertEqual([service.rank[e] for _, _, e in boundary], [0, 0, 0])
                    slots = [u % 12 for u, _, _ in boundary]
                    remote = [v for _, v, _ in boundary]
                    occupied = set(range(12 * cell, 12 * cell + 12)) | set(remote)
                    spare = [v for v in range(screen.n) if v not in occupied][:3]
                    self.assertEqual(rows[cell], [cell, slots, remote, spare])

    def test_bad_wedges(self):
        screen = Screen(4, [(0, 1), (0, 2), (0, 3)])
        good = [[0, 1, 1, 0], [0, 2, 1, 0], [1, 2, 1, 0]]
        for rows, error in (
            (good[1:], "missing wedge"),
            (good + [good[0]], "duplicate wedge"),
            ([[0, 1, 1, 1], *good[1:]], "not exactly one"),
            ([[0, 1, 0, 0], *good[1:]], "not exactly one"),
            ([[0, 1, True, 0], *good[1:]], "invalid transcript row"),
            ([[1, 0, 0, 1], *good[1:]], "nonincident or noncanonical"),
            ([good[0], [0, 2, 0, 1], good[2]], "cyclic"),
        ):
            with self.subTest(error=error), self.assertRaisesRegex(ValueError, error):
                reconstruct(screen, rows)

    def test_foreign_board_or_revision(self):
        for key, value in (("board", "0" * 64), ("pin", "0" * 40), ("level", False), ("format", True)):
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "transcript binding"):
                recover_document(0, self.screen, self.document | {key: value})

    def test_reject_malformed_ranks_and_transcript(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.json"
            for data, error in (([0] * len(self.screen.edges), "incident priority tie"),
                                ([0], "invalid candidate ranks"),
                                ([True] * len(self.screen.edges), "invalid candidate ranks")):
                path.write_text(canonical(data))
                with self.subTest(error=error), self.assertRaisesRegex(ValueError, error):
                    install(0, path)
            path.write_text('{"rows":[],"rows":[]}')
            with self.assertRaisesRegex(ValueError, "duplicate JSON field"):
                read(path)

    def test_coherent_observations_authenticate_only_at_their_own_service(self):
        _, reverse = install(0, "reverse")
        other = observe(0, self.screen, reverse)
        self.assertNotEqual(recover_document(0, self.screen, other), self.relation)
        with self.assertRaisesRegex(ValueError, "false addressed read"):
            authenticate(0, self.screen, self.native, other)

    def test_distinguishing_preparation(self):
        result = self.compare("native", "reverse")
        self.assertFalse(result["equivalent"])
        self.assertEqual(result["experiment"],
                         dict(pair=[0, 1], double=0, vacancies=[1, 5], default_load=1,
                              initial_collars="blank", read=0))
        self.assertEqual(result["addressed_hits"], [[1, 0], [0, 1]])

    def test_cell_reads_the_certified_patch(self):
        result = self.compare("native", "relabeled")
        self.assertTrue(result["equivalent"])
        self.assertEqual(result["classical_bits"], 3)
        self.assertEqual(len(result["certificate"]["order"]), 33)
        self.assertEqual(result["boundary_addresses"], [9, 19, 32])
        self.assertEqual({tuple(w["boundary_hits"]) for w in result["witnesses"]},
                         {(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)})
        self.assertEqual(len({tuple(w["internal_word"]) for w in result["witnesses"]}), 8)

    def test_disjoint_order_inversion(self):
        n = len(self.screen.edges)
        up = [r * (n + 1) + i for i, r in enumerate(self.native.rank)]
        down = [r * (n + 1) - i for i, r in enumerate(self.native.rank)]
        groups = {}
        for i, rank in enumerate(self.native.rank):
            groups.setdefault(rank, []).append(i)
        i, j = next(row[:2] for row in groups.values() if len(row) > 1)
        self.assertTrue(set(self.screen.edges[i]).isdisjoint(self.screen.edges[j]))
        self.assertTrue(up[i] < up[j] and down[i] > down[j])
        with TemporaryDirectory() as tmp:
            paths = [Path(tmp) / name for name in ("up.json", "down.json")]
            for path, ranks in zip(paths, (up, down)):
                path.write_text(canonical(ranks))
            result = self.compare(*paths)
        self.assertTrue(result["equivalent"] and result["native_cell_class"])

    def test_equivalent_reversed_candidates(self):
        result = self.compare("reverse", "reverse")
        self.assertTrue(result["equivalent"])
        self.assertFalse(result["native_cell_class"])
        self.assertNotIn("classical_bits", result)

    def test_unjustified_or_incomplete_closure(self):
        local, boundary = cell_inputs(self.screen, self.relation.order)[0]
        seeds = set(local)
        board = board_hash(0)
        cert = certificate(self.relation, board, seeds)
        removed = [i for i in cert["order"] if i != boundary[0][2]]
        extra = next(i for i in self.relation.order if i not in cert["order"])
        for order, reason in ((removed, "missing or unordered"), (cert["order"][::-1], "missing or unordered"),
                              ([extra, *cert["order"]], "incorrect closure"),
                              (cert["order"] + [True], "invalid address")):
            with self.subTest(reason=reason), self.assertRaisesRegex(ValueError, reason):
                validate(self.relation, board, seeds, dict(cert, order=order))


if __name__ == "__main__":
    unittest.main()
