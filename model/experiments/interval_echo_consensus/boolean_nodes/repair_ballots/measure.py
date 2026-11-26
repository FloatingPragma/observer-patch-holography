"""Read physical overlap ballots and service the first eligible owner."""

import argparse
import hashlib
import json
from itertools import product
from pathlib import Path

import numpy as np

from .. import material

HERE = Path(__file__).resolve().parent


def binary(bits):
    return "".join(str(int(x)) for x in bits)


class Apparatus:
    def __init__(self, n, edges=None):
        self.device = material.Machine(0, n)
        board = self.device.car
        physical = sorted(
            {tuple(sorted((a // 12, b // 12))) for a, b in board.external}
        )
        self.edges = physical if edges is None else edges
        if not all(i < j < n and (i, j) in physical for i, j in self.edges):
            raise ValueError("comparison requires a native owner edge")
        self.neighbors = [[] for _ in range(n)]
        for i, j in self.edges:
            self.neighbors[i].append(j)
            self.neighbors[j].append(i)
        self.routes = {}

    def probe(self, loads, u, v, at, mode="native"):
        """A strict echo, balanced pointer read, inverse echo and exact route undo."""
        edge = self.device.rows[at][3]
        a, b = edge
        key = u, v, at
        if key not in self.routes:
            self.routes[key] = self.device.move(u, a) + self.device.move(
                u if v == a else v, b
            )
        route = self.routes[key]
        for i, j in route:
            loads[i], loads[j] = loads[j], loads[i]
        click = material.balanced_pointer(mode)[loads[a], loads[b]]
        for i, j in reversed(route):
            loads[i], loads[j] = loads[j], loads[i]
        self.probes += 1
        self.exchanges += 2 * len(route)
        return click

    def flip(self, loads, pair):
        a, b = pair
        loads[a], loads[b] = loads[b], loads[a]
        self.exchanges += 1

    def read(self, loads, site, slot, at):
        rail = self.device.rows[site][slot][0]
        standard = self.device.rows[at][2][0]
        return 1 - self.probe(loads, rail, standard, at)

    def local_differences(self, loads, i, mode):
        output = self.device.rows[i][0][0]
        return [
            self.probe(loads, output, self.device.rows[j][0][0], i, mode)
            for j in self.neighbors[i]
        ]

    def run(self, bits, direction, mode):
        n = len(bits)
        self.probes = self.exchanges = 0
        loads = self.device.fresh(bits)
        # Executed only by the independent audit reader.
        audit_before = [self.read(loads, i, 0, n - 1) for i in range(n)]
        initial_audit = self.probes, self.exchanges
        sensor = "parked" if mode == "parked" else "native"
        survey = [self.local_differences(loads, i, sensor) for i in range(n)]
        for i, row in enumerate(survey):
            if row and all(row):
                self.flip(loads, self.device.rows[i][1])
        # Each end measures the same edge independently. Count one endpoint
        # per edge so that the initial rank is in units of mismatches.
        edge_count = sum(survey[i][self.neighbors[i].index(j)] for i, j in self.edges)
        if any(
            survey[i][self.neighbors[i].index(j)]
            != survey[j][self.neighbors[j].index(i)]
            for i, j in self.edges
        ):
            raise ArithmeticError("endpoint reports disagree")
        # The parked arm uses the same independent budget as the native arm.
        audit_rank = sum(audit_before[i] != audit_before[j] for i, j in self.edges)
        if mode != "parked" and edge_count != audit_rank:
            raise ArithmeticError("rank monitor disagrees with reported bits")
        budget = (
            0 if mode == "short" else (audit_rank if mode == "parked" else edge_count)
        )
        attempts = []
        for turn in range(budget):
            order = [(turn + direction * k) % n for k in range(n)]
            scan = []
            chosen = order[0]
            if mode != "blind":
                for i in order:
                    flag = self.read(loads, i, 1, 0)
                    scan.append(flag)
                    if flag:
                        chosen = i
                        break
            observed = self.local_differences(loads, chosen, sensor)
            wrote = bool(observed and all(observed))
            if wrote:
                self.flip(loads, self.device.rows[chosen][0])
            withdraw = (
                [chosen]
                if mode == "stale"
                else ([chosen, *self.neighbors[chosen]] if wrote else [])
            )
            withdrawn = []
            for i in withdraw:
                if self.read(loads, i, 1, i):
                    self.flip(loads, self.device.rows[i][1])
                    withdrawn.append(i)
            attempts.append(
                [chosen, binary(scan), binary(observed), int(wrote), withdrawn]
            )
        work = [self.probes - initial_audit[0], self.exchanges - initial_audit[1]]
        # Independent actual terminal read of the outputs and latched flags.
        end = [self.read(loads, i, 0, n - 1) for i in range(n)]
        flags = [self.read(loads, i, 1, n - 1) for i in range(n)]
        final_differences = [
            self.local_differences(loads, i, "native") for i in range(n)
        ]
        return {
            "input": binary(audit_before),
            "direction": direction,
            "mode": mode,
            "local_survey": [binary(r) for r in survey],
            "budget": budget,
            "attempts": attempts,
            "output": binary(end),
            "latched": binary(flags),
            "final_survey": [binary(r) for r in final_differences],
            "service_probes_exchanges": work,
            "audit_probes_exchanges": [self.probes - work[0], self.exchanges - work[1]],
        }


def jam_census(n, edges):
    """Exhaust every sharp word in bounded chunks using the calibrated XOR law."""
    hist = np.zeros(len(edges) + 1, dtype=np.int64)
    competing_count = 0
    reaches_unanimity = 0
    forced_target = 0
    for start in range(0, 1 << n, 8192):
        words = np.arange(start, min(1 << n, start + 8192), dtype=np.uint32)
        all_different = [np.ones(len(words), dtype=bool) for _ in range(n)]
        can_zero = np.ones(len(words), dtype=bool)
        can_one = np.ones(len(words), dtype=bool)
        seen = set()
        broken = np.zeros(len(words), dtype=np.uint8)
        for i, j in edges:
            left, right = (words >> i) & 1, (words >> j) & 1
            different = left ^ right
            can_zero &= (left & right) == 0
            can_one &= (left | right) == 1
            all_different[i] &= different.astype(bool)
            all_different[j] &= different.astype(bool)
            seen.update((i, j))
            broken += different.astype(np.uint8)
        # An isolated owner can never change its initial report.
        for i in set(range(n)) - seen:
            can_zero &= ((words >> i) & 1) == 0
            can_one &= ((words >> i) & 1) == 1
        reaches_unanimity += int((can_zero | can_one).sum())
        active = np.zeros(len(words), dtype=bool)
        force_zero, force_one = can_zero.copy(), can_one.copy()
        for i in seen:
            active |= all_different[i]
            bit = (words >> i) & 1
            force_zero &= ~(all_different[i] & (bit == 0))
            force_one &= ~(all_different[i] & (bit == 1))
        forced_target += int((force_zero | force_one).sum())
        competing = np.zeros(len(words), dtype=bool)
        for i, j in edges:
            competing |= all_different[i] & all_different[j]
        competing_count += int(competing.sum())
        hist += np.bincount(broken[~active], minlength=len(hist))
    return (
        {str(i): int(v) for i, v in enumerate(hist) if v},
        competing_count,
        reaches_unanimity,
        forced_target,
    )


def report():
    graphs = {}
    for shape, n, edges in (
        ("path", 4, [(0, 1), (1, 2), (2, 3)]),
        ("pentagon", 5, [(0, 1), (1, 2), (2, 3), (3, 4), (0, 4)]),
        ("sphere", 20, None),
    ):
        app = Apparatus(n, edges)
        if shape == "sphere":
            u, v = app.edges[0]
            competing = [0] * n
            competing[v] = 1
            for neighbor in app.neighbors[u]:
                competing[neighbor] = 1
            rng = np.random.default_rng(780126)
            sources = [
                [0] * n,
                [1] * n,
                [1] + [0] * (n - 1),
                competing,
                *rng.integers(2, size=(8, n)).tolist(),
            ]
        else:
            sources = list(product(range(2), repeat=n))
        runs = [
            app.run(bits, direction, mode)
            for bits in sources
            for direction in (1, -1)
            for mode in ("native", "blind", "stale", "parked", "short")
        ]
        quiescent, competing, unanimous, forced = jam_census(n, app.edges)
        graphs[shape] = {
            "edges": app.edges,
            "quiescent_by_residual_mismatches": quiescent,
            "inputs_with_multiple_terminal_words": competing,
            "inputs_with_a_unanimous_terminal_word": unanimous,
            "inputs_forcing_a_fixed_unanimous_word": forced,
            "runs": runs,
        }
    from model.executable_model import carrier
    from model.candidate_laws.interval_echo import candidate

    root = HERE.parents[4]
    files = [
        HERE / "PROTOCOL.md",
        Path(__file__),
        Path(material.__file__),
        Path(candidate.__file__),
        Path(candidate.native_law.__file__),
        Path(carrier.__file__),
        Path(carrier.__file__).parent / "carriers" / "FEDERATION_LEVEL0.npz",
    ]
    return {
        "schema": "Native canonical adaptive ballots 1",
        "attempt_columns": [
            "site",
            "scanned_flags",
            "guard_differences",
            "wrote",
            "withdrawn_flags",
        ],
        "sources": {
            str(p.resolve().relative_to(root)): hashlib.sha256(
                p.read_bytes()
            ).hexdigest()
            for p in files
        },
        "graphs": graphs,
    }


def serialize(receipt):
    # One compact line per complete trial keeps long all-outcome banks readable.
    header = {k: v for k, v in receipt.items() if k != "graphs"}
    graphs = []
    for name, graph in receipt["graphs"].items():
        metadata = {k: v for k, v in graph.items() if k != "runs"}
        rows = ",\n".join(
            json.dumps(run, separators=(",", ":")) for run in graph["runs"]
        )
        graphs.append(
            json.dumps(name)
            + ": "
            + json.dumps(metadata)[:-1]
            + ', "runs": [\n'
            + rows
            + "]}"
        )
    return json.dumps(header)[:-1] + ',\n"graphs": {\n' + ",\n".join(graphs) + "}}\n"


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = report()
    receipt = HERE / "results.json"
    if args.write:
        receipt.write_text(serialize(data))
    if json.loads(receipt.read_bytes()) != json.loads(json.dumps(data)):
        raise ValueError("receipt differs from the native replay")
    for k, r in data["graphs"].items():
        print(k, len(r["runs"]), r["quiescent_by_residual_mismatches"])
