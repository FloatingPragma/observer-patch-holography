"""Native modal record factors on the two finite observer-pair carriers."""

import argparse
import ast
from collections import Counter
from itertools import product
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from ...signed_effects.full_address.measure import local_check
from ...signed_effects.fixed_pointer import copy_word

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
FORMAL = ROOT / "Lean/QFT/SourceOperatorGeneration.lean"
IDS = (86, 88, 247)
SHOTS = 65536
ZERO = (0, 0, 0, 0, 0)


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def formal_source():
    s = FORMAL.read_text()
    result = {}
    for observer in IDS:
        region = s.split(f"def obs{observer} :", 1)[1].split("  cnt :=", 1)[0]
        table, path = region.split("  path := ")
        table = table.split("fun i => ", 1)[1].strip().removesuffix(" i")
        result[observer] = (
            list(map(tuple, ast.literal_eval(table.replace("!", "")))),
            ast.literal_eval(path.strip().removeprefix("!")),
        )
    return result


def packet_bits(label):
    family, checkpoint, stable, sector, load = label
    require(
        family == sector == 0
        and 0 <= checkpoint < 4
        and stable in (0, 1)
        and 0 <= load < 8,
        "packet outside the literal vocabulary",
    )
    return (
        tuple((checkpoint >> b) & 1 for b in range(2))
        + (stable,)
        + tuple((load >> b) & 1 for b in range(3))
    )


def from_ports(ports):
    pairs = np.asarray(ports).reshape(-1, 6, 2)
    left, right = pairs[:, :, 0], pairs[:, :, 1]
    require(((right == 0) | (right == 2)).all(), "invalid pointer reading")
    bits = right // 2
    sums = left + right
    # The two constant fields have separate measurable balance witnesses.
    return np.column_stack(
        (
            abs(sums[:, 0] - 2),
            bits[:, 0] + 2 * bits[:, 1],
            bits[:, 2],
            abs(sums[:, 3] - 2),
            bits[:, 3] + 2 * bits[:, 4] + 4 * bits[:, 5],
        )
    )


def owner_layout():
    board = federation(2)
    owner_graph = nx.Graph([(a // 12, b // 12) for a, b in board.external])
    supports = {
        int(k): set(v)
        for k, v in json.loads((HERE / "layout.json").read_text())["observers"].items()
    }
    occurrence = Counter(x for s in supports.values() for x in s)
    shared = sorted(x for x in occurrence if occurrence[x] > 1)
    seeds = (0, 208, 89)
    distances = {
        s: nx.single_source_shortest_path_length(owner_graph, s) for s in seeds
    }
    zones, mapping, work = {}, {}, {}
    for observer, seed in zip(IDS, (208, 89, 0), strict=True):
        pool = [
            c
            for c in sorted(owner_graph)
            if min(seeds, key=lambda s: distances[s][c]) == seed
        ]
        private = sorted(s for s in supports[observer] if occurrence[s] == 1)
        tree = nx.bfs_tree(owner_graph.subgraph(pool), seed, sort_neighbors=sorted)
        owners = list(tree)[: len(private) + 3]
        require(len(owners) == len(private) + 3, "connected private zone too small")
        # Put the four-port gate on native ports 0,6,10,7 of an eligible owner.
        eligible = [
            c for c in owners if all(board.degree[12 * c + j] == 5 for j in (0, 10))
        ]
        w = min(eligible)
        others = [c for c in owners if c != w]
        zones[observer] = owners
        work[observer] = (w, others[:2])
        mapping.update(zip(private, others[2:], strict=True))
    spare = sorted(set(owner_graph) - {c for z in zones.values() for c in z})
    mapping.update(zip(shared, spare[: len(shared)], strict=True))
    require(len(set(mapping.values())) == len(occurrence), "two source nodes alias")
    return board, supports, shared, zones, mapping, work


def setup():
    board, supports, shared, zones, mapping, work = owner_layout()
    form = formal_source()
    local_check()  # all 81 microscopic outputs; one common phase for all six inputs
    full = nx.Graph(board.seams)
    programs = {}
    for o in IDS:
        labels, _ = form[o]
        w, register_owners = work[o]
        private = sorted(supports[o] - set(shared))
        graph = full.subgraph([12 * c + k for c in zones[o] for k in range(12)])
        require(nx.is_connected(graph), "inaccessible factor address")
        base = 12 * w
        for kind, a, b, _ in copy_word(0):
            require(
                (base + a, base + b) in graph.edges
                if kind == "contact"
                else board.degree[base + a] == 5,
                "copy gate unavailable on board",
            )
        rails = [12 * c + k for c in register_owners for k in range(12)][: len(labels)]
        routes = {
            key: nx.single_source_shortest_path(graph, base + key) for key in (0, 6, 10)
        }
        occupied_bits = [
            np.flatnonzero(packet_bits(label)).tolist() for label in labels
        ]
        # A transposition along a path of l edges takes 2l-1 pi contacts.
        # Three disjoint endpoint swaps take inputs to the workbench; repeat
        # the same transpositions after the local copy. Intermediate nodes restore.
        costs = []
        for i, bs in enumerate(occupied_bits):
            for node in private:
                for b in bs:
                    ports = (
                        rails[i],
                        12 * mapping[node] + 2 * b,
                        12 * mapping[node] + 2 * b + 1,
                    )
                    lengths = [
                        len(routes[r][p]) - 1
                        for r, p in zip((0, 6, 10), ports, strict=True)
                    ]
                    require(min(lengths) >= 1, "workspace overlaps packet or address")
                    costs.append(2 * sum(2 * length - 1 for length in lengths))
        programs[o] = dict(
            labels=labels,
            private=private,
            supports=sorted(supports[o]),
            owners=zones[o],
            workbench=w,
            address=rails,
            calls=len(costs),
            routing_contacts=sum(costs),
        )
    return programs, shared, mapping


def encode(program, index):
    # Exact action of the certified material fanout on all label columns.
    # A data port is 2 for its addressed column and 1 in every other column.
    addresses = np.ones(len(program["labels"]), dtype=int)
    addresses[index] = 2
    bits = np.zeros((len(program["supports"]), 6), dtype=int)
    private_rows = [
        i for i, node in enumerate(program["supports"]) if node in program["private"]
    ]
    for label, raw in zip(program["labels"], addresses, strict=True):
        if raw != 1:
            bits[private_rows] ^= packet_bits(label)
    return np.stack((2 - 2 * bits, 2 * bits), axis=-1)


def mode(packets):
    # np.unique sorts field values, so argmax uses the historical tie rule.
    out = []
    for column in packets.T:
        values, counts = np.unique(column, return_counts=True)
        out.append(int(values[np.argmax(counts)]))
    return tuple(out)


def read(program, index):
    return mode(from_ports(encode(program, index)))


def native_halves():
    exchange, contrast = candidate.native_law.local_term_parts()
    src = candidate.native_law.SRC
    selected = [src.index(x) for x in ((2, 1), (1, 2))]
    kernels = {}
    for name, h in (
        ("native", exchange + candidate.KAPPA * contrast),
        ("exchange_off", candidate.KAPPA * contrast),
    ):
        block = h[np.ix_(selected, selected)]
        np.testing.assert_allclose(
            np.delete(h[:, selected], selected, axis=0), 0, atol=1e-13
        )
        k = expm(-1j * np.pi / 2 * (block - block[0, 0] * np.eye(2)))
        kernels[name] = k
        if name == "exchange_off":
            np.testing.assert_allclose(
                h[:9, :9], np.diag(np.diag(h[:9, :9])), atol=1e-13
            )
    return kernels


def interference(rng, programs):
    kernels = native_halves()
    hist = []
    for o in IDS:
        p = programs[o]
        for arm in ("coherent", "published", "exchange_off"):
            gate = kernels["exchange_off" if arm == "exchange_off" else "native"]
            state = gate[:, 0]
            density = np.outer(state, state.conj())
            if arm == "published":
                density = np.diag(np.diag(density))
            probs = np.diag(gate @ density @ gate.conj().T).real
            exact = np.rint(2 * probs) / 2
            np.testing.assert_allclose(exact, probs, rtol=0, atol=2e-13)
            draws = rng.multinomial(SHOTS, exact)
            counts = [0] * (len(p["labels"]) + 1)
            for index, count in zip((0, len(p["labels"]) - 1), draws, strict=True):
                # With exchange removed both routing and fanout are diagonal.
                packet = ZERO if arm == "exchange_off" else read(p, index)
                idx = p["labels"].index(packet) if packet in p["labels"] else -1
                counts[idx] += int(count)
            hist.append(dict(observer=o, arm=arm, counts=counts))
    for left, right in ((86, 88), (86, 247)):
        dims = (len(programs[left]["labels"]) + 1, len(programs[right]["labels"]) + 1)
        for publish in (False, True):
            g = kernels["native"]
            both = np.kron(g, g)
            rho = np.outer(both[:, 0], both[:, 0].conj()).reshape(2, 2, 2, 2)
            if publish:
                rho[0, :, 1, :] = 0
                rho[1, :, 0, :] = 0
            observed = np.diag(both @ rho.reshape(4, 4) @ both.conj().T).real
            exact = np.rint(2 * observed) / 2
            np.testing.assert_allclose(exact, observed, atol=2e-13, rtol=0)
            counts = np.zeros(dims, dtype=int)
            for (i, j), num in zip(
                product((0, dims[0] - 2), (0, dims[1] - 2)),
                rng.multinomial(SHOTS, exact),
                strict=True,
            ):
                counts[i, j] = num
            other = next(o for o in IDS if o not in (left, right))
            middle = len(programs[other]["labels"]) // 2
            require(
                read(programs[other], middle) == programs[other]["labels"][middle],
                "spectator read changes",
            )
            hist.append(
                dict(
                    pair=[left, right],
                    publish_left=publish,
                    counts=counts.tolist(),
                    spectator=[other, middle],
                )
            )
    return hist


def shift_tags(rng, dimension):
    pending = np.arange(SHOTS)
    tags = np.empty(SHOTS, dtype=int)
    refusals = []
    while len(pending):
        draws = rng.integers(16, size=len(pending))
        accepted = draws < dimension
        tags[pending[accepted]] = draws[accepted]
        pending = pending[~accepted]
        refusals.append(len(pending))
    return tags, refusals


def regional_erasure(programs):
    rng = np.random.default_rng(413329)
    form = formal_source()
    fair = abs(native_halves()["native"][:, 0]) ** 2
    np.testing.assert_allclose(fair, [0.5, 0.5], atol=2e-13, rtol=0)
    banks = []
    for right in (88, 247):
        sides = (86, right)
        modes = ["source", "left", "right", "both"] + (
            ["shared"] if right == 247 else []
        )
        for arm in modes:
            starts = rng.integers(32, size=SHOTS)
            records = [np.array(form[o][1])[starts] for o in sides]
            refusals = {}
            common = None
            for side, o in enumerate(sides):
                which = ("left", "right")[side]
                if arm in ("both", "shared", which):
                    if arm == "shared" and side == 1:
                        draws = common
                    else:
                        draws, retries = shift_tags(rng, len(form[o][0]))
                        refusals[which] = retries
                        if arm == "shared":
                            common = draws
                    indices = (records[side] + draws) % len(form[o][0])
                    reads = [
                        programs[o]["labels"].index(read(programs[o], k))
                        for k in range(len(form[o][0]))
                    ]
                    records[side] = np.array(reads)[indices]
            dims = (len(form[86][0]), len(form[right][0]))
            ids = np.ravel_multi_index(records, dims)
            counts = np.bincount(ids, minlength=dims[0] * dims[1]).reshape(dims)
            banks.append(
                dict(
                    right=right,
                    arm=arm,
                    source_times=np.bincount(starts, minlength=32).tolist(),
                    refusals=refusals,
                    joint=counts.tolist(),
                )
            )
    return banks


def report():
    programs, shared, mapping = setup()
    form = formal_source()
    histories, layout, power = {}, {}, {}
    for o, p in programs.items():
        labels = p["labels"]
        reads = [read(p, i) for i in range(len(labels))]
        require(reads == labels, "wrong modal packet")
        histories[o] = [reads.index(reads[i]) for i in form[o][1]]
        layout[o] = {
            k: v for k, v in p.items() if k not in ("private", "supports", "labels")
        }
        layout[o]["support_size"] = len(p["supports"])
        layout[o]["private_size"] = len(p["private"])
        # Prepare the largest label, overwrite private votes to keep k of them,
        # and switch all shared packets between the baseline and the last label.
        private = [i for i, n in enumerate(p["supports"]) if n in p["private"]]
        common = [i for i, n in enumerate(p["supports"]) if n in shared]
        last = labels[-1]
        panels = []
        for kept in (len(private), 48, 47, 49 - len(common), 48 - len(common)):
            if not 0 <= kept <= len(private):
                continue
            data = from_ports(encode(p, len(labels) - 1))
            data[private[kept:]] = ZERO
            first = mode(data)
            data[common] = last
            panels.append(dict(kept=kept, zero_overlap=first, last_overlap=mode(data)))
        power[o] = panels
    triple = list(product(*(range(len(form[o][0])) for o in IDS)))
    products = {}
    for r in (88, 247):
        dim = len(form[86][0]) * len(form[r][0])
        indices = (IDS.index(86), IDS.index(r))
        require(
            len({(t[indices[0]], t[indices[1]]) for t in triple}) == dim,
            "carrier alias",
        )
        walk = list(zip(histories[86], histories[r], strict=True))
        changes = Counter((a != c, b != d) for (a, b), (c, d) in zip(walk, walk[1:]))
        products[r] = dict(
            dimension=dim,
            units=dim**2,
            joint_walk=walk,
            steps=[changes[key] for key in product((False, True), repeat=2)],
        )
    return dict(
        placement=layout,
        source_to_cell=sorted(mapping.items()),
        shared_nodes=shared,
        triple_dimension=len(triple),
        pairs=products,
        quorum_controls=power,
        shots=SHOTS,
        interferometers=interference(np.random.default_rng(329344), programs),
        expectations=regional_erasure(programs),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = json.dumps(report(), separators=(",", ":"), allow_nan=False) + "\n"
    target = HERE / "records.json"
    if args.write:
        target.write_text(out)
    elif target.read_text() != out:
        raise SystemExit("Modal-slot records differ")
