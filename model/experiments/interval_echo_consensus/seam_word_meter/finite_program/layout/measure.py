"""Compact graph certificates and held sharp material programs."""

import argparse
import json
from itertools import combinations
from pathlib import Path
from types import SimpleNamespace

from model.executable_model import carrier
from model.executable_model.refinement import load_refinement
from ..machine import decode_word
from .mount import Mount, PORTS, TEMPLATE, local_path, station

HERE = Path(__file__).resolve().parent


def all_cells():
    cases = 0
    worst = 0
    for writer in range(12):
        for others in combinations(sorted(PORTS - {writer}), 2):
            gates = (writer, *others)
            head, refs, matching = station(writer, gates)
            frozen = {writer, head, *refs}
            remaining = PORTS - frozen
            assert set().union(*map(set, matching)) == remaining and len(matching) == 4
            assert head not in gates and all(
                r in TEMPLATE[head] and r not in gates for r in refs
            )
            anchor = min(remaining)
            for end in range(12):
                path = local_path(anchor, end, frozen - {anchor, end})
                assert not frozen.intersection(path[1:-1])
                worst = max(worst, len(path) - 1)
            cases += 1
    return {"labelled_writer_gateway_cases": cases, "largest_checked_local_path": worst}


def check_path(mount, source, dest):
    from ...chart import exchange_word

    path = mount.path(source, dest)
    assert len(set(path)) == len(path)
    assert set(path[1:-1]).isdisjoint(mount.protected)
    word = exchange_word(path)
    assert all(mount.has_edge(a, b) and {a, b} != set(mount.writer) for a, b in word)
    # Distinct symbolic labels include any intervening data dimers or flags.
    labels = {port: port for port in path}
    for a, b in word:
        labels[a], labels[b] = labels[b], labels[a]
    assert all(
        labels[p] == (dest if p == source else source if p == dest else p) for p in path
    )
    return dict(source=source, dest=dest, path=list(path), swaps=len(word))


def check_board(board):
    m = Mount(board)
    used = bytearray(12 * m.cells)
    for port in m.protected:
        used[port] = 1
    for j in range(m.capacity + 2):
        a, b = m.dimer(j)
        assert not used[a] and not used[b] and m.has_edge(a, b)
        used[a] = used[b] = 1
    assert all(used)
    assert m.capacity == 6 * m.cells - 6
    # These routes include both read stations, both write endpoints and the
    # last data pair on the whole board. The graph argument covers every pair.
    ships = [*m.writer, *m.dimer(m.capacity + 1)]
    sample = [
        check_path(m, u, v) for u in ships for v in (*m.writer, *m.heads) if u != v
    ]
    bench = m.bench(registers=m.capacity, axes=1, width=1)
    assert len(bench.words) == m.capacity and bench.words[-1][0] == m.dimer(
        m.capacity + 1
    )
    try:
        m.bench(registers=m.capacity + 1, axes=1, width=1)
    except ValueError as ex:
        assert str(ex) == "Insufficient data dimers"
    else:
        raise AssertionError("Accepted one excess data bit")
    return dict(
        cells=m.cells,
        capacity=m.capacity,
        writer=m.writer,
        heads=m.heads,
        references=m.references,
        flags=m.flags,
        owner_tree_edges=m.cells - 1,
        maximum_tree_depth=int(m.depth.max()),
        example_routes=sample,
    )


def program(level, registers, axes, width):
    bench = Mount(carrier.federation(level)).bench(
        registers=registers, axes=axes, width=width
    )
    count = len(bench.words)
    inputs = [(5 * j + 3) % (1 << width) for j in range(count)]
    initial = [(1 << (width - 1)) - 2, 0]
    distant = [(1 << width) - 2, 4]
    for i in (0, 1):
        inputs[i], inputs[count - axes + i] = initial[i], distant[i]
    inputs[2] = (1 << width) - 1
    bench.start(inputs)
    raw_before = bench.loads.copy()
    pairs = [(i, count - axes + i) for i in (0, 1)]
    assert bench.mean_axes(pairs) == "accepted"
    expected = [(a + b) // 2 for a, b in zip(initial, distant)]
    observations = [
        [bench.read(j, station) for j in (0, 1, count - axes, count - axes + 1, 2)]
        for station in (0, 1)
    ]
    assert observations == [[*expected, *expected, inputs[2]]] * 2
    # Each other physical load, including all standards, returned unaltered.
    changed = {
        v for pair in pairs for r in pair for edge in bench.words[r] for v in edge
    }
    assert all(
        a == b
        for j, (a, b) in enumerate(zip(raw_before, bench.loads))
        if j not in changed
    )
    flags = bench.public("high")
    assert all(decode_word(pointer) == 2 for pointer in flags)
    # The high dirty sentinel is odd; the second axis now stores two.
    before = bench.loads.copy()
    status = bench.mean(1, 2)
    assert status == "off_grid" and before == bench.loads
    return dict(
        level=level,
        registers=registers,
        axes=axes,
        width=width,
        dimers=count * width,
        scalar_pairs=pairs,
        input_pairs=list(zip(initial, distant)),
        tapes=bench.tape,
        observations=observations,
        publication=flags,
        invalid_sum=status,
    )


def bundled():
    b = carrier.federation(0)
    disjoint = SimpleNamespace(
        cells=2 * b.cells,
        cell_seams=carrier.CELL_SEAMS,
        external=[*b.external, *[(x + b.n, y + b.n) for x, y in b.external]],
    )
    try:
        Mount(disjoint)
    except ValueError as ex:
        assert str(ex) == "Disconnected owner graph"
        control = str(ex)
    else:
        raise AssertionError("Disconnected fixture accepted")
    return dict(
        local=all_cells(),
        bundled={str(i): check_board(carrier.federation(i)) for i in range(3)},
        programs=[program(0, 19, 2, 3), program(1, 9, 3, 13)],
        disconnected_control=control,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--deposit", type=Path)
    parser.add_argument("--level", type=int, default=3)
    parser.add_argument("--write", action="store_true")
    opts = parser.parse_args()
    if opts.deposit:
        snapshot = load_refinement(opts.deposit, level=opts.level)
        value = dict(
            level=opts.level, sha256=snapshot.sha256, layout=check_board(snapshot)
        )
        out = HERE / f"level{opts.level}.json"
    else:
        value = bundled()
        out = HERE / "records.json"
    payload = json.dumps(value, separators=(",", ":")) + "\n"
    if opts.write:
        out.write_text(payload)
    elif opts.deposit and not out.exists():
        print(
            f"SHA-verified level {opts.level}, {snapshot.cells} owners, {value['layout']['capacity']} available data dimers"
        )
        raise SystemExit(0)
    elif json.loads(payload) != json.loads(out.read_text()):
        raise ValueError("Layout receipt differs")
    print(out.name, "passed")
