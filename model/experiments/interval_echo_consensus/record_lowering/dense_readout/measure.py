"""Acquire complete deterministic receiver transcripts for two finite dense publication layers."""

import argparse
import gzip
import json
from math import pi
from pathlib import Path
import random

import numpy as np

from model.executable_model.carrier import federation
from ...seam_word_meter import apparatus as seam
from ..hardware import calibrated
from .placement import obstruction
from .apparatus import (
    SITES,
    clock_trial,
    compile_layer,
    mount,
    owner_obstruction,
    paths,
)

HERE = Path(__file__).resolve().parent
SEED = 2601174
SIZE = len(SITES)


def physical_columns():
    permutation = np.eye(9)[:, [seam.LAW.SRC.index((b, a)) for a, b in seam.LAW.SRC]]
    np.testing.assert_allclose(seam.pulse(pi), permutation, atol=seam.TOL, rtol=0)
    read, _ = calibrated()
    assert read == {0: 0, 1: 0, 2: 1}
    return read


def read_paths(plan, stores, parked=False, zero="zero"):
    """Interpret every port swap. A token labels a sharp input for all its choices."""
    initial = list(stores)
    tape = []
    for item in plan:
        path = item["path"]
        hops = [] if parked else list(zip(path, path[1:]))
        for a, b in hops:
            stores[a], stores[b] = stores[b], stores[a]
        # The full echo, projection, inverse on (rail,0) returns the original
        # loads. The lamp on a two-unit rail is positive exactly at load two.
        if stores[item["ref"]] != zero:
            raise ValueError("Receiving private standard has changed")
        tape.append(stores[path[-1]])
        for a, b in reversed(hops):
            stores[a], stores[b] = stores[b], stores[a]
    if initial != stores:
        raise ValueError("A source word, cache or transit port changed during a survey")
    return tape


def stores_at_boundary(board, layout):
    # Names of independent bits. Their complements form the balanced second rails.
    stores = [f"transit{p}" for p in range(board.n)]
    for i, (a, b, c, d, h, z) in enumerate(layout):
        stores[a], stores[b] = ("published", i), ("complement-published", i)
        stores[c], stores[d] = ("cached", i), ("complement-cached", i)
        stores[h], stores[z] = "zero", "zero"
    return stores


def appointments(plan):
    lasts = max(e["end"] for e in plan)
    return dict(
        reads=len(plan),
        primitive_slots=sum(e["end"] - e["start"] for e in plan),
        parallel_slots=lasts,
        swaps=sum(2 * (len(x["path"]) - 1) for x in plan),
        local_pointer_writes=len(plan),
        local_projections=len(plan),
        local_inverse_writes=len(plan),
    )


def decode(plan, tape, word, oldword):
    # Completed receiving projections, in the same order as the fixed calendar.
    values = {"published": word, "cached": oldword}
    mask = sum(
        ((values[tag[0]] >> tag[1]) & 1) << j
        for j, tag in enumerate(tape)
        if tag != "zero"
    )
    return format(mask, f"0{(len(plan) + 3) // 4}x")


def settings():
    rng = random.Random(SEED)
    prior, first, second = [rng.getrandbits(SIZE) for _ in range(3)]
    words = [(first, second)]
    words += [(first ^ (1 << bit), second) for bit in range(SIZE)]
    words += [(first, second ^ (1 << bit)) for bit in range(SIZE)]
    words += [(rng.getrandbits(SIZE), rng.getrandbits(SIZE)) for _ in range(8)]
    return prior, words


def transpose_tape(planes, count):
    # A data-blind serialisation. Output k contains only that launch's lamps.
    width = (len(planes) + 3) // 4
    return [
        format(sum(((row >> k) & 1) << j for j, row in enumerate(planes)), f"0{width}x")
        for k in range(count)
    ]


def sharp_replay(board, sites, routes, clean, stale, prior, words):
    """Use one bitplane per physical rail to run all sharp cohorts together."""
    count = len(words)
    ones = (1 << count) - 1
    load = [ones if p % 3 == 1 else 0 for p in range(board.n)]
    for i, (a, b, c, d, head, ref) in enumerate(sites):
        load[a] = ones if (prior >> i) & 1 else 0
        load[b] = ones ^ load[a]
        load[c], load[d], load[head], load[ref] = 0, ones, 0, 0
    # Read the original word locally. Every positive lamp opens the cache swap.
    for i, site in enumerate(sites):
        a, b, c, d, head, ref = site
        result = read_paths([dict(path=routes[i, i], ref=ref)], load, zero=0)[0]
        assert result == load[a]
        load[c], load[d] = result, ones ^ result
    after_copy = load.copy()
    public = {p for site in sites for p in site[:2]}
    layers = []
    for layer in range(2):
        act = []
        before_edit = load.copy()
        for i, (a, b, *_) in enumerate(sites):
            wanted = sum(((w[layer] >> i) & 1) << k for k, w in enumerate(words))
            mask = load[a] ^ wanted
            act.append(mask)
            swap = mask & (load[a] ^ load[b])
            load[a] ^= swap
            load[b] ^= swap
            assert load[a] == wanted and load[b] == (ones ^ wanted)
        assert all(
            load[p] == before_edit[p] == after_copy[p]
            for p in range(board.n)
            if p not in public
        )
        lamps = read_paths(clean, load, zero=0)
        stale_lamps = read_paths(stale, load, zero=0)
        parked_lamps = read_paths(clean, load, parked=True, zero=0)
        native = transpose_tape(lamps, count)
        alternative = transpose_tape(stale_lamps, count)
        parked = transpose_tape(parked_lamps, count)
        assert len(set(alternative)) == len(set(parked)) == 1
        layers.append(
            dict(
                native=native,
                stale=alternative[0],
                parked=parked[0],
                boundary_swaps=transpose_tape(act, count),
            )
        )
    return layers


def direct_pulses(board, sites, routes):
    # Two parallel sharp calibration cohorts per address, source load zero or
    # two. Encoding the lamps as two bits lets the pair be compared directly.
    counts = {}
    for (receiver, source), path in routes.items():
        head, ref = sites[receiver][4:]
        load = [0] * board.n
        # The high-source cohort changes public rails only; totals match.
        a, b = sites[source][:2]
        load[a], load[b] = 2, 1
        before = load.copy()
        hops = list(zip(path, path[1:]))
        for a, b in hops[:-1]:
            load[a], load[b] = load[b], load[a]
        assert load[ref] == 0
        early = load[head]
        a, b = hops[-1]
        load[a], load[b] = load[b], load[a]
        on_time = load[head]
        for a, b in reversed(hops):
            load[a], load[b] = load[b], load[a]
        assert load == before and (early, on_time) == (0, 2)
        key = str(len(hops) + 2)
        counts[key] = counts.get(key, 0) + 1
    return dict(
        by_first_discriminating_slot=counts,
        early_lamps_low_high=[0, 0],
        on_time_lamps_low_high=[0, 1],
        trials_per_address_and_setting=1,
    )


def acquire(level):
    board = federation(level)
    layout = mount(board)
    direct = paths(board, layout)
    clean = compile_layer(layout, direct)
    stale = compile_layer(layout, paths(board, layout, stale=True))
    symbols = stores_at_boundary(board, layout)
    native_tape = read_paths(clean, symbols)
    assert native_tape == [("published", e["source"]) for e in clean]
    stale_tape = read_paths(stale, symbols)
    assert stale_tape == [("cached", e["source"]) for e in stale]
    parked_tape = read_paths(clean, symbols, parked=True)
    prior, words = settings()
    native = [[decode(clean, native_tape, w, prior) for w in pair] for pair in words]
    modes = {
        "stale": decode(stale, stale_tape, 0, prior),
        "parked": decode(clean, parked_tape, 0, prior),
    }
    raw = sharp_replay(board, layout, direct, clean, stale, prior, words)
    for layer in range(2):
        np.testing.assert_array_equal(raw[layer]["native"], [w[layer] for w in native])
        assert {k: raw[layer][k] for k in modes} == modes
    for i in range(SIZE):
        changed = int(native[0][1], 16) ^ int(native[SIZE + i + 1][1], 16)
        expected = sum(1 << j for j, e in enumerate(clean) if e["source"] == i)
        assert changed == expected and native[0][0] == native[SIZE + i + 1][0]
    clock = clock_trial(direct)
    clock["direct_sharp_probes"] = direct_pulses(board, layout, direct)
    clock["any_placement"] = obstruction(clock["ordered_pair_slots"])
    # The local publication knob also belongs to each observer's transcript.
    # Use a remote intervention at fixed prior history and fixed own knob.
    remote_collisions = 0
    for r in range(SIZE):
        for e in clean:
            if e["receiver"] != r or e["source"] == r:
                continue
            k = SIZE + 1 + e["source"]
            own_controls = [
                (int(raw[1]["boundary_swaps"][run], 16) >> r) & 1 for run in [0, k]
            ]
            assert own_controls[0] == own_controls[1]
            remote_collisions += 1
    return dict(
        level=level,
        owners=board.cells,
        layout=layout,
        prior=hex(prior),
        inputs=[[hex(w) for w in pair] for pair in words],
        preparation=[
            dict(
                owner=i,
                local_cache_copy_path=direct[i, i],
                pointer_ref=site[5],
                cache_flip=site[2:4],
                publication_flip=site[:2],
            )
            for i, site in enumerate(layout)
        ],
        clock=clock,
        owner_geometry=owner_obstruction(board),
        native_plan=clean,
        stale_plan=stale,
        layer_cost=appointments(clean),
        boundary_pairs=sum(
            1
            for e in clean
            if sum(
                (a - b) ** 2
                for a, b in zip(SITES[e["source"]], SITES[e["receiver"]], strict=True)
            )
            == 4
        ),
        native=native,
        control_tapes=modes,
        lost_remote_single_bit_interventions_per_control=remote_collisions,
        local_boundary_swaps=[s["boundary_swaps"] for s in raw],
    )


def report():
    physical_columns()
    return dict(
        seed=SEED,
        alphabet="One bit per complete projection, low bit is first read; other observer slots are fixed by the calendar.",
        mounts=[acquire(level) for level in (1, 2)],
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path)
    args = cli.parse_args()
    data = json.dumps(report(), separators=(",", ":")).encode()
    if args.output:
        args.output.write_bytes(gzip.compress(data, mtime=0))
    else:
        assert gzip.decompress((HERE / "records.json.gz").read_bytes()) == data
        print(
            "Verified both native mounts, all source interventions, receiver tapes and deadlines."
        )
