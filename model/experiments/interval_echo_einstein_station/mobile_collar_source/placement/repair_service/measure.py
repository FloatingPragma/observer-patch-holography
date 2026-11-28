"""Stored binary contrast, failed and renewed writes, and matter transport on native matchings."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from model.candidate_laws.interval_echo import observations as obs
from scipy import sparse

from . import compact, local

HERE = Path(__file__).resolve().parent
APPOINTMENTS = (0, 1, 4, 16, 64, 256)
SEED = 117281
COPIES = 512


def hardware(board):
    layers = obs.schedules(board)[0]
    obs.gateway_words(board)
    pairs = np.asarray(board.external) // 12
    lookup = {tuple(e): i for i, e in enumerate(pairs)}
    layers = [(x, y, np.array([lookup[p] for p in zip(x, y)])) for x, y, _, _ in layers]
    a, b = pairs.T
    sizes = board.cells
    graph = sparse.csr_matrix((np.ones(len(a)), (a, b)), shape=(sizes, sizes))
    order, parent = sparse.csgraph.breadth_first_order(graph, 0, directed=False)
    distance = np.zeros(sizes, dtype=int)
    for i in order[1:]:
        distance[i] = distance[parent[i]] + 1
    near = np.argsort(distance, kind="stable")[: sizes // 2]
    patch = np.isin(np.arange(sizes), near).astype(np.int16)
    signed = patch[a] - patch[b]
    incidence = sparse.csr_matrix(
        (
            np.r_[np.ones(len(a), dtype=np.int16), -np.ones(len(a), dtype=np.int16)],
            (np.r_[a, b], np.tile(np.arange(len(a)), 2)),
        ),
        shape=(sizes, len(a)),
    )
    return layers, pairs, patch, signed, incidence


def summarize(x):
    x = np.asarray(x, dtype=float)
    return [float(x.mean()), float(x.std(ddof=1) / np.sqrt(len(x)))]


def fresh_prediction(layers, patch, signed):
    first = patch.astype(float)
    means = np.zeros(len(signed))
    result = {}
    for t in range(257):
        if t in APPOINTMENTS:
            result[str(t)] = dict(
                load_contrast=float(
                    first[patch == 1].mean() - first[patch == 0].mean()
                ),
                directed_unit_head=float(means @ signed / sum(signed != 0)),
            )
        if t < 256:
            for a, b, e in layers:
                means[e] = first[a] - first[b]
                first[a] = first[b] = (first[a] + first[b]) / 2
    return result


def trajectories(board, setting):
    layers, pairs, patch, signed, incidence = hardware(board)
    n, m = board.cells, len(pairs)
    cut = signed != 0
    load = np.broadcast_to(patch, (COPIES, n)).astype(np.uint8).copy()
    mem = np.zeros((COPIES, m), np.uint8)
    births, losses, exports = (np.zeros(COPIES, np.int32) for _ in range(3))
    moved = np.zeros(COPIES, dtype=np.int32)
    last = np.zeros((4, COPIES), dtype=np.int32)
    lookup = local.samplers()[setting]
    rng = np.random.default_rng([SEED, board.cells, local.SETTINGS.index(setting)])
    rows = {}
    for t in range(257):
        if t in APPOINTMENTS:
            f = local.CURRENT[mem]
            active = local.FLAG[mem]
            charge = load + (incidence @ f.T).T
            loss = patch.sum() - load[:, patch == 1].sum(axis=1, dtype=np.int32)
            stored = f @ signed
            count = active.sum(axis=1)
            assert np.array_equal(count, births - losses - exports)
            if setting == "retained":
                np.testing.assert_array_equal(
                    charge, np.broadcast_to(patch, charge.shape)
                )
            np.testing.assert_array_equal(load.sum(axis=1), patch.sum())
            observed = dict(
                load_contrast=load[:, patch == 1].mean(axis=1)
                - load[:, patch == 0].mean(axis=1),
                signed_unit_head=(local.SIGNED[mem] @ signed) / sum(cut),
                cut_current=stored / sum(cut),
                cut_unit_excess=active[:, cut].mean(axis=1)
                - active[:, ~cut].mean(axis=1),
                unit_fraction=count / m,
                stored_unit=count,
                charge_residual_abs=np.abs(charge - patch).sum(axis=1) / n,
                original_patch_loss=loss,
                cumulative_births=births,
                cumulative_returns=losses,
                cumulative_export=exports,
                cumulative_load_moves=moved,
            )
            for label, vec in zip(
                ("birth", "return", "export", "move"), last, strict=True
            ):
                observed[label + "_per_edge_last_cycle"] = vec / m
            rows[str(t)] = {k: summarize(v) for k, v in observed.items()}
            replica = local.SIGNED[mem[::2]] * local.SIGNED[mem[1::2]]
            for label, region in (("cut", cut), ("bulk", ~cut)):
                rows[str(t)]["single_unit_" + label] = summarize(
                    active[:, region].mean(axis=1)
                )
                rows[str(t)]["two_copy_signed_product_" + label] = summarize(
                    replica[:, region].mean(axis=1)
                )
        if t == 256:
            break
        before = np.stack((births, losses, exports, moved))
        for a, b, e in layers:
            old = mem[:, e]
            left = load[:, a]
            right = load[:, b]
            before_flag = local.FLAG[old]
            code = 4 * old + 2 * left + right
            draw = rng.integers(lookup.shape[1], size=code.shape, dtype=np.uint16)
            target = lookup[code, draw]
            new = target // 4
            new_flag = local.FLAG[new]
            if setting == "fresh":
                exports += before_flag.sum(axis=1, dtype=np.int32)
                births += new_flag.sum(axis=1, dtype=np.int32)
            else:
                births += (new_flag > before_flag).sum(axis=1, dtype=np.int32)
                losses += (new_flag < before_flag).sum(axis=1, dtype=np.int32)
            after = (target % 4) // 2
            moved += (left != after).sum(axis=1, dtype=np.int32)
            load[:, a], load[:, b], mem[:, e] = after, target % 2, new
        last = np.stack((births, losses, exports, moved)) - before
    return rows


def control_checks():
    gates = local.operators()[0]
    for name, gate in gates.items():
        for word in (0, 3):
            allowed = [
                word + 4 * j
                for j, mem in enumerate(local.MEM)
                if mem is None or mem[1] == 0
            ]
            assert np.all(gate[allowed][:, np.setdiff1d(np.arange(28), allowed)] == 0)
    # Pure half contacts preserve the blank after any load placement.
    half = local.operators()[2]
    np.testing.assert_array_equal(half[:4, 4:], 0)
    return dict(
        homogeneous_zero_births=0,
        homogeneous_one_births=0,
        half_contact_without_writes_creates_heads=0,
    )


def acquire():
    board_panels = {}
    for level in range(3):
        b = obs.federation(level)
        layers, pairs, patch, signed, _ = hardware(b)
        table = {s: trajectories(b, s) for s in local.SETTINGS}
        board_panels[str(b.cells)] = dict(
            owners=np.flatnonzero(patch).tolist(),
            boundary=int(sum(signed != 0)),
            echo_only_load_contrast_lower=1 - 4 * sum(signed != 0).item() / b.cells,
            edges=len(pairs),
            layer_sizes=[len(a) for a, _, _ in layers],
            fresh_exact=fresh_prediction(layers, patch, signed),
            fresh_fixed_inventory_unit_fraction=b.cells / (2 * (b.cells - 1)),
            trials=table,
        )
        print("complete", b.cells, flush=True)
    return dict(
        protocol_sha256=sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        independent_copies=COPIES,
        seed=SEED,
        appointments=APPOINTMENTS,
        replica_protocol_sha256=sha256((HERE / "REPLICAS.md").read_bytes()).hexdigest(),
        paired_trajectories=COPIES // 2,
        paired_99_radius_per_unit_range=float(
            np.sqrt(np.log(2 * 2048 / 0.01) / (2 * (COPIES // 2)))
        ),
        controls=control_checks(),
        simultaneous_99_radius_per_unit_range=float(
            np.sqrt(np.log(2 * 2048 / 0.01) / (2 * COPIES))
        ),
        tables_quanta=4096,
        local_tables={k: v.tolist() for k, v in local.operators()[0].items()},
        compact={s: compact.solve(obs.federation(0), s) for s in local.SETTINGS},
        boards=board_panels,
    )


def compare(saved, new):
    if isinstance(saved, dict):
        assert saved.keys() == new.keys()
        for k in saved:
            compare(saved[k], new[k])
    elif isinstance(saved, list):
        assert len(saved) == len(new)
        for x, y in zip(saved, new, strict=True):
            compare(x, y)
    elif isinstance(saved, float):
        np.testing.assert_allclose(saved, new, atol=1e-11, rtol=1e-12)
    else:
        assert saved == new


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = acquire()
    path = HERE / "results.json"
    if args.write:
        path.write_text(
            json.dumps(result, separators=(",", ":"), allow_nan=False) + "\n"
        )
    else:
        compare(json.loads(path.read_text()), result)
