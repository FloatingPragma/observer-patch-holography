"""Simultaneous physical energies of the returned clock and a separate native four-port chain."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations
import json

import numpy as np
from scipy.sparse.csgraph import connected_components

from . import measure

HERE = measure.HERE
LOCK = "5925df1b2f9bb1f95af6e073633687c2de4211328ae8c0e2b353be5dd7a29da0"


def direct(board, edges, h0, h8):
    seams = set(board.seams)
    path = next(
        p
        for p in permutations(range(12), 4)
        if all((p[a], p[b]) in seams for a, b in edges)
        and sum(set(s) <= set(p) for s in board.seams) == 3
    )
    cuts = [e for e in board.seams if set(e) & set(path) and not set(e) <= set(path)]
    ix = [
        9 + 9 * measure.LAW.LABELS.index((c, 0)) + i for c in (0, 1) for i in range(9)
    ]
    assert not h0[:, ix].any() and not h8[:, ix].any()
    return dict(
        ports=path,
        physical_seams=[[path[a], path[b]] for a, b in edges],
        touching_dark_zero_contacts=len(cuts),
        dark_private_references=sum(board.degree[v] == 5 for v in path),
    )


def ordinary(board, owners, home):
    # Material can occupy four separate original clock rails. All other ports
    # are ordinary zero and every spectator seam is blank for this comparison.
    words = np.arange(81)
    flat = np.zeros((81, board.n), dtype=int)
    places = 12 * np.array(owners) + home
    flat[:, places] = (words[:, None] // 3 ** np.arange(4)) % 3
    i, j = np.array(board.seams).T
    d = flat[:, i] - flat[:, j]
    # Integer 4H on sharp blank inputs, including the native degree references.
    four = (
        2 * np.count_nonzero(d, axis=1)
        + 32 * (d * d).sum(axis=1)
        + 32 * ((flat**2) @ (6 - board.degree))
    )
    return flat, four


def priced(idxs, codes, baseline, table):
    n = len(idxs)
    live = cross = 0
    for c in codes[:, idxs]:
        counts = np.bincount(c, minlength=99).reshape(11, 9)
        blank = table[0] @ counts.sum(axis=0)
        live += int((table * counts).sum() - blank)
        cross += int(counts.sum(1) @ table @ counts.sum(0) - blank * n)
    bare = F(int(baseline[idxs % 81].sum()), 4 * n)
    diff = F(live, 4 * n)
    swap = F(cross, 4 * n * n) - diff
    return dict(
        full_board=bare + diff,
        matched_blank=bare,
        actual_minus_matched_blank=diff,
        native_swap_work_per_board=swap,
    )


def run():
    assert sha256((HERE / "PARKED.md").read_bytes()).hexdigest() == LOCK
    W, V, h8, h0 = measure.kernels()
    owners, seams, groups, edges, powers, ix, words, codes = measure.carrier()
    board = measure.observations.federation(0)
    home = measure.observations.MATTER
    d = direct(board, edges, h0, h8)
    flat, base = ordinary(board, owners, home)
    livecodes = np.array(
        [
            9 * words[:, 4 + i] + 3 * flat[ix % 81, a] + flat[ix % 81, b]
            for i, (a, b) in enumerate(seams)
        ]
    )
    GW = measure.graph(W, edges, powers, ix, codes)
    GV = measure.graph(V, edges, powers, ix, codes)
    records = {}
    for service, graph in [
        ("exchange_only", GV),
        ("echo", GW),
        ("echo_plus_exchange", GW + GV),
    ]:
        classes = connected_components(graph, directed=False)[1]
        for cl in np.unique(classes[:81]):
            ids = np.flatnonzero(classes == cl)
            summary = priced(ids, livecodes, base, np.diag(h8).reshape(11, 9))
            m2 = (measure.MAG[words[ids, 4:]] == 2).sum(1)
            for s in np.flatnonzero(classes[:81] == cl):
                label = "".join(map(str, words[s, :4]))
                rec = records.setdefault(
                    label, {"original_all_blank": F(int(base[s]), 4)}
                )
                rec[service] = dict(
                    **summary,
                    work_from_original=summary["full_board"]
                    - rec["original_all_blank"],
                    relative_to_same_fuel_exchange_only=summary["full_board"]
                    - (
                        summary["full_board"]
                        if service == "exchange_only"
                        else rec["exchange_only"]["full_board"]
                    ),
                    discarded_fuel_per_matched_blank=F(int(m2.sum()), len(ids)),
                )
    return dict(
        protocol_sha256=LOCK,
        direct_service_free_path=d,
        original_parked_home_ports=(12 * np.array(owners) + home).tolist(),
        original_external_seams=seams,
        original_full_carrier_ports=board.n,
        original_full_carrier_seams=len(board.seams),
        all81_parked_inputs=records,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    w = p.parse_args().write

    def encode(x):
        return str(x) if isinstance(x, F) else x.item()

    out = json.dumps(run(), sort_keys=True, indent=2, default=encode) + "\n"
    path = HERE / "parked.json"
    if w:
        path.write_text(out)
    elif path.read_text() != out:
        raise AssertionError("Parked whole-board replay mismatch")
    print("Verified simultaneous original board and separate direct-path energies.")
