"""Resident populations, fuel and native pairing work in closed ternary repair classes."""

from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import argparse
import json

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from model.candidate_laws.interval_echo import observations
from model.experiments.interval_echo_einstein_station.mobile_collar_source.placement.repair_service import (
    compact,
)
from .. import source

HERE = Path(__file__).resolve().parent
PROTOCOL = "3343849de7eeee26cc0b546e06f6e9e910fbfe0e491f3b5bda0f06b844e8b813"
LAW = source.LAW
MAG = np.array([0, *[abs(d) for _, d in LAW.LABELS]])
LOAD = np.array(LAW.SRC)


def kernels():
    f = LAW._write_pattern()
    r_ = np.kron(LAW._label_ops()[0], np.eye(9, dtype=int))
    q = 2 * np.eye(90, dtype=int) - f @ f.T
    W = np.zeros((99, 99), dtype=int)
    W[9:, :9] = 32 * f**2
    W[:9, 9:] = 32 * f.T**2
    W[9:, 9:] = 16 * q**2
    r2 = np.zeros_like(W)
    r2[:9, :9] = 2 * LAW._load_ops()[0]
    r2[9:, 9:] = 2 * r_ + f @ f.T @ (np.eye(90, dtype=int) - r_)
    V = 32 * np.eye(99, dtype=int) + 8 * r2**2
    identity = 64 * np.eye(99, dtype=int)
    assert not (W @ (W - identity) @ (2 * W + identity)).any()
    assert not (V @ (V - identity) @ (2 * V - identity)).any()
    for mat in (W, V):
        assert np.array_equal(mat, mat.T) and np.all(mat.sum(axis=0) == 64)
    h0 = source.formation.emission.four_h(0)
    h8 = source.formation.emission.four_h(8)
    H = h8 / 4
    e, v = np.linalg.eigh(H)
    half = (v * np.exp(-0.5j * np.pi * e)) @ v.T
    np.testing.assert_allclose(V / 64, abs(half) ** 2, atol=1e-13)
    np.testing.assert_allclose(W / 64, abs(LAW.echo()) ** 2, atol=1e-14)
    mags = MAG[np.arange(99) // 9]
    fuel = np.column_stack(
        [np.sum(LOAD[np.arange(99) % 9] == d, axis=1) + (mags == 2) for d in [0, 2]]
    )
    for mat in (W, V, h0, h8):
        a, b = np.nonzero(mat)
        assert np.array_equal(fuel[a], fuel[b])
    for mat in (V, h0, h8):
        a, b = np.nonzero(mat)
        assert np.array_equal(mags[a], mags[b])
    return W, V, h8, h0


def carrier():
    owners, seams, groups = compact.path(observations.federation(0))
    edges = [(owners.index(a // 12), owners.index(b // 12)) for a, b in seams]
    assert sorted(i for group in groups for i in group) == list(range(3))
    powers = np.array([1, 3, 9, 27, 81, 891, 9801])
    ix = np.arange(81 * 11**3)
    words = (ix[:, None] // powers) % np.array([3, 3, 3, 3, 11, 11, 11])
    codes = np.array(
        [
            9 * words[:, 4 + i] + 3 * words[:, a] + words[:, b]
            for i, (a, b) in enumerate(edges)
        ]
    )
    return owners, seams, groups, edges, powers, ix, words, codes


def graph(mat, edges, powers, ix, codes):
    rows, cols = [], []
    for edge, (a, b) in enumerate(edges):
        code = codes[edge]
        for old, new in zip(*np.nonzero(mat), strict=True):
            if new == old:
                continue
            initial = ix[code == old]
            change = (
                (new // 9 - old // 9) * powers[4 + edge]
                + (new % 9 // 3 - old % 9 // 3) * powers[a]
                + (new % 3 - old % 3) * powers[b]
            )
            rows.append(initial)
            cols.append(initial + change)
    M = coo_matrix(
        (
            np.ones(sum(map(len, rows)), dtype=bool),
            (np.concatenate(rows), np.concatenate(cols)),
        ),
        shape=(len(ix), len(ix)),
    ).tocsr()
    assert (M != M.T).nnz == 0
    return M


def frac(n, d=1):
    return str(Fraction(int(n), int(d)))


def records(idxs, words, codes, h8, h0, W, do_write):
    n = len(idxs)
    r = words[idxs]
    flags = MAG[r[:, 4:]]
    n0 = (r[:, :4] == 0).sum(1)
    n2 = (r[:, :4] == 2).sum(1)
    big = (flags == 2).sum(axis=1)
    assert np.ptp(n0 + big) == 0 and np.ptp(n2 + big) == 0
    out = dict(
        size=n,
        fuel=[int((n0 + big)[0]), int((n2 + big)[0])],
        resident={},
        extrema={},
        repair_flows={},
        energies={},
    )
    for m in (1, 2):
        z = (flags == m).sum(axis=1)
        out["resident"][str(m)] = frac(z.sum(), n)
        out["extrema"][str(m)] = [int(z.min()), int(z.max())]
        active = MAG[np.arange(99) // 9] == m
        born = 0
        lost = 0
        if do_write:
            to = W @ active.astype(int)
            for edge in range(3):
                c = codes[edge, idxs]
                born += to[c[~active[c]]].sum()
                lost += (64 - to[c[active[c]]]).sum()
        assert born == lost
        out["repair_flows"][str(m)] = dict(
            birth_per_cycle=frac(born, 64 * n), return_per_cycle=frac(lost, 64 * n)
        )
    for label, H in [("native", h8), ("exchange_only", h0)]:
        energies = np.diag(H).reshape(11, 9)
        actual = blank = swapped = 0
        edge_rows = []
        for co in codes[:, idxs]:
            table = np.bincount(co, minlength=99).reshape(11, 9)
            bare = energies[0] @ table.sum(axis=0)
            live = (energies * table).sum()
            cross = table.sum(axis=1) @ energies @ table.sum(axis=0)
            actual += int(live)
            blank += int(bare)
            swapped += int(cross)
            edge_rows.append(
                dict(
                    resident=frac(live, 4 * n),
                    matched_blank=frac(bare, 4 * n),
                    independent_swap_work=frac(cross - live * n, 4 * n * n),
                )
            )
        out["energies"][label] = dict(
            resident=frac(actual, 4 * n),
            matched_blank=frac(blank, 4 * n),
            actual_minus_matched_blank=frac(actual - blank, 4 * n),
            swap_work_per_board=frac(swapped - actual * n, 4 * n * n),
            work_for_two_boards=frac(2 * (swapped - actual * n), 4 * n * n),
            edges=edge_rows,
        )
    return out


def run():
    assert sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest() == PROTOCOL
    W, V, h8, h0 = kernels()
    owners, seams, groups, edges, powers, ix, words, codes = carrier()
    GW = graph(W, edges, powers, ix, codes)
    GV = graph(V, edges, powers, ix, codes)
    # Test invariance under the full Hamiltonian using exact graph support,
    # including off-diagonal entries invisible at a single timed appointment.
    GH = graph(abs(h8) + abs(h0), edges, powers, ix, codes)
    hrow, hcol = GH.nonzero()
    starts = ix[np.all(words[:, 4:] == 0, axis=1)]
    panels = {}
    for service, g in [
        ("echo", GW),
        ("echo_plus_exchange", GW + GV),
        ("exchange_only", GV),
    ]:
        _, labels = connected_components(g, directed=False)
        leak = set(labels[hrow[labels[hrow] != labels[hcol]]])
        classes = []
        placement = {}
        for comp in np.unique(labels[starts]):
            idxs = np.flatnonzero(labels == comp)
            rec = records(idxs, words, codes, h8, h0, W, service != "exchange_only")
            rec["stationary_for_direct_path_free_hold"] = int(comp) not in leak
            setup = []
            for s in starts[labels[starts] == comp]:
                key = "".join(map(str, words[s, :4]))
                setup.append(key)
                placement[key] = len(classes)
            rec["blank_inputs"] = setup
            classes.append(rec)
        assert len(placement) == 81
        if service == "echo_plus_exchange":
            assert all(c["stationary_for_direct_path_free_hold"] for c in classes)
        panels[service] = dict(classes=classes, initial_to_class=placement)
    return dict(
        protocol_sha256=PROTOCOL,
        owners=owners,
        physical_seams=seams,
        groups=groups,
        all_states=len(ix),
        full_seam_states=99,
        panels=panels,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    write = parser.parse_args().write
    output = (
        json.dumps(run(), indent=2, sort_keys=True, default=lambda x: x.item()) + "\n"
    )
    target = HERE / "results.json"
    if write:
        target.write_text(output)
    elif target.read_text() != output:
        raise AssertionError("Stationary magnitude-stock replay mismatch")
    print(
        "Verified ternary closed classes, fuel conservation and native stock energies."
    )
