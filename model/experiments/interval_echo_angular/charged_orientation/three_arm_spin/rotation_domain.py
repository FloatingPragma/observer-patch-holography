"""Which proper receiver rotations are reachable using mask-preserving flight pulses?"""

from collections import Counter
import json
from pathlib import Path

import numpy as np

from model.executable_model import cell
from model.executable_model.symmetry import rotations
from .apparatus import LAW, fixture


def sector_check():
    """Build 8h for each registered term exactly, including the blank/written cross block."""
    write = LAW._write_pattern()
    reverse, squares = LAW._label_ops()
    erase = 2 * np.eye(len(write), dtype=int) - write @ write.T
    labels = np.repeat([-1, *(abs(d) for _, d in LAW.LABELS)], len(LAW.SRC))
    mismatch = labels[:, None] != labels[None, :]
    terms = {}
    for name, generator, native in zip(
        ("exchange", "contrast"),
        (np.eye(len(reverse), dtype=int) - reverse, 2 * squares),
        LAW.local_term_parts(),
        strict=True,
    ):
        g2 = np.kron(generator, np.eye(len(LAW.SRC), dtype=int))
        # The two off-diagonal blocks of 8h are sqrt(2) times this integer
        # array and its transpose. The written block has denominator eight.
        cross = write.T @ g2 @ erase
        assert not cross.any()
        h8 = np.zeros((LAW.DIM, LAW.DIM), dtype=int)
        h8[:9, :9] = 2 * write.T @ g2 @ write
        h8[9:, 9:] = erase @ g2 @ erase
        np.testing.assert_allclose(h8, 8 * native, atol=1e-12)
        violations = int(np.count_nonzero(h8[mismatch]))
        assert violations == 0
        terms[name] = dict(
            native_nonzero_entries=int(np.count_nonzero(h8)),
            sector_changing_entries=violations,
        )
    # Counterexample outside the pulse class. An addressed echo writes a
    # blank unequal pair into magnitude one; its inverse restores the word.
    word = LAW.SRC.index((0, 1))
    target = write[:, word]
    assert target @ target == 2
    assert not target[labels[9:] != 1].any()
    column = np.eye(LAW.DIM)[:, word]
    np.testing.assert_allclose(
        LAW.echo() @ column, np.r_[np.zeros(9), target / np.sqrt(2)], atol=1e-13
    )
    np.testing.assert_allclose(LAW.echo() @ (LAW.echo() @ column), column, atol=1e-13)
    return dict(
        dimension=LAW.DIM,
        sector_ranks={str(k): int(np.sum(labels == k)) for k in (-1, 0, 1, 2)},
        terms=terms,
        blank_to_magnitude_one_echo_probability=1,
    )


def rotation_order(r):
    ident = tuple(range(len(r)))
    p = ident
    for n in range(1, 6):
        p = tuple(r[j] for j in p)
        if p == ident:
            return n
    raise ValueError("unexpected proper icosahedral rotation")


def masks():
    lab = fixture()
    owner = lab["tip"] // 12
    hub = lab["tip"] % 12
    ports = {v % 12 for v in lab["ports"] if v // 12 == owner}

    def local(edges):
        return {
            tuple(sorted((a % 12, b % 12)))
            for a, b in edges
            if a // 12 == b // 12 == owner
        }

    live = local(lab["marked"])
    active_blank = local(lab["blank"])
    assert len(live) == 3
    all_edges = set(cell().seams)
    # The compensating dimer is outside this cell. Any untouched seam starts blank.
    assert all(v // 12 != owner for v in lab["donor"])
    silent = {e for e in all_edges if ports.intersection(e)} - live - active_blank
    blank = all_edges - live - silent
    assert len(blank) + len(silent) + len(live) == 30

    def rotated(r, edges):
        return {tuple(sorted((r[a], r[b]))) for a, b in edges}

    rows = []
    for i, r in enumerate(rotations()):
        moved = rotated(r, live)
        rows.append(
            dict(
                rotation=i,
                permutation=list(r),
                order=rotation_order(r),
                hub_after=owner * 12 + r[hub],
                live_overlap=len(moved & live),
                same_live=moved == live,
                same_silent=rotated(r, silent) == silent,
                same_blank=rotated(r, blank) == blank,
                same_port=r[hub] == hub,
            )
        )

    def count(key):
        return sum(row[key] for row in rows)

    half = [row for row in rows if row["order"] == 2]
    assert len(half) == 15 and count("same_port") == 5
    assert count("same_live") == count("same_silent") == 1
    assert not any(row["same_port"] or row["same_live"] for row in half)
    return dict(
        owner=owner,
        global_port=lab["tip"],
        absolute_one_edges=[[12 * owner + a, 12 * owner + b] for a, b in sorted(live)],
        blank_edges=len(blank),
        silent_edges=len(silent),
        mask_stabilizer=[row["rotation"] for row in rows if row["same_live"]],
        port_stabilizer=[row["rotation"] for row in rows if row["same_port"]],
        half_turn_overlap_histogram={
            str(n): v
            for n, v in sorted(Counter(row["live_overlap"] for row in half).items())
        },
        proper_rotations=rows,
    )


def run():
    return dict(native=sector_check(), receiver=masks())


if __name__ == "__main__":
    path = Path(__file__).with_name("rotation_domain.json")
    path.write_text(json.dumps(run(), indent=2) + "\n")
