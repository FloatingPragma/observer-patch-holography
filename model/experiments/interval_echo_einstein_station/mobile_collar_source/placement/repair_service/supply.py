"""Exact replenishment ledger for an explicitly maintained native boundary."""

import argparse
from hashlib import sha256
import json
from math import comb

import numpy as np

from model.candidate_laws.interval_echo import observations as obs
from . import local
from .measure import HERE, hardware, compare


def cost():
    write = local.operators()[1]
    h0, h1 = local.law.local_term_parts()
    ev, vec = np.linalg.eigh(h0 + 8 * h1)
    full = (vec * np.exp(-1j * np.pi * ev)) @ vec.T
    for a in (0, 1):
        for b in (0, 1):
            i = local.law.SRC.index((a, b))
            j = local.law.SRC.index((b, a))
            np.testing.assert_allclose(abs(full[:, i]) ** 2, np.eye(99)[j], atol=1e-12)
            column = write[2 * a + b]
            flip = bool(a != b)
            assert column[local.LEFT != a].sum() == 32 * flip
            assert column @ local.FLAG[local.MEMORY] == 64 * flip
    boards = {}
    for level in range(3):
        graph = obs.federation(level)
        _, _, _, signed, _ = hardware(graph)
        cut = int(np.count_nonzero(signed))
        assert all(
            12 * i + obs.VACANCY in graph.neighbours[12 * i + obs.MATTER]
            for i in range(graph.cells)
        )
        boards[str(graph.cells)] = dict(
            cut=cut,
            count_numerators=[comb(cut, j) for j in range(cut + 1)],
            count_denominator=2**cut,
            mean_pairs=cut / 2,
            native_corrective_swaps_mean=cut,
            signed_head_on_cut=1,
            unit_births_per_cycle=cut,
            held_unit_heads=cut,
            unit_archive_exports_per_cycle_after_first=cut,
        )
    return dict(
        protocol_sha256=sha256((HERE / "SUPPLY.md").read_bytes()).hexdigest(),
        boards=boards,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = cost()
    path = HERE / "supply.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        compare(json.loads(path.read_text()), result)
    for k, b in result["boards"].items():
        print(k, b["mean_pairs"], b["unit_births_per_cycle"])
