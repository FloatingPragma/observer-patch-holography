"""Read the comparison signs with one held load before predicting a free pair."""
from collections import defaultdict

import numpy as np
from scipy import sparse

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.experiments.interval_echo_consensus.record_action import changes, signed_table
from model.experiments.interval_echo_consensus.sign_graph import questions


def read_ports(board):
    """One fresh sharp-zero bank per port. A pi pulse interchanges the two sharp sign states."""
    offers, _ = questions(board, "primitive")
    by_contact = defaultdict(list)
    rows, cols = [], []
    word = signed_table()
    for j, (a, b, address) in enumerate(offers):
        if b is None:
            # Unary reference compares this port with an empty reference port.
            rows.append(a)
            cols.append(j)
        else:
            by_contact[a, b].append((j, board.seams[address]))
    for (a, b), slots in by_contact.items():
        readers = [edge for _, edge in slots]
        loads = {v: 0 for edge in readers + [(a, b)] for v in edge}
        for x, y in ((0, 0), (1, 0), (0, 1), (1, 1)):
            loads[a], loads[b] = x, y
            delta = changes(readers, (a, b), loads, word)
            # The held two-load check tests the predicate at the remaining local input.
            if x == y:
                if any(delta):
                    raise ArithmeticError("An equal-load write drives a comparison sign")
                continue
            for (j, _), d in zip(slots, delta, strict=True):
                if d:
                    rows.append(a if x else b)
                    cols.append(j)
    return sparse.csr_matrix((np.ones(len(cols), dtype=np.int64), (rows, cols)),
                             shape=(board.n, len(offers))), by_contact


def geometry(board):
    """The off-diagonal Gram entries count shared flipped signs, including every nonedge."""
    records, contacts = read_ports(board)
    gram = (records @ records.T).tocsr()
    single = gram.diagonal()
    gram.setdiag(0)
    gram.eliminate_zeros()
    a, b = np.array(board.seams).T
    if (gram.nnz != 2*len(board.seams) or np.min(gram[a, b]) <= 0
            or np.any(single != single[0])):
        raise ArithmeticError("Single-load responses do not cover the board with a flat rest reference")
    return gram, contacts, int(single[0]), records.shape[1]


def pair_source(board, launch):
    gram, contacts, rest, bits = geometry(board)
    pair = PairSector(board)
    common = np.asarray(gram[pair.aa, pair.bb]).ravel()
    tag = contacts[launch]
    width = len(tag)
    if width == 0 or gram[launch] != width:
        raise ArithmeticError("Missing a tagged comparison in the two launch-port readings")
    # Every state uses the same two-single-load rest subtraction, even when
    # the launch bank is sharp. A sector-dependent subtraction changes the sign readout.
    exchange = pair.exchange()
    odd = exchange - sparse.diags(2*common, dtype=float, format="csr")
    either_end = np.isin(pair.aa, launch) ^ np.isin(pair.bb, launch)
    return pair, odd, common, either_end.astype(int), width, rest, bits
