"""Integer block certificates for all-time adjacency on simple echosahedral federations."""
import json
from fractions import Fraction

import numpy as np
from scipy import sparse

from model.executable_model import cell, federation
from .sign_graph import common_width


def validate_structure(board):
    if (board.n < 12 or board.n % 12 or len(set(board.seams)) != len(board.seams)
            or any(not 0 <= a < b < board.n for a, b in board.seams)):
        raise ValueError("The certificate needs complete twelve-port cells with simple contacts")
    internal = [set() for _ in range(board.cells)]
    external, gateways = set(), set()
    for a, b in board.seams:
        x, y = a//12, b//12
        if x == y:
            internal[x].add((a%12, b%12))
        else:
            if (x, y) in external or a in gateways or b in gateways:
                raise ValueError("The block estimate requires a simple intercell graph and an external matching")
            external.add((x, y))
            gateways.update((a, b))
    template = set(cell().seams)
    if any(s != template for s in internal):
        raise ValueError("A primitive cell has a different internal graph")
    if board.cells > 1 and not np.all(np.bincount([v//12 for v in gateways], minlength=board.cells) == 3):
        raise ValueError("Each federated cell needs exactly three gateways")


def certificate(board):
    validate_structure(board)
    adjacent = set(board.seams)
    rows, cols, separated = [], [], {}
    slot = {e: j for j, e in enumerate(board.seams)}
    adjacent_rows, adjacent_cols = [], []
    for j, (a, b) in enumerate(board.seams):
        for moving, fixed in ((a, b), (b, a)):
            for new in board.neighbours[moving]:
                if new == fixed:
                    continue
                target = tuple(sorted((new, fixed)))
                if target in adjacent:
                    adjacent_rows.append(slot[target])
                    adjacent_cols.append(j)
                else:
                    index = separated.setdefault(target, len(separated))
                    rows.append(index)
                    cols.append(j)
    two_b = sparse.csr_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)),
                              shape=(len(separated), len(adjacent)))
    if two_b.data.max() != 1:
        raise ArithmeticError("A pair hop was counted more than once")
    gram = two_b.T @ two_b
    max_row = int(np.asarray(gram.sum(axis=1)).max())
    if max_row > 20:
        raise ArithmeticError("The integer Gram bound exceeds the combinatorial estimate")
    # The triangular-face incidence accounts for every off-diagonal within P.
    faces = board.faces()
    face_columns = [[slot[a, b], slot[a, c], slot[b, c]] for a, b, c in faces]
    F = sparse.csr_matrix((np.ones(3*len(faces), dtype=np.int64),
                          (np.repeat(np.arange(len(faces)), 3), np.ravel(face_columns))),
                         shape=(len(faces), len(slot)))
    diag2 = np.array([len(board.neighbours[a])+len(board.neighbours[b])-2 for a, b in board.seams])
    face_gram = F.T @ F
    face_counts = face_gram.diagonal()
    actual_moves = sparse.csr_matrix((np.ones(len(adjacent_rows), dtype=np.int64),
                                      (adjacent_rows, adjacent_cols)), shape=face_gram.shape)
    if (actual_moves - face_gram + sparse.diags(face_counts, dtype=np.int64)).nnz:
        raise ArithmeticError("A material hop in the adjacent subspace has no unique triangular face")
    limit = int((diag2 + face_counts).max())
    if limit > 12:
        raise ArithmeticError("Adjacent exchange exceeds the spectral estimate")
    width = common_width(board)
    gap = Fraction(2*width)-Fraction(limit,2)
    if gap <= 0:
        raise ValueError("No separated spectral blocks")
    bound = gap**2/(gap**2 + max_row)
    return dict(cells=board.cells, seams=len(adjacent), separated_neighbors=len(separated),
                twice_adjacent_upper=limit, four_cross_norm_squared_upper=max_row,
                minimum_strength=width, all_time_adjacent_lower=str(bound))


if __name__ == "__main__":
    print(json.dumps([certificate(b) for b in (cell(), *(federation(i) for i in range(3)))], indent=2))
