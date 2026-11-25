"""Compare growing native record-cost spectra and blank doubled-load flight spectra."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from math import exp
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from model.candidate_laws.interval_echo.clock import Registers
from model.candidate_laws.interval_echo.projection import Branch
from model.executable_model import federation
from ..flight_confluence import KAPPA, doubled_hamiltonian, native_double


def words(depth):
    return ["1"*k+"0" for k in range(depth)] + ["1"*depth]


def draw(registers, depth):
    word = ""
    while len(word) < depth:
        bit = registers.fair(0, {0})
        word += str(bit)
        if bit == 0:
            break
    return word


def native_counts(depth, copies):
    support = words(depth)
    registers = Registers(federation(0), seed=0)
    assert max(abs(0.5-rate) for rate in registers.rates.values()) < 1e-12
    counts = []
    for starting_rail in (2, 3):
        complete = Counter()
        for branch in product((0, 1), repeat=depth*copies):
            registers.loads[[2, 3]] = (2, 0) if starting_rail == 2 else (0, 2)
            registers.rng = Branch(branch)
            histories = [draw(registers, depth) for _ in range(copies)]
            complete[histories[0] if all(h == histories[0] for h in histories) else "timeout"] += 1
        assert sum(complete.values()) == 2**(depth*copies)
        for h in support:
            assert complete[h] == 2**(copies*(depth-len(h)))
        counts.append(complete)
    assert counts[0] == counts[1]
    return dict(depth=depth, copies=copies, complete_tapes_per_start=2**(depth*copies),
                ordered_records=support, outcomes_per_start=dict(sorted(counts[0].items())))


def flight_panel(level):
    board = federation(level)
    native = doubled_hamiltonian(board, KAPPA)
    laplace = board.laplacian()
    # The stored two-load convention subtracts 12*kappa times the identity.
    assert not (2*native-24*KAPPA*sparse.eye(board.n)-laplace).count_nonzero()
    assert connected_components(laplace, directed=False, return_labels=False) == 1
    assert int(board.degree.max()) == 6 and int(board.degree.min()) == 5
    assert 4*int(board.degree.sum()) == 21*board.n
    panel = dict(cells=board.cells, ports=board.n, degree_sum=int(board.degree.sum()),
                 excited_mean_energy=str(Fraction(21*board.n, 8*(board.n-1))))
    if level <= 1:
        energy = np.linalg.eigvalsh(laplace.toarray()/2)
        assert abs(energy[0]) < 1e-12 and energy[1] > 0 and energy[-1] <= 6+1e-12
        sample = []
        for beta in (1, 4, 16):
            weights = np.exp(-beta*energy)
            prob = weights/weights.sum()
            upper = 1/(1+(board.n-1)*exp(-beta*21*board.n/(8*(board.n-1))))
            assert prob[0] <= upper+1e-13
            sample.append(dict(beta=beta, ground_probability=round(float(prob[0]), 12),
                               entropy_nats=round(float(-np.dot(prob, np.log(prob))), 12),
                               Jensen_ground_upper=round(upper, 12)))
        panel["spectrum_evaluation"] = sample
    return panel


def run():
    native_double()
    return dict(record_panels=[native_counts(d, b) for d, b in ((2, 1), (4, 1), (4, 2), (3, 3))],
                doubled_flight=[flight_panel(level) for level in range(3)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args=parser.parse_args()
    report=json.dumps(run(), indent=2, allow_nan=False)+"\n"
    if args.write:
        Path(__file__).with_name("results.json").write_text(report)
    else:
        print(report, end="")
