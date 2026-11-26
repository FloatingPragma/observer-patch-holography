"""Interrogate a material ray twice on the same two-head 80-cell reference."""
import argparse
import json
from math import log, sqrt
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from ....measure import pairing, ranging

HERE = Path(__file__).resolve().parent
N = 2**26
SOURCE_ERROR = 1.81e-6
OVERLAPS = (1., 0., .5)
SETTING = ("blank", "full", "half")
N_BINS = 13


def probabilities():
    archive = HERE.parents[3] / "full_algebra" / "results.json"
    ref = json.loads(archive.read_text())
    stock = np.array(ref["x_real"])+1j*np.array(ref["x_imag"])
    board = pairing.formation.federation(1)
    heads = pairing.formation.region(board)
    assert board.cells == 80 and heads[0] == 0
    tables = np.zeros((len(OVERLAPS), 2, N_BINS))

    for odd, block in enumerate(stock):
        for charge in (-1, 1):
            H = ranging.generator(board, heads, charge, 8, bool(odd))
            start = np.zeros((H.shape[0], 4), complex)
            start[:4] = np.eye(4)
            U = -.75j*H
            early = expm_multiply(U, start, traceA=0)
            # Keep coherences within and outside the marked four-state port.
            # A and B are separately propagated amplitudes, not distributions.
            split = np.zeros_like(early)
            split[4:8] = early[4:8]
            pair = expm_multiply(U, np.hstack((split, early-split)), traceA=0)
            a, b = np.split(pair.reshape(board.n, 4, 8), 2, axis=2)
            for j, overlap in enumerate(OVERLAPS):
                arms = (b+overlap*a, sqrt(1-overlap**2)*a)
                for k, arm in enumerate(arms):
                    p = np.einsum("xai,ij,xaj->x", arm, block, arm.conj()).real/2
                    tables[j, k, :12] += p[:12]
                    tables[j, k, -1] += sum(p[12:])
    np.testing.assert_allclose(tables.sum(axis=(1, 2)), 1, atol=2e-13, rtol=0)
    assert tables.min() >= 0
    return tables


def intervals(counts, prior, Nprior):
    counts = np.asarray(counts, dtype=np.int64)
    assert counts.shape == (3, 2, N_BINS)
    assert np.all(counts >= 0) and np.all(counts.sum(axis=(1, 2)) == N)
    assert len(prior) == N_BINS and min(prior) >= 0 and sum(prior) == Nprior
    late = counts.sum(axis=1)/N
    # One simultaneous event for 13 late outcomes on each of four banks.
    def radius(n):
        return sqrt(log(2*4*N_BINS/.01)/(2*n)) + SOURCE_ERROR
    e, f = radius(N), radius(Nprior)
    change = late[0]-late[1]
    index = int(np.argmax(abs(change)))
    return dict(terminal_index=index,
                blank_minus_full_at_index=float(change[index]),
                simultaneous_interval=[change[index]-2*e, change[index]+2*e],
                resolved_indices=np.flatnonzero(abs(change) > 2*e).tolist(),
                max_half_minus_midpoint=float(max(abs(late[2]-(late[0]+late[1])/2))),
                half_allowance=2*e,
                max_blank_minus_prior=float(max(abs(late[0]-np.array(prior)/Nprior))),
                prior_allowance=e+f, late_frequency_radius=e)


def run():
    tables = probabilities()
    counts = [np.random.default_rng(211746+i).multinomial(N, table.ravel()).reshape(2, N_BINS)
              for i, table in enumerate(tables)]
    prior = json.loads((HERE.parent.parent/"port_apertures.json").read_text())
    report = intervals(counts, prior["counts_including_outside"], prior["starts"])
    return dict(settings=SETTING, starts_per_setting=N,
                joint_counts=np.array(counts).tolist(), summary=report)


def verify(archive):
    tables = probabilities()
    raw = np.array(archive["joint_counts"])
    prior = json.loads((HERE.parent.parent/"port_apertures.json").read_text())
    s = intervals(raw, prior["counts_including_outside"], prior["starts"])
    for k, v in s.items():
        np.testing.assert_allclose(archive["summary"][k], v, rtol=2e-13, atol=2e-13)
    assert s["resolved_indices"] and s["max_half_minus_midpoint"] < s["half_allowance"]
    assert s["max_blank_minus_prior"] < s["prior_allowance"]
    # Audit the frozen raw joint counts against the native law without
    # requiring two NumPy versions to simulate identical multinomial draws.
    marginal_error = sqrt(log(2*tables.size/.01)/(2*N))
    assert np.max(abs(raw/N-tables)) < marginal_error
    assert not np.any(raw[tables == 0])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    if parser.parse_args().write:
        (HERE/"results.json").write_text(json.dumps(run(), indent=2, allow_nan=False)+"\n")
    else:
        verify(json.loads((HERE/"results.json").read_text()))
        print("Verified all 78 pointer/terminal bins and the independently held exposure.")
