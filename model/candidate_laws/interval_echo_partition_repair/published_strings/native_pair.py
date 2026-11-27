"""Sharp endpoint memories obtained from the original blank Interval Echo seam."""
from dataclasses import dataclass
from functools import lru_cache
import argparse
import json
import math
from pathlib import Path
import random

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate as echo

from .memories import Pair, THETA

LAW = echo.native_law
LOGICAL = ((0, 0), (0, 1), (1, 0), (1, 1))
TICKETS = (0, 0, 0, 0, 0, 1, 2, 3)
HERE = Path(__file__).resolve().parent


@lru_cache(maxsize=1)
def full_generator():
    written_exchange, written_contrast = LAW.generator_parts()
    U = LAW.echo()
    return U @ (written_exchange+echo.KAPPA*written_contrast) @ U.conj().T


@dataclass(frozen=True)
class NativePair:
    """Specify an adjacent load alphabet and the port with reversed logical labels."""
    low: int
    reversed_port: str

    def __post_init__(self):
        if self.low not in (0, 1) or self.reversed_port not in ("left", "right"):
            raise ValueError("choose load alphabet 0/1 or 1/2 and a reversed port")

    def physical(self, a, b):
        if self.reversed_port == "left":
            a = 1-a
        else:
            b = 1-b
        return self.low+a, self.low+b

    @property
    def columns(self):
        return [LAW.SRC.index(self.physical(a, b)) for a, b in LOGICAL]

    @property
    def generator(self):
        full = full_generator()
        columns = self.columns
        leakage = np.delete(full[:, columns], columns, axis=0)
        if np.max(abs(leakage)) > 4e-14:
            raise ArithmeticError("the sharp source leaves the four blank-memory states")
        return full[np.ix_(columns, columns)]

    def gate(self, theta):
        return expm(-2j*theta*self.generator)

    def state(self, theta):
        sharp = np.diag(np.bincount(TICKETS, minlength=4)/len(TICKETS))
        gate = self.gate(theta)
        return gate @ sharp @ gate.conj().T

    def undo(self, state, theta):
        # All allowed physical energies are the integers 0,8,9; wait forward.
        gate = self.gate((-theta) % math.pi)
        return gate @ state @ gate.conj().T


def local_effect(low, digit):
    """Write the port against a known low reference load, read |contrast|, unwrite."""
    if low not in (0, 1) or digit not in (0, 1):
        raise ValueError("expected a binary endpoint form")
    # Either physical orientation of a mounted reference gives the same absolute label.
    labels = [index for index, (_, d) in enumerate(LAW.LABELS) if abs(d) == digit]
    result = []
    for reverse in (False, True):
        columns = [LAW.SRC.index((low, x) if reverse else (x, low)) for x in (low, low+1)]
        written = LAW.echo()[:, columns]
        rows = [LAW.N_SRC*(1+label)+i for label in labels for i in range(LAW.N_SRC)]
        projected = written[rows, :]
        result.append(projected.conj().T @ projected)
    return result


def entropy(state):
    p = np.linalg.eigvalsh(state)
    return float(-p @ np.log(p))


def relative_entropy(state, vacuum):
    eigen, basis = np.linalg.eigh(vacuum)
    diag = np.diag(basis.conj().T @ state @ basis).real
    return float(-entropy(state)-diag @ np.log(eigen))


def endpoint(state, port):
    four = state.reshape((2, 2, 2, 2))
    return np.trace(four, axis1=1, axis2=3) if port == "left" else np.trace(four, axis1=0, axis2=2)


def sharp_counts(state, seed, trials):
    # Both controlled experiments have dyadic probabilities. Recover these
    # tickets from the propagated density before any random draws.
    probabilities = np.diag(state).real
    tickets = np.rint(32*probabilities).astype(int)
    if tickets.sum() != 32 or not np.allclose(tickets/32, probabilities, atol=5e-14, rtol=0):
        raise ValueError("terminal probabilities leave the enrolled dyadic word law")
    cutoffs = np.cumsum(tickets)*2**48
    rng = random.Random(seed)
    counts = [0]*len(cutoffs)
    for _ in range(trials):
        x = rng.getrandbits(53)
        bin_index = next((i for i, stop in enumerate(cutoffs[:-1]) if x < stop), len(counts)-1)
        counts[bin_index] += 1
    return counts


def deposit():
    records = []
    for low in (0, 1):
        for side in ("left", "right"):
            code = NativePair(low, side)
            vacuum = code.state(THETA)
            restored = code.undo(vacuum, THETA)
            marked = code.undo(np.diag(np.diag(vacuum)), THETA)
            index = 2*low + int(side == "right")
            n = 16384
            rec = dict(low_load=low, reversed_port=side, physical=[list(code.physical(*s)) for s in LOGICAL],
                       trials_per_independent_bank=n,
                       restored=sharp_counts(restored, 192843+index*2, n),
                       marked=sharp_counts(marked, 192844+index*2, n))
            rows = []
            for delta in (-1/8, 0., 1/8):
                state = code.state(THETA+delta)
                cut = endpoint(state, "left")
                rows.append(dict(increment=delta, endpoint_one=float(cut[1, 1].real),
                                 whole_KL=relative_entropy(state, vacuum),
                                 cut_KL=relative_entropy(cut, endpoint(vacuum, "left")),
                                 cut_entropy_change=entropy(cut)-entropy(endpoint(vacuum, "left")),
                                 seam_energy=float(np.trace(code.generator @ state).real)))
            rec["exposures"] = rows
            records.append(rec)
    return dict(native_dimension=LAW.DIM, kappa=echo.KAPPA, initial_tickets=list(TICKETS),
                source_basis=[list(s) for s in LOGICAL], records=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = deposit()
    if args.write:
        (HERE/"native_pair.json").write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")
        return
    old = json.loads((HERE/"native_pair.json").read_text())
    for key in ("native_dimension", "kappa", "initial_tickets", "source_basis"):
        if data[key] != old[key]:
            raise ValueError(f"native record differs at {key}")
    for r, v in zip(data["records"], old["records"], strict=True):
        for k in ("low_load", "reversed_port", "physical", "trials_per_independent_bank", "restored", "marked"):
            if r[k] != v[k]:
                raise ValueError(f"native observer counts differ at {k}")
        for i, expected in zip(r["exposures"], v["exposures"], strict=True):
            for k in expected:
                if not math.isclose(i[k], expected[k], abs_tol=5e-13):
                    raise ValueError(f"native matrix replay differs at {k}")
                if k == "cut_entropy_change" and abs(i[k]-Pair(i["increment"]).cut_entropy_change) > 5e-13:
                    raise ValueError("logical entropy disagrees with the free seam")
    print("Both load alphabets and both port orientations reproduce the native pair controls.")


if __name__ == "__main__":
    main()
