"""Acquire port-resolved and coarse passive echoes on the registered one-load boards."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from .. import channel, headline


SHOTS = 2**22
HERE = Path(__file__).resolve().parent


def draw(rng, values):
    rates = np.array(values, dtype=float)
    if not np.isfinite(rates).all() or rates.min() < -2e-13 or abs(rates.sum() - 1) > 2e-12:
        raise ValueError("invalid launch distribution")
    rates = np.maximum(rates, 0)
    return rng.multinomial(SHOTS, rates / rates.sum()).tolist()


def acquire():
    rng = np.random.default_rng(692113)
    panels = []
    for station in channel.SPEC["stations"]:
        board = federation(station["level"])
        h = candidate.generator(board) - 6 * candidate.KAPPA * eye(board.n, format="csr")
        initial = channel.encode(board, station)
        readers = np.arange(12)
        basis = np.zeros((board.n, len(readers)))
        basis[readers, readers] = 1
        for gain in (1, 0.5):
            for t in channel.SPEC["exposures"]:
                step = -1j * gain * t * h
                sent = expm_multiply(step, initial)
                from_readers = expm_multiply(step, basis)
                admitted = np.zeros_like(sent)
                admitted[readers] = sent[readers]
                from_cell = expm_multiply(step, admitted)
                # The fresh one-way challenge uses the raw board incidence and a
                # different polynomial evolution implementation.
                one_way = headline.rebuild(station, gain * t)
                for bit, source in enumerate(station["ports"]):
                    arrived = abs(sent[readers, bit])**2
                    echo = arrived * abs(from_readers[source])**2
                    coarse = abs(from_cell[source, bit])**2
                    quiet = 1 - arrived.sum()
                    resolved = draw(rng, [*echo, *(arrived - echo), quiet])
                    coarse_bank = draw(rng, [coarse, arrived.sum() - coarse, quiet])
                    challenge = draw(rng, one_way[bit])
                    panels.append(dict(cells=board.cells, gain=gain, exposure=t, bit=bit,
                                       resolved=resolved, coarse=coarse_bank, challenge=challenge))
    return dict(shots=SHOTS, panels=panels)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    result = acquire()
    if args.record:
        with args.record.open("x") as dst:
            json.dump(result, dst, separators=(",", ":"))
            dst.write("\n")
    elif result != json.loads((HERE / "records.json").read_text()):
        raise ValueError("passive echo banks differ")
    print(len(result["panels"]), "panels with three full launch banks apiece")
