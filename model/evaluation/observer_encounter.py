"""Two local pair readings after a separated preparation on the 20-cell sphere.

Run from the RER root with python -m model.evaluation.observer_encounter.
The first patch tests one seam. A positive reading prepares a sharp pair there.
The same patch later reads all its internal seams and the complementary outcome.
"""
import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1 import candidate as contrast
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate as echo
from model.executable_model import federation


RECEIPT = Path(__file__).with_suffix(".json")


def measure():
    if echo.native_law is not contrast.law or echo.KAPPA != contrast.KAPPA:
        raise ValueError("The two candidate free sectors must be audited separately.")
    car = federation(0)
    sector = PairSector(car)
    internal = [seam for seam in car.seams if max(seam) < 12]
    target = min(internal)
    starts = [sector.index(a, b) for a, b in combinations(range(12), 2)
              if not sector.adjacent[sector.index(a, b)]]
    target_index = sector.index(*target)
    internal_indices = [sector.index(*edge) for edge in internal]
    seed = np.zeros(sector.P)
    seed[target_index] = 1
    result = {"federation_cells": car.cells, "sector_dimension": sector.P,
              "observer_cells": [0], "first_seam": list(map(int, target)),
              "separated_start_count": len(starts), "local_seams": len(internal),
              "first_wait": 1, "repeat_wait": 4,
              "preparation": "Uniform mixture over distinct nonadjacent ports of the observer cell, with one load on each and no other load.",
              "instrument": "At the first wait one specified seam reports a pair or quiet. After a positive, the same patch reads each of its internal seams or quiet at the second wait. State amplitudes evolve on the whole federation.",
              "omitted_services": "Free native binary sector. The interval filter, gate sequencer and readout duration are idle. The two ideal projective patch reads and their times are supplied.",
              "initial_positive": float(np.mean(np.abs(seed[starts]) ** 2)),
              "laws": {}}
    for name, kappa in (("contrast", contrast.KAPPA), ("exchange_only", 0)):
        h = sector.hamiltonian(kappa)
        if (h - h.T).nnz or np.iscomplexobj(h.data):
            raise ValueError("The reciprocity reduction requires a real symmetric generator.")
        # For real symmetric H, U[target, start] = U[start, target] for every start.
        first = np.abs(expm_multiply(-1j * h, seed)) ** 2
        later = np.abs(expm_multiply(-4j * h, seed)) ** 2
        if not np.allclose([first.sum(), later.sum()], 1, atol=1e-11, rtol=0):
            raise ValueError("Incomplete two-load sector propagation")
        p = first[starts].mean()
        repeat = later[internal_indices].sum()
        values = {"first_positive_per_trial": p, "first_quiet_per_trial": 1-p,
                  "repeat_positive_given_first": repeat, "repeat_quiet_given_first": 1-repeat,
                  "same_seam_given_first": later[target_index],
                  "other_local_seam_given_first": repeat - later[target_index],
                  "both_positive_per_trial": p * repeat}
        result["laws"][name] = {"kappa": kappa, **{key: round(float(val), 12)
                                                  for key, val in values.items()}}
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path, help="write the result to a new file")
    args = parser.parse_args()
    output = measure()
    if args.record:
        with args.record.open("x") as f:
            f.write(json.dumps(output, indent=2) + "\n")
    elif output != json.loads(RECEIPT.read_text()):
        raise SystemExit("Separated-encounter receipt differs from fresh measurement")
    print(json.dumps(output["laws"], indent=2))
