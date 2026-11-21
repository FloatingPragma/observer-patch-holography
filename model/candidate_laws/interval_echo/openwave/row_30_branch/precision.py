"""Full-operator residual check for the two small first shifts on the 80-cell fixture."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import shortest_path

from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from .experiment import HERE, portable
from .refinement import analyze_radius, layout


def run():
    sector = PairSector(federation(1))
    graph, groups, _ = layout(sector.car)
    distances = shortest_path(sector.car.adjacency(), directed=False, unweighted=True)[sector.aa,sector.bb]
    data = analyze_radius(sector,graph,groups,sector.hamiltonian(candidate.KAPPA),distances,10)
    # The detailed primary record contains every line weight. At the tighter radius keep the
    # basis-independent moments, which are the quantities used for the comparison.
    for panel in data["panels"].values():
        for row in panel["classes"]:
            row.pop("low_band_weights")
        panel["residual_pass"] = bool(panel["full_residual_over_first_shift"] < .01)
        panel["relative_shifts_monotone"] = bool(np.all(np.diff([r["shift_from_rest"] for r in panel["classes"]]) >= 0))
    return portable(data)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, help="write a separate receipt to a new path")
    path = cli.parse_args().output
    record = run()
    if path:
        with path.open("x") as f:
            json.dump(record,f,indent=2,allow_nan=False)
            f.write("\n")
    elif record != json.loads((HERE/"precision.json").read_text()):
        raise SystemExit("precision refinement differs from the archive")
    print(record["radius"],record["maximum_full_residual"],
          [(name,p["full_residual_over_first_shift"]) for name,p in record["panels"].items()])
