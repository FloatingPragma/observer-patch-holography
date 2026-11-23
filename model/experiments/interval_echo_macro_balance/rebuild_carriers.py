"""Generate both holdout files with the public OPH geometry builder and authenticate them."""

import argparse
from pathlib import Path

import numpy as np

from model.executable_model.refinement import load_refinement


def build(output):
    from oph_defects import fedgeom

    output.mkdir(parents=True, exist_ok=True)
    for level in (3, 4):
        geometry = fedgeom.build(level, cache=False)
        path = output / f"FEDERATION_LEVEL{level}.npz"
        np.savez_compressed(
            path,
            external=np.column_stack((geometry.glue_a, geometry.glue_b)).astype(np.int32),
            cell_centers=geometry.carrier_xyz.astype(np.float64),
            level=np.array(level),
            cells=np.array(geometry.carriers),
        )
        verified = load_refinement(path, level=level)
        print(f"{verified.cells} cells, SHA-256 {verified.sha256}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    build(parser.parse_args().output)
