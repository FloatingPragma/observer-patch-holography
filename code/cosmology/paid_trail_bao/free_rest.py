"""Hold galaxies out of the calibration of the measured isolated-rest field."""

import argparse
from hashlib import sha256
import json

import numpy as np

from . import kernel_budget as budget
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    isolated_medium as lab,
)


def shell_field(reads, tape):
    board = lab.lab.cold.federation(0)
    marks, _ = lab.lab.cold.read_map(board, tape["spectral_tape"])
    shells = lab.lab.cold.distances(board, marks[0])
    field = np.zeros(board.n)
    for arm, mask, sign in (("blank", [], -1), ("written", marks, 1)):
        p = np.array(reads[arm]["sites"]) / lab.READS
        adj = board.adjacency().tolil()
        for i, j in mask:
            adj[i, j] = adj[j, i] = 0
        field += sign * (p + adj.tocsr() @ p)
    return np.bincount(shells, weights=field) * np.mean(lab.ordinary_gain(tape))


def run():
    origin = lab.HERE / "isolated_medium.json"
    previous = lab.HERE / "mediated_release.json"
    reads = json.loads(origin.read_bytes())["observations"]["rest"]["65536"]
    tape = json.loads(previous.read_bytes())
    shape = shell_field(reads, tape)
    inputs, charts, meta, _ = budget.load_inputs()
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                origin,
                previous,
                budget.HERE / "FREE_REST.md",
                budget.HERE / "source_profile.json",
                budget.local.profiles.maps.HERE / "leroy_table7.dat",
                budget.local.profiles.maps.HERE / "leroy_sample.dat",
                budget.local.profiles.galaxies.sparc.DATA / "table1.dat",
            )
        },
        shell_field=shape.tolist(),
        held=budget.local.evaluate(
            inputs, charts, meta, (np.arange(len(shape)), shape)
        ),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    file = budget.HERE / "free_rest.json"
    if args.write:
        file.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(file.read_bytes()), out)
    print({k: v for k, v in out["held"].items() if k != "galaxies"})
