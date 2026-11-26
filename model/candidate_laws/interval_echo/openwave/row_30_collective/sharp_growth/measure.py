"""Reproduce both sharp-write cohorts, local calibrations and frequency predictions."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from model.experiments.material_inertia.dark_collar import label
from .. import coupling_meter, experiment
from ..confirmation import spatial
from .growth import SEED, commission, instrument

HERE = Path(__file__).resolve().parent


def meter(standard):
    mats, hz, errors = coupling_meter.local_generators()
    law = experiment.candidate.native_law
    exchange, contrast = law.local_term_parts()
    full = exchange + experiment.K * contrast
    parity = 1 - 2 * standard
    memory = label(2, parity)
    for name, neighbour in [("drawn_link", 1), ("vacant_neighbour", 0)]:
        cols = np.stack(
            [
                np.kron(memory, np.eye(law.N_SRC)[:, law.SRC.index((v, neighbour))])
                for v in (0, 1)
            ],
            axis=1,
        )
        mats[name], errors[name] = coupling_meter.restrict(full, cols)
    with (
        patch.object(
            coupling_meter, "local_generators", return_value=(mats, hz, errors)
        ),
        patch.object(coupling_meter, "SEED", SEED + 11 + standard),
    ):
        return coupling_meter.run()


def spectroscopy(folder, plan, level, bit, readings, other):
    captured = experiment.CountStream

    class RecordCounts(captured):
        def keep(self, key, clicks):
            super().keep(key, clicks)
            np.save(
                folder / (key + ".npy"),
                np.asarray(clicks, dtype="<u8"),
                allow_pickle=False,
            )

    # The frozen reader takes its effective connector coefficient from its module.
    # Its zero-contrast foil still has exchange energy one in the odd sign sector.
    physical_j = (4 * experiment.K + bit) / 4
    dual = experiment.dual

    def branch_dual(n, j):
        return dual(n, bit / 4 if j == 0 else j)

    rng = np.random.default_rng([SEED, level, bit, 3])
    with (
        patch.object(experiment, "K", physical_j),
        patch.object(experiment, "CountStream", RecordCounts),
        patch.object(experiment, "dual", branch_dual),
    ):
        obs = experiment.spectroscopy(len(plan["enrolled_route"]), rng)
    checks = spatial.score(dict(route=plan["enrolled_route"], observer=obs), folder)
    result = dict(
        spatial=checks,
        count_sha256=obs["count_sha256"],
        clock=obs["clock"],
        pair_sum_lines=obs["pair_sum_lines"],
        changed_law=obs["same_spatial_source_changed_contrast"],
        held_pair_sum_error=obs["withheld_pair_sum_max_error"],
    )
    if "even_trained_k" not in checks:
        return result
    k = np.asarray(checks["even_trained_k"])
    omega = np.asarray(checks["omega_from_pair_counts"])

    def forecast(bench, gradient):
        return np.sqrt(
            bench["predicted_m"]["value"] ** 2
            + bench["predicted_c_squared"]["value"] * gradient**2
        )

    gradient = 2 * np.sin(k / 2)
    lattice = forecast(readings, gradient)
    swapped = forecast(other, gradient)
    continuum = forecast(readings, k)
    result.update(
        lattice=lattice.tolist(),
        swapped_standard=swapped.tolist(),
        continuum=continuum.tolist(),
        max_lattice_error=float(max(abs(lattice - omega))),
        max_swapped_error=float(max(abs(swapped - omega))),
        spatially_resolved=checks["profile_pass"],
        launch_labels_rejected=checks["launch_label_rejected"],
    )
    return result


def run(directory=None):
    native = instrument()
    meters = [meter(i) for i in (0, 1)]
    shots = np.random.default_rng([SEED, 0, 8]).binomial(
        2**20, native["sharp_pointer_laws"][0][1]
    )
    result = dict(
        law=native,
        first_pointer_counts=[int(2**20 - shots), int(shots)],
        meters=meters,
        boards=[],
    )
    for level in range(3):
        for standard in (0, 1):
            setup = commission(level, standard)
            entry = dict(level=level, first_pointer=standard, setup=setup)
            if setup["complete"]:
                with TemporaryDirectory() as tmp:
                    folder = (
                        (directory / f"level{level}-pointer{standard}")
                        if directory
                        else Path(tmp)
                    )
                    folder.mkdir(parents=True, exist_ok=True)
                    wave = spectroscopy(
                        folder,
                        setup,
                        level,
                        standard,
                        meters[standard]["results"],
                        meters[1 - standard]["results"],
                    )
                entry["wave"] = wave
            result["boards"].append(entry)
    return result


def write_data(obj):
    return json.dumps(
        coupling_meter.canonical(obj),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode()


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--record", action="store_true")
    cli.add_argument("--save-counts", type=Path)
    args = cli.parse_args()
    body = write_data(run(args.save_counts))
    receipt = HERE / "records.json.gz"
    if args.record:
        with receipt.open("xb") as f:
            f.write(gzip.compress(body, mtime=0))
    elif gzip.decompress(receipt.read_bytes()) != body:
        raise SystemExit("Sharp-growth detector or preparation replay differs")
    print(hashlib.sha256(body).hexdigest(), receipt.name)
