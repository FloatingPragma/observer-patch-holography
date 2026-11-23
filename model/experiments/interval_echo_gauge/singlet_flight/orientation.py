"""Transport a reference-tagged orbital qubit with colour filters at the two endpoints."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.interval_echo_angular.charged_orientation import (
    measure as orbital,
)
from model.experiments.interval_echo_gauge.singlet_flight.measure import (
    HERE,
    colour_meter,
)


def run():
    colour, _, _, _ = colour_meter()
    singlet = np.eye(12) - colour
    rows = []
    for level in range(3):
        board = federation(level)
        home, receivers, owners = orbital.design(board)
        launch = orbital.code_columns(owners[home][1], board.n)
        singlet_launch = launch.copy()
        source = slice(12 * home, 12 * (home + 1))
        singlet_launch[source] = singlet @ launch[source]
        starts = np.column_stack((launch, singlet_launch))
        for sign in (-1, 1):
            h = orbital.generator(board, sign, candidate.KAPPA)
            terminal = orbital.step(h, starts, 32)
            for distance, owner in enumerate(receivers, 1):
                patch = slice(12 * owner, 12 * (owner + 1))
                detector = orbital.code_columns(owners[owner][1], board.n)[patch]
                for start, name in enumerate(("direct", "singlet")):
                    basis = slice(2 * start, 2 * (start + 1))
                    local = terminal[patch, basis]
                    selected = (
                        ("all", local),
                        ("singlet", singlet @ local),
                        ("triplet", colour @ local),
                    )
                    for arm, (filter_name, filtered) in enumerate(selected):
                        rng = np.random.default_rng(
                            [852221, level, sign + 1, distance, start, arm]
                        )
                        reading = orbital.tomography(filtered, detector, rng)
                        rows.append(
                            dict(
                                cells=board.cells,
                                sign=sign,
                                from_cell=home,
                                to_cell=owner,
                                steps=distance,
                                launch=name,
                                receive=filter_name,
                                launch_acceptance=float(
                                    np.linalg.norm(starts[:, basis]) ** 2 / 2
                                ),
                                **reading,
                            )
                        )
    return dict(
        protocol_sha256=hashlib.sha256(
            (HERE / "ORIENTATION_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        exposure=32,
        shots_per_setting=orbital.SHOTS,
        negativity_threshold=5 * orbital.EPS,
        outcomes=[
            "code_plus",
            "code_minus",
            "selected_other_orbital",
            "miss_or_filter_rejection",
        ],
        reads=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    with parser.parse_args().output.open("x") as file:
        json.dump(run(), file, allow_nan=False, indent=2)
        file.write("\n")
