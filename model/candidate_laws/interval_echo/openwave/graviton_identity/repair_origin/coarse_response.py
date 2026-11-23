"""Read orbital quadrupoles after a low-energy collective launch on an intact board."""

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.sparse import coo_matrix, eye, kron
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.graviton_identity.repair_origin.measure import angular_frame
from model.executable_model import federation
from model.experiments.interval_echo_space.shell_source import prepare

EXPOSURES = (0, 1, 2, 3, 4, 8)
SHOTS = 65536
SEED = 869112


def compile_wave(carrier, target):
    apparatus = SimpleNamespace(car=carrier, active=carrier.seams)
    return prepare(apparatus, target)[0]


def measure_orbital(state, analyzers, rng):
    # Five separate yes/no instruments, each operated in parallel on disjoint cells.
    amplitudes = state.reshape(-1, 12) @ analyzers.conj()
    probabilities = np.abs(amplitudes) ** 2
    records = []
    for angular_prob in probabilities.T:
        quiet = max(0., 1 - angular_prob.sum())
        distribution = np.append(angular_prob, quiet)
        records.append(rng.multinomial(SHOTS, distribution / distribution.sum()).tolist())
    return {
        "quadrupole_probability": float(probabilities.sum()),
        "counts_by_mode_then_cell_or_complement": records,
        "pooled_reading": sum(SHOTS - row[-1] for row in records) / SHOTS,
    }


def board_run(level, cell_kinetic, analyzers, rng):
    board = federation(level)
    ground = 6 * candidate.KAPPA
    kinetic = candidate.generator(board) - ground * eye(board.n, format="csr")
    inside = kron(eye(board.cells), cell_kinetic, format="csr")
    outside = kinetic - inside
    # Form the independent exterior sum |a-b><a-b|/2 and check every entry.
    a, b = np.array(board.external).T
    incidence = coo_matrix((np.tile([1., -1.], len(a)),
                            (np.repeat(np.arange(len(a)), 2),
                             np.column_stack([a, b]).ravel())),
                           shape=(len(a), board.n)).tocsr()
    residual = outside - incidence.T @ incidence / 2
    np.testing.assert_allclose(residual.data, 0., atol=1e-13)

    # Descendants have consecutive cell indices under this carrier's four-child rule.
    children = 4 ** level
    face = np.arange(board.n) < 12 * children
    boundary = face[a] != face[b]
    assert boundary.sum() == 3 * 2 ** level
    launch = compile_wave(board, face / np.sqrt(12 * children))
    boundary_prob = abs(launch[a[boundary]] - launch[b[boundary]]) ** 2 / 2
    energy = float(boundary_prob.sum())
    np.testing.assert_allclose(energy, 1 / (8 * 2 ** level), atol=1e-13)
    np.testing.assert_allclose(np.vdot(launch, kinetic @ launch), energy, atol=1e-13)
    boundary_records = rng.multinomial(SHOTS, np.append(boundary_prob, 1-energy))
    bound = energy / 3
    runs = {}
    for name, generator in (("registered", kinetic), ("slow_internal", outside + inside / 16)):
        states = expm_multiply(-1j * generator, launch, start=0, stop=8, num=9)
        readings = [dict(exposure=t, **measure_orbital(states[t], analyzers, rng))
                    for t in EXPOSURES]
        peak = max(row["quadrupole_probability"] for row in readings)
        if name == "registered":
            assert peak <= bound + 1e-13
        else:
            assert peak > bound
        runs[name] = readings
    uniform = compile_wave(board, np.full(board.n, 1 / np.sqrt(board.n)))
    np.testing.assert_allclose(kinetic @ uniform, 0., atol=1e-13)
    local_orbital = np.zeros(board.n, complex)
    local_orbital[:12] = analyzers[:, 0]
    checks = {"whole_screen_uniform": measure_orbital(uniform, analyzers, rng),
              "local_quadrupole": measure_orbital(local_orbital, analyzers, rng)}
    np.testing.assert_allclose([r["quadrupole_probability"] for r in checks.values()], [0., 1.], atol=1e-13)
    return {"cells": board.cells, "descendants": children,
            "boundary_seams": np.column_stack([a[boundary], b[boundary]]).tolist(),
            "energy_calibration_counts_then_complement": boundary_records.tolist(),
            "sampled_excess": (SHOTS-int(boundary_records[-1]))/SHOTS,
            "excess_energy": energy, "all_stops_orbital_bound": bound,
            "readings": runs, "controls_at_release": checks}


def acquire():
    cell, frame, _, _ = angular_frame()
    kinetic = candidate.generator(cell) - 6*candidate.KAPPA*eye(cell.n, format="csr")
    analyzers = np.column_stack([compile_wave(cell, column) for column in frame.T])
    projector = analyzers @ analyzers.conj().T
    np.testing.assert_allclose(kinetic @ projector, 3*projector, atol=1e-13)
    assert np.linalg.eigvalsh(kinetic.toarray() - 3*projector).min() > -1e-13
    rng = np.random.default_rng(SEED)
    return {"seed": SEED, "trials_per_instrument": SHOTS, "orbital_gap": 3,
            "boards": [board_run(level, kinetic, analyzers, rng) for level in range(3)]}


def check_replay(saved, fresh):
    """Require exact integer records and allow roundoff only in the numerical audit."""
    if isinstance(fresh, dict):
        assert saved.keys() == fresh.keys()
        for key in fresh:
            check_replay(saved[key], fresh[key])
    elif isinstance(fresh, list):
        for first, second in zip(saved, fresh, strict=True):
            check_replay(first, second)
    elif isinstance(fresh, float):
        np.testing.assert_allclose(saved, fresh, rtol=0, atol=1e-13)
    else:
        assert saved == fresh


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write a new receipt, otherwise check the archive")
    output = parser.parse_args().output
    receipt = acquire()
    if output:
        with output.open("x") as stream:
            json.dump(receipt, stream, separators=(",", ":"), allow_nan=False)
            stream.write("\n")
    else:
        archive = json.loads(Path(__file__).with_name("coarse_response.json").read_text())
        check_replay(archive, receipt)
    for board in receipt["boards"]:
        print(board["cells"], board["all_stops_orbital_bound"],
              {name: max(row["pooled_reading"] for row in rows)
               for name, rows in board["readings"].items()})


if __name__ == "__main__":
    main()
